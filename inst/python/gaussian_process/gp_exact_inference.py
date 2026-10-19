"""REML inference for the covariance actually used by the exact GP.

Zero input weights disable terms in the legacy API. Learned near-zero weights
do not disable them: the structural masks are captured before optimization.
Supplied observation noise stays fixed; GaussianLikelihood has one residual
parameter, irrespective of the number of environment labels.
"""
import numpy as np
import os
import warnings
from dataclasses import replace


def build_exact_inference(model, likelihood, train_x, train_y, *,
        geno_kernel_names, env_labels, env_col_name='Env',
        kernel_term_template='vm(GID,{K})', env_scales=None, nedf=0,
        max_iter=12, tol_loglik=1e-3, tol_theta=1e-4, jitter=1e-6,
        refine_kernel_scales=True, theta0_warm=None, max_iter_warm=None, **unused):
    import torch
    import gpytorch
    from ai_reml import (_exact_dense_precompute, _exact_dense_reml_loglik_from_pre,
        _exact_dense_score_and_AI_from_pre, _ai_reml_optimize)
    from varcomp_asreml import (ParamSpec, DenseVCBackend, build_asreml_varcomp_table,
        build_summary, compute_genetic_residual_covariance, inference_status)

    kernel = model.covar_module
    X = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
    gi, ei = train_x[:, 0].long(), train_x[:, 1].long()
    with torch.no_grad():
        total = kernel(train_x, train_x).to_dense().detach()
    specs, values, derivatives, setters = [], [], [], []
    fixed_covariance = total.clone()
    learned = bool(refine_kernel_scales and kernel.learn_scales)
    absorbed = []
    # Thin QR avoids constructing an n by (n-p) residual-space basis.
    Q, _ = torch.linalg.qr(X, mode='reduced')

    def append_component(matrix, value, attr, index, spec):
        nonlocal fixed_covariance
        fixed_covariance = fixed_covariance - float(value) * matrix
        residual = matrix - Q @ (Q.mT @ matrix)
        residual = residual - (residual @ Q) @ Q.mT
        scale = float(torch.linalg.norm(matrix))
        is_absorbed = scale == 0 or float(torch.linalg.norm(residual)) <= 1e-10 * scale
        if is_absorbed:
            value = 0.
            spec = replace(spec, fixed=True)
            # Masking also removes this term from predictions and row priors.
            mask = getattr(kernel, attr + '_fix')
            with torch.no_grad():
                if index is None: mask.zero_()
                else: mask[index] = 0.
            absorbed.append(spec.component_name())
            warnings.warn(f'{spec.component_name()} absorbed by fixed effects; fixed at zero. '
                          'Set its kernel weight to zero (w_e = 0 for Env), or drop the '
                          'corresponding column (Env) from fixed_effects.', RuntimeWarning)
        specs.append(spec); values.append(float(value)); derivatives.append(matrix)
        setters.append((attr, index))

    scales = kernel.get_effective_scales()
    wg, we = scales[0], scales[-1]
    for i, name in enumerate(geno_kernel_names):
        if bool(kernel.w_g_fix[i] != 0):
            append_component(kernel.G_list[i][gi[:, None], gi[None, :]], wg[i], 'w_g', i,
                ParamSpec('var_kernel', kernel_term_template.format(K=name), sub_label='var', lower=0., fixed=not learned))
    if hasattr(kernel, 'w_ge_fix') and kernel.S_e is not None:
        wge = scales[1]
        S = kernel.S_e[ei[:, None], ei[None, :]]
        for i, name in enumerate(geno_kernel_names):
            if bool(kernel.w_ge_fix[i] != 0):
                append_component(kernel.G_list[i][gi[:, None], gi[None, :]]*S, wge[i], 'w_ge', i,
                    ParamSpec('var_kernel', env_col_name+':'+kernel_term_template.format(K=name), sub_label='var', lower=0., fixed=not learned))
    if bool(kernel.w_e_fix != 0):
        S = ((ei[:, None] == ei[None, :]).to(train_y.dtype) if kernel.S_e is None
             else kernel.S_e[ei[:, None], ei[None, :]])
        append_component(S, we, 'w_e', None,
            ParamSpec('var_random', env_col_name, sub_label='var', lower=0., fixed=not learned))

    fixed_noise = isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood)
    if fixed_noise:
        fixed_covariance += torch.diag(likelihood.noise.detach())
    else:
        noise_floor = float(likelihood.noise_covar.raw_noise_constraint.lower_bound)+1e-12
        specs.append(ParamSpec('resid_env', env_col_name, env_label=str(env_labels[0]) if len(env_labels)==1 else 'common', lower=noise_floor))
        values.append(float(likelihood.noise.detach().item()))
        derivatives.append(torch.eye(len(train_y), device=train_y.device, dtype=train_y.dtype))
        setters.append(('noise', None))
    if not specs:
        raise ValueError('No covariance parameters are available for REML inference')

    cache = {}
    def evaluate(theta):
        key = theta.tobytes()
        if cache.get('key') != key:
            V = fixed_covariance.clone()
            for value, derivative in zip(theta, derivatives):
                V += float(value)*derivative
            with torch.no_grad():
                pre = _exact_dense_precompute(train_y, X, V, jitter=jitter)
                score, ai = _exact_dense_score_and_AI_from_pre(pre,
                    [lambda v, D=D: D@v for D in derivatives], derivatives)
            cache.update(key=key, ll=_exact_dense_reml_loglik_from_pre(pre), score=score, ai=ai)
        return cache
    start = np.array(values)
    if theta0_warm is not None and np.asarray(theta0_warm).shape == start.shape:
        start = np.asarray(theta0_warm, dtype=float).copy()
        # Fixed terms must never be overwritten by an earlier fit's warm start.
        for i, spec in enumerate(specs):
            if spec.fixed: start[i] = values[i]
    backend = DenseVCBackend(specs, lambda t: evaluate(t)['ll'],
        lambda t: evaluate(t)['score'], lambda t: evaluate(t)['ai'])
    def override(name, default, convert):
        value = os.environ.get(name, '').strip()
        if not value: return default
        try:
            parsed = convert(value)
            if not np.isfinite(parsed) or parsed <= 0: raise ValueError()
            return parsed
        except ValueError:
            warnings.warn(f'Ignoring invalid {name}={value!r}; use a positive number.', RuntimeWarning)
            return default
    max_iter = override('PREDICTPRO_GP_MAX_ITER', max_iter, int)
    tol_loglik = override('PREDICTPRO_GP_TOL_LOGLIK', tol_loglik, float)
    tol_theta = override('PREDICTPRO_GP_TOL_THETA', tol_theta, float)
    result = _ai_reml_optimize(backend, start, max_iter=int(max_iter),
        tol_loglik=float(tol_loglik), tol_theta=float(tol_theta), damping=1., max_halvings=8)
    with torch.no_grad():
        for value, (attr, index), spec in zip(result.theta, setters, specs):
            if spec.fixed: continue
            if attr == 'noise':
                likelihood.noise = max(float(value), float(likelihood.noise_covar.raw_noise_constraint.lower_bound)+1e-12)
            else:
                raw = getattr(kernel, attr+'_unconstrained')
                v = torch.as_tensor(max(float(value), 1e-12), device=raw.device, dtype=raw.dtype)
                unconstrained = v if v > 20 else torch.log(torch.expm1(v))
                if index is None: raw.copy_(unconstrained.reshape(raw.shape))
                else: raw[index].copy_(unconstrained)
    table = build_asreml_varcomp_table(specs, result.theta, result.ai, result.bounds,
        result.pct_change, env_scales)
    table['identifiability_reason'] = [
        'absorbed by fixed effects' if s.component_name() in absorbed else '' for s in specs]
    table.loc[table['identifiability_reason'] != '', ['std.error', 'z.ratio']] = np.nan
    summary = build_summary(loglik=result.loglik, nedf=nedf,
        parameters_free=sum(b not in ('B','F') for b in result.bounds),
        n_ai_iter=result.n_iter, converged=result.converged)
    status = inference_status(result.ai, result.bounds, result.converged)
    summary.update(status)
    if result.converged:
        aliased = (table['identifiability_reason'] == '') & ~table['bound'].isin(['B', 'F']) & table['std.error'].isna()
        table.loc[aliased, 'identifiability_reason'] = 'not individually identifiable from active information'
    if not status['inference_available']:
        if not result.converged:
            table[['std.error', 'z.ratio']] = np.nan
        warnings.warn('GP covariance inference ' + status['inference_status'], RuntimeWarning)
    if result.converged:
        cov = compute_genetic_residual_covariance(specs, result.theta, result.ai, result.bounds, env_scales)
        if cov is not None: summary['genetic_residual_covariance'] = cov
    summary.update(theta=result.theta.tolist(), theta_layout='gp_exact_matched_v1',
        theta_shape={'n_params':len(specs)},
        theta_components=[s.component_name() for s in specs],
        residual_structure='fixed_observation' if fixed_noise else 'homogeneous',
        prediction_covariance_source='reported_reml_fit')
    summary['absorbed_components'] = absorbed
    return table, summary, np.asarray(result.ai)
