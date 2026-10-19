"""Average-information and delta-method inference for trait covariances.

Called by the opt-in REML routes, with dense or exact structured solves.
No inverse ridge is used.
Boundary and nonidentified parameterizations are explicitly unavailable.
"""
import numpy as np
import pandas as pd
import warnings


def average_information(theta, unpack, y, X, solve, covariance_matmul, jac_matrices=None):
    """Exact AI in the existing parameter coordinates, using first derivatives only.

    Differentiate the small trait covariances, then apply each V_i to Py.
    Neither an observation-sized covariance Jacobian nor a Hessian is formed.
    """
    import torch
    th = theta.detach()
    if jac_matrices is None:
        jac_matrices = torch.func.jacfwd(lambda t: torch.stack(unpack(t)))(th)
    with torch.no_grad():
        Vy = solve(y[:, None])
        VX = solve(X)
        beta_cov = torch.linalg.inv(X.T @ VX)
        Py = Vy - VX @ beta_cov @ (X.T @ Vy)
        F = torch.cat([covariance_matmul(list(jac_matrices[..., q]), Py)
                       for q in range(th.numel())], dim=1)
        PF = solve(F) - VX @ beta_cov @ (VX.T @ F)
        info = .5 * F.T @ PF
        return .5 * (info + info.T)


def covariance_inference(theta, unpack, covariance_from_matrices, y, X, names, converged,
                         *, solve=None, covariance_matmul=None):
    import torch
    from varcomp_asreml import active_set_inverse
    th = theta.detach()
    matrices = unpack(th)
    T = matrices[0].shape[0]
    tri = torch.tril_indices(T, T, device=th.device)
    def entries(t):
        return torch.cat([m[tri[0], tri[1]] for m in unpack(t)])
    values = entries(th).detach().cpu().numpy()
    labels = [(name, int(i), int(j)) for name in names for i,j in tri.T.tolist()]
    status = 'available'
    # A Cholesky covariance at a singular boundary does not admit ordinary
    # interior Wald inference. A small transformed gradient is not sufficient.
    for matrix in matrices:
        diagonal = matrix.diag().clamp(min=torch.finfo(matrix.dtype).tiny).sqrt()
        corr = matrix / (diagonal[:, None]*diagonal[None, :])
        if float(torch.linalg.eigvalsh(corr).min()) <= 1e-7:
            status = 'unavailable_covariance_boundary'
    if not converged:
        status = 'unavailable_nonconverged'
    C_entries = np.full((len(values),len(values)), np.nan)
    info = None
    if status == 'available':
        # Small covariance Jacobian; do not form an n*n*p autograd tensor.
        jac_matrices = torch.func.jacfwd(lambda t: torch.stack(unpack(t)))(th)
        J = torch.func.jacfwd(entries)(th).detach().cpu().numpy()
        with torch.no_grad():
            if (solve is None) != (covariance_matmul is None):
                raise ValueError('Provide both solve and covariance_matmul, or neither')
            if solve is None:
                V = covariance_from_matrices(matrices)
                L = torch.linalg.cholesky(V)
                solve = lambda rhs: torch.cholesky_solve(rhs, L)
                zero = covariance_from_matrices([torch.zeros_like(m) for m in matrices])
                covariance_matmul = lambda ms, rhs: (covariance_from_matrices(ms)-zero) @ rhs
            info = average_information(th, unpack, y, X, solve, covariance_matmul,
                                       jac_matrices=jac_matrices).cpu().numpy()
        inv = active_set_inverse(info, ['U']*len(th))
        if not np.isfinite(inv).all():
            status = 'unavailable_singular_information'
        else:
            C_entries = J@inv@J.T
            C_entries = .5*(C_entries+C_entries.T)
    if status != 'available':
        warnings.warn('Trait covariance inference ' + status, RuntimeWarning)
    return dict(labels=labels, estimates=values, covariance=C_entries,
        information=info, inference_available=status=='available', inference_status=status,
        method='average_information_delta_method', parameter_scale='standardized',
        covariance_parameter_conditioning='kernel_weights_and_standardization')


def response_inference(payload, trait_names, scales, environment_names=None):
    """Transform all covariance entries and their joint covariance, then r."""
    if payload is None:
        return None
    labels = payload['labels']; estimates = np.asarray(payload['estimates'])
    C = np.asarray(payload['covariance'])
    sd = np.atleast_2d(np.asarray(scales,dtype=float))
    envs = list(environment_names) if environment_names is not None else ['overall']
    if len(envs) != len(sd): raise ValueError('environment scale dimensions do not agree')
    rows, corr_rows = [], []
    index = {tuple(key):i for i,key in enumerate(labels)}
    for env, scale in zip(envs,sd):
        for q,(term,i,j) in enumerate(labels):
            factor = scale[i]*scale[j]
            variance = C[q,q]*factor**2
            se = np.sqrt(variance) if np.isfinite(variance) and variance >= 0 else np.nan
            rows.append(dict(component=f'{term}!{trait_names[i]}:{trait_names[j]}',
                environment=env, trait1=trait_names[i], trait2=trait_names[j],
                estimate=estimates[q]*factor, **{'std.error':se},
                inference_status=payload['inference_status']))
            if i==j: continue
            a=index[(term,i,i)]; b=index[(term,j,j)]
            denom=np.sqrt(estimates[a]*estimates[b])
            r=estimates[q]/denom
            gradient=np.zeros(len(estimates)); gradient[q]=1/denom
            gradient[a]=-r/(2*estimates[a]); gradient[b]=-r/(2*estimates[b])
            rv=float(gradient@C@gradient)
            corr_rows.append(dict(component=f'{term}!cor({trait_names[i]},{trait_names[j]})',
                environment=env, estimate=r, **{'std.error':np.sqrt(rv) if np.isfinite(rv) and rv>=0 else np.nan},
                inference_status=payload['inference_status']))
    return dict(varcomp=pd.DataFrame(rows), correlation_parameters=pd.DataFrame(corr_rows),
        summary={k:payload[k] for k in ['inference_available','inference_status','method','covariance_parameter_conditioning']})
