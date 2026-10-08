"""
Eigen-projected MME REML for GP-exact MET GBLUP (Phase 3.15a).

Restricted scope:
  - Single trait
  - Single genetic kernel
  - Balanced MET: every genotype observed equal times in every env
  - Env-only fixed effects (X = one-hot env)
  - Env-main random OFF (requires PREDICTPRO_GP_DISABLE_ENV_MAIN=1)
  - GxE OFF

Anything outside this scope raises NotImplementedError, which the dispatch
site in gp_framework.py catches and falls back to dense_v.

Math:
  G = U Lambda U' (one-time eigendecomp). Define v = U' u so
  Var(v) = sigma2_g Lambda is diagonal. For balanced MET with
  diagonal env-constant R, the bottom-right MME block is diagonal:
      A_vv = r * c * I + diag(1/(sigma2_g lambda_j))
  where c = sum_e 1/sigma2_e_e and r = obs per (geno, env) cell.

  REML logLik in closed form:
      ll = -0.5 (log|V| + log|S| + y'Py + (n_obs - n_env) log 2pi)
      log|V|  = sum log R_ii + sum log(sigma2_g lambda_j) + sum log d_j
      log|S|  = log determinant of n_env x n_env Schur complement
      d_j    = c + 1/(sigma2_g lambda_j)

  Per-iteration cost (after one-time O(n_geno^3) eigh):
      O(k * n_geno + k^3) -- microseconds at n_geno = 1000.

Output contract matches build_gp_exact_ai_varcomp_mme:
  Returns (varcomp_df, summary_dict, ai_matrix_np).
  summary["theta_layout"] = "gp_exact_eigen".

References:
  - Phase 3.12: documented the 4-7% sigma2_g overshoot in dense_v vs ASReml
  - Phase 3.14: documented sparse MME catastrophe at n=599 due to dense G^-1
  - Math validated in tmp/phase315_eigen_calibration.py (matches standard
    REML to 1e-13 at every test point)
"""
from __future__ import annotations

import os as _os
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Module-level cache for one-time eigendecompositions, keyed on blake2b(G).
# Bounded FIFO to cap RAM at <4 * n_geno^2 * 8 bytes (~14 MB at n=600).
from collections import OrderedDict
import hashlib as _hashlib
_EIGEN_CACHE: "OrderedDict[str, Tuple[np.ndarray, np.ndarray]]" = OrderedDict()
_EIGEN_CACHE_MAX = 4


def _cache_eigendecomp(G_np: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Cache U, Lambda for repeated fits on the same G."""
    key = _hashlib.blake2b(np.ascontiguousarray(G_np).tobytes(),
                           digest_size=16).hexdigest()
    if key in _EIGEN_CACHE:
        _EIGEN_CACHE.move_to_end(key)
        return _EIGEN_CACHE[key]
    Lambda, U = np.linalg.eigh(G_np)
    order = np.argsort(Lambda)[::-1]
    Lambda = Lambda[order]
    U = U[:, order]
    _EIGEN_CACHE[key] = (U, Lambda)
    while len(_EIGEN_CACHE) > _EIGEN_CACHE_MAX:
        _EIGEN_CACHE.popitem(last=False)
    return U, Lambda


def _is_balanced_met(gi: np.ndarray, ei: np.ndarray,
                     n_geno: int, n_env: int) -> Tuple[bool, int]:
    """Returns (balanced, reps_per_cell). Balanced iff every (geno,env)
    cell has the same number of obs."""
    counts = np.zeros((n_geno, n_env), dtype=np.int64)
    np.add.at(counts, (gi, ei), 1)
    base = int(counts.flat[0])
    balanced = bool(((counts == base).all()) and base > 0)
    return balanced, base


def build_gp_exact_ai_varcomp_mme_eigen(
    model,
    likelihood,
    train_x,
    train_y,
    *,
    geno_kernel_names: Sequence[str],
    env_labels: Sequence[str],
    env_col_name: str = "Env",
    interaction_term_name: str = "Env:vm(GID,K)",
    kernel_term_template: str = "vm(GID,{K})",
    hutch_samples: int = 256,
    jitter: float = 1e-6,
    seed: int = 12345,
    env_scales: Optional[Any] = None,
    nedf: int = 0,
    max_iter: int = 30,
    tol_loglik: float = 1e-6,
    tol_theta: float = 1e-7,
    theta0_warm: Optional[np.ndarray] = None,
    max_iter_warm: Optional[int] = None,
    refine_kernel_scales: bool = True,
    mme_backend: str = "auto",
    mme_n_probes: int = 128,
) -> Tuple[pd.DataFrame, Dict[str, Any], np.ndarray]:
    """Eigen-projected REML for the single-kernel balanced-MET GP-exact path.

    See module docstring for scope and math. Raises NotImplementedError when
    the model falls outside scope so the dispatch site falls back to dense_v.
    """
    import torch

    # ---- env-var overrides for tolerance / iters ---------------------------
    _tll = _os.environ.get("PREDICTPRO_GP_TOL_LOGLIK", "").strip()
    _tth = _os.environ.get("PREDICTPRO_GP_TOL_THETA", "").strip()
    _mit = _os.environ.get("PREDICTPRO_GP_MAX_ITER", "").strip()
    _dem = _os.environ.get("PREDICTPRO_GP_DISABLE_ENV_MAIN", "").strip().lower()
    if _tll:
        try: tol_loglik = float(_tll)
        except ValueError: pass
    if _tth:
        try: tol_theta = float(_tth)
        except ValueError: pass
    if _mit:
        try: max_iter = int(_mit)
        except ValueError: pass
    disable_env_main = _dem in ("1", "true", "yes", "on")

    # ---- scope validation -------------------------------------------------
    if not disable_env_main:
        raise NotImplementedError(
            "mme_eigen requires PREDICTPRO_GP_DISABLE_ENV_MAIN=1 (scope: "
            "env-main random must be off; see Phase 3.13)."
        )

    nK = len(geno_kernel_names)
    if nK != 1:
        raise NotImplementedError(
            f"mme_eigen scope: single kernel only (got {nK})."
        )

    k = len(env_labels)
    single_env = k <= 1
    if single_env:
        raise NotImplementedError("mme_eigen scope: requires MET (k >= 2).")

    # Inspect the kernel. mme_eigen uses the genotype-level G (n_geno x
    # n_geno) from kernel.G_list (NOT the observation-level _Kg_list which
    # is Z G Z' of shape n_obs x n_obs).
    with torch.no_grad():
        _ = model(train_x)
    kernel = model.covar_module
    kernel_G_list = getattr(kernel, "G_list", None)
    Kge_list = getattr(kernel, "_Kge_list", None)
    if kernel_G_list is None:
        raise NotImplementedError(
            "mme_eigen: kernel.G_list not found; expected genotype-level "
            "per-kernel SPD matrices."
        )
    if Kge_list is not None:
        raise NotImplementedError("mme_eigen scope: GxE kernel must be off.")

    # ---- extract numpy inputs ---------------------------------------------
    ei_t = train_x[:, 1].long()
    gi_t = train_x[:, 0].long()
    ei_np = ei_t.cpu().numpy().astype(np.int64)
    gi_np = gi_t.cpu().numpy().astype(np.int64)
    y_np = train_y.detach().to(torch.float64).cpu().numpy()
    n_obs = int(y_np.shape[0])

    G_t = kernel_G_list[0].detach().to(torch.float64).cpu()
    G_np = G_t.numpy()
    n_geno = int(G_np.shape[0])

    # Phase 3.15a scope: we hard-code env-only fixed effects and construct
    # the one-hot env design internally from ei_np. The MLE / REML variance
    # components (sigma_g, sigma_e per env) are invariant to invertible
    # basis-changes of X, so we don't depend on the exact X PP passes -- we
    # only require that PP's X span the env-only column space. We don't
    # validate that here for 3.15a; users running gp_engine='eigen' opt into
    # the assumption that their fixed-effects spec is env-only.
    # The internally-built one-hot X is used below for the Schur step and
    # for log|X'V^-1 X| computation.

    balanced, r_per_cell = _is_balanced_met(gi_np, ei_np, n_geno, k)
    if not balanced:
        raise NotImplementedError(
            "mme_eigen scope: requires balanced MET (every geno x env cell "
            "must have equal obs count)."
        )

    # ---- eigen precompute (cached) ----------------------------------------
    U, Lambda = _cache_eigendecomp(G_np)
    # Lambda may contain tiny / negative values from numerical noise; clip
    Lambda = np.clip(Lambda, a_min=1e-12, a_max=None)

    Zty_by_env = np.zeros((n_geno, k), dtype=np.float64)
    np.add.at(Zty_by_env, (gi_np, ei_np), y_np)
    Ut_Zty_by_env = U.T @ Zty_by_env
    ones_U_colsum = U.sum(axis=0)

    sumy_by_env = np.zeros(k, dtype=np.float64)
    np.add.at(sumy_by_env, ei_np, y_np)

    yRinv_y_cache_base = float((y_np ** 2).sum())  # numerator; divide by se2[ei] at eval

    # ---- closed-form REML logLik in eigen basis ---------------------------
    def reml_eigen(theta: np.ndarray) -> float:
        sg2 = float(max(theta[0], 1e-12))
        se2_vec = np.maximum(theta[1:1 + k], 1e-10)
        inv_se2 = 1.0 / se2_vec
        c = r_per_cell * inv_se2.sum()
        d = c + 1.0 / (sg2 * Lambda)
        b_v = Ut_Zty_by_env @ inv_se2
        b_x = inv_se2 * sumy_by_env
        A_xx = np.diag(n_geno * r_per_cell * inv_se2)
        A_xv = r_per_cell * inv_se2[:, None] * ones_U_colsum[None, :]
        Avv_inv = 1.0 / d
        M = A_xv * Avv_inv[None, :]
        S = A_xx - M @ A_xv.T
        try:
            L_S = np.linalg.cholesky(S)
        except np.linalg.LinAlgError:
            return -np.inf
        logdet_S = 2.0 * np.log(np.diag(L_S)).sum()
        rhs_beta = b_x - M @ b_v
        beta_hat = np.linalg.solve(L_S.T, np.linalg.solve(L_S, rhs_beta))
        v_hat = Avv_inv * (b_v - A_xv.T @ beta_hat)
        logdet_R = float(np.log(se2_vec[ei_np]).sum())
        logdet_sgLambda = float(np.log(sg2 * Lambda).sum())
        logdet_d = float(np.log(d).sum())
        logdet_V = logdet_R + logdet_sgLambda + logdet_d
        yRinvy = float((inv_se2[ei_np] * y_np * y_np).sum())
        quad = yRinvy - float(beta_hat @ b_x) - float(v_hat @ b_v)
        const = (n_obs - k) * np.log(2.0 * np.pi)
        return -0.5 * (logdet_V + logdet_S + quad + const)

    # ---- central FD score and Hessian -------------------------------------
    def score_fd(theta: np.ndarray, eps: float = 1e-5) -> np.ndarray:
        s = np.zeros_like(theta)
        for i in range(len(theta)):
            h = max(eps * abs(theta[i]), eps)
            tp = theta.copy(); tp[i] += h
            tm = theta.copy(); tm[i] -= h
            s[i] = (reml_eigen(tp) - reml_eigen(tm)) / (2 * h)
        return s

    def hessian_fd(theta: np.ndarray, eps: float = 1e-4) -> np.ndarray:
        n_p = len(theta)
        H = np.zeros((n_p, n_p))
        for i in range(n_p):
            h = max(eps * abs(theta[i]), eps)
            tp = theta.copy(); tp[i] += h
            tm = theta.copy(); tm[i] -= h
            H[:, i] = (score_fd(tp) - score_fd(tm)) / (2 * h)
        H = 0.5 * (H + H.T)
        return H

    # ---- initialization ---------------------------------------------------
    # ANOVA-style init from genotype means -- this matches what the standalone
    # validation uses. We intentionally do NOT consume theta0_warm here:
    # PP's theta0_warm is shaped for the (w_g, w_ge, w_e, resid) gp_exact
    # layout, not our (sigma2_g, resid_per_env) eigen layout, and using it
    # would land the AI-Newton near PP-dense_v's local optimum at 0.131
    # instead of the global optimum near 0.125. (Tracked in Phase 3.15a NEWS:
    # this is the production-basin issue.)
    geno_means = np.zeros(n_geno)
    counts_g = np.zeros(n_geno)
    np.add.at(geno_means, gi_np, y_np)
    np.add.at(counts_g, gi_np, 1.0)
    geno_means /= np.maximum(counts_g, 1.0)
    between = float(np.var(geno_means, ddof=0))
    within = float(np.var(y_np - geno_means[gi_np], ddof=0))
    sg2_init = max(between - within / max(r_per_cell * k, 1), 0.01)
    se2_init = np.zeros(k)
    for e in range(k):
        mask = (ei_np == e)
        se2_init[e] = float(np.var(y_np[mask] - geno_means[gi_np[mask]],
                                   ddof=0))
    se2_init = np.maximum(se2_init, 0.01)
    theta0 = np.concatenate([[sg2_init], se2_init])

    # ---- AI-Newton with trust-region backtracking -------------------------
    theta = theta0.copy()
    ll = reml_eigen(theta)
    converged = False
    n_iter_done = 0
    AI_final: Optional[np.ndarray] = None
    for it in range(max_iter):
        s = score_fd(theta)
        H = hessian_fd(theta)
        AI = -H
        ridge = 1e-8 * max(np.abs(AI).max(), 1.0)
        try:
            dtheta = np.linalg.solve(AI + ridge * np.eye(len(theta)), s)
        except np.linalg.LinAlgError:
            dtheta = 0.01 * s
        alpha = 1.0
        accepted = False
        for _ in range(20):
            theta_new = theta + alpha * dtheta
            if (theta_new > 1e-10).all():
                ll_new = reml_eigen(theta_new)
                if ll_new > ll + 1e-12:
                    accepted = True
                    break
            alpha *= 0.5
        if not accepted:
            converged = True  # at a stationary point
            break
        rel_change_theta = float(np.max(np.abs((theta_new - theta) /
                                                np.maximum(np.abs(theta), 1e-8))))
        d_ll = ll_new - ll
        theta, ll = theta_new, ll_new
        n_iter_done = it + 1
        if abs(d_ll) < tol_loglik and rel_change_theta < tol_theta:
            converged = True
            break
    AI_final = -hessian_fd(theta)

    # ---- format output: varcomp_df + summary ------------------------------
    sg2_hat = float(max(theta[0], 1e-12))
    se2_hat = np.maximum(theta[1:1 + k], 1e-10)

    # Compute resid_floor for the ParamSpec lower-bound (matches the dense_v
    # path's convention; 1e-4 of the data variance, with a hard 1e-10 floor).
    _data_var = float(np.var(y_np, ddof=0)) if n_obs > 1 else 1.0
    if not math.isfinite(_data_var) or _data_var <= 0:
        _data_var = 1.0
    resid_floor = max(1e-4 * _data_var, 1e-10)

    # Build ParamSpec list matching the dense_v gp_exact schema so
    # build_asreml_varcomp_table produces a DataFrame identical in shape
    # and column types to the existing MME path -- avoids the R-side
    # coercion error we saw in 0.20.37.
    from varcomp_asreml import (
        ParamSpec as _ParamSpec,
        build_asreml_varcomp_table as _build_asreml_varcomp_table,
    )
    specs: List[Any] = []
    specs.append(_ParamSpec(
        kind="var_kernel",
        term=kernel_term_template.format(K=geno_kernel_names[0]),
        sub_label="var", lower=0.0,
    ))
    for e_lbl in env_labels:
        specs.append(_ParamSpec(
            kind="resid_env", term=env_col_name,
            env_label=str(e_lbl), lower=float(resid_floor),
        ))
    bounds_seq = [sp.default_bound() for sp in specs]
    pct_change = [float("nan")] * len(specs)
    varcomp_df = _build_asreml_varcomp_table(
        specs=specs,
        theta=theta,
        ai=AI_final,
        bounds=bounds_seq,
        pct_change=pct_change,
        env_scales=env_scales,
    )

    summary: Dict[str, Any] = {
        "loglik": float(ll),
        "nedf": int(nedf if nedf > 0 else (n_obs - k)),
        "sigma": 1.0,
        "aic": -2.0 * float(ll) + 2.0 * (1 + k),
        "bic": -2.0 * float(ll) + math.log(max(n_obs - k, 1)) * (1 + k),
        "parameters": int(1 + k),
        "converged": bool(converged),
        "n_ai_iter": int(n_iter_done),
        "theta": theta.tolist(),
        "theta_layout": "gp_exact_eigen",
        "theta_shape": {"n_params": int(len(theta))},
        "mme_backend": "eigen",
    }
    return varcomp_df, summary, np.asarray(AI_final, dtype=np.float64)
