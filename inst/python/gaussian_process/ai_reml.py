"""Re-engineered AI-REML pipeline used only when output_level='full_vc'.

This module is the authoritative home for variance-component estimation
logic: REML log-likelihood, dV/dtheta builders, AI matrix assembly,
boundary-aware AI-Newton driver, and the dense VC backend wired for the
framework's gp_exact / gp_icm_fa paths.

Design notes
------------
- Imports are deferred where possible so that merely importing this module
  on the fast path (predict_only / predict_with_se) does not pull in torch
  linear-operator machinery.
- The existing framework-file functions (`reml_loglik_from_V`,
  `compute_AI_Matrix_SEs`, `ai_newton_optimize`) are wrapped here rather
  than duplicated, preserving single-source-of-truth until Tasks #5 and
  later migrate the rewritten bodies.
- Active-set boundary handling (Task #2) is implemented here via
  `compute_ai_matrix_with_bounds`, which consumes the raw AI matrix from
  the framework helper and delegates inversion to
  `varcomp_asreml.active_set_inverse`. This isolates the SE derivation
  in one auditable place.

SE derivation summary (the user's concern)
------------------------------------------
1. AI matrix is computed as
       AI_ij = 0.5 * [ y' P D_i P D_j P y + tr(P D_i P D_j) ]
   with P the REML projector. Trace is a Rademacher-Hutchinson estimate
   using `ai_hutch_samples` probes; probe noise dominates for large
   `ai_hutch_samples`, so SEs are stable only when probe count is adequate
   (recommended >= 256 for reference comparisons).
2. For bound parameters (variances clipped at 0, or fixed), rows/cols are
   removed from AI before inversion (active set). Their SE is NaN, which
   matches ASReml's NA convention.
3. SE_i = sqrt( [AI_active^-1]_ii ). Under per-env standardization, the
   reported SE is back-transformed by the same factor as the estimate
   (linear reparameterization -> same scale for estimate and SE).
4. For FA parameterization (psi_e, lambda_e) the AI is built directly in
   that basis (Task #5), so no delta-method Jacobian is needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import math
import numpy as np

from varcomp_asreml import (
    ParamSpec,
    BoundFlag,
    EnvScales,
    AIRemlResult,
    active_set_inverse,
    apply_boundary_projection,
    ai_reml_optimize as _ai_reml_optimize,
    build_asreml_varcomp_table,
    build_summary,
    DenseVCBackend,
    VCBackend,
)

from mme_reml import (
    mme_precompute,
    reml_loglik_from_assembly,
    ai_and_score,
    select_backend,
)


# ---------------------------------------------------------------------------
# Context object: everything the AI path needs from the caller
# ---------------------------------------------------------------------------

@dataclass
class DenseAIContext:
    """Plain-data container bundling the tensors and callables for AI-REML.

    The caller (gp_exact_with_X / gp_icm_fa_with_X) constructs this and
    hands it to `run_ai_reml`. All torch-specific work is kept behind the
    `v_builder`, `reml_loglik_fn`, and `ai_matrix_fn` callables, so this
    module does not import torch at module scope.
    """
    specs: List[ParamSpec]
    theta0: np.ndarray
    env_scales: Optional[EnvScales]
    reml_loglik_fn: Callable[[np.ndarray], float]       # theta -> loglik
    score_fn: Callable[[np.ndarray], np.ndarray]         # theta -> s (p,)
    ai_matrix_fn: Callable[[np.ndarray], np.ndarray]     # theta -> AI (p, p)
    nedf: int


# ---------------------------------------------------------------------------
# Active-set AI inversion (Task #2 core)
# ---------------------------------------------------------------------------

def ai_covariance_with_bounds(
    ai: np.ndarray,
    bounds: Sequence[BoundFlag],
    jitter: float = 1e-10,
) -> np.ndarray:
    """Return the covariance of theta under REML asymptotics, honoring bounds.

    Thin re-export of `varcomp_asreml.active_set_inverse`; kept here so
    downstream code can import AI machinery from one module.
    """
    return active_set_inverse(ai, bounds, jitter=jitter)


# ---------------------------------------------------------------------------
# AI-REML driver wrapped in boundary-aware optimizer
# ---------------------------------------------------------------------------

def run_ai_reml(
    ctx: DenseAIContext,
    *,
    max_iter: int = 50,
    tol_loglik: float = 1e-4,
    tol_theta: float = 1e-6,
    damping: float = 1.0,
    max_halvings: int = 8,
) -> AIRemlResult:
    """Run boundary-aware AI-Newton on the supplied context.

    Returns an AIRemlResult. The caller formats output via
    `build_varcomp_outputs(...)`.
    """
    backend = DenseVCBackend(
        specs=list(ctx.specs),
        loglik_fn=ctx.reml_loglik_fn,
        score_fn=ctx.score_fn,
        ai_fn=ctx.ai_matrix_fn,
    )
    return _ai_reml_optimize(
        backend,
        ctx.theta0,
        max_iter=max_iter,
        tol_loglik=tol_loglik,
        tol_theta=tol_theta,
        damping=damping,
        max_halvings=max_halvings,
    )


# ---------------------------------------------------------------------------
# Output assembly
# ---------------------------------------------------------------------------

def build_varcomp_outputs(
    ctx: DenseAIContext,
    result: AIRemlResult,
    *,
    extra_summary: Optional[Dict[str, Any]] = None,
) -> Tuple["pd.DataFrame", Dict[str, Any], np.ndarray]:
    """Assemble the ASReml-style varcomp table and summary dict.

    Returns
    -------
    varcomp_df : pandas DataFrame with columns
        ['component','estimate','std.error','z.ratio','bound','pct_change']
    summary : dict
        {'loglik','nedf','sigma','aic','bic','parameters','converged','n_ai_iter', ...}
    ai_matrix : the final AI matrix used for SE derivation (numpy, symmetric)
    """
    varcomp_df = build_asreml_varcomp_table(
        specs=ctx.specs,
        theta=result.theta,
        ai=result.ai,
        bounds=result.bounds,
        pct_change=result.pct_change,
        env_scales=ctx.env_scales,
    )
    parameters_free = int(sum(b not in ("B", "F") for b in result.bounds))
    summary = build_summary(
        loglik=float(result.loglik),
        nedf=int(ctx.nedf),
        parameters_free=parameters_free,
        n_ai_iter=int(result.n_iter),
        converged=bool(result.converged),
        sigma=1.0,
        extra=extra_summary,
    )
    return varcomp_df, summary, np.asarray(result.ai, dtype=np.float64)


# ---------------------------------------------------------------------------
# FA (psi, lambda) dV builders (Task #5 skeleton)
# ---------------------------------------------------------------------------

def make_fa_rank1_dv_builders(
    ei_obs: np.ndarray,
    num_envs: int,
    lambda_vec: np.ndarray,
) -> List[Tuple[str, Callable[[np.ndarray], np.ndarray]]]:
    """Return dV/dpsi_e and dV/dlambda_e matvec callables for rank-1 FA.

    In the env-kernel block, let e_e be the env indicator in R^k. Then
      Sigma_env = diag(psi) + lambda lambda^T
      dSigma/dpsi_e    = e_e e_e^T
      dSigma/dlambda_e = e_e lambda^T + lambda e_e^T

    Each builder applies dSigma to an obs-space vector v by (1) aggregating
    v per env, (2) multiplying by dSigma, (3) scattering back to obs space.
    No torch import at module level; operates on numpy arrays supplied by
    the caller. (Caller may wrap with torch if hot-path performance matters
    once Task #5 wires this in.)

    NOTE: This is the Task #5 skeleton. It is NOT yet wired into gp_icm_fa
    — that step comes when the FA rewrite lands. Kept here now so the
    contract is visible alongside the AI/SE derivation.
    """
    ei = np.asarray(ei_obs, dtype=np.int64)
    k = int(num_envs)
    lam = np.asarray(lambda_vec, dtype=np.float64).reshape(k)

    def _aggregate(v: np.ndarray) -> np.ndarray:
        t = np.zeros(k, dtype=np.float64)
        np.add.at(t, ei, v)
        return t

    def _scatter(u: np.ndarray) -> np.ndarray:
        return u[ei]

    builders: List[Tuple[str, Callable[[np.ndarray], np.ndarray]]] = []

    for e in range(k):
        def make_psi(e=e):
            def DV(v: np.ndarray) -> np.ndarray:
                t = _aggregate(v)
                u = np.zeros(k, dtype=np.float64)
                u[e] = t[e]
                return _scatter(u)
            return DV
        builders.append((f"psi[{e}]", make_psi()))

    for e in range(k):
        def make_lam(e=e):
            def DV(v: np.ndarray) -> np.ndarray:
                t = _aggregate(v)
                u = np.zeros(k, dtype=np.float64)
                u[e] += float(lam @ t)
                u += lam * float(t[e])
                return _scatter(u)
            return DV
        builders.append((f"lambda[{e}]", make_lam()))

    return builders


# ---------------------------------------------------------------------------
# Exact dense AI + score computation (replaces Hutchinson trace)
# ---------------------------------------------------------------------------

def _exact_dense_precompute(
    y: "torch.Tensor",
    X: "torch.Tensor",
    V_dense: "torch.Tensor",
    jitter: float = 1e-6,
):
    """Shared Cholesky precompute for REML loglik AND score/AI.

    Before this helper, ``_exact_dense_reml_loglik`` and
    ``_exact_dense_score_and_AI`` each computed L_V = cholesky(V) and
    L_xvx = cholesky(X' V^-1 X) independently. In the AI-Newton inner loop
    this meant TWO V Choleskys per accepted iteration -- once in the
    trust-region loglik check, once when the optimizer asked for
    score / AI at the same theta.

    Factoring the shared work into a single precompute lets the
    `_ensure` caches in ``build_*_ai_varcomp`` retain a `pre` object
    alongside the cached loglik. When a later `ai()` call lands on the
    same theta, we skip the Cholesky entirely and feed `pre` directly to
    `_exact_dense_score_and_AI_from_pre`. Empirically halves the
    per-iteration cost on a 4-env MET wheat fit.

    Note: Py = V^-1 r where r = y - X beta. By linearity this equals
    Vinv_y - Vinv_X @ beta, so the same vector serves both the REML
    quadratic form (loglik) and the score / AI matvec (Py = P y).

    Returns
    -------
    dict
        Keys: L_V, Vinv_y, Vinv_X, XtVinvX_jit, L_xvx, beta, Py,
        logdet_V, logdet_XtVinvX, quad, n, p, device, dtype.
    """
    import torch

    n = int(y.shape[0])
    p = int(X.shape[1])
    device = y.device
    dtype = y.dtype
    I_n = torch.eye(n, device=device, dtype=dtype)
    I_p = torch.eye(p, device=device, dtype=dtype)

    L_V = torch.linalg.cholesky(V_dense + jitter * I_n)
    logdet_V = 2.0 * torch.log(torch.diag(L_V)).sum()

    Vinv_y = torch.cholesky_solve(y.unsqueeze(-1), L_V).squeeze(-1)
    Vinv_X = torch.cholesky_solve(X, L_V)
    XtVinvX_jit = X.T @ Vinv_X + jitter * I_p
    L_xvx = torch.linalg.cholesky(XtVinvX_jit)
    logdet_XtVinvX = 2.0 * torch.log(torch.diag(L_xvx)).sum()

    XtVinvy = X.T @ Vinv_y
    beta = torch.cholesky_solve(XtVinvy.unsqueeze(-1), L_xvx).squeeze(-1)
    Py = Vinv_y - Vinv_X @ beta
    quad = ((y - X @ beta) * Py).sum()  # r' V^-1 r = r' Py

    return {
        "L_V": L_V,
        "Vinv_y": Vinv_y,
        "Vinv_X": Vinv_X,
        "XtVinvX_jit": XtVinvX_jit,
        "L_xvx": L_xvx,
        "beta": beta,
        "Py": Py,
        "logdet_V": logdet_V,
        "logdet_XtVinvX": logdet_XtVinvX,
        "quad": quad,
        "n": n,
        "p": p,
        "device": device,
        "dtype": dtype,
        "jitter": float(jitter),
    }


def _exact_dense_reml_loglik_from_pre(pre) -> float:
    """REML log-likelihood from a shared precompute dict."""
    import math

    const = (pre["n"] - pre["p"]) * math.log(2.0 * math.pi)
    ll = -0.5 * (pre["logdet_V"] + pre["logdet_XtVinvX"] + pre["quad"] + const)
    return float(ll.detach().cpu().item())


def _exact_dense_reml_loglik(
    y: "torch.Tensor",
    X: "torch.Tensor",
    V_dense: "torch.Tensor",
    jitter: float = 1e-6,
) -> float:
    """Deterministic dense-Cholesky REML log-likelihood.

    Thin wrapper around ``_exact_dense_precompute`` +
    ``_exact_dense_reml_loglik_from_pre`` for callers that don't want to
    handle the shared precompute. The Cholesky cost is identical to the
    pre-0.20.17 implementation; the savings come when a downstream
    ``_exact_dense_score_and_AI_from_pre`` call reuses the same `pre`.

    Replaces ``fw.reml_loglik_from_V`` on the AI-REML hot path so the
    loglik is bit-reproducible at a fixed theta. ``V_op.logdet()`` and
    ``V_op.solve()`` in the framework use Lanczos/CG with stochastic probes,
    which produces noise of order ~0.5 in the returned loglik -- larger than
    any tol_loglik and catastrophic for the trust-region rho test.
    """
    pre = _exact_dense_precompute(y, X, V_dense, jitter=jitter)
    return _exact_dense_reml_loglik_from_pre(pre)


def _exact_dense_score_and_AI_from_pre(
    pre,
    DV_matvecs: list,
    DV_dense_or_diag: list,
):
    """Score + AI from a shared precompute dict (the fused-Cholesky path).

    Identical math to ``_exact_dense_score_and_AI`` but uses the L_V,
    L_xvx, Vinv_X, Py already in `pre`. Saves one full n x n Cholesky
    per call when paired with `_exact_dense_precompute` upstream.
    """
    import torch

    L_V       = pre["L_V"]
    Vinv_X    = pre["Vinv_X"]
    L_xvx     = pre["L_xvx"]
    Py        = pre["Py"]
    n         = pre["n"]
    p         = pre["p"]
    device    = pre["device"]
    dtype     = pre["dtype"]
    q         = len(DV_matvecs)

    I_n = torch.eye(n, device=device, dtype=dtype)
    I_p = torch.eye(p, device=device, dtype=dtype)

    # ---- F columns: F_j = dV/dtheta_j @ Py (one matvec per parameter). The
    #      caller's DV_matvecs encodes how each variance parameter perturbs V,
    #      so a multi-kernel V = sum_k w_k K_k + R is handled the same way as
    #      a single-kernel V: one callable per kernel weight produces K_k @ Py.
    F = torch.stack([dv(Py) for dv in DV_matvecs], dim=1)                # (n, q)

    # ---- PF = V^-1 F - V^-1 X (X' V^-1 X)^-1 X' V^-1 F ----
    # X' V^-1 F = (V^-1 X)' F because V is symmetric, so we can use Vinv_X
    # directly without needing X in the precompute dict.
    Vinv_F = torch.cholesky_solve(F, L_V)                                # (n, q)
    XtVinv_F = Vinv_X.T @ F                                              # (p, q)
    correction = Vinv_X @ torch.cholesky_solve(XtVinv_F, L_xvx)           # (n, q)
    PF = Vinv_F - correction                                             # (n, q)

    # ---- AI = 0.5 * F' @ PF (exact, no stochastic noise) ----
    AI = 0.5 * (F.T @ PF)
    AI = 0.5 * (AI + AI.T)                                                # enforce symmetry

    # ---- Exact score: s_i = 0.5 * (-tr(P D_i) + Py' D_i Py) ----
    Vinv = torch.cholesky_solve(I_n, L_V)
    C    = torch.cholesky_solve(I_p, L_xvx)
    P    = Vinv - Vinv_X @ C @ Vinv_X.T

    score = torch.empty(q, device=device, dtype=dtype)
    for i in range(q):
        quad = float(F[:, i] @ Py)
        D_i = DV_dense_or_diag[i]
        if D_i is not None:
            if D_i.dim() == 1:
                tr_PDi = float((P.diag() * D_i).sum())
            else:
                tr_PDi = float((P * D_i).sum())
        else:
            tr_PDi = 0.0
        score[i] = 0.5 * (-tr_PDi + quad)

    return (
        score.detach().cpu().numpy().astype(np.float64),
        AI.detach().cpu().numpy().astype(np.float64),
    )


def _exact_dense_score_and_AI(
    y: "torch.Tensor",
    X: "torch.Tensor",
    V_dense: "torch.Tensor",
    DV_matvecs: list,
    DV_dense_or_diag: list,
    jitter: float = 1e-6,
):
    """Exact REML score and AI matrix using dense Cholesky.

    Thin wrapper around ``_exact_dense_precompute`` +
    ``_exact_dense_score_and_AI_from_pre`` for callers that don't want to
    handle the shared precompute. The Cholesky cost is identical to the
    pre-0.20.17 implementation; the savings come when a downstream
    ``_exact_dense_reml_loglik_from_pre`` call reuses the same `pre`.

    Replaces the Hutchinson-trace-based ``compute_score_and_AI_operator``.
    The AI matrix is computed as::

        F_j = (dV/dtheta_j) @ Py          (one matvec per parameter)
        PF  = V^-1 F - V^-1 X (X' V^-1 X)^-1 X' V^-1 F
        AI  = 0.5 * F' @ PF               (exact, no stochastic noise)

    The REML score trace ``tr(P dV/dtheta_i)`` is computed exactly from
    the dense REML projector P.

    Parameters
    ----------
    y : (n,) tensor
    X : (n, p) tensor -- fixed-effects design matrix
    V_dense : (n, n) tensor -- full covariance matrix (dense)
    DV_matvecs : list of q callables, each ``v -> dV/dtheta_i @ v``
    DV_dense_or_diag : list of q items, each either:
        - (n, n) tensor: dense derivative matrix dV/dtheta_i
        - (n,) tensor: diagonal of dV/dtheta_i (for per-env residuals)
        - None: use 0 for the trace term (quadratic-only score)
    jitter : float -- Cholesky regularization

    Returns
    -------
    score : (q,) numpy array -- exact REML score
    AI    : (q, q) numpy array -- exact Average Information matrix
    """
    pre = _exact_dense_precompute(y, X, V_dense, jitter=jitter)
    return _exact_dense_score_and_AI_from_pre(pre, DV_matvecs, DV_dense_or_diag)


# ---------------------------------------------------------------------------
# Torch-side FA (psi, lambda) AI builder — the full Task #5 implementation
# ---------------------------------------------------------------------------

def _get_framework_module():
    """Late-bind the framework module without re-executing import."""
    import importlib, sys
    fw_name = "gp_framework"
    if fw_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(fw_name, fw_name + ".py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        sys.modules[fw_name] = mod
    return sys.modules[fw_name]


def build_fa_icm_ai_varcomp(
    model,
    likelihood,
    train_x,
    train_y,
    fa_module,
    *,
    env_labels: Sequence[str],
    geno_kernel_names: Sequence[str],
    term_name: str = "fa(Env,1):vm(GID,GAinv)",
    env_col_name: str = "Env",
    hutch_samples: int = 256,
    jitter: float = 1e-6,
    seed: int = 12345,
    env_scales: Optional["EnvScales"] = None,
    nedf: int = 0,
    max_iter: int = 30,
    tol_loglik: float = 1e-3,
    tol_theta: float = 1e-4,
    theta0_warm: Optional[np.ndarray] = None,
    max_iter_warm: Optional[int] = None,
):
    """Boundary-aware AI-Newton REML on GP-ICM-FA. Emits ASReml-schema output.

    Parameter vector: (psi_0..psi_{k-1}, lambda_{0,0}..lambda_{k-1, rank-1},
    resid_0..resid_{k-1}). The AI-Newton loop refines psi/lambda AND per-env
    residuals jointly so the residual estimates aren't pinned at the Adam-
    converged values.

    Returns: (varcomp_df, summary_dict, ai_matrix_np).
    """
    import torch
    import gpytorch
    from linear_operator.operators import DiagLinearOperator

    device = train_y.device
    dtype_t = train_y.dtype
    k = int(fa_module.num_envs)
    rank = int(fa_module.rank)
    ei_t = train_x[:, 1].long()
    gi_t = train_x[:, 0].long()
    n_obs = int(train_y.shape[0])
    Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
    fw = _get_framework_module()

    # ---- Initial state from converged Adam ---------------------------------
    with torch.no_grad():
        psi0 = (torch.nn.functional.softplus(fa_module.d_unconstrained) + 1e-8).detach().clone()
        L0 = fa_module.L.detach().clone()
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            init_noise = likelihood.noise.detach()
        else:
            init_noise = torch.full((n_obs,), float(likelihood.noise.item()),
                                    dtype=dtype_t, device=device)
        resid0 = torch.zeros(k, dtype=torch.float64, device=device)
        for e in range(k):
            m = (ei_t == e)
            resid0[e] = (init_noise[m].mean() if m.any() else init_noise.mean()).to(torch.float64)

    theta0_np = np.concatenate([
        psi0.cpu().to(torch.float64).numpy(),
        *[L0[:, r].cpu().to(torch.float64).numpy() for r in range(rank)],
        resid0.cpu().numpy(),
    ])

    # ---- Param specs --------------------------------------------------------
    specs: List[ParamSpec] = []
    for e in range(k):
        specs.append(ParamSpec(kind="fa_psi", term=term_name,
                               env_label=str(env_labels[e]), sub_label="var", lower=0.0))
    for r in range(rank):
        sub = f"fa{r+1}"
        for e in range(k):
            specs.append(ParamSpec(kind="fa_lambda", term=term_name,
                                   env_label=str(env_labels[e]), sub_label=sub))
    for e in range(k):
        specs.append(ParamSpec(kind="resid_env", term=env_col_name,
                               env_label=str(env_labels[e]), lower=0.0))

    # ---- State setter: theta (numpy) -> writes into fa_module --------------
    def _softplus_inv_torch(y):
        # robust inverse softplus for y > 0
        y_clamped = torch.clamp(y, min=1e-12)
        return torch.where(y_clamped > 20.0, y_clamped, torch.log(torch.expm1(y_clamped)))

    def _set_state(theta_np: np.ndarray):
        theta_t = torch.as_tensor(theta_np, dtype=torch.float64, device=device)
        idx = 0
        psi = torch.clamp(theta_t[idx:idx+k].to(dtype_t), min=1e-10); idx += k
        d_uncon = _softplus_inv_torch(psi - 1e-8)
        fa_module.d_unconstrained.data.copy_(d_uncon.reshape(fa_module.d_unconstrained.shape))
        for r in range(rank):
            lam = theta_t[idx:idx+k].to(dtype_t); idx += k
            fa_module.L.data[:, r].copy_(lam)
        resid = torch.clamp(theta_t[idx:idx+k].to(dtype_t), min=1e-10); idx += k
        return resid[ei_t]   # per-obs noise diag

    # ---- DV builders parameterized by current model state ------------------
    G_list_torch = model.covar_module.G_list
    # When learn_scales=True, the effective per-kernel scales live in
    # w_fa_unconstrained (softplus-parameterized) and can differ substantially
    # from w_fa_fix (the initial-value buffer frozen at model construction).
    # Using the stale buffer here silently breaks the psi/lambda score on
    # multi-kernel problems — the AI trust region then rejects every step
    # because analytic gradient and true loglik disagree.
    if hasattr(model.covar_module, "get_effective_scales"):
        w_fa_eff, _ = model.covar_module.get_effective_scales()
        w_fa = w_fa_eff.detach().to(dtype=dtype_t, device=device)
    else:
        w_fa = getattr(model.covar_module, "w_fa_fix", None)
        if w_fa is None:
            w_fa = torch.ones(len(geno_kernel_names), dtype=dtype_t, device=device)
        w_fa = w_fa.detach().to(dtype=dtype_t, device=device)

    def _gxe_matvec(dSe: "torch.Tensor", v: "torch.Tensor") -> "torch.Tensor":
        """Apply (sum_g w_fa[g] G_g) (.) (kron) dSe (env block) to obs vector."""
        out = torch.zeros_like(v)
        for r_g, G in enumerate(G_list_torch):
            wg = float(w_fa[r_g]) if r_g < w_fa.numel() else 1.0
            ng = G.shape[0]
            env_geno = [None] * k
            for e_in in range(k):
                m_e = (ei_t == e_in)
                if m_e.any():
                    g_acc = torch.zeros(ng, dtype=dtype_t, device=device)
                    g_acc.scatter_add_(0, gi_t[m_e], v[m_e])
                    env_geno[e_in] = G @ g_acc
            for f_out in range(k):
                m_f = (ei_t == f_out)
                if not m_f.any():
                    continue
                gids_f = gi_t[m_f]
                accum = torch.zeros(int(m_f.sum()), dtype=dtype_t, device=device)
                for e_in in range(k):
                    if env_geno[e_in] is None:
                        continue
                    coef = float(dSe[f_out, e_in])
                    if coef == 0.0:
                        continue
                    accum = accum + (wg * coef) * env_geno[e_in].index_select(0, gids_f)
                out[m_f] = out[m_f] + accum
        return out

    def _build_dv_list():
        DV: List[Callable] = []
        for e in range(k):
            def DVpsi(v, e=e):
                dSe = torch.zeros(k, k, dtype=dtype_t, device=device)
                dSe[e, e] = 1.0
                return _gxe_matvec(dSe, v)
            DV.append(DVpsi)
        for r in range(rank):
            for e in range(k):
                def DVlam(v, e=e, r=r):
                    lam_r = fa_module.L.detach()[:, r]
                    dSe = torch.zeros(k, k, dtype=dtype_t, device=device)
                    dSe[e, :] += lam_r
                    dSe[:, e] += lam_r
                    return _gxe_matvec(dSe, v)
                DV.append(DVlam)
        for e in range(k):
            def DVres(v, e=e):
                m = (ei_t == e).to(dtype_t)
                return v * m
            DV.append(DVres)
        return DV

    # ---- DV dense matrices for FA (exact score traces) --------------------
    def _build_dv_dense_fa():
        """Build dense (n,n) or diagonal (n,) derivative matrices for FA params.

        For psi_e:  dV/dpsi_e = sum_g w_g[g] * (G_g kron e_e e_e^T) in obs space
        For lam_e:  dV/dlam_{e,r} = sum_g w_g[g] * (G_g kron (e_e lam_r^T + lam_r e_e^T))
        For resid:  dV/dresid_e = diag(indicator_e)
        """
        # Precompute per-env obs indices
        env_indices = [(ei_t == e).nonzero(as_tuple=True)[0] for e in range(k)]

        def _build_gxe_dense(dSe):
            """Build n×n dense matrix: sum_g w_g[g] * (G_g kron dSe) in obs space."""
            D = torch.zeros(n_obs, n_obs, dtype=dtype_t, device=device)
            for g, G in enumerate(G_list_torch):
                wg = float(w_fa[g]) if g < w_fa.numel() else 1.0
                for e1 in range(k):
                    idx1 = env_indices[e1]
                    if idx1.numel() == 0:
                        continue
                    gi1 = gi_t[idx1]
                    for e2 in range(k):
                        coef = float(dSe[e1, e2])
                        if abs(coef) < 1e-15:
                            continue
                        idx2 = env_indices[e2]
                        if idx2.numel() == 0:
                            continue
                        gi2 = gi_t[idx2]
                        block = (wg * coef) * G[gi1][:, gi2]
                        D[idx1.unsqueeze(1), idx2.unsqueeze(0)] += block
            return D

        DV_dense = []
        # psi derivatives
        for e in range(k):
            dSe = torch.zeros(k, k, dtype=dtype_t, device=device)
            dSe[e, e] = 1.0
            DV_dense.append(_build_gxe_dense(dSe))

        # lambda derivatives
        lam = fa_module.L.detach()  # (k, rank)
        for r in range(rank):
            for e in range(k):
                dSe = torch.zeros(k, k, dtype=dtype_t, device=device)
                dSe[e, :] += lam[:, r]
                dSe[:, e] += lam[:, r]
                DV_dense.append(_build_gxe_dense(dSe))

        # residual derivatives (diagonal)
        for e in range(k):
            DV_dense.append((ei_t == e).to(dtype_t))  # (n,) diagonal

        return DV_dense

    def _build_V(noise_diag_obs):
        K_op = model.covar_module(train_x, train_x)
        return K_op + DiagLinearOperator(noise_diag_obs)

    # ---- Cached numpy callables --------------------------------------------
    # The cache key is theta-bytes ONLY (pre-0.20.17 the key included
    # need_score_ai, so loglik(theta) and ai(theta) at the same theta
    # missed each other and triggered two independent V Choleskys per
    # accepted iteration). With the precompute we keep `pre` alongside
    # the loglik so a later ai() call at the same theta skips the
    # Cholesky entirely.
    cache = {"key": None, "ll": None, "score": None, "ai": None, "pre": None}

    def _ensure(theta_np: np.ndarray, need_score_ai: bool):
        key = theta_np.tobytes()
        if cache["key"] == key:
            if (not need_score_ai) or (cache["score"] is not None):
                return
            # Same theta, score+AI not yet computed -- reuse the cached pre
            # so we skip the V Cholesky. This is the big saving when the
            # optimizer's trust-region trial loglik is followed by an
            # ai() call at the accepted theta.
            with torch.no_grad():
                DV_matvecs = _build_dv_list()
                DV_dense = _build_dv_dense_fa()
                score_np, AI_np = _exact_dense_score_and_AI_from_pre(
                    cache["pre"], DV_matvecs, DV_dense,
                )
                cache["score"] = score_np
                cache["ai"] = AI_np
            return
        # New theta: full precompute (one V Cholesky), then loglik and
        # optionally score+AI off the shared `pre`.
        with torch.no_grad():
            noise_obs = _set_state(theta_np)
            V_op = _build_V(noise_obs)
            # Deterministic dense-Cholesky path: fw.reml_loglik_from_V uses
            # CG/Lanczos with ~0.5 noise per call which destroys the
            # AI-REML trust-region rho ratio.
            V_dense = V_op.to_dense()
            pre = _exact_dense_precompute(
                train_y, Xtt, V_dense, jitter=float(jitter),
            )
            cache["pre"] = pre
            cache["ll"] = _exact_dense_reml_loglik_from_pre(pre)
            if need_score_ai:
                DV_matvecs = _build_dv_list()
                DV_dense = _build_dv_dense_fa()
                score_np, AI_np = _exact_dense_score_and_AI_from_pre(
                    pre, DV_matvecs, DV_dense,
                )
                cache["score"] = score_np
                cache["ai"] = AI_np
            else:
                cache["score"] = None
                cache["ai"] = None
        cache["key"] = key

    def loglik_fn(theta_np):
        _ensure(theta_np, need_score_ai=False); return cache["ll"]
    def score_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["score"]
    def ai_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["ai"]

    # ---- Run AI-Newton ------------------------------------------------------
    import os as _os
    _verbose = bool(int(_os.environ.get("AIREML_VERBOSE", "0")))
    if _verbose:
        print(f"[fa_icm AI] theta0={np.round(theta0_np, 5).tolist()}  specs={[(s.kind, s.term, s.env_label) for s in specs]}", flush=True)
    backend = DenseVCBackend(specs=specs, loglik_fn=loglik_fn, score_fn=score_fn, ai_fn=ai_fn)

    # Warm-start: when theta0_warm is supplied, replace the Adam-init theta0
    # with the caller's pre-converged theta and tighten max_iter. Same-shape
    # warm-start is the within-CV-fold case (Phase 2.1). The boundary
    # projection in ai_reml_optimize will re-clip if any element drifted
    # past a boundary -- we explicitly pre-project here too so loglik() at
    # iter 0 sees a valid state. SHAPE MISMATCH (Phase 2.2 safety net): if
    # the cached theta is the wrong length (e.g. a 0.20.19 auto-warm-start
    # hit a key whose underlying spec changed between releases), warn and
    # fall back to cold init rather than raising -- the cold path always
    # works, the warm path is best-effort.
    theta_start = theta0_np
    iter_budget = int(max_iter)
    if theta0_warm is not None:
        warm = np.asarray(theta0_warm, dtype=np.float64).reshape(-1)
        if warm.shape[0] != theta0_np.shape[0]:
            import warnings as _wn
            _wn.warn(
                f"theta0_warm length {warm.shape[0]} does not match expected "
                f"{theta0_np.shape[0]} (psi+lambda+resid for k={k}, rank={rank}); "
                f"falling back to cold theta0 (Adam init).",
                RuntimeWarning,
            )
        else:
            theta_start, _bounds_check = apply_boundary_projection(warm, specs)
            if max_iter_warm is not None:
                iter_budget = int(max_iter_warm)
            else:
                # Default: 5 iters is enough for a same-shape warm restart on
                # a neighbouring data subset. The trust-region reject/retry
                # handles cases where the warm theta is too far from the new
                # optimum.
                iter_budget = min(5, int(max_iter))
            if _verbose:
                print(f"[fa_icm AI] WARM-START enabled; iter_budget={iter_budget}", flush=True)
    result = _ai_reml_optimize(
        backend, theta_start,
        max_iter=iter_budget,
        tol_loglik=float(tol_loglik),
        tol_theta=float(tol_theta),
        damping=1.0,
        max_halvings=8,
        verbose=_verbose,
    )
    # ---- Sign-flip identifiability (lesson #9) -----------------------------
    # FA loadings L are determined only up to a column-wise sign: L L^T is
    # invariant under L[:, r] -> -L[:, r]. Force L[r, r] > 0 so output matches
    # the ASReml convention (positive diagonal in the leading rank×rank block).
    # Estimates, loglik, and SEs are unchanged; AI transforms as J AI J with
    # J = diag(signs) and diag(AI) stays fixed, so SE = sqrt(diag(AI^{-1}))
    # is unchanged. Off-diagonal AI entries linking flipped <-> unflipped
    # parameters pick up the sign.
    theta_flipped = result.theta.copy()
    ai_flipped = np.asarray(result.ai, dtype=np.float64).copy()
    sign_vec = np.ones(theta_flipped.shape[0], dtype=np.float64)
    for r in range(rank):
        diag_idx = k + r * k + r
        if theta_flipped[diag_idx] < 0.0:
            lam_slice = slice(k + r * k, k + (r + 1) * k)
            theta_flipped[lam_slice] = -theta_flipped[lam_slice]
            sign_vec[lam_slice] = -1.0
    if not np.all(sign_vec == 1.0):
        ai_flipped = (sign_vec[:, None] * ai_flipped) * sign_vec[None, :]
        result.theta = theta_flipped
        result.ai = ai_flipped

    # Persist the (possibly sign-corrected) state in the model so downstream
    # prediction uses the convention-normalized loadings.
    _set_state(result.theta)

    # ---- Output -------------------------------------------------------------
    pct = np.where(np.isfinite(result.pct_change), result.pct_change, 0.0)
    varcomp_df = build_asreml_varcomp_table(
        specs=specs, theta=result.theta, ai=result.ai, bounds=result.bounds,
        pct_change=pct, env_scales=env_scales,
    )
    parameters_free = int(sum(b not in ("B", "F") for b in result.bounds))
    summary = build_summary(
        loglik=float(result.loglik), nedf=int(nedf),
        parameters_free=parameters_free, n_ai_iter=int(result.n_iter),
        converged=bool(result.converged), sigma=1.0,
    )
    # Delta-method SE for derived sigma2_g and h^2. sigma2_g is NOT a
    # direct REML parameter in FA; it is sum_k lambda_k^2 + psi_e, so SE
    # needs J^T Sigma_theta J on the active set.
    from varcomp_asreml import compute_fa_derived_summary as _compute_fa_derived_summary
    fa_derived = _compute_fa_derived_summary(
        specs=specs, theta=result.theta, ai=result.ai,
        bounds=result.bounds, env_scales=env_scales,
        kernel_weight_total=float(w_fa.detach().sum().cpu().item()),
    )
    if fa_derived is not None:
        summary = dict(summary)
        summary["fa_derived"] = fa_derived
    # Surface the converged variance-component theta and its layout so the
    # R-side warm-start cache can read it (Phase 2.1 Stage B). The Python
    # backend's spec layout for FA is:
    #   [psi[0..k-1], lambda[:,0..rank-1], resid_env[0..k-1]]
    # The R side stores this dict keyed by (model, response, gmatrix-hash)
    # and feeds it back in as theta0_warm on the next compatible fit.
    summary = dict(summary)
    summary["theta"] = result.theta.tolist()
    summary["theta_layout"] = "fa_icm"
    summary["theta_shape"] = {"k": int(k), "rank": int(rank), "n_params": int(len(result.theta))}
    return varcomp_df, summary, np.asarray(result.ai, dtype=np.float64)


# ---------------------------------------------------------------------------
# Torch-side active-set SE (Task #2)
# ---------------------------------------------------------------------------

def compute_ai_se_active_set(
    ai_torch,
    bounds: Sequence[BoundFlag],
    *,
    jitter: float = 1e-10,
) -> Tuple[Any, Any]:
    """Compute SEs from a torch AI tensor honoring bounds.

    Parameters
    ----------
    ai_torch : torch.Tensor of shape (p, p), symmetric.
    bounds   : length-p list of bound flags (P, B, F, U). Rows flagged B or F
               are excluded from the inversion; their SE is NaN.
    jitter   : ridge added to the active sub-block for Cholesky safety.

    Returns
    -------
    se_torch : (p,) tensor. NaN in slots where bound in {B, F}.
    cov_torch : (p, p) tensor. NaN-filled rows/cols where bound in {B, F};
                active sub-block contains the inverse-AI covariance.

    Rationale
    ---------
    Inverting the full AI when a parameter sits at its lower bound pollutes
    adjacent SEs because the Fisher information is singular in the direction
    of the boundary. ASReml handles this by treating bound parameters as
    inactive (fixed) during inversion; their reported SE is NA. This function
    mirrors that behavior on the torch side without copying to numpy.
    """
    import torch

    ai = ai_torch
    if not torch.is_tensor(ai):
        raise TypeError("ai_torch must be a torch.Tensor")
    p = int(ai.shape[0])
    if ai.shape != (p, p):
        raise ValueError(f"ai_torch must be square; got shape {tuple(ai.shape)}")
    if len(bounds) != p:
        raise ValueError(f"bounds length {len(bounds)} != AI size {p}")

    device, dtype = ai.device, ai.dtype
    active_mask = torch.tensor(
        [b not in ("B", "F") for b in bounds], device=device, dtype=torch.bool
    )
    nan = float("nan")
    cov = torch.full((p, p), nan, device=device, dtype=dtype)
    se = torch.full((p,), nan, device=device, dtype=dtype)

    if active_mask.any():
        idx = active_mask.nonzero(as_tuple=True)[0]
        sub = ai.index_select(0, idx).index_select(1, idx)
        sub = 0.5 * (sub + sub.T)
        m = sub.shape[0]
        I = torch.eye(m, device=device, dtype=dtype)

        inv_sub = None
        for mult in (1.0, 10.0, 100.0, 1000.0):
            try:
                L = torch.linalg.cholesky(sub + (jitter * mult) * I)
                inv_sub = torch.cholesky_solve(I, L)
                break
            except Exception:
                inv_sub = None
        if inv_sub is None:
            try:
                evals, evecs = torch.linalg.eigh(sub)
                evals = torch.clamp(evals, min=max(jitter, 1e-12))
                inv_sub = (evecs * evals.reciprocal().unsqueeze(0)) @ evecs.T
            except Exception:
                inv_sub = torch.linalg.pinv(sub)

        # Scatter inv_sub back into the full-size cov matrix
        grid_i = idx.view(-1, 1).expand(m, m)
        grid_j = idx.view(1, -1).expand(m, m)
        cov.index_put_((grid_i.reshape(-1), grid_j.reshape(-1)), inv_sub.reshape(-1))

        se_active = torch.sqrt(torch.clamp(torch.diag(inv_sub), min=0.0))
        se.index_put_((idx,), se_active)

    return se, cov


def bounds_from_specs_and_theta(
    specs: Sequence[ParamSpec],
    theta: np.ndarray,
    *,
    floor_eps: float = 1e-12,
) -> List[BoundFlag]:
    """Classify each (spec, theta_i) pair into a bound flag.

    Thin wrapper around `varcomp_asreml.apply_boundary_projection` that
    returns bounds without projecting (projection is the optimizer's job;
    classification is independently useful, e.g. for final SE assembly).
    """
    _projected, bounds = apply_boundary_projection(
        np.asarray(theta, dtype=np.float64), specs, floor_eps=floor_eps
    )
    return bounds


def build_gp_exact_ai_varcomp(
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
    env_scales: Optional["EnvScales"] = None,
    nedf: int = 0,
    max_iter: int = 12,
    tol_loglik: float = 1e-3,
    tol_theta: float = 1e-4,
    theta0_warm: Optional[np.ndarray] = None,
    max_iter_warm: Optional[int] = None,
    refine_kernel_scales: bool = True,
):
    """Boundary-aware AI-Newton REML for the gp_exact path.

    Parameter vector (in this order):
      - w_g[k]  for each genetic kernel k       -> "vm(GID,K)!var"
      - w_ge[k] for each genetic kernel k       -> "Env:vm(GID,K)!var" (if cached and multi-env)
      - w_e                                     -> "Env!var"           (if cached and multi-env)
      - sigma2_resid_e for each env             -> "<env_col>_<env>!R"

    If ``refine_kernel_scales`` is False, kernel weights are kept fixed and only
    per-env residuals are jointly optimized; this matches the common gp_exact
    default of learn_scales=False but still produces fitted residuals via REML.

    Returns ``(varcomp_df, summary_dict, ai_matrix_np)`` aligned with the
    schema used by build_fa_icm_ai_varcomp.
    """
    import torch
    import gpytorch
    from linear_operator.operators import DiagLinearOperator
    import os as _os
    _tll_env = _os.environ.get("PREDICTPRO_GP_TOL_LOGLIK", "").strip()
    _tth_env = _os.environ.get("PREDICTPRO_GP_TOL_THETA", "").strip()
    _mit_env = _os.environ.get("PREDICTPRO_GP_MAX_ITER", "").strip()
    _dem_env = _os.environ.get("PREDICTPRO_GP_DISABLE_ENV_MAIN", "").strip().lower()
    if _tll_env:
        try: tol_loglik = float(_tll_env)
        except ValueError: pass
    if _tth_env:
        try: tol_theta = float(_tth_env)
        except ValueError: pass
    if _mit_env:
        try: max_iter = int(_mit_env)
        except ValueError: pass
    _disable_env_main = _dem_env in ("1", "true", "yes", "on")

    device = train_y.device
    dtype_t = train_y.dtype
    ei_t = train_x[:, 1].long()
    n_obs = int(train_y.shape[0])
    Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
    fw = _get_framework_module()
    kernel = model.covar_module
    k = len(env_labels)
    nK = len(geno_kernel_names)

    # Cached dense component matrices stored on the kernel during forward.
    # The prediction path may have just evaluated all rows, while AI-REML must
    # use training-row component matrices. Refresh the kernel cache explicitly.
    with torch.no_grad():
        _ = model(train_x)
    Kg_list = getattr(kernel, "_Kg_list", None)
    Kge_list = getattr(kernel, "_Kge_list", None)
    Ke = getattr(kernel, "_Ke", None)
    have_kernel_caches = (Kg_list is not None) and (Ke is not None)
    have_ge = have_kernel_caches and (Kge_list is not None)
    use_kernel_scales = bool(refine_kernel_scales) and have_kernel_caches
    single_env = k <= 1
    optimize_ge = bool(use_kernel_scales and have_ge and not single_env)
    optimize_env_main = bool(use_kernel_scales and not single_env and not _disable_env_main)

    with torch.no_grad():
        data_var = float(torch.var(train_y.detach().to(torch.float64), unbiased=False).item()) if n_obs > 1 else 0.0
    if not math.isfinite(data_var) or data_var <= 0.0:
        data_var = 1.0
    resid_floor = max(1e-4 * data_var, 1e-10)

    # ---- Initial state from converged Adam ---------------------------------
    with torch.no_grad():
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            init_noise = likelihood.noise.detach()
        else:
            init_noise = torch.full((n_obs,), float(likelihood.noise.item()),
                                    dtype=dtype_t, device=device)
        resid0 = torch.zeros(k, dtype=torch.float64, device=device)
        for e in range(k):
            m = (ei_t == e)
            resid0[e] = (init_noise[m].mean() if m.any() else init_noise.mean()).to(torch.float64)

        if use_kernel_scales:
            wg0 = torch.nn.functional.softplus(kernel.w_g_unconstrained.detach()).to(torch.float64)
            we0 = torch.nn.functional.softplus(kernel.w_e_unconstrained.detach()).reshape(()).to(torch.float64)
            if have_ge:
                wge0 = torch.nn.functional.softplus(kernel.w_ge_unconstrained.detach()).to(torch.float64)
            else:
                wge0 = None
        else:
            wg0 = we0 = wge0 = None

    # ---- Specs --------------------------------------------------------------
    specs: List[ParamSpec] = []
    if use_kernel_scales:
        for ki, kn in enumerate(geno_kernel_names):
            specs.append(ParamSpec(kind="var_kernel",
                                   term=kernel_term_template.format(K=kn),
                                   sub_label="var", lower=0.0))
        if optimize_ge:
            for ki, kn in enumerate(geno_kernel_names):
                specs.append(ParamSpec(kind="var_kernel",
                                       term=f"{env_col_name}:{kernel_term_template.format(K=kn)}",
                                       sub_label="var", lower=0.0))
        if optimize_env_main:
            specs.append(ParamSpec(kind="var_kernel", term=env_col_name, sub_label="var", lower=0.0))
    for e in range(k):
        specs.append(ParamSpec(kind="resid_env", term=env_col_name,
                               env_label=str(env_labels[e]), lower=resid_floor))

    n_kernel_params = (
        nK + (nK if optimize_ge else 0) + (1 if optimize_env_main else 0)
    ) if use_kernel_scales else 0

    # ---- theta0 (with MINQUE-0 re-initialization if Adam weights are tiny) --
    pieces = []
    if use_kernel_scales:
        pieces.append(wg0.cpu().numpy())
        if optimize_ge:
            pieces.append(wge0.cpu().numpy())
        if optimize_env_main:
            pieces.append(np.array([float(we0.item())], dtype=np.float64))
    pieces.append(resid0.cpu().numpy())
    theta0_np = np.concatenate(pieces)

    # MINQUE-0 override: if Adam converged to near-zero kernel weights,
    # re-initialize from genotype-mean variance analysis so the AI-REML
    # starts from a basin where genetic variance is not collapsed.
    if use_kernel_scales and nK >= 1:
        wg0_np = wg0.cpu().numpy()
        mean_resid = float(resid0.mean().item())
        if np.all(wg0_np < 0.05 * mean_resid):
            # Estimate genetic variance from genotype means
            with torch.no_grad():
                y_np = train_y.cpu().to(torch.float64).numpy()
                gi_np = train_x[:, 0].long().cpu().numpy()
                ei_np = train_x[:, 1].long().cpu().numpy()
                n_geno_unique = int(gi_np.max()) + 1
                # Genotype means (across envs and reps)
                geno_sums = np.zeros(n_geno_unique, dtype=np.float64)
                geno_counts = np.zeros(n_geno_unique, dtype=np.float64)
                np.add.at(geno_sums, gi_np, y_np)
                np.add.at(geno_counts, gi_np, 1.0)
                geno_means = geno_sums / np.maximum(geno_counts, 1.0)
                mean_reps = float(np.mean(geno_counts[geno_counts > 0]))
                # Between-genotype variance
                grand_mean = float(np.mean(y_np))
                var_between = float(np.var(geno_means[geno_counts > 0], ddof=0))
                # Within-genotype variance
                resid_obs = y_np - geno_means[gi_np]
                var_within = float(np.var(resid_obs, ddof=0))
                # ANOVA h2 estimate
                h2_est = max(0.0, (var_between - var_within / mean_reps)
                             / (var_between + (mean_reps - 1.0) * var_within / mean_reps))
                h2_est = min(h2_est, 0.95)  # cap to avoid degenerate init
                total_var = var_between + var_within
                genetic_var = h2_est * total_var
                if genetic_var > 0.01 * total_var:
                    # Set each kernel weight to genetic_var / nK
                    wg_init = np.full(nK, genetic_var / nK, dtype=np.float64)
                    # Set residuals from within-genotype variance per env
                    resid_init = np.zeros(k, dtype=np.float64)
                    for e in range(k):
                        mask_e = (ei_np == e)
                        if mask_e.any():
                            r_e = y_np[mask_e] - geno_means[gi_np[mask_e]]
                            resid_init[e] = max(float(np.var(r_e, ddof=0)), 1e-6)
                        else:
                            resid_init[e] = var_within
                    # Rebuild theta0
                    pieces_new = [wg_init]
                    if optimize_ge:
                        # Keep GxE at small positive value
                        pieces_new.append(np.full(nK, 0.01 * genetic_var / nK, dtype=np.float64))
                    if optimize_env_main:
                        pieces_new.append(np.array([1e-4], dtype=np.float64))  # w_e near zero
                    pieces_new.append(resid_init)
                    theta0_np = np.concatenate(pieces_new)
                    import warnings as _wn
                    _wn.warn(
                        f"MINQUE-0 init: Adam kernel weights near zero "
                        f"(max={float(wg0_np.max()):.4f}); "
                        f"re-initialized to h2={h2_est:.3f}, "
                        f"w_g={np.round(wg_init,4)}, "
                        f"resid_mean={float(resid_init.mean()):.4f}",
                        stacklevel=2,
                    )

    # ---- State setter -------------------------------------------------------
    def _softplus_inv_torch(y):
        y_clamped = torch.clamp(y, min=1e-12)
        return torch.where(y_clamped > 20.0, y_clamped, torch.log(torch.expm1(y_clamped)))

    def _set_state(theta_np: np.ndarray):
        theta_t = torch.as_tensor(theta_np, dtype=torch.float64, device=device)
        idx = 0
        if use_kernel_scales:
            wg_t = torch.clamp(theta_t[idx:idx+nK].to(dtype_t), min=1e-10); idx += nK
            kernel.w_g_unconstrained.data.copy_(_softplus_inv_torch(wg_t).reshape(kernel.w_g_unconstrained.shape))
            if have_ge:
                if optimize_ge:
                    wge_t = torch.clamp(theta_t[idx:idx+nK].to(dtype_t), min=1e-10); idx += nK
                else:
                    wge_t = torch.full_like(kernel.w_ge_unconstrained, 1e-10)
                kernel.w_ge_unconstrained.data.copy_(_softplus_inv_torch(wge_t).reshape(kernel.w_ge_unconstrained.shape))
            if optimize_env_main:
                we_t = torch.clamp(theta_t[idx:idx+1].to(dtype_t), min=1e-10); idx += 1
            else:
                we_t = torch.full_like(kernel.w_e_unconstrained, 1e-10)
            kernel.w_e_unconstrained.data.copy_(_softplus_inv_torch(we_t).reshape(kernel.w_e_unconstrained.shape))
        resid_t = torch.clamp(theta_t[idx:idx+k].to(dtype_t), min=float(resid_floor)); idx += k
        return resid_t[ei_t]

    # ---- DV builders (matvec) ------------------------------------------------
    def _build_dv_list():
        DV: List[Callable] = []
        if use_kernel_scales:
            for ki in range(nK):
                Km = Kg_list[ki]
                DV.append(lambda v, Km=Km: Km @ v)
            if optimize_ge:
                for ki in range(nK):
                    Km = Kge_list[ki]
                    DV.append(lambda v, Km=Km: Km @ v)
            if optimize_env_main:
                DV.append(lambda v, Km=Ke: Km @ v)
        for e in range(k):
            def DVres(v, e=e):
                m = (ei_t == e).to(dtype_t)
                return v * m
            DV.append(DVres)
        return DV

    # ---- DV dense matrices (for exact score traces) -----------------------
    def _build_dv_dense():
        """Return list of dense (n,n) or diagonal (n,) derivative matrices."""
        DV_dense = []
        if use_kernel_scales:
            for ki in range(nK):
                DV_dense.append(Kg_list[ki])    # (n, n) cached kernel
            if optimize_ge:
                for ki in range(nK):
                    DV_dense.append(Kge_list[ki])  # (n, n) cached GxE kernel
            if optimize_env_main:
                DV_dense.append(Ke)                   # (n, n) cached env kernel
        for e in range(k):
            DV_dense.append((ei_t == e).to(dtype_t))  # (n,) diagonal
        return DV_dense

    def _build_V(noise_diag_obs):
        K_op = kernel(train_x, train_x)
        return K_op + DiagLinearOperator(noise_diag_obs)

    # ---- Cached numpy callables --------------------------------------------
    # Theta-only cache key + shared precompute -- see the FA backend above
    # for the rationale. This handles the multi-kernel case identically to
    # single-kernel: V_dense is built from sum_k w_k K_k + R upstream, and
    # the per-parameter DV_matvecs already encodes the per-kernel
    # derivatives.
    cache = {"key": None, "ll": None, "score": None, "ai": None, "pre": None}

    def _ensure(theta_np: np.ndarray, need_score_ai: bool):
        key = theta_np.tobytes()
        if cache["key"] == key:
            if (not need_score_ai) or (cache["score"] is not None):
                return
            with torch.no_grad():
                DV_matvecs = _build_dv_list()
                DV_dense = _build_dv_dense()
                score_np, AI_np = _exact_dense_score_and_AI_from_pre(
                    cache["pre"], DV_matvecs, DV_dense,
                )
                cache["score"] = score_np
                cache["ai"] = AI_np
            return
        with torch.no_grad():
            noise_obs = _set_state(theta_np)
            V_op = _build_V(noise_obs)
            # Deterministic dense-Cholesky path: fw.reml_loglik_from_V uses
            # CG/Lanczos with stochastic probes whose Monte Carlo noise
            # (~0.5 per call) destroys the trust-region rho ratio.
            V_dense = V_op.to_dense()
            pre = _exact_dense_precompute(
                train_y, Xtt, V_dense, jitter=float(jitter),
            )
            cache["pre"] = pre
            cache["ll"] = _exact_dense_reml_loglik_from_pre(pre)
            if need_score_ai:
                DV_matvecs = _build_dv_list()
                DV_dense = _build_dv_dense()
                score_np, AI_np = _exact_dense_score_and_AI_from_pre(
                    pre, DV_matvecs, DV_dense,
                )
                cache["score"] = score_np
                cache["ai"] = AI_np
            else:
                cache["score"] = None
                cache["ai"] = None
        cache["key"] = key

    def loglik_fn(theta_np):
        _ensure(theta_np, need_score_ai=False); return cache["ll"]
    def score_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["score"]
    def ai_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["ai"]

    # ---- Run AI-Newton ------------------------------------------------------
    import os as _os
    _verbose = bool(int(_os.environ.get("AIREML_VERBOSE", "0")))
    if _verbose:
        print(f"[gp_exact AI] theta0={np.round(theta0_np, 5).tolist()}  specs={[(s.kind, s.term, s.env_label) for s in specs]}", flush=True)
    backend = DenseVCBackend(specs=specs, loglik_fn=loglik_fn, score_fn=score_fn, ai_fn=ai_fn)

    # Warm-start: see build_fa_icm_ai_varcomp for rationale and Phase 2.2
    # shape-mismatch safety net.
    theta_start = theta0_np
    iter_budget = int(max_iter)
    if theta0_warm is not None:
        warm = np.asarray(theta0_warm, dtype=np.float64).reshape(-1)
        if warm.shape[0] != theta0_np.shape[0]:
            import warnings as _wn
            _wn.warn(
                f"theta0_warm length {warm.shape[0]} does not match expected "
                f"{theta0_np.shape[0]} for gp_exact spec layout; falling back "
                f"to cold theta0 (Adam init).",
                RuntimeWarning,
            )
        else:
            theta_start, _bounds_check = apply_boundary_projection(warm, specs)
            if max_iter_warm is not None:
                iter_budget = int(max_iter_warm)
            else:
                iter_budget = min(5, int(max_iter))
            if _verbose:
                print(f"[gp_exact AI] WARM-START enabled; iter_budget={iter_budget}", flush=True)
    result = _ai_reml_optimize(
        backend, theta_start,
        max_iter=iter_budget,
        tol_loglik=float(tol_loglik),
        tol_theta=float(tol_theta),
        damping=1.0,
        max_halvings=8,
        verbose=_verbose,
    )
    _set_state(result.theta)

    # ---- Output -------------------------------------------------------------
    pct = np.where(np.isfinite(result.pct_change), result.pct_change, 0.0)
    varcomp_df = build_asreml_varcomp_table(
        specs=specs, theta=result.theta, ai=result.ai, bounds=result.bounds,
        pct_change=pct, env_scales=env_scales,
    )
    parameters_free = int(sum(b not in ("B", "F") for b in result.bounds))
    summary = build_summary(
        loglik=float(result.loglik), nedf=int(nedf),
        parameters_free=parameters_free, n_ai_iter=int(result.n_iter),
        converged=bool(result.converged), sigma=1.0,
    )
    # Off-diagonal active-set AI inverse for non-FA models. Gives the R-side
    # h2 SE the cross-term it needs:
    #   Var(h2) = (dh2/dg)^2 var_g + (dh2/de)^2 var_e + 2*(dh2/dg)*(dh2/de)*cov_g_e
    # vs the existing 2-variance independence-assuming formula.
    from varcomp_asreml import compute_genetic_residual_covariance as _compute_g_r_cov
    g_r_cov = _compute_g_r_cov(
        specs=specs, theta=result.theta, ai=result.ai,
        bounds=result.bounds, env_scales=env_scales,
    )
    if g_r_cov is not None:
        summary = dict(summary)
        summary["genetic_residual_covariance"] = g_r_cov
    # Surface the converged variance-component theta and its layout for the
    # R-side warm-start cache (Phase 2.1 Stage B). The gp_exact spec layout
    # is dynamic depending on optimize flags; layout key encodes the
    # backend so the R-side cache only feeds it back into a matching fit.
    summary = dict(summary)
    summary["theta"] = result.theta.tolist()
    summary["theta_layout"] = "gp_exact"
    summary["theta_shape"] = {"n_params": int(len(result.theta))}
    return varcomp_df, summary, np.asarray(result.ai, dtype=np.float64)


def build_gp_exact_ai_varcomp_mme(
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
    env_scales: Optional["EnvScales"] = None,
    nedf: int = 0,
    max_iter: int = 12,
    tol_loglik: float = 1e-3,
    tol_theta: float = 1e-4,
    theta0_warm: Optional[np.ndarray] = None,
    max_iter_warm: Optional[int] = None,
    refine_kernel_scales: bool = True,
    mme_backend: str = "auto",
    mme_n_probes: int = 128,
):
    """MME-formulation analogue of build_gp_exact_ai_varcomp.

    Routes through the sparse MME engine in mme_reml.py for the single-trait,
    single- or multi-kernel CS-MET case. Falls back via NotImplementedError when
    GxE-kernel, env-main, or FA parameterizations are requested.
    Returns ``(varcomp_df, summary_dict, ai_matrix_np)``;
    ``summary["theta_layout"] = "gp_exact_mme"``.
    """
    import torch
    import gpytorch
    import scipy.sparse as sp
    import os as _os
    _tll_env = _os.environ.get("PREDICTPRO_GP_TOL_LOGLIK", "").strip()
    _tth_env = _os.environ.get("PREDICTPRO_GP_TOL_THETA", "").strip()
    _mit_env = _os.environ.get("PREDICTPRO_GP_MAX_ITER", "").strip()
    _dem_env = _os.environ.get("PREDICTPRO_GP_DISABLE_ENV_MAIN", "").strip().lower()
    if _tll_env:
        try: tol_loglik = float(_tll_env)
        except ValueError: pass
    if _tth_env:
        try: tol_theta = float(_tth_env)
        except ValueError: pass
    if _mit_env:
        try: max_iter = int(_mit_env)
        except ValueError: pass
    _disable_env_main = _dem_env in ("1", "true", "yes", "on")

    device = train_y.device
    dtype_t = train_y.dtype
    ei_t = train_x[:, 1].long()
    n_obs = int(train_y.shape[0])
    Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
    kernel = model.covar_module
    k = len(env_labels)
    nK = len(geno_kernel_names)

    with torch.no_grad():
        _ = model(train_x)
    Kg_list = getattr(kernel, "_Kg_list", None)
    Kge_list = getattr(kernel, "_Kge_list", None)
    Ke = getattr(kernel, "_Ke", None)
    have_kernel_caches = (Kg_list is not None) and (Ke is not None)
    have_ge = have_kernel_caches and (Kge_list is not None)
    use_kernel_scales = bool(refine_kernel_scales) and have_kernel_caches
    single_env = k <= 1
    optimize_ge = bool(use_kernel_scales and have_ge and not single_env)
    optimize_env_main = bool(use_kernel_scales and not single_env and not _disable_env_main)

    if optimize_ge or optimize_env_main:
        raise NotImplementedError(
            "build_gp_exact_ai_varcomp_mme: GxE and env-main kernels are not "
            "supported in Phase 3.4. Use gp_engine='dense_v' or wait for 3.5."
        )

    with torch.no_grad():
        data_var = float(torch.var(train_y.detach().to(torch.float64), unbiased=False).item()) if n_obs > 1 else 0.0
    if not math.isfinite(data_var) or data_var <= 0.0:
        data_var = 1.0
    resid_floor = max(1e-4 * data_var, 1e-10)

    with torch.no_grad():
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            init_noise = likelihood.noise.detach()
        else:
            init_noise = torch.full((n_obs,), float(likelihood.noise.item()),
                                    dtype=dtype_t, device=device)
        resid0 = torch.zeros(k, dtype=torch.float64, device=device)
        for e in range(k):
            m = (ei_t == e)
            resid0[e] = (init_noise[m].mean() if m.any() else init_noise.mean()).to(torch.float64)

        if getattr(kernel, "w_g_unconstrained", None) is not None:
            wg0 = torch.nn.functional.softplus(kernel.w_g_unconstrained.detach()).to(torch.float64)
        else:
            wg0 = kernel.w_g_fix.detach().to(torch.float64)

    specs: List[ParamSpec] = []
    for ki, kn in enumerate(geno_kernel_names):
        specs.append(ParamSpec(kind="var_kernel",
                               term=kernel_term_template.format(K=kn),
                               sub_label="var", lower=0.0))
    for e in range(k):
        specs.append(ParamSpec(kind="resid_env", term=env_col_name,
                               env_label=str(env_labels[e]), lower=resid_floor))

    theta0_np = np.concatenate([wg0.cpu().numpy().reshape(-1), resid0.cpu().numpy()])

    gi_np = train_x[:, 0].long().cpu().numpy()
    ei_np = train_x[:, 1].long().cpu().numpy()
    train_y_np = train_y.detach().cpu().to(torch.float64).numpy()
    kernel_G_list = getattr(kernel, "G_list", None)
    if kernel_G_list is not None:
        n_geno = int(kernel_G_list[0].shape[0])
    else:
        n_geno = int(gi_np.max()) + 1 if gi_np.size else 0

    if nK >= 1:
        wg0_np = wg0.cpu().numpy().reshape(-1)
        mean_resid = float(resid0.mean().item())
        if np.all(wg0_np < 0.05 * mean_resid):
            y_np = train_y_np
            n_geno_unique = int(gi_np.max()) + 1
            geno_sums = np.zeros(n_geno_unique, dtype=np.float64)
            geno_counts = np.zeros(n_geno_unique, dtype=np.float64)
            np.add.at(geno_sums, gi_np, y_np)
            np.add.at(geno_counts, gi_np, 1.0)
            geno_means = geno_sums / np.maximum(geno_counts, 1.0)
            mean_reps = float(np.mean(geno_counts[geno_counts > 0]))
            var_between = float(np.var(geno_means[geno_counts > 0], ddof=0))
            resid_obs = y_np - geno_means[gi_np]
            var_within = float(np.var(resid_obs, ddof=0))
            h2_est = max(0.0, (var_between - var_within / mean_reps)
                         / (var_between + (mean_reps - 1.0) * var_within / mean_reps))
            h2_est = min(h2_est, 0.95)
            total_var = var_between + var_within
            genetic_var = h2_est * total_var
            if genetic_var > 0.01 * total_var:
                wg_init = np.full(nK, genetic_var / nK, dtype=np.float64)
                resid_init = np.zeros(k, dtype=np.float64)
                for e in range(k):
                    mask_e = (ei_np == e)
                    if mask_e.any():
                        r_e = y_np[mask_e] - geno_means[gi_np[mask_e]]
                        resid_init[e] = max(float(np.var(r_e, ddof=0)), 1e-6)
                    else:
                        resid_init[e] = var_within
                pieces_new = [wg_init, resid_init]
                theta0_np = np.concatenate(pieces_new)
                import warnings as _wn
                _wn.warn(
                    f"MINQUE-0 init (MME): Adam kernel weights near zero "
                    f"(max={float(wg0_np.max()):.4f}); "
                    f"re-initialized to h2={h2_est:.3f}, "
                    f"w_g={np.round(wg_init,4)}, "
                    f"resid_mean={float(resid_init.mean()):.4f}",
                    stacklevel=2,
                )

    G_list_np: List[np.ndarray] = []
    if kernel_G_list is None:
        raise RuntimeError(
            "build_gp_exact_ai_varcomp_mme: kernel.G_list not found; "
            "expected genotype-level per-kernel SPD matrices."
        )
    gi_max = int(gi_np.max()) if gi_np.size else -1
    if gi_max >= n_geno:
        raise RuntimeError(
            f"build_gp_exact_ai_varcomp_mme: gi index {gi_max} out of "
            f"range for kernel.G_list with n_geno={n_geno}."
        )
    for Gk in kernel_G_list[:nK]:
        G_geno = Gk.detach().cpu().to(torch.float64).numpy()
        G_list_np.append(np.ascontiguousarray(G_geno))

    G_inv_list = [sp.csc_matrix(np.linalg.inv(Gk)) for Gk in G_list_np]
    G_logdet_list = [float(np.linalg.slogdet(Gk)[1]) for Gk in G_list_np]

    Z_sp = sp.csr_matrix(
        (np.ones(n_obs, dtype=np.float64), (np.arange(n_obs), gi_np)),
        shape=(n_obs, n_geno),
    )

    X_np = Xtt.detach().cpu().to(torch.float64).numpy()
    if X_np.ndim == 1:
        X_np = X_np.reshape(-1, 1)

    n_obs_per_env = np.bincount(ei_np, minlength=int(k)).astype(np.float64)

    spec_kinds: List[str] = ["var_kernel"] * nK + ["resid_env"] * k
    spec_payloads: List[np.ndarray] = [np.asarray(Gk) for Gk in G_list_np] + \
                                     [(ei_np == e).astype(np.float64) for e in range(k)]

    def _decode_theta(theta_np: np.ndarray):
        sigma2_g_list = [float(v) for v in theta_np[:nK]]
        sigma2_e_per_env = np.asarray(theta_np[nK:nK + k], dtype=np.float64)
        return sigma2_g_list, sigma2_e_per_env

    cache: Dict[str, Any] = {"key": None, "ll": None, "score": None, "ai": None, "asm": None}

    def _ensure(theta_np: np.ndarray, need_score_ai: bool):
        key = theta_np.tobytes()
        if cache["key"] == key:
            if (not need_score_ai) or (cache["score"] is not None):
                return
            score_np, AI_np = ai_and_score(
                asm=cache["asm"], X=X_np, Z=Z_sp,
                spec_kinds=spec_kinds, spec_payloads=spec_payloads,
                K_kernels=nK, n_probes=int(mme_n_probes), seed=int(seed),
            )
            cache["score"] = score_np
            cache["ai"] = AI_np
            return
        sigma2_g_list, sigma2_e_per_env = _decode_theta(theta_np)
        sigma2_g_safe = [max(v, 1e-12) for v in sigma2_g_list]
        sigma2_e_safe = np.clip(sigma2_e_per_env, resid_floor, None)
        asm = mme_precompute(
            X=X_np, Z=Z_sp, G_list=G_list_np, sigma2_g_list=sigma2_g_safe,
            sigma2_e_per_env=sigma2_e_safe, ei=ei_np, y=train_y_np,
            G_inv_list=G_inv_list, backend=str(mme_backend),
        )
        cache["asm"] = asm
        cache["ll"] = reml_loglik_from_assembly(
            asm,
            sigma2_e_per_env=sigma2_e_safe,
            n_obs_per_env=n_obs_per_env,
            sigma2_g_list=sigma2_g_safe,
            G_logdet_list=G_logdet_list,
            q=int(n_geno),
            n_minus_p=int(n_obs - X_np.shape[1]),
        )
        if need_score_ai:
            score_np, AI_np = ai_and_score(
                asm=asm, X=X_np, Z=Z_sp,
                spec_kinds=spec_kinds, spec_payloads=spec_payloads,
                K_kernels=nK, n_probes=int(mme_n_probes), seed=int(seed),
            )
            cache["score"] = score_np
            cache["ai"] = AI_np
        else:
            cache["score"] = None
            cache["ai"] = None
        cache["key"] = key

    def loglik_fn(theta_np):
        _ensure(theta_np, need_score_ai=False); return cache["ll"]
    def score_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["score"]
    def ai_fn(theta_np):
        _ensure(theta_np, need_score_ai=True); return cache["ai"]

    import os as _os
    _verbose = bool(int(_os.environ.get("AIREML_VERBOSE", "0")))
    if _verbose:
        print(f"[gp_exact_mme AI] theta0={np.round(theta0_np, 5).tolist()}  specs={[(s.kind, s.term, s.env_label) for s in specs]}  backend={select_backend(str(mme_backend))}", flush=True)
    backend = DenseVCBackend(specs=specs, loglik_fn=loglik_fn, score_fn=score_fn, ai_fn=ai_fn)

    theta_start = theta0_np
    iter_budget = int(max_iter)
    if theta0_warm is not None:
        warm = np.asarray(theta0_warm, dtype=np.float64).reshape(-1)
        if warm.shape[0] != theta0_np.shape[0]:
            import warnings as _wn
            _wn.warn(
                f"theta0_warm length {warm.shape[0]} does not match expected "
                f"{theta0_np.shape[0]} for gp_exact_mme spec layout; falling "
                f"back to cold theta0 (Adam init).",
                RuntimeWarning,
            )
        else:
            theta_start, _bounds_check = apply_boundary_projection(warm, specs)
            if max_iter_warm is not None:
                iter_budget = int(max_iter_warm)
            else:
                iter_budget = min(5, int(max_iter))
            if _verbose:
                print(f"[gp_exact_mme AI] WARM-START enabled; iter_budget={iter_budget}", flush=True)

    result = _ai_reml_optimize(
        backend, theta_start,
        max_iter=iter_budget,
        tol_loglik=float(tol_loglik),
        tol_theta=float(tol_theta),
        damping=1.0,
        max_halvings=8,
        verbose=_verbose,
    )

    pct = np.where(np.isfinite(result.pct_change), result.pct_change, 0.0)
    varcomp_df = build_asreml_varcomp_table(
        specs=specs, theta=result.theta, ai=result.ai, bounds=result.bounds,
        pct_change=pct, env_scales=env_scales,
    )
    parameters_free = int(sum(b not in ("B", "F") for b in result.bounds))
    summary = build_summary(
        loglik=float(result.loglik), nedf=int(nedf),
        parameters_free=parameters_free, n_ai_iter=int(result.n_iter),
        converged=bool(result.converged), sigma=1.0,
    )
    from varcomp_asreml import compute_genetic_residual_covariance as _compute_g_r_cov
    g_r_cov = _compute_g_r_cov(
        specs=specs, theta=result.theta, ai=result.ai,
        bounds=result.bounds, env_scales=env_scales,
    )
    if g_r_cov is not None:
        summary = dict(summary)
        summary["genetic_residual_covariance"] = g_r_cov
    summary = dict(summary)
    summary["theta"] = result.theta.tolist()
    summary["theta_layout"] = "gp_exact_mme"
    summary["theta_shape"] = {"n_params": int(len(result.theta))}
    summary["mme_backend"] = select_backend(str(mme_backend))
    return varcomp_df, summary, np.asarray(result.ai, dtype=np.float64)


__all__ = [
    "DenseAIContext",
    "ai_covariance_with_bounds",
    "run_ai_reml",
    "build_varcomp_outputs",
    "make_fa_rank1_dv_builders",
    "_exact_dense_score_and_AI",
    "compute_ai_se_active_set",
    "bounds_from_specs_and_theta",
    "build_fa_icm_ai_varcomp",
    "build_gp_exact_ai_varcomp",
    "build_gp_exact_ai_varcomp_mme",
]
