"""Internal math + helpers for fit_multi_trait_gp.

Not a public API. Imported only by mt_gp.py.
"""
from __future__ import annotations
import torch
import numpy as np
import pandas as pd
from typing import Any, Dict, Sequence, Optional

from gp_device import pick_torch_device, torch_to_numpy


def _sanitize_kernel_eigenvalues(
    eigenvalues: "torch.Tensor",
    relative_tolerance: float = 1e-8,
) -> tuple:
    """Clip numerical negative GRM eigenvalues and reject invalid kernels.

    Relationship matrices reaching this backend have already passed the R
    kernel guardrail.  A dense eigendecomposition can nevertheless return
    tiny negative values from roundoff.  Those values are projected to zero;
    a negative eigenvalue beyond the relative tolerance remains a hard input
    error rather than being hidden by covariance jitter.
    """
    if eigenvalues.numel() == 0:
        return eigenvalues, 0, 0.0
    detached = eigenvalues.detach()
    if not bool(torch.isfinite(detached).all().item()):
        raise ValueError("The genomic relationship matrix has non-finite eigenvalues.")
    scale = max(float(torch.max(torch.abs(detached)).item()), 1.0)
    min_eigenvalue = float(torch.min(detached).item())
    tolerance = abs(float(relative_tolerance)) * scale
    if min_eigenvalue < -tolerance:
        raise ValueError(
            "The genomic relationship matrix is not positive semidefinite after "
            f"input validation: minimum eigenvalue={min_eigenvalue:.6e}, "
            f"allowed numerical tolerance={tolerance:.6e}."
        )
    clipped_count = int(torch.sum(detached < 0).item())
    return torch.clamp(eigenvalues, min=0.0), clipped_count, min_eigenvalue


def _is_recoverable_cholesky_error(error: RuntimeError) -> bool:
    """Whether an error is a local non-PD trial-point factorization failure."""
    message = str(error).lower()
    if any(token in message for token in ("cuda", "illegal memory", "device-side")):
        return False
    return "cholesky" in message and any(
        token in message
        for token in ("positive-definite", "positive definite", "not positive")
    )


def _objective_value_or_none(objective, theta: "torch.Tensor") -> Optional[float]:
    """Evaluate a line-search point, rejecting only numerical Cholesky failures."""
    try:
        value = objective(theta)
    except RuntimeError as error:
        if not _is_recoverable_cholesky_error(error):
            raise
        return None
    detached = value.detach()
    if detached.numel() != 1:
        raise RuntimeError("AI-REML objective must return one scalar value.")
    if not bool(torch.isfinite(detached).item()):
        return None
    return float(detached.item())


def _is_balanced_grid(gid_idx, trait_idx, n_gid: int, T: int) -> bool:
    """True iff every (gid, trait) pair appears exactly once.

    Accepts numpy arrays or torch tensors. Used inside `_neg_reml_loglik` to
    decide whether to dispatch to the Kronecker fast path or the dense V_obs
    fallback. The check is `gid_idx.shape[0] == n_gid * T` — sufficient because
    `_build_mt_design` already drops NaN-y rows and disallows duplicate
    (gid, trait) pairs.
    """
    return int(gid_idx.shape[0]) == int(n_gid) * int(T)


def _build_mt_design(
    *,
    pheno_df: pd.DataFrame,
    gid_col: str,
    trait_col: str,
    y_col: str,
    fixed_effects: Sequence[str],
    geno_ids: Sequence[str],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> Dict[str, Any]:
    """Build long-format MT-GP design.

    Returns dict with keys:
      T, n_gid, gid_idx, trait_idx, y_std, X_fixed, train_mask, test_mask,
      trait_levels, trait_means, trait_sds, geno_ids.
    """
    df = pheno_df.copy()
    n_orig = len(df)
    # Filter to rows we care about (train ∪ test)
    keep = np.zeros(n_orig, dtype=bool)
    keep[np.asarray(train_idx, dtype=np.int64)] = True
    keep[np.asarray(test_idx, dtype=np.int64)] = True
    df = df.loc[keep].copy()
    for col in (gid_col, trait_col):
        df[col] = df[col].astype("object").map(lambda x: str(x))
    # Mark which rows are train vs test in the kept set
    is_train_orig = np.zeros(n_orig, dtype=bool)
    is_train_orig[np.asarray(train_idx, dtype=np.int64)] = True
    df["_is_train"] = is_train_orig[keep]
    # Drop missing responses only from training rows. Test rows often carry
    # NaN y by design; they still need design rows so predictions are emitted.
    df = df.loc[(~df["_is_train"]) | df[y_col].notna()].reset_index(drop=True)

    # Trait levels (sorted alphabetically for determinism)
    trait_levels = sorted(df[trait_col].unique().tolist())
    T = len(trait_levels)
    if T < 2:
        raise ValueError(f"MT-GP requires T >= 2 traits; got {T} ({trait_levels!r})")
    trait_to_idx = {t: i for i, t in enumerate(trait_levels)}

    # gid mapping
    gid_to_idx = {str(g): i for i, g in enumerate(geno_ids)}
    n_gid = len(geno_ids)
    df["_trait_idx"] = df[trait_col].map(trait_to_idx).astype(np.int64)
    if not df[gid_col].map(lambda g: g in gid_to_idx).all():
        missing = df.loc[~df[gid_col].isin(gid_to_idx), gid_col].unique()
        raise ValueError(f"gids in pheno_df not found in geno_ids: {missing.tolist()[:5]}")
    df["_gid_idx"] = df[gid_col].map(gid_to_idx).astype(np.int64)

    # Per-trait standardization on training rows
    trait_means = np.zeros(T, dtype=np.float64)
    trait_sds = np.ones(T, dtype=np.float64)
    train_only = df.loc[df["_is_train"]]
    for t, name in enumerate(trait_levels):
        sl = train_only.loc[train_only["_trait_idx"] == t, y_col].to_numpy(dtype=np.float64)
        if sl.size < 2:
            raise ValueError(f"trait {name!r} has <2 training rows; cannot standardize")
        trait_means[t] = float(np.mean(sl))
        sd = float(np.std(sl, ddof=1))
        trait_sds[t] = sd if sd > 1e-12 else 1.0
    y_raw = df[y_col].to_numpy(dtype=np.float64)
    y_std = (y_raw - trait_means[df["_trait_idx"].to_numpy()]
             ) / trait_sds[df["_trait_idx"].to_numpy()]

    # Trait-specific fixed-effect expansion
    # Always start with T intercept columns (one-hot trait indicator)
    trait_idx = df["_trait_idx"].to_numpy(dtype=np.int64)
    n_rows = len(df)
    blocks = [np.eye(T)[trait_idx]]  # (n_rows, T) intercepts
    for col in fixed_effects:
        if col not in df.columns:
            raise ValueError(f"fixed_effects col {col!r} not in pheno_df")
        levels = sorted(df[col].astype(str).unique().tolist())
        if len(levels) < 2:
            continue  # constant column contributes nothing
        # One-hot of the column (no drop_first for fixed effects across traits)
        col_oh = pd.get_dummies(df[col].astype(str), drop_first=False).to_numpy(dtype=np.float64)
        # Trait-specific expansion: each col_oh column × T traits = T expanded cols
        for c in range(col_oh.shape[1]):
            block = np.zeros((n_rows, T), dtype=np.float64)
            for t in range(T):
                m = trait_idx == t
                block[m, t] = col_oh[m, c]
            blocks.append(block)
    X_fixed = np.concatenate(blocks, axis=1)

    train_mask = df["_is_train"].to_numpy(dtype=bool)
    test_mask = ~train_mask

    return {
        "T": T,
        "n_gid": n_gid,
        "gid_idx": df["_gid_idx"].to_numpy(dtype=np.int64),
        "trait_idx": trait_idx,
        "y_std": y_std,
        "X_fixed": X_fixed,
        "train_mask": train_mask,
        "test_mask": test_mask,
        "trait_levels": trait_levels,
        "trait_means": trait_means,
        "trait_sds": trait_sds,
        "geno_ids": list(geno_ids),
    }


def _initialize_params(
    *,
    y_std: np.ndarray,
    X: np.ndarray,
    gid_idx: np.ndarray,
    trait_idx: np.ndarray,
    train_mask: np.ndarray,
    T: int,
    n_gid: int,
    h2_prior: float = 0.5,
) -> tuple:
    """OLS-residual seed for Sigma_eps; diagonal Sigma_G scaled by h2 prior.

    Returns (Sigma_G_init, Sigma_eps_init) as numpy arrays.
    """
    y_tr = y_std[train_mask]
    X_tr = X[train_mask]
    t_tr = trait_idx[train_mask]
    # OLS residuals
    beta_ols, *_ = np.linalg.lstsq(X_tr, y_tr, rcond=None)
    resid = y_tr - X_tr @ beta_ols
    # Per-trait residual variance, organized as a Sigma_eps estimate via outer product of
    # residuals stacked per gid where both traits are present (for off-diagonals).
    Sigma_eps0 = np.zeros((T, T), dtype=np.float64)
    counts = np.zeros((T, T), dtype=np.int64)
    g_tr = gid_idx[train_mask]
    by_gid: Dict[int, np.ndarray] = {}
    for r_idx in range(len(y_tr)):
        gi = int(g_tr[r_idx])
        ti = int(t_tr[r_idx])
        if gi not in by_gid:
            by_gid[gi] = np.full(T, np.nan, dtype=np.float64)
        by_gid[gi][ti] = resid[r_idx]
    for vec in by_gid.values():
        for a in range(T):
            for b in range(T):
                if np.isfinite(vec[a]) and np.isfinite(vec[b]):
                    Sigma_eps0[a, b] += vec[a] * vec[b]
                    counts[a, b] += 1
    counts_clip = np.maximum(counts, 1)
    Sigma_eps0 = Sigma_eps0 / counts_clip
    # PSD-project: clip eigenvalues at small positive
    w, V = np.linalg.eigh(0.5 * (Sigma_eps0 + Sigma_eps0.T))
    w = np.maximum(w, 1e-4)
    Sigma_eps0 = V @ np.diag(w) @ V.T
    # Sigma_G: diagonal seed proportional to h2 prior, scaled per-trait by total resid var
    diag_var = np.diag(Sigma_eps0).copy()
    total = diag_var / max(1.0 - h2_prior, 1e-3)
    Sigma_G0 = np.diag(total * h2_prior)
    Sigma_G0 = Sigma_G0 + 1e-4 * np.eye(T)
    return Sigma_G0, Sigma_eps0


def _fit_mt_ai_reml(
    *,
    K_geno: np.ndarray,
    y_std: np.ndarray,
    X: np.ndarray,
    gid_idx: np.ndarray,
    trait_idx: np.ndarray,
    train_mask: np.ndarray,
    T: int,
    n_gid: int,
    Sigma_G_init: np.ndarray,
    Sigma_eps_init: np.ndarray,
    trait_structure: str = "unstructured",
    trait_fa_rank: Optional[int] = None,
    residual_structure: str = "unstructured",
    residual_fa_rank: Optional[int] = None,
    max_iter: int = 100,
    tol_loglik: float = 1e-4,
    dtype: str = "float64",
    device: Optional[Any] = None,
) -> Dict[str, Any]:
    """Joint AI-REML fit over (Sigma_G, Sigma_eps) using torch.func score + Hessian Newton steps.

    Returns dict with keys: Sigma_G, Sigma_eps, Lambda_G, psi_G, Lambda_eps, psi_eps,
    beta, loglik_history, n_iter, converged.
    """
    torch_dtype = torch.float64 if dtype == "float64" else torch.float32
    device = pick_torch_device(device)
    K_t = torch.tensor(K_geno, dtype=torch_dtype, device=device)
    K_t = 0.5 * (K_t + K_t.T)
    D_K, U_K = torch.linalg.eigh(K_t)
    D_K, kernel_eigen_clip_count, kernel_min_eigenvalue_input = (
        _sanitize_kernel_eigenvalues(D_K)
    )
    y_t = torch.tensor(y_std, dtype=torch_dtype, device=device)
    X_t = torch.tensor(X, dtype=torch_dtype, device=device)
    g_t = torch.tensor(gid_idx, dtype=torch.int64, device=device)
    t_t = torch.tensor(trait_idx, dtype=torch.int64, device=device)
    tm_t = torch.tensor(train_mask, dtype=torch.bool, device=device)

    # --- parameter packing ---
    if trait_structure == "unstructured":
        L_G_init = torch.linalg.cholesky(
            torch.tensor(Sigma_G_init, dtype=torch_dtype, device=device) +
            1e-8 * torch.eye(T, dtype=torch_dtype, device=device)
        )
        theta_G_init = _pack_unstructured(L_G_init)
        n_G = T * (T + 1) // 2
        unpack_G = lambda th: (lambda L: L @ L.T)(_unpack_unstructured(th, T))
    elif trait_structure == "fa":
        if trait_fa_rank is None:
            raise ValueError("trait_fa_rank required when trait_structure='fa'")
        kG = int(trait_fa_rank)
        w, V = torch.linalg.eigh(torch.tensor(Sigma_G_init, dtype=torch_dtype, device=device))
        Lambda0 = (V[:, -kG:] * torch.sqrt(torch.clamp(w[-kG:], min=1e-6))).contiguous()
        for i in range(kG):
            for j in range(i + 1, kG):
                Lambda0[i, j] = 0.0
        psi0 = torch.diagonal(torch.tensor(Sigma_G_init, dtype=torch_dtype, device=device)) - (Lambda0 ** 2).sum(dim=1)
        psi0 = torch.clamp(psi0, min=1e-4)
        log_psi0 = torch.log(psi0)
        theta_G_init = _pack_fa(Lambda0, log_psi0)
        n_G = T * kG - kG * (kG - 1) // 2 + T
        unpack_G = lambda th: (lambda L_psi: L_psi[0] @ L_psi[0].T + torch.diag(torch.exp(L_psi[1])))(_unpack_fa(th, T, kG))
    else:
        raise ValueError(f"unknown trait_structure {trait_structure!r}")

    if residual_structure == "unstructured":
        L_E_init = torch.linalg.cholesky(
            torch.tensor(Sigma_eps_init, dtype=torch_dtype, device=device) +
            1e-8 * torch.eye(T, dtype=torch_dtype, device=device)
        )
        theta_E_init = _pack_unstructured(L_E_init)
        n_E = T * (T + 1) // 2
        unpack_E = lambda th: (lambda L: L @ L.T)(_unpack_unstructured(th, T))
    elif residual_structure == "fa":
        if residual_fa_rank is None:
            raise ValueError("residual_fa_rank required when residual_structure='fa'")
        kE = int(residual_fa_rank)
        w, V = torch.linalg.eigh(torch.tensor(Sigma_eps_init, dtype=torch_dtype, device=device))
        Lambda0 = (V[:, -kE:] * torch.sqrt(torch.clamp(w[-kE:], min=1e-6))).contiguous()
        for i in range(kE):
            for j in range(i + 1, kE):
                Lambda0[i, j] = 0.0
        psi0 = torch.diagonal(torch.tensor(Sigma_eps_init, dtype=torch_dtype, device=device)) - (Lambda0 ** 2).sum(dim=1)
        psi0 = torch.clamp(psi0, min=1e-4)
        log_psi0 = torch.log(psi0)
        theta_E_init = _pack_fa(Lambda0, log_psi0)
        n_E = T * kE - kE * (kE - 1) // 2 + T
        unpack_E = lambda th: (lambda L_psi: L_psi[0] @ L_psi[0].T + torch.diag(torch.exp(L_psi[1])))(_unpack_fa(th, T, kE))
    else:
        raise ValueError(f"unknown residual_structure {residual_structure!r}")

    theta = torch.cat([theta_G_init, theta_E_init]).detach().clone().to(device=device, dtype=torch_dtype)

    def nll_of_theta(theta_v: "torch.Tensor") -> "torch.Tensor":
        Sigma_G = unpack_G(theta_v[:n_G])
        Sigma_eps = unpack_E(theta_v[n_G:])
        return _neg_reml_loglik(
            Sigma_G=Sigma_G, Sigma_eps=Sigma_eps,
            y_std=y_t, X=X_t, gid_idx=g_t, trait_idx=t_t,
            train_mask=tm_t, D=D_K, U=U_K, T=T, n_gid=n_gid,
        )

    loglik_history: list = []
    converged = False
    singular_streak = 0
    line_search_invalid_candidate_rejections = 0
    for it in range(max_iter):
        nll_curr = _objective_value_or_none(nll_of_theta, theta)
        if nll_curr is None:
            raise RuntimeError(
                "AI-REML reached an invalid accepted covariance state before "
                f"iteration {it}; refusing to continue from a non-finite or "
                "non-positive-definite likelihood point."
            )
        loglik_history.append(-nll_curr)
        score = torch.func.grad(nll_of_theta)(theta)
        try:
            H = torch.func.hessian(nll_of_theta)(theta)
            # Levenberg-Marquardt style: ensure H_reg is positive definite before solving.
            # If Hessian has negative eigenvalues, add enough diagonal to make it PD.
            eig_min = torch.linalg.eigvalsh(H).min().item()
            lam = max(1e-6, -eig_min + 1e-3) if eig_min < 1e-3 else 1e-6
            H_reg = H + lam * torch.eye(theta.shape[0], dtype=torch_dtype, device=device)
            step = torch.linalg.solve(H_reg, score)
            singular_streak = 0
        except Exception:
            singular_streak += 1
            if singular_streak >= 2:
                raise RuntimeError(
                    f"AI-REML Hessian singular two iterations in a row at iter {it}. "
                    f"Sigma_G eigvals={np.linalg.eigvalsh(torch_to_numpy(unpack_G(theta[:n_G])))}, "
                    f"Sigma_eps eigvals={np.linalg.eigvalsh(torch_to_numpy(unpack_E(theta[n_G:])))}, "
                    f"loglik_history={loglik_history!r}"
                )
            step = score * 0.1
        # Backtracking line search: try shrinking the step until nll decreases.
        # If no shrink succeeds, skip the update (stay at current theta) to maintain monotonicity.
        step_accepted = False
        for shrink in (1.0, 0.5, 0.25, 0.1, 0.05, 0.01):
            theta_new = theta - shrink * step
            nll_try = _objective_value_or_none(nll_of_theta, theta_new)
            if nll_try is None:
                line_search_invalid_candidate_rejections += 1
                continue
            if nll_try < nll_curr - 1e-12:
                theta = theta_new
                step_accepted = True
                break
        if not step_accepted:
            # Line search exhausted: try a tiny gradient step
            for grad_shrink in (1e-2, 1e-3, 1e-4):
                theta_new = theta - grad_shrink * score
                nll_try = _objective_value_or_none(nll_of_theta, theta_new)
                if nll_try is None:
                    line_search_invalid_candidate_rejections += 1
                    continue
                if nll_try < nll_curr - 1e-12:
                    theta = theta_new
                    step_accepted = True
                    break
        # Convergence: only when we made a step AND log-lik change is small
        if it > 0 and step_accepted and abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik:
            converged = True
            break
        # Hard stop if 3 iterations in a row failed to improve
        if not step_accepted:
            singular_streak += 1
            if singular_streak >= 3:
                converged = False
                break

    Sigma_G_final = torch_to_numpy(unpack_G(theta[:n_G]))
    Sigma_eps_final = torch_to_numpy(unpack_E(theta[n_G:]))
    Sigma_G_t = unpack_G(theta[:n_G])
    Sigma_eps_t = unpack_E(theta[n_G:])
    y_tr = y_t[tm_t]; X_tr = X_t[tm_t]; g_tr = g_t[tm_t]; t_tr = t_t[tm_t]
    K_pairs = K_t[g_tr.unsqueeze(1), g_tr.unsqueeze(0)]
    G_pairs = Sigma_G_t[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]
    E_pairs = Sigma_eps_t[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]
    same_gid = (g_tr.unsqueeze(1) == g_tr.unsqueeze(0)).to(K_pairs.dtype)
    V = G_pairs * K_pairs + E_pairs * same_gid
    V = V + 1e-6 * torch.eye(y_tr.shape[0], dtype=V.dtype, device=V.device)
    L_V = torch.linalg.cholesky(V)
    Vinv_y = torch.cholesky_solve(y_tr.unsqueeze(1), L_V).squeeze(1)
    Vinv_X = torch.cholesky_solve(X_tr, L_V)
    Xt_Vinv_X = X_tr.T @ Vinv_X + 1e-6 * torch.eye(X_tr.shape[1], dtype=V.dtype, device=V.device)
    L_xvx = torch.linalg.cholesky(Xt_Vinv_X)
    beta_hat = torch.cholesky_solve(
        (X_tr.T @ Vinv_y).unsqueeze(1), L_xvx
    ).squeeze(1)

    out: Dict[str, Any] = {
        "Sigma_G": Sigma_G_final,
        "Sigma_eps": Sigma_eps_final,
        "Lambda_G": None, "psi_G": None,
        "Lambda_eps": None, "psi_eps": None,
        "beta": torch_to_numpy(beta_hat),
        "loglik_history": np.asarray(loglik_history, dtype=np.float64),
        "n_iter": len(loglik_history),
        "converged": bool(converged),
        "_K_geno": K_geno,
        "_D_K": torch_to_numpy(D_K),
        "_U_K": torch_to_numpy(U_K),
        "device": str(device),
        "kernel_eigen_clip_count": kernel_eigen_clip_count,
        "kernel_min_eigenvalue_input": kernel_min_eigenvalue_input,
        "line_search_invalid_candidate_rejections": line_search_invalid_candidate_rejections,
    }
    if trait_structure == "fa":
        Lambda_G, log_psi_G = _unpack_fa(theta[:n_G], T, int(trait_fa_rank))
        out["Lambda_G"] = torch_to_numpy(Lambda_G)
        out["psi_G"] = torch_to_numpy(torch.exp(log_psi_G))
    if residual_structure == "fa":
        Lambda_E, log_psi_E = _unpack_fa(theta[n_G:], T, int(residual_fa_rank))
        out["Lambda_eps"] = torch_to_numpy(Lambda_E)
        out["psi_eps"] = torch_to_numpy(torch.exp(log_psi_E))
    return out


def _combine_kernel_bank(kernel_bank: Sequence[np.ndarray],
                         weights: np.ndarray) -> np.ndarray:
    """K(w) = sum_k w_k K_k, symmetrised."""
    K = np.zeros_like(np.asarray(kernel_bank[0], dtype=np.float64))
    for K_k, w_k in zip(kernel_bank, np.asarray(weights, dtype=np.float64).tolist()):
        one = np.asarray(K_k, dtype=np.float64)
        if one.shape != K.shape:
            raise ValueError(
                f"kernel bank shape mismatch: {one.shape} != {K.shape}"
            )
        K = K + float(w_k) * one
    return 0.5 * (K + K.T)


_LOG_WEIGHT_BOUND = 8.0


def _fit_mt_ai_reml_multikernel(
    *,
    kernel_bank: Sequence[np.ndarray],
    y_std: np.ndarray,
    X: np.ndarray,
    gid_idx: np.ndarray,
    trait_idx: np.ndarray,
    train_mask: np.ndarray,
    T: int,
    n_gid: int,
    kernel_weights: Optional[np.ndarray] = None,
    estimate_kernel_weights: bool = False,
    Sigma_G_init: Optional[np.ndarray] = None,
    Sigma_eps_init: Optional[np.ndarray] = None,
    max_outer: int = 3,
    max_weight_evals: int = 80,
    weight_tol_loglik: float = 1e-4,
    **fit_kwargs: Any,
) -> Dict[str, Any]:
    """Joint multi-trait AI-REML over a kernel bank, optionally estimating weights.

    Model:  V = Sigma_G (x) sum_k w_k K_k  +  Sigma_eps (x) I

    `w[0]` is held at 1.0. The overall genetic scale is already carried by
    Sigma_G, so only the ratios w[k]/w[0] are identified; freeing all of them
    would leave the likelihood invariant along `Sigma_G -> c Sigma_G, w -> w/c`.

    With `estimate_kernel_weights=False` this is the historical fixed-weight fit
    with the supplied (or equal) weights, and the returned weights are exactly
    the ones passed in.

    With estimation on, the ratios are chosen by maximising the REML *profile*
    likelihood: each candidate mixture is scored by a full inner AI-REML fit of
    (Sigma_G, Sigma_eps) at that mixture, searched coordinate-wise with a
    bounded Brent line search over log w[1:]. The search is derivative-free --
    the combined kernel's eigendecomposition depends on the weights, and
    differentiating through `eigh` is ill-conditioned when that kernel has
    near-degenerate eigenvalues.

    The profile likelihood in the weights is often very flat: two kernels that
    describe similar relatedness can trade off with almost no likelihood cost.
    `kernel_weight_loglik_gain` reports how much was actually gained over the
    starting weights, so a near-zero gain can be read as "this dataset does not
    identify the mixture" rather than as a confident estimate.

    Returns the inner `_fit_mt_ai_reml` dict plus:
      kernel_weights             final weight vector (w[0] == 1.0 when estimated)
      kernel_weights_estimated   whether the ratios were fitted or supplied
      kernel_weight_outer_iter   completed coordinate sweeps
      kernel_weight_profile_evals inner fits spent scoring candidate mixtures
      kernel_weight_loglik_gain  REML log-lik improvement over the starting weights
    """
    bank = [np.asarray(K, dtype=np.float64) for K in kernel_bank]
    if not bank:
        raise ValueError("kernel_bank must contain at least one kernel")
    n_k = len(bank)

    if kernel_weights is None:
        w = np.ones(n_k, dtype=np.float64)
    else:
        w = np.asarray(kernel_weights, dtype=np.float64).reshape(-1)
        if w.size != n_k:
            raise ValueError(
                f"kernel_weights length {w.size} != number of kernels {n_k}"
            )
        if not np.isfinite(w).all() or (w < 0).any() or not (w > 0).any():
            raise ValueError(
                "kernel_weights must be finite, non-negative, and include a positive value"
            )

    # One kernel carries no estimable ratio.
    estimating = bool(estimate_kernel_weights) and n_k >= 2

    if Sigma_G_init is None or Sigma_eps_init is None:
        Sigma_G0, Sigma_eps0 = _initialize_params(
            y_std=y_std, X=X, gid_idx=gid_idx, trait_idx=trait_idx,
            train_mask=train_mask, T=T, n_gid=n_gid,
        )
        Sigma_G_init = Sigma_G_init if Sigma_G_init is not None else Sigma_G0
        Sigma_eps_init = Sigma_eps_init if Sigma_eps_init is not None else Sigma_eps0

    fit = _fit_mt_ai_reml(
        K_geno=_combine_kernel_bank(bank, w),
        y_std=y_std, X=X, gid_idx=gid_idx, trait_idx=trait_idx,
        train_mask=train_mask, T=T, n_gid=n_gid,
        Sigma_G_init=Sigma_G_init, Sigma_eps_init=Sigma_eps_init,
        **fit_kwargs,
    )

    if not estimating:
        fit["kernel_weights"] = w
        fit["kernel_weights_estimated"] = False
        fit["kernel_weight_at_bound"] = False
        fit["kernel_weight_outer_iter"] = 0
        fit["kernel_weight_profile_evals"] = 0
        fit["kernel_weight_loglik_gain"] = 0.0
        return fit

    start_loglik = float(fit["loglik_history"][-1])
    # Rescale so w[0] == 1: the overall genetic scale is Sigma_G's job.
    w = w / float(w[0]) if w[0] > 0 else np.ones(n_k, dtype=np.float64)

    from scipy.optimize import minimize_scalar

    eval_budget = {"n": 0}
    cache: Dict[tuple, Any] = {}
    best = {"ll": start_loglik, "fit": fit, "log_w": np.log(np.clip(
        w[1:], np.exp(-_LOG_WEIGHT_BOUND), np.exp(_LOG_WEIGHT_BOUND)))}

    def _fit_at(log_vec: np.ndarray):
        """Full inner AI-REML at fixed weights -> (profile loglik, fit).

        This is the exact REML profile likelihood in the weights: every
        evaluation re-maximises (Sigma_G, Sigma_eps) at the candidate mixture,
        so the search optimises the same objective the reported fit uses.
        """
        key = tuple(np.round(np.asarray(log_vec, dtype=np.float64), 6).tolist())
        hit = cache.get(key)
        if hit is not None:
            return hit
        if eval_budget["n"] >= int(max_weight_evals):
            return None
        eval_budget["n"] += 1
        wv = np.concatenate([
            [1.0],
            np.exp(np.clip(np.asarray(key, dtype=np.float64),
                           -_LOG_WEIGHT_BOUND, _LOG_WEIGHT_BOUND)),
        ])
        try:
            cand = _fit_mt_ai_reml(
                K_geno=_combine_kernel_bank(bank, wv),
                y_std=y_std, X=X, gid_idx=gid_idx, trait_idx=trait_idx,
                train_mask=train_mask, T=T, n_gid=n_gid,
                Sigma_G_init=Sigma_G_init, Sigma_eps_init=Sigma_eps_init,
                **fit_kwargs,
            )
        except Exception:
            cache[key] = None
            return None
        ll = float(cand["loglik_history"][-1])
        if not np.isfinite(ll):
            cache[key] = None
            return None
        out = (ll, cand)
        cache[key] = out
        if ll > best["ll"]:
            best["ll"] = ll
            best["fit"] = cand
            best["log_w"] = np.asarray(key, dtype=np.float64)
        return out

    log_w = best["log_w"].copy()
    sweeps_done = 0
    for _sweep in range(int(max_outer)):
        improved = False
        for j in range(log_w.size):
            def _neg(t: float, j=j) -> float:
                trial = log_w.copy()
                trial[j] = float(t)
                got = _fit_at(trial)
                return float(np.inf) if got is None else -got[0]

            lo = max(float(log_w[j]) - 3.0, -_LOG_WEIGHT_BOUND)
            hi = min(float(log_w[j]) + 3.0, _LOG_WEIGHT_BOUND)
            if not hi > lo:
                continue
            before = best["ll"]
            try:
                minimize_scalar(_neg, bounds=(lo, hi), method="bounded",
                                options={"xatol": 1e-2, "maxiter": 30})
            except Exception:
                pass
            # best[] tracks the best point any evaluation reached, so a search
            # that wandered still leaves the fit at its best seen mixture.
            if best["ll"] > before + 1e-9:
                log_w = best["log_w"].copy()
                improved = True
        sweeps_done += 1
        if not improved or eval_budget["n"] >= int(max_weight_evals):
            break

    fit = best["fit"]
    w = np.concatenate([[1.0], np.exp(np.clip(best["log_w"],
                                              -_LOG_WEIGHT_BOUND, _LOG_WEIGHT_BOUND))])
    prev_loglik = best["ll"]
    outer_done = sweeps_done
    fit["kernel_weight_profile_evals"] = int(eval_budget["n"])
    fit["kernel_weights"] = w
    fit["kernel_weights_estimated"] = True
    # A ratio pinned at the search bound is a result, not a glitch: the profile
    # likelihood kept improving as the kernel was driven out of (or made to
    # dominate) the mixture. Surface it so it is not read as a fitted interior
    # optimum.
    fit["kernel_weight_at_bound"] = bool(
        np.any(np.abs(best["log_w"]) >= _LOG_WEIGHT_BOUND - 1e-6)
    )
    fit["kernel_weight_outer_iter"] = int(outer_done)
    fit["kernel_weight_loglik_gain"] = float(prev_loglik - start_loglik)
    return fit


def _predict_mt(
    *,
    K_geno: np.ndarray,
    fit: Dict[str, Any],
    design: Dict[str, Any],
    gid_col: str,
    trait_col: str,
    device: Optional[Any] = None,
) -> tuple:
    """Form long-format predictions DataFrame + (n_test_gid, T, T) cov_traits.

    Posterior:
      g | y_train ~ N(mu_g, C)
      mu_g[i,t] = (Σ_G ⊗ K_{i, train}) V^{-1} (y_train - Xβ̂)
      C = (Σ_G ⊗ K_{ii}) − (Σ_G ⊗ K_{i,train}) V^{-1} (Σ_G ⊗ K_{train,i})

    Predictions are y_pred = X_test β̂ + mu_g, with SE_latent² = diag(C) and
    SE² = SE_latent² + diag(Σ_eps for each trait).
    """
    T = design["T"]
    n_gid = design["n_gid"]
    Sigma_G = fit["Sigma_G"]
    Sigma_eps = fit["Sigma_eps"]
    beta = fit["beta"]
    device = pick_torch_device(device)
    K_t = torch.tensor(K_geno, dtype=torch.float64, device=device)
    K_t = 0.5 * (K_t + K_t.T)
    Sigma_G_t = torch.tensor(Sigma_G, dtype=torch.float64, device=device)
    Sigma_eps_t = torch.tensor(Sigma_eps, dtype=torch.float64, device=device)
    y_std = torch.tensor(design["y_std"], dtype=torch.float64, device=device)
    X = torch.tensor(design["X_fixed"], dtype=torch.float64, device=device)
    gid_idx = torch.tensor(design["gid_idx"], dtype=torch.int64, device=device)
    trait_idx = torch.tensor(design["trait_idx"], dtype=torch.int64, device=device)
    train_mask = torch.tensor(design["train_mask"], dtype=torch.bool, device=device)
    test_mask = torch.tensor(design["test_mask"], dtype=torch.bool, device=device)
    beta_t = torch.tensor(beta, dtype=torch.float64, device=device)

    g_tr = gid_idx[train_mask]
    t_tr = trait_idx[train_mask]
    y_tr = y_std[train_mask]
    X_tr = X[train_mask]
    K_tr_tr = K_t[g_tr.unsqueeze(1), g_tr.unsqueeze(0)]
    G_tr = Sigma_G_t[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]
    E_tr = Sigma_eps_t[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]
    same_gid_tr = (g_tr.unsqueeze(1) == g_tr.unsqueeze(0)).to(K_tr_tr.dtype)
    V = G_tr * K_tr_tr + E_tr * same_gid_tr
    V = V + 1e-6 * torch.eye(V.shape[0], dtype=V.dtype, device=device)
    L_V = torch.linalg.cholesky(V)
    r_tr = y_tr - X_tr @ beta_t
    Vinv_r = torch.cholesky_solve(r_tr.unsqueeze(1), L_V).squeeze(1)

    g_te = gid_idx[test_mask]
    t_te = trait_idx[test_mask]
    X_te = X[test_mask]
    K_te_tr = K_t[g_te.unsqueeze(1), g_tr.unsqueeze(0)]
    G_te_tr = Sigma_G_t[t_te.unsqueeze(1), t_tr.unsqueeze(0)]
    cross = G_te_tr * K_te_tr
    mu_g_std = cross @ Vinv_r
    y_pred_std = X_te @ beta_t + mu_g_std
    Vinv_cross = torch.cholesky_solve(cross.T, L_V).T  # (n_te, n_tr)
    K_te_te_diag = K_t[g_te, g_te]
    G_te_te_diag = Sigma_G_t[t_te, t_te]
    var_lat = G_te_te_diag * K_te_te_diag - (cross * Vinv_cross).sum(dim=1)
    var_lat = torch.clamp(var_lat, min=1e-10)
    E_te_te_diag = Sigma_eps_t[t_te, t_te]
    var_obs = var_lat + E_te_te_diag

    trait_means = torch.tensor(design["trait_means"], dtype=torch.float64, device=device)
    trait_sds = torch.tensor(design["trait_sds"], dtype=torch.float64, device=device)
    y_pred = y_pred_std * trait_sds[t_te] + trait_means[t_te]
    SE_lat = torch.sqrt(var_lat) * trait_sds[t_te]
    SE = torch.sqrt(var_obs) * trait_sds[t_te]
    var_lat_orig = var_lat * trait_sds[t_te] ** 2
    var_obs_orig = var_obs * trait_sds[t_te] ** 2

    test_gid_unique = torch.unique(g_te).tolist()
    n_test_gid = len(test_gid_unique)
    cov_traits = np.zeros((n_test_gid, T, T), dtype=np.float64)
    for s, gi in enumerate(test_gid_unique):
        g_query = torch.tensor([gi] * T, dtype=torch.int64, device=device)
        t_query = torch.arange(T, dtype=torch.int64)
        K_q_tr = K_t[g_query.unsqueeze(1), g_tr.unsqueeze(0)]
        G_q_tr = Sigma_G_t[t_query.unsqueeze(1), t_tr.unsqueeze(0)]
        cross_q = G_q_tr * K_q_tr
        Vinv_cross_q = torch.cholesky_solve(cross_q.T, L_V).T
        K_q_q = K_t[gi, gi]
        Sigma_lat_q = K_q_q * Sigma_G_t - cross_q @ Vinv_cross_q.T
        sd = trait_sds
        Sigma_lat_q_orig = (sd.unsqueeze(1) * sd.unsqueeze(0)) * Sigma_lat_q
        cov_traits[s] = torch_to_numpy(0.5 * (Sigma_lat_q_orig + Sigma_lat_q_orig.T))

    trait_levels = design["trait_levels"]
    geno_ids = design["geno_ids"]
    df_out = pd.DataFrame({
        gid_col: [geno_ids[int(g)] for g in g_te.tolist()],
        trait_col: [trait_levels[int(t)] for t in t_te.tolist()],
        "Prediction": torch_to_numpy(y_pred),
        "SE": torch_to_numpy(SE),
        "SE_latent": torch_to_numpy(SE_lat),
        "PEV": torch_to_numpy(var_lat_orig),
        "Prediction_Var_latent": torch_to_numpy(var_lat_orig),
        "Prediction_Var_observed": torch_to_numpy(var_obs_orig),
    })
    for col in (gid_col, trait_col):
        df_out[col] = pd.Series([str(x) for x in df_out[col].tolist()], index=df_out.index, dtype=object)
    return df_out, cov_traits


def _pack_unstructured(L: "torch.Tensor") -> "torch.Tensor":
    """Pack lower-triangular T×T factor into a flat (T(T+1)/2,) vector.

    Order: row-major over the lower triangle including diagonal.
    """
    T = L.shape[0]
    rows, cols = torch.tril_indices(T, T)
    return L[rows, cols]


def _unpack_unstructured(theta: "torch.Tensor", T: int) -> "torch.Tensor":
    """Inverse of _pack_unstructured. Returns dense T×T lower-Cholesky factor."""
    L = torch.zeros((T, T), dtype=theta.dtype, device=theta.device)
    rows, cols = torch.tril_indices(T, T)
    L[rows, cols] = theta
    return L


def _pack_fa(Lambda: "torch.Tensor", log_psi: "torch.Tensor") -> "torch.Tensor":
    """Pack (Λ, log_ψ) into flat parameter vector.

    Λ is T×k lower-triangular in its first k×k block. We store:
      - Free entries of Λ in row-major order (excluding upper-triangular zeros of the first k×k).
      - log_ψ as the last T entries.
    """
    T, k = Lambda.shape
    parts = []
    for i in range(T):
        for j in range(min(i + 1, k)):
            parts.append(Lambda[i, j])
    parts.append(log_psi)
    return torch.cat([
        torch.stack([p for p in parts[:-1]]),
        parts[-1],
    ])


def _unpack_fa(theta: "torch.Tensor", T: int, k: int) -> tuple:
    """Inverse of _pack_fa. Returns (Lambda T×k, log_psi length-T)."""
    n_lambda = T * k - k * (k - 1) // 2
    assert theta.shape == (n_lambda + T,), \
        f"theta length {theta.shape[0]} != n_lambda+T {n_lambda + T}"
    Lambda = torch.zeros((T, k), dtype=theta.dtype, device=theta.device)
    pos = 0
    for i in range(T):
        for j in range(min(i + 1, k)):
            Lambda[i, j] = theta[pos]
            pos += 1
    log_psi = theta[pos:]
    return Lambda, log_psi


def _neg_reml_loglik_balanced(
    *,
    Sigma_G: "torch.Tensor",      # (T, T) PSD
    Sigma_eps: "torch.Tensor",    # (T, T) PSD
    y_std: "torch.Tensor",        # (N,)
    X: "torch.Tensor",            # (N, p) trait-expanded fixed effects
    gid_idx: "torch.Tensor",      # (N,) int64
    trait_idx: "torch.Tensor",    # (N,) int64
    train_mask: "torch.Tensor",   # (N,) bool
    D: "torch.Tensor",            # (n_gid,) eigenvalues of K
    U: "torch.Tensor",            # (n_gid, n_gid) eigenvectors of K
    T: int,
    n_gid: int,
) -> "torch.Tensor":
    """Negative REML log-likelihood via Kronecker eigendecomp on balanced grid.

    PRECONDITION: caller must have verified _is_balanced_grid(...) is True.
    Math: V = Sigma_G ⊗ K + Sigma_eps ⊗ I is block-diagonal in K's eigenbasis with
    T×T blocks B_i = d_i Sigma_G + Sigma_eps. We work in U-basis throughout.
    """
    # Restrict to training rows
    y_tr = y_std[train_mask]
    X_tr = X[train_mask]
    g_tr = gid_idx[train_mask]
    t_tr = trait_idx[train_mask]
    p = X_tr.shape[1]

    # Reshape long → (n_gid, T) for Y and (n_gid, T, p) for X
    Y_grid = torch.zeros(n_gid, T, dtype=y_std.dtype, device=y_std.device)
    Y_grid[g_tr, t_tr] = y_tr
    X_grid = torch.zeros(n_gid, T, p, dtype=X.dtype, device=X.device)
    X_grid[g_tr, t_tr, :] = X_tr

    # U-basis transform on the gid axis
    Y_tilde = U.T @ Y_grid  # (n_gid, T)
    X_tilde = torch.einsum("gi,itp->gtp", U.T, X_grid)  # (n_gid, T, p)

    # Per-eigenvalue T×T blocks: B_i = d_i Σ_G + Σ_eps  (shape n_gid × T × T)
    B = D.unsqueeze(-1).unsqueeze(-1) * Sigma_G.unsqueeze(0) + Sigma_eps.unsqueeze(0)
    # Match dense path jitter: dense adds 1e-6*I_{n_tr} to V_obs, which in the
    # Kronecker eigenbasis becomes 1e-6*I_T on each B_i block.
    B = B + 1e-6 * torch.eye(T, dtype=B.dtype, device=B.device).unsqueeze(0)
    L_B = torch.linalg.cholesky(B)

    # log|V| = sum_i log|B_i| = 2 sum_i sum_t log L_B[i, t, t]
    log_det_V = 2.0 * torch.log(torch.diagonal(L_B, dim1=-2, dim2=-1)).sum()

    # V^{-1} ỹ in U-basis: per-row B_i^{-1} ỹ_i, batched cholesky_solve
    Vinv_Y = torch.cholesky_solve(Y_tilde.unsqueeze(-1), L_B).squeeze(-1)  # (n_gid, T)
    Vinv_X = torch.cholesky_solve(X_tilde, L_B)  # (n_gid, T, p)

    # X̃^T V^{-1} X̃ = sum_i X̃[i,:,:]^T B_i^{-1} X̃[i,:,:]
    Xt_Vinv_X = torch.einsum("itp,itq->pq", X_tilde, Vinv_X)
    Xt_Vinv_y = torch.einsum("itp,it->p", X_tilde, Vinv_Y)

    # GLS β̂ via Cholesky on (p,p)
    L_xvx = torch.linalg.cholesky(
        Xt_Vinv_X + 1e-6 * torch.eye(p, dtype=B.dtype, device=B.device)
    )
    log_det_xvx = 2.0 * torch.log(torch.diagonal(L_xvx)).sum()
    beta_hat = torch.cholesky_solve(Xt_Vinv_y.unsqueeze(-1), L_xvx).squeeze(-1)

    # Residual r̃ in U-basis
    Xb_tilde = torch.einsum("itp,p->it", X_tilde, beta_hat)  # (n_gid, T)
    r_tilde = Y_tilde - Xb_tilde
    Vinv_r = torch.cholesky_solve(r_tilde.unsqueeze(-1), L_B).squeeze(-1)
    quad = (r_tilde * Vinv_r).sum()

    nll = 0.5 * (log_det_V + log_det_xvx + quad)
    return nll


def _neg_reml_loglik(
    *,
    Sigma_G: "torch.Tensor",      # (T, T) PSD
    Sigma_eps: "torch.Tensor",    # (T, T) PSD
    y_std: "torch.Tensor",        # (N,)
    X: "torch.Tensor",            # (N, p) trait-expanded fixed effects
    gid_idx: "torch.Tensor",      # (N,) int64
    trait_idx: "torch.Tensor",    # (N,) int64
    train_mask: "torch.Tensor",   # (N,) bool
    D: "torch.Tensor",            # (n_gid,) eigenvalues of K
    U: "torch.Tensor",            # (n_gid, n_gid) eigenvectors of K
    T: int,
    n_gid: int,
) -> "torch.Tensor":
    """Negative REML log-likelihood, restricted to training rows.

    Mathematically: REML log-lik = -½(log|V| + log|X^T V^-1 X| + r^T V^-1 r)
    where r = y - X β̂_GLS, β̂_GLS = (X^T V^-1 X)^-1 X^T V^-1 y.

    Auto-routes to the Kronecker fast path (`_neg_reml_loglik_balanced`) when the
    (gid × trait) grid is complete (n_obs == n_gid * T after train_mask).
    Otherwise falls through to the dense V_obs Cholesky body.
    """
    g_tr_check = gid_idx[train_mask]
    t_tr_check = trait_idx[train_mask]
    if _is_balanced_grid(g_tr_check, t_tr_check, n_gid, T):
        return _neg_reml_loglik_balanced(
            Sigma_G=Sigma_G, Sigma_eps=Sigma_eps,
            y_std=y_std, X=X,
            gid_idx=gid_idx, trait_idx=trait_idx,
            train_mask=train_mask,
            D=D, U=U, T=T, n_gid=n_gid,
        )

    # Restrict to training rows
    y_tr = y_std[train_mask]
    X_tr = X[train_mask]
    g_tr = gid_idx[train_mask]
    t_tr = trait_idx[train_mask]
    n_tr = y_tr.shape[0]
    p = X_tr.shape[1]

    # v1 simplicity: build V_obs (n_tr × n_tr) explicitly using Kronecker indexing.
    # V[r, s] = Sigma_G[t_r, t_s] * K[g_r, g_s] + Sigma_eps[t_r, t_s] * δ(g_r, g_s)
    # Cost: O(n_tr²) memory, O(n_tr³) Cholesky. At n_tr=4500 (n_gid=1500, T=3) that's
    # 162 MB and ~4 GFLOPs — under 1s in float64. Kronecker-shortcut left to operator-backend v2.
    # D, U passed for prediction reuse — unused here.
    K_full = U @ torch.diag(D) @ U.T  # (n_gid, n_gid) reconstructed for indexing
    K_pairs = K_full[g_tr.unsqueeze(1), g_tr.unsqueeze(0)]   # (n_tr, n_tr)
    G_pairs = Sigma_G[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]  # (n_tr, n_tr)
    E_pairs = Sigma_eps[t_tr.unsqueeze(1), t_tr.unsqueeze(0)]  # (n_tr, n_tr)
    same_gid = (g_tr.unsqueeze(1) == g_tr.unsqueeze(0)).to(K_pairs.dtype)
    V = G_pairs * K_pairs + E_pairs * same_gid
    V = V + 1e-6 * torch.eye(n_tr, dtype=V.dtype, device=V.device)

    L = torch.linalg.cholesky(V)
    log_det_V = 2.0 * torch.log(torch.diagonal(L)).sum()
    Vinv_y = torch.cholesky_solve(y_tr.unsqueeze(1), L).squeeze(1)
    Vinv_X = torch.cholesky_solve(X_tr, L)
    Xt_Vinv_X = X_tr.T @ Vinv_X
    L_xvx = torch.linalg.cholesky(
        Xt_Vinv_X + 1e-6 * torch.eye(p, dtype=V.dtype, device=V.device)
    )
    log_det_xvx = 2.0 * torch.log(torch.diagonal(L_xvx)).sum()
    beta_hat = torch.cholesky_solve(
        (X_tr.T @ Vinv_y).unsqueeze(1), L_xvx
    ).squeeze(1)
    r = y_tr - X_tr @ beta_hat
    Vinv_r = torch.cholesky_solve(r.unsqueeze(1), L).squeeze(1)
    quad = (r * Vinv_r).sum()

    nll = 0.5 * (log_det_V + log_det_xvx + quad)
    return nll


def _fit_mt_mom_op(
    *,
    design: Dict[str, Any],
    cache: Any,                       # GRMFactorCacheBase
    kernel_weights: Optional[Sequence[float]] = None,
    psd_floor: float = 1e-4,
    pcg_tol: float = 1e-4,
    pcg_max_iter: int = 300,
    return_cov_traits: bool = False,
    force_cov_traits: bool = False,
    force_prediction_se: bool = False,
    device: Optional[Any] = None,
) -> Dict[str, Any]:
    """v2-fast MoM orchestrator: MoM Sigma -> operators -> PCG -> predict.

    Returns dict with: Sigma_G, Sigma_eps, beta, mom_n_full, mom_n_partial,
    mom_psd_clip_count, pcg_n_iter, pcg_residual_norm, predictions, cov_traits.
    """
    import torch
    from mt_gp_mom import _estimate_sigma_mom
    from mt_gp_operators import MTGenoOperator, MTResidualOperator, _MTSumOperator
    from mixed_model_gpu_large_scale_backend import (
        pcg_solve, FoldIndexCache, _unique_inverse,
    )

    device = pick_torch_device(device)
    T = design["T"]
    n_gid = design["n_gid"]
    n_used = cache.n_geno()
    if n_used != n_gid:
        raise ValueError(
            f"cache.n_geno()={n_used} != n_gid={n_gid} from design; "
            f"geno_ids must align with cache row order."
        )
    if kernel_weights is None:
        weights = np.ones(int(cache.num_grms), dtype=np.float64)
    else:
        weights = np.asarray(kernel_weights, dtype=np.float64).reshape(-1)
    if weights.size != int(cache.num_grms):
        raise ValueError(
            f"kernel_weights length {weights.size} != cache.num_grms={int(cache.num_grms)}"
        )
    if not np.isfinite(weights).all() or (weights < 0).any() or not (weights > 0).any():
        raise ValueError("kernel_weights must be finite, non-negative, and include a positive value")
    all_rows = torch.arange(n_gid, dtype=torch.int64)
    phi_blocks = []
    for r, weight in enumerate(weights.tolist()):
        if weight <= 0:
            continue
        phi_r = cache.get_rows(
            r, all_rows, device="cpu", dtype=torch.float64,
        ).detach().cpu().numpy()
        phi_blocks.append(phi_r * np.sqrt(weight))
    Phi_full_rows = np.concatenate(phi_blocks, axis=1)

    mom_out = _estimate_sigma_mom(
        Phi=Phi_full_rows,
        y_std=design["y_std"], X_fixed=design["X_fixed"],
        gid_idx=design["gid_idx"], trait_idx=design["trait_idx"],
        train_mask=design["train_mask"],
        T=T, n_gid=n_gid,
        psd_floor=psd_floor,
    )
    Sigma_G = mom_out["Sigma_G"]
    Sigma_eps = mom_out["Sigma_eps"]

    # Build operators on TRAINING rows
    train_mask = design["train_mask"]
    g_tr = design["gid_idx"][train_mask]
    t_tr = design["trait_idx"][train_mask]
    y_tr_std = design["y_std"][train_mask]
    X_tr = design["X_fixed"][train_mask]
    op_g = MTGenoOperator(
        cache,
        gi=torch.tensor(g_tr, dtype=torch.int64, device=device),
        ti=torch.tensor(t_tr, dtype=torch.int64, device=device),
        Sigma_G=torch.tensor(Sigma_G, dtype=torch.float64, device=device),
        kernel_weights=torch.tensor(weights, dtype=torch.float64, device=device),
    )
    op_e = MTResidualOperator(
        gi=torch.tensor(g_tr, dtype=torch.int64, device=device),
        ti=torch.tensor(t_tr, dtype=torch.int64, device=device),
        Sigma_eps=torch.tensor(Sigma_eps, dtype=torch.float64, device=device),
    )
    V_op = _MTSumOperator(op_g, op_e)
    fold_train, ti_train_fold = op_g.build_fold_cache()

    def A_mv(v):
        return V_op.matvec(v, fold_train, ti_train_fold, dtype_compute=torch.float64)

    # GLS fixed effects under the MoM covariance: beta = (X' V^-1 X)^+ X' V^-1 y.
    # The OLS beta-hat returned by _estimate_sigma_mom is used only for
    # variance-component estimation; the BLUP posterior mean requires GLS so that
    # the prediction is the exact empirical ICM-BLUP rather than an OLS-centred
    # approximation. V^-1 is applied with the same PCG used for alpha. pinv keeps
    # this well-defined when the trait-expanded fixed-effect design is rank
    # deficient (one-hot factors are collinear with the per-trait intercepts).
    y_tr_t = torch.tensor(y_tr_std, dtype=torch.float64, device=device)
    X_tr_t = torch.tensor(X_tr, dtype=torch.float64, device=device)
    Vinv_y = pcg_solve(A_mv, y_tr_t, M_inv_mv=None, tol=pcg_tol, max_iter=pcg_max_iter).x
    p_fixed = int(X_tr_t.shape[1])
    Vinv_X_cols = []
    for j in range(p_fixed):
        sol_j = pcg_solve(A_mv, X_tr_t[:, j].contiguous(), M_inv_mv=None,
                          tol=pcg_tol, max_iter=pcg_max_iter)
        Vinv_X_cols.append(sol_j.x)
    Vinv_X = torch.stack(Vinv_X_cols, dim=1) if p_fixed > 0 else X_tr_t.new_zeros((X_tr_t.shape[0], 0))
    XtVinvX = X_tr_t.transpose(0, 1) @ Vinv_X
    XtVinvy = X_tr_t.transpose(0, 1) @ Vinv_y
    beta_t = torch.linalg.pinv(XtVinvX) @ XtVinvy
    beta = torch_to_numpy(beta_t)

    r = y_tr_t - X_tr_t @ beta_t
    pcg_res = pcg_solve(A_mv, r, M_inv_mv=None, tol=pcg_tol, max_iter=pcg_max_iter)
    alpha = pcg_res.x

    # Predict at test rows
    test_mask = design["test_mask"]
    g_te = design["gid_idx"][test_mask]
    t_te = design["trait_idx"][test_mask]
    X_te = design["X_fixed"][test_mask]
    test_state = {
        "fold_train": fold_train,
        "ti_train":   ti_train_fold,
        "gi_test":    torch.tensor(g_te, dtype=torch.int64, device=device),
        "ti_test":    torch.tensor(t_te, dtype=torch.int64, device=device),
    }
    mu_g_std = torch_to_numpy(op_g.cross_cov_apply(alpha, test_state, dtype_compute=torch.float64))
    y_pred_std = X_te @ beta + mu_g_std

    trait_means = design["trait_means"]
    trait_sds   = design["trait_sds"]
    y_pred = y_pred_std * trait_sds[t_te] + trait_means[t_te]
    trait_levels = design["trait_levels"]
    geno_ids = design["geno_ids"]
    preds_df = pd.DataFrame({
        "gid":        [geno_ids[int(g)] for g in g_te],
        "trait":      [trait_levels[int(t)] for t in t_te],
        "Prediction": y_pred,
    })
    for col in ("gid", "trait"):
        preds_df[col] = pd.Series([str(x) for x in preds_df[col].tolist()], index=preds_df.index, dtype=object)

    cov_traits = None
    compute_cov_traits = False
    if return_cov_traits:
        from mt_gp_reporting import resolve_cov_traits_request

        n_test_gid = int(np.unique(g_te).size)
        decision = resolve_cov_traits_request(
            n_cells=n_test_gid,
            n_traits=T,
            requested=True,
            force_cov_traits=bool(force_cov_traits),
            force_prediction_se=bool(force_prediction_se),
            context="multi-trait GP (single environment; varcomp_mode='reml' "
                    "uses a dense Cholesky path instead)",
        )
        compute_cov_traits = bool(decision["compute"])
    if compute_cov_traits:
        unique_te_gids = np.unique(g_te)
        cov_traits = np.zeros((unique_te_gids.size, T, T), dtype=np.float64)
        Sigma_G_t = torch.tensor(Sigma_G, dtype=torch.float64, device=device)
        for s, gi in enumerate(unique_te_gids):
            query_row = torch.tensor([int(gi)], dtype=torch.int64)
            K_row = torch.zeros(fold_train.ng_used, dtype=torch.float64, device=device)
            K_gi_gi = 0.0
            for r, weight in enumerate(weights.tolist()):
                if weight <= 0:
                    continue
                Phi_row_gi = cache.get_rows(
                    r, query_row, device=device, dtype=torch.float64
                )
                Phi_train = cache.get_rows(
                    r, fold_train.g_unique, device=device, dtype=torch.float64
                )
                K_row = K_row + float(weight) * (
                    Phi_row_gi @ Phi_train.transpose(0, 1)
                ).squeeze(0)
                K_gi_gi += float(weight) * float(
                    (Phi_row_gi @ Phi_row_gi.transpose(0, 1)).item()
                )
            schur = np.zeros((T, T), dtype=np.float64)
            for t in range(T):
                col = (Sigma_G_t[t][torch.tensor(t_tr, dtype=torch.int64, device=device)] *
                       K_row[fold_train.inv]).to(torch.float64)
                pcg_t = pcg_solve(A_mv, col, M_inv_mv=None, tol=pcg_tol, max_iter=pcg_max_iter)
                schur[t, :] = torch_to_numpy(Sigma_G_t[t] * K_gi_gi) - torch_to_numpy(col) @ torch_to_numpy(pcg_t.x)
            schur = 0.5 * (schur + schur.T)
            sd = trait_sds
            cov_traits[s] = (np.outer(sd, sd)) * schur

    return {
        "Sigma_G":             Sigma_G,
        "Sigma_eps":           Sigma_eps,
        "kernel_weights":      weights,
        "beta":                beta,
        "mom_n_full":          mom_out["mom_n_full"],
        "mom_n_partial":       mom_out["mom_n_partial"],
        "mom_psd_clip_count":  mom_out["mom_psd_clip_count"],
        "pcg_n_iter":          int(pcg_res.iters),
        "pcg_residual_norm":   float(pcg_res.rel_resid),
        "device":              str(device),
        "predictions":         preds_df,
        "cov_traits":          cov_traits,
    }
