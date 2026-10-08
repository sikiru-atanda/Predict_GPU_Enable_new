"""Method-of-moments Sigma_G/Sigma_eps estimation + PSD projection for v2-fast path.

Pure-numerical (numpy + light torch). No operator imports.

Public:
  _estimate_sigma_mom — common-matrix MoM on I_full subset.
  _psd_project        — eigenvalue-floor PSD projection.
"""
from __future__ import annotations
from typing import Any, Dict, Tuple
import warnings
import numpy as np


def _estimate_sigma_mom(
    *,
    Phi: np.ndarray,                  # (n_gid, m) low-rank factor; K = Phi @ Phi.T (+ ridge implicit)
    y_std: np.ndarray,                # (N,) standardized training y
    X_fixed: np.ndarray,              # (N, p) trait-expanded fixed-effect matrix
    gid_idx: np.ndarray,              # (N,) int64
    trait_idx: np.ndarray,            # (N,) int64
    train_mask: np.ndarray,           # (N,) bool
    T: int,
    n_gid: int,
    psd_floor: float = 1e-4,
) -> Dict[str, Any]:
    """Method-of-moments Sigma_G, Sigma_eps on the I_full subset (gids with all T traits observed).

    Returns dict with:
      Sigma_G, Sigma_eps        — (T, T) PSD matrices in standardized space.
      beta                      — (p,) OLS beta-hat on I_full restricted X.
      mom_n_full                — |I_full|
      mom_n_partial             — n_gid - |I_full|
      mom_psd_clip_count        — (n_clipped_G, n_clipped_eps)
    """
    # Restrict to training rows
    y_tr = y_std[train_mask]
    X_tr = X_fixed[train_mask]
    g_tr = gid_idx[train_mask]
    t_tr = trait_idx[train_mask]

    # Detect I_full: gids with observations for all T traits
    obs_count = np.zeros(n_gid, dtype=np.int64)
    np.add.at(obs_count, g_tr, 1)
    in_full = obs_count == T
    I_full = np.where(in_full)[0]
    n_full = int(I_full.size)
    n_partial = int(n_gid - n_full)

    if n_full < T + 1:
        raise ValueError(
            f"MoM requires |I_full| >= T+1 = {T+1} fully-observed gids; got {n_full}. "
            f"Use varcomp_mode='reml' for partial-cell-only training."
        )

    # Restrict to I_full rows
    full_row_mask = in_full[g_tr]
    y_full = y_tr[full_row_mask]
    X_full = X_tr[full_row_mask]
    g_full = g_tr[full_row_mask]
    t_full = t_tr[full_row_mask]

    # OLS beta-hat on I_full subset
    beta_hat, *_ = np.linalg.lstsq(X_full, y_full, rcond=None)

    # Residualize
    r = y_full - X_full @ beta_hat

    # Build Z: (n_full, T) residual matrix indexed by gid position in I_full
    gid_to_pos = -np.ones(n_gid, dtype=np.int64)
    gid_to_pos[I_full] = np.arange(n_full)
    Z = np.zeros((n_full, T), dtype=np.float64)
    pos = gid_to_pos[g_full]
    Z[pos, t_full] = r

    # Phi restricted to I_full rows
    Phi_full = Phi[I_full, :]

    # Common system matrix M (2x2)
    # tr(K_full) = ||Phi_full||_F^2 = sum(Phi_full * Phi_full)
    # ||K_full||_F^2 = ||Phi_full^T Phi_full||_F^2 = ||A||_F^2 where A = Phi_full.T @ Phi_full
    A = Phi_full.T @ Phi_full          # (m, m)
    tr_K = float(np.sum(Phi_full * Phi_full))
    tr_K2 = float(np.sum(A * A))
    M = np.array([[tr_K2, tr_K],
                  [tr_K, float(n_full)]], dtype=np.float64)
    M_reg = M + 1e-10 * np.eye(2)

    # Precompute K @ Z = Phi_full @ (Phi_full.T @ Z)
    PhiT_Z = Phi_full.T @ Z            # (m, T)
    KZ = Phi_full @ PhiT_Z             # (n_full, T)

    # Per-pair (a, b) MoM solve
    Sigma_G = np.zeros((T, T), dtype=np.float64)
    Sigma_eps = np.zeros((T, T), dtype=np.float64)
    for a in range(T):
        for b in range(a, T):
            rhs = np.array([
                float(Z[:, a] @ KZ[:, b]),
                float(Z[:, a] @ Z[:, b]),
            ], dtype=np.float64)
            sol = np.linalg.solve(M_reg, rhs)
            Sigma_G[a, b] = sol[0]
            Sigma_G[b, a] = sol[0]
            Sigma_eps[a, b] = sol[1]
            Sigma_eps[b, a] = sol[1]

    # PSD project both matrices
    Sigma_G, n_clip_G = _psd_project(Sigma_G, floor=psd_floor, name="Sigma_G")
    Sigma_eps, n_clip_E = _psd_project(Sigma_eps, floor=psd_floor, name="Sigma_eps")

    return {
        "Sigma_G": Sigma_G,
        "Sigma_eps": Sigma_eps,
        "beta": beta_hat,
        "mom_n_full": n_full,
        "mom_n_partial": n_partial,
        "mom_psd_clip_count": (int(n_clip_G), int(n_clip_E)),
    }


def _psd_project(M: np.ndarray, floor: float = 1e-4, name: str = "matrix",
                 warn_threshold: float = 0.25) -> Tuple[np.ndarray, int]:
    """Project symmetric matrix onto PSD cone with eigenvalue floor.

    Returns (M_projected, n_clipped). Emits UserWarning if more than
    `warn_threshold` fraction of eigenvalues required clipping.
    """
    M_sym = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M_sym)
    n_clipped = int((w < floor).sum())
    w_clipped = np.maximum(w, floor)
    M_proj = V @ np.diag(w_clipped) @ V.T
    M_proj = 0.5 * (M_proj + M_proj.T)
    T = M.shape[0]
    if n_clipped / T > warn_threshold:
        warnings.warn(
            f"{name} has {n_clipped}/{T} eigenvalues at floor {floor}; "
            f"data is sparse for MoM. Consider varcomp_mode='reml' or "
            f"trait_structure='fa' with smaller rank.",
            UserWarning, stacklevel=2,
        )
    return M_proj, n_clipped
