"""Internal math + helpers for fit_multi_trait_gp_met.

Not a public API. Imports v2-fast operators + MoM helpers.
"""
from __future__ import annotations
from typing import Any, Dict, Optional, Sequence, Tuple
import warnings
import numpy as np
import torch

from gp_device import pick_torch_device, torch_to_numpy


def _normalize_kernel_weights(cache: Any, kernel_weights: Optional[Sequence[float]]) -> np.ndarray:
    n = int(cache.num_grms)
    if kernel_weights is None:
        weights = np.ones(n, dtype=np.float64)
    else:
        weights = np.asarray(kernel_weights, dtype=np.float64).reshape(-1)
    if weights.size != n:
        raise ValueError(f"kernel_weights length {weights.size} != grm_factor_cache.num_grms={n}")
    if not np.isfinite(weights).all():
        raise ValueError("kernel_weights must be finite")
    if (weights < 0).any():
        raise ValueError("kernel_weights must be non-negative")
    if float(weights.sum()) <= 0.0:
        raise ValueError("at least one kernel weight must be positive")
    return weights


def _materialize_weighted_phi(cache: Any, n_gid: int, kernel_weights: np.ndarray) -> np.ndarray:
    """Materialize a weighted low-rank root for K = sum_r w_r Phi_r Phi_r'."""
    rows = torch.arange(n_gid, dtype=torch.int64)
    blocks = []
    for r, w in enumerate(kernel_weights.tolist()):
        if w <= 0.0:
            continue
        Phi_r = torch_to_numpy(cache.get_rows(r, rows, device="cpu", dtype=torch.float64))
        blocks.append(Phi_r * np.sqrt(float(w)))
    if not blocks:
        raise ValueError("at least one positive kernel weight is required")
    return np.concatenate(blocks, axis=1) if len(blocks) > 1 else blocks[0]


def _fa_project_covariance(
    M: np.ndarray,
    rank: int,
    *,
    floor: float,
    name: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Principal-factor PSD approximation: Lambda Lambda' + diag(psi)."""
    from mt_gp_mom import _psd_project
    T = int(M.shape[0])
    k = int(rank)
    if k < 1 or k >= T:
        raise ValueError(f"{name} FA rank must satisfy 1 <= rank < T={T}; got {rank}")
    M_psd, _ = _psd_project(M, floor=floor, name=name)
    w, V = np.linalg.eigh(M_psd)
    order = np.argsort(w)[::-1]
    top = order[:k]
    Lambda = V[:, top] * np.sqrt(np.maximum(w[top], floor))
    # Deterministic signs keep output stable across LAPACK sign choices.
    for j in range(k):
        i = int(np.argmax(np.abs(Lambda[:, j])))
        if Lambda[i, j] < 0:
            Lambda[:, j] *= -1.0
    psi = np.diag(M_psd) - np.sum(Lambda * Lambda, axis=1)
    psi = np.maximum(psi, floor)
    M_fa = Lambda @ Lambda.T + np.diag(psi)
    M_fa, _ = _psd_project(M_fa, floor=floor, name=f"{name}_fa")
    return M_fa, Lambda, psi


def _normalize_row_indices(idx: np.ndarray, n_rows: int, name: str) -> np.ndarray:
    """Normalize integer row indices or boolean masks against the original phenotype table."""
    arr = np.asarray(idx)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D integer index array or boolean mask")
    if arr.dtype == np.bool_:
        if arr.size != n_rows:
            raise ValueError(f"{name} boolean mask length {arr.size} != n_rows={n_rows}")
        return np.flatnonzero(arr).astype(np.int64)
    if not np.issubdtype(arr.dtype, np.integer):
        raise ValueError(f"{name} must contain integer row indices or be a boolean mask")
    out = arr.astype(np.int64, copy=False)
    bad = out[(out < 0) | (out >= n_rows)]
    if bad.size:
        raise ValueError(
            f"{name} contains out-of-bounds row indices for n_rows={n_rows}: "
            f"{bad[:5].tolist()}"
        )
    if np.unique(out).size != out.size:
        raise ValueError(f"{name} contains duplicate row indices")
    return out


def _build_mt_met_design(
    *,
    pheno_df,
    gid_col: str, env_col: str, trait_col: str, y_col: str,
    fixed_effects: Sequence[str],
    geno_ids: Sequence[str],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> Dict[str, Any]:
    """Build long-format MET design with per-(env, trait) standardization.

    Returns dict with: T, n_gid, n_env, gid_idx, env_idx, trait_idx, y_std,
    X_fixed, train_mask, test_mask, trait_levels, env_levels, trait_env_means,
    trait_env_sds, geno_ids.
    """
    import pandas as pd
    df = pheno_df.copy()
    n_orig = len(df)
    train_idx = _normalize_row_indices(train_idx, n_orig, "train_idx")
    test_idx = _normalize_row_indices(test_idx, n_orig, "test_idx")
    overlap = np.intersect1d(train_idx, test_idx)
    if overlap.size:
        raise ValueError(
            "train_idx and test_idx must be disjoint; overlapping rows include "
            f"{overlap[:5].tolist()}"
        )
    keep = np.zeros(n_orig, dtype=bool)
    keep[train_idx] = True
    keep[test_idx] = True
    df = df.loc[keep].copy()
    for col in (gid_col, env_col, trait_col):
        df[col] = df[col].astype("object").map(lambda x: str(x))
    is_train_orig = np.zeros(n_orig, dtype=bool)
    is_train_orig[train_idx] = True
    df["_is_train"] = is_train_orig[keep]

    duplicate_mask = df.duplicated(subset=[gid_col, env_col, trait_col], keep=False)
    if bool(duplicate_mask.any()):
        dup = df.loc[duplicate_mask, [gid_col, env_col, trait_col]].head(5)
        raise ValueError(
            "fit_multi_trait_gp_met requires one row per (gid, env, trait) cell; "
            f"duplicate examples: {dup.to_dict(orient='records')}"
        )
    y_all = df[y_col].to_numpy(dtype=np.float64)
    train_y = y_all[df["_is_train"].to_numpy(dtype=bool)]
    if train_y.size == 0:
        raise ValueError("train_idx selects no rows")
    if not np.isfinite(train_y).all():
        raise ValueError(
            "Training rows for fit_multi_trait_gp_met must have finite y values; "
            "missing y is only allowed on prediction/test rows."
        )
    df = df.reset_index(drop=True)

    trait_levels = sorted(df[trait_col].unique().tolist())
    env_levels = sorted(df[env_col].unique().tolist())
    T = len(trait_levels)
    n_env = len(env_levels)
    if T < 2:
        raise ValueError(f"MT-GP MET requires T >= 2 traits; got {T}")
    if n_env < 2:
        raise ValueError(f"MT-GP MET requires n_env >= 2 envs; got {n_env}")
    trait_to_idx = {t: i for i, t in enumerate(trait_levels)}
    env_to_idx = {e: i for i, e in enumerate(env_levels)}
    gid_to_idx = {str(g): i for i, g in enumerate(geno_ids)}
    n_gid = len(geno_ids)

    df["_trait_idx"] = df[trait_col].map(trait_to_idx).astype(np.int64)
    df["_env_idx"] = df[env_col].map(env_to_idx).astype(np.int64)
    if not df[gid_col].map(lambda g: g in gid_to_idx).all():
        missing = df.loc[~df[gid_col].isin(gid_to_idx), gid_col].unique()
        raise ValueError(f"gids in pheno_df not found in geno_ids: {missing.tolist()[:5]}")
    df["_gid_idx"] = df[gid_col].map(gid_to_idx).astype(np.int64)

    # Per-(env, trait) standardization on TRAINING rows. Cells with sparse or
    # no training data fall back to trait-level training scales, which keeps
    # CV0/new-env predictions on the original response scale after inversion.
    train_only = df.loc[df["_is_train"]]
    trained_traits = set(train_only[trait_col].unique().tolist())
    missing_train_traits = [t for t in trait_levels if t not in trained_traits]
    if missing_train_traits:
        raise ValueError(
            "Every trait must have at least one training row; traits missing from "
            f"training: {missing_train_traits[:5]}"
        )
    y_train_all = train_only[y_col].to_numpy(dtype=np.float64)
    finite_train = y_train_all[np.isfinite(y_train_all)]
    global_mean = float(np.mean(finite_train)) if finite_train.size else 0.0
    if finite_train.size >= 2:
        global_sd = float(np.std(finite_train, ddof=1))
        if not np.isfinite(global_sd) or global_sd <= 1e-12:
            global_sd = 1.0
    else:
        global_sd = 1.0

    trait_means = np.full(T, global_mean, dtype=np.float64)
    trait_sds = np.full(T, global_sd, dtype=np.float64)
    for t in range(T):
        sl_t = train_only.loc[train_only["_trait_idx"] == t, y_col].to_numpy(dtype=np.float64)
        sl_t = sl_t[np.isfinite(sl_t)]
        if sl_t.size:
            trait_means[t] = float(np.mean(sl_t))
        if sl_t.size >= 2:
            sd_t = float(np.std(sl_t, ddof=1))
            trait_sds[t] = sd_t if np.isfinite(sd_t) and sd_t > 1e-12 else global_sd

    trait_env_means = np.tile(trait_means[None, :], (n_env, 1))
    trait_env_sds = np.tile(trait_sds[None, :], (n_env, 1))
    for e in range(n_env):
        for t in range(T):
            sl = train_only.loc[
                (train_only["_env_idx"] == e) & (train_only["_trait_idx"] == t),
                y_col,
            ].to_numpy(dtype=np.float64)
            sl = sl[np.isfinite(sl)]
            if sl.size < 2:
                continue
            trait_env_means[e, t] = float(np.mean(sl))
            sd = float(np.std(sl, ddof=1))
            trait_env_sds[e, t] = sd if np.isfinite(sd) and sd > 1e-12 else trait_sds[t]

    env_idx_arr = df["_env_idx"].to_numpy()
    trait_idx_arr = df["_trait_idx"].to_numpy()
    y_std = (df[y_col].to_numpy(dtype=np.float64)
             - trait_env_means[env_idx_arr, trait_idx_arr]
             ) / trait_env_sds[env_idx_arr, trait_idx_arr]

    # X_fixed: user fixed effects only. Per-(env, trait) means are removed by
    # standardization and restored during back-transform, so explicit dense
    # intercept columns are redundant and too expensive for large MET panels.
    n_rows = len(df)
    blocks = []
    for col in fixed_effects:
        if col not in df.columns:
            raise ValueError(f"fixed_effects col {col!r} not in pheno_df")
        levels = sorted(df[col].astype(str).unique().tolist())
        if len(levels) < 2:
            continue
        col_oh = pd.get_dummies(df[col].astype(str)).to_numpy(dtype=np.float64)
        for c in range(col_oh.shape[1]):
            block = np.zeros((n_rows, T), dtype=np.float64)
            for t in range(T):
                m = trait_idx_arr == t
                block[m, t] = col_oh[m, c]
            blocks.append(block)
    X_fixed = np.concatenate(blocks, axis=1) if blocks else np.zeros((n_rows, 0), dtype=np.float64)

    train_mask = df["_is_train"].to_numpy(dtype=bool)
    test_mask = ~train_mask

    return {
        "T": T,
        "n_gid": n_gid,
        "n_env": n_env,
        "gid_idx": df["_gid_idx"].to_numpy(dtype=np.int64),
        "env_idx": env_idx_arr.astype(np.int64),
        "trait_idx": trait_idx_arr.astype(np.int64),
        "y_std": y_std,
        "X_fixed": X_fixed,
        "train_mask": train_mask,
        "test_mask": test_mask,
        "trait_levels": trait_levels,
        "env_levels": env_levels,
        "trait_env_means": trait_env_means,
        "trait_env_sds": trait_env_sds,
        "trait_means": trait_means,
        "trait_sds": trait_sds,
        "geno_ids": list(geno_ids),
    }


def _detect_I_full(
    gid_idx: np.ndarray, env_idx: np.ndarray, trait_idx: np.ndarray,
    n_gid: int, n_env: int, T: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Find the largest balanced sub-grid (gid_subset × env_subset) where every
    (gid, env) pair has all T traits observed in training rows.

    Greedy: start with all (gid, env) pairs that have T traits; iteratively drop
    the gid with fewest envs (or env with fewest gids) until balanced. Tie-break
    by ascending index.

    Returns (gid_subset, env_subset) as int64 arrays.
    """
    presence = np.zeros((n_gid, n_env), dtype=np.int64)
    for r in range(len(gid_idx)):
        presence[gid_idx[r], env_idx[r]] += 1
    cell_full = (presence == T)
    active_gids = set(int(g) for g in np.where(cell_full.any(axis=1))[0])
    active_envs = set(int(e) for e in np.where(cell_full.any(axis=0))[0])

    while True:
        if not active_gids or not active_envs:
            break
        gid_arr = sorted(active_gids); env_arr = sorted(active_envs)
        sub = cell_full[np.ix_(gid_arr, env_arr)]
        gid_counts = sub.sum(axis=1)
        env_counts = sub.sum(axis=0)
        n_active_envs = len(env_arr); n_active_gids = len(gid_arr)
        if (gid_counts == n_active_envs).all() and (env_counts == n_active_gids).all():
            break
        min_gid_count = gid_counts.min(); min_env_count = env_counts.min()
        if min_gid_count <= min_env_count:
            drop_gid = int(gid_arr[np.argmin(gid_counts)])
            active_gids.discard(drop_gid)
        else:
            drop_env = int(env_arr[np.argmin(env_counts)])
            active_envs.discard(drop_env)

    return (np.array(sorted(active_gids), dtype=np.int64),
            np.array(sorted(active_envs), dtype=np.int64))


def _estimate_sigma_G_GE_eps_mom(
    *,
    Phi: np.ndarray,
    K_env: np.ndarray,
    y_std: np.ndarray,
    X_fixed: np.ndarray,
    gid_idx: np.ndarray,
    env_idx: np.ndarray,
    trait_idx: np.ndarray,
    train_mask: np.ndarray,
    T: int,
    n_gid: int,
    n_env: int,
    psd_floor: float = 1e-4,
) -> Dict[str, Any]:
    """Three-moment MoM on a balanced or partially observed MET grid.

    Returns dict: Sigma_G, Sigma_GE, Sigma_eps, beta, mom_n_full_gids,
    mom_n_full_envs, mom_n_full_cells, mom_psd_clip_count, mom_3x3_cond.
    """
    from mt_gp_mom import _psd_project

    g_tr = gid_idx[train_mask]; e_tr = env_idx[train_mask]; t_tr = trait_idx[train_mask]
    y_tr = y_std[train_mask]; X_tr = X_fixed[train_mask]

    gid_subset, env_subset = _detect_I_full(g_tr, e_tr, t_tr, n_gid, n_env, T)
    n_full_gids = int(gid_subset.size); n_full_envs = int(env_subset.size)
    n_full_cells = n_full_gids * n_full_envs

    if n_full_gids < 2 or n_full_envs < 2:
        # CV2 masks genotype-environment cells across all traits. A complete
        # genotype x environment rectangle may therefore not remain even
        # though the three covariance operators are identifiable. Estimate
        # the same moments on all fully observed *cells*, using the exact
        # Frobenius products of the covariance bases restricted to that
        # incomplete grid. This is a genuine observed-grid MoM estimator, not
        # a predictive-variance proxy or a relabelled REML result.
        if X_tr.shape[1]:
            if X_tr.shape[0] <= X_tr.shape[1]:
                raise ValueError(
                    "Observed-grid MoM has too few observations for fixed effects: "
                    f"n_obs={X_tr.shape[0]}, n_fixed={X_tr.shape[1]}. Use bounded "
                    "varcomp_mode='reml' with a sufficiently estimable training fold."
                )
            beta_hat, *_ = np.linalg.lstsq(X_tr, y_tr, rcond=None)
            residual = y_tr - X_tr @ beta_hat
        else:
            beta_hat = np.zeros(0, dtype=np.float64)
            residual = y_tr

        cell_values: Dict[Tuple[int, int], np.ndarray] = {}
        for value, gid_value, env_value, trait_value in zip(
            residual, g_tr, e_tr, t_tr
        ):
            key = (int(gid_value), int(env_value))
            if key not in cell_values:
                cell_values[key] = np.full(T, np.nan, dtype=np.float64)
            cell_values[key][int(trait_value)] = float(value)
        common = [
            (key, values) for key, values in cell_values.items()
            if bool(np.isfinite(values).all())
        ]
        if len(common) < max(3, T + 1):
            raise ValueError(
                "Observed-grid MoM requires at least "
                f"{max(3, T + 1)} genotype-environment cells observed for every "
                f"trait; got {len(common)}. Use bounded varcomp_mode='reml' for "
                "a small partial-trait training panel."
            )

        cell_gid = np.asarray([key[0] for key, _ in common], dtype=np.int64)
        cell_env = np.asarray([key[1] for key, _ in common], dtype=np.int64)
        Z = np.vstack([values for _, values in common])
        # The public design has already removed training environment-by-trait
        # means. Re-center the retained common cells within each environment
        # to avoid a partial-panel mean leaking into the covariance moments.
        for env_value in np.unique(cell_env):
            env_rows = cell_env == env_value
            Z[env_rows, :] -= Z[env_rows, :].mean(axis=0, keepdims=True)

        env_present = np.unique(cell_env)
        S_env = []
        U_env = []
        for env_value in env_present:
            env_rows = cell_env == env_value
            Phi_env = Phi[cell_gid[env_rows], :]
            Z_env = Z[env_rows, :]
            S_env.append(Phi_env.T @ Phi_env)
            U_env.append(Phi_env.T @ Z_env)

        S_total = np.sum(np.stack(S_env, axis=0), axis=0)
        U_total = np.sum(np.stack(U_env, axis=0), axis=0)
        H = np.empty((len(env_present), len(env_present)), dtype=np.float64)
        for left in range(len(env_present)):
            for right in range(left, len(env_present)):
                value = float(np.sum(S_env[left] * S_env[right]))
                H[left, right] = value
                H[right, left] = value
        K_env_obs = K_env[np.ix_(env_present, env_present)]
        norm_G2 = float(np.sum(S_total * S_total))
        inner_G_GE = float(np.sum(K_env_obs * H))
        norm_GE2 = float(np.sum((K_env_obs * K_env_obs) * H))
        phi_diag = np.sum(Phi[cell_gid, :] ** 2, axis=1)
        tr_G = float(phi_diag.sum())
        tr_GE = float(np.sum(phi_diag * K_env[cell_env, cell_env]))
        n_cells_observed = int(len(common))
        M = np.array([
            [norm_G2, inner_G_GE, tr_G],
            [inner_G_GE, norm_GE2, tr_GE],
            [tr_G, tr_GE, float(n_cells_observed)],
        ], dtype=np.float64)
        ridge_scale = max(float(np.max(np.diag(M))), 1.0)
        M_reg = M + (1e-10 * ridge_scale) * np.eye(3)
        cond_M = float(np.linalg.cond(M_reg))
        fallback_no_gxe = (not np.isfinite(cond_M)) or cond_M > 1e10
        if fallback_no_gxe:
            warnings.warn(
                "The observed CV training grid cannot separately identify the "
                "GxE covariance from the main genetic and residual covariance; "
                "fitting the identifiable no-GxE MoM submodel.",
                RuntimeWarning,
                stacklevel=2,
            )
            M2 = np.array([
                [norm_G2, tr_G],
                [tr_G, float(n_cells_observed)],
            ], dtype=np.float64)
            M2_reg = M2 + (1e-10 * max(float(np.max(np.diag(M2))), 1.0)) * np.eye(2)
            cond_M2 = float(np.linalg.cond(M2_reg))
            if (not np.isfinite(cond_M2)) or cond_M2 > 1e12:
                raise ValueError(
                    "Observed-grid MoM cannot identify genetic and residual "
                    "covariance from this training fold."
                )

        Q_G = U_total.T @ U_total
        Q_GE = np.zeros((T, T), dtype=np.float64)
        for left in range(len(env_present)):
            for right in range(len(env_present)):
                Q_GE += K_env_obs[left, right] * (U_env[left].T @ U_env[right])
        Q_E = Z.T @ Z

        Sigma_G = np.zeros((T, T), dtype=np.float64)
        Sigma_GE = np.zeros((T, T), dtype=np.float64)
        Sigma_eps = np.zeros((T, T), dtype=np.float64)
        for left in range(T):
            for right in range(left, T):
                if fallback_no_gxe:
                    solution = np.linalg.solve(
                        M2_reg,
                        np.array([Q_G[left, right], Q_E[left, right]], dtype=np.float64),
                    )
                    values = (solution[0], 0.0, solution[1])
                else:
                    solution = np.linalg.solve(
                        M_reg,
                        np.array([
                            Q_G[left, right],
                            Q_GE[left, right],
                            Q_E[left, right],
                        ], dtype=np.float64),
                    )
                    values = tuple(solution.tolist())
                Sigma_G[left, right] = Sigma_G[right, left] = values[0]
                Sigma_GE[left, right] = Sigma_GE[right, left] = values[1]
                Sigma_eps[left, right] = Sigma_eps[right, left] = values[2]

        Sigma_G, n_clip_G = _psd_project(Sigma_G, floor=psd_floor, name="Sigma_G")
        if fallback_no_gxe:
            Sigma_GE = np.zeros((T, T), dtype=np.float64)
            n_clip_GE = 0
        else:
            Sigma_GE, n_clip_GE = _psd_project(
                Sigma_GE, floor=psd_floor, name="Sigma_GE"
            )
        Sigma_eps, n_clip_E = _psd_project(
            Sigma_eps, floor=psd_floor, name="Sigma_eps"
        )
        return {
            "Sigma_G": Sigma_G,
            "Sigma_GE": Sigma_GE,
            "Sigma_eps": Sigma_eps,
            "beta": beta_hat,
            "mom_n_full_gids": n_full_gids,
            "mom_n_full_envs": n_full_envs,
            "mom_n_full_cells": n_full_cells,
            "mom_n_observed_cells": n_cells_observed,
            "mom_psd_clip_count": (int(n_clip_G), int(n_clip_GE), int(n_clip_E)),
            "mom_3x3_cond": cond_M,
            "mom_fallback": (
                "partial_observed_grid_no_gxe"
                if fallback_no_gxe else "partial_observed_grid_mom"
            ),
        }

    in_gid = np.zeros(n_gid, dtype=bool); in_gid[gid_subset] = True
    in_env = np.zeros(n_env, dtype=bool); in_env[env_subset] = True
    in_full = in_gid[g_tr] & in_env[e_tr]
    g_full = g_tr[in_full]; e_full = e_tr[in_full]; t_full = t_tr[in_full]
    y_full = y_tr[in_full]; X_full = X_tr[in_full]

    if X_full.shape[1]:
        if X_full.shape[0] <= X_full.shape[1]:
            raise ValueError(
                f"MoM balanced sub-grid has too few observations for fixed effects: "
                f"n_obs={X_full.shape[0]}, n_fixed={X_full.shape[1]}"
            )
        beta_hat, *_ = np.linalg.lstsq(X_full, y_full, rcond=None)
        r = y_full - X_full @ beta_hat
    else:
        beta_hat = np.zeros(0, dtype=np.float64)
        r = y_full

    gid_to_pos = -np.ones(n_gid, dtype=np.int64); gid_to_pos[gid_subset] = np.arange(n_full_gids)
    env_to_pos = -np.ones(n_env, dtype=np.int64); env_to_pos[env_subset] = np.arange(n_full_envs)
    Z = np.zeros((n_full_gids, n_full_envs, T), dtype=np.float64)
    for k in range(len(r)):
        Z[gid_to_pos[g_full[k]], env_to_pos[e_full[k]], t_full[k]] = r[k]
    # Equivalent to per-(env, trait) intercept residualization on the balanced
    # MoM grid, without materializing dense intercept columns for the full MET
    # panel. This matters for CV2/partial-cell designs where the I_full subset
    # can have slightly different env-trait means than all training rows.
    Z = Z - Z.mean(axis=0, keepdims=True)

    Phi_full = Phi[gid_subset, :]
    K_env_full = K_env[np.ix_(env_subset, env_subset)]

    A_g = Phi_full.T @ Phi_full
    tr_Kg = float(np.sum(Phi_full * Phi_full))
    norm_Kg2 = float(np.sum(A_g * A_g))
    tr_Ke = float(np.trace(K_env_full))
    norm_Ke2 = float(np.sum(K_env_full * K_env_full))
    sum_Ke = float(K_env_full.sum())
    M = np.array([
        [norm_Kg2 * sum_Ke,    norm_Kg2 * norm_Ke2,   tr_Kg * tr_Ke],
        [norm_Kg2 * n_full_envs, norm_Kg2 * tr_Ke,    tr_Kg * n_full_envs],
        [tr_Kg * n_full_envs,  tr_Kg * tr_Ke,         float(n_full_cells)],
    ], dtype=np.float64)
    M_reg = M + 1e-10 * np.eye(3)
    cond_M = float(np.linalg.cond(M_reg))
    fallback_no_gxe = cond_M > 1e10
    if fallback_no_gxe:
        warnings.warn(
            f"K_env has near-trivial structure (cond(M_3x3) = {cond_M:.2e}); "
            "Sigma_GE is unidentifiable from Sigma_eps. Proceeding with a "
            "no-GxE fallback (Sigma_GE set to zero) for this small or "
            "environment-kernel-limited MT-MET fit. Supply env_covariates or "
            "a richer env_similarity to estimate the GxE covariance.",
            RuntimeWarning,
            stacklevel=2,
        )
        M2 = np.array([
            [norm_Kg2 * n_full_envs, tr_Kg * n_full_envs],
            [tr_Kg * n_full_envs,    float(n_full_cells)],
        ], dtype=np.float64)
        M2_reg = M2 + 1e-10 * np.eye(2)

    Sigma_G = np.zeros((T, T), dtype=np.float64)
    Sigma_GE = np.zeros((T, T), dtype=np.float64)
    Sigma_eps = np.zeros((T, T), dtype=np.float64)

    PhiT_Z = np.einsum("gm, geT -> meT", Phi_full, Z)
    K_geno_Z = np.einsum("gm, meT -> geT", Phi_full, PhiT_Z)
    Ke_Z = np.einsum("ef, gft -> get", K_env_full, Z)
    KK_Z = np.einsum("gm, mef -> gef", Phi_full,
                     np.einsum("gm, gef -> mef", Phi_full, Ke_Z))

    for a in range(T):
        for b in range(a, T):
            m1 = float(np.einsum("ge, ge ->", Z[:, :, a], KK_Z[:, :, b]))
            m2 = float(np.einsum("ge, ge ->", Z[:, :, a], K_geno_Z[:, :, b]))
            m3 = float(np.einsum("ge, ge ->", Z[:, :, a], Z[:, :, b]))
            if fallback_no_gxe:
                rhs = np.array([m2, m3], dtype=np.float64)
                sol = np.linalg.solve(M2_reg, rhs)
                Sigma_G[a, b] = Sigma_G[b, a] = sol[0]
                Sigma_GE[a, b] = Sigma_GE[b, a] = 0.0
                Sigma_eps[a, b] = Sigma_eps[b, a] = sol[1]
            else:
                rhs = np.array([m1, m2, m3], dtype=np.float64)
                sol = np.linalg.solve(M_reg, rhs)
                Sigma_G[a, b] = Sigma_G[b, a] = sol[0]
                Sigma_GE[a, b] = Sigma_GE[b, a] = sol[1]
                Sigma_eps[a, b] = Sigma_eps[b, a] = sol[2]

    Sigma_G, n_clip_G = _psd_project(Sigma_G, floor=psd_floor, name="Sigma_G")
    if fallback_no_gxe:
        Sigma_GE = np.zeros((T, T), dtype=np.float64)
        n_clip_GE = 0
    else:
        Sigma_GE, n_clip_GE = _psd_project(Sigma_GE, floor=psd_floor, name="Sigma_GE")
    Sigma_eps, n_clip_E = _psd_project(Sigma_eps, floor=psd_floor, name="Sigma_eps")

    return {
        "Sigma_G":    Sigma_G,
        "Sigma_GE":   Sigma_GE,
        "Sigma_eps":  Sigma_eps,
        "beta":       beta_hat,
        "mom_n_full_gids":    n_full_gids,
        "mom_n_full_envs":    n_full_envs,
        "mom_n_full_cells":   n_full_cells,
        "mom_psd_clip_count": (int(n_clip_G), int(n_clip_GE), int(n_clip_E)),
        "mom_3x3_cond":       cond_M,
        "mom_fallback":        "no_gxe_near_trivial_env_kernel" if fallback_no_gxe else None,
    }


def _fit_dense_reml_met_variance(
    *,
    K_geno: np.ndarray,
    K_env: np.ndarray,
    y_std: np.ndarray,
    X_fixed: np.ndarray,
    gid_idx: np.ndarray,
    env_idx: np.ndarray,
    trait_idx: np.ndarray,
    train_mask: np.ndarray,
    Sigma_G_init: np.ndarray,
    Sigma_GE_init: np.ndarray,
    Sigma_eps_init: np.ndarray,
    psd_floor: float = 1e-4,
    max_iter: int = 30,
    tol: float = 1e-5,
    device: Optional[Any] = None,
) -> Dict[str, Any]:
    """Small/medium dense REML optimizer for MET variance components.

    This is intentionally bounded by the caller. It estimates variance
    components only; prediction remains in the GP/operator path.
    """
    from mt_gp_internals import _pack_unstructured, _unpack_unstructured

    device = pick_torch_device(device)
    T = int(Sigma_G_init.shape[0])
    n_par = T * (T + 1) // 2
    eye_T = torch.eye(T, dtype=torch.float64, device=device)

    def _pack_init(M: np.ndarray) -> torch.Tensor:
        M = 0.5 * (np.asarray(M, dtype=np.float64) + np.asarray(M, dtype=np.float64).T)
        min_eval = float(np.linalg.eigvalsh(M).min())
        jitter = max(float(psd_floor), -min_eval + float(psd_floor))
        L = torch.linalg.cholesky(torch.tensor(M, dtype=torch.float64, device=device) + jitter * eye_T)
        return _pack_unstructured(L)

    theta = torch.cat([
        _pack_init(Sigma_G_init),
        _pack_init(Sigma_GE_init),
        _pack_init(Sigma_eps_init),
    ]).detach().clone().requires_grad_(True)

    y = torch.tensor(y_std[train_mask], dtype=torch.float64, device=device)
    X = torch.tensor(X_fixed[train_mask], dtype=torch.float64, device=device)
    g = torch.tensor(gid_idx[train_mask], dtype=torch.int64, device=device)
    e = torch.tensor(env_idx[train_mask], dtype=torch.int64, device=device)
    t = torch.tensor(trait_idx[train_mask], dtype=torch.int64, device=device)
    Kg = torch.tensor(K_geno, dtype=torch.float64, device=device)
    Ke = torch.tensor(K_env, dtype=torch.float64, device=device)
    p = int(X.shape[1])
    jitter_I = 1e-6 * torch.eye(int(y.numel()), dtype=torch.float64, device=device)

    def _unpack_all(theta_v: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        Lg = _unpack_unstructured(theta_v[:n_par], T)
        Lge = _unpack_unstructured(theta_v[n_par:2 * n_par], T)
        Le = _unpack_unstructured(theta_v[2 * n_par:], T)
        floor = float(psd_floor)
        return (
            Lg @ Lg.T + floor * eye_T,
            Lge @ Lge.T + floor * eye_T,
            Le @ Le.T + floor * eye_T,
        )

    def _nll(theta_v: torch.Tensor) -> torch.Tensor:
        Sigma_G, Sigma_GE, Sigma_eps = _unpack_all(theta_v)
        Kg_pairs = Kg[g.unsqueeze(1), g.unsqueeze(0)]
        Ke_pairs = Ke[e.unsqueeze(1), e.unsqueeze(0)]
        G_pairs = Sigma_G[t.unsqueeze(1), t.unsqueeze(0)]
        GE_pairs = Sigma_GE[t.unsqueeze(1), t.unsqueeze(0)]
        E_pairs = Sigma_eps[t.unsqueeze(1), t.unsqueeze(0)]
        same_cell = ((g.unsqueeze(1) == g.unsqueeze(0)) &
                     (e.unsqueeze(1) == e.unsqueeze(0))).to(torch.float64)
        V = G_pairs * Kg_pairs + GE_pairs * Kg_pairs * Ke_pairs + E_pairs * same_cell + jitter_I
        L = torch.linalg.cholesky(V)
        log_det_V = 2.0 * torch.log(torch.diagonal(L)).sum()
        if p:
            Vinv_y = torch.cholesky_solve(y.unsqueeze(1), L).squeeze(1)
            Vinv_X = torch.cholesky_solve(X, L)
            Xt_Vinv_X = X.T @ Vinv_X
            L_xvx = torch.linalg.cholesky(
                Xt_Vinv_X + 1e-6 * torch.eye(p, dtype=torch.float64, device=device)
            )
            log_det_xvx = 2.0 * torch.log(torch.diagonal(L_xvx)).sum()
            beta_hat = torch.cholesky_solve((X.T @ Vinv_y).unsqueeze(1), L_xvx).squeeze(1)
            r = y - X @ beta_hat
        else:
            log_det_xvx = torch.zeros((), dtype=torch.float64, device=device)
            r = y
        Vinv_r = torch.cholesky_solve(r.unsqueeze(1), L).squeeze(1)
        quad = (r * Vinv_r).sum()
        return 0.5 * (log_det_V + log_det_xvx + quad)

    history = []
    optimizer = torch.optim.LBFGS(
        [theta],
        max_iter=int(max_iter),
        tolerance_grad=float(tol),
        tolerance_change=float(tol),
        line_search_fn="strong_wolfe",
    )

    def closure():
        optimizer.zero_grad()
        loss = _nll(theta)
        loss.backward()
        history.append(float(-loss.detach().item()))
        return loss

    optimizer.step(closure)
    with torch.no_grad():
        Sigma_G, Sigma_GE, Sigma_eps = _unpack_all(theta)
        nll_final = _nll(theta)
        Kg_pairs = Kg[g.unsqueeze(1), g.unsqueeze(0)]
        Ke_pairs = Ke[e.unsqueeze(1), e.unsqueeze(0)]
        G_pairs = Sigma_G[t.unsqueeze(1), t.unsqueeze(0)]
        GE_pairs = Sigma_GE[t.unsqueeze(1), t.unsqueeze(0)]
        E_pairs = Sigma_eps[t.unsqueeze(1), t.unsqueeze(0)]
        same_cell = ((g.unsqueeze(1) == g.unsqueeze(0)) &
                     (e.unsqueeze(1) == e.unsqueeze(0))).to(torch.float64)
        V = G_pairs * Kg_pairs + GE_pairs * Kg_pairs * Ke_pairs + E_pairs * same_cell + jitter_I
        L = torch.linalg.cholesky(V)
        if p:
            Vinv_y = torch.cholesky_solve(y.unsqueeze(1), L).squeeze(1)
            Vinv_X = torch.cholesky_solve(X, L)
            Xt_Vinv_X = X.T @ Vinv_X
            L_xvx = torch.linalg.cholesky(
                Xt_Vinv_X + 1e-6 * torch.eye(p, dtype=torch.float64, device=device)
            )
            beta_hat = torch.cholesky_solve((X.T @ Vinv_y).unsqueeze(1), L_xvx).squeeze(1)
            beta_np = torch_to_numpy(beta_hat)
        else:
            beta_np = np.zeros(0, dtype=np.float64)

    # `history` holds every objective evaluation, including strong-Wolfe
    # line-search trial points, so its last two entries are not consecutive
    # iterates and cannot decide convergence. Use the optimizer's own state:
    # L-BFGS stopped before max_iter (its gradient/change tolerances were met),
    # or the gradient at the solution is below tolerance.
    hist = np.asarray(history, dtype=np.float64)
    lbfgs_state = optimizer.state.get(theta, {})
    n_iter_done = int(lbfgs_state.get("n_iter", 0))
    theta_check = theta.detach().clone().requires_grad_(True)
    grad = torch.autograd.grad(_nll(theta_check), theta_check)[0]
    grad_max = float(grad.abs().max().item())
    converged = bool(n_iter_done < int(max_iter) or grad_max <= float(tol))
    return {
        "Sigma_G": torch_to_numpy(Sigma_G),
        "Sigma_GE": torch_to_numpy(Sigma_GE),
        "Sigma_eps": torch_to_numpy(Sigma_eps),
        "beta": beta_np,
        "loglik_history": hist,
        "n_iter": n_iter_done,
        "n_function_evals": int(hist.size),
        "grad_max_abs": grad_max,
        "converged": converged,
        "neg_reml": float(nll_final.detach().item()),
        "device": str(device),
    }


def _fit_mt_met_mom_op(
    *,
    design: Dict[str, Any],
    cache: Any,
    K_env: np.ndarray,
    psd_floor: float = 1e-4,
    pcg_tol: float = 1e-4,
    pcg_max_iter: int = 300,
    return_cov_traits: bool = False,
    force_cov_traits: bool = False,
    return_se: bool = False,
    force_prediction_se: bool = False,
    trait_structure: str = "unstructured",
    trait_fa_rank: Optional[int] = None,
    gxe_trait_structure: Optional[str] = None,
    gxe_trait_fa_rank: Optional[int] = None,
    kernel_weights: Optional[Sequence[float]] = None,
    dense_reml: bool = False,
    dense_reml_max_train: int = 3000,
    reml_max_iter: int = 30,
    reml_tol: float = 1e-5,
    device: Optional[Any] = None,
) -> Dict[str, Any]:
    """Multi-env MoM orchestrator: MoM Σs → operators → PCG → predict."""
    import pandas as pd
    from mt_gp_mom import _psd_project
    from mt_gp_operators import MTGenoOperator, MTGxEOperator, MTMETResidualOperator, _MTSumOperator
    from mixed_model_gpu_large_scale_backend import pcg_solve

    device = pick_torch_device(device)
    T = design["T"]; n_gid = design["n_gid"]; n_env = design["n_env"]
    if cache.n_geno() != n_gid:
        raise ValueError(
            f"cache.n_geno()={cache.n_geno()} != n_gid={n_gid}; geno_ids must align with cache."
        )
    kernel_weights_arr = _normalize_kernel_weights(cache, kernel_weights)
    Phi = _materialize_weighted_phi(cache, n_gid, kernel_weights_arr)
    if K_env.shape != (n_env, n_env):
        raise ValueError(f"K_env shape {K_env.shape} != ({n_env}, {n_env})")

    # MoM is the scalable estimator and also supplies dense-REML start values
    # when a balanced training sub-grid exists. CV2 deliberately masks a cell
    # in every genotype, so a balanced sub-grid is not guaranteed. In bounded
    # dense-REML mode only, fall back to a conservative covariance split based
    # on all observed training cells; REML, not this initializer, remains the
    # reported estimator. Never use this fallback for a requested MoM fit.
    try:
        mom_out = _estimate_sigma_G_GE_eps_mom(
            Phi=Phi, K_env=K_env,
            y_std=design["y_std"], X_fixed=design["X_fixed"],
            gid_idx=design["gid_idx"], env_idx=design["env_idx"],
            trait_idx=design["trait_idx"],
            train_mask=design["train_mask"],
            T=T, n_gid=n_gid, n_env=n_env,
            psd_floor=psd_floor,
        )
    except ValueError as exc:
        initializer_shortage = any(
            token in str(exc)
            for token in (
                "MoM requires balanced sub-grid",
                "Observed-grid MoM requires at least",
            )
        )
        if not bool(dense_reml) or not initializer_shortage:
            raise
        train_mask_np = np.asarray(design["train_mask"], dtype=bool)
        y_train = np.asarray(design["y_std"], dtype=np.float64)[train_mask_np]
        X_train = np.asarray(design["X_fixed"], dtype=np.float64)[train_mask_np]
        trait_train = np.asarray(design["trait_idx"], dtype=np.int64)[train_mask_np]
        gid_train = np.asarray(design["gid_idx"], dtype=np.int64)[train_mask_np]
        env_train = np.asarray(design["env_idx"], dtype=np.int64)[train_mask_np]
        beta_init = (
            np.linalg.lstsq(X_train, y_train, rcond=None)[0]
            if X_train.shape[1]
            else np.zeros(0, dtype=np.float64)
        )
        residual = y_train - X_train @ beta_init
        cell_values: Dict[Tuple[int, int], np.ndarray] = {}
        for value, gid_value, env_value, trait_value in zip(
            residual, gid_train, env_train, trait_train
        ):
            key = (int(gid_value), int(env_value))
            if key not in cell_values:
                cell_values[key] = np.full(T, np.nan, dtype=np.float64)
            cell_values[key][int(trait_value)] = float(value)
        cell_matrix = np.vstack(list(cell_values.values()))
        total_covariance = np.zeros((T, T), dtype=np.float64)
        for left in range(T):
            for right in range(left, T):
                paired = np.isfinite(cell_matrix[:, left]) & np.isfinite(cell_matrix[:, right])
                if int(paired.sum()) >= 2:
                    if left == right:
                        value = float(np.var(cell_matrix[paired, left], ddof=1))
                    else:
                        value = float(np.cov(
                            cell_matrix[paired, left],
                            cell_matrix[paired, right],
                            ddof=1,
                        )[0, 1])
                else:
                    value = 1.0 if left == right else 0.0
                total_covariance[left, right] = value
                total_covariance[right, left] = value
        total_covariance, total_clip_count = _psd_project(
            total_covariance,
            floor=psd_floor,
            name="partial-grid total covariance initializer",
        )
        full_gid_subset, full_env_subset = _detect_I_full(
            gid_train, env_train, trait_train, n_gid, n_env, T
        )
        mom_out = {
            "Sigma_G": 0.35 * total_covariance,
            "Sigma_GE": 0.25 * total_covariance,
            "Sigma_eps": 0.40 * total_covariance,
            "beta": beta_init,
            "mom_n_full_gids": int(full_gid_subset.size),
            "mom_n_full_envs": int(full_env_subset.size),
            "mom_n_full_cells": int(full_gid_subset.size * full_env_subset.size),
            "mom_psd_clip_count": (int(total_clip_count), 0, 0),
            "mom_3x3_cond": np.nan,
            "mom_fallback": "dense_reml_partial_grid_initializer",
        }
    Sigma_G = mom_out["Sigma_G"]; Sigma_GE = mom_out["Sigma_GE"]
    Sigma_eps = mom_out["Sigma_eps"]; beta_mom = mom_out["beta"]
    reml_out = None
    beta_reml = None
    if bool(dense_reml):
        n_train = int(np.asarray(design["train_mask"], dtype=bool).sum())
        if n_train > int(dense_reml_max_train):
            raise ValueError(
                f"varcomp_mode='reml' dense MET is capped at dense_reml_max_train="
                f"{int(dense_reml_max_train)} training rows; got {n_train}. "
                "Use varcomp_mode='mom' for scalable MET fitting."
            )
        if str(trait_structure or "unstructured").lower() != "unstructured":
            raise ValueError("dense MET REML currently supports trait_structure='unstructured' only")
        if str(gxe_trait_structure or "unstructured").lower() != "unstructured":
            raise ValueError("dense MET REML currently supports gxe_trait_structure='unstructured' only")
        reml_out = _fit_dense_reml_met_variance(
            K_geno=Phi @ Phi.T,
            K_env=K_env,
            y_std=design["y_std"],
            X_fixed=design["X_fixed"],
            gid_idx=design["gid_idx"],
            env_idx=design["env_idx"],
            trait_idx=design["trait_idx"],
            train_mask=design["train_mask"],
            Sigma_G_init=Sigma_G,
            Sigma_GE_init=Sigma_GE,
            Sigma_eps_init=Sigma_eps,
            psd_floor=psd_floor,
            max_iter=int(reml_max_iter),
            tol=float(reml_tol),
            device=device,
        )
        Sigma_G = reml_out["Sigma_G"]
        Sigma_GE = reml_out["Sigma_GE"]
        Sigma_eps = reml_out["Sigma_eps"]
        beta_reml = reml_out["beta"]
    trait_structure = str(trait_structure or "unstructured").lower()
    gxe_trait_structure = str(gxe_trait_structure or trait_structure).lower()
    Lambda_G = psi_G = Lambda_GE = psi_GE = None
    if trait_structure == "fa":
        if trait_fa_rank is None:
            raise ValueError("trait_fa_rank required when trait_structure='fa'")
        Sigma_G, Lambda_G, psi_G = _fa_project_covariance(
            Sigma_G, int(trait_fa_rank), floor=psd_floor, name="Sigma_G"
        )
    elif trait_structure != "unstructured":
        raise ValueError(f"unknown trait_structure {trait_structure!r}")
    if gxe_trait_structure == "fa":
        rank = gxe_trait_fa_rank if gxe_trait_fa_rank is not None else trait_fa_rank
        if rank is None:
            raise ValueError("gxe_trait_fa_rank required when gxe_trait_structure='fa'")
        Sigma_GE, Lambda_GE, psi_GE = _fa_project_covariance(
            Sigma_GE, int(rank), floor=psd_floor, name="Sigma_GE"
        )
    elif gxe_trait_structure != "unstructured":
        raise ValueError(f"unknown gxe_trait_structure {gxe_trait_structure!r}")
    Sigma_eps, _ = _psd_project(Sigma_eps, floor=psd_floor, name="Sigma_eps")

    # Build operators on TRAINING rows
    train_mask = design["train_mask"]
    g_tr = design["gid_idx"][train_mask]; e_tr = design["env_idx"][train_mask]
    t_tr = design["trait_idx"][train_mask]
    y_tr_std = design["y_std"][train_mask]; X_tr = design["X_fixed"][train_mask]

    op_g = MTGenoOperator(
        cache,
        gi=torch.tensor(g_tr, dtype=torch.int64, device=device),
        ti=torch.tensor(t_tr, dtype=torch.int64, device=device),
        Sigma_G=torch.tensor(Sigma_G, dtype=torch.float64, device=device),
        kernel_weights=torch.tensor(kernel_weights_arr, dtype=torch.float64, device=device),
    )
    op_ge = MTGxEOperator(
        cache,
        gi=torch.tensor(g_tr, dtype=torch.int64, device=device),
        ei=torch.tensor(e_tr, dtype=torch.int64, device=device),
        ti=torch.tensor(t_tr, dtype=torch.int64, device=device),
        Sigma_GE=torch.tensor(Sigma_GE, dtype=torch.float64, device=device),
        K_env=torch.tensor(K_env, dtype=torch.float64, device=device),
        kernel_weights=torch.tensor(kernel_weights_arr, dtype=torch.float64, device=device),
    )
    op_e = MTMETResidualOperator(
        gi=torch.tensor(g_tr, dtype=torch.int64, device=device),
        ei=torch.tensor(e_tr, dtype=torch.int64, device=device),
        ti=torch.tensor(t_tr, dtype=torch.int64, device=device),
        Sigma_eps=torch.tensor(Sigma_eps, dtype=torch.float64, device=device),
    )
    fold_train_g, ti_train_fold = op_g.build_fold_cache()
    ei_train_fold = torch.tensor(e_tr, dtype=torch.int64, device=device)
    fold_train_e, _, _ = op_e.build_fold_cache()
    t_tr_t = torch.tensor(t_tr, dtype=torch.int64, device=device)

    # Adapter for op_ge to fit _MTSumOperator's (v, fold, ti_fold) signature
    class _GxEAdapter:
        def __init__(self, op_ge, ei_fold):
            self._op = op_ge
            self._ei_fold = ei_fold
        def matvec(self, v, fold, ti_fold, dtype_compute=torch.float32):
            return self._op.matvec(v, fold, self._ei_fold, ti_fold, dtype_compute=dtype_compute)

    # Adapter for op_e: ignores the gid-level fold; uses its own cell-level fold
    class _ResidualAdapter:
        def __init__(self, op_e, fold_train_e):
            self._op = op_e
            self._fold_e = fold_train_e
        def matvec(self, v, fold, ti_fold, dtype_compute=torch.float32):
            # fold arg intentionally ignored; cell-level fold fixed at construction
            return self._op.matvec(v, self._fold_e, ti_fold, dtype_compute=dtype_compute)

    V_op = _MTSumOperator(op_g, _GxEAdapter(op_ge, ei_train_fold), _ResidualAdapter(op_e, fold_train_e))

    geno_diag_unique = torch.zeros(fold_train_g.ng_used, dtype=torch.float64, device=device)
    for r, w in enumerate(kernel_weights_arr.tolist()):
        if w <= 0.0:
            continue
        Phi_train_used_r = cache.get_rows(r, fold_train_g.g_unique, device=device, dtype=torch.float64)
        geno_diag_unique = geno_diag_unique + float(w) * (Phi_train_used_r * Phi_train_used_r).sum(dim=1)
    geno_diag = geno_diag_unique.index_select(0, fold_train_g.inv)
    K_env_diag = torch.diag(torch.tensor(K_env, dtype=torch.float64, device=device)).index_select(0, ei_train_fold)
    Sigma_G_diag = torch.diag(torch.tensor(Sigma_G, dtype=torch.float64, device=device)).index_select(0, t_tr_t)
    Sigma_GE_diag = torch.diag(torch.tensor(Sigma_GE, dtype=torch.float64, device=device)).index_select(0, t_tr_t)
    Sigma_eps_diag = torch.diag(torch.tensor(Sigma_eps, dtype=torch.float64, device=device)).index_select(0, t_tr_t)
    diagV = (Sigma_G_diag * geno_diag
             + Sigma_GE_diag * geno_diag * K_env_diag
             + Sigma_eps_diag).clamp_min(torch.tensor(1e-8, dtype=torch.float64, device=device))

    def M_inv_mv(v):
        return v / diagV.to(device=v.device, dtype=v.dtype)

    # Prediction fixed effects are refit on all training rows. The MoM beta above
    # is still the beta used to residualize the balanced sub-grid for variance
    # estimation; this all-train OLS beta avoids throwing away partial training
    # rows in the posterior mean path while preserving the fast MET v1 contract.
    if beta_reml is not None:
        beta = beta_reml
    elif X_tr.shape[1]:
        beta, *_ = np.linalg.lstsq(X_tr, y_tr_std, rcond=None)
    else:
        beta = np.zeros(0, dtype=np.float64)

    # Residualize
    r = torch.tensor(y_tr_std - X_tr @ beta, dtype=torch.float64, device=device)
    def A_mv(v):
        return V_op.matvec(v, fold_train_g, ti_train_fold, dtype_compute=torch.float64)
    pcg_res = pcg_solve(A_mv, r, M_inv_mv=M_inv_mv, tol=pcg_tol, max_iter=pcg_max_iter)
    alpha = pcg_res.x
    if not bool(pcg_res.converged):
        warnings.warn(
            "fit_multi_trait_gp_met PCG did not converge "
            f"(rel_resid={pcg_res.rel_resid:.2e}, max_iter={pcg_max_iter}); "
            "predictions may be unreliable. Increase mom_pcg_max_iter, relax mom_pcg_tol, "
            "or inspect K_env and variance-component conditioning.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Predict at test rows
    test_mask = design["test_mask"]
    g_te = design["gid_idx"][test_mask]; e_te = design["env_idx"][test_mask]
    t_te = design["trait_idx"][test_mask]; X_te = design["X_fixed"][test_mask]
    test_state_g = {
        "fold_train": fold_train_g, "ti_train": ti_train_fold,
        "gi_test": torch.tensor(g_te, dtype=torch.int64, device=device),
        "ti_test": torch.tensor(t_te, dtype=torch.int64, device=device),
    }
    test_state_ge = dict(test_state_g)
    test_state_ge["ei_train"] = ei_train_fold
    test_state_ge["ei_test"] = torch.tensor(e_te, dtype=torch.int64, device=device)

    mu_g_std = torch_to_numpy(op_g.cross_cov_apply(alpha, test_state_g, dtype_compute=torch.float64))
    mu_ge_std = torch_to_numpy(op_ge.cross_cov_apply(alpha, test_state_ge, dtype_compute=torch.float64))
    y_pred_std = X_te @ beta + mu_g_std + mu_ge_std

    # Back-transform per-(env, trait)
    trait_env_means = design["trait_env_means"]; trait_env_sds = design["trait_env_sds"]
    y_pred = (y_pred_std * trait_env_sds[e_te, t_te]
              + trait_env_means[e_te, t_te])
    trait_levels = design["trait_levels"]; env_levels = design["env_levels"]
    geno_ids = design["geno_ids"]
    preds_df = pd.DataFrame({
        "gid":   [geno_ids[int(g)] for g in g_te],
        "env":   [env_levels[int(e)] for e in e_te],
        "trait": [trait_levels[int(t)] for t in t_te],
        "Prediction": y_pred,
    })
    for col in ("gid", "env", "trait"):
        preds_df[col] = pd.Series([str(x) for x in preds_df[col].tolist()], index=preds_df.index, dtype=object)

    # cov_traits opt-in with cap
    cov_traits = None
    cov_trait_cells = None
    cov_traits_downgraded = False
    cov_traits_downgrade_reason = ""
    compute_cov_traits = False
    if return_cov_traits or return_se:
        from mt_gp_reporting import resolve_cov_traits_request

        n_test_cells = int(np.unique(np.column_stack([g_te, e_te]), axis=0).shape[0])
        decision = resolve_cov_traits_request(
            n_cells=n_test_cells,
            n_traits=T,
            requested=True,
            force_cov_traits=bool(force_cov_traits),
            force_prediction_se=bool(force_prediction_se) if return_se else False,
            context="MT-MET GP",
        )
        cov_traits_downgraded = bool(decision["downgraded"])
        cov_traits_downgrade_reason = str(decision["reason"])
        compute_cov_traits = bool(decision["compute"])
    if compute_cov_traits:
        unique_cells = np.unique(np.column_stack([g_te, e_te]), axis=0)
        cov_traits = np.zeros((unique_cells.shape[0], T, T), dtype=np.float64)
        cov_trait_cells = pd.DataFrame({
            "gid": [geno_ids[int(g)] for g in unique_cells[:, 0]],
            "env": [env_levels[int(e)] for e in unique_cells[:, 1]],
        })
        for col in ("gid", "env"):
            cov_trait_cells[col] = pd.Series([str(x) for x in cov_trait_cells[col].tolist()], index=cov_trait_cells.index, dtype=object)
        Sigma_G_t = torch.tensor(Sigma_G, dtype=torch.float64, device=device)
        Sigma_GE_t = torch.tensor(Sigma_GE, dtype=torch.float64, device=device)
        # Hoisted: train-side caches (constant across cells and traits)
        for s, (gi_, ei_) in enumerate(unique_cells.tolist()):
            gi = int(gi_); ei = int(ei_)
            K_geno_row = torch.zeros(fold_train_g.ng_used, dtype=torch.float64, device=device)
            K_gi_gi = 0.0
            for r, w in enumerate(kernel_weights_arr.tolist()):
                if w <= 0.0:
                    continue
                Phi_train_r = cache.get_rows(r, fold_train_g.g_unique, device=device, dtype=torch.float64)
                Phi_row_r = cache.get_rows(r, torch.tensor([gi], dtype=torch.int64),
                                           device=device, dtype=torch.float64)
                K_geno_row = K_geno_row + float(w) * (Phi_row_r @ Phi_train_r.transpose(0, 1)).squeeze(0)
                K_gi_gi += float(w) * float((Phi_row_r @ Phi_row_r.transpose(0, 1)).item())
            K_geno_col = K_geno_row[fold_train_g.inv]                         # (n_train_rows,)
            K_env_col = torch.tensor(K_env[ei, e_tr], dtype=torch.float64, device=device)
            K_env_ee = float(K_env[ei, ei])
            # Build all T columns and solve V α_t = col_t (T PCG solves per cell)
            cols = []
            alphas = []
            for t in range(T):
                col_t = (Sigma_G_t[t][t_tr_t] + Sigma_GE_t[t][t_tr_t] * K_env_col) * K_geno_col
                pcg_t = pcg_solve(A_mv, col_t, M_inv_mv=M_inv_mv, tol=pcg_tol, max_iter=pcg_max_iter)
                cols.append(col_t)
                alphas.append(pcg_t.x)
            # Bilinear Schur: schur[t, t'] = prior[t, t'] - col_t @ alpha_{t'}
            schur = np.zeros((T, T), dtype=np.float64)
            for t in range(T):
                prior_row = torch_to_numpy(
                    Sigma_G_t[t] * K_gi_gi + Sigma_GE_t[t] * K_gi_gi * K_env_ee
                )
                for tp in range(T):
                    schur[t, tp] = prior_row[tp] - float(torch.dot(cols[t], alphas[tp]).item())
            schur = 0.5 * (schur + schur.T)  # numerical symmetrization only
            sd = trait_env_sds[ei, :]   # (T,) sd at this env
            cov_traits[s] = (np.outer(sd, sd)) * schur

    if return_se and cov_traits_downgraded:
        # SE is derived from the per-cell trait covariance, so the downgrade that
        # skipped cov_traits necessarily skips SE too. Leave the SE columns out
        # rather than raising; the caller is told via cov_traits_downgraded.
        return_se = False
    if return_se:
        if cov_traits is None or cov_trait_cells is None:
            raise RuntimeError("internal error: return_se requested but cov_traits was not computed")
        cell_to_pos = {
            (int(g), int(e)): i for i, (g, e) in enumerate(unique_cells.tolist())
        }
        lat_var = np.zeros(g_te.shape[0], dtype=np.float64)
        obs_var = np.zeros(g_te.shape[0], dtype=np.float64)
        for i, (g, e, t) in enumerate(zip(g_te.tolist(), e_te.tolist(), t_te.tolist())):
            pos = cell_to_pos[(int(g), int(e))]
            lv = float(cov_traits[pos, int(t), int(t)])
            lat_var[i] = max(lv, 1e-12)
            obs_var[i] = lat_var[i] + float(Sigma_eps[int(t), int(t)]) * float(trait_env_sds[int(e), int(t)] ** 2)
        preds_df["SE_latent"] = np.sqrt(lat_var)
        preds_df["PEV"] = lat_var
        preds_df["SE"] = np.sqrt(np.maximum(obs_var, 1e-12))
        if not return_cov_traits:
            cov_traits = None
            cov_trait_cells = None

    return {
        "Sigma_G":            Sigma_G,
        "Sigma_GE":           Sigma_GE,
        "Sigma_eps":          Sigma_eps,
        "Lambda_G":           Lambda_G,
        "psi_G":              psi_G,
        "Lambda_GE":          Lambda_GE,
        "psi_GE":             psi_GE,
        "beta":               beta,
        "beta_mom":           beta_mom,
        "beta_method":        "ols_all_train",
        "mom_n_full_gids":    mom_out["mom_n_full_gids"],
        "mom_n_full_envs":    mom_out["mom_n_full_envs"],
        "mom_n_full_cells":   mom_out["mom_n_full_cells"],
        "mom_n_observed_cells": mom_out.get("mom_n_observed_cells"),
        "mom_psd_clip_count": mom_out["mom_psd_clip_count"],
        "mom_3x3_cond":       mom_out["mom_3x3_cond"],
        "mom_fallback":       mom_out.get("mom_fallback"),
        "pcg_n_iter":         int(pcg_res.iters),
        "pcg_residual_norm":  float(pcg_res.rel_resid),
        "pcg_converged":      bool(pcg_res.converged),
        "pcg_preconditioner": "diag",
        "kernel_weights":     kernel_weights_arr,
        "cov_traits_downgraded": bool(cov_traits_downgraded),
        "cov_traits_downgrade_reason": cov_traits_downgrade_reason,
        "device":             str(device),
        "reml_loglik_history": (None if reml_out is None else reml_out["loglik_history"]),
        "reml_n_iter":        (None if reml_out is None else reml_out["n_iter"]),
        "reml_converged":     (None if reml_out is None else reml_out["converged"]),
        "predictions":        preds_df,
        "cov_traits":         cov_traits,
        "cov_trait_cells":    cov_trait_cells,
    }
