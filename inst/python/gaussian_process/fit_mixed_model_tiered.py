"""Connectivity-aware per-tier mixed-model wrapper.

`fit_mixed_model_tiered` classifies each test row by the connectivity of its
environment to the training set and dispatches to the right kernel:

    cv0-tier  (test env has ZERO training rows)   ->  K_env model
    cv1-tier  (test env has training rows)        ->  Sigma_g_mom model

Both models run on the full training set; only the PREDICTION assembly is
tier-aware. This avoids the global-blending compromise that any single
additive kernel suffers when the test set mixes tiers.

Validated on G2F 2025 CV2-mixed split:
    pure K_env          +0.473 overall  RMSE 2.81
    pure Sigma_g_mom    +0.519 overall  RMSE 2.89
    global alpha=0.5    +0.562 overall  RMSE 2.72
    per-tier glue       +0.610 overall  RMSE 2.51   <- +4.8pp, -7.7% RMSE

Design notes
------------
* Under a single-'e' term, any additive/blended kernel K = alpha K_env +
  (1-alpha) Sigma_g_mom COUPLES cv0 and cv1 rows globally through a shared
  latent z_env. No choice of alpha decouples the regimes. The only way to
  decouple is to run two independent models and pick predictions per row.
* The MoM Sigma_g pass is shared work; we run it once and reuse it in the
  cv1-tier fit.
"""
from __future__ import annotations
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import time
import warnings

import numpy as np
import pandas as pd


def _normalize_cov_to_corr(S: np.ndarray, eig_floor: float = 1e-6) -> np.ndarray:
    """Project to nearest PSD correlation matrix (eig-clip + diag-renorm)."""
    S = 0.5 * (S + S.T)
    w, V = np.linalg.eigh(S)
    w = np.maximum(w, eig_floor)
    S = (V * w) @ V.T
    d = np.sqrt(np.maximum(np.diag(S), 1e-12))
    S = S / (d[:, None] * d[None, :])
    S = 0.5 * (S + S.T)
    np.fill_diagonal(S, 1.0)
    return S


def _repair_sigma_g_for_unseen_envs(Sigma: np.ndarray, eps: float = 1e-4) -> np.ndarray:
    """Replace zero-diag rows/cols (unidentified envs in MoM) with identity.

    cv0-tier envs have no training observations so their Sigma_g_mom rows are
    ~0. cov2corr would blow up; we stub them to identity so downstream
    correlation projection is well-defined. These envs will be served by the
    K_env model anyway — this is just to keep the Sigma-path matrix finite.
    """
    S = np.asarray(Sigma, dtype=float).copy()
    bad = np.diag(S) < eps
    if bad.any():
        S[bad, :] = 0.0
        S[:, bad] = 0.0
        S[bad, bad] = 1.0
    return S


def _preflight_validate(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    env_similarity: Optional[np.ndarray],
    env_covariates: Optional[Any],
    geno_ids: Sequence[str],
    grm_factor_cache: Optional[Dict[str, Any]],
) -> None:
    """Fail fast with actionable errors before launching expensive fits."""
    for col in (gid_col, env_col, y_col):
        if col not in pheno_df.columns:
            raise KeyError(
                f"fit_mixed_model_tiered: column {col!r} missing from pheno_df "
                f"(got columns: {list(pheno_df.columns)})"
            )

    n = len(pheno_df)
    tr = np.asarray(train_idx, dtype=np.int64)
    te = np.asarray(test_idx, dtype=np.int64)
    if tr.size == 0:
        raise ValueError("fit_mixed_model_tiered: train_idx is empty")
    if te.size == 0:
        raise ValueError("fit_mixed_model_tiered: test_idx is empty")
    if tr.min() < 0 or tr.max() >= n:
        raise ValueError(
            f"fit_mixed_model_tiered: train_idx out of range [0, {n}) "
            f"(min={tr.min()}, max={tr.max()})"
        )
    if te.min() < 0 or te.max() >= n:
        raise ValueError(
            f"fit_mixed_model_tiered: test_idx out of range [0, {n}) "
            f"(min={te.min()}, max={te.max()})"
        )
    overlap = np.intersect1d(tr, te)
    if overlap.size:
        raise ValueError(
            f"fit_mixed_model_tiered: train_idx and test_idx share "
            f"{overlap.size} rows (expected disjoint). "
            f"First few overlapping indices: {overlap[:5].tolist()}"
        )

    y_train = pd.to_numeric(pheno_df.loc[tr, y_col], errors="coerce").to_numpy()
    n_finite = int(np.isfinite(y_train).sum())
    if n_finite < 10:
        raise ValueError(
            f"fit_mixed_model_tiered: only {n_finite} finite y values in "
            f"train_idx (need >=10). Check that pheno_df[y_col] has values "
            f"at train rows."
        )

    if env_similarity is None and env_covariates is None:
        raise ValueError(
            "fit_mixed_model_tiered: requires env_similarity or env_covariates; "
            "per-tier routing is meaningless without K_env. Pass env-level "
            "features as env_covariates (DataFrame) or a pre-built n_env x "
            "n_env kernel as env_similarity."
        )

    if env_similarity is not None:
        K = np.asarray(env_similarity)
        if K.ndim != 2 or K.shape[0] != K.shape[1]:
            raise ValueError(
                f"fit_mixed_model_tiered: env_similarity must be 2D square; "
                f"got shape {K.shape}"
            )
        n_env_pheno = pheno_df[env_col].astype("category").nunique()
        if K.shape[0] != n_env_pheno:
            raise ValueError(
                f"fit_mixed_model_tiered: env_similarity shape {K.shape} does "
                f"not match n_env={n_env_pheno} inferred from pheno_df[env_col]. "
                f"Kernel rows must correspond to "
                f"pheno_df[env_col].astype('category').cat.categories order."
            )

    geno_ids_set = set(str(g) for g in geno_ids)
    if not geno_ids_set:
        raise ValueError("fit_mixed_model_tiered: geno_ids is empty")
    gids_in_pheno = set(pheno_df[gid_col].astype(str).unique())
    missing = gids_in_pheno - geno_ids_set
    if missing:
        sample = list(missing)[:5]
        raise ValueError(
            f"fit_mixed_model_tiered: {len(missing)} genotypes in pheno_df "
            f"are not in geno_ids (e.g. {sample}). The genotype cohort must "
            f"cover every gid in pheno_df."
        )

    if grm_factor_cache is not None:
        needed = {"type", "paths", "shapes", "dtypes"}
        missing_keys = needed - set(grm_factor_cache.keys())
        if missing_keys:
            raise ValueError(
                f"fit_mixed_model_tiered: grm_factor_cache missing keys "
                f"{missing_keys} (required: {needed})"
            )


def _env_connectivity_mask(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Classify test rows by env-level connectivity.

    A test row is cv0-tier if its env has ZERO training observations with a
    finite y. Otherwise cv1-tier. Returns a numpy array of length len(test_idx)
    with values {"cv0","cv1"} and a summary dict of per-tier counts.
    """
    train_envs = set(
        pheno_df.loc[train_idx][env_col].astype(str).unique().tolist()
    )
    test_envs = pheno_df.loc[test_idx][env_col].astype(str).to_numpy()
    is_cv0 = np.array([e not in train_envs for e in test_envs])
    tier = np.where(is_cv0, "cv0", "cv1")
    counts = {
        "cv0": int(is_cv0.sum()),
        "cv1": int((~is_cv0).sum()),
        "n_cv0_envs": int(len(set(test_envs[is_cv0]))),
        "n_cv1_envs": int(len(set(test_envs[~is_cv0]))),
    }
    return tier, counts


def _build_connectivity_summary(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    K_env: np.ndarray,
    env_levels: Sequence[str],
) -> pd.DataFrame:
    """Per test-env audit: tier, training support, and K_env proximity.

    For each env appearing in the test set, report the dispatch tier and the
    evidence backing the decision. For cv0 envs we also report the highest
    K_env similarity to any training env — a weak K_env connection means the
    reaction-norm prediction for that env rests on sparse evidence.
    """
    env_to_idx = {e: i for i, e in enumerate(env_levels)}
    tr = pheno_df.loc[train_idx]
    train_env_counts = (
        tr[env_col].astype(str).value_counts().to_dict()
    )
    train_env_hybs = (
        tr.groupby(tr[env_col].astype(str))[gid_col].nunique().to_dict()
    )
    train_envs_set = set(train_env_counts.keys())
    train_env_indices = np.array(
        [env_to_idx[e] for e in train_envs_set if e in env_to_idx],
        dtype=np.int64,
    )

    test_envs_unique = sorted(
        pheno_df.loc[test_idx][env_col].astype(str).unique().tolist()
    )
    rows = []
    test_row_env = pheno_df.loc[test_idx][env_col].astype(str).to_numpy()
    for e in test_envs_unique:
        in_train = e in train_envs_set
        tier = "cv1" if in_train else "cv0"
        row = {
            "env": e,
            "tier": tier,
            "n_test_rows": int((test_row_env == e).sum()),
            "n_train_rows": int(train_env_counts.get(e, 0)),
            "n_train_hybrids": int(train_env_hybs.get(e, 0)),
        }
        if (not in_train) and (e in env_to_idx) and train_env_indices.size:
            row_vec = K_env[env_to_idx[e], train_env_indices]
            row["max_Kenv_to_train"] = float(np.max(row_vec))
            row["median_Kenv_to_train"] = float(np.median(row_vec))
        else:
            row["max_Kenv_to_train"] = np.nan
            row["median_Kenv_to_train"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["tier", "env"]).reset_index(drop=True)


def _run_single_fit(
    fw: Any,
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernels: Optional[Dict[str, np.ndarray]],
    geno_ids: Sequence[str],
    env_similarity: Optional[np.ndarray],
    include_components: Sequence[str],
    w_g: Optional[Any],
    w_ge: Optional[Any],
    w_e: float,
    varcomp_mode: str,
    env_structure: Optional[str],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    backend: str,
    method: str,
    dtype: str,
    operator_dtype_compute: str,
    operator_tol: float,
    operator_max_iter: int,
    standardize: str,
    grm_factor_cache: Optional[Dict[str, Any]],
    seed: int,
    return_se: bool = False,
    n_hutchinson_probes: int = 32,
    hutchinson_seed: Optional[int] = None,
    se_pcg_tol: Optional[float] = None,
    se_pcg_max_iter: Optional[int] = None,
    use_nystrom_se_preconditioner: bool = False,
    nystrom_rank: int = 64,
    nystrom_seed: int = 0,
    use_lanczos_se_variance_reduction: bool = False,
    lanczos_iters: int = 80,
    lanczos_seed: int = 0,
    lanczos_theta_floor: float = 1e-8,
    env_resid_robust: str = "none",
    env_resid_shrink_tau: float = 0.0,
) -> Dict[str, Any]:
    extra_kwargs: Dict[str, Any] = {}
    if return_se:
        extra_kwargs["operator_return_se"] = True
        extra_kwargs["operator_n_hutchinson_probes"] = int(n_hutchinson_probes)
        if hutchinson_seed is not None:
            extra_kwargs["operator_hutchinson_seed"] = int(hutchinson_seed)
        if se_pcg_tol is not None:
            extra_kwargs["operator_se_pcg_tol"] = float(se_pcg_tol)
        if se_pcg_max_iter is not None:
            extra_kwargs["operator_se_pcg_max_iter"] = int(se_pcg_max_iter)
        if bool(use_nystrom_se_preconditioner):
            extra_kwargs["operator_use_nystrom_se_preconditioner"] = True
            extra_kwargs["operator_nystrom_rank"] = int(nystrom_rank)
            extra_kwargs["operator_nystrom_seed"] = int(nystrom_seed)
        if bool(use_lanczos_se_variance_reduction):
            extra_kwargs["operator_use_lanczos_se_variance_reduction"] = True
            extra_kwargs["operator_lanczos_iters"] = int(lanczos_iters)
            extra_kwargs["operator_lanczos_seed"] = int(lanczos_seed)
            extra_kwargs["operator_lanczos_theta_floor"] = float(lanczos_theta_floor)
    # Per-env residual variance — forwarded verbatim; framework already accepts both.
    extra_kwargs["env_resid_robust"] = str(env_resid_robust)
    extra_kwargs["env_resid_shrink_tau"] = float(env_resid_shrink_tau)
    output_level = "predict_with_se" if bool(return_se) else "predict_only"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = fw.fit_mixed_model(
            pheno_df=pheno_df,
            gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=list(geno_ids),
            env_similarity=env_similarity,
            include_components=tuple(include_components),
            w_g=w_g, w_ge=w_ge, w_e=w_e,
            fixed_effects=[], method=method,
            env_structure=env_structure, varcomp_mode=varcomp_mode,
            output_level=output_level, point_predictions_only=(not bool(return_se)),
            return_se=bool(return_se),
            prediction_output="test_only", backend=backend,
            grm_factor_cache=grm_factor_cache,
            train_idx=train_idx, test_idx=test_idx,
            standardize=standardize, dtype=dtype,
            operator_dtype_compute=operator_dtype_compute,
            operator_tol=operator_tol, operator_max_iter=operator_max_iter,
            reaction_norm_auto=False, seed=seed,
            **extra_kwargs,
        )
    return res


def _extract_preds(res: Dict[str, Any], gid_col: str, env_col: str) -> pd.DataFrame:
    p = res["result"]["predictions"].copy()
    id_col = "Genotype" if "Genotype" in p.columns else ("Name" if "Name" in p.columns else gid_col)
    p = p.rename(columns={id_col: gid_col, "Env": env_col, "Prediction": "Prediction"})
    keep = [gid_col, env_col, "Prediction"]
    extra = [c for c in p.columns if c not in keep]
    return p[keep + extra]


def fit_mixed_model_tiered(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernels: Optional[Dict[str, np.ndarray]],
    geno_ids: Sequence[str],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    # --- env kernel inputs: supply ONE of these ---
    env_similarity: Optional[np.ndarray] = None,
    env_covariates: Optional[Any] = None,
    feature_qc: bool = True,
    kenv_kernel: str = "matern32",
    kenv_bandwidth: float = 1.0,
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,
    # --- per-env residual variance knobs (forwarded to fit_mixed_model) ---
    env_resid_robust: str = "none",              # "none" | "mad"
    env_resid_shrink_tau: float = 0.0,           # 0=raw per-env, 1=pooled to mean
    # kenv_kernel='auto' embedded mini-CV knobs (forwarded to fit_mixed_model):
    kenv_auto_candidates: Optional[Any] = None,
    kenv_cv_groups: Optional[Any] = None,
    kenv_cv_kfolds: int = 5,
    kenv_cv_ridge: float = 1e-2,
    kenv_cv_seed: int = 0,
    # --- model weights (used for both tier fits unless *_cv0 / *_cv1 overrides given) ---
    w_g: Optional[Any] = None,
    w_ge: Optional[Any] = None,
    w_e: float = 0.5,
    include_components: Sequence[str] = ("g", "ge", "e"),
    # --- MoM Sigma_g aux pass ---
    mom_include_components: Sequence[str] = ("g", "ge"),
    mom_w_e: float = 0.0,
    mom_n_hutch: int = 16,
    mom_h2_prior: float = 0.5,
    # --- compute knobs ---
    method: str = "krr_exact",
    backend: str = "operator",
    dtype: str = "float32",
    operator_dtype_compute: str = "float32",
    operator_tol: float = 1e-4,
    operator_max_iter: int = 300,
    standardize: str = "global",
    grm_factor_cache: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    # --- standard errors (Hutchinson probing) ---
    return_se: bool = False,
    n_hutchinson_probes: int = 32,
    hutchinson_seed: Optional[int] = None,
    se_pcg_tol: Optional[float] = None,
    se_pcg_max_iter: Optional[int] = None,
    use_nystrom_se_preconditioner: bool = False,
    nystrom_rank: int = 64,
    nystrom_seed: int = 0,
    use_lanczos_se_variance_reduction: bool = False,
    lanczos_iters: int = 80,
    lanczos_seed: int = 0,
    lanczos_theta_floor: float = 1e-8,
    # --- diagnostics ---
    verbose: bool = True,
) -> Dict[str, Any]:
    """Connectivity-aware per-tier mixed model.

    Runs two parallel fits with identical architecture but different env
    kernels (K_env for cv0-tier rows, Sigma_g_mom-derived kernel for
    cv1-tier rows), classifies each test row by env connectivity, and glues
    predictions per-tier. Returns a predictions DataFrame with a ``tier``
    column plus diagnostic wall-times and the underlying fit results.

    Parameters
    ----------
    pheno_df, gid_col, env_col, y_col
        Phenotype frame. Training y values present; test y may be NaN.
    geno_kernels, geno_ids, grm_factor_cache
        Forwarded to ``fit_mixed_model``. For n_geno >= 3000 supply
        ``grm_factor_cache`` (see README.md in this folder).
    train_idx, test_idx
        Row-index arrays into ``pheno_df`` defining the fold.
    env_similarity
        Pre-built K_env (n_env x n_env PSD). Aligned to
        ``pheno_df[env_col].astype('category').cat.categories``. If None,
        built from ``env_covariates``.
    env_covariates
        DataFrame/ndarray/dict of env-level features. Passed to
        ``fw._env_covariates_to_aligned_matrix`` then ``_build_kenv_from_Z``
        with kernel/bandwidth determined by ``kenv_kernel`` / ``kenv_bandwidth``.
    feature_qc
        When True and EC features match known G2F patterns, restrict to the
        27 yield-relevant features.
    kenv_kernel, kenv_bandwidth
        Reaction-norm K_env kernel choice. Defaults to ``"matern32"`` with
        ``bandwidth=1.0`` (× median pairwise distance), chosen via 10-fold
        leave-one-year-out validation on G2F 2025 (matern beats linear on
        8/10 folds). Pass ``"linear"`` for the prior classical default.

    Returns
    -------
    dict with keys:
        predictions : DataFrame [gid_col, env_col, Prediction, tier]
        tier_counts : dict (cv0, cv1, n_cv0_envs, n_cv1_envs)
        connectivity : DataFrame, one row per test env with tier, training
            support (n_train_rows, n_train_hybrids), and (for cv0 envs) the
            max / median K_env similarity to any training env. Low max-sim
            flags cv0 envs where the reaction-norm prediction is weak.
        K_env, Sigma_g_mom : n_env x n_env arrays (or None if not computed)
        fit_K_env, fit_Sigma_g_mom, fit_mom_aux : raw fit_mixed_model results
        wall : dict of per-stage seconds
    """
    import importlib
    fw = importlib.import_module(
        "gp_framework"
    )

    _preflight_validate(
        pheno_df=pheno_df, gid_col=gid_col, env_col=env_col, y_col=y_col,
        train_idx=train_idx, test_idx=test_idx,
        env_similarity=env_similarity, env_covariates=env_covariates,
        geno_ids=geno_ids, grm_factor_cache=grm_factor_cache,
    )

    wall: Dict[str, float] = {}
    if verbose:
        print("=" * 72)
        print("fit_mixed_model_tiered: connectivity-aware per-tier dispatch")
        print("=" * 72)

    # ---- [1] tier classification ----
    t0 = time.perf_counter()
    tier, counts = _env_connectivity_mask(
        pheno_df, gid_col, env_col, y_col, train_idx, test_idx
    )
    wall["tier_classify"] = time.perf_counter() - t0
    if verbose:
        print(f"[1] Tier classification: cv0={counts['cv0']} rows "
              f"({counts['n_cv0_envs']} envs), cv1={counts['cv1']} rows "
              f"({counts['n_cv1_envs']} envs)")

    has_cv0 = counts["cv0"] > 0
    has_cv1 = counts["cv1"] > 0

    # ---- [2] build K_env (correlation-scaled) ----
    env_levels = list(
        pheno_df[env_col].astype("category").cat.categories.astype(str)
    )
    _kk_resolved = str(kenv_kernel).lower()
    _kbw_resolved = float(kenv_bandwidth)
    _kkw_resolved: Dict[str, Any] = dict(kenv_kernel_kwargs) if kenv_kernel_kwargs else {}
    kenv_auto_info: Optional[Dict[str, Any]] = None
    if env_similarity is not None:
        K_env = np.asarray(env_similarity, dtype=float)
        if K_env.shape != (len(env_levels), len(env_levels)):
            raise ValueError(
                f"env_similarity shape {K_env.shape} != "
                f"n_env={len(env_levels)} (from pheno_df env categories)"
            )
        K_env = _normalize_cov_to_corr(K_env)
        wall["build_K_env"] = 0.0
    else:
        t0 = time.perf_counter()
        Z, _ = fw._env_covariates_to_aligned_matrix(
            env_covariates, env_levels, env_col=env_col, feature_qc=feature_qc
        )
        if _kk_resolved == "auto":
            y_arr = pheno_df[y_col].to_numpy(dtype=np.float64)
            ei_arr = pheno_df[env_col].astype("category").cat.set_categories(
                env_levels).cat.codes.to_numpy(dtype=np.int64)
            n_env_total = len(env_levels)
            sums = np.zeros(n_env_total, dtype=np.float64)
            counts = np.zeros(n_env_total, dtype=np.int64)
            tr = np.asarray(train_idx, dtype=np.int64)
            tr_finite = tr[np.isfinite(y_arr[tr])]
            np.add.at(sums, ei_arr[tr_finite], y_arr[tr_finite])
            np.add.at(counts, ei_arr[tr_finite], 1)
            env_means = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)
            train_env_mask = counts > 0
            cands = (kenv_auto_candidates if kenv_auto_candidates is not None
                     else (("linear", 1.0), ("matern32", 1.0)))
            cands = fw._merge_kwargs_into_candidates(
                tuple(cands), kenv_kernel_kwargs or {},
            )
            _kk_resolved, _kbw_resolved, kenv_auto_info = fw._autoselect_kenv_kernel(
                Z, env_means, train_env_mask,
                cv_groups=(np.asarray(kenv_cv_groups) if kenv_cv_groups is not None else None),
                candidates=tuple(cands),
                k_folds=int(kenv_cv_kfolds),
                ridge=float(kenv_cv_ridge),
                seed=int(kenv_cv_seed),
            )
            _kkw_resolved = dict(kenv_auto_info.get("kwargs", {}))
            if verbose:
                print(f"[2a] kenv_kernel='auto' chose {_kk_resolved} bw={_kbw_resolved} "
                      f"(cv_design={kenv_auto_info.get('cv_design','n/a')}, "
                      f"summary={kenv_auto_info.get('summary')})")
        K_env = _normalize_cov_to_corr(fw._build_kenv_from_Z(
            Z, kernel=_kk_resolved, bandwidth=_kbw_resolved,
            kernel_kwargs=_kkw_resolved,
        ))
        wall["build_K_env"] = time.perf_counter() - t0
        if verbose:
            print(f"[2] K_env built from env_covariates in "
                  f"{wall['build_K_env']:.1f}s (n_env={len(env_levels)}, "
                  f"kernel={_kk_resolved}, bandwidth={_kbw_resolved})")

    # ---- [3] auxiliary MoM pass to learn Sigma_g (shared by cv1-tier fit) ----
    Sigma_g_mom = None
    K_mom = None
    res_mom_aux = None
    if has_cv1:
        t0 = time.perf_counter()
        if verbose:
            print("[3] Learning Sigma_g via MoM (aux pass)...")
        res_mom_aux = _run_single_fit(
            fw=fw, pheno_df=pheno_df,
            gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=geno_ids,
            env_similarity=None,
            include_components=mom_include_components,
            w_g=None, w_ge=None, w_e=mom_w_e,
            varcomp_mode="mom", env_structure="fa",
            train_idx=train_idx, test_idx=test_idx,
            backend=backend, method=method, dtype=dtype,
            operator_dtype_compute=operator_dtype_compute,
            operator_tol=operator_tol, operator_max_iter=operator_max_iter,
            standardize=standardize, grm_factor_cache=grm_factor_cache,
            seed=seed,
            env_resid_robust=env_resid_robust,
            env_resid_shrink_tau=env_resid_shrink_tau,
        )
        wall["mom_aux"] = time.perf_counter() - t0
        if "env_covariance_fa" not in res_mom_aux:
            raise RuntimeError(
                "MoM aux pass did not return env_covariance_fa; "
                "check env_structure='fa' support in framework."
            )
        Sigma_g_mom = np.asarray(res_mom_aux["env_covariance_fa"], dtype=float)
        Sigma_g_mom = _repair_sigma_g_for_unseen_envs(Sigma_g_mom)
        K_mom = _normalize_cov_to_corr(Sigma_g_mom)
        if verbose:
            print(f"    MoM wall={wall['mom_aux']:.1f}s  "
                  f"Sigma_g diag: min={np.diag(Sigma_g_mom).min():.3f} "
                  f"med={np.median(np.diag(Sigma_g_mom)):.3f} "
                  f"max={np.diag(Sigma_g_mom).max():.3f}")
    else:
        wall["mom_aux"] = 0.0
        if verbose:
            print("[3] Skipping MoM aux pass (no cv1-tier test rows)")

    # ---- [4] fit_A: K_env-driven model (serves cv0-tier rows) ----
    res_K_env = None
    if has_cv0:
        t0 = time.perf_counter()
        if verbose:
            print("[4] Fitting K_env model for cv0-tier rows...")
        res_K_env = _run_single_fit(
            fw=fw, pheno_df=pheno_df,
            gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=geno_ids,
            env_similarity=K_env,
            include_components=include_components,
            w_g=w_g, w_ge=w_ge, w_e=w_e,
            varcomp_mode="reml", env_structure=None,
            train_idx=train_idx, test_idx=test_idx,
            backend=backend, method=method, dtype=dtype,
            operator_dtype_compute=operator_dtype_compute,
            operator_tol=operator_tol, operator_max_iter=operator_max_iter,
            standardize=standardize, grm_factor_cache=grm_factor_cache,
            seed=seed,
            return_se=return_se,
            n_hutchinson_probes=n_hutchinson_probes,
            hutchinson_seed=hutchinson_seed,
            se_pcg_tol=se_pcg_tol,
            se_pcg_max_iter=se_pcg_max_iter,
            use_nystrom_se_preconditioner=use_nystrom_se_preconditioner,
            nystrom_rank=nystrom_rank,
            nystrom_seed=nystrom_seed,
            use_lanczos_se_variance_reduction=use_lanczos_se_variance_reduction,
            lanczos_iters=lanczos_iters,
            lanczos_seed=lanczos_seed,
            lanczos_theta_floor=lanczos_theta_floor,
            env_resid_robust=env_resid_robust,
            env_resid_shrink_tau=env_resid_shrink_tau,
        )
        wall["fit_K_env"] = time.perf_counter() - t0
        if verbose:
            print(f"    K_env fit wall={wall['fit_K_env']:.1f}s")
    else:
        wall["fit_K_env"] = 0.0

    # ---- [5] fit_B: Sigma_g_mom-driven model (serves cv1-tier rows) ----
    res_Sigma_g = None
    if has_cv1:
        t0 = time.perf_counter()
        if verbose:
            print("[5] Fitting Sigma_g_mom model for cv1-tier rows...")
        res_Sigma_g = _run_single_fit(
            fw=fw, pheno_df=pheno_df,
            gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=geno_ids,
            env_similarity=K_mom,
            include_components=include_components,
            w_g=w_g, w_ge=w_ge, w_e=w_e,
            varcomp_mode="reml", env_structure=None,
            train_idx=train_idx, test_idx=test_idx,
            backend=backend, method=method, dtype=dtype,
            operator_dtype_compute=operator_dtype_compute,
            operator_tol=operator_tol, operator_max_iter=operator_max_iter,
            standardize=standardize, grm_factor_cache=grm_factor_cache,
            seed=seed,
            return_se=return_se,
            n_hutchinson_probes=n_hutchinson_probes,
            hutchinson_seed=(None if hutchinson_seed is None else int(hutchinson_seed) + 1),
            se_pcg_tol=se_pcg_tol,
            se_pcg_max_iter=se_pcg_max_iter,
            use_nystrom_se_preconditioner=use_nystrom_se_preconditioner,
            nystrom_rank=nystrom_rank,
            nystrom_seed=(int(nystrom_seed) + 1),
            use_lanczos_se_variance_reduction=use_lanczos_se_variance_reduction,
            lanczos_iters=lanczos_iters,
            lanczos_seed=(int(lanczos_seed) + 1),
            lanczos_theta_floor=lanczos_theta_floor,
            env_resid_robust=env_resid_robust,
            env_resid_shrink_tau=env_resid_shrink_tau,
        )
        wall["fit_Sigma_g_mom"] = time.perf_counter() - t0
        if verbose:
            print(f"    Sigma_g_mom fit wall={wall['fit_Sigma_g_mom']:.1f}s")
    else:
        wall["fit_Sigma_g_mom"] = 0.0

    # ---- [6] per-tier glue ----
    t0 = time.perf_counter()
    test_frame = pheno_df.loc[test_idx, [gid_col, env_col]].reset_index(drop=True).copy()
    test_frame["tier"] = tier
    test_frame["_row"] = np.arange(len(test_frame))

    pred_vec = np.full(len(test_frame), np.nan, dtype=float)
    se_vec = np.full(len(test_frame), np.nan, dtype=float) if return_se else None
    se_lat_vec = np.full(len(test_frame), np.nan, dtype=float) if return_se else None
    var_obs_vec = np.full(len(test_frame), np.nan, dtype=float) if return_se else None
    var_lat_vec = np.full(len(test_frame), np.nan, dtype=float) if return_se else None

    def _merge_tier(res_tier, tier_name: str, is_tier_mask: np.ndarray) -> None:
        p = _extract_preds(res_tier, gid_col, env_col)
        p = p.drop_duplicates(subset=[gid_col, env_col])
        wanted = [gid_col, env_col, "Prediction"]
        if return_se:
            def _first_present(candidates: Sequence[str]) -> Optional[str]:
                for col in candidates:
                    if col in p.columns:
                        return col
                return None

            se_obs_col = _first_present(("SE", "SE_observed", "Prediction_SE_observed", "Prediction_SE"))
            se_lat_col = _first_present(("SE_latent", "Prediction_SE_latent"))
            var_obs_col = _first_present(("Prediction_Var_observed", "Prediction_variance"))
            var_lat_col = _first_present(("PEV", "Prediction_Var_latent"))
            if se_obs_col is None and var_obs_col is not None:
                p["SE"] = np.sqrt(np.maximum(pd.to_numeric(p[var_obs_col], errors="coerce"), 0.0))
                se_obs_col = "SE"
            if se_lat_col is None and var_lat_col is not None:
                p["SE_latent"] = np.sqrt(np.maximum(pd.to_numeric(p[var_lat_col], errors="coerce"), 0.0))
                se_lat_col = "SE_latent"
            if var_obs_col is None and se_obs_col is not None:
                p["Prediction_Var_observed"] = pd.to_numeric(p[se_obs_col], errors="coerce") ** 2
                var_obs_col = "Prediction_Var_observed"
            if var_lat_col is None and se_lat_col is not None:
                p["PEV"] = pd.to_numeric(p[se_lat_col], errors="coerce") ** 2
                var_lat_col = "PEV"
            column_map = {
                se_obs_col: "SE",
                se_lat_col: "SE_latent",
                var_obs_col: "Prediction_Var_observed",
                var_lat_col: "PEV",
            }
            for source, target in column_map.items():
                if source is None:
                    continue
                if source != target:
                    p[target] = p[source]
                if target not in wanted:
                    wanted.append(target)
        m = test_frame.merge(p[wanted], on=[gid_col, env_col], how="left",
                             suffixes=("", f"_{tier_name}"))
        pred_vec[is_tier_mask] = m["Prediction"].to_numpy()[is_tier_mask]
        if return_se and ("SE" in m.columns):
            se_vec[is_tier_mask] = m["SE"].to_numpy()[is_tier_mask]
        if return_se and ("SE_latent" in m.columns):
            se_lat_vec[is_tier_mask] = m["SE_latent"].to_numpy()[is_tier_mask]
        if return_se and ("Prediction_Var_observed" in m.columns):
            var_obs_vec[is_tier_mask] = m["Prediction_Var_observed"].to_numpy()[is_tier_mask]
        if return_se and ("PEV" in m.columns):
            var_lat_vec[is_tier_mask] = m["PEV"].to_numpy()[is_tier_mask]

    if res_K_env is not None:
        is_cv0 = (test_frame["tier"] == "cv0").to_numpy()
        _merge_tier(res_K_env, "kenv", is_cv0)

    if res_Sigma_g is not None:
        is_cv1 = (test_frame["tier"] == "cv1").to_numpy()
        _merge_tier(res_Sigma_g, "mom", is_cv1)

    n_unresolved = int(np.isnan(pred_vec).sum())
    if n_unresolved > 0:
        warnings.warn(
            f"fit_mixed_model_tiered: {n_unresolved} test rows received no "
            f"prediction (likely unknown gid/env in one of the sub-fits)."
        )

    predictions = test_frame.drop(columns=["_row"]).copy()
    predictions["Prediction"] = pred_vec
    if return_se:
        if np.isnan(var_obs_vec).all() and not np.isnan(se_vec).all():
            var_obs_vec = np.square(se_vec)
        if np.isnan(var_lat_vec).all() and not np.isnan(se_lat_vec).all():
            var_lat_vec = np.square(se_lat_vec)
        if np.isnan(se_vec).all() and not np.isnan(var_obs_vec).all():
            se_vec = np.sqrt(np.maximum(var_obs_vec, 0.0))
        if np.isnan(se_lat_vec).all() and not np.isnan(var_lat_vec).all():
            se_lat_vec = np.sqrt(np.maximum(var_lat_vec, 0.0))
        predictions["SE"] = se_vec
        predictions["SE_latent"] = se_lat_vec
        predictions["PEV"] = var_lat_vec
        predictions["Prediction_SE_observed"] = se_vec
        predictions["Prediction_Var_observed"] = var_obs_vec
        predictions["Prediction_SE_latent"] = se_lat_vec
        predictions["Prediction_Var_latent"] = var_lat_vec
    wall["glue"] = time.perf_counter() - t0

    wall["total"] = sum(v for v in wall.values())
    if verbose:
        print(f"[6] Glued predictions: n={len(predictions)} "
              f"(cv0={counts['cv0']} from K_env, cv1={counts['cv1']} from Sigma_g_mom)")
        print(f"    total wall={wall['total']:.1f}s "
              f"(mom_aux={wall['mom_aux']:.1f} + K_env={wall['fit_K_env']:.1f} "
              f"+ Sigma_g={wall['fit_Sigma_g_mom']:.1f})")
        print("=" * 72)

    connectivity = _build_connectivity_summary(
        pheno_df, gid_col, env_col, train_idx, test_idx, K_env, env_levels
    )
    if verbose:
        cv0_rows = connectivity[connectivity["tier"] == "cv0"]
        if len(cv0_rows):
            med_max = float(cv0_rows["max_Kenv_to_train"].median())
            min_max = float(cv0_rows["max_Kenv_to_train"].min())
            print(f"[connectivity] cv0 envs: n={len(cv0_rows)}  "
                  f"max K_env-to-train median={med_max:.3f}  min={min_max:.3f}  "
                  f"(lower => weaker reaction-norm evidence)")

    return {
        "predictions": predictions,
        "tier_counts": counts,
        "connectivity": connectivity,
        "K_env": K_env,
        "Sigma_g_mom": Sigma_g_mom,
        "fit_K_env": res_K_env,
        "fit_Sigma_g_mom": res_Sigma_g,
        "fit_mom_aux": res_mom_aux,
        "wall": wall,
        "kenv_kernel": _kk_resolved,
        "kenv_bandwidth": _kbw_resolved,
        "kenv_kernel_kwargs": _kkw_resolved,
        "kenv_auto": kenv_auto_info,
    }


__all__ = ["fit_mixed_model_tiered"]
