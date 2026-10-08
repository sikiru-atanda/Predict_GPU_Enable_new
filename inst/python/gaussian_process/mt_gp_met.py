"""Multi-environment MT-GP (v1).

Public entry: fit_multi_trait_gp_met(...).
Joint multi-trait GP across environments.
"""
from __future__ import annotations
import time
import warnings
from typing import Any, Dict, Optional, Sequence
import numpy as np
import pandas as pd

from mt_gp_met_internals import _build_mt_met_design, _fit_mt_met_mom_op
from mt_gp_reporting import covariance_by_environment_payload
from gp_device import pick_torch_device


def _validate_env_similarity_matrix(K_env: np.ndarray, n_env: int) -> np.ndarray:
    K_env = np.asarray(K_env, dtype=np.float64)
    if K_env.shape != (n_env, n_env):
        raise ValueError(
            f"env_similarity shape {K_env.shape} != n_env={n_env} (from pheno_df)."
        )
    if not np.isfinite(K_env).all():
        raise ValueError("env_similarity must contain only finite values")
    asym = float(np.max(np.abs(K_env - K_env.T))) if K_env.size else 0.0
    scale = max(1.0, float(np.max(np.abs(K_env))) if K_env.size else 1.0)
    if asym > 1e-8 * scale:
        raise ValueError(f"env_similarity must be symmetric; max asymmetry={asym:.3e}")
    K_env = 0.5 * (K_env + K_env.T)
    eig_min = float(np.linalg.eigvalsh(K_env).min()) if n_env else 0.0
    if eig_min < -1e-8 * scale:
        raise ValueError(
            f"env_similarity must be positive semidefinite; min eigenvalue={eig_min:.3e}"
        )
    return K_env


def _align_env_similarity(
    env_similarity: Any,
    env_levels: Sequence[str],
    env_ids: Optional[Sequence[str]] = None,
) -> np.ndarray:
    """Align a user-supplied environment kernel to the design's sorted env levels."""
    env_levels_str = [str(e) for e in env_levels]
    n_env = len(env_levels_str)
    if isinstance(env_similarity, pd.DataFrame):
        if env_ids is not None:
            warnings.warn(
                "env_ids is ignored when env_similarity is a DataFrame with labels.",
                stacklevel=2,
            )
        K_df = env_similarity.copy()
        K_df.index = K_df.index.map(str)
        K_df.columns = K_df.columns.map(str)
        if not K_df.index.is_unique or not K_df.columns.is_unique:
            raise ValueError("env_similarity DataFrame index/columns must have unique labels")
        missing = [e for e in env_levels_str if e not in K_df.index or e not in K_df.columns]
        if missing:
            raise ValueError(
                "env_similarity DataFrame is missing environments from pheno_df: "
                f"{missing[:5]}"
            )
        return _validate_env_similarity_matrix(
            K_df.loc[env_levels_str, env_levels_str].to_numpy(dtype=np.float64),
            n_env,
        )

    K_env = np.asarray(env_similarity, dtype=np.float64)
    if env_ids is not None:
        if K_env.ndim != 2:
            raise ValueError(f"env_similarity must be a 2D matrix; got ndim={K_env.ndim}")
        ids = [str(e) for e in env_ids]
        if len(ids) != K_env.shape[0] or K_env.shape[0] != K_env.shape[1]:
            raise ValueError(
                "env_ids length must match both dimensions of env_similarity; "
                f"got len(env_ids)={len(ids)}, env_similarity shape={K_env.shape}"
            )
        if len(set(ids)) != len(ids):
            raise ValueError("env_ids contains duplicate environment labels")
        pos = {e: i for i, e in enumerate(ids)}
        missing = [e for e in env_levels_str if e not in pos]
        if missing:
            raise ValueError(f"env_ids is missing environments from pheno_df: {missing[:5]}")
        order = [pos[e] for e in env_levels_str]
        K_env = K_env[np.ix_(order, order)]
    return _validate_env_similarity_matrix(K_env, n_env)


def fit_multi_trait_gp_met(
    *,
    pheno_df: pd.DataFrame,
    gid_col: str, env_col: str, trait_col: str, y_col: str,
    geno_kernels: Dict[str, np.ndarray],
    geno_ids: Sequence[str],
    train_idx: np.ndarray, test_idx: np.ndarray,
    env_similarity: Optional[np.ndarray] = None,
    env_ids: Optional[Sequence[str]] = None,
    env_covariates: Optional[Any] = None,
    reaction_norm_feature_qc: bool = True,
    kenv_kernel: str = "matern32",
    kenv_bandwidth: float = 1.0,
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,
    trait_structure: str = "unstructured",
    trait_fa_rank: Optional[int] = None,
    gxe_trait_structure: Optional[str] = None,
    gxe_trait_fa_rank: Optional[int] = None,
    fixed_effects: Sequence[str] = (),
    kernel_weights: Optional[Sequence[float]] = None,
    varcomp_mode: str = "mom",
    grm_factor_cache: Any = None,
    mom_pcg_tol: float = 1e-4,
    mom_pcg_max_iter: int = 300,
    mom_psd_floor: float = 1e-4,
    dense_reml_max_train: int = 3000,
    reml_max_iter: int = 30,
    reml_tol: float = 1e-5,
    return_se: bool = False,
    prediction_output: str = "test_only",
    force_prediction_se: bool = False,
    return_cov_traits: bool = False,
    force_cov_traits: bool = False,
    return_trait_correlations: bool = False,
    seed: int = 0,
    dtype: str = "float64",
    device: Optional[str] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Multi-environment multi-trait GP.

    `varcomp_mode="mom"` is the scalable default. `varcomp_mode="reml"` runs a
    bounded dense REML variance-estimation pass, then still uses the GP/operator
    path for prediction.

    Important: pheno_df must have one row per (gid, env, trait) cell. If your
    data has multiple replicates per cell (e.g., G2F plot-level data with reps),
    pre-aggregate to cell means via:

        df = df_plot.groupby([gid_col, env_col, trait_col], as_index=False)[y_col].mean()

    The residual operator assumes Cov(eps_r, eps_s) = Σ_eps · δ((g_r,e_r)==(g_s,e_s));
    multiple replicates per cell over-couple the residual structure and prevent
    PCG convergence.

    If env_similarity is a pandas DataFrame, its index/columns are aligned to
    the environment labels in pheno_df. If it is a NumPy array with a custom
    order, pass env_ids to avoid silent environment-kernel misalignment.

    See README.md in this folder.
    """
    t0 = time.perf_counter()
    torch_device = pick_torch_device(device)
    prediction_output = str(prediction_output).strip().lower()
    if prediction_output not in {"test_only", "all"}:
        raise ValueError(
            "prediction_output must be one of {'test_only', 'all'}; "
            f"got {prediction_output!r}"
        )
    mode = str(varcomp_mode).lower()
    if mode not in ("mom", "reml", "reml_dense"):
        raise ValueError(
            "multi-env supports varcomp_mode='mom' for scalable fitting or "
            "varcomp_mode='reml'/'reml_dense' for bounded dense variance estimation."
        )
    dense_reml = mode in ("reml", "reml_dense")
    if grm_factor_cache is None:
        raise ValueError(
            "varcomp_mode='mom' requires grm_factor_cache (memmap/zarr Phi); "
            "see the large-n notes in inst/python/gaussian_process/README.md."
        )
    if (env_similarity is None) == (env_covariates is None):
        # Both or neither
        if env_similarity is None and env_covariates is None:
            raise ValueError(
                "MT-GP MET needs an env kernel; supply env_similarity or env_covariates."
            )
        else:
            raise ValueError("Pass exactly one of env_similarity or env_covariates")

    if (geno_kernels and len(geno_kernels) > 1 and
            len(geno_kernels) != int(grm_factor_cache.num_grms)):
        raise ValueError(
            "geno_kernels and grm_factor_cache must describe the same number "
            f"of kernels; got {len(geno_kernels)} dense kernels and "
            f"{int(grm_factor_cache.num_grms)} factor roots."
        )
    trait_structure = str(trait_structure or "unstructured").lower()
    gxe_trait_structure = str(gxe_trait_structure or trait_structure).lower()
    if trait_structure not in ("unstructured", "fa"):
        raise ValueError("trait_structure must be 'unstructured' or 'fa'")
    if gxe_trait_structure not in ("unstructured", "fa"):
        raise ValueError("gxe_trait_structure must be 'unstructured' or 'fa'")
    if trait_structure == "fa" and trait_fa_rank is None:
        raise ValueError("trait_fa_rank required when trait_structure='fa'")
    if gxe_trait_structure == "fa" and gxe_trait_fa_rank is None and trait_fa_rank is None:
        raise ValueError("gxe_trait_fa_rank required when gxe_trait_structure='fa'")
    if trait_structure == "unstructured" and trait_fa_rank is not None:
        warnings.warn("trait_fa_rank is ignored when trait_structure='unstructured'.", stacklevel=2)
    if gxe_trait_structure == "unstructured" and gxe_trait_fa_rank is not None:
        warnings.warn("gxe_trait_fa_rank is ignored when gxe_trait_structure='unstructured'.", stacklevel=2)
    if int(seed) != 0:
        warnings.warn(f"seed={seed} is ignored in v1 (deterministic MoM path).", stacklevel=2)
    if str(dtype) != "float64":
        warnings.warn(f"dtype={dtype!r} is ignored in v1 (always float64).", stacklevel=2)
    if verbose:
        warnings.warn("verbose=True is ignored in v1.", stacklevel=2)
    if env_similarity is not None and kenv_kernel_kwargs:
        warnings.warn(
            "kenv_kernel_kwargs is ignored when env_similarity is supplied directly.",
            stacklevel=2,
        )
    if env_covariates is not None and env_ids is not None:
        warnings.warn("env_ids is ignored when env_covariates is used.", stacklevel=2)

    # Build design first (need env_levels to size K_env if it's user-supplied)
    design = _build_mt_met_design(
        pheno_df=pheno_df, gid_col=gid_col, env_col=env_col,
        trait_col=trait_col, y_col=y_col,
        fixed_effects=tuple(fixed_effects),
        geno_ids=list(geno_ids),
        train_idx=train_idx,
        test_idx=test_idx,
    )
    if prediction_output == "all":
        design["test_mask"] = np.ones_like(design["test_mask"], dtype=bool)
    n_env = design["n_env"]
    env_levels = design["env_levels"]

    if env_similarity is not None:
        K_env = _align_env_similarity(env_similarity, env_levels, env_ids=env_ids)
        kenv_kernel_used = "user_supplied"
        kenv_bandwidth_used = None
        kenv_kernel_kwargs_used = None
        env_covariates_qc_info = None
    else:
        # Auto-build via single-trait reaction-norm helpers
        # Lazy-import the framework to avoid pulling it in for env_similarity path
        import importlib
        fw = importlib.import_module(
            "gp_framework"
        )
        Z, env_covariates_qc_info = fw._env_covariates_to_aligned_matrix(
            env_covariates, env_levels, env_col=env_col,
            feature_qc=bool(reaction_norm_feature_qc),
        )
        if Z is None:
            raise ValueError("env_covariates produced empty alignment; check column names")
        kenv_kernel_kwargs_used = dict(kenv_kernel_kwargs or {})
        K_env = fw._build_kenv_from_Z(
            Z, kernel=str(kenv_kernel), bandwidth=float(kenv_bandwidth),
            kernel_kwargs=kenv_kernel_kwargs_used,
        )
        K_env = _validate_env_similarity_matrix(K_env, n_env)
        kenv_kernel_used = str(kenv_kernel)
        kenv_bandwidth_used = float(kenv_bandwidth)

    mom_out = _fit_mt_met_mom_op(
        design=design, cache=grm_factor_cache, K_env=K_env,
        psd_floor=float(mom_psd_floor),
        pcg_tol=float(mom_pcg_tol),
        pcg_max_iter=int(mom_pcg_max_iter),
        return_cov_traits=bool(return_cov_traits),
        force_cov_traits=bool(force_cov_traits),
        return_se=bool(return_se),
        force_prediction_se=bool(force_prediction_se),
        trait_structure=trait_structure,
        trait_fa_rank=trait_fa_rank,
        gxe_trait_structure=gxe_trait_structure,
        gxe_trait_fa_rank=gxe_trait_fa_rank,
        kernel_weights=kernel_weights,
        dense_reml=bool(dense_reml),
        dense_reml_max_train=int(dense_reml_max_train),
        reml_max_iter=int(reml_max_iter),
        reml_tol=float(reml_tol),
        device=torch_device,
    )
    preds_df = mom_out["predictions"].rename(
        columns={"gid": gid_col, "env": env_col, "trait": trait_col})
    cov_trait_cells = mom_out.get("cov_trait_cells")
    if cov_trait_cells is not None:
        cov_trait_cells = cov_trait_cells.rename(columns={"gid": gid_col, "env": env_col})
    wall = time.perf_counter() - t0
    result = {
        "predictions": preds_df,
        "cov_traits":  mom_out["cov_traits"],
        "cov_trait_cells": cov_trait_cells,
    }
    # Response-scale covariance matrices are variance-component evidence, not
    # optional correlation diagnostics. Always return them so R never falls
    # back to standardized-scale Sigma matrices when trait correlations were
    # not requested.
    genetic_payload = covariance_by_environment_payload(
        mom_out["Sigma_G"], design["trait_env_sds"],
        design["trait_levels"], design["env_levels"]
    )
    gxe_payload = covariance_by_environment_payload(
        mom_out["Sigma_GE"], design["trait_env_sds"],
        design["trait_levels"], design["env_levels"]
    )
    eps_payload = covariance_by_environment_payload(
        mom_out["Sigma_eps"], design["trait_env_sds"],
        design["trait_levels"], design["env_levels"]
    )
    result["genetic_covariance"] = genetic_payload["mean_covariance"]
    result["genetic_covariance_by_environment"] = genetic_payload["covariance_by_environment"]
    result["gxe_covariance"] = gxe_payload["mean_covariance"]
    result["gxe_covariance_by_environment"] = gxe_payload["covariance_by_environment"]
    result["residual_covariance"] = eps_payload["mean_covariance"]
    result["residual_covariance_by_environment"] = eps_payload["covariance_by_environment"]
    if return_trait_correlations:
        result["genetic_correlation"] = genetic_payload["mean_correlation"]
        result["genetic_correlation_by_environment"] = genetic_payload["correlation_by_environment"]
        result["gxe_correlation"] = gxe_payload["mean_correlation"]
        result["gxe_correlation_by_environment"] = gxe_payload["correlation_by_environment"]
        result["residual_correlation"] = eps_payload["mean_correlation"]
        result["residual_correlation_by_environment"] = eps_payload["correlation_by_environment"]

    return {
        "result": result,
        "fit": {
            "Sigma_G":            mom_out["Sigma_G"],
            "Sigma_GE":           mom_out["Sigma_GE"],
            "Sigma_eps":          mom_out["Sigma_eps"],
            "Sigma_G_response_scale": genetic_payload["mean_covariance"],
            "Sigma_GE_response_scale": gxe_payload["mean_covariance"],
            "Sigma_eps_response_scale": eps_payload["mean_covariance"],
            "Sigma_G_response_scale_by_environment": genetic_payload["covariance_by_environment"],
            "Sigma_GE_response_scale_by_environment": gxe_payload["covariance_by_environment"],
            "Sigma_eps_response_scale_by_environment": eps_payload["covariance_by_environment"],
            "Lambda_G":           mom_out["Lambda_G"],
            "psi_G":              mom_out["psi_G"],
            "Lambda_GE":          mom_out["Lambda_GE"],
            "psi_GE":             mom_out["psi_GE"],
            "beta":               mom_out["beta"],
            "beta_mom":           mom_out["beta_mom"],
            "beta_method":        mom_out["beta_method"],
            "mom_n_full_gids":    mom_out["mom_n_full_gids"],
            "mom_n_full_envs":    mom_out["mom_n_full_envs"],
            "mom_n_full_cells":   mom_out["mom_n_full_cells"],
            "mom_n_observed_cells": mom_out.get("mom_n_observed_cells"),
            "mom_psd_clip_count": mom_out["mom_psd_clip_count"],
            "mom_3x3_cond":       mom_out["mom_3x3_cond"],
            "mom_fallback":       mom_out.get("mom_fallback"),
            "pcg_n_iter":         mom_out["pcg_n_iter"],
            "pcg_residual_norm":  mom_out["pcg_residual_norm"],
            "pcg_converged":      mom_out["pcg_converged"],
            "pcg_preconditioner": mom_out["pcg_preconditioner"],
            "kernel_weights":     mom_out["kernel_weights"],
            "reml_loglik_history": mom_out["reml_loglik_history"],
            "reml_n_iter":         mom_out["reml_n_iter"],
            "reml_converged":      mom_out["reml_converged"],
        },
        "info": {
            "trait_levels":     design["trait_levels"],
            "env_levels":       design["env_levels"],
            "trait_env_means":  design["trait_env_means"],
            "trait_env_sds":    design["trait_env_sds"],
            "kenv_kernel":      kenv_kernel_used,
            "kenv_bandwidth":   kenv_bandwidth_used,
            "kenv_kernel_kwargs": kenv_kernel_kwargs_used,
            "env_covariates_qc": env_covariates_qc_info,
            "wall_s":           float(wall),
            "varcomp_mode":     ("reml_dense" if dense_reml else "mom"),
            "method":           ("mt_gp_met_reml_dense_op_v1" if dense_reml else "mt_gp_met_mom_op_v1"),
            "mom_fallback":     mom_out.get("mom_fallback"),
            "device":           str(torch_device),
            "trait_structure":  trait_structure,
            "gxe_trait_structure": gxe_trait_structure,
            "prediction_se":    bool(return_se),
            "prediction_output": prediction_output,
            "trait_correlations": bool(return_trait_correlations),
            "covariance_scale": "response_by_environment",
        },
    }
