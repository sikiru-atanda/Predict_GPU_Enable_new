"""Multi-trait GP, single environment (v1).

Public entry: fit_multi_trait_gp(...).
Joint multi-trait GP for one environment.
"""
from __future__ import annotations
import os
import time
from typing import Any, Dict, Optional, Sequence

import numpy as np
import pandas as pd

from mt_gp_internals import (
    _fit_mt_ai_reml_multikernel,
    _build_mt_design,
    _initialize_params,
    _fit_mt_ai_reml,
    _predict_mt,
)
from mt_gp_reporting import covariance_correlation_payload, covariance_to_response_scale
from gp_device import pick_torch_device


def _resolve_mt_torch_device(
    device: Optional[Any],
    varcomp_mode: str,
) -> tuple:
    """Resolve the MT-GP device without auto-routing Hessians to CUDA.

    PyTorch's second-order ``torch.func.hessian`` path used by the dense
    AI-REML models can leave a CUDA process unusable after an illegal-memory
    fault.  Auto device selection therefore keeps Hessian-based REML on CPU.
    Scalable MoM fitting retains the ordinary GP device policy, and an
    explicit ``cuda``/``cuda:N`` request is still respected.
    """
    requested_raw = (
        device
        if device is not None
        else os.environ.get("PREDICTPRO_GP_DEVICE", "auto")
    )
    requested = str(requested_raw or "auto").strip().lower()
    mode = str(varcomp_mode or "mom").strip().lower()
    if mode != "mom" and requested in ("", "auto", "gpu"):
        return pick_torch_device("cpu"), "cpu_second_order_autograd_safety"
    reason = "explicit_device" if requested not in ("", "auto", "gpu") else "standard_auto_policy"
    return pick_torch_device(device), reason


def fit_multi_trait_gp(
    *,
    pheno_df: pd.DataFrame,
    gid_col: str,
    trait_col: str,
    y_col: str,
    geno_kernels: Dict[str, np.ndarray],
    geno_ids: Sequence[str],
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    trait_structure: str = "unstructured",
    trait_fa_rank: Optional[int] = None,
    residual_structure: Optional[str] = None,
    residual_fa_rank: Optional[int] = None,
    fixed_effects: Sequence[str] = (),
    kernel_weights: Optional[Sequence[float]] = None,
    estimate_kernel_weights: bool = False,
    max_iter: int = 100,
    tol_loglik: float = 1e-4,
    seed: int = 0,
    dtype: str = "float64",
    device: Optional[str] = None,
    verbose: bool = False,
    # --- v2-fast MoM kwargs ---
    varcomp_mode: str = "mom",          # "mom" (fast, default) | "reml" (v1 AI-REML)
    grm_factor_cache: Any = None,        # GRMFactorCacheBase (memmap or zarr); required when varcomp_mode="mom"
    mom_pcg_tol: float = 1e-4,
    mom_pcg_max_iter: int = 300,
    mom_psd_floor: float = 1e-4,
    return_se: bool = True,
    prediction_output: str = "test_only",
    return_cov_traits: bool = False,
    force_cov_traits: bool = False,
    force_prediction_se: bool = False,
    return_trait_correlations: bool = False,
) -> Dict[str, Any]:
    """Fit single-env multi-trait GP with joint AI-REML.

    See README.md in this folder.
    """
    t0 = time.perf_counter()
    mode = str(varcomp_mode).strip().lower()
    torch_device, device_policy = _resolve_mt_torch_device(device, mode)
    prediction_output = str(prediction_output).strip().lower()
    if prediction_output not in {"test_only", "all"}:
        raise ValueError(
            "prediction_output must be one of {'test_only', 'all'}; "
            f"got {prediction_output!r}"
        )
    if mode == "mom":
        if grm_factor_cache is None:
            raise ValueError(
                "varcomp_mode='mom' requires grm_factor_cache (memmap/zarr Phi); "
                "see the large-n notes in inst/python/gaussian_process/README.md."
            )
        from mt_gp_internals import _fit_mt_mom_op
        design = _build_mt_design(
            pheno_df=pheno_df, gid_col=gid_col, trait_col=trait_col, y_col=y_col,
            fixed_effects=tuple(fixed_effects),
            geno_ids=list(geno_ids),
            train_idx=np.asarray(train_idx, dtype=np.int64),
            test_idx=np.asarray(test_idx, dtype=np.int64),
        )
        if prediction_output == "all":
            design["test_mask"] = np.ones_like(design["test_mask"], dtype=bool)
        mom_out = _fit_mt_mom_op(
            design=design, cache=grm_factor_cache,
            kernel_weights=kernel_weights,
            psd_floor=float(mom_psd_floor),
            pcg_tol=float(mom_pcg_tol),
            pcg_max_iter=int(mom_pcg_max_iter),
            return_cov_traits=bool(return_cov_traits),
            force_cov_traits=bool(force_cov_traits),
            force_prediction_se=bool(force_prediction_se),
            device=torch_device,
        )
        preds_df = mom_out["predictions"].rename(columns={"gid": gid_col, "trait": trait_col})
        result = {
            "predictions": preds_df,
            "cov_traits":  mom_out["cov_traits"],
        }
        if return_trait_correlations:
            Sigma_G_response = covariance_to_response_scale(
                mom_out["Sigma_G"], design["trait_sds"]
            )
            Sigma_eps_response = covariance_to_response_scale(
                mom_out["Sigma_eps"], design["trait_sds"]
            )
            corr_payload = covariance_correlation_payload(
                Sigma_G_response, design["trait_levels"]
            )
            result["genetic_covariance"] = corr_payload["covariance"]
            result["genetic_correlation"] = corr_payload["correlation"]
            eps_payload = covariance_correlation_payload(
                Sigma_eps_response, design["trait_levels"]
            )
            result["residual_covariance"] = eps_payload["covariance"]
            result["residual_correlation"] = eps_payload["correlation"]
        wall = time.perf_counter() - t0
        return {
            "result": result,
            "fit": {
                "Sigma_G":            mom_out["Sigma_G"],
                "Sigma_eps":          mom_out["Sigma_eps"],
                "kernel_weights":     mom_out["kernel_weights"],
                "Sigma_G_response_scale": covariance_to_response_scale(
                    mom_out["Sigma_G"], design["trait_sds"]
                ),
                "Sigma_eps_response_scale": covariance_to_response_scale(
                    mom_out["Sigma_eps"], design["trait_sds"]
                ),
                "beta":               mom_out["beta"],
                "mom_n_full":         mom_out["mom_n_full"],
                "mom_n_partial":      mom_out["mom_n_partial"],
                "mom_psd_clip_count": mom_out["mom_psd_clip_count"],
                "pcg_n_iter":         mom_out["pcg_n_iter"],
                "pcg_residual_norm":  mom_out["pcg_residual_norm"],
            },
            "info": {
                "trait_levels":  design["trait_levels"],
                "trait_means":   design["trait_means"],
                "trait_sds":     design["trait_sds"],
                "wall_s":        float(wall),
                "varcomp_mode":  "mom",
                "method":        "mt_gp_mom_op_v2",
                "device":        str(torch_device),
                "device_policy": device_policy,
                "prediction_se": False,
                "prediction_output": prediction_output,
                "trait_correlations": bool(return_trait_correlations),
                "covariance_scale": "response",
            },
        }
    # Else: fall through to existing v1 REML body (unchanged).
    if "A" not in geno_kernels:
        raise ValueError("geno_kernels must contain key 'A' for the additive K_geno")
    kernel_names = list(geno_kernels)
    if kernel_weights is None:
        weights = np.ones(len(kernel_names), dtype=np.float64)
    else:
        weights = np.asarray(kernel_weights, dtype=np.float64).reshape(-1)
    if weights.size != len(kernel_names):
        raise ValueError(
            f"kernel_weights length {weights.size} != number of geno_kernels={len(kernel_names)}"
        )
    if not np.isfinite(weights).all() or (weights < 0).any() or not (weights > 0).any():
        raise ValueError("kernel_weights must be finite, non-negative, and include a positive value")
    K_geno = np.zeros_like(np.asarray(geno_kernels["A"], dtype=np.float64))
    for name, weight in zip(kernel_names, weights.tolist()):
        one = np.asarray(geno_kernels[name], dtype=np.float64)
        if one.shape != K_geno.shape:
            raise ValueError(f"geno_kernels[{name!r}] shape {one.shape} != {K_geno.shape}")
        K_geno += float(weight) * one
    if K_geno.shape != (len(geno_ids), len(geno_ids)):
        raise ValueError(
            f"K_geno shape {K_geno.shape} != (n_gid={len(geno_ids)}, n_gid={len(geno_ids)})"
        )
    if residual_structure is None:
        residual_structure = trait_structure
    if residual_fa_rank is None and residual_structure == "fa":
        residual_fa_rank = trait_fa_rank

    design = _build_mt_design(
        pheno_df=pheno_df, gid_col=gid_col, trait_col=trait_col, y_col=y_col,
        fixed_effects=tuple(fixed_effects),
        geno_ids=list(geno_ids),
        train_idx=np.asarray(train_idx, dtype=np.int64),
        test_idx=np.asarray(test_idx, dtype=np.int64),
    )
    if prediction_output == "all":
        design["test_mask"] = np.ones_like(design["test_mask"], dtype=bool)
    T = design["T"]
    if trait_structure == "fa":
        if trait_fa_rank is None:
            raise ValueError("trait_fa_rank required when trait_structure='fa'")
        if int(trait_fa_rank) >= T:
            raise ValueError(
                f"trait_fa_rank={trait_fa_rank} must be < T={T} for identifiability"
            )
    if residual_structure == "fa":
        if residual_fa_rank is None:
            raise ValueError("residual_fa_rank required when residual_structure='fa'")
        if int(residual_fa_rank) >= T:
            raise ValueError(
                f"residual_fa_rank={residual_fa_rank} must be < T={T} for identifiability"
            )

    Sigma_G0, Sigma_eps0 = _initialize_params(
        y_std=design["y_std"], X=design["X_fixed"],
        gid_idx=design["gid_idx"], trait_idx=design["trait_idx"],
        train_mask=design["train_mask"], T=T, n_gid=design["n_gid"],
    )
    # The bank is kept separate so the mixture weights can be fitted rather than
    # assumed. With estimation off this reduces to the historical single
    # K_geno = sum_k w_k K_k fit with the supplied (or equal) weights.
    kernel_bank = [np.asarray(geno_kernels[name], dtype=np.float64)
                   for name in kernel_names]
    fit = _fit_mt_ai_reml_multikernel(
        kernel_bank=kernel_bank,
        kernel_weights=weights,
        estimate_kernel_weights=bool(estimate_kernel_weights),
        y_std=design["y_std"], X=design["X_fixed"],
        gid_idx=design["gid_idx"], trait_idx=design["trait_idx"],
        train_mask=design["train_mask"],
        T=T, n_gid=design["n_gid"],
        Sigma_G_init=Sigma_G0, Sigma_eps_init=Sigma_eps0,
        trait_structure=trait_structure, trait_fa_rank=trait_fa_rank,
        residual_structure=residual_structure, residual_fa_rank=residual_fa_rank,
        max_iter=int(max_iter), tol_loglik=float(tol_loglik),
        dtype=str(dtype),
        device=torch_device,
    )
    weights = np.asarray(fit["kernel_weights"], dtype=np.float64)
    # Predictions must use the mixture that was actually fitted.
    K_geno = np.zeros_like(K_geno)
    for name, weight in zip(kernel_names, weights.tolist()):
        K_geno = K_geno + float(weight) * np.asarray(
            geno_kernels[name], dtype=np.float64
        )
    preds_df, cov_traits = _predict_mt(
        K_geno=K_geno, fit=fit, design=design,
        gid_col=gid_col, trait_col=trait_col,
        device=torch_device,
    )
    if not return_se:
        preds_df = preds_df.drop(
            columns=[
                "SE",
                "SE_latent",
                "PEV",
                "Prediction_Var_latent",
                "Prediction_Var_observed",
            ],
            errors="ignore",
        )
    result = {
        "predictions": preds_df,
        "cov_traits": cov_traits,
    }
    if return_trait_correlations:
        Sigma_G_response = covariance_to_response_scale(fit["Sigma_G"], design["trait_sds"])
        Sigma_eps_response = covariance_to_response_scale(fit["Sigma_eps"], design["trait_sds"])
        corr_payload = covariance_correlation_payload(Sigma_G_response, design["trait_levels"])
        result["genetic_covariance"] = corr_payload["covariance"]
        result["genetic_correlation"] = corr_payload["correlation"]
        eps_payload = covariance_correlation_payload(Sigma_eps_response, design["trait_levels"])
        result["residual_covariance"] = eps_payload["covariance"]
        result["residual_correlation"] = eps_payload["correlation"]
    wall = time.perf_counter() - t0
    return {
        "result": result,
        "fit": {
            "Sigma_G":    fit["Sigma_G"],
            "Sigma_eps":  fit["Sigma_eps"],
            "kernel_weights": weights,
            "kernel_names": list(kernel_names),
            "kernel_weights_estimated": bool(fit.get("kernel_weights_estimated", False)),
            "kernel_weight_loglik_gain": float(fit.get("kernel_weight_loglik_gain", 0.0)),
            "kernel_weight_profile_evals": int(fit.get("kernel_weight_profile_evals", 0)),
            "kernel_weight_at_bound": bool(fit.get("kernel_weight_at_bound", False)),
            "Sigma_G_response_scale": covariance_to_response_scale(
                fit["Sigma_G"], design["trait_sds"]
            ),
            "Sigma_eps_response_scale": covariance_to_response_scale(
                fit["Sigma_eps"], design["trait_sds"]
            ),
            "Lambda_G":   fit.get("Lambda_G"),
            "psi_G":      fit.get("psi_G"),
            "Lambda_eps": fit.get("Lambda_eps"),
            "psi_eps":    fit.get("psi_eps"),
            "beta":       fit["beta"],
            "loglik_history": fit["loglik_history"],
            "n_iter":     fit["n_iter"],
            "converged":  fit["converged"],
            "kernel_eigen_clip_count": fit.get("kernel_eigen_clip_count", 0),
            "kernel_min_eigenvalue_input": fit.get("kernel_min_eigenvalue_input"),
            "line_search_invalid_candidate_rejections": fit.get(
                "line_search_invalid_candidate_rejections", 0
            ),
        },
        "info": {
            "trait_levels": design["trait_levels"],
            "trait_means":  design["trait_means"],
            "trait_sds":    design["trait_sds"],
            "wall_s":       float(wall),
            "varcomp_mode": "reml",
            "method":       "mt_gp_dense_v1",
            "device":       str(torch_device),
            "device_policy": device_policy,
            "prediction_se": bool(return_se),
            "prediction_output": prediction_output,
            "trait_correlations": bool(return_trait_correlations),
            "covariance_scale": "response",
        },
    }
