"""Reporting helpers shared by multi-trait GP entry points."""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd


def covariance_to_correlation(cov: np.ndarray) -> np.ndarray:
    """Convert a covariance matrix to a numerically stable correlation matrix."""
    cov_arr = np.asarray(cov, dtype=np.float64)
    if cov_arr.ndim != 2 or cov_arr.shape[0] != cov_arr.shape[1]:
        raise ValueError(f"covariance matrix must be square; got shape {cov_arr.shape}")
    cov_arr = 0.5 * (cov_arr + cov_arr.T)
    diag = np.diag(cov_arr)
    sd = np.sqrt(np.clip(diag, 0.0, None))
    corr = np.full_like(cov_arr, np.nan, dtype=np.float64)
    good = sd > 0
    if np.any(good):
        sub = cov_arr[np.ix_(good, good)] / np.outer(sd[good], sd[good])
        sub = np.clip(sub, -1.0, 1.0)
        corr[np.ix_(good, good)] = 0.5 * (sub + sub.T)
        idx = np.flatnonzero(good)
        corr[idx, idx] = 1.0
    return corr


def labeled_matrix(matrix: np.ndarray, labels: Sequence[str]) -> pd.DataFrame:
    """Return a square matrix as a trait-labeled DataFrame."""
    labels = [str(x) for x in labels]
    arr = np.asarray(matrix, dtype=np.float64)
    if arr.shape != (len(labels), len(labels)):
        raise ValueError(
            f"matrix shape {arr.shape} does not match {len(labels)} labels"
        )
    return pd.DataFrame(arr, index=labels, columns=labels)


def covariance_correlation_payload(cov: np.ndarray, labels: Sequence[str]) -> dict:
    """Return labeled covariance and correlation matrices for a trait covariance."""
    corr = covariance_to_correlation(cov)
    return {
        "covariance": labeled_matrix(cov, labels),
        "correlation": labeled_matrix(corr, labels),
    }


def covariance_to_response_scale(cov: np.ndarray, scales: Sequence[float]) -> np.ndarray:
    """Back-transform a standardized trait covariance with positive trait scales."""
    cov_arr = np.asarray(cov, dtype=np.float64)
    scale_arr = np.asarray(scales, dtype=np.float64).reshape(-1)
    if cov_arr.shape != (scale_arr.size, scale_arr.size):
        raise ValueError(
            f"covariance shape {cov_arr.shape} does not match {scale_arr.size} scales"
        )
    if not np.isfinite(scale_arr).all() or (scale_arr <= 0).any():
        raise ValueError("response scales must be finite and positive")
    out = np.outer(scale_arr, scale_arr) * cov_arr
    return 0.5 * (out + out.T)


def covariance_by_environment_payload(
    cov: np.ndarray,
    scales_by_environment: np.ndarray,
    trait_labels: Sequence[str],
    environment_labels: Sequence[str],
) -> Dict[str, object]:
    """Return response-scale covariance/correlation matrices for every environment."""
    scales = np.asarray(scales_by_environment, dtype=np.float64)
    env_labels = [str(x) for x in environment_labels]
    if scales.shape != (len(env_labels), len(trait_labels)):
        raise ValueError(
            "scales_by_environment shape must equal "
            f"(n_environment, n_trait); got {scales.shape}"
        )
    covariances = {}
    correlations = {}
    covariance_arrays = []
    for env, scale_row in zip(env_labels, scales):
        response_cov = covariance_to_response_scale(cov, scale_row)
        payload = covariance_correlation_payload(response_cov, trait_labels)
        covariances[env] = payload["covariance"]
        correlations[env] = payload["correlation"]
        covariance_arrays.append(response_cov)
    mean_covariance = np.mean(np.stack(covariance_arrays, axis=0), axis=0)
    mean_payload = covariance_correlation_payload(mean_covariance, trait_labels)
    # The mean covariance is a legitimate magnitude summary, but its correlation
    # is not the model's trait correlation. Each environment's response-scale
    # covariance is cov scaled by outer(s_e, s_e), so every per-environment
    # correlation equals correlation(cov) exactly. Averaging the covariances
    # first attenuates the off-diagonals by Cauchy-Schwarz whenever the trait
    # scales are not proportional across environments, which would report a
    # "mean" correlation that matches none of the per-environment values it
    # claims to summarise. Report the scale-invariant correlation instead.
    model_payload = covariance_correlation_payload(
        np.asarray(cov, dtype=np.float64), trait_labels
    )
    return {
        "covariance_by_environment": covariances,
        "correlation_by_environment": correlations,
        "mean_covariance": mean_payload["covariance"],
        "mean_correlation": model_payload["correlation"],
    }


# cov_traits costs one PCG solve per test cell per trait. Past a few hundred
# cells that is hours, so it is capped. Aborting the whole fit is the wrong
# response: predictions and variance components are still valid and usually the
# point of the run. Follow cv_harness.cv_predict -- warn loudly, drop the
# expensive block, and say so in the returned metadata.
COV_TRAITS_SOLVE_BUDGET = 1000


def cov_traits_budget_exceeded(n_cells: int, n_traits: int) -> bool:
    """True when computing cov_traits would exceed the PCG solve budget."""
    return int(n_cells) * int(n_traits) > COV_TRAITS_SOLVE_BUDGET


def resolve_cov_traits_request(
    *,
    n_cells: int,
    n_traits: int,
    requested: bool,
    force_cov_traits: bool = False,
    force_prediction_se: bool = False,
    context: str = "multi-trait GP",
):
    """Decide whether to compute cov_traits, downgrading instead of raising.

    Returns a dict with `compute`, `downgraded` and `reason`. `downgraded` is
    True only when the caller asked for cov_traits and the budget refused it, so
    the caller can record that the block is absent by decision rather than by
    accident.
    """
    import warnings

    if not requested:
        return {"compute": False, "downgraded": False, "reason": ""}
    if bool(force_cov_traits) or bool(force_prediction_se):
        return {"compute": True, "downgraded": False, "reason": ""}
    if not cov_traits_budget_exceeded(n_cells, n_traits):
        return {"compute": True, "downgraded": False, "reason": ""}

    n_solves = int(n_cells) * int(n_traits)
    reason = (
        f"{context}: cov_traits at n_cells={int(n_cells)}, T={int(n_traits)} "
        f"would require {n_solves} PCG solves (budget "
        f"{COV_TRAITS_SOLVE_BUDGET}). Returning predictions without per-cell "
        f"trait covariance, prediction SE or trait correlations. To compute "
        f"them anyway, force it (R: gp_force_prediction_se = TRUE; Python: "
        f"force_cov_traits=True / force_prediction_se=True), or reduce the "
        f"test set."
    )
    warnings.warn(reason, RuntimeWarning, stacklevel=2)
    return {"compute": False, "downgraded": True, "reason": reason}
