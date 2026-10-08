"""K_env (reaction-norm environment kernel) math, extracted from the framework.

Public:
  build_kenv(Z, kernel="matern32", bandwidth=1.0, **kwargs) -> (K, info)
    Returns (kernel_matrix, info_dict). info_dict has keys:
      'kernel', 'effective_kwargs', 'median_dist', 'fallback'.
  validate_kwargs(kernel, kwargs) -> None
    Raises ValueError on unknown keys, type errors, or out-of-range values.

Kernels registered: linear, rbf, matern32, matern52, polynomial, pca_rbf.

Directly-called with uncentered Z, pca_rbf's SVD projection is NOT equivalent
to PCA (it assumes centered input). The framework preprocessor guarantees
centering; direct-module callers must pre-center Z themselves.
"""
from __future__ import annotations
from typing import Any, Callable, Dict, Optional, Set, Tuple
import warnings
import numpy as np


# {kernel_name: {"fn": callable, "allowed": set[str], "defaults": dict}}
_KERNEL_REGISTRY: Dict[str, Dict[str, Any]] = {}


def validate_kwargs(kernel: str, kwargs: Dict[str, Any]) -> None:
    """Strict: raise ValueError on unknown/invalid keys."""
    spec = _KERNEL_REGISTRY.get(kernel)
    if spec is None:
        raise ValueError(
            f"kenv_kernels.validate_kwargs: unknown kernel {kernel!r}; "
            f"registered: {sorted(_KERNEL_REGISTRY.keys())}"
        )
    allowed: Set[str] = spec["allowed"]
    unknown = [k for k in kwargs if k not in allowed]
    if unknown:
        raise ValueError(
            f"kenv_kernel_kwargs: unknown key {unknown[0]!r} for kernel {kernel!r}. "
            f"Valid keys: {sorted(allowed)}. Got kwargs={kwargs}."
        )


def build_kenv(
    Z: np.ndarray,
    kernel: str = "matern32",
    bandwidth: float = 1.0,
    **kwargs: Any,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Build K_env. Returns (K, info). Kernels are registered below."""
    spec = _KERNEL_REGISTRY.get(kernel)
    if spec is None:
        raise ValueError(
            f"kenv_kernels.build_kenv: unknown kernel {kernel!r}; "
            f"registered: {sorted(_KERNEL_REGISTRY.keys())}"
        )
    validate_kwargs(kernel, kwargs)
    Z = np.asarray(Z, dtype=np.float64)
    info: Dict[str, Any] = {
        "kernel": kernel,
        "effective_kwargs": {},
        "median_dist": None,
        "fallback": None,
    }
    K = spec["fn"](Z, bandwidth=float(bandwidth), info=info, **kwargs)
    K = 0.5 * (K + K.T) + 1e-6 * np.eye(K.shape[0])
    return K, info


# -------------------------------------------------- legacy kernels

def _kenv_linear(Z: np.ndarray, bandwidth: float, info: Dict[str, Any]) -> np.ndarray:
    # bandwidth unused: linear kernel has no distance scale
    p = Z.shape[1]
    K = (Z @ Z.T) / float(max(p, 1))
    return K


def _pairwise_sq_dist_and_median(Z: np.ndarray) -> Tuple[np.ndarray, float]:
    """Returns (D2, med_d) where D2 is the pairwise squared Euclidean distance matrix
    and med_d is the median of the upper-triangle Euclidean distances."""
    sq = (Z * Z).sum(axis=1)
    D2 = sq[:, None] + sq[None, :] - 2.0 * (Z @ Z.T)
    D2 = np.clip(D2, 0.0, None)
    n = Z.shape[0]
    iu = np.triu_indices(n, k=1)
    med_d = float(np.median(np.sqrt(D2[iu]))) if iu[0].size else 1.0
    return D2, med_d


def _kenv_rbf(Z: np.ndarray, bandwidth: float, info: Dict[str, Any]) -> np.ndarray:
    D2, med_d = _pairwise_sq_dist_and_median(Z)
    info["median_dist"] = med_d
    h = max(float(bandwidth) * med_d, 1e-8)
    K = np.exp(-D2 / (2.0 * h * h))
    return K


def _kenv_matern32(Z: np.ndarray, bandwidth: float, info: Dict[str, Any]) -> np.ndarray:
    D2, med_d = _pairwise_sq_dist_and_median(Z)
    info["median_dist"] = med_d
    h = max(float(bandwidth) * med_d, 1e-8)
    D = np.sqrt(D2)
    r = D / h
    K = (1.0 + np.sqrt(3.0) * r) * np.exp(-np.sqrt(3.0) * r)
    return K


def _kenv_matern52(Z: np.ndarray, bandwidth: float, info: Dict[str, Any]) -> np.ndarray:
    D2, med_d = _pairwise_sq_dist_and_median(Z)
    info["median_dist"] = med_d
    h = max(float(bandwidth) * med_d, 1e-8)
    D = np.sqrt(D2)
    r = D / h
    K = (1.0 + np.sqrt(5.0) * r + 5.0 * r * r / 3.0) * np.exp(-np.sqrt(5.0) * r)
    return K


_KERNEL_REGISTRY["linear"]   = {"fn": _kenv_linear,   "allowed": set(), "defaults": {}}
_KERNEL_REGISTRY["rbf"]      = {"fn": _kenv_rbf,      "allowed": set(), "defaults": {}}
_KERNEL_REGISTRY["matern32"] = {"fn": _kenv_matern32, "allowed": set(), "defaults": {}}
_KERNEL_REGISTRY["matern52"] = {"fn": _kenv_matern52, "allowed": set(), "defaults": {}}


# -------------------------------------------------- polynomial kernel

def _kenv_polynomial(
    Z: np.ndarray,
    bandwidth: float,
    info: Dict[str, Any],
    degree: int = 2,
    c: float = 1.0,
) -> np.ndarray:
    # bandwidth unused: polynomial kernel has no distance scale
    # Validate types + ranges.
    if not isinstance(degree, (int, np.integer)) or isinstance(degree, bool):
        raise ValueError(f"polynomial: degree must be int; got {type(degree).__name__}")
    if int(degree) <= 0:
        raise ValueError(f"polynomial: degree must be > 0; got {degree}")
    if not isinstance(c, (int, float, np.integer, np.floating)) or isinstance(c, bool):
        raise ValueError(f"polynomial: c must be numeric; got {type(c).__name__}")
    if float(c) < 0:
        raise ValueError(f"polynomial: c must be >= 0; got {c}")
    d = int(degree)
    cc = float(c)
    p = Z.shape[1]
    base = (Z @ Z.T) / float(max(p, 1)) + cc
    K = np.power(base, d)
    info["effective_kwargs"] = {"degree": d, "c": cc}
    return K


_KERNEL_REGISTRY["polynomial"] = {
    "fn": _kenv_polynomial,
    "allowed": {"degree", "c"},
    "defaults": {"degree": 2, "c": 1.0},
}


# -------------------------------------------------- pca_rbf kernel

def _kenv_pca_rbf(
    Z: np.ndarray,
    bandwidth: float,
    info: Dict[str, Any],
    pca_k: Optional[int] = None,
) -> np.ndarray:
    n, p = Z.shape
    if pca_k is None:
        pca_k_req = min(p, max(5, n // 5))
    else:
        if not isinstance(pca_k, (int, np.integer)) or isinstance(pca_k, bool):
            raise ValueError(f"pca_rbf: pca_k must be int; got {type(pca_k).__name__}")
        if int(pca_k) <= 0:
            raise ValueError(f"pca_rbf: pca_k must be > 0; got {pca_k}")
        pca_k_req = int(pca_k)
    pca_k_eff = min(pca_k_req, p)
    if pca_k_eff != pca_k_req:
        info["fallback"] = (
            f"pca_k clipped from {pca_k_req} to {pca_k_eff} (p={p})"
        )
    if pca_k_eff == p and p > 0:
        # Degenerate: projection reconstructs Z fully => equivalent to RBF.
        warnings.warn(
            f"pca_rbf: effective pca_k ({pca_k_eff}) == p ({p}); kernel reduces "
            f"to rbf (no denoising). Pass pca_k < p for real PCA denoising.",
            RuntimeWarning, stacklevel=3,
        )
        prior = info.get("fallback") or ""
        info["fallback"] = (
            (prior + "; " if prior else "")
            + "pca_rbf effectively reduces to rbf (pca_k == p)"
        )
    info["effective_kwargs"] = {"pca_k": pca_k_eff}

    if n < 2:
        # n == 0: return 0x0 below via power path; n == 1: return 1x1 = [[1.0]].
        K = np.ones((n, n), dtype=np.float64)
        return K

    # SVD-based PC projection. Assumes Z is centered (framework preprocessor guarantees).
    U, S, _Vt = np.linalg.svd(Z, full_matrices=False)
    Z_pc = U[:, :pca_k_eff] * S[:pca_k_eff][None, :]

    # RBF on reduced space.
    sq = (Z_pc * Z_pc).sum(axis=1)
    D2 = sq[:, None] + sq[None, :] - 2.0 * (Z_pc @ Z_pc.T)
    D2 = np.clip(D2, 0.0, None)
    iu = np.triu_indices(n, k=1)
    D = np.sqrt(D2)
    med_d = float(np.median(D[iu])) if iu[0].size else 1.0
    info["median_dist"] = med_d
    h = max(float(bandwidth) * med_d, 1e-8)
    K = np.exp(-D2 / (2.0 * h * h))
    return K


_KERNEL_REGISTRY["pca_rbf"] = {
    "fn": _kenv_pca_rbf,
    "allowed": {"pca_k"},
    "defaults": {"pca_k": None},  # None => auto-infer in helper
}
