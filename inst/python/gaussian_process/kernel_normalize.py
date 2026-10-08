"""Kernel normalization for GRM and env-similarity matrices.

Modes
-----
- ``"none"``      : pass through unchanged.
- ``"diag_mean"`` : scale each kernel so mean(diag(K)) = 1. Legacy default.
- ``"asreml"``    : identical to ``"none"`` — kernels are used exactly as the
                    caller supplied them. This is the correct mode when
                    variance components should be reported on the original
                    data scale (matching ``vm(GID, GAinv)`` in ASReml-R).
                    Made the default when ``output_level="full_vc"`` so
                    reported VCs are directly comparable.

The module exposes one function; callers choose mode at call time.
"""
from __future__ import annotations

import warnings
from typing import List, Optional, Tuple

import numpy as np


_KNOWN_MODES = {"none", "diag_mean", "asreml", "identity", "mean_centered"}


def apply_kernel_normalization(
    G_list: List[np.ndarray],
    S_e: Optional[np.ndarray],
    mode: str = "diag_mean",
    *,
    output_level: Optional[str] = None,
) -> Tuple[List[np.ndarray], Optional[np.ndarray], str]:
    """Normalize genetic and env-similarity kernels.

    Parameters
    ----------
    G_list : list of (n_geno, n_geno) arrays (GRMs).
    S_e    : optional (n_env, n_env) env-similarity matrix.
    mode   : normalization mode. See module docstring.
    output_level : if ``"full_vc"`` and mode is the legacy ``"diag_mean"``,
        override to ``"asreml"`` with a warning so VC estimates land on the
        original scale.

    Returns
    -------
    G_out  : normalized GRM list.
    Se_out : normalized env-similarity (or None).
    mode_used : the mode actually applied (may differ from input if overridden).
    """
    mode = (mode or "none").lower()

    if mode not in _KNOWN_MODES:
        raise ValueError(f"Unknown kernel normalization mode {mode!r}; "
                         f"choose from {sorted(_KNOWN_MODES)}")

    if output_level == "full_vc" and mode == "diag_mean":
        warnings.warn(
            "reml_normalize='diag_mean' rescales the genetic kernel; "
            "variance components will NOT be on the original y scale. "
            "Overriding to 'asreml' (no rescale) for full_vc output. "
            "Pass reml_normalize='asreml' explicitly to suppress this warning.",
            stacklevel=2,
        )
        mode = "asreml"

    if mode in ("none", "asreml"):
        return list(G_list), S_e, mode

    if mode == "diag_mean":
        G_out, Se_out = _diag_mean(G_list, S_e)
        return G_out, Se_out, mode

    if mode == "identity":
        G_out, Se_out = _identity(G_list, S_e)
        return G_out, Se_out, mode

    if mode == "mean_centered":
        G_out, Se_out = _mean_centered(G_list, S_e)
        return G_out, Se_out, mode

    return list(G_list), S_e, mode


# -- internal helpers --------------------------------------------------------

def _diag_mean(
    G_list: List[np.ndarray],
    S_e: Optional[np.ndarray],
) -> Tuple[List[np.ndarray], Optional[np.ndarray]]:
    G_out = []
    for K in G_list:
        K = np.asarray(K, dtype=float)
        d = np.trace(K) / max(1, K.shape[0])
        if np.isfinite(d) and d > 0:
            G_out.append(K / d)
        else:
            G_out.append(K)
    Se_out = None
    if S_e is not None:
        Se = np.asarray(S_e, dtype=float)
        d = np.trace(Se) / max(1, Se.shape[0])
        if np.isfinite(d) and d > 0:
            Se_out = Se / d
        else:
            Se_out = Se
    return G_out, Se_out


def _identity(
    G_list: List[np.ndarray],
    S_e: Optional[np.ndarray],
) -> Tuple[List[np.ndarray], Optional[np.ndarray]]:
    G_out = []
    for K in G_list:
        K = np.asarray(K, dtype=float)
        d = np.diag(K)
        safe = np.where(d > 0, d, 1.0)
        scale = 1.0 / np.sqrt(safe)
        G_out.append(K * (scale[:, None] * scale[None, :]))
    Se_out = None
    if S_e is not None:
        Se = np.asarray(S_e, dtype=float)
        d = np.diag(Se)
        safe = np.where(d > 0, d, 1.0)
        scale = 1.0 / np.sqrt(safe)
        Se_out = Se * (scale[:, None] * scale[None, :])
    return G_out, Se_out


def _mean_centered(
    G_list: List[np.ndarray],
    S_e: Optional[np.ndarray],
) -> Tuple[List[np.ndarray], Optional[np.ndarray]]:
    G_out = []
    for K in G_list:
        K = np.asarray(K, dtype=float)
        row_mean = K.mean(axis=1, keepdims=True)
        col_mean = K.mean(axis=0, keepdims=True)
        grand_mean = K.mean()
        K = K - row_mean - col_mean + grand_mean
        G_out.append(K)
    return G_out, S_e


__all__ = ["apply_kernel_normalization"]
