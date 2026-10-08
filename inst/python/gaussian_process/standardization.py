"""Per-env y-standardization with scale export for VC back-transform.

The standardizer z-scores the response per environment using train rows.
Predictions and residual VCs are inverted back to original scale via
``inv_mean`` / ``inv_var``.

For variance-component back-transformation (ASReml-schema output):
  - ``get_env_scales(env_labels)`` returns an ``EnvScales`` object that
    ``build_asreml_varcomp_table`` uses to multiply each VC by the
    correct factor: s_e^2 for variance params, s_e for loading params.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

try:
    from varcomp_asreml import EnvScales
except Exception:
    EnvScales = None


class PerObsEnvStandardizer:
    """z-score within environment using TRAIN rows.

    Modes
    -----
    - ``"none"`` : no scaling.
    - ``"global"`` : one train-set mean/std for all observations.
    - ``"per_env"`` : per-env mean/std from train rows.
    """

    def __init__(self, mode: str = "none"):
        self.mode = (mode or "none").lower()
        self.env_means_: Optional[np.ndarray] = None
        self.env_stds_: Optional[np.ndarray] = None
        self.global_mean_: float = 0.0
        self.global_std_: float = 1.0

    def fit(self, y: np.ndarray, ei: np.ndarray, train_idx: np.ndarray):
        y = np.asarray(y, dtype=float)
        ei = np.asarray(ei, dtype=np.int64)
        tr = np.asarray(train_idx, dtype=np.int64)
        y_tr = y[tr]
        if not np.isfinite(y_tr).any():
            self.global_mean_, self.global_std_ = 0.0, 1.0
            self.env_means_ = None
            self.env_stds_ = None
            return self

        self.global_mean_ = float(np.nanmean(y_tr))
        std = float(np.nanstd(y_tr, ddof=1))
        if not np.isfinite(std) or std <= 0:
            std = 1.0
        self.global_std_ = std

        if self.mode == "per_env":
            n_env = int(np.max(ei)) + 1
            means = np.full(n_env, self.global_mean_, dtype=float)
            stds = np.full(n_env, self.global_std_, dtype=float)
            for e in range(n_env):
                idx = tr[ei[tr] == e]
                if idx.size > 1 and np.isfinite(y[idx]).any():
                    m = float(np.nanmean(y[idx]))
                    s = float(np.nanstd(y[idx], ddof=1))
                    if not np.isfinite(s) or s <= 0:
                        s = self.global_std_
                    means[e] = m
                    stds[e] = s
            self.env_means_, self.env_stds_ = means, stds
        else:
            self.env_means_, self.env_stds_ = None, None
        return self

    def transform(self, y: np.ndarray, ei: np.ndarray) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        ei = np.asarray(ei, dtype=np.int64)
        if self.mode == "none":
            return y.copy()
        if self.mode == "per_env" and self.env_means_ is not None:
            means = self.env_means_[ei]
            stds = self.env_stds_[ei]
        else:
            means = np.full_like(y, self.global_mean_, dtype=float)
            stds = np.full_like(y, self.global_std_, dtype=float)
        z = np.zeros_like(y, dtype=float)
        mask = np.isfinite(y)
        z[mask] = (y[mask] - means[mask]) / stds[mask]
        return z

    def inv_mean(self, m: np.ndarray, ei: np.ndarray) -> np.ndarray:
        m = np.asarray(m, dtype=float)
        ei = np.asarray(ei, dtype=np.int64)
        if self.mode == "none":
            return m.copy()
        if self.mode == "per_env" and self.env_means_ is not None:
            return m * self.env_stds_[ei] + self.env_means_[ei]
        return m * self.global_std_ + self.global_mean_

    def inv_var(self, v: np.ndarray, ei: np.ndarray) -> np.ndarray:
        v = np.asarray(v, dtype=float)
        ei = np.asarray(ei, dtype=np.int64)
        if self.mode == "none":
            return v.copy()
        if self.mode == "per_env" and self.env_stds_ is not None:
            return v * (self.env_stds_[ei] ** 2)
        return v * (self.global_std_ ** 2)

    # ---- NEW: scale export for VC back-transform ----

    def get_env_scales(self, env_labels: Sequence[str]) -> "EnvScales":
        """Return per-env scale factors for VC back-transformation.

        Parameters
        ----------
        env_labels : list of environment level names (same order as env index).

        Returns
        -------
        EnvScales with s[label] = the std-dev multiplier used during
        standardization. If mode is "none", all factors are 1.0.

        Raises
        ------
        ImportError if ``varcomp_asreml`` is not available.
        """
        if EnvScales is None:
            raise ImportError(
                "varcomp_asreml.EnvScales is required for get_env_scales; "
                "ensure varcomp_asreml.py is importable."
            )
        labels = [str(e) for e in env_labels]
        if self.mode == "per_env" and self.env_stds_ is not None:
            s_dict: Dict[str, float] = {}
            for i, lab in enumerate(labels):
                if i < len(self.env_stds_):
                    s_dict[lab] = float(self.env_stds_[i])
                else:
                    s_dict[lab] = float(self.global_std_)
            global_s = float(np.sqrt(np.mean(np.square(list(s_dict.values())))))
            return EnvScales(s=s_dict, global_s=global_s)
        if self.mode == "global":
            return EnvScales(
                s={lab: float(self.global_std_) for lab in labels},
                global_s=float(self.global_std_),
            )
        return EnvScales.identity(labels)

    def get_env_stds_array(self) -> np.ndarray:
        """Return raw per-env std-dev array (n_env,). Falls back to global."""
        if self.mode == "none":
            return np.array([1.0], dtype=np.float64)
        if self.env_stds_ is not None:
            return np.array(self.env_stds_, dtype=np.float64)
        return np.array([self.global_std_], dtype=np.float64)


__all__ = ["PerObsEnvStandardizer"]
