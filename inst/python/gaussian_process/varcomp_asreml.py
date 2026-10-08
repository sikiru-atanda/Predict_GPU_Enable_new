"""ASReml-compatible variance-component output and AI-REML backend protocol. v2

Produces the canonical `varcomp` DataFrame and `summary` dict consumed by the
framework. All VC work is strictly opt-in: callers must explicitly request
`output_level="full_vc"` to trigger any code path that reaches into this module.

Schema matches `summary(asreml_fit)$varcomp`:

    component    str   e.g. "fa(Env,1):vm(GID,GAinv)!E03!var", "...!E03!fa1",
                       "Env_E03!R"
    estimate     float on original y scale (after per-env standardization is
                       inverted)
    std.error    float AI-based asymptotic SE; NaN for bound/fixed rows
    z.ratio      float estimate / std.error; NaN when std.error is NaN
    bound        str   "P" positive, "B" boundary-clipped at 0, "F" fixed, "U"
                       unconstrained
    pct_change   float last AI-step relative change; NaN at convergence floor

The module is import-light (numpy + pandas + typing). No torch import at module
scope, so `from varcomp_asreml import ...` is cheap on the fast path.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Literal, Optional, Protocol, Sequence, Tuple

import math

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Parameter spec
# ---------------------------------------------------------------------------

ParamKind = Literal["fa_psi", "fa_lambda", "resid_env", "var_kernel", "cor", "fixed"]
BoundFlag = Literal["P", "B", "F", "U"]


@dataclass(frozen=True)
class ParamSpec:
    """One scalar VC parameter in the AI vector.

    Fields are sufficient to (a) place the parameter in the AI matrix, (b)
    label it in the output table using ASReml conventions, and (c) apply the
    correct scale back-transform when per-env standardization was used.
    """
    kind: ParamKind
    term: str                    # e.g. "fa(Env,1):vm(GID,GAinv)" or "Env"
    env_label: Optional[str] = None
    sub_label: Optional[str] = None  # "var" | "fa1" | "fa2" | "R" | ...
    lower: float = -math.inf
    upper: float = math.inf
    fixed: bool = False

    def component_name(self) -> str:
        if self.kind == "resid_env":
            return f"{self.term}_{self.env_label}!R"
        if self.env_label is None and self.sub_label is None:
            return self.term
        pieces = [self.term]
        if self.env_label is not None:
            pieces.append(self.env_label)
        if self.sub_label is not None:
            pieces.append(self.sub_label)
        return "!".join(pieces)

    def default_bound(self) -> BoundFlag:
        if self.fixed:
            return "F"
        if self.lower > -math.inf and self.kind in ("fa_psi", "resid_env", "var_kernel"):
            return "P"
        return "U"


# ---------------------------------------------------------------------------
# Scale back-transform
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EnvScales:
    """Per-env scale factors from `_PerObsEnvStandardizer`.

    `s[env_label]` is the y-standard-deviation multiplier applied during fit.
    To invert: multiply variance-like params by s^2, loading-like params by s.
    If standardization was not used, pass `EnvScales.identity(env_labels)`.
    """
    s: Dict[str, float]
    global_s: float = 1.0

    @classmethod
    def identity(cls, env_labels: Sequence[str]) -> "EnvScales":
        return cls(s={str(e): 1.0 for e in env_labels}, global_s=1.0)

    def factor_for(self, spec: ParamSpec) -> float:
        s = (
            float(self.global_s)
            if spec.env_label is None
            else float(self.s.get(str(spec.env_label), self.global_s))
        )
        if spec.kind in ("fa_psi", "resid_env", "var_kernel"):
            return s * s
        if spec.kind == "fa_lambda":
            return s
        return 1.0


# ---------------------------------------------------------------------------
# Active set / AI inversion
# ---------------------------------------------------------------------------

def active_set_inverse(
    ai: np.ndarray,
    bounds: Sequence[BoundFlag],
    jitter: float = 1e-10,
) -> np.ndarray:
    """Invert AI restricted to the active (non-B, non-F) set.

    Parameters
    ----------
    ai : (p, p) AI matrix in natural parameter ordering.
    bounds : length-p flags. "B" and "F" are inactive; SE for those rows is NaN.
    jitter : ridge added before Cholesky if the active block is ill-conditioned.

    Returns
    -------
    cov : (p, p) with NaN in rows/cols for inactive params, inverse AI in the
          active sub-block embedded in its original positions.
    """
    p = int(ai.shape[0])
    active = np.array([b not in ("B", "F") for b in bounds], dtype=bool)
    cov = np.full((p, p), np.nan, dtype=np.float64)
    if not active.any():
        return cov
    idx = np.where(active)[0]
    sub = np.asarray(ai[np.ix_(idx, idx)], dtype=np.float64)
    sub = 0.5 * (sub + sub.T)
    for tries in range(6):
        try:
            L = np.linalg.cholesky(sub + (jitter * (10.0 ** tries)) * np.eye(sub.shape[0]))
            inv = np.linalg.solve(L.T, np.linalg.solve(L, np.eye(sub.shape[0])))
            break
        except np.linalg.LinAlgError:
            inv = None
    if inv is None:
        inv = np.linalg.pinv(sub)
    cov[np.ix_(idx, idx)] = inv
    return cov


def apply_boundary_projection(
    theta: np.ndarray,
    specs: Sequence[ParamSpec],
    floor_eps: float = 1e-12,
) -> Tuple[np.ndarray, List[BoundFlag]]:
    """Clip theta to [lower, upper] per spec; return projected theta and bounds.

    A parameter clipped at its lower bound 0 is tagged "B"; fixed params "F";
    otherwise inherits `default_bound()`. `floor_eps` is the tolerance used to
    decide "at the boundary" vs "interior".
    """
    out = np.array(theta, dtype=np.float64, copy=True)
    bounds: List[BoundFlag] = []
    for i, sp in enumerate(specs):
        if sp.fixed:
            bounds.append("F")
            continue
        if out[i] < sp.lower:
            out[i] = sp.lower
        if out[i] > sp.upper:
            out[i] = sp.upper
        if (sp.lower > -math.inf) and (out[i] - sp.lower <= floor_eps) and sp.kind in (
            "fa_psi", "resid_env", "var_kernel"
        ):
            bounds.append("B")
        else:
            bounds.append(sp.default_bound())
    return out, bounds


def _optimization_active_set(
    bounds: Sequence[BoundFlag],
    score: np.ndarray,
) -> np.ndarray:
    """Return the KKT-aware active set used during constrained optimisation.

    Final covariance/SE reporting still treats boundary parameters as inactive.
    During optimisation, however, a lower-bound parameter with positive score
    must be allowed to re-enter the feasible interior.
    """
    score = np.asarray(score, dtype=np.float64).reshape(-1)
    if score.size != len(bounds):
        raise ValueError("score length must match the variance-component bounds")
    return np.array(
        [
            False if bound == "F" else (
                float(score[j]) > 0.0 if bound == "B" else True
            )
            for j, bound in enumerate(bounds)
        ],
        dtype=bool,
    )


# ---------------------------------------------------------------------------
# Output table assembly
# ---------------------------------------------------------------------------

def build_asreml_varcomp_table(
    specs: Sequence[ParamSpec],
    theta: np.ndarray,
    ai: np.ndarray,
    bounds: Sequence[BoundFlag],
    pct_change: Sequence[float],
    env_scales: Optional[EnvScales] = None,
) -> pd.DataFrame:
    """Assemble the ASReml-style varcomp DataFrame.

    Applies per-env scale back-transform to estimates and SEs (delta method:
    SE scales by the same factor as the estimate since the transform is
    linear in the param).
    """
    p = len(specs)
    if theta.shape[0] != p or ai.shape != (p, p) or len(bounds) != p or len(pct_change) != p:
        raise ValueError("specs / theta / ai / bounds / pct_change length mismatch")

    cov = active_set_inverse(ai, bounds)
    se = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
    se = np.where(np.isfinite(se), se, np.nan)

    scales = np.array(
        [1.0 if env_scales is None else env_scales.factor_for(sp) for sp in specs],
        dtype=np.float64,
    )
    est_out = np.asarray(theta, dtype=np.float64) * scales
    se_out = se * scales  # linear back-transform

    # A near-singular active AI block makes sqrt(diag(AI^-1)) explode, giving a
    # huge, meaningless SE (z-ratio ~ 0). Treat such components as not estimable
    # and report NaN rather than a misleading number. z_floor = 1e-3 means an SE
    # more than ~1000x the estimate is considered non-informative.
    z_floor = 1e-3
    rows = []
    for i, sp in enumerate(specs):
        b = bounds[i]
        est = float(est_out[i])
        s = float(se_out[i]) if (b not in ("B", "F") and np.isfinite(se_out[i])) else float("nan")
        if np.isfinite(s) and s > 0.0 and abs(est) > 0.0 and (abs(est) / s) < z_floor:
            s = float("nan")
        z = (est / s) if (s == s and s > 0.0) else float("nan")
        pc = float(pct_change[i]) if np.isfinite(pct_change[i]) else float("nan")
        rows.append({
            "component": sp.component_name(),
            "estimate": est,
            "std.error": s,
            "z.ratio": z,
            "bound": b,
            "pct_change": pc,
        })
    return pd.DataFrame(rows, columns=["component", "estimate", "std.error", "z.ratio", "bound", "pct_change"])


# ---------------------------------------------------------------------------
# FA-derived variance components: delta-method propagation through (psi, lambda)
# ---------------------------------------------------------------------------

def compute_fa_derived_summary(
    specs: Sequence[ParamSpec],
    theta: np.ndarray,
    ai: np.ndarray,
    bounds: Sequence[BoundFlag],
    env_scales: Optional[EnvScales] = None,
    kernel_weight_total: float = 1.0,
) -> Optional[Dict[str, object]]:
    """Delta-method SE for FA-derived sigma2_g (and h2).

    For an FA(rank=r) covariance, the per-env genetic variance is

        sigma2_g_e = w_K * (sum_k lambda_{e,k}^2 + psi_e)

    where ``w_K`` is the sum of the effective weights on the normalized
    genomic kernels. FA parameters describe the environment covariance
    structure; the genomic-kernel scale is a separate multiplicative part of
    the fitted covariance and must be included in reported genetic variance.

    which is NOT a direct REML parameter. Its asymptotic SE requires
    delta-method propagation through the (fa_psi, fa_lambda) block of the
    AI inverse:

        SE(sigma2_g_e) = sqrt( J_e^T  Sigma_theta  J_e )

    where J_e is the Jacobian vector with entries
        J_e[psi_e]       = 1
        J_e[lambda_{e,k}] = 2 * lambda_{e,k}
        all other entries = 0

    The total sigma2_g (matching `var_components_summary` "genetic_main_total")
    is sum_e sigma2_g_e and uses J_total = sum_e J_e.

    Per-env scale back-transform is applied via `env_scales`: variance-like
    parameters back-transform by s_e^2, loading-like by s_e. Chain rule
    multiplies each Jacobian entry by the corresponding scale.

    Returns None if `specs` contains no `fa_psi` entries (i.e. not an FA fit).

    Notes
    -----
    Pattern follows the canonical delta-method assembly in the
    Mixed_Model_Project reference (api.py:1315-1323): Cholesky-conditioned
    AI inverse + J^T Sigma J. The reference assembles the variance
    structure differently (gamma=sigma_i^2/sigma_e^2); here we propagate
    through the FA (loadings, specific-variances) parameterisation.

    Output schema
    -------------
    {
        "env_labels": [str, ...],          # length n_env
        "sigma2_g_per_env": [float, ...],  # length n_env
        "sigma2_g_se_per_env": [float, ...],
        "sigma2_g_total": float,
        "sigma2_g_total_se": float,
        "h2_per_env": [float, ...],        # NaN if no resid_env match
        "h2_se_per_env": [float, ...],
        "h2_average": float,
        "h2_average_se": float,            # joint delta method
    }
    """
    kernel_weight_total = float(kernel_weight_total)
    if not np.isfinite(kernel_weight_total) or kernel_weight_total < 0.0:
        raise ValueError("kernel_weight_total must be finite and non-negative")

    fa_psi_specs: List[Tuple[int, ParamSpec]] = [
        (i, sp) for i, sp in enumerate(specs) if sp.kind == "fa_psi"
    ]
    if not fa_psi_specs:
        return None

    fa_lambda_specs: List[Tuple[int, ParamSpec]] = [
        (i, sp) for i, sp in enumerate(specs) if sp.kind == "fa_lambda"
    ]
    resid_specs: List[Tuple[int, ParamSpec]] = [
        (i, sp) for i, sp in enumerate(specs) if sp.kind == "resid_env"
    ]

    env_labels = [str(sp.env_label) for _, sp in fa_psi_specs]
    n_env = len(env_labels)
    n_specs = len(specs)

    scales = np.array(
        [1.0 if env_scales is None else env_scales.factor_for(sp) for sp in specs],
        dtype=np.float64,
    )
    theta_np = np.asarray(theta, dtype=np.float64)
    est_scaled = theta_np * scales

    # Active-set AI inverse: rows/cols of bound params are NaN-filled, which
    # propagates to any derived SE that depends on them.
    ai_inv = active_set_inverse(np.asarray(ai, dtype=np.float64), bounds)

    # Index psi and residual rows by env_label for O(1) lookup; lambda indices
    # collected per env (rank may be >= 1).
    psi_idx_by_env = {str(sp.env_label): i for i, sp in fa_psi_specs}
    resid_idx_by_env = {str(sp.env_label): i for i, sp in resid_specs}
    lambda_idxs_by_env: Dict[str, List[int]] = {lbl: [] for lbl in env_labels}
    for i, sp in fa_lambda_specs:
        lbl = str(sp.env_label)
        if lbl in lambda_idxs_by_env:
            lambda_idxs_by_env[lbl].append(i)

    # z-floor for derived SE matches build_asreml_varcomp_table convention:
    # an SE so large that |est|/SE < 1e-3 is statistically meaningless, and
    # the table-level convention is to report NaN rather than the absurd
    # number. The same convention applies to delta-method-derived SE.
    _Z_FLOOR = 1e-3

    def _delta_se(J: np.ndarray, est: float = float("nan")) -> float:
        """SE(derived) = sqrt(J^T Sigma_theta J), NaN if any contributing entry
        of J hits an inactive (NaN) row of the AI inverse OR if the resulting
        z-ratio is below the global z-floor (1e-3).

        Restrict the bilinear form to the contributing index set so that NaN
        entries in `ai_inv` from bound (inactive) parameters do NOT propagate
        through `0 * NaN = NaN` for indices where J is zero.

        `est` is the corresponding point estimate; when supplied and finite,
        the z-floor is applied (matching `build_asreml_varcomp_table`).
        """
        contributing = np.where(J != 0.0)[0]
        if contributing.size == 0:
            return float("nan")
        sub = ai_inv[contributing, :][:, contributing]
        if np.any(np.isnan(sub)):
            return float("nan")
        J_sub = J[contributing]
        V = float(J_sub @ sub @ J_sub)
        if not np.isfinite(V) or V <= 0.0:
            return float("nan")
        se = float(np.sqrt(V))
        if np.isfinite(est) and abs(est) > 0.0 and (abs(est) / se) < _Z_FLOOR:
            return float("nan")
        return se

    sigma2_g_per_env: List[float] = []
    sigma2_g_se_per_env: List[float] = []
    h2_per_env: List[float] = []
    h2_se_per_env: List[float] = []
    h2_jacobians: List[Optional[np.ndarray]] = []
    J_total = np.zeros(n_specs, dtype=np.float64)

    for env_label in env_labels:
        psi_i = psi_idx_by_env[env_label]
        psi_e = float(est_scaled[psi_i])

        lam_is = lambda_idxs_by_env[env_label]
        lam_e_scaled = est_scaled[lam_is] if lam_is else np.zeros(0)
        sigma2_g_e = float(kernel_weight_total * (np.sum(lam_e_scaled ** 2) + psi_e))
        sigma2_g_per_env.append(sigma2_g_e)

        # Jacobian for sigma2_g_e (scaled to natural parameter scale)
        J_e = np.zeros(n_specs, dtype=np.float64)
        # d(scaled_psi_e)/d(theta_psi_e) = scale_psi_e -> chain in scale
        J_e[psi_i] = kernel_weight_total * scales[psi_i]
        for li in lam_is:
            # d((lambda_li * scale_li)^2) / d(theta_lambda_li) = 2 * lambda_li_scaled * scale_li
            J_e[li] = kernel_weight_total * 2.0 * est_scaled[li] * scales[li]
        sigma2_g_se_per_env.append(_delta_se(J_e, est=sigma2_g_e))
        J_total += J_e

        # h2_e = sigma2_g_e / (sigma2_g_e + sigma2_resid_e)
        if env_label in resid_idx_by_env:
            r_i = resid_idx_by_env[env_label]
            sigma2_r_e = float(est_scaled[r_i])
            total_e = sigma2_g_e + sigma2_r_e
            if np.isfinite(total_e) and total_e > 0.0:
                h2_e = sigma2_g_e / total_e
                # Chain rule:
                #   dh2/dtheta = (sigma_r / total^2) * J_e   - (sigma_g / total^2) * J_r
                J_r = np.zeros(n_specs, dtype=np.float64)
                J_r[r_i] = scales[r_i]
                J_h2 = (sigma2_r_e / total_e ** 2) * J_e - (sigma2_g_e / total_e ** 2) * J_r
                h2_se_e = _delta_se(J_h2, est=h2_e)
            else:
                h2_e = float("nan")
                h2_se_e = float("nan")
                J_h2 = None
        else:
            h2_e = float("nan")
            h2_se_e = float("nan")
            J_h2 = None
        h2_per_env.append(h2_e)
        h2_se_per_env.append(h2_se_e)
        h2_jacobians.append(J_h2)

    sigma2_g_total = float(np.sum(sigma2_g_per_env))
    sigma2_g_total_se = _delta_se(J_total, est=sigma2_g_total)
    valid_h2 = [
        i for i, (value, jacobian) in enumerate(zip(h2_per_env, h2_jacobians))
        if np.isfinite(value) and jacobian is not None
    ]
    if valid_h2:
        h2_average = float(np.mean([h2_per_env[i] for i in valid_h2]))
        J_h2_average = np.mean(
            np.stack([h2_jacobians[i] for i in valid_h2]), axis=0
        )
        h2_average_se = _delta_se(J_h2_average, est=h2_average)
    else:
        h2_average = float("nan")
        h2_average_se = float("nan")

    return {
        "env_labels": env_labels,
        "sigma2_g_per_env": sigma2_g_per_env,
        "sigma2_g_se_per_env": sigma2_g_se_per_env,
        "sigma2_g_total": sigma2_g_total,
        "sigma2_g_total_se": sigma2_g_total_se,
        "h2_per_env": h2_per_env,
        "h2_se_per_env": h2_se_per_env,
        "h2_average": h2_average,
        "h2_average_se": h2_average_se,
        "kernel_weight_total": kernel_weight_total,
    }


def compute_genetic_residual_covariance(
    specs: Sequence[ParamSpec],
    theta: np.ndarray,
    ai: np.ndarray,
    bounds: Sequence[BoundFlag],
    env_scales: Optional[EnvScales] = None,
) -> Optional[Dict[str, object]]:
    """Active-set AI inverse off-diagonal for non-FA models.

    For models where sigma2_g and sigma2_e are BOTH direct REML parameters
    (`var_kernel` specs for genetic, `resid_env` specs for residual), the
    standard 2-variance delta-method for h2 SE assumes
    `Cov(sigma2_g, sigma2_e) = 0`. The off-diagonal of the active-set AI
    inverse gives the real covariance; including it in the h2 SE formula
    is straightforward delta-method:

        Var(h2) = (dh2/dg)^2 Var(g) + (dh2/de)^2 Var(e)
                  + 2 (dh2/dg)(dh2/de) Cov(g, e)

    Returns `None` if `specs` does not contain at least one `var_kernel`
    and one `resid_env` (e.g. FA models -- those use
    `compute_fa_derived_summary`).

    Schema
    ------
    {
      "sigma2_g":        float,       # estimate, response scale
      "sigma2_e":        float,       # estimate, response scale (mean across resid_env)
      "var_g":           float,       # ai_inv[g, g] * scale_g^2
      "var_e":           float,       # ai_inv[e, e] * scale_e^2  (mean across resid_env)
      "cov_g_e":         float,       # ai_inv[g, e] * scale_g * scale_e
    }

    For multi-env (multiple resid_env specs), this returns the mean
    residual variance and the average of per-env covariances. Single-env
    is the common case; multi-env results are a defensible aggregate.
    """
    var_kernel_idx = [i for i, sp in enumerate(specs) if sp.kind == "var_kernel"]
    resid_idx = [i for i, sp in enumerate(specs) if sp.kind == "resid_env"]
    if not var_kernel_idx or not resid_idx:
        return None

    scales = np.array(
        [1.0 if env_scales is None else env_scales.factor_for(sp) for sp in specs],
        dtype=np.float64,
    )
    theta_np = np.asarray(theta, dtype=np.float64)
    ai_inv = active_set_inverse(np.asarray(ai, dtype=np.float64), bounds)

    # Genetic: if multiple var_kernel rows (multi-kernel main term), aggregate
    # the variance contributions. For the standard single-kernel case this is
    # just the one row.
    g_indices = np.asarray(var_kernel_idx, dtype=np.int64)
    g_scales = scales[g_indices]
    sigma2_g = float(np.nansum(theta_np[g_indices] * g_scales))
    # Var(sum) = sum_i sum_j scale_i * scale_j * AI_inv[i, j]
    g_sub = ai_inv[g_indices, :][:, g_indices]
    if np.any(np.isnan(g_sub)):
        var_g = float("nan")
    else:
        var_g = float((g_scales[:, None] * g_scales[None, :] * g_sub).sum())
        if not np.isfinite(var_g) or var_g < 0.0:
            var_g = float("nan")

    # Residual: report the mean across env-specific residuals
    r_indices = np.asarray(resid_idx, dtype=np.int64)
    r_scales = scales[r_indices]
    sigma2_e_per = theta_np[r_indices] * r_scales
    sigma2_e = float(np.nanmean(sigma2_e_per))
    # Var(mean) = (1/n^2) * sum_i sum_j scale_i * scale_j * AI_inv[i, j]
    r_sub = ai_inv[r_indices, :][:, r_indices]
    if np.any(np.isnan(r_sub)):
        var_e = float("nan")
    else:
        n_r = float(len(r_indices))
        var_e = float((r_scales[:, None] * r_scales[None, :] * r_sub).sum()) / (n_r ** 2)
        if not np.isfinite(var_e) or var_e < 0.0:
            var_e = float("nan")

    # Cov(sigma2_g_sum, sigma2_e_mean) = (1/n_r) * sum_i sum_j scale_i_g * scale_j_e * AI_inv[i, j]
    cross = ai_inv[g_indices, :][:, r_indices]
    if np.any(np.isnan(cross)):
        cov_g_e = float("nan")
    else:
        n_r = float(len(r_indices))
        cov_g_e = float((g_scales[:, None] * r_scales[None, :] * cross).sum()) / n_r
        if not np.isfinite(cov_g_e):
            cov_g_e = float("nan")

    return {
        "sigma2_g": sigma2_g,
        "sigma2_e": sigma2_e,
        "var_g": var_g,
        "var_e": var_e,
        "cov_g_e": cov_g_e,
    }


def build_summary(
    *,
    loglik: float,
    nedf: int,
    parameters_free: int,
    n_ai_iter: int,
    converged: bool,
    sigma: float = 1.0,
    extra: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """Build the top-level summary dict aligned with ASReml conventions.

    AIC / BIC use the number of *free* (non-B, non-F) parameters.
    """
    k = int(parameters_free)
    n = int(nedf)
    aic = -2.0 * float(loglik) + 2.0 * k
    bic = -2.0 * float(loglik) + math.log(max(n, 1)) * k
    out: Dict[str, object] = {
        "loglik": float(loglik),
        "nedf": n,
        "sigma": float(sigma),
        "aic": float(aic),
        "bic": float(bic),
        "parameters": k,
        "converged": bool(converged),
        "n_ai_iter": int(n_ai_iter),
    }
    if extra:
        out.update(extra)
    return out


# ---------------------------------------------------------------------------
# VC backend protocol (Phase A: Dense; Phase B: Operator+SLQ)
# ---------------------------------------------------------------------------

class VCBackend(Protocol):
    """Minimal interface every VC backend must supply.

    Each method is called only when `output_level="full_vc"`. The protocol
    keeps AI-REML orchestration decoupled from the linear algebra engine.
    """

    def loglik(self, theta: np.ndarray) -> float: ...
    def score(self, theta: np.ndarray) -> np.ndarray: ...
    def ai(self, theta: np.ndarray) -> np.ndarray: ...
    def param_specs(self) -> List[ParamSpec]: ...


@dataclass
class DenseVCBackend:
    """Dense AI-REML backend shim.

    This is a thin adapter. It delegates to callables supplied by the
    framework (`compute_AI_Matrix_SEs`, REML-loglik helpers) so we reuse
    existing kernels without duplicating math. Wiring lives in the
    framework file; this class only enforces the protocol shape.
    """
    specs: List[ParamSpec]
    loglik_fn: Callable[[np.ndarray], float]
    score_fn: Callable[[np.ndarray], np.ndarray]
    ai_fn: Callable[[np.ndarray], np.ndarray]

    def loglik(self, theta: np.ndarray) -> float:
        return float(self.loglik_fn(np.asarray(theta, dtype=np.float64)))

    def score(self, theta: np.ndarray) -> np.ndarray:
        return np.asarray(self.score_fn(np.asarray(theta, dtype=np.float64)), dtype=np.float64)

    def ai(self, theta: np.ndarray) -> np.ndarray:
        return np.asarray(self.ai_fn(np.asarray(theta, dtype=np.float64)), dtype=np.float64)

    def param_specs(self) -> List[ParamSpec]:
        return list(self.specs)


class OperatorVCBackend:
    """Stub for Phase B: SLQ logdet + Hutchinson AI on the large-operator path.

    Intentionally raises NotImplementedError. Phase B implements logdet and
    AI traces via probe vectors against the operator matvec interfaces in
    `mixed_model_gpu_large_scale_backend_*` without changing callers.
    """

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def _notimpl(self):
        raise NotImplementedError(
            "OperatorVCBackend (SLQ+Hutchinson AI) is a Phase B deliverable. "
            "For now, use output_level='predict_only' or 'predict_with_se' on "
            "the large-operator path, or run the dense engine for full_vc."
        )

    def loglik(self, theta: np.ndarray) -> float: self._notimpl()
    def score(self, theta: np.ndarray) -> np.ndarray: self._notimpl()
    def ai(self, theta: np.ndarray) -> np.ndarray: self._notimpl()
    def param_specs(self) -> List[ParamSpec]: self._notimpl()


# ---------------------------------------------------------------------------
# AI-REML driver (backend-agnostic)
# ---------------------------------------------------------------------------

@dataclass
class AIRemlResult:
    theta: np.ndarray
    ai: np.ndarray
    bounds: List[BoundFlag]
    pct_change: np.ndarray
    loglik: float
    n_iter: int
    converged: bool


def ai_reml_optimize(
    backend: VCBackend,
    theta0: np.ndarray,
    *,
    max_iter: int = 50,
    tol_loglik: float = 1e-4,
    tol_theta: float = 1e-6,
    damping: float = 1.0,
    max_halvings: int = 8,
    verbose: bool = False,
) -> AIRemlResult:
    """Boundary-aware AI-Newton with adaptive trust-region damping.

    Uses a Levenberg-Marquardt style rho ratio to adapt the trust region:
      rho = actual_improvement / predicted_improvement
    where predicted_improvement accounts for the damped step size.

    Damping policy (matches ASReml-R / reference Mixed_Model_Project):
      rho > 0.75  -> expand trust region (mu *= 2, capped at 1.0)
      0.25 < rho <= 0.75 -> keep current mu
      0 < rho <= 0.25 -> shrink (mu /= 2)
      rho <= 0 -> reject step, shrink hard (mu /= 4)

    Falls back to steepest-ascent after 10 consecutive rejections to escape
    pathological curvature regions.
    """
    specs = backend.param_specs()
    theta = np.asarray(theta0, dtype=np.float64).copy()
    theta, bounds = apply_boundary_projection(theta, specs)
    last_ll = backend.loglik(theta)
    last_theta = theta.copy()
    pct = np.full_like(theta, np.nan)
    ai_final = np.zeros((theta.size, theta.size), dtype=np.float64)
    converged = False
    mu = float(damping)
    consecutive_rejects = 0

    # Optional one-shot finite-difference check at theta0.
    # Separate env var from AIREML_VERBOSE since it costs 3*q extra loglik
    # calls -- useful when diagnosing score/loglik consistency bugs
    # (e.g. the stochastic-loglik bug fixed 2026-04-18) but otherwise wasted.
    import os as _os
    if verbose and int(_os.environ.get("AIREML_FD_CHECK", "0")):
        try:
            s0 = backend.score(theta)
            ll0 = backend.loglik(theta)
            print(f"  [FD-check] ll0={ll0:.6f}  analytic score={np.round(s0, 3).tolist()}", flush=True)
            for eps in (1e-3, 1e-5, 1e-7):
                fd = np.zeros_like(theta)
                lls = []
                for j in range(len(theta)):
                    tp = theta.copy(); tp[j] += eps
                    tp2, _ = apply_boundary_projection(tp, specs)
                    ll_p = backend.loglik(tp2)
                    lls.append(ll_p)
                    fd[j] = (ll_p - ll0) / eps
                print(f"  [FD-check] eps={eps:.0e} ll_perturbed[0]={lls[0]:.6f} diff[0]={lls[0]-ll0:+.6e}", flush=True)
                print(f"  [FD-check] eps={eps:.0e} fd={np.round(fd, 3).tolist()}", flush=True)
        except Exception as _fd_e:
            print(f"  [FD-check] skipped ({_fd_e})", flush=True)

    for it in range(max_iter):
        ai = backend.ai(theta)
        s = backend.score(theta)
        ai_final = ai.copy()

        # A projected lower-bound parameter is inactive only when its score
        # satisfies the maximisation KKT condition (score <= 0).  A positive
        # score at the lower bound means that moving back into the feasible
        # interior increases the REML log-likelihood, so the component must be
        # allowed to re-enter the Newton system.  Treating every "B" flag as
        # permanently inactive can trap AI-REML at the all-residual solution
        # after one over-sized projected step.
        active = _optimization_active_set(bounds, s)
        if not active.any():
            converged = True
            break

        idx = np.where(active)[0]
        A = ai[np.ix_(idx, idx)]
        A = 0.5 * (A + A.T)
        s_active = s[idx]

        # LM ridge for ill-conditioned AI
        diag_A = np.diag(A).copy()
        pos_diag = np.abs(diag_A[diag_A != 0])
        cond_est = (np.max(pos_diag) / np.min(pos_diag)) if pos_diag.size > 0 else 1.0
        if cond_est > 1e8 or np.any(diag_A <= 0):
            ridge = max(1e-6 * np.max(np.abs(diag_A)), 1e-10)
        else:
            ridge = 1e-10

        # Solve for full Newton step
        try:
            delta_active = np.linalg.solve(A + ridge * np.eye(A.shape[0]), s_active)
        except np.linalg.LinAlgError:
            delta_active = np.linalg.lstsq(A, s_active, rcond=None)[0]

        # Predicted improvement from full Newton step
        predicted_full = float(s_active @ delta_active
                               - 0.5 * delta_active @ A @ delta_active)

        # Apply damped step
        trial = theta.copy()
        trial[idx] = theta[idx] + mu * delta_active
        trial, trial_bounds = apply_boundary_projection(trial, specs)
        trial_ll = backend.loglik(trial)
        actual = trial_ll - last_ll

        # Rho ratio: actual / predicted for the damped step
        predicted_damped = predicted_full * mu * (2.0 - mu)
        if abs(predicted_damped) > 1e-15:
            rho = actual / predicted_damped
        else:
            rho = 1.0 if actual >= 0 else 0.0

        # Trust-region update
        accepted = False
        if rho < 0:
            # Reject step, shrink trust region hard
            mu = max(mu / 4.0, 1e-4)
            consecutive_rejects += 1
        elif rho < 0.25:
            # Accept but shrink trust region
            accepted = True
            mu = max(mu / 2.0, 1e-4)
        else:
            # Accept
            accepted = True
            if rho > 0.75:
                mu = min(mu * 2.0, 1.0)

        if accepted:
            with np.errstate(invalid="ignore", divide="ignore"):
                pct = np.where(
                    np.abs(last_theta) > 1e-12,
                    100.0 * (trial - last_theta) / np.where(
                        np.abs(last_theta) > 1e-12, last_theta, 1.0),
                    np.nan,
                )
            theta = trial
            bounds = trial_bounds
            last_theta = theta.copy()
            last_ll = trial_ll
            consecutive_rejects = 0

        # Steepest-ascent fallback after many consecutive rejections.
        # Lowered threshold from 10 -> 5 so it actually fires when mu sits at
        # the 1e-4 floor with a non-trivial score.
        if consecutive_rejects >= 5:
            max_diag = np.max(np.abs(diag_A))
            if max_diag > 0:
                grad_step = s_active / max_diag
                trial_sa = theta.copy()
                trial_sa[idx] += 0.01 * grad_step
                trial_sa, trial_sa_bounds = apply_boundary_projection(trial_sa, specs)
                trial_sa_ll = backend.loglik(trial_sa)
                if trial_sa_ll > last_ll:
                    theta = trial_sa
                    bounds = trial_sa_bounds
                    last_theta = theta.copy()
                    last_ll = trial_sa_ll
            consecutive_rejects = 0
            mu = 0.1

        if verbose:
            score_norm = float(np.max(np.abs(s_active))) if s_active.size > 0 else 0.0
            max_pct_dbg = float(np.nanmax(np.abs(pct[np.isfinite(pct)]))) if np.isfinite(pct).any() else 0.0
            print(f"  iter={it+1:3d} ll={last_ll:.6f} mu={mu:.4g} "
                  f"rho={rho:+.3f} actual={actual:+.3e} acc={accepted} "
                  f"|s|={score_norm:.3e} max_pct={max_pct_dbg:.3e} "
                  f"theta[:5]={np.round(theta[:5], 5).tolist()}",
                  flush=True)

        # Convergence check (only when step was accepted)
        if accepted:
            finite_pct = pct[np.isfinite(pct)]
            max_pct = float(np.nanmax(np.abs(finite_pct))) if finite_pct.size > 0 else 0.0
            if abs(actual) < tol_loglik and max_pct < 100.0 * tol_theta:
                converged = True
                break

    ai_final = backend.ai(theta)
    return AIRemlResult(
        theta=theta,
        ai=ai_final,
        bounds=bounds,
        pct_change=pct,
        loglik=float(last_ll),
        n_iter=it + 1,
        converged=converged,
    )


__all__ = [
    "ParamKind", "BoundFlag", "ParamSpec",
    "EnvScales",
    "apply_boundary_projection", "active_set_inverse",
    "build_asreml_varcomp_table", "build_summary",
    "VCBackend", "DenseVCBackend", "OperatorVCBackend",
    "AIRemlResult", "ai_reml_optimize",
]
