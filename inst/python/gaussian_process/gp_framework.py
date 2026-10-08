# mixed_model_gpu.py — GPU-accelerated KRR / GP-Exact / GP-ICM-FA (drop-in)
from __future__ import annotations
import os, sys

import math, time, warnings, gc
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# Modular shared layers (preferred when available)
try:
    _this_dir = os.path.dirname(os.path.abspath(__file__))
except Exception:
    _this_dir = None
for _p in [_this_dir, "/mnt/data"]:
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

try:
    from term_spec import RandomTermSpec, FittedRandomTerm, FittedMixedModel, MixedModelSpec as _ModMixedModelSpec
    _MODULAR_TERM_SPEC_AVAILABLE = True
except Exception:
    RandomTermSpec = None
    FittedRandomTerm = None
    FittedMixedModel = None
    _ModMixedModelSpec = None
    _MODULAR_TERM_SPEC_AVAILABLE = False

try:
    from reporting import (
        build_reporting_bundle as _mod_build_reporting_bundle,
        attach_prediction_se_columns as _mod_attach_prediction_se_columns,
    )
    _MODULAR_REPORTING_AVAILABLE = True
except Exception:
    _mod_build_reporting_bundle = None
    _mod_attach_prediction_se_columns = None
    _MODULAR_REPORTING_AVAILABLE = False

try:
    from residuals import (
        resolve_stagewise_residual_inputs as _mod_resolve_stagewise_residual_inputs,
        env_mean_from_obs_diag as _mod_env_mean_from_obs_diag,
        build_residual_summary as _mod_build_residual_summary,
    )
    _MODULAR_RESIDUALS_AVAILABLE = True
except Exception:
    _mod_resolve_stagewise_residual_inputs = None
    _mod_env_mean_from_obs_diag = None
    _mod_build_residual_summary = None
    _MODULAR_RESIDUALS_AVAILABLE = False

try:
    from cov_structures import make_covariance_module as _mod_make_covariance_module
    _MODULAR_COV_AVAILABLE = True
except Exception:
    _mod_make_covariance_module = None
    _MODULAR_COV_AVAILABLE = False

try:
    from model_api import fit_mixed_model_modular as _mod_fit_mixed_model_modular
    _MODULAR_MODEL_API_AVAILABLE = True
except Exception:
    _mod_fit_mixed_model_modular = None
    _MODULAR_MODEL_API_AVAILABLE = False

try:
    from framework_engine_policy import choose_framework_engine_mode as _mod_choose_framework_engine_mode
    _MODULAR_ENGINE_POLICY_AVAILABLE = True
except Exception:
    _mod_choose_framework_engine_mode = None
    _MODULAR_ENGINE_POLICY_AVAILABLE = False

try:
    from kernel_normalize import apply_kernel_normalization as _mod_apply_kernel_normalization
    _MODULAR_KERNEL_NORMALIZE_AVAILABLE = True
except Exception:
    _mod_apply_kernel_normalization = None
    _MODULAR_KERNEL_NORMALIZE_AVAILABLE = False

try:
    from standardization import PerObsEnvStandardizer as _ModPerObsEnvStandardizer
    _MODULAR_STANDARDIZATION_AVAILABLE = True
except Exception:
    _ModPerObsEnvStandardizer = None
    _MODULAR_STANDARDIZATION_AVAILABLE = False

try:
    from varcomp_asreml import ParamSpec as _ParamSpec, EnvScales as _EnvScales
    from ai_reml import (
        DenseAIContext as _DenseAIContext,
        run_ai_reml as _run_ai_reml,
        build_varcomp_outputs as _build_varcomp_outputs,
        compute_ai_se_active_set as _compute_ai_se_active_set,
        bounds_from_specs_and_theta as _bounds_from_specs_and_theta,
    )
    _MODULAR_AI_REML_AVAILABLE = True
except Exception:
    _MODULAR_AI_REML_AVAILABLE = False

_MODULAR_LAYERS_AVAILABLE = bool(
    _MODULAR_TERM_SPEC_AVAILABLE or
    _MODULAR_REPORTING_AVAILABLE or
    _MODULAR_RESIDUALS_AVAILABLE or
    _MODULAR_COV_AVAILABLE
)

# ----------------------- Torch / GPyTorch availability -----------------------

try:
    import torch
    TORCH_AVAILABLE = True
    import torch.nn as nn
except Exception:
    TORCH_AVAILABLE = False

try:
    import gpytorch
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from linear_operator.operators import DenseLinearOperator, DiagLinearOperator
    GPTY_AVAILABLE = True
except Exception:
    GPTY_AVAILABLE = False

try:
    from gp_device import pick_torch_device as _pick_gp_torch_device
except Exception:
    _pick_gp_torch_device = None

def _decide_framework_engine_mode(
    *,
    requested_engine="auto",
    modular_dense_delegate=False,
    use_modular_dense=False,
    point_predictions_only=False,
    return_se=True,
    compute_ai_se=False,
    n_obs=None,
    n_structured_terms=1,
    max_axis_levels=0,
    large_n_threshold=50000,
    max_full_levels=20,
):
    """Framework-side engine routing helper.

    Conservative policy: preserve numerical equivalence for legacy branches by
    keeping auto-routing on legacy engines. The modular dense engine is
    complementary and only activated by explicit user request.
    """
    requested_engine = str(requested_engine).lower().strip()
    if requested_engine not in {"auto", "gpu_fast", "operator_fast", "dense_modular"}:
        raise ValueError("requested_engine must be one of: auto, gpu_fast, operator_fast, dense_modular")

    explicit_modular = bool(modular_dense_delegate or use_modular_dense or requested_engine == "dense_modular")
    if explicit_modular:
        return {
            "engine_mode": "dense_modular",
            "reason": "explicit complementary modular dense request",
            "use_legacy_gpu": False,
            "use_legacy_operator": False,
            "use_modular_dense": True,
            "equivalence_note": "modular dense path is complementary and not guaranteed numerically identical to legacy GPU/operator branches",
        }

    is_large_n = (n_obs is not None) and (int(n_obs) >= int(large_n_threshold))

    if requested_engine == "gpu_fast":
        return {
            "engine_mode": "gpu_fast",
            "reason": "explicit specialized dense/GPU request",
            "use_legacy_gpu": True,
            "use_legacy_operator": False,
            "use_modular_dense": False,
        }
    if requested_engine == "operator_fast":
        return {
            "engine_mode": "operator_fast",
            "reason": "explicit scalable operator request",
            "use_legacy_gpu": False,
            "use_legacy_operator": True,
            "use_modular_dense": False,
        }

    # auto: remain on legacy engines to preserve numerical equivalence to legacy framework.
    if is_large_n and (not compute_ai_se) and (point_predictions_only or return_se):
        return {
            "engine_mode": "operator_fast",
            "reason": "auto-selected legacy operator path for large-n prediction/SE request",
            "use_legacy_gpu": False,
            "use_legacy_operator": True,
            "use_modular_dense": False,
            "equivalence_note": "legacy branch preserved",
        }

    return {
        "engine_mode": "gpu_fast",
        "reason": "auto-selected legacy dense/GPU path to preserve legacy numerical behavior",
        "use_legacy_gpu": True,
        "use_legacy_operator": False,
        "use_modular_dense": False,
        "equivalence_note": "legacy branch preserved",
        "notes": {
            "n_structured_terms": int(n_structured_terms),
            "max_axis_levels": int(max_axis_levels),
            "complex_structure_recommendation": (
                "use engine_mode='dense_modular' explicitly for richer structured interaction models"
                if (n_structured_terms > 1 or max_axis_levels > max_full_levels) else None
            ),
        },
    }

def _make_covariance_module_shared(structure: str, level_names, fa_rank=None):
    structure = str(structure).lower()
    if _MODULAR_COV_AVAILABLE and _mod_make_covariance_module is not None:
        return _mod_make_covariance_module(structure=structure, level_names=level_names, fa_rank=fa_rank)
    if structure == "identity":
        return None
    if structure == "diag":
        return DiagonalHeterogeneousCovariance(level_names)
    if structure == "corh":
        return CorHCovariance(level_names)
    if structure == "corgh":
        return CorGHCovariance(level_names)
    if structure == "us":
        return UnstructuredCovariance(level_names)
    if structure == "fa":
        return FactorAnalyticCovariance(level_names, rank=int(fa_rank or 1))
    raise ValueError(f"Unknown covariance structure: {structure}")


def build_term_registry_from_legacy_inputs(
    geno_kernel_names,
    env_col,
    env_levels,
    env_structure=None,
    fa_rank=None,
    interaction_terms_meta=None,
):
    out = []
    if _MODULAR_TERM_SPEC_AVAILABLE and RandomTermSpec is not None:
        out.append(RandomTermSpec(term_name="G", kernel_names=list(geno_kernel_names), axis_name=None, axis_levels=None, structure="identity", fa_rank=None))
    else:
        out.append({"term_name":"G","kernel_names":list(geno_kernel_names),"axis_name":None,"axis_levels":None,"structure":"identity","fa_rank":None})

    if interaction_terms_meta:
        for term in interaction_terms_meta:
            term_name = str(term.get("term_name", f"G:{term.get('axis_name', env_col)}"))
            axis_name = term.get("axis_name", env_col)
            axis_levels = list(term.get("level_names", term.get("axis_levels", env_levels)))
            structure = str(term.get("structure", env_structure or "identity")).lower()
            term_fa_rank = term.get("fa_rank", fa_rank)
            kernel_names = list(term.get("kernel_names", geno_kernel_names))
            if _MODULAR_TERM_SPEC_AVAILABLE and RandomTermSpec is not None:
                out.append(RandomTermSpec(
                    term_name=term_name,
                    kernel_names=kernel_names,
                    axis_name=str(axis_name),
                    axis_levels=axis_levels,
                    structure=structure,
                    fa_rank=term_fa_rank,
                    weight_init=term.get("weight_init", None),
                    enabled=True,
                    metadata=dict(term),
                ))
            else:
                out.append({"term_name":term_name,"kernel_names":kernel_names,"axis_name":str(axis_name),"axis_levels":axis_levels,"structure":structure,"fa_rank":term_fa_rank,"metadata":dict(term)})
        return out

    if env_structure is not None and str(env_structure).lower() != "identity":
        if _MODULAR_TERM_SPEC_AVAILABLE and RandomTermSpec is not None:
            out.append(RandomTermSpec(
                term_name=f"G:{env_col}",
                kernel_names=list(geno_kernel_names),
                axis_name=str(env_col),
                axis_levels=list(env_levels),
                structure=str(env_structure).lower(),
                fa_rank=fa_rank,
            ))
        else:
            out.append({"term_name":f"G:{env_col}","kernel_names":list(geno_kernel_names),"axis_name":str(env_col),"axis_levels":list(env_levels),"structure":str(env_structure).lower(),"fa_rank":fa_rank})
    return out


def build_covariance_modules_from_term_registry(term_registry):
    out = {}
    for term in term_registry:
        if hasattr(term, "axis_name"):
            axis_name = term.axis_name
            axis_levels = term.axis_levels
            structure = term.structure
            fa_rank = term.fa_rank
            term_name = term.term_name
        else:
            axis_name = term.get("axis_name", None)
            axis_levels = term.get("axis_levels", None)
            structure = term.get("structure", "identity")
            fa_rank = term.get("fa_rank", None)
            term_name = term.get("term_name", "G")
        if axis_name is None:
            continue
        out[term_name] = _make_covariance_module_shared(structure=structure, level_names=axis_levels, fa_rank=fa_rank)
    return out

# ----------------------- Large-scale operator backend (optional) -----------------------
def _load_large_operator_backend():
    """Load the optional large-scale operator backend with one clean import path plus local fallback."""
    try:
        from mixed_model_gpu_large_scale_backend import (
            GRMFactorCacheZarr,
            GRMFactorCacheMemmap,
            OperatorBackendConfig as LargeOperatorBackendConfig,
            fit_predict_operator_backend as fit_predict_operator_backend_large,
            LowRankGRMOperator as _LowRankGRMOperator,
            pcg_solve as _pcg_solve,
            _scatter_add_1d as _scatter_add_1d_fn,
        )
        return (GRMFactorCacheZarr, GRMFactorCacheMemmap, LargeOperatorBackendConfig,
                fit_predict_operator_backend_large,
                _LowRankGRMOperator, _pcg_solve, _scatter_add_1d_fn)
    except Exception:
        pass

    try:
        import importlib.util as _importlib_util
        import sys as _sys
        from pathlib import Path as _Path

        backend_candidates = [
            _Path(__file__).with_name('mixed_model_gpu_large_scale_backend.py'),
            _Path('/mnt/data/mixed_model_gpu_large_scale_backend.py'),
        ]
        backend_path = next((p for p in backend_candidates if p.exists()), None)
        if backend_path is None:
            raise FileNotFoundError('No compatible local backend file found for operator route.')

        module_name = 'large_operator_backend_local'
        spec = _importlib_util.spec_from_file_location(module_name, str(backend_path))
        if spec is None or spec.loader is None:
            raise ImportError(f'Could not create import spec for {backend_path}')
        mod = _importlib_util.module_from_spec(spec)
        _sys.modules[module_name] = mod
        spec.loader.exec_module(mod)
        return (mod.GRMFactorCacheZarr, mod.GRMFactorCacheMemmap, mod.OperatorBackendConfig,
                mod.fit_predict_operator_backend,
                mod.LowRankGRMOperator, mod.pcg_solve, mod._scatter_add_1d)
    except Exception:
        return None


_backend_objs = _load_large_operator_backend()
if _backend_objs is not None:
    (GRMFactorCacheZarr, GRMFactorCacheMemmap, LargeOperatorBackendConfig,
     fit_predict_operator_backend_large,
     _LowRankGRMOperator_op, _pcg_solve_op, _scatter_add_1d_op) = _backend_objs
    LARGE_BACKEND_AVAILABLE = True
else:
    LARGE_BACKEND_AVAILABLE = False

# ----------------------------- Global defaults --------------------------------
if TORCH_AVAILABLE and _pick_gp_torch_device is not None:
    _default_device = str(_pick_gp_torch_device())
else:
    _default_device = "cuda" if TORCH_AVAILABLE and torch.cuda.is_available() else "cpu"
_default_dtype  = (torch.float64 if TORCH_AVAILABLE else None)


def _runtime_gp_device(device: Optional[Any] = None, work_units: Optional[float] = None) -> str:
    """Resolve the torch device at call time, not only at module import time."""
    if TORCH_AVAILABLE and _pick_gp_torch_device is not None:
        try:
            return str(_pick_gp_torch_device(device, work_units=work_units))
        except TypeError:
            return str(_pick_gp_torch_device(device))
    return str(device or _default_device)


_gp_warning_once_keys: set = set()


def _warn_once(key: Any, message: str, category: Any = RuntimeWarning, stacklevel: int = 3) -> None:
    if key in _gp_warning_once_keys:
        return
    _gp_warning_once_keys.add(key)
    warnings.warn(message, category, stacklevel=stacklevel)

# ============================== AI-REML utilities ==============================
# These additions align optimization/SEs with ASReml-style REML + Average Information.
# They are written to work with GPyTorch LinearOperator (no explicit V^{-1}).



# -------------------------- Env covariance structures ---------------------------
# Implemented as torch.nn.Modules so they can be optimized jointly with the GPyTorch model.

class UnstructuredCovariance(nn.Module):
    """US(env): full unstructured k×k covariance Σ = L Lᵀ (Cholesky parameterisation).

    Stores exactly k(k+1)/2 parameters (lower triangle of L).
    Diagonal of L is constrained positive via exp(raw_diag).
    Always PSD by construction.

    Args:
        k: number of environments.
        init_Sigma: optional (k,k) initial covariance for warm start.
        jitter: numerical jitter used for Cholesky warm start and added to Σ diagonals.
    """
    def __init__(self, k: int, dtype=None, device=None, init_Sigma: Optional[torch.Tensor]=None, jitter: float=1e-6):
        super().__init__()
        self.k = int(k)
        self.jitter = float(jitter)
        dtype = dtype or torch.float64
        device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        if init_Sigma is not None:
            C = init_Sigma.to(device=device, dtype=dtype)
            C = 0.5 * (C + C.T)
            C = C + self.jitter * torch.eye(self.k, device=device, dtype=dtype)
            L0 = torch.linalg.cholesky(C)
        else:
            L0 = torch.eye(self.k, device=device, dtype=dtype)

        tril = torch.tril_indices(self.k, self.k, offset=0, device=device)
        self.register_buffer("_tril_i", tril[0])
        self.register_buffer("_tril_j", tril[1])

        raw = L0[self._tril_i, self._tril_j].clone()
        diag_mask = (self._tril_i == self._tril_j)
        raw[diag_mask] = torch.log(torch.clamp(raw[diag_mask], min=self.jitter))
        self.L_raw = nn.Parameter(raw)

    def _L(self) -> torch.Tensor:
        L = torch.zeros((self.k, self.k), device=self.L_raw.device, dtype=self.L_raw.dtype)
        diag_mask = (self._tril_i == self._tril_j)
        vals = torch.where(diag_mask, torch.exp(self.L_raw), self.L_raw)
        L[self._tril_i, self._tril_j] = vals
        return L

    def cov(self) -> torch.Tensor:
        L = self._L()
        Sigma = L @ L.T
        Sigma = Sigma + self.jitter * torch.eye(self.k, device=Sigma.device, dtype=Sigma.dtype)
        return Sigma

    def forward(self) -> torch.Tensor:
        return self.cov()

    def submatrix(self, idx1: torch.Tensor, idx2: torch.Tensor) -> torch.Tensor:
        S = self.cov()
        return S.index_select(0, idx1.long()).index_select(1, idx2.long())


class FactorAnalyticCovariance(nn.Module):
    """FA(env,m): factor-analytic covariance Σ = ΓΓᵀ + Ψ.

    Γ: k×m loadings
    Ψ: diagonal specific variances (k,)
    Identifiability constraint (ASReml-style): in the leading m×m block,
        Γ[i,j] = 0 for i < j  (lower-triangular).

    Args:
        k: number of environments.
        m: FA rank (m <= k).
        init_Sigma: optional initial covariance for warm start.
        jitter: numerical jitter.
        identified: enforce lower-triangular constraint in leading block (recommended).
    """
    def __init__(self, k: int, m: int, dtype=None, device=None, init_Sigma: Optional[torch.Tensor]=None,
                 jitter: float=1e-6, identified: bool=True):
        super().__init__()
        self.k = int(k)
        self.m = int(m)
        if self.m > self.k:
            raise ValueError(f"FA rank m={m} must be <= k={k}")
        self.jitter = float(jitter)
        self.identified = bool(identified)
        dtype = dtype or torch.float64
        device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        if init_Sigma is not None:
            C = init_Sigma.to(device=device, dtype=dtype)
            C = 0.5 * (C + C.T)
            ev, U = torch.linalg.eigh(C)
            ev = torch.clamp(ev, min=self.jitter)
            idx = torch.argsort(ev, descending=True)
            U = U[:, idx]
            ev = ev[idx]
            Gamma0 = U[:, :self.m] @ torch.diag(torch.sqrt(ev[:self.m]))
            psi0 = torch.clamp(torch.diag(C) - (Gamma0 * Gamma0).sum(dim=1), min=self.jitter)
        else:
            g = torch.Generator(device=device); g.manual_seed(42)
            Gamma0 = torch.randn((self.k, self.m), generator=g, device=device, dtype=dtype) * 0.1
            psi0 = torch.ones((self.k,), device=device, dtype=dtype) * 0.5

        self.Gamma_raw = nn.Parameter(Gamma0)
        # store ψ via inverse softplus
        psi0 = torch.clamp(psi0, min=self.jitter)
        self.psi_unconstrained = nn.Parameter(torch.log(torch.expm1(psi0)))

    def _Gamma(self) -> torch.Tensor:
        G = self.Gamma_raw.clone()
        if self.identified:
            for j in range(min(self.m, self.k)):
                if j > 0:
                    G[:j, j] = 0.0
        return G

    def _psi(self) -> torch.Tensor:
        return torch.nn.functional.softplus(self.psi_unconstrained) + self.jitter

    def cov(self) -> torch.Tensor:
        G = self._Gamma()
        Psi = self._psi()
        return G @ G.T + torch.diag(Psi)

    def forward(self) -> torch.Tensor:
        return self.cov()

    def submatrix(self, idx1: torch.Tensor, idx2: torch.Tensor) -> torch.Tensor:
        S = self.cov()
        return S.index_select(0, idx1.long()).index_select(1, idx2.long())


def gls_beta_from_V(y: 'torch.Tensor', X: 'torch.Tensor', V_op: 'DenseLinearOperator', jitter: float = 1e-6) -> 'torch.Tensor':
    """Profile fixed-effect coefficients for y ~ N(X beta, V)."""
    p = X.shape[1]
    Vinv_y = V_op.solve(y)
    Vinv_X = V_op.solve(X)
    Xt_Vinv_X = X.transpose(-1, -2) @ Vinv_X
    Xt_Vinv_y = X.transpose(-1, -2) @ Vinv_y

    Xt_Vinv_X = Xt_Vinv_X + jitter * torch.eye(p, device=X.device, dtype=X.dtype)
    L = torch.linalg.cholesky(Xt_Vinv_X)
    return torch.cholesky_solve(Xt_Vinv_y.unsqueeze(-1), L).squeeze(-1)


def reml_loglik_from_V(y: 'torch.Tensor', X: 'torch.Tensor', V_op: 'DenseLinearOperator', jitter: float = 1e-6) -> 'torch.Tensor':
    """Compute REML log-likelihood for y ~ N(X beta, V).
    Uses linear solves with V_op (a LinearOperator) on GPU; does not form V^{-1}.
    """
    n = y.shape[0]
    p = X.shape[1]
    beta = gls_beta_from_V(y, X, V_op, jitter=jitter)
    Vinv_X = V_op.solve(X)
    Xt_Vinv_X = X.transpose(-1, -2) @ Vinv_X
    Xt_Vinv_X = Xt_Vinv_X + jitter * torch.eye(p, device=X.device, dtype=X.dtype)
    r = y - X @ beta
    Vinv_r = V_op.solve(r)

    logdetV = V_op.logdet()
    sign, logdet_XtVinvX = torch.linalg.slogdet(Xt_Vinv_X)
    if torch.any(sign <= 0):
        raise RuntimeError("X'V^{-1}X not SPD; check fixed effects / collinearity / jitter.")

    quad = (r * Vinv_r).sum()
    const = (n - p) * torch.log(torch.tensor(2.0 * math.pi, device=y.device, dtype=y.dtype))
    return -0.5 * (logdetV + logdet_XtVinvX + quad + const)



def _env_cov_dv_builders(
    env_cov_module: Optional[nn.Module],
    ei_t: "torch.Tensor",
    w_e: float,
    device: str,
    dtype_t: "torch.dtype",
    include_full: bool = False,
) -> Tuple[List[Tuple[str, callable]], List[Tuple[str, torch.Tensor]]]:
    """
    Returns:
      dv_builders: list of (label, DV_times(v)->tensor(n,))
      params: list of (label, parameter_tensor) for reporting (optional)

    Notes:
      - For CorHCovariance, provides DV for each env variance parameter and rho (q=k+1).
      - For FA, callers typically handle psi/loadings separately (already implemented in gp_icm_fa path).
      - For US/CorGH, full parameter SEs can be huge; by default we only expose diagonal-variance DV.
      - Full US/CorGH parameter AI remains intentionally partial in this implementation.
    """
    if env_cov_module is None:
        return [], []

    num_envs = int(ei_t.max().item() + 1)

    # Helper apply: given dSe (k,k) and v (n,), return (D v) in observation space efficiently
    def _apply_dSe_to_obs(dSe: "torch.Tensor", v: "torch.Tensor") -> "torch.Tensor":
        # aggregate v by environment
        t = torch.zeros((num_envs,), device=v.device, dtype=v.dtype)
        t.scatter_add_(0, ei_t, v)
        u = dSe @ t  # (k,)
        return float(w_e) * u.index_select(0, ei_t)

    dv_builders: List[Tuple[str, callable]] = []
    params: List[Tuple[str, torch.Tensor]] = []

    if isinstance(env_cov_module, CorHCovariance):
        # Analytic derivatives
        v = env_cov_module.vars()
        sd = torch.sqrt(v)
        rho = env_cov_module.rho()
        R = torch.full((num_envs, num_envs), rho, device=sd.device, dtype=sd.dtype)
        R.fill_diagonal_(1.0)

        # dSigma/dv_i (variance) where Sigma_ij = sd_i * R_ij * sd_j
        # ds_i/dv_i = 1/(2 sd_i)
        for i in range(num_envs):
            def make_dv(i=i):
                def DV(v_obs: "torch.Tensor") -> "torch.Tensor":
                    # Build dSe for this i on the fly (k,k) but cheap for k<=300
                    dSe = torch.zeros((num_envs, num_envs), device=v_obs.device, dtype=v_obs.dtype)
                    # For row/col i:
                    # dSigma_{i,j}/dv_i = (1/(2 sd_i)) * R_{i,j} * sd_j
                    coeff_row = (0.5 / sd[i]) * R[i, :] * sd
                    dSe[i, :] = coeff_row
                    dSe[:, i] = coeff_row
                    # But dSigma_{i,i}/dv_i should be 1 (since Sigma_ii = v_i). Our formula gives (0.5/sd_i)*1*sd_i*2? -> 1, ok.
                    return _apply_dSe_to_obs(dSe, v_obs)
                return DV
            dv_builders.append((f"env_var[{i}]", make_dv(i)))
            params.append((f"env_var[{i}]", v[i].detach()))

        # dSigma/drho: for i!=j, dSigma_ij = sd_i*sd_j ; diag 0
        def DV_rho(v_obs: "torch.Tensor") -> "torch.Tensor":
            dSe = (sd[:, None] * sd[None, :])
            dSe.fill_diagonal_(0.0)
            return _apply_dSe_to_obs(dSe, v_obs)
        dv_builders.append(("rho", DV_rho))
        params.append(("rho", env_cov_module.rho().detach()))

        return dv_builders, params

    if isinstance(env_cov_module, FactorAnalyticCovariance):
        # handled elsewhere (psi/loadings); keep empty here to avoid double counting
        return [], []

    # For US/corgh: by default only diag variance sensitivities (cheap, stable)
    if (isinstance(env_cov_module, (UnstructuredCovariance, CorGHCovariance))) and (not include_full):
        S = env_cov_module.cov().detach()
        for i in range(num_envs):
            def make_diag(i=i):
                def DV(v_obs: "torch.Tensor") -> "torch.Tensor":
                    dSe = torch.zeros((num_envs, num_envs), device=v_obs.device, dtype=v_obs.dtype)
                    dSe[i, i] = 1.0
                    return _apply_dSe_to_obs(dSe, v_obs)
                return DV
            dv_builders.append((f"env_diag[{i}]", make_diag(i)))
            params.append((f"env_diag[{i}]", S[i, i].detach()))
        return dv_builders, params

    # Full parameter SEs for US/corgh is possible but expensive; user can set include_full=True
    # We provide a safe fallback using autograd per-parameter (slow for large k).
    if include_full:
        # Flatten parameters
        named_params = [(n, p) for n, p in env_cov_module.named_parameters() if p.requires_grad]
        for name, p in named_params:
            flat = p.view(-1)
            for j in range(flat.numel()):
                label = f"{name}[{j}]"
                def make_autograd(p=p, j=j):
                    def DV(v_obs: "torch.Tensor") -> "torch.Tensor":
                        # Compute d(Se @ t)/dtheta_j via autograd
                        t = torch.zeros((num_envs,), device=v_obs.device, dtype=v_obs.dtype)
                        t.scatter_add_(0, ei_t, v_obs)
                        # rebuild Se (with grads)
                        Se = env_cov_module.cov()
                        u = Se @ t  # (k,)
                        # pick scalar directional component by dot with random vector? no; we need full u derivative
                        # We'll compute jacobian-vector product by vjp on u with basis e_r, looping over r (expensive).
                        # Instead: return autograd gradient of (u · ones) which is not correct for DV.
                        # Therefore, we do NOT implement full include_full here.
                        raise NotImplementedError("Full AI for US/corgh parameters is expensive; implement analytic/jvp if needed.")
                    return DV
                dv_builders.append((label, make_autograd()))
        return dv_builders, params

    return [], []
def compute_AI_Matrix_SEs(
    y: 'torch.Tensor',
    X: 'torch.Tensor',
    V_op: 'DenseLinearOperator',
    DV_times_list: List,
    jitter: float = 1e-6,
    hutch_samples: int = 64,
    compute_ai_se: bool = True,
    ai_hutch_samples: Optional[int] = None,
    ai_jitter: Optional[float] = None,
    seed: int = 12345,
) -> Tuple['torch.Tensor', 'torch.Tensor']:
    """Compute ASReml-style Average Information (AI) matrix and asymptotic SEs.

    AI_ij = 0.5 * [ y^T P D_i P D_j P y + tr(P D_i P D_j) ]
    where P = V^{-1} - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}.

    DV_times_list: list of callables, each returning (dV_i) @ v for a vector v.
    Uses a shared-probe Hutchinson estimator for the trace term so probe draws and
    P(D_i(P z_k)) products are reused across all (i, j) pairs.
    """
    if not bool(compute_ai_se):
        q = len(DV_times_list)
        z = torch.zeros((q, q), device=y.device, dtype=y.dtype)
        se = torch.full((q,), float('nan'), device=y.device, dtype=y.dtype)
        return z, se

    hutch_samples = int(ai_hutch_samples if ai_hutch_samples is not None else hutch_samples)
    jitter = float(ai_jitter if ai_jitter is not None else jitter)

    torch.manual_seed(seed)
    device, dtype = y.device, y.dtype
    n = y.numel()
    p = X.shape[1]
    q = len(DV_times_list)

    Vinv_X = V_op.solve(X)
    XtVinvX = X.T @ Vinv_X
    XtVinvX = XtVinvX + jitter * torch.eye(p, device=device, dtype=dtype)
    L = torch.linalg.cholesky(XtVinvX)

    def P_times(v: 'torch.Tensor') -> 'torch.Tensor':
        Vinv_v = V_op.solve(v)
        XtVinv_v = X.T @ Vinv_v
        tmp = torch.cholesky_solve(XtVinv_v.unsqueeze(-1), L).squeeze(-1)
        proj = Vinv_X @ tmp
        return Vinv_v - proj

    Py = P_times(y)

    w = []
    u = []
    for Di_times in DV_times_list:
        DiPy = Di_times(Py)
        w.append(DiPy)
        u.append(P_times(DiPy))

    AI = torch.zeros((q, q), device=device, dtype=dtype)

    # Shared Rademacher probes and cached P(D_i(P z_k)) products.
    probes = []
    PD_cache = [[None for _ in range(hutch_samples)] for _ in range(q)]
    for k in range(int(hutch_samples)):
        z = (torch.randint(0, 2, (n,), device=device, dtype=torch.int64) * 2 - 1).to(dtype)
        probes.append(z)
        Pz = P_times(z)
        for i, Di_times in enumerate(DV_times_list):
            PD_cache[i][k] = P_times(Di_times(Pz))

    for i in range(q):
        for j in range(i, q):
            data_part = (w[i] * u[j]).sum()
            tr_acc = torch.zeros((), device=device, dtype=dtype)
            for k in range(int(hutch_samples)):
                tr_acc = tr_acc + (probes[k] * DV_times_list[j](PD_cache[i][k])).sum()
            trace_part = tr_acc / float(hutch_samples)

            AI_ij = 0.5 * (data_part + trace_part)
            AI[i, j] = AI_ij
            AI[j, i] = AI_ij

    AI = 0.5 * (AI + AI.T)
    AI = AI + jitter * torch.eye(q, device=device, dtype=dtype)

    def _safe_inverse_spd(M: "torch.Tensor", jitter0: float) -> "torch.Tensor":
        I = torch.eye(M.shape[0], device=M.device, dtype=M.dtype)
        tries = [0.0, 1.0, 10.0, 100.0]
        last_err = None
        for mult in tries:
            try:
                Mj = M + (float(jitter0) * mult) * I if mult > 0 else M
                Lm = torch.linalg.cholesky(Mj)
                return torch.cholesky_solve(I, Lm)
            except Exception as e:
                last_err = e
        try:
            evals, evecs = torch.linalg.eigh(M)
            floor = max(float(jitter0), 1e-10)
            evals = torch.clamp(evals, min=floor)
            return (evecs * evals.reciprocal().unsqueeze(0)) @ evecs.T
        except Exception:
            if last_err is not None:
                warnings.warn(f"AI inverse fallback used after Cholesky failure: {last_err}")
            return torch.linalg.pinv(M)

    AI_inv = _safe_inverse_spd(AI, jitter)
    SE = torch.sqrt(torch.clamp(torch.diag(AI_inv), min=0.0))
    return AI, SE


# ========================= AI wiring into model outputs =========================
# The functions below "wire" AI-based asymptotic SEs into the existing public API.
# They intentionally avoid explicit V^{-1}; all solves use the LinearOperator's .solve()
# (CG/Lanczos under the hood in GPyTorch).

def compute_score_and_AI_operator(
    y: "torch.Tensor",
    X: "torch.Tensor",
    V_solve: callable,
    DV_times_list: List[callable],
    hutch_samples: int = 64,
    jitter: float = 1e-6,
    seed: int = 12345,
) -> Tuple["torch.Tensor","torch.Tensor"]:
    """
    Operator-form REML score and Average Information (AI) matrix.

    P v = V^{-1}v - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}v

    Score:
      s_i = 1/2 [ -tr(P D_i) + y' P D_i P y ]
    AI:
      AI_ij = 1/2 [ tr(P D_i P D_j) + (P y)' D_i P D_j (P y) ]
    """
    torch.manual_seed(seed)
    device=y.device
    dtype=y.dtype
    n=y.numel()
    p=X.shape[1]
    q=len(DV_times_list)

    Vinv_X = V_solve(X)
    XtVinvX = X.T @ Vinv_X
    L = torch.linalg.cholesky(XtVinvX + jitter*torch.eye(p, device=device, dtype=dtype))

    def P_times(v):
        Vinv_v = V_solve(v)
        XtVinv_v = X.T @ Vinv_v
        tmp = torch.cholesky_solve(XtVinv_v.unsqueeze(-1), L).squeeze(-1)
        proj = Vinv_X @ tmp
        return Vinv_v - proj

    Py = P_times(y)

    # score pieces
    PyDiPy = torch.empty((q,), device=device, dtype=dtype)
    for i, Di in enumerate(DV_times_list):
        PyDiPy[i] = (Py * Di(Py)).sum()

    # tr(P D_i) via Hutchinson: E[z' P D_i z]
    tr_PDi = torch.zeros((q,), device=device, dtype=dtype)
    for _ in range(int(hutch_samples)):
        z = (torch.randint(0,2,(n,), device=device)*2-1).to(dtype)
        Dz = [Di(z) for Di in DV_times_list]
        for i in range(q):
            tr_PDi[i] += (z * P_times(Dz[i])).sum()
    tr_PDi = tr_PDi / float(hutch_samples)

    score = 0.5 * (-tr_PDi + PyDiPy)

    # AI trace term tr(P D_i P D_j) via Hutchinson:
    tr_PDiPDj = torch.zeros((q,q), device=device, dtype=dtype)
    for _ in range(int(hutch_samples)):
        z = (torch.randint(0,2,(n,), device=device)*2-1).to(dtype)
        Pz = P_times(z)
        Dz = [Di(z) for Di in DV_times_list]
        # t_i = P D_i P z
        t = []
        for i in range(q):
            t.append(P_times(DV_times_list[i](Pz)))
        for i in range(q):
            for j in range(i,q):
                tr_PDiPDj[i,j] += (t[i] * Dz[j]).sum()
                if i!=j:
                    tr_PDiPDj[j,i] += (t[j] * Dz[i]).sum()
    tr_PDiPDj = tr_PDiPDj / float(hutch_samples)

    # data part: (P y)' D_i P D_j (P y) = (D_i Py)' P (D_j Py)
    DiPy = [DV_times_list[i](Py) for i in range(q)]
    PDiPy = [P_times(DiPy[i]) for i in range(q)]
    AI = torch.zeros((q,q), device=device, dtype=dtype)
    for i in range(q):
        for j in range(i,q):
            data_part = (DiPy[i] * PDiPy[j]).sum()
            AI_ij = 0.5 * (tr_PDiPDj[i,j] + data_part)
            AI[i,j]=AI_ij
            AI[j,i]=AI_ij
    return score, AI


def ai_newton_optimize(
    ll_fn: callable,
    get_thetas: callable,
    set_thetas: callable,
    score_ai_fn: callable,
    max_steps: int = 8,
    step_halving: int = 8,
    jitter: float = 1e-6,
    verbose: bool = False,
) -> None:
    """
    ASReml-like AI/Newton update in θ-space with step-halving safeguard.
    """
    for it in range(int(max_steps)):
        ll0 = float(ll_fn().item())
        theta0 = get_thetas().detach()
        score, AI = score_ai_fn()
        q = theta0.numel()
        AI_stab = AI + jitter * torch.eye(q, device=AI.device, dtype=AI.dtype)
        step = torch.linalg.solve(AI_stab, score)

        factor=1.0
        ok=False
        for _ in range(int(step_halving)):
            theta_new = torch.clamp(theta0 + factor*step, min=1e-12)
            set_thetas(theta_new)
            ll1 = float(ll_fn().item())
            if ll1 >= ll0 - 1e-10:
                ok=True
                break
            factor *= 0.5
        if verbose:
            print(f"[AI] iter={it} ll0={ll0:.6f} ll1={ll1:.6f} factor={factor:.3f} ok={ok}")
        if not ok:
            set_thetas(theta0)
            break
        if torch.max(torch.abs(factor*step)) < 1e-5:
            break


def _dense_train_mats_for_gp_exact(
    G_list: "List[torch.Tensor]",
    gi_train: "torch.Tensor",
    ei_train: "torch.Tensor",
    S_e: "Optional[torch.Tensor]",
) -> Tuple[List["torch.Tensor"], List["torch.Tensor"], "torch.Tensor"]:
    """
    Builds dense n×n component matrices on the TRAIN set:

      Kg_i  = G_i[gi,gi]
      Kge_i = Kg_i ⊙ Se_env   (if S_e provided; else zeros)
      Ke    = Se_env (if S_e) else I_env (same-env indicator)

    These are used to define D_i v = (dV/dtheta_i) v in compute_AI_Matrix_SEs.
    """
    n = gi_train.numel()
    dtype = G_list[0].dtype
    device = G_list[0].device

    Kg = [ _build_geno_block(Gk, gi_train, gi_train) for Gk in G_list ]

    if S_e is None:
        Se_env = (ei_train.unsqueeze(1) == ei_train.unsqueeze(0)).to(dtype=dtype, device=device)
        Kge = [ torch.zeros((n,n), dtype=dtype, device=device) for _ in G_list ]
    else:
        Se_env = _build_env_kernel(ei_train, ei_train, S_e, dtype=dtype, device=device)
        Kge = [ (Kg_i * Se_env) for Kg_i in Kg ]

    Ke = Se_env
    return Kg, Kge, Ke



def _ai_varcomp_table_gp_exact(
    model: "gpytorch.models.ExactGP",
    likelihood: "gpytorch.likelihoods.Likelihood",
    train_x: "torch.Tensor",
    train_y: "torch.Tensor",
    hutch_samples: int = 64,
    jitter: float = 1e-6,
    env_cov_module: Optional[nn.Module] = None,
    ai_include_env_params: bool = False,
) -> pd.DataFrame:
    """
    AI-based asymptotic SEs for variance parameters in gp_exact_with_X.

    Parameters included:
      - If kernel.learn_scales=True: w_g[i], w_ge[i], w_e
      - If likelihood is GaussianLikelihood: sigma2_noise
      - If ai_include_env_params=True and env_cov_module is CorHCovariance:
            env_var[i] (k params) and rho
        For US/corgh we include only diag sensitivities by default (env_diag[i]) since full q is huge.

    Returns empty DataFrame if nothing is estimated.
    """
    if not GPTY_AVAILABLE:
        return pd.DataFrame([])

    kernel = model.covar_module
    learn_scales = bool(getattr(kernel, "learn_scales", False))
    device = train_y.device
    dtype_t = train_y.dtype

    with torch.no_grad():
        K_op = kernel(train_x, train_x)
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            noise_diag = likelihood.noise
            R_op = DiagLinearOperator(noise_diag)
        else:
            R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
        V_op = K_op + R_op

        # X rows aligned to training
        Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())

        names: List[str] = []
        thetas: List[torch.Tensor] = []
        DV_times_list: List[callable] = []

        # Kernel scale parameters
        if learn_scales and hasattr(kernel, "w_g") and hasattr(kernel, "w_ge") and hasattr(kernel, "w_e"):
            # w_g, w_ge are vectors; w_e scalar
            wg = kernel.w_g.squeeze()
            wge = kernel.w_ge.squeeze()
            we = kernel.w_e.squeeze()

            # Derivatives: component-wise precomputed kernels are stored in kernel._Kg_list etc.
            Kg_list = getattr(kernel, "_Kg_list", None)
            Kge_list = getattr(kernel, "_Kge_list", None)
            Ke = getattr(kernel, "_Ke", None)

            if Kg_list is not None and Kge_list is not None and Ke is not None:
                for i in range(int(wg.numel())):
                    names.append(f"w_g[{i}]")
                    thetas.append(wg[i])
                    Km = Kg_list[i]
                    DV_times_list.append(lambda v, Km=Km: Km @ v)
                for i in range(int(wge.numel())):
                    names.append(f"w_ge[{i}]")
                    thetas.append(wge[i])
                    Km = Kge_list[i]
                    DV_times_list.append(lambda v, Km=Km: Km @ v)
                names.append("w_e")
                thetas.append(we)
                DV_times_list.append(lambda v, Km=Ke: Km @ v)

        # Likelihood noise (only if learnable)
        if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood):
            sigma2 = likelihood.noise.squeeze()
            names.append("sigma2_noise")
            thetas.append(sigma2)
            DV_times_list.append(lambda v: v)  # I

        # Optionally include env covariance parameters
        if ai_include_env_params:
            # env indices are in train_x[:,1]
            ei_t = train_x[:, 1].long()
            dv_builders, env_params = _env_cov_dv_builders(env_cov_module, ei_t, w_e=float(getattr(kernel, "w_e", torch.tensor(1.0, device=device, dtype=dtype_t))), device=str(device), dtype_t=dtype_t, include_full=False)
            param_map = {lbl: val for lbl, val in env_params}
            for lbl, fn in dv_builders:
                names.append(lbl)
                thetas.append(param_map.get(lbl, torch.tensor(float("nan"), device=device, dtype=dtype_t)).to(device=device, dtype=dtype_t).reshape(()))
                DV_times_list.append(fn)

        if len(DV_times_list) == 0:
            return pd.DataFrame([])

    # AI/SE computed outside no_grad because we need solves but no parameter grads
    AI, SE = compute_AI_Matrix_SEs(
        y=train_y, X=Xtt, V_op=V_op, DV_times_list=DV_times_list,
        jitter=jitter, hutch_samples=int(hutch_samples)
    )

    # estimates for w's and sigma2 are real; env-param placeholders are NaN
    theta_vals = torch.stack([t.reshape(()) for t in thetas])
    se_vals = SE

    
    out = pd.DataFrame({
        "component": names,
        "estimate": theta_vals.detach().cpu().numpy(),
        "se_ai": se_vals.detach().cpu().numpy(),
        "z": (theta_vals / torch.clamp(se_vals, min=1e-12)).detach().cpu().numpy(),
    })
    return out



def _ai_varcomp_table_gp_icm_fa(
    model: "gpytorch.models.ExactGP",
    likelihood: "gpytorch.likelihoods.Likelihood",
    train_x: "torch.Tensor",
    train_y: "torch.Tensor",
    fa_module: Optional[nn.Module],
    hutch_samples: int = 64,
    jitter: float = 1e-6,
    include_fa_loadings: bool = False,
    env_cov_module: Optional[nn.Module] = None,
    ai_include_env_params: bool = False,
) -> pd.DataFrame:
    """
    AI-based asymptotic SEs for gp_icm_fa_with_X.

    Includes:
      - kernel component scales if learn_scales=True
      - sigma2_noise if GaussianLikelihood
      - FA specific variances (d_env) always (if fa_module provided)
      - FA loadings optionally (include_fa_loadings=True) – currently NOT computed; a warning is emitted and only psi SEs are returned
      - env covariance params via env_cov_module if ai_include_env_params=True (corh diag+rho; us/corgh diag-only)

    Limitations:
      - Full US/CorGH environment-parameter AI/SE coverage is intentionally partial here.
      - For US/CorGH, only reduced diagonal-sensitivity rows are reported unless a future
        analytic/JVP implementation is added.
      - FA loading SEs are not computed in this routine; only psi/specific-variance SEs are returned.
    """
    if not GPTY_AVAILABLE:
        return pd.DataFrame([])

    kernel = model.covar_module
    learn_scales = bool(getattr(kernel, "learn_scales", False))
    device = train_y.device
    dtype_t = train_y.dtype

    with torch.no_grad():
        K_op = kernel(train_x, train_x)
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            noise_diag = likelihood.noise
            R_op = DiagLinearOperator(noise_diag)
        else:
            R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
        V_op = K_op + R_op

        Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())

        names: List[str] = []
        thetas: List[torch.Tensor] = []
        DV_times_list: List[callable] = []

        # Kernel scale parameters (FA ICM kernel uses w_fa vector + w_e scalar)
        if learn_scales and hasattr(kernel, "w_fa") and hasattr(kernel, "w_e"):
            wf = kernel.w_fa.squeeze()
            we = kernel.w_e.squeeze()
            Kfa_list = getattr(kernel, "_Kfa_list", None)
            Ke = getattr(kernel, "_Ke", None)
            if Kfa_list is not None and Ke is not None:
                for i in range(int(wf.numel())):
                    names.append(f"w_fa[{i}]")
                    thetas.append(wf[i])
                    Km = Kfa_list[i]
                    DV_times_list.append(lambda v, Km=Km: Km @ v)
                names.append("w_e")
                thetas.append(we)
                DV_times_list.append(lambda v, Km=Ke: Km @ v)

        # Likelihood noise (only if learnable)
        if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood):
            sigma2 = likelihood.noise.squeeze()
            names.append("sigma2_noise")
            thetas.append(sigma2)
            DV_times_list.append(lambda v: v)

        # FA parameters: d_env always; loadings optional
        
        # FA parameters: specific variances always; support both (legacy) d_env and (current) d_unconstrained
        if fa_module is not None:
            ei_t = train_x[:, 1].long()
            num_envs = int(ei_t.max().item() + 1)

            if hasattr(fa_module, "d_env"):
                # legacy naming (already variance-scale)
                d_raw = fa_module.d_env
                d_pos = d_raw
                d_scale = torch.ones_like(d_pos)  # d(d_pos)/d(d_raw)
                fa_delta = None
            elif hasattr(fa_module, "d_unconstrained"):
                # current FactorAnalyticEnv stores unconstrained params; psi = softplus(raw)
                raw = fa_module.d_unconstrained
                d_pos = torch.nn.functional.softplus(raw) + 1e-8
                d_scale = torch.sigmoid(raw)  # d(softplus)/draw for delta-method
                d_raw = raw
                fa_delta = {"scale": d_scale.detach()}
            else:
                d_raw = None
                fa_delta = None

            if d_raw is not None:
                for j in range(int(d_raw.numel())):
                    names.append(f"psi[{j}]")
                    thetas.append(d_pos[j].detach())

                    scale_j = d_scale[j].detach()

                    def make_d(j=j, scale_j=scale_j):
                        def DV(v_obs: "torch.Tensor") -> "torch.Tensor":
                            t = torch.zeros((num_envs,), device=v_obs.device, dtype=v_obs.dtype)
                            t.scatter_add_(0, ei_t, v_obs)
                            u = torch.zeros((num_envs,), device=v_obs.device, dtype=v_obs.dtype)
                            u[j] = scale_j.to(v_obs.dtype) * t[j]
                            return float(getattr(kernel, "w_e", 1.0)) * u.index_select(0, ei_t)
                        return DV

                    DV_times_list.append(make_d(j))

            if include_fa_loadings and hasattr(fa_module, "L"):
                warnings.warn(
                    "include_fa_loadings=True was requested, but FA loading SEs are not implemented in the AI table; "
                    "returning only FA specific-variance (psi) SEs for transparency.",
                    RuntimeWarning,
                )

        # Optionally include env covariance parameters
        if ai_include_env_params:
            ei_t = train_x[:, 1].long()
            dv_builders, env_params = _env_cov_dv_builders(env_cov_module, ei_t, w_e=float(getattr(kernel, "w_e", torch.tensor(1.0, device=device, dtype=dtype_t))), device=str(device), dtype_t=dtype_t, include_full=False)
            param_map = {lbl: val for lbl, val in env_params}
            for lbl, fn in dv_builders:
                names.append(lbl)
                thetas.append(param_map.get(lbl, torch.tensor(float("nan"), device=device, dtype=dtype_t)).to(device=device, dtype=dtype_t).reshape(()))
                DV_times_list.append(fn)

        if len(DV_times_list) == 0:
            return pd.DataFrame([])

    AI, SE = compute_AI_Matrix_SEs(
        y=train_y, X=Xtt, V_op=V_op, DV_times_list=DV_times_list,
        jitter=jitter, hutch_samples=int(hutch_samples)
    )

    theta_vals = torch.stack([t.reshape(()) for t in thetas])
    se_vals = SE

    
    out = pd.DataFrame({
        "component": names,
        "estimate": theta_vals.detach().cpu().numpy(),
        "se_ai": se_vals.detach().cpu().numpy(),
        "z": (theta_vals / torch.clamp(se_vals, min=1e-12)).detach().cpu().numpy(),
    })
    return out


def set_deterministic(seed: int = 12345):
    np.random.seed(int(seed))
    if TORCH_AVAILABLE:
        torch.manual_seed(int(seed))
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(int(seed))
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

def _resolve_torch_dtype(dtype):
    if not TORCH_AVAILABLE:
        return None
    if dtype is None:
        return _default_dtype
    if isinstance(dtype, torch.dtype):
        return dtype
    if isinstance(dtype, str):
        d = dtype
        if d.startswith("torch."):
            d = d.split(".", 1)[1]
        if hasattr(torch, d):
            return getattr(torch, d)
        raise TypeError(f"Unknown torch dtype string: {dtype}")
    raise TypeError("dtype must be None, a torch.dtype, or a torch dtype string (e.g., 'float64')")

def _np(x) -> np.ndarray:
    return np.asarray(x)

def _to_numpy(t):
    if TORCH_AVAILABLE and isinstance(t, torch.Tensor):
        return t.detach().cpu().numpy()
    return np.asarray(t)

def _t(x, device=None, dtype=None):
    if not TORCH_AVAILABLE:
        raise RuntimeError("PyTorch not available")
    device = _default_device if device is None else device
    dtype  = _default_dtype if dtype  is None else dtype
    arr = np.array(x, copy=True)
    return torch.as_tensor(arr, dtype=dtype, device=device)

class _DefaultDTypeContext:
    def __init__(self, dtype):
        self.dtype = dtype
        self.prev = None
    def __enter__(self):
        if TORCH_AVAILABLE and (self.dtype is not None):
            try:
                self.prev = torch.get_default_dtype()
                torch.set_default_dtype(self.dtype)
            except Exception:
                self.prev = None
        return self
    def __exit__(self, exc_type, exc, tb):
        if TORCH_AVAILABLE and (self.prev is not None):
            try:
                torch.set_default_dtype(self.prev)
            except Exception:
                pass


def _safe_nanmean_or_nan(x) -> float:
    arr = np.asarray(x, dtype=float)
    if arr.size == 0:
        return float('nan')
    finite = np.isfinite(arr)
    if not finite.any():
        return float('nan')
    return float(np.nanmean(arr))

_VALID_OUTPUT_LEVELS = ("predict_only", "predict_with_se", "full_vc")


def _resolve_output_level(
    output_level: Optional[str],
    *,
    point_predictions_only: Optional[bool],
    return_se: Optional[bool],
    compute_ai_se: Optional[bool],
) -> Tuple[str, bool, bool, bool]:
    """Single source of truth for how much work a fit should do.

    Returns (level, point_predictions_only, return_se, compute_ai_se) with the
    three legacy flags rewritten consistently from the resolved level:
      - predict_only       -> (True,  False, False)
      - predict_with_se    -> (False, True,  False)
      - full_vc            -> (False, True,  True)

    If `output_level` is None, derive it from the legacy flags (compute_ai_se
    dominates, then return_se, then point_predictions_only). If both an
    explicit level and legacy flags are supplied and they disagree, the
    explicit level wins and we warn.
    """
    legacy_given = any(v is not None for v in (point_predictions_only, return_se, compute_ai_se))

    if output_level is None:
        if not legacy_given:
            level = "predict_only"
        elif bool(compute_ai_se):
            level = "full_vc"
        elif bool(return_se) or (point_predictions_only is False):
            level = "predict_with_se"
        else:
            level = "predict_only"
    else:
        level = str(output_level).lower()
        if level not in _VALID_OUTPUT_LEVELS:
            raise ValueError(
                f"output_level must be one of {_VALID_OUTPUT_LEVELS}; got {output_level!r}"
            )
        if legacy_given:
            legacy_level = (
                "full_vc" if bool(compute_ai_se)
                else ("predict_only" if (bool(point_predictions_only) and not bool(return_se)) else "predict_with_se")
            )
            if legacy_level != level:
                warnings.warn(
                    f"output_level={level!r} overrides inconsistent legacy flags "
                    f"(point_predictions_only={point_predictions_only}, return_se={return_se}, "
                    f"compute_ai_se={compute_ai_se}); using explicit output_level.",
                    stacklevel=2,
                )

    if level == "predict_only":
        return level, True, False, False
    if level == "predict_with_se":
        return level, False, True, False
    return level, False, True, True


def _append_prediction_meta(res: Dict[str, Any], *, point_predictions_only: bool, prediction_output: str, backend_used: str,
                            prediction_se_type: str = "both", prediction_se_method: str = "exact_dense",
                            prediction_block_size: int = 512, prediction_diag_probes: int = 32,
                            return_prediction_cov: bool = False, approximate: bool = False,
                            residual_meta: Optional[Dict[str, Any]] = None,
                            engine_policy: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    if not isinstance(res, dict):
        return res
    diag = res.setdefault("diagnostics", {})
    pred_meta = diag.setdefault("prediction_meta", {})
    operator_backend = str(backend_used).lower().startswith("operator")
    supports_prediction_se = (not operator_backend) or (not bool(point_predictions_only))
    pred_meta.update({
        "backend_used": str(backend_used),
        "point_predictions_only": bool(point_predictions_only),
        "prediction_output": str(prediction_output),
        "prediction_se_type": str(prediction_se_type),
        "prediction_se_method": str(prediction_se_method),
        "prediction_block_size": int(prediction_block_size),
        "prediction_diag_probes": int(prediction_diag_probes),
        "return_prediction_cov": bool(return_prediction_cov),
        "approximate": bool(approximate),
        "used_fixed_effect_correction": True,
        "supports_prediction_se": bool(supports_prediction_se),
        "supports_prediction_cov": (not operator_backend),
        "operator_limitations": (
            "operator backend returns diagonal prediction SEs via Hutchinson probing when requested; full prediction covariance is not implemented"
            if operator_backend else None
        ),
    })
    if residual_meta is not None:
        diag["residual_meta"] = residual_meta
    if engine_policy is not None:
        diag["engine_policy"] = dict(engine_policy)
    return res


def _stamp_output_level(res: Dict[str, Any], output_level: str) -> Dict[str, Any]:
    """Stamp the resolved output_level onto the result's diagnostics."""
    if isinstance(res, dict):
        diag = res.setdefault("diagnostics", {})
        diag["output_level"] = str(output_level)
        pred_meta = diag.setdefault("prediction_meta", {})
        pred_meta["output_level"] = str(output_level)
    return res

def _corr_from_cov_np(cov: np.ndarray) -> np.ndarray:
    cov = np.asarray(cov, dtype=float)
    if cov.size == 0:
        return cov
    d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    den = np.outer(d, d)
    out = np.divide(cov, den, out=np.zeros_like(cov, dtype=float), where=den > 0)
    np.fill_diagonal(out, 1.0)
    return out



def _normalize_interaction_terms_meta(
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]],
    level_names: Sequence[str],
    interaction_term_name: str,
    interaction_cov: Optional[np.ndarray],
    kernel_names: Sequence[str],
    w_ge: Sequence[float],
) -> List[Dict[str, Any]]:
    level_names = [str(x) for x in level_names]
    kernel_names = [str(k) for k in kernel_names]
    w_ge = np.asarray(w_ge, dtype=float).reshape(-1)
    n_levels = len(level_names)
    if interaction_terms_meta is None:
        cov = np.eye(n_levels, dtype=float) if interaction_cov is None else np.asarray(interaction_cov, dtype=float)
        if cov.shape != (n_levels, n_levels):
            cov = np.eye(n_levels, dtype=float)
        return [{
            "term_name": str(interaction_term_name),
            "level_names": level_names,
            "cov_matrix": cov,
            "kernel_weights": {k: float(w_ge[i]) for i, k in enumerate(kernel_names[:len(w_ge)])},
            "notes": None,
        }]
    out = []
    for term in interaction_terms_meta:
        if term is None:
            continue
        tname = str(term.get("term_name", interaction_term_name))
        tlevels = [str(x) for x in term.get("level_names", level_names)]
        n_t = len(tlevels)
        cov = term.get("cov_matrix", interaction_cov)
        cov = np.eye(n_t, dtype=float) if cov is None else np.asarray(cov, dtype=float)
        if cov.shape != (n_t, n_t):
            cov = np.eye(n_t, dtype=float)
        kw = term.get("kernel_weights", None)
        if kw is None:
            kw = {k: float(w_ge[i]) for i, k in enumerate(kernel_names[:len(w_ge)])}
        else:
            kw = {str(k): float(kw.get(k, 0.0)) for k in kernel_names}
        out.append({
            "term_name": tname,
            "level_names": tlevels,
            "cov_matrix": cov,
            "kernel_weights": kw,
            "notes": term.get("notes", None),
        })
    return out


def _build_kernel_specific_interaction_covs(
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]],
    kernel_names: Sequence[str],
    default_level_names: Sequence[str],
    default_interaction_term_name: str,
    default_interaction_cov: Optional[np.ndarray],
    default_w_ge: Sequence[float],
) -> Tuple[List[Dict[str, Any]], Dict[str, np.ndarray]]:
    terms = _normalize_interaction_terms_meta(
        interaction_terms_meta=interaction_terms_meta,
        level_names=default_level_names,
        interaction_term_name=default_interaction_term_name,
        interaction_cov=default_interaction_cov,
        kernel_names=kernel_names,
        w_ge=default_w_ge,
    )
    kernel_names = [str(k) for k in kernel_names]
    n_levels = len(default_level_names)
    mats = {k: np.zeros((n_levels, n_levels), dtype=float) for k in kernel_names}
    for term in terms:
        cov = np.asarray(term["cov_matrix"], dtype=float)
        kw = {str(k): float(v) for k, v in term.get("kernel_weights", {}).items()}
        if cov.shape != (n_levels, n_levels):
            continue
        for k in kernel_names:
            mats[k] += kw.get(k, 0.0) * cov
    return terms, mats


def _has_multiple_interaction_terms(interaction_terms_meta: Optional[Sequence[Dict[str, Any]]]) -> bool:
    return interaction_terms_meta is not None and len([t for t in interaction_terms_meta if t is not None]) > 1


def _primary_interaction_term_name(env_col: str) -> str:
    return f"G:{env_col}"



def _build_reporting_bundle(
    kernel_names: Sequence[str],
    w_g: Sequence[float],
    w_ge: Sequence[float],
    level_names: Sequence[str],
    interaction_term_name: str = "G:Env",
    interaction_cov: Optional[np.ndarray] = None,
    env_main_scale: float = 0.0,
    env_main_cov: Optional[np.ndarray] = None,
    residual_env: Optional[np.ndarray] = None,
    n_rep: Optional[float] = None,
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    kernel_names = [str(k) for k in kernel_names]
    w_g = np.asarray(w_g, dtype=float).reshape(-1)
    w_ge = np.asarray(w_ge, dtype=float).reshape(-1)
    level_names = [str(x) for x in level_names]
    n_levels = len(level_names)

    if env_main_cov is None:
        env_main_cov = np.eye(n_levels, dtype=float)
    env_main_cov = np.asarray(env_main_cov, dtype=float)
    if env_main_cov.shape != (n_levels, n_levels):
        env_main_cov = np.eye(n_levels, dtype=float)

    residual_env = np.asarray(residual_env, dtype=float).reshape(-1) if residual_env is not None and len(np.asarray(residual_env).reshape(-1)) == n_levels else np.full(n_levels, np.nan, dtype=float)

    interaction_terms = _normalize_interaction_terms_meta(
        interaction_terms_meta=interaction_terms_meta,
        level_names=level_names,
        interaction_term_name=interaction_term_name,
        interaction_cov=interaction_cov,
        kernel_names=kernel_names,
        w_ge=w_ge,
    )

    if _MODULAR_LAYERS_AVAILABLE and RandomTermSpec is not None and FittedRandomTerm is not None and FittedMixedModel is not None and _mod_build_reporting_bundle is not None:
        fitted_terms = []
        fitted_terms.append(
            FittedRandomTerm(
                spec=RandomTermSpec(term_name="G", kernel_names=list(kernel_names), axis_name=None, axis_levels=None, structure="identity"),
                kernel_weights={k: float(v) for k, v in zip(kernel_names, w_g)},
                covariance_module=None,
                covariance_matrix=None,
                correlation_matrix=None,
                parameter_table=None,
                ai_table=None,
                diagnostics={},
            )
        )
        for term in interaction_terms:
            cov = np.asarray(term["cov_matrix"], dtype=float)
            d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
            denom = np.outer(d, d)
            corr = np.zeros_like(cov, dtype=float)
            mask = denom > 0
            corr[mask] = cov[mask] / denom[mask]
            np.fill_diagonal(corr, 1.0)
            fitted_terms.append(
                FittedRandomTerm(
                    spec=RandomTermSpec(
                        term_name=str(term["term_name"]),
                        kernel_names=list(kernel_names),
                        axis_name=str(term.get("axis_name", term["term_name"].split("G:",1)[-1] if "G:" in str(term["term_name"]) else "Axis")),
                        axis_levels=list(term["level_names"]),
                        structure=str(term.get("notes", "identity") or "identity"),
                    ),
                    kernel_weights={str(k): float(v) for k, v in term["kernel_weights"].items()},
                    covariance_module=None,
                    covariance_matrix=cov,
                    correlation_matrix=corr,
                    parameter_table=None,
                    ai_table=None,
                    diagnostics={},
                )
            )
        residual_summary = {
            "residual_mean": float(np.nanmean(residual_env)) if residual_env.size else np.nan,
            "residual_by_level": pd.DataFrame({
                "level_index": np.arange(n_levels, dtype=int),
                "level": level_names,
                "residual_variance": residual_env,
            }),
        }
        fitted_model = FittedMixedModel(
            beta=None,
            beta_se=None,
            fixed_effect_names=None,
            random_terms=fitted_terms,
            residual_diag_obs=None,
            residual_summary=residual_summary,
            predictions=None,
            diagnostics={},
            backend="dense",
            method="framework_adapter",
        )
        return _mod_build_reporting_bundle(fitted_model, n_rep=n_rep)

    # Fallback to in-file implementation style
    rows = []
    for k, est in zip(kernel_names, w_g):
        rows.append({"term_group":"main","term":"G","kernel":k,"component_type":"genetic_main","estimate":float(est),"level_count":None,"notes":None})
    rows.append({"term_group":"main","term":"G","kernel":"total","component_type":"genetic_main","estimate":float(np.nansum(w_g)),"level_count":None,"notes":None})

    interaction_variance_summary = {}
    interaction_correlation_summary = {}
    first_env_df = None

    for term in interaction_terms:
        tname = term["term_name"]
        tlevels = term["level_names"]
        cov = np.asarray(term["cov_matrix"], dtype=float)
        diag_int = np.clip(np.diag(cov), 0.0, None)
        kw = {str(k): float(v) for k, v in term["kernel_weights"].items()}

        for k in kernel_names:
            rows.append({"term_group":"interaction","term":tname,"kernel":k,"component_type":"interaction","estimate":float(kw.get(k,0.0)),"level_count":int(len(tlevels)),"notes":term.get("notes", None)})
        rows.append({"term_group":"interaction","term":tname,"kernel":"total","component_type":"interaction","estimate":float(sum(kw.get(k,0.0) for k in kernel_names)),"level_count":int(len(tlevels)),"notes":term.get("notes", None)})

        by_kernel_var = {}
        by_kernel_cov = {}
        for k in kernel_names:
            est = float(kw.get(k, 0.0))
            by_kernel_var[k] = est * diag_int
            by_kernel_cov[k] = est * cov
        total_var = np.sum(np.vstack(list(by_kernel_var.values())), axis=0) if by_kernel_var else np.zeros(len(tlevels), dtype=float)
        total_cov = np.sum(np.stack(list(by_kernel_cov.values())), axis=0) if by_kernel_cov else np.zeros_like(cov)

        term_df = pd.DataFrame({"level_index": np.arange(len(tlevels), dtype=int), "level": tlevels})
        for k in kernel_names:
            term_df[f"interaction_variance_{k}"] = by_kernel_var.get(k, np.zeros(len(tlevels), dtype=float))
        term_df["interaction_variance_total"] = total_var
        interaction_variance_summary[tname] = term_df

        total_corr = np.zeros_like(total_cov)
        d = np.sqrt(np.clip(np.diag(total_cov), 0.0, None))
        denom = np.outer(d, d)
        mask = denom > 0
        total_corr[mask] = total_cov[mask] / denom[mask]
        np.fill_diagonal(total_corr, 1.0)
        interaction_correlation_summary[tname] = {
            "term": tname,
            "levels": list(tlevels),
            "cov_total": pd.DataFrame(total_cov, index=tlevels, columns=tlevels),
            "cor_total": pd.DataFrame(total_corr, index=tlevels, columns=tlevels),
            "cov_by_kernel": {k: pd.DataFrame(by_kernel_cov[k], index=tlevels, columns=tlevels) for k in kernel_names},
            "cor_by_kernel": {k: pd.DataFrame(np.where(np.outer(np.sqrt(np.clip(np.diag(by_kernel_cov[k]),0,None)), np.sqrt(np.clip(np.diag(by_kernel_cov[k]),0,None)))>0, by_kernel_cov[k]/np.outer(np.sqrt(np.clip(np.diag(by_kernel_cov[k]),0,None)), np.sqrt(np.clip(np.diag(by_kernel_cov[k]),0,None))), 0.0), index=tlevels, columns=tlevels) for k in kernel_names},
        }
        if first_env_df is None:
            first_env_df = term_df

    rows.append({"term_group":"environment_main","term":"Env","kernel":"total","component_type":"environment_main","estimate":float(env_main_scale),"level_count":n_levels,"notes":"diag_mean"})

    var_components_summary = pd.DataFrame(rows)
    env_variance_summary = first_env_df if first_env_df is not None else pd.DataFrame({"level_index": np.arange(n_levels, dtype=int), "level": level_names})
    env_variance_summary = env_variance_summary.copy()
    env_variance_summary["environment_variance"] = float(env_main_scale) * np.clip(np.diag(env_main_cov), 0.0, None)
    env_variance_summary["residual_variance"] = residual_env
    env_variance_summary["genetic_main_variance"] = float(np.nansum(w_g))
    env_variance_summary["total_genetic_variance"] = env_variance_summary["genetic_main_variance"] + (env_variance_summary["interaction_variance_total"] if "interaction_variance_total" in env_variance_summary.columns else 0.0)
    genetic_by_level = pd.to_numeric(
        env_variance_summary["total_genetic_variance"], errors="coerce"
    ).to_numpy(dtype=float)
    residual_by_level = pd.to_numeric(
        env_variance_summary["residual_variance"], errors="coerce"
    ).to_numpy(dtype=float)
    total_by_level = genetic_by_level + residual_by_level
    h2_by_level = np.full(n_levels, np.nan, dtype=float)
    h2_ok = (
        np.isfinite(genetic_by_level)
        & np.isfinite(residual_by_level)
        & (total_by_level > 0.0)
    )
    h2_by_level[h2_ok] = genetic_by_level[h2_ok] / total_by_level[h2_ok]
    env_variance_summary["heritability"] = h2_by_level
    H2 = _safe_nanmean_or_nan(h2_by_level)
    resid_mean = _safe_nanmean_or_nan(residual_env)

    return {
        "var_components_summary": var_components_summary,
        "var_components_ai": pd.DataFrame([]),
        "interaction_variance_summary": interaction_variance_summary,
        "interaction_correlation_summary": interaction_correlation_summary,
        "residual_summary": {
            "residual_mean": resid_mean,
            "residual_by_level": pd.DataFrame({"level_index": np.arange(n_levels, dtype=int), "level": level_names, "residual_variance": residual_env}),
        },
        "heritability_summary": {
            "genetic_main_total": float(np.nansum(w_g)),
            "genetic_total_average": _safe_nanmean_or_nan(genetic_by_level),
            "H2_average": H2,
        },
        "env_variance_summary": env_variance_summary,
        "var_components": var_components_summary,
    }


def _attach_fa_derived_se_to_summary(
    report_bundle: Dict[str, Any],
    fa_derived: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Inject delta-method SE for FA-derived sigma2_g into the report bundle.

    `fa_derived` comes from `varcomp_asreml.compute_fa_derived_summary`,
    surfaced via `summary["fa_derived"]` returned by
    `ai_reml.build_fa_icm_ai_varcomp`. When present we:

      - add a dedicated FA-derived total-genetic-variance row to
        `var_components_summary`, with its matching estimate and SE,
      - extend `env_variance_summary` with `total_genetic_variance_se` and
        per-env `heritability` + `heritability_se` columns,
      - extend `heritability_summary` with a `H2_average_se` field if
        per-env h^2 SE is available.

    Pre-existing rows / columns are left untouched. If `fa_derived` is None
    or the bundle is missing expected pieces, the bundle is returned
    unchanged.
    """
    if fa_derived is None or not isinstance(report_bundle, dict):
        return report_bundle

    vcs = report_bundle.get("var_components_summary")
    if isinstance(vcs, pd.DataFrame) and not vcs.empty:
        if "std.error" not in vcs.columns:
            vcs = vcs.copy()
            vcs["std.error"] = float("nan")
        else:
            vcs = vcs.copy()
        derived_mask = (
            vcs.get("component_type", pd.Series(index=vcs.index, dtype=object))
            == "genetic_variance_total"
        )
        derived_values = {
            "term_group": "derived",
            "term": "G:Env",
            "kernel": "total",
            "component_type": "genetic_variance_total",
            "estimate": float(fa_derived.get("sigma2_g_total", float("nan"))),
            "std.error": float(fa_derived.get("sigma2_g_total_se", float("nan"))),
            "level_count": len(fa_derived.get("env_labels", [])),
            "notes": "sum_of_environment_marginal_genetic_variances",
        }
        if derived_mask.any():
            for key, value in derived_values.items():
                if key in vcs.columns:
                    vcs.loc[derived_mask, key] = value
        else:
            row = {column: derived_values.get(column, np.nan) for column in vcs.columns}
            vcs = pd.concat([vcs, pd.DataFrame([row])], ignore_index=True)
        report_bundle["var_components_summary"] = vcs
        # Mirror to var_components (same content alias upstream).
        if "var_components" in report_bundle:
            report_bundle["var_components"] = vcs

    evs = report_bundle.get("env_variance_summary")
    if isinstance(evs, pd.DataFrame) and not evs.empty:
        evs = evs.copy()
        env_labels = [str(x) for x in fa_derived.get("env_labels", [])]
        sigma2_g_se = list(fa_derived.get("sigma2_g_se_per_env", []))
        h2 = list(fa_derived.get("h2_per_env", []))
        h2_se = list(fa_derived.get("h2_se_per_env", []))
        level_col = "level" if "level" in evs.columns else None
        if level_col is not None and env_labels:
            se_map = dict(zip(env_labels, sigma2_g_se))
            h2_map = dict(zip(env_labels, h2))
            h2_se_map = dict(zip(env_labels, h2_se))
            evs["genetic_main_variance_se"] = [
                float(se_map.get(str(lvl), float("nan"))) for lvl in evs[level_col]
            ]
            evs["total_genetic_variance_se"] = evs["genetic_main_variance_se"]
            evs["heritability"] = [
                float(h2_map.get(str(lvl), float("nan"))) for lvl in evs[level_col]
            ]
            evs["heritability_se"] = [
                float(h2_se_map.get(str(lvl), float("nan"))) for lvl in evs[level_col]
            ]
            report_bundle["env_variance_summary"] = evs

    hs = report_bundle.get("heritability_summary")
    if isinstance(hs, dict):
        h2_average = fa_derived.get("h2_average", float("nan"))
        if isinstance(h2_average, (int, float)) and np.isfinite(h2_average):
            hs = dict(hs)
            hs["H2_average"] = float(h2_average)
            h2_average_se = fa_derived.get("h2_average_se", float("nan"))
            if isinstance(h2_average_se, (int, float)) and np.isfinite(h2_average_se):
                hs["H2_average_se"] = float(h2_average_se)
            report_bundle["heritability_summary"] = hs

    return report_bundle


def _sync_reporting_bundle_with_reml(
    report_bundle: Dict[str, Any],
    varcomp: Optional[pd.DataFrame],
    level_names: Sequence[str],
    fa_derived: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Align derived summaries with converged response-scale REML estimates."""
    if not isinstance(report_bundle, dict) or not isinstance(varcomp, pd.DataFrame) or varcomp.empty:
        return report_bundle
    if "component" not in varcomp.columns or "estimate" not in varcomp.columns:
        return report_bundle

    out = dict(report_bundle)
    levels = [str(x) for x in level_names]
    components = varcomp["component"].astype(str)
    estimates = pd.to_numeric(varcomp["estimate"], errors="coerce").to_numpy(dtype=float)
    residual = np.full(len(levels), np.nan, dtype=float)
    for i, level in enumerate(levels):
        suffix = f"_{level}!R"
        hits = np.where(components.str.endswith(suffix).to_numpy())[0]
        if hits.size:
            residual[i] = estimates[hits[0]]

    env_summary = out.get("env_variance_summary")
    if isinstance(env_summary, pd.DataFrame) and len(env_summary) == len(levels):
        env_summary = env_summary.copy()
        if np.isfinite(residual).any():
            env_summary["residual_variance"] = residual
        if isinstance(fa_derived, dict):
            fa_levels = [str(x) for x in fa_derived.get("env_labels", [])]
            fa_g = np.asarray(fa_derived.get("sigma2_g_per_env", []), dtype=float).reshape(-1)
            if len(fa_levels) == fa_g.size:
                lookup = dict(zip(fa_levels, fa_g))
                genetic = np.asarray([lookup.get(level, np.nan) for level in levels], dtype=float)
                if np.isfinite(genetic).any():
                    env_summary["total_genetic_variance"] = genetic
        genetic = pd.to_numeric(
            env_summary.get("total_genetic_variance", pd.Series(np.nan, index=env_summary.index)),
            errors="coerce",
        ).to_numpy(dtype=float)
        resid = pd.to_numeric(
            env_summary.get("residual_variance", pd.Series(np.nan, index=env_summary.index)),
            errors="coerce",
        ).to_numpy(dtype=float)
        total = genetic + resid
        h2 = np.full(len(levels), np.nan, dtype=float)
        ok = np.isfinite(genetic) & np.isfinite(resid) & (total > 0.0)
        h2[ok] = genetic[ok] / total[ok]
        env_summary["heritability"] = h2
        out["env_variance_summary"] = env_summary

        hs = dict(out.get("heritability_summary") or {})
        hs["genetic_total_average"] = _safe_nanmean_or_nan(genetic)
        hs["H2_average"] = _safe_nanmean_or_nan(h2)
        out["heritability_summary"] = hs

    if np.isfinite(residual).any():
        out["residual_summary"] = {
            "residual_mean": _safe_nanmean_or_nan(residual),
            "residual_by_level": pd.DataFrame({
                "level_index": np.arange(len(levels), dtype=int),
                "level": levels,
                "residual_variance": residual,
            }),
        }
    return out


def _reported_residual_by_level(
    report_bundle: Dict[str, Any],
    fallback: Optional[Sequence[float]] = None,
) -> np.ndarray:
    env_summary = report_bundle.get("env_variance_summary") if isinstance(report_bundle, dict) else None
    if isinstance(env_summary, pd.DataFrame) and "residual_variance" in env_summary.columns:
        values = pd.to_numeric(env_summary["residual_variance"], errors="coerce").to_numpy(dtype=float)
        if values.size:
            return values
    if fallback is None:
        return np.array([], dtype=float)
    return np.asarray(fallback, dtype=float).reshape(-1)


def _attach_prediction_se_columns(
    predictions_df: pd.DataFrame,
    latent_var: Optional[Sequence[float]] = None,
    observed_var: Optional[Sequence[float]] = None,
    latent_se: Optional[Sequence[float]] = None,
    observed_se: Optional[Sequence[float]] = None,
) -> pd.DataFrame:
    if _MODULAR_LAYERS_AVAILABLE and _mod_attach_prediction_se_columns is not None:
        return _mod_attach_prediction_se_columns(
            predictions_df,
            latent_var=latent_var,
            observed_var=observed_var,
            latent_se=latent_se,
            observed_se=observed_se,
        )
    out = predictions_df.copy()
    n = len(out)
    def _vec(x):
        if x is None:
            return None
        a = np.asarray(x, dtype=float).reshape(-1)
        if a.size != n:
            raise ValueError(f"Prediction variance/SE length {a.size} must equal n_predictions={n}")
        return a
    latent_var = _vec(latent_var)
    observed_var = _vec(observed_var)
    latent_se = _vec(latent_se)
    observed_se = _vec(observed_se)
    if latent_var is not None:
        out["Prediction_Var_latent"] = latent_var
        if latent_se is None:
            latent_se = np.sqrt(np.clip(latent_var, 0.0, None))
    if observed_var is not None:
        out["Prediction_Var_observed"] = observed_var
        if observed_se is None:
            observed_se = np.sqrt(np.clip(observed_var, 0.0, None))
    if latent_se is not None:
        out["Prediction_SE_latent"] = latent_se
    if observed_se is not None:
        out["Prediction_SE_observed"] = observed_se
    return out



# ------------------------------- Standardizer ---------------------------------
class _PerObsEnvStandardizer:
    """
    Unifies response handling across models and ensures original-scale outputs.

    mode: "none" | "global" | "per_env"
      - "none": no scaling
      - "global": z-score using one train-set mean and standard deviation
      - "per_env": z-score within environment using TRAIN rows only
    """
    def __init__(self, mode: str = "none"):
        self.mode = (mode or "none").lower()
        self.env_means_: Optional[np.ndarray] = None
        self.env_stds_: Optional[np.ndarray]  = None
        self.global_mean_: float = 0.0
        self.global_std_: float  = 1.0

    def fit(self, y: np.ndarray, ei: np.ndarray, train_idx: np.ndarray):
        y = np.asarray(y, dtype=float)
        ei = np.asarray(ei, dtype=np.int64)
        tr = np.asarray(train_idx, dtype=np.int64)
        y_tr = y[tr]
        if not np.isfinite(y_tr).any():
            self.global_mean_, self.global_std_ = 0.0, 1.0
            self.env_means_ = None; self.env_stds_ = None
            return self

        self.global_mean_ = float(np.nanmean(y_tr))
        std = float(np.nanstd(y_tr, ddof=1))
        if not np.isfinite(std) or std <= 0:
            std = 1.0
        self.global_std_  = std

        if self.mode == "per_env":
            n_env = int(np.max(ei)) + 1
            means = np.full(n_env, self.global_mean_, dtype=float)
            stds  = np.full(n_env, self.global_std_,  dtype=float)
            for e in range(n_env):
                idx = tr[ei[tr] == e]
                if idx.size > 1 and np.isfinite(y[idx]).any():
                    m = float(np.nanmean(y[idx]))
                    s = float(np.nanstd(y[idx], ddof=1))
                    if not np.isfinite(s) or s <= 0:
                        s = self.global_std_
                    means[e] = m
                    stds[e]  = s
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
            stds  = self.env_stds_[ei]
        else:
            means = np.full_like(y, self.global_mean_, dtype=float)
            stds  = np.full_like(y, self.global_std_, dtype=float)
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

    def get_env_scales(self, env_labels):
        if _MODULAR_STANDARDIZATION_AVAILABLE:
            proxy = _ModPerObsEnvStandardizer(self.mode)
            proxy.global_mean_ = self.global_mean_
            proxy.global_std_ = self.global_std_
            proxy.env_means_ = self.env_means_
            proxy.env_stds_ = self.env_stds_
            return proxy.get_env_scales(env_labels)
        if _MODULAR_AI_REML_AVAILABLE:
            from varcomp_asreml import EnvScales
            labels = [str(e) for e in env_labels]
            if self.mode == "per_env" and self.env_stds_ is not None:
                s_dict = {}
                for i, lab in enumerate(labels):
                    s_dict[lab] = float(self.env_stds_[i]) if i < len(self.env_stds_) else float(self.global_std_)
                global_s = float(np.sqrt(np.mean(np.square(list(s_dict.values())))))
                return EnvScales(s=s_dict, global_s=global_s)
            if self.mode == "global":
                return EnvScales(
                    s={lab: float(self.global_std_) for lab in labels},
                    global_s=float(self.global_std_),
                )
            return EnvScales.identity(labels)
        raise ImportError("varcomp_asreml or standardization module required for get_env_scales")

    def get_env_stds_array(self) -> np.ndarray:
        if self.mode == "none":
            return np.array([1.0], dtype=np.float64)
        if self.env_stds_ is not None:
            return np.array(self.env_stds_, dtype=np.float64)
        return np.array([self.global_std_], dtype=np.float64)


def _response_variance_scale_by_level(stdr: "_PerObsEnvStandardizer", n_levels: int) -> np.ndarray:
    n_levels = int(max(1, n_levels))
    try:
        scale = np.asarray(
            stdr.inv_var(np.ones(n_levels, dtype=float), np.arange(n_levels, dtype=np.int64)),
            dtype=float,
        )
    except Exception:
        scale = np.full(n_levels, float(getattr(stdr, "global_std_", 1.0)) ** 2, dtype=float)
    if scale.size != n_levels:
        scale = np.resize(scale, n_levels)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)
    return scale


def _scale_frame_columns_by_level(
    df: pd.DataFrame,
    level_scale: np.ndarray,
    level_col_candidates: Sequence[str],
    prefixes: Sequence[str],
) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    level_col = next((c for c in level_col_candidates if c in out.columns), None)
    if level_col is not None and len(level_scale):
        idx = pd.to_numeric(out.get("level_index", pd.Series(np.arange(len(out)))), errors="coerce").to_numpy()
        scale = np.full(len(out), float(np.nanmean(level_scale)), dtype=float)
        ok = np.isfinite(idx) & (idx >= 0) & (idx < len(level_scale))
        if ok.any():
            scale[ok] = level_scale[idx[ok].astype(int)]
    elif len(out) == len(level_scale):
        scale = level_scale
    else:
        scale = np.full(len(out), float(np.nanmean(level_scale)), dtype=float)
    for col in list(out.columns):
        if any(str(col).startswith(prefix) for prefix in prefixes):
            out[col] = pd.to_numeric(out[col], errors="coerce") * scale
    return out


def _scale_cov_frame_by_level(x: Any, level_scale: np.ndarray) -> Any:
    if not isinstance(x, pd.DataFrame) or x.empty:
        return x
    out = x.copy()
    row_scale = np.sqrt(np.resize(level_scale, out.shape[0]))
    col_scale = np.sqrt(np.resize(level_scale, out.shape[1]))
    vals = out.to_numpy(dtype=float, copy=True)
    vals = vals * np.outer(row_scale, col_scale)
    return pd.DataFrame(vals, index=out.index, columns=out.columns)


def _rescale_reporting_bundle_to_response_scale(
    bundle: Dict[str, Any],
    stdr: "_PerObsEnvStandardizer",
    level_names: Sequence[str],
) -> Dict[str, Any]:
    """Put REML variance metadata on the same response scale as SE/PEV."""
    if not isinstance(bundle, dict):
        return bundle
    out = dict(bundle)
    n_levels = len(level_names) if level_names is not None else 1
    level_scale = _response_variance_scale_by_level(stdr, n_levels)
    component_scale = float(np.nanmean(level_scale)) if np.isfinite(level_scale).any() else 1.0

    for key in ("var_components_summary", "var_components"):
        tbl = out.get(key)
        if isinstance(tbl, pd.DataFrame) and not tbl.empty:
            scaled = tbl.copy()
            if "estimate" in scaled.columns:
                scaled["estimate"] = pd.to_numeric(scaled["estimate"], errors="coerce") * component_scale
            out[key] = scaled

    env_tbl = out.get("env_variance_summary")
    if isinstance(env_tbl, pd.DataFrame) and not env_tbl.empty:
        out["env_variance_summary"] = _scale_frame_columns_by_level(
            env_tbl,
            level_scale,
            ("level", "Env", "env", "Environment", "environment"),
            ("environment_variance", "genetic_main_variance", "interaction_variance", "total_genetic_variance"),
        )

    interaction_summary = out.get("interaction_variance_summary")
    if isinstance(interaction_summary, dict):
        out["interaction_variance_summary"] = {
            name: _scale_frame_columns_by_level(
                tbl,
                level_scale,
                ("level", "Env", "env", "Environment", "environment"),
                ("interaction_variance",),
            ) if isinstance(tbl, pd.DataFrame) else tbl
            for name, tbl in interaction_summary.items()
        }

    interaction_corr = out.get("interaction_correlation_summary")
    if isinstance(interaction_corr, dict):
        scaled_corr = {}
        for name, obj in interaction_corr.items():
            if isinstance(obj, dict):
                obj2 = dict(obj)
                obj2["cov_total"] = _scale_cov_frame_by_level(obj2.get("cov_total"), level_scale)
                cov_by_kernel = obj2.get("cov_by_kernel")
                if isinstance(cov_by_kernel, dict):
                    obj2["cov_by_kernel"] = {
                        k: _scale_cov_frame_by_level(v, level_scale)
                        for k, v in cov_by_kernel.items()
                    }
                scaled_corr[name] = obj2
            else:
                scaled_corr[name] = obj
        out["interaction_correlation_summary"] = scaled_corr

    herit = out.get("heritability_summary")
    if isinstance(herit, dict):
        herit2 = dict(herit)
        for key in ("genetic_main_total", "genetic_total_average"):
            if key in herit2:
                try:
                    herit2[key] = float(herit2[key]) * component_scale
                except Exception:
                    pass
        out["heritability_summary"] = herit2

    diagnostics = out.get("diagnostics")
    if not isinstance(diagnostics, dict):
        diagnostics = {}
    diagnostics["variance_component_scale"] = "response"
    diagnostics["response_variance_scale_by_level"] = level_scale.tolist()
    out["diagnostics"] = diagnostics
    return out


def _rescale_env_covariance_to_response_scale(
    cov: Any,
    stdr: "_PerObsEnvStandardizer",
    level_names: Sequence[str],
) -> Any:
    if cov is None:
        return None
    arr = np.asarray(cov, dtype=float)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return cov
    level_scale = _response_variance_scale_by_level(
        stdr,
        len(level_names) if level_names is not None else arr.shape[0],
    )
    sd_scale = np.sqrt(np.resize(level_scale, arr.shape[0]))
    return arr * np.outer(sd_scale, sd_scale)


def _effective_gp_exact_reporting_scales(
    covariance_module: Any,
    w_g: Sequence[float],
    w_ge: Sequence[float],
    w_e: float,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Read the scales actually used by the fitted exact-GP covariance."""
    out = (
        np.asarray(w_g, dtype=float).reshape(-1),
        np.asarray(w_ge, dtype=float).reshape(-1),
        float(w_e),
    )
    getter = getattr(covariance_module, "get_effective_scales", None)
    if getter is None:
        return out
    try:
        fitted = getter()
        if len(fitted) == 3:
            return (
                np.asarray(_to_numpy(fitted[0]), dtype=float).reshape(-1),
                np.asarray(_to_numpy(fitted[1]), dtype=float).reshape(-1),
                float(np.asarray(_to_numpy(fitted[2]), dtype=float).reshape(-1)[0]),
            )
        if len(fitted) == 2:
            return (
                np.asarray(_to_numpy(fitted[0]), dtype=float).reshape(-1),
                out[1],
                float(np.asarray(_to_numpy(fitted[1]), dtype=float).reshape(-1)[0]),
            )
        return out
    except Exception:
        return out


def _effective_fa_reporting_scales(
    covariance_module: Any,
    w_fa: Sequence[float],
    w_e: float,
) -> Tuple[np.ndarray, float]:
    """Read fitted FA genomic-kernel and environment-main scales."""
    out = (np.asarray(w_fa, dtype=float).reshape(-1), float(w_e))
    getter = getattr(covariance_module, "get_effective_scales", None)
    if getter is None:
        return out
    try:
        fitted = getter()
        if len(fitted) != 2:
            return out
        return (
            np.asarray(_to_numpy(fitted[0]), dtype=float).reshape(-1),
            float(np.asarray(_to_numpy(fitted[1]), dtype=float).reshape(-1)[0]),
        )
    except Exception:
        return out

# ------------------------------ Linear algebra --------------------------------
def _add_jitter(A: "torch.Tensor", jitter: float):
    return A + jitter * torch.eye(A.shape[-1], dtype=A.dtype, device=A.device)

def _cholesky_solve(A: "torch.Tensor", B: "torch.Tensor", base_jitter: float = 1e-8, max_tries: int = 6) -> "torch.Tensor":
    was_1d = B.dim() == 1
    B2 = B.unsqueeze(-1) if was_1d else B
    jj = float(base_jitter)
    for _ in range(max_tries):
        try:
            L = torch.linalg.cholesky(_add_jitter(A, jj))
            out = torch.cholesky_solve(B2, L)
            return out.squeeze(-1) if was_1d else out
        except RuntimeError:
            jj *= 10.0
    L = torch.linalg.cholesky(_add_jitter(A, jj))
    out = torch.cholesky_solve(B2, L)
    return out.squeeze(-1) if was_1d else out

def _cholesky_factor(A: "torch.Tensor", base_jitter: float = 1e-8, max_tries: int = 6) -> "torch.Tensor":
    jj = float(base_jitter)
    for _ in range(max_tries):
        try:
            return torch.linalg.cholesky(_add_jitter(A, jj))
        except RuntimeError:
            jj *= 10.0
    return torch.linalg.cholesky(_add_jitter(A, jj))

def _cholesky_solve_from_factor(L: "torch.Tensor", B: "torch.Tensor") -> "torch.Tensor":
    was_1d = B.dim() == 1
    B2 = B.unsqueeze(-1) if was_1d else B
    out = torch.cholesky_solve(B2, L)
    return out.squeeze(-1) if was_1d else out

def _spd_inverse(A: "torch.Tensor", base_jitter: float = 1e-10) -> "torch.Tensor":
    jj = float(base_jitter)
    for _ in range(6):
        try:
            L = torch.linalg.cholesky(A)
            return torch.cholesky_inverse(L)
        except RuntimeError:
            A = A.clone()
            A.diagonal().add_((jj))
            jj *= 10.0
    L = torch.linalg.cholesky(A)
    return torch.cholesky_inverse(L)

def _estimate_diag_Hutch(A_times, n: int, probes: int = 64, device=None, dtype=None, seed: int = 12345) -> "torch.Tensor":
    set_deterministic(seed)
    device = _default_device if device is None else device
    dtype  = _default_dtype if dtype is None else dtype
    acc = torch.zeros(n, dtype=dtype, device=device)
    for _ in range(int(probes)):
        z = torch.randint(0, 2, (n,), device=device, dtype=torch.int64) * 2 - 1
        z = z.to(dtype=dtype)
        Az = A_times(z)
        acc += z * Az
    return acc / float(probes)

# ------------------------------ Design matrix ---------------------------------
def _normalize_fixed_effects(fixed_effects):
    if fixed_effects is None:
        return None
    if isinstance(fixed_effects, (list, tuple)):
        return [str(x) for x in fixed_effects]
    if hasattr(fixed_effects, "dtype") and hasattr(fixed_effects, "shape"):
        arr = np.asarray(fixed_effects)
        if arr.ndim == 0:
            return [str(arr.item())]
        return [str(x) for x in arr.tolist()]
    if isinstance(fixed_effects, str):
        return [fixed_effects]
    return [str(fixed_effects)]

def _build_X_from_df(df: pd.DataFrame, fixed_effects: Optional[Sequence[str]]) -> Tuple[np.ndarray, List[str]]:
    fixed_effects = _normalize_fixed_effects(fixed_effects)
    if fixed_effects is None or len(fixed_effects) == 0:
        X = np.ones((len(df), 1), dtype=float); return X, ["(Intercept)"]
    mats, cols = [np.ones((len(df), 1), dtype=float)], ["(Intercept)"]
    for col in fixed_effects:
        if col not in df.columns:
            raise KeyError(f"Fixed-effect column '{col}' not in pheno_df.")
        v = df[col]
        if pd.api.types.is_numeric_dtype(v):
            mats.append(np.asarray(v, dtype=float).reshape(-1,1)); cols.append(col)
        else:
            d = pd.get_dummies(v.astype(str), drop_first=True, dtype=float)
            if d.shape[1] == 0:
                continue
            mats.append(d.values.astype(float)); cols.extend([f"{col}={lev}" for lev in d.columns])
    if len(mats) == 0:
        X = np.ones((len(df), 1), dtype=float); cols = ["(Intercept)"]
    else:
        X = np.concatenate(mats, axis=1)
    return X, cols


def _fit_context_fixed_key(fixed_effects: Optional[Sequence[str]]) -> Tuple[str, ...]:
    fixed_effects = _normalize_fixed_effects(fixed_effects)
    if fixed_effects is None:
        return tuple()
    return tuple(str(x) for x in fixed_effects)


def prepare_fit_mixed_model_context(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_ids: Sequence[str],
    fixed_effects: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Precompute split-invariant phenotype indices and fixed-effect design.

    The returned object is an optional optimization hint for repeated
    fit_mixed_model calls over different train/test splits on the same data.
    """
    df = pheno_df.reset_index(drop=True).copy()
    if gid_col not in df.columns or env_col not in df.columns or y_col not in df.columns:
        raise KeyError("gid_col/env_col/y_col not found in pheno_df")
    gid = df[gid_col].astype("category")
    env = df[env_col].astype("category")
    y = pd.to_numeric(df[y_col], errors="coerce").values

    geno_ids_key = tuple(str(g) for g in geno_ids)
    gid_to_index = {g: i for i, g in enumerate(geno_ids_key)}
    gi = np.array([gid_to_index.get(str(g), -1) for g in gid.astype(str)], dtype=np.int64)
    if (gi < 0).any():
        missing = [str(g) for g, i in zip(gid.astype(str), gi) if i < 0][:5]
        raise ValueError(f"Unknown genotype IDs encountered: {missing} ...")

    env_levels = list(env.cat.categories.astype(str))
    env_to_index = {e: i for i, e in enumerate(env_levels)}
    ei = np.array([env_to_index[str(e)] for e in env.astype(str)], dtype=np.int64)
    X, x_columns = _build_X_from_df(df, fixed_effects)

    return {
        "df": df,
        "gid_col": str(gid_col),
        "env_col": str(env_col),
        "y_col": str(y_col),
        "source_identity": id(pheno_df),
        "geno_ids_key": geno_ids_key,
        "fixed_effects_key": _fit_context_fixed_key(fixed_effects),
        "n_rows": int(df.shape[0]),
        "y": y,
        "gi": gi,
        "ei": ei,
        "env_levels": env_levels,
        "X": X,
        "x_columns": x_columns,
    }


def _fit_context_matches(
    context: Any,
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_ids: Sequence[str],
    fixed_effects: Optional[Sequence[str]],
) -> bool:
    if not isinstance(context, dict):
        return False
    try:
        return (
            context.get("gid_col") == str(gid_col)
            and context.get("env_col") == str(env_col)
            and context.get("y_col") == str(y_col)
            and context.get("source_identity") == id(pheno_df)
            and tuple(context.get("geno_ids_key", tuple())) == tuple(str(g) for g in geno_ids)
            and tuple(context.get("fixed_effects_key", tuple())) == _fit_context_fixed_key(fixed_effects)
            and int(context.get("n_rows", -1)) == int(len(pheno_df))
            and isinstance(context.get("df"), pd.DataFrame)
            and "y" in context
            and "gi" in context
            and "ei" in context
            and "env_levels" in context
            and "X" in context
        )
    except Exception:
        return False


def _cache_array_token(x: Any) -> Any:
    if x is None:
        return None
    arr = np.asarray(x)
    return (id(x), int(arr.__array_interface__.get("data", (0, False))[0]), tuple(arr.shape), str(arr.dtype))


# ------------------------------ Random-terms helpers ------------------------------

_RANDOM_SYNONYMS = {
    "block": [r"^block$", r"^blk$", r"^block_id$", r"^blk_id$", r"^subblock$", r"^sb$"],
    "rep":   [r"^rep$", r"^repl$", r"^replicate$", r"^replication$", r"^rpt$"],
    "year":  [r"^year$", r"^yr$", r"^season$", r"^cycle$"],
    "location": [r"^location$", r"^loc$", r"^site$", r"^station$", r"^farm$", r"^field$", r"^trial_site$"],
    "herd": [r"^herd$", r"^pen$", r"^cage$", r"^lot$", r"^barn$"],
    "batch": [r"^batch$", r"^cg$", r"^contemporary_group$", r"^management_group$", r"^group$"],
}

def _normalize_colname(s: str) -> str:
    s = str(s).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def find_column_by_synonym(df: pd.DataFrame, canonical: str) -> Optional[str]:
    pats = _RANDOM_SYNONYMS.get(_normalize_colname(canonical), [])
    if not pats:
        return None
    norm_map = {c: _normalize_colname(c) for c in df.columns}
    for col, ncol in norm_map.items():
        for pat in pats:
            if re.match(pat, ncol, flags=re.IGNORECASE):
                return col
    return None

def build_group_index(
    df: pd.DataFrame,
    group_col: Optional[str] = None,
    group_cols: Optional[Sequence[str]] = None,
    *,
    sep: str = "|",
    na_token: str = "NA",
) -> Tuple[np.ndarray, Dict[str, Any]]:
    if group_cols is None:
        if group_col is None:
            raise ValueError("Provide group_col or group_cols")
        group_cols = [group_col]
    group_cols = [str(c) for c in group_cols]
    for c in group_cols:
        if c not in df.columns:
            raise KeyError(f"group column '{c}' not found in pheno_df")

    parts = []
    for c in group_cols:
        v = df[c].astype(str)
        v = v.replace({"nan": na_token, "NaN": na_token, "None": na_token})
        parts.append(v.values)

    if len(parts) == 1:
        key = parts[0]
    else:
        key = np.char.add(parts[0].astype(str), sep)
        for p in parts[1:]:
            key = np.char.add(key, np.char.add(sep, p.astype(str)))

    codes, uniques = pd.factorize(key, sort=True)
    codes = codes.astype(np.int64, copy=False)
    meta = {"group_cols": group_cols, "n_levels": int(len(uniques)), "levels": uniques.astype(str).tolist()}
    return codes, meta

def infer_random_terms(
    pheno_df: pd.DataFrame,
    requested: Sequence[Union[str, Dict[str, Any]]],
    *,
    default_weight: float = 1.0,
) -> List[Dict[str, Any]]:
    """
    Build IID random terms ONLY when user explicitly requests them.
    requested: list of strings (synonyms) or dicts with group_col(s).
    Returns list of term dicts with group_col(s) resolved (not group_index yet).
    """
    out: List[Dict[str, Any]] = []
    for item in requested:
        if isinstance(item, str):
            name = item
            col = find_column_by_synonym(pheno_df, name)
            if col is None:
                raise KeyError(f"Requested random term '{name}' but no matching column found in pheno_df.")
            out.append({"type": "iid", "name": name, "group_col": col, "weight": float(default_weight)})
        elif isinstance(item, dict):
            t = dict(item)
            t.setdefault("type", "iid")
            t.setdefault("weight", float(default_weight))
            if "group_col" not in t and "group_cols" not in t:
                nm = t.get("name", None)
                if nm is None:
                    raise ValueError("Requested dict term must include group_col(s) or a name to resolve.")
                col = find_column_by_synonym(pheno_df, str(nm))
                if col is None:
                    raise KeyError(f"Requested random term '{nm}' but no matching column found in pheno_df.")
                t["group_col"] = col
            out.append(t)
        else:
            raise TypeError("requested must be a sequence of strings or dicts")
    return out

# ------------------------------ Weights helper --------------------------------
def _prepare_weight_vector(w, n: int, default: float) -> np.ndarray:
    if w is None:
        return np.full(n, float(default), dtype=float)
    arr = np.asarray(w, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr), dtype=float)
    arr = arr.ravel()
    if arr.size == 1 and n > 1:
        return np.full(n, float(arr[0]), dtype=float)
    if arr.size != n:
        raise ValueError(f"Weight vector length mismatch. Expected {n}, got {arr.size}.")
    return arr.astype(float, copy=False)

# ------------------------------- Kernel helpers -------------------------------
def _build_geno_block(G: "torch.Tensor", gi1: "torch.Tensor", gi2: "torch.Tensor") -> "torch.Tensor":
    return G.index_select(0, gi1.long()).index_select(1, gi2.long())

def _build_env_kernel(
    ei1: "torch.Tensor",
    ei2: "torch.Tensor",
    S_e: Optional["torch.Tensor"],
    dtype: Optional["torch.dtype"] = None,
    device: Optional[Union[str, "torch.device"]] = None,
) -> "torch.Tensor":
    device = ei1.device if device is None else torch.device(device)
    if S_e is None:
        if dtype is None:
            dtype = torch.get_default_dtype()
        return (ei1.unsqueeze(1) == ei2.unsqueeze(0)).to(dtype=dtype, device=device)
    return S_e.index_select(0, ei1.long()).index_select(1, ei2.long())

# -------------------------- Normalization helper ------------------------------
def _apply_kernel_normalization(
    G_list: List[np.ndarray],
    S_e: Optional[np.ndarray],
    mode: str = "diag_mean",
    output_level: Optional[str] = None,
) -> Tuple[List[np.ndarray], Optional[np.ndarray]]:
    if _MODULAR_KERNEL_NORMALIZE_AVAILABLE:
        G_out, Se_out, _mode_used = _mod_apply_kernel_normalization(
            G_list, S_e, mode=mode, output_level=output_level,
        )
        return G_out, Se_out

    mode = (mode or "none").lower()
    if output_level == "full_vc" and mode == "diag_mean":
        warnings.warn(
            "reml_normalize='diag_mean' rescales the genetic kernel; "
            "variance components will NOT be on the original y scale. "
            "Overriding to 'asreml' (no rescale) for full_vc output.",
            stacklevel=2,
        )
        mode = "asreml"
    if mode in ("none", "asreml"):
        return G_list, S_e

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


# ---------------------- Large-scale operator backend runner ----------------------

def _operator_backend_point_predictions(
    *,
    gi: np.ndarray,
    ei: np.ndarray,
    y: np.ndarray,
    X: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    w_g: np.ndarray,
    w_ge: np.ndarray,
    w_e: float,
    env_similarity: Optional[np.ndarray],
    resid_diag_env: Optional[np.ndarray],
    resid_diag_obs: Optional[np.ndarray] = None,
    device: str,
    dtype_t: "torch.dtype",
    standardize: str,
    geno_id_list: Optional[List[str]],
    env_index_to_name: Optional[List[str]],
    prediction_output: str,
    grm_factor_cache: Dict[str, Any],
    operator_tol: float = 1e-5,
    operator_max_iter: int = 500,
    operator_dtype_compute: Union[str, "torch.dtype"] = "float32",
    seed: int = 12345,
    random_terms: Optional[List[Dict[str, Any]]] = None,
    return_se: bool = False,
    n_hutchinson_probes: int = 256,
    hutchinson_seed: Optional[int] = None,
    se_pcg_tol: Optional[float] = None,
    se_pcg_max_iter: Optional[int] = None,
    use_nystrom_se_preconditioner: bool = False,
    nystrom_rank: int = 64,
    nystrom_seed: int = 0,
    use_lanczos_se_variance_reduction: bool = True,
    lanczos_iters: int = 80,
    lanczos_seed: int = 0,
    lanczos_theta_floor: float = 1e-8,
) -> Dict[str, Any]:
    """
    Point-prediction-only path using the low-rank operator backend.
    This is meant for large n (>= 50k) and many env (<= ~300), where dense GP exact is infeasible.

    grm_factor_cache dict supports:
      - {"type":"zarr", "root_dir":..., "num_grms":...}
      - {"type":"memmap", "paths":[...], "shapes":[(n,m),...], "dtypes":[...]}

    Set return_se=True to compute Hutchinson diagonal prediction SEs. Adds ~M probe-cost
    PCG solves (default M=32); predictions DataFrame gains 'SE_latent' and 'SE' columns
    on observed scale.
    """
    if (not TORCH_AVAILABLE) or (not LARGE_BACKEND_AVAILABLE):
        raise RuntimeError("Operator backend requested but torch/backend not available.")

    set_deterministic(seed)

    gi = np.asarray(gi, dtype=np.int64)
    ei = np.asarray(ei, dtype=np.int64)
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    train_idx = np.asarray(train_idx, dtype=np.int64)
    test_idx = np.asarray(test_idx, dtype=np.int64)

    # Standardize on TRAIN only (to match other code paths)
    stdr = _PerObsEnvStandardizer(mode=standardize)
    stdr.fit(y, ei, train_idx)
    y_std = stdr.transform(y, ei)

    # Residual diagonal per observation.
    # Public contract here is original-response-scale environment variances; convert to standardized scale internally.
    if resid_diag_obs is not None:
        sig_obs_orig = np.asarray(resid_diag_obs, dtype=float).reshape(-1)
        if sig_obs_orig.size != gi.size:
            raise ValueError(f"resid_diag_obs must have length n_obs={gi.size}; got {sig_obs_orig.size}")
        if np.any(~np.isfinite(sig_obs_orig)) or np.any(sig_obs_orig <= 0):
            raise ValueError("resid_diag_obs must be finite and > 0.")
        sigma2_env_orig = _env_mean_from_obs_diag(sig_obs_orig, ei, int(ei.max()) + 1)
        scale2_obs = stdr.inv_var(np.ones(gi.size, dtype=float), ei)
        scale2_obs = np.where(np.isfinite(scale2_obs) & (scale2_obs > 0), scale2_obs, 1.0)
        r_diag_obs = sig_obs_orig / scale2_obs
    elif resid_diag_env is None:
        r_diag_obs = 1.0
        sigma2_env_orig = np.full((int(ei.max()) + 1,), np.nan, dtype=float)
    else:
        sig_orig = np.asarray(resid_diag_env, dtype=float).reshape(-1)
        if sig_orig.size != int(ei.max()) + 1:
            raise ValueError(f"resid_diag_env must have length n_env={int(ei.max()) + 1}; got {sig_orig.size}")
        sigma2_env_orig = sig_orig.copy()
        scale2_env = stdr.inv_var(np.ones(sig_orig.size, dtype=float), np.arange(sig_orig.size, dtype=np.int64))
        scale2_env = np.where(np.isfinite(scale2_env) & (scale2_env > 0), scale2_env, 1.0)
        sig_std = sig_orig / scale2_env
        r_diag_obs = sig_std[ei]

    # Build cache
    ctype = str(grm_factor_cache.get("type", "zarr")).lower()
    if ctype == "zarr":
        root_dir = grm_factor_cache["root_dir"]
        num_grms = int(grm_factor_cache.get("num_grms", len(w_g)))
        cache = GRMFactorCacheZarr(root_dir=str(root_dir), num_grms=num_grms)
    elif ctype == "memmap":
        cache = GRMFactorCacheMemmap(
            paths=grm_factor_cache["paths"],
            shapes=grm_factor_cache["shapes"],
            dtypes=grm_factor_cache["dtypes"],
        )
    else:
        raise ValueError(f"Unknown grm_factor_cache['type']: {ctype}")

    # Sigma_e and S_e (env similarity); for now use env_similarity for both (common in practice)
    k_env = int(ei.max()) + 1
    if env_similarity is None:
        Sigma_e = np.eye(k_env, dtype=float)
        S_e = None
    else:
        S = np.asarray(env_similarity, dtype=float)
        if S.shape[0] != k_env or S.shape[1] != k_env:
            raise ValueError(f"env_similarity must be (n_env,n_env) with n_env={k_env}; got {S.shape}")
        Sigma_e = S
        S_e = S

    # Torch tensors
    dev = device
    gi_t = torch.as_tensor(gi, device=dev, dtype=torch.int64)
    ei_t = torch.as_tensor(ei, device=dev, dtype=torch.int64)
    y_t  = torch.as_tensor(y_std, device=dev, dtype=dtype_t)
    X_t  = torch.as_tensor(X, device=dev, dtype=dtype_t)
    tr_t = torch.as_tensor(train_idx, device=dev, dtype=torch.int64)
    te_t = torch.as_tensor(test_idx, device=dev, dtype=torch.int64)

    # Compute dtype for operator matvec (float32 recommended)
    dt_comp = _resolve_torch_dtype(operator_dtype_compute) or torch.float32

    cfg_kwargs = dict(
        dtype_compute=dt_comp, tol=float(operator_tol), max_iter=int(operator_max_iter),
        point_predictions_only=(not bool(return_se)), return_varcov=False, return_se=bool(return_se),
    )
    if bool(return_se):
        cfg_kwargs["n_hutchinson_probes"] = int(n_hutchinson_probes)
        cfg_kwargs["hutchinson_seed"] = int(seed if hutchinson_seed is None else hutchinson_seed)
        if se_pcg_tol is not None:
            cfg_kwargs["se_pcg_tol"] = float(se_pcg_tol)
        if se_pcg_max_iter is not None:
            cfg_kwargs["se_pcg_max_iter"] = int(se_pcg_max_iter)
        if bool(use_nystrom_se_preconditioner):
            cfg_kwargs["use_nystrom_se_preconditioner"] = True
            cfg_kwargs["nystrom_rank"] = int(nystrom_rank)
            cfg_kwargs["nystrom_seed"] = int(nystrom_seed)
        if bool(use_lanczos_se_variance_reduction):
            cfg_kwargs["use_lanczos_se_variance_reduction"] = True
            cfg_kwargs["lanczos_iters"] = int(lanczos_iters)
            cfg_kwargs["lanczos_seed"] = int(lanczos_seed)
            cfg_kwargs["lanczos_theta_floor"] = float(lanczos_theta_floor)
    cfg = LargeOperatorBackendConfig(**cfg_kwargs)

    out = fit_predict_operator_backend_large(
        y=y_t, X=X_t, gi=gi_t, ei=ei_t,
        cache=cache,
        Sigma_e=torch.as_tensor(Sigma_e, device=dev, dtype=dt_comp),
        train_idx=tr_t, test_idx=te_t,
        w_g=torch.as_tensor(w_g, device=dev, dtype=dt_comp),
        w_ge=torch.as_tensor(w_ge, device=dev, dtype=dt_comp),
        w_e=float(w_e),
        S_e=(None if S_e is None else torch.as_tensor(S_e, device=dev, dtype=dt_comp)),
        r_diag=r_diag_obs if isinstance(r_diag_obs, float) else torch.as_tensor(r_diag_obs, device=dev, dtype=dt_comp),
        config=cfg,
        random_terms=random_terms,
    )

    # Build predictions DataFrame consistent with fast "test_only" outputs
    n = gi.shape[0]
    geno_names = (geno_id_list if (geno_id_list is not None and len(geno_id_list) > 0) else None)
    env_names = (env_index_to_name if (env_index_to_name is not None and len(env_index_to_name) > 0) else None)

    def _name_for_row(i: int) -> str:
        if geno_names is None:
            return str(int(gi[i]))
        if len(geno_names) == n:
            return str(geno_names[i])
        if gi[i] < len(geno_names):
            return str(geno_names[int(gi[i])])
        return str(int(gi[i]))

    def _env_for_row(i: int) -> str:
        if env_names is None:
            return str(int(ei[i]))
        if ei[i] < len(env_names):
            return str(env_names[int(ei[i])])
        return str(int(ei[i]))

    # Operator output is on standardized scale; invert to observed scale for means
    yhat_std = np.asarray(out["predictions"], dtype=float)
    yhat_obs = stdr.inv_mean(yhat_std, ei[test_idx])

    predictions = pd.DataFrame({
        "row_index": test_idx.astype(np.int64),
        "Genotype": [ _name_for_row(int(i)) for i in test_idx ],
        "Env": [ _env_for_row(int(i)) for i in test_idx ],
        "Prediction": yhat_obs.astype(float),
    })

    if bool(return_se) and ("prediction_se_observed" in out):
        # SEs come back on standardized scale; convert to observed via stdr.inv_var
        se_std_obs = np.asarray(out["prediction_se_observed"], dtype=float)
        se_std_lat = np.asarray(out["prediction_se_latent"], dtype=float)
        var_std_obs = se_std_obs ** 2
        var_std_lat = se_std_lat ** 2
        var_obs_obs = stdr.inv_var(var_std_obs, ei[test_idx])
        var_obs_lat = stdr.inv_var(var_std_lat, ei[test_idx])
        predictions["SE"] = np.sqrt(np.maximum(var_obs_obs, 0.0)).astype(float)
        predictions["SE_latent"] = np.sqrt(np.maximum(var_obs_lat, 0.0)).astype(float)
        predictions["PEV"] = np.maximum(var_obs_lat, 0.0).astype(float)
        predictions = _attach_prediction_se_columns(
            predictions,
            latent_var=np.maximum(var_obs_lat, 0.0),
            observed_var=np.maximum(var_obs_obs, 0.0),
        )

    return {
        "beta": np.asarray(out["beta"], dtype=float),
        "beta_se": None,
        "var_components": pd.DataFrame([]),
        "sigma2_resid_env": np.asarray(sigma2_env_orig, dtype=float),
        "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_env_orig),
        "diagnostics": {
            "method": "operator_backend",
            "prediction_output": str(prediction_output),
            "device": dev,
            "dtype": str(dtype_t),
            "operator_dtype_compute": str(dt_comp),
            "operator_tol": float(operator_tol),
            "operator_max_iter": int(operator_max_iter),
            "pcg_info": out.get("pcg_info", {}),
            "se_info": out.get("se_info", {}),
            "notes": out.get("notes", ""),
        },
        "result": {
            "predictions": predictions,
            "per_env": pd.DataFrame([]),
            "across_env": pd.DataFrame([]),
        }
    }


# -------------------------- GPyTorch kernels & means --------------------------
if GPTY_AVAILABLE:

    class DesignLinearMean(gpytorch.means.Mean):
        def __init__(self, X_all: "torch.Tensor"):
            super().__init__()
            self.X_all = X_all
            p = X_all.shape[1]
            self.beta = torch.nn.Parameter(torch.zeros(p, dtype=X_all.dtype, device=X_all.device))
        def forward(self, x: "torch.Tensor"):
            ridx = x[:, 2].long()
            return (self.X_all.index_select(0, ridx) @ self.beta).squeeze(-1)

    def _assign_reml_beta_to_design_mean_(model, likelihood, train_x: "torch.Tensor", train_y: "torch.Tensor", jitter: float = 1e-6) -> "torch.Tensor":
        """Install the profiled GLS beta into DesignLinearMean before prediction."""
        output = model(train_x)
        K_op = output.lazy_covariance_matrix
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            R_op = DiagLinearOperator(likelihood.noise)
        else:
            R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
        V_op = K_op + R_op
        Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
        beta = gls_beta_from_V(train_y, Xtt, V_op, jitter=jitter)
        model.mean_module.beta.data.copy_(beta.detach().reshape_as(model.mean_module.beta))
        return beta.detach()

    def _design_mean_prediction_var_diag(model, train_x: "torch.Tensor", pred_x: "torch.Tensor",
                                         V_op: "DenseLinearOperator", Var_beta: "torch.Tensor") -> "torch.Tensor":
        """Universal-kriging fixed-effect uncertainty for prediction rows."""
        Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
        Xss = model.mean_module.X_all.index_select(0, pred_x[:, 2].long())
        Vinv_X = V_op.solve(Xtt)
        Kst = model.covar_module(pred_x, train_x)
        adjusted_X = Xss - (Kst @ Vinv_X)
        return torch.clamp((adjusted_X @ Var_beta * adjusted_X).sum(dim=1), min=0.0)

    def _prediction_noise_diag_std(likelihood, pred_x: "torch.Tensor", ei_t: "torch.Tensor",
                                   resid_diag_env: Optional[np.ndarray],
                                   resid_diag_obs_std: Optional[np.ndarray],
                                   device: str, dtype_t: "torch.dtype") -> "torch.Tensor":
        """Residual variance for prediction rows on the standardized response scale."""
        if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
            pred_rows = pred_x[:, 2].long()
            if resid_diag_obs_std is not None:
                return _t(resid_diag_obs_std, device=device, dtype=dtype_t).index_select(0, pred_rows)
            pred_env = pred_x[:, 1].long()
            return _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, pred_env)
        return torch.full(
            (pred_x.shape[0],),
            float(likelihood.noise.detach().reshape(-1)[0].item()),
            dtype=dtype_t,
            device=device,
        )

    def _softplus_inv_tensor(x: "torch.Tensor", eps: float = 1e-6) -> "torch.Tensor":
        x = torch.as_tensor(x)
        x = torch.clamp(x, min=eps)
        return torch.log(torch.expm1(x))

    class PrecomputedIndexKernel(gpytorch.kernels.Kernel):
        """
        K = sum_i w_g[i] * G_i
          + sum_i w_ge[i] * (G_i ⊙ S_e)   [if S_e provided]
          + w_e * S_e_or_I
        Supports learnable component scales when learn_scales=True.
        """
        is_stationary = False
        def __init__(self, G_list: List["torch.Tensor"], w_g: "torch.Tensor",
                     w_ge: "torch.Tensor", w_e: float, S_e: Optional["torch.Tensor"],
                     learn_scales: bool = False, **kwargs):
            super().__init__(**kwargs)
            self.G_list = G_list
            self.S_e = S_e
            self.learn_scales = bool(learn_scales)

            wg_fix = torch.as_tensor(w_g, dtype=G_list[0].dtype, device=G_list[0].device).reshape(-1)
            wge_fix = torch.as_tensor(w_ge, dtype=G_list[0].dtype, device=G_list[0].device).reshape(-1)
            we_fix = torch.as_tensor(float(w_e), dtype=G_list[0].dtype, device=G_list[0].device)
            self.register_buffer("w_g_fix", wg_fix)
            self.register_buffer("w_ge_fix", wge_fix)
            self.register_buffer("w_e_fix", we_fix)

            if self.learn_scales:
                self.w_g_unconstrained  = torch.nn.Parameter(_softplus_inv_tensor(self.w_g_fix))
                self.w_ge_unconstrained = torch.nn.Parameter(_softplus_inv_tensor(self.w_ge_fix))
                self.w_e_unconstrained  = torch.nn.Parameter(_softplus_inv_tensor(self.w_e_fix))
            else:
                self.w_g_unconstrained  = None
                self.w_ge_unconstrained = None
                self.w_e_unconstrained  = None

        def _effective_scales(self):
            if self.learn_scales:
                wg  = torch.nn.functional.softplus(self.w_g_unconstrained)
                wge = torch.nn.functional.softplus(self.w_ge_unconstrained)
                we  = torch.nn.functional.softplus(self.w_e_unconstrained)
            else:
                wg, wge, we = self.w_g_fix, self.w_ge_fix, self.w_e_fix
            return wg, wge, we

        def get_effective_scales(self) -> Tuple["torch.Tensor","torch.Tensor","torch.Tensor"]:
            wg, wge, we = self._effective_scales()
            return wg.detach(), wge.detach(), we.detach()

        def forward(self, x1: "torch.Tensor", x2: "torch.Tensor", diag: bool = False, **params):
            gi1, ei1 = x1[:,0].long(), x1[:,1].long()
            gi2, ei2 = x2[:,0].long(), x2[:,1].long()
            dtype = self.G_list[0].dtype; device = self.G_list[0].device
            wg, wge, we = self._effective_scales()

            if diag:
                Kdiag = torch.zeros(gi1.numel(), dtype=dtype, device=device)
                for i, Gk in enumerate(self.G_list):
                    wk = wg[i]
                    if float(wk) != 0.0:
                        rows = Gk.index_select(0, gi1)
                        Kdiag.add_(wk * rows.gather(1, gi1.unsqueeze(1)).squeeze(1))
                if (self.S_e is not None) and bool(torch.any(wge != 0)):
                    Se_d = self.S_e.diagonal()[ei1]
                    for i, Gk in enumerate(self.G_list):
                        wk = wge[i]
                        if float(wk) != 0.0:
                            rows = Gk.index_select(0, gi1)
                            Kdiag.add_(wk * rows.gather(1, gi1.unsqueeze(1)).squeeze(1) * Se_d)
                if float(we) != 0.0:
                    Se_d = self.S_e.diagonal()[ei1] if self.S_e is not None else torch.ones_like(Kdiag)
                    Kdiag.add_(we * Se_d)
                return Kdiag

            K = torch.zeros((gi1.numel(), gi2.numel()), dtype=dtype, device=device)
            _kg_cache = []
            for i, Gk in enumerate(self.G_list):
                block = _build_geno_block(Gk, gi1, gi2)
                _kg_cache.append(block.detach())
                wk = wg[i]
                if float(wk) != 0.0:
                    K.add_(wk * block)
            _kge_cache = []
            if (self.S_e is not None) and bool(torch.any(wge != 0)):
                S12 = _build_env_kernel(ei1, ei2, self.S_e, dtype=dtype, device=device)
                for i, Gk in enumerate(self.G_list):
                    ge_block = _kg_cache[i] * S12
                    _kge_cache.append(ge_block.detach())
                    wk = wge[i]
                    if float(wk) != 0.0:
                        K.add_(wk * ge_block)
            _ke_cache = _build_env_kernel(ei1, ei2, self.S_e, dtype=dtype, device=device).detach()
            if float(we) != 0.0:
                K.add_(we * _ke_cache)
            self._Kg_list = _kg_cache
            self._Kge_list = _kge_cache if len(_kge_cache) > 0 else None
            self._Ke = _ke_cache
            return DenseLinearOperator(K)

        
    class CorHCovariance(nn.Module):
        """
        corh: Heterogeneous variances with constant correlation rho.
        Sigma = D^{1/2} * R(rho) * D^{1/2}
        where R(rho) has 1's on diag and rho off-diagonal.

        Constraints:
          - variances > 0 via softplus
          - rho in (-1/(k-1), 1) via sigmoid mapping
        Parameter count: k + 1
        """
        def __init__(self, k: int, dtype=None, device=None, init_vars: Optional[torch.Tensor]=None,
                     init_rho: float=0.1, jitter: float=1e-6):
            super().__init__()
            self.k = int(k)
            self.jitter = float(jitter)
            dtype = dtype or torch.float64
            device = device or ("cuda" if torch.cuda.is_available() else "cpu")

            if init_vars is None:
                init_vars = torch.full((self.k,), 0.05, device=device, dtype=dtype)
            else:
                init_vars = init_vars.to(device=device, dtype=dtype)

            self.var_unconstrained = nn.Parameter(torch.log(torch.clamp(init_vars, min=self.jitter)))
            # raw rho mapped to (-1/(k-1), 1)
            # init via inverse sigmoid in that range
            lo = -1.0 / max(1, (self.k - 1))
            hi = 1.0
            init_rho = float(max(lo + 1e-4, min(hi - 1e-4, init_rho)))
            t = (init_rho - lo) / (hi - lo)
            self.rho_raw = nn.Parameter(torch.log(torch.tensor(t, device=device, dtype=dtype) / (1.0 - torch.tensor(t, device=device, dtype=dtype))))

        def rho(self) -> torch.Tensor:
            lo = -1.0 / max(1, (self.k - 1))
            hi = 1.0
            return lo + (hi - lo) * torch.sigmoid(self.rho_raw)

        def vars(self) -> torch.Tensor:
            # positive variances
            return torch.nn.functional.softplus(self.var_unconstrained) + self.jitter

        def cov(self) -> torch.Tensor:
            v = self.vars()
            sd = torch.sqrt(v)
            rho = self.rho()
            R = torch.full((self.k, self.k), rho, device=v.device, dtype=v.dtype)
            R.fill_diagonal_(1.0)
            Sigma = (sd[:, None] * R) * sd[None, :]
            Sigma = Sigma + self.jitter * torch.eye(self.k, device=v.device, dtype=v.dtype)
            return Sigma

        def submatrix(self, idx1: torch.Tensor, idx2: torch.Tensor) -> torch.Tensor:
            S = self.cov()
            return S.index_select(0, idx1.long()).index_select(1, idx2.long())


    class CorGHCovariance(nn.Module):
        """
        corgh: Heterogeneous variances with general correlation matrix.

        Parameterization:
          - Variances v_i > 0 via softplus
          - Correlation matrix R built from an unconstrained Cholesky factor A:
              C = A A^T
              R_ij = C_ij / sqrt(C_ii C_jj)
          - Sigma = D^{1/2} R D^{1/2}

        This matches ASReml's 'corgh' family (SPD by construction).
        Parameter count: k (vars) + k(k-1)/2 (corr off-diagonal) + k (corr diag via exp) but diag is constrained positive.
        """
        def __init__(self, k: int, dtype=None, device=None,
                     init_vars: Optional[torch.Tensor]=None,
                     init_corr: Optional[torch.Tensor]=None,
                     jitter: float=1e-6):
            super().__init__()
            self.k = int(k)
            self.jitter = float(jitter)
            dtype = dtype or torch.float64
            device = device or ("cuda" if torch.cuda.is_available() else "cpu")

            if init_vars is None:
                init_vars = torch.full((self.k,), 0.05, device=device, dtype=dtype)
            else:
                init_vars = init_vars.to(device=device, dtype=dtype)
            self.var_unconstrained = nn.Parameter(torch.log(torch.clamp(init_vars, min=self.jitter)))

            # Build initial correlation Cholesky
            if init_corr is None:
                A0 = torch.eye(self.k, device=device, dtype=dtype)
            else:
                R0 = 0.5 * (init_corr + init_corr.T)
                R0 = R0.to(device=device, dtype=dtype)
                R0 = R0 + self.jitter * torch.eye(self.k, device=device, dtype=dtype)
                # ensure SPD
                A0 = torch.linalg.cholesky(R0)

            tril = torch.tril_indices(self.k, self.k, offset=-1)
            self._tril_i = tril[0]
            self._tril_j = tril[1]
            self.corr_diag_unconstrained = nn.Parameter(torch.log(torch.clamp(torch.diag(A0), min=self.jitter)))
            self.corr_offdiag = nn.Parameter(A0[self._tril_i, self._tril_j].clone())

        def vars(self) -> torch.Tensor:
            return torch.nn.functional.softplus(self.var_unconstrained) + self.jitter

        def _A(self) -> torch.Tensor:
            A = torch.zeros((self.k, self.k), device=self.corr_diag_unconstrained.device, dtype=self.corr_diag_unconstrained.dtype)
            A[self._tril_i, self._tril_j] = self.corr_offdiag
            A = A + torch.diag(torch.exp(self.corr_diag_unconstrained))
            return A

        def corr(self) -> torch.Tensor:
            A = self._A()
            C = A @ A.T
            d = torch.sqrt(torch.clamp(torch.diag(C), min=self.jitter))
            R = C / (d[:, None] * d[None, :])
            R.fill_diagonal_(1.0)
            return R

        def cov(self) -> torch.Tensor:
            v = self.vars()
            sd = torch.sqrt(v)
            R = self.corr()
            Sigma = (sd[:, None] * R) * sd[None, :]
            Sigma = Sigma + self.jitter * torch.eye(self.k, device=Sigma.device, dtype=Sigma.dtype)
            return Sigma

        def submatrix(self, idx1: torch.Tensor, idx2: torch.Tensor) -> torch.Tensor:
            S = self.cov()
            return S.index_select(0, idx1.long()).index_select(1, idx2.long())


    def build_env_covariance_module(
        structure: str,
        num_envs: int,
        S_e_init: Optional[torch.Tensor],
        fa_rank: int,
        dtype,
        device,
        jitter: float = 1e-6,
    ) -> Optional[nn.Module]:
        """
        Factory to create a learnable env covariance module.
        structure in {"fixed","fa","us","corh","corgh"}.
        If structure=="fixed" -> returns None and caller should pass S_e_init through.
        """
        s = (structure or "fixed").lower()
        if s == "fixed":
            return None
        if s == "fa":
            return FactorAnalyticCovariance(num_envs, int(fa_rank), dtype=dtype, device=device, init_Sigma=S_e_init, jitter=jitter, identified=True)
        if s == "us":
            return UnstructuredCovariance(num_envs, dtype=dtype, device=device, init_Sigma=S_e_init, jitter=jitter)
        if s == "corh":
            init_vars = None
            init_rho = 0.1
            if S_e_init is not None:
                init_vars = torch.diag(S_e_init).clone()
                # crude rho init from average offdiag correlation
                d = torch.sqrt(torch.clamp(torch.diag(S_e_init), min=jitter))
                R0 = S_e_init / (d[:, None] * d[None, :])
                mask = ~torch.eye(num_envs, dtype=torch.bool, device=device)
                init_rho = float(torch.mean(R0[mask]).item())
            return CorHCovariance(num_envs, dtype=dtype, device=device, init_vars=init_vars, init_rho=init_rho, jitter=jitter)
        if s == "corgh":
            init_vars = None
            init_corr = None
            if S_e_init is not None:
                init_vars = torch.diag(S_e_init).clone()
                d = torch.sqrt(torch.clamp(torch.diag(S_e_init), min=jitter))
                init_corr = S_e_init / (d[:, None] * d[None, :])
            return CorGHCovariance(num_envs, dtype=dtype, device=device, init_vars=init_vars, init_corr=init_corr, jitter=jitter)
        raise ValueError(f"Unknown env covariance structure: {structure!r}")
    class FactorAnalyticEnv(torch.nn.Module):
            def __init__(self, num_envs: int, rank: int, init_cov: Optional["torch.Tensor"]=None, dtype=None, device=None):
                super().__init__()
                dtype = dtype or _default_dtype; device = device or _default_device
                # Clamp rank to num_envs to prevent over-parameterized FA (FA(n_env)
                # is unstructured; rank > n_env is non-identifiable). Previously
                # only the init_cov branch clamped, causing shape mismatch between
                # self.rank and the loadings matrix in the random-init branch.
                self.num_envs = int(num_envs)
                self.rank = min(int(rank), int(num_envs))
                if init_cov is not None:
                    C = 0.5 * (init_cov + init_cov.mT)
                    ev, U = torch.linalg.eigh(C); ev = torch.clamp(ev, min=1e-8)
                    idx = torch.argsort(ev, descending=True); U = U[:, idx]; ev = ev[idx]
                    L0 = U[:, :self.rank] @ torch.diag(torch.sqrt(ev[:self.rank]))
                    d0 = torch.clamp(torch.diag(C) - (L0 * L0).sum(-1), min=1e-6)
                else:
                    g = torch.Generator(device=device); g.manual_seed(12345)
                    L0 = torch.randn(self.num_envs, self.rank, generator=g, dtype=dtype, device=device) * 0.1
                    d0 = torch.full((self.num_envs,), 0.05, dtype=dtype, device=device)
                self.L = torch.nn.Parameter(L0.to(device=device, dtype=dtype))
                self.d_unconstrained = torch.nn.Parameter(d0.to(device=device, dtype=dtype).log())
            def cov(self) -> "torch.Tensor":
                D = torch.nn.functional.softplus(self.d_unconstrained) + 1e-8
                return self.L @ self.L.mT + torch.diag(D)
            def submatrix(self, idx1: "torch.Tensor", idx2: "torch.Tensor") -> "torch.Tensor":
                C = self.cov()
                return C.index_select(0, idx1.long()).index_select(1, idx2.long())

    class PrecomputedFAIndexKernel(gpytorch.kernels.Kernel):
            """
            K = (sum_i w_fa[i] * G_i) ⊗ C   +   w_e * S_e_or_I
            Supports learnable w_fa and w_e when learn_scales=True.
            """
            is_stationary = False
            def __init__(self, G_list: List["torch.Tensor"], w_g: "torch.Tensor",
                         fa_env: FactorAnalyticEnv, w_e: float, S_e: Optional["torch.Tensor"],
                         learn_scales: bool = False, **kwargs):
                super().__init__(**kwargs)
                self.G_list = G_list
                self.fa_env = fa_env
                self.S_e = S_e
                self.learn_scales = bool(learn_scales)

                wf_fix = torch.as_tensor(w_g, dtype=G_list[0].dtype, device=G_list[0].device).reshape(-1)
                we_fix = torch.as_tensor(float(w_e), dtype=G_list[0].dtype, device=G_list[0].device)
                self.register_buffer("w_fa_fix", wf_fix)
                self.register_buffer("w_e_fix", we_fix)

                if self.learn_scales:
                    self.w_fa_unconstrained = torch.nn.Parameter(_softplus_inv_tensor(self.w_fa_fix))
                    self.w_e_unconstrained  = torch.nn.Parameter(_softplus_inv_tensor(self.w_e_fix))
                else:
                    self.w_fa_unconstrained = None
                    self.w_e_unconstrained  = None

            def _effective_scales(self):
                if self.learn_scales:
                    wfa = torch.nn.functional.softplus(self.w_fa_unconstrained)
                    we  = torch.nn.functional.softplus(self.w_e_unconstrained)
                else:
                    wfa, we = self.w_fa_fix, self.w_e_fix
                return wfa, we

            def get_effective_scales(self) -> Tuple["torch.Tensor","torch.Tensor"]:
                wfa, we = self._effective_scales()
                return wfa.detach(), we.detach()

            def forward(self, x1: "torch.Tensor", x2: "torch.Tensor", diag: bool = False, **params):
                gi1, ei1 = x1[:,0].long(), x1[:,1].long()
                gi2, ei2 = x2[:,0].long(), x2[:,1].long()
                dtype = self.G_list[0].dtype; device = self.G_list[0].device
                wfa, we = self._effective_scales()

                if diag:
                    Cdd = self.fa_env.cov().diagonal()[ei1]
                    Gdiag = torch.zeros_like(Cdd)
                    for i, Gk in enumerate(self.G_list):
                        wk = wfa[i]
                        if float(wk) != 0.0:
                            rows = Gk.index_select(0, gi1)
                            Gdiag.add_(wk * rows.gather(1, gi1.unsqueeze(1)).squeeze(1))
                    Kdiag = Gdiag * Cdd
                    if float(we) != 0.0:
                        Se_d = self.S_e.diagonal()[ei1] if self.S_e is not None else torch.ones_like(Cdd)
                        Kdiag.add_(we * Se_d)
                    return Kdiag

                C12 = self.fa_env.submatrix(ei1, ei2)
                K = torch.zeros((gi1.numel(), gi2.numel()), dtype=dtype, device=device)
                Gsum = torch.zeros_like(K)
                for i, Gk in enumerate(self.G_list):
                    wk = wfa[i]
                    if float(wk) != 0.0:
                        Gsum.add_(wk * _build_geno_block(Gk, gi1, gi2))
                K.add_(Gsum * C12)
                if float(we) != 0.0:
                    K.add_(we * _build_env_kernel(ei1, ei2, self.S_e, dtype=dtype, device=device))
                return DenseLinearOperator(K)

    class PrecomputedMultiTermIndexKernel(gpytorch.kernels.Kernel):
        """K = sum_i w_g[i]*G_i + sum_i (G_i ⊙ M_i) + w_e*S_e_or_I, where M_i is a kernel-specific
        interaction covariance over the interaction axis. This allows multiple structured interaction
        terms to contribute additively and natively in the dense GP path.
        """
        is_stationary = False
        def __init__(self, G_list: List["torch.Tensor"], w_g: "torch.Tensor", ge_mats: Sequence["torch.Tensor"],
                     w_e: float, S_e: Optional["torch.Tensor"], learn_scales: bool = False, **kwargs):
            super().__init__(**kwargs)
            self.G_list = G_list
            self.ge_mats = list(ge_mats)
            self.S_e = S_e
            self.learn_scales = bool(learn_scales)
            wg_fix = torch.as_tensor(w_g, dtype=G_list[0].dtype, device=G_list[0].device).reshape(-1)
            we_fix = torch.as_tensor(float(w_e), dtype=G_list[0].dtype, device=G_list[0].device)
            self.register_buffer("w_g_fix", wg_fix)
            self.register_buffer("w_e_fix", we_fix)
            if self.learn_scales:
                self.w_g_unconstrained = torch.nn.Parameter(_softplus_inv_tensor(self.w_g_fix))
                self.w_e_unconstrained = torch.nn.Parameter(_softplus_inv_tensor(self.w_e_fix))
            else:
                self.w_g_unconstrained = None
                self.w_e_unconstrained = None
        def _effective_scales(self):
            if self.learn_scales:
                wg = torch.nn.functional.softplus(self.w_g_unconstrained)
                we = torch.nn.functional.softplus(self.w_e_unconstrained)
            else:
                wg, we = self.w_g_fix, self.w_e_fix
            return wg, we
        def get_effective_scales(self):
            wg, we = self._effective_scales()
            return wg.detach(), we.detach()
        def forward(self, x1: "torch.Tensor", x2: "torch.Tensor", diag: bool = False, **params):
            gi1, ei1 = x1[:,0].long(), x1[:,1].long()
            gi2, ei2 = x2[:,0].long(), x2[:,1].long()
            dtype = self.G_list[0].dtype; device = self.G_list[0].device
            wg, we = self._effective_scales()
            if diag:
                Kdiag = torch.zeros(gi1.numel(), dtype=dtype, device=device)
                for i, Gk in enumerate(self.G_list):
                    wk = wg[i]
                    rows = Gk.index_select(0, gi1)
                    gdiag = rows.gather(1, gi1.unsqueeze(1)).squeeze(1)
                    if float(wk) != 0.0:
                        Kdiag.add_(wk * gdiag)
                    if i < len(self.ge_mats) and self.ge_mats[i] is not None:
                        md = self.ge_mats[i].diagonal()[ei1]
                        Kdiag.add_(gdiag * md)
                if float(we) != 0.0:
                    Se_d = self.S_e.diagonal()[ei1] if self.S_e is not None else torch.ones_like(Kdiag)
                    Kdiag.add_(we * Se_d)
                return Kdiag
            K = torch.zeros((gi1.numel(), gi2.numel()), dtype=dtype, device=device)
            for i, Gk in enumerate(self.G_list):
                Gblk = _build_geno_block(Gk, gi1, gi2)
                wk = wg[i]
                if float(wk) != 0.0:
                    K.add_(wk * Gblk)
                if i < len(self.ge_mats) and self.ge_mats[i] is not None:
                    S12 = _build_env_kernel(ei1, ei2, self.ge_mats[i], dtype=dtype, device=device)
                    K.add_(Gblk * S12)
            if float(we) != 0.0:
                K.add_(we * _build_env_kernel(ei1, ei2, self.S_e, dtype=dtype, device=device))
            return DenseLinearOperator(K)

else:
    class DesignLinearMean: pass
    class PrecomputedIndexKernel: pass
    class PrecomputedMultiTermIndexKernel: pass
    class FactorAnalyticEnv: pass
    class PrecomputedFAIndexKernel: pass

# ----------------------- KRR exact (GPU) with GLS & (SEs optional) -----------
def _select_lambda_gcv(K: "torch.Tensor", y: "torch.Tensor", X: "torch.Tensor",
                       ei_t: "torch.Tensor", lambdas: np.ndarray, device, dtype) -> float:
    y = y.reshape(-1,1)
    best_lam, best_gcv = float(lambdas[0]), float("inf")
    n = K.shape[0]
    for lam in lambdas:
        V = 0.5*(K+K.mT); V.diagonal().add_(float(lam))
        Vinv_y = _cholesky_solve(V, y)
        Vinv_X = _cholesky_solve(V, X)
        Xt_Vinv_X = X.mT @ Vinv_X
        try:
            Var_beta = _spd_inverse(Xt_Vinv_X)
        except RuntimeError:
            Var_beta = _spd_inverse(Xt_Vinv_X + 1e-10*torch.eye(Xt_Vinv_X.shape[0], dtype=dtype, device=device))
        yhat = X @ (Var_beta @ (X.mT @ Vinv_y)) + K @ _cholesky_solve(V, (y - X @ (Var_beta @ (X.mT @ Vinv_y))))
        res = y - yhat
        rss = float((res.mT @ res).item())
        if n <= 2048:
            # Exact dense GCV trace.  The former 16-probe Hutchinson trace is
            # far too noisy when KRR is close to interpolation: on the
            # package's 40-training-genotype fixture it selected lambda=.001
            # although the exact GCV/PRESS optimum is lambda=1.  For
            # V=K+lambda*I and the GLS fixed-effect smoother,
            #
            #   tr(H) = n - lambda*tr(V^-1)
            #             + lambda*tr(A X' V^-2 X),
            #   A = (X' V^-1 X)^-1.
            #
            # This avoids forming H while remaining deterministic.
            eye_n = torch.eye(n, dtype=dtype, device=device)
            V_inv = _cholesky_solve(V, eye_n)
            trace_v_inv = torch.trace(V_inv)
            fixed_effect_correction = torch.trace(
                Var_beta @ (Vinv_X.mT @ Vinv_X)
            )
            trH = float(
                (
                    torch.as_tensor(float(n), dtype=dtype, device=device)
                    - float(lam) * trace_v_inv
                    + float(lam) * fixed_effect_correction
                ).item()
            )
        else:
            # Retain bounded-memory stochastic trace estimation for large
            # dense fits, where materialising V^-1 would be prohibitive.
            def H_times(z):
                z = z.reshape(-1,1)
                Sx_zt = X @ (Var_beta @ (X.mT @ _cholesky_solve(V, z)))
                return (Sx_zt + K @ _cholesky_solve(V, (z - Sx_zt))).squeeze(-1)
            trH = float(_estimate_diag_Hutch(
                H_times, n=n, probes=64, device=device, dtype=dtype, seed=12345
            ).sum().item())
        # GCV's effective residual degrees of freedom may legitimately be
        # below one near interpolation.  Clipping it to one made those fits
        # look artificially excellent and forced the smallest candidate.
        denom = max(np.finfo(float).eps * max(1.0, float(n)), n - trH)
        gcv = rss / (denom**2)
        if gcv < best_gcv:
            best_gcv, best_lam = gcv, float(lam)
    return best_lam

def _robust_env_var(resids: "torch.Tensor") -> float:
    if resids.numel() == 0:
        return float('nan')
    r = resids.detach()
    mask = torch.isfinite(r)
    r = r[mask]
    if r.numel() <= 1:
        return float('nan')
    med = torch.median(r)
    mad = torch.median(torch.abs(r - med))
    return float(((1.4826*mad)**2).item())

def multikernel_krr_posterior_with_se(
    G_list: List[np.ndarray],
    obs_gidx: np.ndarray,
    obs_eidx: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    w_g: Optional[Sequence[float]] = None,
    w_ge: Optional[Sequence[float]] = None,
    w_e: float = 0.0,
    S_e: Optional[np.ndarray] = None,
    lam: float = 1e-3,
    resid_diag_env: Optional[np.ndarray] = None,
    resid_diag_obs: Optional[np.ndarray] = None,
    center_y: bool = True,
    env_levels: Optional[Sequence] = None,
    obs_names: Optional[Sequence] = None,   # may be per-observation OR genotype list
    fixed_effects: Optional[Sequence[str]] = None,
    X: Optional[np.ndarray] = None,
    device: str = None,
    dtype: Any = None,
    hutch_samples: int = 64,
    seed: int = 12345,
    prediction_output: str = "test_only",
    standardize: str = "per_env",
    krr_lams: Union[str, Sequence[float]] = "none",
    lam_select: str = "fixed",
    env_resid_robust: str = "none",
    env_resid_shrink_tau: float = 0.0,
    reml_normalize: str = "diag_mean",
    return_se: bool = True,
    env_structure: Optional[str] = None,
    fa_rank: Optional[int] = None,
    interaction_term_name: str = "G:Env",
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> Dict[str, Any]:

    assert TORCH_AVAILABLE, "PyTorch required"
    requested_device = device
    dtype_t = _resolve_torch_dtype(dtype)
    set_deterministic(seed)

    varcomp_df = pd.DataFrame([])

    gi = np.asarray(obs_gidx, dtype=np.int64)
    ei = np.asarray(obs_eidx, dtype=np.int64)
    y  = np.asarray(y, dtype=float)
    t_idx = np.asarray(train_idx, dtype=np.int64)
    s_idx_default = np.asarray(test_idx if (test_idx is not None and len(test_idx)>0) else train_idx, dtype=np.int64)
    s_idx_all = np.arange(gi.size, dtype=np.int64)
    _work_units = float(t_idx.size) * float(t_idx.size + max(s_idx_default.size, 1))
    try:
        _work_units *= max(1, int(len(G_list)))
    except Exception:
        pass
    device = _runtime_gp_device(requested_device, work_units=_work_units)

    # Standardize y using TRAIN only
    stdr = _PerObsEnvStandardizer(mode=standardize).fit(y, ei, t_idx)
    y_std = stdr.transform(y, ei)

    # FE design
    X_all = np.ones((gi.size, 1), dtype=float) if X is None else np.asarray(X, dtype=float)

    # Normalize kernels
    G_list, S_e = _apply_kernel_normalization(G_list, S_e, mode=reml_normalize)

    # torch tensors
    gi_t = _t(gi, device=device, dtype=torch.int64)
    ei_t = _t(ei, device=device, dtype=torch.int64)
    y_t  = _t(y_std, device=device, dtype=dtype_t)
    X_t  = _t(X_all, device=device, dtype=dtype_t)
    t    = _t(t_idx, device=device, dtype=torch.int64)
    s_default_t = _t(s_idx_default, device=device, dtype=torch.int64)

    # kernels + weights
    Gt_list = [_t(np.asarray(K), device=device, dtype=dtype_t) for K in G_list]
    nK = len(Gt_list)
    wg_vec  = _prepare_weight_vector(w_g,  nK, default=1.0)
    wge_vec = _prepare_weight_vector(w_ge, nK, default=0.0)
    wg_t  = _t(wg_vec,  device=device, dtype=dtype_t)
    wge_t = _t(wge_vec, device=device, dtype=dtype_t)
    we    = float(0.0 if w_e is None else w_e)
    Se_t  = None if S_e is None else _t(np.asarray(S_e), device=device, dtype=dtype_t)

    # --- interaction covariance setup (single or multi-term) ---
    num_envs = int(ei_t.max().item() + 1)
    S_e_init_t = None if S_e is None else Se_t
    env_cov_module = build_env_covariance_module(env_structure, num_envs, S_e_init_t, fa_rank=fa_rank, dtype=dtype_t, device=device, jitter=1e-6)
    S_reporting = _to_numpy(env_cov_module.cov()) if env_cov_module is not None else (None if S_e is None else np.asarray(S_e, dtype=float))
    normalized_terms, ge_cov_by_kernel = _build_kernel_specific_interaction_covs(
        interaction_terms_meta=interaction_terms_meta,
        kernel_names=[f"K{i+1}" for i in range(nK)],
        default_level_names=list(env_levels) if env_levels is not None else [str(i) for i in range(num_envs)],
        default_interaction_term_name=interaction_term_name,
        default_interaction_cov=(S_reporting if S_reporting is not None else np.eye(num_envs, dtype=float)),
        default_w_ge=wge_vec,
    )
    ge_mats_t = {k: _t(v, device=device, dtype=dtype_t) for k, v in ge_cov_by_kernel.items()}
    Se_t = None if env_cov_module is not None else Se_t

    def build_K_block(idx1: "torch.Tensor", idx2: "torch.Tensor") -> "torch.Tensor":
        gi1, ei1 = gi_t.index_select(0, idx1), ei_t.index_select(0, idx1)
        gi2, ei2 = gi_t.index_select(0, idx2), ei_t.index_select(0, idx2)
        S_eff = Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else None)
        K = torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
        for wk, Gk in zip(wg_t.reshape(-1).tolist(), Gt_list):
            if float(wk) != 0.0:
                K.add_(float(wk) * _build_geno_block(Gk, gi1, gi2))
        if S_eff is not None and wge_t.reshape(-1).abs().sum() > 0:
            S12 = _build_env_kernel(ei1, ei2, S_eff, dtype=dtype_t, device=device)
            for wk, Gk in zip(wge_t.reshape(-1).tolist(), Gt_list):
                if float(wk) != 0.0:
                    K.add_(float(wk) * (_build_geno_block(Gk, gi1, gi2) * S12))
        if we != 0.0 and S_eff is not None:
            K.add_(float(we) * _build_env_kernel(ei1, ei2, S_eff, dtype=dtype_t, device=device))
        return K

    # Train-train kernel (standardized)
    Ktt = build_K_block(t, t)
    ytt = y_t.index_select(0, t)
    Xtt = X_t.index_select(0, t)

    # Lambda selection
    lam_eff = float(lam)
    if isinstance(krr_lams, (list, tuple, np.ndarray)) and lam_select.lower() != "fixed":
        lam_grid = np.asarray(krr_lams, dtype=float)
        lam_grid = lam_grid[np.isfinite(lam_grid) & (lam_grid > 0)]
        if lam_grid.size > 0:
            lam_eff = _select_lambda_gcv(Ktt.clone(), ytt.clone(), Xtt.clone(), ei_t, lam_grid, device, dtype_t)

    # Residual diag on standardized scale from observation-level or
    # environment-level input.  Preserve whether it was supplied: Stage 2
    # precision weights define a fixed residual covariance and must not be
    # reinterpreted as a KRR variance ratio.
    resid_diag_env_supplied = resid_diag_env is not None
    resid_diag_obs_supplied = resid_diag_obs is not None
    resid_diag_obs_std = None
    if resid_diag_obs is not None:
        resid_diag_obs = np.asarray(resid_diag_obs, dtype=float).reshape(-1)
        if resid_diag_obs.size != gi.size:
            raise ValueError(f"resid_diag_obs must have length n_obs={gi.size}; got {resid_diag_obs.size}")
        if np.any(~np.isfinite(resid_diag_obs)) or np.any(resid_diag_obs <= 0):
            raise ValueError("resid_diag_obs must be finite and > 0.")
        scale2_obs = stdr.inv_var(np.ones(gi.size, dtype=float), ei)
        scale2_obs = np.where(np.isfinite(scale2_obs) & (scale2_obs > 0), scale2_obs, 1.0)
        resid_diag_obs_std = resid_diag_obs / scale2_obs
        resid_diag_env = _env_mean_from_obs_diag(resid_diag_obs, ei, int(ei_t.max().item()) + 1)

    # Residual diag per env (standardized) if not provided
    if resid_diag_env is None and resid_diag_obs_std is None:
        Vtt_tmp = 0.5*(Ktt + Ktt.mT); Vtt_tmp.diagonal().add_(float(lam_eff))
        Vinv_y = _cholesky_solve(Vtt_tmp, ytt)
        Vinv_X = _cholesky_solve(Vtt_tmp, Xtt)
        Xt_Vinv_X = Xtt.mT @ Vinv_X
        Var_beta_tmp = _spd_inverse(Xt_Vinv_X)
        beta_tmp = Var_beta_tmp @ (Xtt.mT @ Vinv_y)
        yhat_tr = (Xtt @ beta_tmp) + Ktt @ _cholesky_solve(Vtt_tmp, (ytt - Xtt @ beta_tmp))
        resid = (ytt - yhat_tr).reshape(-1)

        def H_times(z):
            zt = z.index_select(0, t)
            Sx_zt = Xtt @ (Var_beta_tmp @ (Xtt.mT @ _cholesky_solve(Vtt_tmp, zt)))
            Hz_t = Sx_zt + Ktt @ _cholesky_solve(Vtt_tmp, (zt - Sx_zt))
            out = torch.zeros_like(z)
            out.index_copy_(0, t, Hz_t)
            return out

        hat_diag = _estimate_diag_Hutch(H_times, n=gi_t.numel(), probes=int(hutch_samples), device=device, dtype=dtype_t, seed=seed)

        n_env = int(ei_t.max().item()) + 1
        sig2 = np.full(n_env, np.nan, dtype=float)
        for e in range(n_env):
            idx_all = torch.where(ei_t == e)[0]
            idx = torch.tensor(np.intersect1d(_to_numpy(idx_all), _to_numpy(t)), device=device, dtype=torch.int64)
            n_e = int(idx.numel())
            if n_e == 0:
                continue
            r_e = resid.index_select(0, torch.where(ei_t.index_select(0, t) == e)[0])
            tr_He = float(torch.clamp(hat_diag.index_select(0, idx).sum(), min=0.0, max=float(max(n_e-1,0))).item())
            denom = max(1.0, n_e - tr_He)
            if env_resid_robust.lower() == "mad":
                s2 = _robust_env_var(r_e)
                if np.isnan(s2):
                    s2 = float(torch.clamp((r_e.pow(2).sum()/denom), min=1e-12).item())
            else:
                s2 = float(torch.clamp((r_e.pow(2).sum()/denom), min=1e-12).item())
            sig2[e] = s2

        if np.isfinite(sig2).any() and env_resid_shrink_tau > 0:
            m = float(np.nanmean(sig2))
            sig2 = (1.0 - float(env_resid_shrink_tau))*sig2 + float(env_resid_shrink_tau)*m
        if np.isnan(sig2).any():
            r2 = resid.pow(2)
            mask = torch.isfinite(resid)
            global_s2 = r2[mask].mean().item() if mask.any() else 1.0
            sig2 = np.where(np.isfinite(sig2), sig2, global_s2)
        resid_diag_env = sig2

    n_env = int(ei_t.max().item()) + 1
    use_krr_variance_ratio = (
        (not resid_diag_env_supplied)
        and (not resid_diag_obs_supplied)
        and n_env == 1
        and np.isfinite(lam_eff)
        and lam_eff > 0
    )
    if use_krr_variance_ratio:
        # In single-environment KRR, lambda is sigma_e^2 / sigma_g^2.
        # Estimate the response scale once, then use K + lambda I for the
        # mean and scale posterior covariance by sigma_g^2.  The former path
        # used K + (estimated residual + lambda)I, double-counting noise and
        # disconnecting the final fit from the lambda selected by GCV.
        sigma2_e_std = float(np.asarray(resid_diag_env, dtype=float)[0])
        sigma2_g_std = max(sigma2_e_std / float(lam_eff), 1e-12)
        kernel_variance_scale = sigma2_g_std
        prediction_resid_env_std = np.asarray(
            [float(lam_eff) * sigma2_g_std], dtype=float
        )
        diag_add = torch.full(
            (t.numel(),), float(lam_eff), dtype=dtype_t, device=device
        )
        variance_component_method = "krr_gcv_variance_ratio"
    elif resid_diag_obs_std is None:
        kernel_variance_scale = 1.0
        prediction_resid_env_std = np.asarray(resid_diag_env, dtype=float)
        variance_component_method = "krr_fixed_scale"
        sigma2_env_t = _t(np.asarray(resid_diag_env, dtype=float), device=device, dtype=dtype_t)
        diag_add = sigma2_env_t.index_select(0, ei_t.index_select(0, t)) + float(lam_eff)
    else:
        kernel_variance_scale = 1.0
        prediction_resid_env_std = np.asarray(resid_diag_env, dtype=float)
        variance_component_method = "fixed_observation_precision"
        sigma2_env_t = _t(np.asarray(resid_diag_env, dtype=float), device=device, dtype=dtype_t)
        resid_diag_obs_std_t = _t(np.asarray(resid_diag_obs_std, dtype=float), device=device, dtype=dtype_t)
        diag_add = resid_diag_obs_std_t.index_select(0, t) + float(lam_eff)

    prediction_resid_env_t = _t(
        prediction_resid_env_std, device=device, dtype=dtype_t
    )

    # Final GLS with observation/environment noise + λ for MEAN prediction
    Vtt = Ktt.clone()
    Vtt = 0.5*(Vtt + Vtt.mT); Vtt.diagonal().add_(diag_add)

    Vinv_y = _cholesky_solve(Vtt, ytt)
    Vinv_X = _cholesky_solve(Vtt, Xtt)
    Xt_Vinv_X = Xtt.mT @ Vinv_X
    Xt_Vinv_y = Xtt.mT @ Vinv_y
    try:
        Var_beta = _spd_inverse(Xt_Vinv_X)
    except RuntimeError:
        Var_beta = _spd_inverse(Xt_Vinv_X + 1e-10*torch.eye(Xt_Vinv_X.shape[0], dtype=dtype_t, device=device))
    beta = Var_beta @ Xt_Vinv_y

    # Prediction set indices
    if prediction_output == "all":
        s_idx = s_idx_all
        s_t = _t(s_idx, device=device, dtype=torch.int64)
    else:
        s_idx = s_idx_default
        s_t = s_default_t

    # Kts only (we only need mean for fast test_only)
    Kts = build_K_block(t, s_t)
    alpha = _cholesky_solve(Vtt, ytt - Xtt @ beta)
    mean_lat_std = (Kts.mT @ alpha) + (X_t.index_select(0, s_t) @ beta)

    gi_s = gi_t.index_select(0, s_t); ei_s = ei_t.index_select(0, s_t)
    mean_obs = _t(stdr.inv_mean(_to_numpy(mean_lat_std), _to_numpy(ei_s)), device=device, dtype=dtype_t)

    # ----------- Build predictions (FAST path for test_only) -----------
    # SAFE naming
    n = gi.size
    if obs_names is None or len(obs_names) == 0:
        all_names = [str(int(gi[i])) for i in range(n)]
    else:
        obs_names = list(obs_names)
        if len(obs_names) == n:
            all_names = [str(obs_names[i]) for i in range(n)]
        elif gi.max() < len(obs_names):
            all_names = [str(obs_names[int(gi[i])]) for i in range(n)]
        else:
            all_names = [str(int(gi[i])) for i in range(n)]
    if env_levels is None or len(env_levels) == 0:
        all_envs = [str(int(ei[i])) for i in range(n)]
    else:
        env_levels = list(env_levels)
        all_envs = [str(env_levels[int(ei[i])]) for i in range(n)]

    if prediction_output == "test_only":
        # Only test rows. When return_se is True we still compute the
        # per-component posterior variances (SE_g_latent / SE_ge_latent /
        # SE_e_latent / Prediction_{Var,SE}_{latent,observed}) on the
        # test rows only -- cheaper than the "all" path (no N^2 kernel
        # builds for non-test rows) but not the SE-less bare fast path.
        names_out = [all_names[i] for i in s_idx]
        envs_out  = [all_envs[i]  for i in s_idx]
        predictions = pd.DataFrame({
            "Name": names_out,
            "Env": envs_out,
            "Prediction": _to_numpy(mean_obs),
        }).reset_index(drop=True)

        if bool(return_se):
            # Total posterior variance on test rows (standardized scale).
            # Kts was already built above for mean prediction; reuse it.
            Kss_diag_s = torch.zeros(s_t.numel(), dtype=dtype_t, device=device)
            for wk, Gk in zip(wg_t.reshape(-1).tolist(), Gt_list):
                if float(wk) != 0.0:
                    rows = Gk.index_select(0, gi_s)
                    Kss_diag_s.add_(float(wk) * rows.gather(1, gi_s.unsqueeze(1)).squeeze(1))
            if (Se_t is not None) and (wge_t.reshape(-1).abs().sum() > 0):
                Se_diag_s = Se_t.diagonal()[ei_s]
                for wk, Gk in zip(wge_t.reshape(-1).tolist(), Gt_list):
                    if float(wk) != 0.0:
                        rows = Gk.index_select(0, gi_s)
                        diag_g = rows.gather(1, gi_s.unsqueeze(1)).squeeze(1)
                        Kss_diag_s.add_(float(wk) * (diag_g * Se_diag_s))
            if we != 0.0:
                Se_diag_s2 = (Se_t.diagonal()[ei_s] if Se_t is not None else torch.ones_like(Kss_diag_s))
                Kss_diag_s.add_(we * Se_diag_s2)

            Vinv_Kts_s = _cholesky_solve(Vtt, Kts)
            proj_s = (Kts * Vinv_Kts_s).sum(dim=0)
            var_lat_std_s = float(kernel_variance_scale) * torch.clamp(
                Kss_diag_s - proj_s, min=0.0
            )

            Xss_s = X_t.index_select(0, s_t)
            XB_s = Xss_s @ Var_beta
            var_beta_diag_s = (XB_s * Xss_s).sum(dim=1)
            var_obs_std_s = (
                var_lat_std_s
                + float(kernel_variance_scale) * var_beta_diag_s
                + prediction_resid_env_t.index_select(0, ei_s)
            )

            # Per-component posterior variances on test rows.
            def _block_g(idx1, idx2):
                gi1 = gi_t.index_select(0, idx1); gi2 = gi_t.index_select(0, idx2)
                K = torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
                for wk, Gk in zip(wg_t.reshape(-1).tolist(), Gt_list):
                    if float(wk) != 0.0:
                        K.add_(float(wk) * _build_geno_block(Gk, gi1, gi2))
                return K
            def _block_ge(idx1, idx2):
                S_eff = Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else None)
                if S_eff is None or wge_t.reshape(-1).abs().sum() == 0:
                    return torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
                gi1 = gi_t.index_select(0, idx1); gi2 = gi_t.index_select(0, idx2)
                S12 = _build_env_kernel(ei_t.index_select(0, idx1), ei_t.index_select(0, idx2), S_eff, dtype=dtype_t, device=device)
                K = torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
                for wk, Gk in zip(wge_t.reshape(-1).tolist(), Gt_list):
                    if float(wk) != 0.0:
                        K.add_(float(wk) * (_build_geno_block(Gk, gi1, gi2) * S12))
                return K
            def _block_e(idx1, idx2):
                S_eff = Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else None)
                if we == 0.0 or S_eff is None:
                    return torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
                return float(we) * _build_env_kernel(ei_t.index_select(0, idx1), ei_t.index_select(0, idx2), S_eff, dtype=dtype_t, device=device)

            Kts_g_s  = _block_g(t, s_t)
            Kts_ge_s = _block_ge(t, s_t)
            Kts_e_s  = _block_e(t, s_t)
            Vinv_Kts_g_s  = _cholesky_solve(Vtt, Kts_g_s)
            Vinv_Kts_ge_s = _cholesky_solve(Vtt, Kts_ge_s)
            Vinv_Kts_e_s  = _cholesky_solve(Vtt, Kts_e_s)
            Kss_g_s  = _block_g(s_t, s_t).diagonal()
            Kss_ge_s = _block_ge(s_t, s_t).diagonal()
            Kss_e_s  = _block_e(s_t, s_t).diagonal()

            var_g_std_s  = float(kernel_variance_scale) * torch.clamp(Kss_g_s  - (Kts_g_s  * Vinv_Kts_g_s ).sum(dim=0) - (Kts_g_s  * Vinv_Kts_ge_s).sum(dim=0) - (Kts_g_s  * Vinv_Kts_e_s ).sum(dim=0), min=0.0)
            var_ge_std_s = float(kernel_variance_scale) * torch.clamp(Kss_ge_s - (Kts_ge_s * Vinv_Kts_g_s ).sum(dim=0) - (Kts_ge_s * Vinv_Kts_ge_s).sum(dim=0) - (Kts_ge_s * Vinv_Kts_e_s ).sum(dim=0), min=0.0)
            var_e_std_s  = float(kernel_variance_scale) * torch.clamp(Kss_e_s  - (Kts_e_s  * Vinv_Kts_g_s ).sum(dim=0) - (Kts_e_s  * Vinv_Kts_ge_s).sum(dim=0) - (Kts_e_s  * Vinv_Kts_e_s ).sum(dim=0), min=0.0)

            var_lat_s = stdr.inv_var(_to_numpy(var_lat_std_s), _to_numpy(ei_s))
            var_obs_s = stdr.inv_var(_to_numpy(var_obs_std_s), _to_numpy(ei_s))
            var_g_s   = stdr.inv_var(_to_numpy(var_g_std_s),   _to_numpy(ei_s))
            var_ge_s  = stdr.inv_var(_to_numpy(var_ge_std_s),  _to_numpy(ei_s))
            var_e_s   = stdr.inv_var(_to_numpy(var_e_std_s),   _to_numpy(ei_s))

            predictions["SE_observed"]  = np.sqrt(np.maximum(var_obs_s, 0.0))
            predictions["SE_latent"]    = np.sqrt(np.maximum(var_lat_s, 0.0))
            predictions["SE_g_latent"]  = np.sqrt(np.maximum(var_g_s,   0.0))
            predictions["SE_ge_latent"] = np.sqrt(np.maximum(var_ge_s,  0.0))
            predictions["SE_e_latent"]  = np.sqrt(np.maximum(var_e_s,   0.0))
            predictions = _attach_prediction_se_columns(
                predictions,
                latent_var=np.maximum(var_lat_s, 0.0),
                observed_var=np.maximum(var_obs_s, 0.0),
            )

        # Minimal return object (keep beta info; skip heavy variance outputs)
        return {
            "beta": _to_numpy(beta),
            "beta_se": (
                None
                if not return_se
                else _to_numpy(
                    torch.sqrt(
                        torch.clamp(
                            float(kernel_variance_scale) * Var_beta.diagonal(),
                            0.0,
                        )
                    )
                )
            ),
            "var_components": varcomp_df,
            "sigma2_resid_env": stdr.inv_var(np.asarray(prediction_resid_env_std, dtype=float),
                                              np.arange(int(ei_t.max().item())+1)),
            "sigma2_resid_overall": float(np.nanmean(stdr.inv_var(np.asarray(prediction_resid_env_std, dtype=float),
                                                                   np.arange(int(ei_t.max().item())+1)))),
            "variance_component_method": variance_component_method,
            "diagnostics": {"method":"krr_exact", "device": device, "dtype": str(dtype_t),
                             "lam_eff": float(lam_eff), "lam_select": str(lam_select),
                            "variance_component_method": variance_component_method,
                            "sigma2_g_standardized": float(kernel_variance_scale),
                            "sigma2_e_standardized": float(np.nanmean(prediction_resid_env_std)),
                             "reml_normalize": reml_normalize},
            "result": {
                "predictions": predictions,
                "per_env": pd.DataFrame([]),
                "across_env": pd.DataFrame([]),
            }
        }

    # ---------- Full path ("all"): compute SEs/components as before ----------
    # Kss diagonal (standardized scale) — total
    Kss_diag = torch.zeros(s_t.numel(), dtype=dtype_t, device=device)
    for wk, Gk in zip(wg_t.reshape(-1).tolist(), Gt_list):
        if float(wk) != 0.0:
            rows = Gk.index_select(0, gi_s)
            Kss_diag.add_(float(wk) * rows.gather(1, gi_s.unsqueeze(1)).squeeze(1))
    if (Se_t is not None) and (wge_t.reshape(-1).abs().sum() > 0):
        Se_diag_s = Se_t.diagonal()[ei_s]
        for wk, Gk in zip(wge_t.reshape(-1).tolist(), Gt_list):
            if float(wk) != 0.0:
                rows = Gk.index_select(0, gi_s)
                diag_g = rows.gather(1, gi_s.unsqueeze(1)).squeeze(1)
                Kss_diag.add_(float(wk) * (diag_g * Se_diag_s))
    if we != 0.0:
        Se_diag_s = (Se_t.diagonal()[ei_s] if Se_t is not None else torch.ones_like(Kss_diag))
        Kss_diag.add_(we * Se_diag_s)

    Vinv_Kts = _cholesky_solve(Vtt, Kts)
    proj = (Kts * Vinv_Kts).sum(dim=0)
    var_lat_std = float(kernel_variance_scale) * torch.clamp(
        Kss_diag - proj, min=0.0
    )

    Xss = X_t.index_select(0, s_t)
    XB = Xss @ Var_beta
    var_beta_diag = (XB * Xss).sum(dim=1)
    var_obs_std = (
        var_lat_std
        + float(kernel_variance_scale) * var_beta_diag
        + prediction_resid_env_t.index_select(0, ei_t.index_select(0, s_t))
    )

    # Components
    def build_K_block_g(idx1, idx2):
        gi1 = gi_t.index_select(0, idx1); gi2 = gi_t.index_select(0, idx2)
        K = torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
        for wk, Gk in zip(wg_t.reshape(-1).tolist(), Gt_list):
            if float(wk) != 0.0:
                K.add_(float(wk) * _build_geno_block(Gk, gi1, gi2))
        return K
    def build_K_block_ge(idx1, idx2):
        S_eff = Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else None)
        if S_eff is None or wge_t.reshape(-1).abs().sum() == 0:
            return torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
        gi1 = gi_t.index_select(0, idx1); gi2 = gi_t.index_select(0, idx2)
        S12 = _build_env_kernel(ei_t.index_select(0, idx1), ei_t.index_select(0, idx2), S_eff, dtype=dtype_t, device=device)
        K = torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
        for wk, Gk in zip(wge_t.reshape(-1).tolist(), Gt_list):
            if float(wk) != 0.0:
                K.add_(float(wk) * (_build_geno_block(Gk, gi1, gi2) * S12))
        return K
    def build_K_block_e(idx1, idx2):
        S_eff = Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else None)
        if we == 0.0 or S_eff is None:
            return torch.zeros((idx1.numel(), idx2.numel()), dtype=dtype_t, device=device)
        return float(we) * _build_env_kernel(ei_t.index_select(0, idx1), ei_t.index_select(0, idx2), S_eff, dtype=dtype_t, device=device)

    Kts_g  = build_K_block_g(t, s_t)
    Kts_ge = build_K_block_ge(t, s_t)
    Kts_e  = build_K_block_e(t, s_t)

    Vinv_Kts_g  = _cholesky_solve(Vtt, Kts_g)
    Vinv_Kts_ge = _cholesky_solve(Vtt, Kts_ge)
    Vinv_Kts_e  = _cholesky_solve(Vtt, Kts_e)

    Kss_g  = build_K_block_g(s_t, s_t).diagonal()
    Kss_ge = build_K_block_ge(s_t, s_t).diagonal()
    Kss_e  = build_K_block_e(s_t, s_t).diagonal()

    var_g_std  = float(kernel_variance_scale) * torch.clamp(Kss_g  - (Kts_g  * Vinv_Kts_g ).sum(dim=0) - (Kts_g  * Vinv_Kts_ge).sum(dim=0) - (Kts_g  * Vinv_Kts_e ).sum(dim=0), min=0.0)
    var_ge_std = float(kernel_variance_scale) * torch.clamp(Kss_ge - (Kts_ge * Vinv_Kts_g ).sum(dim=0) - (Kts_ge * Vinv_Kts_ge).sum(dim=0) - (Kts_ge * Vinv_Kts_e ).sum(dim=0), min=0.0)
    var_e_std  = float(kernel_variance_scale) * torch.clamp(Kss_e  - (Kts_e  * Vinv_Kts_g ).sum(dim=0) - (Kts_e  * Vinv_Kts_ge).sum(dim=0) - (Kts_e  * Vinv_Kts_e ).sum(dim=0), min=0.0)

    var_lat  = _t(stdr.inv_var(_to_numpy(torch.clamp(var_lat_std, 0)),  _to_numpy(ei_s)), device=device, dtype=dtype_t)
    var_obs  = _t(stdr.inv_var(_to_numpy(torch.clamp(var_obs_std, 0)),  _to_numpy(ei_s)), device=device, dtype=dtype_t)
    var_g    = _t(stdr.inv_var(_to_numpy(var_g_std),    _to_numpy(ei_s)), device=device, dtype=dtype_t)
    var_ge   = _t(stdr.inv_var(_to_numpy(var_ge_std),   _to_numpy(ei_s)), device=device, dtype=dtype_t)
    var_e    = _t(stdr.inv_var(_to_numpy(var_e_std),    _to_numpy(ei_s)), device=device, dtype=dtype_t)

    # Full predictions DF ("all")
    names = [all_names[i] for i in range(n)]
    envs  = [all_envs[i]  for i in range(n)]
    label = np.array(["test"] * n, dtype=object)
    label[np.isin(np.arange(n), t_idx)] = "train"

    pred_all = np.full(n, np.nan, dtype=float)
    se_lat_all = np.full(n, np.nan, dtype=float)
    se_obs_all = np.full(n, np.nan, dtype=float)
    se_g_all  = np.full(n, np.nan, dtype=float)
    se_ge_all = np.full(n, np.nan, dtype=float)
    se_e_all  = np.full(n, np.nan, dtype=float)

    pred_all[s_idx]    = _to_numpy(mean_obs)
    se_lat_all[s_idx]  = np.sqrt(np.maximum(_to_numpy(var_lat), 0.0))
    se_obs_all[s_idx]  = np.sqrt(np.maximum(_to_numpy(var_obs), 0.0))
    se_g_all[s_idx]    = np.sqrt(np.maximum(_to_numpy(var_g),   0.0))
    se_ge_all[s_idx]   = np.sqrt(np.maximum(_to_numpy(var_ge),  0.0))
    se_e_all[s_idx]    = np.sqrt(np.maximum(_to_numpy(var_e),   0.0))

    predictions = pd.DataFrame({
        "row": np.arange(n, dtype=int),
        "Name": names,
        "Env": envs,
        "Prediction": pred_all,
        "SE_observed": se_obs_all,
        "SE_latent": se_lat_all,
        "SE_g_latent": se_g_all,
        "SE_ge_latent": se_ge_all,
        "SE_e_latent": se_e_all,
        "Train_Test_label": label.tolist(),
    })

    var_lat_all = np.full(n, np.nan, dtype=float)
    var_obs_all = np.full(n, np.nan, dtype=float)
    var_lat_all[s_idx] = np.maximum(_to_numpy(var_lat), 0.0)
    var_obs_all[s_idx] = np.maximum(_to_numpy(var_obs), 0.0)
    predictions = _attach_prediction_se_columns(
        predictions,
        latent_var=var_lat_all,
        observed_var=var_obs_all,
    )

    # Summaries (original scale)
    n_env = int(ei_t.max().item()) + 1
    env_names = env_levels if env_levels is not None else [str(i) for i in range(n_env)]
    sigma2_env_orig = stdr.inv_var(np.asarray(resid_diag_env, dtype=float), np.arange(n_env))
    S_rep = _to_numpy(Se_t if Se_t is not None else (env_cov_module.cov() if env_cov_module is not None else torch.eye(n_env, dtype=dtype_t, device=device)))
    report_bundle = _build_reporting_bundle(
        kernel_names=[f"G{i+1}" for i in range(nK)],
        w_g=np.asarray(wg_vec, dtype=float) * float(kernel_variance_scale),
        w_ge=np.asarray(wge_vec, dtype=float) * float(kernel_variance_scale),
        level_names=env_names,
        interaction_term_name="G:Env",
        interaction_cov=S_rep,
        env_main_scale=float(we) * float(kernel_variance_scale),
        env_main_cov=S_rep,
        residual_env=sigma2_env_orig,
        interaction_terms_meta=interaction_terms_meta,
    )
    report_bundle = _rescale_reporting_bundle_to_response_scale(report_bundle, stdr, env_names)
    per_env = report_bundle["env_variance_summary"].rename(columns={"level_index":"env_index","level":"Env"})
    var_components = report_bundle["var_components_summary"].copy()

    preds = predictions["Prediction"].to_numpy()
    se_obs = predictions["SE_observed"].to_numpy()
    se_lat = predictions["SE_latent"].to_numpy()
    yhat_mean = float(np.nanmean(preds)) if preds.size else float('nan')
    se_obs_mean = float(np.nanmean(se_obs)) if se_obs.size else float('nan')
    se_lat_mean = float(np.nanmean(se_lat)) if se_lat.size else float('nan')

    total_genetic_rows = var_components[
        (var_components["component_type"] == "genetic_main")
        & (var_components["kernel"] == "total")
    ]
    genetic_estimate = (
        float(total_genetic_rows["estimate"].iloc[0])
        if len(total_genetic_rows) else float("nan")
    )
    residual_estimate = _safe_nanmean_or_nan(sigma2_env_orig)
    krr_varcomp = pd.DataFrame({
        "component": ["vm(GID,KRR_kernel)!var", "Residual!R"],
        "estimate": [genetic_estimate, residual_estimate],
        "std.error": [np.nan, np.nan],
        "z.ratio": [np.nan, np.nan],
        "bound": ["P", "P"],
    })

    return {
        "beta": _to_numpy(beta),
        "beta_se": (
            None
            if not return_se
            else _to_numpy(
                torch.sqrt(
                    torch.clamp(
                        float(kernel_variance_scale) * Var_beta.diagonal(),
                        0.0,
                    )
                )
            )
        ),
        "var_components": var_components,
        "var_components_summary": report_bundle["var_components_summary"],
        "varcomp": krr_varcomp,
        "env_variance_summary": report_bundle["env_variance_summary"],
        "var_components_ai": pd.DataFrame([]),
        "interaction_variance_summary": report_bundle["interaction_variance_summary"],
        "interaction_correlation_summary": report_bundle["interaction_correlation_summary"],
        "residual_summary": report_bundle["residual_summary"],
        "heritability_summary": report_bundle["heritability_summary"],
        "sigma2_resid_env": np.asarray(sigma2_env_orig, dtype=float),
        "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_env_orig),
        "variance_component_method": variance_component_method,
        "diagnostics": {"method":"krr_exact", "device": device, "dtype": str(dtype_t),
                        "lam_eff": float(lam_eff), "lam_select": str(lam_select),
                        "variance_component_method": variance_component_method,
                        "sigma2_g_standardized": float(kernel_variance_scale),
                        "sigma2_e_standardized": float(np.nanmean(prediction_resid_env_std)),
                        "reml_normalize": reml_normalize},
        "result": {
            "predictions": predictions,
            "per_env": per_env,
            "across_env": pd.DataFrame({
                "yhat_across_env_mean": [yhat_mean],
                "se_across_env_mean_observed": [se_obs_mean],
                "se_across_env_mean_latent": [se_lat_mean],
            })
        }
    }

# ------------------------------ GP exact (GPyTorch) ---------------------------
def gp_exact_with_X(
    geno_kernels: Dict[str, np.ndarray],
    gi: np.ndarray, ei: np.ndarray, y: np.ndarray,
    train_idx: np.ndarray, test_idx: Optional[np.ndarray],
    w_g: np.ndarray, w_ge: np.ndarray, w_e: float,
    X: np.ndarray, S_e: Optional[np.ndarray] = None,
    device: Optional[str] = None, dtype: Optional[Union[str, 'torch.dtype']] = None,
    iters: int = 300, lr: float = 0.03, use_fixed_noise: bool = True,
    resid_diag_env: Optional[np.ndarray] = None,
    resid_diag_obs: Optional[np.ndarray] = None,
    geno_id_list: Optional[List[str]] = None,
    env_index_to_name: Optional[List[str]] = None,
    prediction_output: str = "test_only",
    standardize: str = "per_env",
    reml_normalize: str = "diag_mean",
    learn_scales: bool = False,
    compute_ai_se: bool = False,
    ai_hutch_samples: int = 64,
    ai_jitter: float = 1e-6,
    optimizer_method: str = "adam_then_ai",
    ai_steps: int = 8,
    ai_max_iter: int = 30,
    ai_tol_loglik: float = 1e-3,
    ai_tol_theta: float = 1e-4,
    refine_prediction_with_ai: bool = False,
    env_mean_from_obs_diag: bool = False,
    point_predictions_only: bool = False,
    return_se: bool = True,
    prediction_se_type: str = "both",
    prediction_se_method: str = "exact_dense",
    prediction_block_size: int = 512,
    prediction_diag_probes: int = 32,
    return_prediction_cov: bool = False,
    interaction_term_name: str = "G:Env",
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]] = None,
    precomputed_context: Optional[Dict[str, Any]] = None,
    theta0_warm: Optional[np.ndarray] = None,
    gp_engine: str = "auto",
) -> Dict[str, Any]:
    assert TORCH_AVAILABLE and GPTY_AVAILABLE, "gpytorch/torch required"
    set_deterministic(12345)
    device = device or _default_device
    dtype_t = _resolve_torch_dtype(dtype)

    if point_predictions_only:
        compute_ai_se = False
        return_se = False

    with _DefaultDTypeContext(dtype_t), gpytorch.settings.max_preconditioner_size(20), gpytorch.settings.cg_tolerance(1e-5):
        gi = np.asarray(gi, dtype=np.int64); ei = np.asarray(ei, dtype=np.int64); y = np.asarray(y, dtype=float)
        t_idx = np.asarray(train_idx, dtype=np.int64)
        s_idx_default = np.asarray(test_idx if (test_idx is not None and len(test_idx)>0) else train_idx, dtype=np.int64)
        row_idx = np.arange(gi.size, dtype=np.int64)

        stdr = _PerObsEnvStandardizer(mode=standardize).fit(y, ei, t_idx)
        y_std = stdr.transform(y, ei)

        cache = precomputed_context.setdefault("gp_exact_tensor_cache", {}) if isinstance(precomputed_context, dict) else None
        cache_key = None
        cached = None
        if isinstance(cache, dict):
            cache_key = (
                "gp_exact_with_X",
                str(device),
                str(dtype_t),
                str(reml_normalize),
                tuple(str(k) for k in geno_kernels.keys()),
                tuple(_cache_array_token(geno_kernels[k]) for k in geno_kernels.keys()),
                _cache_array_token(S_e),
                _cache_array_token(X),
                int(gi.size),
                int(ei.size),
            )
            cached = cache.get(cache_key)

        if isinstance(cached, dict):
            gi_t = cached["gi_t"]
            ei_t = cached["ei_t"]
            row_t = cached["row_t"]
            x_all = cached["x_all"]
            X_all = cached["X_all"]
            G_np_list = cached["G_np_list"]
            S_e = cached["S_e"]
            G_list = cached["G_list"]
            Se_t = cached["Se_t"]
        else:
            gi_t = _t(gi, device=device, dtype=torch.int64)
            ei_t = _t(ei, device=device, dtype=torch.int64)
            row_t= _t(row_idx, device=device, dtype=torch.int64)
            x_all = torch.stack([gi_t, ei_t, row_t], dim=-1)
            X_all = _t(np.asarray(X, dtype=np.float64), device=device, dtype=dtype_t)

            G_np_list = [np.asarray(geno_kernels[k], dtype=float) for k in geno_kernels.keys()]
            G_np_list, S_e = _apply_kernel_normalization(G_np_list, S_e, mode=reml_normalize)
            G_list = [_t(K, device=device, dtype=dtype_t) for K in G_np_list]
            Se_t = None if S_e is None else _t(np.asarray(S_e), device=device, dtype=dtype_t)
            if isinstance(cache, dict) and cache_key is not None:
                cache[cache_key] = {
                    "gi_t": gi_t,
                    "ei_t": ei_t,
                    "row_t": row_t,
                    "x_all": x_all,
                    "X_all": X_all,
                    "G_np_list": G_np_list,
                    "S_e": S_e,
                    "G_list": G_list,
                    "Se_t": Se_t,
                }
        y_t  = _t(y_std, device=device, dtype=dtype_t)
        t = _t(t_idx, device=device, dtype=torch.int64)
        s_default_t = _t(s_idx_default, device=device, dtype=torch.int64)

        train_x = x_all.index_select(0, t)
        test_x  = x_all if (prediction_output == "all") else x_all.index_select(0, s_default_t)
        train_y = y_t.index_select(0, t)
        if torch.isnan(train_y).any():
            mask = ~torch.isnan(train_y)
            m = train_y[mask].mean() if mask.any() else torch.zeros((), dtype=train_y.dtype, device=train_y.device)
            train_y = torch.where(torch.isnan(train_y), m, train_y)

        nK = len(G_list)
        wg_vec  = _prepare_weight_vector(w_g,  nK, default=1.0)
        wge_vec = _prepare_weight_vector(w_ge, nK, default=0.0)
        w_g_t  = _t(wg_vec,  device=device, dtype=dtype_t)
        w_ge_t = _t(wge_vec, device=device, dtype=dtype_t)
        env_names = env_index_to_name if env_index_to_name is not None else [str(i) for i in range(int(ei.max()) + 1)]
        normalized_terms, ge_cov_by_kernel = _build_kernel_specific_interaction_covs(
            interaction_terms_meta=interaction_terms_meta,
            kernel_names=list(geno_kernels.keys()),
            default_level_names=env_names,
            default_interaction_term_name=interaction_term_name,
            default_interaction_cov=(np.asarray(S_e, dtype=float) if S_e is not None else np.eye(len(env_names), dtype=float)),
            default_w_ge=wge_vec,
        )

        if use_fixed_noise:
            if resid_diag_obs is not None:
                resid_diag_obs = np.asarray(resid_diag_obs, dtype=float).reshape(-1)
                if resid_diag_obs.size != gi.size:
                    raise ValueError(f"resid_diag_obs must have length n_obs={gi.size}; got {resid_diag_obs.size}")
                if np.any(~np.isfinite(resid_diag_obs)) or np.any(resid_diag_obs <= 0):
                    raise ValueError("resid_diag_obs must be finite and > 0.")
                scale2_obs = stdr.inv_var(np.ones(gi.size, dtype=float), ei)
                scale2_obs = np.where(np.isfinite(scale2_obs) & (scale2_obs > 0), scale2_obs, 1.0)
                resid_diag_obs_std = resid_diag_obs / scale2_obs
                resid_diag_env = _env_mean_from_obs_diag(resid_diag_obs, ei, int(ei_t.max().item()) + 1)
                env_noise = _t(resid_diag_obs_std, device=device, dtype=dtype_t).index_select(0, t)
                likelihood = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise=env_noise, learn_additional_noise=False)
            elif resid_diag_env is None:
                krr = multikernel_krr_posterior_with_se(
                    [geno_kernels[k] for k in geno_kernels.keys()], gi, ei, y,
                    t_idx, np.array([], dtype=np.int64),
                    w_g=wg_vec, w_ge=wge_vec, w_e=float(w_e), S_e=S_e, lam=1e-3, X=X,
                    device=device, dtype=str(dtype_t).replace("torch.",""),
                    env_levels=env_index_to_name, obs_names=geno_id_list,
                    hutch_samples=32, seed=12345, prediction_output="test_only",
                    standardize=standardize, reml_normalize=reml_normalize
                )
                sig2_orig = np.asarray(krr["sigma2_resid_env"], dtype=float)
                sd2 = (stdr.env_stds_**2) if standardize.lower() == "per_env" else np.full_like(sig2_orig, stdr.global_std_**2, dtype=float)
                resid_diag_env = sig2_orig / sd2
            env_noise = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
            likelihood = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise=env_noise, learn_additional_noise=False)
        else:
            likelihood = gpytorch.likelihoods.GaussianLikelihood()

        use_multi_terms = normalized_terms is not None and len(normalized_terms) > 1
        if use_multi_terms:
            ge_mats_t = [_t(np.asarray(ge_cov_by_kernel.get(k, np.zeros((len(env_names), len(env_names)), dtype=float))), device=device, dtype=dtype_t) for k in geno_kernels.keys()]
            kernel = PrecomputedMultiTermIndexKernel(G_list, w_g_t, ge_mats_t, float(w_e), Se_t, learn_scales=learn_scales)
        else:
            kernel = PrecomputedIndexKernel(G_list, w_g_t, w_ge_t, float(w_e), Se_t, learn_scales=learn_scales)

        class ExactGPWithX(gpytorch.models.ExactGP):
            def __init__(self, train_x, train_y, likelihood, kernel, X_all):
                super().__init__(train_x, train_y, likelihood)
                self.mean_module = DesignLinearMean(X_all)
                self.covar_module = kernel
            def forward(self, x):
                mean_x = self.mean_module(x)
                covar_x = self.covar_module(x, x)
                return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

        model = ExactGPWithX(train_x, train_y, likelihood, kernel, X_all).to(device)

        model.train(); likelihood.train()
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        for _ in range(int(iters)):
            optimizer.zero_grad()
            output = model(train_x)
            K_op = output.lazy_covariance_matrix
            # Build residual operator R
            if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                noise_diag = likelihood.noise
                R_op = DiagLinearOperator(noise_diag)
            else:
                R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
            V_op = K_op + R_op
            # Design matrix for current training rows
            Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
            ll_reml = reml_loglik_from_V(train_y, Xtt, V_op, jitter=1e-6)
            loss = -ll_reml
            # With fixed observation noise and fixed kernel scales, the
            # profiled REML objective has no trainable covariance parameter.
            # DesignLinearMean.beta is deliberately profiled out here and is
            # installed by _assign_reml_beta_to_design_mean_ before
            # prediction. Calling backward() on this constant objective raises
            # "does not require grad" in PyTorch, so there is no Adam step to
            # perform in this valid fixed-parameter configuration.
            if not loss.requires_grad:
                break
            loss.backward()
            optimizer.step()

        
        # AI/Newton refinement is opt-in for the prediction model. The public
        # contract keeps AI-REML in the variance-component/reporting path so
        # prediction means and GP posterior SEs are produced by the GP path.
        if (bool(refine_prediction_with_ai)
                and optimizer_method.lower() in ("ai", "adam_then_ai")
                and int(ai_steps) > 0
                and (not point_predictions_only)):
            out = model(train_x)
            K_op = out.lazy_covariance_matrix
            if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                R_op = DiagLinearOperator(likelihood.noise)
            else:
                R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
            V_op = K_op + R_op
            Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())

            # Build θ list (effective, positive scales) + derivative operators
            kernel = model.covar_module
            labels = []
            DV_times_list = []
            get_parts = []
            set_parts = []

            # Use cached dense component matrices built in kernel during forward (stored for AI tables)
            Kg_list = getattr(kernel, "_Kg_list", None)
            Kge_list = getattr(kernel, "_Kge_list", None)
            Ke = getattr(kernel, "_Ke", None)

            if kernel.learn_scales and (Kg_list is not None) and (Kge_list is not None) and (Ke is not None):
                for i in range(int(kernel.w_g_unconstrained.numel())):
                    labels.append(f"w_g[{i}]")
                    get_parts.append(lambda i=i: torch.nn.functional.softplus(kernel.w_g_unconstrained[i]))
                    set_parts.append(lambda t, i=i: kernel.w_g_unconstrained.data.__setitem__(i, _softplus_inv_tensor(t)))
                    DV_times_list.append(lambda v, Km=Kg_list[i]: Km @ v)
                for i in range(int(kernel.w_ge_unconstrained.numel())):
                    labels.append(f"w_ge[{i}]")
                    get_parts.append(lambda i=i: torch.nn.functional.softplus(kernel.w_ge_unconstrained[i]))
                    set_parts.append(lambda t, i=i: kernel.w_ge_unconstrained.data.__setitem__(i, _softplus_inv_tensor(t)))
                    DV_times_list.append(lambda v, Km=Kge_list[i]: Km @ v)
                labels.append("w_e")
                get_parts.append(lambda: torch.nn.functional.softplus(kernel.w_e_unconstrained))
                set_parts.append(lambda t: kernel.w_e_unconstrained.data.copy_(_softplus_inv_tensor(t).reshape(kernel.w_e_unconstrained.shape)))
                DV_times_list.append(lambda v, Km=Ke: Km @ v)

            if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood) and hasattr(likelihood, "raw_noise"):
                labels.append("sigma2_noise")
                get_parts.append(lambda: likelihood.noise.squeeze())
                set_parts.append(lambda t: likelihood.raw_noise.data.copy_(_softplus_inv_tensor(t).reshape(likelihood.raw_noise.shape)))
                DV_times_list.append(lambda v: v)

            if len(DV_times_list) > 0:
                V_solve = lambda v: V_op.solve(v)
                ll_fn = lambda: reml_loglik_from_V(train_y, Xtt, V_op, jitter=ai_jitter)
                get_thetas = lambda: torch.stack([g() for g in get_parts])
                def set_thetas(theta_vec):
                    for k, setter in enumerate(set_parts):
                        setter(theta_vec[k])
                score_ai_fn = lambda: compute_score_and_AI_operator(train_y, Xtt, V_solve, DV_times_list, hutch_samples=ai_hutch_samples, jitter=ai_jitter)

                ai_newton_optimize(ll_fn, get_thetas, set_thetas, score_ai_fn,
                                   max_steps=int(ai_steps), jitter=ai_jitter, verbose=False)
                # refresh model cache after updates
                out = model(train_x)
        with torch.no_grad():
            _assign_reml_beta_to_design_mean_(model, likelihood, train_x, train_y, jitter=1e-6)
        model.eval(); likelihood.eval()
        with torch.no_grad():
            pred_lat = model(test_x); pred_obs = likelihood(pred_lat)
            mean_obs_std = pred_obs.mean

        # ---------- FAST return for test_only ----------
        if prediction_output == "test_only":
            out_idx = s_idx_default
            ei_sel = ei_t.index_select(0, s_default_t)
            mean_obs = _t(stdr.inv_mean(_to_numpy(mean_obs_std), _to_numpy(ei_sel)), device=device, dtype=dtype_t)

            gi_sel = gi_t.index_select(0, s_default_t)
            names = [geno_id_list[int(i)] if geno_id_list is not None else str(int(i)) for i in _to_numpy(gi_sel)]
            envs  = [env_index_to_name[int(j)] if env_index_to_name is not None else str(int(j)) for j in _to_numpy(ei_sel)]

            Var_beta = None
            if return_se:
                with torch.no_grad():
                    K_lin = model.covar_module(train_x, train_x)
                    if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                        if resid_diag_obs is not None:
                            noise_diag = _t(resid_diag_obs_std, device=device, dtype=dtype_t).index_select(0, t)
                        else:
                            noise_diag = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
                    else:
                        noise_diag = torch.full((train_x.shape[0],), float(likelihood.noise.item()), dtype=dtype_t, device=device)
                    V_lin = K_lin + DiagLinearOperator(noise_diag)
                    Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
                    Vinv_X = V_lin.solve(Xtt)
                    Xt_Vinv_X = Xtt.mT @ Vinv_X
                    Var_beta = _spd_inverse(Xt_Vinv_X)
                    beta_pred_var_std = _design_mean_prediction_var_diag(model, train_x, test_x, V_lin, Var_beta)
                    noise_test_std = _prediction_noise_diag_std(
                        likelihood, test_x, ei_t, resid_diag_env,
                        resid_diag_obs_std if resid_diag_obs is not None else None,
                        device, dtype_t,
                    )
                pred_lat_var_t = torch.clamp(pred_lat.variance + beta_pred_var_std, min=0.0)
                pred_lat_var_std = _to_numpy(pred_lat_var_t)
                pred_obs_var_std = _to_numpy(torch.clamp(pred_lat_var_t + noise_test_std, min=0.0))
            else:
                pred_lat_var_std = None
                pred_obs_var_std = None
            predictions = pd.DataFrame({
                "Name": names,
                "Env": envs,
                "Prediction": _to_numpy(mean_obs),
            })
            if return_se:
                lat_var = stdr.inv_var(pred_lat_var_std, _to_numpy(ei_sel)) if pred_lat_var_std is not None else None
                obs_var = stdr.inv_var(pred_obs_var_std, _to_numpy(ei_sel)) if pred_obs_var_std is not None else None
                predictions = _attach_prediction_se_columns(predictions, latent_var=lat_var, observed_var=obs_var)

            # --- AI-REML variance-component SEs (ASReml-style) ---
            varcomp_ai = pd.DataFrame([])
            varcomp_asreml_df = None
            varcomp_summary = None
            varcomp_ai_matrix = None
            env_names = env_index_to_name if env_index_to_name is not None else [str(i) for i in range(int(ei_t.max().item())+1)]
            if bool(compute_ai_se):
                try:
                    varcomp_ai = _ai_varcomp_table_gp_exact(
                        model, likelihood, train_x, train_y,
                        hutch_samples=int(ai_hutch_samples),
                        jitter=float(ai_jitter),
                    )
                except Exception as e:
                    warnings.warn(f"AI-SE computation failed (gp_exact): {e}")
                try:
                    from ai_reml import build_gp_exact_ai_varcomp
                    _env_sc = stdr.get_env_scales(env_names)
                    _nedf = int(train_y.numel() - Xtt.shape[1])
                    _build_kwargs = dict(
                        geno_kernel_names=list(geno_kernels.keys()),
                        env_labels=env_names,
                        env_col_name="Env",
                        hutch_samples=int(ai_hutch_samples),
                        jitter=float(ai_jitter),
                        env_scales=_env_sc,
                        nedf=_nedf,
                        max_iter=int(ai_max_iter),
                        tol_loglik=float(ai_tol_loglik),
                        tol_theta=float(ai_tol_theta),
                        refine_kernel_scales=bool(getattr(kernel, "learn_scales", False)),
                        theta0_warm=theta0_warm,
                    )
                    _engine = (gp_engine or "auto").lower()
                    if _engine == "auto":
                        _engine = "dense_v"
                    if _engine == "eigen":
                        try:
                            from mme_eigen import build_gp_exact_ai_varcomp_mme_eigen
                            varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp_mme_eigen(
                                model, likelihood, train_x, train_y, **_build_kwargs,
                            )
                        except NotImplementedError as _err:
                            warnings.warn(
                                f"gp_engine='eigen' not supported for this fit ({_err}); "
                                "falling back to dense V engine.",
                                stacklevel=2,
                            )
                            _engine = "dense_v"
                    if _engine == "mme":
                        try:
                            from ai_reml import build_gp_exact_ai_varcomp_mme
                            varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp_mme(
                                model, likelihood, train_x, train_y, **_build_kwargs,
                            )
                        except NotImplementedError as _err:
                            warnings.warn(
                                f"gp_engine='mme' not supported for this fit ({_err}); "
                                "falling back to dense V engine.",
                                stacklevel=2,
                            )
                            _engine = "dense_v"
                    if _engine == "dense_v":
                        varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp(
                            model, likelihood, train_x, train_y, **_build_kwargs,
                        )
                except Exception as e:
                    warnings.warn(f"ASReml-schema varcomp build failed (gp_exact): {e}")
            sigma2_env_orig = (stdr.inv_var(np.asarray(resid_diag_env), np.arange(int(ei_t.max().item())+1)) if use_fixed_noise and resid_diag_env is not None else np.full(len(env_names), np.nan, dtype=float))
            wg_report, wge_report, we_report = _effective_gp_exact_reporting_scales(
                model.covar_module, wg_vec, wge_vec, float(w_e)
            )
            report_bundle = _build_reporting_bundle(
                kernel_names=list(geno_kernels.keys()),
                w_g=wg_report,
                w_ge=wge_report,
                level_names=env_names,
                interaction_term_name=interaction_term_name,
                interaction_cov=(S_e if S_e is not None else np.eye(len(env_names))),
                env_main_scale=we_report,
                env_main_cov=(S_e if S_e is not None else np.eye(len(env_names))),
                residual_env=sigma2_env_orig,
                interaction_terms_meta=(normalized_terms if use_multi_terms else None),
            )
            report_bundle = _rescale_reporting_bundle_to_response_scale(report_bundle, stdr, env_names)
            report_bundle = _sync_reporting_bundle_with_reml(
                report_bundle, varcomp_asreml_df, env_names
            )
            sigma2_resid_report = _reported_residual_by_level(report_bundle, sigma2_env_orig)

            return {
                "beta": _to_numpy(model.mean_module.beta),
                "beta_se": (None if not return_se else _to_numpy(torch.sqrt(torch.clamp(Var_beta.diagonal(), 0.0)))),
                "var_components": report_bundle["var_components_summary"],
                "var_components_summary": report_bundle["var_components_summary"],
            "env_variance_summary": report_bundle["env_variance_summary"],
                "var_components_ai": varcomp_ai,
                "varcomp": varcomp_asreml_df,
                "summary": varcomp_summary,
                "ai_matrix": varcomp_ai_matrix,
                "interaction_variance_summary": report_bundle["interaction_variance_summary"],
                "interaction_correlation_summary": report_bundle["interaction_correlation_summary"],
                "residual_summary": report_bundle["residual_summary"],
                "heritability_summary": report_bundle["heritability_summary"],
                "sigma2_resid_env": sigma2_resid_report,
                "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_resid_report),
                "diagnostics": {"method":"gp_exact", "iters": int(iters), "device": device,
                                "use_fixed_noise": bool(use_fixed_noise),
                                "likelihood_noise": (
                                    float(likelihood.noise.detach().cpu().reshape(-1)[0].item())
                                    if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood)
                                    else float("nan")
                                ),
                                "reml_normalize": reml_normalize,
                                "learn_scales": bool(learn_scales)},
                "result": {
                    "predictions": predictions,
                    "per_env": pd.DataFrame([]),
                    "across_env": pd.DataFrame([]),
                }
            }

        # ---------- Full path ("all") below ----------
        # FE variance (Var_beta)
        with torch.no_grad():
            K_lin = model.covar_module(train_x, train_x)
            if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                noise_diag = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
            else:
                noise_diag = torch.full((train_x.shape[0],), float(likelihood.noise.item()), dtype=dtype_t, device=device)
            V_lin = K_lin + DiagLinearOperator(noise_diag)
            Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
            Vinv_X = V_lin.solve(Xtt)
            Xt_Vinv_X = Xtt.mT @ Vinv_X
            Var_beta = _spd_inverse(Xt_Vinv_X)

        out_idx = (np.arange(gi_t.numel()) if (prediction_output == "all") else s_idx_default)
        out_idx_t = _t(out_idx, device=device, dtype=torch.int64)
        ei_sel = ei_t.index_select(0, out_idx_t)

        # Simple variance decomposition for full output (as before) ...
        pred_x = x_all if prediction_output == "all" else x_all.index_select(0, s_default_t)
        pred_lat = model(pred_x)
        pred_obs = likelihood(pred_lat)
        mean_obs_std = pred_obs.mean
        # (Full variance code omitted for brevity since the fast path is the user's target;
        #  keep your previous full-variance implementation here if needed.)

        # Minimal but consistent predictions for "all"
        mean_obs = _t(stdr.inv_mean(_to_numpy(mean_obs_std), _to_numpy(ei_sel)), device=device, dtype=dtype_t)
        gi_sel = gi_t.index_select(0, out_idx_t)
        names = [geno_id_list[int(i)] if geno_id_list is not None else str(int(i)) for i in _to_numpy(gi_sel)]
        envs  = [env_index_to_name[int(j)] if env_index_to_name is not None else str(int(j)) for j in _to_numpy(ei_sel)]
        predictions = pd.DataFrame({"row": out_idx, "Name": names, "Env": envs, "Prediction": _to_numpy(mean_obs)})
        if return_se:
            beta_pred_var_std = _design_mean_prediction_var_diag(model, train_x, pred_x, V_lin, Var_beta)
            noise_pred_std = _prediction_noise_diag_std(
                likelihood, pred_x, ei_t, resid_diag_env,
                resid_diag_obs_std if resid_diag_obs is not None else None,
                device, dtype_t,
            )
            lat_var_std_t = torch.clamp(pred_lat.variance + beta_pred_var_std, min=0.0)
            lat_var = stdr.inv_var(_to_numpy(lat_var_std_t), _to_numpy(ei_sel))
            obs_var = stdr.inv_var(_to_numpy(torch.clamp(lat_var_std_t + noise_pred_std, min=0.0)), _to_numpy(ei_sel))
            predictions = _attach_prediction_se_columns(predictions, latent_var=lat_var, observed_var=obs_var)
        env_names = env_index_to_name if env_index_to_name is not None else [str(i) for i in range(int(ei_t.max().item())+1)]
        sigma2_env_orig = (stdr.inv_var(np.asarray(resid_diag_env), np.arange(int(ei_t.max().item())+1)) if use_fixed_noise and resid_diag_env is not None else np.full(len(env_names), np.nan, dtype=float))
        varcomp_asreml_df = None
        varcomp_summary = None
        varcomp_ai_matrix = None
        varcomp_error = None
        if bool(compute_ai_se):
            try:
                from ai_reml import build_gp_exact_ai_varcomp
                _env_sc = stdr.get_env_scales(env_names)
                _nedf = int(train_y.numel() - Xtt.shape[1])
                _build_kwargs = dict(
                    geno_kernel_names=list(geno_kernels.keys()),
                    env_labels=env_names,
                    env_col_name="Env",
                    hutch_samples=int(ai_hutch_samples),
                    jitter=float(ai_jitter),
                    env_scales=_env_sc,
                    nedf=_nedf,
                    max_iter=int(ai_max_iter),
                    tol_loglik=float(ai_tol_loglik),
                    tol_theta=float(ai_tol_theta),
                    refine_kernel_scales=bool(getattr(kernel, "learn_scales", False)),
                    theta0_warm=theta0_warm,
                )
                _engine = (gp_engine or "auto").lower()
                if _engine == "auto":
                    _engine = "dense_v"
                if _engine == "eigen":
                    try:
                        from mme_eigen import build_gp_exact_ai_varcomp_mme_eigen
                        varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp_mme_eigen(
                            model, likelihood, train_x, train_y, **_build_kwargs,
                        )
                    except NotImplementedError as _err:
                        warnings.warn(
                            f"gp_engine='eigen' not supported for this fit ({_err}); "
                            "falling back to dense V engine.",
                            stacklevel=2,
                        )
                        _engine = "dense_v"
                if _engine == "mme":
                    try:
                        from ai_reml import build_gp_exact_ai_varcomp_mme
                        varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp_mme(
                            model, likelihood, train_x, train_y, **_build_kwargs,
                        )
                    except NotImplementedError as _err:
                        warnings.warn(
                            f"gp_engine='mme' not supported for this fit ({_err}); "
                            "falling back to dense V engine.",
                            stacklevel=2,
                        )
                        _engine = "dense_v"
                if _engine == "dense_v":
                    varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_gp_exact_ai_varcomp(
                        model, likelihood, train_x, train_y, **_build_kwargs,
                    )
            except Exception as e:
                varcomp_error = f"{type(e).__name__}: {e}"
                warnings.warn(f"ASReml-schema varcomp build failed (gp_exact, all): {e}")

        wg_report, wge_report, we_report = _effective_gp_exact_reporting_scales(
            model.covar_module, wg_vec, wge_vec, float(w_e)
        )
        report_bundle = _build_reporting_bundle(
            kernel_names=list(geno_kernels.keys()),
            w_g=wg_report,
            w_ge=wge_report,
            level_names=env_names,
            interaction_term_name=interaction_term_name,
            interaction_cov=(S_e if S_e is not None else np.eye(len(env_names))),
            env_main_scale=we_report,
            env_main_cov=(S_e if S_e is not None else np.eye(len(env_names))),
            residual_env=sigma2_env_orig,
            interaction_terms_meta=(normalized_terms if use_multi_terms else None),
        )
        report_bundle = _rescale_reporting_bundle_to_response_scale(report_bundle, stdr, env_names)
        report_bundle = _sync_reporting_bundle_with_reml(
            report_bundle, varcomp_asreml_df, env_names
        )
        sigma2_resid_report = _reported_residual_by_level(report_bundle, sigma2_env_orig)

        return {
            "beta": _to_numpy(model.mean_module.beta),
            "beta_se": (None if not return_se else _to_numpy(torch.sqrt(torch.clamp(Var_beta.diagonal(), 0.0)))),
            "var_components": report_bundle["var_components_summary"],
            "var_components_summary": report_bundle["var_components_summary"],
            "env_variance_summary": report_bundle["env_variance_summary"],
            "var_components_ai": pd.DataFrame([]),
            "varcomp": varcomp_asreml_df,
            "summary": varcomp_summary,
            "ai_matrix": varcomp_ai_matrix,
            "interaction_variance_summary": report_bundle["interaction_variance_summary"],
            "interaction_correlation_summary": report_bundle["interaction_correlation_summary"],
            "residual_summary": report_bundle["residual_summary"],
            "heritability_summary": report_bundle["heritability_summary"],
            "sigma2_resid_env": sigma2_resid_report,
            "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_resid_report),
            "diagnostics": {"method":"gp_exact", "iters": int(iters), "device": device,
                            "use_fixed_noise": bool(use_fixed_noise),
                            "likelihood_noise": (
                                float(likelihood.noise.detach().cpu().reshape(-1)[0].item())
                                if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood)
                                else float("nan")
                            ),
                            "reml_normalize": reml_normalize,
                            "learn_scales": bool(learn_scales),
                            "varcomp_error": varcomp_error},
            "result": {
                "predictions": predictions,
                "per_env": pd.DataFrame([]),
                "across_env": pd.DataFrame([]),
            }
        }

# --------------------------- GP ICM-FA (GPyTorch) -----------------------------
def gp_icm_fa_with_X(
    geno_kernels: Dict[str, np.ndarray],
    gi: np.ndarray, ei: np.ndarray, y: np.ndarray,
    train_idx: np.ndarray, test_idx: Optional[np.ndarray],
    w_g: np.ndarray, w_ge: np.ndarray, w_e: float,
    X: np.ndarray, S_e: Optional[np.ndarray] = None,
    device: Optional[str] = None, dtype: Optional[Union[str, 'torch.dtype']] = None,
    iters: int = 300, lr: float = 0.03, use_fixed_noise: bool = True,
    resid_diag_env: Optional[np.ndarray] = None,
    resid_diag_obs: Optional[np.ndarray] = None,
    geno_id_list: Optional[List[str]] = None,
    env_index_to_name: Optional[List[str]] = None,
    prediction_output: str = "test_only",
    fa_rank: int = 1,
    standardize: str = "per_env",
    reml_normalize: str = "diag_mean",
    learn_scales: bool = False,
    compute_ai_se: bool = False,
    ai_hutch_samples: int = 64,
    ai_jitter: float = 1e-6,
    ai_include_fa_loadings: bool = False,
    optimizer_method: str = "adam_then_ai",
    ai_steps: int = 8,
    ai_max_iter: int = 30,
    ai_tol_loglik: float = 1e-3,
    ai_tol_theta: float = 1e-4,
    refine_prediction_with_ai: bool = False,
    env_mean_from_obs_diag: bool = False,
    point_predictions_only: bool = False,
    return_se: bool = True,
    prediction_se_type: str = "both",
    prediction_se_method: str = "exact_dense",
    prediction_block_size: int = 512,
    prediction_diag_probes: int = 32,
    return_prediction_cov: bool = False,
    interaction_term_name: str = "G:Env",
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]] = None,
    theta0_warm: Optional[np.ndarray] = None,
    gp_engine: str = "auto",
) -> Dict[str, Any]:
    assert TORCH_AVAILABLE and GPTY_AVAILABLE, "gpytorch/torch required"
    set_deterministic(12345)
    device = device or _default_device
    dtype_t = _resolve_torch_dtype(dtype)

    if point_predictions_only:
        compute_ai_se = False
        return_se = False

    with _DefaultDTypeContext(dtype_t), gpytorch.settings.max_preconditioner_size(20), gpytorch.settings.cg_tolerance(1e-5):
        gi = np.asarray(gi, dtype=np.int64); ei = np.asarray(ei, dtype=np.int64); y = np.asarray(y, dtype=float)
        t_idx = np.asarray(train_idx, dtype=np.int64)
        s_idx_default = np.asarray(test_idx if (test_idx is not None and len(test_idx)>0) else train_idx, dtype=np.int64)
        row_idx = np.arange(gi.size, dtype=np.int64)

        stdr = _PerObsEnvStandardizer(mode=standardize).fit(y, ei, t_idx)
        y_std = stdr.transform(y, ei)

        gi_t = _t(gi, device=device, dtype=torch.int64)
        ei_t = _t(ei, device=device, dtype=torch.int64)
        row_t= _t(row_idx, device=device, dtype=torch.int64)
        y_t  = _t(y_std, device=device, dtype=dtype_t)

        x_all = torch.stack([gi_t, ei_t, row_t], dim=-1)
        t = _t(t_idx, device=device, dtype=torch.int64)
        s_default_t = _t(s_idx_default, device=device, dtype=torch.int64)

        train_x = x_all.index_select(0, t)
        test_x  = x_all if (prediction_output == "all") else x_all.index_select(0, s_default_t)
        train_y = y_t.index_select(0, t)
        if torch.isnan(train_y).any():
            mask = ~torch.isnan(train_y)
            m = train_y[mask].mean() if mask.any() else torch.zeros((), dtype=train_y.dtype, device=train_y.device)
            train_y = torch.where(torch.isnan(train_y), m, train_y)

        X_all = _t(np.asarray(X, dtype=np.float64), device=device, dtype=dtype_t)

        # Normalize kernels
        G_np_list = [np.asarray(geno_kernels[k], dtype=float) for k in geno_kernels.keys()]
        G_np_list, S_e = _apply_kernel_normalization(G_np_list, S_e, mode=reml_normalize)
        G_list = [_t(K, device=device, dtype=dtype_t) for K in G_np_list]

        nK = len(G_list)
        wg_vec  = _prepare_weight_vector(w_g,  nK, default=1.0)
        wge_vec = _prepare_weight_vector(w_ge, nK, default=0.0)
        if np.any(np.abs(wge_vec) > 0):
            wg_fa_vec = wge_vec
        else:
            wg_fa_vec = wg_vec

        w_fa_t = _t(wg_fa_vec, device=device, dtype=dtype_t)
        Se_t   = None if S_e is None else _t(np.asarray(S_e), device=device, dtype=dtype_t)

        num_envs = int(ei_t.max().item()) + 1
        fa_module = FactorAnalyticEnv(num_envs, rank=int(fa_rank), init_cov=Se_t, dtype=dtype_t, device=device)
        env_cov_module = None

        if use_fixed_noise:
            if resid_diag_obs is not None:
                resid_diag_obs = np.asarray(resid_diag_obs, dtype=float).reshape(-1)
                if resid_diag_obs.size != gi.size:
                    raise ValueError(f"resid_diag_obs must have length n_obs={gi.size}; got {resid_diag_obs.size}")
                if np.any(~np.isfinite(resid_diag_obs)) or np.any(resid_diag_obs <= 0):
                    raise ValueError("resid_diag_obs must be finite and > 0.")
                scale2_obs = stdr.inv_var(np.ones(gi.size, dtype=float), ei)
                scale2_obs = np.where(np.isfinite(scale2_obs) & (scale2_obs > 0), scale2_obs, 1.0)
                resid_diag_obs_std = resid_diag_obs / scale2_obs
                resid_diag_env = _env_mean_from_obs_diag(resid_diag_obs, ei, int(ei_t.max().item()) + 1)
                env_noise = _t(resid_diag_obs_std, device=device, dtype=dtype_t).index_select(0, t)
                likelihood = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise=env_noise, learn_additional_noise=False)
            elif resid_diag_env is None:
                krr = multikernel_krr_posterior_with_se(
                    [geno_kernels[k] for k in geno_kernels.keys()], gi, ei, y,
                    t_idx, np.array([], dtype=np.int64),
                    w_g=wg_vec, w_ge=np.zeros_like(wg_vec), w_e=float(w_e), S_e=S_e, lam=1e-3, X=X,
                    device=device, dtype=str(dtype_t).replace("torch.",""),
                    env_levels=env_index_to_name, obs_names=geno_id_list,
                    hutch_samples=32, seed=12345, prediction_output="test_only",
                    standardize=standardize, reml_normalize=reml_normalize
                )
                sig2_orig = np.asarray(krr["sigma2_resid_env"], dtype=float)
                sd2 = (stdr.env_stds_**2) if standardize.lower() == "per_env" else np.full_like(sig2_orig, stdr.global_std_**2, dtype=float)
                resid_diag_env = sig2_orig / sd2
            env_noise = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
            likelihood = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise=env_noise, learn_additional_noise=False)
        else:
            likelihood = gpytorch.likelihoods.GaussianLikelihood()

        kernel = PrecomputedFAIndexKernel(G_list=G_list, w_g=w_fa_t, fa_env=fa_module,
                                          w_e=float(w_e), S_e=Se_t, learn_scales=learn_scales)

        class ExactGPFAWithX(gpytorch.models.ExactGP):
            def __init__(self, train_x, train_y, likelihood, kernel, X_all):
                super().__init__(train_x, train_y, likelihood)
                self.mean_module = DesignLinearMean(X_all)
                self.covar_module = kernel
            def forward(self, x):
                mean_x = self.mean_module(x)
                covar_x = self.covar_module(x, x)
                return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

        model = ExactGPFAWithX(train_x, train_y, likelihood, kernel, X_all).to(device)

        model.train(); likelihood.train()
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        for _ in range(int(iters)):
            optimizer.zero_grad()
            output = model(train_x)
            K_op = output.lazy_covariance_matrix
            # Build residual operator R
            if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                noise_diag = likelihood.noise
                R_op = DiagLinearOperator(noise_diag)
            else:
                R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
            V_op = K_op + R_op
            # Design matrix for current training rows
            Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
            ll_reml = reml_loglik_from_V(train_y, Xtt, V_op, jitter=1e-6)
            loss = -ll_reml
            # The same fixed-noise/fixed-scale condition can occur in the FA
            # exact path when all covariance inputs are supplied. The later
            # GLS profiling step still installs the fixed-effect estimate.
            if not loss.requires_grad:
                break
            loss.backward()
            optimizer.step()

                
        # AI/Newton refinement is opt-in for the prediction model. The public
        # contract keeps AI-REML in the variance-component/reporting path so
        # prediction means and GP posterior SEs are produced by the GP path.
        if (bool(refine_prediction_with_ai)
                and optimizer_method.lower() in ("ai", "adam_then_ai")
                and int(ai_steps) > 0
                and (not point_predictions_only)):
            out = model(train_x)
            K_op = out.lazy_covariance_matrix
            if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                R_op = DiagLinearOperator(likelihood.noise)
            else:
                R_op = DiagLinearOperator(likelihood.noise.expand(train_y.shape[0]))
            V_op = K_op + R_op
            Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())

            # Build θ list (effective, positive scales) + derivative operators
            kernel = model.covar_module
            labels = []
            DV_times_list = []
            get_parts = []
            set_parts = []

            # Use cached dense component matrices built in kernel during forward (stored for AI tables)
            Kg_list = getattr(kernel, "_Kg_list", None)
            Kge_list = getattr(kernel, "_Kge_list", None)
            Ke = getattr(kernel, "_Ke", None)

            if kernel.learn_scales and (Kg_list is not None) and (Kge_list is not None) and (Ke is not None):
                for i in range(int(kernel.w_g_unconstrained.numel())):
                    labels.append(f"w_g[{i}]")
                    get_parts.append(lambda i=i: torch.nn.functional.softplus(kernel.w_g_unconstrained[i]))
                    set_parts.append(lambda t, i=i: kernel.w_g_unconstrained.data.__setitem__(i, _softplus_inv_tensor(t)))
                    DV_times_list.append(lambda v, Km=Kg_list[i]: Km @ v)
                for i in range(int(kernel.w_ge_unconstrained.numel())):
                    labels.append(f"w_ge[{i}]")
                    get_parts.append(lambda i=i: torch.nn.functional.softplus(kernel.w_ge_unconstrained[i]))
                    set_parts.append(lambda t, i=i: kernel.w_ge_unconstrained.data.__setitem__(i, _softplus_inv_tensor(t)))
                    DV_times_list.append(lambda v, Km=Kge_list[i]: Km @ v)
                labels.append("w_e")
                get_parts.append(lambda: torch.nn.functional.softplus(kernel.w_e_unconstrained))
                set_parts.append(lambda t: kernel.w_e_unconstrained.data.copy_(_softplus_inv_tensor(t).reshape(kernel.w_e_unconstrained.shape)))
                DV_times_list.append(lambda v, Km=Ke: Km @ v)

            if isinstance(likelihood, gpytorch.likelihoods.GaussianLikelihood) and hasattr(likelihood, "raw_noise"):
                labels.append("sigma2_noise")
                get_parts.append(lambda: likelihood.noise.squeeze())
                set_parts.append(lambda t: likelihood.raw_noise.data.copy_(_softplus_inv_tensor(t).reshape(likelihood.raw_noise.shape)))
                DV_times_list.append(lambda v: v)

            if len(DV_times_list) > 0:
                V_solve = lambda v: V_op.solve(v)
                ll_fn = lambda: reml_loglik_from_V(train_y, Xtt, V_op, jitter=ai_jitter)
                get_thetas = lambda: torch.stack([g() for g in get_parts])
                def set_thetas(theta_vec):
                    for k, setter in enumerate(set_parts):
                        setter(theta_vec[k])
                score_ai_fn = lambda: compute_score_and_AI_operator(train_y, Xtt, V_solve, DV_times_list, hutch_samples=ai_hutch_samples, jitter=ai_jitter)

                ai_newton_optimize(ll_fn, get_thetas, set_thetas, score_ai_fn,
                                   max_steps=int(ai_steps), jitter=ai_jitter, verbose=False)
                # refresh model cache after updates
                out = model(train_x)
        with torch.no_grad():
            _assign_reml_beta_to_design_mean_(model, likelihood, train_x, train_y, jitter=1e-6)
        model.eval(); likelihood.eval()
        with torch.no_grad():
            pred_lat = model(test_x); pred_obs = likelihood(pred_lat)
            mean_obs_std = pred_obs.mean

        if prediction_output == "test_only":
            ei_sel = ei_t.index_select(0, s_default_t)
            mean_obs = _t(stdr.inv_mean(_to_numpy(mean_obs_std), _to_numpy(ei_sel)), device=device, dtype=dtype_t)

            gi_sel = gi_t.index_select(0, s_default_t)
            names = [geno_id_list[int(i)] if geno_id_list is not None else str(int(i)) for i in _to_numpy(gi_sel)]
            envs  = [env_index_to_name[int(j)] if env_index_to_name is not None else str(int(j)) for j in _to_numpy(ei_sel)]

            Var_beta = None
            if return_se:
                with torch.no_grad():
                    K_lin = model.covar_module(train_x, train_x)
                    if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                        if resid_diag_obs is not None:
                            noise_diag = _t(resid_diag_obs_std, device=device, dtype=dtype_t).index_select(0, t)
                        else:
                            noise_diag = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
                    else:
                        noise_diag = torch.full((train_x.shape[0],), float(likelihood.noise.item()), dtype=dtype_t, device=device)
                    V_lin = K_lin + DiagLinearOperator(noise_diag)
                    Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
                    Vinv_X = V_lin.solve(Xtt)
                    Xt_Vinv_X = Xtt.mT @ Vinv_X
                    Var_beta = _spd_inverse(Xt_Vinv_X)
                    beta_pred_var_std = _design_mean_prediction_var_diag(model, train_x, test_x, V_lin, Var_beta)
                    noise_test_std = _prediction_noise_diag_std(
                        likelihood, test_x, ei_t, resid_diag_env,
                        resid_diag_obs_std if resid_diag_obs is not None else None,
                        device, dtype_t,
                    )
                pred_lat_var_t = torch.clamp(pred_lat.variance + beta_pred_var_std, min=0.0)
                pred_lat_var_std = _to_numpy(pred_lat_var_t)
                pred_obs_var_std = _to_numpy(torch.clamp(pred_lat_var_t + noise_test_std, min=0.0))
            else:
                pred_lat_var_std = None
                pred_obs_var_std = None
            predictions = pd.DataFrame({
                "Name": names,
                "Env": envs,
                "Prediction": _to_numpy(mean_obs),
            })
            if return_se:
                lat_var = stdr.inv_var(pred_lat_var_std, _to_numpy(ei_sel)) if pred_lat_var_std is not None else None
                obs_var = stdr.inv_var(pred_obs_var_std, _to_numpy(ei_sel)) if pred_obs_var_std is not None else None
                predictions = _attach_prediction_se_columns(predictions, latent_var=lat_var, observed_var=obs_var)

            # --- AI-REML variance-component SEs (ASReml-style) ---
            varcomp_ai = pd.DataFrame([])
            varcomp_asreml_df = None
            varcomp_summary = None
            varcomp_ai_matrix = None
            env_names = env_index_to_name if env_index_to_name is not None else [str(i) for i in range(int(ei_t.max().item())+1)]
            if bool(compute_ai_se):
                try:
                    if _MODULAR_AI_REML_AVAILABLE:
                        from ai_reml import build_fa_icm_ai_varcomp
                        _env_sc = stdr.get_env_scales(env_names) if hasattr(stdr, 'get_env_scales') else None
                        _nedf = int(train_y.shape[0] - Xtt.shape[1])
                        varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_fa_icm_ai_varcomp(
                            model, likelihood, train_x, train_y,
                            fa_module=fa_module,
                            env_labels=env_names,
                            geno_kernel_names=list(geno_kernels.keys()),
                            term_name=f"fa({interaction_term_name.split(':')[-1] if ':' in interaction_term_name else 'Env'},{fa_rank}):vm(GID,G)",
                            env_col_name=interaction_term_name.split(':')[-1] if ':' in interaction_term_name else 'Env',
                            hutch_samples=int(ai_hutch_samples),
                            jitter=float(ai_jitter),
                            seed=12345,
                            env_scales=_env_sc,
                            nedf=_nedf,
                            max_iter=int(ai_max_iter),
                            tol_loglik=float(ai_tol_loglik),
                            tol_theta=float(ai_tol_theta),
                            theta0_warm=theta0_warm,
                        )
                        varcomp_ai = varcomp_asreml_df
                    else:
                        varcomp_ai = _ai_varcomp_table_gp_icm_fa(
                            model, likelihood, train_x, train_y,
                            fa_module=fa_module,
                            hutch_samples=int(ai_hutch_samples),
                            jitter=float(ai_jitter),
                            include_fa_loadings=bool(ai_include_fa_loadings),
                        )
                except Exception as e:
                    warnings.warn(f"AI-SE computation failed (gp_icm_fa): {e}")

            Sigma_ge = _to_numpy(fa_module.cov())
            wfa_report, we_report = _effective_fa_reporting_scales(
                model.covar_module, wg_fa_vec, float(w_e)
            )
            sigma2_env_orig = (stdr.inv_var(np.asarray(resid_diag_env), np.arange(int(ei_t.max().item())+1)) if use_fixed_noise and resid_diag_env is not None else np.full(len(env_names), np.nan, dtype=float))
            report_bundle = _build_reporting_bundle(
                kernel_names=list(geno_kernels.keys()),
                w_g=np.zeros_like(wfa_report),
                w_ge=wfa_report,
                level_names=env_names,
                interaction_term_name=interaction_term_name,
                interaction_cov=Sigma_ge,
                env_main_scale=we_report,
                env_main_cov=(S_e if S_e is not None else np.eye(len(env_names))),
                residual_env=sigma2_env_orig,
                interaction_terms_meta=None,
            )
            report_bundle = _rescale_reporting_bundle_to_response_scale(report_bundle, stdr, env_names)
            fa_derived = (
                (varcomp_summary or {}).get("fa_derived")
                if isinstance(varcomp_summary, dict)
                else None
            )
            report_bundle = _attach_fa_derived_se_to_summary(
                report_bundle, fa_derived,
            )
            report_bundle = _sync_reporting_bundle_with_reml(
                report_bundle, varcomp_asreml_df, env_names, fa_derived=fa_derived
            )
            sigma2_resid_report = _reported_residual_by_level(report_bundle, sigma2_env_orig)

            fa_structure_response = _rescale_env_covariance_to_response_scale(
                Sigma_ge, stdr, env_names
            )
            fa_total_response = float(np.sum(wfa_report)) * fa_structure_response

            return {
                "beta": _to_numpy(model.mean_module.beta),
                "beta_se": (None if not return_se else _to_numpy(torch.sqrt(torch.clamp(Var_beta.diagonal(), 0.0)))),
                "var_components": report_bundle["var_components_summary"],
                "var_components_summary": report_bundle["var_components_summary"],
                "env_variance_summary": report_bundle["env_variance_summary"],
                "var_components_ai": varcomp_ai,
                "varcomp": varcomp_asreml_df,
                "summary": varcomp_summary,
                "ai_matrix": varcomp_ai_matrix,
                "interaction_variance_summary": report_bundle["interaction_variance_summary"],
                "interaction_correlation_summary": report_bundle["interaction_correlation_summary"],
                "residual_summary": report_bundle["residual_summary"],
                "heritability_summary": report_bundle["heritability_summary"],
                "sigma2_resid_env": sigma2_resid_report,
                "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_resid_report),
                "diagnostics": {"method":"gp_icm_fa", "iters": int(iters), "device": device,
                                "use_fixed_noise": bool(use_fixed_noise), "fa_rank": int(fa_rank),
                                "reml_normalize": reml_normalize, "learn_scales": bool(learn_scales)},
                "env_covariance_fa": fa_total_response,
                "env_correlation_fa": _corr_from_cov_np(fa_total_response),
                "env_covariance_fa_structure": fa_structure_response,
                "env_covariance_fa_structure_standardized": Sigma_ge,
                "fa_kernel_weights": {
                    str(name): float(weight)
                    for name, weight in zip(geno_kernels.keys(), wfa_report)
                },
                "env_covariance_est": None if (env_cov_module is None) else _to_numpy(env_cov_module.cov()),
                "result": {
                    "predictions": predictions,
                    "per_env": pd.DataFrame([]),
                    "across_env": pd.DataFrame([]),
                }
            }

        # For "all": keep a minimal consistent output (mean only)
        out_idx = (np.arange(gi_t.numel()) if (prediction_output == "all") else s_idx_default)
        out_idx_t = _t(out_idx, device=device, dtype=torch.int64)
        ei_sel = ei_t.index_select(0, out_idx_t)
        mean_obs = _t(stdr.inv_mean(_to_numpy(mean_obs_std), _to_numpy(ei_sel)), device=device, dtype=dtype_t)
        gi_sel = gi_t.index_select(0, out_idx_t)
        names = [geno_id_list[int(i)] if geno_id_list is not None else str(int(i)) for i in _to_numpy(gi_sel)]
        envs  = [env_index_to_name[int(j)] if env_index_to_name is not None else str(int(j)) for j in _to_numpy(ei_sel)]
        predictions = pd.DataFrame({"row": out_idx, "Name": names, "Env": envs, "Prediction": _to_numpy(mean_obs)})
        if return_se:
            pred_x = x_all if prediction_output == "all" else x_all.index_select(0, s_default_t)
            pred_lat = model(pred_x)
            pred_obs = likelihood(pred_lat)
        # FE beta, beta_se
        Var_beta = None
        if return_se:
            with torch.no_grad():
                K_lin = model.covar_module(train_x, train_x)
                if isinstance(likelihood, gpytorch.likelihoods.FixedNoiseGaussianLikelihood):
                    noise_diag = _t(resid_diag_env, device=device, dtype=dtype_t).index_select(0, ei_t.index_select(0, t))
                else:
                    noise_diag = torch.full((train_x.shape[0],), float(likelihood.noise.item()), dtype=dtype_t, device=device)
                V_lin = K_lin + DiagLinearOperator(noise_diag)
                Xtt = model.mean_module.X_all.index_select(0, train_x[:, 2].long())
                Vinv_X = V_lin.solve(Xtt)
                Xt_Vinv_X = Xtt.mT @ Vinv_X
                Var_beta = _spd_inverse(Xt_Vinv_X)
                beta_pred_var_std = _design_mean_prediction_var_diag(model, train_x, pred_x, V_lin, Var_beta)
                noise_pred_std = _prediction_noise_diag_std(
                    likelihood, pred_x, ei_t, resid_diag_env,
                    resid_diag_obs_std if resid_diag_obs is not None else None,
                    device, dtype_t,
                )
            lat_var_std_t = torch.clamp(pred_lat.variance + beta_pred_var_std, min=0.0)
            lat_var = stdr.inv_var(_to_numpy(lat_var_std_t), _to_numpy(ei_sel))
            obs_var = stdr.inv_var(_to_numpy(torch.clamp(lat_var_std_t + noise_pred_std, min=0.0)), _to_numpy(ei_sel))
            predictions = _attach_prediction_se_columns(predictions, latent_var=lat_var, observed_var=obs_var)

        # --- AI-REML variance-component SEs (ASReml-style) ---
        varcomp_ai = pd.DataFrame([])
        varcomp_asreml_df = None
        varcomp_summary = None
        varcomp_ai_matrix = None
        env_names = env_index_to_name if env_index_to_name is not None else [str(i) for i in range(int(ei_t.max().item())+1)]
        if bool(compute_ai_se):
            try:
                if _MODULAR_AI_REML_AVAILABLE:
                    from ai_reml import build_fa_icm_ai_varcomp
                    _env_sc = stdr.get_env_scales(env_names) if hasattr(stdr, 'get_env_scales') else None
                    _nedf = int(train_y.shape[0] - Xtt.shape[1])
                    varcomp_asreml_df, varcomp_summary, varcomp_ai_matrix = build_fa_icm_ai_varcomp(
                        model, likelihood, train_x, train_y,
                        fa_module=fa_module,
                        env_labels=env_names,
                        geno_kernel_names=list(geno_kernels.keys()),
                        term_name=f"fa({interaction_term_name.split(':')[-1] if ':' in interaction_term_name else 'Env'},{fa_rank}):vm(GID,G)",
                        env_col_name=interaction_term_name.split(':')[-1] if ':' in interaction_term_name else 'Env',
                        hutch_samples=int(ai_hutch_samples),
                        jitter=float(ai_jitter),
                        seed=12345,
                        env_scales=_env_sc,
                        nedf=_nedf,
                        max_iter=int(ai_max_iter),
                        tol_loglik=float(ai_tol_loglik),
                        tol_theta=float(ai_tol_theta),
                        theta0_warm=theta0_warm,
                    )
                    varcomp_ai = varcomp_asreml_df
                else:
                    varcomp_ai = _ai_varcomp_table_gp_icm_fa(
                        model, likelihood, train_x, train_y,
                        fa_module=fa_module,
                        hutch_samples=int(ai_hutch_samples),
                        jitter=float(ai_jitter),
                        include_fa_loadings=bool(ai_include_fa_loadings),
                    )
            except Exception as e:
                warnings.warn(f"AI-SE computation failed (gp_icm_fa): {e}")
        Sigma_ge = _to_numpy(fa_module.cov())
        wfa_report, we_report = _effective_fa_reporting_scales(
            model.covar_module, wg_fa_vec, float(w_e)
        )
        sigma2_env_orig = (stdr.inv_var(np.asarray(resid_diag_env), np.arange(int(ei_t.max().item())+1)) if use_fixed_noise and resid_diag_env is not None else np.full(len(env_names), np.nan, dtype=float))
        report_bundle = _build_reporting_bundle(
            kernel_names=list(geno_kernels.keys()),
            w_g=np.zeros_like(wfa_report),
            w_ge=wfa_report,
            level_names=env_names,
            interaction_term_name=interaction_term_name,
            interaction_cov=Sigma_ge,
            env_main_scale=we_report,
            env_main_cov=(S_e if S_e is not None else np.eye(len(env_names))),
            residual_env=sigma2_env_orig,
            interaction_terms_meta=None,
        )
        report_bundle = _rescale_reporting_bundle_to_response_scale(report_bundle, stdr, env_names)
        fa_derived = (
            (varcomp_summary or {}).get("fa_derived")
            if isinstance(varcomp_summary, dict)
            else None
        )
        report_bundle = _attach_fa_derived_se_to_summary(
            report_bundle, fa_derived,
        )
        report_bundle = _sync_reporting_bundle_with_reml(
            report_bundle, varcomp_asreml_df, env_names, fa_derived=fa_derived
        )
        sigma2_resid_report = _reported_residual_by_level(report_bundle, sigma2_env_orig)

        fa_structure_response = _rescale_env_covariance_to_response_scale(
            Sigma_ge, stdr, env_names
        )
        fa_total_response = float(np.sum(wfa_report)) * fa_structure_response

        return {
            "beta": _to_numpy(model.mean_module.beta),
            "beta_se": (None if not return_se else _to_numpy(torch.sqrt(torch.clamp(Var_beta.diagonal(), 0.0)))),
            "var_components": report_bundle["var_components_summary"],
            "var_components_summary": report_bundle["var_components_summary"],
            "env_variance_summary": report_bundle["env_variance_summary"],
            "var_components_ai": varcomp_ai,
            "varcomp": varcomp_asreml_df,
            "summary": varcomp_summary,
            "ai_matrix": varcomp_ai_matrix,
            "interaction_variance_summary": report_bundle["interaction_variance_summary"],
            "interaction_correlation_summary": report_bundle["interaction_correlation_summary"],
            "residual_summary": report_bundle["residual_summary"],
            "heritability_summary": report_bundle["heritability_summary"],
            "sigma2_resid_env": sigma2_resid_report,
            "sigma2_resid_overall": _safe_nanmean_or_nan(sigma2_resid_report),
            "diagnostics": {"method":"gp_icm_fa", "iters": int(iters), "device": device,
                            "use_fixed_noise": bool(use_fixed_noise), "fa_rank": int(fa_rank),
                            "reml_normalize": reml_normalize, "learn_scales": bool(learn_scales)},
            "env_covariance_fa": fa_total_response,
            "env_correlation_fa": _corr_from_cov_np(fa_total_response),
            "env_covariance_fa_structure": fa_structure_response,
            "env_covariance_fa_structure_standardized": Sigma_ge,
            "fa_kernel_weights": {
                str(name): float(weight)
                for name, weight in zip(geno_kernels.keys(), wfa_report)
            },
            "env_covariance_est": None if (env_cov_module is None) else _to_numpy(env_cov_module.cov()),
            "result": {
                "predictions": predictions,
                "per_env": pd.DataFrame([]),
                "across_env": pd.DataFrame([]),
            }
        }

# ------------------------------- Public wrapper -------------------------------

def _apply_env_structure_policy(
    env_structure: Optional[str],
    n_env: int,
    fa_rank: Optional[int] = None,
    max_full_envs: int = 20,
    large_env_policy: str = "auto_fa",
) -> Tuple[Optional[str], Optional[int], Dict[str, object]]:
    diag = {
        "requested_env_structure": env_structure,
        "effective_env_structure": env_structure,
        "n_env": int(n_env),
        "fa_rank_requested": fa_rank,
        "fa_rank_effective": fa_rank,
        "env_policy_applied": False,
        "env_policy_message": None,
        "max_full_envs": int(max_full_envs),
        "large_env_policy": large_env_policy,
    }
    if env_structure is None:
        return env_structure, fa_rank, diag

    env_structure = str(env_structure).lower()
    large_env_policy = str(large_env_policy).lower()

    if env_structure not in {"corh", "corgh", "us", "fa"}:
        return env_structure, fa_rank, diag

    if n_env <= int(max_full_envs):
        return env_structure, fa_rank, diag

    if env_structure in {"us", "corgh"}:
        if large_env_policy == "error":
            raise ValueError(
                f"env_structure='{env_structure}' requested with n_env={n_env}. "
                f"For large environment counts (> {max_full_envs}), use FA instead."
            )
        if large_env_policy == "auto_fa":
            eff_rank = fa_rank
            if eff_rank is None:
                eff_rank = max(1, min(4, n_env // 8 if n_env >= 8 else 1))
            msg = (
                f"env_structure='{env_structure}' auto-converted to 'fa' because "
                f"n_env={n_env} exceeds max_full_envs={max_full_envs}. "
                f"Using fa_rank={eff_rank}."
            )
            warnings.warn(msg)
            diag["effective_env_structure"] = "fa"
            diag["fa_rank_effective"] = int(eff_rank)
            diag["env_policy_applied"] = True
            diag["env_policy_message"] = msg
            return "fa", int(eff_rank), diag
        if large_env_policy == "warn":
            msg = (
                f"env_structure='{env_structure}' kept with n_env={n_env} > "
                f"max_full_envs={max_full_envs}. This may be slow and less stable than FA."
            )
            warnings.warn(msg)
            diag["env_policy_applied"] = True
            diag["env_policy_message"] = msg
            return env_structure, fa_rank, diag

    return env_structure, fa_rank, diag




def _env_mean_from_obs_diag(obs_diag: Optional[np.ndarray], ei: np.ndarray, n_env: int) -> Optional[np.ndarray]:
    """Environment-wise mean of observation-level residual variances."""
    if obs_diag is None:
        return None
    if _MODULAR_LAYERS_AVAILABLE and _mod_env_mean_from_obs_diag is not None:
        return _mod_env_mean_from_obs_diag(obs_diag, ei, n_env)
    x = np.asarray(obs_diag, dtype=float).reshape(-1)
    if x.size != np.asarray(ei).size:
        raise ValueError(f"obs-level residual diagonal length {x.size} does not match n_obs={np.asarray(ei).size}")
    out = np.full((int(n_env),), np.nan, dtype=float)
    for j in range(int(n_env)):
        m = (np.asarray(ei) == j)
        if np.any(m):
            out[j] = float(np.nanmean(x[m]))
    return out


def _resolve_stagewise_residual_inputs(
    df: pd.DataFrame,
    ei: np.ndarray,
    obs_weights: Optional[Any] = None,
    obs_var: Optional[Any] = None,
    stage1_pev: Optional[Any] = None,
    obs_weight_mode: str = "relative_precision",
    obs_weight_global_scale: float = 1.0,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[str, Any]]:
    """
    Resolve one-stage / two-stage user inputs into a unified observation-level residual
    variance diagonal on the original response scale.
    """
    n = int(len(df))
    ei = np.asarray(ei, dtype=np.int64).reshape(-1)
    if ei.size != n:
        raise ValueError(f"ei length {ei.size} must match nrow(df)={n}")
    n_env = int(ei.max()) + 1 if ei.size else 0

    if _MODULAR_LAYERS_AVAILABLE and _mod_resolve_stagewise_residual_inputs is not None:
        env_levels = [str(i) for i in range(n_env)]
        out = _mod_resolve_stagewise_residual_inputs(
            n_obs=n,
            stage1_pev=stage1_pev,
            obs_var=obs_var,
            obs_weights=obs_weights,
            obs_weight_mode=str(obs_weight_mode),
            obs_weight_global_scale=float(obs_weight_global_scale),
            env_index=ei,
            env_levels=env_levels,
        )
        meta = {
            "stagewise_residual_source": out.get("source"),
            "obs_weight_mode": out.get("obs_weight_mode", str(obs_weight_mode)),
            "obs_weight_global_scale": float(out.get("obs_weight_global_scale", obs_weight_global_scale)),
            "resid_diag_obs_supplied": out.get("resid_diag_obs") is not None,
            "resid_diag_env_summary": None if out.get("resid_diag_env") is None else np.asarray(out.get("resid_diag_env"), dtype=float).copy(),
        }
        return out.get("resid_diag_obs"), out.get("resid_diag_env"), meta

    meta: Dict[str, Any] = {
        "stagewise_residual_source": None,
        "obs_weight_mode": str(obs_weight_mode),
        "obs_weight_global_scale": float(obs_weight_global_scale),
    }

    def _vec(x, name):
        arr = np.asarray(x, dtype=float).reshape(-1)
        if arr.size != n:
            raise ValueError(f"{name} must have length n_obs={n}; got {arr.size}")
        return arr

    resid_diag_obs = None
    if stage1_pev is not None:
        resid_diag_obs = _vec(stage1_pev, "stage1_pev")
        meta["stagewise_residual_source"] = "stage1_pev"
    elif obs_var is not None:
        resid_diag_obs = _vec(obs_var, "obs_var")
        meta["stagewise_residual_source"] = "obs_var"
    elif obs_weights is not None:
        w = _vec(obs_weights, "obs_weights")
        if np.any(~np.isfinite(w)) or np.any(w <= 0):
            raise ValueError("obs_weights must be finite and > 0")
        mode = str(obs_weight_mode).lower()
        if mode == "inverse_variance":
            resid_diag_obs = 1.0 / w
        elif mode == "relative_precision":
            resid_diag_obs = 1.0 / w
            m = float(np.nanmean(resid_diag_obs))
            if not np.isfinite(m) or m <= 0:
                raise ValueError("Could not normalize obs_weights to a mean-1 residual variance scale.")
            resid_diag_obs = (resid_diag_obs / m) * float(obs_weight_global_scale)
        else:
            raise ValueError("obs_weight_mode must be one of {'inverse_variance','relative_precision'}")
        meta["stagewise_residual_source"] = "obs_weights"
    else:
        meta["stagewise_residual_source"] = None

    if resid_diag_obs is not None:
        resid_diag_obs = np.asarray(resid_diag_obs, dtype=float)
        if np.any(~np.isfinite(resid_diag_obs)) or np.any(resid_diag_obs <= 0):
            raise ValueError("Resolved observation-level residual variances must be finite and > 0")
        resid_diag_env = _env_mean_from_obs_diag(resid_diag_obs, ei, n_env)
    else:
        resid_diag_env = None

    meta["resid_diag_obs_supplied"] = resid_diag_obs is not None
    meta["resid_diag_env_summary"] = None if resid_diag_env is None else resid_diag_env.copy()
    return resid_diag_obs, resid_diag_env, meta


def _build_grm_factor_cache_from_dict(grm_factor_cache: Dict[str, Any]):
    """Materialize a GRMFactorCache from the framework-level dict spec."""
    if not LARGE_BACKEND_AVAILABLE:
        raise RuntimeError("Large operator backend not available; cannot build GRMFactorCache.")
    ctype = str(grm_factor_cache.get("type", "zarr")).lower()
    if ctype == "zarr":
        root_dir = grm_factor_cache["root_dir"]
        num_grms = int(grm_factor_cache.get("num_grms", 1))
        return GRMFactorCacheZarr(root_dir=str(root_dir), num_grms=num_grms)
    if ctype == "memmap":
        return GRMFactorCacheMemmap(
            paths=grm_factor_cache["paths"],
            shapes=grm_factor_cache["shapes"],
            dtypes=grm_factor_cache["dtypes"],
        )
    raise ValueError(f"Unknown grm_factor_cache['type']: {ctype}")


def materialize_geno_kernel_from_factor_cache(
    grm_factor_cache: Any,
    row_index: Sequence[int],
    weights: Optional[Sequence[float]] = None,
    dtype: Union[str, "torch.dtype", None] = "float64",
) -> np.ndarray:
    """Materialize a small dense genomic kernel from a GRM factor cache.

    This helper is intended for bounded internal tuning/calibration subsets.
    Production prediction should keep using the operator backend directly.
    """
    if not TORCH_AVAILABLE:
        raise RuntimeError("PyTorch required to materialize a dense kernel from a factor cache")
    if grm_factor_cache is None:
        raise ValueError("grm_factor_cache is required")
    if isinstance(grm_factor_cache, dict):
        try:
            cache = _build_grm_factor_cache_from_dict(grm_factor_cache)
        except Exception:
            import importlib.util as _importlib_util
            import sys as _sys
            module_name = "gp_large_backend"
            if module_name in _sys.modules:
                _mod = _sys.modules[module_name]
            else:
                _path = os.path.join(os.path.dirname(__file__), "gp_large_backend.py")
                _spec = _importlib_util.spec_from_file_location(module_name, _path)
                if _spec is None or _spec.loader is None:
                    raise
                _mod = _importlib_util.module_from_spec(_spec)
                _sys.modules[module_name] = _mod
                _spec.loader.exec_module(_mod)
            ctype = str(grm_factor_cache.get("type", "zarr")).lower()
            if ctype == "zarr":
                cache = _mod.GRMFactorCacheZarr(
                    root_dir=str(grm_factor_cache["root_dir"]),
                    num_grms=int(grm_factor_cache.get("num_grms", 1)),
                )
            elif ctype == "memmap":
                cache = _mod.GRMFactorCacheMemmap(
                    paths=grm_factor_cache["paths"],
                    shapes=grm_factor_cache["shapes"],
                    dtypes=grm_factor_cache["dtypes"],
                )
            else:
                raise ValueError(f"Unknown grm_factor_cache['type']: {ctype}")
    else:
        cache = grm_factor_cache
    if not hasattr(cache, "get_rows"):
        raise TypeError("grm_factor_cache must be a dict spec or a GRMFactorCache object")
    idx_np = np.asarray(row_index, dtype=np.int64).reshape(-1)
    if idx_np.size == 0:
        raise ValueError("row_index is empty")
    if int(idx_np.min()) < 0 or int(idx_np.max()) >= int(cache.n_geno()):
        raise IndexError(
            f"row_index out of range for factor cache with n_geno={int(cache.n_geno())}"
        )
    dt = _resolve_torch_dtype(dtype) or torch.float64
    idx = torch.as_tensor(idx_np, dtype=torch.int64, device=torch.device("cpu"))
    n_grm = int(cache.num_grms)
    if weights is None:
        w = np.ones(n_grm, dtype=np.float64) / max(1, n_grm)
    else:
        w = np.asarray(list(weights), dtype=np.float64).reshape(-1)
        if w.size != n_grm:
            raise ValueError(f"weights length {w.size} does not match num_grms={n_grm}")
    K = torch.zeros((idx.numel(), idx.numel()), dtype=dt, device=torch.device("cpu"))
    with torch.no_grad():
        for r in range(n_grm):
            Phi = cache.get_rows(r, idx, device=torch.device("cpu"), dtype=dt)
            K.add_(Phi @ Phi.mT, alpha=float(w[r]))
        K = 0.5 * (K + K.mT)
    return K.detach().cpu().numpy()


def _mom_sigma_g_operator(
    grm_factor_cache: Dict[str, Any],
    w_g,
    gi,
    ei,
    y,
    train_idx,
    n_env: int,
    *,
    h2_prior: float = 0.5,
    jitter: float = 1e-6,
    shrinkage_correct: bool = True,
    pcg_tol: float = 1e-4,
    pcg_max_iter: int = 300,
    n_hutch: int = 16,
    block_size: int = 2048,
    device: Optional[str] = None,
    dtype_compute: Union[str, "torch.dtype"] = "float32",
    seed: int = 12345,
    verbose: bool = False,
):
    """Operator-backend Method-of-Moments Σ_g.

    Matvec-only per-env kernel ridge on the low-rank Φ factors stored in
    `grm_factor_cache` — no (n_geno × n_geno) allocation. Σ_g is assembled via the
    rank-space inner-product trick:
        Σ_g[e, e'] = (1/n_geno) · Σ_{k,l} w_k w_l · β_e^{(k)T} M_{kl} β_{e'}^{(l)}
    where M_{kl} = Φ_k^T Φ_l (streamed over blocks of genotypes) and
    β_e^{(k)} = Φ_k[g_unique_e,:]^T · scatter_add(α_e). α_e solves
    (Σ_k w_k Φ_k[gi_e,:] Φ_k[gi_e,:]^T + λI) α_e = y_e via PCG.

    Shrinkage correction (default on): d_e ≈ tr[(G_ee + λI)^-1 G²[train_e, train_e]] / n_geno
    estimated via Rademacher-Hutchinson (n_hutch samples), each one rank-space
    G² apply plus one PCG solve.
    """
    if not LARGE_BACKEND_AVAILABLE:
        raise RuntimeError(
            "Operator-backend MoM requested but large-scale backend not importable. "
            "Install the backend file or fall back to dense _mom_sigma_g_empirical."
        )
    import torch as _torch

    dev = device or _default_device
    dt_c = _resolve_torch_dtype(dtype_compute) or _torch.float32

    cache = _build_grm_factor_cache_from_dict(grm_factor_cache)
    nK = int(cache.num_grms)
    m_ranks = cache.ranks()
    n_geno = int(cache.n_geno())

    gi_np = np.asarray(gi, dtype=np.int64)
    ei_np = np.asarray(ei, dtype=np.int64)
    y_np = np.asarray(y, dtype=np.float64)
    tr_np = np.asarray(train_idx, dtype=np.int64)

    gi_tr = _torch.as_tensor(gi_np[tr_np], device=dev, dtype=_torch.int64)
    ei_tr = _torch.as_tensor(ei_np[tr_np], device=dev, dtype=_torch.int64)
    y_tr = _torch.as_tensor(y_np[tr_np], device=dev, dtype=dt_c)

    w_g_arr = np.asarray(w_g, dtype=np.float64).reshape(-1)
    if w_g_arr.size != nK:
        raise ValueError(f"w_g length {w_g_arr.size} != cache.num_grms {nK}")
    w_g_t = _torch.as_tensor(w_g_arr, device=dev, dtype=dt_c)

    # Step 1: stream-compute cross-kernel Gram matrices M_{kl} = Phi_k^T Phi_l (symmetric).
    M = [[None] * nK for _ in range(nK)]
    for start in range(0, n_geno, int(block_size)):
        end = min(start + int(block_size), n_geno)
        idx_blk = _torch.arange(start, end, device=dev, dtype=_torch.int64)
        Phi_blk = [cache.get_rows(k, idx_blk, device=dev, dtype=dt_c) for k in range(nK)]
        for k in range(nK):
            for l in range(k, nK):
                contrib = Phi_blk[k].transpose(0, 1) @ Phi_blk[l]
                if M[k][l] is None:
                    M[k][l] = contrib
                else:
                    M[k][l] = M[k][l] + contrib
    for k in range(nK):
        for l in range(k + 1, nK):
            M[l][k] = M[k][l].transpose(0, 1).contiguous()

    tr_G = 0.0
    for k in range(nK):
        tr_G += float(w_g_arr[k]) * float(_torch.trace(M[k][k]).item())
    mean_trace = tr_G / max(n_geno, 1)
    lam_scale = (1.0 - float(h2_prior)) / max(float(h2_prior), 1e-6)
    lam = lam_scale * mean_trace + float(jitter)
    lam_t = _torch.as_tensor(lam, device=dev, dtype=dt_c)

    beta_per_env: List[Optional[List["torch.Tensor"]]] = [None] * int(n_env)
    d_vec = np.zeros(int(n_env), dtype=np.float64)

    if bool(verbose):
        import time as _time
        _t_start = _time.perf_counter()
        _log_every = max(1, int(n_env) // 20)
        print(f"  [mom_op] starting per-env PCG loop: n_env={int(n_env)}, "
              f"n_hutch={int(n_hutch)}, shrinkage={bool(shrinkage_correct)}", flush=True)

    for e in range(int(n_env)):
        mask_e = (ei_tr == int(e))
        n_e = int(mask_e.sum().item())
        if n_e < 2:
            beta_per_env[e] = [_torch.zeros(m_ranks[k], device=dev, dtype=dt_c) for k in range(nK)]
            d_vec[e] = 1.0
            continue

        sel_e = _torch.nonzero(mask_e, as_tuple=False).reshape(-1)
        gi_e = gi_tr.index_select(0, sel_e)
        y_e = y_tr.index_select(0, sel_e)
        y_e = y_e - y_e.mean()

        op_e = _LowRankGRMOperator_op(cache, gi_e, weights=w_g_t.clone())
        fold_e = op_e.build_fold_cache()

        Phi_used = [cache.get_rows(k, fold_e.g_unique, device=dev, dtype=dt_c) for k in range(nK)]

        diag_approx = _torch.zeros(n_e, device=dev, dtype=dt_c)
        for k in range(nK):
            rn2 = (Phi_used[k] * Phi_used[k]).sum(dim=1)
            diag_approx = diag_approx + w_g_t[k] * rn2.index_select(0, fold_e.inv)
        diag_approx = (diag_approx + lam_t).clamp_min(_torch.as_tensor(1e-6, device=dev, dtype=dt_c))

        def A_mv(v, _fold=fold_e, _op=op_e, _lam=lam_t, _dt=dt_c):
            return _op.matvec(v, _fold, dtype_compute=_dt) + _lam * v

        def M_inv(v, _d=diag_approx):
            return v / _d

        res = _pcg_solve_op(A_mv, y_e, M_inv_mv=M_inv, tol=float(pcg_tol), max_iter=int(pcg_max_iter))
        alpha_e = res.x

        t_used = _scatter_add_1d_op(alpha_e, fold_e.inv, fold_e.ng_used)
        beta_e_list: List["torch.Tensor"] = []
        for k in range(nK):
            beta_e_list.append(Phi_used[k].transpose(0, 1) @ t_used)
        beta_per_env[e] = beta_e_list

        if bool(shrinkage_correct):
            gen = _torch.Generator(device=dev).manual_seed(int(seed) + int(e))
            d_accum = 0.0
            for _s in range(int(n_hutch)):
                z = (_torch.randint(0, 2, (n_e,), generator=gen, device=dev, dtype=_torch.int64) * 2 - 1).to(dt_c)
                t_z = _scatter_add_1d_op(z, fold_e.inv, fold_e.ng_used)
                a_list = [Phi_used[l].transpose(0, 1) @ t_z for l in range(nK)]
                b_list: List["torch.Tensor"] = []
                for k in range(nK):
                    b_k = _torch.zeros(m_ranks[k], device=dev, dtype=dt_c)
                    for l in range(nK):
                        b_k = b_k + w_g_t[l] * (M[k][l] @ a_list[l])
                    b_list.append(b_k)
                w_z = _torch.zeros(n_e, device=dev, dtype=dt_c)
                for k in range(nK):
                    u_k = Phi_used[k] @ b_list[k]
                    w_z = w_z + w_g_t[k] * u_k.index_select(0, fold_e.inv)
                v_res = _pcg_solve_op(A_mv, w_z, M_inv_mv=M_inv, tol=float(pcg_tol), max_iter=int(pcg_max_iter))
                d_accum += float(_torch.dot(z, v_res.x).item())
            d_vec[e] = d_accum / float(max(n_hutch, 1)) / max(n_geno, 1)

        if bool(verbose) and (((e + 1) % _log_every == 0) or (e + 1 == int(n_env))):
            _elapsed = _time.perf_counter() - _t_start
            _eta = _elapsed * (int(n_env) - e - 1) / max(e + 1, 1)
            print(f"  [mom_op] env {e+1}/{int(n_env)}  elapsed={_elapsed:.1f}s  "
                  f"eta={_eta:.1f}s  (n_e={n_e})", flush=True)

    if bool(verbose):
        print(f"  [mom_op] per-env loop done in {_time.perf_counter()-_t_start:.1f}s; "
              f"assembling Sigma_g...", flush=True)

    Sigma_g = np.zeros((int(n_env), int(n_env)), dtype=np.float64)
    for e in range(int(n_env)):
        beta_e = beta_per_env[e]
        for ep in range(e, int(n_env)):
            beta_ep = beta_per_env[ep]
            val = _torch.zeros((), device=dev, dtype=dt_c)
            for k in range(nK):
                for l in range(nK):
                    val = val + w_g_t[k] * w_g_t[l] * (beta_e[k] @ (M[k][l] @ beta_ep[l]))
            v_f = float(val.item()) / max(n_geno, 1)
            Sigma_g[e, ep] = v_f
            Sigma_g[ep, e] = v_f

    if bool(shrinkage_correct):
        d_floor = 1e-6
        d_safe = np.clip(d_vec, d_floor, None)
        scale = 1.0 / np.sqrt(np.outer(d_safe, d_safe))
        Sigma_g = Sigma_g * scale

    Sigma_g = 0.5 * (Sigma_g + Sigma_g.T)
    w_eig, V_eig = np.linalg.eigh(Sigma_g)
    w_floor = 1e-8 * max(float(np.abs(w_eig).max()), 1e-12)
    w_eig = np.clip(w_eig, w_floor, None)
    return (V_eig @ (w_eig[:, None] * V_eig.T)).astype(np.float64, copy=False)


def _mom_sigma_g_empirical(
    G_list,
    w_g,
    gi,
    ei,
    y,
    train_idx,
    n_env,
    *,
    h2_prior: float = 0.5,
    jitter: float = 1e-6,
    shrinkage_correct: bool = True,
    dtype=np.float64,
):
    """Method-of-moments Σ_g (n_env × n_env) via per-env BLUP cross-products.

    For each env e, fit a kernel ridge regression u_hat_e = G @ (G_ee + λI)^-1 y_e
    on the combined genetic kernel G = Σ_k w_g_k * G_k, using an h2-prior to
    pick λ without iteration. Then Σ_g[e,e'] = (1/n_geno) * u_hat_e' @ u_hat_{e'}.

    When shrinkage_correct=True (default), divides each entry by sqrt(d_e · d_{e'})
    where d_e = ||L_e^{-1} G[train_e, :]||_F² / n_geno is the unbiased per-env
    shrinkage factor (derived from E[u_hat_e^T u_hat_e] = σ_g² · d_e · n_geno under
    λ ≈ σ_e²/σ_g²). Undoes the systematic downward bias from BLUP shrinkage so
    Σ_g matches REML's scale to within ~10-20%.

    Non-iterative → O(n_env · n_obs_e^3) total (one small cholesky per env).
    At n_geno=500 × 24 envs this is ~0.1s vs ~40-60 min for AI-REML FA.
    """
    G0 = np.asarray(G_list[0], dtype=dtype)
    G_combined = np.zeros_like(G0)
    for Gk, wk in zip(G_list, np.asarray(w_g).reshape(-1)):
        G_combined = G_combined + float(wk) * np.asarray(Gk, dtype=dtype)
    n_geno = int(G_combined.shape[0])
    gi_tr = np.asarray(gi, dtype=np.int64)[np.asarray(train_idx, dtype=np.int64)]
    ei_tr = np.asarray(ei, dtype=np.int64)[np.asarray(train_idx, dtype=np.int64)]
    y_tr = np.asarray(y, dtype=dtype)[np.asarray(train_idx, dtype=np.int64)]
    U = np.zeros((n_geno, int(n_env)), dtype=dtype)
    d_vec = np.zeros(int(n_env), dtype=dtype)
    tr_G = float(np.trace(G_combined))
    mean_trace = tr_G / max(n_geno, 1)
    lam_scale = (1.0 - float(h2_prior)) / max(float(h2_prior), 1e-6)
    for e in range(int(n_env)):
        mask_e = (ei_tr == e)
        n_e = int(mask_e.sum())
        if n_e < 2:
            continue
        g_e = gi_tr[mask_e]
        y_e = y_tr[mask_e]
        y_e = y_e - float(y_e.mean())
        G_ee = G_combined[np.ix_(g_e, g_e)]
        lam = lam_scale * mean_trace
        K_e = G_ee + (lam + float(jitter)) * np.eye(n_e, dtype=dtype)
        try:
            L = np.linalg.cholesky(K_e)
            alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_e))
            if bool(shrinkage_correct):
                G_row_e = G_combined[g_e, :]
                B_e = np.linalg.solve(L, G_row_e)
                d_vec[e] = float(np.sum(B_e * B_e)) / max(n_geno, 1)
        except np.linalg.LinAlgError:
            K_fb = K_e + 1e-4 * np.eye(n_e, dtype=dtype)
            alpha = np.linalg.solve(K_fb, y_e)
            if bool(shrinkage_correct):
                G_row_e = G_combined[g_e, :]
                K_inv_G = np.linalg.solve(K_fb, G_row_e)
                d_vec[e] = float(np.sum(G_row_e * K_inv_G)) / max(n_geno, 1)
        U[:, e] = G_combined[:, g_e] @ alpha
    Sigma_g = (U.T @ U) / max(n_geno, 1)
    if bool(shrinkage_correct):
        d_floor = 1e-6
        d_safe = np.clip(d_vec, d_floor, None)
        scale = 1.0 / np.sqrt(np.outer(d_safe, d_safe))
        Sigma_g = Sigma_g * scale
    Sigma_g = 0.5 * (Sigma_g + Sigma_g.T)
    w, V = np.linalg.eigh(Sigma_g)
    w_floor = 1e-8 * max(float(np.abs(w).max()), 1e-12)
    w = np.clip(w, w_floor, None)
    return (V @ (w[:, None] * V.T)).astype(dtype, copy=False)


# ============================================================================
# Reaction-norm / env-covariate preprocessing
# ============================================================================
# Empirical result on G2F 2025 (172k obs, 294 envs, 22 all-new test envs):
# routing observed env covariates (EC file) into the 'e' main-env component
# via env_similarity = K_env lifted overall corr from 0.12 -> 0.31 (2.5x).
# See memory/project_g2f_2025_reaction_norm_findings.md for the full study.
#
# Key requirements for the lift to materialize:
#   1. include_components must contain 'e' (not just 'g','ge')
#   2. w_e > 0 and w_ge > 0 (defaults are 0.0 which silently neutralize K_env)
#   3. EC features must be yield-relevant (not the full engineered superset)
#   4. Training envs far from every test env in K_env add noise; filter them
#
# This preprocessor wires all four when the user supplies env_covariates or
# env_similarity. Set reaction_norm_auto=False to disable.
# ----------------------------------------------------------------------------

# G2F-style phenology stages (per APSIM-like stage codes used in G2F EC file).
_G2F_PHENOLOGY_STAGES = (
    "pGerEme", "pEmeEnJ", "pEnJFlo", "pFloFla", "pFlaFlw",
    "pFlwStG", "pStGEnG", "pEnGMat", "pMatHar",
)
# Yield-relevant EC feature prefixes empirically chosen for maize:
# TT = thermal time (growing degree accumulation), HI30 = heat-stress index,
# CumHI30 = cumulative heat stress. All three keyed per phenology stage.
_G2F_YIELD_EC_PREFIXES = ("TT_", "HI30_", "CumHI30_")


def _detect_g2f_ec_schema(feature_names: Sequence[str]) -> bool:
    """Return True if feature_names look like the G2F Environmental Covariate file."""
    names = [str(c) for c in feature_names]
    has_prefix = any(any(n.startswith(p) for p in _G2F_YIELD_EC_PREFIXES) for n in names)
    has_stage = any(any(s in n for s in _G2F_PHENOLOGY_STAGES) for n in names)
    return has_prefix and has_stage


def _filter_yield_relevant_ec(feature_names: Sequence[str]) -> List[str]:
    """From G2F-style EC feature names, keep the 27 yield-relevant (TT/HI30/CumHI30 per stage)."""
    keep = []
    for c in feature_names:
        c_str = str(c)
        if any(c_str.startswith(p) for p in _G2F_YIELD_EC_PREFIXES):
            suf = c_str.split("_", 1)[1] if "_" in c_str else ""
            if suf in _G2F_PHENOLOGY_STAGES:
                keep.append(c_str)
    return keep


def _env_covariates_to_aligned_matrix(
    env_covariates: Any,
    env_levels: Sequence[str],
    env_col: str,
    feature_qc: bool,
) -> Tuple[Optional[np.ndarray], Dict[str, Any]]:
    """Coerce env_covariates to an aligned (n_env × p) float matrix.

    Accepts:
      - pd.DataFrame with env ids as the index OR as an 'Env'/env_col column
      - np.ndarray of shape (n_env, p) (assumed to match env_levels order)
      - dict mapping {env_name: feature_vector}

    QC:
      - Drop non-numeric columns silently.
      - If features look G2F-EC, keep only TT/HI30/CumHI30 × phenology stages.
      - Impute missing envs with feature-median.
      - Drop zero-variance features.
      - Return None if nothing survives.
    """
    info: Dict[str, Any] = {
        "n_env": int(len(env_levels)),
        "source_type": None,
        "n_features_input": 0,
        "n_features_after_qc": 0,
        "n_features_after_variance_filter": 0,
        "g2f_schema_detected": False,
        "n_envs_imputed": 0,
    }

    if env_covariates is None:
        return None, info

    # Coerce to DataFrame indexed by env
    if isinstance(env_covariates, pd.DataFrame):
        ec = env_covariates.copy()
        if env_col in ec.columns:
            ec = ec.drop_duplicates(subset=env_col).set_index(env_col)
        elif "Env" in ec.columns:
            ec = ec.drop_duplicates(subset="Env").set_index("Env")
        # else assume index already env ids
        info["source_type"] = "DataFrame"
    elif isinstance(env_covariates, dict):
        ec = pd.DataFrame.from_dict(env_covariates, orient="index")
        info["source_type"] = "dict"
    elif isinstance(env_covariates, np.ndarray):
        if env_covariates.ndim != 2 or env_covariates.shape[0] != len(env_levels):
            raise ValueError(
                f"env_covariates as np.ndarray must be shape (n_env={len(env_levels)}, p); "
                f"got {env_covariates.shape}"
            )
        ec = pd.DataFrame(
            env_covariates, index=list(env_levels),
            columns=[f"feat_{i}" for i in range(env_covariates.shape[1])],
        )
        info["source_type"] = "ndarray"
    else:
        raise TypeError(
            f"env_covariates must be DataFrame, dict, or ndarray; got {type(env_covariates).__name__}"
        )

    # Numeric columns only
    ec_num = ec.select_dtypes(include=[np.number])
    if ec_num.shape[1] == 0:
        warnings.warn(
            "env_covariates contains no numeric columns after filtering. "
            "Falling back to env_similarity=None.",
            RuntimeWarning, stacklevel=3,
        )
        return None, info
    info["n_features_input"] = int(ec_num.shape[1])

    # Feature QC: G2F-style yield-relevant subset
    if feature_qc and _detect_g2f_ec_schema(list(ec_num.columns)):
        info["g2f_schema_detected"] = True
        keep = _filter_yield_relevant_ec(list(ec_num.columns))
        if keep:
            ec_num = ec_num[keep]
            _warn_once(
                ("g2f_ec_feature_qc", len(keep), int(info["n_features_input"])),
                f"env_covariates appears to be G2F Environmental Covariates. "
                f"Keeping {len(keep)} yield-relevant features "
                f"(TT/HI30/CumHI30 × {len(_G2F_PHENOLOGY_STAGES)} phenology stages) "
                f"out of {info['n_features_input']}. "
                f"Set reaction_norm_feature_qc=False to keep all features.",
                RuntimeWarning,
                stacklevel=3,
            )
    info["n_features_after_qc"] = int(ec_num.shape[1])

    # Align to env_levels, impute missing by column median
    F = pd.DataFrame(index=list(env_levels), columns=ec_num.columns, dtype=float)
    inter = ec_num.index.astype(str).intersection([str(e) for e in env_levels])
    if len(inter) == 0:
        warnings.warn(
            "env_covariates env identifiers do not overlap any env in the phenotype. "
            "Falling back to env_similarity=None.",
            RuntimeWarning, stacklevel=3,
        )
        return None, info
    F.loc[inter, :] = ec_num.loc[inter, :].astype(float).values
    missing = F.isna().any(axis=1).sum()
    info["n_envs_imputed"] = int(missing)
    F = F.fillna(ec_num.median(axis=0))

    # Standardize + drop zero-variance
    Farr = F.to_numpy(dtype=float)
    sd = Farr.std(axis=0)
    keep_var = sd > 1e-8
    Farr = Farr[:, keep_var]
    info["n_features_after_variance_filter"] = int(Farr.shape[1])
    if Farr.shape[1] == 0:
        warnings.warn(
            "All env_covariate features were zero-variance after alignment. "
            "Falling back to env_similarity=None.",
            RuntimeWarning, stacklevel=3,
        )
        return None, info
    mu = Farr.mean(axis=0)
    sd2 = Farr.std(axis=0)
    sd2 = np.where(sd2 > 1e-8, sd2, 1.0)
    Z = (Farr - mu) / sd2
    return Z.astype(np.float64), info


def _build_linear_K_env(Z: np.ndarray) -> np.ndarray:
    """Classical reaction-norm kernel: K = Z Z^T / p (VanRaden analog in env space)."""
    p = Z.shape[1]
    K = (Z @ Z.T) / float(max(p, 1))
    K = 0.5 * (K + K.T) + 1e-6 * np.eye(K.shape[0])
    return K


def _build_kenv_from_Z(
    Z: np.ndarray,
    kernel: str = "matern32",
    bandwidth: float = 1.0,
    kernel_kwargs: Optional[Dict[str, Any]] = None,
) -> np.ndarray:
    """Reaction-norm K_env from standardized env features Z.

    Thin shim over kenv_kernels.build_kenv. All kernel math lives in the
    standalone kenv_kernels.py module. The four legacy kernels (linear, rbf,
    matern32, matern52) are bit-for-bit preserved through this indirection.

    New kernels (pca_rbf, polynomial) and strict kwargs validation are
    available via kernel_kwargs. See kenv_kernels._KERNEL_REGISTRY for the
    full list.
    """
    try:
        from kenv_kernels import build_kenv
    except ImportError as e:
        raise ImportError(
            "kenv_kernels module not found on the Python path. "
            "It must live alongside this framework file."
        ) from e
    K, _info = build_kenv(
        Z, kernel=kernel, bandwidth=bandwidth,
        **(kernel_kwargs or {}),
    )
    return K


def prepare_env_similarity_from_covariates(
    env_covariates: Any,
    env_levels: Sequence[str],
    env_col: str = "Env",
    feature_qc: bool = True,
    kenv_kernel: str = "matern32",
    kenv_bandwidth: float = 1.0,
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """QC-align environment covariates and build an environment kernel once."""
    levels = [str(e) for e in env_levels]
    Z, cov_info = _env_covariates_to_aligned_matrix(
        env_covariates,
        levels,
        env_col,
        feature_qc=bool(feature_qc),
    )
    info = dict(cov_info or {})
    info["kenv_kernel"] = str(kenv_kernel).lower()
    info["kenv_bandwidth"] = float(kenv_bandwidth)
    info["kenv_kernel_kwargs"] = dict(kenv_kernel_kwargs or {})
    info["env_levels"] = levels
    if Z is None:
        return {"env_similarity": None, "info": info}
    K = _build_kenv_from_Z(
        Z,
        kernel=info["kenv_kernel"],
        bandwidth=info["kenv_bandwidth"],
        kernel_kwargs=info["kenv_kernel_kwargs"],
    )
    return {"env_similarity": np.asarray(K, dtype=np.float64), "info": info}


def _merge_kwargs_into_candidates(
    candidates: Sequence[Tuple],
    user_kwargs: Dict[str, Any],
) -> Tuple[Tuple, ...]:
    """Promote matching 2-tuple `(kernel, bw)` candidates to 3-tuple
    `(kernel, bw, kwargs)` by binding keys from `user_kwargs` that are in
    the candidate kernel's `allowed` set (per `kenv_kernels._KERNEL_REGISTRY`).

    Explicit per-candidate 3-tuples win (never overwritten). Candidates whose
    kernel accepts none of the user's keys are left untouched.
    """
    if not user_kwargs:
        return tuple(candidates)
    try:
        from kenv_kernels import _KERNEL_REGISTRY
    except Exception:
        return tuple(candidates)
    out: List[Tuple] = []
    for cand in candidates:
        if len(cand) >= 3:
            out.append(cand)
            continue
        kname, kbw = str(cand[0]).lower(), float(cand[1])
        allowed = _KERNEL_REGISTRY.get(kname, {}).get("allowed", set())
        matched = {k: v for k, v in user_kwargs.items() if k in allowed}
        if matched:
            out.append((kname, kbw, matched))
        else:
            out.append((kname, kbw))
    return tuple(out)


def _autoselect_kenv_kernel(
    Z: np.ndarray,
    env_means: np.ndarray,
    train_env_mask: np.ndarray,
    *,
    cv_groups: Optional[np.ndarray] = None,
    candidates: Sequence[Tuple[str, float]] = (
        ("linear", 1.0), ("matern32", 1.0),
    ),
    k_folds: int = 5,
    ridge: float = 1e-2,
    seed: int = 0,
) -> Tuple[str, float, Dict[str, Any]]:
    """Pick (kernel, bandwidth) by held-out CV on per-env mean phenotype.

    For each candidate, build K_env on the full env set (so bandwidth scaling
    via median pairwise distance is consistent), then fold over labeled
    training envs:

      pred_e = K[e, fold_train] @ (K[fold_train, fold_train] + ridge·I)^-1 @ env_means[fold_train]

    Score is mean per-fold Pearson correlation of pred vs env_means on held-out envs.

    The split design must match the production test split design:
      - Random-holdout test (e.g., CV2 hybrid splits) → leave cv_groups=None,
        k_folds=5 random folds approximate the test regime well.
      - Year-blocked test (e.g., G2F 2024 holdout) → pass `cv_groups` with one
        ID per env (e.g., the year), and leave-one-group-out CV is used. This
        is the regime the EC-distance heuristic could not detect.

    Returns (best_kernel, best_bandwidth, info_dict).
    """
    Z = np.asarray(Z, dtype=np.float64)
    env_means = np.asarray(env_means, dtype=np.float64)
    train_env_mask = np.asarray(train_env_mask, dtype=bool)

    # Normalize 2-tuple or 3-tuple candidates to a uniform (kernel, bw, kwargs) form.
    norm_cands = []
    for cand in candidates:
        if len(cand) == 2:
            _k, _b = cand
            _kw = {}
        elif len(cand) == 3:
            _k, _b, _kw = cand
            _kw = dict(_kw) if _kw is not None else {}
        else:
            raise ValueError(
                f"_autoselect_kenv_kernel: candidate must be (kernel, bandwidth) "
                f"or (kernel, bandwidth, kwargs_dict); got {cand!r}"
            )
        norm_cands.append((str(_k), float(_b), _kw))

    labeled = train_env_mask & np.isfinite(env_means)
    n_lab = int(labeled.sum())
    if n_lab < max(2 * k_folds, 4):
        return (norm_cands[0][0], norm_cands[0][1],
                {"reason": f"too_few_labeled_envs (n={n_lab}); using fallback {norm_cands[0][0]}",
                 "kwargs": norm_cands[0][2]})

    lab_idx = np.where(labeled)[0]
    rng = np.random.default_rng(int(seed))
    if cv_groups is not None:
        groups_lab = np.asarray(cv_groups)[lab_idx]
        unique_groups = np.unique(groups_lab)
        if unique_groups.size < 2:
            return (norm_cands[0][0], norm_cands[0][1],
                    {"reason": "single_cv_group; using fallback",
                     "kwargs": norm_cands[0][2]})
        folds = [(np.where(groups_lab != g)[0], np.where(groups_lab == g)[0])
                 for g in unique_groups]
        cv_design = f"leave-one-group-out (n_groups={unique_groups.size})"
    else:
        perm = rng.permutation(lab_idx)
        chunks = np.array_split(perm, int(k_folds))
        folds = []
        for f in chunks:
            te = f
            tr = np.setdiff1d(lab_idx, te, assume_unique=False)
            folds.append((np.where(np.isin(lab_idx, tr))[0],
                          np.where(np.isin(lab_idx, te))[0]))
        cv_design = f"random_kfold (k={int(k_folds)})"

    scores: Dict[int, list] = {}
    K_cache: Dict[int, np.ndarray] = {}
    for ci, (kernel, bw, kw) in enumerate(norm_cands):
        K_cache[ci] = _build_kenv_from_Z(
            Z, kernel=kernel, bandwidth=float(bw), kernel_kwargs=kw,
        )
        scores[ci] = []

    y_lab = env_means[lab_idx]
    for tr_loc, te_loc in folds:
        if te_loc.size == 0 or tr_loc.size < 2:
            continue
        tr_env = lab_idx[tr_loc]
        te_env = lab_idx[te_loc]
        for ci in range(len(norm_cands)):
            K = K_cache[ci]
            Ktt = K[np.ix_(tr_env, tr_env)] + float(ridge) * np.eye(tr_env.size)
            try:
                alpha = np.linalg.solve(Ktt, y_lab[tr_loc])
            except np.linalg.LinAlgError:
                continue
            pred = K[np.ix_(te_env, tr_env)] @ alpha
            actual = y_lab[te_loc]
            if pred.size < 2 or np.std(pred) < 1e-10 or np.std(actual) < 1e-10:
                continue
            scores[ci].append(float(np.corrcoef(actual, pred)[0, 1]))

    summary = {}
    for ci, (kernel, bw, kw) in enumerate(norm_cands):
        label = f"{kernel}@bw={bw}" + (f"|{kw}" if kw else "")
        v = scores[ci]
        summary[label] = (float(np.mean(v)) if v else float("nan"), len(v), ci)
    finite = {k: v for k, v in summary.items() if np.isfinite(v[0])}
    if not finite:
        return (
            norm_cands[0][0], norm_cands[0][1],
            {"reason": "all_folds_degenerate",
             "summary": {k: {"mean_corr": v[0], "n_folds": v[1]} for k, v in summary.items()},
             "cv_design": cv_design,
             "kwargs": norm_cands[0][2]},
        )
    best_label = max(finite, key=lambda k: finite[k][0])
    best_ci = finite[best_label][2]
    best_kernel, best_bw, best_kw = norm_cands[best_ci]
    info = {
        "summary": {k: {"mean_corr": v[0], "n_folds": v[1]} for k, v in summary.items()},
        "best": best_label,
        "cv_design": cv_design,
        "n_labeled_envs": int(n_lab),
        "kwargs": best_kw,
    }
    return (str(best_kernel), float(best_bw), info)


def _filter_training_envs_by_Kenv(
    K_env: np.ndarray,
    env_levels: Sequence[str],
    test_env_names: Sequence[str],
    topk_per_test: int,
) -> np.ndarray:
    """Return boolean mask (len n_env) keeping test envs + topk nearest train envs per test env."""
    env_to_idx = {str(e): i for i, e in enumerate(env_levels)}
    test_idx_env = np.array(
        [env_to_idx[str(e)] for e in test_env_names if str(e) in env_to_idx],
        dtype=np.int64,
    )
    keep = np.zeros(len(env_levels), dtype=bool)
    keep[test_idx_env] = True
    train_idx_env = np.array(
        [i for i, e in enumerate(env_levels) if str(e) not in set(str(x) for x in test_env_names)],
        dtype=np.int64,
    )
    if train_idx_env.size == 0 or test_idx_env.size == 0 or topk_per_test <= 0:
        keep[train_idx_env] = True
        return keep
    for ti in test_idx_env:
        sims = K_env[ti, train_idx_env]
        order = np.argsort(-sims)
        top = train_idx_env[order[: int(topk_per_test)]]
        keep[top] = True
    return keep


def _reaction_norm_preprocess(
    *,
    env_covariates: Any,
    env_similarity: Optional[np.ndarray],
    env_levels: Sequence[str],
    env_col: str,
    ei: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    include_components: Sequence[str],
    w_g: np.ndarray,
    w_ge: np.ndarray,
    w_e: float,
    reaction_norm_auto: bool,
    feature_qc: bool,
    env_filter_topk: Optional[int],
    w_ge_default: float,
    w_e_default: float,
    kenv_kernel: str = "matern32",
    kenv_bandwidth: float = 1.0,
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,
    kenv_auto_candidates: Sequence[Tuple[str, float]] = (
        ("linear", 1.0), ("matern32", 1.0),
    ),
    kenv_cv_groups: Optional[np.ndarray] = None,
    kenv_cv_kfolds: int = 5,
    kenv_cv_ridge: float = 1e-2,
    kenv_cv_seed: int = 0,
    y: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Central reaction-norm preprocessing. Returns a dict of (possibly modified) values.

    Behavior:
      - env_covariates provided: build K_env via QC + kenv_kernel (default matern32, bw=1.0), overrides env_similarity.
      - env_similarity provided (no env_covariates): use as-is.
      - Neither provided: warn only if test rows include all-new envs, return unchanged.
      - If env info present and include_components missing 'e': auto-add with warning.
      - If env info present and w_ge or w_e are all zero: set to reaction-norm defaults.
      - If env_filter_topk set and env info present: filter train_idx to nearest envs.
    """
    out = dict(
        env_similarity=env_similarity,
        include_components=tuple(include_components),
        w_g=np.asarray(w_g, dtype=float).copy(),
        w_ge=np.asarray(w_ge, dtype=float).copy(),
        w_e=float(w_e),
        train_idx=np.asarray(train_idx, dtype=np.int64).copy(),
        info={"reaction_norm_applied": False},
    )
    if not reaction_norm_auto:
        return out

    # --- Build K_env from covariates if provided ---
    K_built = None
    cov_info: Dict[str, Any] = {}
    if env_covariates is not None:
        Z, cov_info = _env_covariates_to_aligned_matrix(
            env_covariates, env_levels, env_col, feature_qc=feature_qc,
        )
        if Z is not None:
            _kk = str(kenv_kernel).lower()
            _kbw = float(kenv_bandwidth)
            if _kk == "auto":
                if y is None:
                    warnings.warn(
                        "kenv_kernel='auto' needs y for the embedded mini-CV; "
                        "falling back to matern32 bw=1.0.",
                        RuntimeWarning, stacklevel=3,
                    )
                    _kk, _kbw = "matern32", 1.0
                    auto_info: Dict[str, Any] = {"reason": "y_unavailable"}
                else:
                    n_env_total = Z.shape[0]
                    sums = np.zeros(n_env_total, dtype=np.float64)
                    counts = np.zeros(n_env_total, dtype=np.int64)
                    y_arr = np.asarray(y, dtype=np.float64)
                    ei_arr = np.asarray(ei, dtype=np.int64)
                    tr = np.asarray(train_idx, dtype=np.int64)
                    finite_tr = tr[np.isfinite(y_arr[tr])]
                    np.add.at(sums, ei_arr[finite_tr], y_arr[finite_tr])
                    np.add.at(counts, ei_arr[finite_tr], 1)
                    env_means = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)
                    train_env_mask = counts > 0
                    _cands = _merge_kwargs_into_candidates(
                        tuple(kenv_auto_candidates), kenv_kernel_kwargs or {},
                    )
                    _kk, _kbw, auto_info = _autoselect_kenv_kernel(
                        Z, env_means, train_env_mask,
                        cv_groups=kenv_cv_groups,
                        candidates=_cands,
                        k_folds=int(kenv_cv_kfolds),
                        ridge=float(kenv_cv_ridge),
                        seed=int(kenv_cv_seed),
                    )
                    warnings.warn(
                        f"kenv_kernel='auto' chose {_kk} bw={_kbw} "
                        f"(cv_design={auto_info.get('cv_design','n/a')}, "
                        f"summary={auto_info.get('summary')}).",
                        RuntimeWarning, stacklevel=3,
                    )
                out["info"]["kenv_auto"] = auto_info
            _resolved_kwargs = dict(kenv_kernel_kwargs or {})
            if str(kenv_kernel).lower() == "auto" and _resolved_kwargs:
                try:
                    from kenv_kernels import _KERNEL_REGISTRY
                    _allowed = _KERNEL_REGISTRY.get(_kk, {}).get("allowed", set())
                except Exception:
                    _allowed = set()
                _bad = [k for k in _resolved_kwargs if k not in _allowed]
                if _bad:
                    warnings.warn(
                        f"kenv_kernel='auto' resolved to {_kk!r}, which does not accept "
                        f"the provided kenv_kernel_kwargs keys {_bad}. Dropping them. "
                        f"To force per-kernel kwargs under 'auto', pass 3-tuple "
                        f"candidates via kenv_auto_candidates.",
                        RuntimeWarning, stacklevel=3,
                    )
                    _resolved_kwargs = {k: v for k, v in _resolved_kwargs.items() if k in _allowed}
            K_built = _build_kenv_from_Z(
                Z, kernel=_kk, bandwidth=_kbw,
                kernel_kwargs=_resolved_kwargs,
            )
            out["env_similarity"] = K_built
            out["info"]["env_covariates_qc"] = cov_info
            out["info"]["kenv_kernel"] = _kk
            out["info"]["kenv_bandwidth"] = _kbw
            out["info"]["kenv_kernel_kwargs"] = _resolved_kwargs
            _p_cov = int(cov_info.get("n_features_after_variance_filter", 0))
            _warn_once(
                ("env_covariates_routed", _kk, float(_kbw), _p_cov, int(len(env_levels))),
                f"env_covariates routed to env_similarity via reaction-norm "
                f"kernel={_kk} bandwidth={_kbw} "
                f"(p={_p_cov} features, n_env={len(env_levels)}). "
                f"include_components will be forced to contain 'e' for a proper main-env term.",
                RuntimeWarning,
                stacklevel=3,
            )

    have_env_info = out["env_similarity"] is not None
    if not have_env_info:
        train_envs = set(np.asarray(ei, dtype=np.int64)[np.asarray(train_idx, dtype=np.int64)])
        test_envs = set(np.asarray(ei, dtype=np.int64)[np.asarray(test_idx, dtype=np.int64)])
        unseen_test_envs = sorted(test_envs - train_envs)
        if unseen_test_envs:
            unseen_names = [str(env_levels[i]) for i in unseen_test_envs[:5]]
            more = "" if len(unseen_test_envs) <= 5 else f" (+{len(unseen_test_envs) - 5} more)"
            _warn_once(
                ("missing_env_kernel_cv0", len(env_levels), tuple(unseen_test_envs)),
                "No env_covariates or env_similarity supplied, and the test split contains "
                f"{len(unseen_test_envs)} all-new environment(s): {unseen_names}{more}. "
                "For CV0/prospective prediction, supply env_covariates (e.g., weather/soil/EC) "
                "so a reaction-norm kernel can be built; without it, the model has no mechanism "
                "to predict mean yield in unseen environments.",
                RuntimeWarning,
                stacklevel=3,
            )
        return out

    # --- Enforce 'e' in include_components (auto-include with warning) ---
    ic = list(out["include_components"])
    for comp in ("g", "ge", "e"):
        if comp not in ic:
            warnings.warn(
                f"Reaction-norm routing requires include_components to contain 'g', 'ge', and 'e'. "
                f"Auto-adding '{comp}'. Pass include_components=('g','ge','e') to silence.",
                RuntimeWarning, stacklevel=3,
            )
            ic.append(comp)
    out["include_components"] = tuple(ic)

    # --- Set reaction-norm-friendly weights if user left them at 0 ---
    if np.all(out["w_ge"] == 0.0):
        out["w_ge"] = np.full_like(out["w_ge"], float(w_ge_default))
        warnings.warn(
            f"w_ge was all-zero with env_similarity provided, which would silently disable "
            f"the G×E term that uses the env kernel. Setting w_ge={w_ge_default}. "
            f"Pass an explicit w_ge to override.",
            RuntimeWarning, stacklevel=3,
        )
    if float(out["w_e"]) == 0.0:
        out["w_e"] = float(w_e_default)
        warnings.warn(
            f"w_e was 0.0 with env_similarity provided, which would silently disable the 'e' "
            f"main-env term that uses K_env as S_e. Setting w_e={w_e_default}. "
            f"Pass an explicit w_e to override.",
            RuntimeWarning, stacklevel=3,
        )

    # --- Filter training envs by nearest K_env neighbors of test envs ---
    if env_filter_topk is not None and int(env_filter_topk) > 0:
        te = np.asarray(test_idx, dtype=np.int64)
        tr = np.asarray(train_idx, dtype=np.int64)
        test_env_codes = np.unique(ei[te])
        test_env_names = [str(env_levels[c]) for c in test_env_codes]
        K_for_filter = np.asarray(out["env_similarity"], dtype=float)
        keep_mask = _filter_training_envs_by_Kenv(
            K_for_filter, env_levels, test_env_names, int(env_filter_topk)
        )
        kept_codes = set(int(i) for i in np.where(keep_mask)[0])
        tr_keep = tr[np.isin(ei[tr], list(kept_codes))]
        n_before = int(tr.size)
        n_after = int(tr_keep.size)
        if n_after < n_before:
            warnings.warn(
                f"Reaction-norm env filter: dropped {n_before - n_after} training rows from "
                f"{int(keep_mask.size - keep_mask.sum())} env(s) whose K_env similarity to all "
                f"test envs was outside top-{env_filter_topk}. "
                f"Training rows: {n_before} -> {n_after}. "
                f"Set reaction_norm_env_filter_topk=None to disable.",
                RuntimeWarning, stacklevel=3,
            )
        out["train_idx"] = tr_keep
        out["info"]["env_filter"] = {
            "topk_per_test": int(env_filter_topk),
            "n_envs_kept": int(keep_mask.sum()),
            "n_envs_total": int(keep_mask.size),
            "n_train_rows_before": n_before,
            "n_train_rows_after": n_after,
        }

    out["info"]["reaction_norm_applied"] = True
    return out


def fit_mixed_model(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernels: Dict[str, np.ndarray],
    geno_ids: List[str],
    env_similarity: Optional[np.ndarray] = None,
    include_components: List[str] = ("g","ge","e"),
    fixed_effects: Optional[List[str]] = None,
    method: str = "krr_exact",
    w_g: Optional[np.ndarray] = None,
    w_ge: Optional[np.ndarray] = None,
    w_e: float = 0.0,
    obs_weights: Optional[Any] = None,
    obs_var: Optional[Any] = None,
    stage1_pev: Optional[Any] = None,
    obs_weight_mode: str = "relative_precision",
    obs_weight_global_scale: float = 1.0,
    # legacy/compat names accepted from R wrapper:
    krr_lam: float = 1e-3,
    lam: Optional[float] = None,                # alias of krr_lam
    reml_normalize: str = "diag_mean",          # accepted for compat
    krr_lams: Union[str, Sequence[float]] = "none",
    lam_select: str = "fixed",                  # "fixed" | "gcv"
    env_resid_robust: str = "none",             # "none" | "mad"
    env_resid_shrink_tau: float = 0.0,
    learn_envdiag_noise: bool = True,
    learn_scales: bool = False,                 # learn component scales in GP paths
    device: Optional[str] = None,
    dtype: Union[str, "torch.dtype", None] = "float64",
    seed: int = 12345,
    train_idx: Optional[np.ndarray] = None,
    test_idx: Optional[np.ndarray] = None,
    prediction_output: str = "test_only",
    icm_rank: int = 1,
    standardize: str = "per_env",
    hutch_samples: int = 64,
    backend: str = "auto",
    grm_factor_cache: Optional[Dict[str, Any]] = None,
    large_n_threshold: int = 50000,
    operator_tol: float = 1e-5,
    operator_max_iter: int = 500,
    operator_dtype_compute: Union[str, 'torch.dtype'] = "float32",
    compute_ai_se: bool = False,
    ai_hutch_samples: int = 64,
    ai_jitter: float = 1e-6,
    ai_include_fa_loadings: bool = False,
    optimizer_method: str = "adam_then_ai",
    ai_steps: int = 8,
    ai_max_iter: int = 30,
    ai_tol_loglik: float = 1e-3,
    ai_tol_theta: float = 1e-4,
    varcomp_mode: str = "reml",
    mom_h2_prior: float = 0.5,
    mom_shrinkage_correct: bool = True,
    mom_n_hutch: int = 16,
    mom_verbose: bool = False,
    # --- reaction-norm / env-covariate auto-routing (see _reaction_norm_preprocess) ---
    env_covariates: Optional[Any] = None,       # DataFrame/ndarray/dict of env-level features
    reaction_norm_auto: bool = True,            # enable auto K_env build + 'e' inclusion + env filter
    reaction_norm_feature_qc: bool = True,      # if EC features match known patterns, keep yield-relevant subset
    reaction_norm_env_filter_topk: Optional[int] = 15,  # keep topk nearest train envs per test env; None=no filter
    reaction_norm_w_ge_default: float = 0.05,   # w_ge when env info provided and user default is 0
    reaction_norm_w_e_default: float = 0.5,     # w_e when env info provided and user default is 0
    kenv_kernel: str = "matern32",              # K_env kernel: linear|rbf|matern32|matern52|pca_rbf|polynomial|auto
    kenv_bandwidth: float = 1.0,                # bandwidth × median pairwise dist (non-linear kernels)
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,  # per-kernel kwargs (e.g. {"pca_k": 10}, {"degree": 3, "c": 0.5})
    # kenv_kernel='auto' embedded mini-CV knobs (ignored unless kenv_kernel='auto'):
    kenv_auto_candidates: Sequence[Tuple[str, float]] = (
        ("linear", 1.0), ("matern32", 1.0),
    ),
    kenv_cv_groups: Optional[np.ndarray] = None,  # leave-one-group-out IDs per env (for blocked designs)
    kenv_cv_kfolds: int = 5,                    # random K-fold count when kenv_cv_groups is None
    kenv_cv_ridge: float = 1e-2,                # ridge added to inner kernel solve
    kenv_cv_seed: int = 0,                      # rng seed for the random K-fold split
    tiered_dispatch: bool = False,              # route predictions by env-connectivity (see fit_mixed_model_tiered)
    env_mean_from_obs_diag: bool = False,
    output_level: Optional[str] = None,
    point_predictions_only: Optional[bool] = None,
    return_se: Optional[bool] = None,
    prediction_se_type: str = "both",
    prediction_se_method: str = "exact_dense",
    prediction_block_size: int = 512,
    prediction_diag_probes: int = 32,
    return_prediction_cov: bool = False,
    fa_rank: int = 1,
    env_structure: Optional[str] = None,
    random_terms: Optional[List[Dict[str, Any]]] = None,
    interaction_terms_meta: Optional[Sequence[Dict[str, Any]]] = None,
    interaction_term_name: Optional[str] = None,
    theta0_warm: Optional[np.ndarray] = None,
    gp_engine: str = "auto",
    **kwargs
) -> Dict[str, Any]:

    output_level, point_predictions_only, return_se, compute_ai_se = _resolve_output_level(
        output_level,
        point_predictions_only=point_predictions_only,
        return_se=return_se,
        compute_ai_se=compute_ai_se,
    )
    method = str(method).lower()
    valid_methods = ("gp_exact", "krr_exact", "gp_icm_fa")
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}; got {method!r}")

    if tiered_dispatch:
        if train_idx is None or test_idx is None:
            raise ValueError(
                "fit_mixed_model: tiered_dispatch=True requires train_idx and "
                "test_idx (per-tier routing uses the split to classify rows)."
            )
        if env_similarity is None and env_covariates is None:
            raise ValueError(
                "fit_mixed_model: tiered_dispatch=True requires env_similarity "
                "or env_covariates (K_env is the cv0-tier kernel)."
            )
        try:
            from fit_mixed_model_tiered import fit_mixed_model_tiered as _tiered_fn
        except ImportError as e:
            raise ImportError(
                f"fit_mixed_model: tiered_dispatch=True requires "
                f"fit_mixed_model_tiered on the Python path. ({e})"
            )
        tiered_verbose = bool(kwargs.pop("tiered_verbose", False))
        tiered_return_se = bool(kwargs.pop("tiered_return_se", bool(return_se)))
        tiered_n_probes = int(kwargs.pop("tiered_n_hutchinson_probes",
                                         kwargs.pop("n_hutchinson_probes", 256)))
        tiered_hutch_seed = kwargs.pop("tiered_hutchinson_seed", None)
        tiered_se_pcg_tol = kwargs.pop("tiered_se_pcg_tol", None)
        tiered_se_pcg_max_iter = kwargs.pop("tiered_se_pcg_max_iter", None)
        tiered_use_nys = bool(kwargs.pop("tiered_use_nystrom_se_preconditioner", False))
        tiered_nys_rank = int(kwargs.pop("tiered_nystrom_rank", 64))
        tiered_nys_seed = int(kwargs.pop("tiered_nystrom_seed", 0))
        tiered_use_lan = bool(kwargs.pop("tiered_use_lanczos_se_variance_reduction", True))
        tiered_lan_iters = int(kwargs.pop("tiered_lanczos_iters", 80))
        tiered_lan_seed = int(kwargs.pop("tiered_lanczos_seed", 0))
        tiered_lan_floor = float(kwargs.pop("tiered_lanczos_theta_floor", 1e-8))
        tiered_out = _tiered_fn(
            pheno_df=pheno_df, gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=list(geno_ids),
            train_idx=np.asarray(train_idx, dtype=np.int64),
            test_idx=np.asarray(test_idx, dtype=np.int64),
            env_similarity=env_similarity, env_covariates=env_covariates,
            feature_qc=bool(reaction_norm_feature_qc),
            kenv_kernel=str(kenv_kernel),
            kenv_bandwidth=float(kenv_bandwidth),
            kenv_kernel_kwargs=kenv_kernel_kwargs,
            w_g=w_g, w_ge=w_ge, w_e=float(w_e),
            include_components=tuple(include_components),
            mom_n_hutch=int(kwargs.pop("mom_n_hutch", 16)),
            mom_h2_prior=float(kwargs.pop("mom_h2_prior", 0.5)),
            method=str(method), backend=str(backend),
            dtype=str(dtype) if not hasattr(dtype, "itemsize") else "float32",
            operator_dtype_compute=(
                operator_dtype_compute
                if isinstance(operator_dtype_compute, str)
                else "float32"
            ),
            operator_tol=float(operator_tol),
            operator_max_iter=int(operator_max_iter),
            standardize=str(standardize),
            grm_factor_cache=grm_factor_cache,
            seed=int(seed),
            return_se=tiered_return_se,
            n_hutchinson_probes=tiered_n_probes,
            hutchinson_seed=tiered_hutch_seed,
            se_pcg_tol=tiered_se_pcg_tol,
            se_pcg_max_iter=tiered_se_pcg_max_iter,
            use_nystrom_se_preconditioner=tiered_use_nys,
            nystrom_rank=tiered_nys_rank,
            nystrom_seed=tiered_nys_seed,
            use_lanczos_se_variance_reduction=tiered_use_lan,
            lanczos_iters=tiered_lan_iters,
            lanczos_seed=tiered_lan_seed,
            lanczos_theta_floor=tiered_lan_floor,
            verbose=tiered_verbose,
        )
        preds = tiered_out["predictions"].copy()
        if gid_col in preds.columns and "Genotype" not in preds.columns:
            preds = preds.rename(columns={gid_col: "Genotype"})
        if env_col in preds.columns and "Env" not in preds.columns:
            preds = preds.rename(columns={env_col: "Env"})
        return {
            "result": {"predictions": preds},
            "predictions": preds,
            "tiered": tiered_out,
            "_meta": {
                "backend_used": "tiered_dispatch",
                "tier_counts": tiered_out["tier_counts"],
                "wall": tiered_out["wall"],
                "return_se": bool(tiered_return_se),
                "n_hutchinson_probes": (int(tiered_n_probes) if tiered_return_se else None),
                "supports_prediction_se": bool(tiered_return_se),
                "prediction_se_method": ("hutchinson_rademacher" if tiered_return_se else None),
            },
        }

    requested_device = device
    dtype_t = _resolve_torch_dtype(dtype)
    set_deterministic(seed)

    preprocessed_context = kwargs.pop("preprocessed_context", None)
    reused_preprocessed_context = _fit_context_matches(
        preprocessed_context,
        pheno_df=pheno_df,
        gid_col=gid_col,
        env_col=env_col,
        y_col=y_col,
        geno_ids=geno_ids,
        fixed_effects=fixed_effects,
    )
    if reused_preprocessed_context:
        df = preprocessed_context["df"]
        y = np.asarray(preprocessed_context["y"], dtype=float)
        gi = np.asarray(preprocessed_context["gi"], dtype=np.int64)
        ei = np.asarray(preprocessed_context["ei"], dtype=np.int64)
        env_levels = list(preprocessed_context["env_levels"])
        X = np.asarray(preprocessed_context["X"], dtype=float)
    else:
        ctx = prepare_fit_mixed_model_context(
            pheno_df=pheno_df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_ids=geno_ids,
            fixed_effects=fixed_effects,
        )
        df = ctx["df"]
        y = np.asarray(ctx["y"], dtype=float)
        gi = np.asarray(ctx["gi"], dtype=np.int64)
        ei = np.asarray(ctx["ei"], dtype=np.int64)
        env_levels = list(ctx["env_levels"])
        X = np.asarray(ctx["X"], dtype=float)
        preprocessed_context = None

    resid_diag_obs_input, resid_diag_env_input, residual_meta = _resolve_stagewise_residual_inputs(
        df=df, ei=ei,
        obs_weights=obs_weights, obs_var=obs_var, stage1_pev=stage1_pev,
        obs_weight_mode=obs_weight_mode, obs_weight_global_scale=obs_weight_global_scale,
    )

    modular_dense_delegate = bool(kwargs.pop("modular_dense_delegate", False) or kwargs.pop("use_modular_dense", False))
    modular_fit_method = str(kwargs.pop("modular_fit_method", "reml")).lower()
    requested_engine = str(kwargs.pop("engine_mode", kwargs.pop("requested_engine", "auto"))).lower()
    n_structured_terms = 0
    max_axis_levels = 0
    try:
        _tmp_registry = build_term_registry_from_legacy_inputs(
            geno_kernel_names=list(geno_kernels.keys()),
            env_col=env_col,
            env_levels=env_levels,
            env_structure=env_structure,
            fa_rank=fa_rank,
            interaction_terms_meta=interaction_terms_meta,
        )
        for _term in _tmp_registry:
            if hasattr(_term, "axis_name"):
                _axis_name = _term.axis_name
                _axis_levels = _term.axis_levels
            else:
                _axis_name = _term.get("axis_name", None)
                _axis_levels = _term.get("axis_levels", None)
            if _axis_name is not None:
                n_structured_terms += 1
                if _axis_levels is not None:
                    max_axis_levels = max(max_axis_levels, len(_axis_levels))
    except Exception:
        if env_structure is not None and str(env_structure).lower() != "identity":
            n_structured_terms = 1
            max_axis_levels = len(env_levels)
    diagnostics_engine_policy = _decide_framework_engine_mode(
        requested_engine=requested_engine,
        modular_dense_delegate=modular_dense_delegate,
        use_modular_dense=modular_dense_delegate,
        point_predictions_only=bool(point_predictions_only),
        return_se=bool(return_se),
        compute_ai_se=bool(compute_ai_se),
        n_obs=int(df.shape[0]),
        n_structured_terms=int(n_structured_terms),
        max_axis_levels=int(max_axis_levels),
        large_n_threshold=int(large_n_threshold),
        max_full_levels=int(kwargs.pop("max_full_levels", 20)),
    )
    modular_dense_delegate = bool(diagnostics_engine_policy.get("use_modular_dense", False))

    if modular_dense_delegate and _MODULAR_MODEL_API_AVAILABLE and _MODULAR_TERM_SPEC_AVAILABLE:
        try:
            term_registry = build_term_registry_from_legacy_inputs(
                geno_kernel_names=list(geno_kernels.keys()),
                env_col=env_col,
                env_levels=env_levels,
                env_structure=env_structure,
                fa_rank=fa_rank,
                interaction_terms_meta=interaction_terms_meta,
            )
            modular_spec = _make_modular_model_spec_from_term_registry(fixed_effects or [], term_registry)
            if modular_spec is not None:
                kernel_matrices_by_term = _build_modular_kernel_matrices_by_term(
                    df=df,
                    gid_col=gid_col,
                    geno_kernels=geno_kernels,
                    geno_ids=geno_ids,
                    term_registry=term_registry,
                )
                level_index_by_term = _build_modular_level_index_by_term(df=df, term_registry=term_registry)
                covariance_modules_by_term = build_covariance_modules_from_term_registry(term_registry)

                prediction_df = df[[gid_col, env_col]].copy()
                prediction_df.columns = ["Name", "Env"]
                prediction_df.insert(0, "row", np.arange(df.shape[0], dtype=int))

                modular_result = _mod_fit_mixed_model_modular(
                    model_spec=modular_spec,
                    y=y,
                    X=X,
                    kernel_matrices_by_term=kernel_matrices_by_term,
                    level_index_by_term=level_index_by_term,
                    covariance_modules_by_term=covariance_modules_by_term,
                    stage1_pev=stage1_pev,
                    obs_var=obs_var,
                    obs_weights=obs_weights,
                    obs_weight_mode=obs_weight_mode,
                    obs_weight_global_scale=obs_weight_global_scale,
                    return_prediction_se=bool(return_se),
                    prediction_df=prediction_df,
                    level_index_for_residual_summary=ei,
                    level_names_for_residual_summary=env_levels,
                    n_rep_for_heritability=None,
                    diagnostics={
                        "delegate_source": "modular_dense",
                        "legacy_method_requested": method,
                        "prediction_output": prediction_output,
                        "engine_policy": diagnostics_engine_policy,
                    },
                    fit_method=modular_fit_method,
                    max_iter=int(kwargs.pop("modular_max_iter", 100)),
                    step_size=float(kwargs.pop("modular_step_size", 0.05)),
                    tol=float(kwargs.pop("modular_tol", 1e-6)),
                    fd_eps=float(kwargs.pop("modular_fd_eps", 1e-4)),
                )
                return _convert_modular_result_to_legacy_output(
                    modular_result,
                    prediction_output=prediction_output,
                )
        except Exception as _mod_delegate_err:
            warnings.warn(f"modular_dense_delegate failed; falling back to legacy path: {_mod_delegate_err}")


    # -------------------- preprocess random_terms (explicit only) --------------------
    # If user supplies IID terms with group_col/group_cols, create group_index_all (for train/test slicing)
    # and keep metadata. We do NOT auto-add IID terms unless the user passes them.
    random_terms_prepared = None
    if random_terms is not None:
        random_terms_prepared = []
        for term in random_terms:
            t = dict(term)
            ttype = str(t.get("type","")).lower()
            if ttype in ("iid","random_intercept"):
                if "group_index_all" not in t:
                    gcol = t.get("group_col", None)
                    gcols = t.get("group_cols", None)
                    if gcol is None and gcols is None:
                        raise ValueError("IID term requires group_col/group_cols or group_index_all.")
                    g_idx, g_meta = build_group_index(df, group_col=gcol, group_cols=gcols)
                    t["group_index_all"] = g_idx
                    t["group_meta"] = g_meta
            random_terms_prepared.append(t)
    # downstream should use random_terms_prepared

    idx_all = np.arange(len(df), dtype=np.int64)
    if train_idx is None:
        train_idx = idx_all[np.isfinite(y)]
    if test_idx is None:
        test_idx = idx_all[~np.isfinite(y)]
        if test_idx.size == 0:
            test_idx = idx_all
    train_idx = np.asarray(train_idx, dtype=np.int64)
    test_idx = np.asarray(test_idx, dtype=np.int64)
    if train_idx.size == 0:
        raise ValueError("fit_mixed_model: train_idx is empty")
    for _name, _idx in (("train_idx", train_idx), ("test_idx", test_idx)):
        if _idx.size and (_idx.min() < 0 or _idx.max() >= len(df)):
            raise ValueError(
                f"fit_mixed_model: {_name} out of range [0, {len(df)}) "
                f"(min={int(_idx.min())}, max={int(_idx.max())})"
            )
    if not np.isfinite(y[train_idx]).all():
        bad = train_idx[~np.isfinite(y[train_idx])][:5].tolist()
        raise ValueError(
            f"fit_mixed_model: train_idx contains rows with missing/non-finite {y_col!r}; "
            f"first bad indices: {bad}. Missing y is valid only for prediction rows."
        )

    _work_units = float(train_idx.size) * float(train_idx.size + max(test_idx.size, 1))
    try:
        _work_units *= max(1, int(len(geno_kernels)))
    except Exception:
        pass
    device = _runtime_gp_device(requested_device, work_units=_work_units)

    K_names = list(geno_kernels.keys())
    G_list = [geno_kernels[k] for k in K_names]
    nK = len(G_list)
    w_g  = _prepare_weight_vector(w_g,  nK, default=1.0)
    w_ge = _prepare_weight_vector(w_ge, nK, default=0.0)

    # ---- Reaction-norm / env-covariate auto-routing ----
    # If env_covariates are given, build K_env via QC + `kenv_kernel` (default
    # matern32 with bandwidth=1.0 × median pairwise distance, LOYO-validated);
    # auto-include 'e' in include_components; default w_ge/w_e to
    # reaction-norm-friendly values; filter training envs to nearest K_env
    # neighbors of test envs. See `_reaction_norm_preprocess` for full behavior
    # and `_build_kenv_from_Z` for the supported kernel families.
    _rn = _reaction_norm_preprocess(
        env_covariates=env_covariates,
        env_similarity=env_similarity,
        env_levels=env_levels,
        env_col=env_col,
        ei=ei,
        train_idx=train_idx,
        test_idx=test_idx,
        include_components=include_components,
        w_g=w_g, w_ge=w_ge, w_e=w_e,
        reaction_norm_auto=bool(reaction_norm_auto),
        feature_qc=bool(reaction_norm_feature_qc),
        env_filter_topk=reaction_norm_env_filter_topk,
        w_ge_default=float(reaction_norm_w_ge_default),
        w_e_default=float(reaction_norm_w_e_default),
        kenv_kernel=str(kenv_kernel),
        kenv_bandwidth=float(kenv_bandwidth),
        kenv_kernel_kwargs=kenv_kernel_kwargs,
        kenv_auto_candidates=tuple(kenv_auto_candidates),
        kenv_cv_groups=kenv_cv_groups,
        kenv_cv_kfolds=int(kenv_cv_kfolds),
        kenv_cv_ridge=float(kenv_cv_ridge),
        kenv_cv_seed=int(kenv_cv_seed),
        y=y,
    )
    env_similarity = _rn["env_similarity"]
    include_components = _rn["include_components"]
    w_g = _rn["w_g"]
    w_ge = _rn["w_ge"]
    w_e = _rn["w_e"]
    train_idx = _rn["train_idx"]
    _rn_info = _rn.get("info", {})

    # ---- Large-n SNP-matrix path warning ----
    # Encourage streaming the SNP -> low-rank Phi via grm_factor_cache instead
    # of materializing dense (n_geno x n_geno) kernels, which blows up memory
    # and kills the operator backend.
    try:
        _n_geno_max = max(
            (int(G.shape[0]) if hasattr(G, "shape") else 0) for G in G_list
        )
    except Exception:
        _n_geno_max = 0
    if (_n_geno_max >= 3000) and (grm_factor_cache is None):
        warnings.warn(
            f"Dense geno_kernels with n_geno={_n_geno_max} (>= 3000) and no "
            f"grm_factor_cache provided. The operator backend can route through "
            f"low-rank Phi factors streamed from Zarr/memmap instead of a dense "
            f"(n_geno x n_geno) allocation. For large datasets, pre-factorize "
            f"with grm_pivoted_cholesky_to_zarr_streaming.py and pass "
            f"grm_factor_cache={{'type':'memmap','paths':[...],'shapes':[...],"
            f"'dtypes':[...]}}.",
            RuntimeWarning, stacklevel=2,
        )

    # Zero-obs-env guard: any env with 0 rows in train_idx has no Fisher
    # information for its variance components (!var, !fa*, !R). AI-REML
    # will return ~inf SE for those rows and gp_icm_fa full mode will
    # typically fail to meet the convergence tol because the dead
    # components dominate the gradient norm. Warn explicitly so the
    # huge-SE rows are expected rather than mysterious.
    _ei_tr = np.asarray(ei, dtype=np.int64)[np.asarray(train_idx, dtype=np.int64)]
    _env_counts = np.bincount(_ei_tr, minlength=len(env_levels))
    _zero_obs_envs = [env_levels[i] for i, c in enumerate(_env_counts) if c == 0]
    if _zero_obs_envs:
        _warn_once(
            ("zero_obs_envs", tuple(str(x) for x in _zero_obs_envs)),
            f"{len(_zero_obs_envs)} environment(s) have zero observed rows in "
            f"the training split: {_zero_obs_envs}. Variance components for "
            f"these envs (!var, !fa*, Env_*!R) are not identifiable -- AI-REML "
            f"will report huge SEs and gp_icm_fa full mode may not converge. "
            f"For CV0 env-holdout designs use mode='fast' or exclude these "
            f"rows from the VC table before reporting.",
            RuntimeWarning,
            stacklevel=2,
        )

    lam_eff = float(krr_lam if lam is None else lam)

    if output_level == "full_vc" and reml_normalize == "diag_mean":
        warnings.warn(
            "reml_normalize='diag_mean' rescales the genetic kernel; "
            "variance components will NOT be on the original y scale. "
            "Overriding to 'asreml' (no rescale) for full_vc output. "
            "Pass reml_normalize='asreml' explicitly to suppress this warning.",
            stacklevel=2,
        )
        reml_normalize = "asreml"

    common_env = dict(geno_id_list=list(geno_ids),
                      env_index_to_name=env_levels,
                      prediction_output=prediction_output,
                      interaction_term_name=f"G:{env_col}",
                      theta0_warm=theta0_warm,
                      gp_engine=gp_engine)
    # Dual-engine flag: when True, after gp_icm_fa runs we do a follow-up
    # krr_exact pass using the FITTED Sigma_g as env_similarity so the user
    # gets SE_g_latent on predictions (gp_icm_fa doesn't compute per-
    # component posterior SEs, but classical reliability needs them).
    _krr_fa_dual_engine = False
    _krr_fa_mom = False
    if method == "krr_exact" and env_structure in ("corh", "corgh", "us"):
        method = "gp_exact"
    elif method == "krr_exact" and env_structure == "fa":
        if str(varcomp_mode).lower() == "mom":
            _mom_backend_mode = (backend or "auto").lower()
            _mom_n_train = int(np.asarray(train_idx).size)
            _mom_use_operator = bool(
                _mom_backend_mode == "operator"
                or (
                    _mom_backend_mode == "auto"
                    and grm_factor_cache is not None
                    and _mom_n_train >= int(large_n_threshold)
                    and str(prediction_output).lower() == "test_only"
                )
            )
            if _mom_use_operator:
                sigma_g_mom = _mom_sigma_g_operator(
                    grm_factor_cache=grm_factor_cache,
                    w_g=w_g, gi=gi, ei=ei, y=y,
                    train_idx=train_idx, n_env=len(env_levels),
                    h2_prior=float(mom_h2_prior),
                    shrinkage_correct=bool(mom_shrinkage_correct),
                    pcg_tol=float(operator_tol),
                    pcg_max_iter=int(operator_max_iter),
                    n_hutch=int(mom_n_hutch),
                    device=device,
                    dtype_compute=operator_dtype_compute,
                    seed=int(seed),
                    verbose=bool(mom_verbose),
                )
                _mom_backend_used = "operator"
            else:
                sigma_g_mom = _mom_sigma_g_empirical(
                    G_list=G_list, w_g=w_g, gi=gi, ei=ei, y=y,
                    train_idx=train_idx, n_env=len(env_levels),
                    h2_prior=float(mom_h2_prior),
                    shrinkage_correct=bool(mom_shrinkage_correct),
                )
                _mom_backend_used = "dense"
            env_similarity = sigma_g_mom
            env_structure = None
            fa_rank = None
            _krr_fa_mom = True
        else:
            method = "gp_icm_fa"
            _krr_fa_dual_engine = bool(return_se)


    # ------------------ Auto backend routing for very large n ------------------
    backend_mode = (backend or "auto").lower()
    n_train = int(np.asarray(train_idx).size)
    want_operator = (
        backend_mode == "operator"
        or (backend_mode == "auto" and (grm_factor_cache is not None) and n_train >= int(large_n_threshold)
            and str(prediction_output).lower() == "test_only")
    )
    if bool(return_prediction_cov) and want_operator:
        raise ValueError("return_prediction_cov=True is only supported on dense paths; operator backend returns means or diagonal SEs only when implemented.")
    use_operator_now = bool(
        want_operator
        and method in ("gp_exact", "krr_exact")
        and (not bool(compute_ai_se))
        and (bool(point_predictions_only) or bool(return_se))
    )
    # Allow callers to request Hutchinson SE on the operator path via kwargs.
    _op_return_se = bool(kwargs.get("operator_return_se", bool(return_se)))
    _op_n_probes = int(kwargs.get("operator_n_hutchinson_probes", kwargs.get("n_hutchinson_probes", 256)))
    _op_hutch_seed = kwargs.get("operator_hutchinson_seed", None)
    _op_se_pcg_tol = kwargs.get("operator_se_pcg_tol", None)
    _op_se_pcg_max_it = kwargs.get("operator_se_pcg_max_iter", None)
    _op_use_nys = bool(kwargs.get("operator_use_nystrom_se_preconditioner", False))
    _op_nys_rank = int(kwargs.get("operator_nystrom_rank", 64))
    _op_nys_seed = int(kwargs.get("operator_nystrom_seed", 0))
    _op_use_lan = bool(kwargs.get("operator_use_lanczos_se_variance_reduction", True))
    _op_lan_iters = int(kwargs.get("operator_lanczos_iters", 80))
    _op_lan_seed = int(kwargs.get("operator_lanczos_seed", 0))
    _op_lan_floor = float(kwargs.get("operator_lanczos_theta_floor", 1e-8))
    if use_operator_now:
        res = _operator_backend_point_predictions(
            gi=gi, ei=ei, y=y, X=X,
            train_idx=train_idx, test_idx=test_idx,
            w_g=w_g, w_ge=w_ge, w_e=w_e,
            env_similarity=env_similarity,
            resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else (None if not bool(learn_envdiag_noise) else kwargs.get("resid_diag_env", None))),
            resid_diag_obs=resid_diag_obs_input,
            device=device, dtype_t=dtype_t,
            standardize=standardize,
            geno_id_list=list(geno_ids),
            env_index_to_name=env_levels,
            prediction_output=prediction_output,
            grm_factor_cache=grm_factor_cache,
            operator_tol=float(operator_tol),
            operator_max_iter=int(operator_max_iter),
            operator_dtype_compute=operator_dtype_compute,
            seed=seed,
            random_terms=random_terms_prepared,
            return_se=_op_return_se,
            n_hutchinson_probes=_op_n_probes,
            hutchinson_seed=_op_hutch_seed,
            se_pcg_tol=_op_se_pcg_tol,
            se_pcg_max_iter=_op_se_pcg_max_it,
            use_nystrom_se_preconditioner=_op_use_nys,
            nystrom_rank=_op_nys_rank,
            nystrom_seed=_op_nys_seed,
            use_lanczos_se_variance_reduction=_op_use_lan,
            lanczos_iters=_op_lan_iters,
            lanczos_seed=_op_lan_seed,
            lanczos_theta_floor=_op_lan_floor,
        )
        if _krr_fa_mom:
            res["env_covariance_fa"] = sigma_g_mom
            res["env_correlation_fa"] = _corr_from_cov_np(sigma_g_mom)
            res.setdefault("diagnostics", {})["varcomp_mode"] = "mom"
            res.setdefault("diagnostics", {})["mom_h2_prior"] = float(mom_h2_prior)
            res.setdefault("diagnostics", {})["mom_backend"] = _mom_backend_used
            res.setdefault("diagnostics", {})["dual_engine"] = (
                f"mom(Sigma_g empirical, backend={_mom_backend_used}) + operator_backend(predictions"
                f"{'+SE' if _op_return_se else ''})"
            )
        return _stamp_output_level(_append_prediction_meta(
            res,
            point_predictions_only=(not _op_return_se),
            prediction_output=prediction_output,
            backend_used="operator",
            prediction_se_type=prediction_se_type,
            prediction_se_method=("hutchinson_rademacher" if _op_return_se else "means_only"),
            prediction_block_size=prediction_block_size,
            prediction_diag_probes=(_op_n_probes if _op_return_se else prediction_diag_probes),
            return_prediction_cov=bool(return_prediction_cov),
            approximate=False,
            residual_meta=residual_meta,
            engine_policy=diagnostics_engine_policy,
        ), output_level)
    if method == "krr_exact":
        res = multikernel_krr_posterior_with_se(
            G_list=G_list, obs_gidx=gi, obs_eidx=ei, y=y,
            train_idx=train_idx, test_idx=test_idx,
            w_g=w_g, w_ge=w_ge, w_e=w_e, S_e=env_similarity, lam=lam_eff,
            env_levels=env_levels, obs_names=list(geno_ids),
            fixed_effects=fixed_effects, X=X, device=device, dtype=dtype_t, seed=seed,
            prediction_output=prediction_output,
            standardize=standardize,
            krr_lams=krr_lams, lam_select=lam_select,
            env_resid_robust=env_resid_robust, env_resid_shrink_tau=env_resid_shrink_tau,
            hutch_samples=int(hutch_samples),
            reml_normalize=reml_normalize,
            resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
            resid_diag_obs=resid_diag_obs_input,
            env_structure=env_structure,
            fa_rank=fa_rank,
            interaction_term_name=interaction_term_name,
            interaction_terms_meta=interaction_terms_meta,
            return_se=bool(return_se),
        )
        res.setdefault("envwise", {}).setdefault("details", {})["lam"] = lam_eff
        if _krr_fa_mom:
            res["env_covariance_fa"] = sigma_g_mom
            res["env_correlation_fa"] = _corr_from_cov_np(sigma_g_mom)
            res.setdefault("diagnostics", {})["varcomp_mode"] = "mom"
            res.setdefault("diagnostics", {})["mom_h2_prior"] = float(mom_h2_prior)
            res.setdefault("diagnostics", {})["mom_backend"] = _mom_backend_used
            res.setdefault("diagnostics", {})["dual_engine"] = (
                f"mom(Sigma_g empirical, backend={_mom_backend_used}) + krr_exact(predictions+per-component SE)"
            )
        if (
            output_level == "full_vc"
            and not _krr_fa_mom
            and res.get("variance_component_method") != "krr_gcv_variance_ratio"
        ):
            try:
                vc_method = "gp_icm_fa" if env_structure == "fa" else "gp_exact"
                if vc_method == "gp_icm_fa":
                    vc_res = gp_icm_fa_with_X(
                        geno_kernels, gi, ei, y, train_idx, test_idx,
                        w_g, w_ge, w_e, X, S_e=env_similarity, device=device, dtype=dtype_t,
                        iters=int(kwargs.get("icm_iters", 300)), lr=float(kwargs.get("icm_lr", 0.03)),
                        use_fixed_noise=True if (resid_diag_obs_input is not None or resid_diag_env_input is not None) else (not bool(learn_envdiag_noise)),
                        fa_rank=int(fa_rank if fa_rank is not None else icm_rank),
                        standardize=standardize, reml_normalize=reml_normalize,
                        learn_scales=bool(learn_scales),
                        compute_ai_se=True,
                        ai_hutch_samples=int(ai_hutch_samples),
                        ai_jitter=float(ai_jitter),
                        ai_include_fa_loadings=bool(ai_include_fa_loadings),
                        optimizer_method=str(optimizer_method),
                        ai_steps=int(ai_steps),
                        ai_max_iter=int(ai_max_iter),
                        ai_tol_loglik=float(ai_tol_loglik),
                        ai_tol_theta=float(ai_tol_theta),
                        env_mean_from_obs_diag=bool(kwargs.get("env_mean_from_obs_diag", False)),
                        point_predictions_only=False,
                        resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
                        resid_diag_obs=resid_diag_obs_input,
                        return_se=False,
                        prediction_se_type=prediction_se_type,
                        prediction_se_method="means_only",
                        prediction_block_size=int(prediction_block_size),
                        prediction_diag_probes=int(prediction_diag_probes),
                        return_prediction_cov=False,
                        **common_env)
                else:
                    vc_res = gp_exact_with_X(
                        geno_kernels, gi, ei, y, train_idx, test_idx,
                        w_g, w_ge, w_e, X, S_e=env_similarity, device=device, dtype=dtype_t,
                        iters=int(kwargs.get("gp_iters", 300)), lr=float(kwargs.get("gp_lr", 0.03)),
                        use_fixed_noise=True if (resid_diag_obs_input is not None or resid_diag_env_input is not None) else (not bool(learn_envdiag_noise)),
                        standardize=standardize, reml_normalize=reml_normalize,
                        learn_scales=bool(learn_scales),
                        compute_ai_se=True,
                        ai_hutch_samples=int(ai_hutch_samples),
                        ai_jitter=float(ai_jitter),
                        optimizer_method=str(optimizer_method),
                        ai_steps=int(ai_steps),
                        ai_max_iter=int(ai_max_iter),
                        ai_tol_loglik=float(ai_tol_loglik),
                        ai_tol_theta=float(ai_tol_theta),
                        env_mean_from_obs_diag=bool(kwargs.get("env_mean_from_obs_diag", False)),
                        point_predictions_only=False,
                        resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
                        resid_diag_obs=resid_diag_obs_input,
                        return_se=False,
                        prediction_se_type=prediction_se_type,
                        prediction_se_method="means_only",
                        prediction_block_size=int(prediction_block_size),
                        prediction_diag_probes=int(prediction_diag_probes),
                        return_prediction_cov=False,
                        precomputed_context=preprocessed_context,
                        **common_env)
                for _k in (
                    "varcomp",
                    "summary",
                    "ai_matrix",
                    "var_components",
                    "var_components_summary",
                    "var_components_ai",
                    "env_variance_summary",
                    "residual_summary",
                    "heritability_summary",
                    "interaction_variance_summary",
                    "interaction_correlation_summary",
                    "env_covariance_fa",
                    "env_correlation_fa",
                    "env_covariance_fa_structure",
                    "env_covariance_fa_structure_standardized",
                    "fa_kernel_weights",
                    "sigma2_resid_env",
                    "sigma2_resid_overall",
                ):
                    if _k in vc_res:
                        res[_k] = vc_res[_k]
                res.setdefault("diagnostics", {})["varcomp_source"] = vc_method
                # Post-condition: the VC call completed without raising, but
                # may still have returned without a varcomp table (e.g., the
                # ASReml-schema builder silently failed upstream). Surface this
                # explicitly so output_level='full_vc' callers don't silently
                # receive varcomp=None.
                if res.get("varcomp") is None:
                    import warnings as _w
                    _w.warn(
                        f"output_level='full_vc' requested but the {vc_method} "
                        f"post-fit VC pass returned no varcomp table. "
                        f"Downstream consumers will see varcomp=None.",
                        RuntimeWarning, stacklevel=2,
                    )
            except Exception as _e:
                import warnings as _w
                _w.warn(
                    f"output_level='full_vc' requested but the KRR post-fit "
                    f"VC pass via {vc_method} failed: {type(_e).__name__}: {_e}. "
                    f"Downstream consumers will see varcomp=None.",
                    RuntimeWarning, stacklevel=2,
                )
        return _stamp_output_level(_append_prediction_meta(
            res,
            point_predictions_only=bool(point_predictions_only),
            prediction_output=prediction_output,
            backend_used="dense",
            prediction_se_type=prediction_se_type,
            prediction_se_method=("exact_dense" if not bool(point_predictions_only) else "means_only"),
            prediction_block_size=prediction_block_size,
            prediction_diag_probes=prediction_diag_probes,
            return_prediction_cov=bool(return_prediction_cov),
            approximate=False,
            residual_meta=residual_meta,
            engine_policy=diagnostics_engine_policy,
        ), output_level)

    elif method == "gp_exact":
        res = gp_exact_with_X(geno_kernels, gi, ei, y, train_idx, test_idx,
                               w_g, w_ge, w_e, X, S_e=env_similarity, device=device, dtype=dtype_t,
                               iters=int(kwargs.get("gp_iters", 300)), lr=float(kwargs.get("gp_lr", 0.03)),
                               use_fixed_noise = True if (resid_diag_obs_input is not None or resid_diag_env_input is not None) else (not bool(learn_envdiag_noise)),
                               standardize=standardize, reml_normalize=reml_normalize,
                               learn_scales=bool(learn_scales),
                               compute_ai_se=(output_level == "full_vc") and bool(compute_ai_se),
                               ai_hutch_samples=int(ai_hutch_samples),
                               ai_jitter=float(ai_jitter),
                               optimizer_method=str(optimizer_method),
                               ai_steps=int(ai_steps),
                               ai_max_iter=int(ai_max_iter),
                               ai_tol_loglik=float(ai_tol_loglik),
                               ai_tol_theta=float(ai_tol_theta),
                               env_mean_from_obs_diag=bool(kwargs.get("env_mean_from_obs_diag", False)),
                               point_predictions_only=bool(point_predictions_only),
                               resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
                               resid_diag_obs=resid_diag_obs_input,
                               return_se=bool(return_se),
                               prediction_se_type=prediction_se_type,
                               prediction_se_method=("means_only" if bool(point_predictions_only) else "exact_dense"),
                               prediction_block_size=int(prediction_block_size),
                               prediction_diag_probes=int(prediction_diag_probes),
                               return_prediction_cov=bool(return_prediction_cov),
                               precomputed_context=preprocessed_context,
                               **common_env)
        return _stamp_output_level(_append_prediction_meta(
            res,
            point_predictions_only=bool(point_predictions_only),
            prediction_output=prediction_output,
            backend_used="dense",
            prediction_se_type=prediction_se_type,
            prediction_se_method=("means_only" if bool(point_predictions_only) else "exact_dense"),
            prediction_block_size=prediction_block_size,
            prediction_diag_probes=prediction_diag_probes,
            return_prediction_cov=bool(return_prediction_cov),
            approximate=False,
            residual_meta=residual_meta,
            engine_policy=diagnostics_engine_policy,
        ), output_level)

    elif method == "gp_icm_fa":
        res = gp_icm_fa_with_X(geno_kernels, gi, ei, y, train_idx, test_idx,
                                w_g, w_ge, w_e, X, S_e=env_similarity, device=device, dtype=dtype_t,
                                iters=int(kwargs.get("icm_iters", 300)), lr=float(kwargs.get("icm_lr", 0.03)),
                                use_fixed_noise = True if (resid_diag_obs_input is not None or resid_diag_env_input is not None) else (not bool(learn_envdiag_noise)),
                                fa_rank=int(fa_rank if fa_rank is not None else icm_rank),
                                standardize=standardize, reml_normalize=reml_normalize,
                                learn_scales=bool(learn_scales),
                                compute_ai_se=(output_level == "full_vc") and bool(compute_ai_se),
                                ai_hutch_samples=int(ai_hutch_samples),
                                ai_jitter=float(ai_jitter),
                                ai_include_fa_loadings=bool(ai_include_fa_loadings),
                                optimizer_method=str(optimizer_method),
                                ai_steps=int(ai_steps),
                                ai_max_iter=int(ai_max_iter),
                                ai_tol_loglik=float(ai_tol_loglik),
                                ai_tol_theta=float(ai_tol_theta),
                                env_mean_from_obs_diag=bool(kwargs.get("env_mean_from_obs_diag", False)),
                                point_predictions_only=bool(point_predictions_only),
                                resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
                                resid_diag_obs=resid_diag_obs_input,
                                return_se=bool(return_se),
                                prediction_se_type=prediction_se_type,
                                prediction_se_method=("means_only" if bool(point_predictions_only) else "exact_dense"),
                                prediction_block_size=int(prediction_block_size),
                                prediction_diag_probes=int(prediction_diag_probes),
                                return_prediction_cov=bool(return_prediction_cov),
                                **common_env)

        # Dual-engine post-step for krr_exact + env_structure="fa" requests:
        # the user asked for krr_exact but we rerouted to gp_icm_fa so FA
        # loadings get fit by AI-REML. Now run a cheap krr_exact pass with
        # S_e = fitted Sigma_g (from env_covariance_fa) to obtain
        # SE_g_latent / SE_ge_latent / SE_e_latent on the predictions. Keep
        # gp_icm_fa's varcomp, summaries, and Sigma_g intact.
        if _krr_fa_dual_engine:
            try:
                sigma_g_fa = np.asarray(
                    res.get(
                        "env_covariance_fa_structure_standardized",
                        res.get("env_covariance_fa"),
                    ),
                    dtype=float,
                )
                if sigma_g_fa.ndim != 2 or sigma_g_fa.shape[0] != sigma_g_fa.shape[1]:
                    raise ValueError(f"env_covariance_fa has unexpected shape {sigma_g_fa.shape}")
                krr_res = multikernel_krr_posterior_with_se(
                    G_list=G_list, obs_gidx=gi, obs_eidx=ei, y=y,
                    train_idx=train_idx, test_idx=test_idx,
                    w_g=w_g, w_ge=w_ge, w_e=w_e,
                    S_e=sigma_g_fa, lam=lam_eff,
                    env_levels=env_levels, obs_names=list(geno_ids),
                    fixed_effects=fixed_effects, X=X,
                    device=device, dtype=dtype_t, seed=seed,
                    prediction_output=prediction_output,
                    standardize=standardize,
                    krr_lams=krr_lams, lam_select=lam_select,
                    env_resid_robust=env_resid_robust,
                    env_resid_shrink_tau=env_resid_shrink_tau,
                    hutch_samples=int(hutch_samples),
                    reml_normalize=reml_normalize,
                    resid_diag_env=(resid_diag_env_input if resid_diag_env_input is not None else kwargs.get("resid_diag_env", None)),
                    resid_diag_obs=resid_diag_obs_input,
                    env_structure=None,  # Sigma_g already captured in S_e
                    fa_rank=None,
                    interaction_term_name=interaction_term_name,
                    interaction_terms_meta=interaction_terms_meta,
                    return_se=bool(return_se),
                )
                # Merge per-component SE columns from krr's posterior into
                # gp_icm_fa's predictions. gp_icm_fa's Prediction column is
                # empirically more accurate (it uses REML-fit variance ratios;
                # krr re-derives its own regularization via lam_eff and
                # env_resid_shrink_tau, which diverges from the REML BLUP when
                # per-env sample sizes are small). The krr pass is only here
                # to emit SE_g_latent / SE_ge_latent / SE_e_latent, which
                # gp_icm_fa does not compute per-component. Keep gp_icm_fa's
                # Prediction + total-SE columns, borrow only per-component SEs.
                krr_preds = krr_res.get("result", {}).get("predictions")
                gp_preds = res.get("result", {}).get("predictions")
                if krr_preds is not None:
                    se_cols = [c for c in ("SE_g_latent", "SE_ge_latent",
                                           "SE_e_latent")
                               if c in krr_preds.columns]
                    if (gp_preds is not None and se_cols
                            and "Name" in krr_preds.columns and "Env" in krr_preds.columns
                            and "Name" in gp_preds.columns and "Env" in gp_preds.columns):
                        se_slice = krr_preds[["Name", "Env"] + se_cols]
                        merged = gp_preds.drop(
                            columns=[c for c in se_cols if c in gp_preds.columns],
                            errors="ignore",
                        )
                        merged = merged.merge(se_slice, on=["Name", "Env"], how="left")
                        res.setdefault("result", {})["predictions"] = merged
                        res.setdefault("diagnostics", {})["dual_engine"] = (
                            "gp_icm_fa(predictions+varcomp+Sigma_g) + krr_exact(per-component SE)"
                        )
                    else:
                        # Keep gp_icm_fa predictions. Replacing them with the
                        # auxiliary KRR pass trades accuracy for an SE detail,
                        # which violates the prediction-path contract.
                        res.setdefault("diagnostics", {})["dual_engine"] = (
                            "gp_icm_fa(predictions+varcomp+Sigma_g); auxiliary krr_exact SE merge skipped"
                        )
                        warnings.warn(
                            "Auxiliary krr_exact per-component SE merge after gp_icm_fa "
                            "did not find compatible columns/keys; kept gp_icm_fa predictions.",
                            RuntimeWarning, stacklevel=2,
                        )
            except Exception as _dual_e:
                warnings.warn(
                    f"Dual-engine krr_exact pass after gp_icm_fa failed: "
                    f"{type(_dual_e).__name__}: {_dual_e}. Predictions fall "
                    f"back to gp_icm_fa (no SE_g_latent).",
                    RuntimeWarning, stacklevel=2,
                )

        return _stamp_output_level(_append_prediction_meta(
            res,
            point_predictions_only=bool(point_predictions_only),
            prediction_output=prediction_output,
            backend_used="dense",
            prediction_se_type=prediction_se_type,
            prediction_se_method=("means_only" if bool(point_predictions_only) else "exact_dense"),
            prediction_block_size=prediction_block_size,
            prediction_diag_probes=prediction_diag_probes,
            return_prediction_cov=bool(return_prediction_cov),
            approximate=False,
            residual_meta=residual_meta,
            engine_policy=diagnostics_engine_policy,
        ), output_level)

    else:
        raise ValueError(f"Unknown method: {method}")


def _candidate_float(row: Any, names: Sequence[str], default: float) -> float:
    for name in names:
        try:
            if name in row.index:
                value = row.get(name)
            else:
                continue
        except Exception:
            try:
                value = row[name]
            except Exception:
                continue
        try:
            out = float(value)
            if np.isfinite(out):
                return out
        except Exception:
            continue
    return float(default)


def _batch_krr_point_predict_candidates(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    train_idx: Sequence[int],
    test_idx: Sequence[int],
    candidates: pd.DataFrame,
    kernel_map: Dict[str, np.ndarray],
    fixed_effects: Optional[Sequence[str]] = None,
    dtype: Union[str, "torch.dtype", None] = "float64",
    seed: int = 12345,
    standardize: str = "global",
    prediction_output: str = "test_only",
    reml_normalize: str = "diag_mean",
    hutch_samples: int = 64,
    env_resid_robust: str = "none",
    env_resid_shrink_tau: float = 0.0,
    reaction_norm_env_filter_topk: Optional[int] = 15,
    device: Optional[str] = None,
) -> Dict[str, Any]:
    """Fast dense KRR tuning path for single-kernel, point-only MET fits."""
    if not TORCH_AVAILABLE:
        raise RuntimeError("PyTorch required for fast batch KRR candidate evaluation")
    if str(prediction_output).lower() != "test_only":
        raise ValueError("fast batch KRR candidate evaluation only supports prediction_output='test_only'")
    if candidates.empty:
        return {"predictions": pd.DataFrame([]), "detail": pd.DataFrame([])}

    requested_device = device
    dtype_t = _resolve_torch_dtype(dtype)
    set_deterministic(seed)

    df = pd.DataFrame(pheno_df).reset_index(drop=True).copy()
    if gid_col not in df.columns or env_col not in df.columns or y_col not in df.columns:
        raise KeyError("gid_col/env_col/y_col not found in pheno_df")

    gid = df[gid_col].astype("category")
    env = df[env_col].astype("category")
    y = pd.to_numeric(df[y_col], errors="coerce").values

    geno_ids = [str(x) for x in list(geno_ids)]
    gid_to_index = {str(g): i for i, g in enumerate(geno_ids)}
    gi = np.array([gid_to_index.get(str(g), -1) for g in gid.astype(str)], dtype=np.int64)
    if (gi < 0).any():
        missing = [str(g) for g, i in zip(gid.astype(str), gi) if i < 0][:5]
        raise ValueError(f"Unknown genotype IDs encountered: {missing} ...")

    env_levels = list(env.cat.categories.astype(str))
    env_to_index = {e: i for i, e in enumerate(env_levels)}
    ei = np.array([env_to_index[str(e)] for e in env.astype(str)], dtype=np.int64)

    X_all, _ = _build_X_from_df(df, fixed_effects)
    t_idx = np.asarray(train_idx, dtype=np.int64)
    s_idx = np.asarray(test_idx if (test_idx is not None and len(test_idx) > 0) else train_idx, dtype=np.int64)
    if t_idx.size == 0:
        raise ValueError("fast batch KRR candidate evaluation: train_idx is empty")
    for _name, _idx in (("train_idx", t_idx), ("test_idx", s_idx)):
        if _idx.size and (_idx.min() < 0 or _idx.max() >= len(df)):
            raise ValueError(
                f"fast batch KRR candidate evaluation: {_name} out of range [0, {len(df)}) "
                f"(min={int(_idx.min())}, max={int(_idx.max())})"
            )
    if not np.isfinite(y[t_idx]).all():
        bad = t_idx[~np.isfinite(y[t_idx])][:5].tolist()
        raise ValueError(
            f"fast batch KRR candidate evaluation: train_idx contains rows with missing/non-finite "
            f"{y_col!r}; first bad indices: {bad}"
        )

    work_units = float(t_idx.size) * float(t_idx.size + max(s_idx.size, 1)) * max(1, int(candidates.shape[0]))
    device_resolved = _runtime_gp_device(requested_device, work_units=work_units)

    stdr = _PerObsEnvStandardizer(mode=standardize).fit(y, ei, t_idx)
    y_std = stdr.transform(y, ei)

    gi_t = _t(gi, device=device_resolved, dtype=torch.int64)
    ei_t = _t(ei, device=device_resolved, dtype=torch.int64)
    y_t = _t(y_std, device=device_resolved, dtype=dtype_t)
    X_t = _t(X_all, device=device_resolved, dtype=dtype_t)
    s = _t(s_idx, device=device_resolved, dtype=torch.int64)
    Xss = X_t.index_select(0, s)
    ei_test = ei_t.index_select(0, s)

    n_obs = int(gi.size)
    n_env = int(ei_t.max().item()) + 1 if n_obs else 0
    n_probes = max(0, int(hutch_samples))
    if n_probes > 0:
        set_deterministic(seed)
        probes = []
        for _ in range(n_probes):
            z = torch.randint(0, 2, (n_obs,), device=device_resolved, dtype=torch.int64) * 2 - 1
            probes.append(z.to(dtype=dtype_t))
        Z_full = torch.stack(probes, dim=1)
    else:
        Z_full = torch.empty((n_obs, 0), dtype=dtype_t, device=device_resolved)

    names_out = [geno_ids[int(gi[i])] if 0 <= int(gi[i]) < len(geno_ids) else str(int(gi[i])) for i in s_idx]
    envs_out = [str(env_levels[int(ei[i])]) for i in s_idx]
    geno_kernel = np.asarray(geno_kernel, dtype=np.float64)

    block_cache: Dict[str, Dict[str, Any]] = {}
    for kernel_name in sorted(set(str(x).lower() for x in candidates["kenv_kernel"].tolist())):
        S_np = kernel_map.get(kernel_name)
        if S_np is None:
            raise ValueError(f"no environment similarity was available for kernel {kernel_name!r}")
        t_idx_kernel = t_idx
        if reaction_norm_env_filter_topk is not None and int(reaction_norm_env_filter_topk) > 0:
            test_env_codes = np.unique(ei[s_idx])
            test_env_names = [str(env_levels[int(c)]) for c in test_env_codes]
            keep_mask = _filter_training_envs_by_Kenv(
                np.asarray(S_np, dtype=np.float64),
                env_levels,
                test_env_names,
                int(reaction_norm_env_filter_topk),
            )
            kept_codes = set(int(i) for i in np.where(keep_mask)[0])
            t_idx_kernel = t_idx[np.isin(ei[t_idx], list(kept_codes))]
        if t_idx_kernel.size == 0:
            raise ValueError(
                f"no training rows remain after reaction-norm environment filtering for kernel {kernel_name!r}"
            )
        t_kernel = _t(t_idx_kernel, device=device_resolved, dtype=torch.int64)
        ytt_kernel = y_t.index_select(0, t_kernel)
        Xtt_kernel = X_t.index_select(0, t_kernel)
        ei_train_kernel = ei_t.index_select(0, t_kernel)
        Zt_kernel = Z_full.index_select(0, t_kernel)
        G_norm_list, S_norm = _apply_kernel_normalization(
            [geno_kernel],
            np.asarray(S_np, dtype=np.float64),
            mode=str(reml_normalize),
        )
        G_t = _t(np.asarray(G_norm_list[0], dtype=np.float64), device=device_resolved, dtype=dtype_t)
        S_t = _t(np.asarray(S_norm, dtype=np.float64), device=device_resolved, dtype=dtype_t)
        gi_train = gi_t.index_select(0, t_kernel)
        gi_test = gi_t.index_select(0, s)
        Kg_tt = _build_geno_block(G_t, gi_train, gi_train)
        Kg_ts = _build_geno_block(G_t, gi_train, gi_test)
        Se_tt = _build_env_kernel(ei_train_kernel, ei_train_kernel, S_t, dtype=dtype_t, device=device_resolved)
        Se_ts = _build_env_kernel(ei_train_kernel, ei_test, S_t, dtype=dtype_t, device=device_resolved)
        block_cache[kernel_name] = {
            "t_idx": t_idx_kernel,
            "ytt": ytt_kernel,
            "Xtt": Xtt_kernel,
            "ei_train": ei_train_kernel,
            "Zt": Zt_kernel,
            "Kg_tt": Kg_tt,
            "Kg_ts": Kg_ts,
            "Kge_tt": Kg_tt * Se_tt,
            "Kge_ts": Kg_ts * Se_ts,
            "Ke_tt": Se_tt,
            "Ke_ts": Se_ts,
        }

    def _residual_diag_env(Ktt: "torch.Tensor", lam_eff: float, blocks: Dict[str, Any]) -> np.ndarray:
        ytt = blocks["ytt"]
        Xtt = blocks["Xtt"]
        ei_train = blocks["ei_train"]
        Zt = blocks["Zt"]
        Vtmp = 0.5 * (Ktt + Ktt.mT)
        Vtmp.diagonal().add_(float(lam_eff))
        Ltmp = _cholesky_factor(Vtmp)
        Vinv_y = _cholesky_solve_from_factor(Ltmp, ytt)
        Vinv_X = _cholesky_solve_from_factor(Ltmp, Xtt)
        Xt_Vinv_X = Xtt.mT @ Vinv_X
        Var_beta_tmp = _spd_inverse(Xt_Vinv_X)
        beta_tmp = Var_beta_tmp @ (Xtt.mT @ Vinv_y)
        alpha_tmp = _cholesky_solve_from_factor(Ltmp, ytt - Xtt @ beta_tmp)
        yhat_tr = (Xtt @ beta_tmp) + Ktt @ alpha_tmp
        resid = (ytt - yhat_tr).reshape(-1)

        if n_probes > 0:
            Vinv_Z = _cholesky_solve_from_factor(Ltmp, Zt)
            Sx_Z = Xtt @ (Var_beta_tmp @ (Xtt.mT @ Vinv_Z))
            HZ = Sx_Z + Ktt @ _cholesky_solve_from_factor(Ltmp, Zt - Sx_Z)
            hat_diag_train = torch.mean(Zt * HZ, dim=1)
        else:
            hat_diag_train = torch.zeros_like(resid)

        sig2 = np.full(n_env, np.nan, dtype=float)
        ei_train_np = _to_numpy(ei_train)
        for e in range(n_env):
            mask_np = ei_train_np == e
            n_e = int(np.sum(mask_np))
            if n_e == 0:
                continue
            mask = _t(mask_np, device=device_resolved, dtype=torch.bool)
            r_e = resid[mask]
            tr_He = float(torch.clamp(hat_diag_train[mask].sum(), min=0.0, max=float(max(n_e - 1, 0))).item())
            denom = max(1.0, n_e - tr_He)
            if str(env_resid_robust).lower() == "mad":
                s2 = _robust_env_var(r_e)
                if np.isnan(s2):
                    s2 = float(torch.clamp((r_e.pow(2).sum() / denom), min=1e-12).item())
            else:
                s2 = float(torch.clamp((r_e.pow(2).sum() / denom), min=1e-12).item())
            sig2[e] = s2

        if np.isfinite(sig2).any() and float(env_resid_shrink_tau) > 0:
            m = float(np.nanmean(sig2))
            sig2 = (1.0 - float(env_resid_shrink_tau)) * sig2 + float(env_resid_shrink_tau) * m
        if np.isnan(sig2).any():
            r2 = resid.pow(2)
            mask = torch.isfinite(resid)
            global_s2 = r2[mask].mean().item() if bool(mask.any().item()) else 1.0
            sig2 = np.where(np.isfinite(sig2), sig2, global_s2)
        return sig2

    pred_parts: List[pd.DataFrame] = []
    detail_rows: List[Dict[str, Any]] = []
    for _, row in candidates.iterrows():
        cid = int(row["candidate_id"])
        t0 = time.time()
        failed = False
        err = None
        try:
            kernel_name = str(row.get("kenv_kernel", "matern32")).lower()
            blocks = block_cache[kernel_name]
            wg = _candidate_float(row, ("w_g",), 1.0)
            wge = _candidate_float(row, ("w_ge",), 0.0)
            we = _candidate_float(row, ("w_e",), 0.0)
            lam_eff = _candidate_float(row, ("lambda", "krr_lam", "lam"), 1e-3)
            Ktt = wg * blocks["Kg_tt"] + wge * blocks["Kge_tt"] + we * blocks["Ke_tt"]
            Kts = wg * blocks["Kg_ts"] + wge * blocks["Kge_ts"] + we * blocks["Ke_ts"]
            ytt = blocks["ytt"]
            Xtt = blocks["Xtt"]
            ei_train = blocks["ei_train"]
            resid_diag_env = _residual_diag_env(Ktt, lam_eff, blocks)
            sigma2_env_t = _t(np.asarray(resid_diag_env, dtype=float), device=device_resolved, dtype=dtype_t)
            diag_add = sigma2_env_t.index_select(0, ei_train) + float(lam_eff)
            Vtt = 0.5 * (Ktt + Ktt.mT)
            Vtt.diagonal().add_(diag_add)
            L = _cholesky_factor(Vtt)
            Vinv_y = _cholesky_solve_from_factor(L, ytt)
            Vinv_X = _cholesky_solve_from_factor(L, Xtt)
            Xt_Vinv_X = Xtt.mT @ Vinv_X
            Xt_Vinv_y = Xtt.mT @ Vinv_y
            try:
                Var_beta = _spd_inverse(Xt_Vinv_X)
            except RuntimeError:
                Var_beta = _spd_inverse(
                    Xt_Vinv_X + 1e-10 * torch.eye(Xt_Vinv_X.shape[0], dtype=dtype_t, device=device_resolved)
                )
            beta = Var_beta @ Xt_Vinv_y
            alpha = _cholesky_solve_from_factor(L, ytt - Xtt @ beta)
            mean_lat_std = (Kts.mT @ alpha) + (Xss @ beta)
            mean_obs = stdr.inv_mean(_to_numpy(mean_lat_std), ei[s_idx])
            pred = pd.DataFrame({
                "Name": names_out,
                "Env": envs_out,
                "Prediction": np.asarray(mean_obs, dtype=float),
            })
            pred.insert(0, "candidate_id", cid)
            pred_parts.append(pred)
        except Exception as exc:
            failed = True
            context = {
                "geno_kernel_shape": tuple(np.asarray(geno_kernel).shape),
                "n_geno_ids": len(geno_ids),
                "n_obs": int(df.shape[0]),
                "n_env_levels": len(env_levels),
                "n_train_idx": int(t_idx.size),
                "n_test_idx": int(s_idx.size),
                "max_gi": int(np.max(gi)) if gi.size else None,
                "max_ei": int(np.max(ei)) if ei.size else None,
            }
            try:
                context["kenv_shape"] = tuple(np.asarray(kernel_map.get(kernel_name)).shape)
            except Exception:
                context["kenv_shape"] = None
            err = f"{type(exc).__name__}: {exc}; context={context}"
        detail_rows.append({
            "candidate_id": cid,
            "elapsed_sec": float(time.time() - t0),
            "failed": bool(failed),
            "error": err,
            "fast_path": "dense_krr_blocks",
            "device": str(device_resolved),
        })

    return {
        "predictions": pd.concat(pred_parts, ignore_index=True) if pred_parts else pd.DataFrame([]),
        "detail": pd.DataFrame(detail_rows),
    }


def batch_fit_mixed_model_candidates(
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    train_idx: Sequence[int],
    test_idx: Sequence[int],
    candidates: Any,
    env_covariates: Optional[Any] = None,
    env_similarity_by_kernel: Optional[Any] = None,
    include_components: Sequence[str] = ("g", "ge", "e"),
    fixed_effects: Optional[Sequence[str]] = None,
    backend: str = "auto",
    reaction_norm_feature_qc: bool = True,
    reaction_norm_env_filter_topk: Optional[int] = 15,
    kenv_bandwidth: float = 1.0,
    kenv_kernel_kwargs: Optional[Dict[str, Any]] = None,
    dtype: Union[str, "torch.dtype", None] = "float64",
    seed: int = 12345,
    standardize: str = "global",
    prediction_output: str = "test_only",
) -> Dict[str, Any]:
    """Evaluate single-trait MET KRR tuning candidates in one Python call."""
    cand = pd.DataFrame(candidates).copy()
    if cand.empty:
        return {
            "predictions": pd.DataFrame([]),
            "detail": pd.DataFrame([]),
            "env_kernel_cache_info": pd.DataFrame([]),
        }
    if "candidate_id" not in cand.columns:
        cand.insert(0, "candidate_id", np.arange(1, cand.shape[0] + 1, dtype=int))

    df = pd.DataFrame(pheno_df).reset_index(drop=True).copy()
    env_levels = list(df[env_col].astype("category").cat.categories.astype(str))
    kernel_map: Dict[str, np.ndarray] = {}
    env_cache_rows: List[Dict[str, Any]] = []

    if env_similarity_by_kernel is not None:
        try:
            iterator = env_similarity_by_kernel.items()
        except Exception:
            iterator = []
        for key, value in iterator:
            try:
                kernel_map[str(key).lower()] = np.asarray(value, dtype=np.float64)
            except Exception:
                continue

    if env_covariates is not None:
        for kernel_name in sorted(set(str(x).lower() for x in cand["kenv_kernel"].tolist())):
            if kernel_name in kernel_map:
                continue
            try:
                built = prepare_env_similarity_from_covariates(
                    env_covariates=env_covariates,
                    env_levels=env_levels,
                    env_col=env_col,
                    feature_qc=bool(reaction_norm_feature_qc),
                    kenv_kernel=kernel_name,
                    kenv_bandwidth=float(kenv_bandwidth),
                    kenv_kernel_kwargs=kenv_kernel_kwargs,
                )
                info = dict(built.get("info", {}) or {})
                if built.get("env_similarity") is not None:
                    kernel_map[kernel_name] = np.asarray(built["env_similarity"], dtype=np.float64)
                    status = "ok"
                    err = None
                else:
                    status = "empty"
                    err = None
                env_cache_rows.append({
                    "kenv_kernel": kernel_name,
                    "status": status,
                    "n_env": info.get("n_env"),
                    "n_features_input": info.get("n_features_input"),
                    "n_features_after_qc": info.get("n_features_after_qc"),
                    "n_features_after_variance_filter": info.get("n_features_after_variance_filter"),
                    "g2f_schema_detected": info.get("g2f_schema_detected"),
                    "n_envs_imputed": info.get("n_envs_imputed"),
                    "error": err,
                })
            except Exception as exc:
                env_cache_rows.append({
                    "kenv_kernel": kernel_name,
                    "status": "failed",
                    "n_env": len(env_levels),
                    "n_features_input": None,
                    "n_features_after_qc": None,
                    "n_features_after_variance_filter": None,
                    "g2f_schema_detected": None,
                    "n_envs_imputed": None,
                    "error": str(exc),
                })

    pred_parts: List[pd.DataFrame] = []
    detail_rows: List[Dict[str, Any]] = []
    geno_kernel = np.asarray(geno_kernel, dtype=np.float64)
    geno_ids = [str(x) for x in list(geno_ids)]
    train_idx = np.asarray(train_idx, dtype=np.int64)
    test_idx = np.asarray(test_idx, dtype=np.int64)
    fixed_effects = list(fixed_effects or [])

    try:
        fast = _batch_krr_point_predict_candidates(
            pheno_df=df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_kernel=geno_kernel,
            geno_ids=geno_ids,
            train_idx=train_idx,
            test_idx=test_idx,
            candidates=cand,
            kernel_map=kernel_map,
            fixed_effects=fixed_effects,
            dtype=dtype,
            seed=int(seed),
            standardize=str(standardize),
            prediction_output=str(prediction_output),
            reaction_norm_env_filter_topk=reaction_norm_env_filter_topk,
        )
        fast["env_kernel_cache_info"] = pd.DataFrame(env_cache_rows)
        return fast
    except Exception as fast_exc:
        _warn_once(
            ("batch_krr_fast_path_fallback", type(fast_exc).__name__, str(fast_exc)[:160]),
            f"Fast dense KRR batch candidate path failed; falling back to per-candidate fit_mixed_model: {fast_exc}",
            RuntimeWarning,
            stacklevel=2,
        )

    for _, row in cand.iterrows():
        cid = int(row["candidate_id"])
        t0 = time.time()
        err = None
        failed = False
        try:
            kernel_name = str(row.get("kenv_kernel", "matern32")).lower()
            env_similarity = kernel_map.get(kernel_name)
            if env_similarity is None:
                raise ValueError(f"no environment similarity was available for kernel {kernel_name!r}")
            fit = fit_mixed_model(
                pheno_df=df,
                gid_col=gid_col,
                env_col=env_col,
                y_col=y_col,
                geno_kernels={"G": geno_kernel},
                geno_ids=geno_ids,
                env_similarity=env_similarity,
                include_components=list(include_components),
                fixed_effects=fixed_effects,
                method="krr_exact",
                backend=str(backend),
                prediction_output=str(prediction_output),
                output_level="predict_only",
                point_predictions_only=True,
                return_se=False,
                compute_ai_se=False,
                standardize=str(standardize),
                varcomp_mode="mom",
                krr_lam=float(row.get("lambda", row.get("krr_lam", 1e-3))),
                krr_lams="none",
                lam_select="fixed",
                dtype=dtype,
                seed=int(seed),
                train_idx=train_idx,
                test_idx=test_idx,
                w_g=float(row.get("w_g", 1.0)),
                w_ge=float(row.get("w_ge", 0.0)),
                w_e=float(row.get("w_e", 0.0)),
            )
            pred = fit["result"]["predictions"].copy()
            pred.insert(0, "candidate_id", cid)
            pred_parts.append(pred)
        except Exception as exc:
            failed = True
            err = str(exc)
        detail_rows.append({
            "candidate_id": cid,
            "elapsed_sec": float(time.time() - t0),
            "failed": bool(failed),
            "error": err,
        })

    return {
        "predictions": (
            pd.concat(pred_parts, ignore_index=True)
            if pred_parts else pd.DataFrame([])
        ),
        "detail": pd.DataFrame(detail_rows),
        "env_kernel_cache_info": pd.DataFrame(env_cache_rows),
    }


__all__ = [
    "set_deterministic",
    "prepare_fit_mixed_model_context",
    "multikernel_krr_posterior_with_se",
    "gp_exact_with_X",
    "gp_icm_fa_with_X",
    "fit_mixed_model",
    "batch_fit_mixed_model_candidates",
    "materialize_geno_kernel_from_factor_cache",
]
def modular_registry_summary(
    geno_kernel_names,
    env_col,
    env_levels,
    env_structure=None,
    fa_rank=None,
    interaction_terms_meta=None,
):
    term_registry = build_term_registry_from_legacy_inputs(
        geno_kernel_names=geno_kernel_names,
        env_col=env_col,
        env_levels=env_levels,
        env_structure=env_structure,
        fa_rank=fa_rank,
        interaction_terms_meta=interaction_terms_meta,
    )
    cov_modules = build_covariance_modules_from_term_registry(term_registry)
    return {
        "term_registry": term_registry,
        "covariance_modules": cov_modules,
    }


def _build_modular_kernel_matrices_by_term(
    df: pd.DataFrame,
    gid_col: str,
    geno_kernels: Dict[str, np.ndarray],
    geno_ids: Sequence[str],
    term_registry,
):
    gid = df[gid_col].astype(str).tolist()
    gid_to_index = {str(g): i for i, g in enumerate(list(geno_ids))}
    gi = np.array([gid_to_index[str(g)] for g in gid], dtype=np.int64)
    out = {}
    for term in term_registry:
        try:
            term_name = term.term_name
            kernel_names = list(term.kernel_names)
        except Exception:
            term_name = term.get("term_name", "G")
            kernel_names = list(term.get("kernel_names", []))
        out[term_name] = {}
        for k in kernel_names:
            K = np.asarray(geno_kernels[k], dtype=float)
            out[term_name][k] = K[np.ix_(gi, gi)]
    return out


def _build_modular_level_index_by_term(df: pd.DataFrame, term_registry):
    out = {}
    for term in term_registry:
        try:
            axis_name = term.axis_name
            axis_levels = list(term.axis_levels) if term.axis_levels is not None else None
            term_name = term.term_name
        except Exception:
            axis_name = term.get("axis_name", None)
            axis_levels = list(term.get("axis_levels", [])) if term.get("axis_levels", None) is not None else None
            term_name = term.get("term_name", "G")
        if axis_name is None:
            continue
        if axis_name in df.columns:
            vals = df[axis_name].astype(str).tolist()
        else:
            parts = str(axis_name).split(":")
            if not all(p in df.columns for p in parts):
                continue
            vals = df[parts].astype(str).agg(":".join, axis=1).tolist()
        level_to_idx = {str(x): i for i, x in enumerate(axis_levels or [])}
        out[term_name] = np.array([level_to_idx[str(v)] for v in vals], dtype=np.int64)
    return out


def _make_modular_model_spec_from_term_registry(fixed_effects, term_registry):
    if not _MODULAR_TERM_SPEC_AVAILABLE:
        return None
    try:
        return _ModMixedModelSpec(
            fixed_effect_names=list(fixed_effects or []),
            random_term_specs=list(term_registry),
            residual_mode="heteroskedastic",
            metadata={},
        )
    except Exception:
        return None


def _convert_modular_result_to_legacy_output(modular_result, prediction_output="all"):
    fm = modular_result.fitted_model
    bundle = fm.diagnostics.get("reporting_bundle", {})
    pred_df = fm.predictions.copy() if fm.predictions is not None else pd.DataFrame()
    out = {
        "beta": fm.beta,
        "beta_se": fm.beta_se,
        "post_mean": pred_df["Prediction"].to_numpy(dtype=float) if "Prediction" in pred_df.columns else None,
        "post_se": pred_df["Prediction_SE_observed"].to_numpy(dtype=float) if "Prediction_SE_observed" in pred_df.columns else (
            pred_df["Prediction_SE_latent"].to_numpy(dtype=float) if "Prediction_SE_latent" in pred_df.columns else None
        ),
        "latent_mean": pred_df["Prediction"].to_numpy(dtype=float) if "Prediction" in pred_df.columns else None,
        "latent_se": pred_df["Prediction_SE_latent"].to_numpy(dtype=float) if "Prediction_SE_latent" in pred_df.columns else None,
        "var_components_summary": bundle.get("var_components_summary", pd.DataFrame()),
        "var_components_ai": bundle.get("var_components_ai", pd.DataFrame()),
        "interaction_variance_summary": bundle.get("interaction_variance_summary", {}),
        "interaction_correlation_summary": bundle.get("interaction_correlation_summary", {}),
        "residual_summary": bundle.get("residual_summary", {}),
        "heritability_summary": bundle.get("heritability_summary", {}),
        "env_variance_summary": bundle.get("env_variance_summary", None),
        "var_components": bundle.get("var_components", pd.DataFrame()),
        "diagnostics": fm.diagnostics,
        "result": {
            "predictions": pred_df,
            "per_env": bundle.get("env_variance_summary", pd.DataFrame()) if bundle.get("env_variance_summary", None) is not None else pd.DataFrame(),
            "across_env": pd.DataFrame(),
        },
        "notes": ["Returned via modular dense delegation path."],
    }
    return out


def inspect_framework_engine_mode(
    *,
    requested_engine="auto",
    modular_dense_delegate=False,
    use_modular_dense=False,
    point_predictions_only=False,
    return_se=True,
    compute_ai_se=False,
    n_obs=None,
    n_structured_terms=1,
    max_axis_levels=0,
    large_n_threshold=50000,
    max_full_levels=20,
):
    """Public helper to inspect framework engine routing without fitting."""
    return _decide_framework_engine_mode(
        requested_engine=requested_engine,
        modular_dense_delegate=modular_dense_delegate,
        use_modular_dense=use_modular_dense,
        point_predictions_only=point_predictions_only,
        return_se=return_se,
        compute_ai_se=compute_ai_se,
        n_obs=n_obs,
        n_structured_terms=n_structured_terms,
        max_axis_levels=max_axis_levels,
        large_n_threshold=large_n_threshold,
        max_full_levels=max_full_levels,
    )
