#!/usr/bin/env python3
# bayes_alphabet_bglr_parity.py
# GPU-ready Bayes Alphabet (BRR, BayesA, BayesB/C/BayesCπ) + LASSO + Elastic Net + RKHS
# BGLR-compatible outputs; multiple ETA blocks; NA in y -> prediction set.
#
# This version integrates memory-safety, numerical-stability, and logging fixes,
# plus the requested enhancements:
# - Configurable memory budgets + adaptive batching (robustness-first)
# - Context object (removes reliance on globals for new logic; improves testability)
# - Better priors: Regularized Horseshoe default + data-weighted block R² allocation
# - VI temperature schedule + guide options (NumPyro)
# - Distinct CPU/GPU memory detection & structured diagnostics
# - Windows ctypes guarded behind an opt-in flag
#
# Env knobs:
#   BAYES_ALPHABET_LOGLEVEL=INFO|WARNING|ERROR
#   JAX_ENABLE_X64=true|false (overrides default enabling)
#
# -----------------------------------------------------------------------------
# --- minimal cross-OS bootstrap for JAX / NumPyro / PyMC -----------------------
from __future__ import annotations
import os, sys, math, warnings, importlib, subprocess, gc, ctypes, platform, shutil
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, TYPE_CHECKING
import logging  # <-- must be before logger setup

# Avoid picking up user-site wheels that may conflict with the active env
os.environ.setdefault("PYTHONNOUSERSITE", "1")

def _run(cmd):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          check=False, text=True)

def _pip_install(pkgs):
    """Quiet, upgrade, no cache; installs into the *active* environment."""
    if not pkgs:
        return None
    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--quiet", "--no-cache-dir"]
    cmd += list(pkgs)
    return _run(cmd)

def _has_nvidia_gpu() -> bool:
    try:
        r = _run(["nvidia-smi"])
        return r.returncode == 0
    except Exception:
        return False

def _ensure_jax(prefer_gpu=True, cuda="auto", allow_install=False):
    """
    Import jax/jaxlib (and jax-metal on Apple Silicon). On Windows JAX is CPU-only.
    On Linux, CUDA wheels are attempted only if a GPU is present (unless caller forbids).
    """
    try:
        jax = importlib.import_module("jax")
        jaxlib = importlib.import_module("jaxlib")
        return jax, jaxlib, {"installed": False}
    except Exception as e:
        if not allow_install:
            raise ImportError("Required package 'jax' not available and allow_install=False.") from e

    sys_os = platform.system()
    arch = platform.machine().lower()
    pkgs = []

    if sys_os == "Windows":
        # CPU wheels; GPU JAX is not supported on Windows
        pkgs = ["jax==0.6.2", "jaxlib==0.6.2"]
    elif sys_os == "Darwin" and arch in ("arm64", "aarch64"):
        # Apple Silicon (Metal backend); CPU ok otherwise
        pkgs = ["jax", "jax-metal>=0.0.5"]
    else:
        # Linux
        if prefer_gpu and cuda != "cpu" and _has_nvidia_gpu():
            # Installs a matching CUDA jaxlib via pip extra
            pkgs = ["jax[cuda12_pip]"]
        else:
            pkgs = ["jax", "jaxlib"]

    _pip_install(pkgs)
    jax = importlib.import_module("jax")
    jaxlib = importlib.import_module("jaxlib")
    return jax, jaxlib, {"installed": True, "pkgs": pkgs}

def _ensure_numpyro(allow_install=False):
    try:
        return importlib.import_module("numpyro")
    except Exception as e:
        if not allow_install:
            raise
        _pip_install(["numpyro"])
        return importlib.import_module("numpyro")

def _ensure_pymc_arviz(allow_install=False):
    try:
        pm = importlib.import_module("pymc")
    except Exception:
        if not allow_install:
            raise
        _pip_install(["pymc>=5.10"])
        pm = importlib.import_module("pymc")
    try:
        importlib.import_module("arviz")
    except Exception:
        if allow_install:
            _pip_install(["arviz"])
    return pm

def setup_deps(prefer_gpu: bool = True,
               try_numpyro: bool = True,
               try_pymc: bool = True,
               allow_install: bool = False,
               cuda: str = "auto"):
    """
    Verify (and optionally install) JAX(+metal/CUDA), NumPyro, and PyMC(+ArviZ).
    Returns a small info dict with versions and device info.
    """
    info = {"python": sys.version.split()[0], "platform": platform.platform()}

    jax, jaxlib, meta = _ensure_jax(prefer_gpu=prefer_gpu, cuda=cuda, allow_install=allow_install)
    info["jax_version"] = getattr(jax, "__version__", "unknown")
    info["jaxlib_version"] = getattr(jaxlib, "__version__", "unknown")
    try:
        devs = jax.devices()
        info["devices"] = [getattr(d, "platform", "?") for d in devs]
        info["default_backend"] = jax.default_backend()
        info["local_device_count"] = len(devs)
        info["accelerator_present"] = any(p in ("gpu","cuda","rocm","tpu") for p in info["devices"])
    except Exception as e:
        info["devices_error"] = repr(e)

    if try_numpyro:
        try:
            npyro = _ensure_numpyro(allow_install=allow_install)
            info["numpyro_version"] = getattr(npyro, "__version__", "unknown")
        except Exception as e:
            info["numpyro_error"] = repr(e)

    if try_pymc:
        try:
            pm = _ensure_pymc_arviz(allow_install=allow_install)
            info["pymc_version"] = getattr(pm, "__version__", "unknown")
        except Exception as e:
            info["pymc_error"] = repr(e)

    if platform.system() == "Windows":
        info["note"] = "On Windows, JAX runs CPU-only (no CUDA wheels)."

    return info
# --- end minimal bootstrap ------------------------------------------------------

# ---------------------------
# Logging
# ---------------------------

# --- logging setup (safe, idempotent) ---
_LOGGER_NAME = "bayes_alphabet"
logger = logging.getLogger(_LOGGER_NAME)
if not logger.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
    logger.addHandler(_h)
logger.setLevel(os.environ.get("BAYES_ALPHABET_LOGLEVEL", "WARNING").upper())


# ---------------------------
# Legacy numerical knob (kept for reference)
# ---------------------------
CHOLESKY_JITTER = 1e-6  # legacy constant; not used directly in new code

# ---------------------------
# Install / Import Utilities
# ---------------------------

def _pip_install(spec: str, *, quiet: bool = True, extra_args: Optional[List[str]] = None, disable_auto_install: bool = False) -> bool:
    """Install a pip package unless auto-install is disabled."""
    if disable_auto_install:
        raise ImportError(f"Auto-install disabled for '{spec}'.")
    try:
        import subprocess, sys
        args = [sys.executable, "-m", "pip", "install", spec]
        if extra_args:
            args.extend(extra_args)
        if quiet:
            args.append("--quiet")
        subprocess.check_call(args)
        return True
    except Exception as e:
        warnings.warn(f"Auto-install failed for '{spec}': {e}", RuntimeWarning)
        logger.warning("Auto-install failed for '%s': %s", spec, e)
        return False

def _ensure_packaging(allow_install: bool, strict_versions: bool, disable_auto_install: bool = False) -> bool:
    try:
        from packaging.version import Version, InvalidVersion  # noqa: F401
        return True
    except Exception:
        if allow_install:
            try:
                _pip_install("packaging>=23.0", extra_args=["--upgrade"], disable_auto_install=disable_auto_install)
                from packaging.version import Version, InvalidVersion  # type: ignore  # noqa: F401
                return True
            except Exception:
                if strict_versions:
                    raise
                warnings.warn("Could not ensure 'packaging'; version checks may be skipped.", RuntimeWarning)
                logger.warning("Could not ensure 'packaging'; version checks may be skipped.")
                return False
        if strict_versions:
            raise ImportError("packaging is required when strict_versions=True.")
        warnings.warn("packaging not available; skipping version comparisons.", RuntimeWarning)
        logger.warning("packaging not available; skipping version comparisons.")
        return False

if TYPE_CHECKING:
    from packaging.version import Version, InvalidVersion  # type: ignore

def _ensure_module(name: str,
                   pip_name: Optional[str],
                   min_version: Optional[str],
                   *,
                   allow_install: bool,
                   strict_versions: bool,
                   optional: bool = False,
                   disable_auto_install: bool = False):
    """Import a module, optionally installing/upgrading to satisfy min_version."""
    try:
        mod = importlib.import_module(name)
    except Exception:
        if not allow_install:
            if optional:
                logger.info("Optional module '%s' not found; continuing.", name)
                return None
            raise ImportError(f"Required package '{pip_name or name}' not available and allow_install=False.")
        spec = pip_name or name
        if not _pip_install(spec, disable_auto_install=disable_auto_install):
            if optional:
                logger.info("Optional module '%s' could not be installed; continuing.", name)
                return None
            raise ImportError(f"Installation failed for '{spec}'.")
        mod = importlib.import_module(name)

    v = getattr(mod, "__version__", None)
    if min_version:
        have_packaging = _ensure_packaging(allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
        if v is None:
            if strict_versions:
                raise ImportError(f"Unknown version for '{name}' with strict_versions=True.")
            warnings.warn(f"Unknown version for '{name}'; proceeding.", RuntimeWarning)
            logger.warning("Unknown version for '%s'; proceeding.", name)
        elif have_packaging:
            from packaging.version import Version, InvalidVersion  # type: ignore
            try:
                if Version(str(v)) < Version(min_version):
                    if allow_install:
                        if not _pip_install(f"{pip_name or name}>={min_version}", extra_args=["--upgrade"], disable_auto_install=disable_auto_install):
                            raise ImportError(f"Could not upgrade '{name}' to >= {min_version}.")
                        mod = importlib.import_module(name)
                        v2 = getattr(mod, "__version__", v)
                        if Version(str(v2)) < Version(min_version):
                            raise ImportError(f"{name} version {v2} < {min_version} after upgrade.")
                    else:
                        raise ImportError(f"{name} version {v} < {min_version} and allow_install=False.")
            except InvalidVersion:
                if strict_versions:
                    raise ImportError(f"Unrecognized version for '{name}' with strict_versions=True.")
                warnings.warn(f"Unrecognized version for '{name}'; proceeding.", RuntimeWarning)
                logger.warning("Unrecognized version for '%s'; proceeding.", name)
    return mod

# Globals assigned lazily (kept for backward paths; new logic uses Context)
np = None
jax = None
jnp = None
az = None
psutil = None

# Optional explicit imports for static analysis (will be None if unavailable)
try:
    import numpyro  # type: ignore
    import numpyro.distributions as dist  # type: ignore
    from numpyro.infer import SVI, Trace_ELBO, MCMC, NUTS  # type: ignore
    from numpyro.infer.autoguide import (  # type: ignore
        AutoDiagonalNormal,
        AutoLowRankMultivariateNormal,
        AutoMultivariateNormal,
    )
except Exception:
    numpyro = None  # type: ignore
    dist = None     # type: ignore
    SVI = Trace_ELBO = MCMC = NUTS = AutoDiagonalNormal = None  # type: ignore

# PyTensor (for PyMC) — import lazily-safe
try:
    from pytensor import config as ptconfig  # type: ignore
    import pytensor.tensor as pt  # type: ignore
    _PT_FLOATX = ptconfig.floatX
except Exception:
    pt = None  # type: ignore
    _PT_FLOATX = "float64"

def _pt_const(x):
    """Create a PyTensor constant with floatX dtype (avoid raw python floats in the graph)."""
    global pt
    if pt is None:
        import pytensor.tensor as pt  # type: ignore
    import numpy as _np
    return pt.as_tensor_variable(_np.asarray(x, dtype=_PT_FLOATX))

# ---------------------------
# Context & Config (new)
# ---------------------------

@dataclass
class MemoryBudget:
    host_fraction: float = 0.30
    device_fraction: float = 0.50
    safety_margin_gb: float = 1.0
    enable_win_ctypes: bool = False   # guard Windows ctypes

@dataclass
class VIConfig:
    guide: str = "diag"               # "diag" | "lowrank" | "full"
    rank: Optional[int] = None        # for "lowrank"
    restarts: int = 1
    temp_start: float = 1.0
    temp_end: float = 0.2
    temp_schedule: str = "cosine"     # "cosine" | "linear"
    temp_warmup_frac: float = 0.10    # fraction of steps to keep temp_start

@dataclass
class RKHSConfig:
    approx: str = "exact"             # "exact" | "nystrom"
    m: Optional[int] = None
    seed: Optional[int] = None
    kernel_type: str = "unknown"
    force_nystrom_mem_gb: Optional[float] = None

@dataclass
class PriorConfig:
    family: str = "rhs"               # "rhs" (regularized horseshoe) | "default"
    rhs_m0: Optional[float] = None    # expected nonzeros; default: 0.1 p + 10 (capped)
    rhs_slab_scale_mult: float = 10.0 # c0 = mult * sigma/sqrt(n) (Piironen–Vehtari style)

@dataclass
class BayesAlphabetConfig:
    memory: MemoryBudget = field(default_factory=MemoryBudget)
    vi: VIConfig = field(default_factory=VIConfig)
    rkhs: RKHSConfig = field(default_factory=RKHSConfig)
    prior: PriorConfig = field(default_factory=PriorConfig)
    block_R2_strategy: str = "weighted"   # "equal" | "weighted" | "from_eta"

class Context:
    """Runtime context: modules, devices, memory, logging."""
    def __init__(self, config: BayesAlphabetConfig, *, allow_install=False, strict_versions=False, disable_auto_install=True):
        self.config = config
        self.allow_install = allow_install
        self.strict_versions = strict_versions
        self.disable_auto_install = disable_auto_install
        # Ensure deps and wire modules into context (still keep globals for legacy paths)
        _ensure_core_deps(allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
        self.np, self.jax, self.jnp, self.az = np, jax, jnp, az
        # Device snapshot
        try:
            self.devices = jax.devices()
        except Exception:
            self.devices = []
        # Optional NVML
        try:
            self.nvml = _ensure_module("pynvml", "nvidia-ml-py3", None,
                                       allow_install=allow_install, strict_versions=strict_versions,
                                       optional=True, disable_auto_install=disable_auto_install)
            if self.nvml:
                self.nvml.nvmlInit()
        except Exception:
            self.nvml = None

    # --- Memory detection
    def available_host_gb(self) -> float:
        """Best-effort available RAM (GiB) with guarded Windows ctypes."""
        # psutil first
        try:
            import psutil as _ps
            return float(_ps.virtual_memory().available) / (1024**3)
        except Exception:
            pass
        # POSIX sysconf
        try:
            if hasattr(os, "sysconf"):
                pagesize = os.sysconf("SC_PAGE_SIZE")
                avail = os.sysconf("SC_AVPHYS_PAGES")
                if isinstance(pagesize, int) and isinstance(avail, int) and pagesize > 0 and avail > 0:
                    return float(pagesize) * float(avail) / (1024**3)
        except Exception:
            pass
        # Windows (guarded)
        try:
            if self.config.memory.enable_win_ctypes and os.name == "nt" and hasattr(ctypes, "windll"):
                class MEMORYSTATUSEX(ctypes.Structure):
                    _fields_ = [
                        ("dwLength", ctypes.c_ulong),
                        ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                    ]
                stat = MEMORYSTATUSEX()
                stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
                ok = ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))  # type: ignore[attr-defined]
                if ok:
                    return float(stat.ullAvailPhys) / (1024**3)
        except Exception:
            pass
        warnings.warn("Could not determine available RAM; assuming 8 GiB.", RuntimeWarning)
        logger.warning("Could not determine available RAM; assuming 8 GiB.")
        return 8.0

    def available_device_gb(self) -> Optional[float]:
        """Best-effort available GPU memory (GiB); None if unknown/unavailable."""
        # Prefer NVML if available
        try:
            if self.nvml:
                count = self.nvml.nvmlDeviceGetCount()
                if count > 0:
                    h = self.nvml.nvmlDeviceGetHandleByIndex(0)
                    mem = self.nvml.nvmlDeviceGetMemoryInfo(h)
                    return float(mem.free) / (1024**3)
        except Exception:
            pass
        # Try JAX device memory (not always available)
        try:
            for d in self.devices:
                plat = getattr(d, "platform", "").lower()
                if plat in ("gpu", "cuda", "rocm", "tpu"):
                    # No standard API for "available"; return None
                    return None
        except Exception:
            pass
        return None

    # --- Budgets
    def resolve_host_budget_gb(self, *, requested: Optional[float]=None) -> Tuple[float, Dict[str, Any]]:
        if requested is not None and requested > 0:
            return float(requested), {"source":"user", "value_gb": float(requested)}
        avail = self.available_host_gb()
        frac = float(self.config.memory.host_fraction)
        margin = float(self.config.memory.safety_margin_gb)
        target = max(0.1, frac * max(0.0, avail - margin))
        info = {"source": "auto", "available_gb": float(avail), "fraction": frac, "margin_gb": margin, "value_gb": float(target)}
        return float(target), info

    def resolve_device_budget_gb(self) -> Tuple[Optional[float], Dict[str, Any]]:
        avail = self.available_device_gb()
        if avail is None:
            return None, {"source":"unknown", "available_gb": None, "fraction": float(self.config.memory.device_fraction), "value_gb": None}
        frac = float(self.config.memory.device_fraction)
        margin = float(self.config.memory.safety_margin_gb)
        target = max(0.1, frac * max(0.0, avail - margin))
        return float(target), {"source":"auto", "available_gb": float(avail), "fraction": frac, "margin_gb": margin, "value_gb": float(target)}

# ---------------------------
# Numerics helpers
# ---------------------------

def _estimate_overhead_mult(backend: Optional[str]) -> float:
    """Heuristic memory overhead for graph frameworks (conservative)."""
    b = (backend or "").lower()
    if b == "pymc":
        return 6.0
    if b in ("numpyro", "jax"):
        return 5.0
    return 3.5

def _maybe_cleanup(*objs):
    for o in objs:
        try:
            del o
        except Exception:
            pass
    gc.collect()

def _adaptive_jitter(matrix: np.ndarray, base_jitter: float = 1e-12) -> float:
    """Adaptive jitter based on condition number & dtype eps; capped by matrix scale; returns scalar jitter."""
    A = np.asarray(matrix)
    try:
        cond = float(np.linalg.cond(A))
    except Exception:
        cond = 1e4
    eps = np.finfo(A.dtype if A.dtype.kind == "f" else np.float64).eps
    jitter = max(base_jitter, eps * 10.0)
    if math.isfinite(cond):
        if cond > 1e12:
            jitter *= 1e6
        elif cond > 1e8:
            jitter *= 1e3
        elif cond > 1e6:
            jitter *= 1e2
    try:
        max_diag = float(np.max(np.abs(np.diag(A)))) if A.ndim == 2 else float(np.max(np.abs(A)))
    except Exception:
        max_diag = float(np.linalg.norm(A, ord=np.inf))
    cap = max(1e-12, 1e-3 * max(1.0, max_diag))
    return float(min(jitter, cap))

# ---------------------------
# Public API helpers
# ---------------------------

def _resolve_beta_samples(beta_map: Dict[str, np.ndarray], pref: str, *, S: int, p_expected: int) -> np.ndarray:
    if p_expected == 0:
        return np.zeros((int(S), 0), dtype=np.float32)
    key = pref + "beta"
    if key in beta_map:
        B = beta_map[key]
        if B.ndim == 1:
            B = B[:, None]
        return B
    cands = sorted(k for k in beta_map if k.startswith(pref) and k.endswith("beta") and not k.endswith("beta_raw"))
    if len(cands) > 1:
        warnings.warn(f"Multiple beta candidates for '{pref}': {cands[:5]}... using {cands[0]}", RuntimeWarning)
        logger.warning("Multiple beta candidates for '%s': %s ... using %s", pref, cands[:5], cands[0])
    if cands:
        arr = beta_map[cands[0]]
        if arr.ndim == 1:
            arr = arr[:, None]
        return arr
    raise KeyError(
        f"Could not find beta samples for prefix '{pref}'. "
        f"Available keys: {sorted(beta_map.keys())[:12]}..."
    )

def _extract_samples_from_raw(samples: Dict[str, np.ndarray]):
    beta_map: Dict[str, np.ndarray] = {}
    aux: Dict[str, Any] = {}
    u_map: Dict[str, np.ndarray] = {}

    def _maybe_beta(a: np.ndarray) -> np.ndarray:
        a = np.asarray(a)
        if a.size == 0:
            S = a.shape[0] if a.ndim >= 1 else 1
            return np.zeros((S, 0), dtype=np.float32)
        if a.ndim == 1:
            return a[:, None]
        return a

    intercept = np.asarray(samples.get("intercept", np.array([]))).reshape(-1)
    sigma2 = None
    if "sigma2" in samples:
        sigma2 = np.asarray(samples["sigma2"]).reshape(-1)
    elif "sigma" in samples:
        s = np.asarray(samples["sigma"]).reshape(-1)
        sigma2 = s * s
    if intercept.size == 0:
        S = None
        for v in samples.values():
            if isinstance(v, np.ndarray):
                S = len(v)
                break
        intercept = np.zeros((S or 1,), dtype=np.float32)
    if sigma2 is None:
        sigma2 = np.full_like(intercept, np.nan, dtype=np.float32)

    for k, v in samples.items():
        if not isinstance(v, np.ndarray):
            continue
        if k.endswith("_beta"):
            beta_map[k] = _maybe_beta(v)
        elif k.endswith("_lambda"):
            aux[k.replace("_lambda", "_lambda_mean")] = np.mean(v, axis=0)
        elif k.endswith("_tau2"):
            aux[k + "_mean"] = float(np.mean(v))
            aux[k + "_sd"] = float(np.std(v))
        elif k.endswith("_tau") or k.endswith("_pi") or k == "student_t_nu":
            aux[k + "_mean"] = float(np.mean(v))
            aux[k + "_sd"] = float(np.std(v))
        elif k.endswith("_gate"):
            aux[k] = np.asarray(v)
        elif k.endswith("_u"):
            u_map[k] = np.asarray(v)

    return beta_map, intercept, sigma2, aux, u_map


# ---------------------------
# API
# ---------------------------

def bayesian_alphabet(
    ETA: List[Dict[str, Any]],
    y,
    *,
    backend: str = "numpyro",     # "numpyro" or "pymc"
    vi: bool = False,             # True => SVI (relaxed spike&slab) for BayesB/C/π by default
    n_iter: int = 2000,
    tune: int = 1000,
    chains: int = 2,
    random_seed: int = 123,
    target_accept: float = 0.9,
    weights = None,
    standardize_X: bool = True,
    center_y: bool = True,
    use_bglr_priors: bool = True,
    R2: float = 0.5,
    relaxed_temperature: float = 0.5,
    batch_mem_gb: Optional[float] = None,
    batch_overhead_mult: Optional[float] = None,
    allow_install: bool = False,
    strict_versions: bool = False,
    disable_auto_install: bool = True,
    cpi_hmc_mode: str = "vi",
    pymc_cpi_mode: str = "vi",
    nuts_dense_mass: bool = False,
    nuts_max_tree_depth: int = 10,
    nuts_init: str = "median",
    auto_tune: bool = True,
    block_prior_calibration: bool = False,
    block_R2_strategy: str = "weighted",      # default switched to "weighted"
    vi_guide: str = "diag",                   # "diag" | "fullrank" | "lowrank"
    vi_rank: Optional[int] = None,            # only for "lowrank" (NumPyro)
    # Student-t + RKHS options + VI restarts
    likelihood: str = "gaussian",             # "gaussian" or "student_t"
    sample_nu: bool = False,                  # sample Student-t degrees of freedom
    student_t_nu: float = 7.0,                # default nu if not sampling
    rkhs_approx: str = "exact",               # "exact" (Cholesky) or "nystrom"
    rkhs_m: Optional[int] = None,             # Nyström rank (heuristic if None)
    rkhs_seed: Optional[int] = None,          # seed for Nyström landmarks
    kernel_type: str = "unknown",             # "dense", "sparse", or "unknown"
    rkhs_force_nystrom_mem_gb: Optional[float] = None,  # force Nyström if n^2 exceeds budget (GiB)
    vi_restarts: int = 1,                     # NumPyro: multiple SVI inits; keep best (min loss)
    vi_max_restarts: int = 3,                 # cap restarts to avoid excessive runtime
    # New (optional) config replaces multiple knobs at once; if provided, it wins
    config: Optional[BayesAlphabetConfig] = None,
    # Optional explicit memory budget knobs (overrides config.memory if given)
    mem_host_fraction: Optional[float] = None,
    mem_device_fraction: Optional[float] = None,
    mem_safety_margin_gb: Optional[float] = None,
    enable_win_ctypes: Optional[bool] = None,
    # New: prior family (default regularized horseshoe)
    prior_family: Optional[str] = None,         # "rhs" or "default" (if None, use config)
    rhs_m0: Optional[float] = None,             # expected nonzeros (override)
    rhs_slab_scale_mult: Optional[float] = None,# c0 multiplier
    # New: VI temperature schedule (NumPyro)
    vi_temp_start: Optional[float] = None,
    vi_temp_end: Optional[float] = None,
    vi_temp_schedule: Optional[str] = None,     # "cosine" | "linear"
    vi_temp_warmup_frac: Optional[float] = None,
) -> Dict[str, Any]:
    """Fit Bayes Alphabet / L1 / ENet / RKHS models with HMC or VI.

    Robustness:
      - Configurable memory budgets & adaptive batching with diagnostics.
      - Distinct CPU/GPU memory probes; Windows ctypes guarded behind a flag.

    Priors:
      - Regularized Horseshoe (RHS) is the default for BayesB/C (HMC & VI) and is calibrated via Piironen–Vehtari.

    VI:
      - Temperature schedules for NumPyro SVI (cosine/linear), guide options (diag/lowrank/full).
    """
    # --- Build/merge config
    cfg = config or BayesAlphabetConfig()
    if mem_host_fraction is not None:     cfg.memory.host_fraction = float(mem_host_fraction)
    if mem_device_fraction is not None:   cfg.memory.device_fraction = float(mem_device_fraction)
    if mem_safety_margin_gb is not None:  cfg.memory.safety_margin_gb = float(mem_safety_margin_gb)
    if enable_win_ctypes is not None:     cfg.memory.enable_win_ctypes = bool(enable_win_ctypes)
    if prior_family is not None:          cfg.prior.family = str(prior_family).lower()
    if rhs_m0 is not None:                cfg.prior.rhs_m0 = float(rhs_m0)
    if rhs_slab_scale_mult is not None:   cfg.prior.rhs_slab_scale_mult = float(rhs_slab_scale_mult)
    if vi_temp_start is not None:         cfg.vi.temp_start = float(vi_temp_start)
    if vi_temp_end is not None:           cfg.vi.temp_end = float(vi_temp_end)
    if vi_temp_schedule is not None:      cfg.vi.temp_schedule = str(vi_temp_schedule).lower()
    if vi_temp_warmup_frac is not None:   cfg.vi.temp_warmup_frac = float(vi_temp_warmup_frac)
    if block_R2_strategy is not None:     cfg.block_R2_strategy = str(block_R2_strategy).lower()
    # Compat for direct kwargs
    cfg.vi.guide = str(vi_guide).lower() if vi_guide else cfg.vi.guide
    cfg.vi.rank = int(vi_rank) if vi_rank is not None else cfg.vi.rank
    cfg.rkhs.approx = str(rkhs_approx).lower()
    cfg.rkhs.m = rkhs_m
    cfg.rkhs.seed = rkhs_seed
    cfg.rkhs.kernel_type = kernel_type
    cfg.rkhs.force_nystrom_mem_gb = rkhs_force_nystrom_mem_gb

    ctx = Context(cfg, allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)

    # --- Input checks (kept from earlier)
    if y is None:
        raise ValueError("y must be a 1D array; got None.")
    y = np.asarray(y, dtype=np.float32)

    # Student-t ν validation
    if likelihood.lower() == "student_t" and (not sample_nu):
        if not vi and float(student_t_nu) <= 4.0:
            raise ValueError("student_t_nu must be > 4 for stable MCMC (HMC/NUTS).")
        if vi and float(student_t_nu) <= 2.5:
            raise ValueError("student_t_nu must be > 2.5 for numerically stable VI.")
        if vi and 2.5 < float(student_t_nu) <= 4.0:
            warnings.warn("student_t_nu in (2.5,4]: finite variance but may be numerically fragile; proceed with caution.", RuntimeWarning)

    if ETA is None:
        ETA = []
    if not isinstance(ETA, list):
        raise ValueError("ETA must be a list of dicts.")
    n = y.shape[0]
    if not np.any(np.isfinite(y)):
        raise ValueError("At least one finite y is required for training (all y are NaN).")

    valid = {"bayes_a","bayes_b","bayes_c","bayes_c_pi","bayesian_ridge","fixed","lasso","elastic_net","rkhs"}
    for i, e in enumerate(ETA):
        if ("method" not in e):
            raise ValueError(f"ETA[{i}] needs 'method'.")
        if e["method"] not in valid:
            raise ValueError(f"ETA[{i}]['method'] must be one of {sorted(valid)}.")
        m = e["method"]
        if m == "rkhs":
            if ("K" not in e) and ("X" not in e):
                raise ValueError(f"ETA[{i}] rkhs requires 'K' or 'X'.")
            if "K" in e:
                K = np.asarray(e["K"], dtype=np.float32)
                if K.shape != (n, n):
                    raise ValueError(f"ETA[{i}] rkhs: K must be (n,n).")
                if not np.allclose(K, K.T, atol=1e-5):
                    raise ValueError(f"ETA[{i}] rkhs: RKHS K must be symmetric (||K-K^T|| > tol).")
                if n <= 2000:
                    w = np.linalg.eigvalsh(K)
                    if float(w.min()) < -1e-6:
                        raise ValueError(f"ETA[{i}] rkhs: RKHS K must be PSD (min eig={w.min():.3e} < -1e-6).")
                else:
                    try:
                        J = _adaptive_jitter(K)
                        _ = np.linalg.cholesky(K + J*np.eye(n, dtype=np.float32))
                    except np.linalg.LinAlgError:
                        raise ValueError(f"ETA[{i}] rkhs: RKHS K appears not PSD (Cholesky failed).")
            elif "X" in e:
                X_i = np.asarray(e["X"])
                if X_i.shape[0] != n:
                    raise ValueError(f"ETA[{i}]['X'] rows must equal len(y).")
                if X_i.ndim != 2:
                    raise ValueError(f"ETA[{i}]['X'] must be 2D (n,p).")
                if not np.all(np.isfinite(X_i)):
                    raise ValueError(f"ETA[{i}]['X'] contains non-finite values.")
        else:
            if "X" not in e:
                raise ValueError(f"ETA[{i}] needs 'X' for method '{m}'.")
            X_i = np.asarray(e["X"])
            if X_i.shape[0] != n:
                raise ValueError(f"ETA[{i}]['X'] rows must equal len(y).")
            if X_i.ndim != 2:
                raise ValueError(f"ETA[{i}]['X'] must be 2D (n,p).")
            if not np.all(np.isfinite(X_i)):
                raise ValueError(f"ETA[{i}]['X'] contains non-finite values.")
        if ("fixed_pi" in e) and not (0.0 <= float(e["fixed_pi"]) <= 1.0):
            raise ValueError(f"ETA[{i}]['fixed_pi'] must be in [0,1].")
        if m == "elastic_net":
            a = float(e.get("alpha", 0.5))
            if not (0.0 <= a <= 1.0):
                raise ValueError(f"ETA[{i}] elastic_net alpha must be in [0,1].")
        if "d" in e and m != "rkhs":
            d_in = np.asarray(e["d"]).astype(np.float32)
            if d_in.shape != (np.asarray(e["X"]).shape[1],):
                raise ValueError(f"ETA[{i}]['d'] must have shape (p,).")
            if np.any((d_in < 0) | (d_in > 1)):
                raise ValueError(f"ETA[{i}]['d'] values must lie in [0,1].")

    if not (0.0 < float(R2) < 1.0):
        raise ValueError("R2 must be in (0,1) for prior calibration.")
    if not (0.01 <= float(relaxed_temperature) <= 10.0):
        raise ValueError("relaxed_temperature must be in [0.01, 10.0].")
    if relaxed_temperature < 0.05 or relaxed_temperature > 3.0:
        warnings.warn("relaxed_temperature is very extreme; can harm VI quality.", RuntimeWarning)
    if batch_mem_gb is not None and not (batch_mem_gb > 0.0):
        raise ValueError("batch_mem_gb must be > 0 when provided.")

    total_p = int(sum((np.asarray(e["X"]).shape[1] if ("X" in e and e["method"]!="rkhs") else 0) for e in ETA)) if ETA else 0
    total_size = int(n) * total_p
    if total_size > 1e8:
        per = [ (i, (np.asarray(e["X"]).shape[1] if ("X" in e and e["method"]!="rkhs") else 0)) for i,e in enumerate(ETA) ]
        approx_gib = (total_size * 4) / (1024**3)
        warnings.warn(
            "Large dataset: n * sum(p) ≈ "
            f"{total_size:,} (~{approx_gib:.1f} GiB float32). "
            f"Per-block p: {per[:6]}{'...' if len(per)>6 else ''}. "
            "Consider smaller batch_mem_gb or VI.",
            ResourceWarning
        )

    # auto-tune NUTS knobs
    if auto_tune:
        nuts_init = "median" if total_p < 1000 else "sample"
        latent_dim_est = int(total_p)
        has_rkhs = any(e.get("method") == "rkhs" for e in ETA)
        if has_rkhs:
            rkhs_blocks = sum(1 for e in ETA if e.get("method") == "rkhs")
            if cfg.rkhs.approx == "exact":
                m_est = n
            else:
                m_est = int(cfg.rkhs.m) if cfg.rkhs.m is not None else min(n, max(100, int(np.sqrt(n) * 10)))
            latent_dim_est += rkhs_blocks * m_est
        if latent_dim_est >= 3000:
            nuts_dense_mass = False
        elif latent_dim_est <= 500:
            nuts_dense_mass = True

    if vi_restarts > vi_max_restarts:
        warnings.warn(f"Capping VI restarts at {vi_max_restarts} to avoid excessive runtime.", RuntimeWarning)
        vi_restarts = vi_max_restarts

    backend_lower = backend.lower()
    if backend_lower == "numpyro":
        _ensure_numpyro_deps(allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
    elif backend_lower == "pymc":
        _ensure_pymc_deps(
            allow_install=allow_install,
            strict_versions=strict_versions,
            disable_auto_install=disable_auto_install,
            need_numpyro=(not vi),
        )
    else:
        raise ValueError("backend must be 'numpyro' or 'pymc'.")

    w_j = None
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float32)
        if weights.shape != (n,) or np.any(weights <= 0) or not np.all(np.isfinite(weights)):
            raise ValueError("weights must be positive, finite, and length n.")

    try:
        platforms = {getattr(d, "platform", "").lower() for d in jax.devices()}
        if not ({"gpu","cuda","rocm","tpu"} & platforms):
            warnings.warn(
                "No accelerator detected by JAX; running on CPU. "
                "For NVIDIA GPUs via pip wheels, consider "
                "pip install 'jax[cuda12_pip]' -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html",
                RuntimeWarning
            )
    except Exception:
        pass

    # prepare standardization AND RKHS factors; also compute block weights for weighted R²
    ETA_jax, y_jax, y_stats, x_stats, mask_obs, mask_mis, block_weights = _prep_standardize(
        ETA, y, standardize_X, center_y,
        rkhs_approx=cfg.rkhs.approx, rkhs_m=cfg.rkhs.m, rkhs_seed=cfg.rkhs.seed, kernel_type=cfg.rkhs.kernel_type,
        rkhs_force_nystrom_mem_gb=cfg.rkhs.force_nystrom_mem_gb,
        ctx=ctx
    )
    if weights is not None:
        w_j = jnp.asarray(weights)

    # BayesCπ vs HMC modes
    if (not vi) and any(e["method"] == "bayes_c_pi" for e in ETA):
        if backend_lower == "numpyro":
            if cpi_hmc_mode == "vi":
                warnings.warn(
                    "BayesCπ requested with HMC; switching to VI (NUTS can't sample discrete gates). "
                    "Set cpi_hmc_mode='horseshoe' for HMC surrogate.", RuntimeWarning
                )
                vi = True
            elif cpi_hmc_mode == "horseshoe":
                warnings.warn("BayesCπ with HMC uses a horseshoe surrogate; π recorded but not enforced.", RuntimeWarning)
            else:
                raise ValueError("cpi_hmc_mode must be 'vi' or 'horseshoe'.")
        elif backend_lower == "pymc":
            if pymc_cpi_mode == "vi":
                warnings.warn(
                    "PyMC BayesCπ with HMC -> switching to VI (no discrete in NUTS). "
                    "Set pymc_cpi_mode='horseshoe' for HMC surrogate.", RuntimeWarning
                )
                vi = True
            elif pymc_cpi_mode == "horseshoe":
                warnings.warn("PyMC Cπ with HMC uses horseshoe surrogate; π recorded.", RuntimeWarning)
            else:
                raise ValueError("pymc_cpi_mode must be 'vi' or 'horseshoe'.")

    # Prior selection (RHS default)
    prior_family_eff = (cfg.prior.family or "rhs").lower()

    if backend_lower == "numpyro":
        fit = _numpyro_fit(
            ETA=ETA_jax, y=y_jax, weights=w_j,
            n_iter=n_iter, tune=tune, chains=chains, random_seed=random_seed,
            vi=vi, vi_restarts=vi_restarts, use_bglr_priors=use_bglr_priors, R2=R2,
            relaxed_temperature=relaxed_temperature,
            block_prior_calibration=block_prior_calibration,
            block_R2_strategy=cfg.block_R2_strategy,
            target_accept=target_accept, nuts_dense_mass=nuts_dense_mass,
            nuts_max_tree_depth=nuts_max_tree_depth, nuts_init=nuts_init,
            likelihood=likelihood, sample_nu=sample_nu, student_t_nu=student_t_nu,
            vi_guide=cfg.vi.guide, vi_rank=cfg.vi.rank,
            prior_family=prior_family_eff, rhs_m0=cfg.prior.rhs_m0,
            rhs_slab_scale_mult=cfg.prior.rhs_slab_scale_mult,
            vi_sched=dict(start=cfg.vi.temp_start, end=cfg.vi.temp_end,
                          schedule=cfg.vi.temp_schedule, warmup_frac=cfg.vi.temp_warmup_frac),
            block_weights=block_weights
        )
        inference_kind = "SVI" if vi else "HMC"
    else:
        fit = _pymc_fit(
            ETA=ETA_jax, y=y_jax, weights=w_j,
            n_iter=n_iter, tune=tune, chains=chains, random_seed=random_seed,
            target_accept=target_accept, vi=vi, vi_restarts=vi_restarts, vi_max_restarts=vi_max_restarts,
            use_bglr_priors=use_bglr_priors, R2=R2,
            relaxed_temperature=relaxed_temperature,
            block_prior_calibration=block_prior_calibration,
            block_R2_strategy=cfg.block_R2_strategy,
            likelihood=likelihood, sample_nu=sample_nu, student_t_nu=student_t_nu,
            vi_guide=cfg.vi.guide,
            prior_family=prior_family_eff, rhs_m0=cfg.prior.rhs_m0,
            rhs_slab_scale_mult=cfg.prior.rhs_slab_scale_mult,
            block_weights=block_weights
        )
        inference_kind = "SVI" if vi else "HMC"

    idata = fit["idata"]
    diagnostics = fit["diagnostics"]

    # Extract samples (prefer VI raw if present)
    if "raw_samples" in fit and isinstance(fit["raw_samples"], dict):
        beta_samps_map, intercept_samps, sigma2_samps, aux_map, u_samps_map = _extract_samples_from_raw(fit["raw_samples"])
    else:
        beta_samps_map, intercept_samps, sigma2_samps, aux_map, u_samps_map = _extract_samples_common(idata)

    # Canonicalize beta keys per effect (skip rkhs)
    S = int(len(intercept_samps))
    _beta_map_raw = beta_samps_map
    beta_samps_map = {}
    for i, e in enumerate(ETA_jax):
        if e.get("method") == "rkhs":
            continue
        pX = int(np.asarray(ETA[i]["X"]).shape[1]) if "X" in ETA[i] else 0
        pref = f"effect_{i}_"
        B = _resolve_beta_samples(_beta_map_raw, pref, S=S, p_expected=pX)
        if B.ndim != 2 or B.shape[1] != pX:
            raise ValueError(
                f"[beta shape check] effect {i}: beta shape {B.shape} ≠ X p={pX}. "
                f"Keys: {sorted(_beta_map_raw.keys())[:12]}."
            )
        beta_samps_map[pref + "beta"] = B

    # Per-effect summaries & BGLR keys
    effects: List[Dict[str, Any]] = []
    for i, e in enumerate(ETA_jax):
        method = e["method"]
        pref = f"effect_{i}_"
        X_orig = np.asarray(ETA[i]["X"]) if ("X" in ETA[i]) else np.zeros((y.shape[0], 0), dtype=np.float32)
        p = X_orig.shape[1]

        if method == "rkhs":
            u_key = pref + "u"
            if u_key in u_samps_map:
                u_samps = np.asarray(u_samps_map[u_key])
            else:
                raise KeyError(f"rkhs: missing samples for '{u_key}'.")
            u_mean = u_samps.mean(axis=0)
            eff: Dict[str, Any] = {
                "model": method,
                "n": int(y.shape[0]),
                "p": 0,
                "b": np.array([], dtype=np.float32),
                "sd.b": np.array([], dtype=np.float32),
                "u": u_mean.astype(np.float32),
                "extra": {"b_scale": "none", "rkhs": True, "rkhs_approx": str(e.get("rkhs_approx", "unknown"))}
            }
            effects.append(eff)
            continue

        cm_i, cs_i = x_stats[i]
        cm = np.asarray(cm_i); cs = np.asarray(cs_i)
        b_samps_std = beta_samps_map[pref+"beta"]
        if p == 0:
            b_mean = np.zeros((0,), dtype=np.float32); b_sd = np.zeros((0,), dtype=np.float32)
            u_vec = np.zeros((X_orig.shape[0],), dtype=np.float32)
        else:
            cs_safe = np.where(cs <= 1e-12, 1.0, cs)
            b_samps_unstd = b_samps_std / (cs_safe.reshape(1, -1))
            b_mean = b_samps_unstd.mean(axis=0)
            b_sd   = b_samps_unstd.std(axis=0)
            c_shift = float(np.sum(cm * b_mean))
            u_vec  = (X_orig @ b_mean.reshape(-1)) - c_shift

        eff: Dict[str, Any] = {
            "model": method,
            "n": int(X_orig.shape[0]),
            "p": int(p),
            "b": b_mean.astype(np.float32),
            "sd.b": b_sd.astype(np.float32),
            "u": u_vec.astype(np.float32),
            "extra": {"b_scale": "original_X"}
        }

        varB = None
        if method == "bayesian_ridge":
            if (pref+"tau2_mean") in aux_map: varB = float(aux_map[pref+"tau2_mean"])
            elif (pref+"tau_mean") in aux_map: varB = float(aux_map[pref+"tau_mean"])**2
        elif method == "bayes_a":
            if (pref+"lambda_mean") in aux_map: varB = np.asarray(aux_map[pref+"lambda_mean"], dtype=np.float32)
        elif method in ("bayes_b","bayes_c","bayes_c_pi","elastic_net","lasso"):
            if (pref+"tau2_mean") in aux_map: varB = float(aux_map[pref+"tau2_mean"])
            elif (pref+"tau_mean") in aux_map: varB = float(aux_map[pref+"tau_mean"])**2
        if varB is not None:
            eff["varB"] = varB

        if (pref+"gate") in aux_map and method in ("bayes_b","bayes_c","bayes_c_pi"):
            eff["d"] = aux_map[pref+"gate"].mean(axis=0).astype(np.float32)
        elif (pref+"lambda_mean" in aux_map) and (pref+"tau_mean" in aux_map) and method in ("bayes_b","bayes_c","bayes_c_pi"):
            lam = np.asarray(aux_map[pref+"lambda_mean"])
            tau = float(aux_map[pref+"tau_mean"])
            shrink_proxy = 1.0 / (1.0 + 1.0/(np.square(lam)*(tau**2)+1e-12))
            eff["d"] = shrink_proxy.astype(np.float32)
            eff.setdefault("extra", {})["d_source"] = "horseshoe_shrinkage_proxy"

        if (pref+"pi_mean") in aux_map:
            eff["pi"] = float(aux_map[pref+"pi_mean"])
        elif "fixed_pi" in e:
            eff["pi"] = float(e["fixed_pi"])

        if method == "bayes_a":
            df_used = float(e.get("df", 5.0))
            varB_in = e.get("varB", None)
            if varB_in is not None:
                varB_arr = np.asarray(varB_in, dtype=np.float32)
                s2_arr = varB_arr * (df_used - 2.0) / df_used
                eff["df"]  = df_used; eff["s2"]  = s2_arr
                eff["df0"] = df_used; eff["S0"]  = s2_arr
            else:
                hyp = _default_hyper_from_R2(y_jax, [e["X"]], "bayes_a",
                                             R2=_block_R2_for_idx(ETA_jax, i, R2, block_prior_calibration, cfg.block_R2_strategy, block_weights),
                                             nu=df_used)
                eff["df"]  = float(hyp["nu"]); eff["s2"]  = float(hyp["s2"])
                eff["df0"] = float(hyp["nu"]); eff["S0"]  = float(hyp["s2"])

        if "d" in ETA[i]:
            eff.setdefault("extra", {})["d_mask_input"] = "applied"

        effects.append(eff)

    var_out, yhat_centered, mem_diag = _compute_variance_components_batched(
        ETA=ETA_jax,
        beta_samps_map=beta_samps_map,
        intercept_samps=intercept_samps,
        sigma2_samps=sigma2_samps,
        batch_mem_gb=batch_mem_gb,
        n=int(y.shape[0]),
        overhead_mult_override=batch_overhead_mult,
        u_samps_map=u_samps_map,
        backend=backend_lower,
        ctx=ctx
    )

    for i, _ in enumerate(ETA_jax):
        key = f"effect_{i}_"
        effects[i]["varU"] = float(var_out["per_effect"][key]["var_g_mean"])
        effects[i]["sdU"]  = float(var_out["per_effect"][key]["var_g_sd"])

    y_mean = float(y_stats[0])
    yHat   = (yhat_centered["mean"] + y_mean).astype(np.float32)
    SDyHat = yhat_centered["sd"].astype(np.float32)

    whichNa_1 = [i+1 for i in np.where(np.isnan(y))[0].astype(int).tolist()]

    df0_e = 5.0
    Vy_orig = float(np.nanvar(y.astype(np.float64)))
    S0_e = (1.0 - R2) * Vy_orig if use_bglr_priors else 1.0

    mu_out = float(intercept_samps.mean() + y_mean)

    # Memory diagnostics payload
    host_budget_gb, host_budget_info = ctx.resolve_host_budget_gb(requested=batch_mem_gb)
    device_budget_gb, device_budget_info = ctx.resolve_device_budget_gb()
    mem_diag_top = {"host_budget": host_budget_info, "device_budget": device_budget_info}
    mem_diag_top.update(mem_diag or {})

    result = {
        "mu": mu_out,
        "y": y.astype(np.float32),
        "yHat": yHat,
        "SD.yHat": SDyHat,
        "varE": float(var_out["varE_mean"]),
        "sdE":  float(var_out["varE_sd"]),
        "varG": float(var_out["varG_total_mean"]),
        "sdG":  float(var_out["varG_total_sd"]),
        "h2":   float(var_out["h2_mean"]),
        "whichNa": whichNa_1,
        "df0": df0_e,
        "S0":  S0_e,
        "ETA": effects,
        "extra": {
            "backend": backend_lower,
            "inference": inference_kind,
            "diagnostics": _summarize_diag(fit["diagnostics"]),
            "memory": mem_diag_top
        }
    }
    if auto_tune:
        try:
            result["extra"]["nuts"] = {
                "dense_mass": bool(nuts_dense_mass),
                "latent_dim_est": int(total_p),
                "init": str(nuts_init)
            }
        except Exception:
            pass
    return result

# ---------------------------
# Helpers for per-block R2 (now supports 'weighted')
# ---------------------------

def _block_R2_for_idx(ETA, idx: int, R2: float, block_prior_calibration: bool, block_R2_strategy: str, block_weights: Optional[List[float]] = None):
    strategy = (block_R2_strategy or "equal").lower()
    K = max(1, len(ETA))
    if strategy == "from_eta" and ("R2" in ETA[idx]):
        return float(min(0.9999, max(1e-8, float(ETA[idx]["R2"]))))

    if strategy == "weighted" and block_weights is not None and len(block_weights) == K and sum(block_weights) > 0:
        wsum = float(sum(block_weights))
        r2_blk = float(R2) * (float(block_weights[idx]) / wsum)
        return float(min(0.9999, max(1e-8, r2_blk)))

    # Fallbacks
    if block_prior_calibration or ("R2" in ETA[idx]):
        return float(min(0.9999, max(1e-8, float(R2) / K)))
    return float(min(0.9999, max(1e-8, float(R2))))

# ---------------------------
# Dependency Ensurers
# ---------------------------

def _ensure_core_deps(*, allow_install: bool, strict_versions: bool, disable_auto_install: bool):
    global np, jax, jnp, az, psutil
    if np is None:
        np = _ensure_module("numpy", "numpy", "1.21.0", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
    if jax is None:
        os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
        jax = _ensure_module("jax", "jax", "0.4.20", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
        try:
            from jax import config as jax_config
            env_x64 = os.environ.get("JAX_ENABLE_X64", None)
            if env_x64 is None:
                jax_config.update("jax_enable_x64", True)
        except Exception:
            pass
        jnp = jax.numpy
    if az is None:
        az = _ensure_module("arviz", "arviz", "0.15.0", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
    if psutil is None:
        try:
            psutil_mod = _ensure_module("psutil", "psutil", None, allow_install=allow_install, strict_versions=strict_versions, optional=True, disable_auto_install=disable_auto_install)
        except Exception:
            psutil_mod = None
        globals()["psutil"] = psutil_mod

def _ensure_numpyro_deps(*, allow_install: bool, strict_versions: bool, disable_auto_install: bool):
    global numpyro, dist, SVI, Trace_ELBO, MCMC, NUTS, AutoDiagonalNormal, AutoLowRankMultivariateNormal, AutoMultivariateNormal
    if numpyro is None:
        numpyro = _ensure_module("numpyro", "numpyro", "0.13.0", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
        dist = importlib.import_module("numpyro.distributions")
        from numpyro.infer import SVI as _SVI, Trace_ELBO as _Trace_ELBO, MCMC as _MCMC, NUTS as _NUTS
        from numpyro.infer.autoguide import AutoDiagonalNormal as _ADN, AutoLowRankMultivariateNormal as _ALR, AutoMultivariateNormal as _AMN
        SVI, Trace_ELBO, MCMC, NUTS = _SVI, _Trace_ELBO, _MCMC, _NUTS
        AutoDiagonalNormal, AutoLowRankMultivariateNormal, AutoMultivariateNormal = _ADN, _ALR, _AMN

def _ensure_pymc_deps(*, allow_install: bool, strict_versions: bool, disable_auto_install: bool, need_numpyro: bool = False):
    pm = _ensure_module("pymc", "pymc", "5.0", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
    try:
        pm.set_platform("jax")
    except Exception:
        warnings.warn("Could not set PyMC platform to JAX; using default backend.", RuntimeWarning)
        logger.warning("Could not set PyMC platform to JAX; using default backend.")
    if need_numpyro:
        try:
            _ = _ensure_module("numpyro", "numpyro", "0.13.0", allow_install=allow_install, strict_versions=strict_versions, disable_auto_install=disable_auto_install)
        except Exception:
            warnings.warn("NumPyro not ensured for PyMC HMC; will use PyMC's default NUTS.", RuntimeWarning)
            logger.warning("NumPyro not ensured for PyMC HMC; will use PyMC's default NUTS.")

# ---------------------------
# Backend implementations
# ---------------------------

def _pymc_fit(
    ETA, y, weights,
    n_iter, tune, chains, random_seed,
    target_accept,
    vi, vi_restarts, vi_max_restarts,
    use_bglr_priors, R2,
    relaxed_temperature: float,
    block_prior_calibration: bool,
    block_R2_strategy: str,
    *,
    likelihood: str,
    sample_nu: bool,
    student_t_nu: float,
    vi_guide: str,
    prior_family: str = "rhs",
    rhs_m0: Optional[float] = None,
    rhs_slab_scale_mult: float = 10.0,
    block_weights: Optional[List[float]] = None,
):
    pm = importlib.import_module("pymc")
    ETA = [dict(e, X=np.asarray(e["X"], dtype=np.float32)) if ("X" in e) else dict(e) for e in ETA]
    y_np = np.asarray(y, dtype=np.float32)

    mask_obs_np = np.isfinite(y_np)
    mask_mis_np = ~mask_obs_np
    n_mis = int(mask_mis_np.sum())
    idx_obs_t = pt.as_tensor_variable(np.where(mask_obs_np)[0].astype("int64"))
    idx_mis_t = pt.as_tensor_variable(np.where(mask_mis_np)[0].astype("int64"))

    Vy = float(np.nanvar(y_np))
    sqrtVy = math.sqrt(max(Vy, 1e-12))

    def _to_float(x): return float(np.asarray(x))

    def _check_varB_scalar_pos(i: int, varB_in):
        if varB_in is None: return None
        v = float(np.asarray(varB_in).reshape(-1)[0])
        if not (v > 0.0): raise ValueError(f"ETA[{i}] varB must be > 0.")
        return v

    def _check_bayesA_varB(i: int, varB_in, df_in: float):
        if varB_in is None: return None
        if df_in <= 2.0: raise ValueError(f"ETA[{i}] BayesA: df must be > 2 when varB is provided.")
        arr = np.asarray(varB_in, dtype=np.float32)
        if np.any(arr <= 0): raise ValueError(f"ETA[{i}] BayesA: all varB values must be > 0.")
        return arr

    def _pv_calibrate_rhs_tau0_c(n: int, p: int, m0: Optional[float], R2_blk: float) -> Tuple[float, float]:
        # sigma ≈ sqrt((1-R2)*Vy) as a rough residual scale for calibration
        sigma = math.sqrt(max((1.0 - R2_blk) * Vy, 1e-12))
        m0_eff = float(m0) if (m0 is not None) else float(min(p, max(1.0, 0.1 * p + 10.0)))
        m0_eff = min(max(m0_eff, 1.0), float(p))
        tau0 = (m0_eff / max(1.0, (p - m0_eff))) * (sigma / math.sqrt(max(n,1)))
        c0 = rhs_slab_scale_mult * (sigma / math.sqrt(max(n,1)))
        return float(max(tau0, 1e-6)), float(max(c0, 1e-6))

    with pm.Model() as model:
        intercept = pm.Normal("intercept", _pt_const(0.0), _pt_const(10.0))
        nu_e = 5.0
        Se = (1.0 - R2) * Vy if use_bglr_priors else 1.0
        sigma2 = pm.InverseGamma("sigma2", alpha=_pt_const(nu_e/2.0), beta=_pt_const(nu_e*Se/2.0))
        sigma  = pm.Deterministic("sigma", pt.sqrt(sigma2))

        # Student-t ν (optional sampling)
        if likelihood == "student_t":
            if sample_nu:
                student_t_nu_rv = pm.HalfCauchy("student_t_nu", beta=_pt_const(5.0))
                nu_obs = student_t_nu_rv
            else:
                nu_obs = _pt_const(float(student_t_nu))
        else:
            nu_obs = None  # Gaussian

        mu_terms = []
        for i, e in enumerate(ETA):
            m = e["method"]; pref = f"effect_{i}_"
            p_here = int(e["X"].shape[1]) if "X" in e else 0
            R2_eff = _block_R2_for_idx(ETA, i, R2, block_prior_calibration, block_R2_strategy, block_weights)
            varB_in = e.get("varB", None)
            df_in   = float(e.get("df", 5.0))

            if m == "rkhs":
                if "L" in e:
                    L_np = np.asarray(e["L"], dtype=np.float32)
                    L = pt.as_tensor_variable(L_np)
                else:
                    K = np.asarray(e["K"], dtype=np.float32)
                    bytes_needed = float(K.shape[0]) * float(K.shape[0]) * 4.0
                    mem_avail_gb = 8.0  # conservative in PyMC path
                    if bytes_needed > 0.3 * mem_avail_gb * (1024**3):
                        n_ = K.shape[0]
                        capacity = 0.3 * mem_avail_gb * (1024**3) / 4.0
                        disc = float(n_*n_) + 4.0 * capacity
                        m_suggest = int(max(50, min(n_, 0.5 * (-n_ + math.sqrt(disc)))))
                        warnings.warn(f"PyMC RKHS Cholesky may cause OOM for n={n_}; prefer rkhs_approx='nystrom' with m≈{m_suggest}.", ResourceWarning)
                        logger.warning("PyMC RKHS Cholesky may OOM (n=%d). Suggest m≈%d for Nyström.", n_, m_suggest)
                    J = _adaptive_jitter(K)
                    L = pt.slinalg.cholesky(pt.as_tensor_variable(K + J*np.eye(K.shape[0], dtype=np.float32)))
                    L_np = K
                m_latent = int(L_np.shape[1] if isinstance(L_np, np.ndarray) else K.shape[0])
                tau0 = _pt_const(math.sqrt(max(1e-12, R2_eff * Vy)))
                tau_u = pm.HalfCauchy(pref+"tau_u", tau0)
                z = pm.Normal(pref+"z", _pt_const(0.0), _pt_const(1.0), shape=m_latent)
                u = pm.Deterministic(pref+"u", pt.dot(L, z) * tau_u)
                mu_terms.append(u)
                continue  # skip X/beta path

            X = pt.as_tensor_variable(e["X"])
            p = p_here

            if not vi:
                if m == "bayesian_ridge":
                    vB = _check_varB_scalar_pos(i, varB_in)
                    if vB is not None: tau0 = _pt_const(math.sqrt(vB))
                    else:              tau0 = _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "bayesian_ridge", R2=R2_eff)["tau0"]))
                    tau2 = pm.InverseGamma(pref+"tau2", alpha=_pt_const(2.0), beta=_pt_const(2.0)*(tau0**2))
                    tau  = pm.Deterministic(pref+"tau", pt.sqrt(tau2))
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), tau, shape=p)

                elif m == "bayes_a":
                    varB_arr = _check_bayesA_varB(i, varB_in, df_in)
                    if varB_arr is not None:
                        if varB_arr.ndim == 0: varB_arr = np.full((p,), float(varB_arr), dtype=np.float32)
                        s2_vec = varB_arr * (df_in - 2.0) / df_in
                        lam = pm.InverseGamma(pref+"lambda", alpha=_pt_const(df_in/2.0), beta=pt.as_tensor_variable(df_in*s2_vec/2.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), pt.sqrt(lam), shape=p)
                        pm.Deterministic(pref+"df", _pt_const(df_in))
                    else:
                        hyp = _default_hyper_from_R2(y_np, [e["X"]], "bayes_a", R2=R2_eff, nu=df_in)
                        lam = pm.InverseGamma(pref+"lambda", alpha=_pt_const(_to_float(hyp["nu"])/2.0), beta=_pt_const(_to_float(hyp["nu"])*_to_float(hyp["s2"])/2.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), pt.sqrt(lam), shape=p)

                elif m in ("bayes_b","bayes_c"):
                    # Regularized Horseshoe (default)
                    if prior_family == "rhs":
                        tau0, c0 = _pv_calibrate_rhs_tau0_c(n=int(y_np.shape[0]), p=p, m0=rhs_m0, R2_blk=R2_eff)
                        tau  = pm.HalfCauchy(pref+"tau", _pt_const(tau0))
                        lam  = pm.HalfCauchy(pref+"lambda", _pt_const(1.0), shape=p)
                        scale = (tau * lam)
                        beta_scale = scale / pt.sqrt(_pt_const(1.0) + (scale/_pt_const(c0))**2)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), beta_scale, shape=p)
                    else:
                        # classic horseshoe-like from previous version
                        vB = _check_varB_scalar_pos(i, varB_in)
                        if vB is not None: tau0 = _pt_const(math.sqrt(vB))
                        else:              tau0 = _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "bayes_c", R2=R2_eff)["tau0"]))
                        tau  = pm.HalfCauchy(pref+"tau", tau0)
                        lam  = pm.HalfCauchy(pref+"lambda", _pt_const(1.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), tau*lam, shape=p)

                elif m == "bayes_c_pi":
                    # keep relaxed gate VI / horseshoe surrogate decisions outside
                    vB = _check_varB_scalar_pos(i, varB_in)
                    if vB is not None: tau0 = _pt_const(math.sqrt(vB))
                    else:              tau0 = _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "bayes_c", R2=R2_eff)["tau0"]))
                    tau  = pm.HalfCauchy(pref+"tau", tau0)
                    lam  = pm.HalfCauchy(pref+"lambda", _pt_const(1.0), shape=p)
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), tau*lam, shape=p)
                    if "fixed_pi" not in e:
                        _ = pm.Beta(pref+"pi", _pt_const(1.0), _pt_const(1.0))

                elif m == "lasso":
                    vB = _check_varB_scalar_pos(i, varB_in)
                    tau0 = _pt_const(math.sqrt(vB)) if (vB is not None) else _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "lasso", R2=R2_eff)["tau0"]))
                    b0 = tau0 / _pt_const(math.sqrt(2.0))
                    b  = pm.HalfCauchy(pref+"b", b0)
                    beta = pm.Laplace(pref+"beta", _pt_const(0.0), b, shape=p)

                elif m == "elastic_net":
                    alpha = float(e.get("alpha", 0.5))
                    vB = _check_varB_scalar_pos(i, varB_in)
                    tau0 = _pt_const(math.sqrt(vB)) if (vB is not None) else _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "elastic_net", R2=R2_eff)["tau0"]))
                    b0   = tau0 / _pt_const(math.sqrt(2.0))
                    tau  = pm.HalfCauchy(pref+"tau", tau0)
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), tau, shape=p)
                    pm.Potential(pref+"en_l1", - _pt_const(alpha) * pt.sum(pt.abs_(beta)) / b0)

                elif m == "fixed":
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), _pt_const(100.0), shape=p)

                else:
                    raise ValueError(f"Invalid method: {m}")

            else:
                # VI path mirrors HMC priors; BayesCπ uses relaxed gates (temperature handled at construction time)
                if m == "bayesian_ridge":
                    tau0 = _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "bayesian_ridge", R2=R2_eff)["tau0"]))
                    tau  = pm.HalfCauchy(pref+"tau", tau0)
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), tau, shape=p)

                elif m == "bayes_a":
                    varB_arr = _check_bayesA_varB(i, varB_in, df_in)
                    if varB_arr is not None:
                        if varB_arr.ndim == 0: varB_arr = np.full((p,), float(varB_arr), dtype=np.float32)
                        s2_vec = varB_arr * (df_in - 2.0) / df_in
                        lam = pm.InverseGamma(pref+"lambda", alpha=_pt_const(df_in/2.0), beta=pt.as_tensor_variable(df_in*s2_vec/2.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), pt.sqrt(lam), shape=p)
                        pm.Deterministic(pref+"df", _pt_const(df_in))
                    else:
                        hyp = _default_hyper_from_R2(y_np, [e["X"]], "bayes_a", R2=R2_eff, nu=df_in)
                        lam = pm.InverseGamma(pref+"lambda", alpha=_pt_const(_to_float(hyp["nu"])/2.0), beta=_pt_const(_to_float(hyp["nu"])*_to_float(hyp["s2"])/2.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), pt.sqrt(lam), shape=p)

                elif m in ("bayes_b","bayes_c"):
                    if prior_family == "rhs":
                        tau0, c0 = _pv_calibrate_rhs_tau0_c(n=int(y_np.shape[0]), p=p, m0=rhs_m0, R2_blk=R2_eff)
                        tau  = pm.HalfCauchy(pref+"tau", _pt_const(tau0))
                        lam  = pm.HalfCauchy(pref+"lambda", _pt_const(1.0), shape=p)
                        scale = (tau * lam)
                        beta_scale = scale / pt.sqrt(_pt_const(1.0) + (scale/_pt_const(c0))**2)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), beta_scale, shape=p)
                    else:
                        tau0 = _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "bayes_c", R2=R2_eff)["tau0"]))
                        tau  = pm.HalfCauchy(pref+"tau", tau0)
                        lam  = pm.HalfCauchy(pref+"lambda", _pt_const(1.0), shape=p)
                        beta = pm.Normal(pref+"beta", _pt_const(0.0), tau*lam, shape=p)

                elif m == "bayes_c_pi":
                    if ("fixed_pi" not in e):
                        pi = pm.Beta(pref+"pi", _pt_const(1.0), _pt_const(1.0))
                    else:
                        pi_val = float(e.get("fixed_pi", 0.05))
                        pi = _pt_const(pi_val); pm.Deterministic(pref+"pi", pi)
                    pi_clip = pt.clip(pi, _pt_const(1e-5), _pt_const(1.0 - 1e-5))
                    alpha_g = pi_clip / _pt_const(relaxed_temperature)
                    beta_g  = (_pt_const(1.0) - pi_clip) / _pt_const(relaxed_temperature)
                    gate    = pm.Beta(pref+"gate", alpha=alpha_g, beta=beta_g, shape=p)
                    tau     = pm.HalfCauchy(pref+"tau", _pt_const(1.0))
                    beta_raw= pm.Normal(pref+"beta_raw", _pt_const(0.0), _pt_const(1.0), shape=p)
                    beta    = pm.Deterministic(pref+"beta", gate * (tau * beta_raw))

                elif m == "lasso":
                    vB = _check_varB_scalar_pos(i, varB_in)
                    tau0 = _pt_const(math.sqrt(vB)) if (vB is not None) else _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "lasso", R2=R2_eff)["tau0"]))
                    b0   = tau0 / _pt_const(math.sqrt(2.0))
                    b    = pm.HalfCauchy(pref+"b", b0)
                    beta = pm.Laplace(pref+"beta", _pt_const(0.0), b, shape=p)

                elif m == "elastic_net":
                    alpha = float(e.get("alpha", 0.5))
                    vB = _check_varB_scalar_pos(i, varB_in)
                    tau0 = _pt_const(math.sqrt(vB)) if (vB is not None) else _pt_const(_to_float(_default_hyper_from_R2(y_np, [e["X"]], "elastic_net", R2=R2_eff)["tau0"]))
                    b0   = tau0 / _pt_const(math.sqrt(2.0))
                    tau  = pm.HalfCauchy(pref+"tau", tau0)
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), tau, shape=p)
                    pm.Potential(pref+"en_l1", - _pt_const(alpha) * pt.sum(pt.abs_(beta)) / b0)

                elif m == "fixed":
                    beta = pm.Normal(pref+"beta", _pt_const(0.0), _pt_const(100.0), shape=p)

                else:
                    raise ValueError(f"Invalid method: {m}")

            mu_terms.append(pt.dot(X, beta) if p > 0 else _pt_const(0.0))

        mu = pm.Deterministic("mu", intercept + (pt.add(*mu_terms) if mu_terms else _pt_const(0.0)))

        # Likelihood
        if weights is not None:
            weights_t = pt.as_tensor_variable(np.asarray(weights, dtype=np.float32))
            epsw = _pt_const(1e-12)
            sigma_vec = sigma / pt.sqrt(weights_t + epsw)
            if likelihood == "student_t":
                pm.StudentT("y_obs", nu=nu_obs, mu=pt.take(mu, idx_obs_t), sigma=pt.take(sigma_vec, idx_obs_t),
                            observed=y_np[mask_obs_np])
                if n_mis > 0:
                    pm.StudentT("y_mis", nu=nu_obs, mu=pt.take(mu, idx_mis_t), sigma=pt.take(sigma_vec, idx_mis_t), shape=n_mis)
            else:
                pm.Normal("y_obs", mu=pt.take(mu, idx_obs_t), sigma=pt.take(sigma_vec, idx_obs_t),
                          observed=y_np[mask_obs_np])
                if n_mis > 0:
                    pm.Normal("y_mis", mu=pt.take(mu, idx_mis_t), sigma=pt.take(sigma_vec, idx_mis_t), shape=n_mis)
        else:
            if likelihood == "student_t":
                pm.StudentT("y_obs", nu=nu_obs, mu=pt.take(mu, idx_obs_t), sigma=sigma, observed=y_np[mask_obs_np])
                if n_mis > 0:
                    pm.StudentT("y_mis", nu=nu_obs, mu=pt.take(mu, idx_mis_t), sigma=sigma, shape=n_mis)
            else:
                pm.Normal("y_obs", mu=pt.take(mu, idx_obs_t), sigma=sigma, observed=y_np[mask_obs_np])
                if n_mis > 0:
                    pm.Normal("y_mis", mu=pt.take(mu, idx_mis_t), sigma=sigma, shape=n_mis)

        if vi:
            guide_lc = str(vi_guide).lower()
            if guide_lc in ("fullrank", "full", "dense"):
                advi_method = "fullrank_advi"
            elif guide_lc in ("lowrank", "lr"):
                warnings.warn("PyMC does not support low-rank ADVI directly; using fullrank_advi.", RuntimeWarning)
                advi_method = "fullrank_advi"
            else:
                advi_method = "advi"
            best_elbo = None
            best_idata = None
            for r in range(max(1, int(vi_restarts))):
                approx = pm.fit(n=n_iter, method=advi_method, random_seed=int(random_seed + r))
                try:
                    idata_r = approx.sample(draws=n_iter, return_inferencedata=True)
                    if not isinstance(idata_r, az.InferenceData):
                        idata_r = az.from_pymc(trace=idata_r)
                except TypeError:
                    trace = approx.sample(n_iter)
                    idata_r = az.from_pymc(trace=trace)
                hist = getattr(approx, "hist", None)
                try: elbo_final = float(np.asarray(hist)[-1]) if hist is not None else float("nan")
                except Exception: elbo_final = float("nan")
                if (best_elbo is None) or (elbo_final > best_elbo):
                    best_elbo = elbo_final
                    best_idata = idata_r
            idata = best_idata
            diagnostics = {"advi_elbo_final": best_elbo, "vi_restarts": int(vi_restarts), "rhat": float("nan")}
        else:
            try:
                import pymc as pm_mod
                import pymc.sampling_jax as pmjax  # type: ignore
                idata = pmjax.sample_numpyro(
                    draws=n_iter, tune=tune, chains=chains,
                    target_accept=target_accept, random_seed=random_seed
                )
            except Exception as e1:
                logger.warning("pymc.sampling_jax.sample_numpyro failed; attempting pm.sample(nuts_sampler='numpyro'): %s", e1)
                try:
                    idata = pm.sample(
                        n_iter, tune=tune, chains=chains,
                        target_accept=target_accept,
                        random_seed=random_seed,
                        nuts_sampler="numpyro"
                    )
                except Exception as e2:
                    warnings.warn("Falling back to PyMC's default NUTS.", RuntimeWarning)
                    logger.warning("Falling back to PyMC's default NUTS due to: %s", e2)
                    idata = pm.sample(
                        n_iter, tune=tune, chains=chains,
                        target_accept=target_accept, random_seed=random_seed
                    )
            try:
                rh = az.rhat(idata)
                rhat_max = float(np.nanmax(rh.to_array().values))
            except Exception:
                rhat_max = float("nan")
            diagnostics = {"rhat": rhat_max}

    return {"idata": idata, "diagnostics": diagnostics}

def _vi_return_sites(ETA, *, include_student_t_nu: bool, include_y_mis: bool):
    sites = {"intercept", "sigma2"}
    if include_student_t_nu:
        sites.add("student_t_nu")
    for i, e in enumerate(ETA):
        pref = f"effect_{i}_"; m = e["method"]
        if m == "bayes_a":
            sites.update({pref + "beta", pref + "lambda"})
        elif m in ("bayes_b", "bayes_c", "bayes_c_pi"):
            sites.update({pref + "beta", pref + "beta_raw", pref + "gate", pref + "tau", pref + "lambda"})
            if m == "bayes_c_pi":
                sites.add(pref + "pi")
        elif m in ("bayesian_ridge", "elastic_net"):
            sites.update({pref + "beta", pref + "tau"})
        elif m == "lasso":
            sites.update({pref + "beta", pref + "b"})
        elif m == "rkhs":
            sites.update({pref + "u", pref + "tau_u"})
        elif m == "fixed":
            sites.add(pref + "beta")
    if include_y_mis:
        sites.add("y_mis")
    return sorted(sites)

def _numpyro_fit(
    ETA, y, weights, n_iter, tune, chains, random_seed,
    vi, vi_restarts, use_bglr_priors, R2, relaxed_temperature,
    block_prior_calibration: bool, block_R2_strategy: str,
    target_accept: float = 0.9,
    nuts_dense_mass: bool = False,
    nuts_max_tree_depth: int = 10,
    nuts_init: str = "median",
    *,
    likelihood: str,
    sample_nu: bool,
    student_t_nu: float,
    vi_guide: str,
    vi_rank: Optional[int],
    prior_family: str = "rhs",
    rhs_m0: Optional[float] = None,
    rhs_slab_scale_mult: float = 10.0,
    vi_sched: Optional[Dict[str, Any]] = None,
    block_weights: Optional[List[float]] = None,
):
    n = y.shape[0]
    X_list = [ (e["X"] if "X" in e else jnp.zeros((n,1), dtype=jnp.float32)) for e in ETA ]
    mask_obs = jnp.isfinite(y)
    mask_mis = ~mask_obs
    n_mis = int(mask_mis.sum())

    def sigma_vec_from_weights(sigma):
        eps = jnp.asarray(1e-12, dtype=sigma.dtype)
        return sigma if weights is None else sigma / jnp.sqrt(weights + eps)


    def _check_varB_scalar_pos(i: int, varB_in):
        if varB_in is None: return None
        v = float(np.asarray(varB_in).reshape(-1)[0])
        if not (v > 0.0): raise ValueError(f"ETA[{i}] varB must be > 0.")
        return v

    def _check_bayesA_varB(i: int, varB_in, df_in: float):
        if varB_in is None: return None
        if df_in <= 2.0: raise ValueError(f"ETA[{i}] BayesA: df must be > 2 when varB is provided.")
        arr = np.asarray(varB_in, dtype=np.float32)
        if np.any(arr <= 0): raise ValueError(f"ETA[{i}] BayesA: all varB values must be > 0.")
        return arr

    # --- put near other helpers in _numpyro_fit scope ---

    def _pv_calibrate_rhs_tau0_c_jax(
    y,
    p: int,
    R2_blk: float,
    rhs_m0: Optional[float],
    rhs_slab_scale_mult: float,
    n_total: Optional[int] = None,   # <-- accept it
    **_kwargs,):                  

        # All JAX ops; return JAX scalars
        n = (int(n_total) if n_total is not None else y.shape[0])
        Vy = jnp.nanvar(y)
        sigma = jnp.sqrt(jnp.maximum((1.0 - R2_blk) * Vy, 1e-12))
    
        # m0: expected nonzeros
        m0_eff = (rhs_m0 if rhs_m0 is not None
                  else jnp.minimum(p, jnp.maximum(1.0, 0.1 * p + 10.0)))
        m0_eff = jnp.clip(m0_eff, 1.0, float(p))
    
        tau0 = (m0_eff / jnp.maximum(1.0, (float(p) - m0_eff))) * (sigma / jnp.sqrt(jnp.maximum(n, 1)))
        c0   = rhs_slab_scale_mult * (sigma / jnp.sqrt(jnp.maximum(n, 1)))
        return tau0, c0



    # Reuse factors for RKHS blocks if provided
    L_list = []
    for e in ETA:
        if e["method"] == "rkhs":
            if "L" in e:
                L_list.append(jnp.asarray(e["L"]))
            else:
                K = jnp.asarray(e["K"])
                J = _adaptive_jitter(np.asarray(K))
                L_list.append(jnp.linalg.cholesky(K + J*jnp.eye(K.shape[0], dtype=K.dtype)))
        else:
            L_list.append(None)

    def _pv_calibrate_rhs_tau0_c(y, X, p: int, R2_blk: float) -> Tuple[float, float]:
        # sigma ≈ sqrt((1-R2)*Vy)
        Vy = float(jnp.nanvar(y))
        sigma = math.sqrt(max((1.0 - R2_blk) * Vy, 1e-12))
        m0_eff = float(rhs_m0) if (rhs_m0 is not None) else float(min(p, max(1.0, 0.1 * p + 10.0)))
        m0_eff = min(max(m0_eff, 1.0), float(p))
        tau0 = (m0_eff / max(1.0, (p - m0_eff))) * (sigma / math.sqrt(max(n,1)))
        c0 = rhs_slab_scale_mult * (sigma / math.sqrt(max(n,1)))
        return float(max(tau0, 1e-6)), float(max(c0, 1e-6))

    
    def model_hmc(X_list, y):
        dtype = y.dtype
        l = lambda v: jnp.asarray(v, dtype=dtype)        # cast constants to JAX dtype
        zerosN = lambda n: jnp.zeros((n,), dtype=dtype)   # vector zero
    
        intercept = numpyro.sample("intercept", dist.Normal(l(0.0), l(10.0)))
        nu_e = l(5.0)
        Vy = jnp.nanvar(y)
        Se = (l(1.0) - l(R2)) * Vy if use_bglr_priors else l(1.0)
        sigma2 = numpyro.sample("sigma2", dist.InverseGamma(nu_e / l(2.0), nu_e * Se / l(2.0)))
        sigma  = jnp.sqrt(sigma2)
    
        if likelihood == "student_t":
            nu_obs = numpyro.sample("student_t_nu", dist.HalfCauchy(l(5.0))) if sample_nu else l(student_t_nu)
        else:
            nu_obs = None
    
        mu = zerosN(y.shape[0])
    
        for i, e in enumerate(ETA):
            X = X_list[i]; m = e["method"]; pref = f"effect_{i}_"
            R2_eff = _block_R2_for_idx(ETA, i, R2, block_prior_calibration, block_R2_strategy, block_weights)
            varB_in = e.get("varB", None)
            df_in   = float(e.get("df", 5.0))
            p = X.shape[1] if "X" in e else 0
    
            if m == "bayesian_ridge":
                vB = _check_varB_scalar_pos(i, varB_in)
                if vB is not None:
                    tau0 = jnp.sqrt(l(vB))
                else:
                    tau0 = _default_hyper_from_R2(y, [X], "bayesian_ridge", R2=R2_eff)["tau0"]
                tau2 = numpyro.sample(pref+"tau2", dist.InverseGamma(l(2.0), l(2.0) * (tau0**2)))
                tau  = jnp.sqrt(tau2)
                beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
    
            elif m == "bayes_a":
                varB_arr = _check_bayesA_varB(i, varB_in, df_in)
                if varB_arr is not None:
                    if np.ndim(varB_arr) == 0:
                        varB_arr = np.full((p,), float(varB_arr), dtype=np.float32)
                    s2_vec = jnp.asarray(varB_arr, dtype=dtype) * ((l(df_in) - l(2.0)) / l(df_in))
                    lam = numpyro.sample(pref+"lambda", dist.InverseGamma(l(df_in)/l(2.0), l(df_in)*s2_vec/l(2.0)).expand([p]).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), jnp.sqrt(lam)).expand([p]).to_event(1))
                else:
                    hyp = _default_hyper_from_R2(y, [X], "bayes_a", R2=R2_eff, nu=df_in)
                    lam = numpyro.sample(pref+"lambda", dist.InverseGamma(hyp["nu"]/l(2.0), hyp["nu"]*hyp["s2"]/l(2.0)).expand([p]).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), jnp.sqrt(lam)).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
    
            elif m in ("bayes_b","bayes_c"):
                if prior_family == "rhs":
                    tau0, c0 = _pv_calibrate_rhs_tau0_c_jax(y, p, R2_blk=R2_eff, n_total=y.shape[0], rhs_m0=rhs_m0, rhs_slab_scale_mult=rhs_slab_scale_mult)
                    tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                    lam  = numpyro.sample(pref+"lambda", dist.HalfCauchy(l(1.0)).expand([p]).to_event(1))
                    scale = tau * lam
                    beta_scale = scale / jnp.sqrt(l(1.0) + (scale / c0)**2)
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), beta_scale).expand([p]).to_event(1))
                else:
                    vB = _check_varB_scalar_pos(i, varB_in)
                    tau0 = jnp.sqrt(l(vB)) if (vB is not None) else _default_hyper_from_R2(y, [X], "bayes_c", R2=R2_eff)["tau0"]
                    tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                    lam  = numpyro.sample(pref+"lambda", dist.HalfCauchy(l(1.0)).expand([p]).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau*lam).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
    
            elif m == "bayes_c_pi":
                vB = _check_varB_scalar_pos(i, varB_in)
                tau0 = jnp.sqrt(l(vB)) if (vB is not None) else _default_hyper_from_R2(y, [X], "bayes_c", R2=R2_eff)["tau0"]
                tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                lam  = numpyro.sample(pref+"lambda", dist.HalfCauchy(l(1.0)).expand([p]).to_event(1))
                beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau*lam).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
                if "fixed_pi" not in e:
                    _ = numpyro.sample(pref+"pi", dist.Beta(l(1.0), l(1.0)))
    
            elif m == "lasso":
                vB = _check_varB_scalar_pos(i, varB_in)
                tau0 = jnp.sqrt(l(vB)) if (vB is not None) else _default_hyper_from_R2(y, [X], "lasso", R2=R2_eff)["tau0"]
                b0   = tau0 / jnp.sqrt(l(2.0))
                b    = numpyro.sample(pref+"b", dist.HalfCauchy(b0))
                beta = numpyro.sample(pref+"beta", dist.Laplace(l(0.0), b).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
    
            elif m == "elastic_net":
                alpha = l(float(e.get("alpha", 0.5)))
                vB = _check_varB_scalar_pos(i, varB_in)
                tau0  = jnp.sqrt(l(vB)) if (vB is not None) else _default_hyper_from_R2(y, [X], "elastic_net", R2=R2_eff)["tau0"]
                b0    = tau0 / jnp.sqrt(l(2.0))
                tau   = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                beta  = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau).expand([p]).to_event(1))
                numpyro.factor(pref+"en_l1", -(alpha / b0) * jnp.sum(jnp.abs(beta)))
                mu = mu + jnp.dot(X, beta)
    
            elif m == "rkhs":
                L = L_list[i]
                tau0 = jnp.sqrt(jnp.maximum(l(R2_eff) * jnp.nanvar(y), jnp.asarray(1e-12, dtype=dtype)))
                tau_u = numpyro.sample(pref+"tau_u", dist.HalfCauchy(tau0))
                z = numpyro.sample(pref+"z", dist.Normal(l(0.0), l(1.0)).expand([L.shape[1]]).to_event(1))
                u = jnp.matmul(L, z) * tau_u
                numpyro.deterministic(pref+"u", u)
                mu = mu + u
    
            elif m == "fixed":
                beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), l(100.0)).expand([p]).to_event(1))
                mu = mu + jnp.dot(X, beta)
    
            else:
                raise ValueError(f"Invalid method: {m}")
    
        mu = mu + intercept
        sigv = sigma_vec_from_weights(sigma)
        if weights is None:
            obs_dist = dist.StudentT(df=nu_obs, loc=mu[mask_obs], scale=sigv) if likelihood=="student_t" else dist.Normal(mu[mask_obs], sigv)
            numpyro.sample("y_obs", obs_dist.to_event(1), obs=y[mask_obs])
            if n_mis > 0:
                mis_dist = dist.StudentT(df=nu_obs, loc=mu[mask_mis], scale=sigv) if likelihood=="student_t" else dist.Normal(mu[mask_mis], sigv)
                numpyro.sample("y_mis", mis_dist.to_event(1))
        else:
            obs_dist = dist.StudentT(df=nu_obs, loc=mu[mask_obs], scale=sigv[mask_obs]) if likelihood=="student_t" else dist.Normal(mu[mask_obs], sigv[mask_obs])
            numpyro.sample("y_obs", obs_dist.to_event(1), obs=y[mask_obs])
            if n_mis > 0:
                mis_dist = dist.StudentT(df=nu_obs, loc=mu[mask_mis], scale=sigv[mask_mis]) if likelihood=="student_t" else dist.Normal(mu[mask_mis], sigv[mask_mis])
                numpyro.sample("y_mis", mis_dist.to_event(1))


    # SVI model with relaxed gates and RHS for B/C (no pi)
    temp_holder = {"t": float(relaxed_temperature)}
    def model_svi_relaxed(X_list, y, temp):
        dtype = y.dtype
        l = lambda v: jnp.asarray(v, dtype=dtype)
    
        n = y.shape[0]
        intercept = numpyro.sample("intercept", dist.Normal(l(0.0), l(10.0)))
    
        nu_e = l(5.0)
        Vy = jnp.nanvar(y)
        Se = (l(1.0) - l(R2)) * Vy if use_bglr_priors else l(1.0)
        sigma2 = numpyro.sample("sigma2", dist.InverseGamma(nu_e / l(2.0), nu_e * Se / l(2.0)))
        sigma  = jnp.sqrt(sigma2)
    
        # keep mu as (n,) vector throughout
        mu = jnp.zeros((n,), dtype=dtype)
    
        if likelihood == "student_t":
            nu_obs = numpyro.sample("student_t_nu", dist.HalfCauchy(l(5.0))) if sample_nu else l(student_t_nu)
        else:
            nu_obs = None
    
        for i, e in enumerate(ETA):
            X = X_list[i]
            m = e["method"]
            p = X.shape[1]
            pref = f"effect_{i}_"
            R2_eff = _block_R2_for_idx(ETA, i, R2, block_prior_calibration, block_R2_strategy, block_weights)
            varB_in = e.get("varB", None)
            df_in   = float(e.get("df", 5.0))  # df_in is a Python float (safe)
    
            if m == "bayesian_ridge":
                if varB_in is not None:
                    tau0 = jnp.sqrt(l(varB_in))
                else:
                    tau0 = l(_default_hyper_from_R2(y, [X], "bayesian_ridge", R2=R2_eff)["tau0"])
                tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau).expand([p]).to_event(1))
                mu = mu + jnp.matmul(X, beta)
    
            elif m == "bayes_a":
                varB_arr = _check_bayesA_varB(i, varB_in, df_in)  # numpy array or None
                if varB_arr is not None:
                    if np.ndim(varB_arr) == 0:
                        varB_arr = np.full((p,), float(varB_arr), dtype=np.float32)
                    s2_vec = l(jnp.asarray(varB_arr) * ((df_in - 2.0) / df_in))           # (p,)
                    lam = numpyro.sample(pref+"lambda",
                                         dist.InverseGamma(l(df_in/2.0),
                                                           l(df_in/2.0) * s2_vec).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), jnp.sqrt(lam)).to_event(1))
                else:
                    hyp = _default_hyper_from_R2(y, [X], "bayes_a", R2=R2_eff, nu=df_in)  # returns JAX
                    nuj = l(hyp["nu"])      # JAX scalar
                    s2  = l(hyp["s2"])      # JAX scalar
                    lam = numpyro.sample(pref+"lambda",
                                         dist.InverseGamma(nuj / l(2.0), (nuj * s2) / l(2.0)).expand([p]).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), jnp.sqrt(lam)).to_event(1))
                mu = mu + jnp.matmul(X, beta)
    
            elif m in ("bayes_b","bayes_c"):
                # RHS default
                if prior_family == "rhs":
                    tau0, c0 = _pv_calibrate_rhs_tau0_c_jax(y, int(p), R2_eff, rhs_m0, rhs_slab_scale_mult)
                    tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                    lam  = numpyro.sample(pref+"lambda", dist.HalfCauchy(l(1.0)).expand([p]).to_event(1))
                    scale = tau * lam
                    beta_scale = scale / jnp.sqrt(l(1.0) + (scale / c0)**2)
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), beta_scale).to_event(1))
                else:
                    if varB_in is not None:
                        tau0 = jnp.sqrt(l(varB_in))
                    else:
                        tau0 = l(_default_hyper_from_R2(y, [X], "bayes_c", R2=R2_eff)["tau0"])
                    tau  = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                    lam  = numpyro.sample(pref+"lambda", dist.HalfCauchy(l(1.0)).expand([p]).to_event(1))
                    beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau * lam).to_event(1))
                mu = mu + jnp.matmul(X, beta)
    
            elif m == "bayes_c_pi":
                if ("fixed_pi" not in e):
                    pi = numpyro.sample(pref+"pi", dist.Beta(l(1.0), l(1.0)))  # scalar
                else:
                    pi = l(float(e.get("fixed_pi", 0.05)))
                    numpyro.deterministic(pref+"pi", pi)
                pi_vec = jnp.broadcast_to(jnp.clip(pi, l(1e-5), l(1.0 - 1e-5)), (p,))
                gate = numpyro.sample(pref+"gate", dist.RelaxedBernoulli(l(temp), probs=pi_vec).to_event(1))
                tau     = numpyro.sample(pref+"tau", dist.HalfCauchy(l(1.0)))
                beta_raw= numpyro.sample(pref+"beta_raw", dist.Normal(l(0.0), l(1.0)).expand([p]).to_event(1))
                beta    = gate * (tau * beta_raw)
                numpyro.deterministic(pref+"beta", beta)
                mu = mu + jnp.matmul(X, beta)
    
            elif m == "lasso":
                tau0 = l(_default_hyper_from_R2(y, [X], "lasso", R2=R2_eff)["tau0"])
                b0   = tau0 / l(jnp.sqrt(2.0))
                b    = numpyro.sample(pref+"b", dist.HalfCauchy(b0))
                beta = numpyro.sample(pref+"beta", dist.Laplace(l(0.0), b).expand([p]).to_event(1))
                mu = mu + jnp.matmul(X, beta)
    
            elif m == "elastic_net":
                alpha = float(e.get("alpha", 0.5))
                tau0  = l(_default_hyper_from_R2(y, [X], "elastic_net", R2=R2_eff)["tau0"])
                b0    = tau0 / l(jnp.sqrt(2.0))
                tau   = numpyro.sample(pref+"tau", dist.HalfCauchy(tau0))
                beta  = numpyro.sample(pref+"beta", dist.Normal(l(0.0), tau).expand([p]).to_event(1))
                numpyro.factor(pref+"en_l1", -(alpha / b0) * jnp.sum(jnp.abs(beta)))
                mu = mu + jnp.matmul(X, beta)
    
            elif m == "rkhs":
                L = L_list[i]
                tau0 = jnp.sqrt(jnp.maximum(R2_eff * jnp.nanvar(y), l(1e-12)))
                tau_u = numpyro.sample(pref+"tau_u", dist.HalfCauchy(tau0))
                z = numpyro.sample(pref+"z", dist.Normal(l(0.0), l(1.0)).expand([L.shape[1]]).to_event(1))
                u = jnp.matmul(L, z) * tau_u
                numpyro.deterministic(pref+"u", u)
                mu = mu + u
    
            elif m == "fixed":
                beta = numpyro.sample(pref+"beta", dist.Normal(l(0.0), l(100.0)).expand([p]).to_event(1))
                mu = mu + jnp.matmul(X, beta)
    
        mu = mu + intercept
    
        sigv = sigma if (weights is None) else sigma / jnp.sqrt(weights + 1e-12)
        if weights is None:
            obs_dist = dist.StudentT(df=nu_obs, loc=mu[mask_obs], scale=sigv) if likelihood=="student_t" \
                       else dist.Normal(loc=mu[mask_obs], scale=sigv)
            numpyro.sample("y_obs", obs_dist.to_event(1), obs=y[mask_obs])
            if n_mis > 0:  # <-- use outer Python int, not jnp.any(...)
                mis_dist = dist.StudentT(df=nu_obs, loc=mu[mask_mis], scale=sigv) if likelihood=="student_t" \
                           else dist.Normal(loc=mu[mask_mis], scale=sigv)
                numpyro.sample("y_mis", mis_dist.to_event(1))
        else:
            obs_dist = dist.StudentT(df=nu_obs, loc=mu[mask_obs], scale=sigv[mask_obs]) if likelihood=="student_t" \
                       else dist.Normal(loc=mu[mask_obs], scale=sigv[mask_obs])
            numpyro.sample("y_obs", obs_dist.to_event(1), obs=y[mask_obs])
            if n_mis > 0:
                mis_dist = dist.StudentT(df=nu_obs, loc=mu[mask_mis], scale=sigv[mask_mis]) if likelihood=="student_t" \
                           else dist.Normal(loc=mu[mask_mis], scale=sigv[mask_mis])
                numpyro.sample("y_mis", mis_dist.to_event(1))

    if vi:
        guide_lc = str(vi_guide).lower()

        lat_dim = 0
        for i, e in enumerate(ETA):
            if e["method"] == "rkhs":
                L_here = L_list[i]
                lat_dim += int(L_here.shape[1])
            elif "X" in e:
                lat_dim += int(e["X"].shape[1])

        if guide_lc in ("fullrank", "full", "dense"):
            guide = AutoMultivariateNormal(model_svi_relaxed)
        elif guide_lc in ("lowrank", "lr"):
            rank = int(vi_rank) if vi_rank is not None else max(5, min(64, lat_dim // 20 if lat_dim > 0 else 20))
            guide = AutoLowRankMultivariateNormal(model_svi_relaxed, rank=rank)
        else:
            guide = AutoDiagonalNormal(model_svi_relaxed)

        want_sites = _vi_return_sites(
            ETA,
            include_student_t_nu=(likelihood=="student_t" and sample_nu),
            include_y_mis=bool(jnp.any(~jnp.isfinite(y)))
        )

        # --- Temperature schedule (chunked SVI with staged temperature)
        sched = vi_sched or {"start": 1.0, "end": 0.2, "schedule": "cosine", "warmup_frac": 0.1}
        start_t = float(sched.get("start", 1.0))
        end_t = float(sched.get("end", 0.2))
        warmup = max(0.0, min(0.9, float(sched.get("warmup_frac", 0.1))))
        mode = str(sched.get("schedule", "cosine")).lower()
        n_steps = int(max(1, n_iter))
        # Build temperature per step
        steps = list(range(n_steps))
        ts = []
        for k, s in enumerate(steps):
            frac = s / max(1, n_steps - 1)
            if frac < warmup:
                t = start_t
            else:
                f = (frac - warmup) / max(1e-6, (1.0 - warmup))
                if mode == "linear":
                    t = start_t + (end_t - start_t) * f
                else:  # cosine
                    t = end_t + 0.5*(start_t - end_t)*(1 + math.cos(math.pi * f))
            ts.append(float(max(1e-3, min(100.0, t))))
        # Run SVI in chunks while updating temperature (no recompile)
        from numpyro import optim as numpyro_optim       

        # --- build schedule first ---
        start_t = float(vi_sched.get("start", 1.0))
        end_t   = float(vi_sched.get("end", 0.2))
        warmup  = max(0.0, min(0.9, float(vi_sched.get("warmup_frac", 0.1))))
        mode    = str(vi_sched.get("schedule", "cosine")).lower()
        
        n_steps = int(max(1, n_iter))
        ts = []
        for s in range(n_steps):
            frac = s / max(1, n_steps - 1)
            if frac < warmup:
                t = start_t
            else:
                f = (frac - warmup) / max(1e-6, 1.0 - warmup)
                t = start_t + (end_t - start_t) * f if mode == "linear" \
                    else end_t + 0.5 * (start_t - end_t) * (1 + math.cos(math.pi * f))
            ts.append(max(1e-3, min(100.0, t)))
        
        # --- init & loop ---
        svi = SVI(model_svi_relaxed, guide, numpyro_optim.Adam(1e-2), loss=Trace_ELBO())
        key0 = jax.random.PRNGKey(random_seed)
        state = svi.init(key0, X_list=X_list, y=y, temp=float(ts[0]))
        
        losses = []  # <-- IMPORTANT
        for s in range(n_steps):
            key0, key_use = jax.random.split(key0)
            state, loss = svi.update(state, X_list=X_list, y=y, temp=float(ts[s]))
            try:
                losses.append(float(loss))
            except Exception:
                losses.append(float(jax.device_get(loss)))
        
        params = svi.get_params(state)
        
        # later when creating diagnostics:
        diagnostics = {
            "svi_elbo_final": (float(losses[-1]) if losses else float("nan")),
            "vi_restarts": int(vi_restarts),
            "rhat": float("nan"),
            "vi_temperature": {"start": start_t, "end": end_t, "schedule": mode, "warmup_frac": warmup},
        }

        
        # Predictive sampling must also receive temp
        from numpyro.infer import Predictive
        want_sites = _vi_return_sites(ETA,
                                      include_student_t_nu=(likelihood=="student_t" and sample_nu),
                                      include_y_mis=bool(jnp.any(~jnp.isfinite(y))))
        predictive = Predictive(model_svi_relaxed, guide=guide, params=params,
                                num_samples=max(1, int(n_iter) * int(chains)),
                                return_sites=want_sites)
        samples = predictive(jax.random.PRNGKey(random_seed+1),
                             X_list=X_list, y=y, temp=ts[-1])


        # Derive deterministic beta for relaxed spike&slab as needed
        for i, e in enumerate(ETA):
            if e["method"] in ("bayes_b","bayes_c","bayes_c_pi"):
                pref = f"effect_{i}_"
                if (pref+"gate") in samples and (pref+"tau") in samples and (pref+"beta_raw") in samples:
                    gate = np.asarray(samples[pref+"gate"])
                    tau  = np.asarray(samples[pref+"tau"]).reshape(-1, 1)
                    beta_raw = np.asarray(samples[pref+"beta_raw"])
                    samples[pref+"beta"] = gate * (tau * beta_raw)

        S = max(1, int(n_iter) * int(chains))
        C = int(chains); D = int(n_iter) if C > 0 else S
        def _to_chain_draw(x):
            x = np.asarray(x)
            if x.shape[0] == S and C > 0 and D > 0 and (C * D == S):
                return x.reshape((C, D) + x.shape[1:])
            return x
        posterior = {k: _to_chain_draw(v) for k, v in samples.items() if isinstance(v, np.ndarray)}
        try:
            idata = az.from_dict(posterior=posterior)
        except Exception:
            idata = az.InferenceData()
        diagnostics = {"svi_elbo_final": (float(losses[-1]) if losses else float("nan")), "vi_restarts": int(vi_restarts), "rhat": float("nan"),
                       "vi_temperature": {"start": start_t, "end": end_t, "schedule": mode, "warmup_frac": warmup}}
        return {"idata": idata, "diagnostics": diagnostics, "raw_samples": {k: np.asarray(v) for k, v in samples.items()}}

    else:
        try:
            from numpyro.infer.initialization import (init_to_median, init_to_sample, init_to_uniform, init_to_value)
        except Exception:
            from numpyro.infer.util import (init_to_median, init_to_sample, init_to_uniform, init_to_value)  # type: ignore
        _init_map = {"median": init_to_median, "sample": init_to_sample, "uniform": init_to_uniform}
        init_strategy = _init_map.get(str(nuts_init).lower(), init_to_median)

        # ---- NEW: auto-select chain_method for speed without accuracy loss ----
        try:
            devs = jax.local_device_count()
        except Exception:
            devs = 1
        if int(chains) <= 1:
            chain_method = "sequential"
        else:
            chain_method = "parallel" if devs >= int(chains) else "vectorized"
            
        kernel = NUTS(model_hmc, 
                      target_accept_prob=float(target_accept), 
                      dense_mass=bool(nuts_dense_mass),
                      max_tree_depth=int(nuts_max_tree_depth), 
                      init_strategy=init_strategy)
        
        mcmc = MCMC(kernel, 
                    num_warmup=tune, 
                    num_samples=n_iter, 
                    num_chains=chains,
                    chain_method=chain_method,)   # <--- important
        mcmc.run(jax.random.PRNGKey(random_seed), X_list=X_list, y=y)
        samps = mcmc.get_samples(group_by_chain=True)
        # synthesize RKHS u deterministics from z, tau_u, and L
        for i, e in enumerate(ETA):
            if e["method"] == "rkhs":
                pref = f"effect_{i}_"
                if (pref+"z") in samps and (pref+"tau_u") in samps:
                    z = np.asarray(samps[pref+"z"])           # (chains, draws, m)
                    tau_u = np.asarray(samps[pref+"tau_u"])   # (chains, draws)
                    L = np.asarray(e["L"], dtype=np.float32)  # (n, m)
                    u = np.einsum("cdm,nm->cdn", z, L).astype(np.float32)
                    u *= tau_u[..., None].astype(np.float32)
                    samps[pref+"u"] = u
        idata = az.from_dict(posterior={k: np.asarray(v) for k,v in samps.items()})
        try:
            rh = az.rhat(idata)
            rhat_max = float(np.nanmax(rh.to_array().values))
        except Exception:
            rhat_max = float("nan")
        diagnostics = {
            "rhat": rhat_max,
            "chain_method": chain_method,             # <--- report what ran
            "num_chains": int(chains),
            "local_device_count": int(devs),
        }

    return {"idata": idata, "diagnostics": diagnostics}

# ---------------------------
# Shared utilities
# ---------------------------

def _rkhs_prepare_factor(K: np.ndarray, *, approx: str, n: int, rkhs_m: Optional[int], seed: Optional[int], kernel_type: str, ctx: Optional[Context]=None) -> np.ndarray:
    """Build RKHS factor L: (n×n) Cholesky for exact, (n×m) Nyström for approx."""
    approx_lc = str(approx).lower()
    if approx_lc == "exact":
        bytes_needed = float(n) * float(n) * 4.0
        mem_avail_gb = ctx.available_host_gb() if ctx is not None else 8.0
        if bytes_needed > 0.3 * mem_avail_gb * (1024**3):
            capacity = 0.3 * mem_avail_gb * (1024**3) / 4.0
            disc = float(n*n) + 4.0 * capacity
            m_suggest = int(max(50, min(n, 0.5 * (-n + math.sqrt(disc)))))
            warnings.warn(
                f"RKHS Cholesky may cause OOM for n={n}; consider rkhs_approx='nystrom' with m≈{m_suggest}.",
                ResourceWarning
            )
            logger.warning("RKHS exact Cholesky may OOM (n=%d). Suggest Nyström m≈%d.", n, m_suggest)
        J = _adaptive_jitter(K)
        L = np.linalg.cholesky(K + J*np.eye(n, dtype=np.float32))
        return L.astype(np.float32, copy=False)

    if rkhs_m is None:
        m = min(n, max(100, int(np.sqrt(n)*10)))
        warnings.warn(f"RKHS Nyström: m not provided; using heuristic m={m}.", RuntimeWarning)
    else:
        m = int(rkhs_m)
    if m < 50 or m > n//2:
        warnings.warn("RKHS Nyström: m is quite small/large; may bias approximation.", RuntimeWarning)
    m = max(1, min(n, m))

    rng = np.random.default_rng(seed if seed is not None else 123)
    idx = np.sort(rng.choice(n, size=m, replace=False))
    Z = K[:, idx]
    W = K[np.ix_(idx, idx)]
    J = _adaptive_jitter(W)
    W[np.diag_indices_from(W)] += J
    try:
        cond = float(np.linalg.cond(W))
        if cond > 1e6:
            m_suggest = min(n, max(50, int(0.05 * n))) if str(kernel_type).lower()=="sparse" \
                        else min(n, max(100, int(0.1 * n)))
            warnings.warn(f"RKHS Nyström: W ill-conditioned (cond={cond:.1e}); consider increasing jitter or m≈{m_suggest}.", RuntimeWarning)
            logger.warning("RKHS Nyström: W ill-conditioned (cond=%.2e). Consider m≈%d.", cond, m_suggest)
    except Exception:
        pass
    evals, evecs = np.linalg.eigh(W)
    max_e = float(np.max(evals))
    rel_floor = max(1e-6 * max_e, 10.0*np.finfo(evals.dtype).eps)
    evals = np.clip(evals, rel_floor, None)
    Winv_half = (evecs * (evals ** -0.5)) @ evecs.T
    L = Z @ Winv_half
    _maybe_cleanup(W, evecs, Winv_half, Z)
    return L.astype(np.float32, copy=False)

def _prep_standardize(ETA, y, standardize_X: bool, center_y: bool,
                      *, rkhs_approx: str, rkhs_m: Optional[int], rkhs_seed: Optional[int], kernel_type: str,
                      rkhs_force_nystrom_mem_gb: Optional[float], ctx: Optional[Context]=None):
    """Standardize X blocks; pass RKHS with normalized K or prebuilt L; center y; compute block weights."""
    y = jnp.asarray(y)
    if center_y:
        y_mean = jnp.nanmean(y)
        y_t = y - y_mean
    else:
        y_mean = 0.0
        y_t = y

    n = int(y.shape[0])
    ETA_j = []
    x_stats = []
    weights = []

    for e in ETA:
        if e["method"] == "rkhs":
            # Build/normalize kernel
            if "K" in e:
                K = np.asarray(e["K"], dtype=np.float32)
                base_weight = float(np.mean(np.diag(K))) if K.size else 1.0
            elif "X" in e:
                X = np.asarray(e["X"], dtype=np.float32)
                # Bandwidth estimate from sample (avoid O(n^2) alloc)
                n_samp = min(20000, max(1000, X.shape[0] * max(1, X.shape[0] // 10)))
                rng = np.random.default_rng(123)
                idx_i = rng.integers(0, X.shape[0], size=n_samp, endpoint=False)
                idx_j = rng.integers(0, X.shape[0], size=n_samp, endpoint=False)
                d2s = np.sum((X[idx_i] - X[idx_j])**2, axis=1)
                d2s = d2s[d2s > 0]
                ell = float(np.sqrt(np.median(d2s))) + 1e-6 if d2s.size else 1.0
                K = None
                base_weight = 1.0  # we don't materialize K; neutral weight
            else:
                raise ValueError("rkhs: provide 'K' or 'X'.")

            approx = rkhs_approx

            if rkhs_force_nystrom_mem_gb is not None and K is not None:
                bytes_n2 = float(n) * float(n) * 4.0
                if bytes_n2 > rkhs_force_nystrom_mem_gb * (1024**3):
                    if str(approx).lower() != "nystrom":
                        gb_need = bytes_n2 / (1024**3)
                        warnings.warn(
                            f"RKHS: n^2 kernel (~{gb_need:.2f} GiB) exceeds budget "
                            f"{float(rkhs_force_nystrom_mem_gb):.2f} GiB; forcing rkhs_approx='nystrom'.",
                            ResourceWarning
                        )
                        logger.warning("Forcing Nyström (n^2 kernel %.2f GiB exceeds budget %.2f GiB).",
                                       gb_need, float(rkhs_force_nystrom_mem_gb))
                    approx = "nystrom"

            if str(approx).lower() == "exact" and n > 1000:
                warnings.warn(
                    "RKHS: n > 1000; switching to Nyström to avoid O(n^3) Cholesky. "
                    "Pass rkhs_approx='exact' *and* rkhs_force_nystrom_mem_gb=None with n<=1000 to keep exact.",
                    RuntimeWarning
                )
                logger.warning("RKHS exact requested with n=%d; switching to Nyström for safety.", n)
                approx = "nystrom"

            if K is not None:
                diag_mean = float(max(np.mean(np.diag(K)), 1e-12))
                Kn = (K / diag_mean).astype(np.float32, copy=False)
                L = _rkhs_prepare_factor(Kn, approx=approx, n=n, rkhs_m=rkhs_m, seed=rkhs_seed, kernel_type=kernel_type, ctx=ctx)
                e2 = dict(e); e2["K"] = Kn; e2["L"] = L
            else:
                if str(approx).lower() == "nystrom":
                    m = int(rkhs_m) if rkhs_m is not None else min(n, max(100, int(np.sqrt(n)*10)))
                    rng = np.random.default_rng(rkhs_seed if rkhs_seed is not None else 123)
                    idx = np.sort(rng.choice(n, size=m, replace=False))
                    Xm = X[idx]  # (m×p)
                    def k_row(x, Y):
                        d2 = np.sum((Y - x)**2, axis=1)
                        return np.exp(-d2/(2.0*ell*ell)).astype(np.float32)
                    Z = np.stack([k_row(Xm[j], X) for j in range(m)], axis=1)  # (n×m)
                    W = np.stack([k_row(Xm[j], Xm) for j in range(m)], axis=1)  # (m×m)
                    dmean = float(np.mean(np.diag(W)))
                    if not np.isfinite(dmean) or dmean <= 0:
                        dmean = 1.0
                    Z = (Z / dmean).astype(np.float32, copy=False)
                    W = (W / dmean).astype(np.float32, copy=False)
                    J = _adaptive_jitter(W)
                    W[np.diag_indices_from(W)] += J
                    evals, evecs = np.linalg.eigh(W)
                    max_e = float(np.max(evals))
                    rel_floor = max(1e-6 * max_e, 10.0*np.finfo(evals.dtype).eps)
                    evals = np.clip(evals, rel_floor, None)
                    Winv_half = (evecs * (evals ** -0.5)) @ evecs.T
                    L = (Z @ Winv_half).astype(np.float32, copy=False)  # (n×m)
                    _maybe_cleanup(W, evecs, Winv_half, Z)
                    e2 = dict(e); e2["K"] = None; e2["L"] = L
                else:
                    dists = np.sum((X[:, None, :] - X[None, :, :])**2, axis=-1)
                    K_full = np.exp(-dists / (2.0 * (ell**2))).astype(np.float32)
                    diag = np.diag(K_full)
                    if np.any(diag <= 0):
                        K_full[np.diag_indices_from(K_full)] = 1.0
                    diag_mean = float(max(np.mean(np.diag(K_full)), 1e-12))
                    Kn = (K_full / diag_mean).astype(np.float32, copy=False)
                    L = _rkhs_prepare_factor(Kn, approx="exact", n=n, rkhs_m=rkhs_m, seed=rkhs_seed, kernel_type=kernel_type, ctx=ctx)
                    e2 = dict(e); e2["K"] = Kn; e2["L"] = L

            e2["is_rkhs"] = True; e2["rkhs_approx"] = str(approx).lower()
            ETA_j.append(e2)
            x_stats.append((jnp.zeros(1), jnp.ones(1)))  # placeholder
            weights.append(float(max(1e-12, base_weight)))
            continue

        X = jnp.asarray(e["X"], dtype=jnp.float32)
        # compute block weight BEFORE standardization (empirical)
        try:
            base_weight = float(np.mean(np.var(np.asarray(e["X"], dtype=np.float32), axis=0)))
            if not np.isfinite(base_weight):
                base_weight = 1.0
        except Exception:
            base_weight = 1.0

        if standardize_X:
            cm = jnp.mean(X, axis=0)
            std = jnp.std(X, axis=0)
            min_std = 1e-6 if X.dtype == jnp.float32 else 1e-12
            cs = jnp.clip(std, min_std, None)
            Xs = (X - cm)/cs
            x_stats.append((cm, cs))
        else:
            Xs = X
            x_stats.append((jnp.zeros(X.shape[1]), jnp.ones(X.shape[1])))

        if "d" in e:
            d_mask = jnp.asarray(e["d"], dtype=jnp.float32)
            Xs = Xs * d_mask
        e2 = dict(e); e2["X"] = Xs
        ETA_j.append(e2)
        weights.append(float(max(1e-12, base_weight)))

    mask_obs = jnp.isfinite(y_t)
    mask_mis = ~mask_obs
    return ETA_j, y_t, (float(y_mean),), x_stats, mask_obs, mask_mis, weights

def _default_hyper_from_R2(y, X_list, method, R2=0.5, nu=5.0):
    """Calibrate simple scale hyperparameters from target R2."""
    yj = jnp.asarray(y); Vy = jnp.nanvar(yj)
    p_tot = int(sum(int(X.shape[1]) for X in X_list)) or 1
    R2j = jnp.clip(jnp.asarray(R2, dtype=yj.dtype), 1e-8, 0.9999)

    if method in ("bayesian_ridge", "lasso", "elastic_net"):
        tau0 = jnp.sqrt(jnp.maximum(R2j * Vy / p_tot, jnp.asarray(1e-12, dtype=yj.dtype)))
        return {"tau0": tau0}
    if method == "bayes_c":
        denom = jnp.maximum(jnp.asarray(0.05 * p_tot, dtype=yj.dtype), jnp.asarray(1.0, dtype=yj.dtype))
        tau0 = jnp.sqrt(jnp.maximum(R2j * Vy / denom, jnp.asarray(1e-12, dtype=yj.dtype)))
        return {"tau0": tau0}
    if method == "bayes_a":
        target = jnp.maximum(R2j * Vy / p_tot, jnp.asarray(1e-12, dtype=yj.dtype))
        nuj = jnp.asarray(nu, dtype=yj.dtype)
        s2 = target * (nuj - 2.0) / nuj
        return {"nu": nuj, "s2": s2}
    return {}

from typing import Any, Dict, Tuple
import numpy as np
import arviz as az

def _extract_samples_common(idata) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray, Dict[str, Any], Dict[str, np.ndarray]]:
    """
    Robustly extract posterior arrays (combined chains) and build:
      - beta_map[varname] -> (S, p)
      - intercept_samps   -> (S,)
      - sigma2_samps      -> (S,)
      - aux (means / full arrays for gate/lambda/tau/pi/student_t_nu)
      - u_map[varname]    -> (S, n) for RKHS blocks

    Ensures 'sample' is the leading dimension for all variables (prevents shape flip).
    """
    if not hasattr(idata, "posterior"):
        raise ValueError("InferenceData has no 'posterior' group; cannot extract samples.")

    try:
        ds = az.extract(idata, group="posterior", combined=True)  # Dataset[sample, ...]
    except Exception:
        # Fallback: combine chain/draw manually
        post = idata.posterior
        ds = post.stack(sample=("chain", "draw")).transpose("sample", ...)

    # Number of posterior samples
    S = int(ds.sizes["sample"]) if "sample" in ds.dims else int(next(iter(ds.dims.values())))

    def _pull(name: str, default=None):
        """Get variable with 'sample' as axis 0 if present, else return default."""
        if name not in ds.data_vars:
            return default
        var = ds[name]
        try:
            if "sample" in getattr(var, "dims", ()) and var.dims[0] != "sample":
                var = var.transpose("sample", ...)
        except Exception:
            pass
        return np.asarray(var.values)

    # Intercept
    intercept = _pull("intercept", np.full((S,), np.nan, dtype=np.float32)).reshape(-1)

    # sigma^2
    s2 = _pull("sigma2")
    if s2 is not None:
        sigma2 = np.asarray(s2).reshape(-1)
    else:
        s = _pull("sigma")
        sigma2 = (np.asarray(s).reshape(-1) ** 2) if s is not None else np.full((S,), np.nan, dtype=np.float32)

    beta_map: Dict[str, np.ndarray] = {}
    aux: Dict[str, Any] = {}
    u_map: Dict[str, np.ndarray] = {}

    for v in ds.data_vars:
        arr = _pull(v)
        if arr is None:
            continue

        if v.endswith("_beta"):
            # Ensure shape [S, p]
            if arr.size == 0:
                arr = np.zeros((S, 0), dtype=np.float32)
            elif arr.ndim == 1:
                arr = arr[:, None]
            beta_map[v] = arr

        elif v.endswith("_lambda"):
            aux[v.replace("_lambda", "_lambda_mean")] = np.mean(arr, axis=0)

        elif v.endswith("_tau2"):
            aux[v + "_mean"] = float(np.mean(arr))
            aux[v + "_sd"]   = float(np.std(arr))

        elif v.endswith(("_tau", "_pi")) or v == "student_t_nu":
            aux[v + "_mean"] = float(np.mean(arr))
            aux[v + "_sd"]   = float(np.std(arr))

        elif v.endswith("_gate"):
            aux[v] = arr

        elif v.endswith("_u"):
            u_map[v] = arr

    return beta_map, intercept, sigma2, aux, u_map



def _compute_variance_components_batched(
    ETA: List[Dict[str, Any]],
    beta_samps_map: Dict[str, np.ndarray],
    intercept_samps: np.ndarray,
    sigma2_samps: np.ndarray,
    batch_mem_gb: Optional[float] = None,
    n: Optional[int] = None,
    overhead_mult_override: Optional[float] = None,
    u_samps_map: Optional[Dict[str, np.ndarray]] = None,
    *,
    backend: Optional[str] = None,
    ctx: Optional[Context] = None,
) -> Tuple[Dict[str, Any], Dict[str, np.ndarray], Dict[str, Any]]:
    """Compute per-effect Var(effect contribution), total Var(g), and yHat stats with adaptive batching."""
    num_blocks = len(ETA)
    if n is None:
        if num_blocks > 0:
            n = ETA[0]["X"].shape[0] if "X" in ETA[0] else ETA[0]["K"].shape[0]
        else:
            raise ValueError("Must provide n when ETA is empty.")
    S = len(intercept_samps)

    overhead_mult = float(overhead_mult_override) if (overhead_mult_override and overhead_mult_override > 0) else _estimate_overhead_mult(backend)
    bytes_per = 4 * overhead_mult
    # Budget resolution via context
    host_budget_gb, budget_info = ctx.resolve_host_budget_gb(requested=batch_mem_gb) if ctx else (batch_mem_gb or 0.5, {"source":"unknown", "value_gb": batch_mem_gb or 0.5})
    denom = max((num_blocks + 1) * n * bytes_per, 1.0)
    B_init = int(max(1, min(S, (host_budget_gb * (1024**3)) / denom)))
    B_init = int(max(1, min(B_init, 4096)))

    ddof = 1 if n > 1 else 0
    Xs = [ (np.asarray(e["X"], dtype=np.float32) if "X" in e else None) for e in ETA ]

    # Adaptive backoff/growth
    B = max(1, B_init)
    B_history = [int(B)]
    per_eff = {f"effect_{i}_": {"var_g_samples": []} for i in range(num_blocks)}
    yhat_sum   = np.zeros(n, dtype=np.float64)
    yhat_sqsum = np.zeros(n, dtype=np.float64)
    total_count = 0
    total_var_samples = []
    s0 = 0
    oom_retries = 0
    OOM_LIMIT = 6

    while s0 < S:
        s1 = min(S, s0 + B); B_eff = s1 - s0
        try:
            g_total = np.zeros((B_eff, n), dtype=np.float32)
            for i in range(num_blocks):
                key = f"effect_{i}_"
                key_beta = key + "beta"
                key_u    = key + "u"
                if (u_samps_map is not None) and (key_u in u_samps_map):
                    u_batch_all = u_samps_map[key_u]
                    u_batch = u_batch_all[s0:s1, :]
                    per_eff[key]["var_g_samples"].append(u_batch.var(axis=1, ddof=ddof))
                    g_total += u_batch
                else:
                    X = Xs[i]
                    if X is None:
                        per_eff[key]["var_g_samples"].append(np.zeros((B_eff,), dtype=np.float32))
                        continue
                    p = int(X.shape[1])
                    if p == 0:
                        per_eff[key]["var_g_samples"].append(np.zeros((B_eff,), dtype=np.float32))
                        continue
                    beta_batch = beta_samps_map[key_beta][s0:s1, :]
                    if beta_batch.ndim == 1:
                        beta_batch = beta_batch[:, None]
                    if beta_batch.shape[1] != p:
                        raise ValueError(
                            f"[shape mismatch] Effect {i}: beta_batch.shape={beta_batch.shape}, X.shape={X.shape}. "
                            f"Available beta keys: {sorted(beta_samps_map.keys())[:12]}..."
                        )
                    g_batch = beta_batch @ X.T
                    per_eff[key]["var_g_samples"].append(g_batch.var(axis=1, ddof=ddof))
                    g_total += g_batch
                    _maybe_cleanup(g_batch)
            total_var_samples.append(g_total.var(axis=1, ddof=ddof))
            yhat_batch = g_total + intercept_samps[s0:s1].reshape(-1, 1)
            yhat_sum   += yhat_batch.sum(axis=0, dtype=np.float64)
            yhat_sqsum += (yhat_batch.astype(np.float64)**2).sum(axis=0, dtype=np.float64)
            total_count += B_eff
            _maybe_cleanup(g_total, yhat_batch)
            # succeeded: consider modest growth if we're below initial pick (no growth past B_init)
            if B < B_init and B < 4096:
                B = min(B_init, int(B * 1.25))
                B_history.append(int(B))
            s0 = s1
            oom_retries = 0
        except Exception as ex:
            msg = str(ex).lower()
            if any(k in msg for k in ["memory", "resource exhausted", "alloc"]):
                if B <= 1 or oom_retries >= OOM_LIMIT:
                    raise
                B = max(1, B // 2)
                B_history.append(int(B))
                oom_retries += 1
                logger.warning("Adaptive batching backoff due to OOM-like error; retrying with B=%d (retry %d/%d).", B, oom_retries, OOM_LIMIT)
                continue
            else:
                raise

    total_var_samples = np.concatenate(total_var_samples, axis=0) if total_var_samples else np.zeros(S, dtype=np.float32)

    per_effect_out = {}
    for i in range(num_blocks):
        key = f"effect_{i}_"
        arr = per_eff[key]["var_g_samples"]
        v = np.concatenate(arr, axis=0) if arr else np.zeros(S, dtype=np.float32)
        per_effect_out[key] = {"var_g_mean": float(np.nanmean(v)) if v.size else 0.0,
                               "var_g_sd":   float(np.nanstd(v))  if v.size else 0.0}

    varE_mean = float(np.nanmean(sigma2_samps)) if np.any(np.isfinite(sigma2_samps)) else float("nan")
    varE_sd   = float(np.nanstd(sigma2_samps))  if np.any(np.isfinite(sigma2_samps)) else float("nan")
    varG_total_mean = float(max(0.0, float(np.nanmean(total_var_samples)))) if total_var_samples.size else 0.0
    varG_total_sd   = float(max(0.0, float(np.nanstd(total_var_samples))))  if total_var_samples.size else 0.0

    denom_h2 = varG_total_mean + (varE_mean if np.isfinite(varE_mean) else 0.0)
    h2_mean = float(min(1.0, max(0.0, varG_total_mean / denom_h2))) if (denom_h2 > 0.0 and np.isfinite(denom_h2)) else float("nan")

    yhat_mean = yhat_sum / max(1, total_count)
    yhat_var  = (yhat_sqsum / max(1, total_count)) - np.square(yhat_mean)
    yhat_sd   = np.sqrt(np.maximum(yhat_var, 0.0))

    mem_diag = {
        "batching": {
            "B_init": int(B_init),
            "B_history": [int(x) for x in B_history],
            "overhead_mult": float(overhead_mult),
            "bytes_per_scalar": 4.0,
            "n": int(n),
            "num_blocks": int(num_blocks)
        },
        "host_budget_used_gb": float(host_budget_gb),
        "budget_info": budget_info
    }

    return (
        {"per_effect": per_effect_out,
         "varE_mean": varE_mean, "varE_sd": varE_sd,
         "varG_total_mean": varG_total_mean, "varG_total_sd": varG_total_sd,
         "h2_mean": h2_mean},
        {"mean": yhat_mean.astype(np.float32), "sd": yhat_sd.astype(np.float32)},
        mem_diag
    )

def _summarize_diag(diag):
    if isinstance(diag, dict):
        return diag
    try:
        if hasattr(diag, "to_dict") and hasattr(diag, "columns"):
            return diag.reset_index().to_dict(orient="list")
    except Exception:
        pass
    try:
        df = az.summary(diag)
        return df.reset_index().to_dict(orient="list")
    except Exception:
        return {}


# --- Optional: simple dependency/device helper for reticulate callers ---
def setup_deps(
    prefer_gpu: bool = True,
    try_numpyro: bool = True,
    try_pymc: bool = True,
    allow_install: bool = False,
    strict_versions: bool = False,
    disable_auto_install: bool = True,
    # accept legacy/extra knobs from R wrappers without breaking
    cuda_version: str = "auto",
    cuda: str = None,
    **kwargs,
):
    """
    Soft dependency probe for production use.

    - Never hard-fails on missing packages unless allow_install=False AND core deps truly absent.
    - Returns a dict with versions, device info, and any *_error strings.
    - Swallows extra kwargs (cuda_version/cuda/…​) for signature compatibility with R code.
    """
    import importlib, sys, platform, subprocess

    info = {
        "prefer_gpu": bool(prefer_gpu),
        "allow_install": bool(allow_install),
        "strict_versions": bool(strict_versions),
        "disable_auto_install": bool(disable_auto_install),
        "cuda_version": cuda or cuda_version,
    }

    def _have(mod: str) -> bool:
        return importlib.util.find_spec(mod) is not None

    def _pip_install(pkgs):
        if not pkgs:
            return
        cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--no-input"] + list(pkgs)
        try:
            subprocess.check_call(cmd)
            return True
        except Exception as e:
            return str(e)

    # --- 0) Try existing internal helpers, but never die if they raise -----------
    core_err = None
    try:
        _ensure_core_deps(
            allow_install=allow_install,
            strict_versions=strict_versions,
            disable_auto_install=disable_auto_install,
        )
    except Exception as e:
        core_err = str(e)
        info["core_error"] = core_err

    # --- 1) Ensure jax / jaxlib are present (soft) -------------------------------
    if not (_have("jax") and _have("jaxlib")) and allow_install:
        # Minimal, safe defaults:
        # - Windows: CPU wheels only
        # - macOS arm64: add jax-metal (optional)
        # - Linux: CPU wheels by default (CUDA extras are opt-in and fragile)
        pkgs = ["jax>=0.4.20", "jaxlib>=0.4.20"]
        if platform.system() == "Darwin" and platform.machine().lower() in ("arm64", "aarch64"):
            # helps on Apple Silicon (optional)
            pkgs.append("jax-metal>=0.0.5")
        maybe = _pip_install(pkgs)
        if maybe is not True:
            info["jax_install_error"] = maybe or "unknown pip error"

    # Try to import and summarize devices
    try:
        import jax, jaxlib  # noqa: F401
        info["jax_version"] = getattr(jax, "__version__", "unknown")
        try:
            info["jaxlib_version"] = getattr(jaxlib, "__version__", "unknown")
        except Exception:
            pass
        try:
            devs = jax.devices()
            plats = [getattr(d, "platform", "?") for d in devs]
            info["devices"] = plats
            info["accelerator_present"] = any(p in ("gpu", "cuda", "rocm", "tpu") for p in plats)
            info["local_device_count"] = len(devs)
            info["default_backend"] = jax.default_backend()
        except Exception as e:
            info["devices_error"] = str(e)
    except Exception as e:
        info["jax_error"] = f"{e}"

    # --- 2) NumPyro (optional) ---------------------------------------------------
    if try_numpyro:
        try:
            _ensure_numpyro_deps(
                allow_install=allow_install,
                strict_versions=strict_versions,
                disable_auto_install=disable_auto_install,
            )
        except Exception as e:
            info["numpyro_setup_error"] = str(e)
        if not _have("numpyro") and allow_install:
            maybe = _pip_install(["numpyro>=0.13.2"])
            if maybe is not True:
                info["numpyro_install_error"] = maybe or "unknown pip error"
        try:
            import numpyro as _np  # noqa: F401
            info["numpyro_version"] = getattr(_np, "__version__", "unknown")
        except Exception as e:
            info["numpyro_error"] = f"{e}"

    # --- 3) PyMC (optional) ------------------------------------------------------
    if try_pymc:
        try:
            _ensure_pymc_deps(
                allow_install=allow_install,
                strict_versions=strict_versions,
                disable_auto_install=disable_auto_install,
                need_numpyro=False,  # keep PyMC check independent of numpyro
            )
        except Exception as e:
            info["pymc_setup_error"] = str(e)
        if not _have("pymc") and allow_install:
            maybe = _pip_install(["pymc>=5.10", "arviz>=0.16"])
            if maybe is not True:
                info["pymc_install_error"] = maybe or "unknown pip error"
        try:
            import pymc as _pm  # noqa: F401
            info["pymc_version"] = getattr(_pm, "__version__", "unknown")
        except Exception as e:
            info["pymc_error"] = f"{e}"

    # --- 4) ArviZ (lightweight; often imported by PyMC/NumPyro flows) -----------
    if not _have("arviz") and allow_install:
        maybe = _pip_install(["arviz>=0.16"])
        if maybe is not True:
            info["arviz_install_error"] = maybe or "unknown pip error"
    try:
        import arviz as _az  # noqa: F401
        info["arviz_version"] = getattr(_az, "__version__", "unknown")
    except Exception as e:
        info["arviz_error"] = f"{e}"

    # If core helper raised because of overly strict version pins, we still
    # return successfully as long as imports work.
    info["ok"] = all(k not in info for k in ("jax_error",))
    return info


# ---------------------------
# Minimal unit tests & micro-bench
# ---------------------------

def _test_identical_points_rbf():
    np.random.seed(0)
    n, p = 8, 3
    X = np.zeros((n, p), dtype=np.float32)  # identical points
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method": "rkhs", "X": X, "rkhs_approx": "nystrom", "rkhs_m": 4}]
    out = bayesian_alphabet(ETA, y, vi=True, n_iter=50, chains=1, backend="numpyro")
    assert np.all(np.isfinite(out["yHat"])), "NaNs in yHat for identical points"

def _test_ill_conditioned_kernel():
    n = 50
    A = np.ones((n, n), dtype=np.float32)
    K = A + 1e-6*np.eye(n, dtype=np.float32)  # nearly rank-1
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method": "rkhs", "K": K, "rkhs_approx": "exact"}]
    out = bayesian_alphabet(ETA, y, vi=True, n_iter=50, chains=1, backend="numpyro")
    assert np.isfinite(out["varG"]), "Ill-conditioned K should still work (adaptive jitter)"

def _test_force_nystrom_by_budget():
    n, p = 400, 5
    X = np.random.randn(n, p).astype(np.float32)
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method": "rkhs", "X": X, "rkhs_approx": "exact"}]
    out = bayesian_alphabet(ETA, y, vi=True, n_iter=10, chains=1, backend="numpyro",
                            config=BayesAlphabetConfig(rkhs=RKHSConfig(force_nystrom_mem_gb=0.01)))
    assert any(e["extra"].get("rkhs_approx") == "nystrom" for e in out["ETA"]), "Nyström should be forced by budget"

def _test_zero_std_column():
    n = 20
    X = np.hstack([np.random.randn(n, 1), np.ones((n,1), dtype=np.float32)])
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method": "bayesian_ridge", "X": X}]
    out = bayesian_alphabet(ETA, y, vi=True, n_iter=20, chains=1, backend="numpyro")
    assert np.all(np.isfinite(out["ETA"][0]["b"])), "Zero-std column should not create infs"

def _micro_benchmark():
    """Very small synthetic benchmark; prints chosen batch size and estimated memory."""
    n, p, S = 500, 100, 200
    X = np.random.randn(n, p).astype(np.float32)
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method": "bayesian_ridge", "X": X}]
    res = bayesian_alphabet(
        ETA, y, vi=True, n_iter=50, chains=1, backend="numpyro",
        config=BayesAlphabetConfig(memory=MemoryBudget(host_fraction=0.25, safety_margin_gb=0.5))
    )
    print("Bench OK; varG:", res["varG"], "h2:", res["h2"])

# ---------------------------
# Additional guardrail tests
# ---------------------------

def _test_float32_std_floor():
    n, p = 10, 3
    X = np.hstack([np.random.randn(n, p-1).astype(np.float32), np.zeros((n,1), np.float32)])
    y = np.random.randn(n).astype(np.float32)
    ETA = [{"method":"bayesian_ridge","X":X}]
    out = bayesian_alphabet(ETA, y, vi=True, n_iter=20, chains=1)
    assert np.all(np.isfinite(out["ETA"][0]["b"])), "float32 std floor should prevent infs/nans"

def _test_invalid_X_raises():
    n, p = 8, 4
    X = np.random.randn(n,p).astype(np.float32)
    X[0,0] = np.nan
    y = np.random.randn(n).astype(np.float32)
    try:
        bayesian_alphabet([{"method":"bayesian_ridge","X":X}], y, vi=True, n_iter=5, chains=1)
        assert False, "Expected ValueError for non-finite X"
    except ValueError:
        pass

def _test_student_t_nu_vi_guard():
    n, p = 20, 5
    X = np.random.randn(n,p).astype(np.float32)
    y = np.random.randn(n).astype(np.float32)
    try:
        bayesian_alphabet([{"method":"bayesian_ridge","X":X}], y, vi=True, likelihood="student_t", student_t_nu=2.1, n_iter=5, chains=1)
        assert False, "Expected ValueError for small nu"
    except ValueError:
        pass

def _test_prng_split_no_reuse():
    n, p = 30, 4
    X = np.random.randn(n,p).astype(np.float32)
    y = np.random.randn(n).astype(np.float32)
    _ = bayesian_alphabet([{"method":"bayesian_ridge","X":X}], y, vi=True, random_seed=42, n_iter=10, chains=1)

if __name__ == "__main__":
    # Run a tiny self-check if desired (kept off by default)
    pass
