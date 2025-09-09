# dl_models.py
from __future__ import annotations

# =========================
# Core imports / typing
# =========================
import os, sys, subprocess, importlib, platform, ctypes, shutil, math, random
from typing import Optional, Sequence, Tuple, List, Dict

# Avoid user-site pollution that can reintroduce mismatched wheels
os.environ.setdefault("PYTHONNOUSERSITE", "1")

# =========================
# Pip / system helpers
# =========================
def _run(cmd: Sequence[str]):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, text=True)

def _pip_install(pkgs: Sequence[str], index_url: Optional[str] = None, extra: Sequence[str] = ()):
    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--quiet", "--no-cache-dir"]
    if index_url:
        cmd += ["--index-url", index_url]
    cmd += list(pkgs) + list(extra)
    return _run(cmd)

def _pip_uninstall(pkgs: Sequence[str]):
    if not pkgs:
        return None
    cmd = [sys.executable, "-m", "pip", "uninstall", "-y", "--quiet"] + list(pkgs)
    return _run(cmd)

def _has_nvidia_gpu() -> bool:
    try:
        r = _run(["nvidia-smi"])
        return r.returncode == 0
    except Exception:
        return False

def _lazy_import(name: str, pip_name: Optional[str] = None):
    try:
        return importlib.import_module(name)
    except Exception:
        _pip_install([pip_name or name])
        return importlib.import_module(name)

def _windows_msvc_runtime_ok() -> bool:
    if platform.system() != "Windows":
        return True
    for dll in ("vcruntime140_1.dll", "vcruntime140.dll", "msvcp140.dll"):
        try:
            ctypes.CDLL(dll)
            return True
        except OSError:
            continue
    return False

def _torch_import_probe() -> Tuple[bool, str]:
    try:
        t = importlib.import_module("torch")
        _ = t.__version__
        _ = t.tensor([0.0]).sum().item()  # light op
        return True, ""
    except Exception as e:
        return False, repr(e)

def _pin_extras_to(torch_version: Optional[str], extra_packages: Tuple[str, ...]) -> List[str]:
    out: List[str] = []
    for p in (extra_packages or ()):
        if p in ("torchvision", "torchaudio") and torch_version:
            out.append(f"{p}=={torch_version if p=='torchaudio' else '0.19.0' if torch_version.startswith('2.4') else p}")
        else:
            out.append(p)
    return out

def _install_torch_with_fallback(
    prefer_gpu: bool,
    cuda: str,
    torch_version: Optional[str],
    extra_packages: Tuple[str, ...],
    forced_index: Optional[str],
) -> Dict[str, str]:
    """
    Try GPU wheel (only if GPU detected) then fall back to CPU wheel if import fails.
    Returns meta info including 'index_used', 'flavor', and 'logs'.
    """
    logs: List[str] = []
    sys_os = platform.system()

    # Only try GPU wheels when an NVIDIA GPU is really present, unless index is forced
    gpu_candidate = False
    if not forced_index:
        gpu_candidate = prefer_gpu and (cuda != "cpu") and _has_nvidia_gpu()

    # Build candidate installs: [(desc, index_url, flavor)]
    candidates: List[Tuple[str, Optional[str], str]] = []
    if forced_index:
        candidates.append(("forced", forced_index, "gpu" if "cu" in forced_index else "cpu"))
    else:
        if sys_os == "Darwin":
            # macOS: use default index; MPS handled by standard wheel
            candidates.append(("macos-default", None, "cpu"))
        else:
            if gpu_candidate:
                cu_tag = None
                if isinstance(cuda, str) and cuda.startswith("cu"):
                    cu_tag = cuda.lower()           # e.g., 'cu121'
                elif cuda == "auto":
                    cu_tag = "cu121"                # safe default for recent PyTorch
                if cu_tag:
                    candidates.append((f"pytorch-{cu_tag}", f"https://download.pytorch.org/whl/{cu_tag}", "gpu"))
            # Always queue a CPU fallback
            candidates.append(("pytorch-cpu", "https://download.pytorch.org/whl/cpu", "cpu"))

    def pkg_spec() -> List[str]:
        base = [f"torch=={torch_version}"] if torch_version else ["torch"]
        base += _pin_extras_to(torch_version, extra_packages)
        return base

    last_error = ""
    for desc, index_url, flavor in candidates:
        # Clean slate to avoid ABI mismatches
        _pip_uninstall(["torch", "torchvision", "torchaudio", "torchtext", "functorch"])
        r = _pip_install(pkg_spec(), index_url=index_url)
        logs.append(f"[{desc}] pip stdout:\n{r.stdout or ''}")

        ok, err = _torch_import_probe()
        if ok:
            return {
                "index_used": index_url or "default",
                "flavor": flavor,
                "logs": "\n".join(logs)[-8000:],
                "attempt": desc,
            }

        # Special self-heal for CUDA symbol mixups that can linger on Windows
        bad_cuda_sym = (
            "profiler_allow_cudagraph_cupti_lazy_reinit_cuda12" in err
            or "cupti" in err.lower()
            or "cudagraph" in err.lower()
        )
        if bad_cuda_sym:
            _pip_uninstall(["torch", "torchvision", "torchaudio", "torchtext", "functorch"])
            _run([sys.executable, "-m", "pip", "cache", "purge"])
            r2 = _pip_install(
                ["torch==2.4.0", "torchvision==0.19.0", "torchaudio==2.4.0"],
                index_url="https://download.pytorch.org/whl/cpu",
                extra=["--no-cache-dir"],
            )
            logs.append(f"[cpu-hard-reset] pip stdout:\n{r2.stdout or ''}")
            ok2, err2 = _torch_import_probe()
            if ok2:
                return {
                    "index_used": "https://download.pytorch.org/whl/cpu",
                    "flavor": "cpu-hard-reset",
                    "logs": "\n".join(logs)[-8000:],
                    "attempt": "cpu-hard-reset",
                }
            last_error = err2
            continue  # try next candidate if any

        # DLL issues on Windows? Let the CPU candidate run next.
        last_error = err
        if (("WinError 127" in err) or ("DLL load failed" in err) or ("shm.dll" in err)) and flavor != "cpu":
            continue

    raise RuntimeError(f"Unable to import torch after installs. Last error: {last_error}")

def _ensure_torch(
    prefer_gpu: bool = True,
    cuda: str = "auto",
    torch_version: Optional[str] = None,
    extra_packages: Tuple[str, ...] = ("torchvision", "torchaudio"),
    index_url: Optional[str] = None,
):
    """Ensure torch is importable. If broken, install a suitable wheel and import it."""
    try:
        t = importlib.import_module("torch")
        return t, {"index_used": None, "flavor": "existing", "logs": "", "attempt": "import-existing"}
    except Exception:
        pass

    forced_idx_env = (os.environ.get("PYTORCH_INDEX_URL") or "").strip() or None
    meta = _install_torch_with_fallback(
        prefer_gpu=prefer_gpu,
        cuda=cuda,
        torch_version=torch_version,
        extra_packages=extra_packages,
        forced_index=index_url or forced_idx_env,
    )
    t = importlib.import_module("torch")
    return t, meta

def _ensure_sklearn():
    try:
        return importlib.import_module("sklearn")
    except Exception:
        _pip_install(["scikit-learn"], index_url=None)
        return importlib.import_module("sklearn")

# =========================
# Early bootstrap: numpy / torch / nn
# =========================
try:
    import numpy as np
except Exception:
    np = _lazy_import("numpy", pip_name="numpy")

try:
    import torch
except Exception:
    torch, _ = _ensure_torch()

# MUST exist before any "class X(nn.Module)" executes
try:
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import DataLoader, TensorDataset, random_split
except Exception:
    torch = importlib.import_module("torch")
    nn = importlib.import_module("torch.nn")
    F = importlib.import_module("torch.nn.functional")
    _tud = importlib.import_module("torch.utils.data")
    DataLoader, TensorDataset, random_split = _tud.DataLoader, _tud.TensorDataset, _tud.random_split

try:
    torch.set_default_dtype(torch.float32)
except Exception:
    pass

# =========================
# Public env/setup entrypoint
# =========================
def setup_deps(
    prefer_gpu: bool = True,
    cuda: str = "auto",
    torch_version: Optional[str] = None,
    numpy_spec: str = "numpy>=1.24",
    extra_packages: Tuple[str, ...] = ("torchvision", "torchaudio"),
    index_url: Optional[str] = None,
    install_sklearn: bool = True,
):
    logs = []
    pkgs_installed = []

    # NumPy first (neutral index)
    r = _pip_install([numpy_spec])
    logs.append(r.stdout or "")
    np = _lazy_import("numpy", pip_name="numpy")
    pkgs_installed.append(f"numpy {np.__version__}")

    # Torch (+friends) with GPU->CPU fallback
    torch_mod, torch_meta = _ensure_torch(
        prefer_gpu=prefer_gpu, cuda=cuda, torch_version=torch_version,
        extra_packages=extra_packages, index_url=index_url
    )
    logs.append(torch_meta.get("logs", ""))
    pkgs_installed.append(f"torch {getattr(torch_mod, '__version__', 'unknown')}")
    for _m in ("torchvision", "torchaudio"):
        try:
            _mm = importlib.import_module(_m)  # tolerate failure for torchaudio on Win CPU
            pkgs_installed.append(f"{_m} {getattr(_mm, '__version__', 'unknown')}")
        except Exception:
            pass

    # Info block
    info = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "pytorch_version": getattr(torch_mod, "__version__", "unknown"),
        "cuda_compiled": getattr(getattr(torch_mod, "version", None), "cuda", None),
        "cuda_available": bool(getattr(torch_mod, "cuda", None) and torch_mod.cuda.is_available()),
        "gpu_name": (torch_mod.cuda.get_device_name(0) if getattr(torch_mod, "cuda", None) and torch_mod.cuda.is_available() else None),
        "mps_available": bool(getattr(torch_mod.backends, "mps", None) and torch_mod.backends.mps.is_available()),
        "index_used": torch_meta.get("index_used"),
        "install_flavor": torch_meta.get("flavor"),
        "packages_installed": pkgs_installed,
        "logs": "\n".join(filter(None, logs))[-8000:],
    }

    # cuDNN perf hint (safe try)
    if getattr(torch_mod.backends, "cudnn", None) and torch_mod.cuda.is_available():
        try:
            if not getattr(torch_mod, "are_deterministic_algorithms_enabled", lambda: False)():
                torch_mod.backends.cudnn.benchmark = True
        except Exception:
            torch_mod.backends.cudnn.benchmark = True

    # scikit-learn on default index
    if install_sklearn:
        try:
            sk = _ensure_sklearn()
            info["scikit_learn_version"] = getattr(sk, "__version__", "unknown")
            pkgs_installed.append(f"scikit-learn {info['scikit_learn_version']}")
        except Exception:
            info["scikit_learn_version"] = None

    # Windows VC++ runtime advisory
    if platform.system() == "Windows" and not _windows_msvc_runtime_ok():
        info["windows_hint"] = "Install the Microsoft Visual C++ Redistributable (2015–2022, x64)."

    return info

# =========================
# Utilities: device / determinism / compile
# =========================
def _as_int_list(x: Sequence) -> List[int]:
    return [int(v) for v in list(x)]

def _pick_device(device=None):
    try:
        if device:
            dev = device if isinstance(device, torch.device) else torch.device(device)
            if dev.type == "cuda" and not torch.cuda.is_available():
                return torch.device("cpu")
            return dev
    except Exception:
        return torch.device("cpu")
    if getattr(torch, "cuda", None) and torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")

def _safe_to_device(model: nn.Module, device: torch.device):
    try:
        model = model.to(device)
        return model, device
    except Exception:
        cpu = torch.device("cpu")
        try:
            model = model.to(cpu)
        except Exception:
            pass
        return model, cpu

def _set_deterministic(seed: int = 42):
    random.seed(int(seed))
    np.random.seed(int(seed))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))
    try:
        torch.use_deterministic_algorithms(True)
    except Exception:
        pass
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

def _maybe_compile(model, dev, compile_model=True):
    # Allow users to disable via env
    if not compile_model or not hasattr(torch, "compile"):
        return model
    if os.getenv("TORCH_COMPILE_DISABLE", "").lower() in {"1", "true", "yes"}:
        return model
    if os.getenv("PREDICTDL_COMPILE", "").lower() in {"0", "false", "no"}:
        return model

    # Choose backend safely
    if dev.type == "cuda":
        # On Windows, inductor needs cl.exe; skip if missing
        if platform.system() == "Windows" and shutil.which("cl") is None:
            return model
        backend = "inductor"
    elif dev.type == "mps":
        backend = "aot_eager"  # inductor not useful on MPS today
    else:
        backend = "aot_eager"  # CPU

    try:
        return torch.compile(model, backend=backend, mode="max-autotune")
    except Exception:
        # Don’t crash training; run eager
        try:
            import torch._dynamo as dynamo
            dynamo.config.suppress_errors = True
        except Exception:
            pass
        return model

def _as_int_list(x: Sequence) -> List[int]:
    return [int(v) for v in list(x)]

def _pick_device(device=None):
    try:
        if device:
            dev = device if isinstance(device, torch.device) else torch.device(device)
            if dev.type == "cuda" and not torch.cuda.is_available():
                return torch.device("cpu")
            return dev
    except Exception:
        return torch.device("cpu")
    if getattr(torch, "cuda", None) and torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")

def _safe_to_device(model: nn.Module, device: 'torch.device'):
    try:
        model = model.to(device)
        return model, device
    except Exception:
        cpu = torch.device("cpu")
        try:
            model = model.to(cpu)
        except Exception:
            pass
        return model, cpu


def _infer_task(y, multi_class_max=10):
    y = np.asarray(y)
    uniq = np.unique(y)
    ratio = len(uniq) / max(1, len(y))
    if len(uniq) <= multi_class_max and np.allclose(uniq, np.round(uniq)) and ratio < 0.2:
        return ("binary", 2) if len(uniq) == 2 else ("multiclass", len(uniq))
    return "regression", 1

def _set_deterministic(seed: int = 42):
    random.seed(int(seed))
    np.random.seed(int(seed))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))
    try:
        torch.use_deterministic_algorithms(True)
    except Exception:
        pass
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

# put near other imports
import os, shutil, platform
import torch

def _maybe_compile(model, dev, compile_model=True):
    # Allow users to disable via env
    if not compile_model or not hasattr(torch, "compile"):
        return model
    if os.getenv("TORCH_COMPILE_DISABLE", "").lower() in {"1", "true", "yes"}:
        return model
    if os.getenv("PREDICTDL_COMPILE", "").lower() in {"0", "false", "no"}:
        return model

    # Choose backend safely
    if dev.type == "cuda":
        # On Windows, inductor needs cl.exe; skip if missing
        if platform.system() == "Windows" and shutil.which("cl") is None:
            return model
        backend = "inductor"
    elif dev.type == "mps":
        backend = "aot_eager"  # inductor not useful on MPS today
    else:
        backend = "aot_eager"  # CPU

    try:
        return torch.compile(model, backend=backend, mode="max-autotune")
    except Exception:
        # Don’t crash training; run eager
        try:
            import torch._dynamo as dynamo
            dynamo.config.suppress_errors = True
        except Exception:
            pass
        return model
# ---------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------
class ResidualBlock(nn.Module):
    def __init__(self, in_dim, out_dim, dropout=0.0, batch_norm=True):
        super().__init__()
        in_dim, out_dim = int(in_dim), int(out_dim)
        self.proj = nn.Linear(in_dim, out_dim) if in_dim != out_dim else None
        layers = [nn.Linear(in_dim, out_dim)]
        if batch_norm: layers.append(nn.BatchNorm1d(out_dim))
        layers += [nn.ReLU(inplace=True)]
        if dropout and dropout > 0: layers.append(nn.Dropout(dropout))
        layers.append(nn.Linear(out_dim, out_dim))
        if batch_norm: layers.append(nn.BatchNorm1d(out_dim))
        self.net = nn.Sequential(*layers)
        self.act = nn.ReLU(inplace=True)
    def forward(self, x):
        identity = x if self.proj is None else self.proj(x)
        out = self.net(x)
        out = out + identity
        return self.act(out)

class MLPStack(nn.Module):
    def __init__(self, in_dim, hidden: Sequence[int], dropout=0.0, batch_norm=True):
        super().__init__()
        layers = []
        prev = int(in_dim)
        for h in _as_int_list(hidden):
            layers.append(nn.Linear(prev, h))
            if batch_norm: layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU(inplace=True))
            if dropout and dropout > 0: layers.append(nn.Dropout(dropout))
            prev = h
        self.out_dim = prev
        self.net = nn.Sequential(*layers)
    def forward(self, x): return self.net(x)

class FinalLayerAttention(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.att = nn.Linear(int(dim), int(dim), bias=True)
    def forward(self, h):
        gate = torch.softmax(self.att(h), dim=-1)
        return h * gate

class LayerwiseAttention(nn.Module):
    def __init__(self, layer_dims: Sequence[int], att_dim: Optional[int] = None):
        super().__init__()
        layer_dims = _as_int_list(layer_dims)
        if att_dim is None: att_dim = layer_dims[-1]
        att_dim = int(att_dim)
        self.projections = nn.ModuleList([nn.Linear(int(d), att_dim) for d in layer_dims])
        self.scorer = nn.Linear(att_dim, 1, bias=False)
        self.att_dim = att_dim
    def forward(self, layer_outputs: Sequence[torch.Tensor]):
        toks = [proj(h) for proj, h in zip(self.projections, layer_outputs)]
        H = torch.stack(toks, dim=1)       # (N, L, att_dim)
        scores = self.scorer(torch.tanh(H))
        alpha = torch.softmax(scores, dim=1)
        return (alpha * H).sum(dim=1)      # (N, att_dim)

# ---------------------------------------------------------------------
# Grouping/Tokenization
# ---------------------------------------------------------------------
class GroupTokenizer(nn.Module):
    """
    Map features -> tokens via group_index (F,) and compute per-token means:
        tokens[t] = mean_{i in group t} x_i
    Buffers:
      - group_index: (F,) long, mapping each feature to a token id in [0..T-1]
      - inv_counts:  (T,) float, 1 / (#features in token t)
    """
    def __init__(self, group_index: Optional[torch.LongTensor], group_counts: Optional[torch.Tensor]):
        super().__init__()
        if group_index is None or group_counts is None:
            self.register_buffer("group_index", None, persistent=False)
            self.register_buffer("inv_counts", None, persistent=False)
        else:
            gi = group_index.view(-1).long()
            inv = 1.0 / group_counts.view(-1).clamp_min(1.0).float()
            self.register_buffer("group_index", gi, persistent=True)   # (F,)
            self.register_buffer("inv_counts", inv, persistent=True)   # (T,)

    def has_groups(self) -> bool:
        return (self.group_index is not None) and (self.inv_counts is not None)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        if not self.has_groups():
            return X
        N, F = X.shape
        T = int(self.inv_counts.shape[0])
        idx = self.group_index.unsqueeze(0).expand(N, F)      # (N, F)
        tokens = X.new_zeros((N, T), dtype=X.dtype)
        tokens.scatter_add_(1, idx, X)                        # sum per token
        return tokens * self.inv_counts.to(X.dtype)           # mean per token

def _build_groups_from_chunks(n_features: int, init_group_size: int, max_tokens: int):
    gsize = max(1, int(init_group_size))
    n_tokens = int(min(math.ceil(n_features / gsize), max_tokens))
    group_index = np.repeat(np.arange(n_tokens, dtype=np.int64), gsize)[:n_features]
    counts = np.bincount(group_index, minlength=n_tokens)
    return group_index, counts

def _build_groups_kmeans(X_train: np.ndarray, n_features: int, init_group_size: int,
                         max_tokens: int, kmeans_batch: int, kmeans_iter: int, seed: int):
    try:
        from sklearn.cluster import MiniBatchKMeans
    except Exception:
        return None
    proposed = int(math.ceil(n_features / max(1, init_group_size)))
    n_tokens = int(np.clip(proposed, 64, max_tokens))
    Xf = X_train.astype(np.float32)
    eps = 1e-6
    mu = Xf.mean(axis=0, keepdims=True)
    sd = Xf.std(axis=0, keepdims=True) + eps
    Zt = ((Xf - mu) / sd).T  # (F, N)
    km = MiniBatchKMeans(n_clusters=n_tokens, batch_size=int(kmeans_batch),
                         max_iter=int(kmeans_iter), n_init=3, random_state=int(seed))
    km.fit(Zt)
    labels = km.labels_.astype(np.int64)
    counts = np.bincount(labels, minlength=n_tokens)
    return labels, counts

def _build_feature_groups(
    X_train: np.ndarray,
    n_features: int,
    group_trigger: int = 2048,
    init_group_size: int = 64,
    max_tokens: int = 1024,
    kmeans_batch: int = 4096,
    kmeans_iter: int = 100,
    seed: int = 42,
    method: str = "kmeans"
) -> Optional[List[np.ndarray]]:
    """
    Returns a list of index arrays (groups) or None if no grouping is done.
    """
    if n_features <= int(group_trigger):
        return None
    if method not in ("kmeans", "contiguous"):
        method = "kmeans"

    proposed = int(np.ceil(n_features / max(1, int(init_group_size))))
    n_tokens = int(np.clip(proposed, 64, int(max_tokens)))

    if method == "contiguous":
        groups = []
        step = int(max(1, math.floor(n_features / n_tokens)))
        for start in range(0, n_features, step):
            idxs = np.arange(start, min(n_features, start + step), dtype=np.int64)
            if idxs.size > 0:
                groups.append(idxs)
        return groups

    try:
        from sklearn.cluster import MiniBatchKMeans
        Xf = X_train.astype(np.float32)
        eps = 1e-6
        mu = Xf.mean(axis=0, keepdims=True)
        sd = Xf.std(axis=0, keepdims=True) + eps
        Z = (Xf - mu) / sd
        Zt = Z.T  # (F, N)
        km = MiniBatchKMeans(
            n_clusters=n_tokens, batch_size=int(kmeans_batch),
            max_iter=int(kmeans_iter), n_init=3, random_state=int(seed)
        )
        labels = km.fit_predict(Zt)
        groups = []
        for g in range(n_tokens):
            idxs = np.where(labels == g)[0]
            if idxs.size > 0:
                groups.append(idxs.astype(np.int64))
        if len(groups) == 0:
            return None
        return groups
    except Exception:
        # fallback to contiguous windows
        return _build_feature_groups(
            X_train, n_features, group_trigger, init_group_size, max_tokens,
            kmeans_batch, kmeans_iter, seed, method="contiguous"
        )

import torch.nn.functional as F

class GroupAggregator(nn.Module):
    def __init__(self, groups: List[np.ndarray]):
        super().__init__()
        if len(groups) == 0:
            self.register_buffer("_empty", torch.empty(0, dtype=torch.long), persistent=True)
            self._group_names = []
            self._is_empty = True
            return
        self._is_empty = False
        self._group_names = []
        for i, g in enumerate(groups):
            name = f"group_{i}"
            self._group_names.append(name)
            self.register_buffer(name, torch.as_tensor(g, dtype=torch.long), persistent=True)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self._is_empty:
            return x
        outs = []
        for name in self._group_names:
            idx = getattr(self, name)
            outs.append(x.index_select(1, idx).mean(dim=1, keepdim=True))
        return torch.cat(outs, dim=1)

# ---------------------------------------------------------------------
# Token models (FT-Transformer, SAINT-lite)
# ---------------------------------------------------------------------
class FeatureTokenizer(nn.Module):
    """
    Per-feature tokenizer (continuous only):
      token_i = x_i * W[i] + b[i] + col_emb[i]
    """
    def __init__(self, n_features: int, d_model: int, token_dropout: float = 0.1, use_cls: bool = True):
        super().__init__()
        self.n_features = int(n_features)
        self.d_model = int(d_model)
        self.weight = nn.Parameter(torch.randn(self.n_features, self.d_model) * 0.02)
        self.bias   = nn.Parameter(torch.zeros(self.n_features, self.d_model))
        self.col_emb = nn.Parameter(torch.randn(self.n_features, self.d_model) * 0.02)
        self.use_cls = bool(use_cls)
        self.cls = nn.Parameter(torch.zeros(1, 1, self.d_model)) if self.use_cls else None
        self.drop = nn.Dropout(token_dropout) if token_dropout and token_dropout > 0 else nn.Identity()

    def forward(self, x):  # x: (N, F)
        tokens = x.unsqueeze(-1) * self.weight + self.bias + self.col_emb  # (N, F, d)
        tokens = self.drop(tokens)
        if self.use_cls:
            cls = self.cls.expand(x.size(0), -1, -1)  # (N,1,d)
            tokens = torch.cat([cls, tokens], dim=1)  # (N, 1+F, d)
        return tokens

class ScalarTokenMLP(nn.Module):
    """Shared MLP to lift scalar tokens -> d_model."""
    def __init__(self, d_model: int, dropout: float):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, d_model), nn.GELU(), nn.Dropout(dropout), nn.Linear(d_model, d_model)
        )
    def forward(self, t):  # (N,T,1)
        return self.net(t)

class FTTransformer(nn.Module):
    def __init__(self, n_features: int, out_dim: int, task: str,
                 d_model: int = 192, n_heads: int = 8, n_layers: int = 4,
                 ff_mult: int = 4, dropout: float = 0.1,
                 token_dropout: float = 0.1, use_cls: bool = True,
                 group_index: Optional[torch.LongTensor] = None,
                 group_counts: Optional[torch.Tensor] = None,
                 use_scalar_tokenizer: bool = True):
        super().__init__()
        self.task = task
        self.group_tok = GroupTokenizer(group_index, group_counts)
        self.use_groups = self.group_tok.has_groups()
        self.use_cls = bool(use_cls)
        self.d_model = int(d_model)

        if self.use_groups and use_scalar_tokenizer:
            T = int(self.group_tok.inv_counts.numel())
            self.scalar = ScalarTokenMLP(self.d_model, token_dropout)
            self.col_emb = nn.Embedding(T, self.d_model)  # indices 0..T-1
            if self.use_cls:
                self.cls = nn.Parameter(torch.zeros(1, 1, self.d_model))
            self.token_dropout = nn.Dropout(token_dropout)
        else:
            self.tok = FeatureTokenizer(n_features, self.d_model, token_dropout, use_cls=self.use_cls)

        enc_layer = nn.TransformerEncoderLayer(
            d_model=self.d_model, nhead=int(n_heads),
            dim_feedforward=int(ff_mult) * self.d_model,
            dropout=float(dropout), activation='gelu', batch_first=True
        )
        self.enc = nn.TransformerEncoder(enc_layer, num_layers=int(n_layers))
        self.norm = nn.LayerNorm(self.d_model)
        self.head = nn.Linear(self.d_model, int(out_dim))

    def forward(self, X):
        if self.use_groups and hasattr(self, "scalar"):
            toks = self.group_tok(X)                       # (N, T)
            B, T = toks.shape
            h = self.scalar(toks.unsqueeze(-1))            # (N, T, d)
            ids = torch.arange(T, device=X.device).unsqueeze(0).expand(B, T)
            h = h + self.col_emb(ids)
            h = self.token_dropout(h)
            if self.use_cls:
                cls = self.cls.expand(B, 1, -1)
                z = torch.cat([cls, h], dim=1)             # (N, T+1, d)
            else:
                z = h                                      # (N, T, d)
        else:
            z = self.tok(X)                                # (N, F(+1 if CLS), d)

        z = self.enc(z)
        h = z[:, 0, :] if self.use_cls else z.mean(dim=1)
        h = self.norm(h)
        return self.head(h)

class SAINT(nn.Module):
    def __init__(self, n_features: int, out_dim: int, task: str,
                 d_model: int = 128, n_heads: int = 8, n_layers: int = 4,
                 ff_mult: int = 4, dropout: float = 0.1, token_dropout: float = 0.1,
                 group_index: Optional[torch.LongTensor] = None,
                 group_counts: Optional[torch.Tensor] = None,
                 use_scalar_tokenizer: bool = True,
                 use_cls: bool = True):
        super().__init__()
        self.task = task
        self.d_model = int(d_model)
        self.use_cls = bool(use_cls)

        self.group_tok = GroupTokenizer(group_index, group_counts)
        self.use_groups = self.group_tok.has_groups()

        if self.use_groups and use_scalar_tokenizer:
            T = int(self.group_tok.inv_counts.numel())
            self.scalar = ScalarTokenMLP(self.d_model, token_dropout)
            self.col_emb = nn.Embedding(T, self.d_model)
            if self.use_cls:
                self.cls = nn.Parameter(torch.zeros(1, 1, self.d_model))
            self.token_dropout = nn.Dropout(token_dropout)
        else:
            self.tok = FeatureTokenizer(n_features, self.d_model, token_dropout, use_cls=self.use_cls)

        enc_layer = nn.TransformerEncoderLayer(
            d_model=self.d_model, nhead=int(n_heads),
            dim_feedforward=int(ff_mult) * self.d_model,
            dropout=float(dropout), activation='gelu', batch_first=True
        )
        self.enc = nn.TransformerEncoder(enc_layer, num_layers=int(n_layers))
        self.norm = nn.LayerNorm(self.d_model)
        self.head = nn.Linear(self.d_model, int(out_dim))

    def forward(self, X):
        if self.use_groups and hasattr(self, "scalar"):
            toks = self.group_tok(X)                     # (N, T)
            B, T = toks.shape
            h = self.scalar(toks.unsqueeze(-1))          # (N, T, d)
            ids = torch.arange(T, device=X.device).unsqueeze(0).expand(B, T)
            h = h + self.col_emb(ids)
            h = self.token_dropout(h)
            if self.use_cls:
                cls = self.cls.expand(B, 1, -1)          # (N, 1, d)
                z = torch.cat([cls, h], dim=1)           # (N, T+1, d)
            else:
                z = h                                     # (N, T, d)
        else:
            z = self.tok(X)                               # (N, F(+1 if CLS), d)

        z = self.enc(z)
        h = z[:, 0, :] if self.use_cls else z.mean(dim=1) # CLS or mean pool
        h = self.norm(h)
        return self.head(h)

# ---------------------------------------------------------------------
# TabNet (sparsemax)
# ---------------------------------------------------------------------
class Sparsemax(nn.Module):
    def forward(self, input: torch.Tensor) -> torch.Tensor:
        z = input
        z_sorted, _ = torch.sort(z, descending=True, dim=-1)
        k = torch.arange(1, z.size(-1)+1, device=z.device, dtype=z.dtype).view(1, -1)
        z_cumsum = torch.cumsum(z_sorted, dim=-1)
        support = (1 + k * z_sorted) > z_cumsum
        k_z_int = support.sum(dim=-1, keepdim=True).clamp_min(1)
        tau = (z_cumsum.gather(-1, k_z_int.long()-1) - 1) / k_z_int.to(z.dtype)
        return torch.clamp(z - tau, min=0)

sparsemax = Sparsemax()

class GLU(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.lin = nn.Linear(int(in_dim), int(out_dim) * 2)
        self.bn = nn.BatchNorm1d(int(out_dim) * 2)
    def forward(self, x):
        a, b = self.bn(self.lin(x)).chunk(2, dim=-1)
        return a * torch.sigmoid(b)

class FeatureTransformer(nn.Module):
    def __init__(self, dim, n_glu=2):
        super().__init__()
        layers = []
        for _ in range(n_glu):
            layers += [GLU(dim, dim)]
        self.net = nn.Sequential(*layers)
    def forward(self, x): return self.net(x)

class AttentiveTransformer(nn.Module):
    def __init__(self, feat_dim, input_dim, relax=1.5):
        super().__init__()
        self.lin = nn.Linear(int(feat_dim), int(input_dim))
        self.bn = nn.BatchNorm1d(int(input_dim))
        self.relax = float(relax)
    def forward(self, h_feat, prior):
        scores = self.bn(self.lin(h_feat))
        scores = scores * prior
        mask = sparsemax(scores)
        new_prior = prior * (self.relax - mask)
        return mask, new_prior

class TabNet(nn.Module):
    def __init__(self, in_dim, out_dim, task,
                 steps=5, feature_dim=64, output_dim=64, gamma=1.5,
                 lambda_sparse=1e-4):
        super().__init__()
        self.in_dim = int(in_dim)
        self.initial = nn.Linear(int(in_dim), int(feature_dim))
        self.ft = FeatureTransformer(int(feature_dim), n_glu=2)
        self.atts = nn.ModuleList([AttentiveTransformer(int(feature_dim), int(in_dim), relax=gamma) for _ in range(int(steps))])
        self.dts  = nn.ModuleList([FeatureTransformer(int(feature_dim), n_glu=2) for _ in range(int(steps))])
        self.steps = int(steps)
        self.out = nn.Linear(int(feature_dim), int(out_dim))
        self.task = task
        self.lambda_sparse = float(lambda_sparse)
        self.register_buffer("_last_sparse_loss", torch.tensor(0.0), persistent=False)
    def extra_loss(self):
        return self.lambda_sparse * self._last_sparse_loss
    def forward(self, x):
        h0 = F.relu(self.initial(x))
        h0 = self.ft(h0)
        prior = torch.ones(x.size(0), self.in_dim, device=x.device, dtype=x.dtype)
        agg = torch.zeros_like(h0)
        sparse_loss = x.new_tensor(0.0)
        for att, dt in zip(self.atts, self.dts):
            mask, prior = att(h0, prior)
            x_masked = x * mask
            h_step = dt(F.relu(self.initial(x_masked)))
            agg = agg + h_step
            sparse_loss = sparse_loss + (mask * mask).sum(dim=1).mean()
        self._last_sparse_loss.copy_(sparse_loss.detach())
        return self.out(agg)

# ---------------------------------------------------------------------
# NODE-lite (soft oblivious decision ensembles)
# ---------------------------------------------------------------------
class SoftObliviousEnsemble(nn.Module):
    def __init__(self, in_dim, out_dim, n_trees=8, depth=3):
        super().__init__()
        self.in_dim = int(in_dim)
        self.out_dim = int(out_dim)
        self.n_trees = int(n_trees)
        self.depth = int(depth)
        self.splits = nn.ModuleList([nn.Linear(self.in_dim, self.n_trees) for _ in range(self.depth)])
        self.leaf_values = nn.Parameter(torch.randn(self.n_trees, 2 ** self.depth, self.out_dim) * 0.01)
    def forward(self, x):
        probs = torch.ones(x.size(0), self.n_trees, 1, device=x.device, dtype=x.dtype)
        for lin in self.splits:
            g = torch.sigmoid(lin(x))
            left = probs * g.unsqueeze(-1)
            right = probs * (1.0 - g).unsqueeze(-1)
            probs = torch.cat([left, right], dim=2)
        vals = self.leaf_values.unsqueeze(0)
        out = (probs.unsqueeze(-1) * vals).sum(dim=2).sum(dim=1)
        return out

class NODELite(nn.Module):
    def __init__(self, in_dim, out_dim, task, n_trees=8, depth=3):
        super().__init__()
        self.ensemble = SoftObliviousEnsemble(in_dim, out_dim, n_trees=n_trees, depth=depth)
        self.task = task
    def forward(self, x): return self.ensemble(x)

# ---------------------------------------------------------------------
# DeepFM / DCN
# ---------------------------------------------------------------------
class DeepFM(nn.Module):
    def __init__(self, in_dim, out_dim, task, k=16, deep_hidden=(128, 64), dropout=0.1, batch_norm=True):
        super().__init__()
        in_dim = int(in_dim); out_dim = int(out_dim)
        self.linear = nn.Linear(in_dim, out_dim)
        self.V = nn.Parameter(torch.randn(in_dim, int(k)) * 0.01)
        self.fm_head = nn.Linear(int(k), out_dim)
        self.deep = MLPStack(in_dim, _as_int_list(deep_hidden), dropout=dropout, batch_norm=batch_norm)
        last = self.deep.out_dim if len(_as_int_list(deep_hidden)) > 0 else in_dim
        self.deep_head = nn.Linear(last, out_dim)
        self.task = task
    def forward(self, x):
        lr_out = self.linear(x)
        xv = x @ self.V
        x2v2 = (x * x) @ (self.V * self.V)
        fm_vec = 0.5 * (xv * xv - x2v2)
        fm_out = self.fm_head(fm_vec)
        d_out = self.deep_head(self.deep(x))
        return lr_out + fm_out + d_out

class CrossLayer(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.w = nn.Parameter(torch.randn(int(dim)) * 0.01)
        self.b = nn.Parameter(torch.zeros(int(dim)))
    def forward(self, x0, xl):
        dot = torch.sum(xl * self.w, dim=1, keepdim=True)
        return x0 * dot + self.b + xl

class DCN(nn.Module):
    def __init__(self, in_dim, out_dim, task, n_layers=3, deep_hidden=(256,128), dropout=0.0, batch_norm=True):
        super().__init__()
        self.deep = MLPStack(int(in_dim), _as_int_list(deep_hidden), dropout=dropout, batch_norm=batch_norm)
        self.cross = nn.ModuleList([CrossLayer(int(in_dim)) for _ in range(int(n_layers))])
        last = self.deep.out_dim if len(_as_int_list(deep_hidden)) > 0 else int(in_dim)
        self.head = nn.Linear(int(in_dim) + last, int(out_dim))
        self.task = task
    def forward(self, x):
        x0 = x; xl = x
        for layer in self.cross: xl = layer(x0, xl)
        h = torch.cat([xl, self.deep(x)], dim=1)
        return self.head(h)

# ---------------------------------------------------------------------
# Classic MLP/ResNet/CNN
# ---------------------------------------------------------------------
class ResNetTabular(nn.Module):
    def __init__(self, in_dim, hidden: Sequence[int], dropout=0.0, batch_norm=True, out_dim=1, task="regression"):
        super().__init__()
        in_dim = int(in_dim); hidden = _as_int_list(hidden)
        blocks, prev = [], in_dim
        for h in hidden:
            blocks.append(ResidualBlock(prev, h, dropout=dropout, batch_norm=batch_norm))
            prev = h
        self.backbone = nn.Sequential(*blocks)
        self.head = nn.Linear(prev, int(out_dim))
        self.task = task
    def forward(self, x):
        return self.head(self.backbone(x))

class MLPNet(nn.Module):
    def __init__(self, in_dim, hidden: Sequence[int], dropout=0.0, batch_norm=True, out_dim=1, task="regression", final_attention=False):
        super().__init__()
        in_dim = int(in_dim); hidden = _as_int_list(hidden)
        self.stack = MLPStack(in_dim, hidden, dropout=dropout, batch_norm=batch_norm)
        self.final_att = FinalLayerAttention(self.stack.out_dim) if final_attention else None
        self.head = nn.Linear(self.stack.out_dim, int(out_dim))
        self.task = task
    def forward(self, x):
        h = self.stack(x)
        if self.final_att is not None: h = self.final_att(h)
        return self.head(h)

class MLPWithLayerAttention(nn.Module):
    def __init__(self, in_dim, hidden: Sequence[int], dropout=0.0, batch_norm=True, out_dim=1, task="regression"):
        super().__init__()
        in_dim = int(in_dim); hidden = _as_int_list(hidden)
        self.layers = nn.ModuleList(); self.bns = nn.ModuleList(); self.dropouts = nn.ModuleList()
        prev = in_dim
        for h in hidden:
            self.layers.append(nn.Linear(prev, h))
            self.bns.append(nn.BatchNorm1d(h) if batch_norm else nn.Identity())
            self.dropouts.append(nn.Dropout(dropout) if dropout and dropout > 0 else nn.Identity())
            prev = h
        self.att = LayerwiseAttention(layer_dims=hidden, att_dim=hidden[-1])
        self.head = nn.Linear(hidden[-1], int(out_dim))
        self.task = task
    def forward(self, x):
        outs = []; h = x
        for lin, bn, do in zip(self.layers, self.bns, self.dropouts):
            h = do(torch.relu(bn(lin(h))))
            outs.append(h)
        return self.head(self.att(outs))

class ConvResidualBlock(nn.Module):
    def __init__(self, channels: int, kernel_size: int, batch_norm: bool = True,
                 dropout: float = 0.0, same_supported: bool = True):
        super().__init__()
        pad = "same" if same_supported else (kernel_size // 2)
        def _conv():
            try:
                return nn.Conv1d(channels, channels, kernel_size, padding=pad)
            except TypeError:
                return nn.Conv1d(channels, channels, kernel_size, padding=kernel_size // 2)
        self.conv1 = _conv()
        self.bn1 = nn.BatchNorm1d(channels) if batch_norm else nn.Identity()
        self.conv2 = _conv()
        self.bn2 = nn.BatchNorm1d(channels) if batch_norm else nn.Identity()
        self.dropout = nn.Dropout(dropout) if dropout > 0 else nn.Identity()

    def forward(self, x):
        identity = x
        out = F.relu(self.bn1(self.conv1(x)))
        out = self.dropout(self.bn2(self.conv2(out)))
        return F.relu(out + identity)

# ---- small helpers ---------------------------------------------------------
class _SE1d(nn.Module):
    def __init__(self, channels: int, reduction: int = 16):
        super().__init__()
        r = max(1, int(channels // reduction))
        self.avg = nn.AdaptiveAvgPool1d(1)
        self.fc  = nn.Sequential(
            nn.Linear(channels, r), nn.ReLU(inplace=True),
            nn.Linear(r, channels), nn.Sigmoid()
        )
    def forward(self, x):                 # x: (N, C, L)
        w = self.avg(x).squeeze(-1)       # (N, C)
        w = self.fc(w).unsqueeze(-1)      # (N, C, 1)
        return x * w

def _norm1d(norm: str, c: int):
    norm = (norm or "batch").lower()
    if norm == "batch": return nn.BatchNorm1d(c)
    if norm == "group": return nn.GroupNorm(num_groups=min(32, c), num_channels=c)
    if norm == "layer": return nn.GroupNorm(num_groups=1, num_channels=c)  # LayerNorm over channels for 1D
    return nn.Identity()

def _same_pad_ks(k: int, d: int = 1) -> int:
    return (d * (k - 1)) // 2

class _SepConv1d(nn.Module):
    """Depthwise-separable 1D conv with 'same' padding fallback."""
    def __init__(self, in_c, out_c, k=3, stride=1, dilation=1):
        super().__init__()
        pad = _same_pad_ks(k, dilation)
        self.dw = nn.Conv1d(in_c, in_c, k, stride=stride, padding=pad, dilation=dilation, groups=in_c, bias=False)
        self.pw = nn.Conv1d(in_c, out_c, kernel_size=1, bias=False)
    def forward(self, x):
        return self.pw(self.dw(x))

# ---- upgraded residual block (optionally separable + SE + dilation) --------
class _ConvResBlock1D(nn.Module):
    def __init__(self, channels: int, kernel_size: int, norm_type: str = "batch",
                 dropout: float = 0.0, use_se: bool = True, reduction: int = 16,
                 separable: bool = False, dilation: int = 1):
        super().__init__()
        C = int(channels); K = int(kernel_size); D = int(dilation)
        pad = _same_pad_ks(K, D)

        Conv = _SepConv1d if separable else nn.Conv1d
        self.c1 = Conv(C, C, k=K, stride=1, dilation=D) if separable else nn.Conv1d(C, C, K, padding=pad, dilation=D)
        self.n1 = _norm1d(norm_type, C)
        self.c2 = Conv(C, C, k=K, stride=1, dilation=D) if separable else nn.Conv1d(C, C, K, padding=pad, dilation=D)
        self.n2 = _norm1d(norm_type, C)
        self.se = _SE1d(C, reduction) if use_se else nn.Identity()
        self.drop = nn.Dropout(float(dropout)) if dropout and dropout > 0 else nn.Identity()

    def forward(self, x):
        h = F.gelu(self.n1(self.c1(x)))
        h = self.drop(self.n2(self.c2(h)))
        h = self.se(h)
        return F.gelu(h + x)

# ---- upgraded CNN1DNet -----------------------------------------------------
class CNN1DNet(nn.Module):
    def __init__(
        self, in_len, conv_filters: Sequence[int], kernel_size: int = 3, dropout: float = 0.0,
        batch_norm: bool = True, dense_layers: Sequence[int] = (256, 128, 64),
        out_dim: int = 1, task: str = "regression",
        use_max_pool: bool = False, pool_kernel: int = 2, pool_stride: int = 2, pool_padding: int = 0,

        # NEW knobs (all optional)
        norm_type: str = None,                 # "batch" (default), "group", "layer", "none"
        use_se: bool = True, se_reduction: int = 16,
        separable: bool = True,                # depthwise-separable convs
        dilations: Optional[Sequence[int]] = None,  # e.g., [1,2,4,1,2,...]; default = all 1
        use_positional: bool = True,           # absolute position encoding (breaks shift equivariance)
        pool_type: str = "conv",               # "conv" (anti-alias), "max", "avg"
        use_global_pool: bool = True           # GAP head (robust) vs. flatten+MLP
    ):
        super().__init__()
        in_len = int(in_len)
        conv_filters = _as_int_list(conv_filters)
        kernel_size = int(kernel_size)
        self.task = task
        self.out_dim = int(out_dim)

        # Back-compat: respect 'batch_norm' unless norm_type is set explicitly
        if norm_type is None:
            norm_type = "batch" if batch_norm else "none"
        self._norm_type = norm_type

        # --- Positional embedding (absolute) ---
        self.pos_emb = nn.Parameter(torch.zeros(1, 1, in_len)) if use_positional else None
        if self.pos_emb is not None:
            nn.init.normal_(self.pos_emb, mean=0.0, std=0.02)

        # Dilation schedule
        if dilations is None or len(dilations) == 0:
            dilations = [1] * len(conv_filters)
        dilations = [int(d) for d in (dilations if isinstance(dilations, (list, tuple)) else [dilations])]
        if len(dilations) < len(conv_filters):
            dilations = (dilations * ((len(conv_filters) + len(dilations) - 1) // len(dilations)))[:len(conv_filters)]

        # --- Convolutional trunk ---
        convs = []
        prev_c = 1
        for i, (c_out, d) in enumerate(zip(conv_filters, dilations)):
            c_out = int(c_out)
            K = kernel_size
            pad = _same_pad_ks(K, d)

            # first layer in a stage: expand channels
            if i == 0 or prev_c != c_out:
                if separable and prev_c > 1:
                    # Use separable conv only after we have >=1 channels; first conv must go 1->C via standard conv
                    convs.append(_SepConv1d(prev_c, c_out, k=K, stride=1, dilation=d))
                else:
                    convs.append(nn.Conv1d(prev_c, c_out, kernel_size=K, padding=pad, dilation=d))
                convs.append(_norm1d(self._norm_type, c_out))
                convs.append(nn.GELU())
                if dropout and dropout > 0:
                    convs.append(nn.Dropout(float(dropout)))
            else:
                # residual refinement block with dilation/SE
                convs.append(_ConvResBlock1D(
                    channels=c_out, kernel_size=K, norm_type=self._norm_type, dropout=dropout,
                    use_se=use_se, reduction=int(se_reduction), separable=separable, dilation=d
                ))

            # Downsampling (anti-alias conv-pool or classic)
            if use_max_pool:
                pt = (pool_type or "conv").lower()
                if pt == "conv":
                    # depthwise conv with stride as anti-alias pooling
                    convs.append(nn.Conv1d(c_out, c_out, kernel_size=pool_kernel, stride=pool_stride,
                                           padding=pool_padding, groups=c_out, bias=False))
                elif pt == "avg":
                    convs.append(nn.AvgPool1d(kernel_size=pool_kernel, stride=pool_stride, padding=pool_padding))
                else:
                    convs.append(nn.MaxPool1d(kernel_size=pool_kernel, stride=pool_stride, padding=pool_padding))

            prev_c = c_out

        self.conv = nn.Sequential(*convs)

        # --- Head ---
        self.use_global_pool = bool(use_global_pool)
        self._dense_dropout = float(dropout) if dropout and dropout > 0 else 0.0
        self.dense_layers = _as_int_list(dense_layers)

        # build head shapes with a dry run on CPU (device-agnostic)
        with torch.no_grad():
            z = torch.zeros(1, 1, in_len)
            if self.pos_emb is not None:
                z = z + self.pos_emb
            z = self.conv(z)  # (1, C_last, L_last)
            self._build_head(z)

    def _build_head(self, z: torch.Tensor):
        C_last, L_last = int(z.shape[1]), int(z.shape[2])
        layers = []
        if self.use_global_pool:
            # GAP -> MLP -> out
            flat_dim = C_last
        else:
            # Flatten -> MLP -> out
            flat_dim = C_last * L_last

        prev = flat_dim
        for h in self.dense_layers:
            layers += [nn.Linear(prev, int(h)), nn.GELU()]
            if self._dense_dropout > 0: layers += [nn.Dropout(self._dense_dropout)]
            prev = int(h)

        self.mlp_head = nn.Sequential(*layers) if layers else nn.Identity()
        self.head = nn.Linear(prev, self.out_dim)

    def forward(self, x):
        # x: (N, F) -> (N, 1, F)
        x = x.unsqueeze(1)
        if self.pos_emb is not None and x.shape[-1] == self.pos_emb.shape[-1]:
            x = x + self.pos_emb

        z = self.conv(x)                           # (N, C, L)
        if self.use_global_pool:
            z = z.mean(dim=2)                      # (N, C)
        else:
            z = torch.flatten(z, start_dim=1)      # (N, C*L)

        z = self.mlp_head(z)
        return self.head(z)


# -----------------------------------------------------------------
# NAM / MoE / GP-DKL (+ RFF)
# -----------------------------------------------------------------
class TinyMLP(nn.Module):
    def __init__(self, in_dim=1, hidden=(32, 16), out_dim=1, dropout=0.0, act="relu"):
        super().__init__()
        acts = {"relu": nn.ReLU, "gelu": nn.GELU, "tanh": nn.Tanh, "silu": nn.SiLU}
        A = acts.get(str(act).lower(), nn.ReLU)
        layers, prev = [], int(in_dim)
        for h in _as_int_list(hidden):
            layers += [nn.Linear(prev, h), A()]
            if dropout and float(dropout) > 0:
                layers += [nn.Dropout(float(dropout))]
            prev = h
        layers += [nn.Linear(prev, int(out_dim))]
        self.net = nn.Sequential(*layers)
    def forward(self, x):  # (..., in_dim)
        return self.net(x)

class NAM(nn.Module):
    """
    Per-feature or per-group additive model:
      y_hat = sum_j f_j(g_j(x)) + optional linear term,
    where g_j(x) is either x_j (no grouping) or the mean of a feature group.
    """
    def __init__(self, in_features: int, out_dim: int, task: str,
                 subnet_hidden=(32, 16), activation="relu", dropout=0.0,
                 add_linear=True, l1_lambda=1e-4,
                 feature_groups: Optional[List[np.ndarray]] = None):
        super().__init__()
        self.task = task
        self.out_dim = int(out_dim)
        self.l1_lambda = float(l1_lambda)
        self.agg = GroupAggregator(feature_groups) if feature_groups is not None else None
        G = len(feature_groups) if feature_groups is not None else int(in_features)
        per_out = 1 if task in ("regression", "binary") else int(out_dim)
        self.subnets = nn.ModuleList([TinyMLP(1, subnet_hidden, per_out, dropout, activation) for _ in range(G)])
        self.add_linear = bool(add_linear)
        self.linear = nn.Linear(G, per_out) if self.add_linear else None

    def extra_loss(self):
        if self.l1_lambda <= 0:
            return 0.0
        reg = 0.0
        for m in self.subnets:
            for p in m.parameters():
                reg = reg + p.abs().sum()
        if self.linear is not None:
            for p in self.linear.parameters():
                reg = reg + p.abs().sum()
        return next(self.parameters()).new_tensor(self.l1_lambda) * reg

    def forward(self, x):   # x: (N, F)
        if self.agg is not None:
            x = self.agg(x)  # (N, G)
        contribs = []
        for i, sub in enumerate(self.subnets):
            gi = x[:, i:i+1]
            contribs.append(sub(gi))
        H = torch.stack(contribs, dim=1).sum(dim=1)  # (N, per_out)
        if self.linear is not None:
            H = H + self.linear(x)
        return H

class MoE(nn.Module):
    def __init__(self, in_dim, out_dim, task,
                 n_experts=4, expert_hidden=(128, 64),
                 gate_hidden=128, dropout=0.1, temperature=1.0,
                 sparse_topk: Optional[int] = None,
                 entropy_reg: float = 0.0,
                 feature_groups: Optional[List[np.ndarray]] = None,
                 batch_norm=True):
        super().__init__()
        self.task = task
        self.out_dim = int(out_dim)
        self.temperature = float(max(1e-3, temperature))
        self.sparse_topk = int(sparse_topk) if sparse_topk not in (None, 0) else None
        self.entropy_reg = float(entropy_reg)
        self.agg = GroupAggregator(feature_groups) if feature_groups is not None else None
        eff_in = len(feature_groups) if feature_groups is not None else int(in_dim)

        self.experts = nn.ModuleList([
            MLPStack(eff_in, _as_int_list(expert_hidden), dropout=dropout, batch_norm=batch_norm)
            for _ in range(int(n_experts))
        ])
        last = self.experts[0].out_dim if int(n_experts) > 0 else eff_in
        self.expert_heads = nn.ModuleList([nn.Linear(last, int(out_dim)) for _ in range(int(n_experts))])

        self.gate = nn.Sequential(
            nn.Linear(eff_in, int(gate_hidden)), nn.ReLU(True),
            nn.Dropout(dropout) if dropout and dropout > 0 else nn.Identity(),
            nn.Linear(int(gate_hidden), int(n_experts))
        )
        self._last_gates = None

    def extra_loss(self):
        if self._last_gates is None or self.entropy_reg <= 0:
            return 0.0
        g = self._last_gates
        ent = -(g * (g.clamp_min(1e-8).log())).sum(dim=1).mean()
        return -self.entropy_reg * ent

    def forward(self, x):
        if self.agg is not None:
            x = self.agg(x)
        g = self.gate(x) / self.temperature
        g = torch.softmax(g, dim=-1)  # (N, E)
        if self.sparse_topk is not None and self.sparse_topk < g.size(1):
            topk = torch.topk(g, k=self.sparse_topk, dim=1)
            mask = torch.zeros_like(g).scatter(1, topk.indices, 1.0)
            g = g * mask
            g = g / (g.sum(dim=1, keepdim=True) + 1e-8)
        self._last_gates = g

        outs = [head(ex(x)) for ex, head in zip(self.experts, self.expert_heads)]
        Y = torch.stack(outs, dim=2)                  # (N, out_dim, E)
        return (Y * g.unsqueeze(1)).sum(dim=2)        # (N, out_dim)

def _try_import_gpytorch():
    try:
        import gpytorch  # noqa
        return importlib.import_module("gpytorch")
    except Exception:
        try:
            _pip_install(["gpytorch"])
            return importlib.import_module("gpytorch")
        except Exception:
            return None

class RandomFourierFeatures(nn.Module):
    def __init__(self, in_dim, out_features=1024, lengthscale=1.0, learnable=False):
        super().__init__()
        self.in_dim = int(in_dim)
        self.m = int(out_features)
        W = torch.randn(self.in_dim, self.m) / float(lengthscale)
        b = 2 * math.pi * torch.rand(self.m)
        self.W = nn.Parameter(W, requires_grad=bool(learnable))
        self.b = nn.Parameter(b, requires_grad=bool(learnable))
        self.scale = math.sqrt(2.0 / self.m)
    def forward(self, x):
        proj = x @ self.W + self.b
        return self.scale * torch.cos(proj)

class RFFRegressor(nn.Module):
    def __init__(self, in_dim, out_dim=1, rff_features=1024, lengthscale=1.0, deep_hidden=(128,), dropout=0.0):
        super().__init__()
        deep_hidden = _as_int_list(deep_hidden)
        self.fe = MLPStack(in_dim, deep_hidden, dropout=dropout, batch_norm=True) if len(deep_hidden) > 0 else nn.Identity()
        fe_dim = self.fe.out_dim if isinstance(self.fe, MLPStack) else int(in_dim)
        self.rff = RandomFourierFeatures(fe_dim, rff_features, lengthscale=lengthscale, learnable=False)
        self.head = nn.Linear(int(rff_features), int(out_dim))
        self.task = "regression"
    def forward(self, x):
        h = self.fe(x)
        z = self.rff(h)
        return self.head(z)

class DKLWrapper(nn.Module):
    """Wrapper so predict() works uniformly; training happens in fit_model()."""
    def __init__(self, extractor, gp_model, likelihood):
        super().__init__()
        self.extractor = extractor
        self.gp_model = gp_model
        self.likelihood = likelihood
        self._is_gp = True

    @torch.no_grad()
    def predict_mean(self, x, device=None):
        dev = _pick_device(device)
        self.extractor.to(dev); self.gp_model.to(dev); self.likelihood.to(dev)
        self.extractor.eval(); self.gp_model.eval(); self.likelihood.eval()
        X = np.asarray(x, dtype=np.float32)
        if hasattr(self, "_expected_input_dim") and X.shape[1] != self._expected_input_dim:
            raise ValueError(f"Expected {self._expected_input_dim} features, got {X.shape[1]}.")
        x = torch.as_tensor(X, dtype=torch.float32, device=dev)
        import gpytorch
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            f = self.gp_model(self.extractor(x))
            pred = self.likelihood(f).mean
        return pred.detach().cpu().numpy()

    @torch.no_grad()
    def predict_mean_and_var(self, x, device=None):
        dev = _pick_device(device)
        self.extractor.to(dev); self.gp_model.to(dev); self.likelihood.to(dev)
        self.extractor.eval(); self.gp_model.eval(); self.likelihood.eval()
        X = np.asarray(x, dtype=np.float32)
        if hasattr(self, "_expected_input_dim") and X.shape[1] != self._expected_input_dim:
            raise ValueError(f"Expected {self._expected_input_dim} features, got {X.shape[1]}.")
        x = torch.as_tensor(X, dtype=torch.float32, device=dev)
        import gpytorch
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            f = self.gp_model(self.extractor(x))
            post = self.likelihood(f)
            mean = post.mean
            var  = post.variance
        return mean.detach().cpu().numpy(), var.detach().cpu().numpy()

@torch.no_grad()
def predict_with_uncertainty(model, X, device=None, prefer_bayesian=True):
    X = np.asarray(X, dtype=np.float32)
    if X.ndim == 1:
        X = X.reshape(1, -1)

    # NAM/MoE original feature count check (only if aggregator present)
    if isinstance(model, nn.Module) and hasattr(model, "agg") and isinstance(model.agg, GroupAggregator) and not getattr(model.agg, "_is_empty", False):
        exp = int(sum(getattr(model.agg, name).numel() for name in model.agg._group_names))
        if X.shape[1] != exp:
            raise ValueError(f"Feature count {X.shape[1]} does not match NAM/MoE grouping map ({exp}).")

    if hasattr(model, "_expected_input_dim") and X.shape[1] != model._expected_input_dim:
        raise ValueError(f"Expected {model._expected_input_dim} features, got {X.shape[1]}.")

    # Optional: catch FT/SAINT grouping mismatches early
    if isinstance(model, nn.Module) and hasattr(model, "group_tok") and getattr(model.group_tok, "has_groups", lambda: False)():
        exp_ft = int(model.group_tok.group_index.numel())
        if X.shape[1] != exp_ft:
            raise ValueError(f"Feature count {X.shape[1]} does not match grouping map ({exp_ft}).")

    # GP path (Bayesian)
    if prefer_bayesian and hasattr(model, "_is_gp") and getattr(model, "_is_gp"):
        mu, var = model.predict_mean_and_var(X, device=device)
        return {"mean": mu.reshape(-1), "variance": var.reshape(-1), "kind": "bayesian"}

    # Heteroscedastic regression (non-Bayesian)
    dev = _pick_device(device)
    if isinstance(model, nn.Module) and getattr(model, "_hetero", False):
        model.eval(); model, dev = _safe_to_device(model, dev)
        Xt = torch.as_tensor(X, dtype=torch.float32, device=dev)
        out = model(Xt)
        if isinstance(out, torch.Tensor) and out.shape[-1] >= 2:
            mu = out[:, :1].detach().cpu().numpy().reshape(-1)
            log_var = torch.clamp(out[:, 1:], min=-10.0, max=5.0)
            var = torch.exp(log_var).detach().cpu().numpy().reshape(-1)
            return {"mean": mu, "variance": var, "kind": "heteroscedastic"}

    # Default: mean only
    mu = predict(model, X, device=device)
    if mu.ndim == 2 and mu.shape[1] > 1:
        return {"mean": mu, "variance": None, "kind": "none"}
    return {"mean": mu.reshape(-1), "variance": None, "kind": "none"}


# ---------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------
def build_model(
    model_type: str, input_dim: int, neurons_per_layer: Sequence[int],
    dropout=0.0, batch_norm=True, out_dim=1, task="regression",
    dense_layers_cnn: Sequence[int] = (256, 128, 64), kernel_size=3,
    final_attention=False,
    cnn_use_max_pool: bool = False, cnn_pool_kernel: int = 2, cnn_pool_stride: int = 2, cnn_pool_padding: int = 0,
    norm_type: str = None,  use_se: bool = True, se_reduction: int = 16, separable: bool = True,   
    dilations: Optional[Sequence[int]] = None, use_positional: bool = True,  pool_type: str = "conv", use_global_pool: bool = True,      
    # FT/SAINT
    ft_d_model=192, ft_heads=8, ft_layers=4, ft_ff_mult=4, ft_dropout=0.1, ft_token_dropout=0.1, ft_use_cls=True, saint_use_cls=True,
    saint_d_model=128, saint_heads=8, saint_layers=4, saint_ff_mult=4, saint_dropout=0.1, saint_token_dropout=0.1,
    ft_group_index: Optional[torch.LongTensor] = None, ft_group_counts: Optional[torch.Tensor] = None,
    ft_scalar_tokenizer: bool = True,
    # TabNet
    tabnet_steps=5, tabnet_feature_dim=64, tabnet_output_dim=64, tabnet_gamma=1.5, tabnet_lambda_sparse=1e-4,
    # NODE
    node_trees=8, node_depth=3,
    # DeepFM
    deepfm_k=16, deepfm_hidden=(128, 64),
    # DCN
    dcn_layers=3, dcn_hidden=(256, 128),
    # NAM
    nam_hidden=(32, 16), nam_activation="relu", nam_add_linear=True, nam_l1=1e-4,
    # MoE
    moe_n_experts=4, moe_expert_hidden=(128, 64), moe_gate_hidden=128, moe_temperature=1.0,
    moe_sparse_topk=None, moe_entropy_reg=0.0,
    # Grouping for NAM/MoE
    feature_groups: Optional[List[np.ndarray]] = None,
):
    input_dim = int(input_dim)
    neurons_per_layer = _as_int_list(neurons_per_layer); kernel_size = int(kernel_size)
    dense_layers_cnn = _as_int_list(dense_layers_cnn)

    if model_type == "resnet":
        return ResNetTabular(input_dim, neurons_per_layer, dropout, batch_norm, int(out_dim), task)
    if model_type == "mlp":
        return MLPNet(input_dim, neurons_per_layer, dropout, batch_norm, int(out_dim), task, final_attention=final_attention)
    if model_type == "mlp_with_attention":
        return MLPWithLayerAttention(input_dim, neurons_per_layer, dropout, batch_norm, int(out_dim), task)
    if model_type == "cnn":
        return CNN1DNet(in_len=input_dim, conv_filters=neurons_per_layer, kernel_size=kernel_size,
                        dropout=dropout, batch_norm=batch_norm, dense_layers=dense_layers_cnn, out_dim=int(out_dim), task=task,
                        use_max_pool=cnn_use_max_pool, pool_kernel=int(cnn_pool_kernel), 
                        pool_stride=int(cnn_pool_stride), pool_padding=int(cnn_pool_padding),
                        norm_type = norm_type, use_se= use_se, se_reduction = int(se_reduction),
                        separable = separable, dilations = dilations, use_positional = use_positional,
                        pool_type = pool_type, use_global_pool = use_global_pool)
    if model_type == "ft_transformer":
        return FTTransformer(n_features=input_dim, out_dim=int(out_dim), task=task,
                             d_model=int(ft_d_model), n_heads=int(ft_heads), n_layers=int(ft_layers),
                             ff_mult=int(ft_ff_mult), dropout=float(ft_dropout),
                             token_dropout=float(ft_token_dropout), use_cls=bool(ft_use_cls),
                             group_index=ft_group_index, group_counts=ft_group_counts,
                             use_scalar_tokenizer=bool(ft_scalar_tokenizer))
    if model_type == "saint":
        return SAINT(n_features=input_dim, out_dim=int(out_dim), task=task,
                     d_model=int(saint_d_model), n_heads=int(saint_heads), n_layers=int(saint_layers),
                     ff_mult=int(saint_ff_mult), dropout=float(saint_dropout), token_dropout=float(saint_token_dropout),
                     group_index=ft_group_index, group_counts=ft_group_counts,
                     use_scalar_tokenizer=bool(ft_scalar_tokenizer), use_cls=bool(saint_use_cls))
    if model_type == "tabnet":
        return TabNet(in_dim=input_dim, out_dim=int(out_dim), task=task,
                      steps=int(tabnet_steps), feature_dim=int(tabnet_feature_dim), output_dim=int(tabnet_output_dim),
                      gamma=float(tabnet_gamma), lambda_sparse=float(tabnet_lambda_sparse))
    if model_type == "node":
        return NODELite(in_dim=input_dim, out_dim=int(out_dim), task=task, n_trees=int(node_trees), depth=int(node_depth))
    if model_type == "deepfm":
        return DeepFM(in_dim=input_dim, out_dim=int(out_dim), task=task,
                  k=int(deepfm_k), deep_hidden=_as_int_list(deepfm_hidden),
                  dropout=float(dropout), batch_norm=bool(batch_norm))
    if model_type == "dcnv2":
        return DCN(in_dim=input_dim, out_dim=int(out_dim), task=task, n_layers=int(dcn_layers),
                   deep_hidden=_as_int_list(dcn_hidden), dropout=float(dropout), batch_norm=bool(batch_norm))
    if model_type == "nam":
        return NAM(in_features=input_dim, out_dim=int(out_dim), task=task,
                   subnet_hidden=_as_int_list(nam_hidden), activation=str(nam_activation),
                   dropout=float(dropout), add_linear=bool(nam_add_linear), l1_lambda=float(nam_l1),
                   feature_groups=feature_groups)
    if model_type == "moe":
        return MoE(in_dim=input_dim, out_dim=int(out_dim), task=task,
                   n_experts=int(moe_n_experts), expert_hidden=_as_int_list(moe_expert_hidden),
                   gate_hidden=int(moe_gate_hidden), dropout=float(dropout), temperature=float(moe_temperature),
                   sparse_topk=(None if moe_sparse_topk in (None, 0) else int(moe_sparse_topk)),
                   entropy_reg=float(moe_entropy_reg), feature_groups=feature_groups,
                   batch_norm=bool(batch_norm))
    if model_type == "gp_dkl":
        raise ValueError("gp_dkl is created in fit_model(); do not call build_model() for gp_dkl.")
    raise ValueError(f"Unknown model_type: {model_type}")

# ---------------------------------------------------------------------
# Training / Prediction API
# ---------------------------------------------------------------------
@torch.no_grad()
def predict(model, X, device=None):
    dev = _pick_device(device)
    was_training = getattr(model, "training", False)

    X = np.asarray(X, dtype=np.float32)
    if X.ndim == 1:
        X = X.reshape(1, -1)

    if isinstance(model, nn.Module) and hasattr(model, "agg") and isinstance(model.agg, GroupAggregator) and not getattr(model.agg, "_is_empty", False):
        exp = int(sum(getattr(model.agg, name).numel() for name in model.agg._group_names))
        if X.shape[1] != exp:
            raise ValueError(f"Feature count {X.shape[1]} does not match NAM/MoE grouping map ({exp}).")

    try:
        if hasattr(model, "predict_mean") and callable(getattr(model, "predict_mean")):
            return model.predict_mean(X, device=dev)

        if isinstance(model, nn.Module):
            model.eval()
            model, dev = _safe_to_device(model, dev)

        if hasattr(model, "_expected_input_dim") and X.shape[1] != model._expected_input_dim:
            raise ValueError(f"Expected {model._expected_input_dim} features, got {X.shape[1]}.")
        if hasattr(model, "group_tok") and getattr(model.group_tok, "has_groups", lambda: False)():
            exp_ft = int(model.group_tok.group_index.numel())
            if X.shape[1] != exp_ft:
                raise ValueError(f"Feature count {X.shape[1]} does not match grouping map ({exp_ft}).")

        X_t = torch.as_tensor(X, dtype=torch.float32, device=dev)
        logits = model(X_t)

        if getattr(model, "_hetero", False) and isinstance(logits, torch.Tensor) and logits.shape[-1] >= 2:
            logits = logits[:, :1]

        return logits.detach().cpu().numpy()

    except Exception as e:
        if os.environ.get("DL_VERBOSE_ERRORS") == "1":
            print(f"[WARN] predict() failed on {dev}: {type(e).__name__}: {e}.", flush=True)

        if getattr(dev, "type", "") in ("cuda", "mps"):
            cpu = torch.device("cpu")
            if isinstance(model, nn.Module):
                model, _ = _safe_to_device(model, cpu)
            X_t = torch.as_tensor(X, dtype=torch.float32, device=cpu)
            logits = model(X_t)
            if getattr(model, "_hetero", False) and isinstance(logits, torch.Tensor) and logits.shape[-1] >= 2:
                logits = logits[:, :1]
            return logits.detach().cpu().numpy()
        raise e
    finally:
        if isinstance(model, nn.Module) and was_training:
            model.train()

def fit_model(
    X, y,
    model_type: str = "mlp_with_attention",
    num_hidden_layers: Optional[int] = None,
    neurons_per_layer: Optional[Sequence[int]] = None,
    learning_rate: float = 1e-3,
    epochs: int = 32,
    batch_size: int = 64,
    l2_weight_decay: float = 1e-3,
    dropout: float = 0.5,
    optimizer_name: str = "adam",
    final_attention: bool = True,
    attention_across_multiple_layers: bool = False,
    batch_norm: bool = True,
    validation_split: float = 0.2,
    compile_model: bool = True,
    device: Optional[str] = None,
    # CNN
    kernel_size: int = 3,
    dense_layers_cnn: Sequence[int] = (256, 128, 64),
    cnn_use_max_pool: bool = False, cnn_pool_kernel: int = 2, cnn_pool_stride: int = 2, cnn_pool_padding: int = 0,
    norm_type: str = None,  use_se: bool = True, se_reduction: int = 16, separable: bool = True,   
    dilations: Optional[Sequence[int]] = None, use_positional: bool = True,  pool_type: str = "conv", use_global_pool: bool = True,
    # Repro / AMP / Grad clip
    deterministic: bool = False, random_seed: Optional[int] = None, use_amp: bool = True, max_grad_norm: Optional[float] = 1.0,
    # Heteroscedastic regression
    heteroscedastic: bool = False,
    # NAM
    nam_hidden: Sequence[int] = (32, 16), nam_activation: str = "relu",
    nam_add_linear: bool = True, nam_l1: float = 1e-4,
    # MoE
    moe_n_experts: int = 4, moe_expert_hidden: Sequence[int] = (128, 64),
    moe_gate_hidden: int = 128, moe_temperature: float = 1.0,
    moe_sparse_topk: Optional[int] = None, moe_entropy_reg: float = 0.0,
    # GP-DKL (+ RFF fallback)
    gp_use_variational: bool = True, gp_num_inducing: int = 512, gp_feature_dim: int = 64,
    gp_kernel: str = "rbf", gp_ard: bool = True, gp_lr_mult: float = 0.5,
    rff_features: int = 1024, rff_lengthscale: float = 1.0, rff_deep_hidden: Sequence[int] = (128,),
    # FT / SAINT hyperparams
    ft_d_model: int = 192, ft_heads: int = 8, ft_layers: int = 4, ft_ff_mult: int = 4, ft_dropout: float = 0.1,
    ft_token_dropout: float = 0.1, ft_use_cls: bool = True, saint_use_cls: bool = True,
    saint_d_model: int = 128, saint_heads: int = 8, saint_layers: int = 4, saint_ff_mult: int = 4,
    saint_dropout: float = 0.1, saint_token_dropout: float = 0.1,
    # Grouping controls for FT/SAINT
    use_grouping: bool = True, group_trigger: int = 2048, group_method: str = "auto",
    init_group_size: int = 64, max_tokens: int = 1024, kmeans_batch: int = 4096, kmeans_iter: int = 100,
    # NODE
    node_trees: int = 8, node_depth: int = 3,
    # TabNet
    tabnet_steps: int = 5, tabnet_feature_dim: int = 64, tabnet_output_dim: int = 64, tabnet_gamma: float = 1.5, tabnet_lambda_sparse: float = 1e-4,
    # DeepFM / DCN
    deepfm_k: int = 16, deepfm_hidden: Sequence[int] = (128, 64),
    dcn_layers: int = 3, dcn_hidden: Sequence[int] = (256, 128),
    # Class weighting (optional)
    auto_class_weights: bool = False,
):
    dev = _pick_device(device)

    if deterministic:
        _set_deterministic(42 if random_seed is None else int(random_seed))
    else:
        if hasattr(torch.backends, "cudnn"):
            torch.backends.cudnn.benchmark = (dev.type == "cuda")

    X = np.asarray(X, dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0)
    y = np.asarray(y)

    # BN feasibility decision (must be here, X is known)
    planned_train = X.shape[0] - int(float(X.shape[0]) * float(validation_split))
    if batch_norm and planned_train < 2:
        batch_norm = False

    input_dim = int(X.shape[1])

    task, n_classes = _infer_task(y)
    out_dim = int(n_classes if task == "multiclass" else 1)

    if task == "multiclass":
        classes = np.unique(y)
        remap = {c: i for i, c in enumerate(classes)}
        y = np.vectorize(remap.get)(y)

    hetero = bool(heteroscedastic and task == "regression")
    if hetero:
        out_dim = 2

    prebuilt_model: Optional[nn.Module] = None

    # ---------- Special case: GP-DKL ----------
    if model_type == "gp_dkl":
        gpytorch = _try_import_gpytorch()
        if gpytorch is None or not gp_use_variational or task in ("binary", "multiclass"):
            # Fallback to RFF (regression only)
            prebuilt_model = RFFRegressor(
                input_dim, out_dim=(out_dim if task == "multiclass" else 1),
                rff_features=int(rff_features), lengthscale=float(rff_lengthscale),
                deep_hidden=_as_int_list(rff_deep_hidden), dropout=float(dropout)
            )
            prebuilt_model, dev = _safe_to_device(prebuilt_model, dev)
        else:
            # Build feature extractor
            extractor = MLPStack(input_dim, [int(gp_feature_dim)], dropout=0.0, batch_norm=False) if gp_feature_dim > 0 else nn.Identity()
            extractor, dev = _safe_to_device(extractor, dev)

            # Inducing points
            n_sub = min(int(gp_num_inducing), X.shape[0])
            sel = np.random.default_rng(0 if random_seed is None else int(random_seed)).choice(X.shape[0], size=n_sub, replace=False)
            with torch.no_grad():
                Z = torch.as_tensor(X[sel], dtype=torch.float32, device=dev)
                Zf = extractor(Z)

            class GPModel(gpytorch.models.ApproximateGP):
                def __init__(self, inducing):
                    var_dist = gpytorch.variational.CholeskyVariationalDistribution(inducing.size(0))
                    var_strat = gpytorch.variational.VariationalStrategy(self, inducing, var_dist, learn_inducing_locations=True)
                    super().__init__(var_strat)
                    if str(gp_kernel).lower() == "matern":
                        base = gpytorch.kernels.MaternKernel(nu=2.5, ard_num_dims=inducing.size(1) if gp_ard else None)
                    else:
                        base = gpytorch.kernels.RBFKernel(ard_num_dims=inducing.size(1) if gp_ard else None)
                    self.covar_module = gpytorch.kernels.ScaleKernel(base)
                    self.mean_module  = gpytorch.means.ConstantMean()
                def forward(self, z):
                    mean = self.mean_module(z)
                    cov  = self.covar_module(z)
                    return gpytorch.distributions.MultivariateNormal(mean, cov)

            likelihood = gpytorch.likelihoods.GaussianLikelihood().to(dev)
            gp_model = GPModel(Zf).to(dev)
            wrapper = DKLWrapper(extractor, gp_model, likelihood)
            setattr(wrapper, "_expected_input_dim", input_dim)
            setattr(wrapper, "_hetero", False)

            # Split data
            X_t = torch.as_tensor(X, dtype=torch.float32)
            y_t = torch.as_tensor(y.reshape(-1, 1), dtype=torch.float32) if task != "multiclass" \
                  else torch.as_tensor(y.reshape(-1), dtype=torch.long)
            dataset = TensorDataset(X_t, y_t)
            n_total = len(dataset)
            raw_val = int(float(n_total) * float(validation_split))
            if n_total <= 1 or validation_split <= 0.0:
                n_train, n_val = n_total, 0
            else:
                n_val = max(1, min(n_total - 1, raw_val))
                n_train = n_total - n_val
            gen_seed = 0 if random_seed is None else int(random_seed)
            gen = torch.Generator().manual_seed(gen_seed)
            train_ds, val_ds = random_split(dataset, [n_train, n_val], generator=gen)

            pin = (dev.type == "cuda")
            # --- NEW: choose a safe batch size + drop_last for BN ---
            n_train = len(train_ds)
            any_bn = any(isinstance(m, (nn.BatchNorm1d, nn.BatchNorm2d, nn.BatchNorm3d))
                         for m in extractor.modules())
            safe_batch = max(2, min(int(batch_size), n_train)) if n_train > 1 else 1
            drop_last_train = (any_bn and (n_train % safe_batch == 1))
            val_batch = max(1, min(int(batch_size), len(val_ds)))
            freeze_bn = any_bn and (n_train < 2)    
                     
            train_dl = DataLoader(train_ds, batch_size=safe_batch, shuffle=True, pin_memory=pin, persistent_workers=False, drop_last=drop_last_train)
            val_dl   = DataLoader(val_ds, batch_size=val_batch, shuffle=False, pin_memory=pin, persistent_workers=False, drop_last=False)

            # Optimizer / loss (ELBO)
            params = [
                {"params": extractor.parameters(), "lr": float(learning_rate)},
                {"params": gp_model.parameters(),  "lr": float(learning_rate) * float(gp_lr_mult)},
                {"params": likelihood.parameters(),"lr": float(learning_rate) * float(gp_lr_mult)},
            ]
            opt = torch.optim.Adam(params, weight_decay=float(l2_weight_decay))
            mll = gpytorch.mlls.VariationalELBO(likelihood, gp_model, num_data=n_train)

            history = {"train_loss": [], "val_loss": []}
            best_val = float("inf"); best_state = None; no_improve = 0; patience = 50

            gp_model.train(); likelihood.train()

            if freeze_bn:
                for m in extractor.modules():
                    if isinstance(m, (nn.BatchNorm1d, nn.BatchNorm2d, nn.BatchNorm3d)):
                        m.eval()

            for _ in range(int(epochs)):               
                               
                # Train
                total = 0.0
                for xb, yb in train_dl:
                    if any_bn and xb.size(0) < 2:
                        continue
                    xb = xb.to(dev, non_blocking=pin); yb = yb.to(dev, non_blocking=pin)
                    opt.zero_grad(set_to_none=True)
                    z = extractor(xb)
                    out = gp_model(z)
                    loss = -mll(out, yb.squeeze(-1))
                    loss.backward()
                    opt.step()
                    total += float(loss.detach().item()) * xb.size(0)
                train_loss = total / max(1, len(train_ds))

                # Validate: MSE of posterior mean
                gp_model.eval(); likelihood.eval()
                total = 0.0; cnt = 0
                with torch.no_grad(), gpytorch.settings.fast_pred_var():
                    for xb, yb in val_dl:
                        xb = xb.to(dev, non_blocking=pin); yb = yb.to(dev, non_blocking=pin)
                        post = likelihood(gp_model(extractor(xb)))
                        mu = post.mean
                        mse = F.mse_loss(mu, yb.squeeze(-1))
                        total += float(mse.detach().item()) * xb.size(0); cnt += xb.size(0)
                val_loss = total / max(1, cnt)
                history["train_loss"].append(train_loss); history["val_loss"].append(val_loss)

                if val_loss < best_val - 1e-7:
                    best_val = val_loss; no_improve = 0
                    best_state = {
                        "extractor": {k: v.detach().cpu().clone() for k, v in extractor.state_dict().items()},
                        "gp": {k: v.detach().cpu().clone() for k, v in gp_model.state_dict().items()},
                        "lik": {k: v.detach().cpu().clone() for k, v in likelihood.state_dict().items()},
                    }
                else:
                    no_improve += 1
                    if no_improve >= patience:
                        break
                gp_model.train(); likelihood.train()

            if best_state is not None:
                extractor.load_state_dict(best_state["extractor"])
                gp_model.load_state_dict(best_state["gp"])
                likelihood.load_state_dict(best_state["lik"])

            return wrapper, history
        # else: prebuilt_model is set for RFF fallback

    # ---------- Grouping (before model build so meta can be attached) ----------
    ft_group_index = None; ft_group_counts = None; ft_scalar_tokenizer = True
    feature_groups = None

    if use_grouping and model_type in ("ft_transformer", "saint") and input_dim > int(group_trigger):
        built = None
        if group_method in ("auto", "kmeans"):
            built = _build_groups_kmeans(
                X_train=X, n_features=input_dim, init_group_size=int(init_group_size),
                max_tokens=int(max_tokens), kmeans_batch=int(kmeans_batch), kmeans_iter=int(kmeans_iter),
                seed=(0 if random_seed is None else int(random_seed))
            )
        if built is None:
            labels, counts = _build_groups_from_chunks(input_dim, int(init_group_size), int(max_tokens))
        else:
            labels, counts = built
        ft_group_index = torch.from_numpy(np.asarray(labels, dtype=np.int64))
        ft_group_counts = torch.from_numpy(np.asarray(counts, dtype=np.int64))

    if model_type in ("nam", "moe") and use_grouping and input_dim > int(group_trigger):
        method = "kmeans" if group_method in ("auto", "kmeans") else "contiguous"
        feature_groups = _build_feature_groups(
            X_train=X, n_features=input_dim,
            group_trigger=int(group_trigger), init_group_size=int(init_group_size),
            max_tokens=int(max_tokens), kmeans_batch=int(kmeans_batch), kmeans_iter=int(kmeans_iter),
            seed=(0 if random_seed is None else int(random_seed)), method=method
        )

    # ---------- Build (or reuse prebuilt) ----------
    if prebuilt_model is None:
        if neurons_per_layer is None:
            neurons_per_layer = [256, 128] if model_type not in ("cnn", "node") else ([64, 64, 64] if model_type == "cnn" else [128])
        neurons_per_layer = _as_int_list(neurons_per_layer)
        if num_hidden_layers is not None and int(num_hidden_layers) != len(neurons_per_layer):
            raise ValueError("Mismatch between num_hidden_layers and length(neurons_per_layer).")

        if model_type == "mlp_with_attention":
            if attention_across_multiple_layers:
                pass
            else:
                model_type = "mlp"
                final_attention = True               

        model = build_model(
            model_type=model_type, input_dim=input_dim, neurons_per_layer=neurons_per_layer,
            dropout=dropout, batch_norm=batch_norm, out_dim=out_dim, task=task,
            dense_layers_cnn=dense_layers_cnn, kernel_size=kernel_size, final_attention=final_attention,
            cnn_use_max_pool=cnn_use_max_pool, cnn_pool_kernel=cnn_pool_kernel, cnn_pool_stride=cnn_pool_stride, cnn_pool_padding=cnn_pool_padding,
            norm_type = norm_type, use_se= use_se, se_reduction = int(se_reduction),separable = separable, dilations = dilations, use_positional = use_positional,
            pool_type = pool_type, use_global_pool = use_global_pool,
            ft_d_model=ft_d_model, ft_heads=ft_heads, ft_layers=ft_layers, ft_ff_mult=ft_ff_mult, ft_dropout=ft_dropout,
            ft_token_dropout=ft_token_dropout, ft_use_cls=ft_use_cls, saint_use_cls=saint_use_cls,
            saint_d_model=saint_d_model, saint_heads=saint_heads, saint_layers=saint_layers, saint_ff_mult=saint_ff_mult,
            saint_dropout=saint_dropout, saint_token_dropout=saint_token_dropout,
            ft_group_index=ft_group_index, ft_group_counts=ft_group_counts, ft_scalar_tokenizer=ft_scalar_tokenizer,
            tabnet_steps=tabnet_steps, tabnet_feature_dim=tabnet_feature_dim, tabnet_output_dim=tabnet_output_dim,
            tabnet_gamma=tabnet_gamma, tabnet_lambda_sparse=tabnet_lambda_sparse,
            node_trees=node_trees, node_depth=node_depth,
            deepfm_k=deepfm_k, deepfm_hidden=deepfm_hidden,
            dcn_layers=dcn_layers, dcn_hidden=dcn_hidden,
            nam_hidden=nam_hidden, nam_activation=nam_activation, nam_add_linear=nam_add_linear, nam_l1=nam_l1,
            moe_n_experts=moe_n_experts, moe_expert_hidden=moe_expert_hidden, moe_gate_hidden=moe_gate_hidden,
            moe_temperature=moe_temperature, moe_sparse_topk=moe_sparse_topk, moe_entropy_reg=moe_entropy_reg,
            feature_groups=feature_groups,
        )
    else:
        model = prebuilt_model  # RFF fallback for gp_dkl

    if ft_group_index is not None:
        setattr(model, "_grouping_meta", {
            "type": "ft/saint",
            "group_index_len": int(ft_group_index.numel()),
            "num_tokens": int(ft_group_counts.numel()),
        })
    if feature_groups is not None:
        setattr(model, "_grouping_meta", {
            "type": "nam/moe",
            "num_groups": len(feature_groups),
            "orig_features": input_dim,
        })

    model, dev = _safe_to_device(model, dev)
    setattr(model, "_hetero", hetero)
    setattr(model, "_expected_input_dim", input_dim)

    model = _maybe_compile(model, dev, compile_model=compile_model)

    X_t = torch.as_tensor(X, dtype=torch.float32)
    if task in ("regression", "binary"):
        y_t = torch.as_tensor(y.reshape(-1, 1), dtype=torch.float32)
    else:
        y_t = torch.as_tensor(y.reshape(-1), dtype=torch.long)

    dataset = TensorDataset(X_t, y_t)
    n_total = len(dataset)
    raw_val = int(float(n_total) * float(validation_split))
    if n_total <= 1 or validation_split <= 0.0:
        n_train, n_val = n_total, 0
    else:
        n_val = max(1, min(n_total - 1, raw_val))
        n_train = n_total - n_val
    gen_seed = 0 if random_seed is None else int(random_seed)
    gen = torch.Generator().manual_seed(gen_seed)
    train_ds, val_ds = random_split(dataset, [n_train, n_val], generator=gen)
    n_train = len(train_ds)
    any_bn = any(isinstance(m, (nn.BatchNorm1d, nn.BatchNorm2d, nn.BatchNorm3d))
                 for m in model.modules())
    safe_batch = max(2, min(int(batch_size), n_train)) if n_train > 1 else 1
    drop_last_train = (any_bn and (n_train % safe_batch == 1))
    val_batch = max(1, min(int(batch_size), len(val_ds)))
    pin = (dev.type == "cuda")
    train_dl = DataLoader(train_ds, batch_size=safe_batch, shuffle=True, pin_memory=pin, persistent_workers=False, drop_last=drop_last_train)
    val_dl   = DataLoader(val_ds, batch_size=val_batch, shuffle=False, pin_memory=pin, persistent_workers=False, drop_last=False)

    rebuilt_for_cpu = False

    # Loss
    if task == "regression":
        if hetero:
            def hetero_nll(y_true, out_logits):
                mu = out_logits[:, :1]
                log_var = out_logits[:, 1:]
                log_var = torch.clamp(log_var, min=-10.0, max=5.0)
                inv_var = torch.exp(-log_var)
                return torch.mean(0.5 * (log_var + (y_true - mu) ** 2 * inv_var))
            criterion = hetero_nll
        else:
            criterion = nn.MSELoss()
    elif task == "binary":
        if auto_class_weights:
            y_cpu = y_t.detach().cpu().numpy().reshape(-1)
            pos = float((y_cpu > 0.5).sum()); neg = float(len(y_cpu) - pos)
            pos_weight = torch.tensor((neg + 1e-8) / (pos + 1e-8), device=dev)
            criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        else:
            criterion = nn.BCEWithLogitsLoss()
    else:
        if auto_class_weights:
            y_cpu = y_t.detach().cpu().numpy().reshape(-1).astype(np.int64)
            counts = np.bincount(y_cpu, minlength=out_dim) + 1e-8
            weights = (counts.sum() / counts)
            w = torch.tensor(weights / weights.mean(), device=dev, dtype=torch.float32)
            criterion = nn.CrossEntropyLoss(weight=w)
        else:
            criterion = nn.CrossEntropyLoss()

    # Optimizer
    name = str(optimizer_name).lower()
    if name == "adam":
        opt = torch.optim.AdamW(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
    elif name == "adamax":
        opt = torch.optim.Adamax(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
    elif name == "sgd":
        opt = torch.optim.SGD(model.parameters(), lr=float(learning_rate), momentum=0.9, weight_decay=float(l2_weight_decay))
    elif name == "rmsprop":
        opt = torch.optim.RMSprop(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
    elif name == "adadelta":
        opt = torch.optim.Adadelta(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
    elif name == "nadam":
        opt = torch.optim.NAdam(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
    else:
        raise ValueError("Unsupported optimizer.")

    amp_ok = use_amp and (dev.type == "cuda" and torch.cuda.is_available())
    """
    scaler = torch.cuda.amp.GradScaler("cuda", enabled=amp_ok)
    """
    from torch.amp import autocast, GradScaler
    device_type = "cuda" if (torch.cuda.is_available() and (device != "cpu")) else "cpu"
    scaler = GradScaler(device_type, enabled=amp_ok)

    patience = 50
    best_val = float("inf")
    best_state = None
    no_improve = 0
    history = {"train_loss": [], "val_loss": []}
    use_val = len(val_ds) > 0

    for _ in range(int(epochs)):
        model.train()
        total = 0.0
        for xb, yb in train_dl:
            if any_bn and xb.size(0) < 2:
                continue
            xb = xb.to(dev, non_blocking=pin); yb = yb.to(dev, non_blocking=pin)
            opt.zero_grad(set_to_none=True)
            try:
                with autocast(device_type = device_type, enabled=amp_ok):
                    logits = model(xb)
                    loss = criterion(yb, logits) if hetero else criterion(logits, yb)
                    if hasattr(model, "extra_loss") and callable(model.extra_loss):
                        loss = loss + model.extra_loss()
                scaler.scale(loss).backward()
                if max_grad_norm is not None:
                    scaler.unscale_(opt)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=float(max_grad_norm))
                scaler.step(opt); scaler.update()

            except Exception:
                # Fallback to CPU and rebuild loaders once with pin_memory=False
                model, dev = _safe_to_device(model, torch.device("cpu"))
                amp_ok = False
                pin = False
                """
                scaler = torch.cuda.amp.GradScaler(enabled=False)
                """
                from torch.amp import GradScaler
                scaler = GradScaler("cpu", enabled=False)
                n_train = len(train_ds)
                any_bn = bool(batch_norm)
                safe_batch = max(2, min(int(batch_size), n_train)) if n_train > 1 else 1
                drop_last_train = (any_bn and (n_train % safe_batch == 1))
                val_batch = max(1, min(int(batch_size), len(val_ds)))
                if not rebuilt_for_cpu:
                    train_dl = DataLoader(
                        train_ds, batch_size=safe_batch, shuffle=True,
                        pin_memory=False, persistent_workers=False, drop_last=drop_last_train
                    )
                    val_dl   = DataLoader(
                        val_ds, batch_size=val_batch, shuffle=False,
                        pin_memory=False, persistent_workers=False,drop_last=False
                    )
                    rebuilt_for_cpu = True

                xb = xb.to(dev); yb = yb.to(dev)

                # Recreate optimizer of the same family (and clear grads to be safe)
                if isinstance(opt, torch.optim.SGD):
                    opt = torch.optim.SGD(model.parameters(), lr=float(learning_rate), momentum=0.9, weight_decay=float(l2_weight_decay))
                elif isinstance(opt, torch.optim.AdamW):
                    opt = torch.optim.AdamW(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                elif isinstance(opt, torch.optim.Adamax):
                    opt = torch.optim.Adamax(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                elif isinstance(opt, torch.optim.RMSprop):
                    opt = torch.optim.RMSprop(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                elif isinstance(opt, torch.optim.Adadelta):
                    opt = torch.optim.Adadelta(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                elif isinstance(opt, torch.optim.NAdam):
                    opt = torch.optim.NAdam(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                else:
                    opt = torch.optim.AdamW(model.parameters(), lr=float(learning_rate), weight_decay=float(l2_weight_decay))
                opt.zero_grad(set_to_none=True)

                logits = model(xb)
                loss = criterion(yb, logits) if hetero else criterion(logits, yb)
                if hasattr(model, "extra_loss") and callable(model.extra_loss):
                    loss = loss + model.extra_loss()
                loss.backward(); opt.step()

            total += float(loss.detach().item()) * xb.size(0)

        train_loss = total / max(1, len(train_ds))

        if use_val:
            model.eval()
            total = 0.0
            with torch.no_grad():
                for xb, yb in val_dl:
                    xb = xb.to(dev, non_blocking=pin); yb = yb.to(dev, non_blocking=pin)
                    logits = model(xb)
                    loss = criterion(yb, logits) if hetero else criterion(logits, yb)
                    if hasattr(model, "extra_loss") and callable(model.extra_loss):
                        loss = loss + model.extra_loss()
                    total += float(loss.detach().item()) * xb.size(0)
            val_loss = total / max(1, len(val_ds))
        else:
            val_loss = float(train_loss)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)

        if val_loss < best_val - 1e-7:
            best_val = val_loss
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1
            if use_val and no_improve >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)

    return model, history
