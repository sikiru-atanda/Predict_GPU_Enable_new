from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

os.environ.setdefault("PYTHONNOUSERSITE", "1")


def _run(cmd: Sequence[str]):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, text=True)


def _pip_install(pkgs: Sequence[str], index_url: Optional[str] = None):
    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--quiet", "--no-cache-dir"]
    if index_url:
        cmd += ["--index-url", index_url]
    cmd += list(pkgs)
    return _run(cmd)


def _pip_uninstall(pkgs: Sequence[str]):
    if not pkgs:
        return None
    return _run([sys.executable, "-m", "pip", "uninstall", "-y", "--quiet", *list(pkgs)])


def _conda_binary() -> Optional[str]:
    """Locate the conda executable across OSes.

    Search order: PREDICTPRO_CONDA env var; CONDA_EXE; standard install paths
    (Windows r-miniconda, ~/miniconda3, ~/anaconda3; /opt/conda; /usr/local/conda);
    finally PATH via ``where`` / ``which``. Returns None if conda cannot be found.
    """
    for var in ("PREDICTPRO_CONDA", "CONDA_EXE"):
        p = os.environ.get(var)
        if p and Path(p).exists():
            return p
    home = Path(os.path.expanduser("~"))
    if sys.platform == "win32":
        candidates = [
            Path(sys.prefix) / "Scripts" / "conda.exe",
            home / "AppData" / "Local" / "r-miniconda" / "Scripts" / "conda.exe",
            home / "miniconda3" / "Scripts" / "conda.exe",
            home / "anaconda3" / "Scripts" / "conda.exe",
            Path(r"C:\ProgramData\miniconda3\Scripts\conda.exe"),
            Path(r"C:\ProgramData\Anaconda3\Scripts\conda.exe"),
        ]
        finder = ["where", "conda.exe"]
    else:
        candidates = [
            Path(sys.prefix) / "bin" / "conda",
            home / "miniconda3" / "bin" / "conda",
            home / "anaconda3" / "bin" / "conda",
            Path("/opt/conda/bin/conda"),
            Path("/usr/local/conda/bin/conda"),
            Path("/usr/local/miniconda3/bin/conda"),
        ]
        finder = ["which", "conda"]
    for cand in candidates:
        if cand.exists():
            return str(cand)
    try:
        proc = _run(finder)
        for line in (proc.stdout or "").splitlines():
            ln = line.strip().strip('"')
            if ln and Path(ln).exists():
                return ln
    except Exception:
        pass
    return None


def _conda_env_name() -> Optional[str]:
    """Best-effort detection of the current conda env name.

    Prefers ``CONDA_DEFAULT_ENV`` if set (and not ``base``). Falls back to
    inferring from ``sys.prefix`` when its path looks like ``.../envs/<name>``.
    Returns None when the active Python is not inside a named conda env.
    """
    name = os.environ.get("CONDA_DEFAULT_ENV", "")
    if name and name != "base":
        return name
    parts = Path(sys.prefix).resolve().parts
    if "envs" in parts:
        i = parts.index("envs")
        if i + 1 < len(parts):
            return parts[i + 1]
    return None


def _install_scikit_sparse(quiet: bool = False) -> Optional[str]:
    """Best-effort cross-OS install of scikit-sparse for the MME CHOLMOD backend.

    The PredictProR sparse-MME REML engine in ``mme_reml.py`` auto-selects
    ``sksparse.cholmod`` when available and otherwise falls back to
    ``scipy.sparse.linalg.splu`` (~2x slower on SPD systems but functionally
    equivalent). This helper makes the CHOLMOD path opt-in by attempting the
    install in this order:

    1. If already importable -- return its version (no-op).
    2. Conda-forge into the active named env -- PREFERRED, works on
       Windows / Linux / macOS because conda ships the SuiteSparse C library
       (no separate apt/brew/yum step required).
    3. Pip -- last resort, works on Linux/macOS where SuiteSparse is system-
       installed; on Windows there is no pip wheel so this path fails fast.
    4. Failure -- prints an OS-specific install hint and returns None so the
       caller can record ``mme_factor_backend="splu"``. Never raises.
    """
    def _probe() -> Optional[str]:
        try:
            importlib.import_module("sksparse.cholmod")
        except Exception:
            return None
        try:
            parent = importlib.import_module("sksparse")
            return str(getattr(parent, "__version__", "installed"))
        except Exception:
            return "installed"

    existing = _probe()
    if existing is not None:
        if not quiet:
            print(f"scikit-sparse already installed ({existing}); skipping.",
                  file=sys.stderr, flush=True)
        return existing

    conda_bin = _conda_binary()
    env_name = _conda_env_name()
    if conda_bin and env_name:
        if not quiet:
            print(
                f"Installing scikit-sparse via conda-forge into env "
                f"'{env_name}' (this also pulls SuiteSparse)...",
                file=sys.stderr, flush=True,
            )
        proc = _run([
            conda_bin, "install", "-y", "-n", env_name,
            "-c", "conda-forge", "scikit-sparse",
        ])
        if proc.returncode == 0:
            importlib.invalidate_caches()
            v = _probe()
            if v is not None:
                return v
        if not quiet:
            tail = (proc.stdout or "").splitlines()[-6:]
            print("conda install failed; tail of output:",
                  file=sys.stderr, flush=True)
            for ln in tail:
                print(f"    {ln}", file=sys.stderr, flush=True)

    if not quiet:
        print("Trying pip install scikit-sparse...",
              file=sys.stderr, flush=True)
    proc = _pip_install(["scikit-sparse"])
    if proc.returncode == 0:
        importlib.invalidate_caches()
        v = _probe()
        if v is not None:
            return v
    if not quiet:
        tail = (proc.stdout or "").splitlines()[-6:]
        print("pip install failed; tail of output:",
              file=sys.stderr, flush=True)
        for ln in tail:
            print(f"    {ln}", file=sys.stderr, flush=True)
        plat = sys.platform
        print(
            "\nscikit-sparse is OPTIONAL. The MME REML engine falls back to\n"
            "scipy.sparse.linalg.splu (slower but correct). To install\n"
            "manually:",
            file=sys.stderr, flush=True,
        )
        if plat.startswith("linux"):
            print(
                "  # Debian/Ubuntu:\n"
                "    sudo apt-get install -y libsuitesparse-dev\n"
                "    pip install scikit-sparse\n"
                "  # RHEL/CentOS/Fedora:\n"
                "    sudo yum install -y suitesparse-devel\n"
                "    pip install scikit-sparse\n"
                "  # Or in a conda env (no system package needed):\n"
                "    conda install -n <env> -c conda-forge scikit-sparse",
                file=sys.stderr, flush=True,
            )
        elif plat == "darwin":
            print(
                "  # macOS Homebrew:\n"
                "    brew install suite-sparse\n"
                "    pip install scikit-sparse\n"
                "  # Or in a conda env:\n"
                "    conda install -n <env> -c conda-forge scikit-sparse",
                file=sys.stderr, flush=True,
            )
        elif plat == "win32":
            print(
                "  # Windows -- conda-forge is the only supported path:\n"
                "    conda install -n <env> -c conda-forge scikit-sparse\n"
                "  (no pip wheel available on Windows)",
                file=sys.stderr, flush=True,
            )
    return None


def _has_nvidia_gpu() -> bool:
    try:
        return _run(["nvidia-smi"]).returncode == 0
    except Exception:
        return False


def _torch_import_probe() -> Tuple[bool, str]:
    code = (
        "import json\n"
        "try:\n"
        "    import torch\n"
        "    _ = torch.__version__\n"
        "    _ = torch.tensor([0.0]).sum().item()\n"
        "    print(json.dumps({'ok': True, 'version': getattr(torch, '__version__', None), "
        "'cuda_available': bool(torch.cuda.is_available()), "
        "'cuda_compiled': getattr(getattr(torch, 'version', None), 'cuda', None)}))\n"
        "except Exception as exc:\n"
        "    print(json.dumps({'ok': False, 'error': repr(exc)}))\n"
    )
    try:
        proc = _run([sys.executable, "-c", code])
        payload = json.loads((proc.stdout or "").strip().splitlines()[-1])
        return bool(payload.get("ok")), str(payload.get("error") or "")
    except Exception as exc:
        return False, repr(exc)


def _torch_probe_info() -> Dict[str, Any]:
    code = (
        "import json\n"
        "try:\n"
        "    import torch\n"
        "    print(json.dumps({'ok': True, 'version': getattr(torch, '__version__', None), "
        "'cuda_available': bool(torch.cuda.is_available()), "
        "'cuda_compiled': getattr(getattr(torch, 'version', None), 'cuda', None), "
        "'cuda_device_count': int(torch.cuda.device_count()) if torch.cuda.is_available() else 0}))\n"
        "except Exception as exc:\n"
        "    print(json.dumps({'ok': False, 'error': repr(exc)}))\n"
    )
    try:
        proc = _run([sys.executable, "-c", code])
        return json.loads((proc.stdout or "").strip().splitlines()[-1])
    except Exception as exc:
        return {"ok": False, "error": repr(exc)}


def _torch_cuda_available() -> bool:
    return bool(_torch_probe_info().get("cuda_available"))


def _torch_cuda_compiled() -> Optional[str]:
    value = _torch_probe_info().get("cuda_compiled")
    return None if value in ("", "None") else value


def _install_torch_with_fallback(
    *,
    prefer_gpu: bool,
    cuda: str,
    torch_version: Optional[str],
    forced_index: Optional[str],
) -> Dict[str, str]:
    logs: List[str] = []
    candidates: List[Tuple[str, Optional[str], str]] = []
    if forced_index:
        candidates.append(("forced", forced_index, "gpu" if "cu" in forced_index else "cpu"))
    elif platform.system() == "Darwin":
        candidates.append(("macos-default", None, "cpu"))
    else:
        if bool(prefer_gpu) and cuda != "cpu" and _has_nvidia_gpu():
            cu_tag = cuda.lower() if str(cuda).lower().startswith("cu") else "cu121"
            candidates.append((f"pytorch-{cu_tag}", f"https://download.pytorch.org/whl/{cu_tag}", "gpu"))
        candidates.append(("pytorch-cpu", "https://download.pytorch.org/whl/cpu", "cpu"))

    torch_spec = [f"torch=={torch_version}"] if torch_version else ["torch"]
    last_error = ""
    for desc, index_url, flavor in candidates:
        _pip_uninstall(["torch", "torchvision", "torchaudio", "torchtext", "functorch"])
        result = _pip_install(torch_spec, index_url=index_url)
        logs.append(f"[{desc}] pip stdout:\n{result.stdout or ''}")
        ok, err = _torch_import_probe()
        if ok:
            return {
                "index_used": index_url or "default",
                "flavor": flavor,
                "attempt": desc,
                "logs": "\n".join(logs)[-8000:],
            }
        last_error = err
        if flavor != "cpu" and (("WinError 127" in err) or ("DLL load failed" in err) or ("shm.dll" in err)):
            continue
    raise RuntimeError(f"Unable to import torch after installs. Last error: {last_error}")


def _ensure_torch(prefer_gpu: bool = True,
                  cuda: str = "auto",
                  torch_version: Optional[str] = None,
                  index_url: Optional[str] = None) -> Dict[str, str]:
    forced_index = index_url or (os.environ.get("PYTORCH_INDEX_URL") or "").strip() or None
    ok, _ = _torch_import_probe()
    if ok:
        gpu_requested = bool(prefer_gpu) and str(cuda).lower() != "cpu" and (_has_nvidia_gpu() or forced_index is not None)
        existing_cpu_torch = _torch_cuda_compiled() is None
        if not (gpu_requested and existing_cpu_torch and not _torch_cuda_available()):
            return {
                "index_used": None,
                "flavor": "existing",
                "attempt": "import-existing",
                "logs": "",
            }
    return _install_torch_with_fallback(
        prefer_gpu=bool(prefer_gpu),
        cuda=str(cuda).lower(),
        torch_version=torch_version,
        forced_index=forced_index,
    )


def _requirements_without_torch(project_root: Optional[str]) -> List[str]:
    if not project_root:
        return []
    req_path = os.path.join(project_root, "requirements.txt")
    if not os.path.exists(req_path):
        return []
    out: List[str] = []
    for raw in Path(req_path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        pkg_name = line.split("==", 1)[0].split(">=", 1)[0].split("<=", 1)[0].split("~=", 1)[0].strip().lower()
        if pkg_name == "torch":
            continue
        out.append(line)
    return out


def setup_deps(project_root: Optional[str] = None,
               prefer_gpu: bool = True,
               cuda: str = "auto",
               torch_version: Optional[str] = None,
               index_url: Optional[str] = None,
               install_scikit_sparse: bool = True):
    _pip_install(["numpy", "pandas"])
    torch_meta = _ensure_torch(
        prefer_gpu=bool(prefer_gpu),
        cuda=str(cuda).lower(),
        torch_version=torch_version,
        index_url=index_url,
    )
    reqs = _requirements_without_torch(project_root)
    if not reqs:
        reqs = ["gpytorch", "linear-operator", "zarr", "pytest", "pyreadr"]
    _pip_install(reqs)

    scikit_sparse_version: Optional[str] = None
    if install_scikit_sparse:
        scikit_sparse_version = _install_scikit_sparse(quiet=False)

    out = {"python": sys.executable}
    out.update({
        "torch_index_used": torch_meta.get("index_used"),
        "torch_install_flavor": torch_meta.get("flavor"),
        "torch_install_attempt": torch_meta.get("attempt"),
    })
    for name in ("numpy", "pandas", "torch", "gpytorch", "linear_operator", "zarr"):
        try:
            mod = __import__(name)
            out[f"{name}_version"] = getattr(mod, "__version__", None)
        except Exception:
            out[f"{name}_version"] = None
    out["scikit_sparse_version"] = scikit_sparse_version
    out["mme_factor_backend"] = "cholmod" if scikit_sparse_version else "splu"
    try:
        import torch

        out["cuda_available"] = bool(torch.cuda.is_available())
        out["cuda_compiled"] = getattr(getattr(torch, "version", None), "cuda", None)
        out["cuda_device_count"] = int(torch.cuda.device_count()) if out["cuda_available"] else 0
        out["gpu_name"] = str(torch.cuda.get_device_name(0)) if out["cuda_available"] and out["cuda_device_count"] > 0 else None
    except Exception:
        out["cuda_available"] = False
        out["cuda_device_count"] = 0
        out["gpu_name"] = None
    return out


def _load_framework(project_root: str):
    module_name = "gp_framework"
    py = os.path.join(project_root, module_name + ".py")
    if not os.path.exists(py):
        raise FileNotFoundError(f"Framework file not found: {py}")
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    importlib.invalidate_caches()
    return importlib.import_module(module_name)


def _read_shape(path: str) -> tuple[int, int]:
    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Shape file must contain two lines: {path}")
    return int(lines[0]), int(lines[1])


def _read_matrix_bin(bin_path: str, meta_path: str) -> np.ndarray:
    nrow, ncol = _read_shape(meta_path)
    arr = np.fromfile(bin_path, dtype="<f8")
    if arr.size != nrow * ncol:
        raise ValueError(f"Matrix binary size mismatch for {bin_path}: expected {nrow*ncol}, got {arr.size}")
    return arr.reshape((nrow, ncol), order="F")


def _read_lines(path: str) -> list[str]:
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def _read_idx(path: str) -> np.ndarray:
    vals = _read_lines(path)
    if not vals:
        return np.asarray([], dtype=np.int64)
    return np.asarray([int(v) for v in vals], dtype=np.int64)


def _write_df(df: Any, path: str):
    if df is None:
        return
    if isinstance(df, pd.DataFrame):
        out = df
    else:
        out = pd.DataFrame(df)
    out.to_csv(path, index=False)


def _write_matrix_bin(matrix: Any, bin_path: str, meta_path: str) -> None:
    arr = np.asarray(matrix, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"matrix must be 2D, got shape {arr.shape}")
    np.asarray(arr, dtype="<f8").ravel(order="F").tofile(bin_path)
    Path(meta_path).write_text(f"{arr.shape[0]}\n{arr.shape[1]}\n", encoding="utf-8")


def _json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.ndarray,)):
        return value.tolist()
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _json_sanitize(value: Any):
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if np.isfinite(value) else None
    if isinstance(value, np.ndarray):
        return _json_sanitize(value.tolist())
    if isinstance(value, pd.DataFrame):
        clean = value.replace([np.inf, -np.inf], np.nan)
        clean = clean.where(pd.notnull(clean), None)
        return _json_sanitize(clean.to_dict(orient="records"))
    if isinstance(value, pd.Series):
        clean = value.replace([np.inf, -np.inf], np.nan)
        clean = clean.where(pd.notnull(clean), None)
        return _json_sanitize(clean.tolist())
    if isinstance(value, dict):
        return {str(k): _json_sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_sanitize(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(_json_sanitize(payload), indent=2, default=_json_default, allow_nan=False),
        encoding="utf-8",
    )


def _coalesce(*values: Any, default: Any = None) -> Any:
    for value in values:
        if value is not None:
            return value
    return default


def _spec_get(spec: Dict[str, Any], *paths: str, default: Any = None) -> Any:
    for path in paths:
        cur: Any = spec
        ok = True
        for part in path.split("."):
            if isinstance(cur, dict) and part in cur:
                cur = cur[part]
            else:
                ok = False
                break
        if ok and cur is not None:
            return cur
    return default


def _resolve_path(value: Any, base_dir: Path, *, required: bool = False, label: str = "path") -> Optional[str]:
    if value is None or value == "":
        if required:
            raise ValueError(f"Missing required {label}")
        return None
    path = Path(str(value))
    if not path.is_absolute():
        path = base_dir / path
    return str(path)


def _read_json_file(path: str) -> Dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"JSON file must contain an object: {path}")
    return obj


def _read_matrix_from_spec(obj: Any, base_dir: Path, label: str) -> Optional[np.ndarray]:
    if obj is None:
        return None
    if isinstance(obj, str):
        return np.asarray(pd.read_csv(_resolve_path(obj, base_dir, required=True, label=label), header=None), dtype=float)
    if isinstance(obj, list):
        return np.asarray(obj, dtype=np.float64)
    if not isinstance(obj, dict):
        raise ValueError(f"{label} must be a path, array, or object")
    bin_path = _resolve_path(
        _coalesce(obj.get("bin"), obj.get("bin_path"), obj.get("matrix_bin"), obj.get("path")),
        base_dir,
        required=True,
        label=f"{label}.bin",
    )
    meta_path = _resolve_path(
        _coalesce(obj.get("meta"), obj.get("meta_path"), obj.get("shape"), obj.get("shape_path"), obj.get("matrix_meta")),
        base_dir,
        required=True,
        label=f"{label}.meta",
    )
    return _read_matrix_bin(bin_path, meta_path)


def _read_matrix_mapping_from_spec(obj: Any, base_dir: Path, label: str) -> Optional[Dict[str, np.ndarray]]:
    if obj is None:
        return None
    out: Dict[str, np.ndarray] = {}
    if isinstance(obj, list):
        for i, item in enumerate(obj, start=1):
            if not isinstance(item, dict):
                raise ValueError(f"{label}[{i}] must be an object")
            name = _coalesce(item.get("name"), item.get("key"), item.get("kernel"))
            if name is None:
                raise ValueError(f"{label}[{i}] is missing name/key/kernel")
            matrix_spec = _coalesce(item.get("matrix"), item.get("value"), item)
            out[str(name)] = _read_matrix_from_spec(matrix_spec, base_dir, f"{label}.{name}")
        return out
    if isinstance(obj, dict):
        for name, matrix_spec in obj.items():
            out[str(name)] = _read_matrix_from_spec(matrix_spec, base_dir, f"{label}.{name}")
        return out
    raise ValueError(f"{label} must be an object or list")


def _read_geno_kernels_from_spec(spec: Dict[str, Any], base_dir: Path) -> Dict[str, np.ndarray]:
    kernels_spec = _spec_get(
        spec,
        "geno_kernels",
        "genotype_kernels",
        "inputs.geno_kernels",
        "inputs.genotype_kernels",
        default=None,
    )
    if kernels_spec is not None:
        kernels = _read_matrix_mapping_from_spec(kernels_spec, base_dir, "geno_kernels")
        if not kernels:
            raise ValueError("geno_kernels must contain at least one kernel")
        return kernels

    kernel_spec = _coalesce(
        _spec_get(spec, "geno_kernel", "genotype_kernel", "inputs.geno_kernel", "inputs.genotype_kernel"),
        {
            "bin": _spec_get(spec, "kernel_bin", "geno_kernel_bin", "inputs.kernel_bin", "inputs.geno_kernel_bin"),
            "meta": _spec_get(spec, "kernel_meta", "geno_kernel_meta", "inputs.kernel_meta", "inputs.geno_kernel_meta"),
        },
    )
    geno_kernel = _read_matrix_from_spec(kernel_spec, base_dir, "geno_kernel")
    if geno_kernel is None:
        raise ValueError("geno_kernel is required")
    kernel_name = str(_coalesce(_spec_get(spec, "kernel_name", "geno_kernel_name"), default="G"))
    return {kernel_name: geno_kernel}


def _read_lines_from_spec(value: Any, base_dir: Path, label: str, *, required: bool = False) -> Optional[List[str]]:
    if value is None:
        if required:
            raise ValueError(f"Missing required {label}")
        return None
    if isinstance(value, list):
        return [str(v) for v in value]
    path = _resolve_path(value, base_dir, required=True, label=label)
    return _read_lines(path)


def _read_idx_from_spec(value: Any, base_dir: Path, label: str, *, required: bool = False) -> Optional[np.ndarray]:
    if value is None:
        if required:
            raise ValueError(f"Missing required {label}")
        return None
    if isinstance(value, list):
        return np.asarray([int(v) for v in value], dtype=np.int64)
    path = _resolve_path(value, base_dir, required=True, label=label)
    return _read_idx(path)


def _read_vector_from_spec(
    value: Any,
    base_dir: Path,
    label: str,
    *,
    required: bool = False,
    dtype: Any = np.float64,
) -> Optional[np.ndarray]:
    if value is None:
        if required:
            raise ValueError(f"Missing required {label}")
        return None
    if isinstance(value, list):
        return np.asarray(value, dtype=dtype).reshape(-1)
    path = _resolve_path(value, base_dir, required=True, label=label)
    arr = pd.read_csv(path, header=None).to_numpy().reshape(-1)
    return np.asarray(arr, dtype=dtype).reshape(-1)


def _normalize_factor_cache(cache_spec: Any, base_dir: Path) -> Optional[Dict[str, Any]]:
    if cache_spec is None:
        return None
    if isinstance(cache_spec, str):
        cache_spec = _read_json_file(_resolve_path(cache_spec, base_dir, required=True, label="factor_cache_spec"))
    if not isinstance(cache_spec, dict):
        raise ValueError("factor_cache_spec must be a JSON object or path to one")
    out = dict(cache_spec)
    ctype = str(out.get("type", "zarr")).lower()
    out["type"] = ctype
    if ctype == "zarr":
        root = _coalesce(out.get("root_dir"), out.get("root"), out.get("path"))
        out["root_dir"] = _resolve_path(root, base_dir, required=True, label="factor_cache.root_dir")
        out.pop("root", None)
        out.pop("path", None)
    elif ctype == "memmap":
        paths = out.get("paths")
        if isinstance(paths, dict):
            out["paths"] = {k: _resolve_path(v, base_dir, required=True, label=f"factor_cache.paths.{k}") for k, v in paths.items()}
        elif isinstance(paths, list):
            out["paths"] = [_resolve_path(v, base_dir, required=True, label="factor_cache.paths") for v in paths]
    return out


def _align_env_similarity_if_needed(
    env_similarity: Optional[np.ndarray],
    env_ids: Optional[Sequence[str]],
    pheno_df: pd.DataFrame,
    env_col: str,
) -> Optional[np.ndarray]:
    if env_similarity is None or env_ids is None:
        return env_similarity
    ids = [str(x) for x in env_ids]
    if len(set(ids)) != len(ids):
        raise ValueError("env_ids contains duplicate values")
    K = np.asarray(env_similarity, dtype=np.float64)
    if K.ndim != 2 or K.shape[0] != K.shape[1] or K.shape[0] != len(ids):
        raise ValueError(f"env_similarity shape {K.shape} does not match env_ids length {len(ids)}")
    env_levels = list(pheno_df[env_col].astype("category").cat.categories.astype(str))
    pos = {env_id: i for i, env_id in enumerate(ids)}
    missing = [env for env in env_levels if env not in pos]
    if missing:
        raise ValueError(f"env_similarity/env_ids is missing phenotype environments: {missing[:5]}")
    order = [pos[env] for env in env_levels]
    return K[np.ix_(order, order)]


def _write_fit_outputs(res: Dict[str, Any], out_dir: Path, meta: Dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    result = dict(res.get("result") or {})
    predictions = result.pop("predictions", None)
    if predictions is None:
        predictions = res.get("predictions")
    _write_df(predictions, str(out_dir / "predictions.csv"))
    _write_df(res.get("var_components"), str(out_dir / "var_components.csv"))
    _write_df(res.get("var_components_summary"), str(out_dir / "var_components_summary.csv"))
    for key in (
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
        "fa_kernel_weights",
        "sigma2_resid_env",
        "sigma2_resid_overall",
        "variance_component_method",
        "diagnostics",
    ):
        if key in res and key not in result:
            result[key] = res.get(key)
    if "per_env" not in result and "env_variance_summary" in result:
        result["per_env"] = result.get("env_variance_summary")
    _write_json(out_dir / "result.json", result)
    _write_json(out_dir / "fit.json", res.get("fit") or {})
    _write_json(out_dir / "info.json", res.get("info") or {})
    meta = dict(meta)
    meta["diagnostics_keys"] = sorted(list((res.get("diagnostics") or {}).keys()))
    meta["result_keys"] = sorted(list(result.keys()))
    meta["fit_keys"] = sorted(list((res.get("fit") or {}).keys()))
    meta["info_keys"] = sorted(list((res.get("info") or {}).keys()))
    if isinstance(res.get("_meta"), dict):
        meta["fit_meta"] = res.get("_meta")
    _write_json(out_dir / "meta.json", meta)


def _write_public_result_outputs(res: Dict[str, Any], out_dir: Path, meta: Dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    result = dict(res.get("result") or {})
    for key in (
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
        "fa_kernel_weights",
        "sigma2_resid_env",
        "sigma2_resid_overall",
        "diagnostics",
    ):
        if key in res and key not in result:
            result[key] = res.get(key)
    if "per_env" not in result and "env_variance_summary" in result:
        result["per_env"] = result.get("env_variance_summary")
    predictions = result.pop("predictions", None)
    if predictions is None:
        predictions = res.get("predictions")
    _write_df(predictions, str(out_dir / "predictions.csv"))
    _write_json(out_dir / "result.json", result)
    _write_json(out_dir / "fit.json", res.get("fit") or {})
    _write_json(out_dir / "info.json", res.get("info") or {})
    meta = dict(meta)
    meta["result_keys"] = sorted(list(result.keys()))
    meta["fit_keys"] = sorted(list((res.get("fit") or {}).keys()))
    meta["info_keys"] = sorted(list((res.get("info") or {}).keys()))
    _write_json(out_dir / "meta.json", meta)


def _fit_predict(args) -> int:
    fw = _load_framework(args.project_root)
    return _fit_predict_with_framework(args, fw)


def _fit_predict_spec(args) -> int:
    spec_path = Path(args.spec).resolve()
    spec = _read_json_file(str(spec_path))
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root"), default=None)
    if project_root is None:
        raise ValueError("project_root must be supplied in the spec or with --project-root")
    fw = _load_framework(_resolve_path(project_root, base_dir, required=True, label="project_root"))

    pheno_csv = _resolve_path(
        _spec_get(spec, "pheno_csv", "phenotype_csv", "inputs.pheno_csv", "inputs.phenotype_csv"),
        base_dir,
        required=True,
        label="pheno_csv",
    )
    pheno_df = pd.read_csv(pheno_csv)

    geno_kernels = _read_geno_kernels_from_spec(spec, base_dir)
    geno_ids = _read_lines_from_spec(
        _spec_get(spec, "geno_ids", "genotype_ids", "inputs.geno_ids", "inputs.genotype_ids"),
        base_dir,
        "geno_ids",
        required=True,
    )
    train_idx = _read_idx_from_spec(
        _spec_get(spec, "train_idx", "train_index", "inputs.train_idx", "inputs.train_index"),
        base_dir,
        "train_idx",
        required=True,
    )
    test_idx = _read_idx_from_spec(
        _spec_get(spec, "test_idx", "test_index", "inputs.test_idx", "inputs.test_index"),
        base_dir,
        "test_idx",
        required=True,
    )

    cols = _spec_get(spec, "columns", default={}) or {}
    gid_col = _coalesce(_spec_get(spec, "gid_col"), cols.get("gid"), cols.get("gid_col"), default="GID")
    env_col = _coalesce(_spec_get(spec, "env_col"), cols.get("env"), cols.get("env_col"), default="Env")
    y_col = _coalesce(_spec_get(spec, "y_col"), cols.get("y"), cols.get("y_col"), default="Trait")

    fit_options = dict(_spec_get(spec, "fit", "options", "model", default={}) or {})
    output_level = str(_coalesce(_spec_get(spec, "output_level"), fit_options.pop("output_level", None), default="predict_only"))
    method = str(_coalesce(_spec_get(spec, "method"), fit_options.pop("method", None), default="krr_exact"))
    backend = str(_coalesce(_spec_get(spec, "backend"), fit_options.pop("backend", None), default="auto"))
    prediction_output = str(_coalesce(_spec_get(spec, "prediction_output"), fit_options.pop("prediction_output", None), default="all"))
    varcomp_mode = str(_coalesce(_spec_get(spec, "varcomp_mode"), fit_options.pop("varcomp_mode", None), default="mom"))
    seed = int(_coalesce(_spec_get(spec, "seed"), fit_options.pop("seed", None), default=12345))
    fa_rank = int(_coalesce(_spec_get(spec, "fa_rank"), fit_options.pop("fa_rank", None), default=1))

    env_covariates = None
    env_covariates_csv = _resolve_path(
        _spec_get(spec, "env_covariates_csv", "inputs.env_covariates_csv", "environment_covariates_csv", "env_covariates.csv", "environment_covariates.csv"),
        base_dir,
        label="env_covariates_csv",
    )
    if env_covariates_csv is not None:
        env_covariates = pd.read_csv(env_covariates_csv)

    env_similarity_spec = _coalesce(
        _spec_get(spec, "env_similarity", "environment_similarity", "inputs.env_similarity"),
        {
            "bin": _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin"),
            "meta": _spec_get(spec, "env_similarity_meta", "inputs.env_similarity_meta"),
        }
        if _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin") is not None
        else None,
    )
    env_similarity = _read_matrix_from_spec(env_similarity_spec, base_dir, "env_similarity")
    env_ids = _read_lines_from_spec(
        _spec_get(spec, "env_ids", "environment_ids", "inputs.env_ids", "inputs.environment_ids", "env_similarity.env_ids", "env_similarity.ids", "environment_similarity.env_ids"),
        base_dir,
        "env_ids",
        required=False,
    )
    env_similarity = _align_env_similarity_if_needed(env_similarity, env_ids, pheno_df, str(env_col))

    factor_cache_spec = _coalesce(
        _spec_get(spec, "factor_cache_spec", "grm_factor_cache", "gp_factor_cache", "inputs.factor_cache_spec"),
        _spec_get(spec, "grm_factor_cache_root", "gp_factor_cache_root"),
    )
    if isinstance(factor_cache_spec, str) and not factor_cache_spec.lower().endswith(".json"):
        factor_cache_spec = {"type": "zarr", "root_dir": factor_cache_spec}
    grm_factor_cache = _normalize_factor_cache(factor_cache_spec, base_dir)

    include_components = _coalesce(
        _spec_get(spec, "include_components"),
        fit_options.pop("include_components", None),
        default=(["g", "ge", "e"] if (env_similarity is not None or env_covariates is not None) else ["g"]),
    )
    fixed_effects = _coalesce(_spec_get(spec, "fixed_effects"), fit_options.pop("fixed_effects", None), default=[])

    kw: Dict[str, Any] = {
        "include_components": list(include_components),
        "fixed_effects": list(fixed_effects),
        "method": method,
        "backend": backend,
        "prediction_output": prediction_output,
        "output_level": output_level,
        "point_predictions_only": output_level == "predict_only",
        "return_se": output_level in ("predict_with_se", "full_vc"),
        "compute_ai_se": output_level == "full_vc",
        "standardize": str(_coalesce(_spec_get(spec, "standardize"), fit_options.pop("standardize", None), default="global")),
        "varcomp_mode": varcomp_mode,
        "dtype": str(_coalesce(_spec_get(spec, "dtype"), fit_options.pop("dtype", None), default="float64")),
        "seed": seed,
        "train_idx": train_idx,
        "test_idx": test_idx,
    }
    if method == "gp_icm_fa":
        kw["env_structure"] = str(_coalesce(_spec_get(spec, "env_structure"), fit_options.pop("env_structure", None), default="fa"))
        kw["fa_rank"] = fa_rank
    if env_similarity is not None:
        kw["env_similarity"] = env_similarity
    if env_covariates is not None:
        kw["env_covariates"] = env_covariates
    if grm_factor_cache is not None:
        kw["grm_factor_cache"] = grm_factor_cache
    kw.update(fit_options)

    res = fw.fit_mixed_model(
        pheno_df=pheno_df,
        gid_col=str(gid_col),
        env_col=str(env_col),
        y_col=str(y_col),
        geno_kernels=geno_kernels,
        geno_ids=list(geno_ids or []),
        **kw,
    )

    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))
    _write_fit_outputs(
        res,
        out_dir,
        {
            "python": sys.executable,
            "command": "fit-predict-spec",
            "spec": str(spec_path),
            "method": method,
            "backend": backend,
            "output_level": output_level,
            "prediction_output": prediction_output,
            "varcomp_mode": varcomp_mode,
            "kernel_names": list(geno_kernels.keys()),
            "has_env_covariates": env_covariates is not None,
            "has_env_similarity": env_similarity is not None,
            "has_grm_factor_cache": grm_factor_cache is not None,
        },
    )
    return 0


def _read_cv_idx_from_spec(value: Any, base_dir: Path, label: str) -> np.ndarray:
    if isinstance(value, (int, np.integer)):
        return np.asarray([int(value)], dtype=np.int64)
    return _read_idx_from_spec(value, base_dir, label, required=True)


def _validate_cv_fold_indices(
    *,
    fold_id: str,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    n_rows: int,
    y_values: np.ndarray,
) -> None:
    train_idx = np.asarray(train_idx, dtype=np.int64).reshape(-1)
    test_idx = np.asarray(test_idx, dtype=np.int64).reshape(-1)
    if train_idx.size == 0:
        raise ValueError(f"CV fold {fold_id!r} has no training rows")
    if test_idx.size == 0:
        raise ValueError(f"CV fold {fold_id!r} has no test rows")
    for name, idx in (("train_idx", train_idx), ("test_idx", test_idx)):
        if np.unique(idx).size != idx.size:
            raise ValueError(f"CV fold {fold_id!r} has duplicate {name} values")
        if idx.min() < 0 or idx.max() >= int(n_rows):
            raise ValueError(f"CV fold {fold_id!r} {name} is outside phenotype row bounds")
    if np.intersect1d(train_idx, test_idx).size:
        raise ValueError(f"CV fold {fold_id!r} has overlapping train/test rows")
    if not np.all(np.isfinite(y_values[train_idx])):
        raise ValueError(f"CV fold {fold_id!r} has non-finite training phenotypes")


def _cv_prediction_frame(predictions: Any, fold_id: str, test_idx: np.ndarray) -> pd.DataFrame:
    if predictions is None:
        pred_df = pd.DataFrame()
    elif isinstance(predictions, pd.DataFrame):
        pred_df = predictions.copy()
    else:
        pred_df = pd.DataFrame(predictions)
    test_idx = np.asarray(test_idx, dtype=np.int64).reshape(-1)

    if len(pred_df) != len(test_idx):
        label_col = next((c for c in ("Train_Test_label", "train_test_label", "set", "Set") if c in pred_df.columns), None)
        if label_col is not None:
            labels = pred_df[label_col].astype(str).str.lower()
            keep = labels.isin(("test", "prediction", "pred", "holdout"))
            if keep.any():
                pred_df = pred_df.loc[keep, :].copy()

    framework_row_index = None
    if "row_index" in pred_df.columns:
        framework_row_index = pd.to_numeric(pred_df["row_index"], errors="coerce")
        pred_df = pred_df.rename(columns={"row_index": "framework_row_index"})
    elif "row" in pred_df.columns:
        framework_row_index = pd.to_numeric(pred_df["row"], errors="coerce")

    if len(pred_df) == len(test_idx):
        stable_row_index = test_idx
    elif framework_row_index is not None and len(framework_row_index) == len(pred_df):
        stable_row_index = framework_row_index.to_numpy()
    else:
        stable_row_index = np.full(len(pred_df), np.nan)

    pred_df.insert(0, "set", "test")
    pred_df.insert(0, "row_index", stable_row_index)
    pred_df.insert(0, "cv_row_index", stable_row_index)
    pred_df.insert(0, "fold_id", str(fold_id))
    return pred_df


def _cv_fold_metrics(fold_id: str, pheno_df: pd.DataFrame, y_col: str, pred_df: pd.DataFrame) -> Dict[str, Any]:
    out: Dict[str, Any] = {"fold_id": str(fold_id), "n": 0, "rmse": np.nan, "mae": np.nan, "cor": np.nan}
    if pred_df is None or pred_df.empty or "row_index" not in pred_df.columns:
        return out
    pred_col = next((c for c in ("Prediction", "Predicted_value", "yhat") if c in pred_df.columns), None)
    if pred_col is None:
        return out
    row_index = pd.to_numeric(pred_df["row_index"], errors="coerce").to_numpy()
    yhat = pd.to_numeric(pred_df[pred_col], errors="coerce").to_numpy()
    valid_row = np.isfinite(row_index)
    if not valid_row.any():
        return out
    row_index = row_index[valid_row].astype(np.int64)
    yhat = yhat[valid_row]
    y_true_all = pd.to_numeric(pheno_df[y_col], errors="coerce").to_numpy()
    ok = (row_index >= 0) & (row_index < len(y_true_all)) & np.isfinite(yhat) & np.isfinite(y_true_all[row_index])
    if not ok.any():
        return out
    err = yhat[ok] - y_true_all[row_index[ok]]
    out["n"] = int(ok.sum())
    out["rmse"] = float(np.sqrt(np.mean(err * err)))
    out["mae"] = float(np.mean(np.abs(err)))
    if ok.sum() >= 2 and np.nanstd(yhat[ok]) > 0 and np.nanstd(y_true_all[row_index[ok]]) > 0:
        out["cor"] = float(np.corrcoef(y_true_all[row_index[ok]], yhat[ok])[0, 1])
    return out


def _env_enabled(name: str, default: str = "auto") -> bool:
    value = str(os.environ.get(name, default)).strip().lower()
    return value not in ("0", "false", "no", "n", "off", "disabled", "disable")


def _env_float(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, ""))
        if np.isfinite(value):
            return value
    except Exception:
        pass
    return float(default)


def _env_int(name: str, default: int) -> int:
    try:
        value = int(float(os.environ.get(name, "")))
        if value > 0:
            return value
    except Exception:
        pass
    return int(default)


@contextmanager
def _temporary_envvar(name: str, value: Optional[Any]):
    if value is None or str(value).strip() == "":
        yield
        return
    sentinel = object()
    old = os.environ.get(name, sentinel)
    os.environ[name] = str(value)
    try:
        yield
    finally:
        if old is sentinel:
            os.environ.pop(name, None)
        else:
            os.environ[name] = str(old)


def _none_like(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in ("", "none", "null", "false")
    return False


def _zero_weight(value: Any) -> bool:
    if value is None:
        return True
    try:
        arr = np.asarray(value, dtype=float).reshape(-1)
    except Exception:
        return False
    return arr.size == 0 or bool(np.all(np.isfinite(arr)) and np.max(np.abs(arr)) <= 1e-12)


def _first_numeric_column(df: pd.DataFrame, names: Sequence[str]) -> Optional[np.ndarray]:
    if df is None or df.empty:
        return None
    for name in names:
        if name in df.columns:
            values = pd.to_numeric(df[name], errors="coerce").to_numpy(dtype=float)
            return values
    return None


def _folds_are_delete_blocks(folds: Sequence[Dict[str, Any]], observed_rows: np.ndarray) -> bool:
    observed_sorted = np.sort(np.asarray(observed_rows, dtype=np.int64))
    for fold in folds:
        train_idx = np.sort(np.asarray(fold["train_idx"], dtype=np.int64).reshape(-1))
        test_idx = np.asarray(fold["test_idx"], dtype=np.int64).reshape(-1)
        expected_train = np.setdiff1d(observed_sorted, test_idx, assume_unique=False)
        if train_idx.size != expected_train.size or not np.array_equal(train_idx, expected_train):
            return False
    return True


def _normalized_fixed_effect_names(value: Sequence[Any]) -> List[str]:
    return [str(x).strip() for x in list(value or []) if str(x).strip()]


def _fixed_effects_are_env_only(value: Sequence[Any], env_col: str) -> bool:
    names = [x.lower() for x in _normalized_fixed_effect_names(value)]
    return len(names) == 1 and names[0] == str(env_col).strip().lower()


def _fast_krr_cv_ineligible_reason(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    env_col: str,
    method: str,
    backend: str,
    prediction_output: str,
    output_level: str,
    include_components: Sequence[Any],
    fixed_effects: Sequence[Any],
    env_similarity: Optional[np.ndarray],
    env_covariates: Optional[Any],
    grm_factor_cache: Optional[Any],
    kw_base: Dict[str, Any],
    y_values: np.ndarray,
    folds: Sequence[Dict[str, Any]],
    geno_kernel: np.ndarray,
) -> Optional[str]:
    if not _env_enabled("PREDICTPRO_GP_FAST_CV", "auto"):
        return "disabled"
    if not bool(getattr(fw, "TORCH_AVAILABLE", False)):
        return "torch_unavailable"
    if str(method).lower() != "krr_exact":
        return "method_not_krr_exact"
    if str(prediction_output).lower() != "test_only":
        return "prediction_output_not_test_only"
    output_level_l = str(output_level).lower()
    if output_level_l not in ("predict_only", "predict_with_se"):
        return "output_level_not_supported"
    if str(backend).lower() == "operator" or grm_factor_cache is not None:
        return "operator_or_factor_cache"
    if env_similarity is not None or env_covariates is not None:
        return "environment_kernel_present"
    comps = [str(x).lower() for x in list(include_components or [])]
    if comps and set(comps) != {"g"}:
        return "non_genetic_components_present"
    fixed_effect_names = _normalized_fixed_effect_names(fixed_effects)
    kw_fixed_effect_names = _normalized_fixed_effect_names(kw_base.get("fixed_effects", fixed_effects))
    if [x.lower() for x in fixed_effect_names] != [x.lower() for x in kw_fixed_effect_names]:
        return "fixed_effects_mismatch"
    if any(kw_base.get(name, None) is not None for name in ("obs_weights", "obs_var", "stage1_pev")):
        return "observation_specific_residual_inputs_present"
    if not _zero_weight(kw_base.get("w_ge", None)) or not _zero_weight(kw_base.get("w_e", None)):
        return "nonzero_environment_weights"
    if not _none_like(kw_base.get("krr_lams", "none")) and str(kw_base.get("lam_select", "fixed")).lower() != "fixed":
        return "lambda_selection_not_fixed"
    if str(kw_base.get("lam_select", "fixed")).lower() != "fixed":
        return "lambda_selection_not_fixed"
    if kw_base.get("learn_scales", False):
        return "learn_scales_enabled"
    if kw_base.get("compute_ai_se", False):
        return "ai_se_enabled"
    if kw_base.get("random_terms", None) is not None:
        return "random_terms_present"
    if kw_base.get("interaction_terms_meta", None) is not None:
        return "interaction_terms_present"
    if kw_base.get("env_structure", None) not in (None, "", "identity"):
        return "structured_environment_covariance"
    if str(env_col) not in pheno_df.columns:
        return "env_col_missing"
    observed_rows = np.where(np.isfinite(np.asarray(y_values, dtype=float)))[0].astype(np.int64)
    if observed_rows.size == 0:
        return "no_observed_rows"
    if not _folds_are_delete_blocks(folds, observed_rows):
        return "folds_not_observed_delete_blocks"
    df = pheno_df.reset_index(drop=True)
    env_values = df[str(env_col)].astype(str).to_numpy()
    n_env = int(pd.Series(env_values[observed_rows]).nunique(dropna=False))
    if n_env <= 1:
        if fixed_effect_names:
            return "fixed_effects_present"
    else:
        if not _fixed_effects_are_env_only(fixed_effect_names, env_col):
            return "multi_environment_requires_env_fixed_effect_only"
        if output_level_l != "predict_only":
            return "multi_environment_fast_cv_predict_only_only"
        observed_envs = set(str(x) for x in env_values[observed_rows])
        try:
            X_all, _ = fw._build_X_from_df(df, [str(env_col)])
        except Exception:
            return "fixed_effect_design_failed"
        X_all = np.asarray(X_all, dtype=float)
        for fold in folds:
            train_idx = np.asarray(fold["train_idx"], dtype=np.int64).reshape(-1)
            train_envs = set(str(x) for x in env_values[train_idx])
            if train_envs != observed_envs:
                return "multi_environment_fold_missing_env"
            X_train = X_all[train_idx, :]
            if X_train.shape[0] < X_train.shape[1]:
                return "fixed_effect_design_underdetermined"
            try:
                rank = int(np.linalg.matrix_rank(X_train))
            except Exception:
                return "fixed_effect_design_rank_failed"
            if rank < int(X_train.shape[1]):
                return "fixed_effect_design_rank_deficient"
    n = int(observed_rows.size)
    max_n = _env_int("PREDICTPRO_GP_FAST_CV_MAX_N", 8000)
    if n > max_n:
        return f"n_observed_gt_fast_limit_{max_n}"
    G = np.asarray(geno_kernel)
    if G.ndim != 2 or G.shape[0] != G.shape[1]:
        return "geno_kernel_not_square"
    if not np.all(np.isfinite(G)):
        return "geno_kernel_nonfinite"
    return None


def _resolve_stagewise_residuals_for_fast_cv(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    ei: np.ndarray,
    kw_base: Dict[str, Any],
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    resolver = getattr(fw, "_resolve_stagewise_residual_inputs", None)
    if resolver is None:
        return None, kw_base.get("resid_diag_env", None)
    resid_diag_obs, resid_diag_env, _meta = resolver(
        df=pheno_df,
        ei=ei,
        obs_weights=kw_base.get("obs_weights", None),
        obs_var=kw_base.get("obs_var", None),
        stage1_pev=kw_base.get("stage1_pev", None),
        obs_weight_mode=kw_base.get("obs_weight_mode", "relative_precision"),
        obs_weight_global_scale=float(kw_base.get("obs_weight_global_scale", 1.0)),
    )
    if resid_diag_env is None:
        resid_diag_env = kw_base.get("resid_diag_env", None)
    return resid_diag_obs, resid_diag_env


def _fast_krr_block_cv_predictions_learned_scalar_residual(
    *,
    fw: Any,
    torch: Any,
    df: pd.DataFrame,
    gid_values: np.ndarray,
    env_col: str,
    y_col: str,
    y: np.ndarray,
    ei: np.ndarray,
    observed_rows: np.ndarray,
    row_to_obs_pos: np.ndarray,
    X_obs_np: np.ndarray,
    K: Any,
    lam_eff: float,
    standardize: str,
    folds: Sequence[Dict[str, Any]],
    output_level: str,
    device: str,
    dtype_t: Any,
    kw_base: Dict[str, Any],
) -> Dict[str, Any]:
    n_obs = int(observed_rows.size)
    max_n = _env_int("PREDICTPRO_GP_FAST_CV_LEARNED_RESID_MAX_N", 3000)
    if n_obs > max_n:
        raise ValueError(f"fast KRR learned-residual CV limited to n_observed <= {max_n}; got {n_obs}")

    X_obs_t = fw._t(np.asarray(X_obs_np, dtype=float), device=device, dtype=dtype_t)
    K = 0.5 * (K + K.mT)
    eigvals, eigvecs = torch.linalg.eigh(K)
    eigvals = torch.clamp(eigvals, min=0.0)
    inv_lam = 1.0 / torch.clamp(eigvals + float(lam_eff), min=1e-12)
    C0 = (eigvecs * inv_lam.reshape(1, -1)) @ eigvecs.mT
    UX = eigvecs.mT @ X_obs_t

    return_se = str(output_level).lower() in ("predict_with_se", "full_vc")
    predictions_all: List[pd.DataFrame] = []
    fold_rows: List[Dict[str, Any]] = []
    metric_rows: List[Dict[str, Any]] = []
    started = time.perf_counter()
    sigma_rows: List[float] = []

    for fold in folds:
        fold_id = str(fold["fold_id"])
        fold_number = int(fold["fold_number"])
        train_idx = np.asarray(fold["train_idx"], dtype=np.int64).reshape(-1)
        test_idx = np.asarray(fold["test_idx"], dtype=np.int64).reshape(-1)
        block_pos_np = row_to_obs_pos[test_idx]
        train_pos_np = row_to_obs_pos[train_idx]
        if (block_pos_np < 0).any() or (train_pos_np < 0).any():
            raise ValueError(f"fast KRR CV: fold {fold_id!r} contains non-observed rows")
        t0 = time.perf_counter()

        stdr = fw._PerObsEnvStandardizer(mode=standardize).fit(y, ei, train_idx)
        y_std = stdr.transform(y, ei)
        y_obs_t = fw._t(np.asarray(y_std[observed_rows], dtype=float), device=device, dtype=dtype_t)

        block_pos_t = fw._t(block_pos_np, device=device, dtype=torch.int64)
        train_pos_t = fw._t(train_pos_np, device=device, dtype=torch.int64)
        X_train_t = X_obs_t.index_select(0, train_pos_t)
        y_train_t = y_obs_t.index_select(0, train_pos_t)

        C_BB0 = C0.index_select(0, block_pos_t).index_select(1, block_pos_t)
        C_TB0 = C0.index_select(0, train_pos_t).index_select(1, block_pos_t)
        C_TT0 = C0.index_select(0, train_pos_t).index_select(1, train_pos_t)
        CBB0_inv = fw._spd_inverse(0.5 * (C_BB0 + C_BB0.mT))
        Ainv_T = C_TT0 - C_TB0 @ CBB0_inv @ C_TB0.mT
        Ainv_T = 0.5 * (Ainv_T + Ainv_T.mT)

        Ainv_y = Ainv_T @ y_train_t
        Ainv_X = Ainv_T @ X_train_t
        Var_beta0 = fw._spd_inverse(X_train_t.mT @ Ainv_X)
        beta0 = Var_beta0 @ (X_train_t.mT @ Ainv_y)
        Py_train = Ainv_y - Ainv_X @ beta0
        resid = float(lam_eff) * Py_train.reshape(-1)
        hutch_samples = int(kw_base.get("hutch_samples", 64))
        if hutch_samples > 0:
            K_TT0 = K.index_select(0, train_pos_t).index_select(1, train_pos_t)
            train_idx_full_t = fw._t(train_idx, device=device, dtype=torch.int64)

            def H_times(z):
                zt = z.index_select(0, train_idx_full_t)
                Ainv_z = Ainv_T @ zt
                Sx_zt = X_train_t @ (Var_beta0 @ (X_train_t.mT @ Ainv_z))
                Hz_t = Sx_zt + K_TT0 @ (Ainv_T @ (zt - Sx_zt))
                out = torch.zeros_like(z)
                out.index_copy_(0, train_idx_full_t, Hz_t)
                return out

            hat_diag = fw._estimate_diag_Hutch(
                H_times,
                n=int(df.shape[0]),
                probes=hutch_samples,
                device=device,
                dtype=dtype_t,
                seed=int(kw_base.get("seed", 12345)),
            )
            tr_H = float(
                torch.clamp(
                    hat_diag.index_select(0, train_idx_full_t).sum(),
                    min=0.0,
                    max=float(max(int(train_idx.size) - 1, 0)),
                ).item()
            )
        else:
            tr_H = 0.0
        denom = max(1.0, float(train_idx.size) - tr_H)
        if str(kw_base.get("env_resid_robust", "none")).lower() == "mad" and hasattr(fw, "_robust_env_var"):
            sig2 = float(fw._robust_env_var(resid))
            if not np.isfinite(sig2):
                sig2 = float(torch.clamp((resid.pow(2).sum() / denom), min=1e-12).item())
        else:
            sig2 = float(torch.clamp((resid.pow(2).sum() / denom), min=1e-12).item())
        sigma_rows.append(sig2)

        c_eff = float(lam_eff) + sig2
        inv_c = 1.0 / torch.clamp(eigvals + c_eff, min=1e-12)
        Uy = eigvecs.mT @ y_obs_t
        Cy = eigvecs @ (inv_c * Uy)
        CX = eigvecs @ (inv_c.reshape(-1, 1) * UX)
        Var_beta = fw._spd_inverse(X_obs_t.mT @ CX)
        beta = Var_beta @ (X_obs_t.mT @ Cy)
        Py = Cy - CX @ beta

        U_B = eigvecs.index_select(0, block_pos_t)
        C_BB = (U_B * inv_c.reshape(1, -1)) @ U_B.mT
        CX_B = CX.index_select(0, block_pos_t)
        PBB = C_BB - CX_B @ Var_beta @ CX_B.mT
        PBB = 0.5 * (PBB + PBB.mT)
        Lbb = fw._cholesky_factor(PBB)
        cv_resid_std = fw._cholesky_solve_from_factor(
            Lbb,
            Py.index_select(0, block_pos_t).reshape(-1, 1),
        ).reshape(-1)
        pred_std = y_obs_t.index_select(0, block_pos_t) - cv_resid_std
        pred_obs = stdr.inv_mean(fw._to_numpy(pred_std), ei[test_idx])

        pred_df = pd.DataFrame(
            {
                "Name": [str(gid_values[i]) for i in test_idx],
                "Env": [str(df[str(env_col)].iloc[i]) for i in test_idx],
                "Prediction": np.asarray(pred_obs, dtype=float),
            }
        )
        if return_se:
            U_T = eigvecs.index_select(0, train_pos_t)
            C_TT = (U_T * inv_c.reshape(1, -1)) @ U_T.mT
            C_TB = (U_T * inv_c.reshape(1, -1)) @ U_B.mT
            CBB_inv = fw._spd_inverse(0.5 * (C_BB + C_BB.mT))
            Ainv_T_c = C_TT - C_TB @ CBB_inv @ C_TB.mT
            Ainv_T_c = 0.5 * (Ainv_T_c + Ainv_T_c.mT)
            Ainv_X_c = Ainv_T_c @ X_train_t
            Var_beta_c = fw._spd_inverse(X_train_t.mT @ Ainv_X_c)
            K_TB = K.index_select(0, train_pos_t).index_select(1, block_pos_t)
            Ainv_K = Ainv_T_c @ K_TB
            proj = (K_TB * Ainv_K).sum(dim=0)
            Kss_diag = K.diagonal().index_select(0, block_pos_t)
            var_lat_std = torch.clamp(Kss_diag - proj, min=0.0)
            X_B = X_obs_t.index_select(0, block_pos_t)
            XB = X_B @ Var_beta_c
            var_beta_diag = (XB * X_B).sum(dim=1)
            var_obs_std = torch.clamp(var_lat_std + var_beta_diag + sig2, min=0.0)
            var_obs = stdr.inv_var(fw._to_numpy(var_obs_std), ei[test_idx])
            var_lat = stdr.inv_var(fw._to_numpy(var_lat_std), ei[test_idx])
            pred_df["SE_observed"] = np.sqrt(np.maximum(var_obs, 0.0))
            pred_df["SE_latent"] = np.sqrt(np.maximum(var_lat, 0.0))
            pred_df["SE_g_latent"] = pred_df["SE_latent"]
            pred_df["Prediction_Var_latent"] = np.maximum(var_lat, 0.0)
            pred_df["Prediction_Var_observed"] = np.maximum(var_obs, 0.0)
            pred_df["Prediction_SE_latent"] = pred_df["SE_latent"]
            pred_df["Prediction_SE_observed"] = pred_df["SE_observed"]

        pred_df = _cv_prediction_frame(pred_df, fold_id, test_idx)
        predictions_all.append(pred_df)
        elapsed = time.perf_counter() - t0
        fold_rows.append(
            {
                "fold_id": fold_id,
                "fold_number": fold_number,
                "n_train": int(train_idx.size),
                "n_test": int(test_idx.size),
                "elapsed_sec": float(elapsed),
                "status": "ok",
                "fast_path": "krr_block_delete_learned_residual",
                "sigma2_resid_std": float(sig2),
            }
        )
        metric_rows.append(_cv_fold_metrics(fold_id, df, str(y_col), pred_df))

    combined = pd.concat(predictions_all, axis=0, ignore_index=True) if predictions_all else pd.DataFrame()
    return {
        "predictions": combined,
        "fold_rows": fold_rows,
        "metric_rows": metric_rows,
        "meta": {
            "cv_fast_path": "krr_block_delete_learned_residual",
            "cv_fast_path_device": str(device),
            "cv_fast_path_n_observed": int(n_obs),
            "cv_fast_path_fixed_residual": False,
            "cv_fast_path_sigma2_resid_std_mean": float(np.mean(sigma_rows)) if sigma_rows else np.nan,
            "cv_fast_path_elapsed_sec": float(time.perf_counter() - started),
            "cv_fast_path_self_checked": False,
        },
    }


def _fast_krr_block_cv_predictions_learned_env_residual(
    *,
    fw: Any,
    torch: Any,
    df: pd.DataFrame,
    gid_values: np.ndarray,
    env_col: str,
    y_col: str,
    y: np.ndarray,
    ei: np.ndarray,
    observed_rows: np.ndarray,
    row_to_obs_pos: np.ndarray,
    X_obs_np: np.ndarray,
    K: Any,
    lam_eff: float,
    standardize: str,
    folds: Sequence[Dict[str, Any]],
    output_level: str,
    device: str,
    dtype_t: Any,
    kw_base: Dict[str, Any],
) -> Dict[str, Any]:
    if str(output_level).lower() != "predict_only":
        raise ValueError("fast KRR multi-environment learned-residual CV currently supports predict_only only")

    n_obs = int(observed_rows.size)
    max_n = _env_int("PREDICTPRO_GP_FAST_CV_LEARNED_ENV_RESID_MAX_N", 3000)
    if n_obs > max_n:
        raise ValueError(f"fast KRR learned-env-residual CV limited to n_observed <= {max_n}; got {n_obs}")

    X_obs_t = fw._t(np.asarray(X_obs_np, dtype=float), device=device, dtype=dtype_t)
    K = 0.5 * (K + K.mT)

    eigvals, eigvecs = torch.linalg.eigh(K)
    eigvals = torch.clamp(eigvals, min=0.0)
    inv_lam = 1.0 / torch.clamp(eigvals + float(lam_eff), min=1e-12)
    C0 = (eigvecs * inv_lam.reshape(1, -1)) @ eigvecs.mT

    n_env = int(np.max(ei) + 1) if np.asarray(ei).size else 0
    predictions_all: List[pd.DataFrame] = []
    fold_rows: List[Dict[str, Any]] = []
    metric_rows: List[Dict[str, Any]] = []
    sigma_rows: List[np.ndarray] = []
    started = time.perf_counter()

    for fold in folds:
        fold_id = str(fold["fold_id"])
        fold_number = int(fold["fold_number"])
        train_idx = np.asarray(fold["train_idx"], dtype=np.int64).reshape(-1)
        test_idx = np.asarray(fold["test_idx"], dtype=np.int64).reshape(-1)
        train_pos_np = row_to_obs_pos[train_idx]
        block_pos_np = row_to_obs_pos[test_idx]
        if (train_pos_np < 0).any() or (block_pos_np < 0).any():
            raise ValueError(f"fast KRR CV: fold {fold_id!r} contains non-observed rows")
        t0 = time.perf_counter()

        stdr = fw._PerObsEnvStandardizer(mode=standardize).fit(y, ei, train_idx)
        y_std = stdr.transform(y, ei)
        y_obs_t = fw._t(np.asarray(y_std[observed_rows], dtype=float), device=device, dtype=dtype_t)

        train_pos_t = fw._t(train_pos_np, device=device, dtype=torch.int64)
        block_pos_t = fw._t(block_pos_np, device=device, dtype=torch.int64)
        X_train_t = X_obs_t.index_select(0, train_pos_t)
        X_test_t = X_obs_t.index_select(0, block_pos_t)
        y_train_t = y_obs_t.index_select(0, train_pos_t)

        C_BB0 = C0.index_select(0, block_pos_t).index_select(1, block_pos_t)
        C_TB0 = C0.index_select(0, train_pos_t).index_select(1, block_pos_t)
        C_TT0 = C0.index_select(0, train_pos_t).index_select(1, train_pos_t)
        CBB0_inv = fw._spd_inverse(0.5 * (C_BB0 + C_BB0.mT))
        Ainv_T = C_TT0 - C_TB0 @ CBB0_inv @ C_TB0.mT
        Ainv_T = 0.5 * (Ainv_T + Ainv_T.mT)

        Ainv_y = Ainv_T @ y_train_t
        Ainv_X = Ainv_T @ X_train_t
        Var_beta0 = fw._spd_inverse(X_train_t.mT @ Ainv_X)
        beta0 = Var_beta0 @ (X_train_t.mT @ Ainv_y)
        Py_train = Ainv_y - Ainv_X @ beta0
        resid = float(lam_eff) * Py_train.reshape(-1)

        hutch_samples = int(kw_base.get("hutch_samples", 64))
        if hutch_samples > 0:
            K_TT0 = K.index_select(0, train_pos_t).index_select(1, train_pos_t)
            train_idx_full_t = fw._t(train_idx, device=device, dtype=torch.int64)

            def H_times(z):
                zt = z.index_select(0, train_idx_full_t)
                Ainv_z = Ainv_T @ zt
                Sx_zt = X_train_t @ (Var_beta0 @ (X_train_t.mT @ Ainv_z))
                Hz_t = Sx_zt + K_TT0 @ (Ainv_T @ (zt - Sx_zt))
                out = torch.zeros_like(z)
                out.index_copy_(0, train_idx_full_t, Hz_t)
                return out

            hat_diag = fw._estimate_diag_Hutch(
                H_times,
                n=int(df.shape[0]),
                probes=hutch_samples,
                device=device,
                dtype=dtype_t,
                seed=int(kw_base.get("seed", 12345)),
            )
        else:
            hat_diag = torch.zeros(int(df.shape[0]), dtype=dtype_t, device=device)

        train_env = np.asarray(ei[train_idx], dtype=np.int64)
        sig2 = np.full(n_env, np.nan, dtype=float)
        robust_mode = str(kw_base.get("env_resid_robust", "none")).lower()
        for env_i in range(n_env):
            rel_np = np.where(train_env == env_i)[0].astype(np.int64)
            n_e = int(rel_np.size)
            if n_e == 0:
                continue
            rel_t = fw._t(rel_np, device=device, dtype=torch.int64)
            r_e = resid.index_select(0, rel_t)
            idx_full_t = fw._t(train_idx[rel_np], device=device, dtype=torch.int64)
            tr_He = float(
                torch.clamp(
                    hat_diag.index_select(0, idx_full_t).sum(),
                    min=0.0,
                    max=float(max(n_e - 1, 0)),
                ).item()
            )
            denom = max(1.0, float(n_e) - tr_He)
            if robust_mode == "mad" and hasattr(fw, "_robust_env_var"):
                s2 = float(fw._robust_env_var(r_e))
                if not np.isfinite(s2):
                    s2 = float(torch.clamp((r_e.pow(2).sum() / denom), min=1e-12).item())
            else:
                s2 = float(torch.clamp((r_e.pow(2).sum() / denom), min=1e-12).item())
            sig2[env_i] = s2

        shrink_tau = float(kw_base.get("env_resid_shrink_tau", 0.0) or 0.0)
        if np.isfinite(sig2).any() and shrink_tau > 0:
            m = float(np.nanmean(sig2))
            sig2 = (1.0 - shrink_tau) * sig2 + shrink_tau * m
        if np.isnan(sig2).any():
            resid_np = fw._to_numpy(resid)
            mask = np.isfinite(resid_np)
            global_s2 = float(np.mean(resid_np[mask] ** 2)) if np.any(mask) else 1.0
            sig2 = np.where(np.isfinite(sig2), sig2, global_s2)
        if np.any(~np.isfinite(sig2)) or np.any(sig2 <= 0):
            raise ValueError("fast KRR CV: learned environment residuals must be finite and > 0")
        sigma_rows.append(sig2.copy())

        resid_train_std_t = fw._t(sig2[train_env], device=device, dtype=dtype_t)
        K_TT = K.index_select(0, train_pos_t).index_select(1, train_pos_t)
        Vtt = K_TT.clone()
        Vtt = 0.5 * (Vtt + Vtt.mT)
        Vtt.diagonal().add_(resid_train_std_t + float(lam_eff))
        L = fw._cholesky_factor(Vtt)
        Vinv_y = fw._cholesky_solve_from_factor(L, y_train_t)
        Vinv_X = fw._cholesky_solve_from_factor(L, X_train_t)
        Xt_Vinv_X = X_train_t.mT @ Vinv_X
        try:
            Var_beta = fw._spd_inverse(Xt_Vinv_X)
        except RuntimeError:
            eye = torch.eye(Xt_Vinv_X.shape[0], dtype=dtype_t, device=device)
            Var_beta = fw._spd_inverse(Xt_Vinv_X + 1e-10 * eye)
        beta = Var_beta @ (X_train_t.mT @ Vinv_y)
        alpha = fw._cholesky_solve_from_factor(L, y_train_t - X_train_t @ beta)
        K_TS = K.index_select(0, train_pos_t).index_select(1, block_pos_t)
        mean_lat_std = K_TS.mT @ alpha + X_test_t @ beta
        pred_obs = stdr.inv_mean(fw._to_numpy(mean_lat_std), ei[test_idx])

        pred_df = pd.DataFrame(
            {
                "Name": [str(gid_values[i]) for i in test_idx],
                "Env": [str(df[str(env_col)].iloc[i]) for i in test_idx],
                "Prediction": np.asarray(pred_obs, dtype=float),
            }
        )
        pred_df = _cv_prediction_frame(pred_df, fold_id, test_idx)
        predictions_all.append(pred_df)
        elapsed = time.perf_counter() - t0
        fold_rows.append(
            {
                "fold_id": fold_id,
                "fold_number": fold_number,
                "n_train": int(train_idx.size),
                "n_test": int(test_idx.size),
                "elapsed_sec": float(elapsed),
                "status": "ok",
                "fast_path": "krr_block_delete_learned_env_residual",
                "sigma2_resid_env_std_mean": float(np.mean(sig2)),
                "sigma2_resid_env_std_min": float(np.min(sig2)),
                "sigma2_resid_env_std_max": float(np.max(sig2)),
            }
        )
        metric_rows.append(_cv_fold_metrics(fold_id, df, str(y_col), pred_df))

    combined = pd.concat(predictions_all, axis=0, ignore_index=True) if predictions_all else pd.DataFrame()
    sigma_stack = np.vstack(sigma_rows) if sigma_rows else np.empty((0, 0), dtype=float)
    return {
        "predictions": combined,
        "fold_rows": fold_rows,
        "metric_rows": metric_rows,
        "meta": {
            "cv_fast_path": "krr_block_delete_learned_env_residual",
            "cv_fast_path_device": str(device),
            "cv_fast_path_n_observed": int(n_obs),
            "cv_fast_path_fixed_residual": False,
            "cv_fast_path_residual_mode": "per_environment",
            "cv_fast_path_final_solve": "fold_specific_cholesky",
            "cv_fast_path_sigma2_resid_std_mean": float(np.mean(sigma_stack)) if sigma_stack.size else np.nan,
            "cv_fast_path_elapsed_sec": float(time.perf_counter() - started),
            "cv_fast_path_self_checked": False,
        },
    }


def _fast_krr_block_cv_predictions(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    folds: Sequence[Dict[str, Any]],
    kw_base: Dict[str, Any],
    output_level: str,
) -> Dict[str, Any]:
    torch = fw.torch
    df = pheno_df.reset_index(drop=True).copy()
    y = pd.to_numeric(df[str(y_col)], errors="coerce").to_numpy(dtype=float)
    observed_rows = np.where(np.isfinite(y))[0].astype(np.int64)
    n_obs = int(observed_rows.size)
    row_to_obs_pos = np.full(df.shape[0], -1, dtype=np.int64)
    row_to_obs_pos[observed_rows] = np.arange(n_obs, dtype=np.int64)

    geno_ids = [str(x) for x in list(geno_ids or [])]
    gid_to_index = {str(g): i for i, g in enumerate(geno_ids)}
    gid_values = df[str(gid_col)].astype(str).to_numpy()
    gi = np.asarray([gid_to_index.get(str(g), -1) for g in gid_values], dtype=np.int64)
    if (gi < 0).any():
        missing = [str(g) for g, i in zip(gid_values, gi) if i < 0][:5]
        raise ValueError(f"fast KRR CV: unknown genotype IDs encountered: {missing} ...")

    env = df[str(env_col)].astype("category")
    env_levels = list(env.cat.categories.astype(str))
    env_to_index = {e: i for i, e in enumerate(env_levels)}
    ei = np.asarray([env_to_index[str(e)] for e in env.astype(str)], dtype=np.int64)

    fixed_effects = list(kw_base.get("fixed_effects", []) or [])
    X_all, _ = fw._build_X_from_df(df, fixed_effects)
    X_obs_np = np.asarray(X_all, dtype=float)[observed_rows, :]

    standardize = str(kw_base.get("standardize", "global"))
    stdr = fw._PerObsEnvStandardizer(mode=standardize).fit(y, ei, observed_rows)
    y_std = stdr.transform(y, ei)
    y_obs_np = np.asarray(y_std[observed_rows], dtype=float)

    dtype_t = fw._resolve_torch_dtype(kw_base.get("dtype", "float64"))
    work_units = float(n_obs) * float(n_obs + max(1, n_obs // max(1, len(folds))))
    device = fw._runtime_gp_device(kw_base.get("device", None), work_units=work_units)

    G_list, _ = fw._apply_kernel_normalization(
        [np.asarray(geno_kernel, dtype=np.float64)],
        None,
        mode=str(kw_base.get("reml_normalize", "diag_mean")),
    )
    G_t = fw._t(np.asarray(G_list[0], dtype=np.float64), device=device, dtype=dtype_t)
    gi_obs_t = fw._t(gi[observed_rows], device=device, dtype=torch.int64)
    wg = fw._prepare_weight_vector(kw_base.get("w_g", None), 1, default=1.0)
    K = float(np.asarray(wg, dtype=float).reshape(-1)[0]) * fw._build_geno_block(G_t, gi_obs_t, gi_obs_t)
    K = 0.5 * (K + K.mT)

    y_obs_t = fw._t(y_obs_np, device=device, dtype=dtype_t)
    X_obs_t = fw._t(X_obs_np, device=device, dtype=dtype_t)
    lam_eff = float(_coalesce(kw_base.get("lam", None), kw_base.get("krr_lam", None), default=1e-3))

    resid_diag_obs, resid_diag_env = _resolve_stagewise_residuals_for_fast_cv(
        fw=fw,
        pheno_df=df,
        ei=ei,
        kw_base=kw_base,
    )
    fixed_residual = resid_diag_obs is not None or resid_diag_env is not None

    if not fixed_residual:
        if len(env_levels) > 1:
            if not _fixed_effects_are_env_only(fixed_effects, str(env_col)):
                raise ValueError("fast KRR multi-environment learned-residual CV requires Env as the only fixed effect")
            return _fast_krr_block_cv_predictions_learned_env_residual(
                fw=fw,
                torch=torch,
                df=df,
                gid_values=gid_values,
                env_col=str(env_col),
                y_col=str(y_col),
                y=y,
                ei=ei,
                observed_rows=observed_rows,
                row_to_obs_pos=row_to_obs_pos,
                X_obs_np=X_obs_np,
                K=K,
                lam_eff=lam_eff,
                standardize=standardize,
                folds=folds,
                output_level=output_level,
                device=device,
                dtype_t=dtype_t,
                kw_base=kw_base,
            )
        return _fast_krr_block_cv_predictions_learned_scalar_residual(
            fw=fw,
            torch=torch,
            df=df,
            gid_values=gid_values,
            env_col=str(env_col),
            y_col=str(y_col),
            y=y,
            ei=ei,
            observed_rows=observed_rows,
            row_to_obs_pos=row_to_obs_pos,
            X_obs_np=X_obs_np,
            K=K,
            lam_eff=lam_eff,
            standardize=standardize,
            folds=folds,
            output_level=output_level,
            device=device,
            dtype_t=dtype_t,
            kw_base=kw_base,
        )

    if resid_diag_obs is not None:
        resid_diag_obs = np.asarray(resid_diag_obs, dtype=float).reshape(-1)
        if resid_diag_obs.size != df.shape[0]:
            raise ValueError(f"fast KRR CV: resid_diag_obs length {resid_diag_obs.size} != n_obs_all {df.shape[0]}")
        scale2_obs = stdr.inv_var(np.ones(df.shape[0], dtype=float), ei)
        scale2_obs = np.where(np.isfinite(scale2_obs) & (scale2_obs > 0), scale2_obs, 1.0)
        resid_std_np = resid_diag_obs[observed_rows] / scale2_obs[observed_rows]
    elif resid_diag_env is not None:
        resid_diag_env = np.asarray(resid_diag_env, dtype=float).reshape(-1)
        if resid_diag_env.size != len(env_levels):
            raise ValueError(f"fast KRR CV: resid_diag_env length {resid_diag_env.size} != n_env {len(env_levels)}")
        if np.any(~np.isfinite(resid_diag_env)) or np.any(resid_diag_env <= 0):
            raise ValueError("fast KRR CV: resid_diag_env must be finite and > 0")
        resid_std_np = resid_diag_env[ei[observed_rows]]
    else:
        Vtmp = K.clone()
        Vtmp = 0.5 * (Vtmp + Vtmp.mT)
        Vtmp.diagonal().add_(float(lam_eff))
        Vinv_y_tmp = fw._cholesky_solve(Vtmp, y_obs_t)
        Vinv_X_tmp = fw._cholesky_solve(Vtmp, X_obs_t)
        Xt_Vinv_X_tmp = X_obs_t.mT @ Vinv_X_tmp
        Var_beta_tmp = fw._spd_inverse(Xt_Vinv_X_tmp)
        beta_tmp = Var_beta_tmp @ (X_obs_t.mT @ Vinv_y_tmp)
        alpha_tmp = fw._cholesky_solve(Vtmp, y_obs_t - X_obs_t @ beta_tmp)
        yhat_tmp = (X_obs_t @ beta_tmp) + K @ alpha_tmp
        resid = (y_obs_t - yhat_tmp).reshape(-1)

        hutch_samples = int(kw_base.get("hutch_samples", 64))
        if hutch_samples > 0:
            def H_times(z):
                Vinv_z = fw._cholesky_solve(Vtmp, z)
                Sx_z = X_obs_t @ (Var_beta_tmp @ (X_obs_t.mT @ Vinv_z))
                return Sx_z + K @ fw._cholesky_solve(Vtmp, z - Sx_z)

            hat_diag = fw._estimate_diag_Hutch(
                H_times,
                n=n_obs,
                probes=hutch_samples,
                device=device,
                dtype=dtype_t,
                seed=int(kw_base.get("seed", 12345)),
            )
        else:
            hat_diag = torch.zeros(n_obs, dtype=dtype_t, device=device)

        tr_H = float(torch.clamp(hat_diag.sum(), min=0.0, max=float(max(n_obs - 1, 0))).item())
        denom = max(1.0, float(n_obs) - tr_H)
        if str(kw_base.get("env_resid_robust", "none")).lower() == "mad" and hasattr(fw, "_robust_env_var"):
            sig2 = float(fw._robust_env_var(resid))
            if not np.isfinite(sig2):
                sig2 = float(torch.clamp((resid.pow(2).sum() / denom), min=1e-12).item())
        else:
            sig2 = float(torch.clamp((resid.pow(2).sum() / denom), min=1e-12).item())
        resid_std_np = np.full(n_obs, sig2, dtype=float)

    if np.any(~np.isfinite(resid_std_np)) or np.any(resid_std_np <= 0):
        raise ValueError("fast KRR CV: residual diagonal must be finite and > 0")
    resid_std_t = fw._t(np.asarray(resid_std_np, dtype=float), device=device, dtype=dtype_t)

    V = K.clone()
    V = 0.5 * (V + V.mT)
    V.diagonal().add_(resid_std_t + float(lam_eff))
    L = fw._cholesky_factor(V)
    Vinv_y = fw._cholesky_solve_from_factor(L, y_obs_t)
    Vinv_X = fw._cholesky_solve_from_factor(L, X_obs_t)
    Xt_Vinv_X = X_obs_t.mT @ Vinv_X
    Var_beta = fw._spd_inverse(Xt_Vinv_X)
    beta = Var_beta @ (X_obs_t.mT @ Vinv_y)
    Py = Vinv_y - Vinv_X @ beta

    return_se = str(output_level).lower() in ("predict_with_se", "full_vc")
    predictions_all: List[pd.DataFrame] = []
    fold_rows: List[Dict[str, Any]] = []
    metric_rows: List[Dict[str, Any]] = []
    started = time.perf_counter()

    for fold in folds:
        fold_id = str(fold["fold_id"])
        fold_number = int(fold["fold_number"])
        test_idx = np.asarray(fold["test_idx"], dtype=np.int64).reshape(-1)
        block_pos_np = row_to_obs_pos[test_idx]
        if (block_pos_np < 0).any():
            raise ValueError(f"fast KRR CV: fold {fold_id!r} contains non-observed test rows")
        t0 = time.perf_counter()
        block_pos_t = fw._t(block_pos_np, device=device, dtype=torch.int64)
        b = int(block_pos_np.size)
        E = torch.zeros((n_obs, b), dtype=dtype_t, device=device)
        E[block_pos_t, torch.arange(b, device=device, dtype=torch.int64)] = 1.0
        Vinv_E = fw._cholesky_solve_from_factor(L, E)
        Vinv_BB = Vinv_E.index_select(0, block_pos_t)
        Vinv_X_B = Vinv_X.index_select(0, block_pos_t)
        PBB = Vinv_BB - Vinv_X_B @ Var_beta @ Vinv_X_B.mT
        PBB = 0.5 * (PBB + PBB.mT)
        Lbb = fw._cholesky_factor(PBB)
        cv_resid_std = fw._cholesky_solve_from_factor(
            Lbb,
            Py.index_select(0, block_pos_t).reshape(-1, 1),
        ).reshape(-1)
        pred_std = y_obs_t.index_select(0, block_pos_t) - cv_resid_std
        pred_obs = stdr.inv_mean(fw._to_numpy(pred_std), ei[test_idx])

        pred_df = pd.DataFrame(
            {
                "Name": [str(gid_values[i]) for i in test_idx],
                "Env": [str(df[str(env_col)].iloc[i]) for i in test_idx],
                "Prediction": np.asarray(pred_obs, dtype=float),
            }
        )
        if return_se:
            PBB_inv = torch.cholesky_inverse(Lbb)
            var_obs_std = torch.clamp(PBB_inv.diagonal() - float(lam_eff), min=0.0)
            resid_test_std = resid_std_t.index_select(0, block_pos_t)
            var_lat_std = torch.clamp(var_obs_std - resid_test_std, min=0.0)
            var_obs = stdr.inv_var(fw._to_numpy(var_obs_std), ei[test_idx])
            var_lat = stdr.inv_var(fw._to_numpy(var_lat_std), ei[test_idx])
            pred_df["SE_observed"] = np.sqrt(np.maximum(var_obs, 0.0))
            pred_df["SE_latent"] = np.sqrt(np.maximum(var_lat, 0.0))
            pred_df["Prediction_Var_latent"] = np.maximum(var_lat, 0.0)
            pred_df["Prediction_Var_observed"] = np.maximum(var_obs, 0.0)
            pred_df["Prediction_SE_latent"] = pred_df["SE_latent"]
            pred_df["Prediction_SE_observed"] = pred_df["SE_observed"]

        pred_df = _cv_prediction_frame(pred_df, fold_id, test_idx)
        predictions_all.append(pred_df)
        elapsed = time.perf_counter() - t0
        fold_rows.append(
            {
                "fold_id": fold_id,
                "fold_number": fold_number,
                "n_train": int(np.asarray(fold["train_idx"], dtype=np.int64).size),
                "n_test": int(test_idx.size),
                "elapsed_sec": float(elapsed),
                "status": "ok",
                "fast_path": "krr_block_delete",
            }
        )
        metric_rows.append(_cv_fold_metrics(fold_id, df, str(y_col), pred_df))

    combined = pd.concat(predictions_all, axis=0, ignore_index=True) if predictions_all else pd.DataFrame()
    return {
        "predictions": combined,
        "fold_rows": fold_rows,
        "metric_rows": metric_rows,
        "meta": {
            "cv_fast_path": "krr_block_delete",
            "cv_fast_path_device": str(device),
            "cv_fast_path_n_observed": int(n_obs),
            "cv_fast_path_fixed_residual": bool(fixed_residual),
            "cv_fast_path_elapsed_sec": float(time.perf_counter() - started),
            "cv_fast_path_self_checked": False,
        },
    }


def _fast_cv_self_check_passed(
    *,
    fw: Any,
    fast_result: Dict[str, Any],
    folds: Sequence[Dict[str, Any]],
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    kw_base: Dict[str, Any],
) -> Tuple[bool, str]:
    if not folds:
        return False, "no_folds"
    fold = folds[0]
    kw = dict(kw_base)
    kw["train_idx"] = np.asarray(fold["train_idx"], dtype=np.int64)
    kw["test_idx"] = np.asarray(fold["test_idx"], dtype=np.int64)
    ref = fw.fit_mixed_model(
        pheno_df=pheno_df,
        gid_col=str(gid_col),
        env_col=str(env_col),
        y_col=str(y_col),
        geno_kernels={"G": geno_kernel},
        geno_ids=list(geno_ids or []),
        **kw,
    )
    ref_pred = ref.get("result", {}).get("predictions")
    if ref_pred is None:
        ref_pred = ref.get("predictions")
    ref_df = _cv_prediction_frame(ref_pred, str(fold["fold_id"]), np.asarray(fold["test_idx"], dtype=np.int64))
    fast_df = fast_result["predictions"]
    fast_df = fast_df.loc[fast_df["fold_id"].astype(str) == str(fold["fold_id"]), :].copy()
    if ref_df.empty or fast_df.empty:
        return False, "empty_self_check_predictions"
    ref_df = ref_df.sort_values("row_index").reset_index(drop=True)
    fast_df = fast_df.sort_values("row_index").reset_index(drop=True)
    if len(ref_df) != len(fast_df):
        return False, "self_check_row_count_mismatch"
    if not np.array_equal(
        pd.to_numeric(ref_df["row_index"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(fast_df["row_index"], errors="coerce").to_numpy(dtype=float),
    ):
        return False, "self_check_row_index_mismatch"

    atol = _env_float("PREDICTPRO_GP_FAST_CV_ATOL", 1e-5)
    rtol = _env_float("PREDICTPRO_GP_FAST_CV_RTOL", 1e-5)
    pred_ref = _first_numeric_column(ref_df, ("Prediction", "Predicted_value", "yhat"))
    pred_fast = _first_numeric_column(fast_df, ("Prediction", "Predicted_value", "yhat"))
    if pred_ref is None or pred_fast is None:
        return False, "self_check_missing_prediction"
    pred_ok = np.isfinite(pred_ref) & np.isfinite(pred_fast)
    if not pred_ok.any():
        return False, "self_check_no_finite_predictions"
    pred_diff = np.abs(pred_fast[pred_ok] - pred_ref[pred_ok])
    pred_tol = atol + rtol * np.maximum(1.0, np.abs(pred_ref[pred_ok]))
    pred_check_reason = "ok"
    if bool(np.any(pred_diff > pred_tol)):
        finite_ref = pred_ref[pred_ok][np.isfinite(pred_ref[pred_ok])]
        scale_candidates: List[float] = []
        if finite_ref.size > 1:
            scale_candidates.append(float(np.nanstd(finite_ref)))
            if finite_ref.size >= 4:
                q75, q25 = np.nanpercentile(finite_ref, [75.0, 25.0])
                iqr_scale = float((q75 - q25) / 1.349) if np.isfinite(q75 - q25) else np.nan
                scale_candidates.append(iqr_scale)
        scale_floor = max(0.0, _env_float("PREDICTPRO_GP_FAST_CV_SCALE_FLOOR", 1.0))
        finite_scales = [x for x in scale_candidates if np.isfinite(x) and x > 0]
        pred_scale = max([scale_floor] + finite_scales)
        max_scale_rtol = max(0.0, _env_float("PREDICTPRO_GP_FAST_CV_MAX_SCALE_RTOL", 0.01))
        rmse_scale_rtol = max(0.0, _env_float("PREDICTPRO_GP_FAST_CV_RMSE_SCALE_RTOL", 0.005))
        max_scale_atol = max(0.0, _env_float("PREDICTPRO_GP_FAST_CV_MAX_SCALE_ATOL", 0.0))
        rmse_scale_atol = max(0.0, _env_float("PREDICTPRO_GP_FAST_CV_RMSE_SCALE_ATOL", 0.0))
        max_allowed = max_scale_atol + max_scale_rtol * pred_scale
        rmse_allowed = rmse_scale_atol + rmse_scale_rtol * pred_scale
        pred_max_abs = float(np.nanmax(pred_diff))
        pred_rmse = float(np.sqrt(np.nanmean(pred_diff * pred_diff)))
        if not (
            np.isfinite(pred_max_abs)
            and np.isfinite(pred_rmse)
            and pred_max_abs <= max_allowed
            and pred_rmse <= rmse_allowed
        ):
            return False, f"prediction_mismatch_max_abs_{pred_max_abs:.6g}"
        pred_check_reason = (
            "scale_tolerated_prediction_mismatch"
            f"_max_abs_{pred_max_abs:.6g}"
            f"_rmse_{pred_rmse:.6g}"
        )

    var_ref = _first_numeric_column(ref_df, ("Prediction_Var_observed", "Prediction_error_variance", "PEV"))
    var_fast = _first_numeric_column(fast_df, ("Prediction_Var_observed", "Prediction_error_variance", "PEV"))
    if var_ref is not None and var_fast is not None:
        var_ok = np.isfinite(var_ref) & np.isfinite(var_fast)
        if var_ok.any():
            var_diff = np.abs(var_fast[var_ok] - var_ref[var_ok])
            var_tol = (10.0 * atol) + (10.0 * rtol) * np.maximum(1.0, np.abs(var_ref[var_ok]))
            if bool(np.any(var_diff > var_tol)):
                return False, f"variance_mismatch_max_abs_{float(np.nanmax(var_diff)):.6g}"
    return True, pred_check_reason


def _try_fast_krr_block_cv(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    y_values: np.ndarray,
    folds: Sequence[Dict[str, Any]],
    method: str,
    backend: str,
    prediction_output: str,
    output_level: str,
    include_components: Sequence[Any],
    fixed_effects: Sequence[Any],
    env_similarity: Optional[np.ndarray],
    env_covariates: Optional[Any],
    grm_factor_cache: Optional[Any],
    kw_base: Dict[str, Any],
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    mode = str(os.environ.get("PREDICTPRO_GP_FAST_CV", "auto")).strip().lower()
    force = mode in ("force", "forced", "always")
    reason = _fast_krr_cv_ineligible_reason(
        fw=fw,
        pheno_df=pheno_df,
        env_col=env_col,
        method=method,
        backend=backend,
        prediction_output=prediction_output,
        output_level=output_level,
        include_components=include_components,
        fixed_effects=fixed_effects,
        env_similarity=env_similarity,
        env_covariates=env_covariates,
        grm_factor_cache=grm_factor_cache,
        kw_base=kw_base,
        y_values=y_values,
        folds=folds,
        geno_kernel=geno_kernel,
    )
    if reason is not None and not force:
        return None, reason
    try:
        result = _fast_krr_block_cv_predictions(
            fw=fw,
            pheno_df=pheno_df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_kernel=geno_kernel,
            geno_ids=geno_ids,
            folds=folds,
            kw_base=kw_base,
            output_level=output_level,
        )
    except Exception as exc:
        if force:
            raise
        warnings.warn(f"Fast KRR block CV failed; falling back to per-fold fits: {exc}", RuntimeWarning, stacklevel=2)
        return None, f"fast_path_failed:{type(exc).__name__}"

    # Self-check is on by default for every output level, including predict_only
    # (the common CV case). It validates the block-delete result against an exact
    # per-fold refit on fold 0 and falls back to full per-fold CV on mismatch, so
    # the fast path can never silently diverge from the exact estimator.
    self_check_default = "true"
    self_check = _env_enabled("PREDICTPRO_GP_FAST_CV_SELF_CHECK", self_check_default) and not force
    if self_check:
        ok, check_reason = _fast_cv_self_check_passed(
            fw=fw,
            fast_result=result,
            folds=folds,
            pheno_df=pheno_df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_kernel=geno_kernel,
            geno_ids=geno_ids,
            kw_base=kw_base,
        )
        if not ok:
            warnings.warn(
                f"Fast KRR block CV self-check failed ({check_reason}); falling back to per-fold fits.",
                RuntimeWarning,
                stacklevel=2,
            )
            return None, f"self_check_failed:{check_reason}"
        result.setdefault("meta", {})["cv_fast_path_self_checked"] = True
        result.setdefault("meta", {})["cv_fast_path_self_check_reason"] = check_reason
    return result, None


def _fast_gp_exact_cv_ineligible_reason(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    env_col: str,
    method: str,
    backend: str,
    prediction_output: str,
    output_level: str,
    include_components: Sequence[Any],
    fixed_effects: Sequence[Any],
    env_similarity: Optional[np.ndarray],
    env_covariates: Optional[Any],
    grm_factor_cache: Optional[Any],
    kw_base: Dict[str, Any],
    y_values: np.ndarray,
    folds: Sequence[Dict[str, Any]],
    geno_kernel: np.ndarray,
) -> Optional[str]:
    if not _env_enabled("PREDICTPRO_GP_FAST_CV", "auto"):
        return "disabled"
    if not _env_enabled("PREDICTPRO_GP_EXACT_FAST_CV", "false"):
        return "gp_exact_fast_cv_disabled"
    if not bool(getattr(fw, "TORCH_AVAILABLE", False)) or not bool(getattr(fw, "GPTY_AVAILABLE", False)):
        return "torch_or_gpytorch_unavailable"
    if str(method).lower() != "gp_exact":
        return "method_not_gp_exact"
    if str(prediction_output).lower() != "test_only":
        return "prediction_output_not_test_only"
    output_level_l = str(output_level).lower()
    if output_level_l not in ("predict_only", "predict_with_se"):
        return "output_level_not_supported"
    if output_level_l != "predict_only":
        return "gp_exact_fast_cv_predict_with_se_requires_force"
    if str(backend).lower() == "operator" or grm_factor_cache is not None:
        return "operator_or_factor_cache"
    if env_similarity is not None or env_covariates is not None:
        return "environment_kernel_present"
    comps = [str(x).lower() for x in list(include_components or [])]
    if comps and set(comps) != {"g"}:
        return "non_genetic_components_present"
    if fixed_effects is not None and len(list(fixed_effects)) > 0:
        return "fixed_effects_present"
    if any(kw_base.get(name, None) is not None for name in ("obs_weights", "obs_var", "stage1_pev")):
        return "observation_specific_residual_inputs_present"
    if kw_base.get("learn_scales", False):
        return "learn_scales_enabled"
    if kw_base.get("compute_ai_se", False):
        return "ai_se_enabled"
    if kw_base.get("random_terms", None) is not None:
        return "random_terms_present"
    if kw_base.get("interaction_terms_meta", None) is not None:
        return "interaction_terms_present"
    if kw_base.get("env_structure", None) not in (None, "", "identity"):
        return "structured_environment_covariance"
    if str(env_col) not in pheno_df.columns:
        return "env_col_missing"
    if pheno_df[str(env_col)].astype(str).nunique(dropna=False) != 1:
        return "not_single_environment"
    observed_rows = np.where(np.isfinite(np.asarray(y_values, dtype=float)))[0].astype(np.int64)
    if observed_rows.size == 0:
        return "no_observed_rows"
    if not _folds_are_delete_blocks(folds, observed_rows):
        return "folds_not_observed_delete_blocks"
    n = int(observed_rows.size)
    max_n = _env_int("PREDICTPRO_GP_EXACT_FAST_CV_MAX_N", 3000)
    if n > max_n:
        return f"n_observed_gt_gp_exact_fast_limit_{max_n}"
    G = np.asarray(geno_kernel)
    if G.ndim != 2 or G.shape[0] != G.shape[1]:
        return "geno_kernel_not_square"
    if not np.all(np.isfinite(G)):
        return "geno_kernel_nonfinite"
    return None


def _try_fast_gp_exact_shared_noise_block_cv(
    *,
    fw: Any,
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernel: np.ndarray,
    geno_ids: Sequence[str],
    y_values: np.ndarray,
    folds: Sequence[Dict[str, Any]],
    method: str,
    backend: str,
    prediction_output: str,
    output_level: str,
    include_components: Sequence[Any],
    fixed_effects: Sequence[Any],
    env_similarity: Optional[np.ndarray],
    env_covariates: Optional[Any],
    grm_factor_cache: Optional[Any],
    kw_base: Dict[str, Any],
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    mode = str(os.environ.get("PREDICTPRO_GP_EXACT_FAST_CV", "false")).strip().lower()
    force = mode in ("force", "forced", "always")
    reason = _fast_gp_exact_cv_ineligible_reason(
        fw=fw,
        pheno_df=pheno_df,
        env_col=env_col,
        method=method,
        backend=backend,
        prediction_output=prediction_output,
        output_level=output_level,
        include_components=include_components,
        fixed_effects=fixed_effects,
        env_similarity=env_similarity,
        env_covariates=env_covariates,
        grm_factor_cache=grm_factor_cache,
        kw_base=kw_base,
        y_values=y_values,
        folds=folds,
        geno_kernel=geno_kernel,
    )
    if reason is not None and not force:
        return None, reason

    observed_rows = np.where(np.isfinite(np.asarray(y_values, dtype=float)))[0].astype(np.int64)
    try:
        fit_kw = dict(kw_base)
        fit_kw.update(
            {
                "method": "gp_exact",
                "prediction_output": "test_only",
                "output_level": "predict_only",
                "point_predictions_only": True,
                "return_se": False,
                "compute_ai_se": False,
                "train_idx": observed_rows,
                "test_idx": np.array([], dtype=np.int64),
            }
        )
        full_fit = fw.fit_mixed_model(
            pheno_df=pheno_df,
            gid_col=str(gid_col),
            env_col=str(env_col),
            y_col=str(y_col),
            geno_kernels={"G": geno_kernel},
            geno_ids=list(geno_ids or []),
            **fit_kw,
        )
        diagnostics = dict(full_fit.get("diagnostics", {}) or {})
        noise_variance = float(diagnostics.get("likelihood_noise", np.nan))
        if not np.isfinite(noise_variance) or noise_variance <= 0:
            raise ValueError(f"shared exact-GP noise variance was not finite and positive: {noise_variance!r}")

        block_kw = dict(kw_base)
        block_kw["method"] = "krr_exact"
        block_kw["krr_lam"] = 0.0
        block_kw["lam"] = 0.0
        block_kw["resid_diag_env"] = np.array([noise_variance], dtype=float)
        block_kw["learn_envdiag_noise"] = False
        result = _fast_krr_block_cv_predictions(
            fw=fw,
            pheno_df=pheno_df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_kernel=geno_kernel,
            geno_ids=geno_ids,
            folds=folds,
            kw_base=block_kw,
            output_level=output_level,
        )
    except Exception as exc:
        if force:
            raise
        warnings.warn(
            f"Fast exact-GP shared-noise CV failed; falling back to per-fold fits: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )
        return None, f"fast_gp_exact_path_failed:{type(exc).__name__}"

    result_meta = result.setdefault("meta", {})
    result_meta.update(
        {
            "cv_fast_path": "gp_exact_shared_noise_block_delete",
            "cv_fast_path_model": "gp_exact",
            "cv_fast_path_approximate": True,
            "cv_fast_path_exact_refit": False,
            "cv_fast_path_solver": "krr_block_delete",
            "cv_fast_path_shared_hyperparameters": True,
            "cv_fast_path_likelihood_noise_variance": float(noise_variance),
            # Compatibility alias: the value has always been a variance.
            "cv_fast_path_likelihood_noise_std": float(noise_variance),
            "cv_fast_path_self_checked": False,
        }
    )
    for row in result.get("fold_rows", []) or []:
        if isinstance(row, dict):
            row["fast_path"] = "gp_exact_shared_noise_block_delete"

    self_check = _env_enabled("PREDICTPRO_GP_FAST_CV_SELF_CHECK", "true") and not force
    if self_check:
        ok, check_reason = _fast_cv_self_check_passed(
            fw=fw,
            fast_result=result,
            folds=folds,
            pheno_df=pheno_df,
            gid_col=gid_col,
            env_col=env_col,
            y_col=y_col,
            geno_kernel=geno_kernel,
            geno_ids=geno_ids,
            kw_base=kw_base,
        )
        if not ok:
            warnings.warn(
                f"Fast exact-GP shared-noise CV self-check failed ({check_reason}); falling back to per-fold fits.",
                RuntimeWarning,
                stacklevel=2,
            )
            return None, f"self_check_failed:{check_reason}"
        result_meta["cv_fast_path_self_checked"] = True
        result_meta["cv_fast_path_self_check_reason"] = check_reason
    return result, None


def _fit_predict_cv_spec(args, fw=None) -> int:
    spec_path = Path(args.spec).resolve()
    spec = _read_json_file(str(spec_path))
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root"), default=None)
    if project_root is None:
        raise ValueError("project_root must be supplied in the spec or with --project-root")
    if fw is None:
        fw = _load_framework(_resolve_path(project_root, base_dir, required=True, label="project_root"))

    pheno_csv = _resolve_path(
        _spec_get(spec, "pheno_csv", "phenotype_csv", "inputs.pheno_csv", "inputs.phenotype_csv"),
        base_dir,
        required=True,
        label="pheno_csv",
    )
    pheno_df = pd.read_csv(pheno_csv).reset_index(drop=True)

    geno_kernels = _read_geno_kernels_from_spec(spec, base_dir)
    first_kernel_name = next(iter(geno_kernels))
    geno_kernel = geno_kernels[first_kernel_name]
    geno_ids = _read_lines_from_spec(
        _spec_get(spec, "geno_ids", "genotype_ids", "inputs.geno_ids", "inputs.genotype_ids"),
        base_dir,
        "geno_ids",
        required=True,
    )

    cols = _spec_get(spec, "columns", default={}) or {}
    gid_col = _coalesce(_spec_get(spec, "gid_col"), cols.get("gid"), cols.get("gid_col"), default="GID")
    env_col = _coalesce(_spec_get(spec, "env_col"), cols.get("env"), cols.get("env_col"), default="Env")
    y_col = _coalesce(_spec_get(spec, "y_col"), cols.get("y"), cols.get("y_col"), default="Trait")
    if str(y_col) not in pheno_df.columns:
        raise ValueError(f"y_col {y_col!r} is not present in pheno_csv")
    y_values = pd.to_numeric(pheno_df[str(y_col)], errors="coerce").to_numpy()

    fit_options = dict(_spec_get(spec, "fit", "options", "model", default={}) or {})
    output_level = str(_coalesce(_spec_get(spec, "output_level"), fit_options.pop("output_level", None), default="predict_only"))
    method = str(_coalesce(_spec_get(spec, "method"), fit_options.pop("method", None), default="krr_exact"))
    backend = str(_coalesce(_spec_get(spec, "backend"), fit_options.pop("backend", None), default="auto"))
    prediction_output = str(_coalesce(_spec_get(spec, "prediction_output"), fit_options.pop("prediction_output", None), default="test_only"))
    varcomp_mode = str(_coalesce(_spec_get(spec, "varcomp_mode"), fit_options.pop("varcomp_mode", None), default="mom"))
    seed = int(_coalesce(_spec_get(spec, "seed"), fit_options.pop("seed", None), default=12345))
    fa_rank = int(_coalesce(_spec_get(spec, "fa_rank"), fit_options.pop("fa_rank", None), default=1))
    gp_exact_fast_cv = _coalesce(
        _spec_get(spec, "gp_exact_fast_cv", "exact_fast_cv"),
        fit_options.pop("gp_exact_fast_cv", None),
        default=None,
    )

    env_covariates = None
    env_covariates_csv = _resolve_path(
        _spec_get(spec, "env_covariates_csv", "inputs.env_covariates_csv", "environment_covariates_csv", "env_covariates.csv", "environment_covariates.csv"),
        base_dir,
        label="env_covariates_csv",
    )
    if env_covariates_csv is not None:
        env_covariates = pd.read_csv(env_covariates_csv)

    env_similarity_spec = _coalesce(
        _spec_get(spec, "env_similarity", "environment_similarity", "inputs.env_similarity"),
        {
            "bin": _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin"),
            "meta": _spec_get(spec, "env_similarity_meta", "inputs.env_similarity_meta"),
        }
        if _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin") is not None
        else None,
    )
    env_similarity = _read_matrix_from_spec(env_similarity_spec, base_dir, "env_similarity")
    env_ids = _read_lines_from_spec(
        _spec_get(spec, "env_ids", "environment_ids", "inputs.env_ids", "inputs.environment_ids", "env_similarity.env_ids", "env_similarity.ids", "environment_similarity.env_ids"),
        base_dir,
        "env_ids",
        required=False,
    )
    env_similarity = _align_env_similarity_if_needed(env_similarity, env_ids, pheno_df, str(env_col))

    factor_cache_spec = _coalesce(
        _spec_get(spec, "factor_cache_spec", "grm_factor_cache", "gp_factor_cache", "inputs.factor_cache_spec"),
        _spec_get(spec, "grm_factor_cache_root", "gp_factor_cache_root"),
    )
    if isinstance(factor_cache_spec, str) and not factor_cache_spec.lower().endswith(".json"):
        factor_cache_spec = {"type": "zarr", "root_dir": factor_cache_spec}
    grm_factor_cache = _normalize_factor_cache(factor_cache_spec, base_dir)

    include_components = _coalesce(
        _spec_get(spec, "include_components"),
        fit_options.pop("include_components", None),
        default=(["g", "ge", "e"] if (env_similarity is not None or env_covariates is not None) else ["g"]),
    )
    fixed_effects = _coalesce(_spec_get(spec, "fixed_effects"), fit_options.pop("fixed_effects", None), default=[])

    kw_base: Dict[str, Any] = {
        "include_components": list(include_components),
        "fixed_effects": list(fixed_effects),
        "method": method,
        "backend": backend,
        "prediction_output": prediction_output,
        "output_level": output_level,
        "point_predictions_only": output_level == "predict_only",
        "return_se": output_level in ("predict_with_se", "full_vc"),
        "compute_ai_se": output_level == "full_vc",
        "standardize": str(_coalesce(_spec_get(spec, "standardize"), fit_options.pop("standardize", None), default="global")),
        "varcomp_mode": varcomp_mode,
        "dtype": str(_coalesce(_spec_get(spec, "dtype"), fit_options.pop("dtype", None), default="float64")),
        "seed": seed,
    }
    if method == "gp_icm_fa":
        kw_base["env_structure"] = str(_coalesce(_spec_get(spec, "env_structure"), fit_options.pop("env_structure", None), default="fa"))
        kw_base["fa_rank"] = fa_rank
    if env_similarity is not None:
        kw_base["env_similarity"] = env_similarity
    if env_covariates is not None:
        kw_base["env_covariates"] = env_covariates
    if grm_factor_cache is not None:
        kw_base["grm_factor_cache"] = grm_factor_cache
    kw_base.update(fit_options)

    folds_spec = _spec_get(spec, "folds", "cv_folds", "inputs.folds", default=None)
    if not isinstance(folds_spec, list) or not folds_spec:
        raise ValueError("fit-predict-cv-spec requires a non-empty folds list")

    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))
    out_dir.mkdir(parents=True, exist_ok=True)

    fold_defs: List[Dict[str, Any]] = []
    for fold_number, fold_spec in enumerate(folds_spec, start=1):
        if not isinstance(fold_spec, dict):
            raise ValueError(f"folds[{fold_number}] must be an object")
        fold_id = str(_coalesce(fold_spec.get("fold_id"), fold_spec.get("id"), fold_spec.get("fold"), default=f"fold{fold_number}"))
        train_idx = _read_cv_idx_from_spec(_coalesce(fold_spec.get("train_idx"), fold_spec.get("train_index")), base_dir, f"{fold_id}.train_idx")
        test_idx = _read_cv_idx_from_spec(_coalesce(fold_spec.get("test_idx"), fold_spec.get("test_index")), base_dir, f"{fold_id}.test_idx")
        _validate_cv_fold_indices(
            fold_id=fold_id,
            train_idx=train_idx,
            test_idx=test_idx,
            n_rows=len(pheno_df),
            y_values=y_values,
        )
        fold_defs.append(
            {
                "fold_id": fold_id,
                "fold_number": int(fold_number),
                "train_idx": np.asarray(train_idx, dtype=np.int64).reshape(-1),
                "test_idx": np.asarray(test_idx, dtype=np.int64).reshape(-1),
            }
        )

    started = time.perf_counter()
    if len(geno_kernels) == 1:
        with _temporary_envvar("PREDICTPRO_GP_EXACT_FAST_CV", gp_exact_fast_cv):
            fast_cv, fast_skip_reason = _try_fast_krr_block_cv(
                fw=fw,
                pheno_df=pheno_df,
                gid_col=str(gid_col),
                env_col=str(env_col),
                y_col=str(y_col),
                geno_kernel=geno_kernel,
                geno_ids=list(geno_ids or []),
                y_values=y_values,
                folds=fold_defs,
                method=method,
                backend=backend,
                prediction_output=prediction_output,
                output_level=output_level,
                include_components=include_components,
                fixed_effects=fixed_effects,
                env_similarity=env_similarity,
                env_covariates=env_covariates,
                grm_factor_cache=grm_factor_cache,
                kw_base=kw_base,
            )
            if fast_cv is None and str(method).lower() == "gp_exact":
                fast_cv, fast_skip_reason = _try_fast_gp_exact_shared_noise_block_cv(
                    fw=fw,
                    pheno_df=pheno_df,
                    gid_col=str(gid_col),
                    env_col=str(env_col),
                    y_col=str(y_col),
                    geno_kernel=geno_kernel,
                    geno_ids=list(geno_ids or []),
                    y_values=y_values,
                    folds=fold_defs,
                    method=method,
                    backend=backend,
                    prediction_output=prediction_output,
                    output_level=output_level,
                    include_components=include_components,
                    fixed_effects=fixed_effects,
                    env_similarity=env_similarity,
                    env_covariates=env_covariates,
                    grm_factor_cache=grm_factor_cache,
                    kw_base=kw_base,
                )
    else:
        fast_cv, fast_skip_reason = None, "multi_kernel_fast_cv_not_supported"

    if fast_cv is not None:
        combined = fast_cv["predictions"]
        fold_rows = fast_cv["fold_rows"]
        metric_rows = fast_cv["metric_rows"]
        fast_meta = dict(fast_cv.get("meta", {}) or {})
    else:
        predictions_all: List[pd.DataFrame] = []
        fold_rows: List[Dict[str, Any]] = []
        metric_rows: List[Dict[str, Any]] = []
        fast_meta = {
            "cv_fast_path": "not_used",
            "cv_fast_path_skip_reason": fast_skip_reason or "unknown",
            "cv_preprocessed_context": "not_used",
        }
        shared_context = None
        if hasattr(fw, "prepare_fit_mixed_model_context"):
            try:
                shared_context = fw.prepare_fit_mixed_model_context(
                    pheno_df=pheno_df,
                    gid_col=str(gid_col),
                    env_col=str(env_col),
                    y_col=str(y_col),
                    geno_ids=list(geno_ids or []),
                    fixed_effects=list(fixed_effects or []),
                )
                fast_meta["cv_preprocessed_context"] = "used"
            except Exception as exc:
                shared_context = None
                fast_meta["cv_preprocessed_context"] = "failed"
                fast_meta["cv_preprocessed_context_error"] = f"{type(exc).__name__}: {exc}"
        for fold in fold_defs:
            fold_id = str(fold["fold_id"])
            fold_number = int(fold["fold_number"])
            train_idx = np.asarray(fold["train_idx"], dtype=np.int64)
            test_idx = np.asarray(fold["test_idx"], dtype=np.int64)
            t0 = time.perf_counter()
            kw = dict(kw_base)
            kw["train_idx"] = train_idx
            kw["test_idx"] = test_idx
            if shared_context is not None:
                kw["preprocessed_context"] = shared_context
            res = fw.fit_mixed_model(
                pheno_df=pheno_df,
                gid_col=str(gid_col),
                env_col=str(env_col),
                y_col=str(y_col),
                geno_kernels=geno_kernels,
                geno_ids=list(geno_ids or []),
                **kw,
            )
            pred = res.get("result", {}).get("predictions")
            if pred is None:
                pred = res.get("predictions")
            pred_df = _cv_prediction_frame(pred, fold_id, test_idx)
            predictions_all.append(pred_df)
            elapsed = time.perf_counter() - t0
            fold_rows.append(
                {
                    "fold_id": fold_id,
                    "fold_number": int(fold_number),
                    "n_train": int(train_idx.size),
                    "n_test": int(test_idx.size),
                    "elapsed_sec": float(elapsed),
                    "status": "ok",
                }
            )
            metric_rows.append(_cv_fold_metrics(fold_id, pheno_df, str(y_col), pred_df))
        if shared_context is not None:
            tensor_cache = shared_context.get("gp_exact_tensor_cache", {})
            fast_meta["cv_preprocessed_context_rows"] = int(shared_context.get("n_rows", len(pheno_df)))
            fast_meta["cv_preprocessed_context_tensor_cache_entries"] = int(
                len(tensor_cache) if isinstance(tensor_cache, dict) else 0
            )
        combined = pd.concat(predictions_all, axis=0, ignore_index=True) if predictions_all else pd.DataFrame()

    _write_df(combined, str(out_dir / "predictions.csv"))
    _write_df(pd.DataFrame(fold_rows), str(out_dir / "folds.csv"))
    _write_df(pd.DataFrame(metric_rows), str(out_dir / "metrics.csv"))
    meta = {
        "python": sys.executable,
        "command": "fit-predict-cv-spec",
        "spec": str(spec_path),
        "method": method,
        "backend": backend,
        "output_level": output_level,
        "prediction_output": prediction_output,
        "varcomp_mode": varcomp_mode,
        "n_folds": len(fold_rows),
        "elapsed_sec": float(time.perf_counter() - started),
        "kernel_names": list(geno_kernels.keys()),
        "has_env_covariates": env_covariates is not None,
        "has_env_similarity": env_similarity is not None,
        "has_grm_factor_cache": grm_factor_cache is not None,
        "gp_exact_fast_cv": str(gp_exact_fast_cv) if gp_exact_fast_cv is not None else None,
    }
    meta.update(fast_meta)
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=_json_default), encoding="utf-8")
    return 0


def _fit_batch_candidates_spec(args, fw=None) -> int:
    spec_path = Path(args.spec)
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    if not isinstance(spec, dict):
        raise ValueError("batch-fit-candidates-spec requires a JSON object spec")
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root", "python_project_root"))
    if project_root is None:
        raise ValueError("project_root is required")
    project_root = _resolve_path(project_root, base_dir, required=True, label="project_root")
    if fw is None:
        fw = _load_framework(project_root)

    pheno_csv = _resolve_path(
        _spec_get(spec, "pheno_csv", "phenotype_csv", "inputs.pheno_csv", "inputs.phenotype_csv"),
        base_dir,
        required=True,
        label="pheno_csv",
    )
    pheno_df = pd.read_csv(pheno_csv)
    kernel_spec = _coalesce(
        _spec_get(spec, "geno_kernel", "kernel", "inputs.geno_kernel", "inputs.kernel"),
        {
            "bin": _spec_get(spec, "kernel_bin", "geno_kernel_bin", "inputs.kernel_bin", "inputs.geno_kernel_bin"),
            "meta": _spec_get(spec, "kernel_meta", "geno_kernel_meta", "inputs.kernel_meta", "inputs.geno_kernel_meta"),
        },
    )
    geno_kernel = _read_matrix_from_spec(kernel_spec, base_dir, "geno_kernel")
    if geno_kernel is None:
        raise ValueError("geno_kernel is required")
    geno_ids = _read_lines_from_spec(
        _spec_get(spec, "geno_ids", "genotype_ids", "inputs.geno_ids", "inputs.genotype_ids"),
        base_dir,
        "geno_ids",
        required=True,
    )
    train_idx = _read_idx_from_spec(
        _spec_get(spec, "train_idx", "train_index", "inputs.train_idx", "inputs.train_index"),
        base_dir,
        "train_idx",
        required=True,
    )
    test_idx = _read_idx_from_spec(
        _spec_get(spec, "test_idx", "test_index", "inputs.test_idx", "inputs.test_index"),
        base_dir,
        "test_idx",
        required=True,
    )
    candidates_csv = _resolve_path(
        _spec_get(spec, "candidates_csv", "inputs.candidates_csv"),
        base_dir,
        required=True,
        label="candidates_csv",
    )
    candidates = pd.read_csv(candidates_csv)

    cols = _spec_get(spec, "columns", default={}) or {}
    gid_col = _coalesce(_spec_get(spec, "gid_col"), cols.get("gid"), cols.get("gid_col"), default="GID")
    env_col = _coalesce(_spec_get(spec, "env_col"), cols.get("env"), cols.get("env_col"), default="Env")
    y_col = _coalesce(_spec_get(spec, "y_col"), cols.get("y"), cols.get("y_col"), default="y")

    env_covariates = None
    env_covariates_csv = _resolve_path(
        _spec_get(spec, "env_covariates_csv", "inputs.env_covariates_csv", "environment_covariates_csv"),
        base_dir,
        label="env_covariates_csv",
    )
    if env_covariates_csv is not None:
        env_covariates = pd.read_csv(env_covariates_csv)
    env_similarity_by_kernel = _read_matrix_mapping_from_spec(
        _spec_get(spec, "env_similarity_by_kernel", "environment_similarity_by_kernel", "inputs.env_similarity_by_kernel"),
        base_dir,
        "env_similarity_by_kernel",
    )

    fit_options = dict(_spec_get(spec, "fit", "options", "model", default={}) or {})
    include_components = _coalesce(
        _spec_get(spec, "include_components"),
        fit_options.pop("include_components", None),
        default=["g", "ge", "e"],
    )
    fixed_effects = _coalesce(_spec_get(spec, "fixed_effects"), fit_options.pop("fixed_effects", None), default=[])
    backend = str(_coalesce(_spec_get(spec, "backend"), fit_options.pop("backend", None), default="auto"))
    prediction_output = str(_coalesce(_spec_get(spec, "prediction_output"), fit_options.pop("prediction_output", None), default="test_only"))
    seed = int(_coalesce(_spec_get(spec, "seed"), fit_options.pop("seed", None), default=12345))
    started = time.perf_counter()
    res = fw.batch_fit_mixed_model_candidates(
        pheno_df=pheno_df,
        gid_col=str(gid_col),
        env_col=str(env_col),
        y_col=str(y_col),
        geno_kernel=geno_kernel,
        geno_ids=list(geno_ids or []),
        train_idx=train_idx,
        test_idx=test_idx,
        candidates=candidates,
        env_covariates=env_covariates,
        env_similarity_by_kernel=env_similarity_by_kernel,
        include_components=list(include_components),
        fixed_effects=list(fixed_effects or []),
        backend=backend,
        reaction_norm_feature_qc=bool(_coalesce(
            _spec_get(spec, "reaction_norm_feature_qc"),
            fit_options.pop("reaction_norm_feature_qc", None),
            default=True,
        )),
        reaction_norm_env_filter_topk=_coalesce(
            _spec_get(spec, "reaction_norm_env_filter_topk"),
            fit_options.pop("reaction_norm_env_filter_topk", None),
            default=15,
        ),
        kenv_bandwidth=float(_coalesce(
            _spec_get(spec, "kenv_bandwidth"),
            fit_options.pop("kenv_bandwidth", None),
            default=1.0,
        )),
        kenv_kernel_kwargs=_coalesce(
            _spec_get(spec, "kenv_kernel_kwargs"),
            fit_options.pop("kenv_kernel_kwargs", None),
            default=None,
        ),
        dtype=str(_coalesce(_spec_get(spec, "dtype"), fit_options.pop("dtype", None), default="float64")),
        seed=seed,
        standardize=str(_coalesce(_spec_get(spec, "standardize"), fit_options.pop("standardize", None), default="global")),
        prediction_output=prediction_output,
    )
    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_df(res.get("predictions"), str(out_dir / "predictions.csv"))
    _write_df(res.get("detail"), str(out_dir / "detail.csv"))
    _write_df(res.get("env_kernel_cache_info"), str(out_dir / "env_kernel_cache_info.csv"))
    meta = {
        "python": sys.executable,
        "command": "batch-fit-candidates-spec",
        "spec": str(spec_path),
        "backend": backend,
        "prediction_output": prediction_output,
        "n_candidates": int(candidates.shape[0]),
        "elapsed_sec": float(time.perf_counter() - started),
        "has_env_covariates": env_covariates is not None,
        "has_env_similarity_by_kernel": env_similarity_by_kernel is not None,
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=_json_default), encoding="utf-8")
    return 0


def _load_project_module(project_root: str, module_name: str):
    module_path = os.path.join(project_root, module_name + ".py")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Python module file not found: {module_path}")
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    importlib.invalidate_caches()
    return importlib.import_module(module_name)


def _factor_cache_object(cache_spec: Optional[Dict[str, Any]], project_root: str):
    if cache_spec is None:
        return None
    backend = _load_project_module(project_root, "gp_large_backend")
    ctype = str(cache_spec.get("type", "zarr")).lower()
    if ctype == "zarr":
        return backend.GRMFactorCacheZarr(
            root_dir=str(cache_spec["root_dir"]),
            num_grms=int(cache_spec.get("num_grms", 1)),
        )
    if ctype == "memmap":
        return backend.GRMFactorCacheMemmap(
            paths=cache_spec["paths"],
            shapes=[tuple(int(v) for v in shape) for shape in cache_spec["shapes"]],
            dtypes=[np.dtype(dt) for dt in cache_spec["dtypes"]],
        )
    raise ValueError(f"Unknown grm_factor_cache type: {ctype}")


def _read_common_public_spec(args) -> Tuple[Dict[str, Any], Path, str, pd.DataFrame, Dict[str, np.ndarray], List[str], np.ndarray, np.ndarray, Path]:
    spec_path = Path(args.spec).resolve()
    spec = _read_json_file(str(spec_path))
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root"), default=None)
    if project_root is None:
        raise ValueError("project_root must be supplied in the spec or with --project-root")
    project_root = _resolve_path(project_root, base_dir, required=True, label="project_root")

    pheno_csv = _resolve_path(
        _spec_get(spec, "pheno_csv", "phenotype_csv", "inputs.pheno_csv", "inputs.phenotype_csv"),
        base_dir,
        required=True,
        label="pheno_csv",
    )
    pheno_df = pd.read_csv(pheno_csv)
    geno_kernels = _read_geno_kernels_from_spec(spec, base_dir)
    geno_ids = _read_lines_from_spec(
        _spec_get(spec, "geno_ids", "genotype_ids", "inputs.geno_ids", "inputs.genotype_ids"),
        base_dir,
        "geno_ids",
        required=True,
    )
    train_idx = _read_idx_from_spec(
        _spec_get(spec, "train_idx", "train_index", "inputs.train_idx", "inputs.train_index"),
        base_dir,
        "train_idx",
        required=True,
    )
    test_idx = _read_idx_from_spec(
        _spec_get(spec, "test_idx", "test_index", "inputs.test_idx", "inputs.test_index"),
        base_dir,
        "test_idx",
        required=True,
    )
    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))
    return spec, base_dir, str(project_root), pheno_df, geno_kernels, list(geno_ids or []), train_idx, test_idx, out_dir


def _fit_multi_trait_spec(args) -> int:
    spec, base_dir, project_root, pheno_df, geno_kernels, geno_ids, train_idx, test_idx, out_dir = _read_common_public_spec(args)
    mt_gp = _load_project_module(project_root, "mt_gp")
    cols = _spec_get(spec, "columns", default={}) or {}
    gid_col = _coalesce(_spec_get(spec, "gid_col"), cols.get("gid"), cols.get("gid_col"), default="gid")
    trait_col = _coalesce(_spec_get(spec, "trait_col"), cols.get("trait"), cols.get("trait_col"), default="trait")
    y_col = _coalesce(_spec_get(spec, "y_col"), cols.get("y"), cols.get("y_col"), default="y")
    fit_options = dict(_spec_get(spec, "fit", "fit_args", "options", "model", default={}) or {})
    factor_cache_spec = _coalesce(
        _spec_get(spec, "factor_cache_spec", "grm_factor_cache", "gp_factor_cache", "inputs.factor_cache_spec"),
        _spec_get(spec, "grm_factor_cache_root", "gp_factor_cache_root"),
    )
    if isinstance(factor_cache_spec, str) and not factor_cache_spec.lower().endswith(".json"):
        factor_cache_spec = {"type": "zarr", "root_dir": factor_cache_spec}
    grm_factor_cache = _normalize_factor_cache(factor_cache_spec, base_dir)
    if grm_factor_cache is not None:
        fit_options["grm_factor_cache"] = _factor_cache_object(grm_factor_cache, project_root)

    res = mt_gp.fit_multi_trait_gp(
        pheno_df=pheno_df,
        gid_col=str(gid_col),
        trait_col=str(trait_col),
        y_col=str(y_col),
        geno_kernels=geno_kernels,
        geno_ids=geno_ids,
        train_idx=train_idx,
        test_idx=test_idx,
        **fit_options,
    )
    _write_public_result_outputs(
        res,
        out_dir,
        {
            "python": sys.executable,
            "command": "fit-multi-trait-spec",
            "spec": str(Path(args.spec).resolve()),
            "kernel_names": list(geno_kernels.keys()),
            "has_grm_factor_cache": grm_factor_cache is not None,
        },
    )
    return 0


def _fit_multi_trait_met_spec(args) -> int:
    spec, base_dir, project_root, pheno_df, geno_kernels, geno_ids, train_idx, test_idx, out_dir = _read_common_public_spec(args)
    mt_gp_met = _load_project_module(project_root, "mt_gp_met")
    cols = _spec_get(spec, "columns", default={}) or {}
    gid_col = _coalesce(_spec_get(spec, "gid_col"), cols.get("gid"), cols.get("gid_col"), default="gid")
    env_col = _coalesce(_spec_get(spec, "env_col"), cols.get("env"), cols.get("env_col"), default="env")
    trait_col = _coalesce(_spec_get(spec, "trait_col"), cols.get("trait"), cols.get("trait_col"), default="trait")
    y_col = _coalesce(_spec_get(spec, "y_col"), cols.get("y"), cols.get("y_col"), default="y")
    fit_options = dict(_spec_get(spec, "fit", "fit_args", "options", "model", default={}) or {})

    env_covariates = None
    env_covariates_csv = _resolve_path(
        _spec_get(spec, "env_covariates_csv", "inputs.env_covariates_csv", "environment_covariates_csv", "env_covariates.csv", "environment_covariates.csv"),
        base_dir,
        label="env_covariates_csv",
    )
    if env_covariates_csv is not None:
        env_covariates = pd.read_csv(env_covariates_csv)
    env_similarity_spec = _coalesce(
        _spec_get(spec, "env_similarity", "environment_similarity", "inputs.env_similarity"),
        {
            "bin": _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin"),
            "meta": _spec_get(spec, "env_similarity_meta", "inputs.env_similarity_meta"),
        }
        if _spec_get(spec, "env_similarity_bin", "inputs.env_similarity_bin") is not None
        else None,
    )
    env_similarity = _read_matrix_from_spec(env_similarity_spec, base_dir, "env_similarity")
    env_ids = _read_lines_from_spec(
        _spec_get(spec, "env_ids", "environment_ids", "inputs.env_ids", "inputs.environment_ids", "env_similarity.env_ids", "env_similarity.ids", "environment_similarity.env_ids"),
        base_dir,
        "env_ids",
        required=False,
    )
    factor_cache_spec = _coalesce(
        _spec_get(spec, "factor_cache_spec", "grm_factor_cache", "gp_factor_cache", "inputs.factor_cache_spec"),
        _spec_get(spec, "grm_factor_cache_root", "gp_factor_cache_root"),
    )
    if isinstance(factor_cache_spec, str) and not factor_cache_spec.lower().endswith(".json"):
        factor_cache_spec = {"type": "zarr", "root_dir": factor_cache_spec}
    grm_factor_cache = _normalize_factor_cache(factor_cache_spec, base_dir)
    if grm_factor_cache is not None:
        fit_options["grm_factor_cache"] = _factor_cache_object(grm_factor_cache, project_root)

    res = mt_gp_met.fit_multi_trait_gp_met(
        pheno_df=pheno_df,
        gid_col=str(gid_col),
        env_col=str(env_col),
        trait_col=str(trait_col),
        y_col=str(y_col),
        geno_kernels=geno_kernels,
        geno_ids=geno_ids,
        train_idx=train_idx,
        test_idx=test_idx,
        env_similarity=env_similarity,
        env_ids=env_ids,
        env_covariates=env_covariates,
        **fit_options,
    )
    _write_public_result_outputs(
        res,
        out_dir,
        {
            "python": sys.executable,
            "command": "fit-multi-trait-met-spec",
            "spec": str(Path(args.spec).resolve()),
            "kernel_names": list(geno_kernels.keys()),
            "has_env_covariates": env_covariates is not None,
            "has_env_similarity": env_similarity is not None,
            "has_grm_factor_cache": grm_factor_cache is not None,
        },
    )
    return 0


def _write_hybrid_fit_outputs(fit: Dict[str, Any], out_dir: Path, meta: Dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "fit.json").write_text(json.dumps(fit, indent=2, default=_json_default), encoding="utf-8")
    meta = dict(meta)
    meta["fit_keys"] = sorted(list(fit.keys()))
    if isinstance(fit.get("device_info"), dict):
        meta["device_info"] = fit.get("device_info")
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=_json_default), encoding="utf-8")


def _read_hybrid_kernel(spec: Dict[str, Any], base_dir: Path, name: str, *, required: bool = True) -> Optional[np.ndarray]:
    aliases = [
        f"kernels.{name}",
        f"inputs.kernels.{name}",
        f"K_{name}",
        f"inputs.K_{name}",
        name,
    ]
    kernel_spec = _spec_get(spec, *aliases)
    if kernel_spec is None:
        if required:
            raise ValueError(f"kernels.{name} is required")
        return None
    return _read_matrix_from_spec(kernel_spec, base_dir, f"kernels.{name}")


def _fit_hybrid_additive_spec(args) -> int:
    spec_path = Path(args.spec).resolve()
    spec = _read_json_file(str(spec_path))
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root"), default=None)
    if project_root is None:
        raise ValueError("project_root must be supplied in the spec or with --project-root")
    project_root = _resolve_path(project_root, base_dir, required=True, label="project_root")
    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))

    hybrid_gp = _load_project_module(project_root, "hybrid_gp")
    y = _read_vector_from_spec(_spec_get(spec, "y", "inputs.y"), base_dir, "y", required=True, dtype=np.float64)
    X = _read_matrix_from_spec(_spec_get(spec, "X", "inputs.X", "fixed_effects", "inputs.fixed_effects"), base_dir, "X")
    if X is None:
        raise ValueError("X is required")
    train_idx = _read_idx_from_spec(_spec_get(spec, "train_idx", "inputs.train_idx"), base_dir, "train_idx", required=True)
    k_gxe = _read_hybrid_kernel(spec, base_dir, "gxe", required=False)

    fit = hybrid_gp.fit_hybrid_additive_gp(
        y=y,
        X=X,
        K_female=_read_hybrid_kernel(spec, base_dir, "female", required=True),
        K_male=_read_hybrid_kernel(spec, base_dir, "male", required=True),
        K_sca=_read_hybrid_kernel(spec, base_dir, "sca", required=True),
        K_gxe=k_gxe,
        train_idx=train_idx,
        lambda_value=float(_coalesce(_spec_get(spec, "lambda_value", "lambda"), default=0.1)),
        device=_spec_get(spec, "device", "fit.device"),
        dtype=str(_coalesce(_spec_get(spec, "dtype", "fit.dtype"), default="float64")),
        seed=int(_coalesce(_spec_get(spec, "seed", "fit.seed"), default=12345)),
        min_gpu_work_units=_spec_get(spec, "min_gpu_work_units", "fit.min_gpu_work_units"),
    )
    _write_hybrid_fit_outputs(
        fit,
        out_dir,
        {
            "python": sys.executable,
            "command": "fit-hybrid-additive-spec",
            "spec": str(spec_path),
        },
    )
    return 0


def _fit_hybrid_multi_trait_spec(args) -> int:
    spec_path = Path(args.spec).resolve()
    spec = _read_json_file(str(spec_path))
    base_dir = spec_path.parent
    project_root = _coalesce(args.project_root, _spec_get(spec, "project_root"), default=None)
    if project_root is None:
        raise ValueError("project_root must be supplied in the spec or with --project-root")
    project_root = _resolve_path(project_root, base_dir, required=True, label="project_root")
    out_dir_value = _coalesce(args.out_dir, _spec_get(spec, "out_dir", "output_dir", "outputs.out_dir", "outputs.output_dir"))
    out_dir = Path(_resolve_path(out_dir_value, base_dir, required=True, label="out_dir"))

    hybrid_gp = _load_project_module(project_root, "hybrid_gp")
    k_gxe = _read_hybrid_kernel(spec, base_dir, "gxe", required=False)
    x_joint = _read_matrix_from_spec(_spec_get(spec, "X_joint", "inputs.X_joint"), base_dir, "X_joint")
    trait_cor = _read_matrix_from_spec(_spec_get(spec, "trait_cor", "inputs.trait_cor"), base_dir, "trait_cor")
    if x_joint is None:
        raise ValueError("X_joint is required")
    if trait_cor is None:
        raise ValueError("trait_cor is required")
    fit = hybrid_gp.fit_hybrid_additive_multi_trait_gp(
        y_std=_read_vector_from_spec(_spec_get(spec, "y_std", "y_std_vec", "inputs.y_std"), base_dir, "y_std", required=True, dtype=np.float64),
        X_joint=x_joint,
        K_female=_read_hybrid_kernel(spec, base_dir, "female", required=True),
        K_male=_read_hybrid_kernel(spec, base_dir, "male", required=True),
        K_sca=_read_hybrid_kernel(spec, base_dir, "sca", required=True),
        K_gxe=k_gxe,
        trait_cor=trait_cor,
        gxe_trait_cor=_read_matrix_from_spec(_spec_get(spec, "gxe_trait_cor", "inputs.gxe_trait_cor"), base_dir, "gxe_trait_cor"),
        all_row=_read_idx_from_spec(_spec_get(spec, "all_row", "inputs.all_row"), base_dir, "all_row", required=True),
        all_trait=_read_idx_from_spec(_spec_get(spec, "all_trait", "inputs.all_trait"), base_dir, "all_trait", required=True),
        train_obs_idx=_read_idx_from_spec(_spec_get(spec, "train_obs_idx", "inputs.train_obs_idx"), base_dir, "train_obs_idx", required=True),
        lambda_value=float(_coalesce(_spec_get(spec, "lambda_value", "lambda"), default=0.1)),
        device=_spec_get(spec, "device", "fit.device"),
        dtype=str(_coalesce(_spec_get(spec, "dtype", "fit.dtype"), default="float64")),
        seed=int(_coalesce(_spec_get(spec, "seed", "fit.seed"), default=12345)),
        min_gpu_work_units=_spec_get(spec, "min_gpu_work_units", "fit.min_gpu_work_units"),
    )
    _write_hybrid_fit_outputs(
        fit,
        out_dir,
        {
            "python": sys.executable,
            "command": "fit-hybrid-multi-trait-spec",
            "spec": str(spec_path),
        },
    )
    return 0


def _materialize_factor_cache(args) -> int:
    project_root = _resolve_path(args.project_root, Path.cwd(), required=True, label="project_root")
    fw = _load_framework(project_root)
    cache_spec = _normalize_factor_cache(_read_json_file(args.factor_cache_json), Path(args.factor_cache_json).resolve().parent)
    row_idx = _read_idx(args.row_idx)
    K = fw.materialize_geno_kernel_from_factor_cache(
        grm_factor_cache=cache_spec,
        row_index=row_idx,
        dtype=str(args.dtype),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_matrix_bin(K, str(out_dir / "kernel.bin"), str(out_dir / "kernel_shape.txt"))
    meta = {
        "python": sys.executable,
        "command": "materialize-factor-cache",
        "n_rows": int(np.asarray(K).shape[0]),
        "n_cols": int(np.asarray(K).shape[1]),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=_json_default), encoding="utf-8")
    return 0


def _fit_predict_with_framework(args, fw) -> int:
    pheno_df = pd.read_csv(args.pheno_csv)
    geno_kernel = _read_matrix_bin(args.kernel_bin, args.kernel_meta)
    geno_ids = _read_lines(args.geno_ids)
    train_idx = _read_idx(args.train_idx)
    test_idx = _read_idx(args.test_idx)

    kw: Dict[str, Any] = {
        "include_components": ["g"],
        "fixed_effects": [],
        "method": args.method,
        "backend": args.backend,
        "prediction_output": args.prediction_output,
        "output_level": args.output_level,
        "point_predictions_only": args.output_level == "predict_only",
        "return_se": args.output_level in ("predict_with_se", "full_vc"),
        "compute_ai_se": args.output_level == "full_vc",
        "standardize": "global",
        "varcomp_mode": args.varcomp_mode,
        "dtype": "float64",
        "seed": int(args.seed),
        "train_idx": train_idx,
        "test_idx": test_idx,
    }
    if args.gp_iters is not None and int(args.gp_iters) > 0:
        kw["gp_iters"] = int(args.gp_iters)
    if args.gp_lr is not None and float(args.gp_lr) > 0:
        kw["gp_lr"] = float(args.gp_lr)
    if args.method == "gp_icm_fa":
        kw["env_structure"] = "fa"
        kw["fa_rank"] = int(args.fa_rank)
    if args.grm_factor_cache_root:
        kw["grm_factor_cache"] = {
            "type": "zarr",
            "root": args.grm_factor_cache_root,
            "paths": {},
            "shapes": {},
            "dtypes": {},
        }

    res = fw.fit_mixed_model(
        pheno_df=pheno_df,
        gid_col=args.gid_col,
        env_col=args.env_col,
        y_col=args.y_col,
        geno_kernels={"G": geno_kernel},
        geno_ids=geno_ids,
        **kw,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_df(res.get("result", {}).get("predictions"), str(out_dir / "predictions.csv"))
    _write_df(res.get("var_components"), str(out_dir / "var_components.csv"))
    _write_df(res.get("var_components_summary"), str(out_dir / "var_components_summary.csv"))
    meta = {
        "python": sys.executable,
        "method": args.method,
        "backend": args.backend,
        "output_level": args.output_level,
        "prediction_output": args.prediction_output,
        "diagnostics_keys": sorted(list((res.get("diagnostics") or {}).keys())),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return 0


def _worker_dispatch(cmd: str, ns, fw) -> int:
    if cmd == "fit-predict":
        return _fit_predict_with_framework(ns, fw)
    if cmd == "fit-predict-spec":
        return _fit_predict_spec(ns)
    if cmd == "fit-predict-cv-spec":
        return _fit_predict_cv_spec(ns, fw=fw)
    if cmd == "batch-fit-candidates-spec":
        return _fit_batch_candidates_spec(ns, fw=fw)
    if cmd == "fit-multi-trait-spec":
        return _fit_multi_trait_spec(ns)
    if cmd == "fit-multi-trait-met-spec":
        return _fit_multi_trait_met_spec(ns)
    if cmd == "fit-hybrid-additive-spec":
        return _fit_hybrid_additive_spec(ns)
    if cmd == "fit-hybrid-multi-trait-spec":
        return _fit_hybrid_multi_trait_spec(ns)
    raise ValueError(f"unknown command: {cmd}")


def _serve(args) -> int:
    fw = _load_framework(args.project_root)
    print(json.dumps({"status": "ready"}), flush=True)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        cmd = req.get("cmd")
        if cmd == "shutdown":
            print(json.dumps({"status": "ok"}), flush=True)
            return 0
        if cmd == "ping":
            print(json.dumps({"status": "ok"}), flush=True)
            continue
        try:
            class Request:
                pass
            ns = Request()
            for k, v in req.get("args", {}).items():
                setattr(ns, k.replace("-", "_"), v)
            _worker_dispatch(cmd, ns, fw)
            print(json.dumps({"status": "ok", "out_dir": getattr(ns, "out_dir", None)}), flush=True)
        except Exception as exc:
            print(json.dumps({"status": "error", "message": f"{type(exc).__name__}: {exc}"}), flush=True)


def _serve_socket(args) -> int:
    fw = _load_framework(args.project_root)
    host = args.host
    port = int(args.port)
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as srv:
        srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        srv.bind((host, port))
        srv.listen(5)
        while True:
            conn, _addr = srv.accept()
            with conn:
                conn.settimeout(float(getattr(args, "client_timeout", 30.0)))
                data = b""
                while True:
                    try:
                        chunk = conn.recv(65536)
                    except socket.timeout:
                        break
                    if not chunk:
                        break
                    data += chunk
                    if b"\n" in chunk:
                        break
                if not data:
                    continue
                req = json.loads(data.decode("utf-8").strip())
                cmd = req.get("cmd")
                if cmd == "shutdown":
                    conn.sendall((json.dumps({"status": "ok"}) + "\n").encode("utf-8"))
                    return 0
                if cmd == "ping":
                    conn.sendall((json.dumps({"status": "ok"}) + "\n").encode("utf-8"))
                    continue
                try:
                    class Request:
                        pass
                    ns = Request()
                    for k, v in req.get("args", {}).items():
                        setattr(ns, k.replace("-", "_"), v)
                    _worker_dispatch(cmd, ns, fw)
                    conn.sendall((json.dumps({"status": "ok", "out_dir": getattr(ns, "out_dir", None)}) + "\n").encode("utf-8"))
                except Exception as exc:
                    conn.sendall((json.dumps({"status": "error", "message": f"{type(exc).__name__}: {exc}"}) + "\n").encode("utf-8"))


def _build_factor_cache(args) -> int:
    project_root = args.project_root
    script = os.path.join(project_root, "grm_pivoted_cholesky_to_zarr_streaming.py")
    if not os.path.exists(script):
        raise FileNotFoundError(f"Factor-cache script not found: {script}")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grm = _read_matrix_bin(args.matrix_bin, args.matrix_meta)
    grm_npy = out_dir / "grm.npy"
    np.save(grm_npy, grm)
    cmd = [
        sys.executable,
        script,
        "--grm",
        str(grm_npy),
        "--out",
        str(out_dir),
        "--eps_trace",
        str(args.eps_trace),
        "--max_rank",
        str(int(args.max_rank)),
        "--store_dtype",
        str(args.store_dtype),
    ]
    if args.overwrite:
        cmd.append("--overwrite")
    proc = _run(cmd)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        return int(proc.returncode)
    return 0


def _print_setup_info(info: Dict[str, Any]) -> int:
    for key, value in info.items():
        print(f"{key}={value}")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_setup = sub.add_parser("setup-deps")
    p_setup.add_argument("--project-root", default=None)
    p_setup.add_argument("--prefer-gpu", action=argparse.BooleanOptionalAction, default=True)
    p_setup.add_argument("--cuda", default="auto")
    p_setup.add_argument("--torch-version", default=None)
    p_setup.add_argument("--index-url", default=None)
    p_setup.add_argument(
        "--install-scikit-sparse",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Attempt cross-OS install of scikit-sparse (CHOLMOD) for the "
            "MME REML engine. Default ON; pass --no-install-scikit-sparse "
            "to skip. Failures are non-fatal -- MME falls back to "
            "scipy.sparse.linalg.splu."
        ),
    )

    p_fit = sub.add_parser("fit-predict")
    p_fit.add_argument("--project-root", required=True)
    p_fit.add_argument("--pheno-csv", required=True)
    p_fit.add_argument("--kernel-bin", required=True)
    p_fit.add_argument("--kernel-meta", required=True)
    p_fit.add_argument("--geno-ids", required=True)
    p_fit.add_argument("--train-idx", required=True)
    p_fit.add_argument("--test-idx", required=True)
    p_fit.add_argument("--out-dir", required=True)
    p_fit.add_argument("--gid-col", default="GID")
    p_fit.add_argument("--env-col", default="Env")
    p_fit.add_argument("--y-col", default="Trait")
    p_fit.add_argument("--method", required=True)
    p_fit.add_argument("--backend", default="auto")
    p_fit.add_argument("--output-level", default="predict_only")
    p_fit.add_argument("--prediction-output", default="all")
    p_fit.add_argument("--varcomp-mode", default="mom")
    p_fit.add_argument("--fa-rank", type=int, default=1)
    p_fit.add_argument("--seed", type=int, default=12345)
    p_fit.add_argument("--gp-iters", type=int, default=None)
    p_fit.add_argument("--gp-lr", type=float, default=None)
    p_fit.add_argument("--grm-factor-cache-root", default=None)

    p_fit_spec = sub.add_parser("fit-predict-spec")
    p_fit_spec.add_argument("--spec", required=True)
    p_fit_spec.add_argument("--project-root", default=None)
    p_fit_spec.add_argument("--out-dir", default=None)

    p_fit_cv_spec = sub.add_parser("fit-predict-cv-spec")
    p_fit_cv_spec.add_argument("--spec", required=True)
    p_fit_cv_spec.add_argument("--project-root", default=None)
    p_fit_cv_spec.add_argument("--out-dir", default=None)

    p_batch_candidates_spec = sub.add_parser("batch-fit-candidates-spec")
    p_batch_candidates_spec.add_argument("--spec", required=True)
    p_batch_candidates_spec.add_argument("--project-root", default=None)
    p_batch_candidates_spec.add_argument("--out-dir", default=None)

    p_mt_spec = sub.add_parser("fit-multi-trait-spec")
    p_mt_spec.add_argument("--spec", required=True)
    p_mt_spec.add_argument("--project-root", default=None)
    p_mt_spec.add_argument("--out-dir", default=None)

    p_mt_met_spec = sub.add_parser("fit-multi-trait-met-spec")
    p_mt_met_spec.add_argument("--spec", required=True)
    p_mt_met_spec.add_argument("--project-root", default=None)
    p_mt_met_spec.add_argument("--out-dir", default=None)

    p_hybrid_spec = sub.add_parser("fit-hybrid-additive-spec")
    p_hybrid_spec.add_argument("--spec", required=True)
    p_hybrid_spec.add_argument("--project-root", default=None)
    p_hybrid_spec.add_argument("--out-dir", default=None)

    p_hybrid_mt_spec = sub.add_parser("fit-hybrid-multi-trait-spec")
    p_hybrid_mt_spec.add_argument("--spec", required=True)
    p_hybrid_mt_spec.add_argument("--project-root", default=None)
    p_hybrid_mt_spec.add_argument("--out-dir", default=None)

    p_materialize = sub.add_parser("materialize-factor-cache")
    p_materialize.add_argument("--project-root", required=True)
    p_materialize.add_argument("--factor-cache-json", required=True)
    p_materialize.add_argument("--row-idx", required=True)
    p_materialize.add_argument("--out-dir", required=True)
    p_materialize.add_argument("--dtype", default="float64")

    p_cache = sub.add_parser("build-factor-cache")
    p_cache.add_argument("--project-root", required=True)
    p_cache.add_argument("--matrix-bin", required=True)
    p_cache.add_argument("--matrix-meta", required=True)
    p_cache.add_argument("--out-dir", required=True)
    p_cache.add_argument("--eps-trace", type=float, default=1e-4)
    p_cache.add_argument("--max-rank", type=int, default=4096)
    p_cache.add_argument("--store-dtype", default="float32")
    p_cache.add_argument("--overwrite", action="store_true")

    p_serve = sub.add_parser("serve")
    p_serve.add_argument("--project-root", required=True)

    p_serve_socket = sub.add_parser("serve-socket")
    p_serve_socket.add_argument("--project-root", required=True)
    p_serve_socket.add_argument("--host", default="127.0.0.1")
    p_serve_socket.add_argument("--port", type=int, required=True)
    p_serve_socket.add_argument("--client-timeout", type=float, default=30.0)

    args = parser.parse_args(argv)

    try:
      if args.cmd == "setup-deps":
          return _print_setup_info(setup_deps(
              project_root=args.project_root,
              prefer_gpu=bool(args.prefer_gpu),
              cuda=args.cuda,
              torch_version=args.torch_version,
              index_url=args.index_url,
              install_scikit_sparse=bool(args.install_scikit_sparse),
          ))
      if args.cmd == "fit-predict":
          return _fit_predict(args)
      if args.cmd == "fit-predict-spec":
          return _fit_predict_spec(args)
      if args.cmd == "fit-predict-cv-spec":
          return _fit_predict_cv_spec(args)
      if args.cmd == "batch-fit-candidates-spec":
          return _fit_batch_candidates_spec(args)
      if args.cmd == "fit-multi-trait-spec":
          return _fit_multi_trait_spec(args)
      if args.cmd == "fit-multi-trait-met-spec":
          return _fit_multi_trait_met_spec(args)
      if args.cmd == "fit-hybrid-additive-spec":
          return _fit_hybrid_additive_spec(args)
      if args.cmd == "fit-hybrid-multi-trait-spec":
          return _fit_hybrid_multi_trait_spec(args)
      if args.cmd == "materialize-factor-cache":
          return _materialize_factor_cache(args)
      if args.cmd == "build-factor-cache":
          return _build_factor_cache(args)
      if args.cmd == "serve":
          return _serve(args)
      if args.cmd == "serve-socket":
          return _serve_socket(args)
      raise ValueError(f"Unknown command: {args.cmd}")
    except Exception as exc:
      print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
      raise


if __name__ == "__main__":
    raise SystemExit(main())
