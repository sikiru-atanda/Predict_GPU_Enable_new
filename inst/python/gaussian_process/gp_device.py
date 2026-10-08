"""Shared device selection for PredictProR Gaussian-process backends."""
from __future__ import annotations

import os
import subprocess
import warnings
from typing import Any, Dict, Optional

import torch


def gpu_usage_percent() -> Optional[float]:
    try:
        proc = subprocess.run(
            ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            text=True,
        )
        vals = [float(x.strip()) for x in (proc.stdout or "").splitlines() if x.strip()]
        vals = [v for v in vals if v == v]
        if not vals:
            return None
        return sum(vals) / len(vals)
    except Exception:
        return None


def gpu_busy() -> bool:
    try:
        threshold = float(os.environ.get("PREDICTPRO_GPU_BUSY_THRESHOLD", "85"))
    except Exception:
        threshold = 85.0
    usage = gpu_usage_percent()
    return usage is not None and usage >= max(0.0, min(100.0, threshold))


def _truthy_env(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in ("1", "true", "yes", "on")


def _float_env(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return float(default)


def pick_torch_device(
    device: Optional[Any] = None,
    *,
    work_units: Optional[float] = None,
    min_gpu_work_units: Optional[float] = None,
) -> torch.device:
    """Return the requested torch device.

    ``auto`` keeps the public GPU-first behavior, but avoids sending very small
    dense GP solves to CUDA when launch/transfer overhead is expected to exceed
    compute savings. Set ``PREDICTPRO_GP_GPU_MIN_WORK_UNITS=0`` to disable this
    size gate, or pass ``device='cuda'`` to force CUDA for a call.
    """
    requested = str(device or os.environ.get("PREDICTPRO_GP_DEVICE", "auto")).strip().lower()
    if requested in ("", "auto", "gpu"):
        if torch.cuda.is_available():
            if work_units is not None and not _truthy_env("PREDICTPRO_GP_IGNORE_SMALL_WORKLOAD"):
                threshold = (
                    float(min_gpu_work_units)
                    if min_gpu_work_units is not None
                    else _float_env("PREDICTPRO_GP_GPU_MIN_WORK_UNITS", 1000000.0)
                )
                try:
                    wu = float(work_units)
                except Exception:
                    wu = threshold
                if threshold > 0 and wu < threshold:
                    return torch.device("cpu")
            if gpu_busy() and not _truthy_env("PREDICTPRO_GP_IGNORE_GPU_BUSY"):
                return torch.device("cpu")
            return torch.device("cuda")
        return torch.device("cpu")

    try:
        dev = torch.device(requested)
    except Exception:
        warnings.warn(
            f"Unknown GP torch device {requested!r}; using CPU.",
            RuntimeWarning,
            stacklevel=2,
        )
        return torch.device("cpu")

    if dev.type == "cuda" and not torch.cuda.is_available():
        warnings.warn(
            "CUDA was requested for GP models, but torch.cuda.is_available() is FALSE; using CPU.",
            RuntimeWarning,
            stacklevel=2,
        )
        return torch.device("cpu")
    return dev


def torch_to_numpy(x: Any):
    if torch.is_tensor(x):
        return x.detach().to("cpu").numpy()
    return x


def torch_device_info(device: Optional[Any] = None) -> Dict[str, Any]:
    dev = pick_torch_device(device)
    usage = gpu_usage_percent()
    out: Dict[str, Any] = {
        "device": str(dev),
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
        "gpu_name": None,
        "gpu_usage": usage,
        "gpu_busy": bool(gpu_busy()),
        "torch_version": getattr(torch, "__version__", None),
        "cuda_compiled": getattr(getattr(torch, "version", None), "cuda", None),
    }
    if out["cuda_available"] and out["cuda_device_count"] > 0:
        out["gpu_name"] = str(torch.cuda.get_device_name(0))
    return out
