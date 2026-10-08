"""Torch backend for PredictProR hybrid GP additive kernels."""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

import numpy as np
import torch

try:
    from gp_device import pick_torch_device, torch_device_info
except Exception:  # pragma: no cover - fallback for isolated import tests.
    def pick_torch_device(device=None, *, work_units=None, min_gpu_work_units=None):
        if str(device or "auto").lower() in ("cuda", "gpu") and torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")

    def torch_device_info(device=None):
        dev = pick_torch_device(device)
        return {
            "device": str(dev),
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
            "gpu_name": str(torch.cuda.get_device_name(0)) if torch.cuda.is_available() else None,
            "torch_version": getattr(torch, "__version__", None),
            "cuda_compiled": getattr(getattr(torch, "version", None), "cuda", None),
        }


def _as_numpy(x: Any, dtype=np.float64) -> np.ndarray:
    arr = np.array(x, dtype=dtype, order="C", copy=True)
    if not arr.flags["C_CONTIGUOUS"]:
        arr = np.ascontiguousarray(arr)
    return arr


def _tensor(x: Any, *, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    return torch.as_tensor(_as_numpy(x), device=device, dtype=dtype)


def _solve_spd(A: torch.Tensor, B: torch.Tensor, jitter: float = 1e-8, max_tries: int = 6) -> torch.Tensor:
    n = int(A.shape[0])
    eye = torch.eye(n, device=A.device, dtype=A.dtype)
    B2 = B if B.ndim > 1 else B.reshape(-1, 1)
    for i in range(int(max_tries)):
        add = float(jitter) * (10.0 ** i)
        try:
            L = torch.linalg.cholesky(A + add * eye)
            return torch.cholesky_solve(B2, L)
        except RuntimeError:
            continue
    try:
        return torch.linalg.solve(A + float(jitter) * (10.0 ** max_tries) * eye, B2)
    except RuntimeError:
        return torch.linalg.pinv(A) @ B2


def _to_cpu_numpy(x: torch.Tensor) -> np.ndarray:
    return x.detach().to("cpu").numpy()


def fit_hybrid_additive_gp(
    *,
    y: Sequence[float],
    X: Any,
    K_female: Any,
    K_male: Any,
    K_sca: Any,
    train_idx: Sequence[int],
    lambda_value: float,
    K_gxe: Optional[Any] = None,
    device: Optional[Any] = None,
    dtype: str = "float64",
    seed: int = 12345,
    min_gpu_work_units: Optional[float] = None,
) -> Dict[str, Any]:
    """Fit additive female GCA + male GCA + SCA GP/KRR equations.

    The algebra mirrors ``gp_hybrid_gp_fit_core`` in R. Indices are zero-based.
    """
    torch.manual_seed(int(seed))
    dtype_t = torch.float32 if str(dtype).lower() == "float32" else torch.float64
    y_np = _as_numpy(y, dtype=np.float64).reshape(-1)
    train_np = np.asarray(train_idx, dtype=np.int64).reshape(-1)
    n_all = int(y_np.shape[0])
    n_train = int(train_np.shape[0])
    work_units = float(max(n_train, 1) ** 3 + max(n_all, 1) * max(n_train, 1))
    dev = pick_torch_device(
        device=device,
        work_units=work_units,
        min_gpu_work_units=min_gpu_work_units,
    )

    Kf = _tensor(K_female, device=dev, dtype=dtype_t)
    Km = _tensor(K_male, device=dev, dtype=dtype_t)
    Ks = _tensor(K_sca, device=dev, dtype=dtype_t)
    Kgxe = _tensor(K_gxe if K_gxe is not None else np.zeros_like(_as_numpy(K_sca)), device=dev, dtype=dtype_t)
    X_t = _tensor(X, device=dev, dtype=dtype_t)
    train_t = torch.as_tensor(train_np, device=dev, dtype=torch.long)
    y_train_np = y_np[train_np]
    y_mean = float(np.nanmean(y_train_np))
    y_sd = float(np.nanstd(y_train_np, ddof=1)) if n_train > 1 else 1.0
    if not np.isfinite(y_sd) or y_sd <= 0:
        y_sd = 1.0
    y_std_np = np.full_like(y_np, np.nan, dtype=np.float64)
    finite_y = np.isfinite(y_np)
    y_std_np[finite_y] = (y_np[finite_y] - y_mean) / y_sd
    y_train = torch.as_tensor(y_std_np[train_np], device=dev, dtype=dtype_t).reshape(-1, 1)

    lam = float(lambda_value)
    if not np.isfinite(lam) or lam <= 0:
        lam = 0.1

    K_total = Kf + Km + Ks + Kgxe
    K_train = K_total.index_select(0, train_t).index_select(1, train_t)
    X_train = X_t.index_select(0, train_t)
    V = K_train + lam * torch.eye(n_train, device=dev, dtype=dtype_t)

    Vinv_y = _solve_spd(V, y_train)
    Vinv_X = _solve_spd(V, X_train)
    Xt_Vinv_X = X_train.T @ Vinv_X
    Xt_Vinv_y = X_train.T @ Vinv_y
    beta = _solve_spd(
        Xt_Vinv_X + 1e-8 * torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
        Xt_Vinv_y,
    )
    alpha = Vinv_y - Vinv_X @ beta

    K_all_train = K_total.index_select(1, train_t)
    fixed_std = (X_t @ beta).reshape(-1)
    female_std = (Kf.index_select(1, train_t) @ alpha).reshape(-1)
    male_std = (Km.index_select(1, train_t) @ alpha).reshape(-1)
    sca_std = (Ks.index_select(1, train_t) @ alpha).reshape(-1)
    gxe_std = (Kgxe.index_select(1, train_t) @ alpha).reshape(-1)
    total_std = fixed_std + female_std + male_std + sca_std + gxe_std

    Vinv_K = _solve_spd(V, K_all_train.T)
    k_vinv_k = torch.sum(K_all_train * Vinv_K.T, dim=1)
    base_var = torch.clamp(torch.diagonal(K_total) - k_vinv_k, min=0.0)
    M_inv = _solve_spd(
        Xt_Vinv_X + 1e-8 * torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
        torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
    )
    x_delta = X_t - K_all_train @ Vinv_X
    fixed_var = torch.sum((x_delta @ M_inv) * x_delta, dim=1)
    latent_var_std = torch.clamp(base_var + fixed_var, min=0.0)
    observed_var_std = torch.clamp(latent_var_std + lam, min=0.0)

    scale = float(y_sd)
    info = torch_device_info(dev)
    info["backend"] = "torch"
    info["dtype"] = str(dtype_t).replace("torch.", "")
    info["work_units"] = work_units

    return {
        "lambda": lam,
        "y_mean": y_mean,
        "y_sd": y_sd,
        "beta": _to_cpu_numpy(beta.reshape(-1)),
        "alpha": _to_cpu_numpy(alpha.reshape(-1)),
        "fixed_component": y_mean + scale * _to_cpu_numpy(fixed_std),
        "female_component": scale * _to_cpu_numpy(female_std),
        "male_component": scale * _to_cpu_numpy(male_std),
        "sca_component": scale * _to_cpu_numpy(sca_std),
        "gxe_component": scale * _to_cpu_numpy(gxe_std),
        "prediction": y_mean + scale * _to_cpu_numpy(total_std),
        "latent_var": (scale ** 2) * _to_cpu_numpy(latent_var_std),
        "observed_var": (scale ** 2) * _to_cpu_numpy(observed_var_std),
        "train_idx": train_np + 1,
        "device_info": info,
    }


def fit_hybrid_additive_multi_trait_gp(
    *,
    y_std: Sequence[float],
    X_joint: Any,
    K_female: Any,
    K_male: Any,
    K_sca: Any,
    trait_cor: Any,
    all_row: Sequence[int],
    all_trait: Sequence[int],
    train_obs_idx: Sequence[int],
    lambda_value: float,
    K_gxe: Optional[Any] = None,
    gxe_trait_cor: Optional[Any] = None,
    device: Optional[Any] = None,
    dtype: str = "float64",
    seed: int = 12345,
    min_gpu_work_units: Optional[float] = None,
) -> Dict[str, Any]:
    """Fit joint cross-trait hybrid GP equations.

    The kernel for stacked observation ``(row, trait)`` pairs is
    ``trait_cor[trait_i, trait_j] * K_component[row_i, row_j]`` for each
    additive hybrid component. Indices are zero-based.
    """
    torch.manual_seed(int(seed))
    dtype_t = torch.float32 if str(dtype).lower() == "float32" else torch.float64
    y_np = _as_numpy(y_std, dtype=np.float64).reshape(-1)
    all_row_np = np.asarray(all_row, dtype=np.int64).reshape(-1)
    all_trait_np = np.asarray(all_trait, dtype=np.int64).reshape(-1)
    train_np = np.asarray(train_obs_idx, dtype=np.int64).reshape(-1)
    n_all = int(y_np.shape[0])
    n_train = int(train_np.shape[0])
    n_fixed = int(np.asarray(X_joint).shape[1])
    work_units = float(max(n_train, 1) ** 3 + max(n_all, 1) * max(n_train, 1) + max(n_fixed, 1) ** 3)
    dev = pick_torch_device(
        device=device,
        work_units=work_units,
        min_gpu_work_units=min_gpu_work_units,
    )

    Kf = _tensor(K_female, device=dev, dtype=dtype_t)
    Km = _tensor(K_male, device=dev, dtype=dtype_t)
    Ks = _tensor(K_sca, device=dev, dtype=dtype_t)
    Kgxe = _tensor(K_gxe if K_gxe is not None else np.zeros_like(_as_numpy(K_sca)), device=dev, dtype=dtype_t)
    C = _tensor(trait_cor, device=dev, dtype=dtype_t)
    Cgxe = _tensor(gxe_trait_cor if gxe_trait_cor is not None else trait_cor, device=dev, dtype=dtype_t)
    X_t = _tensor(X_joint, device=dev, dtype=dtype_t)
    y_t = torch.as_tensor(y_np, device=dev, dtype=dtype_t).reshape(-1, 1)
    all_row_t = torch.as_tensor(all_row_np, device=dev, dtype=torch.long)
    all_trait_t = torch.as_tensor(all_trait_np, device=dev, dtype=torch.long)
    train_t = torch.as_tensor(train_np, device=dev, dtype=torch.long)
    train_row_t = all_row_t.index_select(0, train_t)
    train_trait_t = all_trait_t.index_select(0, train_t)

    lam = float(lambda_value)
    if not np.isfinite(lam) or lam <= 0:
        lam = 0.1

    K_base = Kf + Km + Ks
    trait_train = C.index_select(0, train_trait_t).index_select(1, train_trait_t)
    hybrid_train = K_base.index_select(0, train_row_t).index_select(1, train_row_t)
    gxe_trait_train = Cgxe.index_select(0, train_trait_t).index_select(1, train_trait_t)
    gxe_train = Kgxe.index_select(0, train_row_t).index_select(1, train_row_t)
    K_train = trait_train * hybrid_train + gxe_trait_train * gxe_train
    X_train = X_t.index_select(0, train_t)
    y_train = y_t.index_select(0, train_t)
    V = K_train + lam * torch.eye(n_train, device=dev, dtype=dtype_t)

    Vinv_y = _solve_spd(V, y_train)
    Vinv_X = _solve_spd(V, X_train)
    Xt_Vinv_X = X_train.T @ Vinv_X
    Xt_Vinv_y = X_train.T @ Vinv_y
    beta = _solve_spd(
        Xt_Vinv_X + 1e-8 * torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
        Xt_Vinv_y,
    )
    alpha = Vinv_y - Vinv_X @ beta

    def joint_cross(K: torch.Tensor) -> torch.Tensor:
        trait_cross = C.index_select(0, all_trait_t).index_select(1, train_trait_t)
        hybrid_cross = K.index_select(0, all_row_t).index_select(1, train_row_t)
        return trait_cross * hybrid_cross

    def joint_cross_with_trait(K: torch.Tensor, C_use: torch.Tensor) -> torch.Tensor:
        trait_cross = C_use.index_select(0, all_trait_t).index_select(1, train_trait_t)
        hybrid_cross = K.index_select(0, all_row_t).index_select(1, train_row_t)
        return trait_cross * hybrid_cross

    K_base_all_train = joint_cross(K_base)
    K_gxe_all_train = joint_cross_with_trait(Kgxe, Cgxe)
    K_total_all_train = K_base_all_train + K_gxe_all_train
    Kf_all_train = joint_cross(Kf)
    Km_all_train = joint_cross(Km)
    Ks_all_train = joint_cross(Ks)

    fixed_std = (X_t @ beta).reshape(-1)
    female_std = (Kf_all_train @ alpha).reshape(-1)
    male_std = (Km_all_train @ alpha).reshape(-1)
    sca_std = (Ks_all_train @ alpha).reshape(-1)
    gxe_std = (K_gxe_all_train @ alpha).reshape(-1)
    prediction_std = fixed_std + female_std + male_std + sca_std + gxe_std

    Vinv_K = _solve_spd(V, K_total_all_train.T)
    k_vinv_k = torch.sum(K_total_all_train * Vinv_K.T, dim=1)
    diag_total = (
        torch.diagonal(K_base).index_select(0, all_row_t) * C[all_trait_t, all_trait_t]
        + torch.diagonal(Kgxe).index_select(0, all_row_t) * Cgxe[all_trait_t, all_trait_t]
    )
    base_var = torch.clamp(diag_total - k_vinv_k, min=0.0)
    M_inv = _solve_spd(
        Xt_Vinv_X + 1e-8 * torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
        torch.eye(int(Xt_Vinv_X.shape[0]), device=dev, dtype=dtype_t),
    )
    x_delta = X_t - K_total_all_train @ Vinv_X
    fixed_var = torch.sum((x_delta @ M_inv) * x_delta, dim=1)
    latent_var_std = torch.clamp(base_var + fixed_var, min=0.0)
    observed_var_std = torch.clamp(latent_var_std + lam, min=0.0)

    info = torch_device_info(dev)
    info["backend"] = "torch"
    info["dtype"] = str(dtype_t).replace("torch.", "")
    info["work_units"] = work_units
    info["model"] = "joint_hybrid_multi_trait_gp"

    return {
        "lambda": lam,
        "beta": _to_cpu_numpy(beta.reshape(-1)),
        "alpha": _to_cpu_numpy(alpha.reshape(-1)),
        "fixed_std": _to_cpu_numpy(fixed_std),
        "female_std": _to_cpu_numpy(female_std),
        "male_std": _to_cpu_numpy(male_std),
        "sca_std": _to_cpu_numpy(sca_std),
        "gxe_std": _to_cpu_numpy(gxe_std),
        "prediction_std": _to_cpu_numpy(prediction_std),
        "latent_var_std": _to_cpu_numpy(latent_var_std),
        "observed_var_std": _to_cpu_numpy(observed_var_std),
        "train_obs_idx": train_np + 1,
        "device_info": info,
    }
