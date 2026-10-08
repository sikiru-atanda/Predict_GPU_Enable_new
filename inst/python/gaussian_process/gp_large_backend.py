
from __future__ import annotations

"""
Large-scale operator backend for multi-GRM + GxE with env similarity, designed for:
  - n_obs >= 50k
  - n_env ~ 300
  - R (num GRMs) >= 5
  - GRMs provided dense (no markers)
  - Low-rank GRM factors cached to disk (Zarr or memmap)
  - CV-friendly: point_predictions_only=True skips REML logdet / AI / SEs

This module provides:
  - GRMFactorCacheZarr / GRMFactorCacheMemmap
  - LowRankGRMOperator: v -> Z_g Phi Phi^T Z_g^T v
  - LowRankGxEOperator: v -> Z_ge (Phi Phi^T ⊗ Sigma_e) Z_ge^T v
  - SumOperator and DiagResidualOperator
  - PCGSolver (preconditioned conjugate gradient)
  - fit_predict_operator_backend(...) for point predictions at scale

Important:
  - This backend assumes you already factorized each dense GRM G_r into a low-rank root Phi_r (Ngeno x m_r)
    and stored Phi_r to disk once. Recommended factorization: adaptive pivoted Cholesky (APC) with a trace
    tolerance eps_trace (e.g. 1e-4 to 1e-5), then store Phi in Zarr with row-friendly chunking.

  - This backend never allocates Ngeno-length vectors at runtime. It gathers only rows of Phi corresponding to
    genotypes used in a CV fold, and uses scatter_add + feature matmuls.

  - For full REML (logdet) and AI/SEs at this scale, you will need SLQ logdet and Hutchinson traces; the operators
    here are designed for that (matvec-only).
"""

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import os
import warnings
import numpy as np
import torch

try:
    import zarr
except Exception:
    zarr = None


# ----------------------------
# Utilities
# ----------------------------

def _as_torch(x, device, dtype):
    if torch.is_tensor(x):
        return x.to(device=device, dtype=dtype)
    return torch.as_tensor(x, device=device, dtype=dtype)

def _unique_inverse(idx: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    return torch.unique(idx, sorted=True, return_inverse=True)

def _scatter_add_1d(values: torch.Tensor, inv: torch.Tensor, out_len: int) -> torch.Tensor:
    out = torch.zeros(out_len, device=values.device, dtype=values.dtype)
    out.scatter_add_(0, inv, values)
    return out

def _scatter_add_cells(values: torch.Tensor, inv: torch.Tensor, ei: torch.Tensor, k_env: int, ng_used: int) -> torch.Tensor:
    flat = torch.zeros(ng_used * k_env, device=values.device, dtype=values.dtype)
    cell = inv * k_env + ei
    flat.scatter_add_(0, cell, values)
    return flat.view(ng_used, k_env)


# ----------------------------
# Factor caches
# ----------------------------

class GRMFactorCacheBase:
    def __init__(self, num_grms: int):
        self.num_grms = int(num_grms)

    def ranks(self) -> List[int]:
        raise NotImplementedError

    def n_geno(self) -> int:
        raise NotImplementedError

    def get_rows(self, r: int, row_index: torch.Tensor, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
        raise NotImplementedError


class GRMFactorCacheZarr(GRMFactorCacheBase):
    """
    Zarr-backed cache.

    Layout:
      root_dir/
        grm0/Phi  (Ngeno, m0)
        grm1/Phi  (Ngeno, m1)
        ...
    """
    def __init__(self, root_dir: str, num_grms: int):
        super().__init__(num_grms=num_grms)
        if zarr is None:
            raise ImportError("zarr not available; install zarr or use GRMFactorCacheMemmap.")
        self.root_dir = root_dir
        self._phi = []
        for r in range(self.num_grms):
            g = zarr.open_group(os.path.join(root_dir, f"grm{r}"), mode="r")
            self._phi.append(g["Phi"])
        self._n_geno = int(self._phi[0].shape[0])

    def ranks(self) -> List[int]:
        return [int(a.shape[1]) for a in self._phi]

    def n_geno(self) -> int:
        return self._n_geno

    def get_rows(self, r: int, row_index: torch.Tensor, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
        row_cpu = row_index.detach().to("cpu").numpy().astype(np.int64, copy=False)
        Phi_np = np.asarray(self._phi[r].oindex[row_cpu, :])
        return torch.as_tensor(Phi_np, device=device, dtype=dtype)


class GRMFactorCacheMemmap(GRMFactorCacheBase):
    """
    Memmap-backed cache using raw binary files.

    Provide a list of file paths and shapes. Each file stores Phi_r in row-major.
    """
    def __init__(self, paths: Sequence[str], shapes: Sequence[Tuple[int,int]], dtypes: Sequence[np.dtype]):
        super().__init__(num_grms=len(paths))
        assert len(paths) == len(shapes) == len(dtypes)
        # Specs arrive as JSON from R, where a shape can serialise as nested
        # singletons (e.g. list(as.list(c(12, 6))) -> [[12], [6]]); numpy
        # memmap needs a flat tuple of ints.
        shapes = [tuple(int(v) for v in np.asarray(sh, dtype=np.int64).reshape(-1)) for sh in shapes]
        dtypes = [np.asarray(dt).reshape(-1)[0] if isinstance(dt, (list, tuple)) else dt for dt in dtypes]
        for sh in shapes:
            if len(sh) != 2:
                raise ValueError(f"GRM factor cache shape must have two dimensions; got {sh}")
        self.paths = list(paths)
        self.shapes = list(shapes)
        self.dtypes = list(dtypes)
        self._n_geno = int(shapes[0][0])
        self._mm = []
        for p, sh, dt in zip(paths, shapes, dtypes):
            self._mm.append(np.memmap(p, mode="r", dtype=dt, shape=sh, order="C"))

    def ranks(self) -> List[int]:
        return [int(s[1]) for s in self.shapes]

    def n_geno(self) -> int:
        return self._n_geno

    def get_rows(self, r: int, row_index: torch.Tensor, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
        row_cpu = row_index.detach().to("cpu").numpy().astype(np.int64, copy=False)
        Phi_np = np.asarray(self._mm[r][row_cpu, :])
        return torch.as_tensor(Phi_np, device=device, dtype=dtype)


# ----------------------------
# Operators
# ----------------------------

@dataclass
class FoldIndexCache:
    g_unique: torch.Tensor  # (ng_used,) int64
    inv: torch.Tensor       # (n_fold,) int64
    ng_used: int


class LowRankGRMOperator(torch.nn.Module):
    def __init__(self, cache: GRMFactorCacheBase, gi: torch.Tensor, weights: torch.Tensor):
        super().__init__()
        self.cache = cache
        self.register_buffer("gi", gi.long())
        self.weights = torch.nn.Parameter(weights)

    def build_fold_cache(self, obs_index: Optional[torch.Tensor]=None) -> FoldIndexCache:
        gi = self.gi if obs_index is None else self.gi.index_select(0, obs_index.long())
        g_unique, inv = _unique_inverse(gi)
        return FoldIndexCache(g_unique=g_unique, inv=inv, ng_used=int(g_unique.numel()))

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, dtype_compute: torch.dtype=torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        t_used = _scatter_add_1d(v, fold.inv, fold.ng_used)
        out = torch.zeros_like(v)

        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_used = self.cache.get_rows(r, fold.g_unique, device=device, dtype=dtype_compute)
            u = Phi_used.transpose(0, 1) @ t_used
            s_used = Phi_used @ u
            out = out + w * s_used.index_select(0, fold.inv)

        return out


    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        fold_train = test_state["fold_train"]
        gi_test = test_state["gi_test"].long()
        t_used = _scatter_add_1d(alpha_train, fold_train.inv, fold_train.ng_used)
        g_unique_tr = fold_train.g_unique
        pred = torch.zeros((gi_test.numel(),), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_tr = self.cache.get_rows(r, g_unique_tr, device=device, dtype=dtype_compute)
            u = Phi_tr.transpose(0, 1) @ t_used
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            pred = pred + w * (Phi_te @ u)
        return pred

    @torch.no_grad()
    def cross_cov_apply_transpose(self, z_test: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        # K(tr,te) z_te = Σ_r w_r * Z_tr Phi_tr_unique (Phi_te^T z_te)
        device = z_test.device
        z_test = z_test.to(dtype=dtype_compute)
        fold_train = test_state["fold_train"]
        gi_test = test_state["gi_test"].long()
        g_unique_tr = fold_train.g_unique
        n_train = int(fold_train.inv.numel())
        out = torch.zeros((n_train,), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            u = Phi_te.transpose(0, 1) @ z_test
            Phi_tr = self.cache.get_rows(r, g_unique_tr, device=device, dtype=dtype_compute)
            s_used = Phi_tr @ u
            out = out + w * s_used.index_select(0, fold_train.inv)
        return out

    @torch.no_grad()
    def self_var_diag(self, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        # K(te,te)_ii = Σ_r w_r * ||Phi[gi_te_i]||^2
        device = test_state["gi_test"].device
        gi_test = test_state["gi_test"].long()
        diag = torch.zeros((gi_test.numel(),), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            diag = diag + w * (Phi_te * Phi_te).sum(dim=1)
        return diag

class LowRankGxEOperator(torch.nn.Module):
    def __init__(self, cache: GRMFactorCacheBase, gi: torch.Tensor, ei: torch.Tensor, Sigma_e: torch.Tensor, weights: torch.Tensor):
        super().__init__()
        self.cache = cache
        self.register_buffer("gi", gi.long())
        self.register_buffer("ei", ei.long())
        self.weights = torch.nn.Parameter(weights)
        self.register_buffer("Sigma_e", Sigma_e)

    def build_fold_cache(self, obs_index: Optional[torch.Tensor]=None) -> Tuple[FoldIndexCache, torch.Tensor]:
        if obs_index is None:
            gi = self.gi
            ei = self.ei
        else:
            obs_index = obs_index.long()
            gi = self.gi.index_select(0, obs_index)
            ei = self.ei.index_select(0, obs_index)
        g_unique, inv = _unique_inverse(gi)
        return FoldIndexCache(g_unique=g_unique, inv=inv, ng_used=int(g_unique.numel())), ei

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, ei_fold: torch.Tensor,
               dtype_compute: torch.dtype=torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        Sigma_e = self.Sigma_e.to(device=device, dtype=dtype_compute)
        k_env = int(Sigma_e.shape[0])

        T = _scatter_add_cells(v, fold.inv, ei_fold.long(), k_env=k_env, ng_used=fold.ng_used)
        out = torch.zeros_like(v)

        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_used = self.cache.get_rows(r, fold.g_unique, device=device, dtype=dtype_compute)
            A = Phi_used.transpose(0, 1) @ T
            B = Phi_used @ A
            U = B @ Sigma_e.transpose(0, 1)
            out = out + w * U[fold.inv, ei_fold.long()]

        return out



    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        fold_train = test_state["fold_train"]
        ei_train = test_state["ei_train"].long()
        gi_test = test_state["gi_test"].long()
        ei_test = test_state["ei_test"].long()
        k_env = int(self.Sigma_e.shape[0])
        T = _scatter_add_cells(alpha_train, fold_train.inv, ei_train, k_env=k_env, ng_used=fold_train.ng_used)
        Sigma_e_t = self.Sigma_e.to(device=device, dtype=dtype_compute)
        g_unique_tr = fold_train.g_unique
        pred = torch.zeros((gi_test.numel(),), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_tr = self.cache.get_rows(r, g_unique_tr, device=device, dtype=dtype_compute)
            A = Phi_tr.transpose(0, 1) @ T
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            B_te = Phi_te @ A
            U_te = B_te @ Sigma_e_t.transpose(0, 1)
            pred = pred + w * U_te[torch.arange(gi_test.numel(), device=device), ei_test]
        return pred

    @torch.no_grad()
    def cross_cov_apply_transpose(self, z_test: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        # K(tr,te) z_te:
        #   Scatter z_te into (n_test, k_env) by ei_te (one-hot weighted), then B_back = U_back @ Sigma_e,
        #   A_back = Phi_te^T @ B_back, T_back = Phi_tr @ A_back, out[j] = T_back[inv[j], ei_train[j]]
        device = z_test.device
        z_test = z_test.to(dtype=dtype_compute)
        fold_train = test_state["fold_train"]
        ei_train = test_state["ei_train"].long()
        gi_test = test_state["gi_test"].long()
        ei_test = test_state["ei_test"].long()
        Sigma_e_t = self.Sigma_e.to(device=device, dtype=dtype_compute)
        k_env = int(Sigma_e_t.shape[0])
        n_test = int(gi_test.numel())
        n_train = int(fold_train.inv.numel())
        g_unique_tr = fold_train.g_unique
        # U_back[i,e] = z_te[i] * δ(e == ei_test[i])
        U_back = torch.zeros((n_test, k_env), device=device, dtype=dtype_compute)
        U_back[torch.arange(n_test, device=device), ei_test] = z_test
        # B_back[i,:] = U_back[i,:] @ Sigma_e^T  (adjoint of right-mul by Sigma_e^T)
        B_back = U_back @ Sigma_e_t
        out = torch.zeros((n_train,), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            # A_back[r_idx, e] = Σ_i Phi_te[i, r_idx] * B_back[i, e]
            A_back = Phi_te.transpose(0, 1) @ B_back  # (m, k_env)
            Phi_tr = self.cache.get_rows(r, g_unique_tr, device=device, dtype=dtype_compute)
            T_back = Phi_tr @ A_back  # (ng_used, k_env)
            # out[j] = T_back[inv[j], ei_train[j]]
            gathered = T_back.view(-1)[fold_train.inv * k_env + ei_train]
            out = out + w * gathered
        return out

    @torch.no_grad()
    def self_var_diag(self, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        # K(te,te)_ii = Σ_r w_r * ||Phi[gi_te_i]||^2 * Sigma_e[ei_te_i, ei_te_i]
        device = test_state["gi_test"].device
        gi_test = test_state["gi_test"].long()
        ei_test = test_state["ei_test"].long()
        Sigma_e_t = self.Sigma_e.to(device=device, dtype=dtype_compute)
        Se_diag = torch.diag(Sigma_e_t)
        se_te = Se_diag.index_select(0, ei_test)
        diag = torch.zeros((gi_test.numel(),), device=device, dtype=dtype_compute)
        for r in range(self.cache.num_grms):
            w = self.weights[r].to(dtype_compute)
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            diag = diag + w * (Phi_te * Phi_te).sum(dim=1) * se_te
        return diag

class EnvKernelOperator(torch.nn.Module):
    """
    Environment (non-genetic) covariance in observation space:
        v -> w_e * Z_e S_e Z_e^T v
    where Z_e maps observations to environments (one-hot by ei), and S_e is (k,k).

    This costs O(n + k^2) per matvec:
      - aggregate v by env: t[e] = sum_{i:ei_i=e} v_i
      - u = S_e @ t
      - return u[ei]
    """
    def __init__(self, ei: torch.Tensor, S_e: torch.Tensor, weight: Union[float, torch.Tensor] = 1.0):
        super().__init__()
        self.register_buffer("ei", ei.long())
        self.register_buffer("S_e", S_e)
        w = float(weight) if not torch.is_tensor(weight) else None
        if w is None:
            self.weight = torch.nn.Parameter(weight.reshape(()))
        else:
            self.weight = torch.nn.Parameter(torch.as_tensor(w, dtype=torch.float32))

    def build_fold_cache(self, obs_index: Optional[torch.Tensor]=None) -> torch.Tensor:
        return self.ei if obs_index is None else self.ei.index_select(0, obs_index.long())

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, ei_fold: torch.Tensor, dtype_compute: torch.dtype=torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        S_e = self.S_e.to(device=device, dtype=dtype_compute)
        k_env = int(S_e.shape[0])
        t = torch.zeros((k_env,), device=device, dtype=dtype_compute)
        t.scatter_add_(0, ei_fold.long(), v)
        u = S_e @ t
        return self.weight.to(dtype_compute) * u.index_select(0, ei_fold.long())



    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        ei_train = test_state["ei_train"].long()
        ei_test = test_state["ei_test"].long()
        S_e_t = self.S_e.to(device=device, dtype=dtype_compute)
        k_env = int(S_e_t.shape[0])
        t_env = torch.zeros((k_env,), device=device, dtype=dtype_compute)
        t_env.scatter_add_(0, ei_train, alpha_train)
        u_env = S_e_t @ t_env
        return self.weight.to(dtype_compute) * u_env.index_select(0, ei_test)

    @torch.no_grad()
    def cross_cov_apply_transpose(self, z_test: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = z_test.device
        z_test = z_test.to(dtype=dtype_compute)
        ei_train = test_state["ei_train"].long()
        ei_test = test_state["ei_test"].long()
        S_e_t = self.S_e.to(device=device, dtype=dtype_compute)
        k_env = int(S_e_t.shape[0])
        t_te = torch.zeros((k_env,), device=device, dtype=dtype_compute)
        t_te.scatter_add_(0, ei_test, z_test)
        # S_e assumed symmetric; adjoint is S_e @ t_te
        u_tr = S_e_t @ t_te
        return self.weight.to(dtype_compute) * u_tr.index_select(0, ei_train)

    @torch.no_grad()
    def self_var_diag(self, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = test_state["ei_test"].device
        ei_test = test_state["ei_test"].long()
        S_e_t = self.S_e.to(device=device, dtype=dtype_compute)
        Se_diag = torch.diag(S_e_t)
        return self.weight.to(dtype_compute) * Se_diag.index_select(0, ei_test)

# ----------------------------
# Random term interface (extensible)
# ----------------------------

class RandomTerm(torch.nn.Module):
    """
    Base class for random terms in operator backend.

    A term must implement:
      - build_fold_cache(obs_index=None) -> any cache object needed for matvec
      - matvec(v, cache, dtype_compute) -> term contribution in observation space
    """
    term_type: str = "base"

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None):
        raise NotImplementedError

    def matvec(self, v: torch.Tensor, cache, dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        raise NotImplementedError

    def cross_cov_apply(self, alpha_train: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        raise NotImplementedError(f"cross_cov_apply not implemented for term_type={self.term_type}")

    def describe(self) -> Dict[str, object]:
        return {"term_type": getattr(self, "term_type", "base")}


class IIDGroupOperator(RandomTerm):
    """
    IID random intercept by group (Z I Z^T).
    """
    term_type = "iid"

    def __init__(self, group_index: torch.Tensor, weight: Union[float, torch.Tensor] = 1.0):
        super().__init__()
        self.register_buffer("group_index", group_index.long())
        if torch.is_tensor(weight):
            self.weight = torch.nn.Parameter(weight.reshape(()))
        else:
            self.weight = torch.nn.Parameter(torch.as_tensor(float(weight), dtype=torch.float32))

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None) -> torch.Tensor:
        return self.group_index if obs_index is None else self.group_index.index_select(0, obs_index.long())

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, cache: Optional[torch.Tensor], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        g = (self.group_index if cache is None else cache).long()
        n_levels = int(g.max().item() + 1)
        t = torch.zeros((n_levels,), device=device, dtype=dtype_compute)
        t.scatter_add_(0, g, v)
        return self.weight.to(dtype_compute) * t.index_select(0, g)

    def describe(self) -> Dict[str, object]:
        return {"term_type": "iid", "n_levels": int(self.group_index.max().item() + 1)}


    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        g_tr = test_state["group_train"].long()
        g_te = test_state["group_test"].long()
        n_levels = int(torch.max(g_tr).item() + 1)
        t_level = torch.zeros((n_levels,), device=device, dtype=dtype_compute)
        t_level.scatter_add_(0, g_tr, alpha_train)
        return self.weight.to(dtype_compute) * t_level.index_select(0, g_te)

    @torch.no_grad()
    def cross_cov_apply_transpose(self, z_test: torch.Tensor, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = z_test.device
        z_test = z_test.to(dtype=dtype_compute)
        g_tr = test_state["group_train"].long()
        g_te = test_state["group_test"].long()
        n_levels = int(max(torch.max(g_tr).item(), torch.max(g_te).item()) + 1)
        t_level = torch.zeros((n_levels,), device=device, dtype=dtype_compute)
        t_level.scatter_add_(0, g_te, z_test)
        return self.weight.to(dtype_compute) * t_level.index_select(0, g_tr)

    @torch.no_grad()
    def self_var_diag(self, test_state: Dict[str, object], dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = test_state["group_test"].device
        g_te = test_state["group_test"].long()
        return torch.full((g_te.numel(),), float(self.weight.to(dtype_compute).item()), device=device, dtype=dtype_compute)

def build_random_terms_operator_backend(
    *,
    gi: torch.Tensor,
    ei: torch.Tensor,
    cache: "GRMFactorCacheBase",
    Sigma_e: torch.Tensor,
    S_e: Optional[torch.Tensor],
    w_g: torch.Tensor,
    w_ge: torch.Tensor,
    w_e: Union[float, torch.Tensor],
    r_diag: Union[float, torch.Tensor],
    random_terms: Optional[List[Dict[str, object]]] = None,
) -> List[torch.nn.Module]:
    """
    Build operator list from random_terms specification.

    Core model structure is preserved by default:
      g + ge + env_kernel + residual

    Explicit random_terms are additive extensions and/or overrides:
      - iid terms are added on top of the core model
      - grm/gxe/env_kernel/resid entries override the corresponding core component
    """
    if isinstance(r_diag, (float, int)):
        r_diag_obs = torch.full((gi.numel(),), float(r_diag), device=gi.device, dtype=Sigma_e.dtype)
    else:
        r_diag_obs = r_diag

    default_terms: Dict[str, torch.nn.Module] = {
        "grm": LowRankGRMOperator(cache, gi, weights=w_g),
        "gxe": LowRankGxEOperator(cache, gi, ei, Sigma_e=Sigma_e, weights=w_ge),
        "env_kernel": EnvKernelOperator(
            ei,
            S_e if S_e is not None else torch.eye(int(Sigma_e.shape[0]), device=gi.device, dtype=Sigma_e.dtype),
            weight=w_e,
        ),
        "resid": DiagResidualOperator(diag=r_diag_obs),
    }

    if random_terms is None:
        return [
            default_terms["grm"],
            default_terms["gxe"],
            default_terms["env_kernel"],
            default_terms["resid"],
        ]

    extra_ops: List[torch.nn.Module] = []

    for term in random_terms:
        ttype = str(term.get("type", "")).lower()

        if ttype in ("grm", "g", "geno"):
            weights = term.get("weights", w_g)
            default_terms["grm"] = LowRankGRMOperator(
                cache, gi, weights=_as_torch(weights, gi.device, Sigma_e.dtype)
            )

        elif ttype in ("gxe", "ge", "gxe_term", "g_by_e"):
            weights = term.get("weights", w_ge)
            Sigma = term.get("Sigma_e", Sigma_e)
            default_terms["gxe"] = LowRankGxEOperator(
                cache, gi, ei,
                Sigma_e=_as_torch(Sigma, gi.device, Sigma_e.dtype),
                weights=_as_torch(weights, gi.device, Sigma_e.dtype),
            )

        elif ttype in ("env_kernel", "env", "e"):
            Se = term.get("S_e", S_e)
            if Se is None:
                Se = torch.eye(int(Sigma_e.shape[0]), device=gi.device, dtype=Sigma_e.dtype)
            w = term.get("weight", w_e)
            default_terms["env_kernel"] = EnvKernelOperator(
                ei, _as_torch(Se, gi.device, Sigma_e.dtype), weight=w
            )

        elif ttype in ("iid", "random_intercept"):
            grp = term.get("group_index", None)
            if grp is None:
                raise ValueError("IID term requires 'group_index'")
            w = term.get("weight", 1.0)
            op = IIDGroupOperator(_as_torch(grp, gi.device, torch.int64), weight=w)
            op.term_name = str(term.get("name", f"iid_{len(extra_ops)}"))
            extra_ops.append(op)

        elif ttype in ("resid", "residual"):
            diag = term.get("diag", r_diag)
            if isinstance(diag, (float, int)):
                diag_t = torch.full((gi.numel(),), float(diag), device=gi.device, dtype=Sigma_e.dtype)
            else:
                diag_t = _as_torch(diag, gi.device, Sigma_e.dtype)
            default_terms["resid"] = DiagResidualOperator(diag=diag_t)

        else:
            raise ValueError(f"Unsupported random term type: {ttype}")

    return [
        default_terms["grm"],
        default_terms["gxe"],
        default_terms["env_kernel"],
        *extra_ops,
        default_terms["resid"],
    ]
class DiagResidualOperator(torch.nn.Module):
    def __init__(self, diag: torch.Tensor):
        super().__init__()
        self.diag = torch.nn.Parameter(diag)

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, dtype_compute=torch.float32) -> torch.Tensor:
        return self.diag.to(device=v.device, dtype=dtype_compute) * v.to(dtype_compute)


class SumOperator(torch.nn.Module):
    def __init__(self, ops: Sequence[torch.nn.Module]):
        super().__init__()
        self.ops = torch.nn.ModuleList(list(ops))

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, context: dict, dtype_compute=torch.float32) -> torch.Tensor:
        out = torch.zeros_like(v, dtype=dtype_compute)
        for op in self.ops:
            if isinstance(op, LowRankGRMOperator):
                out = out + op.matvec(v, context["fold_g"], dtype_compute=dtype_compute)
            elif isinstance(op, LowRankGxEOperator):
                out = out + op.matvec(v, context["fold_ge"], context["ei_fold"], dtype_compute=dtype_compute)
            elif isinstance(op, EnvKernelOperator):
                out = out + op.matvec(v, context["ei_fold_env"], dtype_compute=dtype_compute)
            elif isinstance(op, IIDGroupOperator):
                # Use fold-specific group index if provided; otherwise operator uses its stored index.
                g_fold = context.get('iid_group_index', None)
                out = out + op.matvec(v, g_fold, dtype_compute=dtype_compute)
            elif isinstance(op, DiagResidualOperator):
                out = out + op.matvec(v, dtype_compute=dtype_compute)
            else:
                raise TypeError(f"Unknown operator type: {type(op)}")
        return out


# ----------------------------
# PCG Solver
# ----------------------------

@dataclass
class PCGResult:
    x: torch.Tensor
    iters: int
    converged: bool
    rel_resid: float


def pcg_solve(
    A_mv: Callable[[torch.Tensor], torch.Tensor],
    b: torch.Tensor,
    M_inv_mv: Optional[Callable[[torch.Tensor], torch.Tensor]] = None,
    x0: Optional[torch.Tensor] = None,
    tol: float = 1e-5,
    max_iter: int = 500,
) -> PCGResult:
    device = b.device
    dtype = b.dtype
    x = torch.zeros_like(b) if x0 is None else x0.clone()
    r = b - A_mv(x)
    b_norm = torch.norm(b).clamp_min(torch.tensor(1e-30, device=device, dtype=dtype))
    rel = (torch.norm(r) / b_norm).item()
    if rel < tol:
        return PCGResult(x=x, iters=0, converged=True, rel_resid=rel)

    z = r if M_inv_mv is None else M_inv_mv(r)
    p = z.clone()
    rz_old = torch.dot(r, z)

    eps = torch.tensor(1e-30, device=device, dtype=dtype)

    for k in range(1, max_iter + 1):
        Ap = A_mv(p)
        alpha = rz_old / torch.dot(p, Ap).clamp_min(eps)
        x = x + alpha * p
        r = r - alpha * Ap
        rel = (torch.norm(r) / b_norm).item()
        if rel < tol:
            return PCGResult(x=x, iters=k, converged=True, rel_resid=rel)
        z = r if M_inv_mv is None else M_inv_mv(r)
        rz_new = torch.dot(r, z)
        beta = rz_new / rz_old.clamp_min(eps)
        p = z + beta * p
        rz_old = rz_new

    return PCGResult(x=x, iters=max_iter, converged=False, rel_resid=rel)


# ----------------------------
# Randomized Nyström preconditioner for V
# ----------------------------

@dataclass
class NystromFactors:
    U: torch.Tensor          # (n, ℓ), approximately orthonormal
    Lambda: torch.Tensor     # (ℓ,), non-negative, approximating top eigenvalues of V
    mu: float                # shift used in P = U Λ U^T + μ (I - U U^T)


@torch.no_grad()
def build_nystrom_preconditioner(
    V_mv: Callable[[torch.Tensor], torch.Tensor],
    n: int,
    rank: int,
    device: torch.device,
    dtype: torch.dtype,
    seed: int = 0,
    mu: Optional[float] = None,
) -> Tuple[NystromFactors, Callable[[torch.Tensor], torch.Tensor]]:
    """
    Randomized Nyström approximation of V, intended as a PCG preconditioner.

    Constructs V̂ = U diag(Λ) U^T + μ (I − U U^T) with
        V̂^{-1} = U diag(1/(Λ+μ)) U^T + (1/μ) (I − U U^T)
    (equivalently, (Λ+μ)^{-1} on the range of U and μ^{-1} on its complement).
    This captures the top ℓ eigenpairs of V so PCG effectively only has to
    handle the small-eigenvalue tail.

    Cost: ℓ matvecs with V plus O(n ℓ² + ℓ³) linear algebra. Apply is O(n ℓ).
    """
    rank = max(1, min(int(rank), int(n)))
    gen = torch.Generator(device=device).manual_seed(int(seed))
    G = torch.randn((n, rank), generator=gen, device=device, dtype=dtype)
    Omega, _ = torch.linalg.qr(G)  # (n, ℓ), orthonormal columns

    # Y = V Omega (ℓ matvecs)
    Y = torch.empty((n, rank), device=device, dtype=dtype)
    for j in range(rank):
        Y[:, j] = V_mv(Omega[:, j])

    eps_mach = torch.finfo(dtype).eps
    y_norm = Y.norm().item()
    nu = float(eps_mach) * max(y_norm, 1.0)
    Y_nu = Y + nu * Omega

    C = Omega.transpose(0, 1) @ Y_nu
    C = 0.5 * (C + C.transpose(0, 1))

    try:
        L = torch.linalg.cholesky(C)
    except Exception:
        diag_max = float(torch.diag(C).abs().max().item())
        jitter = max(1e-6 * diag_max, 1e-10)
        L = torch.linalg.cholesky(
            C + jitter * torch.eye(rank, device=device, dtype=dtype)
        )

    F = torch.linalg.solve_triangular(L, Y_nu.transpose(0, 1), upper=False).transpose(0, 1)
    U, s, _ = torch.linalg.svd(F, full_matrices=False)
    Lambda = (s * s - nu).clamp_min(0.0)

    if mu is None:
        lam_max = float(Lambda.max().item()) if Lambda.numel() > 0 else 1.0
        mu_val = float(Lambda.min().clamp_min(1e-6 * max(lam_max, 1.0)).item())
    else:
        mu_val = float(mu)
    mu_t = torch.as_tensor(mu_val, device=device, dtype=dtype)

    inv_lmu = 1.0 / (Lambda + mu_t)  # (ℓ,)
    inv_mu = 1.0 / mu_t

    def P_inv_mv(v: torch.Tensor) -> torch.Tensor:
        v_c = v.to(dtype=dtype)
        ut_v = U.transpose(0, 1) @ v_c
        return inv_mu * (v_c - U @ ut_v) + U @ (inv_lmu * ut_v)

    return NystromFactors(U=U, Lambda=Lambda, mu=mu_val), P_inv_mv


# ----------------------------
# Lanczos low-rank V^{-1} for Hutchinson variance reduction
# ----------------------------

@dataclass
class LanczosLowRankInv:
    """Low-rank approximation of V^{-1} from Lanczos Ritz pairs.

    V_low := V_vecs diag(1/theta) V_vecs^T is the best rank-k approximation of
    V^{-1} restricted to the Krylov subspace spanned by the starting vector.
    """
    V_vecs: torch.Tensor  # (n, k) orthonormal-ish columns
    theta: torch.Tensor   # (k,) positive Ritz values approximating some eigvals of V


@torch.no_grad()
def build_lanczos_lowrank_inv(
    V_mv: Callable[[torch.Tensor], torch.Tensor],
    n: int,
    k_iter: int,
    device: torch.device,
    dtype: torch.dtype,
    seed: int = 0,
    start_vec: Optional[torch.Tensor] = None,
    reorth: bool = True,
    theta_floor: float = 1e-8,
    breakdown_tol: float = 1e-10,
) -> LanczosLowRankInv:
    """
    Symmetric Lanczos on V (PSD + jitter) with full reorthogonalization.

    Builds a Krylov subspace Q_k = span{q_0, V q_0, V^2 q_0, ..., V^{k-1} q_0}
    and extracts Ritz pairs (Q_ritz, theta) from the tridiagonal T_k.
    The returned LanczosLowRankInv encodes the best rank-k approximation of
    V^{-1} within Q_k, which is V_low = Q_ritz diag(1/theta) Q_ritz^T.

    For SE variance reduction in K V^{-1} K^T, pass start_vec = K^T z
    (Rademacher projected into training space) so the Krylov subspace aligns
    with K's column range — exactly where PEV mass lives.

    Cost: k_iter matvecs with V, plus O(n k_iter) for reorthogonalization and
    O(k_iter^3) for the small tridiagonal eigenproblem.
    """
    k_iter = max(1, min(int(k_iter), int(n)))
    gen = torch.Generator(device=device).manual_seed(int(seed))

    if start_vec is None:
        q = torch.randn((n,), generator=gen, device=device, dtype=dtype)
    else:
        q = start_vec.detach().to(device=device, dtype=dtype).reshape(-1)
        if q.numel() != n:
            raise ValueError(f"start_vec has {q.numel()} entries, expected {n}")
    q_norm = q.norm().item()
    if not np.isfinite(q_norm) or q_norm < 1e-20:
        q = torch.randn((n,), generator=gen, device=device, dtype=dtype)
        q_norm = q.norm().item()
    q = q / q_norm

    Q = torch.empty((n, k_iter), device=device, dtype=dtype)
    alpha = torch.empty((k_iter,), device=device, dtype=dtype)
    beta = torch.empty((k_iter,), device=device, dtype=dtype)

    Q[:, 0] = q
    q_prev = torch.zeros_like(q)
    beta_prev = torch.zeros((), device=device, dtype=dtype)
    k_actual = k_iter
    for j in range(k_iter):
        Vq = V_mv(Q[:, j])
        a = torch.dot(Q[:, j], Vq)
        alpha[j] = a
        r = Vq - a * Q[:, j] - beta_prev * q_prev
        if reorth:
            # two-pass classical Gram-Schmidt against all prior Lanczos vectors
            for _ in range(2):
                proj = Q[:, : j + 1].transpose(0, 1) @ r
                r = r - Q[:, : j + 1] @ proj
        b = r.norm()
        beta[j] = b
        if float(b.item()) < breakdown_tol:
            k_actual = j + 1
            break
        if j + 1 < k_iter:
            q_prev = Q[:, j]
            Q[:, j + 1] = r / b
            beta_prev = b

    k = int(k_actual)
    T = torch.zeros((k, k), device=device, dtype=torch.float64)
    a64 = alpha[:k].to(torch.float64)
    T[torch.arange(k), torch.arange(k)] = a64
    if k > 1:
        off = beta[: k - 1].to(torch.float64)
        T[torch.arange(1, k), torch.arange(k - 1)] = off
        T[torch.arange(k - 1), torch.arange(1, k)] = off
    theta_all, W = torch.linalg.eigh(T)
    theta_all = theta_all.to(dtype)
    W = W.to(dtype)

    # Ritz vectors in the original n-dim space
    Q_ritz = Q[:, :k] @ W  # (n, k)

    # Keep only Ritz values safely above the floor — division by tiny theta would
    # blow up the low-rank apply and defeat variance reduction.
    mask = theta_all > float(theta_floor)
    theta_keep = theta_all[mask]
    Q_keep = Q_ritz[:, mask]

    return LanczosLowRankInv(V_vecs=Q_keep, theta=theta_keep)


# ----------------------------
# High-level backend
# ----------------------------

@dataclass
class OperatorBackendConfig:
    dtype_compute: torch.dtype = torch.float32
    tol: float = 1e-5
    max_iter: int = 500

    point_predictions_only: bool = True
    return_varcov: bool = False
    return_se: bool = False

    # Hutchinson probe settings for diagonal prediction SE.
    # Probe noise scales as 1/sqrt(M) and is the dominant source of SE error at
    # practical M; loose PCG tol is fine because the Monte Carlo floor dominates.
    # Default raised from 32 to 256: 32 probes give unbiased-in-mean but
    # per-genotype-noisy PEV (empirically corr ~0.8 even at 4096 probes), which
    # is unsafe for reporting per-line reliability. Paired with Lanczos variance
    # reduction (enabled by default below).
    n_hutchinson_probes: int = 256
    hutchinson_seed: int = 0
    se_pcg_tol: float = 1e-3
    se_pcg_max_iter: int = 150

    # Nyström preconditioner for the SE PCG solves. Captures top-ℓ eigenpairs of V
    # so each Hutchinson probe converges in far fewer iterations, enabling larger M
    # at the same wall budget. Does NOT reduce per-probe estimator variance on its
    # own — the Monte Carlo variance comes from the off-diagonals of K V^{-1} K^T
    # and is independent of how fast PCG converges.
    use_nystrom_se_preconditioner: bool = False
    nystrom_rank: int = 64
    nystrom_seed: int = 0

    # Lanczos low-rank V^{-1} control variate for SE Hutchinson variance reduction.
    # Builds a rank-`lanczos_iters` Krylov subspace on V (starting from K^T z so
    # the span aligns with K's range), extracts Ritz pairs, and subtracts the
    # deterministic diag(K V_low K^T) from the Hutchinson estimator — the residual
    # Hutchinson term has smaller Frobenius off-diagonals, reducing Monte Carlo
    # variance per probe. Complementary to Nyström preconditioning (which only
    # cuts PCG iters); the Lanczos control variate actually reduces estimator
    # variance. ON by default: it cuts per-genotype PEV variance by ~5 orders of
    # magnitude at negligible cost (see test_lanczos_variance_reduction.py) -- the
    # deterministic low-rank term carries the bulk of diag(K V^-1 K^T), so only a
    # tiny residual is estimated stochastically. Without it, plain Hutchinson is
    # unbiased but hopelessly noisy per-genotype even at 256 probes.
    use_lanczos_se_variance_reduction: bool = True
    lanczos_iters: int = 80
    lanczos_seed: int = 0
    lanczos_theta_floor: float = 1e-8


def _pev_mc_rel_error(probe_sum, probe_sumsq, n_probes, pev, *, rel_floor_frac=0.1):
    """Monte-Carlo relative error of the Hutchinson PEV estimate.

    ``probe_sum`` / ``probe_sumsq`` are the running sum and sum-of-squares of the
    per-probe diagonal contributions ``z * (s - s_low)``. PEV is the Hutchinson
    *mean* of those contributions (subtracted from deterministic terms), so its
    Monte-Carlo standard error is ``sqrt(var_probe / M)`` where ``var_probe`` is
    the unbiased per-probe sample variance; the deterministic self-variance and
    Lanczos pieces contribute no MC noise.

    Returns ``(mc_se, rel_median, rel_max)``: ``mc_se`` is the per-row MC standard
    error of PEV; ``rel_*`` summarise ``mc_se / pev`` over rows whose PEV exceeds
    ``rel_floor_frac * median(positive PEV)`` (near-zero PEV rows are excluded so
    the ratio does not blow up). With ``M < 2`` the error is undefined (NaN).
    """
    M = int(n_probes)
    if M < 2:
        nan = float("nan")
        return torch.full_like(pev, nan), nan, nan
    mean = probe_sum / M
    var = (probe_sumsq - M * mean * mean) / (M - 1)
    var = var.clamp_min(0.0)
    mc_se = torch.sqrt(var / M)
    pos = pev[pev > 0]
    if pos.numel() == 0:
        return mc_se, float("nan"), float("nan")
    scale = float(torch.median(pos))
    mask = pev > (rel_floor_frac * scale)
    if not bool(mask.any()):
        return mc_se, float("nan"), float("nan")
    rel = mc_se[mask] / pev[mask]
    # torch.median returns the lower middle for even counts; use quantile(0.5)
    # for the standard (interpolated) median so the diagnostic is conventional.
    return mc_se, float(torch.quantile(rel, 0.5)), float(torch.max(rel))


def fit_predict_operator_backend(
    y: torch.Tensor,
    X: torch.Tensor,
    gi: torch.Tensor,
    ei: torch.Tensor,
    cache: GRMFactorCacheBase,
    Sigma_e: torch.Tensor,
    train_idx: torch.Tensor,
    test_idx: torch.Tensor,
    w_g: torch.Tensor,
    w_ge: torch.Tensor,
    w_e: Union[float, torch.Tensor],
    S_e: Optional[torch.Tensor],
    r_diag: Union[float, torch.Tensor],
    config: Optional[OperatorBackendConfig] = None,
    random_terms: Optional[List[Dict[str, object]]] = None,
) -> Dict[str, object]:
    """
    Operator backend for GP posterior prediction at scale.

    This backend returns point predictions and, when ``config.return_se`` is
    true, diagonal prediction SEs via Hutchinson probing. Full prediction
    covariance is intentionally not implemented here.

    Returns:
      {
        "predictions": np.ndarray (n_test,),
        "beta": np.ndarray (p,),
        "pcg_info": {...},
        "prediction_meta": {...},
        "notes": ...
      }
    """
    if config is None:
        config = OperatorBackendConfig()

    device = y.device
    dt = config.dtype_compute

    train_idx = train_idx.long()
    test_idx = test_idx.long()

    y_tr = y.index_select(0, train_idx).to(dt)
    X_tr = X.index_select(0, train_idx).to(dt)
    gi_tr = gi.index_select(0, train_idx).long()
    ei_tr = ei.index_select(0, train_idx).long()

    # Build operators from random_terms spec (or default fast path)
    # If random_terms include IID terms specified on full observation set (group_index_all),
    # slice them to training for building V, and keep test indices for prediction cross-cov.
    iid_terms_test = []
    random_terms_train = None
    if random_terms is not None:
        random_terms_train = []
        for term in random_terms:
            t = dict(term)
            ttype = str(t.get("type","")).lower()
            if ttype in ("iid","random_intercept"):
                if "group_index_all" in t:
                    g_all = _as_torch(t["group_index_all"], device, torch.int64)
                    g_tr = g_all.index_select(0, train_idx.long())
                    g_te = g_all.index_select(0, test_idx.long())
                    t["group_index"] = g_tr
                    iid_terms_test.append((g_tr, g_te, t.get("weight", 1.0), t.get("name", None)))
                elif "group_index" in t:
                    # assume already aligned to training; cannot compute cross-cov without test indices
                    pass
            random_terms_train.append(t)
    else:
        random_terms_train = None
    
    # Residual diagonal for training fold
    if isinstance(r_diag, (float, int)):
        r_diag_tr = float(r_diag)
    else:
        r_diag_t = _as_torch(r_diag, device, dt)
        # If r_diag is full-observation length, slice to training; else assume already aligned.
        if r_diag_t.numel() == gi.numel():
            r_diag_tr = r_diag_t.index_select(0, train_idx.long())
        else:
            r_diag_tr = r_diag_t

    ops = build_random_terms_operator_backend(
        gi=gi_tr,
        ei=ei_tr,
        cache=cache,
        Sigma_e=_as_torch(Sigma_e, device, dt),
        S_e=None if S_e is None else _as_torch(S_e, device, dt),
        w_g=_as_torch(w_g, device, dt),
        w_ge=_as_torch(w_ge, device, dt),
        w_e=w_e,
        r_diag=r_diag_tr,
        random_terms=random_terms_train,
    )

    op_g = next((o for o in ops if isinstance(o, LowRankGRMOperator)), None)
    op_ge = next((o for o in ops if isinstance(o, LowRankGxEOperator)), None)
    op_e = next((o for o in ops if isinstance(o, EnvKernelOperator)), None)
    op_r = next((o for o in ops if isinstance(o, DiagResidualOperator)), None)

    V = SumOperator(ops)

    context = {}
    if op_g is not None:
        fold_g = op_g.build_fold_cache()
        context["fold_g"] = fold_g
    if op_ge is not None:
        fold_ge, ei_fold = op_ge.build_fold_cache()
        context["fold_ge"] = fold_ge
        context["ei_fold"] = ei_fold
    if op_e is not None:
        ei_fold_env = op_e.build_fold_cache()
        context["ei_fold_env"] = ei_fold_env

    # Jacobi preconditioner: approximate diag(V)
    @torch.no_grad()
    def approx_diag_V():
        diag = torch.zeros_like(y_tr, dtype=dt)

        if op_g is not None:
            fold_g_local = context["fold_g"]
            for r in range(cache.num_grms):
                Phi_used = cache.get_rows(r, fold_g_local.g_unique, device=device, dtype=dt)
                rownorm2 = (Phi_used * Phi_used).sum(dim=1)
                diag = diag + op_g.weights[r].to(dt) * rownorm2.index_select(0, fold_g_local.inv)

        if op_ge is not None:
            fold_ge_local = context["fold_ge"]
            ei_fold_local = context["ei_fold"]
            Se_diag = torch.diag(op_ge.Sigma_e.to(device=device, dtype=dt))
            for r in range(cache.num_grms):
                Phi_used = cache.get_rows(r, fold_ge_local.g_unique, device=device, dtype=dt)
                rownorm2 = (Phi_used * Phi_used).sum(dim=1)
                diag = diag + op_ge.weights[r].to(dt) * rownorm2.index_select(0, fold_ge_local.inv) * Se_diag.index_select(0, ei_fold_local)

        if op_e is not None:
            ei_fold_env_local = context.get("ei_fold_env", ei_tr)
            diag = diag + op_e.weight.to(dt) * torch.diag(op_e.S_e.to(device=device, dtype=dt)).index_select(0, ei_fold_env_local)


        # IID random intercept terms contribute constant weight to diagonal
        if random_terms_train is not None:
            for term in random_terms_train:
                if str(term.get("type","")).lower() in ("iid","random_intercept"):
                    try:
                        w = float(term.get("weight", 1.0))
                    except Exception:
                        w = 1.0
                    diag = diag + torch.as_tensor(w, device=device, dtype=dt)
        if op_r is not None:
            diag = diag + op_r.diag.to(device=device, dtype=dt)

        return diag.clamp_min(torch.tensor(1e-6, device=device, dtype=dt))

    diagV = approx_diag_V()

    def M_inv_mv(v):
        return v / diagV

    def V_mv(v):
        return V.matvec(v, context=context, dtype_compute=dt)

    # GLS beta: beta = (X'V^-1X)^-1 X'V^-1 y
    Vinv_y = pcg_solve(V_mv, y_tr, M_inv_mv=M_inv_mv, tol=config.tol, max_iter=config.max_iter)
    p = X_tr.shape[1]
    Vinv_X_cols = []
    pcg_iters = [Vinv_y.iters]
    converged_X_cols: List[bool] = []
    rel_resid_X_cols: List[float] = []
    for j in range(p):
        sol = pcg_solve(V_mv, X_tr[:, j], M_inv_mv=M_inv_mv, tol=config.tol, max_iter=config.max_iter)
        Vinv_X_cols.append(sol.x)
        pcg_iters.append(sol.iters)
        converged_X_cols.append(bool(sol.converged))
        rel_resid_X_cols.append(float(sol.rel_resid))
    Vinv_X = torch.stack(Vinv_X_cols, dim=1)

    XtVinvX = X_tr.T @ Vinv_X
    XtVinvY = X_tr.T @ Vinv_y.x

    # Stable small solve in float64
    XtVinvX64 = XtVinvX.to(torch.float64)
    XtVinvY64 = XtVinvY.to(torch.float64)
    L = torch.linalg.cholesky(XtVinvX64 + 1e-8 * torch.eye(p, device=device, dtype=torch.float64))
    beta = torch.cholesky_solve(XtVinvY64.unsqueeze(-1), L).squeeze(-1).to(dt)

    # alpha = V^-1 (y - X beta)
    r_tr = y_tr - X_tr @ beta
    alpha = pcg_solve(V_mv, r_tr, M_inv_mv=M_inv_mv, tol=config.tol, max_iter=config.max_iter)
    pcg_iters.append(alpha.iters)

    # Predict on test: yhat = X_te beta + cov(te,tr) alpha
    X_te = X.index_select(0, test_idx).to(dt)
    gi_te = gi.index_select(0, test_idx).long()
    ei_te = ei.index_select(0, test_idx).long()

    pred_cov = torch.zeros((test_idx.numel(),), device=device, dtype=dt)

    fold_g_train = context.get("fold_g", None)
    fold_ge_train = context.get("fold_ge", None)

    # Build per-op test_states once so SE block can reuse them.
    op_test_states: List[Tuple[torch.nn.Module, Dict[str, object]]] = []
    for op in ops:
        if isinstance(op, LowRankGRMOperator):
            if fold_g_train is None:
                continue
            ts = {"fold_train": fold_g_train, "gi_test": gi_te}
        elif isinstance(op, LowRankGxEOperator):
            if fold_ge_train is None or "ei_fold" not in context:
                continue
            ts = {
                "fold_train": fold_ge_train,
                "ei_train": context["ei_fold"],
                "gi_test": gi_te,
                "ei_test": ei_te,
            }
        elif isinstance(op, EnvKernelOperator):
            ts = {"ei_train": ei_tr.long(), "ei_test": ei_te}
        elif isinstance(op, IIDGroupOperator):
            term_name = getattr(op, "term_name", None)
            matched = None
            if iid_terms_test:
                for item in iid_terms_test:
                    if len(item) == 4:
                        g_tr, g_te, w, nm = item
                        if term_name == nm:
                            matched = (g_tr, g_te)
                            break
                    else:
                        g_tr, g_te, w = item
                        matched = (g_tr, g_te)
                        break
            if matched is None:
                continue
            g_tr, g_te = matched
            ts = {"group_train": g_tr, "group_test": g_te}
        else:
            continue
        op_test_states.append((op, ts))
        pred_cov = pred_cov + op.cross_cov_apply(alpha.x, ts, dtype_compute=dt)

    yhat = X_te @ beta + pred_cov

    # Hutchinson probe-based diagonal prediction SE
    prediction_se_latent_np = None
    prediction_se_observed_np = None
    prediction_pev_np = None
    se_info: Dict[str, object] = {}
    if config.return_se:
        n_te = int(test_idx.numel())
        K_te_te_diag = torch.zeros((n_te,), device=device, dtype=dt)
        for op, ts in op_test_states:
            K_te_te_diag = K_te_te_diag + op.self_var_diag(ts, dtype_compute=dt)

        M_probes = int(config.n_hutchinson_probes)
        gen = torch.Generator(device=device).manual_seed(int(config.hutchinson_seed))
        diag_sum = torch.zeros((n_te,), device=device, dtype=dt)
        diag_sumsq = torch.zeros((n_te,), device=device, dtype=dt)
        se_tol = float(config.se_pcg_tol)
        se_max_it = int(config.se_pcg_max_iter)
        probe_iters: List[int] = []
        probe_converged: List[bool] = []
        probe_rel_resids: List[float] = []

        # Optional Nyström preconditioner for the SE PCG solves.
        se_M_inv_mv = M_inv_mv
        if bool(config.use_nystrom_se_preconditioner):
            _, P_inv_nys = build_nystrom_preconditioner(
                V_mv,
                n=int(y_tr.numel()),
                rank=int(config.nystrom_rank),
                device=device,
                dtype=dt,
                seed=int(config.nystrom_seed),
            )
            se_M_inv_mv = P_inv_nys

        # Optional Lanczos low-rank V^{-1} control variate for variance reduction.
        # Builds Q_k ≈ Krylov(V, K^T z0) once and precomputes F = K Q_k so the
        # deterministic piece diag(K V_low K^T) and the per-probe s_low projection
        # are cheap. Hutchinson then only estimates the residual diag(K (V^{-1} -
        # V_low) K^T), which has smaller off-diagonals and so smaller MC variance.
        lanczos_obj: Optional[LanczosLowRankInv] = None
        F_lan: Optional[torch.Tensor] = None
        det_diag_lan: Optional[torch.Tensor] = None
        inv_theta_lan: Optional[torch.Tensor] = None
        if bool(config.use_lanczos_se_variance_reduction):
            gen_init = torch.Generator(device=device).manual_seed(
                int(config.lanczos_seed) * 2654435761 + 1
            )
            z0_bits = torch.randint(
                0, 2, (n_te,), device=device, generator=gen_init, dtype=torch.int64
            )
            z0 = z0_bits.to(dt) * 2.0 - 1.0
            u0 = torch.zeros_like(y_tr)
            for op, ts in op_test_states:
                u0 = u0 + op.cross_cov_apply_transpose(z0, ts, dtype_compute=dt)
            u0_norm = float(u0.norm().item())
            start = (u0 / u0_norm) if (np.isfinite(u0_norm) and u0_norm > 1e-20) else None
            lanczos_obj = build_lanczos_lowrank_inv(
                V_mv,
                n=int(y_tr.numel()),
                k_iter=int(config.lanczos_iters),
                device=device,
                dtype=dt,
                seed=int(config.lanczos_seed),
                start_vec=start,
                theta_floor=float(config.lanczos_theta_floor),
            )
            if lanczos_obj.theta.numel() > 0:
                V_vecs = lanczos_obj.V_vecs  # (n_tr, k)
                k_lan = int(V_vecs.shape[1])
                F_lan = torch.zeros((n_te, k_lan), device=device, dtype=dt)
                for j in range(k_lan):
                    v = V_vecs[:, j]
                    s_col = torch.zeros((n_te,), device=device, dtype=dt)
                    for op, ts in op_test_states:
                        s_col = s_col + op.cross_cov_apply(v, ts, dtype_compute=dt)
                    F_lan[:, j] = s_col
                inv_theta_lan = (1.0 / lanczos_obj.theta).to(dt)
                det_diag_lan = (F_lan * F_lan) @ inv_theta_lan
            else:
                # Lanczos failed to find positive Ritz values — skip CV quietly.
                lanczos_obj = None

        for m_idx in range(M_probes):
            z_bits = torch.randint(0, 2, (n_te,), device=device, generator=gen, dtype=torch.int64)
            z = z_bits.to(dt) * 2.0 - 1.0
            u = torch.zeros_like(y_tr)
            for op, ts in op_test_states:
                u = u + op.cross_cov_apply_transpose(z, ts, dtype_compute=dt)
            w_sol = pcg_solve(V_mv, u, M_inv_mv=se_M_inv_mv, tol=se_tol, max_iter=se_max_it)
            probe_iters.append(w_sol.iters)
            probe_converged.append(bool(w_sol.converged))
            probe_rel_resids.append(float(w_sol.rel_resid))
            s = torch.zeros((n_te,), device=device, dtype=dt)
            for op, ts in op_test_states:
                s = s + op.cross_cov_apply(w_sol.x, ts, dtype_compute=dt)
            if lanczos_obj is not None:
                # s_low_i = (K V_low u)_i = F (Q^T u) / theta per probe
                qt_u = lanczos_obj.V_vecs.transpose(0, 1) @ u  # (k,)
                s_low = F_lan @ (qt_u * inv_theta_lan)          # (n_te,)
                term = z * (s - s_low)
            else:
                term = z * s
            diag_sum = diag_sum + term
            diag_sumsq = diag_sumsq + term * term
        diag_est = diag_sum / float(max(M_probes, 1))

        if lanczos_obj is not None and det_diag_lan is not None:
            # PEV = K_te_te_diag - det(K V_low K^T) - Hutch(K (V^{-1}-V_low) K^T)
            pev = (K_te_te_diag - det_diag_lan - diag_est).clamp_min(0.0)
        else:
            pev = (K_te_te_diag - diag_est).clamp_min(0.0)
        prediction_pev_np = pev.detach().cpu().numpy()
        prediction_se_latent_np = torch.sqrt(pev).detach().cpu().numpy()

        # Monte-Carlo relative-error diagnostic for the stochastic PEV estimate.
        pev_mc_se_t, pev_mc_rel_error_median, pev_mc_rel_error_max = _pev_mc_rel_error(
            diag_sum, diag_sumsq, M_probes, pev
        )
        if np.isfinite(pev_mc_rel_error_median) and pev_mc_rel_error_median > 0.10:
            warnings.warn(
                "operator backend SE: estimated Monte-Carlo relative error of PEV "
                f"~{pev_mc_rel_error_median:.2f} (median over test rows) exceeds 0.10. "
                "Per-genotype SE/reliability are noisy; raise n_hutchinson_probes or "
                "keep use_lanczos_se_variance_reduction enabled.",
                RuntimeWarning, stacklevel=2,
            )

        if isinstance(r_diag, (float, int)):
            r_diag_te_val = torch.full((n_te,), float(r_diag), device=device, dtype=dt)
        else:
            rd_t = _as_torch(r_diag, device, dt)
            if rd_t.numel() == gi.numel():
                r_diag_te_val = rd_t.index_select(0, test_idx.long())
            else:
                r_diag_te_val = rd_t
        prediction_se_observed_np = torch.sqrt(pev + r_diag_te_val).detach().cpu().numpy()

        n_fail = sum(1 for c in probe_converged if not c)
        if n_fail > 0:
            worst = max(probe_rel_resids) if probe_rel_resids else float("nan")
            warnings.warn(
                f"operator backend SE: {n_fail}/{M_probes} Hutchinson PCG solves did not converge "
                f"(max rel_resid={worst:.2e}). Predictions SEs may be biased; raise se_pcg_max_iter "
                f"or loosen se_pcg_tol.",
                RuntimeWarning, stacklevel=2,
            )
        se_info = {
            "method": "hutchinson_rademacher",
            "n_probes": M_probes,
            "seed": int(config.hutchinson_seed),
            "probe_iters": probe_iters,
            "probe_converged": probe_converged,
            "probe_rel_resids": probe_rel_resids,
            "se_pcg_tol": se_tol,
            "se_pcg_max_iter": se_max_it,
            "preconditioner": (
                "nystrom" if bool(config.use_nystrom_se_preconditioner) else "jacobi"
            ),
            "nystrom_rank": (
                int(config.nystrom_rank) if bool(config.use_nystrom_se_preconditioner) else None
            ),
            "variance_reduction": (
                "lanczos_lowrank_cv"
                if (bool(config.use_lanczos_se_variance_reduction) and lanczos_obj is not None)
                else None
            ),
            "lanczos_iters": (
                int(config.lanczos_iters)
                if bool(config.use_lanczos_se_variance_reduction) else None
            ),
            "lanczos_kept": (
                int(lanczos_obj.theta.numel())
                if (bool(config.use_lanczos_se_variance_reduction) and lanczos_obj is not None)
                else None
            ),
            "lanczos_theta_min": (
                float(lanczos_obj.theta.min().item())
                if (bool(config.use_lanczos_se_variance_reduction)
                    and lanczos_obj is not None
                    and lanczos_obj.theta.numel() > 0)
                else None
            ),
            "lanczos_theta_max": (
                float(lanczos_obj.theta.max().item())
                if (bool(config.use_lanczos_se_variance_reduction)
                    and lanczos_obj is not None
                    and lanczos_obj.theta.numel() > 0)
                else None
            ),
            "pev_mc_rel_error_median": pev_mc_rel_error_median,
            "pev_mc_rel_error_max": pev_mc_rel_error_max,
            "pev_mc_se": pev_mc_se_t.detach().cpu().numpy(),
        }

    # Non-convergence in any of the three PCG calls biases the solution
    # (beta drifts toward OLS when Vinv_y / Vinv_X cols fail; predictions drift
    # toward the prior when alpha fails). Surface this as a warning instead of
    # silently producing biased point predictions.
    failed = []
    if not Vinv_y.converged:
        failed.append(f"V^-1 y (rel_resid={Vinv_y.rel_resid:.2e})")
    n_X_failed = sum(1 for c in converged_X_cols if not c)
    if n_X_failed > 0:
        worst = max(rel_resid_X_cols) if rel_resid_X_cols else float("nan")
        failed.append(f"{n_X_failed}/{len(converged_X_cols)} X-cols (max rel_resid={worst:.2e})")
    if not alpha.converged:
        failed.append(f"alpha (rel_resid={alpha.rel_resid:.2e})")
    if failed:
        warnings.warn(
            "operator backend: PCG did not converge for "
            + "; ".join(failed)
            + f" at tol={config.tol} max_iter={config.max_iter}. "
            "Predictions may be biased; raise max_iter, loosen tol, or improve the preconditioner.",
            RuntimeWarning, stacklevel=2,
        )

    result: Dict[str, object] = {
        "predictions": yhat.detach().cpu().numpy(),
        "beta": beta.detach().cpu().numpy(),
        "pcg_info": {
            "tol": config.tol,
            "max_iter": config.max_iter,
            "iters": pcg_iters,
            "converged_beta_y": Vinv_y.converged,
            "converged_alpha": alpha.converged,
            "converged_X_cols": converged_X_cols,
            "rel_resid_y": Vinv_y.rel_resid,
            "rel_resid_alpha": alpha.rel_resid,
            "rel_resid_X_cols": rel_resid_X_cols,
            "all_converged": (bool(Vinv_y.converged)
                              and bool(alpha.converged)
                              and all(converged_X_cols)),
        },
        "prediction_meta": {
            "backend_used": "operator",
            "point_predictions_only": (not config.return_se),
            "supports_prediction_se": bool(config.return_se),
            "supports_prediction_cov": False,
            "prediction_se_method": ("hutchinson_rademacher" if config.return_se else None),
            "operator_limitations": (
                "operator backend: diagonal prediction SEs via Hutchinson probing; full prediction covariance not implemented"
                if config.return_se else
                "operator backend returned GP posterior point predictions; diagonal SE disabled (set config.return_se=True to enable Hutchinson SE)"
            ),
        },
        "notes": (
            "Operator backend (low-rank GRM roots + gather/scatter). Hutchinson diagonal SE enabled."
            if config.return_se else
            "Operator backend (low-rank GRM roots + gather/scatter). GP posterior point predictions without SE."
        ),
    }
    if config.return_se:
        result["prediction_se_latent"] = prediction_se_latent_np
        result["prediction_se_observed"] = prediction_se_observed_np
        result["prediction_pev"] = prediction_pev_np
        result["se_info"] = se_info
    return result
