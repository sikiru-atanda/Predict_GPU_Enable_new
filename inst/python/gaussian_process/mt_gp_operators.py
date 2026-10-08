"""Multi-trait GP operators for the v2-fast prediction path.

Pure matvec-only operators that compose under SumOperator from the existing
operator backend. Imports primitives from `mixed_model_gpu_large_scale_backend`
(aliased by conftest.py). Does NOT modify the operator backend file.

Public:
  MTGenoOperator    — applies (Sigma_G ⊗ K_geno) · v via Phi low-rank factor.
  MTResidualOperator — applies (Sigma_eps ⊗ I_n) · v with same-gid trait coupling.
"""
from __future__ import annotations
from typing import Optional, Tuple
import torch

from mixed_model_gpu_large_scale_backend import (
    GRMFactorCacheBase,
    FoldIndexCache,
    _scatter_add_cells,
    _unique_inverse,
)


def _kernel_weights_tensor(cache: GRMFactorCacheBase, kernel_weights: Optional[torch.Tensor]) -> torch.Tensor:
    if kernel_weights is None:
        return torch.ones(int(cache.num_grms), dtype=torch.float64)
    w = torch.as_tensor(kernel_weights, dtype=torch.float64).reshape(-1)
    if int(w.numel()) != int(cache.num_grms):
        raise ValueError(
            f"kernel_weights length {int(w.numel())} != cache.num_grms={int(cache.num_grms)}"
        )
    if not bool(torch.isfinite(w).all()):
        raise ValueError("kernel_weights must be finite")
    if bool((w < 0).any()):
        raise ValueError("kernel_weights must be non-negative")
    if float(w.sum().item()) <= 0.0:
        raise ValueError("at least one kernel weight must be positive")
    return w


class MTGenoOperator(torch.nn.Module):
    """Applies (Sigma_G ⊗ K_geno) · v via the Phi low-rank factor.

    Math identical to `LowRankGxEOperator` with `ei → ti` and `Sigma_e → Sigma_G`.
    Supports one or more low-rank GRM roots through `kernel_weights`.
    """

    def __init__(self, cache: GRMFactorCacheBase, gi: torch.Tensor,
                 ti: torch.Tensor, Sigma_G: torch.Tensor,
                 kernel_weights: Optional[torch.Tensor] = None):
        super().__init__()
        self.cache = cache
        self.register_buffer("gi", gi.long())
        self.register_buffer("ti", ti.long())
        self.register_buffer("Sigma_G", Sigma_G)
        self.register_buffer("kernel_weights", _kernel_weights_tensor(cache, kernel_weights))

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None) -> Tuple[FoldIndexCache, torch.Tensor]:
        if obs_index is None:
            gi = self.gi
            ti = self.ti
        else:
            obs_index = obs_index.long()
            gi = self.gi.index_select(0, obs_index)
            ti = self.ti.index_select(0, obs_index)
        g_unique, inv = _unique_inverse(gi)
        return FoldIndexCache(g_unique=g_unique, inv=inv, ng_used=int(g_unique.numel())), ti

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, ti_fold: torch.Tensor,
               dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        Sigma_G = self.Sigma_G.to(device=device, dtype=dtype_compute)
        weights = self.kernel_weights.to(device=device, dtype=dtype_compute)
        T = int(Sigma_G.shape[0])
        T_acc = _scatter_add_cells(v, fold.inv, ti_fold.long(), k_env=T, ng_used=fold.ng_used)
        B = torch.zeros((fold.ng_used, T), dtype=dtype_compute, device=device)
        for r in range(int(self.cache.num_grms)):
            Phi_used = self.cache.get_rows(r, fold.g_unique, device=device, dtype=dtype_compute)
            A = Phi_used.transpose(0, 1) @ T_acc
            B = B + weights[r] * (Phi_used @ A)
        U = B @ Sigma_G.transpose(0, 1)
        return U[fold.inv, ti_fold.long()]

    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor,
                        test_state: dict,
                        dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        """Predict at test rows: Sigma_G[t_test, t_train] · K[g_test, g_train] · alpha_train."""
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        fold_train: FoldIndexCache = test_state["fold_train"]
        ti_train = test_state["ti_train"].long()
        gi_test = test_state["gi_test"].long()
        ti_test = test_state["ti_test"].long()
        Sigma_G = self.Sigma_G.to(device=device, dtype=dtype_compute)
        weights = self.kernel_weights.to(device=device, dtype=dtype_compute)
        T = int(Sigma_G.shape[0])
        T_acc = _scatter_add_cells(alpha_train, fold_train.inv, ti_train,
                                   k_env=T, ng_used=fold_train.ng_used)
        B_te = torch.zeros((gi_test.numel(), T), dtype=dtype_compute, device=device)
        for r in range(int(self.cache.num_grms)):
            Phi_tr = self.cache.get_rows(r, fold_train.g_unique, device=device, dtype=dtype_compute)
            A = Phi_tr.transpose(0, 1) @ T_acc
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            B_te = B_te + weights[r] * (Phi_te @ A)
        U_te = B_te @ Sigma_G.transpose(0, 1)
        return U_te[torch.arange(gi_test.numel(), device=device), ti_test]


class MTResidualOperator(torch.nn.Module):
    """Applies (Sigma_eps ⊗ I_n) · v with same-gid trait coupling.

    Cov(eps_r, eps_s) = Sigma_eps[t_r, t_s] · δ(g_r == g_s).
    """

    def __init__(self, gi: torch.Tensor, ti: torch.Tensor, Sigma_eps: torch.Tensor):
        super().__init__()
        self.register_buffer("gi", gi.long())
        self.register_buffer("ti", ti.long())
        self.register_buffer("Sigma_eps", Sigma_eps)

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None) -> Tuple[FoldIndexCache, torch.Tensor]:
        if obs_index is None:
            gi = self.gi; ti = self.ti
        else:
            obs_index = obs_index.long()
            gi = self.gi.index_select(0, obs_index)
            ti = self.ti.index_select(0, obs_index)
        g_unique, inv = _unique_inverse(gi)
        return FoldIndexCache(g_unique=g_unique, inv=inv, ng_used=int(g_unique.numel())), ti

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, ti_fold: torch.Tensor,
               dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        Sigma_eps = self.Sigma_eps.to(device=device, dtype=dtype_compute)
        T = int(Sigma_eps.shape[0])
        T_acc = _scatter_add_cells(v, fold.inv, ti_fold.long(), k_env=T, ng_used=fold.ng_used)
        out_per_gid = T_acc @ Sigma_eps.transpose(0, 1)
        return out_per_gid[fold.inv, ti_fold.long()]


class MTMETResidualOperator(torch.nn.Module):
    """Applies (Sigma_eps ⊗ I_n_cells) · v with same-(gid, env)-cell trait coupling.

    Cov(eps_r, eps_s) = Sigma_eps[t_r, t_s] · δ((g_r, e_r) == (g_s, e_s)).

    Internally uses a per-(gid, env) cell index in the scatter, not just gid.
    """

    def __init__(self, gi: torch.Tensor, ei: torch.Tensor, ti: torch.Tensor,
                 Sigma_eps: torch.Tensor):
        super().__init__()
        self.register_buffer("gi", gi.long())
        self.register_buffer("ei", ei.long())
        self.register_buffer("ti", ti.long())
        self.register_buffer("Sigma_eps", Sigma_eps)
        # Hoist N_E_max to __init__ so all build_fold_cache calls use the same encoding
        # regardless of which obs_index slice is passed.
        self._N_E_max = int(self.ei.max().item()) + 1 if self.ei.numel() > 0 else 1

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None
                         ) -> Tuple[FoldIndexCache, torch.Tensor, torch.Tensor]:
        """Returns (cell_fold, ei_fold, ti_fold). cell_fold.inv maps each row to its
        (gid, env) cell index. cell_fold.g_unique stores unique cell IDs (encoded as gid*N_E + env).
        """
        if obs_index is None:
            gi = self.gi; ei = self.ei; ti = self.ti
        else:
            obs_index = obs_index.long()
            gi = self.gi.index_select(0, obs_index)
            ei = self.ei.index_select(0, obs_index)
            ti = self.ti.index_select(0, obs_index)
        # Encode (gid, env) into a single int via gid * N_E_max + env
        cell_id = gi * self._N_E_max + ei
        cell_unique, inv = _unique_inverse(cell_id)
        return (FoldIndexCache(g_unique=cell_unique, inv=inv,
                               ng_used=int(cell_unique.numel())),
                ei, ti)

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, ti_fold: torch.Tensor,
               dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        """Same signature as MTResidualOperator (3-arg). fold.inv must come from
        build_fold_cache (per-cell, NOT per-gid)."""
        device = v.device
        v = v.to(dtype=dtype_compute)
        Sigma_eps = self.Sigma_eps.to(device=device, dtype=dtype_compute)
        T = int(Sigma_eps.shape[0])
        T_acc = _scatter_add_cells(v, fold.inv, ti_fold.long(), k_env=T, ng_used=fold.ng_used)
        out_per_cell = T_acc @ Sigma_eps.transpose(0, 1)
        return out_per_cell[fold.inv, ti_fold.long()]


class _MTSumOperator:
    """Composes any number of MT operators into a single matvec callable.

    Each operator must expose `matvec(v, fold, ti_fold, dtype_compute=...)`.
    Lightweight (no torch.nn.Module needed; not optimized over). Avoids modifying
    the backend's SumOperator (which dispatches on hardcoded operator types).
    """

    def __init__(self, *ops):
        if len(ops) == 0:
            raise ValueError("_MTSumOperator needs at least one operator")
        self.ops = list(ops)
        # Backward-compatible attributes (used by some v2-fast call sites)
        self.op_g = ops[0]
        self.op_e = ops[1] if len(ops) >= 2 else None

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache, ti_fold: torch.Tensor,
               dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        out = self.ops[0].matvec(v, fold, ti_fold, dtype_compute=dtype_compute)
        for op in self.ops[1:]:
            out = out + op.matvec(v, fold, ti_fold, dtype_compute=dtype_compute)
        return out


class MTGxEOperator(torch.nn.Module):
    """Applies (Sigma_GE ⊗ K_geno ⊗ K_env) · v via Phi (low-rank K_geno) and dense K_env.

    Math: for row r at (gid g_r, env e_r, trait t_r):
      out[r] = sum_{(g',e',t') in train} Sigma_GE[t_r, t'] * K_geno[g_r, g'] * K_env[e_r, e'] * v_{g',e',t'}

    Matvec: scatter v to 3-tensor (n_gid_used, n_env, T), apply K_env on env axis (dense matmul),
    Phi (Phi^T ·) on gid axis, Sigma_GE on trait axis, then gather back.

    Supports one or more low-rank GRM roots through `kernel_weights`.
    """

    def __init__(self, cache: GRMFactorCacheBase, gi: torch.Tensor, ei: torch.Tensor,
                 ti: torch.Tensor, Sigma_GE: torch.Tensor, K_env: torch.Tensor,
                 kernel_weights: Optional[torch.Tensor] = None):
        super().__init__()
        self.cache = cache
        self.register_buffer("gi", gi.long())
        self.register_buffer("ei", ei.long())
        self.register_buffer("ti", ti.long())
        self.register_buffer("Sigma_GE", Sigma_GE)
        self.register_buffer("K_env", K_env)
        self.register_buffer("kernel_weights", _kernel_weights_tensor(cache, kernel_weights))

    def build_fold_cache(self, obs_index: Optional[torch.Tensor] = None
                         ) -> Tuple[FoldIndexCache, torch.Tensor, torch.Tensor]:
        if obs_index is None:
            gi = self.gi; ei = self.ei; ti = self.ti
        else:
            obs_index = obs_index.long()
            gi = self.gi.index_select(0, obs_index)
            ei = self.ei.index_select(0, obs_index)
            ti = self.ti.index_select(0, obs_index)
        g_unique, inv = _unique_inverse(gi)
        return (FoldIndexCache(g_unique=g_unique, inv=inv, ng_used=int(g_unique.numel())),
                ei, ti)

    @torch.no_grad()
    def matvec(self, v: torch.Tensor, fold: FoldIndexCache,
               ei_fold: torch.Tensor, ti_fold: torch.Tensor,
               dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        device = v.device
        v = v.to(dtype=dtype_compute)
        Sigma_GE = self.Sigma_GE.to(device=device, dtype=dtype_compute)
        K_env = self.K_env.to(device=device, dtype=dtype_compute)
        weights = self.kernel_weights.to(device=device, dtype=dtype_compute)
        T = int(Sigma_GE.shape[0])
        n_env = int(K_env.shape[0])
        ng_used = fold.ng_used

        # Scatter v into 3-tensor (ng_used, n_env, T)
        V_grid = torch.zeros((ng_used, n_env, T), dtype=dtype_compute, device=device)
        flat_idx = (fold.inv * (n_env * T) + ei_fold.long() * T + ti_fold.long())
        V_grid.view(-1).index_add_(0, flat_idx, v)

        # Apply K_env on env axis: V_env[g, e, t] = sum_{e'} K_env[e, e'] * V_grid[g, e', t]
        V_env = torch.einsum("ef,gft->get", K_env, V_grid)

        # Apply K_geno on gid axis via Phi
        V_env_flat = V_env.reshape(ng_used, n_env * T)
        B = torch.zeros((ng_used, n_env * T), dtype=dtype_compute, device=device)
        for r in range(int(self.cache.num_grms)):
            Phi_used = self.cache.get_rows(r, fold.g_unique, device=device, dtype=dtype_compute)
            A = Phi_used.transpose(0, 1) @ V_env_flat  # (m, n_env * T)
            B = B + weights[r] * (Phi_used @ A)         # (ng_used, n_env * T)
        B_grid = B.reshape(ng_used, n_env, T)

        # Apply Sigma_GE on trait axis: U[g, e, t] = sum_{t'} B_grid[g, e, t'] * Sigma_GE[t, t']
        U = torch.einsum("ges,ts->get", B_grid, Sigma_GE)

        # Gather back to long-format
        out = U.view(-1).index_select(0, flat_idx)
        return out

    @torch.no_grad()
    def cross_cov_apply(self, alpha_train: torch.Tensor,
                        test_state: dict,
                        dtype_compute: torch.dtype = torch.float32) -> torch.Tensor:
        """Predict at test rows (including unseen envs):
           Sigma_GE[t_test, t_train] · K_geno[g_test, g_train] · K_env[e_test, e_train] · alpha_train.

        test_state must contain:
          fold_train, ti_train, ei_train (for training)
          gi_test, ei_test, ti_test (for test rows; test envs may be unseen)
        """
        device = alpha_train.device
        alpha_train = alpha_train.to(dtype=dtype_compute)
        fold_train: FoldIndexCache = test_state["fold_train"]
        ei_train = test_state["ei_train"].long()
        ti_train = test_state["ti_train"].long()
        gi_test = test_state["gi_test"].long()
        ei_test = test_state["ei_test"].long()
        ti_test = test_state["ti_test"].long()
        Sigma_GE = self.Sigma_GE.to(device=device, dtype=dtype_compute)
        K_env = self.K_env.to(device=device, dtype=dtype_compute)
        weights = self.kernel_weights.to(device=device, dtype=dtype_compute)
        T = int(Sigma_GE.shape[0])
        n_env = int(K_env.shape[0])
        ng_used = fold_train.ng_used

        # Scatter alpha_train into (ng_used, n_env, T)
        alpha_grid = torch.zeros((ng_used, n_env, T), dtype=dtype_compute, device=device)
        flat_idx_tr = (fold_train.inv * (n_env * T) + ei_train * T + ti_train)
        alpha_grid.view(-1).index_add_(0, flat_idx_tr, alpha_train)

        # Apply K_env on env axis (over all envs; test rows pick the right e_test row later)
        V_env_full = torch.einsum("ef,gft->get", K_env, alpha_grid)

        # Apply Phi on gid axis: project to m, then to test gids
        V_env_flat = V_env_full.reshape(ng_used, n_env * T)
        B_te = torch.zeros((gi_test.numel(), n_env * T), dtype=dtype_compute, device=device)
        for r in range(int(self.cache.num_grms)):
            Phi_tr = self.cache.get_rows(r, fold_train.g_unique, device=device, dtype=dtype_compute)
            A = Phi_tr.transpose(0, 1) @ V_env_flat  # (m, n_env * T)
            Phi_te = self.cache.get_rows(r, gi_test, device=device, dtype=dtype_compute)
            B_te = B_te + weights[r] * (Phi_te @ A)  # (n_test, n_env * T)
        B_te_grid = B_te.reshape(B_te.shape[0], n_env, T)

        # Apply Sigma_GE on trait axis
        U_te = torch.einsum("ges,ts->get", B_te_grid, Sigma_GE)  # (n_test, n_env, T)

        # Gather (e_test[i], t_test[i]) for each test row
        n_test = gi_test.shape[0]
        return U_te[torch.arange(n_test, device=device), ei_test, ti_test]
