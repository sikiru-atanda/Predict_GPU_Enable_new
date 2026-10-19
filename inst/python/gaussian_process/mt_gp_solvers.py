"""Exact covariance solves on the observed single-environment trait grid.

Complete training-genotype grids use independent T-by-T blocks in the
training kernel's eigenbasis. Small missing sets use an exact Schur correction
of the full-grid precision; large missing sets retain dense Cholesky. These
routes preserve the covariance model, jitter and observed-data REML objective.
"""
import os
import math
import torch


def correction_memory_budget():
    """Workspace allowance, not a promise about total process peak memory."""
    mb = float(os.environ.get('PREDICTPRO_GP_MISSING_MAX_MB', '512'))
    if not math.isfinite(mb) or mb <= 0:
        raise ValueError('PREDICTPRO_GP_MISSING_MAX_MB must be finite and positive')
    return int(mb * 2**20)


def correction_workspace_bytes(n_full, n_missing, n_gid, traits, itemsize, block_size):
    # Conservative allowance for stored correction, concatenation/backward
    # intermediates, Wmm factors, block solves and kernel/eigen storage.
    block = min(n_missing, block_size)
    return itemsize * (4*n_full*n_missing + 4*n_missing*n_missing +
                       8*n_full*block + 4*n_gid*n_gid + 4*n_gid*traits*traits)


class MTTrainingLayout:
    def __init__(self, kernel, gid, trait, n_traits, eigen=None, force_dense=False,
                 missing_block_size=32):
        self.gid, self.trait, self.T = gid, trait, int(n_traits)
        self.gids, compact = torch.unique(gid, sorted=True, return_inverse=True)
        self.n_gid = self.gids.numel()
        self.n_obs = gid.numel()
        cells = compact * self.T + trait
        self.n_full = self.n_gid * self.T
        unique_cells = torch.unique(cells).numel() == self.n_obs
        self.balanced = self.n_obs == self.n_full and unique_cells and not force_dense
        n_missing = self.n_full - self.n_obs
        self.missing_block_size = max(1, int(missing_block_size))
        self.correction_workspace_bytes = correction_workspace_bytes(
            self.n_full, max(0,n_missing), self.n_gid, self.T,
            kernel.element_size(), self.missing_block_size)
        # The fraction limit controls correction work; the configurable memory
        # budget replaces fixed element-count/missing-count limits. Exact only.
        self.missing_correction = (not force_dense and unique_cells and
            0 < n_missing <= self.n_full // 4 and
            self.correction_workspace_bytes <= correction_memory_budget())
        self.use_grid = self.balanced or self.missing_correction
        self.cells = cells
        self.compact = compact
        self.kernel = kernel[self.gids[:, None], self.gids[None, :]]
        self.pairs = None
        self.D = self.U = None
        if self.use_grid:
            present = torch.zeros(self.n_full, dtype=torch.bool, device=gid.device)
            present[cells] = True
            self.missing = torch.nonzero(~present).reshape(-1)
            if eigen is not None and self.n_gid == kernel.shape[0]:
                self.D, self.U = eigen
            else:
                self.D, self.U = torch.linalg.eigh(self.kernel)
        else:
            self.pairs = kernel[gid[:, None], gid[None, :]]
            self.same_gid = (gid[:, None] == gid[None, :]).to(kernel.dtype)
        self.method = ('balanced_kronecker' if self.balanced else
                       'missing_cell_correction' if self.missing_correction else 'dense_observed')

    def grid(self, rhs):
        # Zero-fill missing cells, preserving the supplied observation order.
        full = rhs.new_zeros((self.n_full, rhs.shape[1])).index_copy(0, self.cells, rhs)
        return full.reshape(self.n_gid, self.T, -1)

    def rows(self, grid):
        return grid.reshape(self.n_full, -1)[self.cells]

    def dense_covariance(self, matrices):
        """Build the fallback covariance, including the common numerical jitter."""
        g, e = matrices
        t = self.trait
        return (self.pairs * g[t[:, None], t[None, :]] +
                self.same_gid * e[t[:, None], t[None, :]] +
                1e-6 * torch.eye(self.n_obs, dtype=g.dtype, device=g.device))

    def matmul(self, matrices, rhs):
        """Apply the linear covariance components (no numerical jitter)."""
        g, e = matrices
        if self.use_grid:
            grid = self.grid(rhs)
            kg = torch.einsum('ij,jtr->itr', self.kernel, grid)
            return self.rows(torch.einsum('ut,itr->iur', g, kg) +
                             torch.einsum('ut,itr->iur', e, grid))
        t = self.trait
        return (self.pairs * g[t[:, None], t[None, :]] +
                self.same_gid * e[t[:, None], t[None, :]]) @ rhs

    def factor(self, genetic, residual):
        return MTTrainingSolver(self, genetic, residual)


class MTTrainingSolver:
    def __init__(self, layout, genetic, residual):
        self.layout = layout
        if layout.use_grid:
            blocks = layout.D[:, None, None] * genetic + residual
            blocks = blocks + 1e-6 * torch.eye(layout.T, dtype=genetic.dtype, device=genetic.device)
            self.cholesky = torch.linalg.cholesky(blocks)
            if layout.missing_correction:
                # Only a small block of precision columns exists at once;
                # never allocate the full-grid-by-missing indicator matrix.
                observed, missing = [], []
                for indices in layout.missing.split(layout.missing_block_size):
                    basis = genetic.new_zeros((layout.n_full, indices.numel()))
                    basis[indices, torch.arange(indices.numel(),device=indices.device)] = 1
                    wm = self._solve_full(basis.reshape(layout.n_gid, layout.T, -1)).reshape(layout.n_full, -1)
                    observed.append(wm[layout.cells])
                    missing.append(wm[layout.missing])
                self.Wom = torch.cat(observed, dim=1)
                wmm = torch.cat(missing, dim=1)
                self.missing_cholesky = torch.linalg.cholesky(.5 * (wmm + wmm.T))
        else:
            self.cholesky = torch.linalg.cholesky(layout.dense_covariance([genetic, residual]))

    def _solve_full(self, grid):
        transformed = torch.einsum('ji,jtr->itr', self.layout.U, grid)
        solved = torch.cholesky_solve(transformed, self.cholesky)
        return torch.einsum('ji,itr->jtr', self.layout.U, solved)

    def solve(self, rhs):
        if not self.layout.use_grid:
            return torch.cholesky_solve(rhs, self.cholesky)
        full = self._solve_full(self.layout.grid(rhs)).reshape(self.layout.n_full, -1)
        result = full[self.layout.cells]
        if self.layout.missing_correction:
            correction = torch.cholesky_solve(full[self.layout.missing], self.missing_cholesky)
            result = result - self.Wom @ correction
        return result

    def logdet(self):
        value = 2 * torch.log(torch.diagonal(self.cholesky, dim1=-2, dim2=-1)).sum()
        if self.layout.missing_correction:
            value = value + 2 * torch.log(torch.diagonal(self.missing_cholesky)).sum()
        return value

    def reml(self, y, X):
        """Observed-data negative REML log likelihood (constant omitted), GLS."""
        solved = self.solve(torch.cat([y[:, None], X], dim=1))
        vy, vx = solved[:, 0], solved[:, 1:]
        lx = torch.linalg.cholesky(X.T @ vx)
        beta = torch.cholesky_solve((X.T @ vy)[:, None], lx).squeeze(1)
        residual = y - X @ beta
        quad = residual @ (vy - vx @ beta)
        return .5 * (self.logdet() + 2 * torch.log(torch.diagonal(lx)).sum() + quad), beta
