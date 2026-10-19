"""Implicit MET cell means and exact fixed-effect absorption for operator GLS."""
import numpy as np
import torch


class CellMeanDesign:
    """One observed-cell indicator per row; pooled CV0 rows stay implicit.

    Storage is O(rows + observed cells), including unobserved environments.
    Materializing a dense rows-by-cells array is deliberately unsupported.
    """
    def __init__(self, env, trait, train, T, scales):
        self.env = np.asarray(env, dtype=np.int64)
        self.trait = np.asarray(trait, dtype=np.int64)
        codes = self.env*T + self.trait
        levels, counts = np.unique(codes[train], return_counts=True)
        self.cells = np.column_stack([levels//T, levels % T])
        self.counts = counts
        lookup = np.full(scales.shape[0]*T, -1, dtype=np.int64)
        lookup[levels] = np.arange(len(levels))
        self.cell = lookup[codes]
        self.row_sd = scales[self.env, self.trait]
        totals = np.bincount(levels % T, weights=counts, minlength=T)
        if np.any(totals == 0):
            raise ValueError('Every trait needs training observations for its pooled mean')
        self.pool_weights = counts*scales[levels//T, levels % T]/totals[levels % T]
        self.T = T
        self.shape = (len(env), len(levels))
        if len(levels) >= int(np.sum(train)):
            raise ValueError('GP mean design is unidentifiable or leaves no residual degrees of freedom')

    @property
    def storage_nbytes(self):
        return sum(getattr(self, name).nbytes for name in
                   ('env', 'trait', 'cells', 'counts', 'cell', 'row_sd', 'pool_weights'))

    def __getitem__(self, rows):
        out = object.__new__(type(self))
        out.__dict__ = self.__dict__.copy()
        for name in ('env', 'trait', 'cell', 'row_sd'):
            setattr(out, name, np.atleast_1d(getattr(self, name)[rows]))
        out.shape = (len(out.cell), self.shape[1])
        return out

    def __matmul__(self, beta):
        beta = np.asarray(beta)
        out = beta[np.maximum(self.cell, 0)].copy()
        missing = self.cell < 0
        for t in range(self.T):
            rows = missing & (self.trait == t)
            if rows.any():
                cols = self.cells[:, 1] == t
                pooled = self.pool_weights[cols] @ beta[cols]
                out[rows] = pooled / self.row_sd[rows]
        return out

    def row(self, index):
        out = np.zeros(self.shape[1])
        cell = self.cell[index]
        if cell >= 0:
            out[cell] = 1.
        else:
            cols = self.cells[:, 1] == self.trait[index]
            out[cols] = self.pool_weights[cols]/self.row_sd[index]
        return out

    def ols(self, y):
        if (self.cell < 0).any():
            raise ValueError('OLS may only use observed training cells')
        counts = np.bincount(self.cell, minlength=self.shape[1])
        sums = np.bincount(self.cell, weights=y, minlength=self.shape[1])
        # A balanced MoM sub-grid may omit some full-design cells; their
        # minimum-norm coefficients are zero, as in dense least squares.
        return np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)


class AbsorbedGLS:
    """Solve P_V b using Q V Q + (I-Q), with Q removing cell means.

    Q = I-X(X'X)^-1X'. On the residual space this SPD solve equals the
    usual GLS projector P_V, without forming X or solving V against its
    columns. The preconditioner Q diag(V)^-1 Q + (I-Q) is also SPD.
    """
    def __init__(self, design, V_mv, diagonal, solver, *, tol, max_iter):
        if np.any(design.cell < 0):
            raise ValueError('Absorption requires observed training cells')
        self.cell = torch.as_tensor(design.cell, device=diagonal.device, dtype=torch.long)
        self.counts = torch.bincount(self.cell, minlength=design.shape[1]).to(diagonal.dtype)
        if bool((self.counts == 0).any()):
            raise ValueError('Absorbed training mean design must have full column rank')
        self.V_mv, self.diagonal, self.solver = V_mv, diagonal, solver
        self.tol, self.max_iter = tol, max_iter

    def average(self, v):
        sums = torch.zeros_like(self.counts)
        sums.index_add_(0, self.cell, v)
        return sums/self.counts

    def project(self, v):
        return v - self.average(v)[self.cell]

    def solve(self, rhs):
        def A(v):
            q = self.project(v)
            return self.project(self.V_mv(q)) + (v-q)
        def precondition(v):
            q = self.project(v)
            return self.project(q/self.diagonal) + (v-q)
        return self.solver(A, self.project(rhs), M_inv_mv=precondition,
                           tol=self.tol, max_iter=self.max_iter)

    def beta(self, y, alpha):
        return self.average(y - self.V_mv(alpha))

    def mean_weights(self, row):
        """B' l for B=(X'X)^-1 X', and target mean row l."""
        l = torch.as_tensor(row, device=self.diagonal.device, dtype=self.diagonal.dtype)
        return (l/self.counts)[self.cell]
