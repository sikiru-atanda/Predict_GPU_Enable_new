"""Bounded training-only GCV for prediction, without variance estimation."""
import numpy as np


def _feature_gcv(F, y, X, grid):
    # Work in feature space, without an n-by-n kernel or residual projector.
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    rank = int(np.sum(s > max(X.shape)*np.finfo(float).eps*s[0])) if len(s) else 0
    U = U[:, :rank]
    Z = F - U @ (U.T @ F)
    yr = y - U @ (U.T @ y)
    d, W = np.linalg.eigh(Z.T @ Z)
    keep = d > max(Z.shape)*np.finfo(float).eps*max(float(d.max()), 1.)
    d, W = d[keep], W[:, keep]
    H = (Z @ W) / np.sqrt(d)[None, :]
    z = H.T @ yr
    null_rss = float(np.sum((yr - H @ z)**2))
    null_df = len(y) - rank - len(d)
    if null_df < 0 or len(y)-rank < 2:
        raise ValueError('Insufficient residual degrees of freedom for feature GCV.')
    inv = 1. / (d[:, None] + grid[None, :])
    scores = len(y) * (np.sum((z[:, None]*inv)**2, axis=0) + null_rss/grid**2) / (
        np.sum(inv,axis=0) + null_df/grid)**2
    return scores, rank


def select_operator_lambda(cache, gi, ei, y, X, train_idx, w_g, w_ge,
                           w_e, Sigma_e, S_e, lambdas, *, max_rows=512, seed=12345):
    """Exact GCV on a bounded training subset; full-data fit remains matrix-free.

    The fixed-effect space is projected out before diagonalizing the kernel.
    For residual contrasts z and eigenvalues d, GCV is
    n * sum((z/(d+lambda))**2) / sum(1/(d+lambda))**2.
    Subsampling is independent of response values. This is observation-level
    GCV, not a claim of genotype- or environment-blocked validation.
    """
    import torch
    grid = np.asarray(lambdas, dtype=float).reshape(-1)
    if not len(grid) or np.any(~np.isfinite(grid)) or np.any(grid <= 0):
        raise ValueError('Regularization candidates must be finite and positive.')
    if int(max_rows) < 3:
        raise ValueError('operator_gcv_max_rows must be at least 3.')
    idx = np.asarray(train_idx, dtype=np.int64)
    # Repeated genotype observations can have a small exact feature rank even
    # when n is large. Use all training responses in that case. Both feature
    # rank and materialized entries are bounded independently of total n.
    unique_g, inverse_g = np.unique(np.asarray(gi)[idx], return_inverse=True)
    if (len(idx) > max_rows and len(unique_g) <= max_rows
            and len(idx)*len(unique_g) <= 4_000_000
            and not np.any(w_ge) and not float(w_e)):
        Kg = np.zeros((len(unique_g),len(unique_g)))
        for r, wg in enumerate(w_g):
            phi = cache.get_rows(r,torch.as_tensor(unique_g),torch.device('cpu'),torch.float64).numpy()
            Kg += float(wg)*(phi@phi.T)
        d, W = np.linalg.eigh((Kg+Kg.T)*.5)
        if d.min() < -1e-8*max(1.,float(np.max(np.abs(d)))):
            raise ValueError('GCV kernel is not positive semidefinite.')
        F = (W*np.sqrt(np.maximum(d,0.))[None,:])[inverse_g]
        scores, rank = _feature_gcv(F,np.asarray(y,float)[idx],np.asarray(X,float)[idx],grid)
        if not np.all(np.isfinite(scores)):
            raise ValueError('Nonfinite training GCV score.')
        chosen = int(np.argmin(scores))
        return float(grid[chosen]), dict(method='training_feature_gcv',n_train=int(len(idx)),
            n_selection=int(len(idx)),max_rows=int(max_rows),feature_rank_bound=int(len(unique_g)),
            fixed_rank=rank,seed=int(seed),candidates=grid.tolist(),scores=scores.tolist(),
            selected_lambda=float(grid[chosen]),selected_at_grid_edge=chosen in (0,len(grid)-1),
            variance_estimation=False)
    if len(idx) > int(max_rows):
        idx = np.sort(np.random.default_rng(seed).choice(idx, int(max_rows), replace=False))
    xs = np.asarray(X, float)[idx]
    # SVD gives the actual fixed-effect space even when subset columns alias.
    U, s, _ = np.linalg.svd(xs, full_matrices=True)
    rank = int(np.sum(s > (max(xs.shape) * np.finfo(float).eps * s[0]))) if len(s) else 0
    Q = U[:, rank:]
    if Q.shape[1] < 2:
        raise ValueError('GCV requires at least two residual degrees of freedom.')
    g, e = np.asarray(gi)[idx], np.asarray(ei)[idx]
    K = np.zeros((len(idx), len(idx)))
    env_cov = np.asarray(Sigma_e)[np.ix_(e, e)]
    for r, (wg, wge) in enumerate(zip(w_g, w_ge)):
        phi = cache.get_rows(r, torch.as_tensor(g), torch.device('cpu'), torch.float64).numpy()
        K += (phi @ phi.T) * (float(wg) + float(wge) * env_cov)
    if float(w_e):
        K += float(w_e) * (np.eye(int(np.max(ei)) + 1) if S_e is None else S_e)[np.ix_(e, e)]
    projected = Q.T @ K @ Q
    d, Uc = np.linalg.eigh((projected + projected.T) * .5)
    if d.min() < -1e-8 * max(1., float(np.max(np.abs(d)))):
        raise ValueError('GCV kernel is not positive semidefinite.')
    d = np.maximum(d, 0.)
    z = Uc.T @ Q.T @ np.asarray(y, float)[idx]
    inv = 1. / (d[:, None] + grid[None, :])
    scores = len(idx) * np.sum((z[:, None] * inv)**2, axis=0) / np.sum(inv, axis=0)**2
    if not np.all(np.isfinite(scores)):
        raise ValueError('Nonfinite training GCV score.')
    chosen = int(np.argmin(scores))
    return float(grid[chosen]), dict(method='training_subset_gcv' if len(idx) < len(train_idx) else 'training_gcv',
        n_train=int(len(train_idx)), n_selection=int(len(idx)), max_rows=int(max_rows),
        fixed_rank=rank, seed=int(seed), candidates=grid.tolist(), scores=scores.tolist(),
        selected_lambda=float(grid[chosen]), selected_at_grid_edge=chosen in (0, len(grid)-1),
        variance_estimation=False)
