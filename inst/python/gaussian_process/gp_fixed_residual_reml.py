"""One-dimensional REML for V = scale*K + R with known positive diagonal R."""
import numpy as np
from scipy.optimize import minimize_scalar, brentq


def fit_fixed_residual_scale(K, residual, y, X):
    """Fit a nonnegative kernel multiplier using one whitened eigendecomposition.

    The residual covariance is absolute and held fixed. Fixed effects are
    eliminated at each candidate scale; no phenotype-variance substitution or
    added ridge is used. The zero-variance endpoint is compared explicitly.
    """
    K, residual, y, X = [np.asarray(a, dtype=np.float64) for a in (K, residual, y, X)]
    if np.any(~np.isfinite(residual)) or np.any(residual <= 0):
        raise ValueError('Fixed residual variances must be finite and positive')
    if len(y) <= X.shape[1] or (X.shape[1] and np.linalg.matrix_rank(X) != X.shape[1]):
        raise ValueError('Fixed-residual REML requires full-rank means and residual degrees of freedom')
    root = np.sqrt(residual)
    W = (K + K.T)/2 / root[:,None] / root[None,:]
    Q, _ = np.linalg.qr(X/root[:,None],mode='reduced')
    residual_kernel_trace = np.trace(W)-np.trace(Q.T@W@Q)
    if residual_kernel_trace <= 1e-12*max(np.finfo(float).tiny,abs(float(np.trace(W)))):
        raise ValueError('Genetic scale is unidentifiable after removing fixed means')
    eig, U = np.linalg.eigh(W)
    if eig.min() < -1e-8*max(np.finfo(float).tiny,float(np.max(np.abs(eig)))):
        raise ValueError('Fixed-residual REML requires a positive semidefinite kernel')
    eig = np.maximum(eig, 0.)
    yw = U.T@(y/root); Xw = U.T@(X/root[:,None])
    if eig.max() <= 0.:
        raise ValueError('Genetic scale is unidentifiable for a zero kernel')
    # A response-based unit only parameterizes the search; it is not an estimate.
    ols = y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    unit = max(float(ols@ols/len(y)),float(np.mean(residual)))/max(float(np.mean(np.diag(K))),1e-30)

    def objective(scale):
        d = 1. + scale*eig
        VX = Xw/d[:,None]
        information = Xw.T@VX
        sign, logdet = np.linalg.slogdet(information)
        if sign <= 0:
            return np.inf
        beta = np.linalg.solve(information, Xw.T@(yw/d))
        r = yw-Xw@beta
        return .5*(np.log(d).sum()+logdet+np.dot(r/d,r))

    def score(scale):
        d = 1. + scale*eig
        information = Xw.T@(Xw/d[:,None])
        beta = np.linalg.solve(information,Xw.T@(yw/d))
        residual_score = (yw-Xw@beta)/d
        correction = np.trace(np.linalg.solve(information,Xw.T@((eig/d**2)[:,None]*Xw)))
        return .5*(np.sum(eig/d)-correction-np.dot(eig,residual_score**2))

    # Search over many orders of magnitude, then refine each local minimum.
    # Expansion avoids silently accepting the finite search endpoint.
    grid = np.r_[0., unit*np.logspace(-10,6,81)]
    values = np.array([objective(s) for s in grid])
    for _ in range(12):
        if values[-1] > values[-2]:
            break
        grid = np.r_[grid, grid[-1]*10.]
        values = np.r_[values, objective(grid[-1])]
    else:
        raise RuntimeError('Fixed-residual REML scale search did not bracket a minimum')
    candidates = [(values[0],0.)]
    for i in range(1,len(grid)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            result = minimize_scalar(lambda z: objective(unit*np.expm1(z)),
                bounds=(np.log1p(grid[i-1]/unit),np.log1p(grid[i+1]/unit)),
                method='bounded',options={'xatol':1e-12})
            if not result.success:
                raise RuntimeError('Fixed-residual REML scalar minimization failed')
            fitted = unit*np.expm1(result.x)
            lo, hi = grid[i-1], grid[i+1]
            if score(lo) < 0 < score(hi):
                # Objective subtraction loses precision near a flat optimum;
                # finish with the exact REML score at the same local minimum.
                fitted = unit*brentq(lambda z: score(unit*z),lo/unit,hi/unit,
                                     xtol=1e-14,rtol=1e-13)
            candidates.append((objective(fitted),fitted))
    return float(min(candidates,key=lambda pair:pair[0])[1])
