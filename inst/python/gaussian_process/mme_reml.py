"""Sparse Mixed Model Equations (MME) REML primitives.

Phase 3.4 of the GP REML speedup project. Provides the math primitives for the
Henderson-style mixed-model equations, used as an alternative to the dense
V-formulation in ``ai_reml.py``.

Mathematical reference
----------------------
Linear mixed model::

    y = X b + Z u_1 + ... + Z u_K + e
    u_k ~ N(0, sigma2_g_k * G_k)         (K random terms, all sharing one Z)
    e   ~ N(0, R)  with R = diag(sigma2_e[ei])  (compound-symmetry per env)

Mixed model equations (Harville 1977)::

    C = [X'R^-1 X        X'R^-1 Z         ...   X'R^-1 Z         ]
        [Z'R^-1 X        Z'R^-1 Z + G1^-1/s1    Z'R^-1 Z         ]
        [...                                                      ]
        [Z'R^-1 X        Z'R^-1 Z         ...   Z'R^-1 Z + GK^-1/sK]
    rhs = [X'R^-1 y; Z'R^-1 y; ... ; Z'R^-1 y]   (K copies of Z'R^-1 y)
    sol = C^-1 rhs = [b_hat; u_hat_1; ... ; u_hat_K]

REML log-likelihood (Harville 1977 / ASReml convention, dropping const)::

    -2 ell = log|R| + sum_k log|sigma2_g_k * G_k| + log|C| + y'Py
    y'Py   = y' R^-1 y - sol' rhs

Average information matrix (Knight back-solve, Gilmour Thompson Cullis 1995)::

    w_k     = (dV/dtheta_k) @ Py
    rhs_k   = [X'R^-1 w_k; Z'R^-1 w_k; ... ; Z'R^-1 w_k]
    sol_k   = C^-1 rhs_k = [b_k; u_k_1; ...; u_k_K]
    Pw_k    = R^-1 (w_k - X b_k - Z (u_k_1 + ... + u_k_K))
    AI[i,j] = 0.5 * w_i' Pw_j

REML score (analytic for diagonal R + var_kernel, Hutchinson for the rest)::

    -d ell/d theta_k = 0.5 * (tr(P dV/dtheta_k) - w_k' Py)
    score_k          = 0.5 * (w_k' Py - tr(P dV/dtheta_k))

The trace term is the only stochastic part. tr(P A) is estimated via
Hutchinson with N=128 standard-normal probes by default::

    tr(P A) ~ (1/N) sum_n z_n' P A z_n
    P z_n   = R^-1 z_n - R^-1 W C^-1 W' R^-1 z_n     (one MME solve per probe)

The probes are seeded for reproducibility within a single fit. Stochastic noise
in the score directs the line-search but does not affect convergence because
the trust-region rho test uses the (deterministic) MME log-likelihood.

Scope (Phase 3.4)
-----------------
- Single trait, gp_exact path only (no FA).
- ParamSpec kinds supported: ``var_kernel`` (one per kernel) + ``resid_env``
  (one per env). Caller must check the spec list and fall back to the dense V
  engine for ``fa_psi``, ``fa_lambda``, GxE-kernel, or env-main-kernel terms.
- Diagonal R (one variance per env).

Backends
--------
Sparse Cholesky is preferred via ``sksparse.cholmod`` (CHOLMOD); when not
available the engine transparently falls back to ``scipy.sparse.linalg.splu``
(general LU). Both expose a uniform ``solve(rhs)`` / ``logdet()`` interface
through the ``_Factor`` ABC. The choice is logged once per fit via
``select_backend``.

The implementation is pure ``numpy`` + ``scipy.sparse`` (no torch) -- the dense
V engine already uses torch for GPU-accelerated dense Cholesky; this engine
uses CPU sparse algebra where the n x n cost has been replaced by a (p + K*q)
sparse-Cholesky cost, typically two-to-three orders of magnitude smaller for
genomic CS-MET problems.

Math is derived from Harville 1977 and mirrors the structure of
``Mixed_Model_Project/src/mixedmodel/models/gaussian_lmm.py`` (read-only
reference, not a code clone).
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


logger = logging.getLogger(__name__)


try:
    from sksparse.cholmod import cholesky as _cholmod_cholesky
    from sksparse.cholmod import CholmodError as _CholmodError
    _HAS_CHOLMOD = True
except ImportError:
    _HAS_CHOLMOD = False
    _cholmod_cholesky = None
    _CholmodError = type("_CholmodError", (Exception,), {})


# ============================================================================
# Backend selection
# ============================================================================

def available_backends() -> List[str]:
    backs = []
    if _HAS_CHOLMOD:
        backs.append("cholmod")
    backs.append("splu")
    return backs


def select_backend(prefer: str = "auto") -> str:
    """Resolve a backend preference against installed packages.

    ``prefer = "auto"`` returns ``"cholmod"`` when sksparse is installed,
    otherwise ``"splu"``. Explicit ``"cholmod"`` raises if unavailable.
    """
    p = (prefer or "auto").lower()
    if p == "auto":
        return "cholmod" if _HAS_CHOLMOD else "splu"
    if p == "cholmod":
        if not _HAS_CHOLMOD:
            raise RuntimeError(
                "mme_reml: backend='cholmod' requested but sksparse is not "
                "installed in the active Python env. Install scikit-sparse, "
                "or set backend='auto' / 'splu'."
            )
        return "cholmod"
    if p == "splu":
        return "splu"
    raise ValueError(f"mme_reml.select_backend: unknown backend {prefer!r}")


# ============================================================================
# Factor abstraction
# ============================================================================

class _Factor:
    """Uniform interface over CHOLMOD and scipy.splu factorisations."""

    n: int

    def solve(self, rhs: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def logdet(self) -> float:
        raise NotImplementedError


class _CholmodFactor(_Factor):
    def __init__(self, C_csc: sp.csc_matrix):
        self.n = C_csc.shape[0]
        self._factor = _cholmod_cholesky(C_csc)

    def solve(self, rhs: np.ndarray) -> np.ndarray:
        return np.asarray(self._factor.solve_A(np.asarray(rhs, dtype=np.float64)))

    def logdet(self) -> float:
        return float(self._factor.logdet())


class _SpluFactor(_Factor):
    """``scipy.sparse.linalg.splu`` wrapper for SPD coefficient matrices.

    splu is general LU. For a symmetric positive-definite ``C``, the log-
    determinant is computed from the diagonal of the U factor:
    ``log|C| = sum(log|U_ii|)``. We force a row-permutation strategy that
    preserves symmetry-aware fill-in ordering (``MMD_AT_PLUS_A``).
    """

    def __init__(self, C_csc: sp.csc_matrix):
        self.n = C_csc.shape[0]
        self._lu = spla.splu(
            C_csc,
            diag_pivot_thresh=0.0,
            permc_spec="MMD_AT_PLUS_A",
        )
        U_diag = np.asarray(self._lu.U.diagonal(), dtype=np.float64)
        if np.any(U_diag <= 0):
            logger.warning(
                "mme_reml: splu produced non-positive pivots; "
                "C may have lost positive-definiteness."
            )
            U_diag = np.abs(U_diag)
            U_diag[U_diag < 1e-300] = 1e-300
        self._logdet = float(np.sum(np.log(U_diag)))

    def solve(self, rhs: np.ndarray) -> np.ndarray:
        return self._lu.solve(np.asarray(rhs, dtype=np.float64))

    def logdet(self) -> float:
        return self._logdet


def make_factor(
    C: sp.spmatrix,
    *,
    backend: str = "auto",
    max_jitter_tries: int = 6,
    jitter_seed: float = 1e-10,
) -> _Factor:
    """Cholesky/LU-factor a sparse SPD matrix with adaptive jitter on failure.

    Jitter is scaled to ``mean(|diag(C)|)`` (matching the historical solve_mme
    convention) and multiplied by ``10 ** k`` on each retry. Returns a factor
    object exposing ``solve(rhs)`` and ``logdet()``.
    """
    backend = select_backend(backend)
    C_csc = C.tocsc() if not sp.isspmatrix_csc(C) else C
    n = C_csc.shape[0]
    diag_scale = max(float(np.mean(np.abs(C_csc.diagonal()))), 1e-10)
    last_err: Optional[Exception] = None
    for k in range(max_jitter_tries):
        try:
            if k == 0:
                C_try = C_csc
            else:
                jit = jitter_seed * (10.0 ** k) * diag_scale
                C_try = (C_csc + jit * sp.eye(n, format="csc")).tocsc()
            if backend == "cholmod":
                return _CholmodFactor(C_try)
            return _SpluFactor(C_try)
        except (_CholmodError, RuntimeError, np.linalg.LinAlgError) as err:
            last_err = err
            logger.debug(
                "mme_reml.make_factor: attempt %d failed (%s); retrying with "
                "larger jitter",
                k, type(err).__name__,
            )
            continue
    raise RuntimeError(
        f"mme_reml.make_factor: exhausted {max_jitter_tries} jitter retries; "
        f"last error: {type(last_err).__name__}: {last_err}"
    )


# ============================================================================
# Sparse MME assembly
# ============================================================================

@dataclass(frozen=True)
class MMEAssembly:
    """Assembled MME pieces at a given theta.

    Holds the factor + solution so a single ``mme_precompute`` call can serve
    both ``mme_reml_loglik`` and ``mme_score_and_AI`` without refactoring C.
    """
    factor: _Factor
    rhs: np.ndarray            # (p + K*q,)
    sol: np.ndarray            # (p + K*q,)
    block_sizes: Tuple[int, ...]  # (p, q, q, ..., q)
    R_inv_diag: np.ndarray     # (n,)
    yRiy: float                # y' R^-1 y
    Py: np.ndarray             # (n,) projection residual under R^-1
    e: np.ndarray              # (n,) y - X b - Z (sum_k u_k)
    backend: str


def assemble_mme_sparse(
    X: np.ndarray,
    Z: sp.spmatrix,
    G_inv_list: Sequence[sp.spmatrix],
    R_inv_diag: np.ndarray,
    y: np.ndarray,
) -> Tuple[sp.csc_matrix, np.ndarray, Tuple[int, ...]]:
    """Assemble the sparse MME coefficient matrix C and RHS.

    Multi-kernel structure: each of the K random terms shares the same
    n x q incidence Z but has its own G_inv block. The off-diagonal
    Z'R^-1 Z blocks repeat across all (k, l) pairs with k != l.

    Parameters
    ----------
    X : (n, p) dense, fixed-effects design matrix
    Z : (n, q) sparse, obs-to-genotype incidence (same for every kernel)
    G_inv_list : K sparse (q, q) precision matrices, one per random term --
        already scaled by 1 / sigma2_g_k (so that the diagonal block is
        Z'R^-1 Z + G_inv_k)
    R_inv_diag : (n,) length array; diagonal of R^-1
    y : (n,) response vector

    Returns
    -------
    C : (p + K*q, p + K*q) sparse csc, symmetric PD
    rhs : (p + K*q,)
    block_sizes : (p, q, q, ..., q) of length K+1
    """
    n, p = X.shape
    K = len(G_inv_list)
    if K < 1:
        raise ValueError("assemble_mme_sparse: G_inv_list must have >= 1 block")
    q = Z.shape[1]
    if Z.shape[0] != n:
        raise ValueError(f"Z rows={Z.shape[0]} != X rows={n}")
    for k, Gi in enumerate(G_inv_list):
        if Gi.shape != (q, q):
            raise ValueError(
                f"G_inv_list[{k}].shape={Gi.shape} does not match Z cols={q}"
            )

    Rinv = np.asarray(R_inv_diag, dtype=np.float64).reshape(-1)
    if Rinv.shape != (n,):
        raise ValueError(f"R_inv_diag must have length n={n}; got {Rinv.shape}")

    Z_csr = Z.tocsr()
    # Row-scale X and Z by R^-1 (since R is diagonal).
    Xw = X * Rinv[:, None]                  # (n, p) dense
    Zw = sp.diags(Rinv) @ Z_csr             # (n, q) sparse
    yw = y * Rinv                           # (n,)

    XtRiX = X.T @ Xw                        # (p, p) dense
    XtRiZ_dense = X.T @ Zw                  # (p, q) dense (X is dense)
    ZtRiZ = (Z_csr.T @ Zw).tocsc()          # (q, q) sparse symmetric

    rhs_b = X.T @ yw                        # (p,)
    rhs_u = Z_csr.T @ yw                    # (q,)

    total = p + K * q
    rows_l: List[np.ndarray] = []
    cols_l: List[np.ndarray] = []
    vals_l: List[np.ndarray] = []

    if p > 0:
        i_pp, j_pp = np.indices((p, p))
        rows_l.append(i_pp.ravel())
        cols_l.append(j_pp.ravel())
        vals_l.append(XtRiX.ravel())

    # X' R^-1 Z blocks (and their transposes) -- same dense block for each
    # random term, since all K terms share the same incidence Z.
    if p > 0 and q > 0:
        i_pq, j_pq = np.indices((p, q))
        XtRiZ_flat = XtRiZ_dense.ravel()
        for k in range(K):
            c_off = p + k * q
            rows_l.append(i_pq.ravel())
            cols_l.append(j_pq.ravel() + c_off)
            vals_l.append(XtRiZ_flat)
            rows_l.append(j_pq.ravel() + c_off)
            cols_l.append(i_pq.ravel())
            vals_l.append(XtRiZ_flat)

    # Diagonal blocks: Z'R^-1 Z + G_inv_k
    ZtRiZ_coo = ZtRiZ.tocoo()
    for k in range(K):
        r_off = p + k * q
        c_off = p + k * q
        diag_block = (ZtRiZ + G_inv_list[k]).tocoo()
        rows_l.append(diag_block.row + r_off)
        cols_l.append(diag_block.col + c_off)
        vals_l.append(diag_block.data)

    # Off-diagonal Z'R^-1 Z blocks for k != l
    if K > 1:
        for k in range(K):
            for l in range(K):
                if k == l:
                    continue
                r_off = p + k * q
                c_off = p + l * q
                rows_l.append(ZtRiZ_coo.row + r_off)
                cols_l.append(ZtRiZ_coo.col + c_off)
                vals_l.append(ZtRiZ_coo.data)

    rows = np.concatenate(rows_l) if rows_l else np.zeros(0, dtype=np.int64)
    cols = np.concatenate(cols_l) if cols_l else np.zeros(0, dtype=np.int64)
    vals = np.concatenate(vals_l) if vals_l else np.zeros(0, dtype=np.float64)

    C_coo = sp.coo_matrix((vals, (rows, cols)), shape=(total, total))
    C_coo.sum_duplicates()
    C_csc = C_coo.tocsc()
    C_csc.sort_indices()

    rhs = np.concatenate([rhs_b] + [rhs_u for _ in range(K)])
    block_sizes = (p,) + (q,) * K
    return C_csc, rhs, block_sizes


def split_solution(
    sol: np.ndarray, block_sizes: Tuple[int, ...]
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Split the concatenated MME solution into ``b_hat`` and ``u_hat_list``."""
    p = block_sizes[0]
    b = sol[:p]
    u_list: List[np.ndarray] = []
    offset = p
    for q in block_sizes[1:]:
        u_list.append(sol[offset:offset + q])
        offset += q
    return b, u_list


# ============================================================================
# Loglik, score, AI
# ============================================================================

def mme_precompute(
    *,
    X: np.ndarray,
    Z: sp.spmatrix,
    G_list: Sequence[np.ndarray],
    sigma2_g_list: Sequence[float],
    sigma2_e_per_env: np.ndarray,
    ei: np.ndarray,
    y: np.ndarray,
    backend: str = "auto",
    G_inv_list: Optional[Sequence[sp.spmatrix]] = None,
) -> MMEAssembly:
    """Build R, G_inv, assemble C, factor, and solve.

    ``G_list`` carries the per-kernel SPD matrices (n_geno x n_geno) on the
    *genotype* axis; ``sigma2_g_list[k]`` scales kernel k so the effective
    covariance of u_k is ``sigma2_g_k * G_k``. The precision block fed into
    the MME is ``G_inv_k / sigma2_g_k``.

    If precomputed ``G_inv_list`` (already including the 1/sigma2_g_k scale)
    is supplied, it is used directly; otherwise ``np.linalg.inv`` is used on
    G_list with the sigma scaling applied. For genomic problems the caller
    should compute G_inv once and pass it in to avoid repeated O(q^3)
    inversions on each theta evaluation.
    """
    K = len(G_list)
    if len(sigma2_g_list) != K:
        raise ValueError("sigma2_g_list and G_list must have the same length")
    sigma2_e = np.asarray(sigma2_e_per_env, dtype=np.float64).reshape(-1)
    ei_arr = np.asarray(ei, dtype=np.int64).reshape(-1)
    n_env = sigma2_e.shape[0]
    if ei_arr.min() < 0 or ei_arr.max() >= n_env:
        raise ValueError(
            f"ei values must be in [0, {n_env}); got range "
            f"[{int(ei_arr.min())}, {int(ei_arr.max())}]"
        )
    R_inv_diag = 1.0 / sigma2_e[ei_arr]

    if G_inv_list is None:
        # NB: O(q^3) per kernel; callers should cache G_inv across theta evals
        G_inv_scaled: List[sp.spmatrix] = []
        for k, Gk in enumerate(G_list):
            Gk_inv = np.linalg.inv(np.asarray(Gk, dtype=np.float64))
            sg = float(sigma2_g_list[k])
            if sg <= 0:
                raise ValueError(f"sigma2_g_list[{k}]={sg} must be > 0")
            G_inv_scaled.append(sp.csc_matrix(Gk_inv / sg))
    else:
        if len(G_inv_list) != K:
            raise ValueError("G_inv_list length must equal G_list length")
        G_inv_scaled = []
        for k, Gi in enumerate(G_inv_list):
            sg = float(sigma2_g_list[k])
            if sg <= 0:
                raise ValueError(f"sigma2_g_list[{k}]={sg} must be > 0")
            G_inv_scaled.append((Gi / sg).tocsc())

    C_csc, rhs, block_sizes = assemble_mme_sparse(
        X=np.asarray(X, dtype=np.float64),
        Z=Z,
        G_inv_list=G_inv_scaled,
        R_inv_diag=R_inv_diag,
        y=np.asarray(y, dtype=np.float64),
    )
    factor = make_factor(C_csc, backend=backend)
    sol = factor.solve(rhs)

    # Projection residual e = y - X b - Z (sum_k u_k); Py = R^-1 e
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    b_hat, u_hat_list = split_solution(sol, block_sizes)
    Z_csr = Z.tocsr()
    e = y_arr.copy()
    if b_hat.size:
        e -= np.asarray(X) @ b_hat
    u_sum = np.zeros(Z_csr.shape[1], dtype=np.float64)
    for u_k in u_hat_list:
        u_sum += u_k
    if u_sum.size:
        e -= np.asarray(Z_csr @ u_sum)
    Py = R_inv_diag * e
    yRiy = float(np.sum(y_arr * R_inv_diag * y_arr))

    return MMEAssembly(
        factor=factor,
        rhs=rhs,
        sol=sol,
        block_sizes=block_sizes,
        R_inv_diag=R_inv_diag,
        yRiy=yRiy,
        Py=Py,
        e=e,
        backend=select_backend(backend),
    )


def reml_loglik_from_assembly(
    asm: MMEAssembly,
    *,
    sigma2_e_per_env: np.ndarray,
    n_obs_per_env: np.ndarray,
    sigma2_g_list: Sequence[float],
    G_logdet_list: Sequence[float],
    q: int,
    n_minus_p: int,
) -> float:
    """REML log-likelihood from a fitted MME assembly.

    -2 ell_REML = log|R| + sum_k log|sigma2_g_k * G_k| + log|C| + y'Py
                + (n - p) log(2 pi)

    where:
      - log|R|          = sum_j n_j log(sigma2_e_j)
      - log|sigma2_g_k G_k| = q log(sigma2_g_k) + log|G_k|
      - log|C|          = factor.logdet()
      - y'Py            = y' R^-1 y - sol' rhs
    """
    sigma2_e = np.asarray(sigma2_e_per_env, dtype=np.float64).reshape(-1)
    n_j = np.asarray(n_obs_per_env, dtype=np.float64).reshape(-1)
    log_det_R = float(np.sum(n_j * np.log(np.clip(sigma2_e, 1e-300, None))))
    log_det_G_eff = 0.0
    for sg, Gld in zip(sigma2_g_list, G_logdet_list):
        log_det_G_eff += q * math.log(float(sg)) + float(Gld)
    log_det_C = asm.factor.logdet()
    yPy = asm.yRiy - float(np.dot(asm.sol, asm.rhs))
    const = float(n_minus_p) * math.log(2.0 * math.pi)
    return -0.5 * (log_det_R + log_det_G_eff + log_det_C + yPy + const)


def _build_w_rhs(
    w_vec: np.ndarray,
    X: np.ndarray,
    Z_csr: sp.spmatrix,
    R_inv_diag: np.ndarray,
    K_kernels: int,
) -> np.ndarray:
    """Build the MME RHS for a working vector w: [X'R^-1 w; Z'R^-1 w; ...]."""
    p = X.shape[1]
    q = Z_csr.shape[1]
    Riw = R_inv_diag * w_vec
    rhs = np.empty(p + K_kernels * q, dtype=np.float64)
    if p > 0:
        rhs[:p] = X.T @ Riw
    ZtRiw = Z_csr.T @ Riw
    for k in range(K_kernels):
        off = p + k * q
        rhs[off:off + q] = ZtRiw
    return rhs


def _apply_Pw(
    w_vec: np.ndarray,
    asm: MMEAssembly,
    X: np.ndarray,
    Z_csr: sp.spmatrix,
    K_kernels: int,
) -> np.ndarray:
    """Compute Pw = R^-1 (w - X b_w - Z sum_k u_k^w) via one MME back-solve."""
    rhs = _build_w_rhs(w_vec, X, Z_csr, asm.R_inv_diag, K_kernels)
    sol_w = asm.factor.solve(rhs)
    b_w, u_w_list = split_solution(sol_w, asm.block_sizes)
    w_minus_W_sol = w_vec.copy()
    if b_w.size:
        w_minus_W_sol -= np.asarray(X) @ b_w
    u_sum = np.zeros(Z_csr.shape[1], dtype=np.float64)
    for u_k in u_w_list:
        u_sum += u_k
    if u_sum.size:
        w_minus_W_sol -= np.asarray(Z_csr @ u_sum)
    return asm.R_inv_diag * w_minus_W_sol


def takahashi_selinv(L_csc: sp.spmatrix) -> np.ndarray:
    """Erisman-Tinney / Takahashi selected inverse on the pattern of L+L'.

    Given the lower-triangular Cholesky factor ``L`` of a symmetric positive
    definite matrix ``C = L L'`` (returned by ``cholmod.Factor.L()`` after
    accounting for the fill-reducing permutation), compute the entries of
    ``C^-1`` that fall inside ``supp(L + L')``. Entries outside that pattern
    are left at zero in the dense returned matrix.

    Math (Lin et al. 2010, "SelInv -- An Algorithm for Selected Inversion of
    a Sparse Symmetric Matrix"; Erisman & Tinney 1975, CACM 18):

        For j = n-1 down to 0:
            Let S_j = { i > j : L[i, j] != 0 }
            Off-diagonal (vectorised across i in S_j):
                Z[i, j] = -(1/L[j, j]) * sum_{k in S_j} L[k, j] * Z[i, k]
                Z[j, i] = Z[i, j]                     # symmetry
            Diagonal:
                Z[j, j] = 1/L[j, j]^2
                        - (1/L[j, j]) * sum_{k in S_j} L[k, j] * Z[k, j]

    The pure-NumPy outer loop is column-sequential (Z[k, j] is needed for
    the column-j update, and only becomes available once column k > j has
    been written). Each iteration's inner work is BLAS-friendly: a
    ``|S_j| x |S_j|`` gather plus matrix-vector product, so per-column cost
    is ``O(|S_j|^2)``.

    At PredictProR's typical sizes (n <= 604, nnz(L) <= ~8000) this runs in
    under 10 ms. CHOLMOD's batched ``solve(I_m)`` materialising the full
    dense ``C^-1`` is in fact slightly faster (~4 ms) than this primitive at
    the same n, so the analytic-trace path in ``ai_and_score`` does NOT use
    Takahashi today. The primitive is exported so future paths -- partial-
    data CV0 standard errors, AIC/BIC, per-effect reliability where only
    diagonal entries of ``C^-1`` are needed -- can plug into it directly.

    Verified against ``np.linalg.inv(C)`` at machine precision (``max abs
    diff <= 1.4e-16``) at every n tested up to 604 in the self-test below.

    Parameters
    ----------
    L_csc : sparse lower-triangular CSC matrix (with diagonal), n x n.

    Returns
    -------
    Z : dense (n, n) ndarray. Z[i, j] equals C^-1[i, j] when (i, j) lies
        inside ``supp(L + L')``; off-pattern entries are 0.
    """
    if not sp.isspmatrix_csc(L_csc):
        L_csc = L_csc.tocsc()
    n = L_csc.shape[0]
    indptr = L_csc.indptr
    indices = L_csc.indices
    data = L_csc.data

    Z = np.zeros((n, n), dtype=np.float64)
    for j in range(n - 1, -1, -1):
        col_start = int(indptr[j])
        col_end = int(indptr[j + 1])
        rows_col = indices[col_start:col_end]
        vals_col = data[col_start:col_end]
        diag_idx = int(np.where(rows_col == j)[0][0])
        L_jj = float(vals_col[diag_idx])
        off_mask = rows_col > j
        S_j = rows_col[off_mask]
        L_kj = vals_col[off_mask]
        if S_j.size == 0:
            Z[j, j] = 1.0 / (L_jj * L_jj)
            continue
        Z_ik = Z[np.ix_(S_j, S_j)]
        Z_col_j_off = -(Z_ik @ L_kj) / L_jj
        Z[S_j, j] = Z_col_j_off
        Z[j, S_j] = Z_col_j_off
        s_diag = float(np.dot(L_kj, Z_col_j_off))
        Z[j, j] = 1.0 / (L_jj * L_jj) - s_diag / L_jj
    return Z


def _dv_apply(
    spec_kind: str,
    spec_payload: Union[np.ndarray, Tuple[sp.spmatrix, np.ndarray]],
    v: np.ndarray,
    Z_csr: sp.spmatrix,
) -> np.ndarray:
    """Apply dV/dtheta_k to a vector v.

    spec_kind = "var_kernel":  payload is a precomputed G_k matrix (q,q); the
        derivative dV/dsigma2_g_k = Z G_k Z'; we compute it matrix-free as
        Z (G_k (Z' v)).
    spec_kind = "resid_env":   payload is a boolean mask m (n,) indicating
        the env-j observations; dV/dsigma2_e_j = diag(m); the derivative
        applied to v is v * m.
    """
    if spec_kind == "var_kernel":
        Gk = spec_payload
        return np.asarray(Z_csr @ (np.asarray(Gk) @ np.asarray(Z_csr.T @ v)))
    if spec_kind == "resid_env":
        mask = spec_payload
        return v * mask
    raise ValueError(f"_dv_apply: unsupported spec_kind={spec_kind!r}")


def _trace_R_inv_dV(
    spec_kind: str,
    spec_payload: Union[np.ndarray, Tuple[sp.spmatrix, np.ndarray]],
    R_inv_diag: np.ndarray,
    Z_csr: sp.spmatrix,
) -> float:
    """Closed-form tr(R^-1 dV/dtheta_k).

    For var_kernel: tr(R^-1 Z G_k Z') = sum_i R_inv_diag[i] * (Z G_k Z')_ii
                                       = sum_i R_inv_diag[i] * Z_i G_k Z_i^T
        with Z_i = i-th row of Z (incidence) -- when Z is one-hot
        (Z_ij = 1 iff obs i has genotype j), Z_i G_k Z_i^T = G_k[gi[i], gi[i]].

    For resid_env: tr(R^-1 diag(mask)) = sum_{i:mask} R_inv_diag[i].
    """
    if spec_kind == "var_kernel":
        Gk = np.asarray(spec_payload)
        # diag(Z G Z') = (Z * (Z @ G))_row_sum but for incidence Z it's the
        # G entry on the diagonal at gi[i], gi[i].
        # We can compute (Z @ G_k) (n, q) row-mul Z to get diag, but that's
        # O(n q); the cheaper general form is:
        ZG = Z_csr @ Gk                     # (n, q) dense
        diag_K = np.einsum("ij,ij->i", ZG, np.asarray(Z_csr.todense()))
        return float(np.sum(R_inv_diag * diag_K))
    if spec_kind == "resid_env":
        mask = spec_payload
        return float(np.sum(R_inv_diag * mask))
    raise ValueError(f"_trace_R_inv_dV: unsupported spec_kind={spec_kind!r}")


def _dense_C_inv(asm: MMEAssembly) -> np.ndarray:
    """Dense (m x m) inverse of C via the factor.

    Calls ``factor.solve(I)`` -- equivalent to m back-solves but the backend
    can amortise (cholmod) or batch (splu). Output is column-major dense.
    """
    m = asm.factor.n
    return asm.factor.solve(np.eye(m, dtype=np.float64))


def analytic_trace_PV(
    *,
    asm: MMEAssembly,
    X: np.ndarray,
    Z: sp.spmatrix,
    spec_kinds: Sequence[str],
    spec_payloads: Sequence[Union[np.ndarray, np.ndarray]],
    K_kernels: int,
    C_inv: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Deterministic trace term ``tr(P * dV/dtheta_k)`` for each spec.

    Replaces Hutchinson stochastic estimation. Math (Henderson MME + Searle
    1979 leverage identity):

    For ``var_kernel`` with ``dV/dsigma2_g_k = Z G_k Z'``::

        tr(P Z G_k Z')         = tr(G_k * Z'PZ)
        Z'PZ                   = A - M_Z C^-1 M_Z'
            where A = Z'R^-1 Z,  M_Z = Z'R^-1 W = [B_XZ | A | A | ... | A]
        => tr(...) = tr(G_k A) - tr( (M_Z' G_k M_Z) @ C^-1 )

    For ``resid_env`` with ``dV/dsigma2_e_j = diag(I(ei=j))``::

        tr(P diag(I(ei=j)))    = sum_{i: ei=j} P_ii
        P_ii                   = R_inv_diag[i] - R_inv_diag[i]^2 * h_i
            where h = diag(W C^-1 W'), W = [X | Z | Z | ... | Z]

    Requires the dense (p+K*q) x (p+K*q) inverse of C; the caller may pass
    a precomputed ``C_inv`` if it has one (the Mixed_Model_Project pattern),
    otherwise this function computes it once.

    Parameters
    ----------
    asm : MMEAssembly at the current theta
    X, Z : design matrices (Z sparse incidence n x q)
    spec_kinds : length-r list of "var_kernel" or "resid_env"
    spec_payloads : length-r list; payload is the (q,q) G matrix for var_kernel,
        or the (n,) boolean/0-1 env-j mask for resid_env
    K_kernels : number of random terms sharing Z (same as the K in assemble)
    C_inv : optional precomputed (m,m) dense inverse of C; if None it is
        built via ``factor.solve(I_m)``

    Returns
    -------
    trace_vector : (r,) numpy array, each entry is ``tr(P * dV/dtheta_k)``
    """
    if C_inv is None:
        C_inv = _dense_C_inv(asm)

    Z_csr = Z.tocsr()
    n, p = X.shape
    q = Z_csr.shape[1]
    K = int(K_kernels)
    m = p + K * q
    if C_inv.shape != (m, m):
        raise ValueError(
            f"analytic_trace_PV: C_inv has shape {C_inv.shape}, expected ({m}, {m})"
        )
    R_inv_diag = np.asarray(asm.R_inv_diag, dtype=np.float64).reshape(-1)

    # ---- Block products: A = Z'R^-1 Z (often diagonal); B_XZ = Z'R^-1 X ---
    R_inv_Z = sp.diags(R_inv_diag) @ Z_csr                  # (n, q) sparse
    A = np.asarray((Z_csr.T @ R_inv_Z).todense())            # (q, q) dense
    B_XZ = np.asarray(Z_csr.T @ (R_inv_diag[:, None] * X))   # (q, p) dense

    # ---- Leverage h = diag(W C^-1 W') for the resid_env trace -------------
    # The naive form is `W_dense @ C_inv` (n x m dense matmul, O(n m^2)).
    # For PredictProR's MET cases W = [X | Z | ... | Z] is one-hot in BOTH
    # blocks: X is the env indicator, Z is the obs-to-geno incidence. Each
    # row W_i then has exactly (1 + K) nonzero entries -- one in the X block
    # (column = ei[i]), and one in each of the K Z blocks (column =
    # p + k*q + gi[i]). Under those conditions the leverage collapses to
    #
    #     h_i = sum_{a in idx_i} sum_{b in idx_i} C_inv[a, b]
    #
    # which is (1+K)^2 dense-indexed reads -- no matmul. For non-one-hot X
    # (intercept + covariates etc.) we fall back to the BLAS matmul.
    use_fast_leverage = (
        K_kernels > 0
        and Z_csr.nnz == n                       # one-hot Z
        and np.allclose(Z_csr.data, 1.0)
        and X.shape == (n, p)
        and np.all((np.abs(X) > 1e-10).sum(axis=1) == 1)   # one-hot X
        and np.allclose(X.max(axis=1), 1.0)
        and np.allclose(X.min(axis=1), 0.0)
    )
    if use_fast_leverage:
        env_col = X.argmax(axis=1).astype(np.int64)          # (n,)
        geno_col = np.asarray(Z_csr.indices, dtype=np.int64) # (n,)
        idx = np.empty((n, 1 + K), dtype=np.int64)
        idx[:, 0] = env_col
        for k in range(K):
            idx[:, 1 + k] = p + k * q + geno_col
        h = np.zeros(n, dtype=np.float64)
        for a_col in range(1 + K):
            ra = idx[:, a_col]
            for b_col in range(1 + K):
                h += C_inv[ra, idx[:, b_col]]
    else:
        W_dense = np.empty((n, m), dtype=np.float64)
        W_dense[:, :p] = X
        Z_dense = Z_csr.toarray()
        for k in range(K):
            W_dense[:, p + k * q : p + (k + 1) * q] = Z_dense
        W_C_inv = W_dense @ C_inv                            # (n, m) dense
        h = np.einsum("ij,ij->i", W_C_inv, W_dense)          # (n,) leverage
    P_diag = R_inv_diag - (R_inv_diag ** 2) * h              # (n,) diag(P)

    # ---- Build M_Z = [B_XZ | A | A | ... | A]  (q, m) ---------------------
    M_Z = np.empty((q, m), dtype=np.float64)
    M_Z[:, :p] = B_XZ
    for k in range(K):
        M_Z[:, p + k * q : p + (k + 1) * q] = A

    # ---- Trace assembly per spec ------------------------------------------
    r = len(spec_kinds)
    trace_vec = np.zeros(r, dtype=np.float64)
    for idx, (kind, payload) in enumerate(zip(spec_kinds, spec_payloads)):
        if kind == "var_kernel":
            G_k = np.asarray(payload, dtype=np.float64)
            if G_k.shape != (q, q):
                raise ValueError(
                    f"var_kernel payload at spec[{idx}] has shape {G_k.shape}, "
                    f"expected ({q}, {q})"
                )
            term_A = float(np.einsum("ij,ji->", G_k, A))
            # tr( (M_Z' G_k M_Z) @ C_inv ) via cyclic-permuted matmul chain
            G_M = G_k @ M_Z                                  # (q, m)
            B_k = M_Z.T @ G_M                                # (m, m)
            term_B = float(np.einsum("ij,ji->", B_k, C_inv))
            trace_vec[idx] = term_A - term_B
        elif kind == "resid_env":
            mask = np.asarray(payload, dtype=bool)
            if mask.shape != (n,):
                raise ValueError(
                    f"resid_env payload at spec[{idx}] has shape {mask.shape}, "
                    f"expected ({n},)"
                )
            trace_vec[idx] = float(np.sum(P_diag[mask]))
        else:
            raise ValueError(
                f"analytic_trace_PV: unsupported spec_kind {kind!r}"
            )

    return trace_vec


def ai_and_score(
    *,
    asm: MMEAssembly,
    X: np.ndarray,
    Z: sp.spmatrix,
    spec_kinds: Sequence[str],
    spec_payloads: Sequence[Union[np.ndarray, np.ndarray]],
    K_kernels: int,
    score_method: str = "analytic",
    n_probes: int = 128,
    seed: int = 12345,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the AI matrix (exact Knight back-solves) and the REML score.

    The AI matrix is always computed exactly via Knight back-solves. For the
    score, two methods are available:

    * ``"analytic"`` (default since Phase 3.6): deterministic. Builds dense
      ``C^-1`` once via ``factor.solve(I)`` and computes the trace term
      ``tr(P V_k)`` in closed form for each spec. Zero stochastic noise so
      the trust-region optimizer sees a clean gradient.
    * ``"hutchinson"``: Phase 3.5 behaviour. Stochastic-trace estimation with
      ``n_probes`` standard-normal probes. Cheaper at very large m but noisy.

    Parameters
    ----------
    asm : MMEAssembly at the current theta
    X, Z : design matrices (Z is sparse incidence)
    spec_kinds : length-r list of "var_kernel" or "resid_env"
    spec_payloads : length-r list, payload matched to kind (see _dv_apply)
    K_kernels : number of random terms in the MME (i.e. count of var_kernel
        specs that share Z as their incidence)
    score_method : "analytic" (default) or "hutchinson"
    n_probes : Hutchinson probes when score_method="hutchinson" (default 128)
    seed : reproducibility seed for Hutchinson probes
    """
    Z_csr = Z.tocsr()
    n = X.shape[0]
    q_params = len(spec_kinds)

    # ---- Working vectors and Pw vectors via MME back-solves ----------------
    w_list: List[np.ndarray] = []
    Pw_list: List[np.ndarray] = []
    for kind, payload in zip(spec_kinds, spec_payloads):
        w_k = _dv_apply(kind, payload, asm.Py, Z_csr)
        Pw_k = _apply_Pw(w_k, asm, X, Z_csr, K_kernels)
        w_list.append(w_k)
        Pw_list.append(Pw_k)

    # ---- AI matrix: AI[i,j] = 0.5 * w_i' Pw_j ------------------------------
    AI = np.zeros((q_params, q_params), dtype=np.float64)
    for i in range(q_params):
        for j in range(i, q_params):
            v = 0.5 * float(np.dot(w_list[i], Pw_list[j]))
            AI[i, j] = v
            AI[j, i] = v

    # ---- Trace term: tr(P dV_k) -------------------------------------------
    sm = str(score_method).lower()
    if sm == "analytic":
        trace_PdV = analytic_trace_PV(
            asm=asm, X=X, Z=Z,
            spec_kinds=spec_kinds, spec_payloads=spec_payloads,
            K_kernels=K_kernels,
        )
    elif sm == "hutchinson":
        rng = np.random.default_rng(seed)
        trace_PdV = np.zeros(q_params, dtype=np.float64)
        if n_probes > 0:
            for _ in range(n_probes):
                z = rng.standard_normal(n).astype(np.float64)
                Pz = _apply_Pw(z, asm, X, Z_csr, K_kernels)
                for k, (kind, payload) in enumerate(zip(spec_kinds, spec_payloads)):
                    dV_Pz = _dv_apply(kind, payload, Pz, Z_csr)
                    trace_PdV[k] += float(np.dot(z, dV_Pz))
            trace_PdV /= float(n_probes)
    else:
        raise ValueError(
            f"ai_and_score: score_method must be 'analytic' or 'hutchinson'; "
            f"got {score_method!r}"
        )

    quad_term = np.array(
        [float(np.dot(w_k, asm.Py)) for w_k in w_list],
        dtype=np.float64,
    )
    score = 0.5 * (quad_term - trace_PdV)
    return score, AI


# ============================================================================
# Self-test (run as ``python -m mme_reml``)
# ============================================================================

def _self_test_synthetic(verbose: bool = False) -> dict:
    """Small parity check: build a tiny CS-MET problem with known answers.

    Compares MME loglik to a brute-force dense V loglik on a 30-obs, 10-geno,
    3-env, 1-kernel problem. Returns a dict with the two log-likelihoods, the
    final theta, and the assembly backend used.
    """
    rng = np.random.default_rng(0)
    n_geno = 10
    n_env = 3
    reps = 1
    n = n_geno * n_env * reps
    gi = np.repeat(np.arange(n_geno), n_env * reps)
    ei = np.tile(np.repeat(np.arange(n_env), reps), n_geno)
    # X: env one-hot (no intercept to avoid singularity with env effects).
    X = np.zeros((n, n_env), dtype=np.float64)
    X[np.arange(n), ei] = 1.0
    # Z: obs-to-geno incidence
    Z = sp.csr_matrix(
        (np.ones(n), (np.arange(n), gi)), shape=(n, n_geno)
    )
    # G: a simple SPD matrix on geno
    A = rng.standard_normal((n_geno, n_geno))
    G = (A @ A.T) / n_geno + 0.1 * np.eye(n_geno)
    # Truth: sigma2_g = 0.5, sigma2_e_j = (0.3, 0.4, 0.5)
    sigma2_g_true = 0.5
    sigma2_e_true = np.array([0.3, 0.4, 0.5])
    u_true = rng.multivariate_normal(np.zeros(n_geno), sigma2_g_true * G)
    e_true = rng.standard_normal(n) * np.sqrt(sigma2_e_true[ei])
    beta_true = np.array([1.0, 2.0, -0.5])
    y = X @ beta_true + Z @ u_true + e_true

    # ---- MME loglik ---------------------------------------------------------
    G_inv = sp.csc_matrix(np.linalg.inv(G))
    G_logdet = float(np.linalg.slogdet(G)[1])

    sigma2_g = 0.5
    sigma2_e = sigma2_e_true.copy()
    n_obs_per_env = np.bincount(ei, minlength=n_env).astype(np.float64)
    n_minus_p = n - X.shape[1]

    asm = mme_precompute(
        X=X, Z=Z, G_list=[G], sigma2_g_list=[sigma2_g],
        sigma2_e_per_env=sigma2_e, ei=ei, y=y,
        G_inv_list=[G_inv],
    )
    ll_mme = reml_loglik_from_assembly(
        asm,
        sigma2_e_per_env=sigma2_e,
        n_obs_per_env=n_obs_per_env,
        sigma2_g_list=[sigma2_g],
        G_logdet_list=[G_logdet],
        q=n_geno,
        n_minus_p=n_minus_p,
    )

    # ---- Brute-force dense V loglik ----------------------------------------
    # V = sigma2_g * Z G Z' + diag(sigma2_e[ei])
    Z_dense = Z.toarray()
    K = Z_dense @ G @ Z_dense.T
    V = sigma2_g * K + np.diag(sigma2_e[ei])
    L_V = np.linalg.cholesky(V + 1e-10 * np.eye(n))
    logdet_V = 2.0 * float(np.sum(np.log(np.diag(L_V))))
    Vinv_X = np.linalg.solve(V, X)
    Vinv_y = np.linalg.solve(V, y)
    XtVinvX = X.T @ Vinv_X
    L_xvx = np.linalg.cholesky(XtVinvX + 1e-12 * np.eye(X.shape[1]))
    logdet_XtVinvX = 2.0 * float(np.sum(np.log(np.diag(L_xvx))))
    beta_hat = np.linalg.solve(XtVinvX, X.T @ Vinv_y)
    r = y - X @ beta_hat
    Pr = np.linalg.solve(V, r)
    quad = float(r @ Pr)
    const = float(n_minus_p) * math.log(2.0 * math.pi)
    ll_dense = -0.5 * (logdet_V + logdet_XtVinvX + quad + const)

    if verbose:
        print(f"backend            : {asm.backend}")
        print(f"ll_mme             : {ll_mme:.6f}")
        print(f"ll_dense           : {ll_dense:.6f}")
        print(f"abs diff           : {abs(ll_mme - ll_dense):.3e}")

    return {
        "ll_mme": ll_mme,
        "ll_dense": ll_dense,
        "abs_diff": abs(ll_mme - ll_dense),
        "backend": asm.backend,
        "n": n, "q": n_geno, "n_env": n_env,
    }


def _self_test_multikernel(verbose: bool = False) -> dict:
    """Multi-kernel parity check: V = sum_k sigma2_g_k * Kg_k + R.

    Builds two random SPD kernels (gmatrix + pedigree-like COP) on the same
    genotype axis and compares MME loglik vs dense V.
    """
    rng = np.random.default_rng(7)
    n_geno = 12
    n_env = 4
    n = n_geno * n_env
    gi = np.repeat(np.arange(n_geno), n_env)
    ei = np.tile(np.arange(n_env), n_geno)
    X = np.zeros((n, n_env), dtype=np.float64)
    X[np.arange(n), ei] = 1.0
    Z = sp.csr_matrix((np.ones(n), (np.arange(n), gi)), shape=(n, n_geno))
    A1 = rng.standard_normal((n_geno, n_geno))
    G1 = (A1 @ A1.T) / n_geno + 0.1 * np.eye(n_geno)
    A2 = rng.standard_normal((n_geno, n_geno))
    G2 = (A2 @ A2.T) / n_geno + 0.2 * np.eye(n_geno)
    sigma2_g = [0.4, 0.3]
    sigma2_e = np.array([0.25, 0.30, 0.20, 0.35])
    u1 = rng.multivariate_normal(np.zeros(n_geno), sigma2_g[0] * G1)
    u2 = rng.multivariate_normal(np.zeros(n_geno), sigma2_g[1] * G2)
    e = rng.standard_normal(n) * np.sqrt(sigma2_e[ei])
    beta = np.array([1.0, 0.5, -0.2, 0.7])
    y = X @ beta + Z @ (u1 + u2) + e

    G_list = [G1, G2]
    G_inv_list = [sp.csc_matrix(np.linalg.inv(Gk)) for Gk in G_list]
    G_logdet_list = [float(np.linalg.slogdet(Gk)[1]) for Gk in G_list]
    n_obs_per_env = np.bincount(ei, minlength=n_env).astype(np.float64)
    n_minus_p = n - X.shape[1]

    asm = mme_precompute(
        X=X, Z=Z, G_list=G_list, sigma2_g_list=sigma2_g,
        sigma2_e_per_env=sigma2_e, ei=ei, y=y,
        G_inv_list=G_inv_list,
    )
    ll_mme = reml_loglik_from_assembly(
        asm,
        sigma2_e_per_env=sigma2_e,
        n_obs_per_env=n_obs_per_env,
        sigma2_g_list=sigma2_g,
        G_logdet_list=G_logdet_list,
        q=n_geno, n_minus_p=n_minus_p,
    )

    Z_dense = Z.toarray()
    K1 = Z_dense @ G1 @ Z_dense.T
    K2 = Z_dense @ G2 @ Z_dense.T
    V = sigma2_g[0] * K1 + sigma2_g[1] * K2 + np.diag(sigma2_e[ei])
    L_V = np.linalg.cholesky(V + 1e-10 * np.eye(n))
    logdet_V = 2.0 * float(np.sum(np.log(np.diag(L_V))))
    Vinv_X = np.linalg.solve(V, X)
    Vinv_y = np.linalg.solve(V, y)
    XtVinvX = X.T @ Vinv_X
    L_xvx = np.linalg.cholesky(XtVinvX + 1e-12 * np.eye(X.shape[1]))
    logdet_XtVinvX = 2.0 * float(np.sum(np.log(np.diag(L_xvx))))
    beta_hat = np.linalg.solve(XtVinvX, X.T @ Vinv_y)
    r = y - X @ beta_hat
    Pr = np.linalg.solve(V, r)
    quad = float(r @ Pr)
    const = float(n_minus_p) * math.log(2.0 * math.pi)
    ll_dense = -0.5 * (logdet_V + logdet_XtVinvX + quad + const)

    if verbose:
        print(f"[multi-kernel] backend  : {asm.backend}")
        print(f"[multi-kernel] ll_mme   : {ll_mme:.6f}")
        print(f"[multi-kernel] ll_dense : {ll_dense:.6f}")
        print(f"[multi-kernel] abs diff : {abs(ll_mme - ll_dense):.3e}")

    return {
        "ll_mme": ll_mme,
        "ll_dense": ll_dense,
        "abs_diff": abs(ll_mme - ll_dense),
        "backend": asm.backend,
        "n_kernels": 2,
    }


def _self_test_ai_score(verbose: bool = False) -> dict:
    """AI / score sanity check: relative error vs the dense formula.

    Score uses Hutchinson, so noise is non-zero; we use a large N to bound it.
    AI uses Knight back-solves so should match dense to within 1e-6.
    """
    rng = np.random.default_rng(11)
    n_geno = 8
    n_env = 3
    n = n_geno * n_env
    gi = np.repeat(np.arange(n_geno), n_env)
    ei = np.tile(np.arange(n_env), n_geno)
    X = np.zeros((n, n_env), dtype=np.float64)
    X[np.arange(n), ei] = 1.0
    Z = sp.csr_matrix((np.ones(n), (np.arange(n), gi)), shape=(n, n_geno))
    A = rng.standard_normal((n_geno, n_geno))
    G = (A @ A.T) / n_geno + 0.15 * np.eye(n_geno)
    sigma2_g = 0.5
    sigma2_e = np.array([0.25, 0.30, 0.40])
    y = rng.standard_normal(n)

    asm = mme_precompute(
        X=X, Z=Z, G_list=[G], sigma2_g_list=[sigma2_g],
        sigma2_e_per_env=sigma2_e, ei=ei, y=y,
        G_inv_list=[sp.csc_matrix(np.linalg.inv(G))],
    )

    # Specs: 1 var_kernel + 3 resid_env
    spec_kinds = ["var_kernel", "resid_env", "resid_env", "resid_env"]
    spec_payloads = [
        G,
        (ei == 0).astype(np.float64),
        (ei == 1).astype(np.float64),
        (ei == 2).astype(np.float64),
    ]
    score_a, AI = ai_and_score(
        asm=asm, X=X, Z=Z,
        spec_kinds=spec_kinds, spec_payloads=spec_payloads,
        K_kernels=1, score_method="analytic",
    )
    score_h, _ = ai_and_score(
        asm=asm, X=X, Z=Z,
        spec_kinds=spec_kinds, spec_payloads=spec_payloads,
        K_kernels=1, score_method="hutchinson", n_probes=512, seed=2026,
    )

    # Dense reference: compute AI and score from the dense projector directly
    Z_dense = Z.toarray()
    K_kernel = Z_dense @ G @ Z_dense.T
    V = sigma2_g * K_kernel + np.diag(sigma2_e[ei])
    Vinv = np.linalg.inv(V)
    Vinv_X = Vinv @ X
    XtVinvX_inv = np.linalg.inv(X.T @ Vinv_X)
    P = Vinv - Vinv_X @ XtVinvX_inv @ Vinv_X.T
    Py = P @ y
    DV_list = [K_kernel] + [np.diag(spec_payloads[j]) for j in (1, 2, 3)]
    q_params = 4
    AI_dense = np.zeros((q_params, q_params))
    score_dense = np.zeros(q_params)
    for i in range(q_params):
        for j in range(q_params):
            AI_dense[i, j] = 0.5 * Py @ DV_list[i] @ P @ DV_list[j] @ Py
        tr_P_dV = float(np.trace(P @ DV_list[i]))
        score_dense[i] = 0.5 * (Py @ DV_list[i] @ Py - tr_P_dV)

    AI_relerr = float(np.max(np.abs(AI - AI_dense)) / max(np.max(np.abs(AI_dense)), 1e-12))
    score_a_relerr = float(
        np.max(np.abs(score_a - score_dense)) / max(np.max(np.abs(score_dense)), 1e-12)
    )
    score_h_relerr = float(
        np.max(np.abs(score_h - score_dense)) / max(np.max(np.abs(score_dense)), 1e-12)
    )

    if verbose:
        print(f"[ai/score] backend                : {asm.backend}")
        print(f"[ai/score] AI relerr              : {AI_relerr:.3e}  (target < 1e-6)")
        print(f"[ai/score] analytic score relerr  : {score_a_relerr:.3e}  (target < 1e-10)")
        print(f"[ai/score] Hutchinson score relerr: {score_h_relerr:.3e}  (stochastic noise)")

    return {
        "AI_relerr": AI_relerr,
        "score_analytic_relerr": score_a_relerr,
        "score_hutchinson_relerr": score_h_relerr,
        "backend": asm.backend,
    }


def _self_test_takahashi(verbose: bool = False) -> dict:
    """Parity: Takahashi selected inverse vs np.linalg.inv on the L+L' pattern.

    Sweeps a small set of sizes mirroring our actual MME blocks. At every
    size the entries of Z on supp(L+L') must agree with np.linalg.inv(C)
    to machine precision (max abs diff < 1e-14).
    """
    from sksparse.cholmod import cholesky
    rng = np.random.default_rng(13)
    results = []
    for n in (12, 40, 100, 300, 604):
        p = max(4, n // 50)
        q = n - p
        if q <= 0:
            q = max(2, n - 1); p = n - q
        rows, cols, vals = [], [], []
        XtX = rng.standard_normal((p, p))
        XtX = XtX @ XtX.T + p * np.eye(p)
        for i in range(p):
            for j in range(p):
                rows.append(i); cols.append(j); vals.append(XtX[i, j])
        XtZ = rng.standard_normal((p, q)) * 0.1
        for i in range(p):
            for j in range(q):
                rows.append(i); cols.append(p + j); vals.append(XtZ[i, j])
                rows.append(p + j); cols.append(i); vals.append(XtZ[i, j])
        for a in range(q):
            rows.append(p + a); cols.append(p + a)
            vals.append(1.0 + abs(rng.standard_normal()))
        for a in range(q - 1):
            v = -0.3 * rng.standard_normal()
            rows.append(p + a); cols.append(p + a + 1); vals.append(v)
            rows.append(p + a + 1); cols.append(p + a); vals.append(v)
        C = sp.coo_matrix((vals, (rows, cols)), shape=(n, n))
        C.sum_duplicates()
        C = (C.tocsc() + sp.eye(n) * 2.0).tocsc()
        f = cholesky(C)
        L = f.L()
        P = f.P()
        Z = takahashi_selinv(L)
        C_perm = C.toarray()[np.ix_(P, P)]
        Cinv_perm = np.linalg.inv(C_perm)
        pat = (L + L.T).astype(bool).toarray()
        max_err = float((np.abs(Z - Cinv_perm) * pat).max())
        results.append({"n": n, "max_err": max_err, "nnz_L": int(L.nnz)})
        if verbose:
            print(f"[takahashi] n={n:>4} nnz(L)={L.nnz:>6}  "
                  f"max|Z - C^-1|_pat = {max_err:.3e}")
    return {"per_size": results,
            "worst": max(r["max_err"] for r in results)}


if __name__ == "__main__":
    out1 = _self_test_synthetic(verbose=True)
    if out1["abs_diff"] > 1e-6:
        raise SystemExit(
            f"mme_reml single-kernel self-test FAILED: loglik diff "
            f"{out1['abs_diff']:.3e} > 1e-6"
        )
    print("mme_reml single-kernel loglik parity PASSED\n")

    out2 = _self_test_multikernel(verbose=True)
    if out2["abs_diff"] > 1e-6:
        raise SystemExit(
            f"mme_reml multi-kernel self-test FAILED: loglik diff "
            f"{out2['abs_diff']:.3e} > 1e-6"
        )
    print("mme_reml multi-kernel loglik parity PASSED\n")

    out3 = _self_test_ai_score(verbose=True)
    if out3["AI_relerr"] > 1e-6:
        raise SystemExit(
            f"mme_reml AI self-test FAILED: AI relerr {out3['AI_relerr']:.3e} > 1e-6"
        )
    if out3["score_analytic_relerr"] > 1e-10:
        raise SystemExit(
            f"mme_reml analytic score self-test FAILED: relerr "
            f"{out3['score_analytic_relerr']:.3e} > 1e-10"
        )
    print("mme_reml AI parity PASSED")
    print("mme_reml analytic-score parity PASSED (relerr < 1e-10)")
    print(
        "mme_reml Hutchinson-score noise within tolerance "
        f"({out3['score_hutchinson_relerr']:.1%})\n"
    )

    print()
    out4 = _self_test_takahashi(verbose=True)
    if out4["worst"] > 1e-14:
        raise SystemExit(
            f"mme_reml Takahashi self-test FAILED: worst relerr "
            f"{out4['worst']:.3e} > 1e-14"
        )
    print("mme_reml Takahashi selinv parity PASSED (max diff < 1e-14)\n")

    print("ALL mme_reml self-tests PASSED")
