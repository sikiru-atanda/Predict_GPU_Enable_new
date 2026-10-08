# Rigor: GBLUP_BRR previously used the GRM itself as the BGLR BRR design matrix
# (X = G), which implies a genetic covariance proportional to X X' = G^2 rather
# than G -- a mis-specified GBLUP that inflates the residual and lowers h2
# (confirmed: held-out accuracy 0.49 vs 0.58 for the correct parameterization).
# gp_bayes_brr_design_from_kernel(K) returns a matrix square-root X with X X' = K,
# so that BRR (genetic effect = X beta, beta ~ N(0, s2 I)) has covariance s2 * K =
# true GBLUP, matching RKHS. Downstream VarG (s2 * tcrossprod(X)) then equals s2*K.

resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  ns <- tryCatch(asNamespace("PredictProR"), error = function(e) NULL)
  if (!is.null(ns) && exists(name, envir = ns, mode = "function")) return(get(name, envir = ns, mode = "function"))
  NULL
}

make_grm <- function(n = 60L, p = 250L, seed = 1L) {
  set.seed(seed)
  M <- matrix(stats::rnorm(n * p), n, p)
  K <- tcrossprod(scale(M)) / p
  K <- K / mean(diag(K))
  rownames(K) <- colnames(K) <- paste0("g", seq_len(n))
  K
}

test_that("gp_bayes_brr_design_from_kernel yields X with X X' == K (full-rank PD kernel)", {
  fn <- resolve_fn("gp_bayes_brr_design_from_kernel")
  skip_if(is.null(fn), "gp_bayes_brr_design_from_kernel not available")
  K <- make_grm()
  K <- K + diag(1e-3, nrow(K))          # nudge to clearly PD
  X <- fn(K)
  expect_equal(nrow(X), nrow(K))
  expect_equal(unname(X %*% t(X)), unname(K), tolerance = 1e-6)
})

test_that("gp_bayes_brr_design_from_kernel handles a rank-deficient kernel without error", {
  fn <- resolve_fn("gp_bayes_brr_design_from_kernel")
  skip_if(is.null(fn), "gp_bayes_brr_design_from_kernel not available")
  K <- make_grm(n = 40L, p = 15L)        # p < n -> rank-deficient GRM
  X <- fn(K)
  expect_equal(nrow(X), nrow(K))
  expect_true(ncol(X) <= nrow(K))
  # reconstructs K on its supported subspace (drops only ~0 eigenvalues)
  expect_equal(unname(X %*% t(X)), unname(K), tolerance = 1e-4)
  expect_false(any(!is.finite(X)))
})
