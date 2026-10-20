# The hybrid GP fixes the kernel/noise ratio lambda, so the GP scale must be
# estimated (REML profile: s2 = y'Py / (n - p)). Up to 0.20.223 it was taken
# as 1 and the posterior variances followed the kernel diagonal: test-hybrid
# SEs were 2.6 x the phenotypic SD in the scenario suite.

hybrid_fit_fixture <- function(scale = 1) {
  set.seed(11)
  n <- 30L
  Z <- matrix(rnorm(n * 20), n)
  K <- tcrossprod(scale(Z)) / 20 + diag(0.05, n)
  zero <- matrix(0, n, n)
  y <- 10 + drop(t(chol(K)) %*% rnorm(n)) + rnorm(n, sd = 0.5)
  list(y = y, kernels = list(total = scale * K, female = scale * K, male = zero, sca = zero),
       X = matrix(1, n, 1), train_idx = 1:24)
}

test_that("hybrid GP variances do not depend on the kernel's arbitrary scale", {
  a <- hybrid_fit_fixture(1)
  b <- hybrid_fit_fixture(10)
  fa <- PredictProR:::gp_hybrid_gp_fit_core(a$y, a$kernels, a$X, a$train_idx, lambda = 0.3)
  fb <- PredictProR:::gp_hybrid_gp_fit_core(b$y, b$kernels, b$X, b$train_idx, lambda = 3)
  expect_equal(fa$prediction, fb$prediction, tolerance = 1e-8)
  expect_equal(fa$latent_var, fb$latent_var, tolerance = 1e-6)
  expect_equal(fa$observed_var, fb$observed_var, tolerance = 1e-6)
  expect_equal(fa$variance_scale, 10 * fb$variance_scale, tolerance = 1e-6)
})

test_that("the profiled scale is the REML estimate y'Py / (n - p)", {
  a <- hybrid_fit_fixture(1)
  fit <- PredictProR:::gp_hybrid_gp_fit_core(a$y, a$kernels, a$X, a$train_idx, lambda = 0.3)
  ys <- (a$y - mean(a$y[a$train_idx])) / sd(a$y[a$train_idx])
  V <- a$kernels$total[a$train_idx, a$train_idx] + diag(0.3, length(a$train_idx))
  X <- a$X[a$train_idx, , drop = FALSE]
  Vi <- solve(V)
  P <- Vi - Vi %*% X %*% solve(crossprod(X, Vi %*% X), crossprod(X, Vi))
  yt <- ys[a$train_idx]
  expect_equal(fit$variance_scale, drop(crossprod(yt, P %*% yt)) / (length(yt) - 1), tolerance = 1e-6)
  # held-out rows: variance below the phenotypic variance
  test_rows <- setdiff(seq_along(a$y), a$train_idx)
  expect_true(all(sqrt(fit$latent_var[test_rows]) < 2 * sd(a$y)))
})
