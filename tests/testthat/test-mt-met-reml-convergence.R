test_that("MET multi-trait dense REML reports convergence from the optimizer, not line-search evaluations", {
  # _fit_dense_reml_met_variance() recorded the objective at every L-BFGS
  # function evaluation (strong-Wolfe trial points included) and declared
  # convergence from the last two evaluations, which are usually a trial point
  # and an accepted point. It reported not converged regardless of the fit, so
  # MET multi-trait GP always withheld genetic and GxE trait covariances.
  py <- Sys.getenv("PREDICTPRO_GP_PYTHON")
  skip_if(!nzchar(py) || !file.exists(py), "No GP Python runtime (PREDICTPRO_GP_PYTHON)")

  set.seed(11)
  n <- 30L; p <- 60L
  gids <- sprintf("G%02d", seq_len(n)); envs <- c("E1", "E2")
  X <- matrix(stats::rbinom(n * p, 2, 0.5), n, p, dimnames = list(gids, NULL))
  Z <- scale(X, scale = FALSE)
  K <- tcrossprod(Z) / mean(diag(tcrossprod(Z))); K <- K + diag(1e-3, n)
  dimnames(K) <- list(gids, gids)
  b <- matrix(stats::rnorm(p * 2), p, 2) %*% chol(matrix(c(1, 0.5, 0.5, 1), 2))
  g <- scale(Z %*% b)
  pheno <- expand.grid(GID = gids, Env = envs, stringsAsFactors = FALSE)
  gi <- match(pheno$GID, gids)
  env_shift <- ifelse(pheno$Env == "E2", 0.5, 0)
  pheno$Trait1 <- g[gi, 1] + env_shift + stats::rnorm(nrow(pheno), sd = 0.7)
  pheno$Trait2 <- g[gi, 2] - env_shift + stats::rnorm(nrow(pheno), sd = 0.7)
  test_ids <- tail(gids, 3L)

  fit <- PredictProR::gp_multi_trait_met_model(
    pheno_data = pheno, gmatrix = K, response = c("Trait1", "Trait2"), gen_name = "GID",
    heter_groups = "Env", test_set = test_ids,
    env_similarity = matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(envs, envs)),
    return_se = FALSE, return_trait_correlations = TRUE, reml_max_iter = 200L,
    python_bin = py
  )

  expect_true(isTRUE(fit$fit$reml_converged))
  expect_lt(fit$fit$reml_n_iter, 200L)
  expect_identical(fit$result$variance_component_status, "promoted")
  gc <- as.matrix(fit$result$genetic_correlation)
  expect_true(all(is.finite(gc)))
  expect_equal(unname(diag(gc)), c(1, 1), tolerance = 1e-8)
})
