test_that("aliased kernel totals retain their joint delta-method SE in R", {
  # Independent coordinates are genetic total g1+g2 and residual e.
  C <- solve(matrix(c(4, 1, 1, 3), 2))
  x <- list(
    predictions = data.frame(GID = c("g1", "g2"), Predicted_value = c(1, 2)),
    gp_varcomp_mode = "reml", diagnostics = list(method = "gp_exact"),
    varcomp = data.frame(component = c("vm(GID,G1)!var", "vm(GID,G2)!var", "Env!R"),
      estimate = c(1, 2, 4), std.error = c(NA, NA, sqrt(C[2, 2])),
      bound = "P", pct_change = 0),
    summary = list(converged = TRUE, inference_available = FALSE,
      inference_status = "available_estimable_functions_only",
      genetic_residual_covariance = list(sigma2_g = 3, sigma2_e = 4,
        var_g = C[1, 1], var_e = C[2, 2], cov_g_e = C[1, 2]))
  )
  out <- PredictProR:::gp_format_gp_single_env_variance_components(x)
  expect_equal(out["genetic_variance", "Components"], 3)
  expect_equal(out["genetic_variance", "Standard_error"], sqrt(C[1, 1]))
  expect_equal(out["residual_variance", "Standard_error"], sqrt(C[2, 2]))
  gradient <- c(4, -3)/49
  expect_equal(out["heritability", "Standard_error"], sqrt(drop(t(gradient) %*% C %*% gradient)))
  # A covariance for a different aggregate must not be silently reused.
  x$summary$genetic_residual_covariance$sigma2_g <- 10
  mismatch <- PredictProR:::gp_format_gp_single_env_variance_components(x)
  expect_true(is.na(mismatch["genetic_variance", "Standard_error"]))
  expect_true(is.na(mismatch["heritability", "Standard_error"]))
})
