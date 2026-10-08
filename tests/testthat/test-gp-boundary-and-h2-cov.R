# Pinning tests for:
#   A. boundary handling: when only ONE component is at the REML boundary,
#      preserve the well-estimated partner's estimate + SE; report the
#      boundary value as-is with SE = NA (ASReml convention). Previous
#      behaviour wiped both components to NA.
#   B. covariance-aware h^2 SE: when var_components_summary attaches a
#      Cov(sigma2_g, sigma2_e) dict from the AI inverse, the public h^2
#      SE uses the full 3-term delta-method instead of the independence-
#      assuming 2-term approximation.

# ---------------------------------------------------------------------------
# A: boundary preservation
# ---------------------------------------------------------------------------

test_that("residual at zero makes genetic, residual and h2 not estimable", {
  # var_kernel genetic with finite estimate + SE; resid_env clipped at 0.
  varcomp <- data.frame(
    component = c("vm(GID,G)!var", "Env!R"),
    estimate = c(0.65, 0.0),
    std.error = c(0.12, NA_real_),  # resid SE NaN per ASReml at-boundary convention
    z.ratio = c(5.42, NA_real_),
    bound = c("P", "B"),  # resid at boundary
    pct_change = c(0, 0),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  gp_result <- list(
    predictions = predictions,
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact")
  )

  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)

  # Residual estimated at zero: the genetic estimate absorbs the residual, so
  # genetic variance, residual variance and h^2 are all not estimable.
  expect_true(all(is.na(out[c("genetic_variance", "residual_variance", "heritability"), "Components"])))
  expect_true(all(is.na(out[c("genetic_variance", "residual_variance", "heritability"), "Standard_error"])))
  expect_identical(unname(out[, "Boundary"]), rep("not_estimable", 3L))
})

test_that("genetic at zero is reported as 0 with h2 = 0 and a boundary flag", {
  varcomp <- data.frame(
    component = c("vm(GID,G)!var", "Env!R"),
    estimate = c(0.0, 0.30),
    std.error = c(NA_real_, 0.06),
    z.ratio = c(NA_real_, 5.00),
    bound = c("B", "P"),
    pct_change = c(0, 0),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  gp_result <- list(
    predictions = predictions,
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact")
  )

  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)

  # Genetic at zero (residual estimable) is a valid REML estimate: 0, h^2 = 0.
  expect_equal(out["genetic_variance", "Components"], 0)
  expect_true(is.na(out["genetic_variance", "Standard_error"]))
  expect_equal(out["residual_variance", "Components"], 0.30)
  expect_equal(out["residual_variance", "Standard_error"], 0.06)
  expect_equal(out["heritability", "Components"], 0)
  expect_identical(out["genetic_variance", "Boundary"], "bound_at_zero")
})

test_that("one boundary kernel does not erase an estimable multi-kernel genetic total", {
  varcomp <- data.frame(
    component = c("vm(GID,gmatrix)!var", "vm(GID,omic1)!var", "vm(GID,omic2)!var", "Env!R"),
    estimate = c(0.2668689, 0.0, 0.1445661, 0.3022461),
    std.error = c(0.08, NA_real_, 0.06, 0.05),
    z.ratio = c(3.34, NA_real_, 2.41, 6.04),
    bound = c("P", "B", "P", "P"),
    pct_change = c(0, 0, 0, 0),
    stringsAsFactors = FALSE
  )
  gp_result <- list(
    predictions = data.frame(
      GID = c("g1", "g2"),
      Predicted_value = c(1.0, 2.0),
      stringsAsFactors = FALSE
    ),
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact")
  )

  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)
  expected_genetic <- 0.2668689 + 0.1445661
  expected_h2 <- expected_genetic / (expected_genetic + 0.3022461)

  expect_equal(out["genetic_variance", "Components"], expected_genetic)
  expect_equal(out["residual_variance", "Components"], 0.3022461)
  expect_equal(out["heritability", "Components"], expected_h2)
  expect_true(is.finite(out["genetic_variance", "Standard_error"]))
})

# ---------------------------------------------------------------------------
# B: covariance-aware h^2 SE
# ---------------------------------------------------------------------------

test_that("h2 SE uses the full 3-term delta-method when Cov(sigma2_g, sigma2_e) is supplied", {
  varcomp <- data.frame(
    component = c("vm(GID,G)!var", "Env!R"),
    estimate = c(0.65, 0.30),
    std.error = c(0.12, 0.06),
    z.ratio = c(5.42, 5.00),
    bound = c("P", "P"),
    pct_change = c(0, 0),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )

  # AI inverse for (sigma2_g, sigma2_e) with NEGATIVE off-diagonal. This is a
  # very common pattern for REML (e.g. negative correlation between estimated
  # genetic and residual variances). The 3-term formula must subtract twice
  # the negative cross-term -- i.e. ADD a positive contribution -- giving a
  # larger h^2 SE than the 2-term assuming independence.
  g_r_cov <- list(
    sigma2_g = 0.65,
    sigma2_e = 0.30,
    var_g    = 0.12^2,
    var_e    = 0.06^2,
    cov_g_e  = -0.003
  )

  gp_result_with_cov <- list(
    predictions = predictions,
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact"),
    summary = list(genetic_residual_covariance = g_r_cov)
  )
  gp_result_without_cov <- list(
    predictions = predictions,
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact")
  )

  out_full <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result_with_cov)
  out_indep <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result_without_cov)

  total <- 0.65 + 0.30
  d_g <- 0.30 / total^2
  d_r <- -0.65 / total^2
  expected_2term <- sqrt((d_g^2 * 0.12^2) + (d_r^2 * 0.06^2))
  expected_3term <- sqrt(
    (d_g^2 * 0.12^2) + (d_r^2 * 0.06^2) + (2 * d_g * d_r * -0.003)
  )

  expect_equal(out_indep["heritability", "Standard_error"], expected_2term, tolerance = 1e-10)
  expect_equal(out_full["heritability", "Standard_error"], expected_3term, tolerance = 1e-10)
  # Negative cov_g_e with one positive and one negative chain derivative -> the
  # cross-term is POSITIVE -> 3-term SE is larger than the 2-term.
  expect_gt(out_full["heritability", "Standard_error"],
            out_indep["heritability", "Standard_error"])
})

test_that("h2 SE silently falls back to 2-term when Cov has NA entries", {
  varcomp <- data.frame(
    component = c("vm(GID,G)!var", "Env!R"),
    estimate = c(0.65, 0.30),
    std.error = c(0.12, 0.06),
    z.ratio = c(5.42, 5.00),
    bound = c("P", "P"),
    pct_change = c(0, 0),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )
  g_r_cov <- list(
    sigma2_g = 0.65, sigma2_e = 0.30,
    var_g = NA_real_, var_e = NA_real_, cov_g_e = NA_real_
  )
  gp_result <- list(
    predictions = predictions,
    varcomp = varcomp,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact"),
    summary = list(genetic_residual_covariance = g_r_cov)
  )
  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)

  total <- 0.65 + 0.30
  d_g <- 0.30 / total^2
  d_r <- -0.65 / total^2
  expected_2term <- sqrt((d_g^2 * 0.12^2) + (d_r^2 * 0.06^2))
  expect_equal(out["heritability", "Standard_error"], expected_2term, tolerance = 1e-10)
})
