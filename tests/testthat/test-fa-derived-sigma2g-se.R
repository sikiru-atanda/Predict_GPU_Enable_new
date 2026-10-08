# Pinning test for the FA-derived sigma2_g SE delta-method path.
#
# Mathematical context: for factor-analytic mixed models, sigma2_g is NOT a
# direct REML parameter -- it is sum_k lambda_{e,k}^2 + psi_e. Its asymptotic
# SE requires delta-method propagation through (Lambda, Psi) via the AI
# inverse. The Python backend computes this in
# `varcomp_asreml.compute_fa_derived_summary` and surfaces it as a
# `std.error` column on var_components_summary's `genetic_main / total` row
# (via `gp_framework._attach_fa_derived_se_to_summary`).
#
# This test pins the R-side contract: when var_components_summary carries a
# finite std.error for the genetic_main component, the single-env reporter
# uses it for the public table's `genetic_variance` row -- even under
# strict_reml -- bypassing the varcomp aggregation that would otherwise
# return NA (because fa_psi rows are z-floor NaN'd by the Python backend).

test_that("gp_format_gp_single_env_variance_components prefers delta-method SE for FA", {
  # Synthetic FA varcomp: fa_psi rows have NA SE (z-floor flagged); only
  # fa_lambda rows carry an SE.
  varcomp <- data.frame(
    component = c(
      "fa(Env,1):vm(GID,G)!E1!var",
      "fa(Env,1):vm(GID,G)!E2!var",
      "fa(Env,1):vm(GID,G)!E1!fa1",
      "fa(Env,1):vm(GID,G)!E2!fa1",
      "Env_E1!R",
      "Env_E2!R"
    ),
    estimate = c(0.01, 0.20, 0.70, 0.80, 0.30, 0.30),
    std.error = c(NA_real_, NA_real_, 0.15, 0.18, 0.05, 0.06),
    z.ratio = c(NA_real_, NA_real_, 4.67, 4.44, 6.00, 5.00),
    bound = c("P", "P", "U", "U", "P", "P"),
    pct_change = rep(0, 6),
    stringsAsFactors = FALSE
  )

  # Synthetic var_components_summary with the principled delta-method SE
  # on the genetic_main / total row.
  var_components_summary <- data.frame(
    term_group = c("main", "main"),
    term = c("G", "G"),
    kernel = c("gmatrix", "total"),
    component_type = c("genetic_main", "genetic_main"),
    estimate = c(1.20, 1.20),
    std.error = c(NA_real_, 0.42),  # delta-method SE on the total row
    stringsAsFactors = FALSE
  )

  predictions <- data.frame(
    GID = c("g1", "g2", "g3"),
    Predicted_value = c(1.1, 1.5, 2.0),
    stringsAsFactors = FALSE
    # no Env column -> single-env path
  )

  gp_result <- list(
    predictions = predictions,
    varcomp = varcomp,
    var_components_summary = var_components_summary,
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_icm_fa")
  )

  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)

  # genetic_variance row: estimate + SE come from var_components_summary's
  # delta-method values, not from varcomp aggregation.
  expect_equal(out["genetic_variance", "Components"], 1.20)
  expect_equal(out["genetic_variance", "Standard_error"], 0.42)

  # h2_se is computable now (both genetic_se and residual_se are finite).
  expect_true(is.finite(out["heritability", "Standard_error"]))
})

test_that("gp_format_gp_single_env_variance_components falls back to varcomp when summary lacks SE", {
  # var_components_summary has no std.error column at all -> existing
  # varcomp-based path should still apply (non-FA case, REML kernel).
  varcomp <- data.frame(
    component = c("vm(GID,G)!var", "Env!R"),
    estimate = c(0.65, 0.30),
    std.error = c(0.12, 0.05),
    z.ratio = c(5.42, 6.00),
    bound = c("P", "P"),
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
    var_components_summary = data.frame(
      term_group = "main",
      term = "G",
      kernel = "total",
      component_type = "genetic_main",
      estimate = 0.650003,
      stringsAsFactors = FALSE
    ),
    gp_varcomp_mode = "reml",
    diagnostics = list(method = "gp_exact")
  )

  out <- PredictProR:::gp_format_gp_single_env_variance_components(gp_result)
  expect_equal(out["genetic_variance", "Components"], 0.65)
  expect_equal(out["genetic_variance", "Standard_error"], 0.12)
  expect_equal(
    unique(PredictProR:::gp_extract_gp_genetic_variance(
      gp_result, n = 2L, env_values = rep("ENV1", 2L)
    )),
    out["genetic_variance", "Components"]
  )
})
