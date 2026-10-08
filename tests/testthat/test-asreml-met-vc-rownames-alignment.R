# Pins the ASReml MET variance_components rowname pattern: must lead with
# the three explicitly labelled across-environment summary rows followed by
# per-env (genetic_variance_<env> /
# residual_variance_<env> / heritability_<env>) rows -- the same pattern
# the Bayesian kernel-MET path emits from
# bayes_multitrait_heter_resid.R::bayes_multitrait_env_heter_fit.
#
# Background: the cross-family user probe at 0.20.8 found the Bayes MET path
# emitted (summary rows) then (per-env rows), while ASReml corgh emitted
# only per-env rows. This test pins the now-aligned ASReml output so a
# regression in either engine surfaces visibly.

test_that("asreml_variance_components MET no-omics emits summary rows then per-env rows", {
  # Three environments, single genomic kernel (no omics).
  envs <- c("B2IR", "F5I", "B5I")
  varG_per_omics <- matrix(c(0.20, 0.30, 0.18),
                           nrow = 1L, ncol = 3L,
                           dimnames = list("Genomic", envs))
  total_varG <- varG_per_omics
  var_residual <- matrix(c(0.18, 0.10, 0.20), nrow = 3L, ncol = 1L,
                         dimnames = list(envs, "Residual"))
  hertiability <- matrix(c(0.526, 0.750, 0.474), nrow = 1L, ncol = 3L,
                         dimnames = list("h2", envs))

  res <- list(
    varG_per_omics = varG_per_omics,
    Total_genetic_var = total_varG,
    Residual_Var = var_residual,
    Heritability = hertiability,
    Heritability_SE = matrix(0.05, nrow = 1L, ncol = 3L,
                             dimnames = list("h2_se", envs))
  )

  out <- asreml_variance_components(res_var_cov_h_ve = res)

  # Summary rows must come first and must state that they are aggregates.
  expect_identical(rownames(out)[1:3],
                   c("mean_genetic_variance_across_environments",
                     "mean_residual_variance_across_environments",
                     "heritability_from_mean_variances"))

  # Per-env rows must still be present after the summary rows.
  per_env_expected <- c(
    "genetic_variance_B2IR", "residual_variance_B2IR", "heritability_B2IR",
    "genetic_variance_F5I",  "residual_variance_F5I",  "heritability_F5I",
    "genetic_variance_B5I",  "residual_variance_B5I",  "heritability_B5I"
  )
  expect_true(all(per_env_expected %in% rownames(out)),
              info = paste("Missing rows:",
                           paste(setdiff(per_env_expected, rownames(out)),
                                 collapse = ", ")))

  # Summary values must match mean(per-env) for Vg and Ve, and the
  # ratio-of-means for h2.
  expected_mean_vg <- mean(as.numeric(varG_per_omics[1L, ]))
  expected_mean_ve <- mean(as.numeric(var_residual))
  expected_h2     <- expected_mean_vg / (expected_mean_vg + expected_mean_ve)
  expect_equal(out["mean_genetic_variance_across_environments", "Components"], expected_mean_vg,
               tolerance = 1e-10)
  expect_equal(out["mean_residual_variance_across_environments", "Components"], expected_mean_ve,
               tolerance = 1e-10)
  expect_equal(out["heritability_from_mean_variances", "Components"], expected_h2,
               tolerance = 1e-10)

  # Summary SE column is NA (matches Bayes pattern).
  expect_true(is.na(out["mean_genetic_variance_across_environments", "Standard_error"]))
  expect_true(is.na(out["mean_residual_variance_across_environments", "Standard_error"]))
  expect_true(is.na(out["heritability_from_mean_variances", "Standard_error"]))

  # Per-env h2 row values are unchanged (sanity check the existing
  # extractor still works on the new layout).
  expect_equal(out["heritability_B2IR", "Components"], 0.526, tolerance = 1e-6)
  expect_equal(out["heritability_F5I",  "Components"], 0.750, tolerance = 1e-6)
  expect_equal(out["heritability_B5I",  "Components"], 0.474, tolerance = 1e-6)
})

test_that("asreml_variance_components MET single-env path is unchanged", {
  # Single environment, no omics. The prepend-summary code path is gated on
  # ncol(varG_per_omics) > 1, so single-env runs must keep the original
  # (genetic_variance / residual_variance / heritability) rownames.
  varG_per_omics <- matrix(0.25, nrow = 1L, ncol = 1L,
                           dimnames = list("Genomic", "E1"))
  var_residual <- matrix(0.15, nrow = 1L, ncol = 1L,
                         dimnames = list("E1", "Residual"))
  hertiability <- matrix(0.625, nrow = 1L, ncol = 1L,
                         dimnames = list("h2", "E1"))
  res <- list(
    varG_per_omics = varG_per_omics,
    Total_genetic_var = varG_per_omics,
    Residual_Var = var_residual,
    Heritability = hertiability,
    Heritability_SE = matrix(0.05, nrow = 1L, ncol = 1L)
  )
  out <- asreml_variance_components(res_var_cov_h_ve = res)
  expect_identical(rownames(out),
                   c("genetic_variance", "residual_variance", "heritability"))
})
