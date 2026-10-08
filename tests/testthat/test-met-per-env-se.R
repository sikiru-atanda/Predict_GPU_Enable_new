# Pinning the per-env SE plumbing in the MET reporter.
#
# Before this fix `gp_variance_components_add_per_env` hard-coded
# `Standard_error <- NA_real_` on every per-env row, so the public MET
# variance_components table showed estimates per env but no SE -- even though
# the underlying varcomp rows (`vm(GID,gmatrix)!var`, `Env_<e>!R`) had finite
# REML SEs. This test pins:
#
#  1. Compound-symmetry layout (Kernel-GBLUP style): the single sigma2_g SE
#     is shared across env rows; each env's residual SE is taken from the
#     matching `Env_<e>!R` row; h2 SE is computed via 2-term delta-method.
#  2. FA layout: per-env genetic SE and h2 SE come from the FA delta-method
#     summary (`fa_derived`) attached by `build_fa_icm_ai_varcomp`. Per-env
#     residual SE still comes from the `Env_<e>!R` varcomp row.

test_that("gp_met_varcomp_augment_per_env propagates SE for compound-symmetry layout", {
  vc <- data.frame(
    Component       = c("vm(GID,gmatrix)!var", "Env_E1!R", "Env_E2!R"),
    Components      = c(0.30, 1.20, 0.60),
    Standard_error  = c(0.08, 0.20, 0.10),
    stringsAsFactors = FALSE
  )
  rownames(vc) <- vc$Component
  out <- PredictProR:::gp_met_varcomp_augment_per_env(vc)

  # Estimates per env (sigma2_g shared, sigma2_e per env).
  expect_equal(out["genetic_variance_E1",  "Components"], 0.30)
  expect_equal(out["genetic_variance_E2",  "Components"], 0.30)
  expect_equal(out["residual_variance_E1", "Components"], 1.20)
  expect_equal(out["residual_variance_E2", "Components"], 0.60)

  # SE per env: sigma2_g SE shared; sigma2_e SE from Env_<e>!R rows.
  expect_equal(out["genetic_variance_E1",  "Standard_error"], 0.08)
  expect_equal(out["genetic_variance_E2",  "Standard_error"], 0.08)
  expect_equal(out["residual_variance_E1", "Standard_error"], 0.20)
  expect_equal(out["residual_variance_E2", "Standard_error"], 0.10)

  # h2 SE: 2-term delta-method.
  total_E1 <- 0.30 + 1.20
  d_g_E1 <- 1.20 / total_E1^2
  d_r_E1 <- -0.30 / total_E1^2
  h2_se_E1 <- sqrt((d_g_E1^2) * 0.08^2 + (d_r_E1^2) * 0.20^2)
  expect_equal(out["heritability_E1", "Standard_error"], h2_se_E1, tolerance = 1e-10)
})

test_that("gp_met_varcomp_augment_per_env uses fa_derived per-env SEs for the FA layout", {
  vc <- data.frame(
    Component = c(
      "fa(Env,1):vm(GID,G)!E1!var", "fa(Env,1):vm(GID,G)!E2!var",
      "fa(Env,1):vm(GID,G)!E1!fa1", "fa(Env,1):vm(GID,G)!E2!fa1",
      "Env_E1!R", "Env_E2!R"
    ),
    Components     = c(0.10, 0.20, 0.70, 0.80, 0.30, 0.25),
    Standard_error = c(NA_real_, NA_real_, 0.15, 0.18, 0.05, 0.04),
    stringsAsFactors = FALSE
  )
  rownames(vc) <- vc$Component

  fa_derived <- list(
    env_labels         = c("E1", "E2"),
    # Includes an effective total genomic-kernel weight of 2. The raw FA
    # psi/loading rows describe the environment structure only.
    sigma2_g_per_env   = c(2 * (0.10 + 0.70^2), 2 * (0.20 + 0.80^2)),
    sigma2_g_se_per_env = c(0.44, 0.54),
    h2_per_env         = c(1.18 / (1.18 + 0.30), 1.68 / (1.68 + 0.25)),
    h2_se_per_env      = c(0.17, 0.19)
  )

  out <- PredictProR:::gp_met_varcomp_augment_per_env(vc, fa_derived = fa_derived)

  # Per-env genetic variance includes the effective genomic-kernel scale.
  expect_equal(out["genetic_variance_E1", "Components"], 1.18, tolerance = 1e-10)
  expect_equal(out["genetic_variance_E2", "Components"], 1.68, tolerance = 1e-10)

  # Per-env genetic SE comes from fa_derived (NOT NA).
  expect_equal(out["genetic_variance_E1", "Standard_error"], 0.44)
  expect_equal(out["genetic_variance_E2", "Standard_error"], 0.54)

  # h2 SE comes from fa_derived directly (delta-method already done in Python).
  expect_equal(out["heritability_E1", "Standard_error"], 0.17)
  expect_equal(out["heritability_E2", "Standard_error"], 0.19)

  # Per-env residual SE still comes from varcomp's Env_<e>!R rows.
  expect_equal(out["residual_variance_E1", "Standard_error"], 0.05)
  expect_equal(out["residual_variance_E2", "Standard_error"], 0.04)
})

test_that("gp_met_varcomp_augment_per_env falls back to 2-term h2 SE when fa_derived absent", {
  # Same FA varcomp shape but NO fa_derived -> compound-symmetry-style 2-term
  # h2 SE; per-env genetic SE is NA because FA estimates without a delta-method
  # summary cannot have a meaningful SE.
  vc <- data.frame(
    Component = c(
      "fa(Env,1):vm(GID,G)!E1!var",
      "fa(Env,1):vm(GID,G)!E1!fa1",
      "Env_E1!R"
    ),
    Components     = c(0.10, 0.70, 0.30),
    Standard_error = c(NA_real_, 0.15, 0.05),
    stringsAsFactors = FALSE
  )
  rownames(vc) <- vc$Component
  out <- PredictProR:::gp_met_varcomp_augment_per_env(vc, fa_derived = NULL)

  expect_equal(out["genetic_variance_E1",  "Components"], 0.59, tolerance = 1e-10)
  expect_true(is.na(out["genetic_variance_E1", "Standard_error"]))
  expect_equal(out["residual_variance_E1", "Standard_error"], 0.05)
  # h2 SE without genetic SE -> NA (the 2-term formula needs both).
  expect_true(is.na(out["heritability_E1", "Standard_error"]))
})
