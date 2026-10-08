# Regression: Bayesian kernel-MET (RKHS / GBLUP_BRR) variance_components must
# carry per-environment rows, matching the ASReml / GP engines.
#
# Background (user probe): the BGLR MET path emitted only per-ETA-term
# `genetic_variance_<i>` rows (i counts random ETA terms -- e.g. 2 kernels x
# {GID, GID:env} = 4, NOT environments) plus a single pooled heritability, so a
# 5-environment run showed `genetic_variance_1..4` with no per-env breakdown.
# bayes_variance_components_add_per_env() appends
# genetic_variance_<env> / residual_variance_<env> / heritability_<env> using
# the fitted model's marginal genetic covariance and the BGLR residual.

test_that("bayes_variance_components_add_per_env appends per-env Vg/Ve/h2", {
  base_vc <- data.frame(
    Components = c(0.20, 0.30, 0.50, 0.40, 0.5556),
    Standard_error = c(0.05, 0.06, 0.07, 0.04, 0.01),
    row.names = c("genetic_variance_1", "genetic_variance_2",
                  "total_genetic_variance", "residual_variance", "heritability"),
    stringsAsFactors = FALSE
  )

  envs <- c("E1", "E2", "E3")
  env_labels <- rep(envs, each = 5L)
  # Observation-level diag(G) values; their per-env means are the marginal Vg.
  set_e1 <- c(0.8, 0.9, 1.0, 1.1, 1.2)
  set_e2 <- c(0.3, 0.4, 0.5, 0.6, 0.7)
  set_e3 <- c(1.6, 1.8, 2.0, 2.2, 2.4)
  marginal_vg <- c(set_e1, set_e2, set_e3)
  varE <- 0.4

  out <- bayes_variance_components_add_per_env(
    variance_components = base_vc,
    env_labels = env_labels,
    marginal_genetic_variance_per_obs = marginal_vg,
    varE = varE
  )

  # Original rows preserved.
  expect_true(all(rownames(base_vc) %in% rownames(out)))

  # Per-env rows present.
  for (e in envs) {
    expect_true(paste0("genetic_variance_", e) %in% rownames(out))
    expect_true(paste0("residual_variance_", e) %in% rownames(out))
    expect_true(paste0("heritability_", e) %in% rownames(out))
  }

  # Per-env genetic variance = mean marginal variance from diag(G).
  expect_equal(out["genetic_variance_E1", "Components"], mean(set_e1), tolerance = 1e-8)
  expect_equal(out["genetic_variance_E2", "Components"], mean(set_e2), tolerance = 1e-8)
  expect_equal(out["genetic_variance_E3", "Components"], mean(set_e3), tolerance = 1e-8)

  # Homogeneous residual -> constant per-env residual = varE.
  expect_equal(out["residual_variance_E1", "Components"], varE, tolerance = 1e-8)
  expect_equal(out["residual_variance_E3", "Components"], varE, tolerance = 1e-8)

  # Heritability = Vg_env / (Vg_env + Ve_env), and therefore env-varying.
  vg3 <- mean(set_e3)
  expect_equal(out["heritability_E3", "Components"], vg3 / (vg3 + varE), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(out["heritability_E1", "Components"],
                                out["heritability_E3", "Components"])))
})

test_that("bayes_variance_components_add_per_env is a no-op on non-MET / mismatched input", {
  base_vc <- data.frame(
    Components = c(0.25, 0.15, 0.625),
    Standard_error = c(0.05, 0.04, 0.01),
    row.names = c("genetic_variance", "residual_variance", "heritability"),
    stringsAsFactors = FALSE
  )
  # Genetic-values length mismatch with env labels -> unchanged.
  out <- bayes_variance_components_add_per_env(
    variance_components = base_vc,
    env_labels = rep(c("E1", "E2"), each = 3L),
    marginal_genetic_variance_per_obs = c(1, 2, 3),  # length 3 != 6
    varE = 0.15
  )
  expect_identical(out, base_vc)
})
