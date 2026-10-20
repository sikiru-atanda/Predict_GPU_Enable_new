# Rule agreed 2026-10-07 (0.20.185): where the residual variance is estimated
# at zero, genetic and residual variance are not separable, so reliability is
# not estimable for that trait (single environment) or environment (MET).
# 0.20.221 dropped the rule from the GP prediction table; 0.20.224 restores it.

fa_varcomp <- function() {
  data.frame(
    component = c("fa(Env,1):vm(GID,G)!E1!var", "fa(Env,1):vm(GID,G)!E2!var",
                  "fa(Env,1):vm(GID,G)!E1!fa1", "fa(Env,1):vm(GID,G)!E2!fa1",
                  "Env_E1!R", "Env_E2!R"),
    estimate = c(0.36, 0, 0.38, 0.93, 0.31, 0),
    std.error = c(0.35, NA, 0.12, 0.10, 0.64, NA),
    bound = c("P", "B", "U", "U", "P", "B"),
    stringsAsFactors = FALSE
  )
}

test_that("rows of an environment whose residual is at zero are flagged", {
  gp_result <- list(varcomp = fa_varcomp())
  env <- rep(c("E1", "E2"), each = 3)
  expect_identical(PredictProR:::gp_residual_at_zero_rows(gp_result, env),
                   rep(c(FALSE, TRUE), each = 3))
  single <- list(varcomp = data.frame(component = c("vm(GID,G)!var", "Env_ENV1!R"),
                                      estimate = c(1.2, 0.0001), std.error = c(0.6, NA),
                                      bound = c("P", "B"), stringsAsFactors = FALSE))
  expect_true(all(PredictProR:::gp_residual_at_zero_rows(single, rep("ENV1", 4))))
  ok <- list(varcomp = transform(fa_varcomp(), bound = "P", estimate = estimate + 0.1))
  expect_false(any(PredictProR:::gp_residual_at_zero_rows(ok, env)))
})

test_that("the GP prediction table reports no reliability where the residual is at zero", {
  ids <- paste0("L", 1:3)
  pheno <- data.frame(GID = rep(ids, 2), Env = rep(c("E1", "E2"), each = 3),
                      Yield = c(10, 11, NA, 12, 13, NA), stringsAsFactors = FALSE)
  preds <- data.frame(GID = pheno$GID, Env = pheno$Env,
                      Predicted_value = c(10.1, 10.9, 10.5, 12.2, 12.8, 12.4),
                      Standard_error = c(0.4, 0.4, 0.9, 0.01, 0.01, 0.9),
                      Reliability_variance_input = c(0.2, 0.2, 0.8, 0.01, 0.01, 0.8),
                      Reliability_reference_variance = c(1, 1, 1, 1.7, 1.7, 1.7),
                      Reliability_basis = "genetic_plus_gxe_prior_variance",
                      stringsAsFactors = FALSE)
  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = preds, pheno_data = pheno, response = "Yield", gen_name = "GID",
    heter_groups = "Env", gp_result = list(varcomp = fa_varcomp())
  )
  e2 <- out$Env == "E2"
  expect_true(all(is.na(out$Reliability[e2])))
  expect_true(all(grepl("_not_estimable$", out$Reliability_basis[e2])))
  expect_true(all(grepl("residual variance at zero", out$Reliability_remarks[e2])))
  expect_true(all(is.finite(out$Reliability[!e2])))
})

asreml_fake_model <- function(bound_e2 = "B") {
  vc <- data.frame(
    component = c(0.36, 0, 0.38, 0.93, 0.31, 0.0000009, 1.47),
    std.error = c(0.35, NA, 0.12, 0.10, 0.64, NA, 0.34),
    z.ratio = NA_real_,
    bound = c("P", "B", "U", "U", "P", bound_e2, "P"),
    stringsAsFactors = FALSE,
    row.names = c("fa(Env, 1):vm(GID, G)!E1!var", "fa(Env, 1):vm(GID, G)!E2!var",
                  "fa(Env, 1):vm(GID, G)!E1!fa1", "fa(Env, 1):vm(GID, G)!E2!fa1",
                  "Env_E1!R", "Env_E2!R", "Env_E3!R")
  )
  structure(list(varcomp = vc), class = "pp_fake_asreml")
}

test_that("ASReml follows the same rule in its variance table and reliability", {
  registerS3method("summary", "pp_fake_asreml", function(object, ...) object,
                   envir = asNamespace("PredictProR"))
  envs <- c("E1", "E2", "E3")
  vc <- data.frame(
    Components = c(0.5, 0.4, 0.556, 0.51, 0.31, 0.62, 0.86, 0.0000009, 0.99999, 0.12, 1.47, 0.0755),
    Standard_error = NA_real_,
    row.names = c("mean_genetic_variance_across_environments",
                  "mean_residual_variance_across_environments",
                  "heritability_from_mean_variances",
                  as.vector(t(outer(c("genetic_variance", "residual_variance", "heritability"),
                                    envs, paste, sep = "_")))[order(rep(1:3, 3))])
  )
  rownames(vc)[4:12] <- c(rbind(paste0("genetic_variance_", envs), paste0("residual_variance_", envs),
                                paste0("heritability_", envs)))
  out <- PredictProR:::asreml_flag_residual_at_zero_vc(vc, asreml_fake_model())
  e2 <- rownames(out) %in% c("genetic_variance_E2", "residual_variance_E2", "heritability_E2")
  expect_true(all(is.na(out$Components[e2])))
  expect_true(all(out$Boundary[e2] == "not_estimable"))
  expect_true(all(out$Boundary[!e2] == ""))
  expect_equal(out["mean_genetic_variance_across_environments", "Components"],
               mean(out[c("genetic_variance_E1", "genetic_variance_E3"), "Components"]))
  expect_identical(PredictProR:::asreml_residual_at_zero_rows(asreml_fake_model(), rep(envs, each = 2)),
                   rep(c(FALSE, TRUE, FALSE), each = 2))
  untouched <- PredictProR:::asreml_flag_residual_at_zero_vc(vc, asreml_fake_model(bound_e2 = "P"))
  expect_equal(untouched$Components, vc$Components)
})
