# A GP MET fit with one residual for every environment (heter_resid = FALSE)
# labels it `Env_common!R` (`Residual!R` on the KRR route). 0.20.221-0.20.223
# reported per-"environment" rows named genetic_variance_common /
# genetic_variance_Residual, dropped the G x Env row estimated at zero, and so
# lost the environment covariance/correlation tables.

common_vc <- function(resid_label = "Env_common!R", gen_label = "vm(GID,gmatrix)!var") {
  data.frame(
    Component = c(gen_label, resid_label),
    Components = c(0.27, 1.06),
    Standard_error = c(0.11, 0.17),
    Estimation_method = "REML",
    row.names = c(gen_label, resid_label),
    stringsAsFactors = FALSE
  )
}

test_that("a common residual is reported for every environment, with the GxE row", {
  out <- PredictProR:::gp_met_varcomp_augment_per_env(
    common_vc(), env_levels = c("E1", "E2", "E3"), gxe_variance = 0
  )
  rn <- rownames(out)
  expect_false(any(grepl("_common$|_Residual$", rn)))
  expect_true(all(paste0("heritability_", c("E1", "E2", "E3")) %in% rn))
  expect_equal(out["residual_variance_E2", "Components"], 1.06)
  expect_equal(out["genetic_variance_E3", "Components"], 0.27)
  expect_equal(out["heritability_E1", "Components"], 0.27 / (0.27 + 1.06))
  expect_true("Env!var" %in% rn)
  expect_identical(out["Env!var", "Boundary"], "bound_at_zero")
  expect_equal(out["Env!var", "Components"], 0)
})

test_that("the KRR combined kernel gets per-environment rows but no GxE row", {
  out <- PredictProR:::gp_met_varcomp_augment_per_env(
    common_vc("Residual!R", "vm(GID,KRR_kernel)!var"),
    env_levels = c("E1", "E2"), gxe_variance = 0
  )
  expect_false("Env!var" %in% rownames(out))
  expect_true(all(c("genetic_variance_E1", "residual_variance_E2") %in% rownames(out)))
  expect_false(any(grepl("_Residual$", rownames(out))))
})

test_that("per-environment residuals are unchanged", {
  vc <- data.frame(
    Component = c("vm(GID,gmatrix)!var", "Env!var", "Env_E1!R", "Env_E2!R"),
    Components = c(0.3, 0.1, 0.8, 1.2), Standard_error = c(0.1, 0.05, 0.2, 0.3),
    row.names = c("vm(GID,gmatrix)!var", "Env!var", "Env_E1!R", "Env_E2!R"),
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::gp_met_varcomp_augment_per_env(vc, env_levels = c("E1", "E2"), gxe_variance = 0.1)
  expect_equal(sum(rownames(out) == "Env!var"), 1L)
  expect_equal(out["genetic_variance_E1", "Components"], 0.4)
  expect_equal(out["residual_variance_E2", "Components"], 1.2)
})

test_that("environment levels and the GxE variance are read from the GP result", {
  x <- list(gp_result = list(
    env_variance_summary = data.frame(level_index = 0:2, level = c("E1", "E2", "E3")),
    var_components_summary = data.frame(
      term_group = c("main", "interaction"), term = c("G", "G:Env"),
      kernel = c("total", "total"), component_type = c("genetic_main", "interaction"),
      estimate = c(0.27, 0)
    )
  ))
  info <- PredictProR:::gp_met_varcomp_env_info(x, "Env")
  expect_identical(info$env_levels, c("E1", "E2", "E3"))
  expect_identical(info$gxe_variance, 0)
})
