# Regression: GP MET variance_components must flag a REML boundary-clipped
# kernel ("bound at zero") so it is not mistaken for an omitted kernel.
#
# Background (user probe, Kernel-GBLUP / Scalable-GBLUP with gmatrix + omics):
# REML clips the redundant omics kernel weight to ~0 (boundary "B"). Phase 3.28
# keeps the value visible (vm(GID,omic1_kernel)!var = 0) but NA's its SE; a user
# could still read "0" as "kernel dropped". The optional `Boundary` column makes
# the distinction explicit, and is only present when a boundary component exists
# (so the pinned 3-column contract for boundary-free tables is unchanged).

make_gp_met_result <- function(omic_bound = "B", omic_se = NA_real_, omic_val = 0.0) {
  envs <- c("E1", "E2")
  pred <- data.frame(
    GID = rep(paste0("g", 1:6), times = 2L),
    Env = rep(envs, each = 6L),
    Predicted_value = rnorm_stub <- c(seq(-1, 1, length.out = 6L), seq(-1.2, 1.2, length.out = 6L)),
    Train_Test_Label = "Train",
    stringsAsFactors = FALSE
  )
  varcomp <- data.frame(
    Component = c("vm(GID,gmatrix)!var", "vm(GID,omic1_kernel)!var",
                  "Env_E1!R", "Env_E2!R"),
    Components = c(0.50, omic_val, 0.30, 0.40),
    Standard_error = c(0.10, omic_se, 0.05, 0.06),
    bound = c("P", omic_bound, "P", "P"),
    stringsAsFactors = FALSE
  )
  list(
    gp_result = TRUE,
    gp_varcomp_mode = "reml",
    predicted_values = pred,
    varcomp = varcomp
  )
}

test_that("GP MET varcomp flags a boundary-clipped kernel with a Boundary column", {
  x <- make_gp_met_result(omic_bound = "B", omic_se = NA_real_, omic_val = 0.0)
  out <- gp_format_gp_multi_env_variance_components(x, heter_groups = "Env")

  expect_true("Boundary" %in% names(out))
  omic_row <- out[["Component"]] == "vm(GID,omic1_kernel)!var"
  gmat_row <- out[["Component"]] == "vm(GID,gmatrix)!var"
  expect_identical(out[["Boundary"]][omic_row], "bound_at_zero")
  expect_identical(out[["Boundary"]][gmat_row], "")

  # Phase 3.28 behaviour preserved: boundary kernel keeps its value, NA SE.
  expect_equal(out[["Components"]][omic_row], 0.0, tolerance = 1e-12)
  expect_true(is.na(out[["Standard_error"]][omic_row]))
  # The well-estimated kernel keeps its SE.
  expect_equal(out[["Standard_error"]][gmat_row], 0.10, tolerance = 1e-12)
})

test_that("GP MET varcomp without a boundary keeps the pinned 3-column contract", {
  x <- make_gp_met_result(omic_bound = "P", omic_se = 0.02, omic_val = 0.08)
  out <- gp_format_gp_multi_env_variance_components(x, heter_groups = "Env")
  expect_false("Boundary" %in% names(out))
  expect_identical(names(out), gp_variance_component_columns())
})

test_that("GP MET varcomp preserves non-boundary optimizer statuses", {
  x <- make_gp_met_result(omic_bound = "F", omic_se = NA_real_, omic_val = 0.08)
  out <- gp_format_gp_multi_env_variance_components(x, heter_groups = "Env")

  omic_row <- out[["Component"]] == "vm(GID,omic1_kernel)!var"
  expect_identical(out[["Boundary"]][omic_row], "fixed")
  expect_false(identical(out[["Boundary"]][omic_row], "bound_at_zero"))
})

test_that("explicit GP boundary status survives a small nonzero clipped estimate", {
  x <- make_gp_met_result(omic_bound = "B", omic_se = NA_real_, omic_val = 1e-6)
  out <- gp_format_gp_multi_env_variance_components(x, heter_groups = "Env")

  omic_row <- out[["Component"]] == "vm(GID,omic1_kernel)!var"
  expect_identical(out[["Boundary"]][omic_row], "bound_at_zero")
  expect_equal(out[["Components"]][omic_row], 1e-6, tolerance = 1e-12)
})

test_that("GP MET: residual at zero is not estimable per environment; GxE / FA specific at zero stay 0", {
  envs <- c("E1", "E2")
  pred <- data.frame(GID = rep(paste0("g", 1:6), times = 2L), Env = rep(envs, each = 6L),
                     Predicted_value = c(seq(-1, 1, length.out = 6L), seq(-1.2, 1.2, length.out = 6L)),
                     Train_Test_Label = "Train", stringsAsFactors = FALSE)
  # lower-case raw column names, as the GP backend returns them
  fa <- list(gp_result = TRUE, gp_varcomp_mode = "reml", predicted_values = pred,
             varcomp = data.frame(
               component = c("fa(Env,1):vm(GID,G)!E1!var", "fa(Env,1):vm(GID,G)!E2!var",
                             "fa(Env,1):vm(GID,G)!E1!fa1", "fa(Env,1):vm(GID,G)!E2!fa1",
                             "Env_E1!R", "Env_E2!R"),
               estimate = c(0.40, 0.0, 0.50, 0.90, 0.30, 0.0),
               std.error = c(0.10, NA, 0.10, 0.10, 0.05, NA),
               bound = c("P", "B", "U", "U", "P", "B"), stringsAsFactors = FALSE))
  out <- gp_format_gp_multi_env_variance_components(fa, heter_groups = "Env")
  out <- gp_met_varcomp_augment_per_env(out)
  rn <- rownames(out)
  val <- stats::setNames(out$Components, rn)
  flag <- stats::setNames(out$Boundary, rn)
  # FA specific variance at zero: a valid 0
  expect_equal(unname(val["fa(Env,1):vm(GID,G)!E2!var"]), 0)
  expect_identical(unname(flag["fa(Env,1):vm(GID,G)!E2!var"]), "bound_at_zero")
  # residual at zero in E2: residual, genetic and h2 of E2 not estimable
  for (r in c("Env_E2!R", "genetic_variance_E2", "residual_variance_E2", "heritability_E2")) {
    expect_true(is.na(val[[r]]), info = r)
    expect_identical(unname(flag[r]), "not_estimable", info = r)
  }
  # E1 fully estimable
  expect_equal(unname(val["genetic_variance_E1"]), 0.50^2 + 0.40, tolerance = 1e-12)
  expect_equal(unname(val["residual_variance_E1"]), 0.30)
  expect_true(is.finite(val[["heritability_E1"]]))

  # compound symmetry: GxE (Env!var) at zero is a valid 0, genetic stays estimable
  cs <- fa
  cs$varcomp <- data.frame(component = c("vm(GID,G)!var", "Env!var", "Env_E1!R", "Env_E2!R"),
                           estimate = c(0.6, 0.0, 0.3, 0.4), std.error = c(0.1, NA, 0.05, 0.06),
                           bound = c("P", "B", "P", "P"), stringsAsFactors = FALSE)
  out <- gp_met_varcomp_augment_per_env(gp_format_gp_multi_env_variance_components(cs, heter_groups = "Env"))
  val <- stats::setNames(out$Components, rownames(out))
  expect_equal(unname(val["Env!var"]), 0)
  expect_equal(unname(val[c("genetic_variance_E1", "genetic_variance_E2")]), c(0.6, 0.6))
  expect_true(all(is.finite(val[c("heritability_E1", "heritability_E2")])))
})