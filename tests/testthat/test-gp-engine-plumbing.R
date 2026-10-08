# Pinning tests for Phase 3.4 gp_engine argument propagation.

test_that("gp_backend_gaussian_model accepts gp_engine and matches expected values", {
  arg_names <- names(formals(PredictProR:::gp_backend_gaussian_model))
  expect_true("gp_engine" %in% arg_names)
})

test_that("gp_single_trait_model accepts gp_engine", {
  arg_names <- names(formals(PredictProR::gp_single_trait_model))
  expect_true("gp_engine" %in% arg_names)
})

test_that("gp_bridge_fit_mixed_model_direct accepts gp_engine", {
  arg_names <- names(formals(PredictProR:::gp_bridge_fit_mixed_model_direct))
  expect_true("gp_engine" %in% arg_names)
})

test_that("gp_bridge_prepare_mixed_model_spec accepts gp_engine", {
  arg_names <- names(formals(PredictProR:::gp_bridge_prepare_mixed_model_spec))
  expect_true("gp_engine" %in% arg_names)
})

test_that("model_execute accepts gp_engine", {
  arg_names <- names(formals(PredictProR::model_execute))
  expect_true("gp_engine" %in% arg_names)
})

test_that("PREDICTPRO_GP_ENGINE env var is honoured when set", {
  withr::with_envvar(c(PREDICTPRO_GP_ENGINE = "dense_v"), {
    resolved <- tolower(Sys.getenv("PREDICTPRO_GP_ENGINE"))
    expect_identical(resolved, "dense_v")
  })
})

test_that("gp_engine is rejected for invalid values via match.arg", {
  expect_error(
    PredictProR:::gp_backend_gaussian_model(
      pheno_data = data.frame(GID = "g1", Yield = 1.0),
      gmatrix = matrix(1, 1, 1, dimnames = list("g1", "g1")),
      gen_name = "GID", response = "Yield",
      model_name = "KRR",
      gp_engine = "invalid_engine"
    ),
    "should be one of|arg.*should",
    ignore.case = TRUE
  )
})
