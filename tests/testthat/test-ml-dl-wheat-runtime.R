if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

ml_dl_wheat_slice <- function(n = 24L, p = 60L) {
  data("wheat", package = "BGLR", envir = environment())
  X <- as.matrix(get("wheat.X", envir = environment()))[seq_len(n), seq_len(p), drop = FALSE]
  Y <- as.matrix(get("wheat.Y", envir = environment()))[seq_len(n), 1:2, drop = FALSE]
  rownames(X) <- paste0("w", seq_len(n))
  colnames(X) <- paste0("m", seq_len(ncol(X)))
  pheno <- data.frame(
    GID = rownames(X),
    Trait1 = Y[, 1L],
    Trait2 = Y[, 2L],
    stringsAsFactors = FALSE
  )
  pheno$Trait1[c(3L, 9L, 17L)] <- NA_real_
  pheno$Trait2[c(5L, 12L, 21L)] <- NA_real_
  list(pheno = pheno, geno = X)
}

expect_wheat_predictive_contract <- function(out) {
  pred <- out$model_results$predicted_values
  expect_identical(
    names(pred),
    PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
  )
  expect_true(all(is.finite(pred$Standard_error)))
  expect_equal(pred$PEV, pred$Standard_error^2)
  expect_true(all(is.finite(pred$Reliability)))
  expect_equal(
    pred$Reliability,
    pred$Reliability_reference_variance /
      (pred$Reliability_reference_variance + pred$PEV)
  )
  expect_equal(pred$Reliability_variance_input, pred$PEV)
  expect_true(all(is.finite(pred$Prediction_stability)))
  expect_true(all(grepl("not genetic reliability", pred$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("predictive precision index", pred$Reliability_basis, fixed = TRUE)))
  # Uncertainty is cross-fitted and calibrated to held-out predictive error;
  # the method is recorded in Prediction_uncertainty_source and PEV_basis
  # describes the calibration (and that it is not mixed-model PEV).
  expect_true(all(grepl("cross_fitted", pred$Prediction_uncertainty_source, fixed = TRUE)))
  expect_true(all(grepl("held-out predictive MSE", pred$PEV_basis, fixed = TRUE)))
  expect_true(all(grepl("not mixed-model PEV", pred$PEV_basis, fixed = TRUE)))
  expect_true(is.matrix(out$model_results$Prediction_covariance_traits))
  expect_true(is.matrix(out$model_results$Prediction_error_covariance_traits))
  expect_true(is.matrix(out$model_results$Prediction_covariance_traits_SE))
  expect_false("Genetic_covariance_traits" %in% names(out$model_results))
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "phenotypic_reference_variance",
    "prediction_variance",
    "estimated_prediction_mse"
  ) %in% out$model_results$variance_components$Component))
}

test_that("BGLR wheat multi-trait Ridge uses the predictive uncertainty contract", {
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  skip_if(is.null(py) || !file.exists(py), "No preferred ML Python runtime found")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")
  dat <- ml_dl_wheat_slice()

  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    GS_model = "Ridge_Regression",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    lambda_rr = 0.1,
    internal_cv_nfolds = 2L,
    internal_cv_replication = 2L,
    random_seed = 17L
  ))

  expect_wheat_predictive_contract(out)
})

test_that("BGLR wheat multi-trait MLP uses the predictive uncertainty contract", {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  skip_if(is.null(py) || !file.exists(py), "No preferred DL Python runtime found")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py)
  dat <- ml_dl_wheat_slice()

  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    GS_model = "mlp",
    response_family = "gaussian",
    multi_trait_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 17L,
    validation_split = 0.2,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    epochs = 2L,
    batch_size = 8L,
    internal_cv_nfolds = 2L,
    internal_cv_replication = 2L,
    mlp_neurons_per_layer = as.integer(c(16L, 8L)),
    mlp_learning_rate = 1e-3,
    dropout = 0
  ))

  expect_wheat_predictive_contract(out)
})
