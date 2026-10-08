test_that("multi-trait DL predicts rows with all traits missing", {
  ids <- paste0("g", 1:6)
  pheno <- data.frame(
    GID = ids,
    Trait1 = c(1, 2, 3, 4, NA, NA),
    Trait2 = c(2, 3, 4, 5, NA, NA),
    stringsAsFactors = FALSE
  )
  features <- matrix(seq_len(18), nrow = 6,
                     dimnames = list(ids, paste0("m", 1:3)))
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    gp_dl_bridge_fit_predict = function(model_type, X_train, y_train, X_test,
                                        response_family, dl_args) {
      captured$n_train <- nrow(X_train)
      captured$n_test <- nrow(X_test)
      matrix(rep(colMeans(y_train), each = nrow(X_test)),
             nrow = nrow(X_test), byrow = FALSE)
    },
    gp_multitrait_dl_gaussian_cv = function(...) {
      args <- list(...)
      captured$cv_pheno <- args$pheno_object
      captured$cv_features <- args$geno_omic_object
      NULL
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_multitrait_dl_gaussian_model(
    model_type = "mlp",
    pheno_object = pheno,
    response = c("Trait1", "Trait2"),
    geno_omic_object = features,
    gen_name = "GID",
    scaling = FALSE,
    centering = FALSE
  )

  expect_equal(captured$n_train, 4L)
  expect_equal(captured$n_test, 6L)
  expect_equal(nrow(captured$cv_pheno), 6L)
  expect_equal(nrow(captured$cv_features), 6L)
  expect_setequal(as.character(captured$cv_pheno$GID), ids)
  expect_equal(nrow(out$predicted_values), 12L)
  expect_equal(sum(out$predicted_values$Train_Test_Label == "Test"), 4L)
  expect_setequal(
    unique(out$predicted_values$GID[out$predicted_values$Train_Test_Label == "Test"]),
    c("g5", "g6")
  )
})

test_that("multi-trait DL resampling retries one transient backend failure", {
  attempts <- 0L
  out <- NULL
  expect_warning(
    out <- PredictProR:::gp_multitrait_dl_retry_fit(
      fit_call = function() {
        attempts <<- attempts + 1L
        if (attempts == 1L) stop("transient subprocess exit")
        c(1, 2, 3)
      },
      max_attempts = 2L,
      context = "test resampling fit"
    ),
    "retrying once"
  )
  expect_equal(attempts, 2L)
  expect_equal(out, c(1, 2, 3))
})

test_that("multi-trait DL can skip internal uncertainty calibration", {
  ids <- paste0("g", 1:6)
  pheno <- data.frame(
    GID = ids,
    Trait1 = c(1, 2, 3, 4, NA, NA),
    Trait2 = c(2, 3, 4, 5, NA, NA)
  )
  features <- matrix(seq_len(18), nrow = 6,
                     dimnames = list(ids, paste0("m", 1:3)))
  local_mocked_bindings(
    gp_dl_bridge_fit_predict = function(model_type, X_train, y_train, X_test,
                                        response_family, dl_args) {
      matrix(rep(colMeans(y_train), each = nrow(X_test)),
             nrow = nrow(X_test))
    },
    gp_multitrait_dl_gaussian_cv = function(...) {
      stop("internal calibration must not run")
    },
    .package = "PredictProR"
  )
  out <- PredictProR:::gp_multitrait_dl_gaussian_model(
    model_type = "mlp", pheno_object = pheno,
    response = c("Trait1", "Trait2"), geno_omic_object = features,
    gen_name = "GID", scaling = FALSE, centering = FALSE,
    dl_internal_calibration = FALSE
  )
  pred <- out$predicted_values
  expect_equal(nrow(pred), 12L)
  expect_true(all(is.finite(pred$Predicted_value)))
  expect_true(all(is.na(pred$Standard_error)))
  expect_true(all(is.na(pred$PEV)))
  expect_true(all(pred$Prediction_uncertainty_source ==
    "unavailable_internal_calibration_disabled_by_user"))
  expect_identical(out$model_parameters$summary[
    out$model_parameters$stat == "dl_internal_calibration"
  ], "FALSE")
})
