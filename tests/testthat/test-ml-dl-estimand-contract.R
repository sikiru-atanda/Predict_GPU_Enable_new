test_that("ML/DL predictive covariance is not labeled genetic covariance", {
  pred <- expand.grid(
    GID = paste0("g", 1:6),
    Trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  base <- rep(seq_len(6), times = 2)
  pred$Predicted_value <- ifelse(pred$Trait == "T1", base, 2 * base + 1)
  pred$Observed_value <- pred$Predicted_value +
    ifelse(pred$Trait == "T1", rep(c(-1, 1), 6), rep(c(1, -1), 6))

  out <- PredictProR:::gp_ml_predictive_matrix_bundle(
    pred = pred,
    id_col = "GID",
    dimension_col = "Trait",
    error_pred = pred,
    n_bootstrap = 50L,
    random_seed = 7L
  )

  expect_true(all(c(
    "Prediction_covariance_traits",
    "Prediction_correlation_traits",
    "Prediction_covariance_traits_SE",
    "Prediction_error_covariance_traits",
    "Prediction_error_correlation_traits"
  ) %in% names(out)))
  expect_equal(out$Prediction_correlation_traits["T1", "T2"], 1)
  expect_equal(
    out$Prediction_covariance_traits,
    stats::cov(cbind(T1 = 1:6, T2 = 2 * (1:6) + 1))
  )
  expect_true(all(is.finite(out$Prediction_covariance_traits_SE)))
  expect_false(any(grepl("Genetic", names(out), fixed = TRUE)))
})

test_that("constant multi-trait ML/DL target uncertainty is unavailable", {
  pred <- expand.grid(
    GID = paste0("g", 1:4),
    Trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred$Observed_value <- rep(c(1, 2, 3, 4), times = 2)
  pred$Predicted_value <- pred$Observed_value + 0.1
  pred$Train_Test_Label <- "Train"

  cv <- pred[rep(seq_len(nrow(pred)), each = 2), , drop = FALSE]
  cv$rep <- rep(1:2, times = nrow(pred))
  cv$Predicted_value <- cv$Observed_value + rep(c(-0.2, 0.3), times = nrow(pred))

  out <- PredictProR:::gp_ml_apply_cv_uncertainty(
    pred = pred,
    cv_pred = cv,
    id_col = "GID",
    group_col = "Trait"
  )

  expect_true(all(is.na(out$Standard_error)))
  expect_true(all(is.na(out$PEV)))
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(is.na(out$Reliability_variance_input)))
  expect_true(all(is.finite(out$Prediction_stability)))
  expect_true(all(
    out$Prediction_uncertainty_source ==
      "unavailable_no_target_specific_predictive_uncertainty"
  ))
  expect_true(all(grepl("Reliability unavailable", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("target-specific", out$PEV_basis, fixed = TRUE)))
  expect_true(all(is.finite(out$lower_bound)))
  expect_true(all(is.finite(out$upper_bound)))
  expect_true(all(out$Prediction_interval_nominal_coverage == 0.95))
})

test_that("varying multi-trait repeated predictions yield target-specific PEV", {
  pred <- expand.grid(
    GID = paste0("g", 1:4),
    Trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred$Observed_value <- rep(c(1, 2, 3, 4), times = 2)
  pred$Predicted_value <- pred$Observed_value + 0.2
  pred$Train_Test_Label <- "Train"

  cv <- pred[rep(seq_len(nrow(pred)), each = 2), , drop = FALSE]
  cv$rep <- rep(1:2, times = nrow(pred))
  amplitudes <- rep(c(0.05, 0.1, 0.2, 0.4), times = 2)
  cv$Predicted_value <- cv$Observed_value + unlist(lapply(
    amplitudes,
    function(amplitude) c(0.2 - amplitude, 0.2 + amplitude)
  ))

  out <- PredictProR:::gp_ml_apply_cv_uncertainty(
    pred = pred,
    cv_pred = cv,
    id_col = "GID",
    group_col = "Trait"
  )

  expect_true(all(is.finite(out$Standard_error)))
  expect_equal(out$PEV, out$Standard_error^2)
  expect_gt(diff(range(out$Standard_error[out$Trait == "T1"])), 0)
  expect_equal(
    as.numeric(tapply(out$PEV, out$Trait, mean)),
    rep(0.2^2, 2L),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(out$Reliability)))
  expect_equal(
    out$Reliability,
    out$Reliability_reference_variance /
      (out$Reliability_reference_variance + out$PEV)
  )
  expect_true(all(
    out$Prediction_uncertainty_source ==
      "heldout_cross_fitted_target_resampling_calibrated"
  ))
  expect_true(all(grepl("not mixed-model PEV", out$PEV_basis, fixed = TRUE)))
})

test_that("multi-trait ML/DL reliability includes labeled Train and Test predictions", {
  pred <- expand.grid(
    GID = paste0("g", 1:5),
    Trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred$Observed_value <- rep(c(1, 2, 3, NA, NA), times = 2)
  pred$Predicted_value <- c(1, 2, 3, 10, 12, 2, 4, 6, 15, 18)
  pred$Train_Test_Label <- rep(c("Train", "Train", "Train", "Test", "Test"), times = 2)

  cv <- pred[pred$Train_Test_Label == "Train", , drop = FALSE]
  cv <- cv[rep(seq_len(nrow(cv)), each = 2), , drop = FALSE]
  cv$rep <- rep(1:2, times = nrow(cv) / 2)
  cv$Predicted_value <- cv$Observed_value + rep(c(-0.1, 0.1), times = nrow(cv) / 2)

  out <- PredictProR:::gp_ml_apply_cv_uncertainty(
    pred = pred,
    cv_pred = cv,
    id_col = "GID",
    group_col = "Trait"
  )

  expected <- tapply(pred$Observed_value, pred$Trait, stats::var, na.rm = TRUE)
  expect_identical(out$Train_Test_Label, pred$Train_Test_Label)
  expect_equal(
    out$Reliability_reference_variance,
    as.numeric(expected[out$Trait])
  )
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(grepl(
    "Reliability unavailable",
    out$Reliability_basis,
    fixed = TRUE
  )))
})

test_that("all-target resampling gives unseen multi-trait Test GIDs target-specific uncertainty", {
  pred <- expand.grid(
    GID = paste0("g", 1:5),
    Trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred$Observed_value <- rep(c(1, 2, 3, NA, NA), times = 2)
  pred$Predicted_value <- c(1, 2, 3, 4, 5, 2, 4, 6, 8, 10)
  pred$Train_Test_Label <- rep(c("Train", "Train", "Train", "Test", "Test"), times = 2)

  cv <- pred[pred$Train_Test_Label == "Train", , drop = FALSE]
  cv <- cv[rep(seq_len(nrow(cv)), each = 2), , drop = FALSE]
  cv$Predicted_value <- cv$Observed_value + rep(c(-0.2, 0.3), times = nrow(cv) / 2)

  target <- pred[rep(seq_len(nrow(pred)), each = 3), , drop = FALSE]
  amplitudes <- rep(c(0.04, 0.08, 0.16, 0.32, 0.64), times = 2)
  target$Predicted_value <- unlist(Map(
    function(center, amplitude) center + c(-amplitude, 0, amplitude),
    pred$Predicted_value,
    amplitudes
  ))

  out <- PredictProR:::gp_ml_apply_cv_uncertainty(
    pred = pred,
    cv_pred = cv,
    target_pred = target,
    id_col = "GID",
    group_col = "Trait"
  )

  test_rows <- out$Train_Test_Label == "Test"
  expect_true(all(is.finite(out$Standard_error[test_rows])))
  expect_true(all(is.finite(out$PEV[test_rows])))
  expect_equal(out$PEV, out$Standard_error^2)
  expect_gt(length(unique(round(out$Standard_error[test_rows], 8))), 1L)
  expect_true(all(is.finite(out$Reliability[test_rows])))
  expect_true(all(grepl("not genetic reliability", out$Reliability_basis[test_rows], fixed = TRUE)))
})

test_that("prediction-table variance summaries make genetic limits explicit", {
  pred <- data.frame(
    GID = paste0("g", 1:5),
    Predicted_value = c(1.0, 1.4, 1.9, 2.1, 2.8),
    Observed_value = c(0.9, 1.5, 2.0, 2.0, 2.7),
    Train_Test_Label = "Train",
    PEV = c(0.04, 0.05, 0.06, 0.05, 0.07),
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::gp_prediction_table_variance_components(pred)

  unavailable <- grepl("not_identifiable$", out$Component)
  expect_true(all(is.na(out$Components[unavailable])))
  expect_true(all(c(
    "phenotypic_reference_variance",
    "prediction_variance",
    "estimated_prediction_mse"
  ) %in% out$Component))
  expect_equal(
    out$Components[out$Component == "prediction_variance"],
    stats::var(pred$Predicted_value)
  )
  expect_equal(
    out$Components[out$Component == "estimated_prediction_mse"],
    mean(pred$PEV)
  )
})
