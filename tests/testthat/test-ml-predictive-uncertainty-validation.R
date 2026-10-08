test_that("finite-sample residual intervals attain independent marginal coverage", {
  set.seed(812)
  calibration_residuals <- stats::rnorm(1200, sd = 1.5)
  heldout <- list(
    heldout_residuals = calibration_residuals,
    heldout_observed = calibration_residuals,
    heldout_predictions = rep(0, length(calibration_residuals)),
    calibration_source = "heldout_internal_cv"
  )

  target_prediction <- rep(0, 5000)
  target_instability <- seq(0.005, 0.05, length.out = length(target_prediction))
  calibrated <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = target_prediction,
    pred_se = sqrt(target_instability),
    pred_variances = target_instability,
    heldout_calibration = heldout,
    confidence_level = 0.95,
    require_heldout_calibration = TRUE
  )
  future_y <- stats::rnorm(length(target_prediction), sd = 1.5)
  coverage <- mean(
    future_y >= calibrated$calibrated_lower_bound &
      future_y <= calibrated$calibrated_upper_bound
  )

  expect_equal(calibrated$calibration_source, "heldout_cross_fitted_residuals")
  expect_gt(coverage, 0.935)
  expect_lt(coverage, 0.965)
  expect_gt(length(unique(round(calibrated$calibrated_standard_error, 12))), 1L)
  expect_equal(
    mean(calibrated$calibrated_prediction_error_var),
    mean(calibration_residuals^2),
    tolerance = 1e-10
  )
})

test_that("marginal held-out MSE is not repeated as target-specific SE or PEV", {
  heldout <- list(
    heldout_residuals = c(-2, -1, 0.5, 1, 2),
    heldout_observed = c(-2, -1, 0.5, 1, 2),
    heldout_predictions = rep(0, 5),
    risk_low_cut = 0.25,
    risk_high_cut = 0.75,
    calibration_source = "heldout_internal_cv"
  )
  out <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = 1:4,
    pred_se = rep(0.1, 4),
    pred_variances = rep(0.01, 4),
    heldout_calibration = heldout,
    require_heldout_calibration = TRUE
  )

  expect_true(all(is.na(out$calibrated_standard_error)))
  expect_true(all(is.na(out$calibrated_prediction_error_var)))
  expect_true(all(is.finite(out$calibrated_lower_bound)))
  expect_true(all(grepl(
    "marginal cross-fitted MSE is not repeated",
    out$prediction_error_estimand,
    fixed = TRUE
  )))
})

test_that("effectively constant target SE and PEV are reported as unavailable", {
  heldout <- list(
    heldout_residuals = c(-2, -1, 0.5, 1, 2),
    heldout_observed = c(-2, -1, 0.5, 1, 2),
    heldout_predictions = rep(0, 5),
    risk_low_cut = 0.25,
    risk_high_cut = 0.75,
    calibration_source = "heldout_internal_cv"
  )
  factor_net_like_instability <- c(
    0.000020, 0.000031, 0.000044, 0.000037, 0.000050
  )
  out <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = seq_len(5),
    pred_se = sqrt(factor_net_like_instability),
    pred_variances = factor_net_like_instability,
    heldout_calibration = heldout,
    require_heldout_calibration = TRUE
  )

  expect_true(PredictProR:::gp_ml_has_substantive_target_variation(
    factor_net_like_instability
  ))
  expect_true(all(is.na(out$calibrated_standard_error)))
  expect_true(all(is.na(out$calibrated_prediction_error_var)))
  expect_true(all(is.finite(out$calibrated_lower_bound)))
  expect_match(out$prediction_error_estimand, "1% minimum relative")
  expect_match(
    out$prediction_error_estimand,
    "marginal cross-fitted MSE is not repeated"
  )
  expect_true(all(
    PredictProR:::gp_ml_gaussian_prediction_output(
      ids = paste0("g", seq_len(5)),
      predicted_value = seq_len(5),
      train_test_label = rep("Train", 5),
      pred_se = sqrt(factor_net_like_instability),
      pred_variances = factor_net_like_instability,
      lower_bound = seq_len(5) - 0.2,
      upper_bound = seq_len(5) + 0.2,
      uncertainty = rep(0.4, 5),
      uncertainty_remarks = rep("bootstrap", 5),
      reference_variance = 1,
      observed_y = seq_len(5),
      risk_calibration = heldout,
      gen_name = "GID",
      require_heldout_calibration = TRUE
    )$Prediction_risk_basis == "Unavailable_No_Target_Specific_PEV"
  ))
})

test_that("SE variation threshold cannot be passed only by squaring into PEV", {
  n <- 48L
  catboost_like_pev <- seq(0.7804, 0.7929, length.out = n)
  catboost_like_se <- sqrt(catboost_like_pev)
  expect_true(PredictProR:::gp_ml_has_substantive_target_variation(
    catboost_like_pev
  ))
  expect_false(PredictProR:::gp_ml_has_substantive_target_variation(
    catboost_like_se
  ))

  heldout <- list(
    heldout_residuals = c(-1, -0.5, 0.5, 1),
    heldout_observed = c(-1, -0.5, 0.5, 1),
    heldout_predictions = rep(0, 4),
    heldout_fold_id = 1:4,
    target_fold_predictions = matrix(0, nrow = 4, ncol = n),
    calibration_source = "heldout_internal_cv"
  )
  local_mocked_bindings(
    gp_ml_cross_fitted_target_uncertainty = function(...) {
      list(
        prediction_mse = catboost_like_pev,
        lower = rep(-1, n),
        upper = rep(1, n),
        calibration_n = 4L
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = rep(0, n),
    pred_se = catboost_like_se,
    pred_variances = catboost_like_pev,
    observed_y = rep(0, n),
    train_test_label = rep("Train", n),
    heldout_calibration = heldout,
    require_heldout_calibration = TRUE
  )

  expect_true(all(is.na(out$calibrated_standard_error)))
  expect_true(all(is.na(out$calibrated_prediction_error_var)))
  expect_match(out$prediction_error_estimand, "1% minimum relative")
})

test_that("unavailable target SE still yields calibrated diagnostics without Rplots leak", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  heldout <- list(
    heldout_residuals = c(-2, -1, 0.5, 1, 2),
    heldout_observed = c(-2, -1, 0.5, 1, 2),
    heldout_predictions = rep(0, 5),
    calibration_source = "heldout_internal_cv"
  )
  predictions <- seq_len(5)
  instability <- c(0.000020, 0.000031, 0.000044, 0.000037, 0.000050)
  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)
  unlink("Rplots.pdf")

  plot <- PredictProR:::diagnostic_plot_true_prediction(
    predictions = predictions,
    standard_errors = sqrt(instability),
    prediction_error_var = instability,
    lower_bound = predictions - 0.2,
    upper_bound = predictions + 0.2,
    observed_y = predictions,
    train_test_label = rep("Train", length(predictions)),
    reference_variance = stats::var(predictions),
    risk_calibration = heldout,
    require_heldout_calibration = TRUE,
    model_for_CI_cal = "ML",
    system_database = FALSE
  )

  expect_s3_class(plot, "gtable")
  expect_false(file.exists("Rplots.pdf"))

  plot_list <- PredictProR:::diagnostic_plot_true_prediction(
    predictions = predictions,
    standard_errors = sqrt(instability),
    prediction_error_var = instability,
    lower_bound = predictions - 0.2,
    upper_bound = predictions + 0.2,
    observed_y = predictions,
    train_test_label = rep("Train", length(predictions)),
    reference_variance = stats::var(predictions),
    risk_calibration = heldout,
    require_heldout_calibration = TRUE,
    model_for_CI_cal = "ML",
    system_database = TRUE
  )
  expect_identical(
    plot_list$prediction_inter_vs_reliability$labels$title,
    "Prediction Intervals and Stability"
  )
  expect_match(
    plot_list$prediction_inter_vs_reliability$labels$subtitle,
    "Target-specific SE/PEV unavailable"
  )
  expect_lt(
    nchar(plot_list$prediction_inter_vs_reliability$labels$subtitle),
    120L
  )
})

test_that("Gaussian ML/DL output retains observed training responses", {
  predicted <- c(1.1, 1.8, 2.9, 4.2)
  heldout <- list(
    heldout_residuals = c(-0.4, 0.2, 0.3, -0.1, 0.5, -0.2),
    heldout_observed = 1:6,
    heldout_predictions = 1:6 - c(-0.4, 0.2, 0.3, -0.1, 0.5, -0.2),
    risk_low_cut = 0.25,
    risk_high_cut = 0.75,
    calibration_source = "heldout_internal_cv"
  )
  out <- PredictProR:::gp_ml_gaussian_prediction_output(
    ids = paste0("g", 1:4),
    predicted_value = predicted,
    train_test_label = c("Train", "Train", "Train", "Test"),
    pred_se = sqrt(c(0.01, 0.03, 0.06, 0.09)),
    pred_variances = c(0.01, 0.03, 0.06, 0.09),
    lower_bound = predicted - 0.2,
    upper_bound = predicted + 0.2,
    uncertainty = rep(0.4, 4),
    uncertainty_remarks = rep("bootstrap", 4),
    reference_variance = 1,
    observed_y = c(1, 2, 3),
    risk_calibration = heldout,
    gen_name = "GID",
    require_heldout_calibration = TRUE
  )

  expect_equal(out$Observed_value, c(1, 2, 3, NA_real_))
  expect_gt(length(unique(round(out$Standard_error, 12))), 1L)
  expect_equal(out$PEV, out$Standard_error^2)
  expect_equal(mean(out$PEV[1:3]), mean(heldout$heldout_residuals^2))
})

test_that("strict ML/DL uncertainty is unavailable without held-out calibration", {
  out <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = c(1, 2, 3),
    pred_se = rep(0.1, 3),
    pred_variances = rep(0.01, 3),
    observed_y = c(1, 2, 3),
    train_test_label = rep("Train", 3),
    require_heldout_calibration = TRUE
  )

  expect_identical(out$calibration_source, "unavailable_no_heldout_calibration")
  expect_true(all(is.na(out$calibrated_prediction_error_var)))
  expect_true(all(is.na(out$calibrated_lower_bound)))
})

test_that("ML/DL Reliability requires target-specific PEV and stays separate from bootstrap stability", {
  final_prediction <- seq_len(6)
  fold_shift <- rbind(
    final_prediction - c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
    final_prediction,
    final_prediction + c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
  )
  observed <- seq_len(12)
  residual <- rep(c(-0.3, 0.2, 0.1, -0.2), 3)
  fold_id <- rep(seq_len(3), each = 4)
  heldout <- PredictProR:::gp_ml_cv_rank_risk_calibration(
    predicted_cv = observed - residual,
    observed_y = observed,
    fold_id = fold_id,
    target_fold_predictions = fold_shift
  )
  out <- PredictProR:::gp_ml_gaussian_prediction_output(
    ids = paste0("g", 1:6),
    predicted_value = final_prediction,
    train_test_label = rep("Train", 6),
    pred_se = sqrt(seq(0.01, 0.06, length.out = 6)),
    pred_variances = seq(0.01, 0.06, length.out = 6),
    lower_bound = final_prediction - 0.2,
    upper_bound = final_prediction + 0.2,
    uncertainty = rep(0.4, 6),
    uncertainty_remarks = rep("bootstrap", 6),
    reference_variance = stats::var(final_prediction),
    observed_y = final_prediction,
    risk_calibration = heldout,
    gen_name = "GID",
    require_heldout_calibration = TRUE
  )

  expect_true(all(is.finite(out$Prediction_stability)))
  expect_true(all(is.finite(out$PEV)))
  expect_gt(length(unique(round(out$PEV, 8))), 1L)
  expect_equal(out$Reliability_variance_input, out$PEV)
  expect_equal(
    unique(out$Reliability_reference_variance),
    stats::var(final_prediction)
  )
  expect_equal(
    out$Reliability,
    out$Reliability_reference_variance /
      (out$Reliability_reference_variance + out$PEV)
  )
  expect_false(isTRUE(all.equal(out$Reliability, out$Prediction_stability)))
  expect_true(all(grepl("not genetic reliability", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("V_y/(V_y + PEV_i)", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("cross-fitted", out$PEV_basis, fixed = TRUE)))
  expect_true(all(out$Prediction_interval_nominal_coverage == 0.95))

  formatted <- PredictProR:::gp_format_gaussian_prediction_table(out, gen_name = "GID")
  expect_equal(formatted$Reliability, out$Reliability)
  expect_equal(formatted$Reliability_variance_input, out$Reliability_variance_input)
  expect_equal(formatted$Prediction_stability, out$Prediction_stability)
  expect_equal(formatted$PEV_basis, out$PEV_basis)

  public <- list(
    model_parameters = data.frame(stat = "model", summary = "mock"),
    predicted_values = formatted,
    diagnostic_plots = list(),
    variance_components = data.frame(
      Component = "genetic_variance_not_identifiable",
      Components = NA_real_,
      Standard_error = NA_real_
    )
  )
  expect_true(validate_prediction_output(public)$valid)
  public$predicted_values$Reliability[[1L]] <-
    public$predicted_values$Reliability[[1L]] + 0.05
  expect_error(
    validate_prediction_output(public),
    "V_y/\\(V_y \\+ PEV_i\\)"
  )
})

test_that("local cross-fitted Xgboost uncertainty separates conditional residual risk", {
  heldout_ids <- paste0("g", 1:12)
  target_ids <- c(heldout_ids, "test_low", "test_high")
  heldout_x <- cbind(
    marker_1 = c(seq(0, 0.5, length.out = 6), seq(5, 5.5, length.out = 6)),
    marker_2 = c(seq(0.2, 0.7, length.out = 6), seq(4.8, 5.3, length.out = 6))
  )
  rownames(heldout_x) <- heldout_ids
  target_x <- rbind(heldout_x, test_low = c(0.1, 0.3), test_high = c(5.1, 4.9))
  predicted_cv <- seq_len(12)
  residuals <- c(0.10, -0.15, 0.12, -0.20, 0.18, -0.11,
                 0.90, -1.10, 1.25, -1.35, 1.05, -1.20)
  observed <- predicted_cv + residuals
  final_prediction <- c(seq_len(12), 1.5, 10.5)
  instability_scale <- seq(0.02, 0.16, length.out = length(target_ids))
  target_fold_predictions <- rbind(
    final_prediction - instability_scale,
    final_prediction,
    final_prediction + instability_scale
  )

  calibration <- PredictProR:::gp_ml_cv_rank_risk_calibration(
    predicted_cv = predicted_cv,
    observed_y = observed,
    fold_id = rep(1:3, each = 4),
    target_fold_predictions = target_fold_predictions,
    heldout_predictors = heldout_x,
    target_predictors = target_x,
    heldout_ids = heldout_ids,
    target_ids = target_ids
  )
  out <- PredictProR:::gp_ml_cross_fitted_target_uncertainty(
    calibration = calibration,
    final_prediction = final_prediction
  )

  expect_identical(
    out$source,
    "local_cross_fitted_residual_risk_plus_fold_instability"
  )
  expect_equal(
    mean(out$prediction_mse[seq_along(heldout_ids)]),
    mean(residuals^2),
    tolerance = 1e-12
  )
  expect_gt(
    mean(out$prediction_mse[7:12]),
    mean(out$prediction_mse[1:6])
  )
  expect_gt(out$prediction_mse[14], out$prediction_mse[13])
  expect_gt(length(unique(round(sqrt(out$prediction_mse), 8))), 2L)
  expect_true(all(is.finite(out$lower)))
  expect_true(all(is.finite(out$upper)))
  expect_true(is.finite(out$empirical_coverage))
  expect_true(out$empirical_coverage >= 0 && out$empirical_coverage <= 1)
  expect_true(all(grepl("locally weighted", out$error_estimand, fixed = TRUE)))
})

test_that("every registered DL backend receives local target uncertainty", {
  expected <- sort(c(
    "cnn", "mlp", "mlp_with_attention", "ft_transformer", "resnet",
    "saint", "tabnet", "node", "deepfm", "dcnv2", "nam", "moe",
    "gp_dkl"
  ))
  registered <- PredictProR:::gp_dl_registered_model_types()

  expect_identical(registered, expected)
  expect_true(all(vapply(
    registered,
    PredictProR:::gp_dl_supports_local_target_uncertainty,
    logical(1)
  )))
  expect_false(PredictProR:::gp_dl_supports_local_target_uncertainty("not_a_dl_model"))
})

test_that("all three boosted-tree models receive local target uncertainty", {
  expected <- sort(c("xgboost", "catboost", "lightgbm"))

  expect_identical(
    PredictProR:::gp_py_boosted_tree_local_uncertainty_models(),
    expected
  )
  expect_true(all(vapply(
    expected,
    PredictProR:::gp_py_uses_local_target_uncertainty,
    logical(1)
  )))
  expect_false(PredictProR:::gp_py_uses_local_target_uncertainty("randomforest"))
})

test_that("marker-adjustment reliability uses the original clipped formula without a 0.7 floor", {
  out <- PredictProR:::gp_ml_gaussian_stability(
    prediction_error_var = c(0.05, 0.3, 1.2),
    reference_variance = 1
  )

  expect_equal(out$stability, c(0.95, 0.7, 0))
  expect_lt(out$stability[[3L]], 0.7)
})

test_that("ML/DL reliability uses training phenotype variance and calibrated PEV", {
  predicted <- c(1, 2, 3, 10, 12)
  labels <- c("Train", "Train", "Train", "Test", "Test")
  out <- PredictProR:::gp_ml_gaussian_prediction_output(
    ids = paste0("g", seq_along(predicted)),
    predicted_value = predicted,
    train_test_label = labels,
    pred_se = rep(0.1, length(predicted)),
    pred_variances = rep(0.01, length(predicted)),
    lower_bound = predicted - 0.2,
    upper_bound = predicted + 0.2,
    uncertainty = rep(0.4, length(predicted)),
    uncertainty_remarks = rep("bootstrap", length(predicted)),
    reference_variance = stats::var(predicted[labels == "Train"]),
    observed_y = predicted[labels == "Train"],
    gen_name = "GID"
  )

  expect_identical(out$Train_Test_Label, labels)
  expect_equal(
    unique(out$Reliability_reference_variance),
    stats::var(predicted[labels == "Train"])
  )
  expect_equal(
    out$Reliability,
    out$Reliability_reference_variance /
      (out$Reliability_reference_variance + out$PEV)
  )
  expect_true(all(grepl(
    "training phenotype reference variance",
    out$Reliability_basis,
    fixed = TRUE
  )))
})

test_that("external uncertainty validation is held-out and reports calibration", {
  set.seed(91)
  n <- 800
  observed <- stats::rnorm(n)
  pred <- data.frame(
    Observed_value = observed,
    Predicted_value = rep(0, n),
    Train_Test_Label = "Test",
    PEV = 1,
    lower_bound = -1.96,
    upper_bound = 1.96,
    Prediction_stability = stats::runif(n),
    stringsAsFactors = FALSE
  )
  audit <- validate_predictive_uncertainty(pred)

  expect_s3_class(audit, "PredictProR_uncertainty_validation")
  expect_equal(audit$overall$pev_to_mse_ratio, 1, tolerance = 0.15)
  expect_gt(audit$overall$empirical_coverage, 0.925)
  expect_lt(audit$overall$empirical_coverage, 0.975)
  expect_true(
    audit$overall$coverage_ci_lower <= 0.95 &&
      audit$overall$coverage_ci_upper >= 0.95
  )
})

test_that("external uncertainty validation refuses in-sample claims by default", {
  pred <- data.frame(
    Observed_value = 1:5,
    Predicted_value = 1:5,
    Train_Test_Label = "Train"
  )
  expect_error(
    validate_predictive_uncertainty(pred),
    "No held-out evaluation rows"
  )
})

test_that("MET ML/DL does not promote in-sample residuals to uncertainty", {
  pred <- data.frame(
    GID = paste0("g", 1:8),
    Env = rep(c("E1", "E2"), each = 4),
    Predicted_value = seq_len(8),
    Observed_value = seq_len(8) + 0.1,
    Train_Test_Label = "Train",
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::gp_met_add_gaussian_uncertainty(
    pred_obs = pred,
    gen_name = "GID",
    heter_groups = "Env",
    ml_dl_estimand = TRUE
  )

  expect_true(all(is.na(out$PEV)))
  expect_true(all(is.na(out$lower_bound)))
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(is.na(out$Reliability_variance_input)))
  expect_true(all(is.finite(out$Prediction_stability)))
  expect_true(all(grepl("Reliability unavailable", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("not genetic reliability", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(
    out$Prediction_uncertainty_source ==
      "unavailable_no_heldout_or_model_based_uncertainty"
  ))
})
