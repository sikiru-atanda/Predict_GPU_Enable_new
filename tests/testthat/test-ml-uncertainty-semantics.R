test_that("ML bootstrap target metrics use training response variance as the reference scale", {
  boot_mat <- rbind(
    c(1, 2, 10),
    c(1, 4, 12),
    c(1, 6, 14)
  )
  labels <- c("Train", "Train", "Test")
  observed_y <- c(2, 4)

  out <- PredictProR:::gp_ml_bootstrap_target_metrics(
    boot_matrix = boot_mat,
    train_test_label = labels,
    observed_y = observed_y
  )

  expect_equal(out$predicted_mean, c(1, 4, 12))
  expect_equal(out$prediction_error_var, c(0, 4, 4))
  expect_equal(out$standard_error, c(0, 2, 2))
  expect_equal(out$reference_variance, stats::var(observed_y))
})

test_that("Gaussian ML stability uses the original marker-adjustment scale", {
  out <- PredictProR:::gp_ml_gaussian_stability(
    prediction_error_var = c(0.1, 0.2),
    reference_variance = 0.5,
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5
  )

  expect_equal(out$stability, 1 - c(0.1, 0.2) / 0.5)
  expect_identical(out$remarks, c("Moderately Stable", "Moderately Stable"))
  expect_identical(out$stability_percentage, 100)
})

test_that("Gaussian ML stability falls back to zero when the reference variance is invalid", {
  out <- PredictProR:::gp_ml_gaussian_stability(
    prediction_error_var = c(0.1, 0.2),
    reference_variance = 0,
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5
  )

  expect_identical(out$stability, c(0, 0))
  expect_identical(out$remarks, c("Unstable", "Unstable"))
  expect_identical(out$stability_percentage, 0)
})

test_that("Gaussian ML uncertainty calibration rescales PEV and intervals to observed error", {
  out <- PredictProR:::gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = c(10, 20, 30),
    pred_se = c(2, 2, 2),
    pred_variances = c(4, 4, 4),
    lower_bound = c(8, 18, 28),
    upper_bound = c(12, 22, 32),
    observed_y = c(11, 19, 30),
    train_test_label = c("Train", "Train", "Test")
  )

  expected_factor <- sqrt(mean(c(1, 1)^2) / 4)
  expect_equal(out$calibration_factor, expected_factor)
  expect_equal(out$calibrated_standard_error[1:2], c(1, 1))
  expect_equal(out$calibrated_prediction_error_var[1:2], c(1, 1))
  expect_equal(out$calibrated_lower_bound[1:2], c(9, 19))
  expect_equal(out$calibrated_upper_bound[1:2], c(11, 21))
  expect_equal(out$empirical_coverage, 1)
  expect_identical(out$calibration_source, "observed_rows")
})

test_that("Gaussian ML local rank confidence drops when neighbors are too close relative to SE", {
  out <- PredictProR:::gp_ml_local_rank_confidence(
    predicted_value = c(10, 10.1, 15),
    pred_se = c(0.5, 0.5, 0.1),
    groups = c("Test", "Test", "Test")
  )

  expect_true(out$rank_confidence[1] < out$rank_confidence[3])
  expect_true(out$rank_confidence[2] < out$rank_confidence[3])
  expect_true(is.finite(out$nearest_rank_margin[1]))
})

test_that("Gaussian ML CV rank risk calibration returns usable thresholds", {
  out <- PredictProR:::gp_ml_cv_rank_risk_calibration(
    predicted_cv = c(1.0, 2.0, 3.1, 4.2, 5.0, 6.4, 7.2, 8.1, 9.4, 10.0, 11.3, 12.2),
    observed_y = c(1.1, 2.4, 2.7, 6.0, 4.1, 7.2, 5.5, 8.4, 10.3, 9.0, 12.5, 10.8)
  )

  expect_true(is.list(out))
  expect_true(is.finite(out$risk_low_cut))
  expect_true(is.finite(out$risk_high_cut))
  expect_true(out$risk_low_cut < out$risk_high_cut)
})

test_that("Gaussian ML predictor unfamiliarity flags feature-space outliers", {
  x <- rbind(
    c(0, 0),
    c(0.1, 0.05),
    c(-0.05, 0.02),
    c(0.08, -0.04),
    c(8, 8)
  )
  out <- PredictProR:::gp_ml_predictor_unfamiliarity(
    predictor_matrix = x,
    train_test_label = c("Train", "Train", "Train", "Train", "Test")
  )

  expect_true(is.list(out))
  expect_length(out$unfamiliarity, 5L)
  expect_true(out$unfamiliarity[5] > out$unfamiliarity[1])
  expect_identical(out$remarks[5], "Highly Unfamiliar")
})

test_that("Gaussian ML output separates marker stability from calibrated confidence", {
  out <- PredictProR:::gp_ml_gaussian_prediction_output(
    ids = c("g1", "g2"),
    predicted_value = c(10, 20),
    train_test_label = c("Train", "Train"),
    pred_se = c(2, 2),
    pred_variances = c(4, 4),
    lower_bound = c(8, 18),
    upper_bound = c(12, 22),
    uncertainty = c(4, 4),
    uncertainty_remarks = c("high", "high"),
    reference_variance = 10,
    observed_y = c(11, 19),
    gen_name = "GID"
  )

  expect_true(all(c(
    "Standard_error", "PEV", "Prediction_stability",
    "Prediction_confidence", "Prediction_risk_score",
    "Rank_instability_risk", "Rank_instability_risk_remarks",
    "Prediction_risk_basis", "Rank_instability_risk_basis"
  ) %in% names(out)))
  expect_equal(out$Standard_error, c(1, 1))
  expect_equal(out$PEV, c(1, 1))
  expect_equal(out$Prediction_stability, rep(0.92, 2))
  expect_equal(out$Prediction_confidence, rep(10 / 11, 2))
  expect_equal(out$Prediction_risk_score, rep(1 / 11, 2))
  expect_equal(out$Rank_instability_risk, out$Prediction_risk_score)
  expect_identical(out$Prediction_stability_remarks, c("Stable", "Stable"))
  expect_identical(out$Prediction_confidence_remarks, c("High Confidence", "High Confidence"))
  expect_identical(out$Prediction_risk_remarks, c("Low Risk", "Low Risk"))
  expect_identical(out$Rank_instability_risk_remarks, out$Prediction_risk_remarks)
  expect_identical(out$Prediction_risk_basis, c("Calibrated_Final_Prediction", "Calibrated_Final_Prediction"))
  expect_identical(out$Rank_instability_risk_basis, out$Prediction_risk_basis)

  cal <- attr(out, "uncertainty_calibration")
  expect_true(is.list(cal))
  expect_equal(cal$calibration_factor, 0.5)
})

test_that("Gaussian ML prediction output adds unfamiliarity and raises risk for novel rows", {
  x <- rbind(
    c(0, 0),
    c(0.1, 0.05),
    c(-0.05, 0.02),
    c(0.08, -0.04),
    c(8, 8)
  )
  out <- PredictProR:::gp_ml_gaussian_prediction_output(
    ids = paste0("g", 1:5),
    predicted_value = c(1, 1.1, 0.9, 1.05, 1.2),
    train_test_label = c("Train", "Train", "Train", "Train", "Test"),
    pred_se = rep(0.2, 5),
    pred_variances = rep(0.04, 5),
    lower_bound = rep(0.8, 5),
    upper_bound = rep(1.2, 5),
    uncertainty = rep(0.04, 5),
    uncertainty_remarks = rep("moderate", 5),
    reference_variance = 1,
    observed_y = c(1, 1.2, 0.8, 1.1),
    predictor_matrix = x,
    gen_name = "GID"
  )

  expect_true(all(c(
    "Prediction_unfamiliarity",
    "Prediction_unfamiliarity_remarks",
    "Prediction_unfamiliarity_basis"
  ) %in% names(out)))
  expect_identical(out$Prediction_unfamiliarity_basis[1], "predictor_space_novelty")
  expect_true(out$Prediction_unfamiliarity[5] > out$Prediction_unfamiliarity[1])
  expect_true(out$Prediction_risk_score[5] >= out$Rank_instability_risk[5])
  expect_true(any(out$Prediction_risk_score != out$Rank_instability_risk))
})

test_that("Gaussian ML prediction output marks training rows as held-out OOB risk when bootstrap OOB is available", {
  skip_if_not_installed("boot")

  set.seed(1)
  x <- cbind(v1 = rnorm(12), v2 = rnorm(12))
  y <- x[, 1] + rnorm(12, sd = 0.2)
  stat_fun <- function(data, indices) {
    fit <- stats::lm(y ~ ., data = data.frame(y = data[, 3], v1 = data[, 1], v2 = data[, 2])[indices, , drop = FALSE])
    stats::predict(fit, newdata = data.frame(v1 = data[, 1], v2 = data[, 2]))
  }
  dat <- cbind(x, y)
  boot_res <- boot::boot(dat, statistic = stat_fun, R = 8)

  pred_mean <- colMeans(boot_res$t)
  pred_var <- apply(boot_res$t, 2, stats::var)
  out <- suppressWarnings(PredictProR:::gp_ml_gaussian_prediction_output(
    ids = paste0("g", seq_along(y)),
    predicted_value = pred_mean,
    train_test_label = rep("Train", length(y)),
    pred_se = sqrt(pred_var),
    pred_variances = pred_var,
    lower_bound = apply(boot_res$t, 2, stats::quantile, probs = 0.05),
    upper_bound = apply(boot_res$t, 2, stats::quantile, probs = 0.95),
    uncertainty = pred_var,
    uncertainty_remarks = rep("high", length(y)),
    reference_variance = stats::var(y),
    observed_y = y,
    boot_results = boot_res,
    gen_name = "GID"
  ))

  expect_true(all(out$Prediction_risk_basis == "Heldout_Training_OOB"))
  expect_identical(out$Rank_instability_risk_basis, out$Prediction_risk_basis)
})

test_that("Gaussian ML summary statistics expose train/test risk basis and risk distribution", {
  predicted <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Predicted_value = c(1.1, 1.8, 2.2, 2.7),
    Train_Test_Label = c("Train", "Train", "Test", "Test"),
    Standard_error = c(0.1, 0.2, 0.15, 0.18),
    PEV = c(0.01, 0.04, 0.0225, 0.0324),
    Prediction_stability = c(0.9, 0.8, 0.85, 0.82),
    Prediction_confidence = c(0.7, 0.6, 0.5, 0.4),
    Prediction_risk_score = c(0.3, 0.4, 0.5, 0.6),
    Prediction_risk_remarks = c("Low Risk", "Moderate Risk", "Moderate Risk", "High Risk"),
    Prediction_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB", "Calibrated_Final_Prediction", "Calibrated_Final_Prediction"),
    Prediction_unfamiliarity = c(0.1, 0.3, 0.5, 0.8),
    Prediction_unfamiliarity_remarks = c("Familiar", "Moderately Unfamiliar", "Moderately Unfamiliar", "Highly Unfamiliar"),
    Prediction_unfamiliarity_basis = c("predictor_space_novelty", "predictor_space_novelty", "predictor_space_novelty", "predictor_space_novelty"),
    Rank_instability_risk = c(0.3, 0.4, 0.5, 0.6),
    Rank_instability_risk_remarks = c("Low Risk", "Moderate Risk", "Moderate Risk", "High Risk"),
    Rank_instability_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB", "Calibrated_Final_Prediction", "Calibrated_Final_Prediction"),
    stringsAsFactors = FALSE
  )
  pheno <- data.frame(GID = c("g1", "g2"), Trait = c(1, 2), stringsAsFactors = FALSE)
  geno <- matrix(c(1, 0, 0, 1), nrow = 2)

  out <- PredictProR:::summary_statistics_AI(
    predicted_object = predicted,
    pheno_object = pheno,
    response = "Trait",
    geno_omic_object = geno,
    model_parameters = data.frame(
      stat = c("model_type", "epochs", "batch_size"),
      summary = c("deepfm", "1", "16"),
      stringsAsFactors = FALSE
    ),
    response_family = "gaussian",
    GS_model = "Ridge_Regression"
  )$summary_statistics

  expect_identical(out$summary[out$stat == "Number_TrainingSet"], "2")
  expect_identical(out$summary[out$stat == "Number_TestingSet"], "2")
  expect_identical(out$summary[out$stat == "model_type"], "deepfm")
  expect_identical(out$summary[out$stat == "epochs"], "1")
  expect_identical(out$summary[out$stat == "batch_size"], "16")
  expect_identical(out$summary[out$stat == "Train_Risk_Basis"], "Heldout_Training_OOB")
  expect_identical(out$summary[out$stat == "Test_Risk_Basis"], "Calibrated_Final_Prediction")
  expect_identical(out$summary[out$stat == "Train_Low_Risk_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Train_Moderate_Risk_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Train_High_Risk_Percentage"], "0.00")
  expect_identical(out$summary[out$stat == "Test_Low_Risk_Percentage"], "0.00")
  expect_identical(out$summary[out$stat == "Test_Moderate_Risk_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Test_High_Risk_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Train_Prediction_Unfamiliarity_Basis"], "predictor_space_novelty")
  expect_identical(out$summary[out$stat == "Test_Prediction_Unfamiliarity_Basis"], "predictor_space_novelty")
  expect_identical(out$summary[out$stat == "Train_Mean_Prediction_Unfamiliarity"], "0.2000")
  expect_identical(out$summary[out$stat == "Test_Mean_Prediction_Unfamiliarity"], "0.6500")
  expect_identical(out$summary[out$stat == "Train_Familiar_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Train_Moderately_Unfamiliar_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Train_Highly_Unfamiliar_Percentage"], "0.00")
  expect_identical(out$summary[out$stat == "Test_Familiar_Percentage"], "0.00")
  expect_identical(out$summary[out$stat == "Test_Moderately_Unfamiliar_Percentage"], "50.00")
  expect_identical(out$summary[out$stat == "Test_Highly_Unfamiliar_Percentage"], "50.00")
})

test_that("DL model parameters identify requested public model and effective backend", {
  parameters <- data.frame(
    stat = c("model_type", "epochs", "batch_size"),
    summary = c("deepfm", "1", "16"),
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::gp_annotate_requested_model_parameters(
    parameters,
    model = "deepfm"
  )

  expect_identical(out$summary[out$stat == "model_type"], "deepfm")
  expect_identical(out$summary[out$stat == "requested_model"], "FactorNet")
})

test_that("ML variance report retains marginal MSE when target PEV is unavailable", {
  out <- PredictProR:::gp_ml_gaussian_variance_components(
    boot_matrix = matrix(seq_len(30), nrow = 6),
    train_test_label = rep("Train", 5),
    observed_y = seq_len(5),
    predicted_value = seq_len(5),
    reference_variance = 2.5,
    prediction_error_var = rep(NA_real_, 5),
    marginal_prediction_mse = 1.75
  )

  expect_equal(
    out$Components[out$Component == "estimated_prediction_mse"],
    1.75
  )
})
