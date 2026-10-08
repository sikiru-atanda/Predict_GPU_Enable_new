test_that("ML stability labels map from legacy reliability categories", {
  out <- PredictProR:::gp_ml_stability_labels(
    metric = c(0.9, 0.5, 0.1),
    remarks = c("Reliable", "Acceptable", "Unreliable"),
    percentage = 75
  )

  expect_identical(out$Prediction_stability, c(0.9, 0.5, 0.1))
  expect_identical(
    out$Prediction_stability_remarks,
    c("Stable", "Moderately Stable", "Unstable")
  )
  expect_identical(out$Prediction_stability_percentage, 75)
})

test_that("legacy ML bootstrap diagnostic requires explicit opt-in and uses stability labels", {
  plots <- PredictProR::diagnostic_plot_true_prediction(
    boot_results = structure(
      list(
        t = matrix(
          c(
            1.0, 1.1, 0.9,
            2.0, 2.4, 1.8,
            3.0, 3.8, 2.2
          ),
          nrow = 3,
          byrow = TRUE
        )
      ),
      class = "boot"
    ),
    GID_names = c("g1", "g2", "g3"),
    reference_variance = 0.5,
    model_for_CI_cal = "ML",
    require_heldout_calibration = FALSE,
    system_database = TRUE
  )

  expect_true(is.list(plots))
  expect_true("predicted_vs_reliability" %in% names(plots))
  p <- plots$predicted_vs_reliability
  expect_identical(p$labels$x, "Prediction Stability")
  expect_identical(p$labels$title, "Prediction Stability Analysis")
  expect_true(is.character(p$labels$subtitle))
  expect_match(
    plots$prediction_inter_vs_reliability$labels$title,
    "Prediction Intervals and Stability",
    fixed = TRUE
  )
  expect_match(
    plots$prediction_inter_vs_reliability$labels$title,
    "nominal prediction-interval coverage is not guaranteed",
    fixed = TRUE
  )
})

test_that("Gaussian ML diagnostic refuses uncalibrated bootstrap and standard-error inputs", {
  expect_error(
    PredictProR::diagnostic_plot_true_prediction(
      predictions = c(1.1, 1.5, 2.0),
      standard_errors = c(0.1, 0.2, 0.15),
      prediction_error_var = c(0.01, 0.04, 0.0225),
      boot_results = structure(
        list(t = matrix(c(1, 1.1, 0.9, 2, 2.2, 1.8, 3, 3.1, 2.9), nrow = 3, byrow = TRUE)),
        class = "boot"
      ),
      model_for_CI_cal = "ML",
      response_family = "gaussian",
      system_database = TRUE
    ),
    "require held-out or cross-fitted calibration",
    fixed = TRUE
  )
})

test_that("ML gaussian diagnostic uses the same held-out estimand as prediction output", {
  risk_calibration <- list(
    heldout_residuals = c(-0.2, 0.1, -0.1, 0.3),
    heldout_observed = c(1.0, 1.4, 2.0, 2.2),
    heldout_predictions = c(1.2, 1.3, 2.1, 1.9),
    calibration_source = "heldout_internal_cv"
  )
  plots <- PredictProR::diagnostic_plot_true_prediction(
    predictions = c(1.1, 1.5, 2.0, 2.2),
    standard_errors = c(0.1, 0.2, 0.15, 0.18),
    prediction_error_var = c(0.01, 0.04, 0.0225, 0.0324),
    lower_bound = c(0.9, 1.1, 1.7, 1.8),
    upper_bound = c(1.3, 1.9, 2.3, 2.6),
    observed_y = c(1.0, 1.4, NA, NA),
    train_test_label = c("Train", "Train", "Test", "Test"),
    reference_variance = 0.5,
    risk_calibration = risk_calibration,
    model_for_CI_cal = "ML",
    response_family = "gaussian",
    system_database = TRUE
  )

  p <- plots$predicted_vs_reliability
  expect_true(grepl("row-specific bootstrap prediction instability", p$labels$subtitle, fixed = TRUE))
  expect_true(grepl("calibration n = 4", p$labels$subtitle, fixed = TRUE))
  interval_plot <- plots$prediction_inter_vs_reliability
  expect_match(interval_plot$labels$title, "Prediction Intervals and Stability", fixed = TRUE)
  expect_equal(interval_plot$data$lower_bound, c(0.8, 1.2, 1.7, 1.9))
  expect_equal(interval_plot$data$upper_bound, c(1.4, 1.8, 2.3, 2.5))
})
