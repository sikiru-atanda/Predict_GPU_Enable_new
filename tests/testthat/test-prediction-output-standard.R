make_public_gaussian_result <- function() {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.2, 1.8),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.1, NA_real_),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    PEV_basis = "model-reported breeding-value PEV",
    Prediction_uncertainty_source = "model_reported",
    Prediction_interval_method = "model_based_gaussian",
    Prediction_interval_nominal_coverage = 0.95,
    Prediction_interval_calibration_n = NA_real_,
    lower_bound = c(1.004, 1.408),
    upper_bound = c(1.396, 2.192),
    Uncertainty = c(0.392, 0.784),
    Uncertainty_remarks = "model_based",
    Prediction_stability = NA_real_,
    Prediction_stability_remarks = NA_character_,
    Prediction_stability_reference_variance = NA_real_,
    Reliability = c(0.9, 0.8),
    Reliability_remarks = "Reliable",
    Reliability_variance_input = c(0.01, 0.04),
    Reliability_reference_variance = c(1, 1),
    Reliability_basis = "model_reported_genetic_reliability",
    stringsAsFactors = FALSE
  )
  list(
    model_parameters = data.frame(stat = "model", summary = "mock"),
    predicted_values = pred,
    diagnostic_plots = list(),
    variance_components = data.frame(
      Component = c("genetic_variance", "residual_variance"),
      Components = c(1, 0.2),
      Standard_error = c(0.1, 0.05)
    )
  )
}

test_that("prediction output standard describes stable family schemas", {
  gaussian <- prediction_output_standard("gaussian", "multi_trait_multi_environment")
  expect_identical(
    gaussian$top_level_required,
    c("model_parameters", "predicted_values", "diagnostic_plots", "variance_components")
  )
  expect_identical(gaussian$prediction_columns[1:3], c("GID", "Env", "Trait"))

  classification <- prediction_output_standard("classification", "single_trait")
  expect_true("Probability_<class>" %in% classification$prediction_columns)

  hybrid <- prediction_output_standard("gaussian", "hybrid_multi_environment")
  expect_identical(hybrid$prediction_columns[1:4], c("Hybrid_ID", "Female", "Male", "Env"))
})

test_that("prediction output validator checks direct and indexed results", {
  result <- make_public_gaussian_result()
  direct <- validate_prediction_output(result)
  expect_true(direct$valid)
  expect_identical(direct$response_family, "gaussian")

  indexed <- list(
    model_results_by_model = list(
      RandomForest = list(Yield = result),
      Ridge_Regression = list(Yield = result)
    )
  )
  report <- validate_prediction_output(indexed)
  expect_equal(nrow(report), 2L)
  expect_true(all(report$valid))
  expect_identical(report$path, c("RandomForest/Yield", "Ridge_Regression/Yield"))
})

test_that("prediction output validator fails on uncertainty drift", {
  result <- make_public_gaussian_result()
  result$predicted_values$PEV[[2L]] <- 0.4

  report <- validate_prediction_output(result, strict = FALSE)
  expect_false(report$valid)
  expect_match(report$issues, "PEV must equal")
  expect_error(validate_prediction_output(result), "PEV must equal")
})

test_that("public standardizer returns a validated canonical class without legacy aliases", {
  input <- make_public_gaussian_result()
  input$Predicted_value <- input$predicted_values
  input$Total_Predicted_value <- input$predicted_values

  result <- PredictProR:::gp_standardize_public_model_result(input, gen_name = "GID")
  expect_s3_class(result, "predictpror_model_result")
  expect_true("across_environment_predicted_values" %in% names(result))
  expect_false("Total_Predicted_value" %in% names(result))
  expect_true(validate_prediction_output(result)$valid)
})

test_that("empty covariance tables are omitted from simple model results", {
  result <- PredictProR:::gp_standardize_public_model_result(
    make_public_gaussian_result(),
    gen_name = "GID"
  )

  expect_identical(
    names(result),
    c("model_parameters", "predicted_values", "diagnostic_plots", "variance_components")
  )
})
