# Response-based feature selection adds a trailing Feature_selection_remark
# column; validate_prediction_output() rejected those tables up to 0.20.223.
test_that("a trailing feature-selection remark is part of the prediction contract", {
  std <- PredictProR::prediction_output_standard(response_family = "gaussian", task = "single_trait",
                                                 gen_name = "GID")
  expect_identical(std$optional_trailing_prediction_columns, "Feature_selection_remark")
  cols <- std$prediction_columns
  pred <- as.data.frame(stats::setNames(replicate(length(cols), NA, simplify = FALSE), cols))
  pred <- pred[rep(1, 3), , drop = FALSE]
  pred$GID <- paste0("g", 1:3)
  pred$Predicted_value <- c(1, 2, 3)
  pred$Train_Test_Label <- "Train"
  base <- list(model_parameters = data.frame(stat = "m", summary = "x"), predicted_values = pred,
               diagnostic_plots = list(),
               variance_components = data.frame(Component = "genetic_variance", Components = 1,
                                                Standard_error = NA_real_))
  with_remark <- base
  with_remark$predicted_values$Feature_selection_remark <- "selected markers"
  r1 <- PredictProR::validate_prediction_output(base, strict = FALSE)
  r2 <- PredictProR::validate_prediction_output(with_remark, strict = FALSE)
  order_issue <- "prediction columns do not match the canonical order"
  expect_false(grepl(order_issue, r2$issues, fixed = TRUE))
  expect_identical(grepl(order_issue, r1$issues, fixed = TRUE), FALSE)
  moved <- with_remark
  moved$predicted_values <- moved$predicted_values[, c("Feature_selection_remark", cols)]
  r3 <- PredictProR::validate_prediction_output(moved, strict = FALSE)
  expect_true(grepl(order_issue, r3$issues, fixed = TRUE))
})
