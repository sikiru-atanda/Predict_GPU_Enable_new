test_that("hybrid model-implied CV summaries retain their model identity", {
  pred <- data.frame(
    HybridID = c("H1", "H2", "H3"),
    Observed_value = c(1.0, 2.0, 3.0),
    Predicted_value = c(1.1, 1.8, 2.9),
    cv_scenario = rep("Known_Parents_New_Hybrid", 3L),
    fold = c(1L, 1L, 2L),
    rep = 1L,
    stringsAsFactors = FALSE
  )
  eval <- data.frame(
    cv_scenario = c("Known_Parents_New_Hybrid", "Known_Parents_New_Hybrid"),
    fold = c(1L, 2L),
    rep = c(1L, 1L),
    root_mean_squared_error = c(0.15, 0.10),
    stringsAsFactors = FALSE
  )

  asreml_summary <- PredictProR:::gp_hybrid_asreml_cv_process(
    pred_df = pred,
    eval_df = eval,
    response = "Yield",
    model_type = "GBLUP"
  )$hybrid_metric_summary
  bayes_summary <- PredictProR:::gp_hybrid_bayes_cv_process(
    pred_df = pred,
    eval_df = eval,
    response = "Yield",
    model_type = "RKHS"
  )$hybrid_metric_summary
  gp_summary <- PredictProR:::gp_hybrid_gp_cv_process(
    pred_df = pred,
    eval_df = eval,
    response = "Yield",
    model_type = "GP"
  )$hybrid_metric_summary

  expect_identical(names(asreml_summary)[1:2], c("cv_scenario", "model"))
  expect_identical(unique(asreml_summary$model), "GBLUP")
  expect_identical(unique(bayes_summary$model), "RKHS")
  expect_identical(unique(gp_summary$model), "GP")
})
