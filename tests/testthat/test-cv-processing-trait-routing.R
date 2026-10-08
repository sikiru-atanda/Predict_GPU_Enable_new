test_that("single-location CV processing ranks available models per trait without requiring a full trait grid", {
  skip_if_not_installed("plotly")

  cv_results_data <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L,
        accuracy = 0.80,
        mean_squared_error = 0.20
      )
    ),
    list(
      trait = "Yield",
      model = "ModelB",
      rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L,
        accuracy = 0.90,
        mean_squared_error = 0.10
      )
    ),
    list(
      trait = "Height",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L,
        accuracy = 0.70,
        mean_squared_error = 0.30
      )
    )
  )

  res <- PredictProR:::cv_single_loc_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = c("accuracy", "mean_squared_error")
  )

  best_acc <- res$best_models_list$accuracy
  expect_identical(sort(unique(best_acc$trait)), c("Height", "Yield"))
  expect_identical(best_acc$model[best_acc$trait == "Yield"], "ModelB")
  expect_identical(best_acc$model[best_acc$trait == "Height"], "ModelA")
})

test_that("single-location CV processing can skip interactive plot rendering", {
  cv_results_data <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L,
        accuracy = 0.80,
        mean_squared_error = 0.20
      )
    ),
    list(
      trait = "Yield",
      model = "ModelB",
      rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L,
        accuracy = 0.90,
        mean_squared_error = 0.10
      )
    )
  )

  res <- PredictProR:::cv_single_loc_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = "accuracy",
    render_interactive = FALSE
  )

  expect_null(res$plot_reps_list$accuracy$plotly_boxplot_reps)
})

test_that("true prediction task table uses one CV-selected model per trait", {
  cv_results_processed <- list(
    best_models_list = list(
      accuracy = data.frame(
        trait = c("HDDT", "HT"),
        model = c("BRR", "RKHS"),
        accuracy = c(0.82, 0.79),
        stringsAsFactors = FALSE
      )
    ),
    aggregated_data_list = list(
      aggregated_across_reps = data.frame(
        trait = c("HDDT", "HDDT", "HT", "HT"),
        model = c("BRR", "RKHS", "BRR", "RKHS"),
        accuracy = c(0.82, 0.70, 0.60, 0.79),
        stringsAsFactors = FALSE
      )
    )
  )

  tasks <- PredictProR:::gp_build_true_prediction_task_table(
    best_models = NULL,
    GS_model = NULL,
    GS_model_cv = c("BRR", "RKHS"),
    response = c("HDDT", "HT"),
    cv_results_processed = cv_results_processed,
    metric_for_ranking = "accuracy",
    cross_validation = TRUE
  )

  expect_identical(tasks$trait, c("HDDT", "HT"))
  expect_identical(tasks$model, c("BRR", "RKHS"))
  expect_true(all(tasks$is_cv_best))
})

test_that("true prediction task table reduces full CV rankings to one best model per trait", {
  full_cv_rankings <- data.frame(
    trait = c("HDDT", "HDDT", "HT", "HT", "YIELD", "YIELD"),
    model = c("BRR", "BayesA", "BRR", "BayesA", "BayesB", "BayesC"),
    accuracy = c(0.41, 0.88, 0.76, 0.52, 0.61, 0.93),
    stringsAsFactors = FALSE
  )

  tasks <- PredictProR:::gp_build_true_prediction_task_table(
    best_models = full_cv_rankings,
    GS_model = NULL,
    GS_model_cv = c("BRR", "BayesA", "BayesB", "BayesC"),
    response = c("HDDT", "HT", "YIELD"),
    cv_results_processed = NULL,
    metric_for_ranking = "accuracy",
    cross_validation = TRUE
  )

  expect_identical(tasks$trait, c("HDDT", "HT", "YIELD"))
  expect_identical(tasks$model, c("BayesA", "BRR", "BayesC"))
  expect_equal(nrow(tasks), 3)
  expect_false(any(tasks$trait == "HDDT" & tasks$model == "BRR"))
  expect_false(any(tasks$trait == "HT" & tasks$model == "BayesA"))
  expect_true(all(tasks$is_cv_best))
})

test_that("true prediction task table pairs explicit no-CV models with traits by position", {
  tasks <- PredictProR:::gp_build_true_prediction_task_table(
    best_models = NULL,
    GS_model = c("BRR", "RKHS"),
    GS_model_cv = NULL,
    response = c("HDDT", "HT"),
    cv_results_processed = NULL,
    metric_for_ranking = "accuracy",
    cross_validation = FALSE
  )

  expect_identical(tasks$trait, c("HDDT", "HT"))
  expect_identical(tasks$model, c("BRR", "RKHS"))
})

test_that("predicted-vs-observed CV diagnostics keep trait-local models when another trait/model combo is missing", {
  local_mocked_bindings(
    predicted_vs_observed_ranking_plot = function(...) list(plot = "ok"),
    .package = "PredictProR"
  )

  results <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        y = c(1, 2),
        yhat = c(1.1, 1.9)
      )
    ),
    list(
      trait = "Yield",
      model = "ModelB",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        y = c(1, 2),
        yhat = c(0.9, 2.1)
      )
    ),
    list(
      trait = "Height",
      model = "ModelA",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        y = c(5, 6),
        yhat = c(5.2, 5.8)
      )
    )
  )

  res <- PredictProR:::single_predicted_vs_observed_result_plots_process(
    results = results,
    pheno_data = data.frame(
      GID = c("g1", "g2"),
      Yield = c(1, 2),
      Height = c(5, 6),
      stringsAsFactors = FALSE
    ),
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05
  )

  expect_identical(
    sort(names(res$predicted_vs_observed_plots$Yield$predicted_vs_observed_plots)),
    c("ModelA", "ModelB")
  )
  expect_identical(
    names(res$predicted_vs_observed_plots$Height$predicted_vs_observed_plots),
    "ModelA"
  )
  expect_true("Yield" %in% names(res$mod_res_per_trait_per_model$ModelB))
  expect_false("Height" %in% names(res$mod_res_per_trait_per_model$ModelB))
})

test_that("predicted-vs-observed CV diagnostics average repeated predictions by row identity", {
  local_mocked_bindings(
    predicted_vs_observed_ranking_plot = function(...) list(plot = "ok"),
    .package = "PredictProR"
  )

  results <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L, 3L),
        y = c(10, 20, 30),
        yhat = c(11, NA, 31),
        cv_role = c("test", "train", "test"),
        stringsAsFactors = FALSE
      )
    ),
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 2L,
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L, 3L),
        y = c(10, 20, 30),
        yhat = c(13, 22, NA),
        cv_role = c("test", "test", "train"),
        stringsAsFactors = FALSE
      )
    )
  )

  res <- PredictProR:::single_predicted_vs_observed_result_plots_process(
    results = results,
    pheno_data = data.frame(
      GID = c("g1", "g2", "g3"),
      Yield = c(10, 20, 30),
      stringsAsFactors = FALSE
    ),
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05
  )

  stats <- res$mod_res_per_trait_per_model$ModelA$Yield
  expect_equal(stats$y_mean, c(10, 20, 30))
  expect_equal(stats$pred_mean, c(12, 22, 31))
})

test_that("CV processing returns aggregated classification probability summaries", {
  cv_results_data <- list(
    list(
      trait = "Class",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(Rep = 1L, accuracy = 0.8, log_loss = 0.4),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c("A", "B"),
        yhat = c("A", "B"),
        cv_role = c("test", "test"),
        Env = c("E1", "E2"),
        stringsAsFactors = FALSE
      ),
      yprob_cv_Reps_all = data.frame(
        Prob_A = c(0.8, 0.2),
        Prob_B = c(0.1, 0.7),
        Prob_C = c(0.1, 0.1)
      )
    ),
    list(
      trait = "Class",
      model = "ModelA",
      rep = 2L,
      eval_metrics_reps = data.frame(Rep = 2L, accuracy = 0.9, log_loss = 0.3),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c("A", "B"),
        yhat = c("A", "B"),
        cv_role = c("test", "test"),
        Env = c("E1", "E2"),
        stringsAsFactors = FALSE
      ),
      yprob_cv_Reps_all = data.frame(
        Prob_A = c(0.7, 0.3),
        Prob_B = c(0.2, 0.6),
        Prob_C = c(0.1, 0.1)
      )
    )
  )

  res <- PredictProR:::cv_single_loc_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = c("accuracy", "log_loss")
  )

  probs <- res$classification_probability_summaries$aggregated_probabilities
  prob_cols <- res$classification_probability_summaries$probability_columns
  expect_true(is.data.frame(probs))
  expect_true(all(c("Prob_A", "Prob_B", "Prob_C") %in% names(probs)))
  expect_true("Env" %in% names(probs))
  expect_equal(prob_cols, c("Prob_A", "Prob_B", "Prob_C"))
  expect_true(any(is.finite(as.matrix(probs[, prob_cols, drop = FALSE]))))
  expect_equal(probs$Prob_A, c(0.75, 0.25))
  expect_equal(probs$Prob_B, c(0.15, 0.65))
})

test_that("CV processing returns aggregated gaussian uncertainty summaries", {
  cv_results_data <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(Rep = 1L, root_mean_squared_error = 0.3),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c(10, 20),
        yhat = c(11, 19),
        cv_role = c("test", "test"),
        Standard_error = c(0.5, 0.7),
        PEV = c(0.25, 0.49),
        lower_bound = c(9, 18),
        upper_bound = c(13, 21),
        stringsAsFactors = FALSE
      )
    ),
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 2L,
      eval_metrics_reps = data.frame(Rep = 2L, root_mean_squared_error = 0.4),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c(10, 20),
        yhat = c(12, 20),
        cv_role = c("test", "test"),
        Standard_error = c(0.6, 0.8),
        PEV = c(0.36, 0.64),
        lower_bound = c(9.5, 18.5),
        upper_bound = c(13.5, 21.5),
        stringsAsFactors = FALSE
      )
    )
  )

  res <- PredictProR:::cv_single_loc_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = c("root_mean_squared_error")
  )

  unc <- res$gaussian_uncertainty_summaries$aggregated_uncertainty
  expect_true(is.data.frame(unc))
  expect_true(all(c(
    "mean_absolute_error",
    "root_mean_squared_error",
    "mean_standard_error",
    "mean_pev",
    "mean_interval_width",
    "empirical_coverage"
  ) %in% names(unc)))
  expect_equal(unc$mean_absolute_error, 1)
  expect_equal(unc$root_mean_squared_error, sqrt(1.5), tolerance = 1e-8)
  expect_equal(unc$mean_standard_error, 0.65, tolerance = 1e-8)
  expect_equal(unc$mean_pev, 0.435, tolerance = 1e-8)
  expect_equal(unc$mean_interval_width, 3.5, tolerance = 1e-8)
  expect_equal(unc$empirical_coverage, 1)
})

test_that("CV processing returns aggregated gaussian risk summaries", {
  cv_results_data <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(Rep = 1L, root_mean_squared_error = 0.3),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c(10, 20),
        yhat = c(11, 19),
        cv_role = c("test", "test"),
        Prediction_confidence = c(0.7, 0.5),
        Prediction_risk_score = c(0.3, 0.5),
        Prediction_risk_remarks = c("Low Risk", "Moderate Risk"),
        Prediction_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        Rank_instability_risk = c(0.3, 0.5),
        Rank_instability_risk_remarks = c("Low Risk", "Moderate Risk"),
        Rank_instability_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        stringsAsFactors = FALSE
      )
    ),
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 2L,
      eval_metrics_reps = data.frame(Rep = 2L, root_mean_squared_error = 0.4),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c(10, 20),
        yhat = c(12, 20),
        cv_role = c("test", "test"),
        Prediction_confidence = c(0.6, 0.4),
        Prediction_unfamiliarity = c(0.1, 0.3),
        Prediction_unfamiliarity_remarks = c("Familiar", "Moderately Unfamiliar"),
        Prediction_unfamiliarity_basis = c("predictor_space_novelty", "predictor_space_novelty"),
        Prediction_risk_score = c(0.4, 0.6),
        Prediction_risk_remarks = c("Moderate Risk", "High Risk"),
        Prediction_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        Rank_instability_risk = c(0.4, 0.6),
        Rank_instability_risk_remarks = c("Moderate Risk", "High Risk"),
        Rank_instability_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        stringsAsFactors = FALSE
      )
    )
  )

  res <- PredictProR:::cv_single_loc_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = c("root_mean_squared_error")
  )

  risk <- res$gaussian_risk_summaries$aggregated_risk
  expect_true(is.data.frame(risk))
  expect_true(all(c(
    "mean_prediction_confidence",
    "mean_prediction_unfamiliarity",
    "mean_prediction_risk_score",
    "mean_rank_instability_risk",
    "low_risk_percentage",
    "moderate_risk_percentage",
    "high_risk_percentage",
    "low_unfamiliarity_percentage",
    "moderate_unfamiliarity_percentage",
    "high_unfamiliarity_percentage",
    "prediction_unfamiliarity_basis",
    "prediction_risk_basis",
    "rank_instability_risk_basis"
  ) %in% names(risk)))
  expect_equal(risk$mean_prediction_confidence, 0.55, tolerance = 1e-8)
  expect_equal(risk$mean_prediction_unfamiliarity, 0.20, tolerance = 1e-8)
  expect_equal(risk$mean_prediction_risk_score, 0.45, tolerance = 1e-8)
  expect_equal(risk$mean_rank_instability_risk, 0.45, tolerance = 1e-8)
  expect_equal(risk$low_risk_percentage, 25, tolerance = 1e-8)
  expect_equal(risk$moderate_risk_percentage, 50, tolerance = 1e-8)
  expect_equal(risk$high_risk_percentage, 25, tolerance = 1e-8)
  expect_equal(risk$low_unfamiliarity_percentage, 50, tolerance = 1e-8)
  expect_equal(risk$moderate_unfamiliarity_percentage, 50, tolerance = 1e-8)
  expect_equal(risk$high_unfamiliarity_percentage, 0, tolerance = 1e-8)
  expect_identical(risk$prediction_unfamiliarity_basis, "predictor_space_novelty")
  expect_identical(risk$prediction_risk_basis, "Heldout_Training_OOB")
  expect_identical(risk$rank_instability_risk_basis, "Heldout_Training_OOB")

  expect_true(is.list(res$gaussian_risk_plots))
  expect_true(all(c("mean_prediction_unfamiliarity", "mean_rank_instability_risk", "risk_distribution") %in% names(res$gaussian_risk_plots)))
  expect_true(inherits(res$gaussian_risk_plots$mean_prediction_unfamiliarity, "ggplot"))
  expect_true(inherits(res$gaussian_risk_plots$mean_rank_instability_risk, "ggplot"))
  expect_true(inherits(res$gaussian_risk_plots$risk_distribution, "ggplot"))
})

test_that("multi-environment CV processing returns gaussian risk plots with environment splits", {
  cv_results_data <- list(
    list(
      trait = "Yield",
      model = "ModelA",
      rep = 1L,
      eval_metrics_reps = data.frame(Rep = 1L, root_mean_squared_error = 0.3, Env = "E1"),
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L),
        y = c(10, 20),
        yhat = c(11, 19),
        cv_role = c("test", "test"),
        Env = c("E1", "E1"),
        Prediction_confidence = c(0.7, 0.5),
        Prediction_unfamiliarity = c(0.2, 0.4),
        Prediction_unfamiliarity_remarks = c("Familiar", "Moderately Unfamiliar"),
        Prediction_unfamiliarity_basis = c("predictor_space_novelty", "predictor_space_novelty"),
        Prediction_risk_score = c(0.3, 0.5),
        Prediction_risk_remarks = c("Low Risk", "Moderate Risk"),
        Prediction_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        Rank_instability_risk = c(0.3, 0.5),
        Rank_instability_risk_remarks = c("Low Risk", "Moderate Risk"),
        Rank_instability_risk_basis = c("Heldout_Training_OOB", "Heldout_Training_OOB"),
        stringsAsFactors = FALSE
      )
    ),
    list(
      trait = "Yield",
      model = "ModelB",
      rep = 1L,
      eval_metrics_reps = data.frame(Rep = 1L, root_mean_squared_error = 0.4, Env = "E2"),
      ypred_cv_Reps_all = data.frame(
        row_id = c(3L, 4L),
        y = c(12, 22),
        yhat = c(10, 24),
        cv_role = c("test", "test"),
        Env = c("E2", "E2"),
        Prediction_confidence = c(0.4, 0.3),
        Prediction_unfamiliarity = c(0.5, 0.8),
        Prediction_unfamiliarity_remarks = c("Highly Unfamiliar", "Highly Unfamiliar"),
        Prediction_unfamiliarity_basis = c("predictor_space_novelty", "predictor_space_novelty"),
        Prediction_risk_score = c(0.6, 0.7),
        Prediction_risk_remarks = c("High Risk", "High Risk"),
        Prediction_risk_basis = c("Calibrated_Final_Prediction", "Calibrated_Final_Prediction"),
        Rank_instability_risk = c(0.6, 0.7),
        Rank_instability_risk_remarks = c("High Risk", "High Risk"),
        Rank_instability_risk_basis = c("Calibrated_Final_Prediction", "Calibrated_Final_Prediction"),
        stringsAsFactors = FALSE
      )
    )
  )

  res <- PredictProR:::cv1_cv2_and_across_env_result_plot_process(
    cv_results_data = cv_results_data,
    eval_metrics = c("root_mean_squared_error"),
    heter_groups = "Env"
  )

  expect_true(is.list(res$gaussian_risk_plots))
  expect_true(inherits(res$gaussian_risk_plots$mean_prediction_unfamiliarity, "ggplot"))
  expect_true(inherits(res$gaussian_risk_plots$mean_rank_instability_risk, "ggplot"))
  expect_true(inherits(res$gaussian_risk_plots$risk_distribution, "ggplot"))
})

test_that("classification CV diagnostics build plots from stored probabilities", {
  results <- list(
    list(
      trait = "Class",
      model = "ModelA",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L, 3L),
        y = c("A", "B", "A"),
        yhat = c("A", "B", "B"),
        cv_role = c("test", "test", "test"),
        stringsAsFactors = FALSE
      ),
      yprob_cv_Reps_all = data.frame(
        Prob_A = c(0.8, 0.2, 0.4),
        Prob_B = c(0.2, 0.8, 0.6)
      )
    )
  )

  res <- PredictProR:::single_predicted_vs_observed_result_plots_process(
    results = results,
    pheno_data = data.frame(
      GID = c("g1", "g2", "g3"),
      Class = c("A", "B", "A"),
      stringsAsFactors = FALSE
    ),
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05
  )

  expect_true("Class" %in% names(res$predicted_vs_observed_plots))
  expect_true("ModelA" %in% names(res$predicted_vs_observed_plots$Class$predicted_vs_observed_plots))
  expect_true(inherits(res$predicted_vs_observed_plots$Class$predicted_vs_observed_plots$ModelA, "gtable"))
  expect_true(all(c("Prob_A", "Prob_B") %in% names(res$mod_res_per_trait_per_model$ModelA$Class)))
  expect_true("ModelA__Class" %in% names(res$classification_diagnostic_summaries))
  expect_true(is.data.frame(res$classification_diagnostic_summaries$ModelA__Class$confusion_matrix))
  expect_true(is.data.frame(res$classification_diagnostic_summaries$ModelA__Class$calibration_summary))
})

test_that("multi-environment classification diagnostics build per-environment plots", {
  results <- list(
    list(
      trait = "Class",
      model = "ModelA",
      rep = 1L,
      ypred_cv_Reps_all = data.frame(
        row_id = c(1L, 2L, 3L, 4L),
        y = c("A", "B", "A", "B"),
        yhat = c("A", "B", "B", "B"),
        cv_role = c("test", "test", "test", "test"),
        Env = c("E1", "E1", "E2", "E2"),
        stringsAsFactors = FALSE
      ),
      yprob_cv_Reps_all = data.frame(
        Prob_A = c(0.8, 0.2, 0.4, 0.3),
        Prob_B = c(0.2, 0.8, 0.6, 0.7)
      )
    )
  )

  res <- PredictProR:::multi_env_predicted_vs_observed_result_plots_process(
    results = results,
    pheno_data = data.frame(
      GID = c("g1", "g2", "g3", "g4"),
      Class = c("A", "B", "A", "B"),
      Env = c("E1", "E1", "E2", "E2"),
      stringsAsFactors = FALSE
    ),
    heter_groups = "Env",
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05
  )

  expect_true("Class" %in% names(res$predicted_vs_observed_plots))
  expect_true(all(c("ModelA_E1", "ModelA_E2") %in% names(res$predicted_vs_observed_plots$Class$predicted_vs_observed_plots)))
  expect_true(inherits(res$predicted_vs_observed_plots$Class$predicted_vs_observed_plots$ModelA_E1, "gtable"))
  expect_true("Class" %in% names(res$mod_res_per_trait_per_model$ModelA_E1))
  expect_true(all(c("Prob_A", "Prob_B") %in% names(res$mod_res_per_trait_per_model$ModelA_E2$Class)))
  expect_true("ModelA_E1__Class" %in% names(res$classification_diagnostic_summaries))
  expect_true(is.data.frame(res$classification_diagnostic_summaries$ModelA_E1__Class$confusion_matrix))
})

test_that("predicted_vs_observed_ranking_plot avoids fragile loess fitting on tiny samples", {
  expect_no_warning(
    res <- PredictProR:::predicted_vs_observed_ranking_plot(
      observed_value = c(1, 2),
      predicted_value = c(1.1, 1.9),
      abs_very_close_threshold = 0.01,
      abs_close_threshold = 0.05
    )
  )

  expect_true(inherits(res, "gtable"))
})
