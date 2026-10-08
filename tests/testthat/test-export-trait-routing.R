test_that("gp_finalize_model_execute_results keeps per-trait exports separated", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-export-routing-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  results <- list(
    Yield = list(
      GS_model = "BayesA",
      res_model_output = list(
        bayes_result = list(
          Predicted_value = data.frame(
            GID = c("g1", "g2", "g3"),
            Predicted_value = c(1, 2, 3),
            Train_Test_Label = c("Train", "Test", "Train"),
            stringsAsFactors = FALSE
          )
        )
      ),
      res_summary_stat = list(
        summary_statistics = data.frame(stat = "trait", summary = "Yield", stringsAsFactors = FALSE)
      )
    ),
    Height = list(
      GS_model = "BayesA",
      res_model_output = list(
        bayes_result = list(
          Predicted_value = data.frame(
            GID = c("g1", "g2", "g3"),
            Predicted_value = c(4, 5, 6),
            Train_Test_Label = c("Test", "Train", "Train"),
            stringsAsFactors = FALSE
          )
        )
      ),
      res_summary_stat = list(
        summary_statistics = data.frame(stat = "trait", summary = "Height", stringsAsFactors = FALSE)
      )
    )
  )

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = results,
    GS_model = "GS_model",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(
      list(
        trait = "Yield",
        model = "BayesA",
        rep = 1L,
        ypred_cv_Reps_all = data.frame(y = 1, yhat = 1)
      ),
      list(
        trait = "Height",
        model = "BayesA",
        rep = 2L,
        ypred_cv_Reps_all = data.frame(y = 2, yhat = 2)
      )
    ),
    run_metadata = data.frame(
      key = c("context", "policy_backend"),
      value = c("model_execute", "sequential"),
      stringsAsFactors = FALSE
    ),
    system_database = FALSE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_true(all(c(
    "cv_results_raw",
    "cv_results_processed",
    "model_results_by_trait",
    "model_results_by_model",
    "run_metadata"
  ) %in% names(out)))
  expect_false(any(c("Yield", "Height") %in% names(out)))
  expect_identical(out$model_results_by_trait$Yield$BayesA$predicted_values$Train_Test_Label, c("Train", "Test", "Train"))
  expect_identical(out$model_results_by_trait$Height$BayesA$predicted_values$Train_Test_Label, c("Test", "Train", "Train"))
  expect_identical(out$run_metadata$key, c("context", "policy_backend"))
  expect_identical(out$run_metadata$value, c("model_execute", "sequential"))
  expect_identical(names(out$cv_results_raw), c("Yield_BayesA_rep1", "Height_BayesA_rep2"))
  expect_identical(names(out$model_results_by_trait), c("Yield", "Height"))
  expect_true("BayesA" %in% names(out$model_results_by_model))
  expect_true("Yield" %in% names(out$model_results_by_model$BayesA))

  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 2)

  yield_dir <- created_dirs[basename(created_dirs) %in% basename(created_dirs)[grepl("^Yield_", basename(created_dirs))]]
  height_dir <- created_dirs[basename(created_dirs) %in% basename(created_dirs)[grepl("^Height_", basename(created_dirs))]]
  expect_length(yield_dir, 1)
  expect_length(height_dir, 1)

  yield_csv <- utils::read.csv(file.path(yield_dir, "Predicted_Value.csv"), stringsAsFactors = FALSE)
  height_csv <- utils::read.csv(file.path(height_dir, "Predicted_Value.csv"), stringsAsFactors = FALSE)
  yield_meta <- utils::read.csv(file.path(yield_dir, "Run_metadata.csv"), stringsAsFactors = FALSE)
  height_meta <- utils::read.csv(file.path(height_dir, "Run_metadata.csv"), stringsAsFactors = FALSE)

  expect_identical(yield_csv$Train_Test_Label, c("Train", "Test", "Train"))
  expect_identical(height_csv$Train_Test_Label, c("Test", "Train", "Train"))
  expect_identical(names(yield_csv), PredictProR:::gp_gaussian_prediction_columns())
  expect_identical(names(height_csv), PredictProR:::gp_gaussian_prediction_columns())
  expect_identical(yield_meta$key, c("context", "policy_backend"))
  expect_identical(height_meta$value, c("model_execute", "sequential"))
})

test_that("CV-only outputs expose predicted-vs-observed diagnostic plots", {
  cv_plot <- ggplot2::ggplot(
    data.frame(observed = c(1, 2), predicted = c(1.1, 1.9)),
    ggplot2::aes(observed, predicted)
  ) + ggplot2::geom_point()
  cv_plot_bundle <- list(
    Yield = list(
      predicted_vs_observed_plots = list(KRR_E1 = cv_plot)
    )
  )
  cv_model_results <- list(
    KRR_E1 = list(
      Yield = data.frame(y = c(1, 2), yhat = c(1.1, 1.9))
    )
  )

  out <- PredictProR:::results_handling(
    GS_model = NULL,
    res_model_output = NULL,
    res_summary_stat = NULL,
    res_plot = NULL,
    res_plot_mean = NULL,
    res_plot_result_diagnostic = NULL,
    test_diagonistic_plots = NULL,
    res_plot_result_diagnostic_cv_only = cv_plot_bundle,
    res_mod_results_cv_per_trait_model = cv_model_results,
    cv_results_processed = list(best_models_list = list()),
    cv_results_raw = list(
      list(
        trait = "Yield",
        model = "KRR",
        rep = 1L,
        ypred_cv_Reps_all = data.frame(y = c(1, 2), yhat = c(1.1, 1.9))
      )
    ),
    system_database = TRUE,
    plot_filename = "CV_results",
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_identical(out$res_plot_result_diagnostic_cv_only, cv_plot_bundle)
  expect_identical(
    out$cv_results_predicted_vs_observed$predicted_vs_observed_plots,
    cv_plot_bundle
  )
  expect_identical(
    out$cv_results_predicted_vs_observed$mod_res_per_trait_per_model,
    cv_model_results
  )
})

test_that("gp_finalize_model_execute_results keeps all final models nested by trait and model", {
  make_result <- function(trait, model, values) {
    list(
      trait = trait,
      GS_model = model,
      res_model_output = list(
        predicted_values = data.frame(
          GID = c("g1", "g2"),
          Predicted_value = values,
          Train_Test_Label = c("Train", "Test"),
          Standard_error = c(0.1, 0.2),
          PEV = c(0.01, 0.04),
          Reliability = c(0.9, 0.8),
          Observed_value = c(values[[1]], NA),
          stringsAsFactors = FALSE
        )
      ),
      res_summary_stat = list(
        summary_statistics = data.frame(
          stat = c("trait", "model"),
          summary = c(trait, model),
          stringsAsFactors = FALSE
        )
      )
    )
  }

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(
      Yield_Kernel = make_result("Yield", "Kernel-GBLUP", c(1.1, 1.2)),
      Yield_Scalable = make_result("Yield", "Scalable-GBLUP", c(1.3, 1.4)),
      Height_Kernel = make_result("Height", "Kernel-GBLUP", c(2.1, 2.2)),
      Height_Scalable = make_result("Height", "Scalable-GBLUP", c(2.3, 2.4))
    ),
    GS_model = "GS_model",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = NULL,
    system_database = TRUE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  # Phase 3.26 (revert of 3.25 dual-storage): the user-supplied label is
  # the SOLE list key. R's [[name]] and backtick-`name` accessors handle
  # dashes, so no make.names() alias is needed. Duplicate keys in
  # names(...) are a regression, not a feature.
  expect_setequal(names(out$model_results_by_model),
                  c("Kernel-GBLUP", "Scalable-GBLUP"))
  expect_setequal(names(out$model_results_by_model[["Kernel-GBLUP"]]), c("Yield", "Height"))
  expect_setequal(names(out$model_results_by_model[["Scalable-GBLUP"]]), c("Yield", "Height"))
  expect_setequal(names(out$model_results_by_trait), c("Yield", "Height"))
  expect_setequal(names(out$model_results_by_trait$Yield),
                  c("Kernel-GBLUP", "Scalable-GBLUP"))
  expect_setequal(names(out$model_results_by_trait$Height),
                  c("Kernel-GBLUP", "Scalable-GBLUP"))
  expect_false(any(c("Yield", "Height") %in% names(out)))
  expect_true("predicted_values" %in% names(out$model_results_by_trait$Yield[["Kernel-GBLUP"]]))
  expect_identical(
    names(out$model_results_by_trait$Yield[["Kernel-GBLUP"]]$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns()
  )
  expect_identical(
    intersect(c("Predicted_value", "diagnostic_tst_plot", "Variance_components"),
              names(out$model_results_by_trait$Yield[["Kernel-GBLUP"]])),
    character()
  )
  expect_true(all(c("model_parameters", "predicted_values", "diagnostic_plots", "variance_components") %in%
                    names(out$model_results_by_trait$Yield[["Kernel-GBLUP"]])))
  expect_false("raw_python_result" %in% names(out$Yield$models[["Kernel-GBLUP"]]))
  expect_false("gp_info" %in% names(out$Yield$models[["Kernel-GBLUP"]]))
})

test_that("public finalizer gives ML DL Bayes GBLUP and GP models the same user-facing result contract", {
  make_result <- function(model, payload) {
    list(
      trait = "Yield",
      GS_model = model,
      res_model_output = payload,
      res_summary_stat = list(
        summary_statistics = data.frame(stat = "model", summary = model, stringsAsFactors = FALSE)
      )
    )
  }
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 1.2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.0, NA),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    Reliability = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(
      make_result("RandomForest", list(predicted_values = pred, raw_python_result = list(debug = TRUE))),
      make_result("mlp", list(predicted_values = pred, raw_python_result = list(debug = TRUE))),
      make_result("BayesB", list(bayes_result = list(Predicted_value = pred))),
      make_result("GBLUP", list(Predicted_value = pred, Asreml_model = list())),
      make_result("Kernel-GBLUP", list(predicted_values = pred, gp_info = list(debug = TRUE)))
    ),
    GS_model = "GS_model",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = NULL,
    system_database = TRUE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expected_names <- c(
    "model_parameters",
    "predicted_values",
    "diagnostic_plots",
    "variance_components"
  )
  expected_cols <- PredictProR:::gp_gaussian_prediction_columns()
  for (model_key in c("RandomForest", "mlp", "BayesB", "GBLUP", "Kernel-GBLUP")) {
    expect_true(model_key %in% names(out$model_results_by_model))
    res <- out$model_results_by_model[[model_key]]$Yield
    expect_identical(names(res), expected_names)
    expect_identical(names(res$predicted_values), expected_cols)
    expect_false("raw_python_result" %in% names(res))
    expect_false("gp_info" %in% names(res))
    expect_false("covariance_components" %in% names(res))
    expect_false("correlation_components" %in% names(res))
  }
})

test_that("single-environment model_results_by_model predictions expose only the common gaussian columns", {
  pred <- data.frame(
    GID = c("g1", "g2", "g3"),
    Env = c(NA_character_, NA_character_, NA_character_),
    Trait = c(NA_character_, NA_character_, NA_character_),
    Predicted_value = c(1.1, 2.2, 3.3),
    Train_Test_Label = c("Train", "Test", "Train"),
    Observed_value = c(1, NA, 3),
    Standard_error = c(0.1, 0.2, 0.1),
    PEV = c(0.01, 0.04, 0.01),
    Prediction_error_variance = c(0.01, 0.04, 0.01),
    Genetic_variance = c(1, 1, 1),
    Reliability_reference_variance = c(1, 1, 1),
    Reliability = c(0.9, 0.8, 0.85),
    Reliability_remarks = c("Reliable", "Acceptable", "Acceptable"),
    Reliability_percentage = c(100, 100, 100),
    Reliability_basis = c("internal", "internal", "internal"),
    Train_Test = c("Train", "Test", "Train"),
    Prediction_stability = c(0.9, 0.8, 0.85),
    Prediction_stability_remarks = c("Stable", "Stable", "Stable"),
    Predicted_class = NA_character_,
    Prediction_confidence = NA_real_,
    Classification_uncertainty = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(
      list(
        trait = "HDDT",
        GS_model = "BRR",
        res_model_output = list(predicted_values = pred),
        res_summary_stat = list(summary_statistics = data.frame(stat = "trait", summary = "HDDT"))
      )
    ),
    GS_model = "BRR",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = NULL,
    system_database = TRUE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_identical(
    names(out$model_results_by_model$BRR$HDDT$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns()
  )
  expect_false("Env" %in% names(out$model_results_by_model$BRR$HDDT$predicted_values))
  expect_false("Trait" %in% names(out$model_results_by_model$BRR$HDDT$predicted_values))
  expect_false("Prediction_error_variance" %in% names(out$model_results_by_model$BRR$HDDT$predicted_values))
  expect_false("Predicted_class" %in% names(out$model_results_by_model$BRR$HDDT$predicted_values))
})

test_that("finalizer carries custom genotype ID name into public prediction table", {
  pred <- data.frame(
    NAME = c("g1", "g2", "g3"),
    Predicted_value = c(1.1, 2.2, 3.3),
    Train_Test_Label = c("Train", "Test", "Train"),
    Observed_value = c(1, NA, 3),
    Standard_error = c(0.1, 0.2, 0.1),
    PEV = c(0.01, 0.04, 0.01),
    Reliability = c(0.9, 0.8, 0.85),
    Reliability_remarks = c("Reliable", "Acceptable", "Reliable"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(
      list(
        trait = "HDDT",
        response = "HDDT",
        GS_model = "BRR",
        res_model_output = list(predicted_values = pred),
        res_summary_stat = list(summary_statistics = data.frame(stat = "trait", summary = "HDDT"))
      )
    ),
    GS_model = "BRR",
    gen_name = "NAME",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = NULL,
    system_database = TRUE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  returned <- out$model_results_by_model$BRR$HDDT$predicted_values
  # The finalizer must carry the user-supplied gen_name through to the public
  # output column name. Here the caller passed gen_name = "NAME" (above) so
  # the returned table's first column is "NAME", not a hardcoded "GID".
  expect_identical(names(returned), PredictProR:::gp_gaussian_prediction_columns(gen_name = "NAME"))
  expect_identical(returned$NAME, c("g1", "g2", "g3"))
})

test_that("nested Bayesian predictions are compact in returned objects and exported Predicted_Value", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-bayes-compact-output-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pred <- data.frame(
    GID = c("g1", "g2"),
    Env = c(NA_character_, NA_character_),
    Trait = c(NA_character_, NA_character_),
    Predicted_value = c(1.1, 2.2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1, NA),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    Prediction_error_variance = c(0.01, 0.04),
    Genetic_variance = c(1, 1),
    Reliability_reference_variance = c(1, 1),
    Reliability = c(0.9, 0.8),
    Reliability_remarks = c("Reliable", "Acceptable"),
    Reliability_percentage = c(100, 100),
    Reliability_basis = c("internal", "internal"),
    Train_Test = c("Train", "Test"),
    Prediction_stability = c(0.9, 0.8),
    Prediction_stability_remarks = c("Stable", "Stable"),
    Predicted_class = NA_character_,
    Prediction_confidence = NA_real_,
    Classification_uncertainty = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(
      list(
        trait = "HDDT",
        GS_model = "BRR",
        res_model_output = list(bayes_result = list(Predicted_value = pred)),
        res_summary_stat = list(summary_statistics = data.frame(stat = "trait", summary = "HDDT"))
      )
    ),
    GS_model = "BRR",
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = NULL,
    system_database = FALSE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expected_cols <- PredictProR:::gp_gaussian_prediction_columns()
  returned <- out$model_results_by_model$BRR$HDDT$predicted_values
  expect_identical(names(returned), expected_cols)
  expect_false("Env" %in% names(returned))
  expect_false("Trait" %in% names(returned))
  expect_false("Prediction_error_variance" %in% names(returned))
  expect_false("Predicted_class" %in% names(returned))

  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)
  exported <- utils::read.csv(file.path(created_dirs, "Predicted_Value.csv"), stringsAsFactors = FALSE)
  expect_identical(names(exported), expected_cols)
})

test_that("Predicted_Value exports keep the common cross-model schema", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-predicted-value-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  out <- PredictProR:::results_handling(
    GS_model = "Kernel-GBLUP",
    res_model_output = list(
      predicted_values = data.frame(
        GID = c("g1", "g2"),
        Predicted_value = c(1.5, 2.5),
        Train_Test_Label = c("Train", "Test"),
        Standard_error = c(0.1, 0.2),
        PEV = c(0.01, 0.04),
        Genetic_variance = c(3, 3),
        Reliability = c(0.9, 0.8),
        Observed_value = c(1.4, NA),
        stringsAsFactors = FALSE
      )
    ),
    system_database = FALSE,
    plot_filename = "Yield"
  )

  expect_false(is.null(out))
  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)
  pred_csv <- utils::read.csv(file.path(created_dirs[[1]], "Predicted_Value.csv"), stringsAsFactors = FALSE)
  expect_identical(
    names(pred_csv),
    PredictProR:::gp_gaussian_prediction_columns()
  )
})

test_that("results_handling exports run metadata for cv-only outputs", {
  skip_if_not_installed("ggplot2")

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-cv-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  run_metadata <- data.frame(
    key = c("context", "policy_backend"),
    value = c("cross_validation", "sequential"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::results_handling(
    GS_model = NULL,
    res_model_output = NULL,
    res_summary_stat = NULL,
    res_plot = ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
    res_plot_mean = NULL,
    test_diagonistic_plots = NULL,
    res_plot_result_diagnostic = NULL,
    res_plot_result_diagnostic_cv_only = NULL,
    res_mod_results_cv_per_trait_model = list(),
    cv_results_processed = list(
      best_models_list = list(),
      run_metadata = run_metadata,
      gaussian_uncertainty_summaries = list(
        aggregated_uncertainty = data.frame(
          trait = "Yield",
          model = "mlp",
          n = 2L,
          mean_observed = 1.5,
          mean_predicted = 1.4,
          mean_absolute_error = 0.2,
          mean_squared_error = 0.05,
          root_mean_squared_error = 0.2236068,
          residual_sd = 0.1,
          mean_standard_error = 0.3,
          mean_pev = 0.09,
          mean_interval_width = 0.8,
          empirical_coverage = 1,
          stringsAsFactors = FALSE
        )
      ),
      gaussian_risk_summaries = list(
        aggregated_risk = data.frame(
          trait = "Yield",
          model = "mlp",
          n = 2L,
          mean_prediction_confidence = 0.55,
          mean_prediction_unfamiliarity = 0.20,
          mean_prediction_risk_score = 0.45,
          mean_rank_instability_risk = 0.45,
          low_risk_percentage = 25,
          moderate_risk_percentage = 50,
          high_risk_percentage = 25,
          low_unfamiliarity_percentage = 25,
          moderate_unfamiliarity_percentage = 50,
          high_unfamiliarity_percentage = 0,
          prediction_unfamiliarity_basis = "predictor_space_novelty",
          prediction_risk_basis = "Heldout_Training_OOB",
          rank_instability_risk_basis = "Heldout_Training_OOB",
          stringsAsFactors = FALSE
        )
      ),
      gaussian_risk_plots = list(
        mean_prediction_unfamiliarity = ggplot2::ggplot(
          data.frame(model = "mlp", mean_prediction_unfamiliarity = 0.2),
          ggplot2::aes(model, mean_prediction_unfamiliarity)
        ) + ggplot2::geom_col(),
        mean_rank_instability_risk = ggplot2::ggplot(
          data.frame(model = "mlp", mean_rank_instability_risk = 0.45),
          ggplot2::aes(model, mean_rank_instability_risk)
        ) + ggplot2::geom_col(),
        risk_distribution = ggplot2::ggplot(
          data.frame(model = "mlp", risk_band = "Moderate Risk", percentage = 50),
          ggplot2::aes(model, percentage, fill = risk_band)
        ) + ggplot2::geom_col()
      ),
      classification_probability_summaries = list(
        aggregated_probabilities = data.frame(
          trait = "Yield",
          model = "mlp",
          row_id = 1L,
          Prob_no = 0.2,
          Prob_yes = 0.8,
          stringsAsFactors = FALSE
        )
      ),
      classification_diagnostic_summaries = list(
        `mlp__Yield` = list(
          confusion_matrix = data.frame(
            trait = "Yield",
            model = "mlp",
            cv_role = "test",
            Observed = "yes",
            Predicted = "yes",
            Freq = 1L,
            stringsAsFactors = FALSE
          ),
          calibration_summary = data.frame(
            trait = "Yield",
            model = "mlp",
            cv_role = "test",
            bin = "[0.8,0.9]",
            n = 1L,
            mean_confidence = 0.8,
            observed_rate = 1,
            stringsAsFactors = FALSE
          )
        )
      )
    ),
    cv_results_raw = list(
      list(
        response = "Yield",
        model = "mlp",
        rep = 1L,
        eval_metrics_reps = data.frame(Rep = 1L, mean_squared_error = 0.1),
        ypred_cv_Reps_all = data.frame(y = 1, yhat = 1),
        yprob_cv_Reps_all = data.frame(no = 0.2, yes = 0.8)
      )
    ),
    run_metadata = run_metadata,
    system_database = FALSE,
    plot_filename = "CV_results",
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_identical(out$run_metadata$key, c("context", "policy_backend"))
  expect_identical(out$cv_results_processed$run_metadata$value, c("cross_validation", "sequential"))

  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)

  meta <- utils::read.csv(file.path(created_dirs[[1]], "Run_metadata.csv"), stringsAsFactors = FALSE)
  yprob <- utils::read.csv(
    file.path(created_dirs[[1]], "cv_results_raw", "Yield_mlp_rep1_yprob_cv_Reps_all.csv"),
    stringsAsFactors = FALSE
  )
  yprob_summary <- utils::read.csv(
    file.path(created_dirs[[1]], "cv_results_processed_classification_probability_summaries_aggregated_probabilities.csv"),
    stringsAsFactors = FALSE
  )
  gaussian_summary <- utils::read.csv(
    file.path(created_dirs[[1]], "cv_results_processed_gaussian_uncertainty_summaries_aggregated_uncertainty.csv"),
    stringsAsFactors = FALSE
  )
  gaussian_risk <- utils::read.csv(
    file.path(created_dirs[[1]], "cv_results_processed_gaussian_risk_summaries_aggregated_risk.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(file.exists(file.path(created_dirs[[1]], "CV_gaussian_mean_prediction_unfamiliarity.pdf")))
  expect_true(file.exists(file.path(created_dirs[[1]], "CV_gaussian_mean_rank_instability_risk.pdf")))
  expect_true(file.exists(file.path(created_dirs[[1]], "CV_gaussian_risk_distribution.pdf")))
  conf_summary <- utils::read.csv(
    file.path(created_dirs[[1]], "cv_results_processed_classification_diagnostic_summaries_mlp__Yield_confusion_matrix.csv"),
    stringsAsFactors = FALSE
  )
  expect_identical(meta$key, c("context", "policy_backend"))
  expect_identical(meta$value, c("cross_validation", "sequential"))
  expect_identical(names(yprob), c("no", "yes"))
  expect_equal(as.numeric(yprob[1, ]), c(0.2, 0.8))
  expect_true(all(c("mean_standard_error", "mean_interval_width", "empirical_coverage") %in% names(gaussian_summary)))
  expect_equal(as.numeric(gaussian_summary[1, c("mean_standard_error", "mean_interval_width", "empirical_coverage")]), c(0.3, 0.8, 1))
  expect_true(all(c(
    "mean_prediction_confidence",
    "mean_prediction_unfamiliarity",
    "mean_prediction_risk_score",
    "mean_rank_instability_risk",
    "prediction_unfamiliarity_basis",
    "prediction_risk_basis",
    "rank_instability_risk_basis"
  ) %in% names(gaussian_risk)))
  expect_equal(as.numeric(gaussian_risk[1, c("mean_prediction_confidence", "mean_prediction_unfamiliarity", "mean_prediction_risk_score", "mean_rank_instability_risk")]), c(0.55, 0.2, 0.45, 0.45))
  expect_identical(gaussian_risk$prediction_unfamiliarity_basis, "predictor_space_novelty")
  expect_identical(gaussian_risk$prediction_risk_basis, "Heldout_Training_OOB")
  expect_true(all(c("Prob_no", "Prob_yes") %in% names(yprob_summary)))
  expect_equal(as.numeric(yprob_summary[1, c("Prob_no", "Prob_yes")]), c(0.2, 0.8))
  expect_true(all(c("Observed", "Predicted", "Freq") %in% names(conf_summary)))
  expect_identical(conf_summary$Observed, "yes")
})

test_that("results_handling exports nested GP metadata without coercion errors", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-gp-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 2.2),
    Train_Test_Label = c("Train", "Test"),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    stringsAsFactors = FALSE
  )

  expect_no_error(
    out <- PredictProR:::results_handling(
      GS_model = "Kernel-GBLUP",
      res_model_output = list(
        predicted_values = pred,
        model_parameters = data.frame(stat = "gp_model", summary = "KRR"),
        gp_info = list(
          backend = "r",
          empty_branch = list(),
          nested = list(iterations = integer(), converged = TRUE)
        ),
        gp_result = list(predictions = pred, diagnostics = list())
      ),
      res_summary_stat = list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
      system_database = FALSE,
      plot_filename = "gp_export",
      plot_extension = "pdf"
    )
  )

  expect_identical(out$export_status, "Successful")
  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)
  expect_true(file.exists(file.path(created_dirs, "Predicted_Value.csv")))
  expect_true(file.exists(file.path(created_dirs, "model_parameters.csv")))
  expect_false(file.exists(file.path(created_dirs, "Prediction_Details.csv")))
})

test_that("nested Bayesian model parameters survive the public export contract", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-bayes-parameter-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 2.2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.0, NA_real_),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    stringsAsFactors = FALSE
  )
  params <- data.frame(
    stat = c("stage2_observation_weighting", "bglr_native_weight_transform"),
    summary = c("precision_via_BGLR_inverse_squared_weight", "sqrt(user_precision)"),
    stringsAsFactors = FALSE
  )

  standardized <- PredictProR:::gp_standardize_public_model_result(
    list(bayes_result = list(Predicted_value = pred, model_parameters = params)),
    gen_name = "GID"
  )
  expect_identical(standardized$model_parameters, params)

  out <- PredictProR:::results_handling(
    GS_model = "GBLUP_BRR",
    res_model_output = standardized,
    system_database = FALSE,
    plot_filename = "bayes_parameter_export"
  )
  expect_identical(out$export_status, "Successful")
  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1L)
  native <- utils::read.csv(
    file.path(created_dirs, "model_parameters.csv"),
    stringsAsFactors = FALSE
  )
  expect_identical(native, params)
})

test_that("results_handling exports one common prediction table for ML outputs", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-ml-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  out <- PredictProR:::results_handling(
    GS_model = "Ridge_Regression",
    res_model_output = list(
      model_parameters = data.frame(stat = "lambda", summary = 1, stringsAsFactors = FALSE),
      predicted_values = data.frame(
        GID = c("g1", "g2"),
        Predicted_value = c(1.2, 2.3),
        Train_Test_Label = c("Train", "Test"),
        Standard_error = c(0.1, 0.2),
        PEV = c(0.01, 0.04),
        Prediction_stability = c(0.95, 0.80),
        Prediction_stability_remarks = c("Stable", "Moderately Stable"),
        Prediction_stability_percentage = c(100, 100),
        Prediction_confidence = c(0.7, 0.6),
        Prediction_confidence_remarks = c("High Confidence", "Moderate Confidence"),
        Prediction_risk_score = c(0.3, 0.4),
        Prediction_risk_remarks = c("Low Risk", "Moderate Risk"),
        Prediction_risk_basis = c("Heldout_Training_OOB", "Calibrated_Final_Prediction"),
        Rank_instability_risk = c(0.3, 0.4),
        Rank_instability_risk_remarks = c("Low Risk", "Moderate Risk"),
        Rank_instability_risk_basis = c("Heldout_Training_OOB", "Calibrated_Final_Prediction"),
        stringsAsFactors = FALSE
      )
    ),
    res_summary_stat = list(
      summary_statistics = data.frame(stat = "trait", summary = "Yield", stringsAsFactors = FALSE)
    ),
    system_database = FALSE,
    plot_filename = "Yield",
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_false(is.null(out$model_results))

  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)

  pred_common <- utils::read.csv(file.path(created_dirs[[1]], "Predicted_Value.csv"), stringsAsFactors = FALSE)
  expect_identical(names(pred_common), PredictProR:::gp_gaussian_prediction_columns())
  expect_false(file.exists(file.path(created_dirs[[1]], "Prediction_Details.csv")))
})

test_that("results_handling preserves hybrid identifiers in native prediction CSV", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-hybrid-id-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pred_df <- data.frame(
    HybridID = c("H001", "H002", "H003"),
    Female = c("F1", "F1", "F2"),
    Male = c("M1", "M2", "M1"),
    Observed_value = c(10, 8, NA_real_),
    Predicted_value = c(9.8, 8.1, 9.0),
    Train_Test_Label = c("Train", "Train", "Test"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::results_handling(
    GS_model = "RandomForest",
    res_model_output = list(predicted_values = pred_df),
    res_summary_stat = list(
      summary_statistics = data.frame(
        stat = "trait", summary = "Yield", stringsAsFactors = FALSE
      )
    ),
    system_database = FALSE,
    plot_filename = "hybrid_id_export"
  )

  expect_identical(out$export_status, "Successful")
  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1L)
  native <- utils::read.csv(
    file.path(created_dirs[[1L]], "Predicted_Value.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_identical(names(native), PredictProR:::gp_hybrid_prediction_columns())
  expect_identical(native$Hybrid_ID, pred_df$HybridID)
  expect_identical(native$Female, pred_df$Female)
  expect_identical(native$Male, pred_df$Male)
  expect_false("GID" %in% names(native))
  expect_false(anyNA(native[, c("Hybrid_ID", "Female", "Male")]))
})

test_that("results_handling exports hybrid summary tables and hybrid CV plots", {
  skip_if_not_installed("ggplot2")

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-hybrid-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pred_df <- data.frame(
    HybridID = c("H1", "H2", "H3"),
    Female = c("F1", "F1", "F2"),
    Male = c("M1", "M2", "M1"),
    Observed_value = c(10, NA, 8),
    Predicted_value = c(9.7, 9.1, 8.2),
    Train_Test_Label = c("Train", "Test", "Train"),
    Prediction_intercept = c(9, 9, 9),
    Female_additive_contribution = c(0.6, 0.6, -0.3),
    Male_additive_contribution = c(0.2, -0.1, 0.2),
    Hybrid_interaction_contribution = c(-0.1, 0.3, 0.1),
    stringsAsFactors = FALSE
  )
  summary_stat <- PredictProR:::summary_statistics_hybrid(
    predicted_object = pred_df,
    response = "Yield",
    female_parent = "Female",
    male_parent = "Male",
    mode = "hybrid_ml",
    model_type = "RandomForest"
  )

  out <- PredictProR:::results_handling(
    GS_model = "RandomForest",
    res_model_output = list(
      predicted_values = pred_df,
      diagnostic_plots = ggplot2::ggplot(
        pred_df[!is.na(pred_df$Observed_value), ],
        ggplot2::aes(Observed_value, Predicted_value)
      ) + ggplot2::geom_point()
    ),
    res_summary_stat = summary_stat,
    res_plot = NULL,
    res_plot_mean = NULL,
    test_diagonistic_plots = ggplot2::ggplot(
      pred_df[!is.na(pred_df$Observed_value), ],
      ggplot2::aes(Observed_value, Predicted_value)
    ) + ggplot2::geom_point(),
    res_plot_result_diagnostic = NULL,
    res_plot_result_diagnostic_cv_only = NULL,
    res_mod_results_cv_per_trait_model = NULL,
    cv_results_processed = list(
      hybrid_metric_summary = data.frame(
        cv_scenario = "One_New_Parent",
        root_mean_squared_error = 1.1,
        mean_absolute_error = 0.9,
        stringsAsFactors = FALSE
      ),
      hybrid_prediction_counts = data.frame(
        cv_scenario = "One_New_Parent",
        n_test_hybrids = 6L,
        stringsAsFactors = FALSE
      ),
      hybrid_cv_predictions = data.frame(
        HybridID = "H4",
        Observed_value = 10,
        Predicted_value = 9.2,
        cv_scenario = "One_New_Parent",
        stringsAsFactors = FALSE
      ),
      hybrid_cv_plot = ggplot2::ggplot(
        data.frame(Observed_value = 10, Predicted_value = 9.2, cv_scenario = "One_New_Parent"),
        ggplot2::aes(Observed_value, Predicted_value)
      ) + ggplot2::geom_point()
    ),
    cv_results_raw = NULL,
    run_metadata = data.frame(
      key = c("context", "policy_backend"),
      value = c("hybrid_ml", "sequential"),
      stringsAsFactors = FALSE
    ),
    system_database = FALSE,
    plot_filename = "hybrid_export",
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )

  expect_true(all(c("summary_statistics", "hybrid_component_summary", "hybrid_train_test_summary") %in% names(out$summary_statistic)))

  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1)
  export_dir <- created_dirs[[1]]

  expect_true(file.exists(file.path(export_dir, "summary_statistics.csv")))
  expect_true(file.exists(file.path(export_dir, "hybrid_component_summary.csv")))
  expect_true(file.exists(file.path(export_dir, "hybrid_train_test_summary.csv")))
  expect_true(file.exists(file.path(export_dir, "cv_results_processed_hybrid_metric_summary.csv")))
  expect_true(file.exists(file.path(export_dir, "cv_results_processed_hybrid_prediction_counts.csv")))
  expect_true(file.exists(file.path(export_dir, "CV_hybrid_observed_vs_predicted.pdf")))
  expect_true(file.exists(file.path(export_dir, "test_diagonistic_plots_GS_model.pdf")))
  pred_native <- utils::read.csv(
    file.path(export_dir, "Predicted_Value.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_identical(names(pred_native), PredictProR:::gp_hybrid_prediction_columns())
  expect_identical(pred_native$Hybrid_ID, pred_df$HybridID)
  expect_identical(pred_native$Female, pred_df$Female)
  expect_identical(pred_native$Male, pred_df$Male)
})
