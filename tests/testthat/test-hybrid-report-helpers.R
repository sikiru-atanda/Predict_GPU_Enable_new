hybrid_report_call <- function(...) {
  if (exists("hybrid_benchmark_report", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("hybrid_benchmark_report", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::hybrid_benchmark_report(...))
  }
  stop("hybrid_benchmark_report is unavailable for this test.", call. = FALSE)
}

hybrid_process_call <- function(...) {
  if (exists("hybrid_benchmark_process", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("hybrid_benchmark_process", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::hybrid_benchmark_process(...))
  }
  stop("hybrid_benchmark_process is unavailable for this test.", call. = FALSE)
}

test_that("hybrid_benchmark_report returns ranked tables and plots", {
  skip_if_not_installed("ggplot2")

  summary_tab <- data.frame(
    family = c("hybrid_ml", "hybrid_ml", "hybrid_ml", "hybrid_ml"),
    model = c("Ridge_Regression", "SupportVectorMachine", "Ridge_Regression", "SupportVectorMachine"),
    split = c("masked_true_prediction", "masked_true_prediction", "Hybrid_One_New_Parent", "Hybrid_One_New_Parent"),
    rmse = c(1.2, 1.0, 0.9, 1.1),
    mae = c(1.0, 0.8, 0.7, 0.9),
    stringsAsFactors = FALSE
  )

  out <- hybrid_report_call(
    summary_table = summary_tab,
    masked_split_name = "masked_true_prediction",
    cv_split_name = "Hybrid_One_New_Parent",
    export = FALSE
  )

  expect_true(is.data.frame(out$summary))
  expect_true(is.data.frame(out$ranked_masked_true_prediction))
  expect_true(is.data.frame(out$ranked_cv))
  expect_true(is.data.frame(out$recommendations))
  expect_true(inherits(out$plots$rmse_by_model, "ggplot"))
  expect_true(inherits(out$plots$mae_by_model, "ggplot"))
  expect_identical(out$ranked_masked_true_prediction$model[1], "SupportVectorMachine")
  expect_identical(out$ranked_cv$model[1], "Ridge_Regression")
})

test_that("hybrid_benchmark_report exports CSV and PDF artifacts", {
  skip_if_not_installed("ggplot2")

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-hybrid-benchmark-report-")
  dir.create(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  summary_tab <- data.frame(
    family = c("hybrid_ml", "hybrid_ml"),
    model = c("RandomForest", "Ridge_Regression"),
    split = c("masked_true_prediction", "Hybrid_One_New_Parent"),
    rmse = c(1.1, 0.9),
    mae = c(0.95, 0.8),
    stringsAsFactors = FALSE
  )

  out <- hybrid_report_call(
    summary_table = summary_tab,
    masked_split_name = "masked_true_prediction",
    cv_split_name = "Hybrid_One_New_Parent",
    output_dir = tmp_dir,
    export = TRUE
  )

  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_summary.csv")))
  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_ranked_masked_true_prediction.csv")))
  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_ranked_cv.csv")))
  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_recommendations.csv")))
  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_rmse_by_model.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "hybrid_shortlist_mae_by_model.pdf")))
  expect_true(inherits(out$plots$rmse_by_model, "ggplot"))
})

test_that("hybrid_benchmark_table_from_runs assembles table from run objects", {
  true_run <- list(
    model_results = list(
      predicted_values = data.frame(
        HybridID = c("H1", "H2", "H3"),
        Predicted_value = c(10.1, 8.9, 7.8),
        Train_Test_Label = c("Train", "Test", "Test"),
        stringsAsFactors = FALSE
      )
    ),
    summary_statistic = list(
      summary_statistics = data.frame(
        stat = c("Mode", "Model_Type"),
        summary = c("hybrid_ml", "RandomForest"),
        stringsAsFactors = FALSE
      )
    )
  )

  cv_run <- list(
    cv_results_processed = list(
      hybrid_metric_summary = data.frame(
        cv_scenario = c("Hybrid_One_New_Parent", "Hybrid_Both_New_Parents"),
        root_mean_squared_error = c(0.9, 1.2),
        mean_absolute_error = c(0.7, 1.0),
        stringsAsFactors = FALSE
      )
    ),
    summary_statistic = list(
      summary_statistics = data.frame(
        stat = c("Mode", "Model_Type"),
        summary = c("hybrid_ml", "RandomForest"),
        stringsAsFactors = FALSE
      )
    )
  )

  truth_lookup <- data.frame(
    HybridID = c("H2", "H3"),
    Observed_truth = c(9.4, 7.6),
    stringsAsFactors = FALSE
  )

  tab <- hybrid_benchmark_table_from_runs(
    run_entries = list(
      rf = list(
        true_prediction = true_run,
        cv = cv_run,
        cv_split_name = "Hybrid_One_New_Parent"
      )
    ),
    truth_lookup = truth_lookup
  )

  expect_true(all(c("family", "model", "split", "rmse", "mae") %in% names(tab)))
  expect_equal(nrow(tab), 2)
  expect_true(any(tab$split == "masked_true_prediction"))
  expect_true(any(tab$split == "Hybrid_One_New_Parent"))
  expect_true(all(tab$family == "hybrid_ml"))
  expect_true(all(tab$model == "RandomForest"))
})

test_that("hybrid_benchmark_process returns processed hybrid comparison bundle", {
  skip_if_not_installed("ggplot2")

  true_run <- list(
    model_results = list(
      predicted_values = data.frame(
        HybridID = c("H1", "H2", "H3"),
        Predicted_value = c(10.1, 8.9, 7.8),
        Train_Test_Label = c("Train", "Test", "Test"),
        stringsAsFactors = FALSE
      )
    ),
    summary_statistic = list(
      summary_statistics = data.frame(
        stat = c("Mode", "Model_Type"),
        summary = c("hybrid_ml", "RandomForest"),
        stringsAsFactors = FALSE
      )
    )
  )

  cv_run <- list(
    cv_results_processed = list(
      hybrid_metric_summary = data.frame(
        cv_scenario = "Hybrid_One_New_Parent",
        root_mean_squared_error = 0.9,
        mean_absolute_error = 0.7,
        stringsAsFactors = FALSE
      )
    ),
    summary_statistic = list(
      summary_statistics = data.frame(
        stat = c("Mode", "Model_Type"),
        summary = c("hybrid_ml", "RandomForest"),
        stringsAsFactors = FALSE
      )
    )
  )

  truth_lookup <- data.frame(
    HybridID = c("H2", "H3"),
    Observed_truth = c(9.4, 7.6),
    stringsAsFactors = FALSE
  )

  out <- hybrid_process_call(
    run_entries = list(
      rf = list(
        true_prediction = true_run,
        cv = cv_run,
        cv_split_name = "Hybrid_One_New_Parent"
      )
    ),
    truth_lookup = truth_lookup
  )

  expect_true(is.data.frame(out$hybrid_benchmark_table))
  expect_true(is.data.frame(out$hybrid_benchmark_summary))
  expect_true(is.data.frame(out$hybrid_benchmark_ranked_masked_true_prediction))
  expect_true(is.data.frame(out$hybrid_benchmark_ranked_cv))
  expect_true(is.data.frame(out$hybrid_benchmark_recommendations))
  expect_true(inherits(out$hybrid_benchmark_plots$rmse_by_model, "ggplot"))
  expect_true(inherits(out$hybrid_benchmark_plots$mae_by_model, "ggplot"))
})
