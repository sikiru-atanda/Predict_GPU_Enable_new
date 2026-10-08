if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

mk_multitrait_ml_runtime_data <- function(n = 20L, p = 8L) {
  set.seed(42)
  gids <- paste0("g", seq_len(n))
  # Traits come from a latent continuous signal; geno_data carries it as 0/1/2
  # homozygous 0/2 dosages: genotype QC rejects values above 2 without ploidy and,
  # for line data, removes markers with more than 10% heterozygotes.
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- ifelse(Z > 0, 2, 0)   # homozygous 0/2 coding, like inbred lines
  storage.mode(X) <- "double"
  rownames(X) <- gids
  colnames(X) <- paste0("m", seq_len(p))
  base <- Z[, 1] * 0.7 - Z[, 2] * 0.3
  t1 <- base + rnorm(n, sd = 0.15)
  t2 <- base * 0.6 + Z[, 3] * 0.4 + rnorm(n, sd = 0.15)
  t3 <- -base * 0.5 + Z[, 4] * 0.6 + rnorm(n, sd = 0.15)
  t1[c(3, 7, 12)] <- NA
  t2[c(5, 9, 15)] <- NA
  t3[c(2, 11, 18)] <- NA
  list(
    pheno = data.frame(
      GID = gids,
      Trait1 = t1,
      Trait2 = t2,
      Trait3 = t3,
      stringsAsFactors = FALSE
    ),
    geno = X
  )
}

configure_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")
  invisible(py)
}

test_that("multi-trait gaussian ML runtime returns long and wide outputs", {
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model = "RandomForest",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    ntree = 100,
    internal_cv_nfolds = 3L,
    internal_cv_replication = 2L
  )

  preds <- out$model_results$predicted_values
  wide <- out$model_results$multitrait_prediction_wide
  counts <- out$model_results$multitrait_trait_counts

  expect_identical(
    names(preds),
    PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
  )
  expect_equal(sort(unique(preds$Trait)), c("Trait1", "Trait2", "Trait3"))
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_equal(nrow(wide), nrow(dat$pheno))
  expect_true(all(c("GID", "Trait1", "Trait2", "Trait3") %in% names(wide)))
  expect_equal(sort(counts$trait), c("Trait1", "Trait2", "Trait3"))
  expect_true(all(counts$predicted_count > 0))
  expect_true(any(vapply(
    out$model_results$diagnostic_plots,
    inherits,
    logical(1L),
    what = "ggplot"
  )))
  expect_true(all(is.finite(preds$Standard_error)))
  expect_equal(preds$PEV, preds$Standard_error^2)
  expect_true(all(grepl("not genetic reliability", preds$Reliability_basis, fixed = TRUE)))
  expect_true(is.matrix(out$model_results$Prediction_covariance_traits))
  expect_true(is.matrix(out$model_results$Prediction_error_covariance_traits))
  expect_true(is.matrix(out$model_results$Prediction_covariance_traits_SE))
  expect_false("Genetic_covariance_traits" %in% names(out$model_results))
})

test_that("multi-trait gaussian Ridge, PLS, Lasso, SVM, Xgboost, and CatBoost runtime return long and wide outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  for (mdl in c("Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost")) {
    out <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = c("Trait1", "Trait2", "Trait3"),
      gen_name = "GID",
      GS_model = mdl,
      response_family = "gaussian",
      multi_trait_ml = TRUE,
      system_database = TRUE,
      message = FALSE,
      lambda_rr = if (identical(mdl, "Ridge_Regression")) 0.1 else NULL,
      ncomp = if (identical(mdl, "PartialLeastSquare")) 2L else 3L,
      svm_kernel = if (identical(mdl, "SupportVectorMachine")) "Linear" else "Gaussian",
      C_value = if (identical(mdl, "SupportVectorMachine")) 1 else 1,
      iteration = if (identical(mdl, "Xgboost")) 25L else 100L,
      learning_rate = if (identical(mdl, "Xgboost")) 0.05 else 0.01,
      max_depth = if (identical(mdl, "Xgboost")) 3L else 6L,
      subsample = if (identical(mdl, "Xgboost")) 0.8 else 0.7,
      xgb_booster = if (identical(mdl, "Xgboost")) "gbtree" else "dart",
      colsample_bytree = if (identical(mdl, "Xgboost")) 0.8 else 0.7,
      catboost_iterations = if (identical(mdl, "CatBoost")) 50L else 500L,
      catboost_depth = if (identical(mdl, "CatBoost")) 4L else 6L,
      catboost_learning_rate = if (identical(mdl, "CatBoost")) 0.05 else 0.03,
      internal_cv_nfolds = 3L,
      internal_cv_replication = 2L
    )

    preds <- out$model_results$predicted_values
    wide <- out$model_results$multitrait_prediction_wide

    expect_identical(
      names(preds),
      PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
    )
    expect_equal(sort(unique(preds$Trait)), c("Trait1", "Trait2", "Trait3"))
    expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
    expect_equal(nrow(wide), nrow(dat$pheno))
    expect_true(any(vapply(
      out$model_results$diagnostic_plots,
      inherits,
      logical(1L),
      what = "ggplot"
    )))
    expect_true(all(is.finite(preds$Standard_error)))
    expect_equal(preds$PEV, preds$Standard_error^2)
    expect_true(all(grepl("not genetic reliability", preds$Reliability_basis, fixed = TRUE)))
    expect_true(is.matrix(out$model_results$Prediction_covariance_traits))
    expect_true(is.matrix(out$model_results$Prediction_error_covariance_traits))
    expect_true(is.matrix(out$model_results$Prediction_covariance_traits_SE))
    expect_false("Genetic_covariance_traits" %in% names(out$model_results))
  }
})

test_that("multi-trait ML rejects non-gaussian requests", {
  dat <- mk_multitrait_ml_runtime_data()

  expect_error(
    suppressWarnings(PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      GS_model = "RandomForest",
      response_family = "binary",
      multi_trait_ml = TRUE,
      system_database = TRUE,
      message = FALSE
    )),
    "support gaussian traits only"
  )
})

test_that("multi-trait gaussian ML CV returns held-out fold outputs", {
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "RandomForest",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    ntree = 100
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(c("multitrait_metric_summary", "multitrait_prediction_counts", "multitrait_cv_plots", "multitrait_uncertainty_summaries", "multitrait_rank_risk_summaries", "multitrait_rank_risk_plot") %in% names(proc)))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
  expect_true(all(c("combined_uncertainty", "aggregated_uncertainty") %in% names(proc$multitrait_uncertainty_summaries)))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error", "relative_rmse") %in% names(proc$multitrait_uncertainty_summaries$aggregated_uncertainty)))
  expect_true(all(c("combined_risk", "aggregated_risk") %in% names(proc$multitrait_rank_risk_summaries)))
  expect_true(all(c("mean_rank_shift", "mean_rank_instability_risk", "high_risk_percentage") %in% names(proc$multitrait_rank_risk_summaries$aggregated_risk)))
  expect_true(inherits(proc$multitrait_rank_risk_plot, "ggplot"))
})

test_that("multi-trait gaussian Ridge CV returns held-out fold outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "Ridge_Regression",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    lambda_rr = 0.1
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
})

test_that("multi-trait gaussian Lasso CV returns held-out fold outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "Lasso",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
})

test_that("multi-trait gaussian SVM CV returns held-out fold outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "SupportVectorMachine",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    svm_kernel = "Linear",
    C_value = 1
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
})

test_that("multi-trait gaussian Xgboost CV returns held-out fold outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "Xgboost",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    iteration = 25L,
    learning_rate = 0.05,
    max_depth = 3L,
    subsample = 0.8,
    xgb_booster = "gbtree",
    colsample_bytree = 0.8
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
})

test_that("multi-trait gaussian CatBoost CV returns held-out fold outputs", {
  configure_test_python_runtime()
  dat <- mk_multitrait_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "CatBoost",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 50L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.05
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed

  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c("GID", "Trait", "Observed_value", "Predicted_value", "fold", "rep", "cv_role") %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_true(all(sort(unique(proc$multitrait_metric_summary$trait)) == c("Trait1", "Trait2", "Trait3")))
  expect_true(all(c("mean_absolute_error", "root_mean_squared_error") %in% proc$multitrait_metric_summary$metric))
  expect_true(inherits(proc$multitrait_cv_plots, "ggplot"))
})

test_that("multi-trait ML exports true-prediction and CV artifacts", {
  skip_if_not_installed("withr")

  dat <- mk_multitrait_ml_runtime_data()
  tmp_dir <- tempfile("predictpror-mt-ml-export-")
  dir.create(tmp_dir)

  withr::with_dir(tmp_dir, {
    out_true <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = c("Trait1", "Trait2", "Trait3"),
      gen_name = "GID",
      GS_model = "RandomForest",
      response_family = "gaussian",
      multi_trait_ml = TRUE,
      system_database = FALSE,
      message = FALSE,
      ntree = 100,
      internal_cv_nfolds = 3L,
      internal_cv_replication = 2L
    )
    expect_identical(out_true$export_status, "Successful")

    out_cv <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = c("Trait1", "Trait2", "Trait3"),
      gen_name = "GID",
      GS_model_cv = "RandomForest",
      response_family = "gaussian",
      multi_trait_ml = TRUE,
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      cross_validation_meth = "K-Folds",
      nfolds = 3L,
      replication = 1L,
      eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
      system_database = FALSE,
      message = FALSE,
      ntree = 100
    )
    expect_identical(out_cv$export_status, "Successful")

    created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
    true_dir <- created_dirs[grepl("^multi_trait_ml_", basename(created_dirs))]
    cv_dir <- created_dirs[grepl("^CV_results_multi_trait_ml_", basename(created_dirs))]
    expect_length(true_dir, 1)
    expect_length(cv_dir, 1)

    pred_csv <- utils::read.csv(file.path(true_dir, "Predicted_Value.csv"), stringsAsFactors = FALSE)
    expect_identical(
      names(pred_csv),
      PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
    )

    expect_true(file.exists(file.path(cv_dir, "cv_results_processed_multitrait_metric_summary.csv")))
    expect_true(file.exists(file.path(cv_dir, "cv_results_processed_multitrait_uncertainty_summaries_aggregated_uncertainty.csv")))
    expect_true(file.exists(file.path(cv_dir, "cv_results_processed_multitrait_rank_risk_summaries_aggregated_risk.csv")))
    expect_true(file.exists(file.path(cv_dir, "CV_multitrait_observed_vs_predicted.pdf")))
    expect_true(file.exists(file.path(cv_dir, "CV_multitrait_mean_rank_instability_risk.pdf")))
  })
})
