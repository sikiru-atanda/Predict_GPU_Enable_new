if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

configure_met_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")
  invisible(py)
}

met_result_block <- function(out) {
  if (is.list(out) && "model_results" %in% names(out)) {
    return(out)
  }
  if (is.list(out) && length(out) > 0L && is.list(out[[1L]])) {
    return(out[[1L]])
  }
  out
}

met_model_results <- function(out) {
  met_result_block(out)$model_results
}

met_summary_statistic <- function(out) {
  met_result_block(out)$summary_statistic
}

met_test_diagnostic_plots <- function(out) {
  block <- met_result_block(out)
  plots <- block$test_diagonistic_plots
  if (is.null(plots) && is.list(block$model_results)) {
    plots <- block$model_results$diagnostic_plots
  }
  plots
}

expect_met_diagnostic_plots <- function(plots) {
  has_plot <- inherits(plots, "gtable") || inherits(plots, "ggplot") ||
    (is.list(plots) && any(vapply(plots, function(x) !is.null(x), logical(1L))))
  expect_true(has_plot)
}

met_probability_columns <- function(pred) {
  grep("^(Prob_|Probability_)", names(pred), value = TRUE)
}

mk_met_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  set.seed(20260424 + n_gid + length(envs))
  gids <- paste0("g", seq_len(n_gid))
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  geno <- matrix(rnorm(n_gid * 12L), nrow = n_gid, dimnames = list(gids, paste0("m", seq_len(12L))))
  omic1 <- matrix(rnorm(n_gid * 10L), nrow = n_gid, dimnames = list(gids, paste0("o", seq_len(10L))))

  g_eff <- setNames(rnorm(n_gid, 0, 0.45), gids)
  o_eff <- setNames(scale(rowMeans(omic1))[, 1], gids)
  e_eff <- setNames(seq(-0.25, 0.35, length.out = length(envs)), envs)
  pheno$Yield <- 5 + g_eff[pheno$GID] + 0.2 * o_eff[pheno$GID] + e_eff[pheno$Env] + rnorm(nrow(pheno), 0, 0.18)

  grm <- tcrossprod(scale(geno, center = TRUE, scale = FALSE))
  grm <- grm / mean(diag(grm))
  cop <- tcrossprod(scale(omic1, center = TRUE, scale = FALSE))
  cop <- cop / mean(diag(cop))

  list(
    pheno = pheno,
    grm = grm,
    cop = cop
  )
}

mk_met_binary_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_runtime_data(n_gid = n_gid, envs = envs)
  cutoff <- stats::median(dat$pheno$Yield, na.rm = TRUE)
  dat$pheno$Trait <- factor(ifelse(dat$pheno$Yield > cutoff, "yes", "no"), levels = c("no", "yes"))
  dat
}

mk_met_multiclass_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_runtime_data(n_gid = n_gid, envs = envs)
  gids <- unique(dat$pheno$GID)
  gid_class <- setNames(rep(c("low", "mid", "high"), length.out = length(gids)), gids)
  dat$pheno$Trait3 <- factor(gid_class[dat$pheno$GID], levels = c("low", "mid", "high"))
  dat
}

met_model_runtime_args <- function(model) {
  switch(
    model,
    CatBoost = list(
      catboost_iterations = 30L,
      catboost_depth = 4L,
      catboost_learning_rate = 0.1,
      catboost_l2_leaf_reg = 3L,
      catboost_thread_count = 1L,
      n_bootstrap = 3L,
      internal_cv_nfolds = 2L
    ),
    LightGBM = list(
      lightgbm_nrounds = 30L,
      lightgbm_learning_rate = 0.1,
      lightgbm_num_leaves = 15L,
      lightgbm_feature_fraction = 0.8,
      lightgbm_bagging_fraction = 0.8,
      lightgbm_min_data_in_leaf = 3L,
      lightgbm_lambda_l1 = 0,
      lightgbm_lambda_l2 = 0,
      lightgbm_nthread = 1L,
      n_bootstrap = 3L,
      internal_cv_nfolds = 2L
    ),
    Xgboost = list(
      iteration = 30L,
      learning_rate = 0.1,
      max_depth = 4L,
      subsample = 0.8,
      xgb_gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      xgb_alpha = 0,
      xgb_lambda = 1,
      xgb_booster = "gbtree",
      xgb_nthread = 1L,
      n_bootstrap = 3L,
      internal_cv_nfolds = 2L
    ),
    RandomForest = list(
      ntree = 100L,
      mtry = 3L,
      maxnodes = 20L,
      nodesize = 1L,
      n_bootstrap = 3L,
      internal_cv_nfolds = 2L
    ),
    stop("Unsupported MET runtime model in test: ", model, call. = FALSE)
  )
}

test_that("public MET CatBoost true prediction supports multi-kernel per-environment output", {
  configure_met_test_python_runtime()

  dat <- mk_met_runtime_data(n_gid = 6L)
  pheno <- dat$pheno
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Yield[hold_rows] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "CatBoost",
    response_family = "gaussian",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L,
    n_bootstrap = 3L,
    internal_cv_nfolds = 2L
  )

  expect_true(is.list(out))
  expect_true(length(out) >= 1L)
  expect_true("run_metadata" %in% names(met_result_block(out)))

  model_results <- met_model_results(out)
  pred <- model_results$predicted_values
  expect_true(is.data.frame(pred))
  expect_equal(nrow(pred), nrow(pheno))
  expect_true(all(c(
    "GID", "Env", "Predicted_value", "Train_Test_Label",
    "Standard_error", "PEV", "Reliability"
  ) %in% names(pred)))
  expect_true(all(is.finite(pred$Standard_error)))
  expect_true(all(is.finite(pred$PEV)))
  expect_equal(pred$PEV, pred$Standard_error^2)
  expect_gt(length(unique(round(pred$Standard_error, 10))), 1L)
  expect_true(all(is.finite(pred$Reliability)))
  expect_true(all(is.finite(pred$Observed_value[pred$Train_Test_Label == "Train"])))
  expect_true(all(is.na(pred$Observed_value[pred$Train_Test_Label == "Test"])))
  expect_true(all(is.finite(pred$Prediction_stability)))
  expect_true(all(grepl("cross_fitted", pred$Prediction_uncertainty_source)))
  expect_true(all(grepl("not genetic reliability", pred$Reliability_basis, fixed = TRUE)))
  expect_setequal(unique(pred$Train_Test_Label), c("Train", "Test"))
  expect_equal(
    sum(pred$Train_Test_Label == "Test"),
    sum(hold_rows)
  )
  total_pred <- model_results$across_environment_predicted_values
  if (is.null(total_pred)) {
    total_pred <- model_results$Total_Predicted_value
  }
  expect_true(is.data.frame(total_pred))
  expect_equal(nrow(total_pred), length(unique(pheno$GID)))
  expect_true(all(c(
    "GID", "Predicted_value", "Train_Test_Label",
    "Standard_error", "PEV", "Reliability"
  ) %in% names(total_pred)))
  expect_true(all(is.finite(total_pred$Standard_error)))
  expect_true(all(is.finite(total_pred$PEV)))
  expect_equal(total_pred$PEV, total_pred$Standard_error^2)
  expect_true(all(is.finite(total_pred$Reliability)))
  expect_true(all(total_pred$Prediction_uncertainty_source ==
    "across_environment_independence_approximation"))
  expect_true(all(grepl("not genetic reliability", total_pred$Reliability_basis, fixed = TRUE)))
  expect_true(all(c(
    "summary_statistics",
    "MET_pooled_metrics",
    "MET_metrics_by_environment",
    "MET_prediction_counts_by_environment"
  ) %in% names(met_summary_statistic(out))))
  expect_true("predicted_values" %in% names(model_results))
  expect_true(is.data.frame(model_results$predicted_values))
  expect_met_diagnostic_plots(met_test_diagnostic_plots(out))
})

test_that("public MET python tabular models support true prediction with multi-kernel input", {
  configure_met_test_python_runtime()

  dat <- mk_met_runtime_data(n_gid = 6L)
  pheno <- dat$pheno
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Yield[hold_rows] <- NA_real_

  for (model in c("LightGBM", "Xgboost", "RandomForest")) {
    args <- c(list(
      pheno_data = pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = model,
      response_family = "gaussian",
      met_ml_dl = TRUE,
      system_database = TRUE,
      message = FALSE
    ), met_model_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)

    expect_true(is.list(out))
    pred <- met_model_results(out)$predicted_values
    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), nrow(pheno))
    expect_true(all(c("GID", "Env", "Predicted_value", "Train_Test_Label") %in% names(pred)))
    expect_setequal(unique(pred$Train_Test_Label), c("Train", "Test"))
  }
})

test_that("public MET CatBoost export writes MET summaries and plot files", {
  configure_met_test_python_runtime()

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-met-catboost-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  dat <- mk_met_runtime_data(n_gid = 6L)
  pheno <- dat$pheno
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Yield[hold_rows] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "CatBoost",
    response_family = "gaussian",
    met_ml_dl = TRUE,
    system_database = FALSE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L
  )

  expect_identical(met_result_block(out)$export_status, "Successful")
  created_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_length(created_dirs, 1L)

  expect_true(file.exists(file.path(created_dirs[[1]], "MET_pooled_metrics.csv")))
  expect_true(file.exists(file.path(created_dirs[[1]], "MET_metrics_by_environment.csv")))
  expect_true(file.exists(file.path(created_dirs[[1]], "MET_prediction_counts_by_environment.csv")))
  expect_true(file.exists(file.path(created_dirs[[1]], "Predicted_Value.csv")))
  expect_true(file.exists(file.path(created_dirs[[1]], "across_environment_predicted_values.csv")))
  expect_true(any(grepl("^test_diagonistic_plots_GS_model_.*\\.pdf$", list.files(created_dirs[[1]]))))
})

test_that("public MET CatBoost binary true prediction returns probabilities and confidence outputs", {
  configure_met_test_python_runtime()

  dat <- mk_met_binary_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait[hold_rows] <- NA

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "CatBoost",
    response_family = "binary",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L
  )

  pred <- met_model_results(out)$predicted_values
  expect_true(all(c(
    "Predicted_class",
    "Prediction_confidence",
    "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(met_probability_columns(pred)), 2L)
  expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
  expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(met_summary_statistic(out))))
  expect_met_diagnostic_plots(met_test_diagnostic_plots(out))
})

test_that("public MET python tabular models support binary true prediction on the shared path", {
  configure_met_test_python_runtime()

  dat <- mk_met_binary_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait[hold_rows] <- NA

  for (model in c("LightGBM", "Xgboost", "RandomForest")) {
    args <- c(list(
      pheno_data = pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Trait",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = model,
      response_family = "binary",
      met_ml_dl = TRUE,
      system_database = TRUE,
      message = FALSE
    ), met_model_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)
    pred <- met_model_results(out)$predicted_values

    expect_true(all(c(
      "Predicted_class",
      "Prediction_confidence",
      "Classification_uncertainty"
    ) %in% names(pred)))
    expect_gte(length(met_probability_columns(pred)), 2L)
    expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
    expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(met_summary_statistic(out))))
    expect_met_diagnostic_plots(met_test_diagnostic_plots(out))
  }
})

test_that("public MET CatBoost CV runtime covers CV0 CV1 and CV2", {
  configure_met_test_python_runtime()

  dat <- mk_met_runtime_data(n_gid = 8L)

  check_cv_case <- function(cv_method) {
    out <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model_cv = "CatBoost",
      response_family = "gaussian",
      met_ml_dl = TRUE,
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      cross_validation_meth = cv_method,
      nfolds = 2L,
      replication = 1L,
      eval_metrics = c("root_mean_squared_error", "kendalls_tau"),
      system_database = TRUE,
      message = FALSE,
      catboost_iterations = 30L,
      catboost_depth = 4L,
      catboost_learning_rate = 0.1,
      catboost_l2_leaf_reg = 3L,
      catboost_thread_count = 1L
    )

    expect_true(is.list(out))
    expect_true("cv_results_raw" %in% names(out))
    expect_true("cv_results_processed" %in% names(out))
    expect_length(out$cv_results_raw, 1L)
    expect_true(is.list(out$cv_results_processed))
    expect_true("run_metadata" %in% names(out$cv_results_processed))

    raw <- out$cv_results_raw[[1L]]
    expect_true(is.data.frame(raw$ypred_cv_Reps_all))
    expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
    expect_true(is.data.frame(raw$eval_metrics_reps))
    expect_gt(nrow(raw$eval_metrics_reps), 0L)
    expect_true(any(is.finite(raw$ypred_cv_Reps_all$yhat)))
    expect_true(all(c("root_mean_squared_error", "kendalls_tau") %in% names(raw$eval_metrics_reps)))
    expect_true(any(is.finite(raw$eval_metrics_reps$root_mean_squared_error)))
    expect_true(any(is.finite(raw$eval_metrics_reps$kendalls_tau)))

    proc <- out$cv_results_processed
    expect_true(all(c(
      "aggregated_data_list",
      "plot_reps_list",
      "plot_mean_list",
      "best_models_list"
    ) %in% names(proc)))
  }

  lapply(c("CV0", "CV1", "CV2"), check_cv_case)
})

test_that("public MET python tabular models support CV runtime on the shared path", {
  configure_met_test_python_runtime()

  dat <- mk_met_runtime_data(n_gid = 8L)

  for (model in c("LightGBM", "Xgboost", "RandomForest")) {
    args <- c(list(
      pheno_data = dat$pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model_cv = model,
      response_family = "gaussian",
      met_ml_dl = TRUE,
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      cross_validation_meth = "CV2",
      nfolds = 2L,
      replication = 1L,
      eval_metrics = c("root_mean_squared_error", "kendalls_tau"),
      system_database = TRUE,
      message = FALSE
    ), met_model_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)

    expect_true(is.list(out))
    expect_true("cv_results_raw" %in% names(out))
    expect_length(out$cv_results_raw, 1L)
    raw <- out$cv_results_raw[[1L]]
    expect_true(is.data.frame(raw$ypred_cv_Reps_all))
    expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
    expect_true(any(is.finite(raw$ypred_cv_Reps_all$yhat)))
    expect_true(is.data.frame(raw$eval_metrics_reps))
    expect_true(any(is.finite(raw$eval_metrics_reps$root_mean_squared_error)))
  }
})

test_that("public MET CatBoost binary CV runtime exposes probabilities and classification summaries", {
  configure_met_test_python_runtime()

  dat <- mk_met_binary_runtime_data(n_gid = 8L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = "CatBoost",
    response_family = "binary",
    met_ml_dl = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = c("accuracy", "log_loss", "brier_score", "ece"),
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L
  )

  expect_true(is.list(out))
  expect_true("cv_results_raw" %in% names(out))
  expect_true("cv_results_processed" %in% names(out))
  expect_length(out$cv_results_raw, 1L)

  raw <- out$cv_results_raw[[1L]]
  expect_true(is.data.frame(raw$ypred_cv_Reps_all))
  expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
  expect_true(is.data.frame(raw$yprob_cv_Reps_all))
  prob_cols <- setdiff(names(raw$yprob_cv_Reps_all), c("y", "yhat", "row_id", "rep", "trait", "model", "cv_role"))
  expect_gte(length(prob_cols), 2L)
  expect_true(any(is.finite(as.matrix(raw$yprob_cv_Reps_all[, prob_cols, drop = FALSE]))))
  expect_true(is.data.frame(raw$eval_metrics_reps))
  expect_true(any(is.finite(raw$eval_metrics_reps$accuracy)))
  expect_true(any(is.finite(raw$eval_metrics_reps$log_loss)))

  proc <- out$cv_results_processed
  expect_true("classification_probability_summaries" %in% names(proc))
  probs <- proc$classification_probability_summaries$aggregated_probabilities
  expect_true(is.data.frame(probs))
  expect_true(all(c("trait", "model", "row_id") %in% names(probs)))
  agg_prob_cols <- proc$classification_probability_summaries$probability_columns
  expect_true(all(agg_prob_cols %in% names(probs)))
  expect_gte(length(agg_prob_cols), 2L)
  expect_true(any(is.finite(as.matrix(probs[, agg_prob_cols, drop = FALSE]))))
  expect_true(!is.null(out$res_plot_result_diagnostic) || !is.null(out$res_mod_results_cv_per_trait_model))
})

test_that("public MET CatBoost multiclass true prediction returns probabilities and confidence outputs", {
  configure_met_test_python_runtime()

  dat <- mk_met_multiclass_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait3[hold_rows] <- NA

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait3",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "CatBoost",
    response_family = "multiclass",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L
  )

  pred <- met_model_results(out)$predicted_values
  expect_true(all(c(
    "Predicted_class",
    "Prediction_confidence",
    "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(met_probability_columns(pred)), 3L)
  expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
  expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(met_summary_statistic(out))))
  expect_met_diagnostic_plots(met_test_diagnostic_plots(out))
})

test_that("public MET python tabular models support multiclass true prediction on the shared path", {
  configure_met_test_python_runtime()

  dat <- mk_met_multiclass_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait3[hold_rows] <- NA

  for (model in c("LightGBM", "Xgboost", "RandomForest")) {
    args <- c(list(
      pheno_data = pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Trait3",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = model,
      response_family = "multiclass",
      met_ml_dl = TRUE,
      system_database = TRUE,
      message = FALSE
    ), met_model_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)
    pred <- met_model_results(out)$predicted_values

    expect_true(all(c(
      "Predicted_class",
      "Prediction_confidence",
      "Classification_uncertainty"
    ) %in% names(pred)))
    expect_gte(length(met_probability_columns(pred)), 3L)
    expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
    expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(met_summary_statistic(out))))
    expect_met_diagnostic_plots(met_test_diagnostic_plots(out))
  }
})

test_that("public MET CatBoost multiclass CV runtime exposes probabilities and classification summaries", {
  configure_met_test_python_runtime()

  dat <- mk_met_multiclass_runtime_data(n_gid = 8L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait3",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = "CatBoost",
    response_family = "multiclass",
    met_ml_dl = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = c("accuracy", "log_loss", "macro_f1", "ece"),
    system_database = TRUE,
    message = FALSE,
    catboost_iterations = 30L,
    catboost_depth = 4L,
    catboost_learning_rate = 0.1,
    catboost_l2_leaf_reg = 3L,
    catboost_thread_count = 1L
  )

  raw <- out$cv_results_raw[[1L]]
  expect_true(is.data.frame(raw$ypred_cv_Reps_all))
  expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
  expect_true(is.data.frame(raw$yprob_cv_Reps_all))
  prob_cols <- setdiff(names(raw$yprob_cv_Reps_all), c("y", "yhat", "row_id", "rep", "trait", "model", "cv_role"))
  expect_gte(length(prob_cols), 3L)
  expect_true(any(is.finite(as.matrix(raw$yprob_cv_Reps_all[, prob_cols, drop = FALSE]))))

  proc <- out$cv_results_processed
  expect_true("classification_probability_summaries" %in% names(proc))
  probs <- proc$classification_probability_summaries$aggregated_probabilities
  expect_true(is.data.frame(probs))
  agg_prob_cols <- proc$classification_probability_summaries$probability_columns
  expect_true(all(agg_prob_cols %in% names(probs)))
  expect_gte(length(agg_prob_cols), 3L)
  expect_true(any(is.finite(as.matrix(probs[, agg_prob_cols, drop = FALSE]))))
})

test_that("public MET python tabular models support multiclass CV runtime on the shared path", {
  configure_met_test_python_runtime()

  dat <- mk_met_multiclass_runtime_data(n_gid = 8L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]

  for (model in c("LightGBM", "Xgboost", "RandomForest")) {
    args <- c(list(
      pheno_data = pheno,
      gmatrix = dat$grm,
      omic1_kernel = dat$cop,
      omics_kernel_label = list(omic1_kernel = "COP"),
      response = "Trait3",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model_cv = model,
      response_family = "multiclass",
      met_ml_dl = TRUE,
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      cross_validation_meth = "CV2",
      nfolds = 2L,
      replication = 1L,
      eval_metrics = c("accuracy", "log_loss", "macro_f1", "ece"),
      system_database = TRUE,
      message = FALSE
    ), met_model_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)
    raw <- out$cv_results_raw[[1L]]
    expect_true(is.data.frame(raw$ypred_cv_Reps_all))
    expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
    expect_true(is.data.frame(raw$yprob_cv_Reps_all))
    prob_cols <- setdiff(names(raw$yprob_cv_Reps_all), c("y", "yhat", "row_id", "rep", "trait", "model", "cv_role"))
    expect_gte(length(prob_cols), 3L)
    expect_true(any(is.finite(as.matrix(raw$yprob_cv_Reps_all[, prob_cols, drop = FALSE]))))

    proc <- out$cv_results_processed
    expect_true("classification_probability_summaries" %in% names(proc))
    probs <- proc$classification_probability_summaries$aggregated_probabilities
    expect_true(is.data.frame(probs))
    agg_prob_cols <- proc$classification_probability_summaries$probability_columns
    expect_true(all(agg_prob_cols %in% names(probs)))
    expect_gte(length(agg_prob_cols), 3L)
    expect_true(any(is.finite(as.matrix(probs[, agg_prob_cols, drop = FALSE]))))
  }
})
