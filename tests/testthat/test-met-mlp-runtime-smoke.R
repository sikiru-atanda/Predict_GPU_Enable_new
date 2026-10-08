if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

configure_met_dl_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_PYTHON = py)
  invisible(py)
}

mk_met_dl_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  set.seed(20260425 + n_gid + length(envs))
  gids <- paste0("g", seq_len(n_gid))
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  geno <- matrix(rnorm(n_gid * 10L), nrow = n_gid, dimnames = list(gids, paste0("m", seq_len(10L))))
  omic1 <- matrix(rnorm(n_gid * 8L), nrow = n_gid, dimnames = list(gids, paste0("o", seq_len(8L))))

  g_eff <- setNames(rnorm(n_gid, 0, 0.4), gids)
  o_eff <- setNames(scale(rowMeans(omic1))[, 1], gids)
  e_eff <- setNames(seq(-0.3, 0.3, length.out = length(envs)), envs)
  pheno$Yield <- 4.5 + g_eff[pheno$GID] + 0.15 * o_eff[pheno$GID] + e_eff[pheno$Env] + rnorm(nrow(pheno), 0, 0.2)

  grm <- tcrossprod(scale(geno, center = TRUE, scale = FALSE))
  grm <- grm / mean(diag(grm))
  cop <- tcrossprod(scale(omic1, center = TRUE, scale = FALSE))
  cop <- cop / mean(diag(cop))

  list(pheno = pheno, grm = grm, cop = cop)
}

mk_met_dl_binary_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_dl_runtime_data(n_gid = n_gid, envs = envs)
  cutoff <- stats::median(dat$pheno$Yield, na.rm = TRUE)
  dat$pheno$Trait <- factor(ifelse(dat$pheno$Yield > cutoff, "yes", "no"), levels = c("no", "yes"))
  dat
}

mk_met_dl_multiclass_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_dl_runtime_data(n_gid = n_gid, envs = envs)
  gids <- unique(dat$pheno$GID)
  gid_class <- setNames(rep(c("low", "mid", "high"), length.out = length(gids)), gids)
  dat$pheno$Trait3 <- factor(gid_class[dat$pheno$GID], levels = c("low", "mid", "high"))
  dat
}

met_mlp_runtime_args <- list(
  compile_model = FALSE,
  deterministic = TRUE,
  random_seed = 1L,
  validation_split = 0.2,
  device = "cpu",
  use_amp = FALSE,
  batch_norm = FALSE,
  epochs = 1L,
  batch_size = 4L,
  mlp_neurons_per_layer = as.integer(c(16L, 8L)),
  mlp_learning_rate = 1e-3,
  dropout = 0.1,
  n_bootstrap = 2L
)

test_that("public MET mlp true prediction supports multi-kernel per-environment output", {
  configure_met_dl_test_python_runtime()

  dat <- mk_met_dl_runtime_data(n_gid = 6L)
  pheno <- dat$pheno
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Yield[hold_rows] <- NA_real_

  args <- c(list(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "mlp",
    response_family = "gaussian",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE
  ), met_mlp_runtime_args)

  out <- do.call(PredictProR::model_execute, args)

  expect_true(is.list(out))
  pred <- out$model_results$predicted_values
  expect_true(is.data.frame(pred))
  expect_equal(nrow(pred), nrow(pheno))
  expect_true(all(c("GID", "Env", "Predicted_value", "Train_Test_Label") %in% names(pred)))
  expect_setequal(unique(pred$Train_Test_Label), c("Train", "Test"))
  expect_true(all(c(
    "summary_statistics",
    "MET_pooled_metrics",
    "MET_metrics_by_environment",
    "MET_prediction_counts_by_environment"
  ) %in% names(out$summary_statistic)))
  expect_true(is.matrix(out$model_results$Prediction_covariance_environments))
  expect_true(is.matrix(out$model_results$Prediction_error_covariance_environments))
  expect_true(is.matrix(out$model_results$Prediction_covariance_environments_SE))
  expect_false("Genetic_covariance_environments" %in% names(out$model_results))
  expect_true(all(grepl("not genetic reliability", pred$Reliability_basis, fixed = TRUE)))
  expect_true(any(vapply(
    out$model_results$diagnostic_plots,
    inherits,
    logical(1L),
    what = "ggplot"
  )))
})

test_that("public MET mlp CV2 runtime uses the shared long-format kernel route", {
  configure_met_dl_test_python_runtime()

  dat <- mk_met_dl_runtime_data(n_gid = 8L)

  args <- c(list(
    pheno_data = dat$pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = "mlp",
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
  ), met_mlp_runtime_args)

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
})

test_that("public MET mlp binary true prediction returns probabilities and confidence outputs", {
  configure_met_dl_test_python_runtime()

  dat <- mk_met_dl_binary_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait[hold_rows] <- NA

  args <- c(list(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "mlp",
    response_family = "binary",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE
  ), met_mlp_runtime_args)

  out <- do.call(PredictProR::model_execute, args)

  pred <- out$model_results$predicted_values
  expect_true(all(c(
    "Predicted_class",
    "Prediction_confidence",
    "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(grep("^Probability_", names(pred), value = TRUE)), 2L)
  expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
  expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(out$summary_statistic)))
  expect_true(any(vapply(
    out$model_results$diagnostic_plots,
    function(plot) inherits(plot, c("ggplot", "gtable")),
    logical(1L)
  )))
})

test_that("public MET mlp binary CV runtime exposes probabilities and classification summaries", {
  configure_met_dl_test_python_runtime()

  dat <- mk_met_dl_binary_runtime_data(n_gid = 8L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]

  args <- c(list(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = "mlp",
    response_family = "binary",
    met_ml_dl = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = c("accuracy", "log_loss", "brier_score", "ece"),
    system_database = TRUE,
    message = FALSE
  ), met_mlp_runtime_args)

  out <- do.call(PredictProR::model_execute, args)

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

test_that("public MET mlp multiclass true prediction returns probabilities and confidence outputs", {
  configure_met_dl_test_python_runtime()

  dat <- mk_met_dl_multiclass_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait3[hold_rows] <- NA

  args <- c(list(
    pheno_data = pheno,
    gmatrix = dat$grm,
    omic1_kernel = dat$cop,
    omics_kernel_label = list(omic1_kernel = "COP"),
    response = "Trait3",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "mlp",
    response_family = "multiclass",
    met_ml_dl = TRUE,
    system_database = TRUE,
    message = FALSE
  ), met_mlp_runtime_args)

  out <- do.call(PredictProR::model_execute, args)

  pred <- out$model_results$predicted_values
  expect_true(all(c(
    "Predicted_class",
    "Prediction_confidence",
    "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(grep("^Probability_", names(pred), value = TRUE)), 3L)
  expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
  expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(out$summary_statistic)))
  expect_true(any(vapply(
    out$model_results$diagnostic_plots,
    function(plot) inherits(plot, c("ggplot", "gtable")),
    logical(1L)
  )))
})
