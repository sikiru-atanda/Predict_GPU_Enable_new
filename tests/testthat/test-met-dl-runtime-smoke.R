if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

configure_met_tabular_dl_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_PYTHON = py)
  invisible(py)
}

mk_met_tabular_dl_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  set.seed(20260426 + n_gid + length(envs))
  gids <- paste0("g", seq_len(n_gid))
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  geno <- matrix(rnorm(n_gid * 10L), nrow = n_gid, dimnames = list(gids, paste0("m", seq_len(10L))))
  omic1 <- matrix(rnorm(n_gid * 8L), nrow = n_gid, dimnames = list(gids, paste0("o", seq_len(8L))))

  g_eff <- setNames(rnorm(n_gid, 0, 0.35), gids)
  o_eff <- setNames(scale(rowMeans(omic1))[, 1], gids)
  e_eff <- setNames(seq(-0.25, 0.25, length.out = length(envs)), envs)
  pheno$Yield <- 4.7 + g_eff[pheno$GID] + 0.12 * o_eff[pheno$GID] + e_eff[pheno$Env] + rnorm(nrow(pheno), 0, 0.18)

  grm <- tcrossprod(scale(geno, center = TRUE, scale = FALSE))
  grm <- grm / mean(diag(grm))
  cop <- tcrossprod(scale(omic1, center = TRUE, scale = FALSE))
  cop <- cop / mean(diag(cop))

  list(pheno = pheno, grm = grm, cop = cop)
}

mk_met_tabular_dl_binary_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_tabular_dl_runtime_data(n_gid = n_gid, envs = envs)
  cutoff <- stats::median(dat$pheno$Yield, na.rm = TRUE)
  dat$pheno$Trait <- factor(ifelse(dat$pheno$Yield > cutoff, "yes", "no"), levels = c("no", "yes"))
  dat
}

mk_met_tabular_dl_multiclass_runtime_data <- function(n_gid = 8L, envs = c("E1", "E2", "E3")) {
  dat <- mk_met_tabular_dl_runtime_data(n_gid = n_gid, envs = envs)
  gids <- unique(dat$pheno$GID)
  gid_class <- setNames(rep(c("low", "mid", "high"), length.out = length(gids)), gids)
  dat$pheno$Trait3 <- factor(gid_class[dat$pheno$GID], levels = c("low", "mid", "high"))
  dat
}

met_tabular_dl_runtime_args <- function(model) {
  base <- list(
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    validation_split = 0.2,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    epochs = 1L,
    batch_size = 4L,
    dropout = 0.1,
    l2_weight_decay = 1e-4,
    n_bootstrap = 2L,
    internal_cv_nfolds = 2L
  )
  switch(
    model,
    ft_transformer = c(base, list(
      ft_d_model = 32L,
      ft_heads = 4L,
      ft_layers = 2L,
      ft_ff_mult = 2L,
      ft_dropout = 0.1,
      ft_token_dropout = 0.1,
      ft_use_cls = TRUE
    )),
    saint = c(base, list(
      saint_d_model = 32L,
      saint_heads = 4L,
      saint_layers = 2L,
      saint_ff_mult = 2L,
      saint_dropout = 0.1,
      saint_token_dropout = 0.1,
      saint_use_cls = TRUE
    )),
    tabnet = c(base, list(
      tabnet_steps = 3L,
      tabnet_feature_dim = 16L,
      tabnet_output_dim = 16L,
      tabnet_gamma = 1.3,
      tabnet_lambda_sparse = 1e-4
    )),
    moe = c(base, list(
      moe_n_experts = 3L,
      moe_expert_hidden = as.integer(c(32L, 16L)),
      moe_gate_hidden = 32L,
      moe_temperature = 1.0,
      moe_sparse_topk = NA,
      moe_entropy_reg = 0.0
    )),
    stop("Unsupported test model: ", model, call. = FALSE)
  )
}

expect_met_dl_classification_contract <- function(pred, class_levels) {
  probability_columns <- paste0(
    "Probability_",
    make.names(class_levels, unique = TRUE)
  )
  expect_identical(
    names(pred),
    c(
      PredictProR:::gp_public_prediction_columns(
        "multiclass",
        include_env = TRUE
      ),
      probability_columns
    )
  )
  expect_true(all(pred$Train_Test_Label %in% c("Train", "Test")))
  expect_true(all(is.finite(pred$Prediction_confidence)))
  expect_true(all(is.finite(pred$Classification_uncertainty)))
  expect_true(all(is.finite(pred$Reliability)))
  probability <- as.matrix(pred[, probability_columns, drop = FALSE])
  expect_true(all(is.finite(probability)))
  expect_true(all(probability >= 0 & probability <= 1))
  expect_equal(rowSums(probability), rep(1, nrow(pred)), tolerance = 1e-6)
  confidence <- apply(probability, 1L, max)
  expect_equal(pred$Prediction_confidence, confidence, tolerance = 1e-6)
  expect_equal(pred$Classification_uncertainty, 1 - confidence, tolerance = 1e-6)
  expect_equal(pred$Reliability, confidence, tolerance = 1e-6)
  expect_identical(
    pred$Predicted_class,
    class_levels[max.col(probability, ties.method = "first")]
  )
}

test_that("public MET tabular DL models support true prediction with multi-kernel input", {
  configure_met_tabular_dl_test_python_runtime()

  dat <- mk_met_tabular_dl_runtime_data(n_gid = 6L)
  pheno <- dat$pheno
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Yield[hold_rows] <- NA_real_

  for (model in c("ft_transformer", "saint", "tabnet", "moe")) {
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
    ), met_tabular_dl_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)

    expect_true(is.list(out))
    pred <- out$model_results$predicted_values
    expect_true(is.data.frame(pred))
    expect_equal(nrow(pred), nrow(pheno))
    expect_true(all(c("GID", "Env", "Predicted_value", "Train_Test_Label") %in% names(pred)))
    expect_setequal(unique(pred$Train_Test_Label), c("Train", "Test"))
    expect_true(all(c(
      "Standard_error", "PEV", "Prediction_stability", "Reliability",
      "Prediction_uncertainty_source", "Reliability_basis"
    ) %in% names(pred)))
    expect_true(all(is.finite(pred$Standard_error)))
    expect_true(all(is.finite(pred$PEV)))
    expect_equal(pred$PEV, pred$Standard_error^2)
    expect_gt(length(unique(round(pred$Standard_error, 10))), 1L)
    expect_true(all(is.finite(pred$Prediction_stability)))
    expect_true(all(is.finite(pred$Reliability)))
    expect_true(all(grepl("cross_fitted", pred$Prediction_uncertainty_source)))
    expect_true(all(grepl("not genetic reliability", pred$Reliability_basis, fixed = TRUE)))
    expect_true(all(c(
      "summary_statistics",
      "MET_pooled_metrics",
      "MET_metrics_by_environment",
      "MET_prediction_counts_by_environment"
    ) %in% names(out$summary_statistic)))
    expect_true("predicted_values" %in% names(out$model_results))
    expect_true(any(vapply(
      out$model_results$diagnostic_plots,
      function(plot) inherits(plot, c("ggplot", "gtable")),
      logical(1L)
    )))
  }
})

test_that("public MET tabular DL models support CV2 runtime", {
  configure_met_tabular_dl_test_python_runtime()

  dat <- mk_met_tabular_dl_runtime_data(n_gid = 8L)

  for (model in c("ft_transformer", "saint", "tabnet", "moe")) {
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
    ), met_tabular_dl_runtime_args(model))

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

test_that("public MET tabular DL models support binary true prediction", {
  configure_met_tabular_dl_test_python_runtime()

  dat <- mk_met_tabular_dl_binary_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait[hold_rows] <- NA

  for (model in c("ft_transformer", "saint", "tabnet", "moe")) {
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
    ), met_tabular_dl_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)
    pred <- out$model_results$predicted_values

    expect_met_dl_classification_contract(pred, c("no", "yes"))
    expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
    expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(out$summary_statistic)))
    expect_true(any(vapply(
      out$model_results$diagnostic_plots,
      function(plot) inherits(plot, c("ggplot", "gtable")),
      logical(1L)
    )))
  }
})

test_that("public MET tabular DL models support multiclass true prediction", {
  configure_met_tabular_dl_test_python_runtime()

  dat <- mk_met_tabular_dl_multiclass_runtime_data(n_gid = 6L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]
  hold_rows <- pheno$GID %in% c("g2", "g5") & pheno$Env %in% c("E2", "E3")
  pheno$Trait3[hold_rows] <- NA

  for (model in c("ft_transformer", "saint", "tabnet", "moe")) {
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
    ), met_tabular_dl_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)
    pred <- out$model_results$predicted_values

    expect_met_dl_classification_contract(pred, c("low", "mid", "high"))
    expect_false(any(c("Standard_error", "PEV", "Prediction_stability") %in% names(pred)))
    expect_true(all(c("MET_pooled_metrics", "MET_metrics_by_environment") %in% names(out$summary_statistic)))
    expect_true(any(vapply(
      out$model_results$diagnostic_plots,
      function(plot) inherits(plot, c("ggplot", "gtable")),
      logical(1L)
    )))
  }
})

test_that("public MET tabular DL models support multiclass CV2 runtime", {
  configure_met_tabular_dl_test_python_runtime()

  dat <- mk_met_tabular_dl_multiclass_runtime_data(n_gid = 8L)
  pheno <- dat$pheno[, c("GID", "Env", "Trait3"), drop = FALSE]

  for (model in c("ft_transformer", "saint", "tabnet", "moe")) {
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
    ), met_tabular_dl_runtime_args(model))

    out <- do.call(PredictProR::model_execute, args)

    expect_true(is.list(out))
    expect_true("cv_results_raw" %in% names(out))
    expect_length(out$cv_results_raw, 1L)
    raw <- out$cv_results_raw[[1L]]
    expect_true(is.data.frame(raw$ypred_cv_Reps_all))
    expect_gt(nrow(raw$ypred_cv_Reps_all), 0L)
    expect_true(is.data.frame(raw$yprob_cv_Reps_all))
    prob_cols <- setdiff(names(raw$yprob_cv_Reps_all), c("y", "yhat", "row_id", "rep", "trait", "model", "cv_role"))
    expect_gte(length(prob_cols), 3L)
    expect_true(any(is.finite(as.matrix(raw$yprob_cv_Reps_all[, prob_cols, drop = FALSE]))))
    expect_true(is.data.frame(raw$eval_metrics_reps))
    expect_true(any(is.finite(raw$eval_metrics_reps$accuracy)))
    expect_true(any(is.finite(raw$eval_metrics_reps$log_loss)))

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
