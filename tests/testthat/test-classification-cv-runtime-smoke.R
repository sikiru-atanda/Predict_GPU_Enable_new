if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

configure_cv_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")
  invisible(py)
}

mk_cv_runtime_data <- function(kind = c("binary", "ordinal", "multiclass"), n = 30L, p = 8L) {
  kind <- match.arg(kind)
  set.seed(20260425 + match(kind, c("binary", "ordinal", "multiclass")))
  latent_x <- matrix(rnorm(n * p), nrow = n)
  X <- matrix(
    as.numeric(latent_x > -0.45) + as.numeric(latent_x > 0.45),
    nrow = n, ncol = p
  )
  rownames(X) <- paste0(substr(kind, 1, 1), seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))

  y <- switch(
    kind,
    binary = factor(ifelse(runif(n) < plogis(X[, 1] - X[, 2]), "yes", "no"), levels = c("no", "yes")),
    ordinal = {
      score <- X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3] + rnorm(n, sd = 0.15)
      cuts <- stats::quantile(score, probs = c(1 / 3, 2 / 3), names = FALSE)
      ordered(
        ifelse(score <= cuts[[1L]], "low", ifelse(score <= cuts[[2L]], "medium", "high")),
        levels = c("low", "medium", "high")
      )
    },
    multiclass = {
      sc <- cbind(A = X[, 1], B = X[, 2], C = -X[, 1] - X[, 2]) + matrix(rnorm(n * 3, sd = 0.2), nrow = n)
      factor(colnames(sc)[max.col(sc, ties.method = "first")], levels = c("A", "B", "C"))
    }
  )

  list(
    pheno = data.frame(GID = rownames(X), Trait = y, stringsAsFactors = FALSE),
    geno = X
  )
}

test_that("binary classification CV runtime exposes probability summaries and diagnostics", {
  configure_cv_test_python_runtime()

  dat <- mk_cv_runtime_data("binary", n = 24L, p = 6L)

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "binary",
    ploidy = 2L,
    qc_filtering = FALSE,
    impute = FALSE,
    ld_prunning_qc = FALSE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hold_Out",
    test_size = 0.25,
    random_state = 123L,
    replication = 1L,
    GS_model_cv = c("CatBoost"),
    eval_metrics = c("accuracy", "log_loss", "ece"),
    para_tunning = TRUE,
    catboost_paras_tunning = list(
      catboost_iterations = c(25L),
      catboost_depth = c(4L),
      catboost_learning_rate = c(0.05),
      catboost_l2_leaf_reg = c(3L)
    ),
    system_database = TRUE,
    message = FALSE,
    scaling = FALSE,
    centering = FALSE
  )

  expect_true(is.list(out))
  expect_true("cv_results_processed" %in% names(out))

  proc <- out$cv_results_processed
  expect_true("classification_probability_summaries" %in% names(proc))
  expect_true("run_metadata" %in% names(proc))

  probs <- proc$classification_probability_summaries$aggregated_probabilities
  expect_true(is.data.frame(probs))
  expect_true(all(c("trait", "model", "row_id") %in% names(probs)))
  expect_true(all(c("no", "yes") %in% names(probs)))
  expect_true(any(is.finite(probs$yes)))
  expect_true(is.list(proc$classification_metrics))
  expect_true(is.list(proc$classification_reliability))
  expect_true(all(is.finite(proc$classification_reliability$out_of_fold_predictions$Reliability)))
  expect_true(all(proc$classification_reliability$out_of_fold_predictions$Reliability >= 0 &
                    proc$classification_reliability$out_of_fold_predictions$Reliability <= 1))
  expect_true(is.finite(proc$classification_reliability$summary$ece[[1L]]))
  expect_true(!is.null(out$res_plot_result_diagnostic) || !is.null(out$res_mod_results_cv_per_trait_model))

  raw_probs <- out$cv_results_raw[[1L]]$yprob_cv_Reps_all
  expect_true(is.data.frame(raw_probs))
  expect_true(all(c("no", "yes") %in% names(raw_probs)))
})

test_that("multiclass classification CV runtime exposes probability summaries and diagnostics", {
  configure_cv_test_python_runtime()

  dat <- mk_cv_runtime_data("multiclass", n = 27L, p = 6L)

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "multiclass",
    ploidy = 2L,
    qc_filtering = FALSE,
    impute = FALSE,
    ld_prunning_qc = FALSE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hold_Out",
    test_size = 0.25,
    random_state = 321L,
    replication = 1L,
    GS_model_cv = c("mlp"),
    eval_metrics = c("accuracy", "log_loss", "macro_f1"),
    para_tunning = TRUE,
    dpl_paras_tunning = list(
      epochs = c(1L),
      batch_size = c(4L),
      mlp_learning_rate = c(1e-3)
    ),
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
    system_database = TRUE,
    message = FALSE,
    scaling = FALSE,
    centering = FALSE
  )

  expect_true(is.list(out))
  expect_true("cv_results_processed" %in% names(out))

  proc <- out$cv_results_processed
  probs <- proc$classification_probability_summaries$aggregated_probabilities
  expect_true(is.data.frame(probs))
  expect_true(all(c("trait", "model", "row_id") %in% names(probs)))
  expect_true(all(c("A", "B", "C") %in% names(probs)))
  expect_true(any(is.finite(probs$A)))
  expect_true(is.list(proc$classification_metrics))
  expect_true(is.list(proc$classification_reliability))
  expect_true(all(is.finite(proc$classification_reliability$out_of_fold_predictions$Reliability)))
  expect_true(is.finite(proc$classification_reliability$summary$ece[[1L]]))
  expect_true(!is.null(out$res_plot_result_diagnostic) || !is.null(out$res_mod_results_cv_per_trait_model))

  raw_probs <- out$cv_results_raw[[1L]]$yprob_cv_Reps_all
  expect_true(is.data.frame(raw_probs))
  expect_true(all(c("A", "B", "C") %in% names(raw_probs)))
})

test_that("ordinal Bayesian CV exposes ordinal metrics and calibrated reliability", {
  skip_if_not_installed("BGLR")
  dat <- mk_cv_runtime_data("ordinal", n = 24L, p = 6L)

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Trait",
    gen_name = "GID",
    fixed = ~ 1,
    random = ~ GID,
    response_family = "ordinal",
    ploidy = 2L,
    qc_filtering = FALSE,
    impute = FALSE,
    ld_prunning_qc = FALSE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 2L,
    sampling_method = "stratified",
    random_state = 456L,
    replication = 1L,
    GS_model_cv = "BayesA",
    eval_metrics = c(
      "accuracy", "balanced_accuracy", "mean_absolute_error_class",
      "within_one_class_accuracy", "quadratic_weighted_kappa", "log_loss", "ece"
    ),
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE,
    scaling = FALSE,
    centering = FALSE
  )

  proc <- out$cv_results_processed
  metrics <- proc$classification_metrics$aggregated
  expect_true(all(c(
    "accuracy", "balanced_accuracy", "mean_absolute_error_class",
    "within_one_class_accuracy", "quadratic_weighted_kappa", "log_loss", "ece"
  ) %in% names(metrics)))
  expect_true(all(is.finite(as.numeric(metrics[1L, setdiff(names(metrics), c("trait", "model"))]))))

  point <- proc$classification_reliability$out_of_fold_predictions
  expect_true(all(c("Probability_low", "Probability_medium", "Probability_high") %in% names(point)))
  expect_equal(rowSums(point[, c("Probability_low", "Probability_medium", "Probability_high")]),
               rep(1, nrow(point)), tolerance = 1e-6)
  expect_true(all(is.finite(point$Reliability)))
})
