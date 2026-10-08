if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

local_predictdl_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_PYTHON = py)
  invisible(py)
}

mk_multitrait_runtime_data <- function(n = 20L, p = 8L) {
  set.seed(42)
  gids <- paste0("g", seq_len(n))
  # Traits from a latent signal; geno_data carries it as homozygous 0/2
  # dosages (genotype QC rejects values above 2 without ploidy and removes
  # markers with more than 10% heterozygotes).
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- ifelse(Z > 0, 2, 0)
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

test_that("multi-trait gaussian DL runtime returns long and wide outputs", {
  local_predictdl_runtime()

  dat <- mk_multitrait_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model = "mlp",
    response_family = "gaussian",
    multi_trait_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    validation_split = 0.2,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    epochs = 8L,
    batch_size = 8L,
    internal_cv_nfolds = 3L,
    internal_cv_replication = 2L,
    mlp_neurons_per_layer = as.integer(c(32L, 16L)),
    mlp_learning_rate = 1e-3,
    dropout = 0.0
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

test_that("multi-trait DL rejects non-gaussian requests", {
  dat <- mk_multitrait_runtime_data()

  expect_error(
    suppressWarnings(PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      GS_model = "mlp",
      response_family = "binary",
      multi_trait_dl = TRUE,
      system_database = TRUE,
      message = FALSE
    )),
    "support gaussian traits only"
  )

})

test_that("multi-trait gaussian DL CV returns held-out fold outputs", {
  local_predictdl_runtime()

  dat <- mk_multitrait_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model_cv = "mlp",
    response_family = "gaussian",
    multi_trait_dl = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    validation_split = 0.2,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    epochs = 8L,
    batch_size = 8L,
    mlp_neurons_per_layer = as.integer(c(32L, 16L)),
    mlp_learning_rate = 1e-3,
    dropout = 0.0
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
