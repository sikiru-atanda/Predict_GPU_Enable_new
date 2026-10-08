if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

mk_hybrid_ml_runtime_data <- function() {
  set.seed(99)
  ph <- data.frame(
    HybridID = c(
      "A_M1", "A_M2", "B_M1", "B_M2",
      "C_M1", "C_M2", "D_M1", "D_M2",
      "E_M5", "F_M3"
    ),
    Female = c("A", "A", "B", "B", "C", "C", "D", "D", "E", "F"),
    Male = c("M1", "M2", "M1", "M2", "M1", "M2", "M1", "M2", "M5", "M3"),
    stringsAsFactors = FALSE
  )
  female_eff <- c(A = 0.8, B = 0.4, C = -0.3, D = 0.2, E = 1.0, F = -0.5)
  male_eff <- c(M1 = 0.5, M2 = -0.2, M3 = 0.7, M5 = 1.1)
  y <- female_eff[ph$Female] + male_eff[ph$Male] + rnorm(nrow(ph), sd = 0.05)
  y[c(2, 10)] <- NA_real_
  ph$Yield <- as.numeric(y)

  # Hybrid-level genotypes must be dosages (0/1/2 for diploids); genotype
  # QC rejects values above 2 without ploidy metadata. Columns 1-2 carry the
  # female/male effects as dosage classes and column 3 their interaction.
  dose <- function(v) findInterval(v, stats::quantile(v, c(1 / 3, 2 / 3)))
  X <- matrix(sample(0:2, nrow(ph) * 12, replace = TRUE), nrow = nrow(ph), ncol = 12)
  storage.mode(X) <- "double"
  rownames(X) <- ph$HybridID
  colnames(X) <- paste0("m", seq_len(ncol(X)))
  X[, 1] <- dose(female_eff[ph$Female])
  X[, 2] <- dose(male_eff[ph$Male])
  X[, 3] <- (X[, 1] * X[, 2]) %/% 2

  list(pheno = ph, geno = X)
}

mk_hybrid_ml_met_runtime_data <- function() {
  dat <- mk_hybrid_ml_runtime_data()
  base <- subset(dat$pheno, !is.na(HybridID))
  ph1 <- base
  ph1$Env <- "E1"
  ph2 <- base
  ph2$Env <- "E2"
  ph2$Yield <- ph2$Yield + 0.6
  ph <- rbind(ph1, ph2)
  ph$Yield[c(2, 13)] <- NA_real_
  list(pheno = ph, geno = dat$geno)
}

configure_test_python_runtime <- function() {
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  skip_if(is.null(py) || !file.exists(py), "No preferred Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")
  invisible(py)
}

expect_hybrid_ml_contract <- function(preds, include_env = FALSE, info = NULL) {
  expect_identical(
    names(preds),
    PredictProR:::gp_hybrid_prediction_columns(include_env = include_env),
    info = info
  )
  expect_true(all(is.finite(preds$Predicted_value)), info = info)
  expect_true(all(is.finite(preds$Standard_error)), info = info)
  expect_true(all(is.finite(preds$PEV)), info = info)
  expect_equal(preds$PEV, preds$Standard_error^2, tolerance = 1e-8, info = info)
  # expect_gt() has no `info` argument in current testthat
  expect_true(diff(range(preds$Standard_error)) / median(preds$Standard_error) > 0.01, info = info)
  expect_true(all(is.finite(preds$lower_bound)), info = info)
  expect_true(all(is.finite(preds$upper_bound)), info = info)
  expect_true(all(preds$lower_bound <= preds$Predicted_value), info = info)
  expect_true(all(preds$upper_bound >= preds$Predicted_value), info = info)
  expect_true(all(is.finite(preds$Reliability)), info = info)
  expect_true(all(preds$Reliability >= 0 & preds$Reliability <= 1), info = info)
  expect_true(all(is.finite(preds$Reliability_variance_input)), info = info)
  expect_true(all(is.finite(preds$Prediction_stability)), info = info)
  expect_true(
    all(grepl("not genetic reliability", preds$Reliability_basis, fixed = TRUE)),
    info = info
  )
  expect_true(
    all(grepl("predictive", preds$PEV_basis, ignore.case = TRUE)),
    info = info
  )
  expect_true(all(grepl("heldout|cross.fitted|local", preds$Prediction_uncertainty_source, ignore.case = TRUE)), info = info)
}

test_that("hybrid gaussian ML true prediction returns expected outputs", {
  dat <- mk_hybrid_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    ntree = 100
  )

  preds <- out$model_results$predicted_values
  expect_hybrid_ml_contract(preds)
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable"
  ) %in% out$model_results$variance_components$Component))
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_equal(sum(preds$Train_Test_Label == "Test"), sum(is.na(dat$pheno$Yield)))
  expect_identical(
    names(out$model_results$diagnostic_plots),
    c("predicted_vs_observed", "predicted_vs_reliability", "prediction_interval", "reliability_summary")
  )
  expect_true(all(c("summary_statistics", "hybrid_component_summary", "hybrid_train_test_summary") %in% names(out$summary_statistic)))
  expect_true(is.data.frame(out$summary_statistic$hybrid_component_summary))
  expect_true(is.data.frame(out$summary_statistic$hybrid_train_test_summary))
})

test_that("hybrid ML keeps heterozygous hybrid markers by default", {
  # Hybrid-level genotypes are heterozygous by design. The inbred-line
  # heterozygosity filter (het_threshold = 0.1) used to drop every marker here.
  dat <- mk_hybrid_ml_runtime_data()
  dat$geno[] <- 1
  dat$geno[cbind(seq_len(nrow(dat$geno)), (seq_len(nrow(dat$geno)) %% ncol(dat$geno)) + 1L)] <- 0
  dat$geno[cbind(seq_len(nrow(dat$geno)), ((seq_len(nrow(dat$geno)) + 5L) %% ncol(dat$geno)) + 1L)] <- 2
  expect_true(all(colMeans(dat$geno == 1) > 0.1))

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    ntree = 50
  )
  expect_true(all(is.finite(out$model_results$predicted_values$Predicted_value)))
})

test_that("hybrid ML builds F1 features from QC'd inbred parent genotypes", {
  dat <- mk_hybrid_ml_runtime_data()
  set.seed(7)
  fem <- unique(dat$pheno$Female)
  mal <- unique(dat$pheno$Male)
  parents <- matrix(sample(c(0, 2), (length(fem) + length(mal)) * 15, replace = TRUE),
                    ncol = 15, dimnames = list(c(fem, mal), paste0("p", 1:15)))
  parents[fem[1:3], "p1"] <- 1  # residual heterozygosity in the female pool

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = parents[fem, ],
    male_geno_data = parents[mal, ],
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    ntree = 50
  )
  expect_hybrid_ml_contract(out$model_results$predicted_values)
  qc <- setNames(out$hybrid_parent_qc$value, out$hybrid_parent_qc$item)
  expect_equal(qc[["hybrid_feature_source"]], "expected_F1_from_parents")
  expect_equal(qc[["markers_removed_parent_heterozygosity"]], "1")
})

test_that("hybrid gaussian ML CV returns scenario-aware held-out outputs", {
  dat <- mk_hybrid_ml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_One_New_Parent",
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    ntree = 100
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed
  expect_true(all(c("trait", "model", "rep", "eval_metrics_reps", "ypred_cv_Reps_all") %in% names(raw)))
  expect_true(all(c(
    "HybridID", "Female", "Male", "Observed_value", "Predicted_value",
    "cv_scenario", "fold", "rep", "cv_role",
    "Prediction_intercept",
    "Female_additive_contribution",
    "Male_additive_contribution",
    "Hybrid_interaction_contribution",
    "Female_GCA_like_contribution",
    "Male_GCA_like_contribution",
    "SCA_like_contribution"
  ) %in% names(raw$ypred_cv_Reps_all)))
  expect_true(all(raw$ypred_cv_Reps_all$Train_Test_Label == "Test"))
  expect_equal(
    raw$ypred_cv_Reps_all$Predicted_value,
    raw$ypred_cv_Reps_all$Prediction_intercept +
      raw$ypred_cv_Reps_all$Female_additive_contribution +
      raw$ypred_cv_Reps_all$Male_additive_contribution +
      raw$ypred_cv_Reps_all$Hybrid_interaction_contribution,
    tolerance = 1e-8
  )
  expect_true(all(c("hybrid_metric_summary", "hybrid_prediction_counts", "hybrid_cv_predictions", "hybrid_cv_plot") %in% names(proc)))
  expect_true(any(proc$hybrid_cv_predictions$cv_scenario == "One_New_Parent"))
  expect_true(inherits(proc$hybrid_cv_plot, "ggplot"))
})

test_that("hybrid gaussian Ridge ML runtime and CV work", {
  configure_test_python_runtime()
  dat <- mk_hybrid_ml_runtime_data()

  out_true <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "Ridge_Regression",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    lambda_rr = 0.1
  )
  expect_hybrid_ml_contract(out_true$model_results$predicted_values)

  out_cv <- PredictProR::model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "Ridge_Regression",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_Both_New_Parents",
    replication = 1L,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    lambda_rr = 0.1
  )
  expect_true(nrow(out_cv$cv_results_raw[[1]]$ypred_cv_Reps_all) > 0)
  expect_true(any(out_cv$cv_results_processed$hybrid_cv_predictions$cv_scenario == "Both_New_Parents"))
})

test_that("all remaining supported hybrid gaussian ML models honor the public contract", {
  configure_test_python_runtime()
  dat <- mk_hybrid_ml_runtime_data()
  models <- setdiff(
    PredictProR:::gp_hybrid_ml_supported_models(),
    c("RandomForest", "Ridge_Regression")
  )

  for (model in models) {
    out <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = "Yield",
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      GS_model = model,
      response_family = "gaussian",
      hybrid_ml = TRUE,
      system_database = TRUE,
      message = FALSE,
      n_bootstrap = 5L,
      ncomp = 2L,
      iteration = 20L,
      catboost_iterations = 20L,
      catboost_depth = 3L,
      xgb_nthread = 1L,
      catboost_thread_count = 1L
    )
    expect_hybrid_ml_contract(out$model_results$predicted_values, info = model)
  }
})

test_that("hybrid gaussian ML supports multi-environment rows", {
  dat <- mk_hybrid_ml_met_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    heter_groups = "Env",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    ntree = 100
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_ml_contract(preds, include_env = TRUE)
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))

  out_cv <- PredictProR::model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    heter_groups = "Env",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    # Multi-environment hybrid CV uses CV0/CV1/CV2; parent-novelty scenarios
    # are single-environment only.
    cross_validation_meth = "CV1",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    ntree = 100
  )
  expect_true("Env" %in% names(out_cv$cv_results_processed$hybrid_cv_predictions))
  expect_true(nrow(out_cv$cv_results_processed$hybrid_cv_predictions) > 0)
})
