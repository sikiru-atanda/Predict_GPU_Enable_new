if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

mk_hybrid_bayes_runtime_data <- function() {
  set.seed(314)
  female <- paste0("F", 1:3)
  male <- paste0("M", 1:3)
  parents <- c(female, male)
  X <- matrix(sample(0:2, length(parents) * 18L, replace = TRUE), nrow = length(parents), ncol = 18L)
  rownames(X) <- parents
  colnames(X) <- paste0("m", seq_len(ncol(X)))

  crosses <- expand.grid(
    Female = female,
    Male = male,
    stringsAsFactors = FALSE
  )
  female_signal <- c(0.9, 0.1, -0.5)
  male_signal <- c(0.4, -0.2, 0.8)
  sca_signal <- c(0.15, -0.08, 0.05, 0.22, -0.10, 0.04, -0.18, 0.09, -0.06)
  y <- female_signal[match(crosses$Female, female)] +
    male_signal[match(crosses$Male, male)] +
    sca_signal +
    rnorm(nrow(crosses), sd = 0.03)
  y[c(2, 8)] <- NA_real_

  pheno <- data.frame(
    HybridID = paste(crosses$Female, crosses$Male, sep = "_x_"),
    Female = crosses$Female,
    Male = crosses$Male,
    Yield = y,
    stringsAsFactors = FALSE
  )

  list(
    pheno = pheno,
    geno = X,
    female_geno = X[female, , drop = FALSE],
    male_geno = X[male, , drop = FALSE]
  )
}

mk_hybrid_bayes_met_runtime_data <- function() {
  dat <- mk_hybrid_bayes_runtime_data()
  obs <- subset(dat$pheno, !is.na(Yield))
  ph1 <- obs
  ph1$Env <- "E1"
  ph2 <- obs
  ph2$Env <- "E2"
  ph2$Yield <- ph2$Yield + 0.4
  ph <- rbind(ph1, ph2)
  ph$Yield[c(2, 11)] <- NA_real_
  list(
    pheno = ph,
    geno = dat$geno,
    female_geno = dat$female_geno,
    male_geno = dat$male_geno
  )
}

hybrid_bayes_model_execute <- function(...) {
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::model_execute(...))
  }
  get("model_execute", envir = .GlobalEnv)(...)
}

expect_hybrid_bayes_contract <- function(preds, include_env = FALSE) {
  expect_identical(
    names(preds),
    PredictProR:::gp_hybrid_prediction_columns(include_env = include_env)
  )
  expect_true(all(is.finite(preds$Predicted_value)))
  expect_true(all(is.finite(preds$Standard_error)))
  expect_true(all(is.finite(preds$PEV)))
  expect_equal(preds$Standard_error^2, preds$PEV, tolerance = 1e-7)
  expect_true(all(is.finite(preds$lower_bound)))
  expect_true(all(is.finite(preds$upper_bound)))
  expect_true(all(preds$lower_bound <= preds$Predicted_value))
  expect_true(all(preds$Predicted_value <= preds$upper_bound))
  expect_true(all(is.finite(preds$Reliability)))
  expect_true(all(is.finite(preds$Reliability_variance_input)))
  expect_true(all(is.finite(preds$Reliability_reference_variance)))
  expect_equal(preds$Reliability_variance_input, preds$PEV, tolerance = 1e-10)
  ratio_fallback <- grepl(
    "ratio-form fallback used",
    ifelse(is.na(preds$Reliability_remarks), "", preds$Reliability_remarks),
    fixed = TRUE
  )
  standard_reliability <- pmax(0, pmin(
    1,
    1 - preds$Reliability_variance_input / preds$Reliability_reference_variance
  ))
  ratio_reliability <- preds$Reliability_reference_variance /
    (preds$Reliability_reference_variance + preds$Reliability_variance_input)
  expect_equal(
    preds$Reliability,
    ifelse(ratio_fallback, ratio_reliability, standard_reliability),
    tolerance = 1e-7
  )
  expect_true(all(grepl("BGLR model-implied hybrid reliability", preds$Reliability_basis, fixed = TRUE)))
  expect_false(any(grepl("ML/DL marker-adjustment", preds$Reliability_basis, fixed = TRUE)))
}

expect_hybrid_bayes_variance_components <- function(vc) {
  expected <- c(
    "female_gca_variance", "male_gca_variance", "sca_variance",
    "total_genetic_variance", "residual_variance", "heritability"
  )
  expect_true(all(expected %in% vc$Component))
  idx <- match(expected, vc$Component)
  expect_true(all(is.finite(vc$Components[idx])))
  expect_true(all(is.finite(vc$Standard_error[idx])))
  expect_true(all(vc$Components[idx[seq_len(5L)]] >= 0))
  expect_true(vc$Components[idx[[6L]]] >= 0 && vc$Components[idx[[6L]]] <= 1)
}

test_that("hybrid gaussian GBLUP_BRR true prediction returns decomposition outputs", {
  dat <- mk_hybrid_bayes_runtime_data()

  out <- hybrid_bayes_model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "GBLUP_BRR",
    response_family = "gaussian",
    hybrid_bayes = TRUE,
    system_database = TRUE,
    message = FALSE,
    nIter = 120,
    burnIn = 40,
    thin = 2
  )

  preds <- out$model_results$predicted_values
  expect_hybrid_bayes_contract(preds)
  expect_hybrid_bayes_variance_components(out$model_results$variance_components)
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_identical(
    names(out$model_results$diagnostic_plots),
    c("predicted_vs_observed", "predicted_vs_reliability", "prediction_interval", "reliability_summary")
  )
})

test_that("hybrid gaussian RKHS true prediction accepts separate parent genotype tables", {
  dat <- mk_hybrid_bayes_runtime_data()

  out <- hybrid_bayes_model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "RKHS",
    response_family = "gaussian",
    hybrid_bayes = TRUE,
    system_database = TRUE,
    message = FALSE,
    nIter = 120,
    burnIn = 40,
    thin = 2
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_bayes_contract(preds)
  expect_hybrid_bayes_variance_components(out$model_results$variance_components)
})

test_that("hybrid Bayesian known-parents CV returns scenario summaries", {
  dat <- mk_hybrid_bayes_runtime_data()

  out <- hybrid_bayes_model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "GBLUP_BRR",
    response_family = "gaussian",
    hybrid_bayes = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_Known_Parents",
    nfolds = 3,
    replication = 1,
    eval_metrics = c("mean_absolute_error", "root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    nIter = 120,
    burnIn = 40,
    thin = 2
  )

  proc <- out$cv_results_processed
  expect_true(is.data.frame(proc$hybrid_metric_summary))
  expect_true(is.data.frame(proc$hybrid_cv_predictions))
  expect_true(any(proc$hybrid_cv_predictions$cv_scenario == "Known_Parents_New_Hybrid"))
  expect_equal(
    proc$hybrid_cv_predictions$Predicted_value,
    proc$hybrid_cv_predictions$Prediction_intercept +
      proc$hybrid_cv_predictions$Female_additive_contribution +
      proc$hybrid_cv_predictions$Male_additive_contribution +
      proc$hybrid_cv_predictions$Hybrid_interaction_contribution,
    tolerance = 1e-6
  )
})

test_that("hybrid Bayesian one-new-parent CV works for RKHS", {
  dat <- mk_hybrid_bayes_runtime_data()

  out <- hybrid_bayes_model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "RKHS",
    response_family = "gaussian",
    hybrid_bayes = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_One_New_Parent",
    nfolds = 3,
    replication = 1,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    nIter = 120,
    burnIn = 40,
    thin = 2
  )

  preds <- out$cv_results_processed$hybrid_cv_predictions
  expect_true(nrow(preds) > 0)
  expect_true(all(preds$cv_scenario == "One_New_Parent"))
})

test_that("hybrid Bayesian models support multi-environment rows", {
  dat <- mk_hybrid_bayes_met_runtime_data()

  out <- hybrid_bayes_model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    heter_groups = "Env",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "GBLUP_BRR",
    response_family = "gaussian",
    hybrid_bayes = TRUE,
    system_database = TRUE,
    message = FALSE,
    nIter = 120,
    burnIn = 40,
    thin = 2
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_bayes_contract(preds, include_env = TRUE)
  expect_hybrid_bayes_variance_components(out$model_results$variance_components)
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
})
