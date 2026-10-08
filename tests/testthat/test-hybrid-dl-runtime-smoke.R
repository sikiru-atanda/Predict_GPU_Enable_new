if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

mk_hybrid_dl_runtime_data <- function() {
  set.seed(199)
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

  p <- 12L
  # Hybrid-level genotypes must be dosages (0/1/2 for diploids); genotype
  # QC rejects values above 2 without ploidy metadata. Columns 1-2 carry the
  # female/male effects as dosage classes and column 3 their interaction.
  dose <- function(v) findInterval(v, stats::quantile(v, c(1 / 3, 2 / 3)))
  X <- matrix(sample(0:2, nrow(ph) * p, replace = TRUE), nrow = nrow(ph), ncol = p)
  storage.mode(X) <- "double"
  rownames(X) <- ph$HybridID
  colnames(X) <- paste0("m", seq_len(p))
  X[, 1] <- dose(female_eff[ph$Female])
  X[, 2] <- dose(male_eff[ph$Male])
  X[, 3] <- (X[, 1] * X[, 2]) %/% 2

  list(pheno = ph, geno = X)
}

mk_hybrid_dl_met_runtime_data <- function() {
  dat <- mk_hybrid_dl_runtime_data()
  ph1 <- dat$pheno
  ph1$Env <- "E1"
  ph2 <- dat$pheno
  ph2$Env <- "E2"
  ph2$Yield <- ph2$Yield + 0.5
  ph <- rbind(ph1, ph2)
  ph$Yield[c(2, 13)] <- NA_real_
  list(pheno = ph, geno = dat$geno)
}

configure_test_dl_python_runtime <- function() {
  py <- if (exists("gp_preferred_python", envir = .GlobalEnv, inherits = FALSE)) {
    get("gp_preferred_python", envir = .GlobalEnv)(purpose = "dl")
  } else {
    PredictProR:::gp_preferred_python(purpose = "dl")
  }
  skip_if(is.null(py) || !file.exists(py), "No preferred DL Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_PYTHON = py)
  invisible(py)
}

expect_hybrid_dl_contract <- function(preds, include_env = FALSE) {
  expect_identical(
    names(preds),
    PredictProR:::gp_hybrid_prediction_columns(include_env = include_env)
  )
  expect_true(all(is.finite(preds$Predicted_value)))
  expect_true(all(is.finite(preds$Standard_error)))
  expect_true(all(is.finite(preds$PEV)))
  expect_equal(preds$PEV, preds$Standard_error^2, tolerance = 1e-8)
  expect_gt(diff(range(preds$Standard_error)) / median(preds$Standard_error), 0.01)
  expect_true(all(is.finite(preds$lower_bound)))
  expect_true(all(is.finite(preds$upper_bound)))
  expect_true(all(preds$lower_bound <= preds$Predicted_value))
  expect_true(all(preds$upper_bound >= preds$Predicted_value))
  expect_true(all(is.finite(preds$Reliability)))
  expect_true(all(preds$Reliability >= 0 & preds$Reliability <= 1))
  expect_true(all(is.finite(preds$Reliability_variance_input)))
  expect_true(all(is.finite(preds$Prediction_stability)))
  expect_true(all(grepl("not genetic reliability", preds$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl(
    "predictive",
    preds$PEV_basis,
    ignore.case = TRUE
  )))
  expect_true(all(grepl("heldout|cross.fitted|local", preds$Prediction_uncertainty_source, ignore.case = TRUE)))
}

test_that("hybrid gaussian DL true prediction returns decomposition outputs", {
  configure_test_dl_python_runtime()
  dat <- mk_hybrid_dl_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "mlp",
    response_family = "gaussian",
    hybrid_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    epochs = 4L,
    batch_size = 8L,
    mlp_neurons_per_layer = as.integer(c(32L, 16L)),
    validation_split = 0.1
  )

  preds <- out$model_results$predicted_values
  expect_hybrid_dl_contract(preds)
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable"
  ) %in% out$model_results$variance_components$Component))
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
})

test_that("hybrid gaussian DL CV returns scenario-aware held-out outputs", {
  configure_test_dl_python_runtime()
  dat <- mk_hybrid_dl_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = subset(dat$pheno, !is.na(Yield)),
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "mlp",
    response_family = "gaussian",
    hybrid_dl = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_One_New_Parent",
    replication = 1L,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    system_database = TRUE,
    message = FALSE,
    epochs = 4L,
    batch_size = 8L,
    mlp_neurons_per_layer = as.integer(c(32L, 16L)),
    validation_split = 0.1
  )

  raw <- out$cv_results_raw[[1]]
  proc <- out$cv_results_processed
  expect_true(all(c(
    "HybridID", "Female", "Male", "Observed_value", "Predicted_value",
    "cv_scenario", "fold", "rep", "cv_role",
    "Prediction_intercept",
    "Female_additive_contribution",
    "Male_additive_contribution",
    "Hybrid_interaction_contribution"
  ) %in% names(raw$ypred_cv_Reps_all)))
  expect_equal(
    raw$ypred_cv_Reps_all$Predicted_value,
    raw$ypred_cv_Reps_all$Prediction_intercept +
      raw$ypred_cv_Reps_all$Female_additive_contribution +
      raw$ypred_cv_Reps_all$Male_additive_contribution +
      raw$ypred_cv_Reps_all$Hybrid_interaction_contribution,
    tolerance = 1e-8
  )
  expect_true(any(proc$hybrid_cv_predictions$cv_scenario == "One_New_Parent"))
})

test_that("hybrid gaussian ft_transformer DL true prediction works", {
  configure_test_dl_python_runtime()
  dat <- mk_hybrid_dl_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "ft_transformer",
    response_family = "gaussian",
    hybrid_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    epochs = 3L,
    batch_size = 8L,
    ft_d_model = 32L,
    ft_heads = 4L,
    ft_layers = 2L,
    ft_ff_mult = 2L,
    validation_split = 0.1
  )

  expect_hybrid_dl_contract(out$model_results$predicted_values)
})

test_that("hybrid gaussian DL supports multi-environment rows", {
  configure_test_dl_python_runtime()
  dat <- mk_hybrid_dl_met_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "HybridID",
    heter_groups = "Env",
    female_parent = "Female",
    male_parent = "Male",
    GS_model = "mlp",
    response_family = "gaussian",
    hybrid_dl = TRUE,
    system_database = TRUE,
    message = FALSE,
    n_bootstrap = 5L,
    epochs = 3L,
    batch_size = 8L,
    mlp_neurons_per_layer = as.integer(c(32L, 16L)),
    validation_split = 0.1
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_dl_contract(preds, include_env = TRUE)
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
})
