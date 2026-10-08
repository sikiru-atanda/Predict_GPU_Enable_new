mk_hybrid_gp_runtime_data <- function() {
  set.seed(20260503)
  female <- paste0("F", 1:4)
  male <- paste0("M", 1:4)
  parents <- c(female, male)
  X <- matrix(sample(0:2, length(parents) * 18L, replace = TRUE), nrow = length(parents), ncol = 18L)
  rownames(X) <- parents
  colnames(X) <- paste0("m", seq_len(ncol(X)))

  crosses <- expand.grid(
    Female = female,
    Male = male,
    stringsAsFactors = FALSE
  )
  crosses <- crosses[c(1, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 16), , drop = FALSE]
  female_signal <- c(0.9, 0.4, -0.3, -0.7)
  male_signal <- c(0.5, -0.2, 0.6, -0.4)
  sca_signal <- c(0.15, -0.05, 0.10, 0.20, -0.10, 0.05, -0.20, 0.12, -0.08, 0.04, -0.03, 0.09)
  y <- female_signal[match(crosses$Female, female)] +
    male_signal[match(crosses$Male, male)] +
    sca_signal +
    stats::rnorm(nrow(crosses), sd = 0.03)
  moisture <- 0.6 * female_signal[match(crosses$Female, female)] -
    0.4 * male_signal[match(crosses$Male, male)] +
    rev(sca_signal) +
    stats::rnorm(nrow(crosses), sd = 0.03)
  y[c(4, 10)] <- NA_real_
  moisture[c(3, 11)] <- NA_real_

  list(
    pheno = data.frame(
      HybridID = paste(crosses$Female, crosses$Male, sep = "_x_"),
      Female = crosses$Female,
      Male = crosses$Male,
      Yield = y,
      Moisture = moisture,
      stringsAsFactors = FALSE
    ),
    geno = X,
    female_geno = X[female, , drop = FALSE],
    male_geno = X[male, , drop = FALSE]
  )
}

mk_hybrid_gp_met_runtime_data <- function() {
  dat <- mk_hybrid_gp_runtime_data()
  ph <- dat$pheno
  envs <- c("E1", "E2")
  met <- do.call(rbind, lapply(seq_along(envs), function(i) {
    out <- ph
    out$Env <- envs[[i]]
    out$HybridID <- paste(out$HybridID, envs[[i]], sep = "_")
    out$Yield <- out$Yield + c(0.35, -0.20)[[i]]
    out$Moisture <- out$Moisture + c(-0.15, 0.25)[[i]]
    out
  }))
  rownames(met) <- NULL
  dat$pheno <- met
  dat
}

mk_hybrid_gp_env_similarity <- function() {
  K <- matrix(c(1, 0.45, 0.45, 1), nrow = 2L)
  rownames(K) <- colnames(K) <- c("E1", "E2")
  K
}

mk_hybrid_gp_env_covariates <- function() {
  data.frame(
    Env = c("E1", "E2"),
    heat_index = c(0.25, 1.30),
    rainfall = c(1.10, 0.35),
    stringsAsFactors = FALSE
  )
}

expect_hybrid_public_contract <- function(preds, include_env = FALSE, include_trait = FALSE) {
  expect_identical(
    names(preds),
    PredictProR:::gp_hybrid_prediction_columns(
      include_env = include_env,
      include_trait = include_trait
    )
  )
  expect_true(all(is.finite(preds$Predicted_value)))
  expect_true(all(is.finite(preds$Standard_error)))
  expect_true(all(is.finite(preds$PEV)))
  expect_equal(preds$Standard_error^2, preds$PEV, tolerance = 1e-7)
  expect_true(all(is.finite(preds$lower_bound)))
  expect_true(all(is.finite(preds$upper_bound)))
  expect_true(all(is.finite(preds$Reliability)))
  expect_true(all(is.finite(preds$Reliability_variance_input)))
  expect_true(all(is.finite(preds$Reliability_reference_variance)))
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
}

expect_public_plot_contract <- function(plots) {
  expect_true(is.list(plots))
  expect_identical(
    names(plots),
    c(
      "predicted_vs_observed",
      "predicted_vs_reliability",
      "prediction_interval",
      "reliability_summary"
    )
  )
}

test_that("hybrid GP true prediction returns public hybrid uncertainty outputs", {
  dat <- mk_hybrid_gp_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_public_contract(preds)
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_public_plot_contract(out$model_results$diagnostic_plots)
  expect_true(is.data.frame(out$model_results$variance_components))
  expect_identical(
    names(out$model_results$variance_components),
    PredictProR:::gp_hybrid_variance_component_columns(include_trait = FALSE, include_env = FALSE)
  )
})

test_that("all supported hybrid GP model names honor the public uncertainty contract", {
  dat <- mk_hybrid_gp_runtime_data()
  models <- PredictProR:::gp_hybrid_gp_supported_models()

  expect_identical(models, "GP")

  for (model in models) {
    out <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      gmatrix_method = "VanRaden",
      response = "Yield",
      gen_name = "HybridID",
      GS_model = model,
      response_family = "gaussian",
      hybrid_gp = TRUE,
      female_parent = "Female",
      male_parent = "Male",
      system_database = TRUE,
      message = FALSE
    )
    expect_hybrid_public_contract(out$model_results$predicted_values)
  }
})

test_that("hybrid GP true prediction supports single-trait multi-environment data", {
  dat <- mk_hybrid_gp_met_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    heter_groups = "Env",
    env_similarity = mk_hybrid_gp_env_similarity(),
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_hybrid_public_contract(preds, include_env = TRUE)
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
  expect_true(any(out$model_results$variance_components$Component %in% c("gxe", "gxe_variance")))
})

test_that("hybrid GP derives multi-environment kernel from environment covariates", {
  dat <- mk_hybrid_gp_met_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    heter_groups = "Env",
    env_covariates = mk_hybrid_gp_env_covariates(),
    kenv_kernel = "rbf",
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_hybrid_public_contract(preds, include_env = TRUE)
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
})

test_that("hybrid GP true prediction supports multi-trait single-environment data", {
  dat <- mk_hybrid_gp_runtime_data()

  out <- NULL
  expect_warning(
    out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = c("Yield", "Moisture"),
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
    ),
    "trait-specific inferred testing sets"
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), 2L * nrow(dat$pheno))
  expect_hybrid_public_contract(preds, include_trait = TRUE)
  expect_equal(sort(unique(preds$Trait)), c("Moisture", "Yield"))
  expect_true("Trait" %in% names(out$model_results$variance_components))
  expect_true(is.matrix(out$model_results$Genetic_correlation_traits))
  expect_equal(dim(out$model_results$Genetic_correlation_traits), c(2L, 2L))
  expect_equal(unname(diag(out$model_results$Genetic_correlation_traits)), c(1, 1), tolerance = 1e-7)
})

test_that("hybrid GP true prediction supports multi-trait multi-environment data", {
  dat <- mk_hybrid_gp_met_runtime_data()

  out <- NULL
  expect_warning(
    out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = c("Yield", "Moisture"),
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    heter_groups = "Env",
    env_covariates = mk_hybrid_gp_env_covariates(),
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
    ),
    "trait-specific inferred testing sets"
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), 2L * nrow(dat$pheno))
  expect_hybrid_public_contract(preds, include_env = TRUE, include_trait = TRUE)
  expect_equal(sort(unique(preds$Trait)), c("Moisture", "Yield"))
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
  expect_true(is.matrix(out$model_results$Genetic_correlation_traits))
  expect_true(is.matrix(out$model_results$Genetic_correlation_trait_environment))
})

test_that("hybrid GP CV returns parent-aware scenario summaries", {
  dat <- mk_hybrid_gp_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model_cv = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    cross_validation = TRUE,
    cross_validation_meth = "Hybrid_Known_Parents",
    nfolds = 3,
    replication = 1,
    eval_metrics = c("mean_absolute_error", "root_mean_squared_error"),
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  expect_true(is.list(out$cv_results_processed))
  expect_true(is.data.frame(out$cv_results_processed$hybrid_metric_summary))
  expect_true(is.data.frame(out$cv_results_processed$hybrid_cv_predictions))
  expect_true(any(out$cv_results_processed$hybrid_cv_predictions$cv_scenario == "Known_Parents_New_Hybrid"))
  expect_true(all(c("Female_GCA", "Male_GCA", "SCA_effect") %in% names(out$cv_results_processed$hybrid_cv_predictions)))
})

test_that("hybrid GP CV supports multi-trait multi-environment CV0 predictions", {
  dat <- mk_hybrid_gp_met_runtime_data()

  out <- NULL
  expect_warning(
    out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = c("Yield", "Moisture"),
    gen_name = "HybridID",
    GS_model_cv = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    heter_groups = "Env",
    env_covariates = mk_hybrid_gp_env_covariates(),
    cross_validation = TRUE,
    cross_validation_meth = "CV0",
    nfolds = 2,
    replication = 1,
    eval_metrics = c("mean_absolute_error", "root_mean_squared_error"),
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
    ),
    "trait-specific inferred testing sets"
  )

  preds <- out$cv_results_processed$hybrid_cv_predictions
  expect_true(is.data.frame(preds))
  expect_true(all(c("Trait", "Env", "Standard_error", "Prediction_error_variance") %in% names(preds)))
  expect_equal(sort(unique(preds$Trait)), c("Moisture", "Yield"))
  expect_true(all(is.finite(preds$Predicted_value)))
  expect_true(all(is.finite(preds$GxE_effect)))
  expect_equal(preds$Standard_error^2, preds$Prediction_error_variance, tolerance = 1e-7)
  expect_true("Trait" %in% names(out$cv_results_processed$hybrid_metric_summary))
})

test_that("hybrid GP can force the R backend", {
  dat <- mk_hybrid_gp_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    gp_backend = "r",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_hybrid_public_contract(preds)
})

test_that("hybrid GP Python/Torch backend matches the R solver when available", {
  py <- tryCatch(PredictProR:::gp_detect_gp_python(), error = function(e) NULL)
  testthat::skip_if(is.null(py), "No configured GP Python runtime.")

  dat <- mk_hybrid_gp_runtime_data()
  base_args <- list(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    hybrid_gp_lambda = 0.1,
    system_database = TRUE,
    message = FALSE
  )
  out_r <- do.call(PredictProR::model_execute, c(base_args, list(gp_backend = "r")))
  out_py <- tryCatch(
    do.call(PredictProR::model_execute, c(base_args, list(gp_backend = "python"))),
    error = function(e) testthat::skip(paste("Hybrid GP Python/Torch backend unavailable:", conditionMessage(e)))
  )

  pred_r <- out_r$model_results$predicted_values
  pred_py <- out_py$model_results$predicted_values
  expect_equal(pred_py$Predicted_value, pred_r$Predicted_value, tolerance = 1e-5)
  expect_equal(pred_py$PEV, pred_r$PEV, tolerance = 1e-5)
})

test_that("joint multi-trait hybrid GP Python/Torch backend matches the R solver when available", {
  py <- tryCatch(PredictProR:::gp_detect_gp_python(), error = function(e) NULL)
  testthat::skip_if(is.null(py), "No configured GP Python runtime.")

  dat <- mk_hybrid_gp_runtime_data()
  dat$pheno$Moisture[is.na(dat$pheno$Moisture)] <- mean(dat$pheno$Moisture, na.rm = TRUE)
  dat$pheno$Moisture[is.na(dat$pheno$Yield)] <- NA_real_
  base_args <- list(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = c("Yield", "Moisture"),
    gen_name = "HybridID",
    GS_model = "GP",
    response_family = "gaussian",
    hybrid_gp = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    hybrid_gp_lambda = 0.1,
    system_database = TRUE,
    message = FALSE
  )
  out_r <- do.call(PredictProR::model_execute, c(base_args, list(gp_backend = "r")))
  out_py <- tryCatch(
    do.call(PredictProR::model_execute, c(base_args, list(gp_backend = "python"))),
    error = function(e) testthat::skip(paste("Joint hybrid GP Python/Torch backend unavailable:", conditionMessage(e)))
  )

  pred_r <- out_r$model_results$predicted_values
  pred_py <- out_py$model_results$predicted_values
  expect_equal(pred_py$Predicted_value, pred_r$Predicted_value, tolerance = 1e-5)
  expect_equal(pred_py$PEV, pred_r$PEV, tolerance = 1e-5)
  expect_equal(
    out_py$model_results$Genetic_correlation_traits,
    out_r$model_results$Genetic_correlation_traits,
    tolerance = 1e-7
  )
})
