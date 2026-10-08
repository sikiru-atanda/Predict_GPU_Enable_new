general_standard_call <- function(...) {
  if (exists("general_prediction_data_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("general_prediction_data_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::general_prediction_data_standard(...))
  }
  stop("general_prediction_data_standard is unavailable for this test.", call. = FALSE)
}

general_validate_call <- function(...) {
  if (exists("validate_general_input_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("validate_general_input_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::validate_general_input_standard(...))
  }
  stop("validate_general_input_standard is unavailable for this test.", call. = FALSE)
}

general_model_execute_call <- function(...) {
  if (exists("model_execute", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("model_execute", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::model_execute(...))
  }
  stop("model_execute is unavailable for this test.", call. = FALSE)
}

mk_general_contract_data <- function() {
  ph_single <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1.1, 2.2, 3.3),
    stringsAsFactors = FALSE
  )
  ph_multi_env <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1.1, 1.2, 2.1, 2.3),
    stringsAsFactors = FALSE
  )
  X <- matrix(c(
    0, 1,
    1, 0,
    1, 1
  ), nrow = 3, byrow = TRUE)
  rownames(X) <- ph_single$GID
  colnames(X) <- c("m1", "m2")
  G <- tcrossprod(X)
  rownames(G) <- colnames(G) <- ph_single$GID
  list(ph_single = ph_single, ph_multi_env = ph_multi_env, X = X, G = G)
}

test_that("general_prediction_data_standard returns the package contract", {
  spec <- general_standard_call()
  expect_true(is.list(spec))
  expect_true(all(c("phenotype", "model_families", "genomic_contracts", "multi_environment_notes") %in% names(spec)))
})

test_that("validate_general_input_standard accepts a valid single-environment ML contract", {
  dat <- mk_general_contract_data()
  out <- general_validate_call(
    pheno_data = dat$ph_single,
    geno_data = dat$X,
    response = "Yield",
    gen_name = "GID",
    GS_model = "Ridge_Regression",
    AI_valid_models = c("Ridge_Regression"),
    bayes_valid_models = c("BRR"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    asreml_model = "GBLUP"
  )
  expect_true(isTRUE(out$valid))
})

test_that("validate_general_input_standard rejects ML without feature data", {
  dat <- mk_general_contract_data()
  expect_error(
    general_validate_call(
      pheno_data = dat$ph_single,
      response = "Yield",
      gen_name = "GID",
      GS_model = "Ridge_Regression",
      AI_valid_models = c("Ridge_Regression"),
      bayes_valid_models = c("BRR"),
      bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
      asreml_model = "GBLUP"
    ),
    regexp = "Genomic standard for single-environment ML/DL and Bayesian marker models"
  )
})

test_that("validate_general_input_standard rejects repeated-row phenotype without heter_groups", {
  dat <- mk_general_contract_data()
  expect_error(
    general_validate_call(
      pheno_data = dat$ph_multi_env,
      geno_data = dat$X[1:2, , drop = FALSE],
      response = "Yield",
      gen_name = "GID",
      GS_model = "RandomForest",
      AI_valid_models = c("RandomForest"),
      bayes_valid_models = c("BRR"),
      bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
      asreml_model = "GBLUP"
    ),
    regexp = "Repeated genotype rows were detected"
  )
})

test_that("validate_general_input_standard treats GP models as kernel relationship models", {
  dat <- mk_general_contract_data()
  gp_models <- PredictProR:::gp_lowrank_supported_models()
  expect_error(
    general_validate_call(
      pheno_data = dat$ph_multi_env,
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = "KRR",
      AI_valid_models = c("RandomForest"),
      bayes_valid_models = c("BRR"),
      bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
      gp_valid_models = gp_models,
      asreml_model = "GBLUP"
    ),
    regexp = "selected kernel/relationship models require"
  )

  out <- general_validate_call(
    pheno_data = dat$ph_multi_env,
    gmatrix = dat$G,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "KRR",
    AI_valid_models = c("RandomForest"),
    bayes_valid_models = c("BRR"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    gp_valid_models = gp_models,
    asreml_model = "GBLUP"
  )
  expect_true(isTRUE(out$valid))
})

test_that("model_execute fails early for general ML input without genomic features", {
  dat <- mk_general_contract_data()
  expect_error(
    general_model_execute_call(
      pheno_data = dat$ph_single,
      response = "Yield",
      gen_name = "GID",
      response_family = "gaussian",
      GS_model = "Ridge_Regression",
      system_database = TRUE,
      message = FALSE
    ),
    regexp = "Genomic standard for single-environment ML/DL and Bayesian marker models"
  )
})

test_that("model_execute fails early for repeated-row phenotype without heter_groups", {
  dat <- mk_general_contract_data()
  expect_error(
    general_model_execute_call(
      pheno_data = dat$ph_multi_env,
      geno_data = dat$X[1:2, , drop = FALSE],
      response = "Yield",
      gen_name = "GID",
      response_family = "gaussian",
      GS_model = "RandomForest",
      system_database = TRUE,
      message = FALSE
    ),
    regexp = "Repeated genotype rows were detected"
  )
})

test_that("model_execute MET marker-Bayes guardrail points users to GP-capable models", {
  dat <- mk_general_contract_data()
  err <- tryCatch(
    general_model_execute_call(
      pheno_data = dat$ph_multi_env,
      geno_data = dat$X[1:2, , drop = FALSE],
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      random = ~GID,
      response_family = "gaussian",
      GS_model = "BRR",
      system_database = TRUE,
      message = FALSE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Bayesian marker-regression")
  expect_match(conditionMessage(err), "Kernel-GBLUP")
  expect_match(conditionMessage(err), "Gaussian-Process-GBLUP")
  expect_false(grepl("thus, use GBLUP, RKHS, or GBLUP_BRR", conditionMessage(err), fixed = TRUE))
})
