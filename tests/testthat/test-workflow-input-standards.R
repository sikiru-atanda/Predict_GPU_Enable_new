multi_trait_standard_call <- function(...) {
  if (exists("multi_trait_data_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("multi_trait_data_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::multi_trait_data_standard(...))
  }
  stop("multi_trait_data_standard is unavailable for this test.", call. = FALSE)
}

multi_trait_validate_call <- function(...) {
  if (exists("validate_multi_trait_input_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("validate_multi_trait_input_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::validate_multi_trait_input_standard(...))
  }
  stop("validate_multi_trait_input_standard is unavailable for this test.", call. = FALSE)
}

met_standard_call <- function(...) {
  if (exists("met_data_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("met_data_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::met_data_standard(...))
  }
  stop("met_data_standard is unavailable for this test.", call. = FALSE)
}

met_validate_call <- function(...) {
  if (exists("validate_met_input_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("validate_met_input_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::validate_met_input_standard(...))
  }
  stop("validate_met_input_standard is unavailable for this test.", call. = FALSE)
}

workflow_model_execute_call <- function(...) {
  if (exists("model_execute", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("model_execute", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::model_execute(...))
  }
  stop("model_execute is unavailable for this test.", call. = FALSE)
}

mk_multitrait_contract_data <- function() {
  ph <- data.frame(
    GID = paste0("g", 1:4),
    Trait1 = c(1.0, 2.0, 3.0, NA),
    Trait2 = c(1.5, NA, 2.5, 4.0),
    stringsAsFactors = FALSE
  )
  X <- matrix(c(
    0, 1, 2,
    1, 1, 0,
    2, 0, 1,
    1, 2, 1
  ), nrow = 4, byrow = TRUE)
  rownames(X) <- ph$GID
  colnames(X) <- paste0("m", 1:3)
  list(ph = ph, X = X)
}

mk_met_contract_data <- function() {
  ph <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1.2, 1.4, 2.1, 2.3),
    stringsAsFactors = FALSE
  )
  X <- matrix(c(
    0, 1,
    1, 0
  ), nrow = 2, byrow = TRUE)
  rownames(X) <- c("g1", "g2")
  colnames(X) <- c("m1", "m2")
  list(ph = ph, X = X)
}

test_that("multi_trait_data_standard returns the package contract", {
  spec <- multi_trait_standard_call()
  expect_true(is.list(spec))
  expect_true(all(c("phenotype", "genomics", "current_scope") %in% names(spec)))
})

test_that("validate_multi_trait_input_standard accepts a valid ML contract", {
  dat <- mk_multitrait_contract_data()
  out <- multi_trait_validate_call(
    pheno_data = dat$ph,
    geno_data = dat$X,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    response_family = "gaussian",
    multi_trait_ml = TRUE
  )
  expect_true(isTRUE(out$valid))
  expect_identical(out$mode, "multi_trait_ml")
})

test_that("validate_multi_trait_input_standard rejects too few responses with standard guidance", {
  dat <- mk_multitrait_contract_data()
  expect_error(
    multi_trait_validate_call(
      pheno_data = dat$ph,
      geno_data = dat$X,
      response = "Trait1",
      gen_name = "GID",
      response_family = "gaussian",
      multi_trait_ml = TRUE
    ),
    regexp = "Multi-trait phenotype standard"
  )
})

test_that("model_execute fails early for malformed multi-trait genomic input", {
  dat <- mk_multitrait_contract_data()
  bad_x <- dat$X
  rownames(bad_x) <- paste0("bad_", rownames(bad_x))
  expect_error(
    workflow_model_execute_call(
      pheno_data = dat$ph,
      geno_data = bad_x,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      response_family = "gaussian",
      multi_trait_ml = TRUE,
      GS_model = "Ridge_Regression",
      system_database = TRUE,
      message = FALSE
    ),
    regexp = "Current genomic standards"
  )
})

test_that("met_data_standard returns the package contract", {
  spec <- met_standard_call()
  expect_true(is.list(spec))
  expect_true(all(c("phenotype", "genomics", "cv_scenarios", "current_models", "current_scope") %in% names(spec)))
  expect_identical(
    spec$current_scope$gp_models,
    c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP")
  )
  expect_true(all(c("GBLUP_BRR", "RKHS") %in% spec$current_scope$kernel_bayesian_models))
  expect_identical(
    names(spec$current_scope$rkhs_environment_modes),
    c("heterogeneous_environment", "weighted_common_variance", "incompatible_combination")
  )
  expect_match(
    spec$current_scope$rkhs_environment_modes[["heterogeneous_environment"]],
    "weights = NULL",
    fixed = TRUE
  )
  expect_match(
    spec$current_scope$rkhs_environment_modes[["incompatible_combination"]],
    "rejected before fitting",
    fixed = TRUE
  )
  expect_true(all(c("Xgboost", "DenseNeuralNet", "TabTransformer") %in%
                    spec$current_scope$predictive_ml_dl_models))
  expect_false("Scalable-GBLUP" %in% spec$current_models)
})

test_that("validate_met_input_standard accepts a valid MET contract", {
  dat <- mk_met_contract_data()
  out <- met_validate_call(
    pheno_data = dat$ph,
    geno_data = dat$X,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    met_ml_dl = TRUE,
    GS_model = "RandomForest"
  )
  expect_true(isTRUE(out$valid))
})

test_that("validate_met_input_standard rejects missing heter_groups with standard guidance", {
  dat <- mk_met_contract_data()
  expect_error(
    met_validate_call(
      pheno_data = dat$ph,
      geno_data = dat$X,
      response = "Yield",
      gen_name = "GID",
      met_ml_dl = TRUE,
      GS_model = "RandomForest"
    ),
    regexp = "MET phenotype standard"
  )
})

test_that("validate_met_input_standard rejects a mixed supported and unsupported MET request", {
  dat <- mk_met_contract_data()
  expect_error(
    met_validate_call(
      pheno_data = dat$ph,
      geno_data = dat$X,
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      met_ml_dl = TRUE,
      GS_model = c("Xgboost", "Lasso")
    ),
    regexp = "Lasso The request was stopped before fitting"
  )
})

test_that("model_execute fails early for malformed MET input", {
  dat <- mk_met_contract_data()
  expect_error(
    workflow_model_execute_call(
      pheno_data = dat$ph,
      geno_data = dat$X,
      response = "Yield",
      gen_name = "GID",
      response_family = "gaussian",
      met_ml_dl = TRUE,
      GS_model = "RandomForest",
      system_database = TRUE,
      message = FALSE
    ),
    regexp = "MET phenotype standard"
  )
})
