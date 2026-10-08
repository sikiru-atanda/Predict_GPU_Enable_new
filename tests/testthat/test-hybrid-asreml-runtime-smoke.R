if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

local_hybrid_asreml_runtime <- function() {
  major_minor <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  candidate_libs <- c(
    Sys.getenv("R_LIBS_USER", unset = ""),
    file.path(Sys.getenv("LOCALAPPDATA", unset = ""), "R", "win-library", major_minor),
    file.path(Sys.getenv("USERPROFILE", unset = ""), "AppData", "Local", "R", "win-library", major_minor)
  )
  candidate_libs <- candidate_libs[nzchar(candidate_libs) & dir.exists(candidate_libs)]
  if (length(candidate_libs)) {
    .libPaths(unique(c(.libPaths(), candidate_libs)))
  }
  if (!requireNamespace("asreml", quietly = TRUE)) {
    skip("{asreml} is not installed")
  }
  status <- tryCatch(asreml::asreml.license.status(), error = function(e) NULL)
  if (is.null(status) || is.null(status$status) || is.na(status$status) || as.integer(status$status) < 0L) {
    msg <- if (!is.null(status$statusMessage) && nzchar(status$statusMessage)) status$statusMessage else "ASReml-R license unavailable"
    skip(msg)
  }
}

mk_hybrid_asreml_runtime_data <- function() {
  set.seed(42)
  parents <- paste0("P", seq_len(6))
  X <- matrix(sample(0:2, length(parents) * 12L, replace = TRUE), nrow = length(parents), ncol = 12L)
  rownames(X) <- parents
  colnames(X) <- paste0("m", seq_len(ncol(X)))

  crosses <- expand.grid(
    Female = parents[1:3],
    Male = parents[4:6],
    stringsAsFactors = FALSE
  )
  crosses <- crosses[c(1, 2, 3, 4, 5, 7, 8, 9), , drop = FALSE]
  female_signal <- c(1.0, 0.3, -0.6)
  male_signal <- c(0.4, -0.2, 0.8)
  sca_signal <- c(0.2, -0.1, 0.0, 0.3, -0.2, 0.15, -0.25, 0.1)
  y <- female_signal[match(crosses$Female, parents[1:3])] +
    male_signal[match(crosses$Male, parents[4:6])] +
    sca_signal +
    rnorm(nrow(crosses), sd = 0.05)
  y[c(2, 6)] <- NA

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
    female_geno = X[parents[1:3], , drop = FALSE],
    male_geno = X[parents[4:6], , drop = FALSE]
  )
}

mk_hybrid_asreml_cv_data <- function() {
  set.seed(77)
  female <- paste0("F", 1:3)
  male <- paste0("M", 1:3)
  parents <- c(female, male)
  X <- matrix(sample(0:2, length(parents) * 15L, replace = TRUE), nrow = length(parents), ncol = 15L)
  rownames(X) <- parents
  colnames(X) <- paste0("m", seq_len(ncol(X)))

  crosses <- expand.grid(
    Female = female,
    Male = male,
    stringsAsFactors = FALSE
  )
  female_signal <- c(0.9, 0.2, -0.5)
  male_signal <- c(0.4, -0.1, 0.7)
  sca_signal <- c(0.15, -0.05, 0.10, 0.20, -0.10, 0.05, -0.20, 0.12, -0.08)
  y <- female_signal[match(crosses$Female, female)] +
    male_signal[match(crosses$Male, male)] +
    sca_signal +
    rnorm(nrow(crosses), sd = 0.03)

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

mk_hybrid_asreml_met_runtime_data <- function() {
  dat <- mk_hybrid_asreml_runtime_data()
  obs <- subset(dat$pheno, !is.na(Yield))
  ph1 <- obs
  ph1$Env <- "E1"
  ph2 <- obs
  ph2$Env <- "E2"
  ph2$Yield <- ph2$Yield + 0.5
  ph <- rbind(ph1, ph2)
  ph$Yield[c(2, 11)] <- NA_real_
  list(
    pheno = ph,
    geno = dat$geno,
    female_geno = dat$female_geno,
    male_geno = dat$male_geno
  )
}

test_that("hybrid gaussian ASReml-R runtime returns total and component outputs", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_identical(names(preds), PredictProR:::gp_hybrid_prediction_columns())
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_true(is.data.frame(out$model_results$variance_components))
  expect_identical(
    names(out$model_results$diagnostic_plots),
    c("predicted_vs_observed", "predicted_vs_reliability", "prediction_interval", "reliability_summary")
  )
})

test_that("hybrid gaussian ASReml-R runtime accepts separate female and male parent genotype tables", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_identical(names(preds), PredictProR:::gp_hybrid_prediction_columns())
})

test_that("hybrid ASReml-R requires parent columns", {
  dat <- mk_hybrid_asreml_runtime_data()

  expect_error(
    suppressWarnings(PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      gmatrix_method = "VanRaden",
      response = "Yield",
      gen_name = "HybridID",
      GS_model = "GBLUP",
      engine = "asreml",
      response_family = "gaussian",
      hybrid_asreml = TRUE,
      system_database = TRUE,
      message = FALSE
    )),
    "Provide valid female_parent and male_parent"
  )
})

test_that("hybrid known-parents CV returns scenario summaries", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_cv_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model_cv = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
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
})

test_that("hybrid one-new-parent CV accepts separate parent genotype tables", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_cv_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    female_geno_data = dat$female_geno,
    male_geno_data = dat$male_geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model_cv = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
    cross_validation = TRUE,
    cross_validation_meth = "Hybrid_One_New_Parent",
    nfolds = 3,
    replication = 1,
    eval_metrics = c("mean_absolute_error", "root_mean_squared_error"),
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$cv_results_processed$hybrid_cv_predictions
  expect_true(nrow(preds) > 0)
  expect_true(all(preds$cv_scenario == "One_New_Parent"))
})

test_that("hybrid both-new-parents CV returns both-new-parent rows", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_cv_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    GS_model_cv = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
    cross_validation = TRUE,
    cross_validation_meth = "Hybrid_Both_New_Parents",
    nfolds = 3,
    replication = 1,
    eval_metrics = c("mean_absolute_error", "root_mean_squared_error"),
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$cv_results_processed$hybrid_cv_predictions
  expect_true(nrow(preds) > 0)
  expect_true(all(preds$cv_scenario == "Both_New_Parents"))
})

test_that("hybrid ASReml-R supports multi-environment rows", {
  local_hybrid_asreml_runtime()

  dat <- mk_hybrid_asreml_met_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gmatrix_method = "VanRaden",
    response = "Yield",
    gen_name = "HybridID",
    heter_groups = "Env",
    GS_model = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    hybrid_asreml = TRUE,
    female_parent = "Female",
    male_parent = "Male",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  expect_equal(nrow(preds), nrow(dat$pheno))
  expect_identical(names(preds), PredictProR:::gp_hybrid_prediction_columns(include_env = TRUE))
  expect_true(all(c("E1", "E2") %in% unique(preds$Env)))
})
