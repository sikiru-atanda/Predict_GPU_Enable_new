if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

mk_bayesian_cv_runtime_data <- function(n = 12L, p = 6L) {
  set.seed(20260503)
  X <- matrix(sample(c(0L, 2L), n * p, replace = TRUE), nrow = n)
  for (j in seq_len(ncol(X))) {
    if (length(unique(X[, j])) < 2L) {
      X[, j] <- 2L * ((seq_len(n) + j) %% 2L)
    }
  }
  rownames(X) <- paste0("g", seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))
  pheno <- data.frame(
    GID = rownames(X),
    Yield = drop(0.5 * X[, 1] - 0.2 * X[, 2] + rnorm(n, sd = 0.1)),
    stringsAsFactors = FALSE
  )
  list(pheno = pheno, geno = X)
}

mk_bayesian_met_cv_runtime_data <- function(n_gid = 6L, n_marker = 6L) {
  set.seed(20260504)
  gids <- paste0("g", seq_len(n_gid))
  envs <- paste0("E", seq_len(3L))
  X <- matrix(sample(c(0L, 2L), n_gid * n_marker, replace = TRUE), nrow = n_gid)
  for (j in seq_len(ncol(X))) {
    if (length(unique(X[, j])) < 2L) {
      X[, j] <- 2L * ((seq_len(n_gid) + j) %% 2L)
    }
  }
  rownames(X) <- gids
  colnames(X) <- paste0("m", seq_len(n_marker))
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pheno$Yield <- 1 +
    0.2 * match(pheno$Env, envs) +
    0.5 * X[pheno$GID, 1] +
    rnorm(nrow(pheno), sd = 0.2)
  list(pheno = pheno, geno = X)
}

expect_bayesian_cv_output <- function(out, model, token = NULL) {
  expect_true(is.list(out))
  expect_true("cv_results_raw" %in% names(out))
  expect_gt(length(out$cv_results_raw), 0L)
  raw <- out$cv_results_raw[[1L]]
  expect_identical(raw$model, model)
  if (!is.null(token)) {
    expect_identical(raw$cv_info$token, token)
  }
  expect_s3_class(raw$ypred_cv_Reps_all, "data.frame")
  expect_true(any(raw$ypred_cv_Reps_all$cv_role == "test"))
  expect_true(any(is.finite(as.numeric(raw$ypred_cv_Reps_all$yhat))))
}

test_that("Bayesian single-environment CV runtime covers stratified and unstratified folds", {
  skip_if_not_installed("BGLR")
  dat <- mk_bayesian_cv_runtime_data()
  set.seed(20260907)
  omic1 <- matrix(rnorm(nrow(dat$geno) * 3L), nrow = nrow(dat$geno), ncol = 3L)
  omic2 <- matrix(rnorm(nrow(dat$geno) * 2L), nrow = nrow(dat$geno), ncol = 2L)
  rownames(omic1) <- rownames(omic2) <- rownames(dat$geno)
  colnames(omic1) <- paste0("transcript", seq_len(ncol(omic1)))
  colnames(omic2) <- paste0("metabolite", seq_len(ncol(omic2)))

  bayesa <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    omic1_data = omic1,
    omic2_data = omic2,
    response = "Yield",
    gen_name = "GID",
    GS_model_cv = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Stratified_K-Folds",
    nfolds = 2,
    random_state = 1,
    replication = 1,
    eval_metrics = "mean_squared_error",
    system_database = TRUE,
    message = FALSE,
    nIter = 40L,
    burnIn = 10L,
    thin = 2L
  )
  expect_bayesian_cv_output(bayesa, "BayesA", "stratified_k_folds")

  gblup_brr <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Yield",
    gen_name = "GID",
    GS_model_cv = "GBLUP_BRR",
    gmatrix_method = "VanRaden",
    fixed = ~ 1,
    random = ~ GID,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 2,
    random_state = 1,
    replication = 1,
    eval_metrics = "mean_squared_error",
    system_database = TRUE,
    message = FALSE,
    nIter = 40L,
    burnIn = 10L,
    thin = 2L
  )
  expect_bayesian_cv_output(gblup_brr, "GBLUP_BRR", "k_folds")
})

test_that("Bayesian multi-environment GBLUP_BRR runtime covers CV0 CV1 and CV2", {
  skip_if_not_installed("BGLR")
  dat <- mk_bayesian_met_cv_runtime_data()

  for (cv_method in c("CV0", "CV1", "CV2")) {
    out <- PredictProR::model_execute(
      pheno_data = dat$pheno,
      geno_data = dat$geno,
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      GS_model_cv = "GBLUP_BRR",
      gmatrix_method = "VanRaden",
      fixed = ~ Env,
      random = ~ GID + GID:Env,
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      cross_validation_meth = cv_method,
      nfolds = 2,
      random_state = 1,
      replication = 1,
      eval_metrics = "mean_squared_error",
      system_database = TRUE,
      message = FALSE,
      nIter = 40L,
      burnIn = 10L,
      thin = 2L
    )
    expect_bayesian_cv_output(out, "GBLUP_BRR", tolower(cv_method))
  }
}
)
