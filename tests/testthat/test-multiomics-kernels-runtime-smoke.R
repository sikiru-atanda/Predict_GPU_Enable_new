if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

toy_multiomics_data <- function(n = 18L) {
  safe_set_seed(20260517)
  ids <- sprintf("g%02d", seq_len(n))
  make_block <- function(prefix, p) {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    rownames(x) <- ids
    colnames(x) <- paste0(prefix, seq_len(p))
    x
  }
  omic1 <- make_block("rna", 5L)
  omic2 <- make_block("met", 4L)
  omic3 <- make_block("prot", 3L)
  y <- 2 + 0.8 * omic1[, 1] - 0.5 * omic2[, 2] + 0.3 * omic3[, 1] +
    stats::rnorm(n, sd = 0.05)
  y[c(4L, 11L, 17L)] <- NA_real_
  pheno <- data.frame(GID = ids, Yield = y, stringsAsFactors = FALSE)
  list(pheno = pheno, omic1 = omic1, omic2 = omic2, omic3 = omic3)
}

toy_kernel_data <- function(n = 16L) {
  safe_set_seed(20260518)
  ids <- sprintf("k%02d", seq_len(n))
  make_kernel <- function(prefix, p) {
    z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    k <- tcrossprod(scale(z, center = TRUE, scale = FALSE))
    k <- k + diag(0.05, n)
    k <- k / mean(diag(k))
    rownames(k) <- ids
    colnames(k) <- ids
    k
  }
  k0 <- make_kernel("g", 6L)
  k1 <- make_kernel("o1", 5L)
  k2 <- make_kernel("o2", 4L)
  k3 <- make_kernel("o3", 3L)
  y <- 1.5 + stats::rnorm(n, sd = 0.4)
  y[c(3L, 9L, 14L)] <- NA_real_
  pheno <- data.frame(GID = ids, Yield = y, stringsAsFactors = FALSE)
  list(pheno = pheno, gmatrix = k0, omic1_kernel = k1, omic2_kernel = k2, omic3_kernel = k3)
}

safe_set_seed <- function(seed) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) &&
      !is.integer(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed)
}

find_prediction_table <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.list(x)) {
    if (is.data.frame(x[["predicted_values"]])) {
      return(x[["predicted_values"]])
    }
    if (is.data.frame(x[["Predicted_value"]])) {
      return(x[["Predicted_value"]])
    }
    for (nm in names(x)) {
      out <- find_prediction_table(x[[nm]])
      if (is.data.frame(out)) {
        return(out)
      }
    }
  }
  NULL
}

expect_gaussian_public_contract <- function(pred, test_n) {
  expect_s3_class(pred, "data.frame")
  required_columns <- c(
    "GID", "Predicted_value", "Train_Test_Label", "Observed_value",
    "Standard_error", "PEV", "lower_bound", "upper_bound",
    "Uncertainty", "Uncertainty_remarks", "Reliability", "Reliability_remarks"
  )
  expect_true(all(required_columns %in% names(pred)))
  expect_equal(sum(pred$Train_Test_Label == "Test"), test_n)
  expect_true(any(is.finite(pred$Predicted_value)))
}

test_that("user API handles multiple omics feature matrices", {
  dat <- toy_multiomics_data()
  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    omic1_data = dat$omic1,
    omic2_data = dat$omic2,
    omic3_data = dat$omic3,
    response = "Yield",
    gen_name = "GID",
    random = ~ GID,
    GS_model = "Ridge_Regression",
    response_family = "gaussian",
    para_tunning = FALSE,
    n_bootstrap = 2L,
    system_database = TRUE,
    message = FALSE
  )
  pred <- find_prediction_table(out)
  expect_gaussian_public_contract(pred, sum(is.na(dat$pheno$Yield)))
})

test_that("user API handles multiple user-provided kernels", {
  dat <- toy_kernel_data()
  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    gmatrix = dat$gmatrix,
    omic1_kernel = dat$omic1_kernel,
    omic2_kernel = dat$omic2_kernel,
    omic3_kernel = dat$omic3_kernel,
    response = "Yield",
    gen_name = "GID",
    random = ~ GID,
    GS_model = "Kernel-GBLUP",
    response_family = "gaussian",
    para_tunning = FALSE,
    system_database = TRUE,
    message = FALSE,
    gp_full_vc = TRUE,
    gp_return_se = TRUE
  )
  pred <- find_prediction_table(out)
  expect_gaussian_public_contract(pred, sum(is.na(dat$pheno$Yield)))
})
