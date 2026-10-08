if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

local_asreml_runtime <- function() {
  skip_if_not_installed("asreml")
  status <- tryCatch(asreml::asreml.license.status(), error = function(e) NULL)
  if (is.null(status) || is.null(status$status) || is.na(status$status) || as.integer(status$status) < 0L) {
    msg <- if (!is.null(status$statusMessage) && nzchar(status$statusMessage)) status$statusMessage else "ASReml-R license unavailable"
    skip(msg)
  }
}

mk_multitrait_asreml_runtime_data <- function(n = 18L, p = 8L) {
  set.seed(123)
  gids <- paste0("g", seq_len(n))
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- gids
  colnames(X) <- paste0("m", seq_len(p))
  gmat <- tcrossprod(scale(X))
  gmat <- gmat / mean(diag(gmat), na.rm = TRUE)
  diag(gmat) <- diag(gmat) + 0.1
  base <- X[, 1] * 0.8 - X[, 2] * 0.4
  t1 <- base + rnorm(n, sd = 0.2)
  t2 <- base * 0.5 + X[, 3] * 0.5 + rnorm(n, sd = 0.2)
  t3 <- -base * 0.4 + X[, 4] * 0.6 + rnorm(n, sd = 0.2)
  t1[c(2, 7, 15)] <- NA
  t2[c(4, 9, 12)] <- NA
  t3[c(3, 11, 16)] <- NA
  list(
    pheno = data.frame(
      GID = gids,
      Trait1 = t1,
      Trait2 = t2,
      Trait3 = t3,
      stringsAsFactors = FALSE
    ),
    gmatrix = gmat
  )
}

test_that("multi-trait gaussian ASReml-R runtime returns long and wide outputs", {
  local_asreml_runtime()

  dat <- mk_multitrait_asreml_runtime_data()

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    gmatrix = dat$gmatrix,
    response = c("Trait1", "Trait2", "Trait3"),
    gen_name = "GID",
    GS_model = "GBLUP",
    engine = "asreml",
    response_family = "gaussian",
    multi_trait_asreml = TRUE,
    var_cov_str = "diag",
    system_database = TRUE,
    message = FALSE
  )

  preds <- out$model_results$predicted_values
  wide <- out$model_results$multitrait_prediction_wide
  counts <- out$model_results$multitrait_trait_counts

  expected_columns <- PredictProR::prediction_output_standard(
    response_family = "gaussian",
    task = "multi_trait"
  )$prediction_columns
  expect_identical(names(preds), expected_columns)
  expect_true(PredictProR::validate_prediction_output(out$model_results)$valid)
  expect_equal(sort(unique(preds$Trait)), c("Trait1", "Trait2", "Trait3"))
  expect_true(all(c("Train", "Test") %in% unique(preds$Train_Test_Label)))
  expect_equal(nrow(wide), nrow(dat$pheno))
  expect_true(all(c("GID", "Trait1", "Trait2", "Trait3") %in% names(wide)))
  expect_equal(sort(counts$trait), c("Trait1", "Trait2", "Trait3"))
  expect_true(all(counts$predicted_count > 0))
  expect_true(is.list(out$model_results$diagnostic_plots))
  expect_true(any(vapply(
    out$model_results$diagnostic_plots,
    inherits,
    logical(1L),
    what = "ggplot"
  )))
})

test_that("multi-trait ASReml-R rejects non-gaussian requests", {
  dat <- mk_multitrait_asreml_runtime_data()

  expect_error(
    suppressWarnings(PredictProR::model_execute(
      pheno_data = dat$pheno,
      gmatrix = dat$gmatrix,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      GS_model = "GBLUP",
      engine = "asreml",
      response_family = "binary",
      multi_trait_asreml = TRUE,
      system_database = TRUE,
      message = FALSE
    )),
    "support gaussian traits only"
  )
})
