if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

mk_bayes_runtime_data <- function(kind = c("binary", "ordinal"), n = 24L, p = 8L) {
  kind <- match.arg(kind)
  set.seed(20260426 + match(kind, c("binary", "ordinal")))
  X <- matrix(
    sample(0:2, n * p, replace = TRUE, prob = c(0.4, 0.4, 0.2)),
    nrow = n
  )
  rownames(X) <- paste0(substr(kind, 1, 1), seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))

  y <- switch(
    kind,
    binary = {
      eta <- 1.0 * X[, 1] - 0.7 * X[, 2] + 0.4 * X[, 3]
      factor(ifelse(runif(n) < plogis(eta), "yes", "no"), levels = c("no", "yes"))
    },
    ordinal = {
      score <- X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3] + rnorm(n, sd = 0.25)
      cuts <- stats::quantile(score, probs = c(1 / 3, 2 / 3))
      ordered(
        ifelse(score <= cuts[[1]], "low", ifelse(score <= cuts[[2]], "medium", "high")),
        levels = c("low", "medium", "high")
      )
    }
  )

  pheno <- data.frame(GID = rownames(X), Trait = y, stringsAsFactors = FALSE)
  tst <- sample(seq_len(n), ceiling(n * 0.25))
  pheno$Trait[tst] <- NA
  list(pheno = pheno, geno = X, tst = tst)
}

extract_single_trait_result <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.null(x$model_results)) return(x)
  x[[1L]]
}

extract_bayes_predicted_value <- function(res) {
  if (is.null(res) || is.null(res$model_results)) {
    return(NULL)
  }
  res$model_results$predicted_values %||% res$model_results$Predicted_value
}

find_summary_statistics <- function(res) {
  if (is.null(res) || is.null(res$summary_statistic)) {
    return(NULL)
  }
  res$summary_statistic$summary_statistics
}

probability_columns <- function(x) {
  grep("^(Prob|Probability)_", names(x), value = TRUE)
}

expect_bayesian_classification_contract <- function(pred, class_levels) {
  prob_cols <- paste0("Probability_", make.names(class_levels, unique = TRUE))
  expect_identical(
    names(pred),
    c(PredictProR:::gp_public_prediction_columns("ordinal"), prob_cols)
  )
  expect_true(all(pred$Train_Test_Label %in% c("Train", "Test")))
  expect_true(all(is.finite(pred$Prediction_confidence)))
  expect_true(all(is.finite(pred$Classification_uncertainty)))
  expect_true(all(is.finite(pred$Reliability)))
  expect_true(all(!is.na(pred$Reliability_remarks) & nzchar(pred$Reliability_remarks)))
  prob <- as.matrix(pred[, prob_cols, drop = FALSE])
  expect_true(all(is.finite(prob)))
  expect_true(all(prob >= 0 & prob <= 1))
  expect_equal(rowSums(prob), rep(1, nrow(pred)), tolerance = 1e-6)
  confidence <- apply(prob, 1L, max)
  expect_equal(pred$Prediction_confidence, confidence, tolerance = 1e-6)
  expect_equal(pred$Classification_uncertainty, 1 - confidence, tolerance = 1e-6)
  expect_equal(pred$Reliability, confidence, tolerance = 1e-6)
  expect_identical(
    pred$Predicted_class,
    class_levels[max.col(prob, ties.method = "first")]
  )
}

find_new_export_dir <- function(root, before_dirs) {
  after_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  new_dirs <- setdiff(after_dirs, before_dirs)
  if (!length(new_dirs)) {
    return(NULL)
  }
  new_dirs[[which.max(file.info(new_dirs)$mtime)]]
}

find_predicted_value_file <- function(export_dir) {
  files <- list.files(
    export_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  pred_files <- files[grepl("^predicted(_value|_values)?\\.csv$", basename(files), ignore.case = TRUE)]
  if (!length(pred_files)) {
    return(NULL)
  }
  pred_files[[1L]]
}

test_that("Bayesian binary runtime smoke returns posterior classification outputs", {
  skip_if_not_installed("BGLR")

  dat <- mk_bayes_runtime_data("binary", n = 20L, p = 6L)

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Trait",
    gen_name = "GID",
    GS_model = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    ld_pruning = FALSE,
    response_family = "binary",
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )

  res <- extract_single_trait_result(out)
  pred <- extract_bayes_predicted_value(res)
  prob_cols <- probability_columns(pred)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(dat$pheno))
  expect_true(all(c(
    "GID", "Predicted_class", "Train_Test_Label",
    "Prediction_confidence", "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(prob_cols), 2L)
  expect_true(any(pred$Train_Test_Label == "Test"))
  expect_bayesian_classification_contract(pred, c("no", "yes"))

  summ <- find_summary_statistics(res)
  expect_s3_class(summ, "data.frame")
  expect_true("Response_Family" %in% summ$stat)
  expect_true("Prediction_Confidence_Definition" %in% summ$stat)
})

test_that("Bayesian ordinal runtime smoke returns posterior classification outputs", {
  skip_if_not_installed("BGLR")

  dat <- mk_bayes_runtime_data("ordinal", n = 21L, p = 6L)

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    response = "Trait",
    gen_name = "GID",
    GS_model = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    ld_pruning = FALSE,
    response_family = "ordinal",
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )

  res <- extract_single_trait_result(out)
  pred <- extract_bayes_predicted_value(res)
  prob_cols <- probability_columns(pred)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(dat$pheno))
  expect_true(all(c(
    "GID", "Predicted_class", "Train_Test_Label",
    "Prediction_confidence", "Classification_uncertainty"
  ) %in% names(pred)))
  expect_gte(length(prob_cols), 3L)
  expect_true(any(pred$Train_Test_Label == "Test"))
  expect_bayesian_classification_contract(pred, c("low", "medium", "high"))

  summ <- find_summary_statistics(res)
  expect_s3_class(summ, "data.frame")
  expect_true("Response_Family" %in% summ$stat)
  expect_true("Prediction_Confidence_Definition" %in% summ$stat)
})

test_that("Bayes Alphabet multi-omics classification reports liability variance", {
  skip_if_not_installed("BGLR")

  dat <- mk_bayes_runtime_data("binary", n = 20L, p = 6L)
  set.seed(20260907)
  omic1 <- matrix(rnorm(nrow(dat$geno) * 4L), nrow = nrow(dat$geno), ncol = 4L)
  rownames(omic1) <- rownames(dat$geno)
  colnames(omic1) <- paste0("omic", seq_len(ncol(omic1)))

  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    omic1_data = omic1,
    response = "Trait",
    gen_name = "GID",
    GS_model = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    ld_pruning = FALSE,
    response_family = "binary",
    nIter = 120L,
    burnIn = 20L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )

  res <- extract_single_trait_result(out)
  vc <- res$model_results$variance_components %||%
    res$model_results$Variance_components
  expect_s3_class(vc, "data.frame")
  expect_true(all(c(
    "liability_scale_genetic_variance_geno_data",
    "liability_scale_genetic_variance_omic1_data",
    "liability_scale_total_genetic_variance",
    "liability_scale_residual_variance_fixed",
    "liability_scale_heritability"
  ) %in% vc$Component))
  expect_true(all(is.finite(vc$Components)))
  expect_equal(
    vc$Components[vc$Component == "liability_scale_residual_variance_fixed"],
    1
  )
})

test_that("Bayesian classification exports write posterior probability columns", {
  skip_if_not_installed("BGLR")

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-bayes-class-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  before_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  bdat <- mk_bayes_runtime_data("binary", n = 18L, p = 6L)
  out_bin <- PredictProR::model_execute(
    pheno_data = bdat$pheno,
    geno_data = bdat$geno,
    response = "Trait",
    gen_name = "GID",
    GS_model = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    ld_pruning = FALSE,
    response_family = "binary",
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    system_database = FALSE,
    message = FALSE
  )
  expect_true(length(out_bin) >= 1L)
  bin_dir <- find_new_export_dir(tmp_dir, before_dirs)
  expect_false(is.null(bin_dir))
  pred_bin_file <- find_predicted_value_file(bin_dir)
  expect_false(is.null(pred_bin_file))
  pred_bin <- utils::read.csv(pred_bin_file, stringsAsFactors = FALSE)
  expect_true("Prediction_confidence" %in% names(pred_bin))
  expect_gte(length(probability_columns(pred_bin)), 2L)
  expect_bayesian_classification_contract(pred_bin, c("no", "yes"))

  before_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  odat <- mk_bayes_runtime_data("ordinal", n = 18L, p = 6L)
  out_ord <- PredictProR::model_execute(
    pheno_data = odat$pheno,
    geno_data = odat$geno,
    response = "Trait",
    gen_name = "GID",
    GS_model = "BayesA",
    fixed = ~ 1,
    random = ~ GID,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    ld_pruning = FALSE,
    response_family = "ordinal",
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    system_database = FALSE,
    message = FALSE
  )
  expect_true(length(out_ord) >= 1L)
  ord_dir <- find_new_export_dir(tmp_dir, before_dirs)
  expect_false(is.null(ord_dir))
  pred_ord_file <- find_predicted_value_file(ord_dir)
  expect_false(is.null(pred_ord_file))
  pred_ord <- utils::read.csv(pred_ord_file, stringsAsFactors = FALSE)
  expect_true("Prediction_confidence" %in% names(pred_ord))
  expect_gte(length(probability_columns(pred_ord)), 3L)
  expect_bayesian_classification_contract(pred_ord, c("low", "medium", "high"))
})
