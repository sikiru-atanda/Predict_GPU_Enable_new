# Cross-validation through model_execute() with mirai workers (each fold is a
# task in a daemon process): the out-of-fold predictions must be the
# sequential ones, belong to the right lines, and track the simulated
# genetic values; the refit of the best model must predict the lines without
# phenotype, and Run_metadata must show the CV's backend.

sim_cv_data <- function(n = 150L, p = 300L, n_qtl = 20L, n_hidden = 20L) {
  set.seed(11)
  freq <- stats::runif(p, 0.1, 0.9)
  ids <- sprintf("L%03d", seq_len(n))
  G <- matrix(2L * stats::rbinom(n * p, 1, rep(freq, each = n)), n, p,
              dimnames = list(ids, sprintf("m%04d", seq_len(p))))
  G <- G[, apply(G, 2, stats::var) > 0, drop = FALSE]
  qtl <- sample(ncol(G), n_qtl)
  g <- as.numeric(scale(G[, qtl, drop = FALSE] %*% stats::rnorm(n_qtl)))
  y <- 10 + g + stats::rnorm(n, sd = sqrt(0.3 / 0.7))
  hidden <- sample(ids, n_hidden)
  pheno <- data.frame(GID = ids, Yield = y, stringsAsFactors = FALSE)
  pheno$Yield[pheno$GID %in% hidden] <- NA
  list(pheno = pheno, geno = G, g = stats::setNames(g, ids), hidden = hidden)
}

run_cv <- function(dat, mode, evaluation_only = TRUE) {
  withr::local_dir(withr::local_tempdir())
  suppressWarnings(PredictProR::model_execute(
    pheno_data = dat$pheno, geno_data = dat$geno, gen_name = "GID", response = "Yield",
    response_family = "gaussian", fixed = ~1, random = ~GID, gmatrix_method = "VanRaden",
    cross_validation = TRUE, cross_validation_meth = "K-Folds", nfolds = 5L, replication = 1,
    random_state = 20261010L, cv_evaluation_only = evaluation_only, cv_generate_plots = FALSE,
    nIter = 600L, burnIn = 100L, thin = 5L, GS_model_cv = c("GBLUP_BRR", "BayesB"),
    num_cores = 2L, parallel_mode = mode, parallel_backend_prefer_fork = FALSE,
    system_database = TRUE, message = FALSE, verbose = FALSE
  ))
}

oof_table <- function(res, model) {
  el <- Filter(function(r) identical(r$model, model), res$cv_results_raw)[[1L]]
  el$ypred_cv_Reps_all
}

test_that("mirai CV predictions equal the sequential ones, track the true genetic values, and Run_metadata shows the CV backend", {
  skip_on_cran()
  skip_if_not_installed("mirai", "2.7.0")
  skip_if_not_installed("purrr")
  skip_if_not_installed("BGLR")

  dat <- sim_cv_data()
  phenotyped <- dat$pheno$GID[!is.na(dat$pheno$Yield)]
  seq_res <- run_cv(dat, "sequential")
  mirai_res <- run_cv(dat, "mirai", evaluation_only = FALSE)

  for (model in c("GBLUP_BRR", "BayesB")) {
    seq_oof <- oof_table(seq_res, model)
    mirai_oof <- oof_table(mirai_res, model)
    expect_identical(mirai_oof, seq_oof, info = model)
    # one out-of-fold prediction per phenotyped line, next to its own phenotype
    oof <- mirai_oof[!is.na(mirai_oof$y), ]
    expect_false(anyDuplicated(mirai_oof$GID) > 0L, info = model)
    expect_setequal(oof$GID, phenotyped)
    expect_identical(oof$y, dat$pheno$Yield[match(oof$GID, dat$pheno$GID)], info = model)
    expect_true(all(is.finite(oof$yhat)), info = model)
    expect_true(all(is.na(mirai_oof$yhat[mirai_oof$GID %in% dat$hidden])), info = model)
    # out-of-fold: clearly above chance, clearly below a fit that saw the line
    expect_gt(stats::cor(oof$yhat, dat$g[oof$GID]), 0.35)
    expect_lt(stats::cor(oof$yhat, oof$y), 0.95)
  }

  # the refit of the best model predicts the lines without phenotype
  pv <- mirai_res$model_results$predicted_values
  hidden_pred <- pv[pv$GID %in% dat$hidden, ]
  expect_setequal(hidden_pred$GID, dat$hidden)
  expect_true(all(is.finite(hidden_pred$Predicted_value)))
  expect_gt(stats::cor(hidden_pred$Predicted_value, dat$g[hidden_pred$GID]), 0.35)

  # Run_metadata: policy_* for the refit, cv_policy_* for the CV on the daemons
  meta <- mirai_res$run_metadata
  val <- function(k) meta$value[meta$key == k]
  expect_identical(val("policy_stage"), "best_model_refit")
  expect_identical(val("cv_policy_backend"), "mirai")
  expect_identical(val("cv_policy_workers"), "2")
  expect_identical(val("cv_policy_decision_reason"), "user_forced_mirai")
  expect_gt(as.integer(val("cv_policy_parallel_count")), 1L)
  expect_identical(val("cv_context"), "cross_validation")
  expect_false(anyDuplicated(meta$key) > 0L)
  cv_meta <- mirai_res$cv_results_processed$run_metadata
  expect_identical(cv_meta$value[cv_meta$key == "policy_backend"], "mirai")
})
