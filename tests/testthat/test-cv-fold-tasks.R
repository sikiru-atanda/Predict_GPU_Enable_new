mk_fold_task_data <- function(n = 60L, p = 120L, envs = NULL) {
  set.seed(42)
  ids <- sprintf("L%02d", seq_len(n))
  X <- matrix(sample(c(0, 2), n * p, replace = TRUE), n, p, dimnames = list(ids, sprintf("m%03d", seq_len(p))))
  u <- as.numeric(scale(X[, 1:10]) %*% rnorm(10))
  if (is.null(envs)) {
    ph <- data.frame(GID = ids, Yield = 5 + u + rnorm(n, sd = 0.5), stringsAsFactors = FALSE)
  } else {
    ph <- expand.grid(GID = ids, Env = envs, stringsAsFactors = FALSE)
    ph$Yield <- 5 + u[match(ph$GID, ids)] + rnorm(nrow(ph), sd = 0.5)
  }
  list(pheno = ph, geno = X)
}

run_fold_task_cv <- function(dat, mode, ..., models = c("BRR", "GBLUP_BRR"), random = ~GID) {
  withr::local_dir(withr::local_tempdir())  # exported results stay out of tests/
  suppressWarnings(PredictProR::model_execute(
    pheno_data = dat$pheno, geno_data = dat$geno, response = "Yield", gen_name = "GID",
    fixed = ~1, random = random, gmatrix_method = "VanRaden",
    GS_model_cv = models, cross_validation = TRUE, cv_evaluation_only = TRUE,
    nIter = 300L, burnIn = 100L, thin = 5L, eval_metrics = c("accuracy", "root_mean_squared_error"),
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = mode, system_database = FALSE, message = FALSE, verbose = FALSE, ...
  ))
}

oof_table <- function(out) {
  do.call(rbind, lapply(out$cv_results_raw, function(r) {
    data.frame(model = r$model, rep = r$rep, row = r$ypred_cv_Reps_all$row_id,
               yhat = r$ypred_cv_Reps_all$yhat, role = r$ypred_cv_Reps_all$cv_role)
  }))
}

test_that("fold-level CV tasks reproduce unsplit CV exactly (K-Folds)", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 3L, "needs more than one usable core for fold tasks")
  withr::local_envvar(GP_CV_FOLD_TASKS = "true", GP_PAR_WORKER_OVERHEAD_GB = "0.5")
  dat <- mk_fold_task_data()
  seq_out <- run_fold_task_cv(dat, "sequential", cross_validation_meth = "K-Folds", nfolds = 3L, replication = 2L)
  # "auto" splits into fold tasks; on this tiny data the cost model runs them
  # in-process, so this checks the split + merge path itself.
  auto_out <- run_fold_task_cv(dat, "auto", cross_validation_meth = "K-Folds", nfolds = 3L, replication = 2L)
  expect_identical(length(auto_out$cv_results_raw), length(seq_out$cv_results_raw))
  a <- oof_table(seq_out); b <- oof_table(auto_out)
  m <- merge(a, b, by = c("model", "rep", "row"))
  expect_equal(nrow(m), nrow(a))
  expect_identical(m$role.x, m$role.y)
  expect_equal(m$yhat.x, m$yhat.y, tolerance = 1e-10)
  agg_a <- seq_out$cv_results_processed$aggregated_data_list$aggregated_across_reps
  agg_b <- auto_out$cv_results_processed$aggregated_data_list$aggregated_across_reps
  agg_b <- agg_b[match(paste(agg_a$model, agg_a$Rep), paste(agg_b$model, agg_b$Rep)), ]
  expect_equal(agg_b$accuracy, agg_a$accuracy, tolerance = 1e-10)
  meta <- auto_out$cv_results_processed$run_metadata
  expect_match(meta$value[meta$key == "policy_decision_reason"], "cost_model_")
})

test_that("fold-level CV tasks reproduce unsplit MET CV1 exactly", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 3L, "needs more than one usable core for fold tasks")
  withr::local_envvar(GP_CV_FOLD_TASKS = "true", GP_PAR_WORKER_OVERHEAD_GB = "0.5")
  dat <- mk_fold_task_data(n = 40L, envs = c("E2", "E1"))
  args <- list(heter_groups = "Env", cross_validation_meth = "CV1", nfolds = 3L, replication = 1L,
               random = ~ GID + GID:Env, models = c("GBLUP_BRR", "RKHS"), kernel_method = "Gaussian")
  seq_out <- do.call(run_fold_task_cv, c(list(dat, "sequential"), args))
  auto_out <- do.call(run_fold_task_cv, c(list(dat, "auto"), args))
  m <- merge(oof_table(seq_out), oof_table(auto_out), by = c("model", "rep", "row"))
  expect_equal(m$yhat.x, m$yhat.y, tolerance = 1e-10)
})

test_that("the cost model keeps small jobs sequential and parallelises heavy ones", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    gp_parallel_memory_budget_gb = function() 32,
    .package = "PredictProR"
  )
  withr::local_envvar(GP_PAR_WORKER_OVERHEAD_GB = "0.5")
  decide <- function(models, cost) {
    PredictProR:::sp_decide_policy(
      models = models, model_params_list = rep(list(list()), length(models)),
      globals = list(params = list(x = matrix(0, 10, 10))), n_tasks = length(models),
      user_mode = "auto", num_cores = 8L, verbose = FALSE, sys_name = "Windows", task_cost = cost
    )
  }
  heavy <- list(n_train = 225, n_markers = 8068, nIter = 3000, ntree = 300, iterations = 300, folds_per_task = 1)
  d <- decide(rep(c("GBLUP_BRR", "BayesB", "RandomForest", "Ridge_Regression"), 10), heavy)
  expect_identical(d$backend, "base_parallel")
  expect_match(d$decision_reason, "^cost_model_parallel")
  tiny <- list(n_train = 50, n_markers = 200, nIter = 500, ntree = 50, iterations = 50, folds_per_task = 1)
  d2 <- decide(rep(c("BRR", "GBLUP_BRR"), 3), tiny)
  expect_identical(d2$backend, "sequential")
  expect_match(d2$decision_reason, "^cost_model_sequential")
  withr::local_envvar(GP_PAR_COST_MODEL = "false")
  d3 <- decide(rep(c("GBLUP_BRR", "BayesB", "RandomForest", "Ridge_Regression"), 10), heavy)
  expect_false(grepl("cost_model", d3$decision_reason))
})
