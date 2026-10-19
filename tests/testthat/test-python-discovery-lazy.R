# Finding a Python interpreter starts every candidate on the machine (13 on a
# Windows test machine: 13.8 s per fresh session, most of a small CV). Runs of
# R-only models must not do it, and run metadata never does.

count_python_discovery <- function(env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  local_mocked_bindings(
    gp_python_candidates = function() { calls$n <- calls$n + 1L; character() },
    gp_python_probe = function(...) { calls$n <- calls$n + 1L; NULL },
    .package = "PredictProR",
    .env = env
  )
  withr::local_envvar(PREDICTPRO_GP_PYTHON = NA, RETICULATE_PYTHON = NA, .local_envir = env)
  # a choice cached by an earlier test would hide the discovery
  cache <- PredictProR:::PredictProR_runtime_cache
  keys <- grep("^(preferred_python|ml_python|python_probe)::", ls(cache), value = TRUE)
  saved <- mget(keys, envir = cache)
  rm(list = keys, envir = cache)
  withr::defer(list2env(saved, envir = cache), envir = env)
  calls
}

small_cv_data <- function(n = 60, p = 120) {
  set.seed(11)
  ids <- sprintf("L%03d", seq_len(n))
  G <- matrix(2L * stats::rbinom(n * p, 1, 0.5), n, p, dimnames = list(ids, sprintf("m%03d", seq_len(p))))
  list(pheno = data.frame(GID = ids, Yield = as.numeric(G[, 1:5] %*% stats::rnorm(5)) + stats::rnorm(n)),
       geno = G)
}

test_that("run metadata does not start Python interpreters to describe them", {
  calls <- count_python_discovery()
  md <- PredictProR:::gp_runtime_metadata(context = "test")
  expect_identical(calls$n, 0L)
  expect_true(is.data.frame(md) || is.list(md))
})

test_that("run metadata reports an interpreter set by environment variable without starting it", {
  calls <- count_python_discovery()
  py <- withr::local_tempfile(fileext = ".exe")
  writeLines("", py)
  withr::local_envvar(RETICULATE_PYTHON = py)
  md <- PredictProR:::gp_runtime_metadata(context = "test")
  expect_identical(calls$n, 0L)
  expect_identical(md$value[md$key == "active_python"], py)
})

# Without nvidia-smi, detect_num_gpus() asks torch, which means finding an
# interpreter and starting embedded Python: only GPU-capable models need that.
test_that("the parallel policy for R-only models looks for no Python on a machine without nvidia-smi", {
  calls <- count_python_discovery()
  local_mocked_bindings(detect_physical_num_gpus = function() NA_integer_, .package = "PredictProR")
  res <- PredictProR:::sp_decide_policy(
    models = c("BayesB", "Lasso"), model_params_list = list(list(), list()),
    globals = list(pheno_data = data.frame(y = 1:50)), n_tasks = 2L, user_mode = "auto",
    verbose = FALSE, num_cores = 2L, worker_memory_gb = 0.5,
    task_cost = list(n_train = 50, n_markers = 100, folds_per_task = 1)
  )
  expect_identical(calls$n, 0L)
  expect_identical(res$num_gpus, 0L)
})

test_that("a cross-validation of R-only models looks for no Python", {
  skip_on_cran()
  calls <- count_python_discovery()
  d <- small_cv_data()
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(model_execute(
    pheno_data = d$pheno, geno_data = d$geno, gen_name = "GID", response = "Yield", random = ~GID,
    GS_model_cv = c("BayesB", "Lasso"), cross_validation = TRUE, cross_validation_meth = "K-Folds",
    nfolds = 3L, replication = 1, cv_evaluation_only = TRUE, cv_generate_plots = FALSE,
    nIter = 200L, burnIn = 50L, thin = 2L, qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE, verbose = FALSE
  ))
  expect_true(length(res$cv_results_raw) >= 1L)
  expect_identical(calls$n, 0L)
})

test_that("a true prediction with an R-only model looks for no Python", {
  skip_on_cran()
  calls <- count_python_discovery()
  d <- small_cv_data()
  d$pheno$Yield[1:10] <- NA
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(model_execute(
    pheno_data = d$pheno, geno_data = d$geno, gen_name = "GID", response = "Yield", random = ~GID,
    GS_model = "BRR", nIter = 200L, burnIn = 50L, thin = 2L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE, verbose = FALSE
  ))
  expect_false(is.null(res$run_metadata))
  expect_identical(calls$n, 0L)
})
