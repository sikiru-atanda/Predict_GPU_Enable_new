# A model_execute() run, and the helpers it calls, leave the caller's session
# as they found it. KNN imputation and the PLINK-engine VCF QC set a future
# plan for their workers and then set "sequential", and model_execute() set
# "sequential" when it ended, so a plan the user had chosen was lost. The
# single-trait ASReml routes put their vm() inverses (gmatrix_inv, omic1_inv
# ...) in the global environment and left them there, over a user's objects
# of the same names. Their fits also attached asreml (and Matrix with it) to
# the search path, and kept a copy of every fit's data in a package store for
# the rest of the session, under a key drawn from the caller's random-number
# stream. A cross-validation whose parallel tasks the policy ran on the
# sequential backend pinned BLAS/OpenMP to one thread in the caller's session
# (OMP_NUM_THREADS = 1 and the rest), and later Python fits inherited it. A
# true-prediction run left the random_state option at its task's seed, and
# the ML seeds of later runs came from it.

# The user's plan is multicore: setting it starts no process (on Windows its
# futures run in the session). Every future version records the call that set
# a plan, which tells the user's plan apart from one a helper sets (multicore
# too, on Unix).
.plan_call <- function() paste(deparse(attr(future::plan(), "call")), collapse = "")
.local_user_plan <- function(env = parent.frame()) {
  old <- future::plan(future::multicore, workers = 2L)
  withr::defer(future::plan(old), envir = env)
  call <- .plan_call()
  stopifnot(grepl("multicore", call, fixed = TRUE))
  call
}

test_that("KNN imputation puts the caller's future plan back", {
  user_plan <- .local_user_plan()
  withr::local_envvar(PREDICTPROR_KNN_IMPUTE_BACKEND = "r")
  parallel_plan <- logical()
  local_mocked_bindings(
    impute_chunk = function(chunk_data, k = 5) { chunk_data[is.na(chunk_data)] <- 0; chunk_data },
    gp_future_lapply_session = function(X, FUN, ...) {
      parallel_plan <<- c(parallel_plan, inherits(future::plan(), "multiprocess"))
      lapply(X, FUN)
    },
    .package = "PredictProR"
  )
  x <- matrix(c(0, 1, NA, 2, 1, 0, NA, 2), 2, 4)

  out <- PredictProR:::handle_large_scale_knn(x, k = 2, chunk_size = 50, num_cores = 2)  # one chunk
  expect_false(anyNA(out))
  expect_identical(.plan_call(), user_plan)

  out <- PredictProR:::handle_large_scale_knn(x, k = 2, chunk_size = 2, num_cores = 2)   # two chunks
  expect_false(anyNA(out))
  expect_identical(parallel_plan, c(FALSE, TRUE))
  expect_identical(.plan_call(), user_plan)
})

test_that("KNN imputation puts the caller's future plan back when it fails", {
  user_plan <- .local_user_plan()
  withr::local_envvar(PREDICTPROR_KNN_IMPUTE_BACKEND = "r")
  local_mocked_bindings(
    impute_chunk = function(chunk_data, k = 5) stop("imputation failed"),
    gp_future_lapply_session = function(X, FUN, ...) lapply(X, FUN),
    .package = "PredictProR"
  )
  expect_error(PredictProR:::handle_large_scale_knn(matrix(c(0, NA, 1, 2), 2), k = 2, chunk_size = 50, num_cores = 2),
               "imputation failed")
  expect_identical(.plan_call(), user_plan)
})

test_that("the batched VCF QC puts the caller's future plan back", {
  user_plan <- .local_user_plan()
  parallel_plan <- logical()
  local_mocked_bindings(
    gp_future_lapply_session = function(X, FUN, ...) {
      parallel_plan <<- c(parallel_plan, inherits(future::plan(), "multiprocess"))
      lapply(X, FUN)
    },
    .package = "PredictProR"
  )
  gt <- c("0/0", "0/1", "1/1", "0/0")
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sprintf("S%d", 1:4)), collapse = "\t"),
    vapply(1:5, function(i) {
      paste(c("1", 100 * i, sprintf("rs%d", i), "A", "G", ".", "PASS", ".", "GT", gt[(i + 0:3) %% 4 + 1]), collapse = "\t")
    }, character(1))
  ), vcf)

  res <- PredictProR:::process_vcf_in_batches_parallel(vcf, batch_size = 2, het_threshold = NULL,
                                                      ind_call_rate_threshold = 0.9, recode_format = "0,1,2",
                                                      num_cores = 2)
  expect_identical(res$vcf_processed_data$ID, sprintf("rs%d", 1:5))
  expect_identical(parallel_plan, TRUE)
  expect_identical(.plan_call(), user_plan)
})

test_that("model_execute() puts the caller's future plan back", {
  skip_if_not_installed("BGLR")
  user_plan <- .local_user_plan()
  set.seed(3)
  G <- matrix(rbinom(20 * 15, 2, 0.4), 20, dimnames = list(sprintf("g%02d", 1:20), sprintf("m%02d", 1:15)))
  ph <- data.frame(GID = rownames(G), Yield = rnorm(20))
  ph$Yield[1:3] <- NA
  withr::local_dir(withr::local_tempdir())
  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "GBLUP_BRR", nIter = 60L, burnIn = 20L, thin = 2L, qc_filtering = FALSE, impute = FALSE,
    ld_prunning_qc = FALSE, gmatrix_method = "VanRaden", parallel_mode = "sequential",
    system_database = TRUE, message = FALSE
  ))
  expect_equal(nrow(out$model_results$predicted_values), 20L)
  expect_identical(.plan_call(), user_plan)
})

# The user's thread settings: the variables unset (as in a fresh session) and
# two OpenMP threads, which no pin to one thread leaves in place.
.thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS", "PREDICTPRO_IN_PARALLEL_WORKER")
.thread_state <- function() {
  list(env = Sys.getenv(.thread_vars, unset = NA_character_), omp = RhpcBLASctl::omp_get_max_threads())
}
.local_user_threads <- function(env = parent.frame()) {
  withr::local_envvar(stats::setNames(rep(NA_character_, length(.thread_vars)), .thread_vars), .local_envir = env)
  old_omp <- RhpcBLASctl::omp_get_max_threads()
  RhpcBLASctl::omp_set_num_threads(2L)
  withr::defer(RhpcBLASctl::omp_set_num_threads(old_omp), envir = env)
  state <- .thread_state()
  stopifnot(identical(state$omp, 2L))
  state
}

test_that("sp_apply() puts the caller's thread settings back after running tasks in this process", {
  user_threads <- .local_user_threads()
  # the sequential backend runs the initializer here; one that pins threads
  # (as the CV's worker initializer did) must not outlast the call
  out <- PredictProR:::sp_apply(1:2, function(i) i * 10L, decision = list(backend = "sequential", workers = 1L),
                                seed = FALSE, initializer = PredictProR:::init_single_thread_blas)
  expect_identical(out, list(10L, 20L))
  expect_identical(.thread_state(), user_threads)
  # a future backend with a sequential plan runs the worker set-up here
  # (future warns that the sequential plan takes no `workers`)
  out <- suppressWarnings(PredictProR:::sp_apply(1:2, function(i) i * 10L,
                                                 decision = list(backend = "future", plan = future::sequential, workers = 1L),
                                                 seed = FALSE))
  expect_identical(out, list(10L, 20L))
  expect_identical(.thread_state(), user_threads)
})

test_that("CV tasks the policy leaves to the sequential backend keep the caller's thread settings", {
  user_threads <- .local_user_threads()
  # the cost model's choice when start-up outweighs the work: the tasks stay
  # parallel-eligible, on the sequential backend
  local_mocked_bindings(
    sp_decide_policy = function(models, n_tasks, ...) {
      list(backend = "sequential", workers = 4L, decision_reason = "cost_model_sequential(est_seq=1s,est_par=9s)",
           internal_flags = rep(FALSE, n_tasks), num_gpus = 0L)
    },
    .package = "PredictProR"
  )
  seen <- PredictProR:::gp_execute_partitioned_tasks(
    model_ids = c("BayesB", "GBLUP_BRR"),
    run_task = function(i) Sys.getenv("OMP_NUM_THREADS", unset = "<unset>"),
    decision_args = list(payload_gb = 0),
    seed = FALSE
  )
  # tasks run in this process use the user's threads...
  expect_identical(unlist(seen), c("<unset>", "<unset>"))
  # ...and leave them as they were
  expect_identical(.thread_state(), user_threads)
})

test_that("a cross-validation run leaves the caller's thread settings as it found them", {
  skip_if_not_installed("BGLR")
  user_threads <- .local_user_threads()
  set.seed(3)
  G <- matrix(rbinom(40 * 30, 2, 0.4), 40, dimnames = list(sprintf("g%02d", 1:40), sprintf("m%02d", 1:30)))
  ph <- data.frame(GID = rownames(G), Yield = rnorm(40))
  withr::local_dir(withr::local_tempdir())
  # with RKHS (internally threaded, kept sequential) next to BayesB and
  # GBLUP_BRR, the Windows test machine's policy put the other two on the
  # sequential backend
  out <- suppressWarnings(suppressMessages(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model_cv = c("BayesB", "GBLUP_BRR", "RKHS"), nIter = 60L, burnIn = 20L, thin = 2L,
    cross_validation = TRUE, cross_validation_meth = "Hold_Out", nfolds = 3, replication = 1, random_state = 1L,
    cv_generate_plots = FALSE, cv_evaluation_only = TRUE, eval_metrics = "accuracy",
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, gmatrix_method = "VanRaden",
    parallel_mode = "auto", system_database = FALSE, message = FALSE, verbose = FALSE
  )))
  expect_true(is.list(out))
  expect_identical(.thread_state(), user_threads)
})

# The ML seeds below the routes (uncertainty folds, the folds Lasso chooses
# its penalty on, bootstraps) come from getOption("random_state"). The
# true-prediction refit set it to its task's seed and left it; later runs
# that do not set it (hybrid ML, the CV) read that, or 123 when unset, and
# not their own random_state.
.lasso_data <- function() {
  set.seed(5)
  G <- matrix(rbinom(80 * 40, 2, 0.4), 80, dimnames = list(sprintf("g%02d", 1:80), sprintf("m%02d", 1:40)))
  ph <- data.frame(GID = rownames(G), Yield = drop(G[, 1:8] %*% c(2, -2, 1, 1, -1, 1, -1, 0.5)) + rnorm(80, sd = 0.5))
  list(G = G, ph = ph)
}
.hybrid_lasso_data <- function() {
  set.seed(11)
  fem <- sprintf("F%d", 1:6)
  mal <- sprintf("M%d", 1:5)
  fg <- matrix(sample(c(0, 2), 6 * 30, TRUE), 6, dimnames = list(fem, sprintf("m%02d", 1:30)))
  mg <- matrix(sample(c(0, 2), 5 * 30, TRUE), 5, dimnames = list(mal, sprintf("m%02d", 1:30)))
  ph <- expand.grid(Female = fem, Male = mal, stringsAsFactors = FALSE)
  ph$HybridID <- paste(ph$Female, ph$Male, sep = "_")
  X <- (fg[ph$Female, ] + mg[ph$Male, ]) / 2
  rownames(X) <- ph$HybridID
  ph$Yield <- drop(X[, 1:6] %*% c(2, -1, 1, -1, 1, 0.5)) + rnorm(nrow(ph), sd = 0.5)
  ph$Yield[c(3, 11, 19, 27)] <- NA
  list(X = X, ph = ph)
}

test_that("a true-prediction run leaves the random_state option as it found it", {
  skip_if_not_installed("BGLR")
  set.seed(3)
  G <- matrix(rbinom(20 * 15, 2, 0.4), 20, dimnames = list(sprintf("g%02d", 1:20), sprintf("m%02d", 1:15)))
  ph <- data.frame(GID = rownames(G), Yield = rnorm(20))
  ph$Yield[1:3] <- NA
  withr::local_dir(withr::local_tempdir())
  tp <- function() suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "GBLUP_BRR", nIter = 60L, burnIn = 20L, thin = 2L, qc_filtering = FALSE, impute = FALSE,
    ld_prunning_qc = FALSE, gmatrix_method = "VanRaden", parallel_mode = "sequential", random_state = 7L,
    system_database = FALSE, message = FALSE
  ))
  withr::local_options(random_state = NULL)
  tp()
  expect_null(getOption("random_state"))
  options(random_state = 77L)   # an option of the user's own
  tp()
  expect_identical(getOption("random_state"), 77L)
})

test_that("a cross-validation's ML seeds come from its random_state, not from the session", {
  skip_if_not_installed("glmnet")
  d <- .lasso_data()
  withr::local_dir(withr::local_tempdir())
  withr::local_options(random_state = NULL)
  cv <- function(session_option) {
    options(random_state = session_option)
    out <- suppressWarnings(suppressMessages(PredictProR::model_execute(
      pheno_data = d$ph, geno_data = d$G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
      GS_model_cv = "Lasso", cross_validation = TRUE, cross_validation_meth = "Hold_Out", nfolds = 3,
      replication = 1, random_state = 7L, cv_generate_plots = FALSE, cv_evaluation_only = TRUE,
      eval_metrics = "accuracy", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
      gmatrix_method = "VanRaden", parallel_mode = "sequential", system_database = FALSE, message = FALSE,
      verbose = FALSE
    )))
    expect_identical(getOption("random_state"), session_option)
    out$cv_results_raw[[1]]$ypred_cv_Reps_all
  }
  first <- cv(1L)
  expect_true(any(is.finite(first$yhat)))
  expect_identical(cv(2L), first)
  expect_identical(cv(NULL), first)
})

# Workers do not have the session's options: a CV Lasso there took 123 while
# the same task in this process took the run's random_state (99 and 123 choose
# different penalties on these data; 7 and 123 the same one).
test_that("a cross-validation's ML seeds are the same on parallel workers as in this process", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  skip_if(parallel::detectCores() < 3L, "needs more than one usable core")
  withr::local_envvar(GP_PAR_WORKER_OVERHEAD_GB = "0.5")
  d <- .lasso_data()
  withr::local_dir(withr::local_tempdir())
  withr::local_options(random_state = 5L)   # an option of the user's own
  cv <- function(mode) {
    out <- suppressWarnings(suppressMessages(PredictProR::model_execute(
      pheno_data = d$ph, geno_data = d$G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
      GS_model_cv = c("Lasso", "Ridge_Regression"), cross_validation = TRUE,
      cross_validation_meth = "Hold_Out", nfolds = 3, replication = 1, random_state = 99L,
      cv_generate_plots = FALSE, cv_evaluation_only = TRUE, eval_metrics = "accuracy",
      qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, gmatrix_method = "VanRaden",
      parallel_mode = mode, num_cores = 2L, system_database = FALSE, message = FALSE, verbose = FALSE
    )))
    expect_identical(getOption("random_state"), 5L)
    out
  }
  lasso <- function(out) {
    raw <- out$cv_results_raw
    raw[[which(vapply(raw, function(r) identical(r$model, "Lasso"), logical(1)))]]$ypred_cv_Reps_all
  }
  seq_out <- cv("sequential")
  par_out <- cv("base_parallel")
  meta <- par_out$run_metadata
  expect_identical(meta$value[meta$key == "policy_backend"], "base_parallel")
  expect_identical(meta$value[meta$key == "policy_parallel_count"], "2")
  expect_true(any(is.finite(lasso(seq_out)$yhat)))
  expect_identical(lasso(par_out), lasso(seq_out))
})

test_that("a hybrid ML run's seeds come from its random_state, not from the session", {
  skip_if_not_installed("glmnet")
  d <- .hybrid_lasso_data()
  withr::local_dir(withr::local_tempdir())
  withr::local_options(random_state = NULL)
  hyb <- function(session_option, random_state = 7L) {
    options(random_state = session_option)
    out <- suppressWarnings(suppressMessages(PredictProR::model_execute(
      pheno_data = d$ph, geno_data = d$X, response = "Yield", gen_name = "HybridID",
      female_parent = "Female", male_parent = "Male", GS_model = "Lasso", response_family = "gaussian",
      hybrid_ml = TRUE, random_state = random_state, system_database = FALSE, message = FALSE, verbose = FALSE
    )))
    expect_identical(getOption("random_state"), session_option)
    out$model_results$predicted_values
  }
  first <- hyb(1L)
  expect_true(all(is.finite(first$PEV)))
  expect_identical(hyb(2L), first)
  expect_identical(hyb(NULL), first)
  # and the run's own random_state does reach them
  expect_false(identical(hyb(NULL, random_state = 8L), first))
})

test_that("objects assigned for ASReml during a run are put back when it ends", {
  assign("gmatrix_inv", "the user's object", envir = .GlobalEnv)
  withr::defer(suppressWarnings(rm("gmatrix_inv", "omic1_inv", envir = .GlobalEnv)))
  close_run <- PredictProR:::gp_asreml_globals_scope()
  PredictProR:::gp_asreml_assign_run_globals(list(gmatrix_inv = "fold 1", omic1_inv = "fold 1"))
  PredictProR:::gp_asreml_assign_run_globals(list(gmatrix_inv = "fold 2", omic1_inv = "fold 2"))
  expect_identical(get("gmatrix_inv", envir = .GlobalEnv), "fold 2")  # ASReml finds the current fit's
  # a nested run does not close the outer one
  PredictProR:::gp_asreml_globals_scope()()
  expect_identical(get("omic1_inv", envir = .GlobalEnv), "fold 2")
  close_run()
  expect_identical(get("gmatrix_inv", envir = .GlobalEnv), "the user's object")
  expect_false(exists("omic1_inv", envir = .GlobalEnv, inherits = FALSE))
})

test_that("outside a run the ASReml objects stay, for the caller's predict()", {
  withr::defer(suppressWarnings(rm("omic1_inv", envir = .GlobalEnv)))
  PredictProR:::gp_asreml_assign_run_globals(list(omic1_inv = "direct call"))
  expect_identical(get("omic1_inv", envir = .GlobalEnv), "direct call")
})

test_that("a single-trait ASReml run with cross-validation leaves the session as it was", {
  skip_if_not_installed("asreml")
  assign("gmatrix_inv", "the user's object", envir = .GlobalEnv)
  withr::defer(if (exists("gmatrix_inv", envir = .GlobalEnv, inherits = FALSE)) rm("gmatrix_inv", envir = .GlobalEnv))
  before <- ls(.GlobalEnv, all.names = TRUE)
  # ASReml works from its namespace: its licence and engine are set up when it
  # loads (.onLoad); attaching it only prints its version
  if ("package:asreml" %in% search()) {
    detach("package:asreml")
    withr::defer(suppressPackageStartupMessages(library(asreml)))
  }
  search_before <- search()
  registry <- PredictProR:::PredictProR_asreml_data_registry
  registry_before <- ls(registry)
  in_store <- integer()
  register <- PredictProR:::asreml_model_data_register
  local_mocked_bindings(asreml_model_data_register = function(pheno_data) {
    key <- register(pheno_data)
    in_store <<- c(in_store, length(setdiff(ls(registry), registry_before)))
    key
  }, .package = "PredictProR")
  # sim_single_env(n = 60, p = 300) of inst/examples/scenarios/00_setup.R
  set.seed(1)
  n <- 60L; p <- 300L; ids <- sprintf("L%03d", seq_len(n))
  freq <- stats::runif(p, 0.1, 0.9)
  X <- matrix(2L * stats::rbinom(n * p, 1, rep(freq, each = n)), n, p, dimnames = list(ids, sprintf("m%04d", seq_len(p))))
  X <- X[, apply(X, 2, stats::var) > 0, drop = FALSE]
  set.seed(2)
  qtl <- sample(ncol(X), 30)
  g <- as.numeric(scale(X[, qtl, drop = FALSE] %*% stats::rnorm(30)))
  set.seed(3)
  ph <- data.frame(GID = ids, Yield = 10 + g + stats::rnorm(n, sd = sqrt(0.3 / 0.7)))
  ph$Yield[sample(n, 12)] <- NA
  withr::local_dir(withr::local_tempdir())
  out <- tryCatch(suppressWarnings(suppressMessages(PredictProR::model_execute(
    pheno_data = ph, geno_data = X, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model_cv = "GBLUP", engine = "asreml", gmatrix_method = "VanRaden",
    cross_validation = TRUE, cross_validation_meth = "K-Folds", nfolds = 3, replication = 1,
    random_state = 1L, cv_generate_plots = FALSE, eval_metrics = c("accuracy", "root_mean_squared_error"),
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = "sequential", system_database = FALSE, message = FALSE, verbose = FALSE
  ))), error = function(e) e)
  if (inherits(out, "error") && grepl("License|licen", conditionMessage(out))) skip("ASReml license unavailable")
  expect_false(inherits(out, "error"), info = if (inherits(out, "error")) conditionMessage(out) else "")
  expect_true(all(is.finite(out$model_results$predicted_values$Predicted_value)))
  expect_identical(get("gmatrix_inv", envir = .GlobalEnv), "the user's object")
  expect_identical(setdiff(ls(.GlobalEnv, all.names = TRUE), c(before, ".Random.seed")), character())
  expect_identical(search(), search_before)
  # three folds, then the refit: a fold's data is released once the fold is predicted
  expect_identical(in_store, c(1L, 1L, 1L, 1L))
  expect_identical(ls(registry), registry_before)
})

test_that("registering ASReml data leaves the random-number stream alone", {
  withr::local_seed(7)
  seed_before <- get(".Random.seed", envir = .GlobalEnv)
  df <- data.frame(GID = c("a", "b"), y = 1:2)
  keys <- c(PredictProR:::asreml_model_data_register(df), PredictProR:::asreml_model_data_register(df))
  withr::defer(rm(list = keys, envir = PredictProR:::PredictProR_asreml_data_registry))
  expect_identical(get(".Random.seed", envir = .GlobalEnv), seed_before)
  expect_false(identical(keys[1], keys[2]))
  expect_identical(PredictProR:::asreml_model_data_lookup(keys[2]), df)
})

test_that("a new ASReml data key never replaces data still in the store", {
  # test-python-runtime-policy.R empties the runtime cache, the count with it
  registry <- PredictProR:::PredictProR_asreml_data_registry
  cache <- PredictProR:::PredictProR_runtime_cache
  suppressWarnings(rm("asreml_data_count", envir = cache))
  first <- PredictProR:::asreml_model_data_register(data.frame(y = 1))
  withr::defer(rm(list = intersect(c(first, second), ls(registry)), envir = registry))
  rm("asreml_data_count", envir = cache)
  second <- PredictProR:::asreml_model_data_register(data.frame(y = 2))
  expect_false(identical(second, first))
  expect_identical(PredictProR:::asreml_model_data_lookup(first), data.frame(y = 1))
})

test_that("ASReml data registered during a run is released when it ends; outside a run it stays", {
  registry <- PredictProR:::PredictProR_asreml_data_registry
  close_run <- PredictProR:::gp_asreml_globals_scope()
  in_run <- PredictProR:::asreml_model_data_register(data.frame(y = 1))
  PredictProR:::gp_asreml_globals_scope()()  # a nested run does not release it
  expect_true(exists(in_run, envir = registry, inherits = FALSE))
  close_run()
  expect_false(exists(in_run, envir = registry, inherits = FALSE))
  direct <- PredictProR:::asreml_model_data_register(data.frame(y = 2))
  withr::defer(rm(list = direct, envir = registry))
  expect_true(exists(direct, envir = registry, inherits = FALSE))
})
