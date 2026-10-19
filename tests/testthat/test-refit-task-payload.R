# The refit's task is a closure made inside model_execute(), so every worker
# received that whole frame: 528 objects and 83.5 MB for a 3-trait refit on the
# Windows test machine, of which the task reads 15.4 MB. The policy could not
# measure that frame (reading its empty `...` fails) and fell back to the data.

test_that("a trimmed task keeps only the names it and its helpers read", {
  make <- function() {
    unused <- stats::runif(2e5)
    data <- 1:10
    found <- NULL
    lookups <- 0L
    find <- function() {
      if (is.null(found)) {
        lookups <<- lookups + 1L
        found <<- "python"
      }
      found
    }
    init <- function(model = NULL) find()
    task <- function(i) {
      data <- data
      list(total = sum(data) + i, py = init())
    }
    list(task = task, init = init, find = find)
  }
  # as inside model_execute(), whose frame's parent is the namespace (a test
  # environment as the parent would be sent with the copy)
  environment(make) <- asNamespace("PredictProR")
  f <- make()
  # looked up before the trim, as the refit does in the main session
  f$find()
  out <- PredictProR:::gp_trim_closure_frame(f[c("task", "init")])
  env <- environment(out$task)
  expect_setequal(ls(env, all.names = TRUE), c("data", "init", "find", "found", "lookups"))
  expect_identical(environment(out$init), env)
  expect_identical(environment(get("find", envir = env)), env)
  # `unused` alone is 1.6 MB (the rest is the test file's source references)
  expect_lt(length(serialize(out$task, NULL)), 5e5)
  expect_identical(out$task(1), list(total = 56, py = "python"))
  expect_identical(out$init(), "python")
  # the copy kept what was found and did not look again
  expect_identical(get("lookups", envir = env), 1L)
  # the original frame is untouched
  expect_true(exists("unused", envir = environment(f$task), inherits = FALSE))
})

test_that("a trimmed task keeps the fields it is asked to, and what they partially match", {
  make <- function() {
    alpha <- 1
    threshold_value <- 2
    other <- 3
    function(i) i
  }
  out <- PredictProR:::gp_trim_closure_frame(list(make()), fields = c("alpha", "threshold"))
  # `ctx$threshold` on the whole frame would return threshold_value
  expect_setequal(ls(environment(out[[1L]]), all.names = TRUE), c("alpha", "threshold_value"))
})

test_that("functions from the global or a namespace environment are not trimmed", {
  fun <- PredictProR:::gp_merge_task_context
  expect_identical(PredictProR:::gp_trim_closure_frame(list(fun))[[1L]], fun)
  glob <- local(function(i) i, envir = globalenv())
  expect_identical(environment(PredictProR:::gp_trim_closure_frame(list(glob))[[1L]]), globalenv())
})

# The fields come from reading gp_run_best_model_task()'s code; any other use
# of its context (passing it on whole, a computed name) would be missed.
test_that("the refit task reads its context only as ctx$name or ctx[[\"name\"]]", {
  bad <- character()
  walk <- function(e) {
    if (is.name(e) && identical(as.character(e), "ctx")) {
      bad <<- c(bad, "ctx")
      return(invisible())
    }
    if (!is.call(e)) return(invisible())
    f <- e[[1L]]
    if ((identical(f, as.name("$")) || identical(f, as.name("[["))) &&
        identical(e[[2L]], as.name("ctx")) && length(e) == 3L &&
        (is.name(e[[3L]]) || is.character(e[[3L]]))) {
      return(invisible())
    }
    for (x in as.list(e)) if (!missing(x)) walk(x)
  }
  walk(body(PredictProR:::gp_run_best_model_task))
  expect_identical(bad, character())
  fields <- PredictProR:::gp_best_model_task_fields()
  expect_true(all(c("pheno_clean", "ml_dat_res", "gen_name", "var_cov_str", "nIter") %in% fields))
})

test_that("the refit's initializer shares the task's trimmed frame and the interpreter found first", {
  make <- function() {
    unused <- stats::runif(2e5)
    py <- NULL
    resolved <- FALSE
    lookups <- 0L
    resolve <- function() {
      if (!resolved) {
        lookups <<- lookups + 1L
        py <<- "python"
        resolved <<- TRUE
      }
      py
    }
    init_py <- function(model = NULL) resolve()
    run_one_best <- function(i) {
      init_py("PyOnlyModel")
      list(ok = i)
    }
    list(run = run_one_best, init = init_py, resolve = resolve)
  }
  f <- make()
  seen <- new.env()
  local_mocked_bindings(
    gp_execute_partitioned_tasks = function(model_ids, run_task, initializer = NULL, ...) {
      seen$run_task <- run_task
      seen$initializer <- initializer
      list(list(ok = TRUE))
    },
    .package = "PredictProR"
  )
  PredictProR:::gp_execute_best_model_tasks(
    best_models = data.frame(model = "PyOnlyModel", trait = "Yield", stringsAsFactors = FALSE),
    run_one_best = f$run,
    canonical_names = character(),
    friendly_names = character(),
    parallel_mode = "auto",
    num_cores = 2L,
    globals_max_GB = 1,
    verbose = FALSE,
    parallel_backend_prefer_fork = FALSE,
    pheno_clean = list(pheno_clean_data = data.frame(GID = c("g1", "g2"), Yield = c(1, 2))),
    ml_dat_res = list(),
    gmatrix_kernel_model_ready_list = list(),
    geno_omic_model_ready_list = list(),
    response = "Yield",
    init_py = f$init,
    resolve_python = f$resolve
  )
  env <- environment(seen$run_task)
  expect_false(exists("unused", envir = env, inherits = FALSE))
  expect_identical(environment(seen$initializer), env)
  expect_true(get("resolved", envir = env))
  seen$initializer()
  seen$run_task(1L)
  expect_identical(get("lookups", envir = env), 1L)
})

test_that("the refit sends its task without the objects it never reads", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(11)
  n <- 60
  p <- 120
  ids <- sprintf("L%03d", seq_len(n))
  G <- matrix(2L * stats::rbinom(n * p, 1, 0.5), n, p, dimnames = list(ids, sprintf("m%03d", seq_len(p))))
  pheno <- data.frame(GID = ids, Yield = as.numeric(G[, 1:5] %*% stats::rnorm(5)) + stats::rnorm(n))
  pheno$Yield[1:10] <- NA
  seen <- new.env()
  real_payload_gb <- PredictProR:::gp_parallel_closure_payload_gb
  local_mocked_bindings(
    gp_parallel_closure_payload_gb = function(fun) {
      seen$fun <- fun
      seen$payload_gb <- tryCatch(real_payload_gb(fun), error = function(e) NA_real_)
      seen$payload_gb
    },
    gp_python_prepare_runtimes = function(...) invisible(character()),
    .package = "PredictProR"
  )
  withr::local_dir(withr::local_tempdir())
  # a formula carries the environment it was written in; this one's would hold
  # the test's objects, `res` included, where a user's is the global one
  random <- local(~GID, envir = globalenv())
  res <- suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = G, gen_name = "GID", response = "Yield", random = random,
    GS_model = "GBLUP_BRR", gmatrix_method = "VanRaden", nIter = 200L, burnIn = 50L, thin = 2L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE, verbose = FALSE
  ))
  expect_true(is.function(seen$fun))
  env <- environment(seen$fun)
  expect_false(exists("input_guardrail_ctx", envir = env, inherits = FALSE))
  mb <- function(x) length(serialize(x, NULL)) / 1024^2
  data_mb <- mb(mget(c("pheno_clean", "ml_dat_res", "geno_omic_model_ready_list",
                       "gmatrix_kernel_model_ready_list"), envir = env))
  # as an installed package sends it: under load_all() each function keeps its
  # file's text (1.5 MB for model_execute_para.R)
  sent <- local({
    vals <- mget(ls(env, all.names = TRUE), envir = env)
    copy <- new.env(parent = parent.env(env))
    for (i in seq_along(vals)) {
      if (is.function(vals[[i]])) {
        moved <- identical(environment(vals[[i]]), env)
        vals[[i]] <- utils::removeSource(vals[[i]])
        if (moved) environment(vals[[i]]) <- copy
      }
    }
    list2env(vals, envir = copy)
    fun <- utils::removeSource(seen$fun)
    environment(fun) <- copy
    fun
  })
  expect_lt(mb(sent), 1.5 * data_mb)
  # nothing is sent twice: two names bound to one object are serialized twice
  # (the kernel bank went as outer_kernel_bank and gmatrix_kernel_model_ready_list)
  vals <- mget(ls(env, all.names = TRUE), envir = env)
  big <- Filter(function(x) !is.function(x) && mb(x) > 0.01, vals)
  pairs <- if (length(big) > 1L) utils::combn(names(big), 2, simplify = FALSE) else list()
  twice <- Filter(function(p) identical(big[[p[1]]], big[[p[2]]]), pairs)
  expect_identical(vapply(twice, paste, "", collapse = " = "), character())
  # the policy measures what is sent (model_execute()'s own frame failed: its
  # empty `...` cannot be read)
  expect_true(is.finite(seen$payload_gb) && seen$payload_gb > 0)
  expect_false(is.null(res$run_metadata))
})
