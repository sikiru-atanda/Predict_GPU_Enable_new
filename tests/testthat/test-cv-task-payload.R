# The CV's task function is a closure in models_execute_crossval()'s frame.
# It carried all of it (127 objects, 52.3 MB in a 3-trait run), and params
# with five of its elements held under their own names too. Workers now get a
# copy with only what the task reads, each object once.

# every object of at least min_kb in a function's frames and the lists in them
.cv_leaves <- function(fun, min_kb = 10) {
  out <- list()
  seen_env <- list()
  walk <- function(x, path) {
    if (is.environment(x) || is.function(x)) return(invisible())
    if (is.list(x) && !is.data.frame(x) && is.null(attr(x, "class", exact = TRUE))) {
      for (i in seq_along(x)) walk(x[[i]], paste0(path, "[[", i, "]]"))
    } else if (length(serialize(x, NULL)) >= min_kb * 1024) {
      out[[path]] <<- x
    }
  }
  env <- environment(fun)
  while (is.environment(env) && !isNamespace(env) && !identical(env, globalenv()) &&
         !identical(env, emptyenv()) && !any(vapply(seen_env, identical, logical(1), env))) {
    seen_env[[length(seen_env) + 1L]] <- env
    for (nm in ls(env, all.names = TRUE)) {
      val <- get(nm, envir = env)
      if (is.function(val) && !isNamespace(environment(val) %||% emptyenv()) &&
          !any(vapply(seen_env, identical, logical(1), environment(val)))) {
        inner <- .cv_leaves(val, min_kb)
        out[paste0(nm, "/", names(inner))] <- inner
      } else {
        walk(val, nm)
      }
    }
    env <- parent.env(env)
  }
  out
}

test_that("CV workers get the task with each object once, and the same predictions", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  lib <- PredictProR:::gp_predictpror_loaded_lib()
  skip_if(is.null(lib), "PredictProR not loaded from an installed library")
  set.seed(5)
  n <- 60
  p <- 150
  ids <- sprintf("L%03d", seq_len(n))
  G <- matrix(2L * stats::rbinom(n * p, 1, 0.5), n, p, dimnames = list(ids, sprintf("m%03d", seq_len(p))))
  pheno <- data.frame(GID = ids, Yield = as.numeric(G[, 1:5] %*% stats::rnorm(5)) + stats::rnorm(n))
  run <- function(mode) suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = G, gen_name = "GID", response = "Yield", random = ~GID,
    GS_model_cv = c("BRR", "GBLUP_BRR"), gmatrix_method = "VanRaden", nIter = 200L, burnIn = 50L, thin = 2L,
    cross_validation = TRUE, cross_validation_meth = "K-Folds", nfolds = 2L, replication = 1,
    cv_evaluation_only = TRUE, cv_generate_plots = FALSE, random_state = 7L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = mode, num_cores = 2L, system_database = TRUE, message = FALSE, verbose = FALSE
  ))
  withr::local_dir(withr::local_tempdir())
  seen <- new.env()
  real_sp_apply <- PredictProR:::sp_apply
  local_mocked_bindings(
    sp_apply = function(X, FUN, ...) {
      if (is.null(seen$fun)) seen$fun <- FUN
      real_sp_apply(X, FUN, ...)
    },
    gp_python_prepare_runtimes = function(...) invisible(character()),
    .package = "PredictProR"
  )
  par <- run("base_parallel")
  skip_if(is.null(seen$fun), "the CV ran nothing in parallel")
  leaves <- .cv_leaves(seen$fun)
  expect_gt(length(leaves), 0L)
  pairs <- if (length(leaves) > 1L) utils::combn(names(leaves), 2, simplify = FALSE) else list()
  twice <- Filter(function(pr) identical(leaves[[pr[1]]], leaves[[pr[2]]]), pairs)
  expect_identical(vapply(twice, paste, "", collapse = " = "), character())

  seq_res <- run("sequential")
  yhat <- function(res) lapply(res$cv_results_raw, function(el) el$ypred_cv_Reps_all$yhat)
  expect_identical(yhat(par), yhat(seq_res))
})
