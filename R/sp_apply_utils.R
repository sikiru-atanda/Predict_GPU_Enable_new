
# Library directory of the installed PredictProR this session loaded; NULL for
# a source tree (devtools::load_all), which workers cannot load from.
gp_predictpror_loaded_lib <- function() {
  path <- tryCatch(getNamespaceInfo("PredictProR", "path"), error = function(e) "")
  if (!nzchar(path) || !file.exists(file.path(path, "Meta", "package.rds"))) return(NULL)
  normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
}

# Identity of the PredictProR loaded in this process: version plus the
# install's build stamp ("Built"), so a reinstall that kept the version number
# but changed the code is still told apart. A source tree (load_all) has no
# build stamp and is identified by its version only.
gp_predictpror_identity <- function() {
  pkg_dir <- tryCatch(getNamespaceInfo("PredictProR", "path"), error = function(e) "")
  version <- tryCatch(as.character(getNamespaceVersion("PredictProR")), error = function(e) "")
  built <- if (nzchar(pkg_dir) && file.exists(file.path(pkg_dir, "Meta", "package.rds"))) {
    tryCatch(utils::packageDescription("PredictProR", lib.loc = dirname(pkg_dir))[["Built"]] %||% "",
             error = function(e) "")
  } else ""
  list(version = unname(version), built = built, path = pkg_dir)
}

gp_parallel_runtime_env <- function() {
  # Interpreters already resolved in this session are handed to the workers
  # so each worker does not re-probe every candidate Python (each probe
  # starts a process and may import torch).
  cached_path <- function(key) {
    if (is.na(key) || !exists(key, envir = PredictProR_runtime_cache, inherits = FALSE)) return(NULL)
    path <- get(key, envir = PredictProR_runtime_cache, inherits = FALSE)
    if (is.character(path) && length(path) == 1L && file.exists(path)) path else NULL
  }
  resolved <- list(
    PREDICTPRO_ML_PYTHON = cached_path(grep("^ml_python::", ls(PredictProR_runtime_cache), value = TRUE)[1L]),
    PREDICTPRO_GP_PYTHON = cached_path("preferred_python::gp"),
    PREDICTPRO_DL_PYTHON = cached_path("preferred_python::dl")
  )
  env_keys <- c(
    "PREDICTPRO_GP_PYTHON",
    "PREDICTPRO_GP_DEVICE",
    "PREDICTPRO_PYTHON",
    "PREDICTPRO_ML_PYTHON",
    "PREDICTPRO_DL_PYTHON",
    "PREDICTPRO_PYTHON_DL",
    "RETICULATE_PYTHON",
    "CUDA_VISIBLE_DEVICES",
    "PYTHONPATH",
    "R_LIBS",
    "R_LIBS_USER",
    "R_LIBS_SITE",
    "PREDICTPRO_R_LIBPATHS"
  )
  vals <- Sys.getenv(env_keys, unset = NA_character_)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  for (env_name in names(resolved)) {
    if (!is.null(resolved[[env_name]]) && !env_name %in% names(vals)) {
      vals[[env_name]] <- resolved[[env_name]]
    }
  }
  # Workers must load the PredictProR this session is running. The library it
  # was loaded from goes first: with library(PredictProR, lib.loc = ...) (or
  # renv/project libraries) it is not on .libPaths(), and workers silently
  # loaded an older copy from the default library instead.
  lib_paths <- normalizePath(c(gp_predictpror_loaded_lib(), .libPaths()), winslash = "/", mustWork = FALSE)
  lib_paths <- unique(lib_paths[nzchar(lib_paths)])
  # Every worker compares its own PredictProR with these and refuses to run
  # on a mismatch (gp_parallel_check_worker_version). A load_all() session is
  # checked by version: its workers load an installed copy, which must match.
  identity <- gp_predictpror_identity()
  if (nzchar(identity$version)) {
    vals[["PREDICTPRO_EXPECTED_VERSION"]] <- identity$version
  }
  if (nzchar(identity$built)) {
    vals[["PREDICTPRO_EXPECTED_BUILD"]] <- identity$built
  }
  if (length(lib_paths)) {
    lib_env <- paste(lib_paths, collapse = .Platform$path.sep)
    vals[["R_LIBS"]] <- lib_env
    vals[["PREDICTPRO_R_LIBPATHS"]] <- lib_env
  }
  vals
}

gp_env_map_value <- function(env_map, key) {
  if (key %in% names(env_map)) as.character(env_map[[key]]) else NULL
}

gp_parallel_apply_libpaths <- function(env_map) {
  lib_env <- gp_env_map_value(env_map, "PREDICTPRO_R_LIBPATHS")
  if (is.null(lib_env) || !nzchar(lib_env)) {
    lib_env <- gp_env_map_value(env_map, "R_LIBS")
  }
  if (is.null(lib_env) || !nzchar(lib_env)) {
    return(invisible(NULL))
  }
  lib_paths <- unlist(strsplit(lib_env, .Platform$path.sep, fixed = TRUE), use.names = FALSE)
  lib_paths <- normalizePath(lib_paths[nzchar(lib_paths)], winslash = "/", mustWork = FALSE)
  lib_paths <- unique(lib_paths[dir.exists(lib_paths)])
  if (length(lib_paths)) {
    .libPaths(unique(c(lib_paths, .libPaths())))
  }
  invisible(NULL)
}

gp_parallel_set_runtime_env <- function(env_map) {
  if (is.null(env_map) || !length(env_map)) {
    return(invisible(NULL))
  }
  do.call(Sys.setenv, as.list(env_map))
  gp_parallel_apply_libpaths(env_map)
  invisible(NULL)
}

gp_parallel_restore_runtime_env <- function(old_env) {
  if (is.null(old_env) || !length(old_env)) {
    return(invisible(NULL))
  }
  missing <- names(old_env)[is.na(old_env)]
  present <- old_env[!is.na(old_env)]
  if (length(missing)) {
    Sys.unsetenv(missing)
  }
  if (length(present)) {
    do.call(Sys.setenv, as.list(present))
  }
  invisible(NULL)
}

gp_psock_port_candidates <- function(pid = Sys.getpid(),
                                     count = 64L,
                                     lower = 20000L,
                                     upper = 59999L) {
  pid <- suppressWarnings(as.double(pid)[1L])
  count <- suppressWarnings(as.integer(count)[1L])
  lower <- suppressWarnings(as.integer(lower)[1L])
  upper <- suppressWarnings(as.integer(upper)[1L])
  if (!is.finite(pid) || !is.finite(count) || count < 1L ||
      !is.finite(lower) || !is.finite(upper) || lower < 1024L ||
      upper > 65535L || lower > upper) {
    stop("Invalid PSOCK port-candidate configuration.", call. = FALSE)
  }
  span <- as.double(upper - lower + 1L)
  count <- min(count, as.integer(span))
  # 7919 is coprime to the 40,000-port default span. Adjacent process IDs
  # therefore start at distinct candidates while still walking the full span.
  start <- (pid * 7919) %% span
  offsets <- (as.double(seq_len(count) - 1L) * 7919) %% span
  as.integer(lower + ((start + offsets) %% span))
}

gp_psock_cluster_port <- function() {
  configured_raw <- trimws(Sys.getenv("R_PARALLEL_PORT", unset = ""))
  if (nzchar(configured_raw)) {
    configured <- suppressWarnings(as.integer(configured_raw))
    if (length(configured) == 1L && is.finite(configured) &&
        configured >= 1024L && configured <= 65535L) {
      return(configured)
    }
    warning(
      "Ignoring invalid R_PARALLEL_PORT; selecting a process-distinct port.",
      call. = FALSE
    )
  }

  candidates <- gp_psock_port_candidates()
  if (requireNamespace("parallelly", quietly = TRUE)) {
    return(as.integer(parallelly::freePort(
      ports = candidates,
      default = candidates[[1L]],
      randomize = FALSE
    )))
  }
  candidates[[1L]]
}

gp_make_psock_cluster <- function(workers) {
  workers <- max(1L, as.integer(workers)[1L])
  port <- gp_psock_cluster_port()
  logger::log_info("Starting PSOCK cluster on process-distinct port {port}")
  parallel::makePSOCKcluster(workers, port = port)
}

gp_parallel_init_runtime <- function(env_map = NULL) {
  if (!is.null(env_map) && length(env_map)) {
    do.call(Sys.setenv, as.list(env_map))
    lib_env <- gp_env_map_value(env_map, "PREDICTPRO_R_LIBPATHS")
    if (is.null(lib_env) || !nzchar(lib_env)) {
      lib_env <- gp_env_map_value(env_map, "R_LIBS")
    }
    if (!is.null(lib_env) && nzchar(lib_env)) {
      lib_paths <- unlist(strsplit(lib_env, .Platform$path.sep, fixed = TRUE), use.names = FALSE)
      lib_paths <- normalizePath(lib_paths[nzchar(lib_paths)], winslash = "/", mustWork = FALSE)
      lib_paths <- unique(lib_paths[dir.exists(lib_paths)])
      if (length(lib_paths)) {
        .libPaths(unique(c(lib_paths, .libPaths())))
      }
    }
  }
  py_bin <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (!nzchar(py_bin)) {
    py_bin <- Sys.getenv("PREDICTPRO_PYTHON", unset = "")
  }
  if (nzchar(py_bin) && exists("gp_init_python_once", mode = "function")) {
    try(gp_init_python_once(py_bin), silent = TRUE)
  }
  invisible(NULL)
}

# A worker running a different PredictProR than the session would compute
# fold results with other code (it happened with a stale copy in the default
# library). Refuse instead of mixing versions silently.
gp_parallel_check_worker_version <- function(expected = Sys.getenv("PREDICTPRO_EXPECTED_VERSION", unset = ""),
                                             expected_build = Sys.getenv("PREDICTPRO_EXPECTED_BUILD", unset = "")) {
  loaded <- gp_predictpror_identity()
  if (nzchar(expected) && !identical(loaded$version, expected)) {
    stop(sprintf(
      "A parallel worker loaded PredictProR %s from %s, but this session runs %s. Install %s into that library, remove the old copy, or use parallel_mode = \"sequential\".",
      loaded$version, dirname(loaded$path), expected, expected
    ), call. = FALSE)
  }
  if (nzchar(expected_build) && nzchar(loaded$built) && !identical(loaded$built, expected_build)) {
    stop(sprintf(
      "A parallel worker loaded a different install of PredictProR %s from %s (built %s; this session's was built %s). Reinstall PredictProR in one library only, or use parallel_mode = \"sequential\".",
      loaded$version, dirname(loaded$path), loaded$built, expected_build
    ), call. = FALSE)
  }
  invisible(loaded$version)
}

# Runs future.apply::future_lapply on workers that load this session's
# PredictProR: the loaded library goes first on .libPaths() (future workers
# copy the session's paths) and every task checks the worker's version/build
# before running. For helpers that set their own future plan.
gp_future_lapply_session <- function(X, FUN, ...) {
  runtime_env <- gp_parallel_runtime_env()
  loaded_lib <- gp_predictpror_loaded_lib()
  if (!is.null(loaded_lib) && !loaded_lib %in% normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)) {
    old_lib_paths <- .libPaths()
    .libPaths(c(loaded_lib, old_lib_paths))
    on.exit(.libPaths(old_lib_paths), add = TRUE)
  }
  future.apply::future_lapply(X, gp_session_checked_fun(FUN, runtime_env, Sys.getpid()), ...)
}

# Built in its own small frame so the wrapper ships only these three objects
# to the workers (not the caller's X). In-process tasks (plan sequential) skip
# the env setup, which would otherwise leak into the session.
gp_session_checked_fun <- function(task_fun, runtime_env, session_pid) {
  force(task_fun); force(runtime_env); force(session_pid)
  function(x, ...) {
    if (!identical(Sys.getpid(), session_pid)) {
      gp_parallel_init_runtime(runtime_env)
      gp_parallel_check_worker_version()
    }
    task_fun(x, ...)
  }
}

gp_parallel_init_worker_runtime <- function(env_map = NULL) {
  gp_parallel_init_runtime(env_map)
  gp_parallel_check_worker_version()
  init_single_thread_blas()
  # Marks this process as a parallel worker: nested process pools (e.g. the
  # ML bootstrap) stay single-process here to avoid oversubscription.
  Sys.setenv(PREDICTPRO_IN_PARALLEL_WORKER = "1")
  invisible(NULL)
}

# Processes for the Python ML bootstrap (prediction uncertainty refits).
# PREDICTPRO_BOOTSTRAP_JOBS overrides; inside a parallel worker it is 1;
# otherwise all but one logical core, capped by the memory budget (one Python
# process with its data copy per job) and by the number of replicates.
gp_ml_bootstrap_jobs <- function(n_bootstrap, data_gb = 0) {
  n_bootstrap <- max(1L, as.integer(n_bootstrap %||% 1L))
  user <- suppressWarnings(as.integer(Sys.getenv("PREDICTPRO_BOOTSTRAP_JOBS", "")))
  if (is.finite(user) && user >= 1L) {
    return(min(user, n_bootstrap))
  }
  cores <- tryCatch(gp_detect_cores(logical = TRUE, reserve = 1L), error = function(e) 1L)
  budget <- tryCatch(gp_parallel_memory_budget_gb(), error = function(e) Inf)
  if (identical(Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER", ""), "1")) {
    # Inside a parallel worker: use this worker's share of the machine, so a
    # long bootstrap is not left single-process while other cores idle.
    n_workers <- suppressWarnings(as.integer(Sys.getenv("PREDICTPRO_PARALLEL_WORKERS", "")))
    if (!is.finite(n_workers) || n_workers < 1L) {
      return(1L)
    }
    cores <- max(1L, as.integer(floor(cores / n_workers)))
    if (is.finite(budget)) budget <- budget / n_workers
  }
  # A pool process holds the pickled start-up data and a resampled copy of
  # the training rows: count the data twice (Ridge pool workers measured
  # 0.69 GB at 0.27 GB of data). With 1x, the larger pools allowed by the
  # 0.7 budget fraction failed to start RF workers on mice (Errno 22).
  per_job <- 0.5 + 2 * max(as.numeric(data_gb), 0)
  mem_cap <- if (is.finite(budget)) max(1L, as.integer(floor(budget / per_job))) else .Machine$integer.max
  max(1L, min(cores, mem_cap, n_bootstrap))
}

gp_parallel_seed_value <- function(seed) {
  if (is.null(seed) || identical(seed, FALSE)) {
    return(NULL)
  }
  value <- if (identical(seed, TRUE)) {
    suppressWarnings(as.numeric(Sys.getenv("GP_RANDOM_SEED", "123")))
  } else {
    suppressWarnings(as.numeric(seed)[1L])
  }
  if (!is.finite(value)) {
    stop("seed must be FALSE, TRUE, or one finite numeric value.", call. = FALSE)
  }
  as.integer(value %% .Machine$integer.max)
}

gp_parallel_seed_streams <- function(n, seed = TRUE) {
  n <- as.integer(n)
  base_seed <- gp_parallel_seed_value(seed)
  if (is.null(base_seed) || !is.finite(n) || n <= 0L) {
    return(NULL)
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  restore_seed <- had_seed && is.integer(old_seed) && length(old_seed) > 1L
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (restore_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  gp_set_seed(base_seed)
  streams <- vector("list", n)
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (n > 1L) {
    for (i in 2:n) {
      streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
    }
  }
  streams
}

gp_mirai_required_versions <- function() {
  c(mirai = "2.7.0", nanonext = "1.9.0")
}

gp_mirai_check_runtime <- function() {
  required <- gp_mirai_required_versions()
  missing <- names(required)[!vapply(names(required), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf(
        "The mirai backend requires %s. Install the current stable async stack before using backend='mirai'.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  outdated <- names(required)[vapply(names(required), function(pkg) {
    utils::packageVersion(pkg) < package_version(required[[pkg]])
  }, logical(1))]
  if (length(outdated)) {
    details <- paste(
      sprintf("%s >= %s", outdated, unname(required[outdated])),
      collapse = ", "
    )
    stop(
      sprintf(
        "The mirai backend requires current production async packages: %s.",
        details
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

gp_mirai_queue_memory_mb <- function(workers = 1L, n_tasks = 1L) {
  configured <- suppressWarnings(as.numeric(Sys.getenv("GP_MIRAI_QUEUE_MEMORY_MB", "")))
  if (is.finite(configured) && configured > 0) {
    return(as.numeric(ceiling(configured)))
  }
  workers <- max(1L, as.integer(workers %||% 1L))
  n_tasks <- max(1L, as.integer(n_tasks %||% 1L))
  active_fanout <- min(workers, n_tasks)
  as.numeric(max(256L, min(4096L, 128L * active_fanout)))
}

gp_mirai_queue_memory_bytes <- function(workers = 1L, n_tasks = 1L) {
  gp_mirai_queue_memory_mb(workers = workers, n_tasks = n_tasks) * 1024^2
}

gp_mirai_submit_timeout_sec <- function() {
  timeout <- suppressWarnings(as.numeric(Sys.getenv("GP_MIRAI_SUBMIT_TIMEOUT_SEC", "600")))
  if (!is.finite(timeout) || timeout <= 0) {
    return(600)
  }
  timeout
}

gp_mirai_submit_sleep_sec <- function() {
  sleep <- suppressWarnings(as.numeric(Sys.getenv("GP_MIRAI_SUBMIT_SLEEP_SEC", "0.02")))
  if (!is.finite(sleep) || sleep < 0) {
    return(0.02)
  }
  sleep
}

gp_mirai_start_daemons <- function(n, .compute = NULL, .expr = NULL, sync = FALSE, memory = NULL) {
  mirai::daemons(
    n = n,
    .compute = .compute,
    .expr = .expr,
    sync = sync,
    memory = memory
  )
}

gp_mirai_submit <- function(.expr, .args = list(), .compute = NULL) {
  expr <- substitute(.expr)
  do.call(
    mirai::try_mirai,
    list(
      .expr = expr,
      .args = .args,
      .compute = .compute
    )
  )
}

gp_mirai_collect <- function(job) {
  mirai::collect_mirai(job)
}

gp_mirai_is_error_value <- function(x) {
  isTRUE(mirai::is_error_value(x))
}

gp_mirai_error_message <- function(x) {
  tryCatch(conditionMessage(x), error = function(e) paste(as.character(x), collapse = " "))
}

gp_mirai_queue_depth <- function() {
  tryCatch(mirai::status()$daemons$tasks, error = function(e) NA_integer_)
}

gp_future_globals_limit_bytes <- function(decision,
                                          current_option = getOption(
                                            "future.globals.maxSize", NULL
                                          )) {
  if (!is.null(current_option)) {
    value <- suppressWarnings(as.numeric(current_option)[1L])
    return(list(
      bytes = value,
      source = "user_option",
      temporary = FALSE
    ))
  }

  budget_gb <- suppressWarnings(as.numeric(decision$memory_budget_gb)[1L])
  workers <- suppressWarnings(as.integer(decision$workers)[1L])
  if (!is.finite(budget_gb) || budget_gb <= 0 ||
      !is.finite(workers) || workers < 1L) {
    return(list(
      bytes = 500 * 1024^2,
      source = "future_default",
      temporary = FALSE
    ))
  }

  # future's default 500 MiB scan is independent of PredictProR's memory
  # policy. Use the already-budgeted memory available to each worker so a
  # closure that exposes shared objects more than once to future's recursive
  # global scan is not rejected before execution. Explicit user options always
  # win and are never changed here.
  list(
    bytes = max(1024^2, budget_gb * 1024^3 / workers),
    source = "predictpror_per_worker_memory_budget",
    temporary = TRUE
  )
}

sp_apply <- function(X, FUN, decision, packages = NULL, seed = TRUE, initializer = NULL, fail_fast = FALSE) {
  backend <- decision$backend
  packages <- unique(c(packages %||% character(), "reticulate", "logger"))
  runtime_env <- gp_parallel_runtime_env()
  # workers learn how many siblings share the machine (nested bootstrap pools)
  runtime_env[["PREDICTPRO_PARALLEL_WORKERS"]] <- as.character(max(1L, as.integer(decision$workers %||% 1L)))
  # ...and the run's memory budget, so a worker's nested pool sizes itself
  # from the budget this run was given (not from memory free at that moment,
  # which the sibling workers are already using).
  run_budget <- suppressWarnings(as.numeric(decision$memory_budget_gb)[1L])
  if (length(run_budget) && is.finite(run_budget) && run_budget > 0) {
    runtime_env[["GP_PAR_MEMORY_BUDGET_GB"]] <- format(run_budget, scientific = FALSE)
  }
  old_runtime_env <- Sys.getenv(names(runtime_env), unset = NA_character_)
  gp_parallel_set_runtime_env(runtime_env)
  on.exit(gp_parallel_restore_runtime_env(old_runtime_env), add = TRUE)
  # future/parallelly workers re-prepend this session's .libPaths(); include the
  # library PredictProR was loaded from for the duration of this call.
  loaded_lib <- gp_predictpror_loaded_lib()
  if (!is.null(loaded_lib) && !loaded_lib %in% normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)) {
    old_lib_paths <- .libPaths()
    .libPaths(c(loaded_lib, old_lib_paths))
    on.exit(.libPaths(old_lib_paths), add = TRUE)
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # In-process backends run tasks in the caller's session, and each task
  # installs an L'Ecuyer-CMRG stream seed; restore the caller's RNG kind and
  # seed afterwards so model fitting does not change the user's random stream.
  caller_rng_kind <- RNGkind()
  caller_had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  caller_seed <- if (caller_had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(caller_rng_kind))
    if (caller_had_seed) {
      assign(".Random.seed", caller_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  task_seed_streams <- gp_parallel_seed_streams(length(X), seed = seed)
  if (!is.null(task_seed_streams)) {
    task_fun <- FUN
    X <- lapply(seq_along(X), function(i) {
      list(value = X[[i]], rng_seed = task_seed_streams[[i]])
    })
    FUN <- function(task) {
      assign(".Random.seed", task$rng_seed, envir = .GlobalEnv)
      task_fun(task$value)
    }
  }

  if (identical(backend, "sequential")) {
    logger::log_info("Running {length(X)} tasks sequentially")
    gp_parallel_init_runtime(runtime_env)
    if (is.function(initializer)) initializer()
    return(lapply(X, FUN))
  }

  if (identical(backend, "base_parallel")) {
    sys_name <- Sys.info()[["sysname"]]
    use_fork <- sys_name != "Windows" &&
      .Platform$OS.type == "unix" &&
      !identical(decision$prefer_fork, FALSE)
    # Task costs differ a lot (one RandomForest fold can cost 100 Ridge
    # folds). Run the most expensive tasks first and hand tasks out
    # dynamically; static contiguous chunks left workers idle.
    task_order <- seq_along(X)
    task_seconds <- suppressWarnings(as.numeric(decision$task_seconds))
    if (length(task_seconds) == length(X) && all(is.finite(task_seconds))) {
      task_order <- order(-task_seconds)
    }
    X_run <- X[task_order]
    restore_order <- function(res) {
      out <- vector("list", length(X))
      out[task_order] <- res
      out
    }
    if (isTRUE(use_fork)) {
      logger::log_info("Running {length(X)} tasks with load-balanced mclapply on {decision$workers} workers")
      return(restore_order(parallel::mclapply(X_run, function(x) {
        gp_parallel_init_worker_runtime(runtime_env)
        if (is.function(initializer)) initializer()
        FUN(x)
      }, mc.cores = max(1L, decision$workers), mc.preschedule = FALSE,
      mc.set.seed = !is.null(task_seed_streams))))
    } else {
      logger::log_info("Running {length(X)} tasks with load-balanced parLapply on {decision$workers} workers")
      cl <- gp_make_psock_cluster(max(1L, decision$workers))
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      parallel::clusterCall(cl, gp_parallel_init_worker_runtime, runtime_env)
      if (length(packages)) {
        parallel::clusterCall(cl, function(pkgs) lapply(pkgs, requireNamespace, quietly = TRUE), packages)
      }
      if (is.function(initializer)) parallel::clusterCall(cl, initializer)
      # Ship the task function (and the shared data it closes over) once per
      # worker; each task then sends only a tiny wrapper and its index.
      parallel::clusterCall(cl, function(f) {
        assign(".gp_sp_task_fun", f, envir = globalenv())
        invisible(NULL)
      }, FUN)
      task_wrapper <- function(x) get(".gp_sp_task_fun", envir = globalenv())(x)
      environment(task_wrapper) <- globalenv()
      return(restore_order(parallel::parLapplyLB(cl, X_run, task_wrapper, chunk.size = 1L)))
    }
  }

  if (identical(backend, "foreach")) {
    if (!requireNamespace("foreach", quietly = TRUE)) {
      logger::log_error("foreach package required for backend='foreach'")
      stop("foreach package required")
    }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      logger::log_error("doParallel package required for backend='foreach'")
      stop("doParallel package required")
    }
    logger::log_info("Running {length(X)} tasks with foreach on {decision$workers} workers")
    cl <- gp_make_psock_cluster(max(1L, decision$workers))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterCall(cl, gp_parallel_init_worker_runtime, runtime_env)
    doParallel::registerDoParallel(cl)
    if (is.function(initializer)) parallel::clusterCall(cl, initializer)
    `%dopar%` <- foreach::`%dopar%`
    return(
      foreach::foreach(i = X, .packages = packages) %dopar% {
        FUN(i)
      }
    )
  }

  if (identical(backend, "mirai")) {
    if (!requireNamespace("mirai", quietly = TRUE)) {
      logger::log_error("mirai package required for backend='mirai'")
      stop("mirai package required")
    }
    gp_mirai_check_runtime()
    if (!requireNamespace("purrr", quietly = TRUE)) {
      logger::log_error("purrr (>=1.1.0) package required for backend='mirai'")
      stop("purrr package required")
    }
    has_progressr <- requireNamespace("progressr", quietly = TRUE)

    # Start metrics endpoint
    metrics_server <- start_metrics_endpoint()
    on.exit({
      stop_metrics_endpoint(metrics_server)
      stop_daemons()
    }, add = TRUE)

    # Reset daemons
    stop_daemons()

    # Initialize metrics
    metrics_registry$initialize()

    # Setup daemons
    sync_startup <- !(tolower(Sys.getenv("GP_MIRAI_SYNC_STARTUP", "true")) %in% c("false", "0", "no"))
    num_gpus <- decision$num_gpus
    use_gpu <- isTRUE(decision$gpu_parallel_enabled) && num_gpus > 1L && decision$workers > 1L
    profiles <- "default"
    daemons <- list()
    queue_memory <- gp_mirai_queue_memory_bytes(workers = decision$workers, n_tasks = length(X))

    if (use_gpu) {
      profiles <- sprintf("gpu%d", seq_along(decision$gpu_ids))
      for (i in seq_along(profiles)) {
        load_expr <- substitute({
          gp_parallel_init_worker_runtime(env_map)
          Sys.setenv(CUDA_VISIBLE_DEVICES = cuda_id)
          lapply(pkgs, requireNamespace, quietly = TRUE)
        }, list(pkgs = packages, cuda_id = Sys.getenv("GP_CUDA_VISIBLE_DEVICES", decision$gpu_ids[i]), env_map = runtime_env))
        daemons[[profiles[i]]] <- tryCatch({
          gp_mirai_start_daemons(n = 1L, .compute = profiles[i], .expr = load_expr, sync = sync_startup, memory = queue_memory)
          logger::log_info("Started daemon for GPU {decision$gpu_ids[i]}")
          1L
        }, error = function(e) {
          logger::log_error("Failed to start daemon for GPU {decision$gpu_ids[i]}: {conditionMessage(e)}")
          stop("Daemon startup failed")
        })
      }
    } else {
      load_expr <- substitute({
        gp_parallel_init_worker_runtime(env_map)
        lapply(pkgs, requireNamespace, quietly = TRUE)
      }, list(pkgs = packages, env_map = runtime_env))
      daemons[["default"]] <- tryCatch({
        gp_mirai_start_daemons(n = max(1L, decision$workers), .expr = load_expr, sync = sync_startup, memory = queue_memory)
        logger::log_info("Started {decision$workers} default daemons with bounded queue memory {round(queue_memory / 1024^2)} MB")
        decision$workers
      }, error = function(e) {
        logger::log_error("Failed to start default daemons: {conditionMessage(e)}")
        stop("Daemon startup failed")
      })
    }

    # Wrap FUN with initializer
    innerFUN <- if (is.function(initializer)) {
      function(x) {
        initializer()
        FUN(x)
      }
    } else FUN

    # Assign profiles round-robin
    task_profiles <- profiles[(seq_along(X) - 1L) %% length(profiles) + 1L]

    submit_mirai_one <- function(task_x, task_profile, task_index) {
      deadline <- Sys.time() + gp_mirai_submit_timeout_sec()
      sleep_base <- gp_mirai_submit_sleep_sec()
      attempts <- 0L
      last_error <- NULL
      repeat {
        attempts <- attempts + 1L
        job <- gp_mirai_submit(
          innerFUN(task_x),
          .args = list(innerFUN = innerFUN, task_x = task_x),
          .compute = task_profile
        )
        if (!gp_mirai_is_error_value(job)) {
          return(job)
        }
        last_error <- gp_mirai_error_message(job)
        if (Sys.time() >= deadline) {
          stop(
            sprintf(
              "mirai submission timed out for task %d after %d attempts; last dispatcher response: %s",
              task_index,
              attempts,
              last_error
            ),
            call. = FALSE
          )
        }
        Sys.sleep(min(1, sleep_base * max(1, attempts)))
      }
    }

    submit_mirai_jobs <- function() {
      submit_times <- vector("list", length(X))
      jobs <- purrr::map(seq_along(X), function(i) {
        task_x <- X[[i]]
        submit_times[[i]] <<- Sys.time()
        submit_mirai_one(task_x = task_x, task_profile = task_profiles[i], task_index = i)
      })
      list(jobs = jobs, submit_times = submit_times)
    }

    collect_mirai_job <- function(i, submitted, progressor = NULL) {
      result <- tryCatch({
        res <- gp_mirai_collect(submitted$jobs[[i]])
        if (gp_mirai_is_error_value(res)) {
          err_msg <- gp_mirai_error_message(res)
          stop(sprintf("mirai task returned error value: %s", err_msg), call. = FALSE)
        }
        if (is.function(progressor)) {
          progressor(message = sprintf("Task %d completed", i))
        }
        res
      }, error = function(e) {
        logger::log_error("Task {i} failed: {conditionMessage(e)}")
        if (fail_fast) stop("Task failed and fail_fast is TRUE")
        NULL
      })
      start_time <- submitted$submit_times[[i]] %||% Sys.time()
      metrics_registry$update(
        queue_depth = gp_mirai_queue_depth(),
        task_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        cpu = if (requireNamespace("ps", quietly = TRUE)) tryCatch(sum(ps::ps_cpu_times()$percent), error = function(e) NA_real_) else NA_real_,
        gpu = if (use_gpu) tryCatch(get_gpu_usage(), error = function(e) NA_real_) else NA_real_,
        mem = if (requireNamespace("ps", quietly = TRUE)) tryCatch(sum(ps::ps_memory_info()$rss) / 1024^2, error = function(e) NA_real_) else NA_real_
      )
      result
    }

    # Submit all work before collecting so the mirai worker pool can actually fan out.
    logger::log_info("Running {length(X)} tasks with mirai on {decision$workers} workers")
    res <- if (has_progressr) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(X))
        submitted <- submit_mirai_jobs()
        purrr::map(seq_along(X), function(i) {
          collect_mirai_job(i, submitted, progressor = p)
        })
      })
    } else {
      submitted <- submit_mirai_jobs()
      purrr::map(seq_along(X), function(i) {
        collect_mirai_job(i, submitted)
      })
    }

    return(res)
  }

  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    logger::log_error("future/future.apply packages required for backend='future'")
    stop("future/future.apply packages required")
  }
  logger::log_info("Running {length(X)} tasks with future on {decision$workers} workers")
  future_globals_limit <- gp_future_globals_limit_bytes(decision)
  if (isTRUE(future_globals_limit$temporary)) {
    options(future.globals.maxSize = future_globals_limit$bytes)
    on.exit(options(future.globals.maxSize = NULL), add = TRUE)
    logger::log_debug(
      "Temporary future globals limit set from the PredictProR per-worker memory budget: {round(future_globals_limit$bytes / 1024^2, 1)} MiB"
    )
  }
  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
  future::plan(decision$plan, workers = max(1L, decision$workers))
  innerFUN <- if (is.function(initializer)) {
    function(x) { gp_parallel_init_worker_runtime(runtime_env); initializer(); FUN(x) }
  } else FUN
  if (!is.function(initializer)) {
    innerFUN <- function(x) { gp_parallel_init_worker_runtime(runtime_env); FUN(x) }
  }
  future.apply::future_lapply(
    X, innerFUN,
    future.seed = !is.null(task_seed_streams),
    future.chunk.size = decision$chunk_size,
    future.packages = packages
  )
}

gp_sp_prepare_globals <- function(globals, decision) {
  if (is.null(globals)) {
    return(list(globals = globals, mori_summary = NULL))
  }
  share_res <- gp_mori_share_globals(
    globals = globals,
    backend = decision$backend %||% NULL,
    plan = decision$plan %||% NULL,
    workers = decision$workers %||% 1L,
    n_tasks = decision$n_tasks %||% decision$workers %||% 1L
  )
  gp_mori_log_summary(share_res$summary)
  list(globals = share_res$globals, mori_summary = share_res$summary)
}

gp_sp_bind_shared_globals <- function(fun, shared_globals) {
  if (!is.function(fun) || is.null(shared_globals) || !is.list(shared_globals) || !length(shared_globals)) {
    return(fun)
  }
  env <- environment(fun)
  if (is.null(env)) {
    return(fun)
  }
  for (nm in names(shared_globals)) {
    assign(nm, shared_globals[[nm]], envir = env)
  }
  fun
}
