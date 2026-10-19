gp_normalize_internal_flags <- function(internal_flags, n_tasks) {
  if (is.null(internal_flags) || length(internal_flags) == 0L) {
    return(rep(TRUE, n_tasks))
  }
  as.logical(rep_len(internal_flags, n_tasks))
}

gp_partition_task_indices <- function(model_ids, sequential_models = NULL) {
  all_idx <- seq_along(model_ids)
  seq_models <- unique(sequential_models %||% character())
  idx_seq_user <- if (length(seq_models) > 0L) which(model_ids %in% seq_models) else integer(0)
  idx_rest <- setdiff(all_idx, idx_seq_user)
  list(
    idx_seq_user = idx_seq_user,
    idx_rest = idx_rest
  )
}

gp_task_device_routes <- function(model_ids, decision, task_indices) {
  routes <- rep(NA_character_, length(model_ids))
  caps <- decision$model_capabilities
  if (!is.data.frame(caps) || !("device_route" %in% names(caps)) || !length(task_indices)) {
    return(routes)
  }
  n <- min(length(task_indices), nrow(caps))
  routes[task_indices[seq_len(n)]] <- as.character(caps$device_route[seq_len(n)])
  routes
}

gp_with_task_device_route <- function(task_index, device_routes, fun) {
  route <- device_routes[[task_index]] %||% NA_character_
  route <- as.character(route)[1L]
  if (is.na(route) || !nzchar(route)) {
    return(fun(task_index))
  }

  old <- Sys.getenv("PREDICTPRO_GP_DEVICE", unset = NA_character_)
  Sys.setenv(PREDICTPRO_GP_DEVICE = route)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("PREDICTPRO_GP_DEVICE")
    } else {
      Sys.setenv(PREDICTPRO_GP_DEVICE = old)
    }
  }, add = TRUE)
  fun(task_index)
}

# The task function workers receive. Built here, not inside
# gp_execute_partitioned_tasks(): a closure created there carries that frame,
# with the policy's `globals` list and its mori-prepared copy -- the run's data
# twice more (heavy CV: 314 MB serialised per worker, of which ~100 MB the
# data; sent to each PSOCK worker in turn, measured 2026-10-10).
gp_routed_task <- function(run_task, device_routes) {
  force(run_task)
  force(device_routes)
  function(i) gp_with_task_device_route(i, device_routes, run_task)
}

# A closure sends its frame to each PSOCK worker, and R writes an object once
# for every name or list element that holds it: the CV's frame held params
# and five of its elements under their own names too (52.3 MB, 17 MB of it
# twice; 3-trait run). Returns `funs` (a named list; NULLs pass through) with
# a copy of their frame in which every repeat of an object of `min_bytes` or
# more is a reference to the place first holding it (breadth first, so a
# name before a list element). Each returned function fills the references
# in on its first call. Make the copy right before sending it: a call in
# this process fills it in too.
gp_share_repeated_objects <- function(funs, min_bytes = 10 * 1024) {
  idx <- which(vapply(funs, is.function, logical(1L)))
  if (!length(idx)) return(funs)
  src <- environment(funs[[idx[[1L]]]])
  if (!is.environment(src) || isNamespace(src) || identical(src, globalenv()) ||
      identical(src, baseenv()) || identical(src, emptyenv())) {
    return(funs)
  }
  # every name: what the task needs was decided by the caller (the refit's
  # copy keeps the context fields it reads by name)
  copied <- gp_trim_closure_frame(funs[idx], fields = ls(src, all.names = TRUE))
  env <- environment(copied[[1L]])
  refs <- gp_mark_repeated_objects(env, min_bytes)
  if (length(refs)) {
    assign(".gp_shared_refs", refs, envir = env)
    copied <- lapply(copied, function(f) {
      if (identical(environment(f), env)) gp_shared_objects_entry(f) else f
    })
  }
  funs[idx] <- copied
  funs
}

# Replaces, in `env`, each repeat of an object by a marker and returns where
# each marker is and where its object stays: list(at = path, from = path), a
# path being a name followed by list positions.
gp_mark_repeated_objects <- function(env, min_bytes) {
  plain_list <- function(x) is.list(x) && !is.data.frame(x) && is.null(attr(x, "class", exact = TRUE))
  queue <- lapply(sort(ls(env, all.names = TRUE)), function(nm) list(path = list(nm), value = get(nm, envir = env)))
  kept <- list()        # by type and length: the objects seen first, with their paths
  refs <- list()
  while (length(queue)) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    x <- node$value
    if (is.null(x) || is.function(x) || is.environment(x) ||
        as.numeric(utils::object.size(x)) < min_bytes) next
    key <- paste(typeof(x), length(x))
    hit <- NULL
    for (k in kept[[key]]) if (identical(k$value, x)) { hit <- k$path; break }
    if (!is.null(hit)) {
      refs[[length(refs) + 1L]] <- list(at = node$path, from = hit)
      next
    }
    kept[[key]] <- c(kept[[key]], list(list(path = node$path, value = x)))
    if (plain_list(x)) {
      queue <- c(queue, lapply(seq_along(x), function(i) list(path = c(node$path, list(i)), value = x[[i]])))
    }
  }
  marker <- structure(list(), class = "gp_shared_ref")
  for (r in refs) gp_set_path(env, r$at, marker)
  refs
}

gp_get_path <- function(env, path) {
  x <- get(path[[1L]], envir = env, inherits = FALSE)
  if (length(path) > 1L) x <- x[[unlist(path[-1L])]]
  x
}

gp_set_path <- function(env, path, value) {
  if (length(path) == 1L) {
    assign(path[[1L]], value, envir = env)
    return(invisible())
  }
  x <- get(path[[1L]], envir = env, inherits = FALSE)
  x[[unlist(path[-1L])]] <- value
  assign(path[[1L]], x, envir = env)
  invisible()
}

# Fills the references gp_share_repeated_objects() left in `env`, once. A
# reference waits while its source still holds a marker: `z <- a`, with a
# repeat of `b` inside `a`, is found before that repeat.
gp_fill_shared_objects <- function(env) {
  pending <- get0(".gp_shared_refs", envir = env, inherits = FALSE)
  if (is.null(pending)) return(invisible())
  within <- function(path, inside) length(inside) >= length(path) && identical(inside[seq_along(path)], path)
  while (length(pending)) {
    ready <- vapply(pending, function(r) {
      !any(vapply(pending, function(o) within(r$from, o$at), logical(1L)))
    }, logical(1L))
    if (!any(ready)) stop("Shared task objects refer to each other.", call. = FALSE)
    for (r in pending[ready]) gp_set_path(env, r$at, gp_get_path(env, r$from))
    pending <- pending[!ready]
  }
  rm(".gp_shared_refs", envir = env)
  invisible()
}

# Built here, so the function holds only `f`.
gp_shared_objects_entry <- function(f) {
  force(f)
  function(...) {
    gp_fill_shared_objects(environment(f))
    f(...)
  }
}

gp_execute_partitioned_tasks <- function(model_ids,
                                         run_task,
                                         sequential_models = NULL,
                                         decision_args = list(),
                                         packages = NULL,
                                         seed = TRUE,
                                         initializer = NULL,
                                         fail_fast = FALSE,
                                         unsplit = NULL,
                                         log_messages = list(
                                           seq_user = "Running {n} sequential tasks",
                                           seq_policy = "Running {n} additional sequential tasks",
                                           parallel = "Running {n} parallel tasks"
                                         )) {
  # Every task gets the initializer as given; workers pin BLAS/OpenMP to one
  # thread themselves when they start (gp_parallel_init_worker_runtime()).
  main_initializer <- initializer

  parts <- gp_partition_task_indices(
    model_ids = model_ids,
    sequential_models = sequential_models
  )
  task_device_routes <- rep(NA_character_, length(model_ids))
  routed_run_task <- function(i) {
    gp_with_task_device_route(i, task_device_routes, run_task)
  }
  safe_run_task <- function(i) {
    if (is.function(main_initializer)) {
      main_initializer()
    }
    routed_run_task(i)
  }
  # Preserve the caller's task-table order across every execution partition.
  # Appending the user-sequential, policy-sequential, and parallel subsets
  # reorders interleaved task indices and lets the caller attach model/trait
  # labels to the wrong fitted result.  Preallocate by the original task index
  # and assign each partition back into its source positions.
  results <- vector("list", length(model_ids))

  # `unsplit` (optional): the same work as fewer, larger tasks, list(model_ids,
  # run_task). Cross-validation splits each model's task into one task per
  # fold so workers can share the folds; when nothing runs in parallel the
  # unsplit tasks run instead, since every fold task repeats its model's set-up
  # (10 fold tasks: 18.2 s vs 11.5 s unsplit, measured 2026-10-10).
  run_unsplit <- function(policy, routes = rep(NA_character_, length(model_ids))) {
    logger::log_info("Running {length(unsplit$model_ids)} sequential tasks (one per model; nothing runs in parallel)")
    unsplit_routes <- routes[match(unsplit$model_ids, model_ids)]
    out <- lapply(seq_along(unsplit$model_ids), function(i) {
      if (is.function(main_initializer)) main_initializer()
      gp_with_task_device_route(i, unsplit_routes, unsplit$run_task)
    })
    policy$decision_reason <- paste0(policy$decision_reason %||% "", ";fold_tasks_merged")
    policy$model_ids <- unsplit$model_ids
    policy$internal_flags <- rep(TRUE, length(unsplit$model_ids))
    policy$internal_sequential_count <- length(unsplit$model_ids)
    policy$parallel_count <- 0L
    attr(out, "gp_execution_policy") <- policy
    attr(out, "gp_unsplit") <- TRUE
    out
  }

  run_seq_user <- function() {
    if (length(parts$idx_seq_user) > 0L) {
      logger::log_info(glue::glue(log_messages$seq_user, n = length(parts$idx_seq_user)))
      results[parts$idx_seq_user] <<- lapply(parts$idx_seq_user, safe_run_task)
    }
  }
  # with `unsplit`, only after the policy decided (they may run unsplit)
  if (!is.list(unsplit)) run_seq_user()

  if (length(parts$idx_rest) == 0L) {
    num_gpus <- tryCatch(gp_parallel_num_gpus(model_ids), error = function(e) 0L)
    if (is.list(unsplit)) {
      return(run_unsplit(list(
        backend = "sequential", workers = 1L, plan = NULL, chunk_size = NULL, backend_score = Inf,
        decision_reason = "all_tasks_user_sequential", num_gpus = num_gpus,
        user_sequential_models = unique(sequential_models %||% character())
      )))
    }
    attr(results, "gp_execution_policy") <- list(
      backend = "sequential",
      workers = 1L,
      plan = NULL,
      chunk_size = NULL,
      backend_score = Inf,
      decision_reason = "all_tasks_user_sequential",
      num_gpus = num_gpus,
      model_ids = model_ids,
      user_sequential_models = unique(sequential_models %||% character()),
      internal_flags = rep(TRUE, length(model_ids)),
      internal_sequential_count = length(model_ids),
      parallel_count = 0L
    )
    return(results)
  }

  # Size what workers really receive: run_task's enclosing frames, each
  # object once (gp_share_repeated_objects(), as sent below).
  if (is.null(decision_args$payload_gb)) {
    decision_args$payload_gb <- tryCatch(
      gp_parallel_closure_payload_gb(gp_share_repeated_objects(list(run = run_task, init = main_initializer))$run),
      error = function(e) NA_real_
    )
  }
  decision <- do.call(
    sp_decide_policy,
    c(
      list(
        models = model_ids[parts$idx_rest],
        n_tasks = length(parts$idx_rest)
      ),
      decision_args
    )
  )

  shared_globals_info <- NULL
  if (!is.null(decision_args$globals)) {
    shared_globals_info <- gp_sp_prepare_globals(decision_args$globals, decision = decision)
    run_task <- gp_sp_bind_shared_globals(run_task, shared_globals_info$globals)
  }
  task_device_routes <- gp_task_device_routes(
    model_ids = model_ids,
    decision = decision,
    task_indices = parts$idx_rest
  )
  routed_run_task <- function(i) {
    gp_with_task_device_route(i, task_device_routes, run_task)
  }
  safe_run_task <- function(i) {
    if (is.function(main_initializer)) {
      main_initializer()
    }
    routed_run_task(i)
  }

  int_flags <- gp_normalize_internal_flags(
    internal_flags = decision$internal_flags,
    n_tasks = length(parts$idx_rest)
  )
  idx_seq_policy <- which(int_flags)
  idx_parallel <- which(!int_flags)

  if (is.list(unsplit)) {
    # nothing in parallel: no parallel tasks, or the sequential backend (the
    # cost model's choice when start-up outweighs the work)
    if (!length(idx_parallel) || identical(decision$backend, "sequential")) {
      return(run_unsplit(c(decision, list(
        mori_summary = shared_globals_info$mori_summary %||% NULL,
        user_sequential_models = unique(sequential_models %||% character())
      )), routes = task_device_routes))
    }
    run_seq_user()
  }

  if (length(idx_seq_policy) > 0L) {
    logger::log_info(glue::glue(log_messages$seq_policy, n = length(idx_seq_policy)))
    seq_task_indices <- parts$idx_rest[idx_seq_policy]
    results[seq_task_indices] <- lapply(seq_task_indices, safe_run_task)
  }

  if (length(idx_parallel) > 0L) {
    logger::log_info(glue::glue(log_messages$parallel, n = length(idx_parallel)))
    parallel_task_indices <- parts$idx_rest[idx_parallel]
    # Workers get each object of the task's frame once. Made here, after the
    # tasks this process ran: a call fills the copy in.
    sent <- gp_share_repeated_objects(list(run = run_task, init = main_initializer))
    # BLAS/OpenMP pinning prevents oversubscription in parallel workers only,
    # which pin themselves when they start (gp_parallel_init_worker_runtime(),
    # every backend). It is not added to the initializer: sp_apply() runs that
    # in this process when it picks the sequential backend, and the pin stayed
    # in the user's session (OMP_NUM_THREADS = 1 and the rest).
    worker_initializer <- sent$init
    parallel_results <- sp_apply(
      X = parallel_task_indices,
      FUN = gp_routed_task(sent$run, task_device_routes),
      decision = decision,
      packages = packages,
      seed = seed,
      initializer = worker_initializer,
      fail_fast = fail_fast
    )
    if (length(parallel_results) != length(parallel_task_indices)) {
      stop(
        "Parallel task execution returned a different number of results than task indices.",
        call. = FALSE
      )
    }
    results[parallel_task_indices] <- parallel_results
  }

  aligned_internal_flags <- rep(FALSE, length(model_ids))
  aligned_internal_flags[parts$idx_seq_user] <- TRUE
  aligned_internal_flags[parts$idx_rest] <- as.logical(
    rep_len(int_flags, length(parts$idx_rest))
  )

  attr(results, "gp_execution_policy") <- c(
    decision,
    list(
      mori_summary = shared_globals_info$mori_summary %||% NULL,
      model_ids = model_ids,
      user_sequential_models = unique(sequential_models %||% character()),
      internal_flags = aligned_internal_flags,
      internal_sequential_count = length(parts$idx_seq_user) + length(idx_seq_policy),
      parallel_count = length(idx_parallel)
    )
  )

  results
}

gp_filter_crossval_results <- function(results) {
  # Identify two categories of dropped tasks so users see them instead of a
  # silent "Cross-validation completed with N valid results" with N < input
  # count:
  #   * collapsed   -> all per-fold predictions came back NA (every fold
  #                    failed; typically a missing-kernel/feature setup that
  #                    bailed per-fold).
  #   * nonfinite   -> at least one fold succeeded but every requested metric
  #                    is non-finite (e.g. only one observed value left, all
  #                    predictions degenerate). One mathematically undefined
  #                    classification metric must not discard otherwise valid
  #                    accuracy/probability metrics and OOF predictions.
  describe <- function(res) {
    paste0(
      res$model_canonical %||% res$model %||% "?",
      "/", res$trait %||% "?",
      if (!is.null(res$rep)) paste0(" rep=", res$rep) else ""
    )
  }
  collapsed <- character()
  # the first fold error of each collapsed task (recorded by the CV engine), so
  # the warning and the final "no result" error can name the cause
  collapsed_errors <- list()
  note_collapsed <- function(res) {
    collapsed <<- c(collapsed, describe(res))
    err <- unique(as.character(res$fold_errors %||% character()))
    if (length(err)) collapsed_errors[[describe(res)]] <<- err[[1L]]
  }
  filtered_results <- Filter(function(res) {
    ok <- !is.null(res$eval_metrics_reps) && !is.null(res$ypred_cv_Reps_all)
    if (!ok) {
      note_collapsed(res)
      return(FALSE)
    }
    df <- res$ypred_cv_Reps_all
    has_any <- "yhat" %in% names(df) && any(!is.na(df$yhat))
    if (!has_any) note_collapsed(res)
    has_any
  }, results)

  nonfinite <- character()
  filtered_results <- Filter(function(res) {
    if (is.null(res$eval_metrics_reps)) {
      nonfinite <<- c(nonfinite, describe(res))
      return(FALSE)
    }
    metric_cols <- setdiff(
      names(res$eval_metrics_reps)[vapply(res$eval_metrics_reps, is.numeric, logical(1))],
      c("Rep", "feature_k")
    )
    if (!length(metric_cols)) {
      nonfinite <<- c(nonfinite, describe(res))
      return(FALSE)
    }
    keep <- any(is.finite(unlist(res$eval_metrics_reps[metric_cols], use.names = FALSE)))
    if (!keep) nonfinite <<- c(nonfinite, describe(res))
    keep
  }, filtered_results)

  if (length(collapsed)) {
    dropped_text <- vapply(unique(collapsed), function(task) {
      err <- collapsed_errors[[task]]
      if (is.null(err)) task else paste0(task, " (fold error: ", err, ")")
    }, character(1))
    warning(
      "Cross-validation dropped ", length(collapsed),
      " task(s) because all folds produced NA predictions ",
      "(check input contract: e.g. MET ML/DL requires a kernel; ",
      "see logger output for the underlying fold errors). Dropped: ",
      paste(dropped_text, collapse = "; "),
      call. = FALSE
    )
  }
  if (length(nonfinite)) {
    warning(
      "Cross-validation dropped ", length(nonfinite),
      " task(s) because every requested aggregated metric was non-finite. Dropped: ",
      paste(unique(nonfinite), collapse = "; "),
      call. = FALSE
    )
  }
  attr(filtered_results, "dropped_tasks") <- list(
    collapsed_all_na_predictions = unique(collapsed),
    nonfinite_metric = unique(nonfinite),
    fold_errors = collapsed_errors
  )
  filtered_results
}
