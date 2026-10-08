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

gp_execute_partitioned_tasks <- function(model_ids,
                                         run_task,
                                         sequential_models = NULL,
                                         decision_args = list(),
                                         packages = NULL,
                                         seed = TRUE,
                                         initializer = NULL,
                                         fail_fast = FALSE,
                                         log_messages = list(
                                           seq_user = "Running {n} sequential tasks",
                                           seq_policy = "Running {n} additional sequential tasks",
                                           parallel = "Running {n} parallel tasks"
                                         )) {
  # BLAS/OpenMP pinning prevents oversubscription in PARALLEL WORKERS only.
  # Tasks run in the main process keep the user's thread settings: pinning
  # there leaked OMP_NUM_THREADS = 1 into the session and throttled models
  # kept sequential precisely so they can use internal threads.
  main_initializer <- initializer
  if (exists("compose_initializers", mode = "function") &&
      exists("init_single_thread_blas", mode = "function")) {
    initializer <- compose_initializers(init_single_thread_blas, initializer)
  }

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

  if (length(parts$idx_seq_user) > 0L) {
    logger::log_info(glue::glue(log_messages$seq_user, n = length(parts$idx_seq_user)))
    results[parts$idx_seq_user] <- lapply(parts$idx_seq_user, safe_run_task)
  }

  if (length(parts$idx_rest) == 0L) {
    num_gpus <- if (exists("detect_num_gpus", mode = "function")) {
      tryCatch(detect_num_gpus(), error = function(e) 0L)
    } else {
      0L
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

  # Size what workers really receive: run_task's enclosing frames.
  if (is.null(decision_args$payload_gb)) {
    decision_args$payload_gb <- tryCatch(gp_parallel_closure_payload_gb(run_task), error = function(e) NA_real_)
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

  if (length(idx_seq_policy) > 0L) {
    logger::log_info(glue::glue(log_messages$seq_policy, n = length(idx_seq_policy)))
    seq_task_indices <- parts$idx_rest[idx_seq_policy]
    results[seq_task_indices] <- lapply(seq_task_indices, safe_run_task)
  }

  if (length(idx_parallel) > 0L) {
    logger::log_info(glue::glue(log_messages$parallel, n = length(idx_parallel)))
    parallel_task_indices <- parts$idx_rest[idx_parallel]
    parallel_results <- sp_apply(
      X = parallel_task_indices,
      FUN = routed_run_task,
      decision = decision,
      packages = packages,
      seed = seed,
      initializer = initializer,
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
  filtered_results <- Filter(function(res) {
    ok <- !is.null(res$eval_metrics_reps) && !is.null(res$ypred_cv_Reps_all)
    if (!ok) {
      collapsed <<- c(collapsed, describe(res))
      return(FALSE)
    }
    df <- res$ypred_cv_Reps_all
    has_any <- "yhat" %in% names(df) && any(!is.na(df$yhat))
    if (!has_any) collapsed <<- c(collapsed, describe(res))
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
    warning(
      "Cross-validation dropped ", length(collapsed),
      " task(s) because all folds produced NA predictions ",
      "(check input contract: e.g. MET ML/DL requires a kernel; ",
      "see logger output for the underlying fold errors). Dropped: ",
      paste(unique(collapsed), collapse = "; "),
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
    nonfinite_metric = unique(nonfinite)
  )
  filtered_results
}
