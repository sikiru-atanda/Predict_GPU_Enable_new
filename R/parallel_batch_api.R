# Public shared-data and multi-job capacity planning APIs.

gp_parallel_shared_require_mori <- function() {
  if (!requireNamespace("mori", quietly = TRUE) ||
      utils::packageVersion("mori") < package_version("0.2.0")) {
    stop("Shared job data requires mori >= 0.2.0.", call. = FALSE)
  }
  invisible(TRUE)
}

gp_parallel_shared_reference_fields <- function(reference) {
  required <- c(
    "name", "bytes", "size_mib", "created_pid", "created_at",
    "immutable", "direct_final_worker_map", "schema_version"
  )
  if (!inherits(reference, "predictpror_shared_reference") ||
      !is.list(reference) || !all(required %in% names(reference))) {
    stop(
      "reference must be created by parallel_shared_reference().",
      call. = FALSE
    )
  }
  if (length(reference$name) != 1L || is.na(reference$name) ||
      !nzchar(reference$name)) {
    stop("The shared-data reference has no valid region name.", call. = FALSE)
  }
  bytes <- suppressWarnings(as.numeric(reference$bytes)[1L])
  if (!is.finite(bytes) || bytes < 0) {
    stop("The shared-data reference has an invalid byte size.", call. = FALSE)
  }
  reference
}

#' Create one owned shared-memory data region
#'
#' Converts an immutable R object to a `mori` shared-memory object and returns
#' an owner that must remain reachable until every consuming job has finished.
#' Send only the result of [parallel_shared_reference()] to other R processes;
#' never serialize the owner or an already mapped data object through another
#' Windows PSOCK boundary.
#'
#' @param x An atomic vector, matrix, data frame, or supported list that will
#'   remain immutable while jobs are using it.
#'
#' @return A `predictpror_shared_owner` list containing concise metadata and the
#'   owning shared object. Keep this object alive in the supervisor.
#' @rdname parallel_shared_data
#' @export
parallel_shared_create <- function(x) {
  gp_parallel_shared_require_mori()
  if (is.null(x)) stop("x cannot be NULL.", call. = FALSE)

  bytes <- as.numeric(utils::object.size(x))
  shared <- tryCatch(mori::share(x), error = function(e) e)
  if (inherits(shared, "error")) {
    stop("Could not create shared data: ", conditionMessage(shared), call. = FALSE)
  }
  if (!isTRUE(mori::is_shared(shared))) {
    stop(
      "x is not a data type supported by mori shared memory.",
      call. = FALSE
    )
  }

  reference <- structure(list(
    name = as.character(mori::shared_name(shared)),
    bytes = bytes,
    size_mib = bytes / 1024^2,
    created_pid = Sys.getpid(),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
    immutable = TRUE,
    direct_final_worker_map = TRUE,
    schema_version = 1L
  ), class = c("predictpror_shared_reference", "list"))

  structure(list(
    reference = reference,
    data = shared
  ), class = c("predictpror_shared_owner", "list"))
}

#' Create a lightweight cross-job shared-data reference
#'
#' Removes the owning data object and returns only the metadata required for a
#' final worker to map the supervisor-owned region. The returned reference is
#' intentionally small and safe to serialize.
#'
#' @param owner A `predictpror_shared_owner` returned by
#'   [parallel_shared_create()], or an existing shared reference.
#'
#' @return A `predictpror_shared_reference`. The reference does not keep the
#'   region alive; its owner must remain reachable in the supervisor.
#' @rdname parallel_shared_data
#' @export
parallel_shared_reference <- function(owner) {
  if (inherits(owner, "predictpror_shared_reference")) {
    return(gp_parallel_shared_reference_fields(owner))
  }
  if (!inherits(owner, "predictpror_shared_owner") ||
      !is.list(owner) || is.null(owner$reference) || is.null(owner$data)) {
    stop(
      "owner must be created by parallel_shared_create().",
      call. = FALSE
    )
  }
  gp_parallel_shared_require_mori()
  reference <- gp_parallel_shared_reference_fields(owner$reference)
  if (!isTRUE(mori::is_shared(owner$data)) ||
      !identical(as.character(mori::shared_name(owner$data)), reference$name)) {
    stop("The owner no longer matches its shared-data reference.", call. = FALSE)
  }
  reference
}

#' Map supervisor-owned data in the final worker
#'
#' Opens the shared-memory region described by a lightweight reference. Call
#' this function in the process that will actually read the data. Mapping in an
#' intermediate job and forwarding the mapped object through a second Windows
#' PSOCK boundary can materialize a full copy and is deliberately rejected when
#' an owner is supplied.
#'
#' @param reference A `predictpror_shared_reference` created by
#'   [parallel_shared_reference()].
#'
#' @return The shared, immutable R object backed by the supervisor's region.
#' @rdname parallel_shared_data
#' @export
parallel_shared_map <- function(reference) {
  if (inherits(reference, "predictpror_shared_owner")) {
    stop(
      "Do not send the shared owner to a worker. Send parallel_shared_reference(owner) and map it only in the final worker.",
      call. = FALSE
    )
  }
  gp_parallel_shared_require_mori()
  reference <- gp_parallel_shared_reference_fields(reference)
  mapped <- tryCatch(mori::map_shared(reference$name), error = function(e) e)
  if (inherits(mapped, "error")) {
    stop(
      "The shared-data region is unavailable. Keep its supervisor owner alive until all jobs finish: ",
      conditionMessage(mapped),
      call. = FALSE
    )
  }
  if (is.null(mapped) || !isTRUE(mori::is_shared(mapped)) ||
      !identical(as.character(mori::shared_name(mapped)), reference$name)) {
    stop(
      "The shared-data region could not be mapped. Keep its supervisor owner alive until all jobs finish.",
      call. = FALSE
    )
  }
  mapped
}

#' Inspect a shared-data owner or reference
#'
#' @param x A shared-data owner or reference.
#'
#' @return A one-row data frame describing accessibility and lifecycle rules.
#' @rdname parallel_shared_data
#' @export
parallel_shared_status <- function(x) {
  gp_parallel_shared_require_mori()
  is_owner <- inherits(x, "predictpror_shared_owner")
  reference <- if (is_owner) parallel_shared_reference(x) else
    gp_parallel_shared_reference_fields(x)
  mapped <- if (is_owner) {
    x$data
  } else {
    tryCatch(mori::map_shared(reference$name), error = function(e) NULL)
  }
  accessible <- !is.null(mapped) && isTRUE(mori::is_shared(mapped)) &&
    identical(as.character(mori::shared_name(mapped)), reference$name)
  data.frame(
    name = reference$name,
    bytes = as.numeric(reference$bytes),
    size_mib = as.numeric(reference$size_mib),
    accessible = accessible,
    owner_held_here = is_owner,
    owner_pid = as.integer(reference$created_pid),
    current_pid = Sys.getpid(),
    immutable = isTRUE(reference$immutable),
    direct_final_worker_map = isTRUE(reference$direct_final_worker_map),
    schema_version = as.integer(reference$schema_version),
    stringsAsFactors = FALSE
  )
}

gp_parallel_batch_positive_integer <- function(x, name) {
  value <- suppressWarnings(as.numeric(x)[1L])
  if (!is.finite(value) || value < 1 || value != floor(value)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

gp_parallel_batch_object_bytes <- function(x) {
  if (inherits(x, "predictpror_shared_owner")) {
    return(as.numeric(parallel_shared_reference(x)$bytes))
  }
  if (inherits(x, "predictpror_shared_reference")) {
    return(as.numeric(gp_parallel_shared_reference_fields(x)$bytes))
  }
  as.numeric(utils::object.size(x))
}

gp_parallel_batch_components <- function(payload_gb,
                                         workers_per_job,
                                         data_mode,
                                         payload_multiplier,
                                         worker_overhead_gb) {
  process_overhead <- worker_overhead_gb
  if (identical(data_mode, "independent")) {
    return(list(
      fixed_gb = 0,
      per_job_gb = payload_gb +
        workers_per_job * payload_gb * payload_multiplier +
        (workers_per_job + 1L) * process_overhead
    ))
  }
  if (identical(data_mode, "per_job_shared")) {
    return(list(
      fixed_gb = 0,
      per_job_gb = payload_gb * (1 + payload_multiplier) +
        (workers_per_job + 1L) * process_overhead
    ))
  }
  list(
    fixed_gb = payload_gb * (1 + payload_multiplier) + process_overhead,
    per_job_gb = (workers_per_job + 1L) * process_overhead
  )
}

#' Plan explicit capacity for concurrent PredictProR jobs
#'
#' Computes a conservative, inspectable admission plan before jobs are
#' launched. The caller must explicitly choose whether each job owns ordinary
#' data, owns one per-job shared region, or maps one supervisor-owned region.
#' The function does not launch work or mutate environment variables.
#'
#' @param objects Named list of large phenotype, genotype, kernel, omic, cache,
#'   shared-owner, or shared-reference objects used by every job.
#' @param jobs Number of top-level jobs requested concurrently.
#' @param workers_per_job Worker processes requested by each admitted job.
#' @param tasks_per_job Number of tasks in each job. Worker count is capped by
#'   this value.
#' @param data_mode Required explicit choice: `"independent"`,
#'   `"per_job_shared"`, or `"supervisor_shared"`.
#' @param memory_budget_gb Total batch memory budget. Defaults to the current
#'   PredictProR parallel-memory budget.
#' @param total_worker_limit Maximum worker slots for the entire batch.
#'   Defaults to detected logical cores minus one.
#' @param payload_multiplier Allowance for transient payload copies and model
#'   transformations. Defaults to `GP_PAR_PAYLOAD_MULTIPLIER` or 4.
#' @param worker_overhead_gb Memory allowance for each R supervisor, job, or
#'   worker process. Defaults to `GP_PAR_WORKER_OVERHEAD_GB` or 0.25 GiB.
#'
#' @return A `predictpror_batch_plan` list with object inventory, one-row
#'   decision table, and per-admitted-job environment settings. Jobs beyond
#'   `jobs_admitted` must remain queued.
#' @export
parallel_batch_plan <- function(objects,
                                jobs,
                                workers_per_job,
                                tasks_per_job = workers_per_job,
                                data_mode = NULL,
                                memory_budget_gb = NULL,
                                total_worker_limit = NULL,
                                payload_multiplier = NULL,
                                worker_overhead_gb = NULL) {
  if (!is.list(objects)) stop("objects must be a named list.", call. = FALSE)
  if (is.null(names(objects))) names(objects) <- paste0("object_", seq_along(objects))
  if (is.null(data_mode) || length(data_mode) != 1L || is.na(data_mode)) {
    stop(
      "data_mode must be explicitly selected: 'independent', 'per_job_shared', or 'supervisor_shared'.",
      call. = FALSE
    )
  }
  data_mode <- match.arg(
    as.character(data_mode),
    c("independent", "per_job_shared", "supervisor_shared")
  )
  jobs_requested <- gp_parallel_batch_positive_integer(jobs, "jobs")
  workers_requested <- gp_parallel_batch_positive_integer(
    workers_per_job, "workers_per_job"
  )
  tasks_per_job <- gp_parallel_batch_positive_integer(tasks_per_job, "tasks_per_job")

  detected <- parallel::detectCores(logical = TRUE)
  if (!is.finite(detected) || detected < 1L) detected <- 1L
  worker_limit <- if (is.null(total_worker_limit)) {
    max(1L, as.integer(detected) - 1L)
  } else {
    gp_parallel_batch_positive_integer(total_worker_limit, "total_worker_limit")
  }
  budget <- if (is.null(memory_budget_gb)) gp_parallel_memory_budget_gb() else
    suppressWarnings(as.numeric(memory_budget_gb)[1L])
  if (!is.finite(budget) || budget <= 0) {
    stop("memory_budget_gb must be finite and positive.", call. = FALSE)
  }
  payload_multiplier <- payload_multiplier %||%
    suppressWarnings(as.numeric(Sys.getenv("GP_PAR_PAYLOAD_MULTIPLIER", "4")))
  worker_overhead_gb <- worker_overhead_gb %||%
    suppressWarnings(as.numeric(Sys.getenv("GP_PAR_WORKER_OVERHEAD_GB", "0.25")))
  if (!is.finite(payload_multiplier) || payload_multiplier < 1) {
    stop("payload_multiplier must be finite and at least 1.", call. = FALSE)
  }
  if (!is.finite(worker_overhead_gb) || worker_overhead_gb < 0) {
    stop("worker_overhead_gb must be finite and non-negative.", call. = FALSE)
  }
  if (!identical(data_mode, "independent")) gp_parallel_shared_require_mori()

  sizes <- vapply(objects, gp_parallel_batch_object_bytes, numeric(1L))
  inventory <- data.frame(
    name = names(objects),
    class = vapply(objects, function(x) paste(class(x), collapse = "/"),
                   character(1L)),
    bytes = sizes,
    size_mib = sizes / 1024^2,
    stringsAsFactors = FALSE
  )
  inventory <- inventory[order(inventory$bytes, decreasing = TRUE), , drop = FALSE]
  rownames(inventory) <- NULL
  payload_gb <- sum(sizes) / 1024^3

  workers_selected <- min(workers_requested, tasks_per_job, worker_limit)
  workers_before_memory <- workers_selected
  repeat {
    components <- gp_parallel_batch_components(
      payload_gb, workers_selected, data_mode,
      payload_multiplier, worker_overhead_gb
    )
    if (components$fixed_gb + components$per_job_gb <= budget + 1e-12 ||
        workers_selected <= 1L) break
    workers_selected <- workers_selected - 1L
  }
  minimum_required_gb <- components$fixed_gb + components$per_job_gb
  memory_job_cap <- if (minimum_required_gb > budget + 1e-12) {
    0L
  } else if (components$per_job_gb <= 0) {
    .Machine$integer.max
  } else {
    max(0L, as.integer(floor(
      (budget - components$fixed_gb + 1e-12) / components$per_job_gb
    )))
  }
  core_job_cap <- max(1L, as.integer(floor(worker_limit / workers_selected)))
  jobs_admitted <- min(jobs_requested, core_job_cap, memory_job_cap)
  jobs_queued <- jobs_requested - jobs_admitted
  total_slots <- jobs_admitted * workers_selected
  projected_gb <- if (jobs_admitted > 0L) {
    components$fixed_gb + jobs_admitted * components$per_job_gb
  } else {
    0
  }
  per_job_budget <- if (jobs_admitted > 0L) {
    max(0, (budget - components$fixed_gb) / jobs_admitted)
  } else {
    0
  }

  reasons <- character()
  if (workers_requested > tasks_per_job) reasons <- c(reasons, "workers_capped_by_tasks")
  if (min(workers_requested, tasks_per_job) > worker_limit) {
    reasons <- c(reasons, "workers_capped_by_total_limit")
  }
  if (workers_selected < workers_before_memory) {
    reasons <- c(reasons, "workers_capped_by_memory")
  }
  if (jobs_requested > core_job_cap) reasons <- c(reasons, "jobs_capped_by_total_limit")
  if (jobs_requested > memory_job_cap) reasons <- c(reasons, "jobs_capped_by_memory")
  if (!length(reasons)) reasons <- "all_requested_capacity_admitted"

  decision <- data.frame(
    data_mode = data_mode,
    jobs_requested = jobs_requested,
    jobs_admitted = jobs_admitted,
    jobs_queued = jobs_queued,
    workers_per_job_requested = workers_requested,
    workers_per_job_selected = workers_selected,
    tasks_per_job = tasks_per_job,
    total_worker_slots_requested = jobs_requested * workers_requested,
    total_worker_slots_admitted = total_slots,
    total_worker_limit = worker_limit,
    oversubscription_ratio = total_slots / worker_limit,
    payload_gb = payload_gb,
    memory_budget_gb = budget,
    fixed_shared_memory_gb = components$fixed_gb,
    incremental_memory_per_job_gb = components$per_job_gb,
    projected_batch_memory_gb = projected_gb,
    per_job_memory_budget_gb = per_job_budget,
    can_run = jobs_admitted > 0L,
    decision_reason = paste(unique(reasons), collapse = ";"),
    stringsAsFactors = FALSE
  )
  handoff <- switch(
    data_mode,
    independent = "ordinary_value_per_job",
    per_job_shared = "one_owner_per_job",
    supervisor_shared = "reference_only_final_worker_maps"
  )
  job_environment <- if (jobs_admitted > 0L) {
    data.frame(
      job_id = sprintf("job_%03d", seq_len(jobs_admitted)),
      GP_PAR_MEMORY_BUDGET_GB = format(per_job_budget, digits = 15,
                                       scientific = FALSE, trim = TRUE),
      GP_BLAS_THREADS = "1",
      GP_USE_MORI = if (identical(data_mode, "independent")) "false" else "true",
      data_handoff = handoff,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      job_id = character(), GP_PAR_MEMORY_BUDGET_GB = character(),
      GP_BLAS_THREADS = character(), GP_USE_MORI = character(),
      data_handoff = character(), stringsAsFactors = FALSE
    )
  }

  structure(list(
    inventory = inventory,
    decision = decision,
    job_environment = job_environment
  ), class = c("predictpror_batch_plan", "list"))
}

#' @export
print.predictpror_shared_owner <- function(x, ...) {
  ref <- x$reference
  cat("PredictProR shared-data owner\n")
  cat("  name: ", ref$name, "\n", sep = "")
  cat("  size: ", format(round(ref$size_mib, 3), nsmall = 3), " MiB\n", sep = "")
  cat("  owner pid: ", ref$created_pid, "\n", sep = "")
  cat("  rule: keep owner alive; send only parallel_shared_reference(owner)\n")
  invisible(x)
}

#' @export
print.predictpror_shared_reference <- function(x, ...) {
  cat("PredictProR shared-data reference\n")
  cat("  name: ", x$name, "\n", sep = "")
  cat("  size: ", format(round(x$size_mib, 3), nsmall = 3), " MiB\n", sep = "")
  cat("  rule: map only in the final worker\n")
  invisible(x)
}

#' @export
print.predictpror_batch_plan <- function(x, ...) {
  decision <- x$decision[1L, , drop = FALSE]
  cat("PredictProR parallel batch plan\n")
  cat("  data mode: ", decision$data_mode, "\n", sep = "")
  cat("  jobs: ", decision$jobs_admitted, " admitted / ",
      decision$jobs_requested, " requested; ", decision$jobs_queued,
      " queued\n", sep = "")
  cat("  workers per job: ", decision$workers_per_job_selected,
      "; total slots: ", decision$total_worker_slots_admitted, "/",
      decision$total_worker_limit, "\n", sep = "")
  cat("  projected memory: ",
      format(round(decision$projected_batch_memory_gb, 3), nsmall = 3),
      "/", format(round(decision$memory_budget_gb, 3), nsmall = 3),
      " GiB\n", sep = "")
  cat("  reason: ", decision$decision_reason, "\n", sep = "")
  invisible(x)
}
