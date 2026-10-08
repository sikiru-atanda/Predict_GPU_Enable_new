# Parallel-policy tuning knobs and decision-reason glossary.
#
# Pure-additive helpers introduced to address the many environment variables
# and options that tune the parallel policy, almost none documented,
# gap surfaced by the policy audit. Both functions return data.frames
# describing the policy code.

#' Tuning knobs for the parallel orchestration policy
#'
#' Returns a data.frame describing every environment variable and
#' `options()` key the parallel-policy engine reads, with default values
#' and a short effect summary. The list is the authoritative reference
#' used by `inst/parallel-policy-knobs.md`.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{name}{The environment variable or option key.}
#'     \item{type}{`"env"` for `Sys.getenv` lookups, `"option"` for
#'       `getOption()` lookups.}
#'     \item{default}{The default value used when the knob is unset.}
#'     \item{scope}{One of `"blas"`, `"scoring"`, `"sequential_gate"`,
#'       `"mori"`, `"mirai"`, `"gpu"`, `"reproducibility"`.}
#'     \item{effect}{One-line description of what the knob does.}
#'   }
#'
#' @details
#' The policy is **score-based**: each parallel backend (`sequential`,
#' `future`, `mirai`, `base_parallel`, `foreach`) gets a numeric score
#' built from object size, task fanout, OS, GPU availability, and
#' workload shape; the highest-scoring eligible backend wins. The
#' scoring weights, sequential-gate thresholds, mori share thresholds,
#' and mirai daemon timing are all tunable through the knobs in this
#' table.
#'
#' Companion helper: \code{\link{gp_policy_decision_reason_glossary}}.
#'
#' @seealso \code{\link{gp_policy_decision_reason_glossary}}
#' @export
gp_parallel_policy_knobs <- function() {
  data.frame(
    stringsAsFactors = FALSE,
    name = c(
      # BLAS / OpenMP thread pinning (applied in every worker)
      "GP_BLAS_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
      "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
      # Scoring weights (sp_decide_policy)
      "GP_PAR_SCORE_TASKS", "GP_PAR_SCORE_WORKERS", "GP_PAR_SCORE_GLOBALS",
      "GP_PAR_SCORE_WINDOWS_FUTURE", "GP_PAR_SCORE_WINDOWS_FUTURE_STARTUP",
      "GP_PAR_SCORE_WINDOWS_MIRAI_PENALTY", "GP_PAR_SCORE_MIRAI_WORKER_POOL",
      "GP_PAR_SCORE_GPU_MIRAI", "GP_PAR_SCORE_GPU_FUTURE_PENALTY",
      "GP_PAR_SCORE_HIDDEN_THREADS", "GP_PAR_SCORE_MORI",
      "GP_PAR_SCORE_TINY_SEQ", "GP_PAR_SCORE_THREADED_SEQ",
      "GP_PAR_SCORE_OVERHEAD", "GP_PAR_SCORE_WORKER_STARTUP",
      "GP_PAR_SCORE_FORK_WORKER_STARTUP", "GP_PAR_SCORE_WINDOWS_WORKER_STARTUP",
      "GP_PAR_SCORE_BASE_PARALLEL_UNIX", "GP_PAR_SCORE_FOREACH_PENALTY",
      # Sequential gates
      "GP_GLOBALS_MAX_GB", "GP_PAR_MEMORY_BUDGET_GB", "GP_PAR_MEMORY_FRACTION",
      "GP_PAR_PAYLOAD_MULTIPLIER", "GP_PAR_WORKER_OVERHEAD_GB",
      "gp.rkhs.mem.threshold.GB", "gp.always.seq.RKHS",
      # Mori shared-globals
      "GP_USE_MORI", "GP_MORI_MIN_MB", "GP_MORI_SCORE_THRESHOLD",
      "GP_MORI_SCORE_SIZE", "GP_MORI_SCORE_FANOUT", "GP_MORI_SCORE_REUSE",
      "GP_MORI_SCORE_COPYVOL", "GP_MORI_SCORE_OVERHEAD",
      # Mirai daemon
      "GP_MIRAI_QUEUE_MEMORY_MB", "GP_MIRAI_SUBMIT_TIMEOUT_SEC",
      "GP_MIRAI_SUBMIT_SLEEP_SEC", "GP_MIRAI_SYNC_STARTUP",
      # GPU
      "GP_CUDA_VISIBLE_DEVICES", "PREDICTPRO_GPU_BUSY_THRESHOLD",
      # Reproducibility
      "GP_RANDOM_SEED"
    ),
    type = c(
      rep("env", 6),    # BLAS
      rep("env", 19),   # scoring weights
      rep("env", 5), "option", "option",  # sequential gates
      rep("env", 8),    # mori
      rep("env", 4),    # mirai
      rep("env", 2),    # GPU
      "env"             # reproducibility
    ),
    default = c(
      "1", "(= GP_BLAS_THREADS)", "(= GP_BLAS_THREADS)", "(= GP_BLAS_THREADS)",
      "(= GP_BLAS_THREADS)", "(= GP_BLAS_THREADS)",
      "0.9", "0.5", "0.35", "0.9", "2.0", "0.45", "0.8", "3.0", "1.75", "0.4",
      "0.6", "2.2", "1.8", "0.45", "3.0", "1.2", "3.0", "0.9", "0.35",
      "4", "(auto: available memory x fraction)", "0.7", "4", "0.75", "1.0", "TRUE",
      "auto", "128", "1.0", "1.0", "0.75", "0.5", "0.35", "1.0",
      "(auto: 128 MB * fanout, capped 256-4096)", "600", "0.02", "true",
      "(0:max(0, num_gpus - 1))", "85",
      "123"
    ),
    scope = c(
      rep("blas", 6),
      rep("scoring", 19),
      rep("sequential_gate", 7),
      rep("mori", 8),
      rep("mirai", 4),
      rep("gpu", 2),
      "reproducibility"
    ),
    effect = c(
      # BLAS
      "Pins each worker's BLAS/OpenMP thread count. Default 1 prevents BLAS oversubscription when many workers run in parallel; see also init_single_thread_blas().",
      "Mirrors GP_BLAS_THREADS into OMP_NUM_THREADS so OpenMP-based BLAS libraries respect the cap.",
      "Mirrors GP_BLAS_THREADS into OPENBLAS_NUM_THREADS for OpenBLAS builds.",
      "Mirrors GP_BLAS_THREADS into MKL_NUM_THREADS for Intel MKL builds.",
      "Mirrors GP_BLAS_THREADS into BLIS_NUM_THREADS for BLIS builds.",
      "Mirrors GP_BLAS_THREADS into VECLIB_MAXIMUM_THREADS for Apple vecLib.",
      # Scoring
      "Weight on log1p(n_parallel_tasks). Larger -> more aggressive parallelism for big task queues.",
      "Weight on log1p(n_workers). Larger -> stronger preference for parallel when many cores are available.",
      "Weight on log1p(total_globals_GB). Larger -> stronger preference for parallel when globals are big (amortises copy cost).",
      "Future-backend bonus on Windows. Larger -> Windows future preferred more often.",
      "Penalty applied to future on Windows for PSOCK startup time. Larger -> stronger Windows-future avoidance.",
      "Penalty applied to mirai on Windows (socket overhead). Larger -> mirai less preferred on Windows.",
      "Bonus to mirai for reusable worker pool (queues > 2 tasks, no GPU). Larger -> stronger mirai preference for steady streams.",
      "Bonus to mirai (parallel path) and sequential (single-GPU path) when num_gpus > 1. Drives GPU-aware scheduling.",
      "Penalty to future when num_gpus > 1 (process-per-GPU overhead). Larger -> push GPU work to mirai instead.",
      "Penalty to parallel scores when internal-threaded models (xgboost/catboost/lightgbm/randomforest) are present.",
      "Bonus to parallel scores when a global object qualifies for mori sharing (size signal x weight).",
      "Bonus to sequential for tiny workloads (tasks <= 2 or tasks_per_worker <= 1). Larger -> more aggressive sequential on small queues.",
      "Bonus to sequential when thread-sensitive models would otherwise contend in parallel.",
      "Penalty applied to every parallel score for general parallelization overhead.",
      "Generic PSOCK worker startup penalty (applied to parallel backends with PSOCK transport).",
      "Fork-cluster startup penalty (Unix base_parallel). Lower than PSOCK because fork is cheap.",
      "Windows-specific worker startup penalty (PSOCK only on Windows; no fork).",
      "Bonus to base_parallel on Unix (cheap fork).",
      "Generic foreach penalty (extra dispatch overhead).",
      # Sequential gates
      "Total globals (GB) threshold beyond which the policy will not attempt parallelism (avoids OOM on parallel copies). Default 4 GB.",
      "Explicit memory budget (GB) available for additional worker processes on every OS. When unset, the policy derives a budget from currently available physical memory.",
      "Fraction of currently available physical memory reserved for worker processes when GP_PAR_MEMORY_BUDGET_GB is unset. Default 0.5.",
      "Multiplier applied to the serialized payload per worker to allow for transient copies and model transformations. Must be at least 1; default 1.25.",
      "Fixed GiB allowance per R worker for the interpreter, loaded packages, and runtime state. Default 0.25.",
      "RKHS kernel size (GB) threshold beyond which RKHS is forced sequential to avoid duplicated kernel copies across workers.",
      "When TRUE (default), every RKHS task is forced sequential regardless of size/thread heuristics. Set FALSE to let scoring decide.",
      # Mori
      "Master switch: 'auto'/'true' enables mori shared-object transport when the package is available; 'false'/'0'/'no'/'off' disables.",
      "Minimum object size (MB) before mori sharing is even considered. Smaller objects copy normally.",
      "Numeric threshold for mori scoring; objects scoring below this are skipped even when eligible.",
      "Weight on log1p(size_mb / min_mb) in the mori share score.",
      "Weight on log1p(max(fanout - 1, 0)) -- favours sharing when many workers fan out.",
      "Weight on log1p(max(reuse_ratio - 1, 0)) -- favours sharing when the same object is reused across many tasks.",
      "Weight on log1p(copy_volume_mb / min_mb) -- favours sharing when total copy bytes would be large.",
      "Penalty (1 / log1p(size_mb)) -- discourages sharing tiny objects where overhead dominates.",
      # Mirai
      "Per-daemon memory budget (MB). Empty (auto): max(256, min(4096, 128 * fanout)). Override for very large kernels.",
      "Timeout (seconds) waiting for mirai task submission before failing the fold/task.",
      "Sleep (seconds) between mirai submission attempts while waiting on backpressure.",
      "When TRUE (default), block until daemons report ready before submitting tasks. 'false'/'0'/'no' skips the sync.",
      # GPU
      "Colon-separated list of CUDA device IDs visible to workers. Default covers all detected GPUs.",
      "GPU utilisation percent threshold (0-100) above which the policy treats a GPU as busy and routes new GP tasks to CPU.",
      # Reproducibility
      "Base RNG seed used to derive deterministic L'Ecuyer-CMRG streams per task across sequential, future, mirai, base_parallel and foreach backends."
    )
  )
}

#' Glossary of parallel-policy decision-reason codes
#'
#' Returns the lookup table for the `policy_decision_reason` strings
#' that appear in `Run_metadata.csv` and the policy log line printed at
#' CV start (e.g. `reason=rkhs_hidden_threads;single_gpu_model`). Each
#' code maps to a short human-readable description.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{reason}{The literal string emitted by the policy.}
#'     \item{category}{Grouping (`user_override`, `sequential_driver`,
#'       `backend_choice`, `routing`, `bonus_penalty`, `mori`,
#'       `dispatch`).}
#'     \item{meaning}{One-line explanation of when the reason fires.}
#'   }
#'
#' @details
#' Reason codes are appended (and joined with `;`) as the policy makes
#' decisions, so a single run can carry several codes that together
#' explain the chosen backend, worker count, and any forced sequential
#' subset. This glossary covers every literal emitted from
#' `R/parallel_policy.R`, `R/cv_execution_helpers.R`, and
#' `R/main_crossvalidation_execution_logic.R`.
#'
#' Companion helper: \code{\link{gp_parallel_policy_knobs}}.
#'
#' @seealso \code{\link{gp_parallel_policy_knobs}}
#' @export
gp_policy_decision_reason_glossary <- function() {
  data.frame(
    stringsAsFactors = FALSE,
    reason = c(
      # User override
      "user_forced_sequential", "user_forced_future", "user_forced_mirai",
      "user_forced_base_parallel", "user_forced_foreach",
      # High-level sequential drivers
      "no_parallel_eligible_tasks", "single_gpu_exclusive_deep_models",
      "multiple_sequential_constraints", "all_tasks_user_sequential",
      "single_direct_gp_cv_task", "backend_unavailable", "memory_limited_workers",
      # Backend choice (top-level summary)
      "general_psock_backend", "low_overhead_worker_pool",
      "reusable_worker_pool", "base_parallel_fallback", "foreach_fallback",
      # Sequential tunings
      "tiny_workload", "thread_sensitive_parallel", "small_cpu_queue",
      "small_windows_cpu_queue", "single_gpu_model",
      # Penalty signals (folded into the chosen backend's score)
      "small_cpu_queue_penalty", "small_windows_cpu_queue_penalty",
      "windows_future_startup_penalty", "windows_psock",
      "windows_psock_penalty", "windows_socket_penalty",
      "small_queue_penalty", "gpu_penalty", "thread_sensitive_penalty",
      "memory_capped_workers",
      # Bonus signals
      "multi_gpu_model_bonus", "unix_fork_bonus", "mori_shared_globals",
      # Routing diagnostics
      "deep_gpu_selected", "deep_cpu_or_no_gpu",
      "gp_user_cpu", "gp_cpu_no_gpu", "gp_cpu_gpu_busy", "gp_gpu_selected",
      "xgboost_gpu_selected",
      # Routing constraints (force sequential at task-level)
      "xgboost_internal_threads", "catboost_internal_threads",
      "lightgbm_internal_threads", "randomforest_internal_jobs",
      "rkhs_hidden_threads", "rkhs_memory_heavy", "rkhs_forced_seq",
      "asreml_prepared", "single_gpu_exclusive"
    ),
    category = c(
      rep("user_override", 5),
      rep("sequential_driver", 7),
      rep("backend_choice", 5),
      rep("sequential_driver", 5),
      rep("bonus_penalty", 10),
      rep("bonus_penalty", 3),
      rep("routing", 7),
      rep("routing", 9)
    ),
    meaning = c(
      # User overrides
      "User passed parallel_mode = 'sequential'.",
      "User passed parallel_mode = 'future'.",
      "User passed parallel_mode = 'mirai'.",
      "User passed parallel_mode = 'base_parallel'.",
      "User passed parallel_mode = 'foreach'.",
      # Sequential drivers
      "All tasks were filtered out (e.g. all forced sequential by their internal_flags); no parallel work to do.",
      "All tasks request deep-learning models that are GPU-exclusive on a single-GPU host; can't fan out.",
      "Multiple distinct model classes each force sequential (no single 'sole blocker' to call out).",
      "Every task in the queue is a user-specified sequential model.",
      "Single GP CV task routed directly through the gp backend bridge, bypassing the partitioned scheduler.",
      "Requested backend's required package isn't installed; sequential is the only viable option.",
      "Available-memory fanout permits fewer than two workers, so the queue runs sequentially.",
      # Backend choice
      "future (PSOCK or multicore) was chosen as the parallel backend.",
      "mirai chosen for its reusable worker pool with low per-task overhead.",
      "mirai chosen; the queue is large enough that reusing the daemon pool pays off.",
      "base_parallel (mclapply / parLapply) chosen as the parallel backend.",
      "foreach chosen as the parallel backend (typically when nothing else was eligible).",
      # Sequential tunings
      "tasks <= 2 or tasks_per_worker <= 1; sequential wins on small queues.",
      "Models with internal threading (xgboost/catboost/lightgbm/rf with n_threads/n_jobs > 1) push toward sequential.",
      "CPU queue is small relative to startup cost; sequential wins.",
      "Same as small_cpu_queue but with the Windows-specific startup cost.",
      "Task uses a GPU-exclusive model and only one GPU is available -> sequential.",
      # Penalties
      "Startup penalty appended to a parallel backend's score when the CPU queue is small.",
      "Windows-specific startup penalty added to a parallel backend's score.",
      "Windows future startup penalty (PSOCK cluster spin-up).",
      "future on Windows uses PSOCK (no fork available).",
      "base_parallel on Windows uses PSOCK with the PSOCK penalty.",
      "mirai on Windows pays a socket-overhead penalty vs. Unix.",
      "Penalty applied to mirai when only 2 or fewer tasks are queued.",
      "future penalised when num_gpus > 1 (process-per-GPU overhead).",
      "future or sequential score adjustment when thread-sensitive models are present.",
      "Requested worker count was reduced so process globals fit the configured or detected memory budget.",
      # Bonuses
      "mirai bonus when num_gpus > 1 (daemon-per-GPU profile available).",
      "base_parallel bonus on Unix (fork is cheap).",
      "Bonus applied (per backend) when a global qualifies for mori shared transport.",
      # Routing diagnostics
      "Deep-learning model routed to GPU.",
      "Deep-learning model routed to CPU (no GPU or device='cpu').",
      "Gaussian-process model routed to CPU because user set device='cpu'.",
      "GP routed to CPU because no GPU was detected.",
      "GP routed to CPU because the GPU is busy (utilisation above PREDICTPRO_GPU_BUSY_THRESHOLD).",
      "GP routed to GPU.",
      "Xgboost routed to GPU (booster != 'gblinear' and device != 'cpu').",
      # Routing constraints
      "Xgboost forced sequential because xgb_nthread > 1.",
      "CatBoost forced sequential because catboost_thread_count > 1.",
      "LightGBM forced sequential because lightgbm_nthread > 1.",
      "RandomForest forced sequential because rf_n_jobs > 1.",
      "RKHS forced sequential because detected hidden BLAS threads > 1.",
      "RKHS forced sequential because the kernel exceeds gp.rkhs.mem.threshold.GB.",
      "RKHS forced sequential by the gp.always.seq.RKHS option (TRUE by default).",
      "ASReml prep detected in globals; forced sequential because asreml is not parallel-safe inside R workers.",
      "GPU-exclusive model on a single-GPU host -> forced sequential for that task."
    )
  )
}
