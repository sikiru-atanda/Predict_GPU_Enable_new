# Parallel-policy tuning knobs and decision-reason codes

This is the user-facing reference for tuning PredictProR's parallel
orchestration engine and for interpreting the `policy_decision_reason`
strings that land in `Run_metadata.csv`.

Two helpers expose the same information programmatically:

```r
# Every env var / option key the policy reads, with defaults and effect.
gp_parallel_policy_knobs()

# Every literal reason code the policy emits, grouped by category.
gp_policy_decision_reason_glossary()
```

The policy is **score-based**: each backend (`sequential`, `future`,
`mirai`, `base_parallel`, `foreach`) gets a numeric score built from
object size, task fanout, OS, GPU availability, and workload shape;
the highest-scoring *eligible* backend wins. Sequential is always
eligible. Tuning knobs adjust the weights, gates, and timeouts that
feed that score.

## Where to set things

- **Environment variables**: read once via `Sys.getenv()` in the same
  R session that calls `model_execute()`. Set them via `Sys.setenv()`
  in R, or in the shell / `.Renviron` before launching R.
- **Options**: read via `getOption()`. Set them via `options()` in R
  or in `.Rprofile`. Two options exist:
  `gp.rkhs.mem.threshold.GB`, `gp.always.seq.RKHS`.

## BLAS / OpenMP thread pinning

`init_single_thread_blas()` runs in every parallel worker and sets the
BLAS thread cap so each worker runs single-threaded BLAS. This avoids the
classic `n_workers * BLAS_threads` oversubscription where everything
slows down because threads contend for the same physical cores.

| Knob | Default | Effect |
|---|---|---|
| `GP_BLAS_THREADS` | `1` | Master cap; every other BLAS env var below mirrors it. |
| `OMP_NUM_THREADS` | (= `GP_BLAS_THREADS`) | OpenMP cap (most BLAS libs respect this). |
| `OPENBLAS_NUM_THREADS` | (= `GP_BLAS_THREADS`) | OpenBLAS-specific cap. |
| `MKL_NUM_THREADS` | (= `GP_BLAS_THREADS`) | Intel MKL cap. |
| `BLIS_NUM_THREADS` | (= `GP_BLAS_THREADS`) | BLIS cap. |
| `VECLIB_MAXIMUM_THREADS` | (= `GP_BLAS_THREADS`) | Apple vecLib cap. |

Note: thread-sensitive models (RKHS, Xgboost with `xgb_nthread > 1`,
CatBoost with `catboost_thread_count > 1`, LightGBM with
`lightgbm_nthread > 1`, RandomForest with `rf_n_jobs > 1`) are forced
sequential so they can use their internal threading without
contending with worker parallelism.

On Unix-like systems, `parallel_backend_prefer_fork = FALSE` selects
process-spawned workers instead of forked workers. PredictProR also disables
forking automatically when a model batch requires per-process runtime
initialization, such as an embedded Python session.

## Scoring weights

These weights drive the score computed in `parallel_policy.R` for
each backend. All defaults are tuned for typical genomic prediction
workloads on developer-class hardware (8-16 logical cores, 32-64 GB
RAM, single GPU optional). The score itself is logged at INFO when
verbose is on; reason codes for each contribution land in
`policy_decision_reason`.

| Knob | Default | Backend(s) | Effect |
|---|---|---|---|
| `GP_PAR_SCORE_TASKS` | `0.9` | future/mirai/base_parallel/foreach | Weight on `log1p(n_parallel_tasks)`. |
| `GP_PAR_SCORE_WORKERS` | `0.5` | parallel | Weight on `log1p(n_workers)`. |
| `GP_PAR_SCORE_GLOBALS` | `0.35` | parallel | Weight on `log1p(total_globals_GB)`. |
| `GP_PAR_SCORE_WINDOWS_FUTURE` | `0.9` | future | Windows bonus to future. |
| `GP_PAR_SCORE_WINDOWS_FUTURE_STARTUP` | `2.0` | future | Windows future startup penalty. |
| `GP_PAR_SCORE_WINDOWS_MIRAI_PENALTY` | `0.45` | mirai | Windows socket-overhead penalty. |
| `GP_PAR_SCORE_MIRAI_WORKER_POOL` | `0.8` | mirai | Reusable-pool bonus on larger queues. |
| `GP_PAR_SCORE_GPU_MIRAI` | `3.0` | mirai (multi-GPU); seq (single-GPU) | GPU-aware bonus. |
| `GP_PAR_SCORE_GPU_FUTURE_PENALTY` | `1.75` | future | Multi-GPU process-per-GPU overhead. |
| `GP_PAR_SCORE_HIDDEN_THREADS` | `0.4` | parallel | Penalty when internal-thread models present. |
| `GP_PAR_SCORE_MORI` | `0.6` | parallel | Bonus when mori sharing applies. |
| `GP_PAR_SCORE_TINY_SEQ` | `2.2` | sequential | Tiny-workload sequential bonus. |
| `GP_PAR_SCORE_THREADED_SEQ` | `1.8` | sequential | Thread-sensitive sequential bonus. |
| `GP_PAR_SCORE_OVERHEAD` | `0.45` | parallel | Generic parallelization overhead penalty. |
| `GP_PAR_SCORE_WORKER_STARTUP` | `3.0` | PSOCK-style | PSOCK startup penalty. |
| `GP_PAR_SCORE_FORK_WORKER_STARTUP` | `1.2` | fork-style | Fork startup penalty (cheaper than PSOCK). |
| `GP_PAR_SCORE_WINDOWS_WORKER_STARTUP` | `3.0` | Windows parallel | Windows-specific startup penalty. |
| `GP_PAR_SCORE_BASE_PARALLEL_UNIX` | `0.9` | base_parallel | Unix fork bonus. |
| `GP_PAR_SCORE_FOREACH_PENALTY` | `0.35` | foreach | Generic foreach overhead. |

## Sequential gates

Conditions that force sequential regardless of score.

| Knob | Default | Effect |
|---|---|---|
| `GP_GLOBALS_MAX_GB` | `4` | Total globals (GB) threshold beyond which parallel is not attempted (OOM safety on parallel copies). |
| `GP_PAR_MEMORY_BUDGET_GB` | auto | Explicit GB budget for additional worker processes on Windows, Linux, and macOS. Auto uses currently available physical memory times `GP_PAR_MEMORY_FRACTION`. |
| `GP_PAR_MEMORY_FRACTION` | `0.7` | Fraction of currently available physical memory available for worker processes. |
| `GP_PAR_PAYLOAD_MULTIPLIER` | `4` | Per-worker memory = `GP_PAR_WORKER_OVERHEAD_GB` (or the measured idle worker) + this x the payload each worker receives (its task closure). Calibrated on system-level memory for mice CV (BGLR workers ~3.8x, mixed BGLR/Python ~3.5x including Python children). |
| `GP_PAR_WORKER_OVERHEAD_GB` | `0.75` | Fixed per-worker allowance when the idle worker footprint cannot be measured. |
| `PREDICTPRO_ML_AUTO_THREADS` | `true` | Tree models (RandomForest, XGBoost, CatBoost, LightGBM) with an unset/1-thread setting use this process's share of the cores and therefore run sequentially, not in a worker pool. `false` keeps them single-threaded (and pool-eligible). |
| `PREDICTPRO_SVM_PRECOMPUTE_MAX_N` | `15000` | SVM fits use a BLAS-precomputed kernel up to this many training lines (n x n kernel matrix); libsvm's own kernel above it. |

Workers load PredictProR from the library the session loaded it from and stop
if their version differs from the session's (no silent mixing of versions).
| `gp.rkhs.mem.threshold.GB` (`option`) | `1.0` | RKHS kernel size (GB) beyond which RKHS is forced sequential. |
| `gp.always.seq.RKHS` (`option`) | `TRUE` | When TRUE, every RKHS task is sequential regardless of size; set `FALSE` to let scoring decide. |

## Mori (shared-globals transport)

Mori shares large read-only objects (kernels, gmatrix, etc.) across
workers via shared memory, so a 2 GB kernel isn't copied 16 times for
16 workers. **Shared objects must be treated as immutable** — the
runtime does not enforce this.

| Knob | Default | Effect |
|---|---|---|
| `GP_USE_MORI` | `auto` | Master switch. `auto`/`true` enables when the `mori` package is installed; `false`/`0`/`no`/`off` disables. |
| `GP_MORI_MIN_MB` | `128` | Minimum object size (MB) to even consider sharing. |
| `GP_MORI_SCORE_THRESHOLD` | `1.0` | Score threshold below which sharing is skipped. |
| `GP_MORI_SCORE_SIZE` | `1.0` | Weight on `log1p(size_mb / min_mb)`. |
| `GP_MORI_SCORE_FANOUT` | `0.75` | Weight on `log1p(max(fanout - 1, 0))`. |
| `GP_MORI_SCORE_REUSE` | `0.5` | Weight on `log1p(max(reuse_ratio - 1, 0))`. |
| `GP_MORI_SCORE_COPYVOL` | `0.35` | Weight on `log1p(copy_volume_mb / min_mb)`. |
| `GP_MORI_SCORE_OVERHEAD` | `1.0` | Penalty `1 / log1p(size_mb)` — discourages sharing tiny objects. |

See also `inst/mori-parallel-policy-notes.md` for the longer-form
mori design discussion.

## Mirai daemon

| Knob | Default | Effect |
|---|---|---|
| `GP_MIRAI_QUEUE_MEMORY_MB` | auto | Per-daemon memory budget. Auto: `max(256, min(4096, 128 * fanout))` MB. Override for large kernels. |
| `GP_MIRAI_SUBMIT_TIMEOUT_SEC` | `600` | Submission timeout. |
| `GP_MIRAI_SUBMIT_SLEEP_SEC` | `0.02` | Submission backpressure sleep. |
| `GP_MIRAI_SYNC_STARTUP` | `true` | Wait for daemons to report ready before submitting. `false`/`0`/`no` skips. |

## GPU

| Knob | Default | Effect |
|---|---|---|
| `GP_CUDA_VISIBLE_DEVICES` | `(0:max(0, num_gpus - 1))` | Colon-separated CUDA device IDs visible to workers. |
| `PREDICTPRO_GPU_BUSY_THRESHOLD` | `85` | GPU utilisation % above which the policy treats the GPU as busy and routes new GP tasks to CPU. |

## Reproducibility

| Knob | Default | Effect |
|---|---|---|
| `GP_RANDOM_SEED` | `123` | Base seed used to derive deterministic L'Ecuyer-CMRG streams per task across every backend. |

---

## Decision-reason codes in `Run_metadata.csv`

The full list with descriptions is in
`gp_policy_decision_reason_glossary()`. Codes are joined with `;`
when multiple apply (e.g. `windows_psock;windows_future_startup_penalty`).
Categories: `user_override`, `sequential_driver`, `backend_choice`,
`routing`, `bonus_penalty`, `mori`, `dispatch`.

The two most common questions a `policy_decision_reason` answers:

1. **Why did the policy go sequential here?** Look for any of:
   `tiny_workload`, `single_gpu_model`, `single_gpu_exclusive_deep_models`,
   `multiple_sequential_constraints`, `all_tasks_user_sequential`,
   `no_parallel_eligible_tasks`, `memory_limited_workers`, `rkhs_forced_seq`, `asreml_prepared`,
   the `*_internal_threads` family.

2. **Why did the policy pick *this* backend?** The terminal entry
   (`general_psock_backend` / `reusable_worker_pool` /
   `base_parallel_fallback` / `foreach_fallback`) names the winning
   backend; the surrounding codes name the bonuses and penalties that
   tilted the score toward it.

For programmatic interpretation:

```r
reason <- "windows_psock;windows_future_startup_penalty;general_psock_backend"
codes  <- strsplit(reason, ";", fixed = TRUE)[[1]]
glossary <- gp_policy_decision_reason_glossary()
glossary[glossary$reason %in% codes, c("reason", "category", "meaning")]
```
