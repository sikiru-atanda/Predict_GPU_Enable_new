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

With `parallel_mode = "auto"` the policy answers one question: run the
tasks sequentially or in parallel? When parallel, it uses `base_parallel`,
the fastest parallel backend measured on Windows and Linux, or `mirai` when
GPU models run on several GPUs (the only backend that gives each worker its
own GPU). Cross-validation and the best-model refit answer the question
with a cost model (expected work against worker start-up); anything else
scores `sequential` against that one backend with the weights below.
`future`, `mirai` and `foreach` remain available through `parallel_mode`.
Tuning knobs adjust the weights and gates that feed the decision.

## How the operating system enters the decision

Nothing has to be set per operating system: `parallel_mode = "auto"`
detects it (`Sys.info()[["sysname"]]`) and uses what was measured fastest
there.

- **How workers start.** `base_parallel` forks its workers on Linux and
  macOS, once per worker (a fork cluster), and starts separate R processes
  (PSOCK) on Windows, which cannot fork. Future uses multicore (fork) or
  multisession the same way.
- **Which backend.** `base_parallel` on every OS, the fastest parallel
  backend measured on both systems (forked on Linux and macOS, PSOCK on
  Windows). The OS changes how its workers start, not which backend runs.
  The score path (cost model off) gives `base_parallel` the bonus
  `GP_PAR_SCORE_BASE_PARALLEL_UNIX` on Linux and macOS, where its workers
  start cheaply, and charges `GP_PAR_SCORE_WINDOWS_WORKER_STARTUP` on Windows.
- **Start-up cost.** The cost model counts more start-up time on Windows
  (`6 + 1.5 x workers` s) than on Linux and macOS (`1 + 0.2 x workers` s).
  Its sequential-or-parallel choices matched the measured faster option on
  all six real cross-validations (three per system), and again on Windows
  after 0.20.202 cut the overhead around the tasks (small CV sequential
  8.6 s vs parallel 10.1 s; medium parallel 11.9 vs 15.7 s; heavy parallel
  30.2 vs 65.3 s, 4 workers).
- **Session state.** Once mirai daemons have run in a session, Linux and
  macOS stop forking for the rest of that session and use PSOCK workers:
  forked workers crashed after mirai.
- **GPUs.** When GPU models run on several GPUs, mirai is used on every OS
  (reason `multi_gpu_mirai`), because it is the backend that gives each
  worker its own GPU; without mirai installed, `base_parallel`. Multi-GPU
  runs have not been timed (the test machine has one GPU).

The decision is recorded in `Run_metadata.csv` (`policy_backend`,
`policy_decision_reason`, e.g. `base_parallel_preferred`, `unix_fork_bonus`).
The timings come from Windows 11 and from Ubuntu under WSL; macOS follows
the Linux rules but was not measured.

### Where a Windows parallel CV spent its time (0.20.201 → 0.20.202)

Heavy CV (600 lines x 2,000 markers, four models, 5 folds, 4 workers),
fresh R session, timestamps from an instrumented build:

| Stage | 0.20.201 | 0.20.202 |
|---|---|---|
| finding a Python interpreter (13 on the test machine), before any task | 12.9 s | 0.0 s |
| worker-memory probe (an R worker; plus a Python child when Lasso counted as Python) | 6.3 s | 2.4 s |
| starting and setting up 4 PSOCK workers (mostly loading packages) | 3.5 s | 3.8 s |
| task function each worker receives and holds | 314 MB | 104 MB |
| the 20 fold tasks on 4 workers | 23.2 s | 23.8 s |
| whole CV stage | 46.4 s | 30.9 s |

Sequential runs paid the interpreter search too (77.7 s → 65.4 s). The
smaller task function does not shorten the set-up; each worker holds the
run's data once instead of three times.

### Worker start-up (0.20.203 → 0.20.204)

CV followed by the refit, 3 traits (500 lines x 1,500 markers; BayesB,
GBLUP_BRR, Lasso; 5 folds; `base_parallel`, 4 workers), the two versions
run alternately in fresh sessions, three rounds:

| Stage | 0.20.203 | 0.20.204 |
|---|---|---|
| worker-memory check | 2.3–5.3 s in every session | 2.4 s in the first session after an install, then 0.03 s |
| CV workers: start and set-up | 2.6–4.2 s | 1.7–1.8 s |
| refit workers: start and set-up | 1.7–2.4 s | 1.7–1.9 s |
| whole run (rounds 2 and 3) | 78.5 / 79.0 s | 72.0 / 71.8 s |

Where a PSOCK worker's set-up time goes, measured on fresh workers:

- **Packages, loaded by every worker at the same time:** dplyr 0.23 s and
  logger 0.04 s. reticulate (0.58 s, 0.43 s of it Matrix) is no longer
  loaded; BGLR tasks never used it.
- **The task function, sent to one worker after another:** 0.13 s per
  worker for the CV's 52 MB, 0.16 s for the refit's 83.5 MB.
- **Each model family's first task on a worker** loads that family's
  packages: BGLR 0.02 s; Lasso and the other R ML models loaded caret
  (1.33 s, 61 namespaces) for centring, scaling and folds. Since 0.20.206
  PredictProR does those itself (`gp_center_scale_fit()`,
  `gp_create_folds()`, caret's results exactly); the first scaler call
  takes 0.00–0.02 s.

The refit starts its own workers: `sp_apply()` starts and stops a cluster
per call. Reusing the CV's workers would save about 0.6–0.9 s on Windows
(the start and the package loading; the task function must still be
sent); not done. On Unix the workers are forked, which costs little.

### Refit task size (0.20.204 → 0.20.205)

The refit's task function was a closure in `model_execute()`'s frame, and
carried all of it: 83.5 MB in the run above, mostly objects the refit
does not use (the input-check contexts, the CV's artifacts, the raw
genotypes). It now carries a copy holding only what it reads
(`gp_trim_closure_frame()`): 15.4 MB. Same run, same alternation:

| Stage | 0.20.204 | 0.20.205 |
|---|---|---|
| refit workers: start and set-up | 1.67–1.69 s | 0.79–0.90 s |
| refit tasks (6, on 3 workers) | 13.6–13.8 s | 12.3–12.6 s |
| whole run (rounds 2 and 3) | 76.4 / 76.7 s | 73.3 / 74.3 s |

Up to 0.20.204 the policy never had the refit's task size:
`gp_parallel_closure_payload_gb()` failed on `model_execute()`'s frame
(an empty `...` cannot be read), and it sized the refit by its data
(`globals`). It now measures the copy. In this run it chose the same
workers (3).

### Without caret (0.20.205 → 0.20.206)

Same run, same alternation:

| Stage | 0.20.205 | 0.20.206 |
|---|---|---|
| worker-memory check | 2.40–2.51 s | 2.06–2.07 s |
| CV tasks (9, on 4 workers) | 50.3–51.0 s | 49.2–49.3 s |
| whole run | 74.5–75.3 s | 73.0–73.3 s |
| measured worker memory (with 25% headroom) | 0.457 GB | 0.499 GB |

The CV's Lasso folds loaded caret on all 4 workers at once (about 1.3 s);
the refit fitted BayesB only. The memory check's reading depends on the
order it loads PredictProR's imports in: on Windows the same 127
namespaces read 376 MB of private bytes with caret loaded third (its place
in 0.20.205's DESCRIPTION) and 417 MB with it last, and 409 MB without
caret. R's live objects are the same in both orders (206 MB); private
bytes follow the heap's high-water mark (268 or 305 MB), which depends on
when collections ran during the loading. The reading repeats to 0.5 MB
within a version but can move by about 40 MB either way when the import
list changes; the 25% headroom covers that. Collecting after every package
load does not steady it (372–469 MB over the same three orders) and makes
the check 2–3 s slower.

### What PSOCK workers receive (0.20.208 → 0.20.209)

Whole `model_execute()` runs without tracing, alternating (CV of 5 folds
and the refit, `base_parallel`, 4 workers):

| | 0.20.208 | 0.20.209 |
|---|---|---|
| 1,000 x 6,000, GBLUP_BRR + Lasso: each CV worker receives | 827 MB | 61.3 MB |
| same run, its CV (4 rounds) | 25.9–27.0 s | 12.4–14.5 s |
| same run, whole | 74.3–83.1 s | 66.7–69.1 s |
| 3 traits, 500 x 1,500 (3 rounds) | 71.4–73.7 s | 63.9–70.7 s |

`sp_apply()`'s seeded-task wrapper was a closure made in its own frame,
where `fail_fast` was still a promise holding the caller's frame
(`gp_execute_partitioned_tasks()`, with the policy's `globals`): every
PSOCK worker received that frame. Now `gp_seeded_task()` builds it and
`sp_apply()` forces its arguments. Measure what workers receive inside a
run, by serialising what `parallel:::sendCall()` sends: the task measured at
`sp_apply()`'s entry, or saved and sent again, did not show this. Each
worker deserialises its message before it starts, and the messages go to
one worker after another.

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

These weights score `sequential` against auto mode's parallel backend
(`base_parallel`, or `mirai` for GPU models on several GPUs) in
`parallel_policy.R`. The score itself is logged at INFO when verbose is
on; reason codes for each contribution land in `policy_decision_reason`.

| Knob | Default | Backend(s) | Effect |
|---|---|---|---|
| `GP_PAR_SCORE_TASKS` | `0.9` | parallel | Weight on `log1p(n_parallel_tasks)`. |
| `GP_PAR_SCORE_WORKERS` | `0.5` | parallel | Weight on `log1p(n_workers)`. |
| `GP_PAR_SCORE_GLOBALS` | `0.35` | parallel | Weight on `log1p(total_globals_GB)`. |
| `GP_PAR_SCORE_GPU_MIRAI` | `3.0` | mirai (multi-GPU); seq (single-GPU) | GPU-aware bonus. |
| `GP_PAR_SCORE_HIDDEN_THREADS` | `0.4` | parallel | Penalty when internal-thread models present. |
| `GP_PAR_SCORE_MORI` | `0.6` | parallel | Bonus when mori sharing applies. |
| `GP_PAR_SCORE_TINY_SEQ` | `2.2` | sequential | Tiny-workload sequential bonus. |
| `GP_PAR_SCORE_THREADED_SEQ` | `1.8` | sequential | Thread-sensitive sequential bonus. |
| `GP_PAR_SCORE_OVERHEAD` | `0.45` | sequential, mirai | Lowers sequential as the queue grows; penalises mirai for 2 or fewer tasks. |
| `GP_PAR_SCORE_WORKER_STARTUP` | `3.0` | base_parallel (Unix PSOCK) | Small-queue startup penalty for new R processes. |
| `GP_PAR_SCORE_FORK_WORKER_STARTUP` | `1.2` | base_parallel (forked) | Fork startup penalty (cheaper than PSOCK), when the workers will fork. |
| `GP_PAR_SCORE_WINDOWS_WORKER_STARTUP` | `3.0` | Windows parallel | Windows-specific startup penalty. |
| `GP_PAR_SCORE_BASE_PARALLEL_UNIX` | `0.9` | base_parallel | Unix bonus (forked, or PSOCK after mirai ran in the session). |

These scores decide only when the cost model is off
(`GP_PAR_COST_MODEL=false`) or has no task costs; cross-validation and the
best-model refit use the cost model (next section).

Retired in 0.20.201, no longer read: `GP_PAR_SCORE_WINDOWS_FUTURE`,
`GP_PAR_SCORE_WINDOWS_FUTURE_STARTUP`, `GP_PAR_SCORE_WINDOWS_MIRAI_PENALTY`,
`GP_PAR_SCORE_MIRAI_WORKER_POOL`, `GP_PAR_SCORE_GPU_FUTURE_PENALTY` and
`GP_PAR_SCORE_FOREACH_PENALTY`. They ranked future, mirai and foreach
against base_parallel, which auto mode no longer does.

### Measured backend timings (0.20.199)

`tools/benchmark_parallel_backends.R` times `sp_apply()` per backend on CPU
tasks of a fixed size (4, 16 or 64 tasks of 0.05 s or about 1 s, with or
without 0.25 GB the tasks close over, 4 workers);
`tools/analyze_parallel_backends.R` summarises the runs and
`tools/tune_parallel_score_weights.R` replays auto mode's score path against
them.
On an i7-12700H laptop (Windows 11, R 4.3.3; Ubuntu under WSL, R 4.4.1):

| | Windows | Linux |
|---|---|---|
| Fastest parallel backend | base_parallel (PSOCK) in every shape | base_parallel (forked) in 11 of 12 shapes |
| 64 short tasks (0.05 s) | base_parallel 2.1 s, foreach 2.3, future 2.6, sequential 2.9, mirai 3.0 | fork 1.0 s, future multicore 1.7, PSOCK 2.1, mirai 2.5, sequential 3.0 |
| 64 longer tasks (0.7 s Windows, 0.9 s Linux) | base_parallel 16.5 s, foreach 18.0, future 20.1, mirai 22.3, sequential 38.8 | mirai 21.0 s, fork 21.4, PSOCK 21.5, foreach 26.8, sequential 57.2 |
| Start-up of 2 workers | base_parallel 1.2 s, foreach 1.0, future 1.3, mirai 1.8 | fork 0.1 s, future multicore 0.1, PSOCK 1.0, mirai 1.2 |

With the weights before 0.20.199 the score path picked mirai on Windows for
every larger queue; over the measured shapes its picks lost 33 s against the
fastest backend. Since 0.20.201 auto mode does not choose among parallel
backends at all, and base_parallel alone would have lost 0 s on Windows and
0.5 s on Linux had each shape been run sequentially or in parallel,
whichever was faster. The score path still loses 14.8 s, all of it on
sequential versus parallel: it does not know how long a task takes, so no
weights separate 4 x 0.05 s tasks from 4 x 1 s tasks. The cost model knows.

## Sequential gates

Conditions that force sequential regardless of score.

| Knob | Default | Effect |
|---|---|---|
| `GP_GLOBALS_MAX_GB` | `4` | Total globals (GB) threshold beyond which parallel is not attempted (OOM safety on parallel copies). |
| `GP_PAR_MEMORY_BUDGET_GB` | auto | Explicit GB budget for additional worker processes on Windows, Linux, and macOS. Auto uses currently available physical memory times `GP_PAR_MEMORY_FRACTION`. |
| `GP_PAR_MEMORY_FRACTION` | `0.7` | Fraction of currently available physical memory available for worker processes. |
| `GP_PAR_PAYLOAD_MULTIPLIER` | `4` | Per-worker memory = `GP_PAR_WORKER_OVERHEAD_GB` (or the measured idle worker) + this x the payload each worker receives (its task closure). Calibrated on system-level memory for mice CV (BGLR workers ~3.8x, mixed BGLR/Python ~3.5x including Python children). |
| `GP_PAR_WORKER_OVERHEAD_GB` | `0.75` | Fixed per-worker allowance when the idle worker footprint cannot be measured. Otherwise the policy starts one worker with PredictProR's imports (plus a Python child with numpy/sklearn when Python models run), adds 25%, and keeps the result in `python-discovery.rds` (`tools::R_user_dir("PredictProR", "cache")`) until PredictProR or an import is reinstalled, R or the ML interpreter changes, for at most 30 days; `PREDICTPRO_PROBE_CACHE=false` measures in every session. Setting this variable skips the measurement. |
| `PREDICTPRO_ML_AUTO_THREADS` | `true` | Tree models (RandomForest, XGBoost, CatBoost, LightGBM) with an unset/1-thread setting use this process's share of the cores and therefore run sequentially, not in a worker pool. `false` keeps them single-threaded (and pool-eligible). |
| `PREDICTPRO_SVM_PRECOMPUTE_MAX_N` | `15000` | SVM fits use a BLAS-precomputed kernel up to this many training lines (n x n kernel matrix); libsvm's own kernel above it. |
| `gp.rkhs.mem.threshold.GB` (`option`) | `1.0` | RKHS kernel size (GB) beyond which RKHS is forced sequential. |
| `gp.always.seq.RKHS` (`option`) | `TRUE` | When TRUE, every RKHS task is sequential regardless of size; set `FALSE` to let scoring decide. |

Workers load PredictProR from the library the session loaded it from and stop
if their version differs from the session's (no silent mixing of versions).

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
| `GP_MIRAI_QUEUE_MEMORY_MB` | auto | Cap on the task payloads queued at the dispatcher in the main session; a full queue makes the next submission wait. Auto: `max(256, min(4096, 128 * fanout))` MB. |

`sp_apply(backend = "mirai")` starts real daemon processes. Each one gets
the session's runtime environment, its GPU and the task function with the
data it closes over once (`mirai::everywhere()`); a task then sends only
its own element. A daemon that fails its set-up (for example a different
PredictProR in its library) stops the run. Before 0.20.199 the daemons were
started with `sync = TRUE`, which in mirai 2.x runs every task in the main
session, so `parallel_mode = "mirai"` ran sequentially.
`GP_MIRAI_SYNC_STARTUP`, `GP_MIRAI_SUBMIT_TIMEOUT_SEC` and
`GP_MIRAI_SUBMIT_SLEEP_SEC` were removed; they had no effect.

## GPU

| Knob | Default | Effect |
|---|---|---|
| `GP_CUDA_VISIBLE_DEVICES` | every GPU (`0:1:...`) | Colon-separated CUDA device IDs for the GPU workers (e.g. `2:3`). The mirai backend starts one daemon per id, each seeing only its own GPU. |
| `PREDICTPRO_GPU_BUSY_THRESHOLD` | `85` | GPU utilisation % above which the policy treats the GPU as busy and routes new GP tasks to CPU. |

## Reproducibility

| Knob | Default | Effect |
|---|---|---|
| `GP_RANDOM_SEED` | `123` | Base seed used to derive deterministic L'Ecuyer-CMRG streams per task across every backend. |

---

## Decision-reason codes in `Run_metadata.csv`

The full list with descriptions is in
`gp_policy_decision_reason_glossary()`. Codes are joined with `;`
when multiple apply (e.g. `cost_model_parallel(est_seq=600s,est_par=301s);base_parallel_preferred`).
Categories: `user_override`, `sequential_driver`, `backend_choice`,
`routing`, `bonus_penalty`, `mori`, `dispatch`.

A run that cross-validates and then refits the best models
(`cv_evaluation_only = FALSE`) makes two decisions, and `Run_metadata.csv`
holds both:

- `policy_*` describe the refit; `policy_stage = best_model_refit` says so.
  With one trait it is often a single sequential task
  (`no_parallel_eligible_tasks`).
- `cv_policy_*` describe the cross-validation: `cv_policy_backend`,
  `cv_policy_workers`, `cv_policy_decision_reason`,
  `cv_policy_parallel_count`, and so on. The CV's other fields follow,
  prefixed the same way: `cv_task_unit`, `cv_metric_for_ranking` and the
  `cv_profile_*` timings.
- The CV stage's own table is also written as
  `cv_results_processed_run_metadata.csv`.

A CV-only run (`cv_evaluation_only = TRUE`) writes the CV table as
`Run_metadata.csv`, with its decision in `policy_*`.

The two most common questions a `policy_decision_reason` answers:

1. **Why did the policy go sequential here?** Look for any of:
   `tiny_workload`, `single_gpu_model`, `single_gpu_exclusive_deep_models`,
   `multiple_sequential_constraints`, `all_tasks_user_sequential`,
   `no_parallel_eligible_tasks`, `memory_limited_workers`, `rkhs_forced_seq`, `asreml_prepared`,
   the `*_internal_threads` family.

2. **Why did the policy pick *this* backend?** In auto mode it is
   `base_parallel_preferred` or, for GPU models on several GPUs,
   `multi_gpu_mirai`; `cost_model_parallel(...)` gives the cost model's
   estimates, and the other codes name the bonuses and penalties of the
   score path.

For programmatic interpretation:

```r
reason <- "base_parallel_preferred;unix_fork_bonus"
codes  <- strsplit(reason, ";", fixed = TRUE)[[1]]
glossary <- gp_policy_decision_reason_glossary()
glossary[glossary$reason %in% codes, c("reason", "category", "meaning")]
```
