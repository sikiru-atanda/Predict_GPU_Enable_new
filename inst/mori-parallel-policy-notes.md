# mori parallel memory policy

PredictProR can optionally use the `mori` package to reduce repeated memory
copying of large read-only R objects across local worker processes.

> Companion reference: `inst/parallel-policy-knobs.md` lists every env
> var / option key the parallel policy reads (including the
> `GP_MORI_*` and `GP_PAR_SCORE_*` weights below), and
> `gp_parallel_policy_knobs()` / `gp_policy_decision_reason_glossary()`
> return the same content as data.frames for programmatic inspection.

## Current first-wave scope

The first-wave integration is at the backend/orchestration layer, not inside
individual models.

PredictProR now uses the same orchestration style more broadly for backend
selection itself. In `parallel_mode = "auto"`, the engine does not simply
default to one library. It scores candidate backends and records the decision
reason in runtime metadata.

Eligible backends:

- `mirai`
- `future` with `multisession`
- `base_parallel`
- `foreach`

Currently targeted object types:

- phenotype tables
- genotype matrices
- GRMs and kernels
- processed omics feature tables
- large MET and multi-trait derived feature objects

## How to enable it

Set:

- `GP_USE_MORI=true`

Optional threshold:

- `GP_MORI_MIN_MB=128`

Objects smaller than the threshold are left alone.

The current policy is score-based rather than a single hard cutoff. The
decision considers:

- object size
- backend eligibility
- worker fanout
- task reuse across workers
- estimated copy volume

For the broader backend orchestration, the current `auto` decision also
considers:

- number of parallel-eligible tasks
- worker count
- hidden thread pressure
- total global-object size
- GPU availability for deep models
- Windows vs non-Windows execution
- `mori` suitability as a transport bonus

The size threshold is still used as a minimum gate to avoid obvious
small-object overhead.

## Runtime reason codes

The runtime metadata exported as `Run_metadata.csv` includes a
`policy_decision_reason` field. The most important current values are:

- `single_gpu_exclusive_deep_models`
  - deep-learning tasks were routed to sequential execution because only one
    CUDA GPU was available and outer parallelism would contend for that device
- `rkhs_hidden_threads`
  - `RKHS` was kept sequential because hidden BLAS/thread pressure made outer
    parallel execution unsafe or wasteful
- `xgboost_internal_threads`
  - `Xgboost` was kept sequential because its own thread configuration would
    otherwise collide with outer parallel workers
- `catboost_internal_threads`
  - `CatBoost` was kept sequential because its own thread configuration would
    otherwise collide with outer parallel workers
- `lightgbm_internal_threads`
  - `LightGBM` was kept sequential because its own thread configuration would
    otherwise collide with outer parallel workers
- `multiple_sequential_constraints`
  - more than one overload-safety reason collapsed the queue to sequential
- `no_parallel_eligible_tasks`
  - no task remained worth parallelizing after workload-size and safety checks
- `all_tasks_user_sequential`
  - all tasks were explicitly forced into the user sequential set
- `user_forced_sequential`
  - the user selected `parallel_mode = "sequential"`
- `user_forced_future`
  - the user selected `parallel_mode = "future"`
- `user_forced_mirai`
  - the user selected `parallel_mode = "mirai"`

## What it does not do

- does not affect sequential execution
- does not optimize Python-side memory
- does not optimize GPU tensors
- does not change model results

## Safety rules

Only large immutable objects should be shared. Worker code should treat shared
objects as read-only.

If `mori` is unavailable or an object type is unsupported, PredictProR falls
back automatically to normal copying behavior.

If a backend package is unavailable, the orchestration engine keeps that
candidate in the score table as ineligible and chooses another backend rather
than failing immediately.

## Validation status

The current implementation has been validated on:

- synthetic policy tests
- live `future` multisession synthetic benchmark with a large matrix

In that benchmark, the shared-object path reduced elapsed time while preserving
the returned task values.
