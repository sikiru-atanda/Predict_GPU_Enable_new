# PredictProR Refactor Baseline

Date: 2026-04-21

## Scope of this baseline

This document captures the package state before the architectural refactor starts. It is meant to give a stable reference point for runtime consolidation, dependency cleanup, and later performance work.

## Current package shape

- R source files: 118
- Approximate R code size: 1.96 MB
- Function definitions detected: 473
- Main complexity concentration:
  - `R/model_execute_para.R`
  - `R/models_execute_crossval.R`
  - `R/models_execute_crossval_mirai_purri.R`
  - `R/AI_deep_learning.R`
  - `R/ETA_compiler_bayes_GBLUP.R`

## Baseline findings

### 1. Runtime duplication

The package currently defines overlapping execution and scheduling logic in multiple places.

- Cross-validation entry points exist in:
  - `R/main_crossvalidation_execution_logic.R`
  - `R/models_execute_crossval.R`
  - `R/models_execute_crossval_mirai_purri.R`
- Parallel policy and task application logic are also duplicated across:
  - `R/parallel_policy.R`
  - `R/sp_apply_utils.R`
  - `R/models_execute_crossval.R`
  - `R/models_execute_crossval_mirai_purri.R`

This creates drift risk and makes performance tuning unreliable because runtime behavior depends on which copy is actually used.

### 2. Package boundary problems

- `NAMESPACE` exported almost everything through `exportPattern("^[[:alpha:]]+")`.
- Public and internal helpers are not clearly separated.
- `DESCRIPTION` was missing several runtime dependencies used by the package boundary and scheduler stack.

### 3. Load-time overhead and side effects

- Several package source files used top-level `library(...)` calls.
- `R/zzz.R` defined `.onAttach()` twice, causing silent override of startup behavior.
- Startup messaging existed alongside commented installation/setup logic, increasing noise around package load responsibilities.

### 4. Scheduler/resource issues

- The package tries to support multiple orchestration layers at once:
  - base `parallel`
  - `foreach`
  - `future`
  - `mirai`
  - library-internal threading
  - Python/PyTorch GPU execution
- CPU, BLAS/OpenMP, XGBoost, and GPU resource handling is only partially centralized.
- Large prepared objects are still captured by task closures, which is especially expensive for Windows worker models.

## Verification checklist for refactor safety

These checks should continue to pass after each refactor phase.

1. Package load succeeds without triggering installation side effects.
2. Startup messaging is deterministic and emitted once.
3. `model_execute()` remains callable with existing arguments.
4. `models_execute_crossval()` still returns results for a minimal hold-out workflow.
5. `predict_with_model()` continues to dispatch for Bayesian, asreml, and AI model families.
6. Runtime helpers can choose sequential mode when models already parallelize internally.
7. Worker initialization does not require repeated package attach side effects.
8. Python-backed code paths still initialize through `reticulate` only when needed.

## Step 1 deliverables completed in this pass

- Baseline documented.
- Public API export surface made explicit in `NAMESPACE`.
- Core runtime dependencies declared more clearly in `DESCRIPTION`.
- Startup hooks consolidated in `R/zzz.R`.
- Top-level `library(...)` calls removed from the core runtime helper files.

## Deferred to later phases

- Runtime consolidation into a single canonical cross-validation engine.
- Model registry introduction.
- Serialization-light task payload design.
- Python bridge isolation.
- Scheduler simplification to one canonical backend policy.
