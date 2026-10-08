# PredictProR architecture and statistical audit

Audit date: 2026-08-12

## Contract boundary

`model_execute()` remains the single fitting entry point. Its model results use
four required fields in every family:

- `model_parameters`
- `predicted_values`
- `diagnostic_plots`
- `variance_components`

Environment totals, covariance/correlation tables, and covariance matrices are
optional named additions. Legacy input names are accepted internally, but the
public model result uses only canonical field names. Use
`prediction_output_standard()` to inspect the schema and
`validate_prediction_output()` to audit a returned object.

The unified input facade is `prediction_input_standard()` plus
`validate_prediction_input()`. It dispatches to the existing general, hybrid,
multi-trait, and multi-environment validators, so the same validation rules are
used by documentation, preflight checks, and model execution.

## Statistical model inventory

| Family | Modes in scope | Covariance/kernel representation | Validation boundary |
|---|---|---|---|
| ASReml GBLUP | single trait, MET, multi-trait, hybrid | inverse relationship matrices in `vm()` terms; CS, US, CORGH/CORH/CORV, FA and RR environment structures | inverse/name alignment tests; checked covariance reconstruction; per-environment variance and heritability tests |
| Bayesian marker | BRR, BayesA/B/C, Bayesian Lasso | marker design matrix with BGLR priors | design/response routing, MCMC defaults, shared-fit and classification smoke tests |
| Bayesian kernel | GBLUP_BRR and RKHS; single/MET/multi-trait | one or more named PSD relationship kernels | kernel alignment, per-environment variance, heteroscedastic residual, CV equivalence tests |
| Gaussian process | KRR, GP, GP_FA, LowRankGP | named kernel bank; exact or low-rank factorization; optional environment covariance | boundary, environment-label, sparse-mask, warm-cache and generalization tests |
| Classical ML | RF, ridge, lasso, PLS, SVM, KNN, XGBoost, CatBoost, LightGBM | marker/omic features; MET kernels become independent PCA feature blocks | no-leakage preprocessing, classification, missingness and runtime smoke tests |
| Deep learning | CNN, transformer/attention, TabNet, SAINT, NODE, DeepFM, DCNv2, NAM, MoE, MLP, ResNet, DKL | marker/omic features; MET kernels become independent PCA feature blocks | prediction-scale, deterministic seed, device routing and runtime smoke tests |
| Hybrid | ASReml, Bayesian, GP, ML, DL | female GCA, male GCA and optional SCA components | parent/cross ID contract, component output, prediction consistency and runtime smoke tests |

## Covariance findings and decisions

1. Every kernel in a multi-kernel ASReml MET model now receives the requested
   environment covariance structure. The former behavior applied FA/US/COR to
   only the first kernel and silently reduced later kernels to `idv`, which did
   not represent the requested multi-kernel model.
2. FA covariance is reconstructed as `Lambda %*% t(Lambda) + Psi` through one
   checked helper. Loading dimensions, finite values, and non-negative specific
   variances are enforced.
3. US and correlation structures validate the fitted parameter count before
   reconstruction. Correlation matrices must have a unit diagonal, be
   symmetric, and be PSD. This also fixes the previous CORH diagonal error that
   multiplied each reported genetic variance by the common correlation.
4. Multiple independent genetic/omic random terms contribute additively to
   total genetic variance. Per-environment heritability uses the sum of their
   covariance diagonals divided by that sum plus the corresponding residual
   variance.
5. Multi-trait ASReml data are stacked in observational-unit-major,
   trait-minor order, matching `units:<structure>(Trait)`. The public route
   accepts `us`, `corgh`, or `diag` trait covariance instead of silently
   forcing an unstructured model.

## Licensed ASReml reference validation

ASReml-R 4.2.0.392 was exercised against a locally configured network license
on Windows with R 4.3 on 2026-08-12. The reference run covered direct
`vpredict()` agreement for delta-method heritability standard errors, CORGH
per-environment reconstruction, multi-trait diagonal covariance, hybrid GCA/SCA
and prediction routes, and a converged two-kernel FA1 MET. For the FA1 model,
both kernels had explicit FA terms and both reconstructed covariance matrices
were finite, symmetric, and positive semidefinite.

Use `tools/validate_asreml_reference.R` for the quick licensed check and add
`--full` for the multi-trait and hybrid routes. License configuration remains
external to the package; see `inst/asreml-local-validation.md`.

## Kernel findings and decisions

- VanRaden, Yang, epistatic, Vitezica dominance, and Su heterozygosity formulas
  are pinned against independent R reference calculations and native C++ parity
  tests.
- Weighted VanRaden is explicitly defined as
  `Z W t(Z) / sum(2 p (1-p))`. Weights must be finite, non-negative, and contain
  at least one positive value so the result remains PSD under this definition.
- Gaussian `theta` must be positive. Composite-kernel `alpha` must be in
  `[0, 1]`, preserving a convex combination of PSD kernels.
- `relationship_matrix_diagnostics()` checks dimensions, finite values,
  symmetry, sample-name identity, eigenvalue tolerance, rank, and conditioning.
  `validate_relationship_matrix()` turns any violation into an early error.

## Parallel and memory policy

The canonical CV scheduler is `models_execute_crossval()` with
`gp_execute_partitioned_tasks()` and `sp_apply()`. The files named `legacy_*`
remain compatibility implementations and are not the public scheduler.

Worker selection now inventories the complete CV `params` payload, including
large kernel/marker/omic objects and caches. The memory cap uses:

`per_worker = serialized_payload * GP_PAR_PAYLOAD_MULTIPLIER + GP_PAR_WORKER_OVERHEAD_GB`

Defaults are 1.25 and 0.25 GiB. The worker budget is explicitly configured with
`GP_PAR_MEMORY_BUDGET_GB`, or derived on Windows, Linux, and macOS from available
physical memory times `GP_PAR_MEMORY_FRACTION` (default 0.7). BLAS/OpenMP threads
are pinned per worker to avoid nested oversubscription. GPU-exclusive and
internally threaded models remain sequential unless the policy can route them
without contention.

Use `parallel_memory_plan()` before a large run. The local development probe
`tools/benchmark_parallel_memory.R` records both policy estimates and observed
PSOCK worker RSS without uploading data.

Local Windows probe (2,500 x 2,500 identity kernel, 2026-08-11): the serialized
kernel was 47.99 MiB and the full listed payload was 48.33 MiB. Two PSOCK
workers used 107.82 MiB mean RSS each (215.64 MiB total). The policy estimated
316.42 MiB per worker after its 1.25 payload multiplier and 0.25 GiB fixed
allowance, confirming that the default cap is conservative for this workload.
The observed value is a calibration point, not a portable constant; model fit
allocations and loaded backends can increase it.

## Claim boundaries

- Algebraic and mocked variance-component tests validate all supported
  reconstruction branches without a commercial ASReml installation. The
  licensed reference run above covers selected high-risk branches, not
  numerical parity for every covariance structure, dataset, or ASReml version.
- Runtime smoke tests establish routing and output contracts; they are not
  evidence that every ML/DL architecture is optimally calibrated for a given
  breeding population.
- Dense relationship matrices remain O(n^2) storage. The worker policy prevents
  unsafe fan-out but does not make dense kernels sparse; use low-rank GP or a
  deliberately sparse representation when population size makes dense storage
  infeasible.
