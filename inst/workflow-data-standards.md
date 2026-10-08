# Workflow data standards

PredictProR now has package-side input contracts for the main workflow
families:

- general single-trait prediction
- hybrid prediction
- multi-trait prediction
- MET ML/DL

The unified entry points are `prediction_input_standard()` and
`validate_prediction_input()`. The workflow-specific helpers below remain
available and are the validators used by the unified facade.

Exported helpers:

- `PredictProR::general_prediction_data_standard()`
- `PredictProR::validate_general_input_standard()`
- `PredictProR::hybrid_data_standard()`
- `PredictProR::validate_hybrid_input_standard()`
- `PredictProR::multi_trait_data_standard()`
- `PredictProR::validate_multi_trait_input_standard()`
- `PredictProR::met_data_standard()`
- `PredictProR::validate_met_input_standard()`

## General single-trait standard

Required:

- `gen_name`
- one or more response columns

Current expectations:

- single-environment workflows should use one row per genotype
- if genotype rows repeat across environments, provide `heter_groups`

Genomic contracts:

- ML and Bayesian marker models:
  - `geno_data` and/or omics data keyed by `gen_name`
- kernel or relationship models:
  - `gmatrix`, `gkernel`, or omic kernels keyed by `gen_name`
  - or feature data together with relationship/kernel construction inputs

## Multi-trait standard

Required:

- `gen_name`
- at least two response columns
- one row per genotype

Current genomic contracts:

- `multi_trait_ml`
  - `geno_data` and/or omics keyed by `gen_name`
- `multi_trait_dl`
  - `geno_data` and/or omics keyed by `gen_name`
- `multi_trait_asreml`
  - `gmatrix` keyed by `gen_name`

Current scope:

- gaussian only

Cross-validation:

- `multi_trait_ml` — `cv_evaluation_only = TRUE` with K-Folds
- `multi_trait_dl` — `cv_evaluation_only = TRUE` with K-Folds
- `multi_trait_asreml` — `cv_evaluation_only = TRUE` with K-Folds. Each
  fold masks all trait responses for the held-out genotypes, refits the
  joint model, and routes the predict step through the same
  `gp_multitrait_asreml_predict_or_extract()` fallback chain the
  true-prediction path uses, so `asreml::predict()` failures degrade
  to coefficient-based BLUP extraction instead of stopping the run.
  Fold-level fit failures emit a warning, leave NA predictions for that
  fold, and (if all folds collapse) surface through the dropped-tasks
  warning in `gp_filter_crossval_results()`.

ASReml predict-failure fallback coverage (only these scopes are supported;
no general claim is made for other routes):

- single-trait true-prediction + CV, single-env and MET: handled by
  `asreml_predict_or_extract()` → `extract_asreml_prediction()`.
  Verified against `asreml::predict()` (0.20.166, wheat data) for single-env
  and MET `us`, `corgh`, `corh`, `corv`, `fa1`–`fa3`, `rr2`. Fallback rows
  carry point predictions only and are flagged
  `Prediction_uncertainty_source = "unavailable_asreml_predict_failed"`
  (SE, PEV and reliability are NA).
- multi-trait true-prediction + CV (single-env only): handled by
  `gp_multitrait_asreml_predict_or_extract()` →
  `extract_asreml_multitrait_prediction()` (verified for `diag`, `us`).
- multi-trait MET (grouped MT-MET, `fa` structures only): NO coefficient
  fallback. Only workspace errors are remedied (pworkspace retry, then
  per-group `predict(levels = ...)`, which is exact); any other predict
  error stops the run with a message.
- hybrid (asreml + female/male/SCA): does not call `asreml::predict()`;
  reads random effects directly via
  `gp_hybrid_asreml_random_coefs(model)`. Each of the five
  silent-degradation branches in the extractor (summary error, NULL
  coefs, zero-row coefs, missing solution/effect column, no
  matching `vm()` rows) now emits a warning identifying which
  component (Female_GCA / Male_GCA / SCA_effect) collapsed to NA, so
  the intercept-only fallback isn't silent.

## MET ML/DL standard

Required:

- `gen_name`
- one response column
- repeated genotype rows across environments
- an environment column passed through `heter_groups`

Current genomic contracts (any one is sufficient):

- `geno_data` keyed by `gen_name` — a GRM is auto-built; default
  `gmatrix_method = "Yang"`, override (e.g. `"VanRaden"`) to taste
- omic data keyed by `gen_name`
- pre-computed kernels keyed by `gen_name`: any of `gmatrix`, `gkernel`,
  `omic1_kernel`, `omic2_kernel`, `omic3_kernel`
- `kernel_list` — named list of additional N×N PSD kernels keyed by
  `gen_name`. Each kernel (including the GRM and omic kernels above)
  is eigen-decomposed independently and the top components (default
  `var_explained = 0.95`) become a feature block. Blocks are
  concatenated before the ML/DL fit, so multi-kernel inputs are
  supported natively.

Notes on the kernel pipeline:

- The N×N kernel matrix is consumed at decomposition time only; the
  ML/DL model itself sees the per-kernel PCA feature blocks plus a
  one-hot environment block.
- When only `geno_data` is supplied (no kernel, no explicit
  `gmatrix_method`), the input guardrail auto-sets
  `gmatrix_method = "Yang"` and emits a one-line message so the user
  sees what was built.
- CV tasks where every fold returned NA (e.g. mis-keyed kernel rows
  that escape the upfront contract, or a degenerate fold setup) now
  surface as a `warning()` from `gp_filter_crossval_results()`
  identifying the dropped `model/trait/rep` tuples, instead of the
  earlier silent "completed with N valid results".

Current MET CV scenarios:

- `CV0`
- `CV1`
- `CV2`
- `Repeated_CV0`
- `Repeated_CV1`
- `Repeated_CV2`

When uploaded data does not satisfy one of these contracts, `model_execute()`
now stops early and prints the expected standard instead of proceeding into a
later workflow-specific failure.

Recommended user flow:

1. inspect `prediction_input_standard()` before fitting
2. validate the uploaded data with `validate_prediction_input()`
3. only then call `model_execute()`

After fitting, use `prediction_output_standard()` and
`validate_prediction_output()` to verify the common result schema.
