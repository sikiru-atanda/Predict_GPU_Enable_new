# Gaussian Process Mixed-Model Framework

This project exposes three public model routes:

- `gp_exact`: exact GP posterior prediction with optional GP posterior SEs.
- `krr_exact`: kernel-ridge / GP-equivalent posterior prediction for fast dense problems.
- `gp_icm_fa`: exact GP posterior prediction with factor-analytic environment covariance.

## Output Levels

Use `output_level` to control work:

- `predict_only`: return predictions only.
- `predict_with_se`: return predictions plus GP posterior prediction SEs.
- `full_vc`: return predictions plus GP posterior SEs, then run AI-REML only for variance-component reporting.

Prediction means and prediction SEs must come from the selected GP/KRR posterior path. AI-REML must not mutate the prediction model by default; it is reserved for variance-component estimates, standard errors, and ASReml-style reporting.

## Scalable Backend

For large `n`, pass `grm_factor_cache` and use `backend="auto"` or `backend="operator"`. The operator backend uses low-rank GRM roots and supports:

- point predictions;
- diagonal prediction SEs via Hutchinson probing when requested;
- no full prediction covariance.

The canonical import path is `mixed_model_gpu_large_scale_backend.py`; the dated backend implementation remains as the implementation file.

## Multi-Trait MET Status

`fit_multi_trait_gp_met` is a MET-specific GP/operator path, not a full `gp_icm_fa` replacement. It supports scalable multi-trait, multi-environment posterior mean prediction with `varcomp_mode="mom"`, optional exact posterior SE columns via `return_se=True`, optional FA projection for trait covariance matrices, and multi-root GRM caches through `kernel_weights`.

For small/medium problems, `varcomp_mode="reml"` runs a bounded dense REML variance-estimation pass and then returns to the GP/operator prediction path. It refuses large training sets through `dense_reml_max_train` instead of silently doing cubic work.

User-supplied `env_similarity` should be passed as a labeled `pandas.DataFrame` or with `env_ids` so the environment kernel is aligned correctly.

The MET path requires one pre-aggregated row per `(gid, env, trait)` cell. Prediction rows may have missing `y`; training rows must have finite `y`.

## Core Dependencies

Runtime: `numpy`, `pandas`, `torch`, `gpytorch`, `linear_operator`, `zarr`.

Tests/benchmarks additionally use `pytest` and optional data readers such as `pyreadr`.
