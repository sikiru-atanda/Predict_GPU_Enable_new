# Multi-seed deep-learning true prediction

PredictProR exposes DL training seeds separately from bootstrap resampling.
This matters because the two dimensions answer different questions:

- `n_bootstrap` varies the sampled training observations and supports the
  existing true-prediction uncertainty calculation.
- `dl_n_seeds` or `dl_seeds` varies neural-network initialization and training
  randomness while holding each bootstrap sample fixed.

## Preview the work

```r
dl_seed_plan(
  n_bootstrap = 100,
  dl_n_seeds = 10,
  random_seed = 123
)
```

The returned seed manifest is the exact schedule used by the model. The fit
count includes bootstrap and five-fold Gaussian risk-calibration fits. It does
not include optional hyperparameter tuning, a bounded retry after a subprocess
failure, or model fits performed elsewhere in a larger workflow.

## Run derived or explicit seeds

```r
fit <- model_execute(
  pheno_data = phenotype,
  omic1_data = predictors,
  response = "Yield",
  gen_name = "GID",
  GS_model = "DenseNeuralNet",
  n_bootstrap = 100,
  dl_n_seeds = 10,
  random_seed = 123
)

fit_explicit <- model_execute(
  pheno_data = phenotype,
  omic1_data = predictors,
  response = "Yield",
  gen_name = "GID",
  GS_model = "FactorNet",
  n_bootstrap = 100,
  dl_seeds = c(11, 31, 53, 79)
)
```

`random_seed` remains the bootstrap master seed. When `dl_seeds` is omitted,
it is also the first derived training seed. Supplying explicit `dl_seeds`
never changes or replaces those values.

## Aggregation and outputs

For every bootstrap sample, PredictProR fits all requested training seeds on
that same sample and averages their predictions. It then applies the existing
bootstrap uncertainty calculation to those seed-averaged bootstrap
predictions. Unknown true-prediction responses are never used to select a
seed, and there is no best-seed option.

The model result retains:

- `predicted_values`: the standard point-prediction and bootstrap-uncertainty
  table;
- `dl_seed_manifest`: exact training seeds and whether they were explicit or
  derived;
- `dl_computation_plan`: bootstrap, calibration, and total planned fit counts;
- `dl_seed_predictions`: each seed's prediction averaged across bootstrap
  samples;
- `dl_seed_variability`: within-bootstrap seed variance and range, averaged
  across bootstrap samples.

These seed tables are separate from `PEV`, `Standard_error`, and reliability.
Seed variation is optimizer sensitivity; it is not genetic variance and is not
silently relabeled as bootstrap uncertainty.

The multi-seed arguments currently apply to single-trait true prediction with
`cross_validation = FALSE`. PredictProR stops with an explicit scope error if
they are supplied to MET, hybrid, multi-trait DL, or outer cross-validation
routes.
