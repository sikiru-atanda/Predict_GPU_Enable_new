# MET ML/DL Kernel Design

This document records the intended first implementation path for multi-environment
machine learning and deep learning prediction in PredictProR.

For user-facing public examples, see:

- `inst/met-ml-dl-examples.md`

## Scope

The current version targets:

- per-environment prediction
- across-environment predicted values for Gaussian true prediction
- `CV0`, `CV1`, and `CV2`
- multiple kernels through separate eigendecomposition plus feature concatenation

Current verified model scope:

- Tabular ML
  - `CatBoost`
  - `LightGBM`
  - `Xgboost`
  - `RandomForest`
- Deep learning
  - `mlp`
  - `ft_transformer`
  - `saint`
  - `tabnet`
  - `moe`

Current verified response-family scope:

- `gaussian`
  - all models listed above
- `binary`
  - all models listed above
- `multiclass`
  - all models listed above

For Gaussian true prediction, this path returns empirical
residual-calibrated uncertainty fields for both per-environment and
across-environment predictions:

- `Standard_error`
- `PEV`
- `lower_bound`
- `upper_bound`
- `Uncertainty`
- `Uncertainty_remarks`
- `Reliability`
- `Reliability_remarks`

These ML/DL uncertainty fields are not ASReml/BGLR variance-component
estimates. The across-environment table is an aggregate prediction contract
for breeder-facing ranking and plotting, not a mixed-model total breeding
value decomposition.

## Representation

1. Convert genomic / omic sources to `GRM` or kernels
2. Eigendecompose each kernel
3. Retain principal components up to a variance threshold
4. `cbind()` retained component blocks
5. Expand to long format with one row per `GID x Env`

## Why this path

- supports sparse `CV2` because training uses only observed rows
- borrows information across environments by fitting all observed `GID x Env`
  rows jointly
- keeps source identity for multiple kernels
- reduces dimensionality for ML/DL training

## Current helper functions

- `gp_kernel_to_eigenfeatures()`
- `gp_prepare_met_kernel_features()`
- `gp_build_met_long_data()`
- `gp_met_cv_assignments()`
- `gp_met_summary_statistics()`
- `gp_met_true_prediction_plot()`

## Current user-facing outputs

- `predicted_values.csv`
  - one row per `GID x Env`
  - gaussian: includes `Predicted_value`, `Train_Test_Label`,
    `Standard_error`, `PEV`, and `Reliability`
  - binary/multiclass: includes `Predicted_value`, `Predicted_class`,
    `Prediction_confidence`, `Classification_uncertainty`, and `Prob_*`
- in-memory `model_results$across_environment_predicted_values`
  - gaussian true prediction only
  - one row per `GID`
  - also exposed through the legacy `model_results$Total_Predicted_value`
- `MET_pooled_metrics.csv`
  - gaussian: pooled row-level regression metrics across observed rows
  - binary: pooled classification metrics such as `accuracy`, `balanced_accuracy`,
    `log_loss`, `brier_score`, `ece`
  - multiclass: pooled classification metrics such as `accuracy`,
    `balanced_accuracy`, `macro_precision`, `macro_recall`, `macro_f1`,
    `log_loss`, `ece`
- `MET_metrics_by_environment.csv`
  - per-environment regression or classification metrics
- `MET_prediction_counts_by_environment.csv`
  - train/test/observed row counts by environment
- `MET_prediction_table.csv`
  - joined prediction/observed table used for MET summaries
- `MET_test_diagonistic_plots_GS_model.pdf`
  - gaussian: observed-vs-predicted plus per-environment prediction distribution
  - binary: calibration-style plus confidence distribution panels
  - multiclass: confusion-matrix plus confidence panels
