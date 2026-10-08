# MET ML / DL Example Workflows

These examples show the current public `model_execute()` workflow for
multi-environment machine-learning and deep-learning prediction in PredictProR.

The current MET ML/DL path:

- uses long-format `GID x Env` rows
- uses kernel eigenfeatures as the predictor representation
- supports `CV0`, `CV1`, and `CV2`
- returns per-environment predictions and, for Gaussian true prediction,
  across-environment predicted values

For Gaussian MET true prediction, the output includes empirical
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

These ML/DL uncertainty fields are calibrated from model residual behavior.
They are not ASReml/BGLR variance-component estimates.

## Example 1: Gaussian MET with CatBoost

Use this pattern when the trait is continuous and the goal is per-environment
prediction.

```r
out <- PredictProR::model_execute(
  pheno_data = pheno_met,
  gmatrix = grm,
  omic1_kernel = cop,
  omics_kernel_label = list(omic1_kernel = "COP"),
  response = "Yield",
  gen_name = "GID",
  heter_groups = "Env",
  GS_model = "CatBoost",
  response_family = "gaussian",
  met_ml_dl = TRUE,
  system_database = FALSE,
  message = FALSE,
  catboost_iterations = 300L,
  catboost_depth = 6L,
  catboost_learning_rate = 0.05,
  catboost_l2_leaf_reg = 3L,
  catboost_thread_count = 1L
)
```

Expected outputs:

- `predicted_values.csv`
  - one row per `GID x Env`
  - columns include `Predicted_value`, `Train_Test_Label`, `Standard_error`,
    `PEV`, and `Reliability`
- in-memory `model_results$across_environment_predicted_values`
  - also exposed through the legacy `model_results$Total_Predicted_value`
  - one row per `GID`
  - columns include `Predicted_value`, `Train_Test_Label`, `Standard_error`,
    `PEV`, and `Reliability`
- `MET_pooled_metrics.csv`
- `MET_metrics_by_environment.csv`
- `MET_prediction_counts_by_environment.csv`
- `MET_prediction_table.csv`
- `MET_test_diagonistic_plots_GS_model.pdf`

For Gaussian ML/DL runs, the current predictive-uncertainty columns are:

- `Standard_error`
- `PEV`
- `lower_bound`
- `upper_bound`
- `Uncertainty`
- `Uncertainty_remarks`
- `Reliability`
- `Reliability_remarks`

Interpretation:

- `Standard_error`: empirical prediction standard error calibrated from
  observed residual behavior
- `PEV`: prediction error variance, equal to `Standard_error^2`
- `lower_bound` / `upper_bound`: approximate prediction interval bounds
- `Reliability`: bounded 0-1 reliability score derived from prediction
  error variance relative to observed trait variance
- `Uncertainty` and `Uncertainty_remarks`: breeder-facing uncertainty
  category and label

For Gaussian ML/DL outputs:

- `Train` rows have observed phenotypes available during the fitted run
- `Test` rows are masked or user-provided prediction targets

## Example 2: Binary MET with mlp

Use this pattern for a binary trait such as resistant / susceptible,
survived / failed, or yes / no.

```r
out <- PredictProR::model_execute(
  pheno_data = pheno_met,
  gmatrix = grm,
  omic1_kernel = cop,
  omics_kernel_label = list(omic1_kernel = "COP"),
  response = "Trait",
  gen_name = "GID",
  heter_groups = "Env",
  GS_model = "mlp",
  response_family = "binary",
  met_ml_dl = TRUE,
  system_database = FALSE,
  message = FALSE,
  compile_model = FALSE,
  deterministic = TRUE,
  random_seed = 1L,
  validation_split = 0.2,
  device = "cpu",
  use_amp = FALSE,
  batch_norm = FALSE,
  epochs = 10L,
  batch_size = 32L,
  mlp_neurons_per_layer = as.integer(c(128L, 64L)),
  mlp_learning_rate = 1e-3,
  dropout = 0.1
)
```

Expected row-level prediction columns:

- `Predicted_value`
- `Predicted_class`
- `Prediction_confidence`
- `Classification_uncertainty`
- `Prob_*`

Expected summary files:

- `MET_pooled_metrics.csv`
  - `accuracy`
  - `balanced_accuracy`
  - `log_loss`
  - `brier_score`
  - `ece`
- `MET_metrics_by_environment.csv`
- `MET_test_diagonistic_plots_GS_model.pdf`

## Example 3: Multiclass MET with CV2

This example shows the public CV path for a multiclass trait.

```r
cv_out <- PredictProR::model_execute(
  pheno_data = pheno_met,
  gmatrix = grm,
  omic1_kernel = cop,
  omics_kernel_label = list(omic1_kernel = "COP"),
  response = "Trait3",
  gen_name = "GID",
  heter_groups = "Env",
  GS_model_cv = "ft_transformer",
  response_family = "multiclass",
  met_ml_dl = TRUE,
  cross_validation = TRUE,
  cv_evaluation_only = TRUE,
  cross_validation_meth = "CV2",
  nfolds = 5L,
  replication = 1L,
  eval_metrics = c("accuracy", "log_loss", "macro_f1", "ece"),
  system_database = TRUE,
  message = FALSE,
  compile_model = FALSE,
  deterministic = TRUE,
  random_seed = 1L,
  validation_split = 0.2,
  device = "cpu",
  use_amp = FALSE,
  batch_norm = FALSE,
  epochs = 10L,
  batch_size = 32L,
  ft_d_model = 64L,
  ft_heads = 4L,
  ft_layers = 2L,
  ft_ff_mult = 2L,
  ft_dropout = 0.1,
  ft_token_dropout = 0.1,
  ft_use_cls = TRUE
)
```

Expected CV outputs:

- `cv_results_raw[[1]]$ypred_cv_Reps_all`
- `cv_results_raw[[1]]$yprob_cv_Reps_all`
- `cv_results_processed$classification_probability_summaries`

Recommended multiclass metrics:

- `accuracy`
- `balanced_accuracy`
- `macro_precision`
- `macro_recall`
- `macro_f1`
- `log_loss`
- `ece`

## Supported MET models

Current verified MET ML/DL models:

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

Current verified response families:

- `gaussian`
- `binary`
- `multiclass`

See also:

- `inst/gaussian-ml-dl-risk-report-example.md`
  - example of comparing Gaussian ML/DL models on accuracy, uncertainty, and ranking risk
