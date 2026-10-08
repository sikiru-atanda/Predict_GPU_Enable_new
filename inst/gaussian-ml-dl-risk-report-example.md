# Gaussian ML / DL model-comparison report example

This note shows how to read the processed cross-validation outputs for
Gaussian machine-learning and deep-learning models in PredictProR.

Use this view when you want to compare models on:

- prediction accuracy
- predictive uncertainty
- ranking risk

## Typical workflow

```r
cv_out <- PredictProR::model_execute(
  pheno_data = pheno,
  geno_data = geno,
  response = "Yield",
  gen_name = "GID",
  GS_model_cv = c("Ridge_Regression", "CatBoost", "mlp"),
  response_family = "gaussian",
  cross_validation = TRUE,
  cv_evaluation_only = TRUE,
  nfolds = 5L,
  replication = 1L,
  eval_metrics = c("RMSE", "MAE", "cor"),
  system_database = TRUE,
  message = FALSE
)
```

## Key processed outputs

### Accuracy summaries

Use the standard CV metric tables first.

- `cv_out$cv_results_processed$accuracy_table`
- `cv_out$cv_results_processed$model_rankings`

These answer:

- which model predicts best on average?
- which model is best by the chosen evaluation metric?

### Gaussian uncertainty summaries

Use:

- `cv_out$cv_results_processed$gaussian_uncertainty_summaries$aggregated_uncertainty`

Useful columns include:

- `mean_absolute_error`
- `root_mean_squared_error`
- `mean_standard_error`
- `mean_pev`
- `mean_interval_width`
- `empirical_coverage`

These answer:

- how uncertain are the predictions on average?
- are the intervals too narrow or too wide?

### Gaussian ranking-risk summaries

Use:

- `cv_out$cv_results_processed$gaussian_risk_summaries$aggregated_risk`

Useful columns include:

- `mean_prediction_confidence`
- `mean_prediction_unfamiliarity`
- `mean_prediction_risk_score`
- `mean_rank_instability_risk`
- `low_risk_percentage`
- `moderate_risk_percentage`
- `high_risk_percentage`
- `prediction_risk_basis`
- `rank_instability_risk_basis`

These answer:

- which model produces safer generic predictions?
- which model looks more unfamiliar to the training feature distribution?
- which model produces safer ranking decisions?
- how often does each model fall into low / moderate / high ranking risk?

## Plot outputs

Processed Gaussian risk plots are available at:

- `cv_out$cv_results_processed$gaussian_risk_plots$mean_prediction_unfamiliarity`
- `cv_out$cv_results_processed$gaussian_risk_plots$mean_rank_instability_risk`
- `cv_out$cv_results_processed$gaussian_risk_plots$risk_distribution`

Exported PDF files:

- `CV_gaussian_mean_prediction_unfamiliarity.pdf`
- `CV_gaussian_mean_rank_instability_risk.pdf`
- `CV_gaussian_risk_distribution.pdf`

## How to interpret the report

Read the comparison in this order:

1. Accuracy:
   choose the models that predict well.
2. Uncertainty:
   remove models with poorly calibrated or very wide uncertainty.
3. Unfamiliarity:
   check whether a model is making predictions far from the training feature
   distribution.
4. Ranking risk:
   prefer models with lower `mean_rank_instability_risk` and a larger
   `low_risk_percentage`, especially if the downstream use is genotype ranking
   or selection.

## Important interpretation note

For Gaussian ML/DL models:

- `Observed_value` contains every known training response and remains `NA` only
  for genuinely unobserved prediction targets
- `Standard_error` and `PEV = Standard_error^2` are target-specific predictive
  quantities only when target-level resampling variation is available; a
  marginal CV RMSE/MSE is never repeated in every row as individual uncertainty
- when target-specific uncertainty cannot be estimated, `Standard_error` and
  `PEV` are `NA`, while a separately labelled marginal prediction interval may
  still be available from held-out residual calibration
- `Prediction_stability` is the clipped marker-adjustment surrogate
  `1 - SE_i^2 / var(yhat_all)`, where `yhat_all` includes final `Train` and
  `Test` predictions
- `Reliability` is the same non-genetic score retained for plotting compatibility
- `Reliability_variance_input`, `Reliability_reference_variance`, and
  `Reliability_basis` make the calculation and uncertainty source explicit
- `Prediction_confidence` is calibrated predictive trust
- `Prediction_unfamiliarity` is predictor-space novelty
- `Prediction_risk_score` is combined generic prediction risk
- `Rank_instability_risk` is the breeder-facing ranking-risk score

These are not Bayesian / ASReml reliability estimates.
