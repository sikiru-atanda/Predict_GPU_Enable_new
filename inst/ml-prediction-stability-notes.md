# Gaussian ML / DL uncertainty and ranking-risk notes

For Gaussian machine-learning and deep-learning outputs in PredictProR,
the user-facing uncertainty and risk columns now have separate meanings.

## Row-level columns

- `Standard_error`
  - square root of held-out cross-fitted predictive MSE when available
- `PEV`
  - held-out cross-fitted predictive MSE; not mixed-model PEV
- `lower_bound`, `upper_bound`
  - interval bounds for the predicted phenotype
- `Prediction_stability`
  - bounded marker-adjustment score `max(0, min(1, 1 - U_i / V_yhat))`
- `Reliability`
  - plotting-compatible alias of `Prediction_stability` for ML/DL
- `Reliability_variance_input`
  - `U_i`, normally prediction-resampling `SE_i^2`
- `Reliability_reference_variance`
  - fitted prediction variance `V_yhat = var(yhat)` across all final
    prediction rows (`Train` and `Test` when both are present)
- `Reliability_basis`
  - exact uncertainty source and explicit non-genetic interpretation
- `Prediction_confidence`
  - calibrated predictive trust score
- `Prediction_risk_score`
  - calibrated generic prediction risk
  - now combines held-out calibrated prediction risk with
    `Prediction_unfamiliarity`
- `Rank_instability_risk`
  - breeder-facing ranking-risk score
  - remains on the held-out ranking-risk path and is kept separate from the
    unfamiliarity-adjusted generic prediction risk
- `Prediction_unfamiliarity`
  - predictor-space novelty / out-of-distribution style signal

## What these columns are not

For Gaussian ML/DL outputs, these columns should not be interpreted as:

- mixed-model or breeding-value genetic reliability
- Bayesian / ASReml `PEV` semantics for latent breeding values
- model-based `Genetic_variance`

Those quantitative-genetic reliability semantics are still reserved for
Bayesian and mixed-model outputs that actually estimate them.

For ML/DL, `Reliability_reference_variance` holds the **fitted prediction
variance across the complete final prediction table**, not phenotypic or
genetic variance. `Train_Test_Label` identifies each row but is not used to
subset this denominator.
`Genetic_variance` remains unavailable. The exact ML/DL calculation is
`Reliability = max(0, min(1, 1 - U_i / V_yhat))`, and both inputs are returned.
Normally `U_i` is prediction-resampling `SE_i^2`; when repeated predictions are
unavailable it may be the legacy uncalibrated in-sample residual-dispersion
proxy used only to keep plots available. Neither input is promoted to
inferential PEV.

## Train vs test risk basis

Risk semantics are different for rows with observed response values and rows
without them.

- `Prediction_risk_basis`
- `Rank_instability_risk_basis`

Typical values:

- `Heldout_Training_OOB`
  - training rows scored from held-out OOB behavior
- `Calibrated_Final_Prediction`
  - final test rows scored with the deployed calibrated risk model

This distinction matters because only training rows have observed responses
available for held-out calibration.

## Interpretation

- `Prediction_stability`
  - answers: how small is the selected predictive-uncertainty input relative
    to observed-response variance?
- `Prediction_confidence`
  - answers: how trustworthy is this prediction after training-only calibration?
- `Prediction_risk_score`
  - answers: how risky is this prediction in a generic calibrated sense,
    including unfamiliarity?
- `Rank_instability_risk`
  - answers: how likely is the prediction to be unstable for ranking decisions?

## For comparison across models

When comparing Gaussian ML/DL models, separate:

- accuracy metrics
  - e.g. `RMSE`, `MAE`, correlation
- uncertainty summaries
  - e.g. `mean_standard_error`, `mean_pev`, interval coverage
- ranking-risk summaries
  - e.g. `mean_rank_instability_risk`
  - percentage `Low/Moderate/High Risk`

PredictProR now exposes processed Gaussian CV comparison outputs for this:

- `cv_results_processed$gaussian_uncertainty_summaries`
- `cv_results_processed$gaussian_risk_summaries`
- `cv_results_processed$gaussian_risk_plots`

Exported Gaussian risk PDFs:

- `CV_gaussian_mean_prediction_unfamiliarity.pdf`
- `CV_gaussian_mean_rank_instability_risk.pdf`
- `CV_gaussian_risk_distribution.pdf`

## Summary report rows

The Gaussian `summary_statistics.csv` report now also exposes
train/test unfamiliarity summaries directly. Useful row names include:

- `Train_Prediction_Unfamiliarity_Basis`
- `Test_Prediction_Unfamiliarity_Basis`
- `Train_Mean_Prediction_Unfamiliarity`
- `Test_Mean_Prediction_Unfamiliarity`
- `Train_Familiar_Percentage`
- `Train_Moderately_Unfamiliar_Percentage`
- `Train_Highly_Unfamiliar_Percentage`
- `Test_Familiar_Percentage`
- `Test_Moderately_Unfamiliar_Percentage`
- `Test_Highly_Unfamiliar_Percentage`

These rows let users see, at the report-summary level:

- whether the unfamiliarity signal is being computed from the predictor-space
  novelty path
- whether the test set is systematically more unfamiliar than the training set
- how much of the train/test output falls into familiar vs unfamiliar bands
