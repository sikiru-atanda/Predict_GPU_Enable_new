# Response Family Support Audit

Date: 2026-04-23

## Current state

PredictProR is still organized primarily around continuous-response regression.

Observed code paths:
- `R/evaluation_metrics_new.R` currently only implements regression-style metrics:
  - `accuracy` as correlation
  - `mean_squared_error`
  - `bias`
  - `root_mean_squared_error`
  - `relative_squared_error`
  - `mean_absolute_error`
  - `mean_absolute_percent_error`
  - `kendalls_tau`
- `R/summary_statistics_AI.R`, `R/summary_statistics_asreml.R`, and `R/summary_statistics_plot.R`
  assume numeric observed and predicted values.
- ML true-prediction diagnostics in `R/diagnostic_plot_true_prediction.R`
  assume continuous predictions.
- Bayesian and ASReml reliability logic is still continuous / latent-value oriented.

## Current internal model signals

### BGLR

Based on current CRAN documentation:
- `gaussian`: supported
- `ordinal`: supported
- binary: supported as the two-class case of ordinal probit
- unordered nominal multiclass: not documented as supported

Sources:
- https://cran.r-project.org/web/packages/BGLR/index.html
- https://cran.r-project.org/web/packages/BGLR/BGLR.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4196607/

### Classical ML

Current public wrappers are coded as regression wrappers:
- `R/AI_extreme_gradient_boosting.R`
- `R/AI_random_forestnew.R`
- `R/AI_Partial_least_square.R`
- `R/AI_support_vector_machinesnew.R`
- `R/AI_K-Nearest_Neighbors.R`
- `R/AI_RidgeRegression_Lasso.R`

They currently scale `y` numerically and return continuous predictions, so public
classification support is not yet cleanly exposed.

### Deep learning

The Python deep-learning stack already contains explicit classification logic:
- `inst/python/dl_models.py`
  - infers `binary` vs `multiclass`
  - uses `BCEWithLogitsLoss` for binary
  - uses multiclass heads / softmax logic for multiclass
- `R/AI_deep_learning_helpers.R`
  - has internal `infer_task()` and `encode_y()` logic
- `R/AI_deep_learning_model_utility.R`
  - already distinguishes `binary`, `multiclass`, and `regression` in output handling

So DL is ahead of the public R API.

## Recommended public response families

Add an explicit argument:

- `response_family = "gaussian" | "binary" | "ordinal" | "multiclass"`

Recommended meaning:
- `gaussian`: numeric continuous trait
- `binary`: unordered two-class outcome
- `ordinal`: ordered categorical trait with 3+ levels
- `multiclass`: unordered nominal trait with 3+ classes

## Recommended support matrix

### Bayesian / BGLR

- `BayesA`: gaussian only
- `BayesB`: gaussian only
- `BayesC`: gaussian only
- `BL`: gaussian only
- `BRR` marker model: gaussian only
- `GBLUP_BRR`: gaussian, binary, ordinal if implemented through documented `BGLR` response types
- `RKHS`: gaussian, binary, ordinal if implemented through documented `BGLR` response types
- block nominal multiclass for all BGLR-backed models

### ASReml

Not yet audited here for non-Gaussian family execution. Treat as continuous-only
until the family-specific fit/predict path is reviewed explicitly.

### Classical ML

Recommended near-term classification candidates:
- `Xgboost`: binary, multiclass
- `RandomForest`: binary, multiclass
- `SupportVectorMachine`: binary, multiclass
- `K-NearestNeighbors`: binary, multiclass

Recommended continuous-only for now:
- `PartialLeastSquare`
- `Ridge_Regression`
- `Lasso`

### Deep learning

Recommended support:
- all current DL families: gaussian, binary, multiclass
- ordinal should be treated as a separate planned head/loss, not silently mapped to nominal multiclass

## Required metric split

### Continuous / gaussian

Keep current regression metrics:
- correlation-based `accuracy` if you keep the current name for backward compatibility
- `mean_squared_error`
- `root_mean_squared_error`
- `mean_absolute_error`
- `bias`
- `relative_squared_error`
- `kendalls_tau`

### Binary

Recommended:
- `roc_auc`
- `pr_auc`
- `log_loss`
- `brier_score`
- `accuracy`
- `balanced_accuracy`
- `f1`
- `precision`
- `recall`
- `specificity`
- `mcc`

### Nominal multiclass

Recommended:
- `accuracy`
- `balanced_accuracy`
- `macro_f1`
- `macro_precision`
- `macro_recall`
- `log_loss`
- optional one-vs-rest multiclass AUC where probabilities are available

### Ordinal

Recommended:
- `accuracy`
- `balanced_accuracy`
- `mean_absolute_error_class`
- `quadratic_weighted_kappa`
- `within_one_class_accuracy`
- ranked probability score when cumulative / class probabilities are available

## Required plot split

### Binary

- confusion matrix
- ROC curve
- precision-recall curve
- calibration plot
- probability histogram by true class

### Multiclass nominal

- confusion matrix heatmap
- class probability distribution
- one-vs-rest ROC / PR where available
- calibration-by-class

### Ordinal

- ordered confusion matrix
- observed vs predicted class distribution
- class-distance error plot
- cumulative probability / threshold diagnostics

## Uncertainty / reliability guidance

Do not reuse continuous `Reliability` semantics for classification.

Recommended naming:
- continuous latent-value models:
  - keep `Reliability` where it is model-based and defensible
- ML / DL classification:
  - use `Prediction_confidence`
  - use `Classification_uncertainty`
  - use calibration summaries separately
- BGLR binary / ordinal:
  - base uncertainty on posterior class probabilities, margins, entropy, or liability uncertainty

Do not label plain class-probability confidence as quantitative-genetic reliability.

## Concrete implementation order

1. Add `response_family` validation at the top-level public API
2. Split `evaluation_metrics()` into family-specific metric sets
3. Split diagnostic plots by response family
4. Add model-family validation so unsupported combinations fail early
5. Expose binary / multiclass DL publicly
6. Add BGLR binary / ordinal support using documented response types
7. Add ordinal-specific uncertainty and plots
8. Add nominal multiclass uncertainty/calibration outputs

## Immediate code hotspots

- `R/model_execute_para.R`
- `R/evaluation_metrics_new.R`
- `R/summary_statistics_AI.R`
- `R/summary_statistics_bayes.R`
- `R/summary_statistics_asreml.R`
- `R/diagnostic_plot_true_prediction.R`
- `R/predict_with_model.R`
- `R/AI_deep_learning_helpers.R`
- `R/AI_deep_learning_model_utility.R`
- `inst/python/dl_models.py`

