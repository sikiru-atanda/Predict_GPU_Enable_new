# Classification output contract

PredictProR uses one standardized public prediction table for binary, ordinal,
nominal, and multiclass responses. `nominal`, `multinomial`, `multi-class`, and
`multi_class` are aliases for the canonical unordered `multiclass` family.
Ordered responses use the canonical `ordinal` family and retain the declared
factor-level order.

For observation `i` and class `k`, let `p_hat[i,k]` be the fitted predictive
class probability. Ordinal Bayesian models use BGLR posterior class
probabilities. Nominal/multiclass ML and DL models use the predictive class
probabilities returned by their fitted classifier.

The public quantities are:

- `Predicted_class = argmax_k p_hat[i,k]`;
- `Prediction_confidence = max_k p_hat[i,k]`;
- `Classification_uncertainty = 1 - Prediction_confidence`;
- `Reliability = Prediction_confidence`.

Thus classification reliability is predictive class confidence. It is not
Gaussian PEV reliability, genetic reliability, or heritability. Row-level
classification prediction tables deliberately omit `Standard_error`, `PEV`,
genetic/residual variance, and Gaussian prediction-interval columns. Each supported backend must provide
finite probabilities in `[0, 1]` whose row sum is one. A single confidence value
can reconstruct both probabilities only for a binary outcome; PredictProR does
not fabricate unidentified probabilities for outcomes with three or more
classes.

The separate `variance_components` table is model-family aware. ML/DL
classifiers report `genetic_variance_not_identifiable`,
`residual_variance_not_identifiable`, and `heritability_not_identifiable` with
`NA` estimates. They may additionally report descriptive
`mean_prediction_confidence`, `prediction_confidence_variance`, and
`mean_classification_uncertainty`; none is a genetic or residual variance
component. The same rule applies in single-environment and MET output and is
grouped by trait/environment whenever those dimensions are present.

Bayesian categorical models do not use this ML/DL fallback. When their fitted
latent model identifies variance, their liability-scale genetic, fixed
residual, and heritability summaries remain authoritative. If any
variance-capable model fails to return its fitted component table, PredictProR
reports `model_variance_components_not_returned` rather than relabelling
prediction dispersion as genetic variance.

For binary and ordinal true prediction, RKHS and single-environment
`GBLUP_BRR` additionally expose model-implied variance components outside the
row-level prediction table. These are on BGLR's probit latent-liability scale:
the genetic variance is obtained from retained posterior variance draws and
the fitted covariance basis, while the latent residual variance is fixed at 1
by the probit model and is not estimated from the observed classes. The public
result contains `variance_components` and `variance_component_intervals`; with
`system_database = FALSE`, PredictProR also writes these as package-native CSV
files. Posterior standard deviations and credible limits describe MCMC
uncertainty in those liability-scale summaries. They do not change the
definition of row-level `Reliability` above.

MET RKHS binary and ordinal results separate genotype-main and
genotype-by-environment liability variance and report their total and implied
heritability overall and by environment. Unordered multiclass classifiers do
not receive genetic variance components because those backends do not fit this
latent genetic model.

Every table begins with the same base columns:

1. genotype ID (and environment/trait identifiers when applicable);
2. `Predicted_class`;
3. `Train_Test_Label`;
4. `Observed_class`;
5. `Prediction_confidence`;
6. `Classification_uncertainty`;
7. `Reliability`;
8. `Reliability_remarks`;
9. one `Probability_<class>` column per class.

Run the deterministic contract and routing audit from the package source:

```r
Rscript tools/audit_classification_output_consistency.R
```

Run the real backend smoke tests with:

```r
Sys.setenv(PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS = "1")
testthat::test_file("tests/testthat/test-bayesian-classification-runtime-smoke.R")
testthat::test_file("tests/testthat/test-python-ml-dl-runtime-smoke.R")
```

The Bayesian smoke chains are intentionally short integration checks and must
not be interpreted as recommended inferential MCMC settings.

## Classification cross-validation

Classification cross-validation reports metrics appropriate to the declared
response family. Binary outcomes support accuracy, balanced accuracy,
precision, recall, specificity, F1, MCC, Brier score, log loss, and ECE.
Ordinal outcomes support accuracy, balanced accuracy, mean absolute class
error, within-one-class accuracy, quadratic weighted kappa, log loss, and ECE.
Multiclass outcomes support accuracy, balanced accuracy, macro precision,
macro recall, macro F1, log loss, and ECE. Accuracy-like metrics are maximized;
class error, probability loss, and calibration error are minimized.

Cross-validation is the model-comparison stage. Supply two or more candidate
models in `GS_model_cv` and use `cv_evaluation_only = TRUE` when no final refit
is wanted. `metric_for_ranking = "auto"` uses balanced accuracy for binary
responses, macro F1 for multiclass responses, and quadratic-weighted kappa for
ordinal responses. Exact primary-metric ties are resolved by log loss for
binary/multiclass responses and by mean absolute class error followed by log
loss for ordinal responses. Users can override the primary metric and ordered
`ranking_tie_breakers`.

For binary responses, set `positive_class` explicitly whenever precision,
recall, specificity, F1, MCC, Brier score, or a one-column probability vector
is used. When omitted, the second factor level is used and recorded.

For an in-memory run `out`, users can access:

```r
out$cv_results_processed$classification_metrics$per_replication
out$cv_results_processed$classification_metrics$aggregated
out$cv_results_processed$classification_metrics$metric_directions
out$cv_results_processed$classification_metrics$confusion_matrices$per_replication
out$cv_results_processed$classification_metrics$confusion_matrices$aggregated

out$cv_results_processed$classification_model_selection$policy
out$cv_results_processed$classification_model_selection$models_compared
out$cv_results_processed$classification_model_selection$candidate_metrics
out$cv_results_processed$classification_model_selection$selected_models

out$cv_results_processed$classification_reliability$out_of_fold_predictions
out$cv_results_processed$classification_reliability$calibration_by_bin
out$cv_results_processed$classification_reliability$summary
```

The out-of-fold table contains genotype ID, observed and predicted class,
correctness, all class probabilities, `Prediction_confidence`,
`Classification_uncertainty`, and `Reliability`. Reliability remains the
maximum predicted class probability. ECE is the frequency-weighted absolute
gap between mean confidence and empirical accuracy across confidence bins.
It evaluates calibration of out-of-fold predictions; it is not a PEV-based or
genetic reliability estimate.

With `system_database = FALSE`, the same objects are written by PredictProR as
native CSV files whose names begin with
`cv_results_processed_classification_metrics_` and
`cv_results_processed_classification_reliability_`. Raw fold predictions and
probabilities remain available under `cv_results_raw/`.

## Single-trait multi-environment classification

MET classification is exposed only for models with a real environment-aware
route. Binary responses support CatBoost, LightGBM, Xgboost, RandomForest,
DenseNeuralNet, TabTransformer, TabAttention, TabNet, MixtureOfExperts,
and RKHS. Multiclass responses support the nine environment-aware
ML/DL models. Ordinal responses support RKHS; unordered ML/DL
classifiers are not relabelled as ordinal models. Kernel-Bayesian binary and
ordinal calls require an environment fixed effect and genotype main plus
genotype-by-environment random terms.

GBLUP_BRR is not exposed as a second categorical MET candidate because its
kernel eigen-square-root BRR parameterization induces the same covariance as
the RKHS route. Keeping both would let numerical/MCMC noise decide between two
equivalent model labels. Consequently ordinal MET CV currently evaluates one
unique model and records `cv_compares_models = FALSE`; it is not presented as
a model-selection comparison until another genuinely distinct ordinal MET
implementation is available.

Run CV1 and CV2 as separate model-comparison analyses. CV1 holds complete
genotypes out across environments. CV2 holds genotype-by-environment records
out while retaining those genotypes in other environments. In either case,
model selection uses the pooled out-of-fold metric recorded in the `ALL` row,
averaged across replications. The `ALL` row is not averaged again with the
individual environment rows. Per-environment metrics, confusion matrices, and
calibration summaries are diagnostic outputs and do not change the selected
model.

MET-specific diagnostic objects include:

```r
out$cv_results_processed$classification_metrics$confusion_matrices$per_environment_per_replication
out$cv_results_processed$classification_metrics$confusion_matrices$per_environment_aggregated
out$cv_results_processed$classification_reliability$calibration_by_environment_bin
out$cv_results_processed$classification_reliability$summary_by_environment
```

The MET out-of-fold prediction and probability tables retain the environment
column. The selection policy records the CV scenario, environment column,
pooled selection scope, and `ALL_rows_mean_across_replications` aggregation
rule.
