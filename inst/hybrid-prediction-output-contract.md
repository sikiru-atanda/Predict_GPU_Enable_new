# Hybrid prediction output contract

PredictProR hybrid prediction is currently a Gaussian-response workflow. A
hybrid row is identified by `Hybrid_ID`, `Female`, and `Male`, with `Env` and
`Trait` added when applicable. Classification, ordinal, nominal, and
multiclass hybrid requests are outside the supported scope and fail explicitly.

The supported engines and model names are:

- ASReml-R: `GBLUP`;
- Bayesian BGLR: `GBLUP_BRR`, `RKHS`;
- Gaussian process: `GP`, `KRR`, `LowRankGP`;
- classical ML: `RandomForest`, `Ridge_Regression`, `PartialLeastSquare`,
  `Lasso`, `SupportVectorMachine`, `Xgboost`, `CatBoost`;
- deep learning: `mlp`, `ft_transformer`.

Every true-prediction table uses this ordered base schema:

1. hybrid and parent identifiers, followed by optional environment and trait;
2. `Predicted_value`, `Train_Test_Label`, and `Observed_value`;
3. `Standard_error`, `PEV`, and their provenance fields;
4. prediction-interval and uncertainty fields;
5. prediction-stability fields;
6. `Reliability`, its remarks, exact variance input, reference variance, and
   basis.

## Model-implied uncertainty

ASReml, Bayesian, and GP hybrid models may provide model-implied prediction
uncertainty. Wherever it is identified:

```text
PEV_i = Standard_error_i^2
Reliability_i = max(0, min(1, 1 - PEV_i / Vg_i))
```

`Vg_i` is the model-implied genetic reference variance reported in
`Reliability_reference_variance`. If the bounded standard form degenerates to
the same zero for a group while PEV still varies, PredictProR uses
`Vg_i / (Vg_i + PEV_i)` for that group and writes `ratio-form fallback used`
in `Reliability_remarks`. The fallback is therefore visible rather than silent.
Reliability is never clipped upward to 0.7.

For hybrid BGLR `GBLUP_BRR`, prediction PEV and intervals use posterior target
draws. Female-GCA, male-GCA, and SCA variances are posterior means of the
training-row variance of their respective `X beta` draws; total genetic
variance uses the variance of the summed genetic target draws. For hybrid
BGLR `RKHS`, each component is the posterior `varU` draw multiplied by the
mean fitted-kernel diagonal on observed rows. Residual variance uses BGLR's
`varE` chain. The public component standard errors are posterior standard
deviations.

## ML/DL predictive uncertainty and precision

Hybrid ML/DL has no model-implied genetic reliability. True prediction reuses
the package's batched target bootstrap and held-out internal-CV calibration.
When the data identify target-specific uncertainty, PredictProR reports:

```text
PEV_i = estimated target-specific predictive MSE_i
Standard_error_i = sqrt(PEV_i)
Reliability_i = V_y / (V_y + PEV_i)
```

`V_y` is the training-phenotype reference variance. This `Reliability` is a
predictive precision index, not breeding-value or genetic reliability; that
boundary is stated in `Reliability_basis`. `Prediction_stability` remains a
separate bootstrap descriptive index. Genetic variance, residual variance,
heritability, and genetic reliability remain explicitly non-identifiable for
ML/DL models.

This predictive decomposition is selected from the declared model family, not
from the shape of the prediction table. Hybrid ASReml, Bayesian, and GP models
must retain their fitted variance components. If such a model does not return
its component table, the public result reports
`model_variance_components_not_returned`; it never falls back to the ML/DL
non-identifiability or prediction-dispersion rows.

The held-out residuals calibrate the predictive-error scale; bootstrap and
fold-to-fold target predictions preserve target-to-target differences. A
single marginal residual variance is never repeated as every hybrid's PEV. If
held-out or target-specific uncertainty cannot be identified, PEV, standard
error, intervals, and predictive precision remain unavailable rather than
being fabricated from in-sample residual dispersion.

Run the deterministic contract and routing audit from package source with:

```r
Rscript tools/audit_hybrid_prediction_consistency.R
```
