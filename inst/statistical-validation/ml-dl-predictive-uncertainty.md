# ML/DL predictive-uncertainty validation

## Claim boundary

For Gaussian machine-learning and deep-learning predictions, PredictProR does
not estimate genetic variance, residual variance, heritability, or
breeding-value reliability from fitted predictions. Those quantities remain
unavailable unless a model with an explicit genetic/random-effects estimand
provides them.

The compatibility column `PEV` is an estimated predictive mean squared error
(MSPE), not mixed-model prediction-error variance. `PEV_basis` records this
distinction. `Prediction_stability` is a bounded descriptive index, and the
ML/DL `Reliability` field is its plotting-compatible alias.

Beginning with 0.20.60, the public `Reliability` column is populated again so
existing user plots continue to work. For ML/DL it is exactly a compatibility
alias of `Prediction_stability`, not breeding-value or genetic reliability.
Version 0.20.61 restores the original marker-adjustment calculation:

```text
Reliability_i = max(0, min(1, 1 - U_i / V_yhat))
U_i           = SE_resampling,i^2
V_yhat        = var(yhat) across all final prediction rows (Train and Test).
```

`Reliability_variance_input` reports `U_i` and
`Reliability_reference_variance` reports `V_yhat`, so the score can be
recomputed from every output row. `Reliability_basis` states which uncertainty
source was used and repeats the non-genetic claim boundary.

## Estimation

Single-trait ML uses cross-fitted residuals. When fold predictions at each
target are available, signed held-out residuals are combined with fold target
predictions to estimate target-specific predictive MSE. Otherwise a marginal
cross-fitted MSE is reported for all targets. Multi-trait ML/DL estimates a
separate cross-fitted MSE within each trait.

Intervals use held-out absolute residuals and the finite-sample corrected order
statistic. They are prediction intervals for future outcomes, not confidence
intervals for a latent genetic value. Their supported claim is empirical
marginal coverage under exchangeability. They do not promise conditional
coverage for every genotype, family, environment, or uncertainty stratum.

No ML/DL PEV or interval is produced when neither held-out calibration nor a
model-based uncertainty distribution is available. In-sample residuals are not
silently promoted to PEV.

For the plotting-only reliability alias, the original ML/DL marker-adjustment
rule is kept separate from predictive PEV:

1. `U_i` is the row-level variance of resampled predictions, equivalently the
   squared resampling standard error used by the original script.
   `V_yhat` is the variance of fitted predictions across the complete final
   prediction table. Both `Train` and `Test` rows contribute when both are
   present; `Train_Test_Label` identifies rows but does not subset `V_yhat`.
2. When repeated prediction draws are unavailable, the package uses the legacy
   descriptive fallback. Within
   each environment, it starts from in-sample residual mean square, scales it
   by the clipped prediction-distance factor
   `1 + (prediction - median_prediction)^2 / prediction_variance`, and, for
   rows with an observation, averages that result with the row's squared
   residual. Fallbacks based on residual variance or observed/predicted
   variance are used only when the preceding quantity is unavailable.

Both rules define an uncalibrated visualization heuristic. Neither input is
silently written to inferential `PEV`, neither creates a prediction interval,
and the score must not be interpreted as a probability of correct selection,
genetic reliability, or a coverage guarantee. There is no imposed training
floor: training values can be low when resampled fits are unstable or the model
produces little fitted-value spread.

Gaussian ML/DL diagnostic plots consume the same held-out calibration object
as the prediction table. The exported diagnostic helper requires that
calibration by default. Its legacy bootstrap/standard-error display is
available only through an explicit opt-in and does not claim nominal coverage.

## Reproducible simulation gate

Run from the package root:

```powershell
Rscript tools\validate_ml_dl_predictive_uncertainty.R 100
```

The committed run used 100 independently generated training/test datasets per
scenario. Each training dataset used five-fold cross-fitting, and every
replicate was evaluated on 400 newly generated test outcomes.

| Scenario | Coverage | Monte Carlo SE | Mean PEV / test MSE | Mean squared standardized error |
|---|---:|---:|---:|---:|
| Exchangeable homoskedastic | 0.95480 | 0.00186 | 1.05185 | 0.96944 |
| Exchangeable heteroscedastic, pooled | 0.95322 | 0.00164 | 1.03083 | 0.99455 |
| Group-calibrated high-noise trait | 0.94840 | 0.00245 | 0.99047 | 1.03247 |
| Group-calibrated low-noise trait | 0.95173 | 0.00193 | 1.01160 | 1.00487 |

The pooled heteroscedastic scenario demonstrates the limitation of marginal
coverage: coverage was 0.99635 in the lower-noise half and 0.91010 in the
higher-noise half. Trait-specific calibration restored coverage close to 0.95
in both simulated traits. Scientists should therefore validate by family,
population, site, year, environment, trait, or other design-relevant group
whenever those groups define the intended deployment population.

## Backend integration checks

The helpers were also exercised through the actual Python-backed public model
paths on fixed-seed Gaussian simulations:

| Backend | Training/Test rows | Test coverage | Wilson 95% interval | Mean PEV / test MSE |
|---|---:|---:|---:|---:|
| Ridge | 120 / 1,000 | 0.9390 | 0.9224--0.9522 | 0.9257 |
| MLP, two-epoch smoke configuration | 80 / 400 | 0.9425 | 0.9152--0.9614 | 0.8832 |

Both confidence intervals contain the nominal 0.95 target. These fixed-seed
checks establish integration and scale coherence; the repeated simulation gate
above, rather than either single run, supports the marginal-coverage claim. The
multi-trait Ridge and MLP BGLR-wheat runtime checks also passed with finite
trait-wise cross-fitted MSE and prediction intervals. Their finite
`Reliability` values are the explicitly non-genetic prediction-stability alias
defined above.

## User-side validation

`validate_predictive_uncertainty()` evaluates independent Test/Validation rows
and reports:

- predictive MSE and RMSE;
- mean PEV, PEV-to-MSE ratio, and mean squared standardized error;
- empirical interval coverage with Wilson binomial confidence limits;
- interval width;
- stability-versus-squared-error rank association; and
- the same diagnostics within requested groups.

The function fails by default when held-out evaluation rows cannot be
identified. Passing training predictions as validation data would invalidate
the assessment.

## Methodological references

- Barber RF, Candes EJ, Ramdas A, Tibshirani RJ. Predictive inference with the
  jackknife+. *Annals of Statistics* 49(1), 2021. DOI: 10.1214/20-AOS1965.
- Romano Y, Patterson E, Candes EJ. Conformalized Quantile Regression. NeurIPS,
  2019.
- Bengio Y, Grandvalet Y. No Unbiased Estimator of the Variance of K-Fold
  Cross-Validation. *JMLR* 5, 2004.
