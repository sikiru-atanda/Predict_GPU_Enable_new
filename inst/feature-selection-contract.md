# Feature-selection contract

PredictProR scores and selects predictors separately for every response trait.
The scores are predictive feature-ranking quantities. They are not genetic
variance components, heritabilities, or causal effects.

## Choosing the scoring policy

Use `feature_scoring_cv = "fold_internal"` when cross-validation performance is
the target. The feature scorer is refitted inside every outer training fold and
never sees that fold's held-out responses. The returned
`feature_selection_metadata` records the trait, fold, replication, training
sample count, training-ID hash, score, rank, and selected flag.

Use `feature_scoring_cv = "fixed"` for a ranking established outside the
assessment data, for exploratory candidate curves, or for the final refit on all
available training observations. If a fixed ranking is estimated from the same
responses used in cross-validation, its CV metric is not leakage-free.

For standard single-trait and MET routes, `feature_scoring_cv = "both"` returns
fixed and fold-internal tasks separately. Their metrics are not pooled, and the
fold-internal result is used for automatic model choice. Specialized hybrid and
joint multi-trait routes require separate runs for the two policies.

On standard single-trait and MET routes every CV grid also contains the full
predictor count (k = all predictors, no selection), so the comparison and the
automatic model choice can conclude that selection does not help. That
candidate is not rescored inside the folds, because keeping every predictor
cannot change the set.

If several `feature_k_grid` values are assessed, each candidate's fold-internal
ranking is honest for that pre-specified k. Choosing the best k and quoting the
same outer-CV value as an unbiased final performance estimate is still a tuning
selection and is not nested CV. Specialized routes therefore require one
pre-specified k per run.

## Reusing a selected predictor set

Use `feature_selected` when the predictors were selected previously and should
be reused without rescoring. It is an exact set, not a ranking or an instruction
to choose k. Accepted forms are:

- a character vector for one common set;
- a list keyed by trait for different sets per trait;
- a list keyed by `geno_data`, `omic1_data`, `omic2_data`, or `omic3_data`;
- a nested trait/source list in either orientation; or
- a data frame with `predictor` and optional `trait`, `source_block`, and
  logical `selected` columns.

For example:

```r
selected <- list(
  Yield = list(geno_data = c("M12", "M44"), omic1_data = "RNA7"),
  Protein = list(geno_data = c("M3", "M91"))
)

fit <- model_execute(
  pheno_data = phenotype,
  geno_data = markers,
  omic1_data = transcriptome,
  response = c("Yield", "Protein"),
  gen_name = "ID",
  GS_model_cv = c("RandomForest", "BayesA", "GBLUP_BRR"),
  cross_validation = TRUE,
  feature_selected = selected
)
```

Separate trait fits use their own sets. A genuinely joint multi-trait fit uses
the union of the declared trait sets because it must have one shared predictor
space. In a source-keyed set, a source that is not named is excluded; list every
source that should remain active.

Named raw marker/omics matrices let PredictProR rebuild marker designs and
kernels from the exact set. A precomputed kernel has no feature-to-column map.
If it is supplied with `feature_selected`, PredictProR assumes the user rebuilt
that kernel externally from the declared predictors; this cannot be verified
from the kernel alone. Explicit runs report `selection_mode = "explicit"` and
`scoring_model = "user_supplied"` in selection provenance.

## Trait and model integration

- Single-trait models use that trait's selected predictors.
- Separate multi-trait fits retain separate selected sets for every trait.
- A genuinely joint multi-trait fit must use one predictor space, so it uses the
  union of the trait-specific top-k sets. Metadata contains each trait's set and
  the joint-union size.
- Classical ML and DL models receive the selected predictor columns directly.
- Marker-regression Bayesian models rebuild their design terms from the selected
  raw columns.
- GBLUP, kernel Bayesian, GP, and kernel-based MET models rebuild relationship
  matrices from the selected raw marker/omics columns for every applicable
  trait and fold.
- Hybrid ML/DL models can use automatic selection when their feature matrix is
  keyed by hybrid ID.
- Parent-kernel hybrid models cannot infer an outcome-based ranking directly:
  phenotype rows represent hybrids while marker rows represent parents. They
  accept `feature_selected` or externally calculated fixed score metadata from
  an explicitly documented hybrid feature encoding; fold-internal automatic
  scoring fails instead of silently inventing an encoding.

Precomputed kernels do not retain a feature-to-column map and cannot be subset
by predictor name. Automatic selection therefore requires the corresponding raw
marker/omics matrices and rebuilds the kernels. A run that combines automatic
selection with an incompatible user-supplied precomputed kernel stops with an
actionable error.

## Response families

Ridge and BayesB feature scoring are defined only for Gaussian responses.
Binary, ordinal, nominal, and multiclass traits must use
`feature_scoring_model = "RandomForest"`, which uses a classifier and avoids
rankings that depend on arbitrary factor-level integer codes.

## Example

```r
fit <- model_execute(
  pheno_data = phenotype,
  geno_data = markers,
  response = c("Yield", "Protein"),
  gen_name = "ID",
  GS_model_cv = c("RandomForest", "BayesB", "GP"),
  cross_validation = TRUE,
  cross_validation_meth = "K-Folds",
  feature_scoring = TRUE,
  feature_scoring_model = "Ridge_Regression",
  feature_scoring_cv = "fold_internal",
  feature_k_grid = c(100L, 500L)
)

fit$cv_results_processed$feature_selection_metadata
```

For a nominal or ordinal response, change the scorer to `"RandomForest"`.
