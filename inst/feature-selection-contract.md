# Feature-selection contract

Ordinary multi-response runs score each response separately. Joint multi-trait
routes score the numeric Gaussian responses together and use one shared set.
The scores are predictive feature-ranking quantities. They are not genetic
variance components, heritabilities, or causal effects.

## Choosing the scoring policy

The default, `feature_scoring_cv = "auto"`, resolves before any route runs: to
`"fold_internal"` when PredictProR scores the features (`feature_scoring =
TRUE`), and to `"fixed"` when you supply `feature_score_metadata`, an external
ranking the package cannot recompute per fold (recorded as `ranking_source =
"user"`). The resolved policy is recorded in the CV output.

Use `feature_scoring_cv = "fold_internal"` when cross-validation performance is
the target. The feature scorer is refitted inside every outer training fold and
never sees that fold's held-out responses: it scores exactly the fold's training
records (in multi-environment CV0/CV2 a line's held-out records are not scored,
and a hybrid scenario that keeps further rows out of training, such as
`Hybrid_Both_New_Parents`, scores only its training rows). The returned
`feature_selection_metadata` records the trait, fold, replication, training
sample count, training-ID hash, score, rank, and selected flag. For standard separate-trait routes, each fold is scored once per trait and
replication; its models and `feature_k_grid` values reuse that ranking. Joint
ML has a separate scoring problem for each target trait and fold, reflecting
which auxiliary traits are known; ASReml and DL score once per fold.

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

## How predictors are scored

- **Design effects are adjusted for.** Scoring treats the terms of `fixed` and
  the environments in `heter_groups`, reducing rankings driven by design
  confounding:
  Ridge residualizes the trait and the markers on them (exactly ridge with
  unpenalized covariate effects), BayesB adds them as a fixed-effect term (with
  the markers standardized within them, as for ridge), and RandomForest scores
  the residualized trait (Gaussian) or uses them as unranked predictors.
  Specialised routes (`hybrid_*`, `multi_trait_*`) adjust for environments
  only. Records with a missing covariate are left out and counted
  (`n_records_dropped`); `covariates_adjusted` lists the terms adjusted for. A
  `fixed` term base R cannot evaluate, such as ASReml's `at()`, is left out of
  the adjustment with a warning. With `fixed = ~1` and one environment, nothing
  changes.
- **RandomForest scores are the actual impurity reduction (AIR).** Each
  predictor's impurity decrease minus that of a shadow copy whose rows share
  one permutation (as ranger's `impurity_corrected`). Plain impurity importance
  puts a continuous predictor near the top whether or not it matters; AIR
  removes that bias (its mean is about 0), but a continuous predictor's AIR is
  more variable than a marker dosage's, so it still reaches the top ranks by
  chance more often (in simulation, the top 4 of 31 in about 22% of data sets,
  against about 3% for a null marker). AIR scores can be negative; ranks follow
  the raw score.
- **Ties at the cut-off** are broken by predictor name, recorded as
  `tied_at_cutoff`, and the first one in a run warns.

## Joint multi-trait scoring and selection

`feature_score_predictors(..., joint = TRUE)` requires two or more distinct
numeric Gaussian traits. Trait names cannot contain `+`, which joins the joint
label (for example `Yield+Protein`). The default `joint = FALSE` keeps the
separate-trait behavior. Public joint Bayes, GP, ASReml, ML and DL routes request
joint scoring automatically; single-response hybrid ML/DL continue to score one
trait. A parent `multi_trait` comparison forwards the same rules to its children.

`Ridge_Regression` uses BGLR Multitrait BRR with estimated marker covariance and
shrinkage; `ridge_lambda` is ignored in this joint engine. `BayesB` uses
Multitrait SpikeSlab with trait-specific inclusion, correlated Gaussian slabs,
and inclusion prior probIn = 0.5, counts = 2. It does not impose group inclusion.
This differs from the separate-trait BGLR BayesB prior/model; their common public
label does not make the two scoring procedures matched controls.
Both fit unstructured marker and residual covariance. The implementation follows
the available BRR and SpikeSlab models described by
[Perez-Rodriguez and de los Campos (2022)](https://doi.org/10.1093/genetics/iyac112).

Each trait is scaled by its visible design-adjusted phenotypic SD, and markers
by their visible adjusted SD. The Bayesian score is the row norm of posterior
mean standardized coefficients: sqrt(sum_t E(beta_jt | data)^2). This differs
from E(sqrt(sum_t beta_jt^2) | data). Changing measurement units does not change
the intended ranking; arbitrary linear transformations or duplicating correlated
traits can change it. Higher-heritability traits may contribute more under this
scaling, without a marker-by-marker guarantee. Under LD, genetic variance contains
cross-marker terms, beta' Cov(X) beta. Squared scores are not additive explained
variance and cannot be summed to obtain heritability.

Prior scales are explicitly calibrated from the visible, scaled responses after
adjusting design means; likelihood responses retain their design effects and
FIXED terms. Prior modes use half the response covariance for markers and half for
residuals, after dividing the marker mode by MSx. A pairwise calibration correlation with minimum eigenvalue below
1e-6 is shrunk toward identity to that floor. This warns and reports the shrinkage
and original minimum eigenvalue; it does not constrain posterior correlations.
The marker inverse-Wishart df0 is T+1 with S0 = 0.5*S_pd/MSx*(df0+T+1); the
residual df0 is max(5,T) with S0 = 0.5*S_pd*(df0+T+1). These depend on the visible
panel and constitute an empirical prior-mode calibration, not a specification of
prior means (which need not exist) or fixed posterior heritability. Mean-imputed startup
response covariance that is singular or numerically singular stops instead of
adding artificial observations. FIXED effects have diffuse Gaussian variance
1e6, so arbitrary environment-offset invariance is not claimed.

BGLR retains a record with any observed response and samples missing traits.
Scoring assumes conditional independence of records given markers and design
means, with no extra line effect or GxE term. Repeated records do not create
independent genetic information. Ignorable missingness requires MAR conditional
on the modeled observed information, distinct missingness parameters, appropriate
prior factorization and a suitable response model; MAR alone does not guarantee
unbiased results. Without co-observation of a trait pair there is no direct
residual cross-covariance likelihood contribution; covariance priors can still
couple posterior inference. `min_pair_coobserved` reports the minimum visible
pair count.

A zero standardized marker column is integrated at centered prior mean zero:
score 0, `effect_sd = NA`, `scoring_status = "no_predictor_variation"`. All-missing
columns and a panel without informative predictors stop. Constant columns still
count in p for requested k. For estimated Bayesian effects, `effect_sd` is
sqrt(trace posterior covariance of the coefficient row), the square root of the
sum of marginal variances. It is not the posterior SD or credible interval of
the nonlinear score, its Monte Carlo standard error, or a convergence diagnostic.

Joint RandomForest scores only complete-response records. It warns with the
incomplete count and records dropped/unmatched/all-missing record counts. Each
trait is scaled and residualized on the design matrix, and the response matrix
is fitted by one forest. sklearn's squared-error
[criterion](https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/tree/_criterion.pyx)
averages over outputs; a sum rescales every candidate by the same factor and
leaves split choices and ranking unchanged. AIR uses one shared permutation for
shadow predictors, returns signed importance and reduces impurity bias; it does
not remove LD proxy competition, establish causality or provide a valid p-value.
`effect_sd` is NA. Complete-case scoring targets the retained subset unless the
retention is representative; finite-sample unbiasedness is not promised. The
scoring forest is a separate fit from any downstream prediction forest.

Both Bayesian joint engines use `feature_bayes_nIter`, `feature_bayes_burnIn` and
`feature_bayes_thin` (defaults 1500/500/5). Defaults and a retained-draw count do
not certify convergence or reproducible cut-off ranks. Compare multiple seeds
and diagnose mixing; more iterations can help, while thinning is not a guarantee
of effective sample size. Three-seed Jaccard overlap is a descriptive stability
diagnostic; 0.9 is a heuristic rather than a statistical certification threshold.

On one prespecified n300/p3000 LD panel (top100, three seeds), default marker
Jaccards were .379–.429 for joint BRR and .379–.493 for joint SpikeSlab. At
12000/2000/5 they were .709–.739 and .325–.770 respectively. Neither met .9, and
the longer SpikeSlab chains did not consistently stabilize selection. These are
dependent diagnostics on one panel, not convergence or population-level
stability guarantees. Treat short default rankings as exploratory and examine
multiple seeds, chains, cutoff sensitivity and LD-block overlap for the data at hand.


A supplied ranking must have `scoring_scope = "joint"`, one label matching the
response set (either trait order), unique predictor names, finite raw scores,
and unique positive integer ranks ordered by score with name-based ties. It must
rank every column of the post-QC model matrix. Per-trait rankings, partial rankings
and mismatched labels stop with an actionable message. Explicit shared sets carry
`selection_mode = "explicit"` and `scoring_model = "user_supplied"`; rank and
raw_score are NA because no scoring was done. Omitted sources in a source list
are excluded. Positional masks or indices are resolved against their input columns.

Selected true-prediction model payloads (`fit$model_results` and their by-model/
by-trait views) expose `feature_selection_metadata`, `feature_selected` and integer
`feature_k`. The results envelope also retains `feature_selected`. The selected name vector is passed to
native `feature_selected.RData` saving. Native exports also retain the full joint
ranking as `feature_selection_metadata.csv` and the scoring table as
`feature_score_metadata.RData` when scores were supplied or computed. Obsolete by-trait/union selection fields
are removed. Source-specific raw columns rebuild their corresponding kernels;
precomputed kernels alone cannot verify their originating selected predictors.

## Joint CV visibility and supported routes

For ML, each target-trait fold hides that trait on the test lines while retaining
observed auxiliary traits on those lines. The scorer sees all visible traits;
the target-specific fit keeps its existing training rows and training-only
auxiliary imputation. Marker k excludes auxiliary response features. This is
prediction with known auxiliary traits, distinct from new-line prediction where
all target traits are unavailable.

For single-environment ASReml and DL, every response on a test line is masked
before scoring; ASReml fits the same masked phenotype table, and DL learns target
scaling from its training rows. Selected marker columns are applied in the same
order to training and prediction targets. Missing responses are never filled
from held-out observations for scoring. The training-record hash identifies
retained records and can be the same across ML folds because auxiliary traits
keep every ID visible; it is not a fingerprint of trait-cell visibility.

Joint specialized CV accepts one positive whole k or `"all"` and one policy per
run. `fold_internal` masks according to the route; `fixed` reuses the supplied
ranking and is optimistic if that ranking was learned from the assessment
phenotypes. Selecting k or a model with the same CV score remains a tuning
selection and needs an independent or nested assessment for a final performance
claim. Scored joint GP/Bayes CV and grouped MT-MET ASReml CV selection remain
unsupported and stop. True prediction supports joint Bayes/GP/ASReml/ML/DL,
including adjusted MT-MET scoring for the kernel routes.

## Estimates after selection

Markers chosen with the same phenotypes the model is then fitted to can inflate
genetic variance and heritability estimates and make prediction-error variances,
standard errors and reliabilities optimistic. A fit on selected markers (k below
the number of predictors, or a `feature_selected` set) therefore carries a
`Feature_selection_remark` column on its variance-component and prediction
tables and on its covariance and correlation tables (a matrix carries it as an
attribute), in the returned list and in the exported CSV files. The wording
depends on where the selection came from: package scoring of these phenotypes,
a supplied `feature_score_metadata` ranking, or a supplied `feature_selected`
set (for the last two the caution is conditional on whether these phenotypes
were used). Runs without selection get no column. CV tables (`cv_*`) are not
marked.

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

Separate trait fits use their own sets. A joint multi-trait fit requires one
shared set, or identical sets for every response (ignoring order but preserving
source identity). Differing or missing per-trait sets stop before fitting. In a source-keyed set, a source that is not named is excluded; list every
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
- Joint multi-trait fits use one joint ranking and exactly min(k,p) shared
  predictors after QC, with all p retained for `feature_k = "all"`. They do not
  take the union of separate rankings.
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
  Multi-response parent-kernel hybrid GP requires a shared `feature_selected`
  set (or identical per-trait sets). Supplied score tables stop on that route;
  the first response's ranking is never silently reused for other responses.

RandomForest rankings put zero-variation predictors last, retaining their finite
zero AIR scores and `no_predictor_variation` status. Informative predictors retain
their signed AIR scores, including negatives. New ties use radix name order;
valid supplied tied ranks retain their order across locales. This eligibility
rule applies to both single-trait and joint scoring and can change selected sets.

A supplied joint ranking must cover exactly the fitted post-QC predictor matrix.
Mismatch errors show both counts and at most ten example names. Score the post-QC
matrix, or explicitly disable genotype QC when intentionally fitting the raw matrix.
Joint selection tables can be passed back as `feature_selected`; only selected
rows are used, and the joint trait label may list responses in either order.

Fold-internal joint selection still computes rankings when k retains all columns.
This preserves score metadata and the established score-based column order, which
can affect downstream stochastic fits. Skipping scoring is a separate optimization
requiring an explicit ordering/output contract and compatibility validation.

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
