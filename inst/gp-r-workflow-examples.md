# Gaussian Process R Workflows

This note shows the public R workflow for the Gaussian-process backends in
PredictProR. It covers:

- single-environment, single-trait GP
- multi-environment, single-trait GP
- multi-trait, single-environment GP
- multi-trait, multi-environment GP

For default-tuning and release validation policy, see
`inst/gp-validation-policy.md`.

Run the GP setup once before using these functions:

```r
PredictProR::setup_predictgp_env()
```

The examples assume:

- `pheno` is a data frame with columns `GID`, `Env`, `Trait1`, and `Trait2`
- `G` is a genomic relationship matrix with row and column names matching `GID`
- `K_env` is an environment relationship matrix with row and column names matching `Env`
- `test_ids` is the character vector of genotype IDs to predict

For prediction rows, you can either pass `test_set = test_ids` or set the
response values to `NA`.

## 1. Single-Environment, Single-Trait GP

Use this when one trait is measured in one environment.

```r
pheno_se <- pheno[pheno$Env == "E1", ]

fit_se <- PredictProR::gp_single_trait_model(
  pheno_data = pheno_se,
  gmatrix = G,
  response = "Trait1",
  gen_name = "GID",
  test_set = test_ids,
  gp_return_se = TRUE
)
```

Main output:

```r
pred <- fit_se$predictions
head(pred)
```

Useful columns:

- `Prediction`
- `Prediction_SE_observed`
- `Prediction_Var_observed`
- `Prediction_SE_latent`
- `Prediction_Var_latent`

For these outputs, the variance columns should match the squared standard
errors:

```r
all.equal(pred$Prediction_SE_observed^2, pred$Prediction_Var_observed)
```

## 2. Multi-Environment, Single-Trait GP

Use this when one trait is measured across environments. Supply
`heter_groups` and either an environment relationship matrix or environment
covariates.

```r
fit_met <- PredictProR::gp_single_trait_model(
  pheno_data = pheno,
  gmatrix = G,
  response = "Trait1",
  gen_name = "GID",
  heter_groups = "Env",
  test_set = test_ids,
  env_similarity = K_env,
  include_components = c("g", "ge", "e"),
  w_g = c(0.7),
  w_ge = c(0.2),
  w_e = 0.3,
  gp_return_se = TRUE
)
```

Main output:

```r
pred_met <- fit_met$predictions
table(pred_met$Env)
```

Expected prediction rows are one row per predicted `GID x Env` cell. Useful
uncertainty columns are the same as the single-environment single-trait model:

- `Prediction_SE_observed`
- `Prediction_Var_observed`
- `Prediction_SE_latent`
- `Prediction_Var_latent`

For large MET with more than five environments, pass exactly one of
`env_similarity` or `env_covariates`. If `env_covariates` is supplied, the
backend QC-aligns the covariates to the phenotype environments, drops
non-numeric and zero-variance features, applies G2F EC feature QC by default,
standardizes the retained features, and converts them to an environment kernel
with `kenv_kernel`. That covariate-derived kernel replaces environment FA or a
separate unstructured environment covariance.

For prospective MET prediction where the deployment data may contain new
environments, use the blocked-environment tuner. It tunes only on finite
training rows, masks validation responses, and stores the selected settings in
`fit_met_tuned$info$gp_tuning`.

```r
fit_met_tuned <- PredictProR::gp_single_trait_model(
  pheno_data = pheno,
  gmatrix = G,
  response = "Trait1",
  gen_name = "GID",
  heter_groups = "Env",
  test_set = test_ids,
  env_covariates = env_covariates,
  tune_gp = TRUE,
  tuning_strategy = "blocked_environment",
  tuning_objective = "rmse",
  mean_adjustment = "auto"
)
```

For production runs, leave `gp_auto_policy = TRUE` unless you need a fully
manual benchmark. The policy records its decisions in
`fit_met_tuned$info$gp_auto_policy`, uses GPU automatically when available and
not busy, routes large factor-cache MET jobs to the operator backend, turns on
tiered dispatch for prospective environment-covariate MET prediction, and
enables `tune_gp = "auto"` only when historical validation can be done without
looking at target-year truth. With `tuning_strategy = "auto"`, known-environment
genotype prediction uses blocked-genotype validation, while prospective
new-environment prediction uses blocked-environment validation. Dense `gmatrix`
fits tune on the selected split directly. Factor-cache-only production runs
materialize a bounded dense kernel for a sampled calibration/validation subset,
then use the selected settings in the full operator fit. The default subset
limits can be overridden with
`tuning_max_calibration_rows`, `tuning_max_validation_rows`, and
`tuning_max_genotypes`, or the corresponding environment variables
`PREDICTPRO_GP_TUNING_MAX_CALIBRATION_ROWS`,
`PREDICTPRO_GP_TUNING_MAX_VALIDATION_ROWS`, and
`PREDICTPRO_GP_TUNING_MAX_GENOTYPES`.

For genotype-ranking use cases, set `tuning_objective = "within_env_spearman"`
or `tuning_objective = "centered_rmse"` so validation removes environment mean
differences before selecting kernel weights and lambda.

## 3. Multi-Trait, Single-Environment GP

Use this when two or more traits are measured in one environment.

```r
fit_mt <- PredictProR::gp_multi_trait_model(
  pheno_data = pheno_se,
  gmatrix = G,
  response = c("Trait1", "Trait2"),
  gen_name = "GID",
  test_set = test_ids,
  return_se = TRUE,
  return_trait_correlations = TRUE,
  varcomp_mode = "reml"
)
```

For large MET with more than five environments, supply `env_similarity` from
weather, soil, management covariates, FA scores, or another defensible
environment-similarity construction. Once that covariate-derived environment
kernel is supplied, do not also try to model a separate FA or unstructured
environment covariance for the same signal.

Main prediction output:

```r
pred_mt <- fit_mt$predictions
head(pred_mt)
```

Useful columns:

- `Prediction`
- `SE`
- `SE_latent`
- `PEV`
- `trait`

For multi-trait predictions, `PEV` is the latent prediction error variance:

```r
all.equal(pred_mt$SE_latent^2, pred_mt$PEV)
```

Trait covariance and correlation outputs:

```r
fit_mt$result$genetic_covariance
fit_mt$result$genetic_correlation
fit_mt$result$residual_covariance
fit_mt$result$residual_correlation
```

Use `fit_mt$result$genetic_correlation` for the estimated genetic correlation
between traits.

## 4. Multi-Trait, Multi-Environment GP

Use this when two or more traits are measured across environments.

```r
fit_mt_met <- PredictProR::gp_multi_trait_met_model(
  pheno_data = pheno,
  gmatrix = G,
  response = c("Trait1", "Trait2"),
  gen_name = "GID",
  heter_groups = "Env",
  test_set = test_ids,
  env_similarity = K_env,
  return_se = TRUE,
  return_trait_correlations = TRUE,
  varcomp_mode = "mom"
)
```

Here, pass exactly one of `env_similarity` or `env_covariates`; that input
defines the environment kernel. When `env_covariates` is provided, the backend
uses the same QC-align-standardize-convert process described above before
building `K_env`. That covariate-derived environment kernel is the environment
covariance signal; do not add a second environment FA or unstructured
environment covariance on top of it. The `trait_structure` and
`gxe_trait_structure` arguments control cross-trait covariance only.

Main prediction output:

```r
pred_mt_met <- fit_mt_met$predictions
head(pred_mt_met)
```

Expected prediction rows are one row per predicted `GID x Env x trait` cell.
Useful columns include:

- `Prediction`
- `SE`
- `SE_latent`
- `PEV`
- `gid`
- `env`
- `trait`

Correlation outputs:

```r
fit_mt_met$result$genetic_correlation
fit_mt_met$result$gxe_correlation
fit_mt_met$result$residual_correlation
```

Use:

- `genetic_correlation` for cross-trait genetic correlation
- `gxe_correlation` for cross-trait genotype-by-environment correlation
- `residual_correlation` for cross-trait residual correlation

## Optional SoyNAM Data Pattern

The CRAN `SoyNAM` package can be used to build a small real-data example. The
same four workflows above apply after creating a phenotype table and genomic
relationship matrix.

```r
install.packages("SoyNAM")
library(SoyNAM)
data(soynam)

pheno_raw <- data.line
markers <- gen.raw

traits <- c("yield", "protein")
envs <- names(sort(table(pheno_raw$environ), decreasing = TRUE))[1:3]

pheno_soy <- pheno_raw[pheno_raw$environ %in% envs, c("strain", "environ", traits)]
pheno_soy <- pheno_soy[complete.cases(pheno_soy), ]
pheno_soy <- stats::aggregate(
  pheno_soy[traits],
  pheno_soy[c("strain", "environ")],
  mean,
  na.rm = TRUE
)

gids <- intersect(unique(pheno_soy$strain), rownames(markers))
gids <- head(gids, 50)
pheno_soy <- pheno_soy[pheno_soy$strain %in% gids, ]

X <- markers[gids, seq_len(500), drop = FALSE]
X <- scale(X, center = TRUE, scale = FALSE)
G_soy <- tcrossprod(X)
G_soy <- G_soy / mean(diag(G_soy))
rownames(G_soy) <- colnames(G_soy) <- gids

K_env_soy <- diag(length(envs))
rownames(K_env_soy) <- colnames(K_env_soy) <- envs

test_ids <- tail(gids, 10)

fit_soy_mt_met <- PredictProR::gp_multi_trait_met_model(
  pheno_data = pheno_soy,
  gmatrix = G_soy,
  response = traits,
  gen_name = "strain",
  heter_groups = "environ",
  test_set = test_ids,
  env_similarity = K_env_soy,
  return_se = TRUE,
  return_trait_correlations = TRUE
)
```

The prediction table is in:

```r
fit_soy_mt_met$predictions
```

The multi-trait, multi-environment genetic and GxE correlations are in:

```r
fit_soy_mt_met$result$genetic_correlation
fit_soy_mt_met$result$gxe_correlation
```

## Real-Data Accuracy Benchmark

For a reproducible real-data accuracy smoke on both BGLR wheat and SoyNAM, run:

```r
source("tools/gp_real_data_accuracy_benchmark.R")
```

The benchmark fits all four GP scenarios and writes:

```r
tools/tmp_gp_real_data_accuracy_results.csv
```

The output columns include:

- `dataset`
- `scenario`
- `trait`
- `source_url`
- `n_test`
- `rmse`
- `mae`
- `pearson_cor`
- `relative_rmse`
- `mean_se`
- `mean_pev`
- `genetic_corr_12`
- `gxe_corr_12`

Note: BGLR wheat has one yield trait measured in multiple environments. For
the two multi-trait wheat checks, the benchmark recasts environment columns as
pseudo-traits. SoyNAM is the true multi-trait check because it includes yield,
protein, and oil phenotypes.

Use this benchmark as the first public-data generalization gate. The current
open validation tiers are:

- BGLR wheat: CIMMYT wheat lines with marker data and yield measured in four
  environments, available through the CRAN `BGLR` package.
- SoyNAM: soybean genomic and multi-environmental phenotype data, available
  through the CRAN `SoyNAM` package.
- G2F competition data: maize genotype-by-environment prediction competition
  data and full public yearly releases from the Genomes to Fields resources.
- learnMET, optional: GitHub package with G2F-derived maize MET toy data plus
  rice MET data and environmental covariates. Install this only for an
  expanded validation run because it pulls GitHub dependencies.

## Repeated CV Comparison With GBLUP

For a repeated grouped cross-validation comparison between GP and a
GBLUP-equivalent kernel BLUP baseline, run:

```r
source("tools/gp_vs_gblup_cv_benchmark.R")
```

By default this runs 5 folds x 5 repetitions for:

- BGLR wheat single-environment yield
- BGLR wheat multi-environment yield
- SoyNAM single-environment yield
- SoyNAM multi-environment yield

The GBLUP baseline uses the same genomic relationship matrix and the same
genotype folds as GP. Its residual-to-genetic shrinkage parameter is selected
inside each training fold by a small inner CV grid, so the GP and GBLUP scores
are compared on the same held-out genotypes.

The benchmark writes:

```r
tools/tmp_gp_vs_gblup_cv_details.csv
tools/tmp_gp_vs_gblup_cv_summary.csv
```

The repeated-CV outputs include `source_url`, `gp_auto_policy`,
`gp_tuning_status`, and `gp_tuning_reason` columns. For BGLR wheat and the
basic SoyNAM grouped-genotype folds, the package auto-policy usually records
`manual_or_default_settings_kept` and `gp_tuning_status = "not_requested"`
because these folds use fixed `env_similarity` or single-environment genomic
prediction. Environment-covariate auto tuning is validated by:

```r
source("tools/gp_public_envcov_auto_cv_benchmark.R")
```

This benchmark writes:

```r
tools/tmp_gp_public_envcov_auto_cv/public_envcov_auto_cv_detail.csv
tools/tmp_gp_public_envcov_auto_cv/public_envcov_auto_cv_summary.csv
```

The public env-covariate benchmark currently runs SoyNAM and a G2F historical
subset. Its output includes `gp_tuning_strategy` so it is clear whether the
policy used blocked-genotype validation for known environments or
blocked-environment validation for prospective environments. The full 2024 G2F
competition-style script remains the prospective, factor-cache validation.

You can reduce runtime while testing the script by setting:

```r
Sys.setenv(PREDICTPROR_GP_CV_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_CV_REPS = "1")
```

## Deeper Kernel And Weight Tuning

For a wider GP tuning pass that compares linear, RBF, polynomial, blended, and
low-rank genomic kernels, expands the lambda grid, and retunes
multi-environment G/GxE/environment weights, run:

```r
source("tools/gp_deep_tuning_benchmark.R")
```

By default this runs 5 folds x 3 repetitions with 3-fold inner CV. It writes:

```r
tools/tmp_gp_deep_tuning_details.csv
tools/tmp_gp_deep_tuning_summary.csv
tools/tmp_gp_deep_tuning_selection_counts.csv
```

Useful runtime controls:

```r
Sys.setenv(PREDICTPROR_GP_DEEP_CV_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_DEEP_CV_REPS = "1")
Sys.setenv(PREDICTPROR_GP_DEEP_INNER_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_DEEP_MIN_KERNEL_GAIN = "0.02")
```

For same-environment multi-environment prediction, `gp_single_trait_model()`
uses the environment column as a fixed effect when every environment is present
in the training rows. For CV0-style prediction into completely unseen
environments, pass explicit `fixed_effects = character(0)` and supply
environment covariates or a defensible environment-similarity matrix.

## Multi-Trait GP Benchmark

For repeated CV on multi-trait single-environment and multi-trait
multi-environment GP models, run:

```r
source("tools/gp_mt_deep_tuning_benchmark.R")
```

The benchmark compares:

- trait mean or trait-by-environment mean
- single-trait GBLUP per trait
- default multi-trait GP
- kernel-tuned multi-trait GP

It reports per-trait RMSE, MAE, Pearson correlation, relative RMSE, mean SE,
mean PEV, SE/RMSE ratio, genetic correlation, GxE correlation where available,
selected kernel, and runtime. By default it runs 3 folds x 2 repetitions with
3-fold inner tuning and writes:

```r
tools/tmp_gp_mt_deep_tuning_details.csv
tools/tmp_gp_mt_deep_tuning_summary.csv
tools/tmp_gp_mt_deep_tuning_selection_counts.csv
```

Useful runtime controls:

```r
Sys.setenv(PREDICTPROR_GP_MT_CV_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_MT_CV_REPS = "1")
Sys.setenv(PREDICTPROR_GP_MT_INNER_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_MT_KERNELS = "linear_full,rbf_bw0.5_full")
```

## Multi-Trait REML, FA, and Calibration Diagnostic

To check multi-trait uncertainty calibration and compare unstructured versus
factor-analytic covariance fits, run:

```r
source("tools/gp_mt_reml_fa_calibration_diagnostic.R")
```

The diagnostic runs multi-trait single-environment and multi-trait
multi-environment folds on BGLR wheat and SoyNAM. It reports RMSE, MAE, bias,
Pearson correlation, relative RMSE, mean SE, mean PEV, SE/RMSE ratio, 68% and
95% interval coverage, normalized-error diagnostics, genetic correlation, GxE
correlation, runtime, and any fitting error. It writes:

```r
tools/tmp_gp_mt_reml_fa_calibration_details.csv
tools/tmp_gp_mt_reml_fa_calibration_summary.csv
```

Useful runtime controls:

```r
Sys.setenv(PREDICTPROR_GP_MT_DIAG_CV_FOLDS = "2")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_CV_REPS = "1")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_WHEAT_N = "40")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_SOYNAM_N = "48")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_SOYNAM_MARKERS = "128")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_SINGLE_MAX_ITER = "16")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_MET_REML_MAX_ITER = "12")
Sys.setenv(PREDICTPROR_GP_MT_DIAG_MET_MOM_PCG_MAX_ITER = "160")
```

The practical diagnostic run used 3 folds x 1 repetition with 60 wheat lines,
72 SoyNAM lines, and 192 SoyNAM markers. The package default is therefore:
multi-trait multi-environment GP uses the scalable MoM path with unstructured
cross-trait covariance and an explicit environment kernel supplied through
`env_similarity` or `env_covariates`. For SoyNAM multi-trait
single-environment, REML unstructured completed all folds and was better than
FA1, while FA1 failed in one fold with a non-positive-definite Cholesky
factorization. For SoyNAM multi-trait multi-environment, MoM unstructured and
MoM FA1 were effectively tied for protein, REML unstructured was slightly
better for yield but worse for protein.

Use dense REML unstructured as a small-data diagnostic or refinement path, not
as the default fast path:

```r
fit_reml_us <- gp_multi_trait_met_model(
  pheno_data = pheno,
  gmatrix = G,
  response = c("yield", "protein"),
  gen_name = "gid",
  heter_groups = "environ",
  test_set = holdout_gids,
  env_similarity = env_similarity,
  varcomp_mode = "reml",
  trait_structure = "unstructured",
  gxe_trait_structure = "unstructured"
)
```

Use cross-trait FA only with `varcomp_mode = "mom"` and only when covariance
clipping or unstable cross-trait correlations are observed and a local
benchmark shows an advantage. The R wrapper rejects dense MET REML with FA
cross-trait covariance because that backend currently supports only
unstructured trait and GxE covariance. For large MET with more than five
environments, first build a covariate-derived environment kernel and pass it as
`env_similarity` or pass the covariates as `env_covariates` so the backend can
QC-align them and convert them to `K_env`; that replaces the need for
environment FA or unstructured environment covariance.
