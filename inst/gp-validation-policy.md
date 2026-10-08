# Gaussian Process Validation Policy

This note records how PredictProR should validate Gaussian-process defaults and
tuning policy before changing them. The goal is not to make GP win every local
benchmark. The goal is to keep GP accurate, fast, uncertainty-aware, and
general enough for new crops, locations, countries, and trial designs.

## Principles

- Do not tune defaults to one dataset, one crop, or one competition split.
- Compare GP against both a simple baseline and a linear genomic baseline
  such as GBLUP whenever both are available.
- Treat external datasets that disagree with each other as useful evidence,
  not as failures to be optimized away.
- Keep variance and standard-error output under regression checks when it is
  requested and reported.
- For known-environment genotype holdout, preserve the GP signal by default;
  do not automatically blend predictions back to an environment mean unless
  validation shows a clear baseline-dominant gain.
- For prospective all-new environments, require environment covariates or an
  environment similarity kernel for proper environment transfer.
- When `env_covariates` are supplied, route them through QC into an
  environment similarity kernel. Do not add FA or unstructured environment
  covariance on top of that path just because the MET is large.
- For large MET with more than about five environments and no covariate-derived
  kernel, avoid unstructured covariance as a default because it is fragile and
  expensive. Prefer scalable low-rank or factor-analytic alternatives.
- Prefer GPU execution when an available GPU is not busy; fall back to CPU when
  no GPU is available or the GPU is above the configured busy threshold.
- Keep exact-GP cross-validation exact by default. The shared-noise
  block-delete shortcut is allowed only as an explicit option because it reuses
  full-data hyperparameters and is therefore approximate. With
  `gp_exact_fast_cv = TRUE`, the backend must self-check against a real refit
  and fall back on mismatch. Benchmark-only force modes must report
  `cv_fast_path_approximate = TRUE` and `cv_fast_path_exact_refit = FALSE`.

## Current Reference Evidence

The current validation evidence is built from:

- G2F public maize data and the G2F 2025 competition-style 2024 prediction
  exercise: https://www.genomes2fields.org/resources/
- SoyNAM soybean data from the CRAN `SoyNAM` package:
  https://search.r-project.org/CRAN/refmans/SoyNAM/html/SoyNAM-package.html
- SFSI wheat HTP data from the CRAN source archive:
  https://cran.r-project.org/src/contrib/Archive/SFSI/
- Zenodo 3GS wheat processed data:
  https://doi.org/10.5281/zenodo.19217089

The Dryad wheat DOI `https://doi.org/10.5061/dryad.vx0k6dk3p` is useful as an
external source, but the direct download path can be browser-gated. When the
real archive is available locally, it can be added to the same report/check
workflow.

## Validation Workflow

Run the benchmark scripts that are feasible for the local machine and data
availability, then consolidate the outputs:

```r
source("tools/gp_generalization_validation_report.R")
```

This writes:

- `tools/tmp_gp_generalization_report/gp_generalization_validation_long.csv`
- `tools/tmp_gp_generalization_report/gp_generalization_validation_comparison.csv`
- `tools/tmp_gp_generalization_report/gp_generalization_validation_dataset_summary.csv`

Then run the lightweight regression gate:

```r
source("tools/check_gp_generalization_regression.R")
```

The gate intentionally uses broad rules:

- G2F reference checks must beat simple environment/location baselines.
- External probes may be mixed, but GP must not collapse relative to simple
  baselines or paired GBLUP.
- Across all available rows, GP must beat the simple baseline and GBLUP often
  enough to show broad value.
- Reported GP SE and PEV must be positive and internally consistent.

Thresholds are configurable through environment variables:

- `PREDICTPROR_GP_MAX_BASELINE_RMSE_LOSS_PCT`, default `30`
- `PREDICTPROR_GP_MAX_GBLUP_RMSE_LOSS_PCT`, default `20`
- `PREDICTPROR_GP_MIN_BASELINE_WIN_RATE`, default `0.50`
- `PREDICTPROR_GP_MIN_GBLUP_WIN_RATE`, default `0.50`
- `PREDICTPROR_GP_MAX_SE_PEV_LOG_RATIO`, default `log(2)`
- `PREDICTPROR_GP_REQUIRE_REFERENCE_ROWS`, default `1`

If a proposed default improves one benchmark but fails this gate, keep the old
default and add the new behavior as an explicit user-controlled option.

For exact-GP CV speed work, run the synthetic fast-path benchmark before
changing defaults:

```bash
python tools/benchmark_gp_exact_fast_cv.py --n 80 --n-markers 120 --folds 5 --iters 3
python tools/benchmark_gp_exact_fast_cv.py --n 36 --n-markers 80 --folds 4 --iters 3 --output-level predict_with_se
```

The benchmark compares cached exact refits, self-checked fast CV, and forced
fast CV. A forced speedup alone is not enough to change defaults; the
self-checked route must pass or fall back cleanly.

Then repeat the same decision check on real single-environment slices from G2F
and SoyNAM:

```bash
Rscript tools/benchmark_gp_exact_fast_cv_real_data.R
```

Useful local overrides:

```bash
PREDICTPROR_G2F2025_SINGLE_ENV=NEH1_2019 \
PREDICTPROR_SOYNAM_SINGLE_ENV=IA_2012 \
PREDICTPROR_GP_EXACT_FAST_CV_REAL_N_GID=72 \
PREDICTPROR_GP_EXACT_FAST_CV_REAL_MARKERS=160 \
Rscript tools/benchmark_gp_exact_fast_cv_real_data.R
```

For uncertainty checks, set
`PREDICTPROR_GP_EXACT_FAST_CV_REAL_RETURN_SE=true`. Forced fast-CV changes in
SE/variance are evidence to keep it explicit rather than default.
