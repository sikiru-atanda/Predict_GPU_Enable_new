# Hybrid focused workflow

PredictProR's current user-facing hybrid prediction focus is:

- `hybrid_asreml`
- `hybrid_ml` with:
  - `Ridge_Regression`
  - `SupportVectorMachine`
  - `RandomForest`

The next-stage hybrid kernel additions are now also available:

- `hybrid_bayes` with:
  - `GBLUP_BRR`
  - `RKHS`

## Current recommendation

The current practical recommendation is to keep the hybrid user-facing default
set narrow:

- `hybrid_asreml`
- `hybrid_ml` with:
  - `Ridge_Regression`
  - `SupportVectorMachine`
  - `RandomForest`

That recommendation is no longer based on a single local slice only. It now
has package-side breadth support from:

- local G2F `MNH1_2018`
- local G2F `MNH1_2019`
- public canola `CX`
- public canola `HZ`

Observed pattern across those breadth checks:

- `RandomForest` remains the strongest masked true-prediction baseline in most
  evaluated settings
- `Ridge_Regression` remains one of the strongest and most stable parent
  generalization baselines
- `SupportVectorMachine` remains competitive and, in canola `CX`, was the best
  model on both masked true prediction and `One_New_Parent`

So the practical message is:

- keep all three ML models available and promoted for hybrid work
- do not force one universal winner across all crops and environments
- use the package-side benchmark process to decide which of the three should be
  favored for a given dataset

Before running any hybrid path, use the package-side input contract if your raw
files are not already standardized:

- `PredictProR::hybrid_data_standard()`
- `PredictProR::validate_hybrid_input_standard()`
- [hybrid-data-standard.md](hybrid-data-standard.md)

## Current interpretation

- `hybrid_asreml`
  - use when you need explicit `Female_GCA + Male_GCA + SCA`
  - requires parent-resolved genotype or relationship matrices
  - now supports repeated hybrid rows across environments when
    `heter_groups` is supplied
- `hybrid_ml`
  - use when you have hybrid-level genotype rows
  - current promoted models are the focused shortlist above
  - now supports multi-environment hybrid rows through `heter_groups`,
    with environment features appended to the hybrid design matrix
- `hybrid_bayes`
  - use when you want parent-resolved hybrid kernels without the `ASReml-R`
    dependency
  - uses the same female GCA, male GCA, and SCA hybrid definition as
    `hybrid_asreml`
  - now supports multi-environment hybrid rows when `heter_groups` is supplied
  - supports:
    - true prediction
    - `Hybrid_Known_Parents`
    - `Hybrid_One_New_Parent`
    - `Hybrid_Both_New_Parents`

These Bayesian hybrid kernels are available now, but they are not yet the
default promoted shortlist in the local G2F benchmark runner.

## Local G2F workflow

If you already have the local G2F 2025 files used in this repository:

1. build the phenotype-only bundle

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\prepare_hybrid_g2f2025_local.R
```

2. run the focused hybrid shortlist benchmark

```powershell
$env:PREDICTPRO_HYBRID_MAX_N="40"
$env:PREDICTPRO_HYBRID_CV_METHOD="Hybrid_One_New_Parent"
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\benchmark_hybrid_shortlist_g2f2025.R
```

Default promoted models in that script are:

- `Ridge_Regression`
- `SupportVectorMachine`
- `RandomForest`

Hybrid Bayesian kernels are currently separate from that shortlist runner:

- `hybrid_bayes` with `GBLUP_BRR`
- `hybrid_bayes` with `RKHS`

## Key outputs

The shortlist benchmark writes:

- `hybrid_shortlist_summary.csv`
  - raw summary for the promoted models
- `hybrid_shortlist_ranked_masked_true_prediction.csv`
  - ranking by masked true prediction
- `hybrid_shortlist_ranked_cv.csv`
  - ranking by the chosen hybrid CV scenario
- `hybrid_shortlist_recommendations.csv`
  - one-line current best model by split
- `hybrid_shortlist_rmse_by_model.pdf`
  - side-by-side RMSE comparison across the evaluated split types
- `hybrid_shortlist_mae_by_model.pdf`
  - side-by-side MAE comparison across the evaluated split types

If you already have hybrid run objects in memory, the preferred package-side
comparison path is now:

- `PredictProR::hybrid_benchmark_process()`

That helper:

- builds the benchmark table from hybrid run objects
- ranks masked true prediction and hybrid CV scenarios
- returns comparison plots
- can export the same CSV/PDF outputs directly

If you already have a model-comparison table in memory, use the lower-level:

- `PredictProR::hybrid_benchmark_report()`

For a real installed-package example on the local G2F hybrid data, see:

- [hybrid-benchmark-process-g2f2025-example.md](hybrid-benchmark-process-g2f2025-example.md)
- [tools/hybrid_benchmark_process_g2f2025_example.R](../tools/hybrid_benchmark_process_g2f2025_example.R)

For the broader installed-package real-data validation used to support the
current shortlist recommendation, see:

- [tools/hybrid_benchmark_process_g2f2025_multi_env.R](../tools/hybrid_benchmark_process_g2f2025_multi_env.R)
- [tools/hybrid_benchmark_process_canola_public_example.R](../tools/hybrid_benchmark_process_canola_public_example.R)
- [hybrid_benchmark_process_g2f2025_multi_env_summary.csv](../tools/tmp_hybrid_benchmark_process_g2f2025_multi_env/hybrid_benchmark_process_g2f2025_multi_env_summary.csv)
- [hybrid_benchmark_process_canola_public_summary.csv](../tools/tmp_hybrid_benchmark_process_canola_public/hybrid_benchmark_process_canola_public_summary.csv)

If you want an end-to-end package example, including `hybrid_bayes` with both
`GBLUP_BRR` and `RKHS`, see:

- [hybrid-benchmark-report-example.md](hybrid-benchmark-report-example.md)
- [tools/hybrid_benchmark_report_example.R](../tools/hybrid_benchmark_report_example.R)

For the matching `hybrid_asreml` package-side example, see:

- [hybrid-asreml-benchmark-report-example.md](hybrid-asreml-benchmark-report-example.md)
- [tools/hybrid_asreml_benchmark_report_example.R](../tools/hybrid_asreml_benchmark_report_example.R)

All of those are written under:

- `tools/tmp_hybrid_shortlist_g2f2025`

For direct hybrid prediction runs from `model_execute()`, the exported report
folder now also writes:

- `summary_statistics.csv`
  - hybrid mode, model, parent-count, and train/test prediction summary
- `hybrid_component_summary.csv`
  - mean and spread of the three hybrid contribution channels
- `hybrid_train_test_summary.csv`
  - train/test row counts and prediction-error summary
- `CV_hybrid_observed_vs_predicted.pdf`
  - hybrid CV observed-vs-predicted plot when a hybrid CV path is run

The exact three component columns depend on the path:

- `hybrid_asreml` and `hybrid_bayes`
  - `Female_GCA`
  - `Male_GCA`
  - `SCA_effect`
- `hybrid_ml`
  - `Female_additive_contribution`
  - `Male_additive_contribution`
  - `Hybrid_interaction_contribution`

## Current boundary

The local `g2f_2025/geno_numerical.txt` file is hybrid-level genotype.
That means:

- it is suitable for `hybrid_ml`
- it is not sufficient for parent-resolved `hybrid_asreml` or `hybrid_bayes`

For `hybrid_asreml` and `hybrid_bayes`, you still need parent genotype or
relationship matrices keyed by the female and male parent IDs.

## Multi-environment hybrid note

All current hybrid families now accept repeated hybrid rows across environments
through `heter_groups`, with these semantics:

- CV scenarios still hold out hybrid IDs, not single environment rows
- when a hybrid is held out, all of its environment rows are held out together
- `hybrid_ml` and `hybrid_dl` add environment features to the hybrid-level
  design matrix
- `hybrid_asreml` and `hybrid_bayes` keep the same GCA/SCA kernel structure and
  add environment-aware fixed effects
