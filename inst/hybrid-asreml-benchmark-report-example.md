# Hybrid ASReml-R Benchmark Report Example

This note mirrors the hybrid Bayesian benchmark-report example, but uses the
explicit `ASReml-R` hybrid path.

Use this when you want:

- model-native `Female_GCA`
- model-native `Male_GCA`
- model-native `SCA_effect`

and then want to convert the resulting run objects into one ranked hybrid
benchmark report.

## Workflow

1. run a masked true-prediction `hybrid_asreml` model
2. run one hybrid CV scenario
3. collect the run objects
4. build a benchmark table with:
   - `PredictProR::hybrid_benchmark_table_from_runs()`
5. export ranked summaries and plots with:
   - `PredictProR::hybrid_benchmark_report()`

## Example structure

```r
truth_lookup <- data.frame(
  HybridID = c("P1_x_P5", "P2_x_P4"),
  Observed_truth = c(0.20, -0.10),
  stringsAsFactors = FALSE
)

benchmark_tab <- PredictProR::hybrid_benchmark_table_from_runs(
  run_entries = list(
    asreml = list(
      true_prediction = asreml_true,
      cv = asreml_cv,
      cv_split_name = "One_New_Parent"
    )
  ),
  truth_lookup = truth_lookup,
  hybrid_id_col = "HybridID",
  truth_col = "Observed_truth"
)

report <- PredictProR::hybrid_benchmark_report(
  summary_table = benchmark_tab,
  masked_split_name = "masked_true_prediction",
  cv_split_name = "One_New_Parent",
  output_dir = "hybrid_asreml_report_outputs",
  export = TRUE
)
```

## Runnable example script

The repository now includes a runnable example script:

- [tools/hybrid_asreml_benchmark_report_example.R](../tools/hybrid_asreml_benchmark_report_example.R)

It:

- creates a small synthetic parent-resolved hybrid panel
- runs:
  - true prediction
  - `Hybrid_One_New_Parent` CV
- builds a benchmark table
- exports the ranked report artifacts under:
  - `tools/tmp_hybrid_asreml_report_example`

## Expected exported files

- `hybrid_shortlist_summary.csv`
- `hybrid_shortlist_ranked_masked_true_prediction.csv`
- `hybrid_shortlist_ranked_cv.csv`
- `hybrid_shortlist_recommendations.csv`
- `hybrid_shortlist_rmse_by_model.pdf`
- `hybrid_shortlist_mae_by_model.pdf`

## Important note

This example requires:

- `asreml` installed
- a valid `ASReml-R` license in the active R runtime

If either is unavailable, the example cannot run end to end.
