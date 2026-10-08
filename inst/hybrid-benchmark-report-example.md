# Hybrid Benchmark Report Example

This note shows the intended package-side workflow for turning hybrid
`model_execute()` outputs into one ranked hybrid benchmark report.

The same pattern works for:

- `hybrid_asreml`
- `hybrid_bayes` with:
  - `GBLUP_BRR`
  - `RKHS`
- `hybrid_ml`

## Core idea

1. run one or more hybrid models
2. keep the returned run objects
3. build a standardized benchmark table with:
   - `PredictProR::hybrid_benchmark_table_from_runs()`
4. rank and export it with:
   - `PredictProR::hybrid_benchmark_report()`

## Example layout

Assume you have:

- one masked true-prediction run
- one hybrid CV run for the same model
- a truth lookup table for the masked hybrids

```r
truth_lookup <- data.frame(
  HybridID = c("H2", "H3"),
  Observed_truth = c(9.4, 7.6),
  stringsAsFactors = FALSE
)

benchmark_tab <- PredictProR::hybrid_benchmark_table_from_runs(
  run_entries = list(
    rf = list(
      true_prediction = rf_true,
      cv = rf_cv,
      cv_split_name = "Hybrid_One_New_Parent"
    ),
    brr = list(
      true_prediction = brr_true,
      cv = brr_cv,
      cv_split_name = "Hybrid_One_New_Parent"
    ),
    rkhs = list(
      true_prediction = rkhs_true,
      cv = rkhs_cv,
      cv_split_name = "Hybrid_One_New_Parent"
    )
  ),
  truth_lookup = truth_lookup,
  hybrid_id_col = "HybridID",
  truth_col = "Observed_truth"
)

hybrid_report <- PredictProR::hybrid_benchmark_report(
  summary_table = benchmark_tab,
  masked_split_name = "masked_true_prediction",
  cv_split_name = "Hybrid_One_New_Parent",
  output_dir = "hybrid_report_outputs",
  export = TRUE
)
```

The exported folder will contain:

- `hybrid_shortlist_summary.csv`
- `hybrid_shortlist_ranked_masked_true_prediction.csv`
- `hybrid_shortlist_ranked_cv.csv`
- `hybrid_shortlist_recommendations.csv`
- `hybrid_shortlist_rmse_by_model.pdf`
- `hybrid_shortlist_mae_by_model.pdf`

## Hybrid Bayesian examples

### `GBLUP_BRR`

```r
brr_true <- PredictProR::model_execute(
  pheno_data = hybrid_pheno,
  geno_data = parent_geno,
  gmatrix_method = "VanRaden",
  response = "Yield",
  gen_name = "HybridID",
  female_parent = "Female",
  male_parent = "Male",
  GS_model = "GBLUP_BRR",
  response_family = "gaussian",
  hybrid_bayes = TRUE,
  system_database = TRUE,
  message = FALSE,
  nIter = 1200,
  burnIn = 400,
  thin = 5
)

brr_cv <- PredictProR::model_execute(
  pheno_data = subset(hybrid_pheno, !is.na(Yield)),
  geno_data = parent_geno,
  gmatrix_method = "VanRaden",
  response = "Yield",
  gen_name = "HybridID",
  female_parent = "Female",
  male_parent = "Male",
  GS_model_cv = "GBLUP_BRR",
  response_family = "gaussian",
  hybrid_bayes = TRUE,
  cross_validation = TRUE,
  cv_evaluation_only = TRUE,
  cross_validation_meth = "Hybrid_One_New_Parent",
  eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
  system_database = TRUE,
  message = FALSE,
  nIter = 1200,
  burnIn = 400,
  thin = 5
)
```

### `RKHS`

```r
rkhs_true <- PredictProR::model_execute(
  pheno_data = hybrid_pheno,
  geno_data = parent_geno,
  gmatrix_method = "VanRaden",
  response = "Yield",
  gen_name = "HybridID",
  female_parent = "Female",
  male_parent = "Male",
  GS_model = "RKHS",
  response_family = "gaussian",
  hybrid_bayes = TRUE,
  system_database = TRUE,
  message = FALSE,
  nIter = 1200,
  burnIn = 400,
  thin = 5
)

rkhs_cv <- PredictProR::model_execute(
  pheno_data = subset(hybrid_pheno, !is.na(Yield)),
  geno_data = parent_geno,
  gmatrix_method = "VanRaden",
  response = "Yield",
  gen_name = "HybridID",
  female_parent = "Female",
  male_parent = "Male",
  GS_model_cv = "RKHS",
  response_family = "gaussian",
  hybrid_bayes = TRUE,
  cross_validation = TRUE,
  cv_evaluation_only = TRUE,
  cross_validation_meth = "Hybrid_One_New_Parent",
  eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
  system_database = TRUE,
  message = FALSE,
  nIter = 1200,
  burnIn = 400,
  thin = 5
)
```

## Notes

- `hybrid_bayes` and `hybrid_asreml` need parent-resolved genotype or
  relationship information.
- `hybrid_ml` can work with hybrid-level genotype rows.
- the report helper does not care which hybrid family produced the runs, as
  long as the run objects follow the standard PredictProR hybrid output shape.
