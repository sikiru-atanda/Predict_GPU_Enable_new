# Hybrid benchmark process G2F 2025 example

This example shows the preferred package-side hybrid comparison workflow on a
small real local G2F hybrid slice.

It uses:

- installed `PredictProR`
- direct `PredictProR::model_execute()` calls
- `PredictProR::hybrid_benchmark_process()`

Current example scope:

- `hybrid_ml`
- `Ridge_Regression`
- `SupportVectorMachine`
- `RandomForest`
- `Hybrid_One_New_Parent`

This local installed-package example is currently **CV-centered**. On the
hybrid-level G2F path used here, the stable package-side comparison output is
the `Hybrid_One_New_Parent` ranking. The masked true-prediction ranking file is
still written, but it can be empty on this path when explicit `Test` rows are
not materialized in the returned hybrid ML prediction table.

It assumes these local files already exist:

- [tools/tmp_g2f_2025_hybrid_pheno_only/hybrid_pheno_only_bundle.rds](../tools/tmp_g2f_2025_hybrid_pheno_only/hybrid_pheno_only_bundle.rds)
- `<g2f_2025 data dir>/geno_numerical.txt`

Run it with:

```powershell
$env:PREDICTPRO_HYBRID_MAX_N="30"
$env:PREDICTPRO_HYBRID_CV_METHOD="Hybrid_One_New_Parent"
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\hybrid_benchmark_process_g2f2025_example.R
```

The example writes:

- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_summary.csv`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_ranked_masked_true_prediction.csv`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_ranked_cv.csv`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_recommendations.csv`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_rmse_by_model.pdf`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_shortlist_mae_by_model.pdf`
- `tools/tmp_hybrid_benchmark_process_g2f2025/hybrid_benchmark_process_g2f2025_example.rds`

This is the preferred package-side comparison pattern when you already have a
small hybrid benchmark set and want:

- one standardized benchmark table
- ranked masked-true-prediction output
- ranked hybrid CV output
- recommendation table
- RMSE/MAE comparison plots
