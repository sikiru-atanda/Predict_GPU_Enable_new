# PredictProR scenario examples

Plain R scripts that simulate small data in the PredictProR input format and
run every workflow with full, commented `model_execute()` calls, grouped the
way you would write them for your own data:

```r
model_runtime <- system.time({
  res <- PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,
    ...
    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    ...
    # Candidate models
    GS_model = m,
    ...
  )
})
```

Each script runs one full call per workflow; where only the model changes, the
same call runs inside a loop over a model vector. Use them to learn the
package, as templates for your own data, or to check an installation.

## Run them

```r
library(PredictProR)
folder <- system.file("examples", "scenarios", package = "PredictProR")
list.files(folder)

# one script
source(file.path(folder, "03_multi_environment.R"))

# everything, with a PASS / SKIP / FAIL table at the end (about 20-30 minutes)
source(file.path(folder, "run_all_scenarios.R"))

# only some scripts
only <- c("01", "09")
source(file.path(folder, "run_all_scenarios.R"))
```

To adapt a script, copy it and replace the simulated `pheno` / `Geno.data`
with your own data:
`file.copy(file.path(folder, "01_single_env_true_prediction.R"), ".")`.

## What each script covers

| Script | Full `model_execute()` calls |
|---|---|
| `01_single_env_true_prediction.R` | One environment, one trait, predict untested lines - one call looped over every compatible model (see below). This call lists **every model setting** (tuning grids, XGBoost, random forest, LightGBM, Bayesian, GP and deep-learning controls) at the package defaults, so you can see what exists |
| `02_single_env_cross_validation.R` | One call looped over K-Folds, stratified, repeated and hold-out CV; one call comparing all installed one-environment models, then refitting the best |
| `03_multi_environment.R` | Unbalanced multi-environment trials: Bayesian (incl. environment-specific variances), GP, ML and DL; ASReml looped over us / corgh / corh / corv / fa1 / fa2 / rr1 / rr2; CV0 / CV1 / CV2 and repeated versions; ASReml CV |
| `04_multi_trait.R` | Separate fits and joint Bayesian / ASReml (us, corgh, diag) / GP / ML / DL; multi-trait CV; automatic model choice; multi-trait multi-environment (Bayesian, GP, ASReml fa1) |
| `05_hybrid.R` | Hybrids: Bayesian GCA + SCA, ASReml, GP, ML and DL; known-parent, one-new-parent and both-new-parent CV |
| `06_classification.R` | Binary, ordinal and multiclass traits with every compatible model; classification CV; multi-environment binary, multiclass and ordinal traits |
| `07_feature_selection.R` | Standalone scoring + top-k; selection inside CV with each scoring model (Ridge_Regression, BayesB, RandomForest); selection for true prediction |
| `08_multi_omics_kernels.R` | Markers + raw omics layers (Bayesian marker, ML and DL models); raw omics -> kernel, precomputed omics kernel (COP) and several kernels (kernel_list) with every kernel model |
| `09_genotype_files.R` | VCF, HapMap (file and table) and CSV/TXT input with QC and imputation, including Beagle; standalone recoding and Beagle functions |
| `10_parallel_and_output.R` | Each parallel backend; saving results to a folder |
| `11_input_standards.R` | Printing the expected input format; validating data; messages for common mistakes |
| `12_polyploid.R` | Polyploid crops: tetraploid with every compatible model, tetraploid VCF / HapMap / CSV input, CV, multi-environment, multi-trait, hexaploid; the messages for common polyploid mistakes |
| `run_all_scenarios.R` | Runs all scripts and prints a PASS / SKIP / FAIL table |
| `00_setup.R` | Data simulators, optional-software checks and `log_result()` (used by every script) |

## Which models fit which scenario

Only models that fit and work well in a scenario are used in its script.

| Scenario | Models |
|---|---|
| One environment, one trait (Gaussian) | Bayesian: BayesA, BayesB, BayesC, BL, BRR, GBLUP_BRR, RKHS; ASReml: GBLUP; ML: Lasso, Ridge_Regression, PartialLeastSquare, SupportVectorMachine, K-NearestNeighbors, RandomForest, Xgboost, CatBoost, LightGBM; DL: DenseNeuralNet, TabAttention, TabNet, LightTreeNet, FactorNet, CrossNet, MixtureOfExperts, GPNet, DenseAttentionNet, Conv1DNet, ResNet (NeuralAdditive works but is slow; TabTransformer works but is unstable and less accurate than GBLUP on very small panels, so the examples leave it out there); GP: Kernel-GBLUP, Gaussian-Process-GBLUP |
| Multi-environment | Bayesian: GBLUP_BRR, RKHS; ASReml: GBLUP (`var_cov_str` us, corgh, corh, corv, faK, rrK); GP: Kernel-GBLUP, Gaussian-Process-GBLUP, FA-GBLUP; ML (`met_ml_dl = TRUE`): RandomForest, Xgboost, CatBoost, LightGBM; DL: DenseNeuralNet, TabTransformer, TabAttention, TabNet, MixtureOfExperts |
| Multi-trait | Bayesian: GBLUP_BRR, RKHS; ASReml: GBLUP (us, corgh, diag); GP: Gaussian-Process-GBLUP, FA-GBLUP, Scalable-GBLUP; ML: RandomForest, Ridge_Regression, PartialLeastSquare, Lasso, SupportVectorMachine, Xgboost, CatBoost; DL: DenseNeuralNet, TabTransformer |
| Multi-trait, multi-environment | GBLUP_BRR, RKHS, GBLUP (ASReml), Gaussian-Process-GBLUP, FA-GBLUP, Scalable-GBLUP |
| Hybrid | Bayesian: GBLUP_BRR, RKHS; ASReml: GBLUP; GP: Gaussian-Process-GBLUP; ML (hybrid genotypes): RandomForest, Ridge_Regression, PartialLeastSquare, Lasso, SupportVectorMachine, Xgboost, CatBoost; DL: DenseNeuralNet, TabTransformer |
| Binary trait | the 7 Bayesian models; SupportVectorMachine, K-NearestNeighbors, RandomForest, Xgboost, CatBoost, LightGBM; all deep-learning models (NeuralAdditive works but is left out of the example: slow) |
| Ordinal trait | the 7 Bayesian models |
| Multiclass trait | SupportVectorMachine, K-NearestNeighbors, RandomForest, Xgboost, CatBoost, LightGBM; all deep-learning models (NeuralAdditive works but is left out of the example: slow) |
| Multi-environment classification | binary: RKHS + RandomForest, Xgboost, CatBoost, LightGBM and the 5 MET DL models; multiclass: the 9 MET ML/DL models; ordinal: RKHS |
| Multi-omics | raw layers: Bayesian marker models, ML, DL; omics kernels: GBLUP_BRR, RKHS, GBLUP (ASReml), Kernel-GBLUP, Gaussian-Process-GBLUP; several kernels from the same markers (kernel_list): GBLUP_BRR, RKHS, Kernel-GBLUP, Gaussian-Process-GBLUP |
| Polyploid (`ploidy = 4L`, `gmatrix_method = "VanRaden"`, `het_threshold = NULL`) | one environment: the same models as diploid except GPNet (not supported for polyploids: stops with a message; NeuralAdditive works but is slow); multi-environment: GBLUP_BRR, RKHS, Kernel-GBLUP, Gaussian-Process-GBLUP, FA-GBLUP, ASReml GBLUP, the 9 MET ML/DL models; multi-trait: Bayesian, ASReml, GP, ML, DL. Not for polyploids: GPNet, the Yang and dominance relationship matrices, Beagle imputation (each stops with a message) |

ASReml GBLUP, Lasso, Ridge_Regression, PartialLeastSquare and the Gaussian
process models are for Gaussian (continuous) traits only.
## Input format in brief

* **Phenotypes** (`pheno_data`): a data frame with a line ID column (`gen_name`)
  and trait column(s). Lines to predict have `NA`. Multi-environment data are in
  long format (one row per line x environment) with an environment column given
  as `heter_groups`. Hybrids have `HybridID`, `Female`, `Male` columns.
* **Genotypes** (`geno_data`): numeric matrix, lines in rows (row names = IDs),
  markers in columns, coded 0/1/2. Or a VCF / HapMap / CSV file. Or a
  precomputed relationship matrix (`gmatrix`, square, row and column names = IDs).
  Every genotyped line is predicted, also lines without a phenotype row.
* `prediction_input_standard()` prints the full description.

## Optional software

Calls that need optional software are **skipped** (not failed) when it is
missing: **ASReml-R** (licensed) for `engine = "asreml"`; **Python** for classical
ML (`PREDICTPRO_PYTHON`), Gaussian process (`PREDICTPRO_GP_PYTHON`) and deep
learning (`PREDICTPRO_DL_PYTHON`) - `setup_predictgp_env()` and
`setup_predictdl_env()` can create them; **Java 8+** for Beagle.

## Notes

* QC, imputation and LD pruning are switched on as in real use. The examples
  use short MCMC chains and few trees (marked in the scripts); deep-learning models use 100 epochs because the attention models need that many to learn so
  they finish quickly; use the suggested values for real analyses.
* Results stay in memory (`system_database = TRUE`). Set it to `FALSE` to also
  write a results folder; choose its location with
  `options(PredictProR.output_dir = "...")`.
