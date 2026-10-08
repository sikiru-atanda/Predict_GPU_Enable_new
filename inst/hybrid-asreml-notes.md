# Hybrid ASReml-R notes

PredictProR now has a first-pass hybrid Gaussian `ASReml-R` path for:

- female GCA
- male GCA
- SCA

Current public path:

- `hybrid_asreml = TRUE`
- `GS_model = "GBLUP"`
- `engine = "asreml"`
- `response_family = "gaussian"`
- single- or multi-environment hybrid rows via `heter_groups`

Package-side hybrid input validation is now available before model execution:

- `PredictProR::hybrid_data_standard()`
- `PredictProR::validate_hybrid_input_standard()`
- [hybrid-data-standard.md](hybrid-data-standard.md)

Required phenotype columns:

- hybrid ID
- female parent ID
- male parent ID
- response
- and `heter_groups` when hybrids are repeated across environments

Current supported parent-data patterns:

1. shared parent matrix
- `geno_data` plus `gmatrix_method`
- or one shared `gmatrix`

2. separate female and male parent matrices
- `female_geno_data`
- `male_geno_data`
- plus `gmatrix_method`

3. separate female and male relationship matrices
- `female_gmatrix`
- `male_gmatrix`

Returned hybrid components:

- `Predicted_value`
- `Female_GCA`
- `Male_GCA`
- `SCA_effect`

Current multi-environment behavior:

- repeated hybrid rows are allowed when `heter_groups` is supplied
- the same hybrid GCA/SCA kernels are reused across environments
- environment effects are handled on the fixed side
- hybrid CV scenarios still hold out hybrid IDs across all their environment
  rows together

User-facing hybrid report exports now also include:

- `summary_statistics.csv`
- `hybrid_component_summary.csv`
- `hybrid_train_test_summary.csv`
- `CV_hybrid_observed_vs_predicted.pdf` for hybrid CV runs

## Current hybrid shortlist

PredictProR does not need to promote every ML/DL model for hybrid prediction.
The current recommended shortlist is:

- `hybrid_asreml`
  - `GS_model = "GBLUP"`
  - `engine = "asreml"`
  - strongest option when explicit `Female_GCA + Male_GCA + SCA` decomposition
    is required
- `hybrid_ml` with:
  - `Ridge_Regression`
  - `SupportVectorMachine`
  - `RandomForest`

These are the hybrid models that should be surfaced first in examples and
benchmarks. Other hybrid ML/DL models can remain technically available for
experimentation, but they are not currently the promoted default set.

### Why this shortlist

- `ASReml-R` is the only path with explicit mixed-model `GCA/SCA`
  interpretation.
- `Ridge_Regression` is a strong and fast lean baseline.
- `SupportVectorMachine` has shown robust hybrid prediction behavior and
  remains competitive on parent-generalization CV.
- `RandomForest` is the strongest current masked direct hybrid prediction
  baseline on the local G2F benchmark slice.

Current local G2F evidence on the `MNH1_2018` 40-hybrid slice:

- masked true prediction:
  - `RandomForest` was best
  - `SupportVectorMachine` was second
  - `Ridge_Regression` followed
- `Hybrid_One_New_Parent`:
  - `Ridge_Regression` was best
  - `SupportVectorMachine` was second
  - `RandomForest` followed

That is the practical reason these three ML models are the current promoted
hybrid set.

### Breadth validation note

That shortlist is now supported by broader package-side real-data validation,
not only the original local `MNH1_2018` slice.

Installed-package comparison runs were completed on:

- local G2F `MNH1_2018`
- local G2F `MNH1_2019`
- public canola `CX`
- public canola `HZ`

Practical takeaway from those runs:

- `RandomForest` is the strongest repeated masked true-prediction performer
- `Ridge_Regression` remains one of the strongest and most stable
  `One_New_Parent` baselines
- `SupportVectorMachine` remains robust, and on canola `CX` it was best for
  both masked true prediction and `One_New_Parent`

So PredictProR should not claim one universal hybrid ML winner. The defensible
position is:

- promote `Ridge_Regression`, `SupportVectorMachine`, and `RandomForest`
- use package-side benchmarking to decide which of the three is preferred for a
  given crop, environment set, and validation objective

Reference outputs:

- [hybrid_benchmark_process_g2f2025_multi_env_summary.csv](../tools/tmp_hybrid_benchmark_process_g2f2025_multi_env/hybrid_benchmark_process_g2f2025_multi_env_summary.csv)
- [hybrid_benchmark_process_canola_public_summary.csv](../tools/tmp_hybrid_benchmark_process_canola_public/hybrid_benchmark_process_canola_public_summary.csv)

### Hybrid Bayesian kernel additions

The next-stage hybrid kernel methods are now available through:

- `hybrid_bayes = TRUE`
- `GS_model = "GBLUP_BRR"` or `GS_model = "RKHS"`

Current scope:

- gaussian only
- explicit female GCA, male GCA, and SCA hybrid kernel construction
- true prediction
- hybrid CV scenarios:
  - `Hybrid_Known_Parents`
  - `Hybrid_One_New_Parent`
  - `Hybrid_Both_New_Parents`

These Bayesian hybrid kernels are now implemented, but they are not yet the
default promoted shortlist in the local G2F benchmark runner. The current
default-facing shortlist remains:

- `Ridge_Regression`
- `SupportVectorMachine`
- `RandomForest`

### Current local G2F benchmarking scripts

Using the local hybrid-level G2F 2025 bundle now available in the workspace:

- [tools/benchmark_hybrid_ml_g2f2025.R](../tools/benchmark_hybrid_ml_g2f2025.R)
- [tools/benchmark_hybrid_dl_g2f2025.R](../tools/benchmark_hybrid_dl_g2f2025.R)
- [tools/compare_hybrid_ml_dl_g2f2025.R](../tools/compare_hybrid_ml_dl_g2f2025.R)
- [tools/benchmark_hybrid_shortlist_g2f2025.R](../tools/benchmark_hybrid_shortlist_g2f2025.R)
- [inst/hybrid-focused-workflow.md](hybrid-focused-workflow.md)
- [inst/hybrid-benchmark-report-example.md](hybrid-benchmark-report-example.md)
- [tools/hybrid_benchmark_report_example.R](../tools/hybrid_benchmark_report_example.R)
- [inst/hybrid-asreml-benchmark-report-example.md](hybrid-asreml-benchmark-report-example.md)
- [tools/hybrid_asreml_benchmark_report_example.R](../tools/hybrid_asreml_benchmark_report_example.R)

Important boundary:

- the local `g2f_2025/geno_numerical.txt` file is hybrid-level genotype, so it
  is suitable for `hybrid_ml` and `hybrid_dl`
- it is not sufficient for parent-resolved `hybrid_asreml`
- `hybrid_asreml` still needs parent genotype or parent relationship matrices
  keyed by the female and male parent IDs

For the current promoted shortlist, the default shortlist runner is focused on:

- `Ridge_Regression`
- `SupportVectorMachine`
- `RandomForest`

with no DL models promoted by default in that script.

## Real public data candidates

### 1. Maize Genomes to Fields (G2F) - best first benchmark

This is the strongest immediate public candidate for the first real PredictProR
hybrid benchmark.

Why it fits:

- public hybrid phenotypes
- public inbred parent genotypes
- multiple years and environments
- established maize hybrid prediction use case
- explicit documentation that the released inbred genotypes correspond to the
  parents used to produce the featured hybrids

Good entry points:

- 2014-2015 release note:
  - https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-018-3508-1
- 2014-2017 data release note:
  - https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-020-4922-8
- 2018-2019 data release note:
  - https://bmcgenomdata.biomedcentral.com/articles/10.1186/s12863-023-01129-2
- CyVerse 2018 dataset:
  - https://doi.org/10.25739/anqq-sg86
- CyVerse 2019 dataset:
  - https://doi.org/10.25739/t651-yy97

Key practical evidence from the release notes:

- 2014-2015 note states the project released DNA sequences of inbreds,
  including the inbreds used to produce featured hybrids, plus phenotype
  measurements for inbreds and hybrids.
- 2018-2019 note states 1153 publicly available hybrids were evaluated and the
  inbred parents of the tested hybrids were genotyped using the Practical
  Haplotype Graph (PHG).

Practical note:

- G2F releases hybrid phenotypes and inbred parental genotypes.
- The main preparation step is reconstructing a clean hybrid table with:
  - hybrid ID
  - female parent
  - male parent
  - response trait
- This is the best source for validating the dedicated `Female_GCA + Male_GCA +
  SCA` implementation because the public notes explicitly connect hybrids to
  their inbred parents.

### 2. MaizeGEP - strong secondary benchmark

This looks promising as a modern hybrid-prediction benchmark because it was
presented specifically as a maize hybrids dataset with genotype, phenotype, and
envirotype for genomic selection.

Reference:

- https://pubmed.ncbi.nlm.nih.gov/41485082/
- https://academic.oup.com/gpb/advance-article/doi/10.1093/gpbjnl/qzaf140/8413486

Practical note:

- MaizeGEP contains 260 hybrid maize varieties, 12,233 selected tag SNPs,
  phenotypes for 11 traits across 2382 year-county locations, and daily
  meteorological records.
- It is excellent for multi-environment hybrid benchmarking and for comparing
  hybrid ML/DL against hybrid mixed models.
- However, from the public summary it is less immediately clear than G2F
  whether the released genotype is explicitly decomposable into separate female
  and male parent tables for GCA/SCA construction.
- So it is very useful, but still second to G2F for the first parent-resolved
  `ASReml-R` hybrid validation.

### 3. Sunflower hybrid panels - useful if parent tables can be recovered

Sunflower is a strong domain fit for PredictProR hybrid prediction, but public
dataset usability depends on whether the released materials include explicit
female and male parent identities or at least a recoverable crossing design.

Useful references:

- Genomic prediction of sunflower hybrid oil content:
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC5613134/
- Sunflower hybrid breeding review:
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC5776114/

Why it matters:

- the oil-content study used an incomplete factorial design with 36 CMS lines
  and 36 restorer lines to characterize 452 sunflower hybrids
- this is conceptually very close to the `Female_GCA + Male_GCA + SCA` path now
  implemented in PredictProR

Practical note:

- sunflower is highly relevant scientifically
- but the public benchmark is only immediately actionable if we can obtain:
  - hybrid phenotypes
  - female line IDs
  - male line IDs
  - parent genotypes or a parent relationship matrix

### 4. Canola hybrid panels - very strong practical benchmark

Canola is also a strong fit, especially because the published hybrid prediction
studies explicitly use parent-derived predictor information and large hybrid
trial panels.

Useful references:

- Testcross performance in canola:
  - https://pubmed.ncbi.nlm.nih.gov/26824924/
- Multi-omics hybrid prediction in canola:
  - https://pubmed.ncbi.nlm.nih.gov/33523261/

Why it matters:

- the 2021 canola study reports:
  - 950 F1 hybrids
  - 477 parental lines
  - 13,201 SNP markers
  - seven agronomic traits in multi-location field trials
- that is a very strong downstream benchmark once parent-resolved data access is
  available

Practical note:

- canola may become one of the best non-maize hybrid benchmarks for PredictProR
- but G2F still remains the clearest first public source for a parent-resolved
  end-to-end benchmark in the package today

## Recommended next benchmark path

For the first real hybrid benchmark in PredictProR:

1. use one G2F release year first
2. build a clean hybrid table with parent IDs
3. fit:
   - female GCA
   - male GCA
   - SCA
4. then define hybrid CV scenarios:
   - new hybrid from known parents
   - one new parent
   - both parents new

After that:

5. use MaizeGEP as a second benchmark for:
   - multi-environment hybrid prediction
   - hybrid ML/DL comparison
   - envirotype-aware hybrid prediction
6. if public parent-resolved files are available, add:
   - one sunflower benchmark
   - one canola benchmark
