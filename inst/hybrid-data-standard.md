# Hybrid data standard

PredictProR now validates hybrid prediction inputs before entering the model
paths. When the uploaded data is not in a supported hybrid format, the package
stops early and prints the required contract instead of proceeding into a later
model failure.

Package-side helpers:

- `PredictProR::hybrid_data_standard()`
- `PredictProR::validate_hybrid_input_standard()`

## Phenotype standard

Required columns:

- `HybridID`
- `Female`
- `Male`
- one selected response column

Optional columns:

- `Env`
- `Rep`
- `Year`
- `Location`

For repeated hybrid rows across environments:

- include the environment column in `pheno_data`
- pass that column name through `heter_groups`

## Parent-resolved genomic standard

Used by:

- `hybrid_asreml`
- `hybrid_bayes`

Accepted inputs:

- `female_geno_data` + `male_geno_data`
- `female_gmatrix` + `male_gmatrix`
- shared `geno_data` keyed by parent IDs
- shared `gmatrix` keyed by parent IDs

Requirement:

- row names must match the parent IDs appearing in `Female` and `Male`

## Hybrid-level genomic standard

Used by:

- `hybrid_ml`
- `hybrid_dl`

Accepted input:

- `geno_data` keyed by `HybridID`

Requirement:

- row names must match the `HybridID` values in `pheno_data`

## Hybrid CV scenarios

Supported standard scenarios:

- `Hybrid_Known_Parents`
- `Hybrid_One_New_Parent`
- `Hybrid_Both_New_Parents`
