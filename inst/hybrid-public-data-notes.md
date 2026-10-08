## Public Hybrid Data Downloads

The project now includes locally downloaded public hybrid data under:

- `external_data/hybrid_public/canola_plosgen_2021`
- `external_data/hybrid_public/sunflower_hesva0`

### Canola

Source:
- PLOS Genetics supplementary tables from the canola hybrid paper

Downloaded:
- `pgen.1009879.s008.xlsx`
- `pgen.1009879.s009.xlsx`
- `pgen.1009879.s010.xlsx`

What it gives:
- parental line information
- hybrid phenotype observations across environments
- explicit female and male parent mapping for each hybrid

What it does not yet give:
- a parent marker matrix ready for `hybrid_asreml` or `hybrid_bayes`

Prepared output:
- `tools/tmp_canola_plosgen2021_hybrid_bundle/canola_hybrid_bundle.rds`

Build command:
```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\prepare_canola_plosgen2021_hybrid_bundle.R
```

### Sunflower

Source:
- INRAE / recherche.data.gouv HESVA0 dataset

Downloaded:
- `15EX05_field.tab`
- `15EX05_genotype.tab`
- `15EX05_expression_noNA.tab`
- `15EX05_geneList.tab`
- `15EX05_geneList_173genes.xls`

What it gives:
- hybrid IDs
- hybrid-level molecular predictors
- field metadata

Current limitation:
- the downloaded field table does not include the target response trait needed for the current hybrid prediction benchmark path

Prepared output:
- `tools/tmp_sunflower_hesva0_hybrid_omics/sunflower_hybrid_omics_bundle.rds`

Build command:
```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\prepare_sunflower_hesva0_hybrid_omics_bundle.R
```

### Practical Status

Current direct benchmark candidate:
- canola phenotype bundle

Current staged future-use dataset:
- sunflower hybrid omics bundle

For parent-resolved `Female_GCA + Male_GCA + SCA` benchmarking, both crops still need parent genotype matrices aligned to the named female and male parents.

### Canola Parent Marker Recovery

The project now also includes:
- `external_data/hybrid_public/BnaSNPDB`

This is the processed public rapeseed SNP database used to test whether the canola paper's hybrid parents can be matched directly to a released accession panel.

Audit script:
```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\audit_canola_parent_marker_recovery.R
```

Audit outputs:
- `tools/tmp_canola_parent_marker_audit/canola_parent_direct_bnasnpdb_matches.csv`
- `tools/tmp_canola_parent_marker_audit/canola_parent_pedigree_token_hits.csv`
- `tools/tmp_canola_parent_marker_audit/canola_parent_marker_recovery_summary.csv`

Current finding:
- the 60 canola hybrid parents do not map directly by released name to the 1,007-accession BnaSNPDB panel
- so a full parent SNP matrix for `hybrid_asreml` or `hybrid_bayes` still requires either:
  - raw-read processing from the study BioProject `PRJNA664250`, or
  - another exact parent-level SNP release for those `Ogu-CMS` and `Ogu-R` lines
