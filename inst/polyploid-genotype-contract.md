# Polyploid genotype contract

PredictProR uses one uniform crop ploidy per analysis. Genotypes are represented
internally as ALT-allele dosage from 0 through `ploidy`; fractional values are
allowed after mean or KNN imputation. This contract applies before feature
selection, kernel construction, relationship-matrix construction, and model
fitting.

## Supplying ploidy

- Raw biallelic VCF files can infer ploidy from the number of allele fields in
  `GT`. Mixed-ploidy calls fail loudly.
- Numeric matrices with observed values above two must set `ploidy` explicitly
  or carry `attr(x, "ploidy")`. The maximum observed dosage is not treated as a
  ploidy estimator because some dosage classes may be absent.
- HapMap polyploid heterozygotes must give every allele copy. At ploidy four,
  `AAAA`, `AAAG`, `AAGG`, `AGGG`, and `GGGG` encode dosages 0 through 4 when A
  is the first declared allele and G is the second. An IUPAC heterozygote such
  as `R` is dosage-ambiguous above diploidy and fails by default.

```r
qc <- vcf_qc_recode("potato.vcf.gz", ploidy = 4, impute = TRUE)

x <- matrix(c(0, 1, 2, 3, 4, NA), nrow = 2)
imp <- impute_genotypes_native(x, ploidy = 4, method = "mean")
```

## QC, recoding, and native imputation

For marker `j`, PredictProR estimates ALT frequency as
`p_j = mean(dosage_j) / ploidy`, minor-allele frequency as
`min(p_j, 1 - p_j)`, and heterozygosity as the observed proportion satisfying
`0 < dosage < ploidy`. Call rates use the observed/non-missing indicator and do
not depend on ploidy.

Public recoding choices are:

- `alt_dosage`: `0, ..., ploidy`;
- `centered_dosage`: `dosage - ploidy / 2`;
- `allele_frequency`: `dosage / ploidy`.

Native mean, median, mode, and KNN imputation operate on dosage and enforce the
closed interval `[0, ploidy]`. They do not phase haplotypes and do not infer
untyped variants from a reference panel.

## Beagle boundary

The Beagle 5.4 integration is deliberately restricted to diploid GT calls.
PredictProR rejects `ploidy != 2` before running Beagle. Polyploid data should
use native dosage imputation unless the user performs phasing or reference-panel
imputation with a separately validated external polyploid program before import.

## Relationship matrices and models

The additive VanRaden matrix centers marker dosage by its sample mean and uses
the denominator `sum(ploidy * p * (1 - p))`. Weighted VanRaden and additive
epistasis follow the same ploidy-aware additive base. The currently implemented
Yang and dominance parameterizations are diploid-only and fail loudly for
polyploid input. Generic linear, Gaussian, polynomial, Matern, Laplacian, and
rational-quadratic kernels accept numeric polyploid dosage and retain the
`ploidy` attribute.

`model_execute(ploidy = ...)` propagates the resolved value through genotype
QC, trait-specific feature selection, rebuilt marker kernels/GRMs, ordinary and
multi-trait workflows, multi-environment workflows, multi-omics workflows, and
supported hybrid prediction routes.

## Development-time public reference benchmark

Public polyploid packages are validation references, not PredictProR runtime
dependencies. Maintainers can run:

```r
Rscript tools/benchmark_polyploid_reference.R
```

The benchmark checks P2, P4, and P6 VanRaden matrices against AGHmatrix; P2,
P4, and P6 matrices against polyBreedR's VR1 implementation; a real tetraploid
potato matrix; missing-value mean imputation; and P2, P3, P4, and P6 VCF GT
dosage conversion. The benchmark must report all comparisons within its stated
numeric tolerance before a release claim is made.
