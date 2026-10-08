test_that("genotype QC runs with a filter disabled via a NULL threshold", {
  # A NULL threshold switches that QC filter off (e.g. het_threshold = NULL
  # for heterozygous hybrid genotypes), but the QC summary table was built
  # with data.frame(het_threshold = NULL, ...) and failed with
  # "arguments imply differing number of rows".
  set.seed(5)
  geno <- matrix(sample(0:2, 40 * 30, replace = TRUE), 40, 30,
                 dimnames = list(sprintf("H%02d", 1:40), sprintf("m%02d", 1:30)))
  storage.mode(geno) <- "double"

  out <- PredictProR:::geno_precheck(
    object_geno = geno, qc_filtering = TRUE, het_threshold = NULL,
    maf_threshold = 0.01, snp_call_rate_threshold = 0.9, ind_call_rate_threshold = 0.9,
    impute = FALSE, ld_prunning_qc = FALSE, ploidy = 2L, message = FALSE
  )
  kept <- out$object_geno %||% out[[1]]
  # heterozygous markers are kept when the filter is off
  expect_equal(ncol(kept), ncol(geno))
})

test_that("genotype QC stops with a clear message when it removes every marker", {
  # Heterozygous material with the default inbred het_threshold used to fail
  # later with "Feature matrices must have non-missing predictor names".
  geno <- matrix(1, 20, 8, dimnames = list(sprintf("H%02d", 1:20), sprintf("m%02d", 1:8)))
  geno[cbind(1:20, (1:20 %% 8) + 1L)] <- 0
  geno[cbind(1:20, ((1:20 + 3L) %% 8) + 1L)] <- 2
  expect_error(
    PredictProR:::geno_precheck(
      object_geno = geno, qc_filtering = TRUE, het_threshold = 0.1,
      maf_threshold = 0.01, snp_call_rate_threshold = 0.9, ind_call_rate_threshold = 0.9,
      impute = FALSE, ld_prunning_qc = FALSE, ploidy = 2L, message = FALSE
    ),
    "removed every marker.*heterozygosity 8"
  )
})