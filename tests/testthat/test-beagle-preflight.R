write_test_vcf <- function(path, pos, gt) {
  # gt: markers x samples character matrix of GT calls
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(gt)), collapse = "\t")
  )
  body <- vapply(seq_along(pos), function(i) {
    paste(c("1", pos[i], paste0("m", i), "A", "G", ".", "PASS", ".", "GT", gt[i, ]), collapse = "\t")
  }, character(1))
  writeLines(c(header, body), path)
  path
}

test_that("GT ploidy extraction is vectorised and keeps the per-call contract", {
  expect_identical(
    PredictProR:::gp_extract_gt_ploidy(c("0/1", ".", "./.", "0|1|1|0", NA, "1", "")),
    c(2L, NA, 2L, 4L, NA, 1L, NA)
  )
  gt <- rep(c("0/0", "0/1", "./."), length.out = 3e5)
  expect_lt(system.time(PredictProR:::gp_resolve_gt_ploidy(gt, "auto"))[["elapsed"]], 5)
})

test_that("Beagle preflight rejects unsorted positions with a clear message", {
  gt <- matrix("0/1", 3, 2, dimnames = list(NULL, c("s1", "s2")))
  ok <- write_test_vcf(tempfile(fileext = ".vcf"), c(100, 200, 300), gt)
  bad <- write_test_vcf(tempfile(fileext = ".vcf"), c(100, 300, 200), gt)
  expect_true(PredictProR:::gp_beagle_check_sorted_positions(ok))
  expect_error(PredictProR:::gp_beagle_check_sorted_positions(bad), "not sorted within chromosome '1'.*300.*200")
})

test_that("call-rate QC before Beagle uses observed calls (markers first, then samples)", {
  gt <- matrix("0/1", 10, 4, dimnames = list(NULL, c("s1", "s2", "s3", "s4")))
  gt[1:6, "s4"] <- "./."       # s4 60% missing
  gt[2, c("s1", "s2")] <- "./." # marker 2 has 3/4 missing -> removed first
  vcf <- write_test_vcf(tempfile(fileext = ".vcf"), seq(100, 1000, by = 100), gt)
  out <- PredictProR:::gp_beagle_call_rate_exclusions(
    vcf, snp_call_rate_threshold = 0.7, ind_call_rate_threshold = 0.7,
    output_prefix = tempfile()
  )
  expect_equal(out$markers_removed, 1L)
  expect_identical(readLines(out$excludemarkers), "1:200")
  expect_equal(out$samples_removed, 1L)
  expect_identical(readLines(out$excludesamples), "s4")
  # thresholds < 0.5 are maximum missing rates, as in vcf_qc_recode()
  out2 <- PredictProR:::gp_beagle_call_rate_exclusions(
    vcf, snp_call_rate_threshold = 0.3, ind_call_rate_threshold = NULL, output_prefix = tempfile()
  )
  expect_equal(out2$markers_removed, 1L)
  expect_null(out2$excludesamples)
  none <- PredictProR:::gp_beagle_call_rate_exclusions(vcf, NULL, NULL, output_prefix = tempfile())
  expect_null(none$excludemarkers)
})

test_that("Beagle map-file validation rejects headers and accepts text chromosome names", {
  mf <- tempfile(fileext = ".map")
  writeLines(c("chr\tid\tcM\tbp", "chr1A\tm1\t0.1\t100"), mf)
  expect_error(PredictProR:::validate_map_file(mf), "must be numeric")
  writeLines(c("chr1A\tm1\t0.1\t100", "chr1A\tm2\t0.2\t200"), mf)
  expect_message(PredictProR:::validate_map_file(mf), "validated")
})
