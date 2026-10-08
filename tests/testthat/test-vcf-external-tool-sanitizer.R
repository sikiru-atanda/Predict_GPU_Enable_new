make_external_tool_vcf <- function(path, lines) {
  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

test_that("VCF sanitizer removes rows that commonly break PLINK and Beagle", {
  vcf <- make_external_tool_vcf(
    tempfile("predictpror_bad_", fileext = ".vcf"),
    c(
      "##fileformat=VCFv4.2",
      "##source=predictpror-test",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
      "1\t100\trs_valid\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t0\trs_bad_pos\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
      ".\t105\trs_bad_chrom\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t110\trs_missing_alt\tA\t.\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t120\trs_equal\tC\tC\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t130\trs_multiallelic\tC\tT,G\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t140\trs_symbolic\tC\t<DEL>\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t150\trs_bad_gt\tA\tT\t.\tPASS\t.\tGT\t0/2\t0/1",
      "2\t100\trs_valid_same_pos_other_chr\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
      "1\t160\trs_duplicate_first\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t160\trs_duplicate_second\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t170\trs_conflict_a\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t170\trs_conflict_b\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
    )
  )

  res <- PredictProR:::sanitize_vcf_for_external_tools(vcf, mode = "plink")
  expect_true(res$changed)
  expect_equal(res$metrics$input_variants, 13L)
  expect_equal(res$metrics$output_variants, 3L)
  expect_equal(res$metrics$removed_bad_chromosome, 1L)
  expect_equal(res$metrics$removed_bad_position, 1L)
  expect_equal(res$metrics$removed_missing_allele, 1L)
  expect_equal(res$metrics$removed_ref_alt_equal, 1L)
  expect_equal(res$metrics$removed_non_acgt_or_multiallelic, 2L)
  expect_equal(res$metrics$removed_unsupported_genotype, 1L)
  expect_equal(res$metrics$removed_conflicting_alleles, 2L)
  expect_equal(res$metrics$removed_duplicate_marker, 1L)

  cleaned <- PredictProR:::vcf_read_table(res$path)
  expect_identical(
    cleaned$ID,
    c("rs_valid", "rs_valid_same_pos_other_chr", "rs_duplicate_first")
  )
  expect_true(all(c("##fileformat=VCFv4.2", "##source=predictpror-test") %in% PredictProR:::vcf_read_all_lines(res$path)))
})

test_that("VCF sanitizer reads gzipped VCF without forcing a rewrite when clean", {
  vcf <- tempfile("predictpror_clean_", fileext = ".vcf.gz")
  con <- gzfile(vcf, open = "wt")
  closed <- FALSE
  on.exit(if (!closed) close(con), add = TRUE)
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
      "chr1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:12"
    ),
    con,
    useBytes = TRUE
  )
  close(con)
  closed <- TRUE

  res <- PredictProR:::sanitize_vcf_for_external_tools(vcf, mode = "beagle")
  expect_false(res$changed)
  expect_identical(res$path, normalizePath(vcf, winslash = "/", mustWork = TRUE))
  expect_equal(res$metrics$output_variants, 1L)
})

test_that("Beagle jar path can be supplied explicitly without download", {
  jar <- tempfile("beagle.", fileext = ".jar")
  writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04)), jar)

  old <- Sys.getenv("PREDICTPRO_BEAGLE_JAR", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("PREDICTPRO_BEAGLE_JAR")
    } else {
      Sys.setenv(PREDICTPRO_BEAGLE_JAR = old)
    }
  }, add = TRUE)
  Sys.setenv(PREDICTPRO_BEAGLE_JAR = jar)

  expect_identical(
    PredictProR:::download_beagle(output_dir = tempfile("beagle-cache-")),
    normalizePath(jar, winslash = "/", mustWork = TRUE)
  )
})
