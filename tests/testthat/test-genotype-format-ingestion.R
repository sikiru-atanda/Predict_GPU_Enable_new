make_hapmap_fixture <- function() {
  data.frame(
    `rs#` = c("m_major_alt", "m_iupac", "m_missing"),
    alleles = c("A/G", "C/T", "G/A"),
    chrom = c("1", "1", "2"),
    pos = c(10L, 20L, 30L),
    strand = "+",
    `assembly#` = NA,
    center = NA,
    protLSID = NA,
    assayLSID = NA,
    panelLSID = NA,
    QCcode = NA,
    s1 = c("GG", "Y", "G/G"),
    s2 = c("G", "T/T", "G/A"),
    s3 = c("A/G", "C", "N"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

no_filter_hapmap <- function(hapmap, recode_format = "0,1,2") {
  hmp_qc_recode(
    hapmap = hapmap,
    maf_threshold = 0,
    het_threshold = NULL,
    ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0,
    impute = FALSE,
    recode_format = recode_format,
    message = FALSE
  )
}

test_that("plain and gzip HapMap files are read with the standard columns", {
  hapmap <- make_hapmap_fixture()
  plain <- tempfile(fileext = ".hmp.txt")
  gz <- tempfile(fileext = ".hmp.txt.gz")
  data.table::fwrite(hapmap, plain, sep = "\t", na = "NA")
  con <- gzfile(gz, "wt")
  writeLines(readLines(plain), con)
  close(con)

  plain_read <- read_hapmap_file(plain)
  gz_read <- read_hapmap_file(gz)
  expect_identical(names(plain_read), names(hapmap))
  expect_equal(as.character(plain_read$chrom), hapmap$chrom)
  expect_equal(as.data.frame(gz_read), as.data.frame(plain_read))
})

test_that("a zip archive containing one HapMap file is read deterministically", {
  skip_if(!nzchar(Sys.which("zip")), "A zip utility is not available")
  hapmap <- make_hapmap_fixture()
  source_dir <- tempfile("hapmap-zip-source-")
  dir.create(source_dir)
  plain <- file.path(source_dir, "fixture.hmp.txt")
  archive <- tempfile(fileext = ".zip")
  data.table::fwrite(hapmap, plain, sep = "\t", na = "NA")
  old_dir <- setwd(source_dir)
  on.exit(setwd(old_dir), add = TRUE)
  utils::zip(archive, files = basename(plain))

  observed <- read_hapmap_file(archive)
  expect_equal(as.data.frame(observed), as.data.frame(read_hapmap_file(plain)))
})

test_that("common GAPIT HapMap header aliases are normalized by position", {
  hapmap <- make_hapmap_fixture()
  names(hapmap)[c(1L, 6L, 10L)] <- c("rs", "assembly", "panel")
  file <- tempfile(fileext = ".hmp.txt")
  data.table::fwrite(hapmap, file, sep = "\t", na = "NA")
  read <- read_hapmap_file(file)
  expect_identical(
    names(read)[1:11],
    c(
      "rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
      "protLSID", "assayLSID", "panelLSID", "QCcode"
    )
  )
  out <- no_filter_hapmap(read)
  expect_named(
    out$format_report$normalized_header_aliases,
    c("rs#", "assembly#", "panelLSID")
  )
})

test_that("HapMap alternate-allele coding is stable and both schemes agree", {
  hapmap <- make_hapmap_fixture()
  dosage <- no_filter_hapmap(hapmap, "0,1,2")
  centred <- no_filter_hapmap(hapmap, "-1,0,1")

  expect_identical(dosage$coding, "alternate_allele_dosage")
  expect_equal(
    unclass(dosage$snps_matrix),
    rbind(s1 = c(2, 1, 0), s2 = c(2, 2, 1), s3 = c(1, 0, NA)),
    ignore_attr = TRUE
  )
  expect_equal(
    unclass(centred$snps_matrix),
    unclass(dosage$snps_matrix) - 1,
    ignore_attr = TRUE
  )
})

test_that("HapMap to VCF preserves declared REF and ALT and round-trips dosage", {
  hapmap <- make_hapmap_fixture()
  direct <- no_filter_hapmap(hapmap)
  vcf <- tempfile(fileext = ".vcf")
  converted <- convert_hapmap_to_vcf(hapmap, vcf)
  lines <- readLines(vcf)
  expect_true("##fileformat=VCFv4.3" %in% lines)

  table <- PredictProR:::vcf_read_table(vcf)
  expect_identical(table$REF, c("A", "C", "G"))
  expect_identical(table$ALT, c("G", "T", "A"))
  expect_equal(attr(converted, "conversion_report")$retained_markers, 3L)

  roundtrip <- vcf_qc_recode(
    basename(vcf), dirname(vcf), maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0, impute = FALSE, message = FALSE
  )
  expect_equal(
    unclass(roundtrip$snps_matrix),
    unclass(direct$snps_matrix),
    ignore_attr = TRUE
  )
})

test_that("invalid HapMap calls fail loudly and non-biallelic policy is reported", {
  hapmap <- make_hapmap_fixture()
  hapmap$s1[[1L]] <- "C"
  expect_error(no_filter_hapmap(hapmap), "inconsistent with the declared allele pair")

  hapmap <- make_hapmap_fixture()
  hapmap$alleles[[2L]] <- "A/C/G"
  out <- no_filter_hapmap(hapmap)
  expect_equal(out$format_report$dropped_non_biallelic, 1L)
  expect_error(
    hmp_qc_recode(
      hapmap = hapmap, non_biallelic = "error", maf_threshold = 0,
      het_threshold = NULL, ind_call_rate_threshold = 0,
      snp_call_rate_threshold = 0, message = FALSE
    ),
    "exactly two"
  )
})

make_vcf_format_fixture <- function(path, version = "4.3", format = "GT") {
  samples <- switch(
    format,
    GT = c("0|0", "1/0", "1/1"),
    `GT:DS` = c("0|0:0", "1/0:1", "1/1:2"),
    `AD:GT:DP` = c("9,0:0|0:9", "4,5:1/0:9", "0,8:1/1:8")
  )
  writeLines(
    c(
      paste0("##fileformat=VCFv", version),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
      paste(c("1", "100", "m1", "A", "G", ".", "PASS", ".", format, samples), collapse = "\t")
    ),
    path
  )
  path
}

test_that("VCF 4.1 through 4.5 and common FORMAT layouts give the same dosage", {
  observed <- list()
  k <- 0L
  for (version in c("4.1", "4.2", "4.3", "4.4", "4.5")) {
    for (format in c("GT", "GT:DS", "AD:GT:DP")) {
      k <- k + 1L
      vcf <- make_vcf_format_fixture(tempfile(fileext = ".vcf"), version, format)
      out <- vcf_qc_recode(
        basename(vcf), dirname(vcf), maf_threshold = 0,
        het_threshold = NULL, ind_call_rate_threshold = 0,
        snp_call_rate_threshold = 0, message = FALSE
      )
      observed[[k]] <- unname(out$snps_matrix[, 1L])
    }
  }
  expect_true(all(vapply(observed, identical, logical(1L), c(0, 1, 2))))
})

test_that("VCF fileformat declaration is mandatory and bounded", {
  missing_version <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
      "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/1"
    ),
    missing_version
  )
  expect_error(PredictProR:::vcf_read_header(missing_version), "fileformat")

  unsupported <- make_vcf_format_fixture(tempfile(fileext = ".vcf"), "4.0", "GT")
  expect_error(PredictProR:::vcf_read_header(unsupported), "VCFv4.1 through VCFv4.5")
})

test_that("VCF call-rate thresholds are minimum called fractions", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
      "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1",
      "1\t20\tm2\tC\tT\t.\tPASS\t.\tGT\t0/0\t./.\t./.\t./.",
      "1\t30\tm3\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./."
    ),
    vcf
  )
  marker_filtered <- vcf_qc_recode(
    basename(vcf), dirname(vcf), maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0.9, message = FALSE
  )
  expect_identical(colnames(marker_filtered$snps_matrix), "m1")

  sample_filtered <- vcf_qc_recode(
    basename(vcf), dirname(vcf), maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0.9,
    snp_call_rate_threshold = 0, message = FALSE
  )
  expect_identical(rownames(sample_filtered$snps_matrix), c("s1", "s2", "s3"))

  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
      "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t./.\t./.\t1/1",
      "1\t20\tm2\tC\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1\t./.",
      "1\t30\tm3\tG\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t./.\t./."
    ),
    vcf
  )
  strict_sample_filtered <- vcf_qc_recode(
    basename(vcf), dirname(vcf), maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0.9,
    snp_call_rate_threshold = 0, message = FALSE
  )
  expect_identical(rownames(strict_sample_filtered$snps_matrix), "s1")
})

test_that("VCF format preflight rejects haploid and multiallelic dosage rows", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
      "1\t10\tvalid\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
      "1\t20\thaploid\tC\tT\t.\tPASS\t.\tGT\t0\t1",
      "1\t30\tmulti\tG\tA,T\t.\tPASS\t.\tGT\t0/2\t1/2"
    ),
    vcf
  )
  clean <- PredictProR:::sanitize_vcf_for_external_tools(vcf, mode = "beagle")
  expect_equal(clean$metrics$output_variants, 1L)
  expect_equal(clean$metrics$removed_unsupported_genotype, 1L)
  expect_equal(clean$metrics$removed_non_acgt_or_multiallelic, 1L)
})
