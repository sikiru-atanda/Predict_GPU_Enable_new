make_beagle_test_vcf <- function(path) {
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
      "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "1\t20\tm2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0"
    ),
    path
  )
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

make_fake_jar <- function(path = tempfile(fileext = ".jar")) {
  writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04)), path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# shQuote() uses double quotes on Windows and single quotes on Unix.
unquote_beagle_arg <- function(x) {
  if (grepl("^'.*'$", x)) {
    return(gsub("'\\''", "'", substr(x, 2L, nchar(x) - 1L), fixed = TRUE))
  }
  sub('"$', "", sub('^"', "", x))
}

test_that("Beagle command execution preserves paths and returns provenance", {
  vcf <- make_beagle_test_vcf(tempfile(fileext = ".vcf"))
  jar <- make_fake_jar()
  java <- tempfile("fake java ", fileext = ".exe")
  file.create(java)
  output_dir <- tempfile("beagle output ")
  dir.create(output_dir)
  captured <- NULL

  testthat::local_mocked_bindings(
    beagle_system2 = function(command, args, stdout = TRUE, stderr = TRUE) {
      captured <<- list(command = command, args = vapply(args, unquote_beagle_arg, character(1L)))
      gt <- sub("^gt=", "", captured$args[startsWith(captured$args, "gt=")])
      out <- sub("^out=", "", captured$args[startsWith(captured$args, "out=")])
      con <- gzfile(paste0(out, ".vcf.gz"), "wt")
      writeLines(readLines(gt), con)
      close(con)
      writeLines("fake Beagle log", paste0(out, ".log"))
      "fake Beagle success"
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::run_beagle(
    vcf_file = vcf,
    output_prefix = file.path(output_dir, "result with spaces"),
    beagle_jar = jar,
    generate_map = FALSE,
    java_bin = java,
    java_memory = "1g",
    nthreads = 2,
    seed = 42L
  )
  expect_equal(out$status, 0L)
  expect_true(file.exists(out$output_vcf))
  expect_true(file.exists(out$log_file))
  expect_true(file.exists(out$beagle_log_file))
  expect_true(any(captured$args == paste0("gt=", vcf)))
  expect_true(any(grepl("out=.*result with spaces$", captured$args)))
  expect_true("seed=42" %in% captured$args)
  expect_true("em=true" %in% captured$args)
  expect_equal(out$input_preflight$output_variants, 2L)
})

test_that("Beagle nonzero status and missing output fail loudly with a wrapper log", {
  vcf <- make_beagle_test_vcf(tempfile(fileext = ".vcf"))
  jar <- make_fake_jar()
  java <- tempfile(fileext = ".exe")
  file.create(java)
  prefix <- file.path(tempfile("beagle-failure-"), "result")

  testthat::local_mocked_bindings(
    beagle_system2 = function(command, args, stdout = TRUE, stderr = TRUE) {
      structure("simulated failure", status = 9L)
    },
    .package = "PredictProR"
  )
  expect_error(
    PredictProR:::run_beagle(
      vcf, prefix, jar, generate_map = FALSE, java_bin = java
    ),
    "Beagle failed"
  )
  expect_true(file.exists(paste0(prefix, "_log.txt")))

  testthat::local_mocked_bindings(
    beagle_system2 = function(command, args, stdout = TRUE, stderr = TRUE) "status zero",
    .package = "PredictProR"
  )
  expect_error(
    PredictProR:::run_beagle(
      vcf, paste0(prefix, "-missing"), jar,
      generate_map = FALSE, java_bin = java
    ),
    "did not create the expected VCF"
  )
})

test_that("the public HapMap Beagle API converts, runs, and recodes consistently", {
  hapmap <- data.frame(
    `rs#` = c("m1", "m2"), alleles = c("A/G", "C/T"),
    chrom = "1", pos = c(10L, 20L), strand = "+", `assembly#` = NA,
    center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA,
    QCcode = NA, s1 = c("A", "Y"), s2 = c("R", "T"),
    s3 = c("G", "C"), check.names = FALSE, stringsAsFactors = FALSE
  )
  jar <- make_fake_jar()
  java <- tempfile(fileext = ".exe")
  file.create(java)

  testthat::local_mocked_bindings(
    beagle_system2 = function(command, args, stdout = TRUE, stderr = TRUE) {
      args <- vapply(args, unquote_beagle_arg, character(1L))
      gt <- sub("^gt=", "", args[startsWith(args, "gt=")])
      out <- sub("^out=", "", args[startsWith(args, "out=")])
      con <- gzfile(paste0(out, ".vcf.gz"), "wt")
      writeLines(readLines(gt), con)
      close(con)
      writeLines("fake", paste0(out, ".log"))
      "ok"
    },
    .package = "PredictProR"
  )

  result <- impute_genotypes_with_beagle(
    input = hapmap,
    input_format = "hapmap",
    output_prefix = file.path(tempfile("beagle-hapmap-"), "out"),
    beagle_jar = jar,
    java_installed = TRUE,
    java_bin = java,
    generate_map = FALSE,
    return_genotypes = TRUE,
    qc_args = list(
      maf_threshold = 0, het_threshold = NULL,
      ind_call_rate_threshold = 0, snp_call_rate_threshold = 0
    )
  )
  expect_s3_class(result, "predictpror_beagle_result")
  expect_equal(result$conversion_report$output_format, "VCFv4.3")
  direct <- hmp_qc_recode(
    hapmap = hapmap, maf_threshold = 0, het_threshold = NULL,
    ind_call_rate_threshold = 0, snp_call_rate_threshold = 0,
    impute = FALSE, message = FALSE
  )
  expect_equal(
    unclass(result$genotypes$snps_matrix),
    unclass(direct$snps_matrix),
    ignore_attr = TRUE
  )
})

test_that("invalid jars and removed Beagle 4-era options are rejected", {
  bad_jar <- tempfile(fileext = ".jar")
  writeLines("not a jar", bad_jar)
  expect_error(PredictProR:::validate_beagle_jar(bad_jar), "ZIP/JAR signature")

  vcf <- make_beagle_test_vcf(tempfile(fileext = ".vcf"))
  java <- tempfile(fileext = ".exe")
  file.create(java)
  expect_error(
    PredictProR:::run_beagle(
      vcf, tempfile(), make_fake_jar(), generate_map = FALSE,
      java_bin = java, markers_file = tempfile()
    ),
    "not Beagle 5.4 arguments"
  )
})

test_that("model_execute exposes Beagle preprocessing and fails early on invalid use", {
  expect_true("beagle_options" %in% names(formals(model_execute)))
  expect_error(
    model_execute(imputation_method = "beagle", impute = FALSE),
    "requires impute = TRUE"
  )
  expect_error(
    model_execute(imputation_method = "beagle", impute = TRUE),
    "exactly one raw genotype source"
  )
  expect_error(
    model_execute(
      imputation_method = "beagle", impute = TRUE,
      beagle_options = list("not named")
    ),
    "must be a named list"
  )
})

test_that("model_execute keeps Beagle intermediates out of the working directory", {
  vcf <- make_beagle_test_vcf(tempfile(fileext = ".vcf"))
  seen <- list()
  testthat::local_mocked_bindings(
    impute_genotypes_with_beagle = function(input, output_prefix = NULL, ...) {
      seen[[length(seen) + 1L]] <<- output_prefix
      stop("captured")
    },
    .package = "PredictProR"
  )
  run <- function(opts) {
    expect_error(
      model_execute(
        vcf_file_name = basename(vcf), vcf_file_path = dirname(vcf),
        imputation_method = "beagle", impute = TRUE, beagle_options = opts
      ),
      "captured"
    )
  }
  run(list())
  tmp <- normalizePath(tempdir(), winslash = "/")
  expect_true(startsWith(normalizePath(dirname(dirname(seen[[1L]])), winslash = "/", mustWork = FALSE), tmp))
  user_prefix <- file.path(tempdir(), "my_beagle_out")
  run(list(output_prefix = user_prefix))
  expect_identical(seen[[2L]], user_prefix)
})
