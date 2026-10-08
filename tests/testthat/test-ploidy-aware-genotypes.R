make_polyploid_hapmap <- function() {
  data.frame(
    `rs#` = c("p4_m1", "p4_m2"),
    alleles = c("A/G", "C/T"), chrom = "1", pos = c(10L, 20L),
    strand = "+", `assembly#` = NA, center = NA, protLSID = NA,
    assayLSID = NA, panelLSID = NA, QCcode = NA,
    s0 = c("AAAA", "CCCC"),
    s1 = c("AAAG", "CCCT"),
    s2 = c("AAGG", "CCTT"),
    s3 = c("AGGG", "CTTT"),
    s4 = c("GGGG", "NNNN"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

poly_no_filter <- function(hapmap, recode_format = "alt_dosage", impute = FALSE,
                           imputation_method = "mean") {
  hmp_qc_recode(
    hapmap = hapmap, ploidy = 4L, maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0, impute = impute,
    imputation_method = imputation_method, recode_format = recode_format,
    message = FALSE
  )
}

test_that("tetraploid HapMap dosage, centered coding, and frequency coding are exact", {
  hapmap <- make_polyploid_hapmap()
  dosage <- poly_no_filter(hapmap)
  centered <- poly_no_filter(hapmap, "centered_dosage")
  frequency <- poly_no_filter(hapmap, "allele_frequency")

  expect_equal(unname(dosage$snps_matrix[, "p4_m1"]), 0:4)
  expect_equal(unname(dosage$snps_matrix[, "p4_m2"]), c(0:3, NA))
  expect_equal(unname(centered$snps_matrix[, "p4_m1"]), -2:2)
  expect_equal(unname(frequency$snps_matrix[, "p4_m1"]), (0:4) / 4)
  expect_identical(dosage$ploidy, 4L)
  expect_identical(attr(dosage$snps_matrix, "ploidy"), 4L)
  expect_identical(dosage$coding_contract, "alt_dosage")
})

test_that("polyploid HapMap IUPAC heterozygotes fail because dosage is ambiguous", {
  hapmap <- make_polyploid_hapmap()
  hapmap$s1[[1L]] <- "R"
  expect_error(poly_no_filter(hapmap), "dosage-ambiguous for ploidy 4")
  accepted <- hmp_qc_recode(
    hapmap = hapmap, ploidy = 4L, invalid_call = "missing",
    maf_threshold = 0, het_threshold = NULL,
    ind_call_rate_threshold = 0, snp_call_rate_threshold = 0,
    impute = FALSE, message = FALSE
  )
  expect_equal(accepted$format_report$ambiguous_polyploid_calls, 1L)
  expect_true(is.na(accepted$snps_matrix["s1", "p4_m1"]))
})

test_that("tetraploid HapMap converts to VCF and round-trips all dosage classes", {
  hapmap <- make_polyploid_hapmap()
  direct <- poly_no_filter(hapmap)
  vcf <- tempfile(fileext = ".vcf")
  converted <- convert_hapmap_to_vcf(hapmap, vcf, ploidy = 4L)
  lines <- readLines(vcf)
  expect_true(any(grepl("0/0/0/0", lines, fixed = TRUE)))
  expect_true(any(grepl("0/0/0/1", lines, fixed = TRUE)))
  expect_true(any(grepl("0/0/1/1", lines, fixed = TRUE)))
  expect_true(any(grepl("0/1/1/1", lines, fixed = TRUE)))
  expect_true(any(grepl("1/1/1/1", lines, fixed = TRUE)))
  expect_identical(attr(converted, "conversion_report")$ploidy, 4L)

  roundtrip <- vcf_qc_recode(
    basename(vcf), dirname(vcf), ploidy = "auto", maf_threshold = 0,
    het_threshold = NULL, ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0, impute = FALSE, message = FALSE
  )
  expect_identical(roundtrip$ploidy, 4L)
  expect_equal(
    unclass(roundtrip$snps_matrix), unclass(direct$snps_matrix),
    ignore_attr = TRUE
  )
})

test_that("triploid and hexaploid VCF GT calls infer ploidy and dosage", {
  for (ploidy in c(3L, 6L)) {
    samples <- paste0("s", 0:ploidy)
    calls <- vapply(0:ploidy, function(d) {
      paste(c(rep("0", ploidy - d), rep("1", d)), collapse = "/")
    }, character(1L))
    vcf <- tempfile(fileext = ".vcf")
    writeLines(
      c(
        "##fileformat=VCFv4.3",
        paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"),
        paste(c("1", "10", "m1", "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
      ),
      vcf
    )
    out <- vcf_qc_recode(
      basename(vcf), dirname(vcf), ploidy = "auto", maf_threshold = 0,
      het_threshold = NULL, ind_call_rate_threshold = 0,
      snp_call_rate_threshold = 0, impute = FALSE, message = FALSE
    )
    expect_identical(out$ploidy, ploidy)
    expect_equal(unname(out$snps_matrix[, "m1"]), 0:ploidy)
  }
})

test_that("native polyploid imputation respects dosage bounds", {
  mean_imputed <- poly_no_filter(make_polyploid_hapmap(), impute = TRUE, imputation_method = "mean")
  expect_equal(mean_imputed$snps_matrix["s4", "p4_m2"], 1.5)
  expect_false(anyNA(mean_imputed$snps_matrix))
  expect_true(all(mean_imputed$snps_matrix >= 0 & mean_imputed$snps_matrix <= 4))

  mode_imputed <- poly_no_filter(make_polyploid_hapmap(), impute = TRUE, imputation_method = "mode")
  expect_equal(mode_imputed$snps_matrix["s4", "p4_m2"], 0)
})

test_that("mixed VCF ploidy and polyploid Beagle requests fail loudly", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.3",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
      "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0/1/1"
    ),
    vcf
  )
  expect_error(
    vcf_qc_recode(basename(vcf), dirname(vcf), ploidy = "auto", message = FALSE),
    "mixed GT ploidies"
  )
  expect_error(
    impute_genotypes_with_beagle(
      make_polyploid_hapmap(), input_format = "hapmap", ploidy = 4L
    ),
    "diploid genotypes only.*imputation_method = \"knn\""
  )
})

test_that("Beagle stops for polyploid HapMap/CSV with ploidy = 'auto' and suggests KNN", {
  knn_msg <- "diploid genotypes only.*imputation_method = \"knn\""
  # HapMap calls with four allele copies imply ploidy 4
  expect_error(
    impute_genotypes_with_beagle(make_polyploid_hapmap(), input_format = "hapmap", ploidy = "auto"),
    knn_msg
  )
  # CSV dosage above 2
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(CHROM = 1L, POS = 10L, ID = "m1", REF = "A", ALT = "G", QUAL = ".",
                              FILTER = "PASS", INFO = ".", FORMAT = "GT", s0 = 0, s1 = 4, s2 = 2,
                              check.names = FALSE), csv, row.names = FALSE)
  expect_error(impute_genotypes_with_beagle(csv, input_format = "csv", ploidy = "auto"), knn_msg)
  # diploid inputs are not flagged
  expect_null(PredictProR:::gp_beagle_hapmap_call_ploidy(data.frame(matrix("AG", 2, 13))))
  csv2 <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(CHROM = 1L, POS = 10L, ID = "m1", REF = "A", ALT = "G", QUAL = ".",
                              FILTER = "PASS", INFO = ".", FORMAT = "GT", s0 = 0, s1 = 1, s2 = 2,
                              check.names = FALSE), csv2, row.names = FALSE)
  expect_silent(PredictProR:::gp_beagle_check_csv_dosage(csv2))
})

test_that("model_execute stops Beagle on a polyploid dosage matrix and suggests KNN", {
  x <- matrix(c(0, 1, 2, 3, 4, 2, 1, 3), nrow = 4, ncol = 6,
              dimnames = list(paste0("g", 1:4), paste0("m", 1:6)))
  ph <- data.frame(GID = rownames(x), Yield = c(1, 2, NA, 3))
  expect_error(
    model_execute(pheno_data = ph, geno_data = x, gen_name = "GID", response = "Yield",
                  GS_model = "GBLUP_BRR", ploidy = 4L, impute = TRUE, imputation_method = "beagle",
                  message = FALSE),
    "diploid genotypes only.*imputation_method = \"knn\""
  )
})

test_that("GPNet is rejected for polyploid data with a clear message", {
  msg <- "GPNet is not supported for polyploid data \\(ploidy 4\\)"
  expect_error(PredictProR:::gp_reject_polyploid_unsupported_models("GPNet", 4L), msg)
  expect_error(PredictProR:::gp_reject_polyploid_unsupported_models(list(NULL, c("BayesB", "gp_dkl")), 4L), msg)
  expect_true(PredictProR:::gp_reject_polyploid_unsupported_models("GPNet", 2L))
  expect_true(PredictProR:::gp_reject_polyploid_unsupported_models(c("DenseNeuralNet", "TabNet"), 6L))
  expect_true(PredictProR:::gp_reject_polyploid_unsupported_models("GPNet", "auto"))

  x <- matrix(c(0, 1, 2, 3, 4, 2, 1, 3), nrow = 4, ncol = 6,
              dimnames = list(paste0("g", 1:4), paste0("m", 1:6)))
  ph <- data.frame(GID = rownames(x), Yield = c(1, 2, NA, 3))
  expect_error(
    model_execute(pheno_data = ph, geno_data = x, gen_name = "GID", response = "Yield",
                  GS_model = "GPNet", ploidy = 4L, message = FALSE),
    msg
  )
  # ploidy read from tetraploid VCF GT calls
  vcf <- tempfile(fileext = ".vcf")
  gt <- function(d) paste(c(rep("0", 4 - d), rep("1", d)), collapse = "/")
  writeLines(c(
    "##fileformat=VCFv4.3",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", rownames(x)), collapse = "\t"),
    vapply(seq_len(ncol(x)), function(j) paste(c("1", j * 10, colnames(x)[j], "A", "G", ".", "PASS", ".", "GT",
                                                 vapply(x[, j], gt, "")), collapse = "\t"), "")
  ), vcf)
  expect_error(
    model_execute(pheno_data = ph, vcf_file_name = basename(vcf), vcf_file_path = dirname(vcf),
                  gen_name = "GID", response = "Yield", GS_model = "GPNet", ploidy = "auto",
                  maf_threshold = 0, het_threshold = NULL, message = FALSE),
    msg
  )
})

test_that("native imputation API preserves dimensions, names, and ploidy", {
  x <- matrix(
    c(0, 1, 2, 3, 4, NA, 2, 2, 1, 0),
    nrow = 5L,
    dimnames = list(paste0("id", 1:5), c("m1", "m2"))
  )
  out <- impute_genotypes_native(x, ploidy = 4L, method = "mean")
  expect_identical(dimnames(out$snps_matrix), dimnames(x))
  expect_identical(out$ploidy, 4L)
  expect_identical(out$imputed_cells, 1L)
  expect_false(anyNA(out$snps_matrix))
  expect_identical(attr(out$snps_matrix, "ploidy"), 4L)
  expect_error(impute_genotypes_native(x, ploidy = "auto"), "Supply ploidy explicitly")
})

test_that("polyploid QC metrics use dosage divided by ploidy", {
  x <- matrix(
    c(0, 1, 2, 3, 4, 0, 0, 4, 4, NA),
    nrow = 5L,
    dimnames = list(paste0("id", 1:5), c("balanced", "polarized"))
  )
  metrics <- PredictProR:::geno_qc_metrics(x, ploidy = 4L)
  expect_equal(unname(metrics$maf[["balanced"]]), 0.5)
  expect_equal(unname(metrics$heterozygosity[["balanced"]]), 0.6)
  expect_equal(unname(metrics$maf[["polarized"]]), 0.5)
  expect_equal(unname(metrics$heterozygosity[["polarized"]]), 0)
})

test_that("autopolyploid VanRaden equals its closed-form dosage definition", {
  x <- matrix(
    c(
      0, 1, 4,
      1, 2, 3,
      2, 3, 2,
      3, 4, 1,
      4, 0, 0
    ),
    nrow = 5L, byrow = TRUE,
    dimnames = list(paste0("id", 1:5), paste0("m", 1:3))
  )
  p <- colMeans(x) / 4
  z <- sweep(x, 2L, colMeans(x), "-")
  expected <- tcrossprod(z) / sum(4 * p * (1 - p))
  attr(expected, "ploidy") <- 4L
  observed <- PredictProR:::grm_calculation(
    x, method = "VanRaden", ploidy = 4L, backend = "r"
  )
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(attr(observed, "ploidy"), 4L)
  expect_error(
    PredictProR:::grm_calculation(x, method = "Yang", ploidy = 4L),
    "diploid-only"
  )
  expect_error(
    PredictProR:::grm_calculation(x, method = "Dominance", ploidy = 4L),
    "diploid-only"
  )
})

test_that("feature-selected polyploid GRM is rebuilt with the same ploidy", {
  x <- matrix(
    c(0, 1, 4, 1, 2, 3, 2, 3, 2, 3, 4, 1, 4, 0, 0),
    nrow = 5L, byrow = TRUE,
    dimnames = list(paste0("id", 1:5), paste0("m", 1:3))
  )
  attr(x, "ploidy") <- 4L
  selected <- PredictProR:::gp_feature_subset_source_matrices(
    source_matrices = list(geno_data = x),
    selected_by_source = list(geno_data = c("m1", "m3"))
  )
  bank <- PredictProR:::gp_feature_build_kernel_bank(
    selected_sources = selected,
    gmatrix_method = "VanRaden",
    ploidy = "auto"
  )
  direct <- PredictProR:::grm_calculation(
    x[, c("m1", "m3"), drop = FALSE],
    method = "VanRaden", ploidy = 4L
  )
  expect_equal(bank$gmatrix_model_ready, direct, tolerance = 1e-12)
  expect_identical(attr(bank$gmatrix_model_ready, "ploidy"), 4L)
})

test_that("tetraploid MET ML keeps the user's VanRaden GRM (was replaced by diploid-only Yang)", {
  if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
    testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
  }
  py <- PredictProR:::gp_preferred_python(purpose = "ml")
  testthat::skip_if(is.null(py) || !file.exists(py), "No ML Python runtime found")
  withr::local_dir(withr::local_tempdir())
  set.seed(61)
  ids <- sprintf("P%02d", 1:40)
  G <- matrix(stats::rbinom(40 * 120, 4, 0.4), 40, 120, dimnames = list(ids, paste0("m", 1:120)))
  G <- G[, apply(G, 2, stats::var) > 0]
  ph <- expand.grid(GID = ids, Env = c("E1", "E2"), stringsAsFactors = FALSE)
  ph$Yield <- as.numeric(scale(G[ph$GID, 1:10] %*% rnorm(10))) + stats::rnorm(nrow(ph))
  ph$Yield[ph$GID %in% ids[1:6]] <- NA
  for (method in list("VanRaden", NULL)) {
    res <- model_execute(
      pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield",
      heter_groups = "Env", met_ml_dl = TRUE, GS_model = "RandomForest",
      ploidy = 4L, gmatrix_method = method, het_threshold = NULL,
      ntree = 20L, n_bootstrap = 3, message = FALSE, verbose = FALSE
    )
    pred <- res$model_results$predicted_values
    expect_true(is.data.frame(pred) && any(is.finite(pred$Predicted_value)))
  }
})

test_that("model_execute exposes a package-wide ploidy argument", {
  expect_true("ploidy" %in% names(formals(model_execute)))
  expect_identical(formals(model_execute)$ploidy, "auto")
})

test_that("CSV conversion writes ploidy-length GT calls from explicit dosage coding", {
  csv <- tempfile(fileext = ".csv")
  input <- data.frame(
    CHROM = 1L, POS = 10L, ID = "p4_m1", REF = "A", ALT = "G",
    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
    s0 = 0, s1 = 1, s2 = 2, s3 = 3, s4 = 4,
    check.names = FALSE
  )
  utils::write.csv(input, csv, row.names = FALSE, na = ".")
  lines <- PredictProR:::convert_csv_to_vcf(csv, ploidy = 4L)
  record <- strsplit(lines[[length(lines)]], "\t", fixed = TRUE)[[1L]]
  expect_identical(
    record[10:14],
    c("0/0/0/0", "0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1")
  )

  input$s0 <- -2
  input$s1 <- -1
  input$s2 <- 0
  input$s3 <- 1
  input$s4 <- 2
  utils::write.csv(input, csv, row.names = FALSE, na = ".")
  centered <- PredictProR:::convert_csv_to_vcf(
    csv, ploidy = 4L, input_coding = "centered_dosage"
  )
  centered_record <- strsplit(
    centered[[length(centered)]], "\t", fixed = TRUE
  )[[1L]]
  expect_identical(centered_record[10:14], record[10:14])
})

test_that("coding detection requires a defensible ploidy contract", {
  x <- matrix(0:4, ncol = 1L)
  expect_identical(
    detect_genomic_coding(x),
    "Potential polyploid ALT dosage; supply ploidy"
  )
  expect_identical(
    detect_genomic_coding(x, ploidy = 4L),
    "ALT dosage (0..4), ploidy 4"
  )
  attr(x, "ploidy") <- 4L
  expect_identical(
    detect_genomic_coding(x),
    "ALT dosage (0..4), ploidy 4"
  )
})

test_that("monomorphic removal validates and preserves polyploid metadata", {
  x <- cbind(variable = 0:4, monomorphic = rep(4, 5), missing = rep(NA_real_, 5))
  rownames(x) <- paste0("id", seq_len(nrow(x)))
  out <- Remove_NA_Mono_SNP(x, ploidy = 4L, message = FALSE)
  expect_identical(colnames(out), "variable")
  expect_identical(attr(out, "ploidy"), 4L)
  expect_error(
    Remove_NA_Mono_SNP(x, ploidy = 3L, message = FALSE),
    "dosage.*\\[0, 3\\]"
  )
})
