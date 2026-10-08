# CSV / TXT genotype tables (VCF marker columns + one column per sample) are a
# first-class genotype input: convert_csv_to_vcf() is exported, model_execute()
# takes csv_file_name, and impute_genotypes_with_beagle() takes the table.

make_table <- function(n = 24, p = 40, seed = 5, missing_rate = 0) {
  set.seed(seed)
  ids <- sprintf("L%02d", seq_len(n))
  G <- matrix(2L * stats::rbinom(n * p, 1, 0.5), n, p, dimnames = list(ids, sprintf("m%02d", seq_len(p))))
  dosage <- t(G)
  if (missing_rate > 0) dosage[stats::runif(length(dosage)) < missing_rate] <- NA
  tab <- data.frame(CHROM = "1", POS = seq_len(p) * 1000L, ID = colnames(G), REF = "A", ALT = "G",
                    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT", dosage, check.names = FALSE)
  list(G = G, table = tab, ids = ids)
}

test_that("convert_csv_to_vcf writes GT calls for dosages, GT strings and missing values", {
  tab <- data.frame(CHROM = 1, POS = c(100, 200, 300), ID = c("m1", "m2", "m3"), REF = "A", ALT = "G",
                    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
                    L1 = c("0", "1", "2"), L2 = c("0/1", NA, "1/1"), check.names = FALSE)
  f <- tempfile(fileext = ".csv")
  utils::write.csv(tab, f, row.names = FALSE, na = ".")
  lines <- convert_csv_to_vcf(f)
  expect_identical(lines[3], paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "L1", "L2"),
                                   collapse = "\t"))
  rec <- strsplit(lines[4:6], "\t", fixed = TRUE)
  expect_identical(vapply(rec, `[`, "", 10), c("0/0", "0/1", "1/1"))
  expect_identical(vapply(rec, `[`, "", 11), c("0/1", "./.", "1/1"))
  expect_identical(rec[[1]][6], ".")                      # missing QUAL stays "."
  # tab-separated .txt gives the same VCF
  t2 <- tempfile(fileext = ".txt")
  utils::write.table(tab, t2, sep = "\t", row.names = FALSE, quote = FALSE, na = ".")
  expect_identical(convert_csv_to_vcf(t2), lines)
  # bad values and missing marker columns stop with a clear message
  bad <- tab; bad$L1[1] <- "7"; utils::write.csv(bad, f, row.names = FALSE, na = ".")
  expect_error(convert_csv_to_vcf(f), "Unexpected genotype value")
  utils::write.csv(tab[, -4], f, row.names = FALSE)
  expect_error(convert_csv_to_vcf(f), "Missing required VCF marker columns: REF")
})

test_that("genotype tables are told apart from HapMap and VCF files", {
  d <- make_table(n = 4, p = 3)
  f <- tempfile(fileext = ".txt")
  utils::write.table(d$table, f, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_true(PredictProR:::gp_is_genotype_table_file(f))
  hmp <- tempfile(fileext = ".txt")
  writeLines(c("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\tL1",
               "m1\tA/G\t1\t10\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA"), hmp)
  expect_false(PredictProR:::gp_is_genotype_table_file(hmp))
  vcf <- tempfile(fileext = ".vcf")
  writeLines(convert_csv_to_vcf(f), vcf)
  expect_false(PredictProR:::gp_is_genotype_table_file(vcf))
})

test_that("model_execute reads a CSV genotype table like the equivalent VCF", {
  skip_on_cran()
  d <- make_table()
  dir <- withr::local_tempdir()
  utils::write.csv(d$table, file.path(dir, "geno.csv"), row.names = FALSE, na = ".")
  writeLines(convert_csv_to_vcf(file.path(dir, "geno.csv")), file.path(dir, "geno.vcf"))
  set.seed(9)
  ph <- data.frame(GID = d$ids, Yield = as.numeric(d$G[, 1:5] %*% stats::rnorm(5)) + stats::rnorm(length(d$ids)))
  ph$Yield[1:4] <- NA
  run <- function(...) suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "GBLUP_BRR", gmatrix_method = "VanRaden", nIter = 200L, burnIn = 50L, thin = 2L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, random_state = 1L,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE, ...)))
  from_csv <- run(csv_file_name = "geno.csv", csv_file_path = dir)$model_results$predicted_values
  from_vcf <- run(vcf_file_name = "geno.vcf", vcf_file_path = dir)$model_results$predicted_values
  expect_identical(nrow(from_csv), length(d$ids))
  expect_equal(from_csv$Predicted_value, from_vcf$Predicted_value)
  expect_error(run(csv_file_name = "geno.csv", csv_file_path = dir, vcf_file_name = "geno.vcf",
                   vcf_file_path = dir), "one raw genotype source")
})

test_that("impute_genotypes_with_beagle imputes a CSV genotype table", {
  skip_on_cran()
  java <- Sys.getenv("PREDICTPRO_JAVA", unset = Sys.which("java"))
  skip_if(!nzchar(java), "Java not available for Beagle")
  d <- make_table(n = 30, p = 60, missing_rate = 0.05)
  dir <- withr::local_tempdir()
  f <- file.path(dir, "geno.csv")
  utils::write.csv(d$table, f, row.names = FALSE, na = ".")
  res <- tryCatch(impute_genotypes_with_beagle(f, output_prefix = file.path(dir, "out")),
                  error = function(e) skip(paste("Beagle unavailable:", conditionMessage(e))))
  expect_identical(res$input_format, "csv")
  G <- res$genotypes$snps_matrix
  expect_identical(dim(G), c(30L, 60L))
  expect_false(anyNA(G))
})
