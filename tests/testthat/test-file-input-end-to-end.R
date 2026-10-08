# model_execute() with genotype FILES (not only matrices) must run end to end.
# The VCF/HapMap recoders tag their matrix with an extra "genotype" class, and
# an exact-class check in geno_to_model() used to reject every file input.

make_e2e_vcf <- function(dir, n = 24L, p = 40L) {
  set.seed(11)
  ids <- sprintf("S%02d", seq_len(n))
  dos <- matrix(sample(0:2, n * p, TRUE, prob = c(.45, .35, .2)), p)
  gt <- matrix(c("0/0", "0/1", "1/1")[dos + 1L], p)
  lines <- c(
    "##fileformat=VCFv4.3",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids), collapse = "\t"),
    vapply(seq_len(p), function(j) {
      paste(c("1", j * 1000L, paste0("m", j), "A", "G", ".", "PASS", ".", "GT", gt[j, ]), collapse = "\t")
    }, character(1))
  )
  path <- file.path(dir, "e2e.vcf")
  writeLines(lines, path)
  geno <- t(dos)
  dimnames(geno) <- list(ids, paste0("m", seq_len(p)))
  list(path = path, ids = ids, geno = geno)
}

test_that("model_execute runs GBLUP_BRR from a VCF file and matches the matrix input", {
  skip_if_not_installed("BGLR")
  withr::local_dir(withr::local_tempdir())
  dat <- make_e2e_vcf(getwd())
  ph <- data.frame(GID = dat$ids, Yield = as.numeric(scale(dat$geno[, 1:5] %*% rep(1, 5))) + rnorm(24, sd = 0.5))
  ph$Yield[1:4] <- NA
  run <- function(...) {
    suppressWarnings(PredictProR::model_execute(
      pheno_data = ph, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
      GS_model = "GBLUP_BRR", nIter = 200L, burnIn = 50L, thin = 2L, random_state = 3L,
      qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, gmatrix_method = "VanRaden",
      parallel_mode = "sequential", system_database = TRUE, message = FALSE, ...
    ))
  }
  from_vcf <- run(vcf_file_name = basename(dat$path), vcf_file_path = dirname(dat$path))
  from_matrix <- run(geno_data = dat$geno)
  pv_vcf <- from_vcf$model_results$predicted_values
  pv_mat <- from_matrix$model_results$predicted_values
  expect_equal(nrow(pv_vcf), 24L)
  expect_true(all(is.finite(pv_vcf$Predicted_value)))
  expect_equal(
    pv_vcf$Predicted_value[match(dat$ids, pv_vcf$GID)],
    pv_mat$Predicted_value[match(dat$ids, pv_mat$GID)],
    tolerance = 1e-8
  )
})

test_that("recoded VCF data is QC-filtered once and LD-pruned at most once", {
  withr::local_dir(withr::local_tempdir())
  dat <- make_e2e_vcf(getwd())
  ph <- data.frame(GID = dat$ids, Yield = rnorm(24))
  seen <- list()
  local_mocked_bindings(
    vcf_qc_recode = function(...) {
      m <- dat$geno
      class(m) <- c("matrix", "array", "genotype")
      list(snps_matrix = m, qc_metrics_and_summary_stat = data.frame())
    },
    geno_to_model = function(geno_data = NULL, qc_filtering = NULL, ld_prunning_qc = NULL, ...) {
      seen[[length(seen) + 1L]] <<- list(qc_filtering = qc_filtering, ld_prunning_qc = ld_prunning_qc)
      stop("captured")
    },
    .package = "PredictProR"
  )
  run <- function(ld_pruning, pattern) {
    expect_error(suppressWarnings(PredictProR::model_execute(
      pheno_data = ph, vcf_file_name = basename(dat$path), vcf_file_path = dirname(dat$path),
      gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID, GS_model = "GBLUP_BRR",
      qc_filtering = TRUE, ld_pruning = ld_pruning, ld_prunning_qc = TRUE,
      gmatrix_method = "VanRaden", system_database = TRUE, message = FALSE
    )), pattern)
  }
  # The recoder filtered the file once; the matrix stage must not filter again
  # and is the single LD-pruning step.
  run(FALSE, "captured")
  expect_identical(seen[[1L]], list(qc_filtering = FALSE, ld_prunning_qc = TRUE))
  # Recoder (PLINK) LD pruning is not run by model_execute(): fail clearly
  # instead of pruning twice or naming an argument the user cannot set.
  run(TRUE, "needs the PLINK QC engine")
  expect_length(seen, 1L)
})

test_that("geno_to_model accepts a QC-passed matrix carrying the recoder's genotype class", {
  m <- matrix(c(0, 1, 2, 1, 0, 2, 2, 1, 0, 1, 1, 0), 4,
              dimnames = list(paste0("g", 1:4), paste0("m", 1:3)))
  class(m) <- c("matrix", "array", "genotype")
  out <- PredictProR::geno_to_model(geno_data = m, qc_filtering = FALSE, impute = FALSE,
                                    ld_prunning_qc = FALSE, message = FALSE)
  expect_true(is.matrix(out[[1]]))
  expect_false(inherits(out[[1]], "genotype"))
  expect_identical(attr(out[[1]], "cleared"), "for_model_fit")
})
