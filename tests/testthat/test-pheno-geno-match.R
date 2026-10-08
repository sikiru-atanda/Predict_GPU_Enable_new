make_match_pheno <- function(ids = c("g1", "g2", "g3")) {
  data.frame(
    GID = ids,
    Yield = seq_along(ids),
    stringsAsFactors = FALSE
  )
}

make_match_matrix <- function(ids = c("g1", "g2", "g3", "g4")) {
  out <- outer(seq_along(ids), seq_len(3), function(i, j) (i + j) %% 3)
  dimnames(out) <- list(ids, paste0("m", 1:3))
  storage.mode(out) <- "double"
  out
}

test_that("pheno_geno_match treats genomic-only individuals as test set", {
  pheno <- make_match_pheno(c("g1", "g2", "g3"))
  geno <- make_match_matrix(c("g1", "g2", "g3", "g4", "g5"))

  out <- suppressMessages(PredictProR:::pheno_geno_match(
    object_geno = geno,
    object_pheno = pheno,
    gen_name = "GID",
    message = FALSE
  ))

  expect_identical(out$test_data, c("g4", "g5"))
  expect_identical(rownames(out$geno_pheno_match_data), c("g1", "g2", "g3", "g4", "g5"))
})

test_that("pheno_geno_match preserves explicit test set when all IDs have phenotype", {
  pheno <- make_match_pheno(c("g1", "g2", "g3"))
  geno <- make_match_matrix(c("g1", "g2", "g3"))

  out <- PredictProR:::pheno_geno_match(
    object_geno = geno,
    object_pheno = pheno,
    gen_name = "GID",
    test_set = "g2",
    message = FALSE
  )

  expect_identical(out$test_data, "g2")
  expect_identical(rownames(out$geno_pheno_match_data), c("g1", "g2", "g3"))
})

test_that("pheno_geno_match combines explicit and genomic-only test IDs", {
  pheno <- make_match_pheno(c("g1", "g2", "g3"))
  geno <- make_match_matrix(c("g1", "g2", "g3", "g4"))

  out <- suppressMessages(PredictProR:::pheno_geno_match(
    object_geno = geno,
    object_pheno = pheno,
    gen_name = "GID",
    test_set = "g2",
    message = FALSE
  ))

  expect_identical(out$test_data, c("g2", "g4"))
  expect_identical(rownames(out$geno_pheno_match_data), c("g1", "g2", "g3", "g4"))
})

test_that("pheno_geno_match lets explicit test set override genomic-only inference", {
  pheno <- make_match_pheno(c("g1", "g2", "g3"))
  geno <- make_match_matrix(c("g1", "g2", "g3", "g4", "g5"))

  out <- suppressMessages(PredictProR:::pheno_geno_match(
    object_geno = geno,
    object_pheno = pheno,
    gen_name = "GID",
    test_set = "g2",
    test_set_overrides_inferred = TRUE,
    message = FALSE
  ))

  expect_identical(out$test_data, "g2")
  expect_identical(rownames(out$geno_pheno_match_data), c("g1", "g2", "g3"))
})

test_that("pheno_geno_match derives test set from train set", {
  pheno <- make_match_pheno(c("g1", "g2", "g3"))
  geno <- make_match_matrix(c("g1", "g2", "g3"))

  out <- PredictProR:::pheno_geno_match(
    object_geno = geno,
    object_pheno = pheno,
    gen_name = "GID",
    train_set = c("g1", "g3"),
    message = FALSE
  )

  expect_identical(out$test_data, "g2")
})

test_that("process_geno_data propagates genomic-only test IDs", {
  pheno_clean <- list(pheno_clean_data = make_match_pheno(c("g1", "g2", "g3")))
  geno <- make_match_matrix(c("g1", "g2", "g3", "g4"))

  out <- suppressMessages(PredictProR::process_geno_data(
    geno_data = geno,
    pheno_clean_list = pheno_clean,
    gen_name = "GID",
    maf_threshold = 0.01,
    het_threshold = 0.95,
    ind_call_rate_threshold = 0.95,
    snp_call_rate_threshold = 0.95,
    impute = FALSE,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    message = FALSE
  ))

  expect_identical(out$test_set, "g4")
  expect_identical(rownames(out$geno_model_ready), c("g1", "g2", "g3", "g4"))
})

test_that("process_geno_data honors explicit test set override for genomic-only IDs", {
  pheno_clean <- list(
    pheno_clean_data = make_match_pheno(c("g1", "g2", "g3")),
    test_set = "g2",
    test_set_source = "explicit"
  )
  geno <- make_match_matrix(c("g1", "g2", "g3", "g4"))

  out <- suppressMessages(PredictProR::process_geno_data(
    geno_data = geno,
    pheno_clean_list = pheno_clean,
    gen_name = "GID",
    test_set = "g2",
    maf_threshold = 0.01,
    het_threshold = 0.95,
    ind_call_rate_threshold = 0.95,
    snp_call_rate_threshold = 0.95,
    impute = FALSE,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    test_set_overrides_inferred = TRUE,
    message = FALSE
  ))

  expect_identical(out$test_set, "g2")
  expect_identical(rownames(out$geno_model_ready), c("g1", "g2", "g3"))
})

test_that("process_geno_data removes post-QC failed individuals from phenotype matching", {
  pheno_clean <- list(pheno_clean_data = make_match_pheno(c("g1", "g2", "g3", "g4")))
  geno <- matrix(
    c(
      0, 1, 2, 0,
      NA, NA, 2, 1,
      2, 1, 0, 2,
      1, 0, 1, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3", "g4"), paste0("m", 1:4))
  )

  out <- suppressMessages(PredictProR::process_geno_data(
    geno_data = geno,
    pheno_clean_list = pheno_clean,
    gen_name = "GID",
    test_set = c("g2", "g4"),
    maf_threshold = 0,
    het_threshold = 1,
    ind_call_rate_threshold = 0.75,
    snp_call_rate_threshold = 0.75,
    impute = FALSE,
    qc_filtering = TRUE,
    ld_prunning_qc = FALSE,
    message = FALSE
  ))

  expect_identical(out$low_call_rate_inds_removed, "g2")
  expect_identical(out$test_set, "g4")
  expect_false("g2" %in% rownames(out$geno_model_ready))
  expect_identical(rownames(out$geno_model_ready), c("g1", "g3", "g4"))
})
