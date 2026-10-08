test_that("native genotype QC metrics match existing R semantics", {
  ns <- asNamespace("PredictProR")
  expect_true(exists("geno_qc_cpp_available", envir = ns, inherits = FALSE))
  skip_if_not(PredictProR:::geno_qc_cpp_available())

  geno <- matrix(
    c(
      0, 0, 1, NA, 2, 0,
      0, 1, 1, NA, 2, 1,
      0, 2, 1, NA, NA, 2,
      0, NA, 1, NA, 0, 2,
      0, 1, 1, NA, 0, NA
    ),
    nrow = 5L,
    byrow = TRUE
  )
  colnames(geno) <- paste0("m", seq_len(ncol(geno)))
  rownames(geno) <- paste0("g", seq_len(nrow(geno)))

  fast <- PredictProR:::geno_qc_metrics_cpp(geno)

  ref_monomorphic <- which(apply(geno, 2L, function(x) length(table(x)) <= 1L))
  ref_marker_missing <- colMeans(is.na(geno))
  ref_ind_missing <- rowMeans(is.na(geno))
  phat <- colMeans(geno, na.rm = TRUE) / 2
  ref_maf <- ifelse(phat < 0.5, phat, 1 - phat)
  ref_het <- apply(geno, 2L, function(x) sum(x == 1, na.rm = TRUE) / length(x))

  expect_identical(fast$monomorphic, as.integer(ref_monomorphic))
  expect_equal(fast$marker_missing_rate, ref_marker_missing)
  expect_equal(fast$individual_missing_rate, ref_ind_missing)
  expect_equal(fast$maf, ref_maf)
  expect_equal(fast$heterozygosity, ref_het)
})

test_that("native mean imputation preserves matrix shape and observed values", {
  ns <- asNamespace("PredictProR")
  expect_true(exists("geno_impute_summary_cpp_available", envir = ns, inherits = FALSE))
  skip_if_not(PredictProR:::geno_impute_summary_cpp_available())

  geno <- matrix(
    c(
      0, NA, 2,
      1, 1, NA,
      NA, 2, 2,
      2, NA, 0
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:4), paste0("m", 1:3))
  )

  fast <- PredictProR:::geno_impute_summary_cpp(geno, method = "mean")
  ref <- geno
  for (j in seq_len(ncol(ref))) {
    tmp <- ref[, j]
    ref[, j] <- ifelse(is.na(tmp), round(mean(tmp, na.rm = TRUE)), tmp)
  }

  expect_true(is.matrix(fast))
  expect_identical(dimnames(fast), dimnames(geno))
  expect_equal(fast, ref)
  expect_equal(fast[!is.na(geno)], geno[!is.na(geno)])
})

test_that("native KNN imputation fills missing genotypes and keeps observed calls", {
  ns <- asNamespace("PredictProR")
  expect_true(exists("geno_impute_knn_cpp_available", envir = ns, inherits = FALSE))
  skip_if_not(PredictProR:::geno_impute_knn_cpp_available())

  geno <- matrix(
    c(
      0, 0, 2, 2,
      0, NA, 2, 2,
      2, 2, 0, NA,
      2, 2, 0, 0
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:4), paste0("m", 1:4))
  )

  fast <- PredictProR:::geno_impute_knn_cpp(geno, k = 1L)

  expect_false(anyNA(fast))
  expect_identical(dimnames(fast), dimnames(geno))
  expect_equal(fast[!is.na(geno)], geno[!is.na(geno)])
  expect_equal(fast[2L, 2L], 0)
  expect_equal(fast[3L, 4L], 0)
})

test_that("native KNN imputation exposes a donor cap for large data", {
  ns <- asNamespace("PredictProR")
  expect_true(exists("geno_impute_knn_cpp_available", envir = ns, inherits = FALSE))
  skip_if_not(PredictProR:::geno_impute_knn_cpp_available())
  expect_true("max_donors" %in% names(formals(PredictProR:::geno_impute_knn_cpp)))

  set.seed(20260512)
  geno <- matrix(
    sample(c(0, 1, 2, NA), 80L * 40L, replace = TRUE, prob = c(0.3, 0.3, 0.35, 0.05)),
    nrow = 80L
  )
  rownames(geno) <- paste0("g", seq_len(nrow(geno)))
  colnames(geno) <- paste0("m", seq_len(ncol(geno)))

  fast <- PredictProR:::geno_impute_knn_cpp(geno, k = 3L, max_donors = 16L)

  expect_false(anyNA(fast))
  expect_equal(fast[!is.na(geno)], geno[!is.na(geno)])
})
