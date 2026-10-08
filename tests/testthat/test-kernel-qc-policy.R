test_that("block duplicate diagnostics find high-correlation kernel pairs", {
  K <- diag(4)
  rownames(K) <- colnames(K) <- paste0("g", 1:4)
  K[1, 2] <- K[2, 1] <- 0.98
  K[3, 4] <- K[4, 3] <- 0.97

  res <- PredictProR:::grm_kernel_diagnostic_fix(
    K,
    duplicate_cut_off = 0.95,
    duplicate_scan = "full",
    duplicate_block_size = 2
  )

  expect_true("potential_duplicates" %in% names(res))
  expect_equal(nrow(res$potential_duplicates), 2L)
  expect_equal(res$kernel_qc$duplicate_scan, "full")
  expect_true(res$kernel_qc$duplicate_scan_complete)
})

test_that("C++ duplicate scanner matches R block scanner", {
  skip_if_not(PredictProR:::gp_kernel_cpp_duplicate_scan_available())

  K <- diag(6)
  rownames(K) <- colnames(K) <- paste0("g", 1:6)
  K[1, 2] <- K[2, 1] <- 0.98
  K[1, 3] <- K[3, 1] <- 0.97
  K[4, 6] <- K[6, 4] <- 0.99

  r_res <- PredictProR:::gp_kernel_block_duplicate_pairs(
    K,
    threshold = 0.95,
    block_size = 2L,
    max_pairs = Inf,
    row_index = c(6L, 4L, 3L, 2L, 1L)
  )
  cpp_res <- PredictProR:::gp_kernel_block_duplicate_pairs_cpp(
    K,
    threshold = 0.95,
    block_size = 2L,
    max_pairs = Inf,
    row_index = c(6L, 4L, 3L, 2L, 1L)
  )

  expect_equal(cpp_res, r_res, tolerance = 1e-12)
})

test_that("C++ duplicate scanner respects max_pairs cap", {
  skip_if_not(PredictProR:::gp_kernel_cpp_duplicate_scan_available())

  K <- matrix(0.99, nrow = 5, ncol = 5)
  diag(K) <- 1
  rownames(K) <- colnames(K) <- paste0("g", 1:5)

  res <- PredictProR:::gp_kernel_block_duplicate_pairs_cpp(
    K,
    threshold = 0.95,
    block_size = 2L,
    max_pairs = 3L
  )

  expect_equal(nrow(res), 3L)
  expect_equal(names(res), c("Row", "Col", "Corr", "RowName", "ColName"))
})

test_that("auto duplicate scanner uses effective scan size threshold", {
  skip_if_not(PredictProR:::gp_kernel_cpp_duplicate_scan_available())

  old_cpp <- Sys.getenv("PREDICTPRO_KERNEL_DUP_CPP", unset = NA_character_)
  old_min <- Sys.getenv("PREDICTPRO_KERNEL_DUP_CPP_MIN_N", unset = NA_character_)
  on.exit({
    if (is.na(old_cpp)) {
      Sys.unsetenv("PREDICTPRO_KERNEL_DUP_CPP")
    } else {
      Sys.setenv(PREDICTPRO_KERNEL_DUP_CPP = old_cpp)
    }
    if (is.na(old_min)) {
      Sys.unsetenv("PREDICTPRO_KERNEL_DUP_CPP_MIN_N")
    } else {
      Sys.setenv(PREDICTPRO_KERNEL_DUP_CPP_MIN_N = old_min)
    }
  }, add = TRUE)

  Sys.setenv(
    PREDICTPRO_KERNEL_DUP_CPP = "auto",
    PREDICTPRO_KERNEL_DUP_CPP_MIN_N = "3000"
  )

  expect_false(PredictProR:::gp_kernel_use_cpp_duplicate_scan(n = 5000L, selected_n = 2000L))
  expect_true(PredictProR:::gp_kernel_use_cpp_duplicate_scan(n = 3000L, selected_n = 3000L))
  expect_true(PredictProR:::gp_kernel_use_cpp_duplicate_scan(n = 5899L, selected_n = 5899L))
})

test_that("duplicate optimization keeps one representative per duplicate cluster", {
  K <- matrix(0.99, nrow = 3, ncol = 3)
  diag(K) <- c(1.00, 1.02, 0.98)
  rownames(K) <- colnames(K) <- paste0("g", 1:3)

  res <- PredictProR:::grm_kernel_diagnostic_fix(
    K,
    duplicate_cut_off = 0.95,
    optimize_duplicate = TRUE,
    duplicate_scan = "full"
  )

  expect_equal(nrow(res$clean_matrix), 1L)
  expect_identical(rownames(res$clean_matrix), "g1")
  expect_true("duplicate_clusters" %in% names(res))
  expect_equal(length(res$removed_ids), 2L)
})

test_that("large-kernel auto mode uses bounded duplicate sampling", {
  K <- diag(10)
  rownames(K) <- colnames(K) <- paste0("g", seq_len(10))
  K[1, 10] <- K[10, 1] <- 0.99

  res <- PredictProR:::grm_kernel_precheck(
    K,
    kernel_large_n_threshold = 5L,
    duplicate_scan = "auto",
    duplicate_sample_size = 4L,
    duplicate_block_size = 2L,
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_identical(qc$check_level, "light")
  expect_identical(qc$duplicate_scan, "sample")
  expect_identical(qc$pd_check, "sample")
  expect_identical(attr(res, "cleared"), "for_model_fit")
})

test_that("large-kernel auto mode applies bounded prefit sanitizer", {
  K <- diag(10)
  rownames(K) <- colnames(K) <- paste0("g", seq_len(10))
  K[2, 8] <- K[8, 2] <- 0.30

  res <- PredictProR:::grm_kernel_precheck(
    K,
    kernel_large_n_threshold = 5L,
    duplicate_scan = "auto",
    duplicate_sample_size = 4L,
    duplicate_block_size = 2L,
    kernel_sanitize = "auto",
    kernel_sanitize_value = 0.02,
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_true(qc$sanitizer_applied)
  expect_identical(qc$sanitizer_reason, "large_kernel_light_diagnostics")
  expect_equal(res[2, 8], 0.30 * 0.98)
})

test_that("large-kernel sanitizer can be disabled explicitly", {
  K <- diag(10)
  rownames(K) <- colnames(K) <- paste0("g", seq_len(10))
  K[2, 8] <- K[8, 2] <- 0.30

  res <- PredictProR:::grm_kernel_precheck(
    K,
    kernel_large_n_threshold = 5L,
    duplicate_scan = "auto",
    duplicate_sample_size = 4L,
    duplicate_block_size = 2L,
    kernel_sanitize = "none",
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_false(qc$sanitizer_applied)
  expect_equal(res[2, 8], 0.30)
})

test_that("kernel precheck stabilizes near-duplicate kernels with ridge blending", {
  K <- diag(4)
  rownames(K) <- colnames(K) <- paste0("g", 1:4)
  K[1, 2] <- K[2, 1] <- 0.999999

  res <- PredictProR:::grm_kernel_precheck(
    K,
    bending = TRUE,
    bend_value = 0.01,
    blending_value = 0.05,
    duplicate_cut_off = 0.95,
    kernel_check_level = "full",
    duplicate_scan = "full",
    kernel_fix_method = "ridge",
    show_message = FALSE
  )

  expect_true(PredictProR:::matrix_diagonistic_check(res, "is_positive_definite"))
  expect_lt(res[1, 2] / sqrt(res[1, 1] * res[2, 2]), 0.95)
})
