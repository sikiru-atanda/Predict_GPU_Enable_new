test_that("C++ spectral repair floors indefinite kernels", {
  skip_if_not(PredictProR:::gp_kernel_cpp_repair_available())

  K <- matrix(c(1, 1.8, 1.8, 1), nrow = 2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  res <- PredictProR:::gp_kernel_spd_repair_cpp(K, min_eigen = 0.01, keep_diag = TRUE)

  expect_identical(res$method, "nearPD_cpp")
  expect_gt(res$adjusted_count, 0L)
  expect_equal(unname(diag(res$matrix)), c(1, 1), tolerance = 1e-8)
  expect_true(PredictProR:::matrix_diagonistic_check(res$matrix, "is_positive_definite"))
  expect_identical(rownames(res$matrix), rownames(K))
})

test_that("kernel precheck can use C++ nearPD repair explicitly", {
  skip_if_not(PredictProR:::gp_kernel_cpp_repair_available())

  K <- matrix(c(1, 1.8, 1.8, 1), nrow = 2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  res <- PredictProR:::grm_kernel_precheck(
    K,
    bending = TRUE,
    bend_value = 0.01,
    duplicate_scan = "none",
    kernel_check_level = "full",
    kernel_fix_method = "nearPD_cpp",
    kernel_sanitize = "none",
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_identical(qc$pd_fix_method, "nearPD_cpp")
  expect_true(PredictProR:::matrix_diagonistic_check(res, "is_positive_definite"))
  expect_true(is.list(qc$pd_fix_cpp))
  expect_gt(qc$pd_fix_cpp$adjusted_count, 0L)
})

test_that("auto speed policy prefers C++ repair after ridge fails", {
  skip_if_not(PredictProR:::gp_kernel_cpp_repair_available())

  K <- matrix(c(1, 10, 10, 1), nrow = 2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  res <- PredictProR:::grm_kernel_precheck(
    K,
    bending = TRUE,
    bend_value = 0.01,
    duplicate_scan = "none",
    kernel_check_level = "full",
    kernel_fix_method = "auto",
    kernel_repair_priority = "speed",
    kernel_sanitize = "none",
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_identical(qc$pd_fix_method, "nearPD_cpp")
  expect_identical(qc$pd_fix_repair_priority, "speed")
  expect_true(is.list(qc$pd_fix_cpp))
  expect_true(PredictProR:::matrix_diagonistic_check(res, "is_positive_definite"))
})

test_that("auto repair priority can prefer C++ repair for structure", {
  skip_if_not(PredictProR:::gp_kernel_cpp_repair_available())

  K <- matrix(c(1, 10, 10, 1), nrow = 2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  res <- PredictProR:::grm_kernel_precheck(
    K,
    bending = TRUE,
    bend_value = 0.01,
    duplicate_scan = "none",
    kernel_check_level = "full",
    kernel_fix_method = "auto",
    kernel_repair_priority = "structure",
    kernel_sanitize = "none",
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_identical(qc$pd_fix_method, "nearPD_cpp")
  expect_identical(qc$pd_fix_repair_priority, "structure")
  expect_true(is.list(qc$pd_fix_cpp))
  expect_true(PredictProR:::matrix_diagonistic_check(res, "is_positive_definite"))
})

test_that("auto keeps Matrix nearPD fallback when C++ repair is outside limit", {
  K <- matrix(c(1, 10, 10, 1), nrow = 2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  res <- PredictProR:::grm_kernel_precheck(
    K,
    bending = TRUE,
    bend_value = 0.01,
    duplicate_scan = "none",
    kernel_check_level = "full",
    kernel_fix_method = "auto",
    kernel_repair_priority = "speed",
    kernel_cpp_repair_size_limit = 1L,
    kernel_nearpd_size_limit = 10L,
    kernel_sanitize = "none",
    show_message = FALSE
  )

  qc <- attr(res, "kernel_qc")
  expect_identical(qc$pd_fix_method, "nearPD")
  expect_identical(qc$pd_fix_repair_priority, "speed")
  expect_true(PredictProR:::matrix_diagonistic_check(res, "is_positive_definite"))
})
