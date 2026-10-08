kernel_collision_fixture <- function(n = 5L) {
  ids <- paste0("g", seq_len(n))
  kernels <- lapply(seq_len(4L), function(i) {
    K <- outer(seq_len(n), seq_len(n), function(a, b) (i / 5)^abs(a - b))
    dimnames(K) <- list(ids, ids)
    K
  })
  names(kernels) <- c("genomic", "RNA-A", "RNA_A", "RNA_A_1")
  kernels
}

test_that("prepared kernel banks retain matrices whose names sanitize identically", {
  original <- kernel_collision_fixture()
  for (key in c("gmatrix_model_ready", "omic1_kernel_model_ready", "kernel_list_model_ready")) {
    prepared <- gp_add_kernel_inputs(list(), key, original)
    expect_length(prepared, length(original))
    expect_identical(anyDuplicated(names(prepared)), 0L)
    expect_identical(
      vapply(prepared, gp_kernel_source_name, character(1L), fallback = "missing"),
      stats::setNames(names(original), names(prepared))
    )
    expect_identical(lapply(unname(prepared), as.vector), lapply(unname(original), as.vector))
  }
})

test_that("kernel collection retains repeated names and source tags in input order", {
  original <- kernel_collision_fixture()
  tagged <- lapply(original, function(K) {
    attr(K, "predictpror_source_kernel") <- "shared"
    K
  })
  bank <- gp_collect_kernel_inputs(
    gmatrix = tagged[[1L]],
    gkernel = list(shared = tagged[[2L]], shared = tagged[[3L]]),
    kernel_list = list(shared = original[[4L]], shared = original[[1L]])
  )
  expect_length(bank, 5L)
  expect_identical(anyDuplicated(names(bank)), 0L)
  expected <- original[c(1L, 2L, 3L, 4L, 1L)]
  expect_identical(lapply(unname(bank), as.vector), lapply(unname(expected), as.vector))

  unnamed <- gp_collect_kernel_inputs(kernel_list = unname(original))
  expect_length(unnamed, length(original))
  expect_identical(lapply(unname(unnamed), as.vector), lapply(unname(original), as.vector))
})

test_that("eigenfeatures reconstruct every kernel despite colliding source names", {
  original <- kernel_collision_fixture()
  prepared <- gp_add_kernel_inputs(list(), "kernel_list_model_ready", original)
  for (bank in list(original, prepared)) {
    features <- gp_prepare_met_kernel_features(kernel_list = bank, var_explained = 1)
    expect_length(features$feature_blocks, length(original))
    expect_identical(anyDuplicated(features$feature_summary$source), 0L)
    expect_identical(anyDuplicated(colnames(features$feature_table)), 0L)
    for (i in seq_along(features$feature_blocks)) {
      expect_equal(
        tcrossprod(as.matrix(features$feature_blocks[[i]])),
        original[[i]], tolerance = 1e-10
      )
    }
  }
})

test_that("eigenfeature labels cannot hide explicit kernels or duplicate feature columns", {
  original <- kernel_collision_fixture()
  features <- gp_prepare_met_kernel_features(
    gmatrix = original[[1L]],
    omic1_kernel = original[[2L]],
    omic2_kernel = original[[3L]],
    kernel_list = list(GRM = original[[4L]]),
    omics_kernel_label = list(omic1_kernel = "GRM", omic2_kernel = "GRM"),
    var_explained = 1
  )
  expect_length(features$feature_blocks, 4L)
  expect_identical(anyDuplicated(colnames(features$feature_table)), 0L)
  for (i in seq_along(features$feature_blocks)) {
    expect_equal(tcrossprod(as.matrix(features$feature_blocks[[i]])), original[[i]], tolerance = 1e-10)
  }

  # Prepared canonical entries mirror explicit arguments and must not be counted twice.
  prepared <- gp_add_kernel_inputs(list(), "gmatrix_model_ready", original[[1L]])
  mirrored <- gp_prepare_met_kernel_features(gmatrix = original[[1L]], kernel_list = prepared)
  expect_length(mirrored$feature_blocks, 1L)
})
