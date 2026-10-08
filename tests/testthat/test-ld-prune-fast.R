test_that("native LD graph pruning matches R reference on representative data", {
  skip_if_not(PredictProR:::ld_prune_cpp_available())

  set.seed(20260507)
  geno <- matrix(sample(c(0, 1, 2, NA), 36L * 10L, replace = TRUE, prob = c(0.3, 0.35, 0.3, 0.05)),
                 nrow = 36L)
  geno[, 2] <- geno[, 1]
  geno[, 5] <- geno[, 4]
  geno[, 9] <- 0
  geno[1:2, 9] <- 1
  colnames(geno) <- paste0("m", seq_len(ncol(geno)))

  ref <- suppressMessages(PredictProR:::ld_prune_graph_r(
    geno,
    window = 3L,
    R2_threshold = 0.8,
    maf_thresh = 0.03
  ))
  fast <- suppressMessages(PredictProR:::ld_prune_graph_cpp(
    geno,
    window = 3L,
    R2_threshold = 0.8,
    maf_thresh = 0.03
  ))

  expect_identical(fast, ref)
  expect_false("m9" %in% fast)
})

test_that("LD graph pruning keeps public API and supports unnamed integer matrices", {
  skip_if_not(PredictProR:::ld_prune_cpp_available())

  old_backend <- Sys.getenv("PREDICTPRO_LD_PRUNE_BACKEND", unset = NA_character_)
  on.exit({
    if (is.na(old_backend)) {
      Sys.unsetenv("PREDICTPRO_LD_PRUNE_BACKEND")
    } else {
      Sys.setenv(PREDICTPRO_LD_PRUNE_BACKEND = old_backend)
    }
  }, add = TRUE)
  Sys.setenv(PREDICTPRO_LD_PRUNE_BACKEND = "cpp")

  geno <- matrix(
    as.integer(c(
      0, 1, 2, 0,
      0, 1, 2, 0,
      2, 1, 0, 2,
      2, 1, 0, 2,
      0, 0, 0, 0
    )),
    nrow = 5L,
    byrow = TRUE
  )
  out <- suppressMessages(PredictProR::ld_prune_graph(
    geno,
    window = 2L,
    R2_threshold = 0.5,
    maf_thresh = 0.01
  ))

  expect_type(out, "integer")
  expect_true(length(out) >= 1L)
})

test_that("native polyploid LD pruning matches dosage-correlation reference", {
  skip_if_not(PredictProR:::ld_prune_cpp_available())

  set.seed(20260815)
  geno <- matrix(sample(c(0:4, NA), 48L * 8L, replace = TRUE), nrow = 48L)
  geno[, 2L] <- geno[, 1L]
  geno[, 8L] <- 0
  colnames(geno) <- paste0("p4_m", seq_len(ncol(geno)))
  attr(geno, "ploidy") <- 4L

  ref <- suppressMessages(PredictProR:::ld_prune_graph_r(
    geno, window = 3L, R2_threshold = 0.8, maf_thresh = 0.03,
    ploidy = 4L
  ))
  fast <- suppressMessages(PredictProR:::ld_prune_graph_cpp(
    geno, window = 3L, R2_threshold = 0.8, maf_thresh = 0.03,
    ploidy = 4L
  ))
  expect_identical(fast, ref)
  expect_false("p4_m8" %in% fast)

  pair <- ld_pair_vec(geno[, 1L], geno[, 3L], ploidy = 4L)
  expected_r2 <- stats::cor(
    geno[, 1L], geno[, 3L], use = "pairwise.complete.obs"
  )^2
  expect_equal(unname(pair[[1L]]), unname(expected_r2), tolerance = 1e-12)
  expect_true(is.na(pair[[2L]]))
})
