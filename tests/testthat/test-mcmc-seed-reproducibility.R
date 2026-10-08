test_that("gp_with_pinned_seed is reproducible and restores the caller's RNG", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  before <- .Random.seed

  a <- PredictProR:::gp_with_pinned_seed(7L, stats::rnorm(5))
  b <- PredictProR:::gp_with_pinned_seed(7L, stats::rnorm(5))
  expect_identical(a, b)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, before)

  # NULL seed leaves the caller's stream in charge
  set.seed(3); x <- PredictProR:::gp_with_pinned_seed(NULL, stats::runif(1))
  set.seed(3); y <- stats::runif(1)
  expect_identical(x, y)
})

test_that("joint multi-trait Bayesian fit is reproducible for a fixed random_state", {
  skip_if_not_installed("BGLR")
  ids <- paste0("g", seq_len(6))
  pheno <- data.frame(
    GID = ids,
    T1 = c(1.1, 1.8, NA, 2.7, 3.0, 3.4),
    T2 = c(2.2, NA, 2.8, 3.3, 3.8, 4.1),
    stringsAsFactors = FALSE
  )
  K <- diag(seq(0.7, 1.3, length.out = length(ids)))
  rownames(K) <- colnames(K) <- ids
  fit <- function(seed) {
    out <- PredictProR:::bayes_multitrait_joint_fit(
      pheno_data = pheno,
      response = c("T1", "T2"),
      gen_name = "GID",
      kernels = list(gmatrix = K),
      GS_model = "RKHS",
      bayes_para = list(nIter = 45L, burnIn = 15L, thin = 2L),
      verbose = FALSE,
      random_state = seed
    )
    out$bayes_result$Predicted_value$Predicted_value
  }
  first <- fit(11L)
  expect_identical(fit(11L), first)
  expect_false(identical(fit(12L), first))
})
