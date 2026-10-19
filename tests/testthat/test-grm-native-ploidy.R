# The native GRM backend computed every relationship with its own triple loop
# (143 s for 1,000 lines x 10,000 markers where R's BLAS took 9 s) and refused
# polyploid dosage. It now builds the marker covariates in one pass and forms
# the product with BLAS dsyrk, for every ploidy the method allows, in the
# order of operations of the R backend: the two give the same numbers, so a
# run's results do not depend on which one built the matrix.

.dosage <- function(n, p, ploidy, seed) {
  set.seed(seed)
  freq <- stats::runif(p, 0.1, 0.9)
  x <- matrix(stats::rbinom(n * p, ploidy, rep(freq, each = n)), n, p,
              dimnames = list(sprintf("L%02d", seq_len(n)), sprintf("m%03d", seq_len(p))))
  storage.mode(x) <- "double"
  attr(x, "ploidy") <- ploidy
  x
}

test_that("the native GRM matches the R backend for every ploidy and method", {
  for (ploidy in c(2L, 4L, 6L)) {
    x <- .dosage(30, 120, ploidy, seed = ploidy)
    w <- matrix(stats::runif(ncol(x), 0, 2), ncol = 1, dimnames = list(colnames(x), NULL))
    methods <- c("VanRaden", "Weighted_VanRaden", "Epistasis")
    if (ploidy == 2L) methods <- c(methods, "Yang", "Dominance_Vitezica", "Dominance_Su")
    for (m in methods) {
      wt <- if (m == "Weighted_VanRaden") w else NULL
      r <- grm_calculation(x, method = m, weight = wt, backend = "r", ploidy = ploidy)
      cpp <- grm_calculation(x, method = m, weight = wt, backend = "cpp", ploidy = ploidy)
      expect_identical(unclass(cpp), unclass(r), info = paste("ploidy", ploidy, m))
      expect_identical(attr(cpp, "ploidy"), ploidy, info = paste("ploidy", ploidy, m))
      expect_true(isSymmetric(unname(unclass(cpp))), info = paste("ploidy", ploidy, m))
    }
  }
})

test_that("the native GRM drops uninformative markers as the R backend does", {
  x <- .dosage(25, 60, 4L, seed = 9)
  x[, 3] <- 2           # constant
  x[, 7] <- 0           # constant
  r <- grm_calculation(x, method = "VanRaden", backend = "r", ploidy = 4L)
  cpp <- grm_calculation(x, method = "VanRaden", backend = "cpp", ploidy = 4L)
  expect_identical(unclass(cpp), unclass(r))
})

test_that("the native VanRaden GRM is the matrix earlier versions built in R", {
  x <- .dosage(40, 200, 4L, seed = 11)
  freq <- colMeans(x) / 4L
  before <- tcrossprod(scale(x, center = TRUE, scale = FALSE)) / sum(4L * freq * (1 - freq))
  cpp <- grm_calculation(x, method = "VanRaden", backend = "cpp", ploidy = 4L)
  expect_identical(unname(unclass(cpp))[, ], unname(before)[, ])
  expect_identical(unname(unclass(grm_calculation(x, method = "Epistasis", backend = "cpp", ploidy = 4L)))[, ],
                   unname(before * before)[, ])
})

test_that("dominance and Yang stay diploid-only on the native backend too", {
  x <- .dosage(20, 50, 4L, seed = 4)
  for (m in c("Yang", "Dominance_Vitezica", "Dominance_Su")) {
    expect_error(grm_calculation(x, method = m, backend = "cpp", ploidy = 4L), "diploid-only")
  }
})
