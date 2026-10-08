# BGLR (1.1.0 through at least 1.1.4) initialises ordinal/binary thresholds as
# qnorm(cumsum(table(z)) / n) with n counting NA rows. The top class then gets a
# finite upper bound that is never updated, the latent scale is compressed and
# NA rows inflate the top class. PredictProR passes the lines to predict as NA,
# so every Bayesian classification fit was affected.

simulate_ordinal <- function(n = 300L, p = 60L, missing = 0.3, seed = 3) {
  set.seed(seed)
  X <- matrix(stats::rbinom(n * p, 2, 0.5), n, p)
  Z <- scale(X, scale = FALSE)
  g <- as.vector(Z %*% stats::rnorm(p))
  g <- (g - mean(g)) / stats::sd(g)
  liability <- 2 * g + stats::rnorm(n, sd = 0.7)
  lv <- c("low", "medium", "high")
  truth <- ordered(lv[findInterval(liability, stats::quantile(liability, c(1/3, 2/3))) + 1L], levels = lv)
  y <- truth
  y[sample(n, round(missing * n))] <- NA
  K <- tcrossprod(Z) / mean(diag(tcrossprod(Z)))
  list(Z = Z, K = K, y = y, truth = truth, lv = lv)
}

fit_classification <- function(sim, eta, prefix) {
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  set.seed(1)
  PredictProR:::gp_bayes_bglr_fit(
    fam = "ordinal",
    y = PredictProR:::gp_bayes_prepare_response(sim$y, "ordinal"),
    response_type = "ordinal", ETA = eta,
    nIter = 3000L, burnIn = 1000L, thin = 5L, verbose = FALSE, saveAt = prefix
  )
}

check_heldout <- function(fit, sim) {
  test <- is.na(sim$y)
  expect_identical(dim(fit$probs), c(length(sim$y), length(sim$lv)))
  expect_identical(colnames(fit$probs), sim$lv)
  expect_true(all(is.na(fit$y[test])))
  expect_equal(unname(rowSums(fit$probs)), rep(1, length(sim$y)), tolerance = 1e-6)
  held <- fit$probs[test, , drop = FALSE]
  # classes are balanced; NA rows must not be pushed into the top class
  expect_lt(max(colMeans(held)), 0.45)
  predicted <- max.col(held)
  expect_gte(length(unique(predicted)), 3L)
  truth_code <- as.integer(sim$truth[test])
  expect_gt(stats::cor(as.vector(held %*% seq_along(sim$lv)), truth_code), 0.6)
}

test_that("kernel (RKHS) classification predicts NA lines without BGLR's threshold bias", {
  skip_if_not_installed("BGLR")
  sim <- simulate_ordinal()
  fit <- fit_classification(sim, list(list(K = sim$K, model = "RKHS", saveEffects = TRUE)), "cls_rkhs_")
  check_heldout(fit, sim)
})

test_that("marker (BayesB) classification predicts NA lines without BGLR's threshold bias", {
  skip_if_not_installed("BGLR")
  sim <- simulate_ordinal()
  fit <- fit_classification(sim, list(list(X = sim$Z, model = "BayesB", saveEffects = TRUE)), "cls_bayesb_")
  check_heldout(fit, sim)
})

test_that("classification fits without NA rows are left to BGLR unchanged", {
  skip_if_not_installed("BGLR")
  sim <- simulate_ordinal(missing = 0)
  fit <- fit_classification(sim, list(list(K = sim$K, model = "RKHS")), "cls_full_")
  expect_identical(colnames(fit$probs), sim$lv)
  expect_null(fit$heldout_prediction)
})
