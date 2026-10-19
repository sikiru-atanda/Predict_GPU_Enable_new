# Part A: joint posterior-mean feature ranking and visible-record input contract.
.joint_fixture <- function(n = 90L, p = 12L) {
  gp_with_pinned_seed(20261014L, {
    X <- matrix(rnorm(n * p), n, p,
                dimnames = list(paste0("L", seq_len(n)), paste0("m", seq_len(p))))
    ph <- data.frame(ID = rownames(X), Env = rep(c("E1", "E2", "E3"), length.out = n),
                     T1 = 4 * X[, 1] + rnorm(n, sd = 0.4),
                     T2 = 3 * X[, 1] + rnorm(n, sd = 0.4))
    list(X = X, ph = ph)
  })
}
.joint_score <- function(d, scoring_model = "Ridge_Regression", ...) {
  feature_score_predictors(d$X, d$ph, c("T1", "T2"), "ID",
                          scoring_model = scoring_model, joint = TRUE,
                          bayes_nIter = 300L, bayes_burnIn = 150L,
                          bayes_thin = 3L, seed = 17L, ...)
}

test_that("joint BRR and SpikeSlab return one row per predictor and shared signal", {
  d <- .joint_fixture()
  for (method in c("Ridge_Regression", "BayesB")) {
    m <- .joint_score(d, method)
    expect_equal(nrow(m), ncol(d$X))
    expect_identical(unique(m$scoring_scope), "joint")
    expect_identical(unique(m$trait), "T1+T2")
    expect_identical(unique(m$scoring_engine),
                     paste("BGLR Multitrait", if (method == "BayesB") "SpikeSlab" else "BRR"))
    expect_identical(m$predictor[1L], "m1")
    expect_true(all(is.finite(m$raw_score)))
    expect_true(all(is.finite(m$effect_sd)))
    expect_identical(unique(m$retained_draws), 50L)
  }
})

test_that("joint BRR score and coefficient spread equal a direct same-prior fit", {
  d <- .joint_fixture()
  z <- gp_feature_joint_inputs(d$X, d$ph, c("T1", "T2"), "ID")
  prior <- gp_feature_joint_prior_setup(z$x, z$y, z$C)
  prefix <- tempfile("direct_feature_joint_")
  on.exit(unlink(Sys.glob(paste0(prefix, "*"))), add = TRUE)
  invisible(capture.output(
    direct <- gp_with_pinned_seed(17L, BGLR::Multitrait(
      y = z$y, ETA = list(markers = list(X = z$x, model = "BRR", Cov = prior$marker)),
      resCov = prior$residual, nIter = 300L, burnIn = 150L, thin = 3L,
      saveAt = prefix, verbose = FALSE
    ))
  ))
  m <- .joint_score(d)
  idx <- match(colnames(d$X), m$predictor)
  expect_equal(m$raw_score[idx], sqrt(rowSums(direct$ETA$markers$beta^2)), tolerance = 1e-12)
  expect_equal(m$effect_sd[idx], sqrt(rowSums(direct$ETA$markers$SD.beta^2)), tolerance = 1e-12)
})

test_that("joint prior calibration has independent analytic intercept and MET references", {
  d <- .joint_fixture()
  z <- gp_feature_joint_inputs(d$X, d$ph, c("T1", "T2"), "ID")
  prior <- gp_feature_joint_prior_setup(z$x, z$y, z$C)
  S <- cov(z$y)
  ms <- sum(colSums(z$x^2)) / nrow(z$x) - sum(colMeans(z$x)^2)
  expect_equal(prior$marker$S0, 0.5 * S / ms * 6, tolerance = 1e-12)
  expect_equal(prior$residual$S0, 0.5 * S * 8, tolerance = 1e-12)
  expect_equal(prior$shrinkage, 0)
  z <- gp_feature_joint_inputs(d$X, d$ph, c("T1", "T2"), "ID", "Env")
  prior <- gp_feature_joint_prior_setup(z$x, z$y, z$C)
  yr <- vapply(seq_len(ncol(z$y)), function(j) lm.fit(z$C, z$y[, j])$residuals,
               numeric(nrow(z$y)))
  colnames(yr) <- colnames(z$y)
  ms <- sum(colSums(z$x^2)) / nrow(z$x) - sum(colMeans(z$x)^2)
  expect_equal(prior$marker$S0, 0.5 * cov(yr) / ms * 6, tolerance = 1e-12)
  shifted <- z$y + z$C %*% matrix(seq_len(ncol(z$C) * 2), ncol(z$C), 2)
  expect_equal(gp_feature_joint_prior_setup(z$x, shifted, z$C)$covariance,
               prior$covariance, tolerance = 1e-12)
})

test_that("joint inputs scale each trait on its visible records and retain partial records", {
  d <- .joint_fixture()
  d$ph$T1[1:10] <- NA_real_
  d$ph$T2[11:20] <- NA_real_
  z <- gp_feature_joint_inputs(d$X, d$ph, c("T1", "T2"), "ID")
  expect_equal(z$n_train, 90L)
  expect_equal(z$min_pair_coobserved, 70L)
  expect_equal(z$trait_sd["T1"], c(T1 = sd(d$ph$T1, na.rm = TRUE)))
  expect_equal(z$trait_sd["T2"], c(T2 = sd(d$ph$T2, na.rm = TRUE)))
  m <- .joint_score(d)
  expect_identical(unique(m$n_train), 90L)
})

test_that("joint BGLR runs without a co-observed pair or a complete three-trait record", {
  d <- .joint_fixture()
  d$ph$T1[1:45] <- NA_real_
  d$ph$T2[46:90] <- NA_real_
  m <- .joint_score(d)
  expect_identical(unique(m$min_pair_coobserved), 0L)
  expect_identical(unique(m$n_train), 90L)
  d <- .joint_fixture()
  d$ph$T3 <- gp_with_pinned_seed(11L, rnorm(90))
  d$ph$T1[1:30] <- NA_real_
  d$ph$T2[31:60] <- NA_real_
  d$ph$T3[61:90] <- NA_real_
  m <- suppressWarnings(feature_score_predictors(
    d$X, d$ph, c("T1", "T2", "T3"), "ID", joint = TRUE,
    bayes_nIter = 200L, bayes_burnIn = 100L, bayes_thin = 2L))
  expect_identical(unique(m$n_train), 90L)
  expect_identical(unique(m$min_pair_coobserved), 30L)
})

test_that("joint inputs adjust MET records and allow a trait confined to one environment", {
  d <- .joint_fixture()
  ph <- d$ph
  d$ph <- rbind(transform(ph, Env = "E1"), transform(ph, Env = "E2", T1 = T1 + 8, T2 = T2 - 6))
  d$ph$T1[d$ph$Env == "E2"] <- NA_real_
  d$ph$Env[1L] <- NA_character_
  z <- gp_feature_joint_inputs(d$X, d$ph, c("T1", "T2"), "ID", "Env")
  expect_equal(z$n_train, 179L)
  expect_equal(z$n_dropped, 1L)
  expect_identical(z$covariates_adjusted, "Env")
  expect_equal(z$trait_sd["T1"], c(T1 = sd(d$ph$T1[2:90])))
  m <- .joint_score(d, heter_groups = "Env")
  expect_true(all(is.finite(m$raw_score)))
  expect_identical(unique(m$n_records_dropped), 1L)
})

test_that("joint scores do not depend on trait units or ridge lambda", {
  d <- .joint_fixture()
  a <- .joint_score(d, ridge_lambda = 1)
  d$ph$T1 <- d$ph$T1 * 1e-9
  d$ph$T2 <- d$ph$T2 * 100
  b <- .joint_score(d, ridge_lambda = 1000)
  expect_equal(a$raw_score, b$raw_score, tolerance = 1e-10)
  expect_identical(a$predictor, b$predictor)
})

test_that("uninformative columns are integrated analytically and invalid panels stop", {
  d <- .joint_fixture()
  d$X[, "m12"] <- 1
  m <- .joint_score(d)
  j <- match("m12", m$predictor)
  expect_identical(m$raw_score[j], 0)
  expect_true(is.na(m$effect_sd[j]))
  expect_identical(m$scoring_status[j], "no_predictor_variation")
  d$X[,] <- 1
  expect_error(.joint_score(d), "informative")
  d <- .joint_fixture(); d$X[, "m12"] <- NA_real_
  expect_error(.joint_score(d), "m12")
  d <- .joint_fixture(); d$ph$T2 <- d$ph$T1
  expect_error(.joint_score(d), "singular")
})

test_that("indefinite pairwise covariance gets a declared positive-definite prior", {
  d <- .joint_fixture()
  gp_with_pinned_seed(19L, {
    ph <- data.frame(ID = rownames(d$X), T1 = NA_real_, T2 = NA_real_, T3 = NA_real_)
    v <- rep(seq(-2, 2, length.out = 30), 3)
    ph$T1[1:60] <- v[1:60]
    ph$T2[c(1:30, 61:90)] <- v[c(1:30, 61:90)]
    ph$T3[31:60] <- v[31:60]; ph$T3[61:90] <- -v[61:90]
    ph[, c("T1", "T2", "T3")] <- ph[, c("T1", "T2", "T3")] + matrix(rnorm(270, sd = .03), 90, 3)
    z <- gp_feature_joint_inputs(d$X, ph, c("T1", "T2", "T3"), "ID")
    expect_warning(pr <- gp_feature_joint_prior_setup(z$x, z$y, z$C), "prior covariance")
    expect_gt(pr$shrinkage, 0)
    expect_lt(pr$minimum_eigenvalue, 0)
    expect_gt(min(eigen(pr$covariance, symmetric = TRUE)$values), 0)
    expect_warning(m <- feature_score_predictors(
      d$X, ph, c("T1", "T2", "T3"), "ID", joint = TRUE,
      bayes_nIter = 200L, bayes_burnIn = 100L, bayes_thin = 2L), "prior covariance")
    expect_true(all(is.finite(m$raw_score)))
    expect_gt(unique(m$prior_covariance_shrinkage), 0)
  })
})

test_that("joint scorer validates Gaussian scope, trait visibility, chain and names", {
  d <- .joint_fixture()
  expect_error(feature_score_predictors(d$X, d$ph, "T1", "ID", joint = TRUE), "two")
  expect_error(.joint_score(d, response_family = "binary"), "Gaussian|gaussian")
  expect_error(.joint_score(d, fixed = ~ Env), "fixed")
  expect_error(.joint_score(d, heter_groups = "missing"), "environment")
  d$ph$T2 <- NA_real_; expect_error(.joint_score(d), "T2.*visible")
  d <- .joint_fixture(); d$ph$T2 <- 1; expect_error(.joint_score(d), "T2.*flat|T2.*constant")
  d <- .joint_fixture(); d$ph$T2[1] <- Inf; expect_error(.joint_score(d), "T2")
  d <- .joint_fixture(); d$X[1,1] <- Inf; expect_error(.joint_score(d), "Infinite predictor")
  d <- .joint_fixture()
  expect_error(gp_feature_score_joint_bglr(d$X, as.matrix(d$ph[,c("T1","T2")]),
                                          matrix(1,90,1), nIter=5, burnIn=4, thin=5), "retained")
  expect_error(gp_feature_score_joint_bglr(d$X, as.matrix(d$ph[,c("T1","T2")]),
                                          matrix(1,90,1), thin=0), "thin")
  colnames(d$X)[2] <- colnames(d$X)[1]; expect_error(.joint_score(d), "unique")
})

test_that("metadata extras remain aligned with predictors after ranking", {
  meta <- gp_feature_metadata_table(
    trait="T1+T2", predictors=c("a","b","c"), source_block=c(a="geno_data",b="geno_data",c="geno_data"),
    raw_score=c(1,3,2), scoring_model="Ridge_Regression", seed=1L, n_train=10L,
    scoring_scope="joint", extra=list(effect_sd=c(.1,.2,.3), chain_nIter=300L))
  expect_identical(meta$predictor, c("b","c","a"))
  expect_equal(meta$effect_sd, c(.2,.3,.1))
  expect_error(gp_feature_metadata_table(
    trait="T1+T2", predictors=c("a","b","c"), source_block=c(a="geno_data",b="geno_data",c="geno_data"),
    raw_score=c(1,3,2), scoring_model="Ridge_Regression", seed=1L, n_train=10L,
    extra=list(effect_sd=c(.1,.2))), "length")
})

test_that("joint scoring preserves RNG state and cleans task-owned sampler files", {
  d <- .joint_fixture()
  set.seed(901L); before <- .Random.seed; kind <- RNGkind()
  old_files <- Sys.glob(file.path(tempdir(), "feature_joint_*"))
  .joint_score(d)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_setequal(Sys.glob(file.path(tempdir(), "feature_joint_*")), old_files)
  rm(".Random.seed", envir=.GlobalEnv)
  .joint_score(d)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits=FALSE))
  assign(".Random.seed", before, .GlobalEnv)
  owned_connections <- list()
  local_mocked_bindings(Multitrait=function(...) {
    args <- list(...)
    con <- file(paste0(args$saveAt, "failure.dat"), "w")
    owned_connections[[length(owned_connections) + 1L]] <<- con
    writeLines("partial sampler output", con)
    runif(1)
    stop("sampler failure")
  }, .package="BGLR")
  expect_error(.joint_score(d), "sampler failure")
  expect_setequal(Sys.glob(file.path(tempdir(), "feature_joint_*")), old_files)
  expect_error(isOpen(owned_connections[[1L]]), "invalid connection")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  saved <- .Random.seed
  rm(".Random.seed", envir=.GlobalEnv)
  on.exit(assign(".Random.seed", saved, .GlobalEnv), add=TRUE)
  expect_error(.joint_score(d), "sampler failure")
  expect_false(exists(".Random.seed", .GlobalEnv, inherits=FALSE))
})

test_that("joint covariance calibration uses a proper residual shape for six traits", {
  d <- .joint_fixture()
  Y <- gp_with_pinned_seed(7L, matrix(rnorm(90 * 6), 90, 6))
  C <- matrix(1, 90, 1)
  pr <- gp_feature_joint_prior_setup(gp_feature_scale_matrix(d$X), Y, C)
  expect_identical(pr$residual$df0, 6L)
  expect_identical(pr$marker$df0, 7L)
  expect_gt(min(eigen(pr$residual$S0, symmetric=TRUE)$values), 0)
})
