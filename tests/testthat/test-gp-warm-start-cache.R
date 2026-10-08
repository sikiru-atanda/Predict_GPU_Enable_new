# Pinning tests for the Phase 2.1 Stage C session theta cache.
#
# The cache stores the converged variance-component theta returned by the
# Python AI-REML backend, keyed by
# (canonical_model_name, response, gen_name, heter_groups, gp_varcomp_mode,
#  gmatrix-digest). On a subsequent fit with the same key the cached theta
# is supplied as theta0_warm so AI-Newton converges in 1-2 iterations
# instead of 8-12.

test_that("gp_warm_start_cache_key is deterministic for the same inputs", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  k1 <- PredictProR:::gp_warm_start_cache_key(
    model_name      = "KRR", response = "Yield", gen_name = "GID",
    heter_groups    = "Env", gp_varcomp_mode = "reml", gmatrix = G
  )
  k2 <- PredictProR:::gp_warm_start_cache_key(
    model_name      = "KRR", response = "Yield", gen_name = "GID",
    heter_groups    = "Env", gp_varcomp_mode = "reml", gmatrix = G
  )
  expect_identical(k1, k2)
})

test_that("gp_warm_start_cache_key differs when any component changes", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  base <- PredictProR:::gp_warm_start_cache_key(
    model_name = "KRR", response = "Yield", gen_name = "GID",
    heter_groups = "Env", gp_varcomp_mode = "reml", gmatrix = G
  )
  # Model differs
  expect_false(identical(base, PredictProR:::gp_warm_start_cache_key(
    model_name = "GP_FA", response = "Yield", gen_name = "GID",
    heter_groups = "Env", gp_varcomp_mode = "reml", gmatrix = G
  )))
  # Response differs
  expect_false(identical(base, PredictProR:::gp_warm_start_cache_key(
    model_name = "KRR", response = "Height", gen_name = "GID",
    heter_groups = "Env", gp_varcomp_mode = "reml", gmatrix = G
  )))
  # heter_groups differs (single-env vs MET)
  expect_false(identical(base, PredictProR:::gp_warm_start_cache_key(
    model_name = "KRR", response = "Yield", gen_name = "GID",
    heter_groups = NULL, gp_varcomp_mode = "reml", gmatrix = G
  )))
  # gmatrix differs
  G2 <- G; G2[1, 1] <- 2.0
  expect_false(identical(base, PredictProR:::gp_warm_start_cache_key(
    model_name = "KRR", response = "Yield", gen_name = "GID",
    heter_groups = "Env", gp_varcomp_mode = "reml", gmatrix = G2
  )))
})

test_that("save/lookup cycle returns the same theta vector and tolerates NULL", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  key <- PredictProR:::gp_warm_start_cache_key(
    "KRR", "Yield", "GID", "Env", "reml", G
  )
  theta <- c(0.6, 0.5, 0.4, 0.3)
  on.exit(PredictProR:::gp_warm_start_clear(), add = TRUE)

  PredictProR:::gp_warm_start_clear()
  expect_null(PredictProR:::gp_warm_start_lookup(key))

  PredictProR:::gp_warm_start_save(key, theta, layout = "gp_exact")
  got <- PredictProR:::gp_warm_start_lookup(key)
  expect_equal(got, theta)

  # NULL / non-finite refused
  PredictProR:::gp_warm_start_save(key = "k2", theta = NULL)
  expect_null(PredictProR:::gp_warm_start_lookup("k2"))
  PredictProR:::gp_warm_start_save(key = "k3", theta = c(0.5, NA))
  expect_null(PredictProR:::gp_warm_start_lookup("k3"))
})

test_that("gp_warm_start_clear empties the session cache", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  k <- PredictProR:::gp_warm_start_cache_key(
    "KRR", "Yield", "GID", "Env", "reml", G
  )
  PredictProR:::gp_warm_start_save(k, c(0.5, 0.3))
  expect_false(is.null(PredictProR:::gp_warm_start_lookup(k)))
  PredictProR:::gp_warm_start_clear()
  expect_null(PredictProR:::gp_warm_start_lookup(k))
})

test_that("lookup is keyed by gmatrix content (different content -> separate slot)", {
  G1 <- diag(10); rownames(G1) <- colnames(G1) <- paste0("g", 1:10)
  G2 <- G1; G2[2, 2] <- 1.5
  k1 <- PredictProR:::gp_warm_start_cache_key("KRR","Y","GID","Env","reml", G1)
  k2 <- PredictProR:::gp_warm_start_cache_key("KRR","Y","GID","Env","reml", G2)
  on.exit(PredictProR:::gp_warm_start_clear(), add = TRUE)
  PredictProR:::gp_warm_start_clear()
  PredictProR:::gp_warm_start_save(k1, c(0.7, 0.4))
  PredictProR:::gp_warm_start_save(k2, c(0.65, 0.45))
  expect_equal(PredictProR:::gp_warm_start_lookup(k1), c(0.7, 0.4))
  expect_equal(PredictProR:::gp_warm_start_lookup(k2), c(0.65, 0.45))
})
