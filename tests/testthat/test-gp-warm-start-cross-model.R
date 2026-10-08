# Pinning tests for Phase 2.2 cross-model warm-start sharing.
#
# The 0.20.19 cache key uses the canonical model name resolved via
# gp_backend_resolve_model(). LowRankGP (the canonical name behind the
# user-facing "Scalable-GBLUP" label) resolves to "KRR" -- so Scalable
# and Kernel-GBLUP fits on the same data SHARE the warm-start cache by
# construction. This test pins that behaviour so future refactors of
# the resolver don't silently break the sharing.

test_that("Scalable-GBLUP canonical name resolves to KRR (same as Kernel-GBLUP)", {
  expect_identical(PredictProR:::gp_backend_resolve_model("KRR"),        "KRR")
  expect_identical(PredictProR:::gp_backend_resolve_model("LowRankGP"),  "KRR")
  expect_identical(PredictProR:::gp_backend_resolve_model("GP"),         "GP")
  expect_identical(PredictProR:::gp_backend_resolve_model("GP_FA"),      "GP_FA")
})

test_that("Kernel-GBLUP and Scalable-GBLUP produce the SAME cache key", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  k_kernel <- PredictProR:::gp_warm_start_cache_key(
    model_name      = PredictProR:::gp_backend_resolve_model("KRR"),
    response        = "Yield", gen_name = "GID",
    heter_groups    = "Env",   gp_varcomp_mode = "reml", gmatrix = G
  )
  k_scalable <- PredictProR:::gp_warm_start_cache_key(
    model_name      = PredictProR:::gp_backend_resolve_model("LowRankGP"),
    response        = "Yield", gen_name = "GID",
    heter_groups    = "Env",   gp_varcomp_mode = "reml", gmatrix = G
  )
  expect_identical(k_kernel, k_scalable)
})

test_that("GP-GBLUP gets a SEPARATE cache key from KRR-equivalents", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  k_kernel <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("KRR"),
    "Yield", "GID", "Env", "reml", G
  )
  k_gp <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("GP"),
    "Yield", "GID", "Env", "reml", G
  )
  expect_false(identical(k_kernel, k_gp))
})

test_that("GP_FA gets a SEPARATE cache key (different param layout: psi+lambda+resid)", {
  G <- diag(10); rownames(G) <- colnames(G) <- paste0("g", 1:10)
  k_krr <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("KRR"),
    "Yield", "GID", "Env", "reml", G
  )
  k_fa <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("GP_FA"),
    "Yield", "GID", "Env", "reml", G
  )
  expect_false(identical(k_krr, k_fa))
})

test_that("save under Kernel-GBLUP, lookup as Scalable-GBLUP -> hit (same canonical name)", {
  G <- diag(8); rownames(G) <- colnames(G) <- paste0("g", 1:8)
  on.exit(PredictProR:::gp_warm_start_clear(), add = TRUE)
  PredictProR:::gp_warm_start_clear()

  k_kernel <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("KRR"),
    "Yield", "GID", "Env", "reml", G
  )
  theta_kernel <- c(0.63, 0.55, 0.42, 0.30)
  PredictProR:::gp_warm_start_save(k_kernel, theta_kernel, layout = "gp_exact")

  # Look up using "LowRankGP" (Scalable-GBLUP canonical) -> same key -> hit.
  k_scalable <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("LowRankGP"),
    "Yield", "GID", "Env", "reml", G
  )
  got <- PredictProR:::gp_warm_start_lookup(k_scalable)
  expect_equal(got, theta_kernel)
})

test_that("save under KRR, lookup as GP -> MISS (different canonical names)", {
  G <- diag(8); rownames(G) <- colnames(G) <- paste0("g", 1:8)
  on.exit(PredictProR:::gp_warm_start_clear(), add = TRUE)
  PredictProR:::gp_warm_start_clear()

  k_kernel <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("KRR"),
    "Yield", "GID", "Env", "reml", G
  )
  PredictProR:::gp_warm_start_save(k_kernel, c(0.6, 0.4), layout = "gp_exact")

  # GP has a different canonical name so the key differs and lookup misses.
  k_gp <- PredictProR:::gp_warm_start_cache_key(
    PredictProR:::gp_backend_resolve_model("GP"),
    "Yield", "GID", "Env", "reml", G
  )
  expect_null(PredictProR:::gp_warm_start_lookup(k_gp))
})
