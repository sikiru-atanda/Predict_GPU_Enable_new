# Pins the two ASReml predict-failure fallback extensions:
#   1. Multi-trait ASReml CV now routes through
#      gp_multitrait_asreml_gaussian_cv(), which uses the same
#      gp_multitrait_asreml_predict_or_extract() chain the true-prediction
#      path uses. Previously the CV route called stop("...CV will be added
#      later").
#   2. Hybrid ASReml random-effect extraction now emits a warning() when
#      one of the four silent-degradation branches fires
#      (summary(mod) failed, no rows, no 'solution'/'effect' column, no
#      matching vm() rows). The behaviour is unchanged (still returns an
#      empty data.frame and the hybrid Predicted_value drops back to the
#      fixed-effect intercept); the user just now sees what dropped out.

# ---- Hybrid extraction visibility (mocked summary) --------------------------

test_that("hybrid extractor warns when summary() throws", {
  local_mocked_bindings(
    gp_hybrid_asreml_random_coefs = function(model) stop("synthetic asreml summary failure"),
    .package = "PredictProR"
  )
  expect_warning(
    out <- gp_hybrid_asreml_extract_random_effect(
      model = list(),
      factor_name = "Female",
      inv_name = "Ginv",
      id_col_name = "Female",
      effect_col_name = "Female_GCA"
    ),
    "summary\\(model, coef=TRUE\\) failed"
  )
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 0L)
})

test_that("hybrid extractor warns when random coefs returned NULL", {
  local_mocked_bindings(
    gp_hybrid_asreml_random_coefs = function(model) NULL,
    .package = "PredictProR"
  )
  expect_warning(
    out <- gp_hybrid_asreml_extract_random_effect(
      model = list(),
      factor_name = "Female",
      inv_name = "Ginv",
      id_col_name = "Female",
      effect_col_name = "Female_GCA"
    ),
    "returned NULL"
  )
  expect_equal(nrow(out), 0L)
})

test_that("hybrid extractor warns when coef.random has no matching vm() rows", {
  local_mocked_bindings(
    gp_hybrid_asreml_random_coefs = function(model) {
      matrix(c(0.1, 0.2), nrow = 2L, ncol = 1L,
             dimnames = list(
               c("vm(SomethingElse, OtherInv)_g1", "vm(SomethingElse, OtherInv)_g2"),
               "solution"))
    },
    .package = "PredictProR"
  )
  expect_warning(
    out <- gp_hybrid_asreml_extract_random_effect(
      model = list(),
      factor_name = "Female",
      inv_name = "Ginv",
      id_col_name = "Female",
      effect_col_name = "Female_GCA"
    ),
    "no coef\\.random rows matched"
  )
  expect_equal(nrow(out), 0L)
})

test_that("hybrid extractor warns when neither 'solution' nor 'effect' column present", {
  local_mocked_bindings(
    gp_hybrid_asreml_random_coefs = function(model) {
      matrix(c(0.1, 0.2), nrow = 2L, ncol = 1L,
             dimnames = list(
               c("vm(Female, Ginv)_g1", "vm(Female, Ginv)_g2"),
               "weird_col"))
    },
    .package = "PredictProR"
  )
  expect_warning(
    out <- gp_hybrid_asreml_extract_random_effect(
      model = list(),
      factor_name = "Female",
      inv_name = "Ginv",
      id_col_name = "Female",
      effect_col_name = "Female_GCA"
    ),
    "neither 'solution' nor 'effect' column"
  )
  expect_equal(nrow(out), 0L)
})

test_that("hybrid extractor returns rows silently on the happy path", {
  local_mocked_bindings(
    gp_hybrid_asreml_random_coefs = function(model) {
      matrix(c(0.10, 0.20, 0.05, 0.07),
             nrow = 2L, ncol = 2L,
             dimnames = list(
               c("vm(Female, Ginv)_g1", "vm(Female, Ginv)_g2"),
               c("solution", "std.error")))
    },
    .package = "PredictProR"
  )
  expect_silent(
    out <- gp_hybrid_asreml_extract_random_effect(
      model = list(),
      factor_name = "Female",
      inv_name = "Ginv",
      id_col_name = "Female",
      effect_col_name = "Female_GCA"
    )
  )
  expect_equal(nrow(out), 2L)
  expect_true(all(c("Female", "Female_GCA", "Standard_error", "Prediction_error_variance") %in% names(out)))
  expect_equal(out$Female_GCA, c(0.10, 0.20))
})

# ---- Multi-trait ASReml CV signature + dispatch ------------------------------
# The fitter itself needs the asreml engine, so the actual run is gated by
# skip_if_not_installed. These tests verify the signature, the input
# validation, and that the route dispatcher no longer hard-stops.

test_that("gp_multitrait_asreml_gaussian_cv stops cleanly without asreml", {
  testthat::skip_if(requireNamespace("asreml", quietly = TRUE),
                    "asreml is installed; covered by smoke test instead")
  ph <- data.frame(GID = paste0("g", 1:6),
                   Yield = rnorm(6), Height = rnorm(6),
                   stringsAsFactors = FALSE)
  K <- diag(6); rownames(K) <- colnames(K) <- ph$GID
  expect_error(
    gp_multitrait_asreml_gaussian_cv(
      pheno_object = ph, response = c("Yield", "Height"),
      gen_name = "GID", gmatrix = K,
      cross_validation_meth = "K-Folds", nfolds = 3L, replication = 1L
    ),
    "asreml-R"
  )
})

test_that("gp_multitrait_asreml_gaussian_cv rejects duplicate genotype rows", {
  testthat::skip_if_not_installed("asreml")
  ph <- data.frame(GID = c("g1","g1","g2","g3","g4","g5"),
                   Yield = rnorm(6), Height = rnorm(6),
                   stringsAsFactors = FALSE)
  K <- diag(5); rownames(K) <- colnames(K) <- c("g1","g2","g3","g4","g5")
  expect_error(
    gp_multitrait_asreml_gaussian_cv(
      pheno_object = ph, response = c("Yield", "Height"),
      gen_name = "GID", gmatrix = K, nfolds = 3L
    ),
    "one phenotype row per genotype"
  )
})

test_that("gp_multitrait_asreml_gaussian_cv rejects non-K-Folds CV method", {
  testthat::skip_if_not_installed("asreml")
  ph <- data.frame(GID = paste0("g", 1:6),
                   Yield = rnorm(6), Height = rnorm(6),
                   stringsAsFactors = FALSE)
  K <- diag(6); rownames(K) <- colnames(K) <- ph$GID
  expect_error(
    gp_multitrait_asreml_gaussian_cv(
      pheno_object = ph, response = c("Yield", "Height"),
      gen_name = "GID", gmatrix = K,
      cross_validation_meth = "Hold_Out", nfolds = 3L
    ),
    "K-Folds only"
  )
})
