# FA-GBLUP on single-env data is mathematically degenerate:
#   sigma2_g = sum_k lambda_{e,k}^2 + psi_e
# is not identifiable when there is only one env -- only the total scalar
# is identified, not the (lambda, psi) split. Empirically on real data
# the REML optimiser drives psi to absorb all genetic variance, producing
# sigma2_g 100x+ larger than Kernel-GBLUP/GP-GBLUP/Scalable-GBLUP on the
# same data (real benchmark: FA reported sigma2_g = 73.6 vs the other
# three models all returning sigma2_g = 0.598 on barley YIELD).
#
# The dispatcher in gp_backend_gaussian_model rejects single-env FA requests.
# It must not fit KRR (Kernel-GBLUP) under an FA-GBLUP label.
#
# This test verifies the fail-loud decision via a mocked backend, without
# running an actual GP fit (which would need Python + RAM). The mock
# captures the model_name that reaches the backend and the pheno_data
# n_env signal that triggered the decision.

test_that("GP_FA on heter_groups = NULL fails before fitting", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1.0, 2.0, NA, 3.0),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- pheno$GID
  seen <- list(called = FALSE)
  testthat::local_mocked_bindings(
    gp_single_trait_model = function(model_name, ...) {
      seen$called <<- TRUE
      seen$model_name <<- model_name
      list(
        predictions = data.frame(
          GID = pheno$GID,
          Predicted_value = c(1.0, 2.0, 2.0, 3.0),
          Standard_error = c(0.1, 0.1, 0.1, 0.1),
          PEV = c(0.01, 0.01, 0.01, 0.01),
          stringsAsFactors = FALSE
        ),
        result = list(),
        info = list(backend = "auto")
      )
    },
    .package = "PredictProR"
  )

  expect_error(
    PredictProR:::gp_backend_gaussian_model(
      model_name      = "GP_FA",
      pheno_data      = pheno,
      response        = "Yield",
      gen_name        = "GID",
      gmatrix         = K,
      heter_groups    = NULL,
      gp_backend      = "auto",
      gp_output_level = "predict_only"
    ),
    "FA-GBLUP requires multi-environment data"
  )
  expect_false(seen$called)
})

test_that("GP_FA with heter_groups carrying one unique env fails before fitting", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Env = c("E1", "E1", "E1", "E1"),       # single unique env
    Yield = c(1.0, 2.0, NA, 3.0),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- pheno$GID
  seen <- list(called = FALSE)
  testthat::local_mocked_bindings(
    gp_single_trait_model = function(model_name, ...) {
      seen$called <<- TRUE
      seen$model_name <<- model_name
      list(
        predictions = data.frame(
          GID = pheno$GID,
          Predicted_value = c(1.0, 2.0, 2.0, 3.0),
          Standard_error = c(0.1, 0.1, 0.1, 0.1),
          PEV = c(0.01, 0.01, 0.01, 0.01),
          stringsAsFactors = FALSE
        ),
        result = list(),
        info = list(backend = "auto")
      )
    },
    .package = "PredictProR"
  )

  expect_error(
    PredictProR:::gp_backend_gaussian_model(
      model_name      = "GP_FA",
      pheno_data      = pheno,
      response        = "Yield",
      gen_name        = "GID",
      gmatrix         = K,
      heter_groups    = "Env",
      gp_backend      = "auto",
      gp_output_level = "predict_only"
    ),
    "FA-GBLUP requires multi-environment data"
  )
  expect_false(seen$called)
})

test_that("GP_FA with >=2 distinct envs reaches the FA backend", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1.0, 1.5, 2.0, 2.5),
    stringsAsFactors = FALSE
  )
  K <- diag(2); rownames(K) <- colnames(K) <- c("g1", "g2")
  seen <- list()
  testthat::local_mocked_bindings(
    gp_single_trait_model = function(model_name, ...) {
      seen$model_name <<- model_name
      list(
        predictions = data.frame(
          GID = pheno$GID,
          Predicted_value = pheno$Yield,
          Standard_error = c(0.1, 0.1, 0.1, 0.1),
          PEV = c(0.01, 0.01, 0.01, 0.01),
          stringsAsFactors = FALSE
        ),
        result = list(),
        info = list(backend = "auto")
      )
    },
    .package = "PredictProR"
  )

  # No warning expected.
  out <- PredictProR:::gp_backend_gaussian_model(
    model_name      = "GP_FA",
    pheno_data      = pheno,
    response        = "Yield",
    gen_name        = "GID",
    gmatrix         = K,
    heter_groups    = "Env",
    gp_backend      = "auto",
    gp_output_level = "predict_only"
  )
  expect_identical(seen$model_name, "GP_FA")
  expect_false("fa_auto_routed_from" %in% out$model_parameters$stat)
})

test_that("the public general-data scope distinguishes FA structured fits", {
  scope <- general_prediction_data_standard()$model_families

  expect_setequal(
    scope$single_response_single_environment_kernel_or_relationship_models,
    c(
      "GBLUP_BRR", "RKHS", "GBLUP", "Kernel-GBLUP",
      "Gaussian-Process-GBLUP"
    )
  )
  expect_identical(scope$structured_multi_response_gp_models, "FA-GBLUP")
})

test_that("non-FA models are never auto-routed regardless of single-env", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1.0, 2.0, NA, 3.0),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- pheno$GID
  seen <- list()
  testthat::local_mocked_bindings(
    gp_single_trait_model = function(model_name, ...) {
      seen$model_name <<- model_name
      list(
        predictions = data.frame(
          GID = pheno$GID,
          Predicted_value = c(1.0, 2.0, 2.0, 3.0),
          Standard_error = c(0.1, 0.1, 0.1, 0.1),
          PEV = c(0.01, 0.01, 0.01, 0.01),
          stringsAsFactors = FALSE
        ),
        result = list(),
        info = list(backend = "auto")
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_backend_gaussian_model(
    model_name      = "KRR",
    pheno_data      = pheno,
    response        = "Yield",
    gen_name        = "GID",
    gmatrix         = K,
    heter_groups    = NULL,
    gp_backend      = "auto",
    gp_output_level = "predict_only"
  )
  expect_identical(seen$model_name, "KRR")
  expect_false("fa_auto_routed_from" %in% out$model_parameters$stat)
})
