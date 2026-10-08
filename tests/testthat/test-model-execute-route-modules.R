test_that("model_execute route modules expose all major model families", {
  ns <- asNamespace("PredictProR")
  route_functions <- c(
    "gp_model_execute_route_families",
    "gp_route_pre_input_specialized_models",
    "gp_route_cross_validation_specialized_models",
    "gp_route_true_prediction_specialized_models"
  )
  expect_true(all(vapply(route_functions, exists, logical(1), envir = ns, inherits = FALSE)))

  routes <- get("gp_model_execute_route_families", envir = ns)()
  expect_true(all(c("route", "family", "stage") %in% names(routes)))
  expect_true(all(c(
    "hybrid_asreml",
    "hybrid_bayes",
    "hybrid_gp",
    "hybrid_ml",
    "hybrid_dl",
    "multi_trait_asreml",
    "multi_trait_ml",
    "multi_trait_dl",
    "multi_trait_gp",
    "met_ml_dl",
    "cross_validation",
    "general_true_prediction"
  ) %in% routes$route))
  expect_true(all(c(
    "asreml",
    "bayesian",
    "gp",
    "classical_ml",
    "deep_learning",
    "met_ml_dl",
    "general"
  ) %in% routes$family))
})

test_that("model_execute delegates specialized route blocks to route modules", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  source_file <- file.path(root, "R", "model_execute_para.R")
  if (!file.exists(source_file)) {
    testthat::skip("raw R sources are not available in installed-package test context")
  }
  src <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  expect_true(grepl("gp_route_pre_input_specialized_models\\(as.list\\(environment\\(\\)\\)\\)", src))
  expect_true(grepl("gp_route_cross_validation_specialized_models\\(as.list\\(environment\\(\\)\\)\\)", src))
  expect_true(grepl("gp_route_true_prediction_specialized_models\\(as.list\\(environment\\(\\)\\)\\)", src))
})

test_that("model_execute route context import gives every declared context name a local binding", {
  ns <- asNamespace("PredictProR")
  import_context <- get("gp_model_execute_import_context", envir = ns)
  context_names <- get("gp_model_execute_route_context_variables", envir = ns)()

  exists_after_import <- local({
    import_context(list())
    import_env <- environment()
    vapply(
      context_names,
      exists,
      logical(1),
      envir = import_env,
      inherits = FALSE
    )
  })

  if (!all(exists_after_import)) {
    testthat::fail(
      paste(
        "Missing imported context bindings:",
        paste(names(exists_after_import)[!exists_after_import], collapse = ", ")
      )
    )
  }
  testthat::succeed()
})

test_that("model_execute route context import preserves explicit NULL caller values", {
  ns <- asNamespace("PredictProR")
  import_context <- get("gp_model_execute_import_context", envir = ns)

  null_context_values <- local({
    import_context(list(fixed = NULL, random = NULL, cova = NULL))
    import_env <- environment()
    vapply(
      c("fixed", "random", "cova"),
      exists,
      logical(1),
      envir = import_env,
      inherits = FALSE
    )
  })

  if (!all(null_context_values)) {
    testthat::fail(
      paste(
        "Explicit NULL context values were dropped:",
        paste(names(null_context_values)[!null_context_values], collapse = ", ")
      )
    )
  }
  testthat::succeed()
})

test_that("model_execute route modules tolerate inactive partial contexts", {
  ns <- asNamespace("PredictProR")
  expect_null(get("gp_route_pre_input_specialized_models", envir = ns)(list()))
  expect_null(get("gp_route_cross_validation_specialized_models", envir = ns)(list()))
  expect_null(get("gp_route_true_prediction_specialized_models", envir = ns)(list()))
})

test_that("model_execute route modules use only declared context variables", {
  ns <- asNamespace("PredictProR")
  context_names <- get("gp_model_execute_route_context_variables", envir = ns)()
  allowed <- c(context_names, ls(ns, all.names = TRUE), "is.null")
  route_functions <- c(
    "gp_route_pre_input_specialized_models",
    "gp_route_cross_validation_specialized_models",
    "gp_route_true_prediction_specialized_models"
  )

  for (route_function in route_functions) {
    globals <- codetools::findGlobals(
      get(route_function, envir = ns),
      merge = FALSE
    )$variables
    expect_equal(
      setdiff(globals, allowed),
      character(),
      info = route_function
    )
  }
})
