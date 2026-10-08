# Rigor patch (ASReml finding #2): the REML retry loops in asreml_utilis_new.R,
# asreml_cv_model.R and multi_trait_asreml_helpers.R retry update.asreml up to 3x
# but historically returned the model with NO signal when it still failed to
# converge. A non-converged REML fit yields unreliable variance components, EBVs
# and PEVs, so the user must be warned. This tests the shared warning helper
# (pure R; needs no asreml license to verify).

resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  src <- file.path("..", "..", "R", "asreml_utilis_new.R")
  if (file.exists(src)) {
    env <- new.env(parent = globalenv())
    suppressWarnings(try(sys.source(src, envir = env), silent = TRUE))
    if (exists(name, envir = env, mode = "function")) return(get(name, envir = env, mode = "function"))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    f <- tryCatch(getFromNamespace(name, "PredictProR"), error = function(e) NULL)
    if (!is.null(f)) return(f)
  }
  NULL
}

test_that("gp_asreml_warn_if_not_converged warns only when converge is explicitly FALSE", {
  fn <- resolve_fn("gp_asreml_warn_if_not_converged")
  if (is.null(fn)) skip("gp_asreml_warn_if_not_converged not available")

  # Not converged after retries -> must warn, with an informative message.
  expect_warning(fn(list(converge = FALSE), context = "single-location GBLUP"),
                 "did not converge")

  # Converged -> silent (no false alarm).
  expect_silent(fn(list(converge = TRUE)))

  # Unknown convergence state (NULL / missing field, e.g. a non-asreml object or
  # a failed fit returning NULL) -> silent, so we don't raise spurious warnings.
  expect_silent(fn(list(converge = NULL)))
  expect_silent(fn(list()))
  expect_silent(fn(NULL))

  # The model object is returned invisibly so the helper can be used inline.
  mod <- list(converge = TRUE, x = 1)
  expect_identical(suppressWarnings(fn(mod)), mod)
})

test_that("gp_asreml_warn_if_not_converged includes the context in the warning", {
  fn <- resolve_fn("gp_asreml_warn_if_not_converged")
  if (is.null(fn)) skip("gp_asreml_warn_if_not_converged not available")
  expect_warning(fn(list(converge = FALSE), context = "multi-trait GBLUP"),
                 "multi-trait GBLUP")
})
