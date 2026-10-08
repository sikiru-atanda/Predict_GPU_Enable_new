find_hybrid_gp_consistency_root <- function() {
  candidates <- unique(c(
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = FALSE)
  ))
  hits <- candidates[file.exists(file.path(candidates, "tools", "validate_hybrid_gp_performance_consistency.R"))]
  testthat::skip_if_not(length(hits) > 0L, "Hybrid GP performance consistency script was not found.")
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("hybrid GP performance consistency gate passes", {
  skip_if_not(
    identical(Sys.getenv("PREDICTPROR_RUN_GP_PERFORMANCE_TESTS"), "1"),
    "Set PREDICTPROR_RUN_GP_PERFORMANCE_TESTS=1 to run hybrid GP performance consistency tests."
  )

  root <- find_hybrid_gp_consistency_root()
  script <- file.path(root, "tools", "validate_hybrid_gp_performance_consistency.R")
  out_dir <- file.path(root, "tools", "tmp_hybrid_gp_performance_consistency")
  checks_file <- file.path(out_dir, "hybrid_gp_consistency_checks.csv")
  env_keys <- c(
    "PREDICTPROR_HYBRID_GP_CONSISTENCY_SEEDS",
    "PREDICTPROR_HYBRID_GP_CONSISTENCY_BACKEND",
    "PREDICTPROR_HYBRID_GP_CONSISTENCY_STRICT"
  )
  old_env <- Sys.getenv(env_keys, unset = NA_character_)
  on.exit({
    for (nm in env_keys) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        Sys.setenv(structure(old_env[[nm]], names = nm))
      }
    }
  }, add = TRUE)
  Sys.setenv(
    PREDICTPROR_HYBRID_GP_CONSISTENCY_SEEDS = "7101",
    PREDICTPROR_HYBRID_GP_CONSISTENCY_BACKEND = "r",
    PREDICTPROR_HYBRID_GP_CONSISTENCY_STRICT = "1"
  )

  child_expr <- sprintf(
    "setwd(%s); source(%s)",
    deparse(root),
    deparse(script)
  )
  out <- system2(
    command = "Rscript",
    args = c("--vanilla", "-e", shQuote(child_expr)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(as.integer(status), 0L)) {
    testthat::fail(paste(c("Hybrid GP performance consistency script failed:", out), collapse = "\n"))
  }

  expect_true(file.exists(checks_file))
  checks <- utils::read.csv(checks_file, stringsAsFactors = FALSE)
  expect_true(all(c("check", "status", "detail") %in% names(checks)))
  expect_false(any(checks$status == "fail"))
  expect_true(all(c("finite_uncertainty", "accuracy_floor", "deterministic_r_backend") %in% checks$check))
})
