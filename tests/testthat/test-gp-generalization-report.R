test_that("GP generalization report and regression gate accept manifest inputs", {
  testthat::skip_if_not_installed("data.table")
  testthat::skip_if_not_installed("withr")

  find_repo_root <- function() {
    candidates <- unique(c(
      getwd(),
      normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
      normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE),
      normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = FALSE)
    ))
    hits <- candidates[file.exists(file.path(candidates, "tools", "gp_generalization_validation_report.R"))]
    testthat::skip_if_not(length(hits) > 0L, "GP generalization report script was not found.")
    normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
  }

  root <- find_repo_root()
  tmp <- tempfile("gp-generalization-test-")
  dir.create(tmp, recursive = TRUE)
  out_dir <- file.path(tmp, "report")
  dir.create(out_dir, recursive = TRUE)

  g2f_summary <- data.table::data.table(
    dataset = "G2F_2025",
    scenario = "historical_multi_environment_cv",
    model = c("EnvMean", "GBLUP_env_kernel_tuned", "GP_env_covariates_tuned"),
    runs = 3L,
    total_test = 40L,
    mean_rmse = c(2.0, 1.8, 1.7),
    mean_pearson_cor = c(0.50, 0.60, 0.65),
    mean_se = c(NA_real_, NA_real_, 1.10),
    mean_pev = c(NA_real_, NA_real_, 1.21)
  )
  zenodo_summary <- data.table::data.table(
    dataset = "Zenodo_3GS_wheat",
    scenario = "genotype_holdout_single_environment_cv",
    model = c("EnvMean", "GBLUP_same_lambda", "GP_auto"),
    folds = 3L,
    total_scored = 30L,
    mean_rmse = c(0.50, 0.60, 0.62),
    mean_pearson_cor = c(NA_real_, 0.18, 0.17),
    mean_se = c(NA_real_, NA_real_, 0.20),
    mean_pev = c(NA_real_, NA_real_, 0.04)
  )

  g2f_file <- file.path(tmp, "g2f_summary.csv")
  zenodo_file <- file.path(tmp, "zenodo_summary.csv")
  manifest_file <- file.path(tmp, "manifest.csv")
  data.table::fwrite(g2f_summary, g2f_file)
  data.table::fwrite(zenodo_summary, zenodo_file)
  data.table::fwrite(
    data.table::data.table(
      run_group = c("g2f_synthetic", "zenodo_synthetic"),
      path = c(normalizePath(g2f_file, winslash = "/", mustWork = TRUE),
               normalizePath(zenodo_file, winslash = "/", mustWork = TRUE))
    ),
    manifest_file
  )

  withr::local_envvar(c(
    PREDICTPROR_GP_GENERALIZATION_OUT_DIR = normalizePath(out_dir, winslash = "/", mustWork = TRUE),
    PREDICTPROR_GP_GENERALIZATION_SUMMARY_MANIFEST = normalizePath(manifest_file, winslash = "/", mustWork = TRUE),
    PREDICTPROR_GP_REQUIRE_REFERENCE_ROWS = "0"
  ))

  source(file.path(root, "tools", "gp_generalization_validation_report.R"), local = new.env(parent = globalenv()))
  comparison <- data.table::fread(file.path(out_dir, "gp_generalization_validation_comparison.csv"))
  expect_equal(nrow(comparison), 2L)
  expect_true(all(c(
    "rmse_improvement_pct_vs_baseline",
    "beats_baseline_rmse",
    "beats_gblup_rmse"
  ) %in% names(comparison)))
  g2f_row <- comparison[comparison[["dataset"]] == "G2F_2025", ]
  zenodo_row <- comparison[comparison[["dataset"]] == "Zenodo_3GS_wheat", ]
  expect_equal(nrow(g2f_row), 1L)
  expect_equal(nrow(zenodo_row), 1L)
  expect_true(isTRUE(g2f_row[["beats_baseline_rmse"]]))
  expect_false(isTRUE(zenodo_row[["beats_baseline_rmse"]]))

  source(file.path(root, "tools", "check_gp_generalization_regression.R"), local = new.env(parent = globalenv()))
  check_summary <- data.table::fread(file.path(out_dir, "gp_generalization_regression_check_summary.csv"))
  expect_equal(check_summary$failures, 0L)
  expect_equal(check_summary$gp_rows, 2L)
})
