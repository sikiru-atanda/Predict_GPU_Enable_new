test_that("platform helpers normalize supported OS families", {
  expect_identical(PredictProR:::gp_platform_family("Windows", "windows"), "windows")
  expect_identical(PredictProR:::gp_platform_family("Linux", "unix"), "linux")
  expect_identical(PredictProR:::gp_platform_family("Darwin", "unix"), "macos")
  expect_identical(PredictProR:::gp_platform_family("FreeBSD", "unix"), "linux")
})

test_that("Python executable layouts cover Windows and Unix virtualenvs", {
  win <- PredictProR:::gp_virtualenv_python_path("root", "env", os_type = "windows")
  unix <- PredictProR:::gp_virtualenv_python_path("root", "env", os_type = "unix")

  expect_match(gsub("\\\\", "/", win), "root/env/Scripts/python[.]exe$")
  expect_match(gsub("\\\\", "/", unix), "root/env/bin/python$")
  expect_identical(
    PredictProR:::gp_system_python_commands("unix"),
    c("python3", "python")
  )
})

test_that("external command arguments survive paths and values with spaces", {
  spaced_dir <- file.path(tempdir(), "PredictProR portability space")
  dir.create(spaced_dir, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(spaced_dir, "echo args.R")
  writeLines("cat(commandArgs(trailingOnly = TRUE), sep = '\\n')", script)

  rscript <- file.path(
    R.home("bin"),
    if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript"
  )
  out <- system2(
    rscript,
    args = PredictProR:::gp_quote_system_args(c("--vanilla", script, "value with spaces")),
    stdout = TRUE,
    stderr = TRUE
  )

  # Some configured Windows R installations emit locale-startup warnings even
  # under --vanilla. The portability assertion concerns the trailing argument,
  # not unrelated interpreter startup diagnostics.
  expect_identical(tail(out[nzchar(out)], 1L), "value with spaces")
})

test_that("core detection fails closed on unavailable host information", {
  expect_identical(PredictProR:::gp_detect_cores(detected = NA_integer_), 1L)
  expect_identical(PredictProR:::gp_detect_cores(detected = 8L, reserve = 1L), 7L)
  expect_identical(PredictProR:::gp_detect_cores(detected = 1L, reserve = 1L), 1L)
})

test_that("forking is disabled when workers require process initialization", {
  expect_true(PredictProR:::gp_parallel_safe_fork_preference(TRUE, FALSE))
  expect_false(PredictProR:::gp_parallel_safe_fork_preference(TRUE, TRUE))
  expect_false(PredictProR:::gp_parallel_safe_fork_preference(FALSE, FALSE))
})

test_that("shipped R runtime has no developer-specific absolute path", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
  source <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)

  # any Windows user-profile path (C:/Users/<name>/ or C:\Users\<name>\)
  expect_false(any(grepl("[A-Za-z]:[/\\\\]Users[/\\\\][^/\\\\\"']+[/\\\\]", source)))
})

# Shipped docs and Python modules too (benchmark/test/scratch scripts are
# .Rbuildignored and excluded here the same way).
test_that("shipped docs and Python modules have no developer-specific absolute path", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  inst <- file.path(root, "inst")
  skip_if_not(dir.exists(file.path(inst, "python")), "source tree not available")
  files <- c(list.files(inst, pattern = "[.]md$", full.names = TRUE),
             list.files(file.path(inst, "python"), pattern = "[.]py$", full.names = TRUE, recursive = TRUE))
  ignored <- "(^|/)(bench_|test_|smoke_|probe_|export_|diag_)|old|OLD|legacy|experimental"
  files <- files[!grepl(ignored, sub(paste0("^", inst, "/"), "", files))]
  hits <- Filter(function(f) any(grepl("[A-Za-z]:[/\\\\]+Users[/\\\\]+[^/\\\\\"'<> ]+[/\\\\]|D:/PredictProR",
                                       readLines(f, warn = FALSE))), files)
  expect_identical(basename(hits), character())
})

# Exported runs used to leave the session in the output folder when
# PredictProR.output_dir / PREDICTPRO_OUTPUT_DIR pointed elsewhere.
test_that("an exported run writes to the output folder and keeps the working directory", {
  skip_on_cran()
  set.seed(3)
  ids <- sprintf("L%02d", 1:30)
  G <- matrix(2L * stats::rbinom(30 * 60, 1, 0.5), 30, dimnames = list(ids, sprintf("m%02d", 1:60)))
  ph <- data.frame(GID = ids, Yield = as.numeric(G[, 1:5] %*% stats::rnorm(5)) + stats::rnorm(30))
  ph$Yield[1:5] <- NA
  start <- withr::local_tempdir()
  withr::local_dir(start)
  out_dir <- file.path(withr::local_tempdir(), "exports")
  withr::local_options(PredictProR.output_dir = out_dir)
  res <- suppressWarnings(model_execute(
    pheno_data = ph, geno_data = G, response = "Yield", gen_name = "GID", random = ~GID, GS_model = "BRR",
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, nIter = 120L, burnIn = 20L, thin = 2L,
    parallel_mode = "sequential", system_database = FALSE, message = FALSE, verbose = FALSE
  ))
  expect_identical(normalizePath(getwd(), winslash = "/"), normalizePath(start, winslash = "/"))
  expect_true(length(list.files(out_dir, pattern = "Predicted_Value\\.csv$", recursive = TRUE)) >= 1L)
})

test_that("results export folder follows the option, then the env var, then the working directory", {
  withr::local_options(PredictProR.output_dir = NULL)
  withr::local_envvar(PREDICTPRO_OUTPUT_DIR = NA)
  expect_identical(PredictProR:::gp_output_root(), getwd())

  env_dir <- file.path(withr::local_tempdir(), "from-env")
  withr::local_envvar(PREDICTPRO_OUTPUT_DIR = env_dir)
  expect_identical(PredictProR:::gp_output_root(create = TRUE), normalizePath(env_dir, winslash = "/"))
  expect_true(dir.exists(env_dir))

  opt_dir <- withr::local_tempdir()
  withr::local_options(PredictProR.output_dir = opt_dir)
  expect_identical(PredictProR:::gp_output_root(), normalizePath(opt_dir, winslash = "/"))
})

test_that("base_parallel honors PSOCK opt-out from Unix forking", {
  skip_if_not(identical(.Platform$OS.type, "unix"), "Unix-specific PSOCK behavior")

  parent_pid <- Sys.getpid()
  out <- PredictProR:::sp_apply(
    X = 1:2,
    FUN = function(i) Sys.getpid(),
    decision = list(
      backend = "base_parallel",
      workers = 1L,
      prefer_fork = FALSE
    ),
    packages = character(),
    seed = FALSE
  )

  expect_true(all(as.integer(unlist(out)) != parent_pid))
})
