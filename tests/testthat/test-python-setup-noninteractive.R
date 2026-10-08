test_that("setup_predictgp_env provisions through subprocess without switching reticulate", {
  skip_if_not_installed("reticulate")

  conda_envs <- data.frame(
    name = "predictgp",
    python = "C:/conda/envs/predictgp/python.exe",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    conda_binary = function(...) "conda",
    conda_list = function(...) conda_envs,
    use_condaenv = function(...) stop("setup_predictgp_env should not switch reticulate Python"),
    .package = "reticulate"
  )
  local_mocked_bindings(
    gp_bridge_run_cli = function(command, args, python_bin = NULL) {
      expect_identical(command, "setup-deps")
      expect_identical(python_bin, "C:/conda/envs/predictgp/python.exe")
      c(
        "python=C:/conda/envs/predictgp/python.exe",
        "torch_version=2.4.0",
        "gpytorch_version=1.12",
        "cuda_available=false",
        "cuda_device_count=0",
        "torch_install_flavor=cpu"
      )
    },
    .package = "PredictProR"
  )

  out <- setup_predictgp_env(quiet = TRUE, prefer_gpu = FALSE, cuda = "cpu")
  expect_identical(out$python, "C:/conda/envs/predictgp/python.exe")
  expect_false(out$cuda_available)
  expect_identical(Sys.getenv("PREDICTPRO_GP_PYTHON"), "C:/conda/envs/predictgp/python.exe")
})

test_that("setup_predictdl_env provisions through subprocess without switching reticulate", {
  skip_if_not_installed("reticulate")

  conda_envs <- data.frame(
    name = "predictdl",
    python = "C:/conda/envs/predictdl/python.exe",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    conda_binary = function(...) "conda",
    conda_list = function(...) conda_envs,
    use_condaenv = function(...) stop("setup_predictdl_env should not switch reticulate Python"),
    .package = "reticulate"
  )
  local_mocked_bindings(
    gp_dl_run_cli = function(command, args, python_bin = NULL) {
      expect_identical(command, "setup-deps")
      expect_identical(python_bin, "C:/conda/envs/predictdl/python.exe")
      c(
        "python=C:/conda/envs/predictdl/python.exe",
        "pytorch_version=2.4.0",
        "cuda_compiled=12.1",
        "cuda_available=true",
        "cuda_device_count=1",
        "gpu_name=Mock GPU",
        "mps_available=false"
      )
    },
    .package = "PredictProR"
  )

  out <- setup_predictdl_env(quiet = TRUE, prefer_gpu = TRUE, cuda = "auto")
  expect_identical(out$python, "C:/conda/envs/predictdl/python.exe")
  expect_true(out$cuda_available)
  expect_identical(out$num_gpus, 1L)
  expect_identical(Sys.getenv("PREDICTPRO_DL_PYTHON"), "C:/conda/envs/predictdl/python.exe")
})

test_that("startup message describes setup helpers as optional provisioning", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  source_file <- file.path(root, "R", "zzz.R")
  if (!file.exists(source_file)) {
    testthat::skip("raw R sources are not available in installed-package test context")
  }

  src <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  expect_true(grepl("optional one-time provisioning helpers", src, fixed = TRUE))
  expect_true(grepl("PREDICTPRO_DL_PYTHON", src, fixed = TRUE))
  expect_true(grepl("PREDICTPRO_GP_PYTHON", src, fixed = TRUE))
  expect_false(grepl("Run `setup_predictdl_env(prefer_gpu = TRUE, cuda = 'auto')` only when you need", src, fixed = TRUE))
})
