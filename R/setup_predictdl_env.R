gp_conda_env_python <- function(env_name) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(NULL)
  }

  envs <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.null(envs) || !nrow(envs)) {
    return(NULL)
  }

  name_col <- intersect(c("name", "envname"), names(envs))
  py_col <- intersect(c("python", "python_path"), names(envs))
  if (!length(name_col) || !length(py_col)) {
    return(NULL)
  }

  matched <- envs[envs[[name_col[[1L]]]] %in% env_name, , drop = FALSE]
  if (!nrow(matched)) {
    return(NULL)
  }

  py <- matched[[py_col[[1L]]]][1]
  if (is.null(py) || !nzchar(py)) {
    return(NULL)
  }

  normalizePath(py, winslash = "/", mustWork = FALSE)
}

#' Create or update the package-owned Python environment for deep learning
#'
#' @param env_name Conda environment name. Defaults to `predictdl`.
#' @param python_version Python version to install when creating the environment.
#' @param prefer_gpu Logical; if `TRUE`, install a GPU-capable Torch build when possible.
#' @param cuda Torch build selector. One of `"auto"`, `"cpu"`, `"cu121"`, `"cu118"`.
#' @param torch_version Optional Torch version passed through to the Python installer.
#' @param index_url Optional custom Python package index URL for Torch wheels.
#' @param quiet Logical; if `TRUE`, suppress routine status messages.
#'
#' @return A named list describing the configured environment and Torch runtime.
#' @export
setup_predictdl_env <- function(
    env_name = "predictdl",
    python_version = "3.11",
    prefer_gpu = TRUE,
    cuda = c("auto", "cpu", "cu121", "cu118"),
    torch_version = NULL,
    index_url = NULL,
    quiet = FALSE
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to setup the predictdl environment.", call. = FALSE)
  }

  cuda <- match.arg(cuda)

  has_conda <- tryCatch({
    nzchar(reticulate::conda_binary())
  }, error = function(e) FALSE)
  if (!has_conda) {
    if (!quiet) {
      message("Conda not found. Installing Miniconda.")
    }
    reticulate::install_miniconda()
  }

  envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
  if (!(env_name %in% envs)) {
    if (!quiet) {
      message("Creating conda env '", env_name, "' with Python ", python_version, ".")
    }
    reticulate::conda_create(envname = env_name, packages = paste0("python=", python_version))
  } else {
    if (!quiet) {
      message("Using existing conda env '", env_name, "'.")
    }
  }

  env_python <- gp_conda_env_python(env_name)
  if (!is.null(env_python) && nzchar(env_python)) {
    Sys.setenv(PREDICTPRO_DL_PYTHON = env_python, PREDICTPRO_PYTHON_DL = env_python)
  }

  setup_args <- c(
    if (isTRUE(prefer_gpu)) "--prefer-gpu" else "--no-prefer-gpu",
    "--cuda", cuda
  )
  if (!is.null(torch_version) && nzchar(torch_version)) {
    setup_args <- c(setup_args, "--torch-version", as.character(torch_version)[1L])
  }
  if (!is.null(index_url) && nzchar(index_url)) {
    setup_args <- c(setup_args, "--index-url", as.character(index_url)[1L])
  }

  bridge_out <- gp_dl_run_cli("setup-deps", setup_args, python_bin = env_python)
  info <- gp_parse_cli_key_values(bridge_out)

  info_python <- info[["python"]] %||% NULL
  if (!is.null(info_python) && nzchar(info_python) && file.exists(info_python)) {
    env_python <- info_python
  }
  if (!is.null(env_python) && nzchar(env_python)) {
    env_python <- normalizePath(env_python, winslash = "/", mustWork = FALSE)
    Sys.setenv(PREDICTPRO_DL_PYTHON = env_python, PREDICTPRO_PYTHON_DL = env_python)
  }

  out <- list(
    env_name = env_name,
    python = env_python,
    pytorch_version = info[["pytorch_version"]] %||% NULL,
    cuda_compiled = info[["cuda_compiled"]] %||% NULL,
    cuda_available = gp_cli_bool(info[["cuda_available"]]),
    gpu_name = info[["gpu_name"]] %||% NULL,
    mps_available = gp_cli_bool(info[["mps_available"]]),
    reticulate_initialized = FALSE,
    torch_available = !is.null(info[["pytorch_version"]]) && nzchar(info[["pytorch_version"]]),
    num_gpus = gp_cli_int(info[["cuda_device_count"]] %||% info[["num_gpus"]])
  )

  if (!quiet) {
    cat(sprintf(
      "\nEnvironment ready\n  - env: %s\n  - Python: %s\n  - torch: %s (compiled CUDA: %s)\n  - CUDA available: %s | GPUs: %s | GPU: %s | MPS: %s\n",
      out$env_name,
      out$python %||% "NA",
      out$pytorch_version %||% "NA",
      out$cuda_compiled %||% "NA",
      out$cuda_available,
      out$num_gpus %||% 0L,
      out$gpu_name %||% "None",
      out$mps_available
    ))
  }

  invisible(out)
}
