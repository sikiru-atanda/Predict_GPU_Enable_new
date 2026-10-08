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

#' Create or update the package-owned Python environment for Gaussian-process models
#'
#' @param env_name Conda environment name. Defaults to `predictgp`.
#' @param python_version Python version to install when creating the environment.
#' @param prefer_gpu Logical; if `TRUE`, install/use a CUDA-capable PyTorch
#' build when an NVIDIA GPU is available.
#' @param cuda Torch CUDA wheel selector. One of `"auto"`, `"cpu"`,
#' `"cu121"`, or `"cu118"`.
#' @param torch_version Optional Torch version passed through to the Python
#' installer.
#' @param index_url Optional custom Python package index URL for Torch wheels.
#' @param install_scikit_sparse Logical; if `TRUE` (default), additionally
#' attempt a cross-OS install of `scikit-sparse` so the new sparse Mixed Model
#' Equations REML engine (`gp_engine = "mme"`) can use the CHOLMOD backend
#' (~2x faster than the `scipy.splu` fallback). The install is best-effort: it
#' tries conda-forge first (the only supported path on Windows, also works on
#' Linux/macOS) and falls back to pip on Linux/macOS where SuiteSparse is
#' system-installed. Failures are non-fatal -- the MME engine transparently
#' falls back to `scipy.sparse.linalg.splu`.
#' @param quiet Logical; if `TRUE`, suppress routine status messages.
#'
#' @return A named list describing the configured environment and GP runtime,
#' including the resolved MME factor backend (`mme_factor_backend` is either
#' `"cholmod"` or `"splu"`) and the `scikit-sparse` version if installed.
#' @export
setup_predictgp_env <- function(
    env_name = "predictgp",
    python_version = "3.11",
    prefer_gpu = TRUE,
    cuda = c("auto", "cpu", "cu121", "cu118"),
    torch_version = NULL,
    index_url = NULL,
    install_scikit_sparse = TRUE,
    quiet = FALSE
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to setup the predictgp environment.", call. = FALSE)
  }

  cuda <- match.arg(cuda)

  has_conda <- tryCatch(nzchar(reticulate::conda_binary()), error = function(e) FALSE)
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
  } else if (!quiet) {
    message("Using existing conda env '", env_name, "'.")
  }

  env_python <- gp_conda_env_python(env_name)
  if (!is.null(env_python) && nzchar(env_python)) {
    Sys.setenv(PREDICTPRO_GP_PYTHON = env_python)
  }

  setup_args <- c(
    "--project-root", gp_project_root(),
    if (isTRUE(prefer_gpu)) "--prefer-gpu" else "--no-prefer-gpu",
    "--cuda", cuda,
    if (isTRUE(install_scikit_sparse)) "--install-scikit-sparse" else "--no-install-scikit-sparse"
  )
  if (!is.null(torch_version) && nzchar(torch_version)) {
    setup_args <- c(setup_args, "--torch-version", as.character(torch_version)[1L])
  }
  if (!is.null(index_url) && nzchar(index_url)) {
    setup_args <- c(setup_args, "--index-url", as.character(index_url)[1L])
  }
  bridge_out <- gp_bridge_run_cli("setup-deps", setup_args, python_bin = env_python)
  info <- gp_parse_cli_key_values(bridge_out)
  if (!is.null(info[["python"]]) && nzchar(info[["python"]]) && file.exists(info[["python"]])) {
    env_python <- normalizePath(info[["python"]], winslash = "/", mustWork = FALSE)
    Sys.setenv(PREDICTPRO_GP_PYTHON = env_python)
  }

  out <- list(
    env_name = env_name,
    python = env_python,
    prefer_gpu = isTRUE(prefer_gpu),
    cuda_requested = cuda,
    torch_install_flavor = info[["torch_install_flavor"]] %||% NULL,
    torch_index_used = info[["torch_index_used"]] %||% NULL,
    cuda_available = gp_cli_bool(info[["cuda_available"]]),
    cuda_compiled = info[["cuda_compiled"]] %||% NULL,
    num_gpus = gp_cli_int(info[["cuda_device_count"]]),
    gpu_name = info[["gpu_name"]] %||% NULL,
    torch_version = info[["torch_version"]] %||% NULL,
    gpytorch_version = info[["gpytorch_version"]] %||% NULL,
    scikit_sparse_version = info[["scikit_sparse_version"]] %||% NULL,
    mme_factor_backend = info[["mme_factor_backend"]] %||% "splu"
  )

  if (!quiet) {
    cat(sprintf(
      "\nEnvironment ready\n  - env: %s\n  - Python: %s\n  - torch: %s (compiled CUDA: %s; install: %s)\n  - gpytorch: %s\n  - CUDA available: %s | GPUs: %s | GPU: %s\n  - MME factor backend: %s (scikit-sparse: %s)\n",
      out$env_name,
      out$python %||% "NA",
      out$torch_version %||% "NA",
      out$cuda_compiled %||% "NA",
      out$torch_install_flavor %||% "NA",
      out$gpytorch_version %||% "NA",
      out$cuda_available,
      out$num_gpus %||% 0L,
      out$gpu_name %||% "None",
      out$mme_factor_backend,
      out$scikit_sparse_version %||% "not installed"
    ))
  }

  invisible(out)
}
