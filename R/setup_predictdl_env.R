# Create a single env with both stacks
setup_predictpro_env <- function(
    #env_name = "predictpro",
    env_name = "predictdl",
    python_version = "3.10",
    prefer_gpu = TRUE,
    add_metal_on_macos = TRUE,
    add_cuda_on_linux = FALSE
) {
  # 1) Ensure conda
  has_conda <- tryCatch({ reticulate::conda_binary(); TRUE }, error = function(e) FALSE)
  if (!has_conda) reticulate::install_miniconda()

  # 2) Create env if missing
  envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
  if (!(env_name %in% envs)) {
    message("Creating conda env '", env_name, "' with Python ", python_version, " …")
    reticulate::conda_create(envname = env_name, packages = paste0("python=", python_version))
  } else {
    message("Using existing conda env '", env_name, "'.")
  }

  # 3) Activate for this session
  reticulate::use_condaenv(env_name, required = TRUE)

  # 4) Install both stacks (CPU-safe everywhere)
  base_pkgs <- c(
    # numerics
    "numpy>=1.23,<2.0", "scipy>=1.10", "pandas>=1.5",
    # JAX stack
    "jax>=0.4.20", "jaxlib", "numpyro>=0.13.2", "pymc>=5.10", "arviz>=0.16",
    "xarray>=2024.1", "xarray-einstats>=0.6",
    # DL stack
    "torch", "torchvision", "torchaudio"
  )

  # Optional accelerators (best-effort; safe if wheels are missing)
  sysname <- Sys.info()[["sysname"]]
  arch <- R.version$arch
  if (identical(sysname, "Darwin") && grepl("arm64", arch) && isTRUE(add_metal_on_macos)) {
    base_pkgs <- c(base_pkgs, "jax-metal>=0.0.5")
  }
  if (identical(sysname, "Linux") && isTRUE(add_cuda_on_linux)) {
    # You can install CUDA JAX via pip extra, but reticulate::py_install won't accept extras cleanly per-item.
    # Leave JAX CPU for reliability; use CUDA in power users’ machines via your Python setup_deps if desired.
  }

  reticulate::py_install(base_pkgs, pip = TRUE)
  invisible(reticulate::py_config())
}



# setup_predictdl_env <- function(
#     env_name = "predictdl",
#     python_version = "3.11",
#     prefer_gpu = TRUE,
#     cuda = c("auto", "cpu", "cu121", "cu118"),
#     torch_version = NULL,
#     index_url = NULL
# ) {
#   cuda <- match.arg(cuda)
#
#   # 1) Ensure Conda
#   has_conda <- tryCatch({ reticulate::conda_binary(); TRUE }, error = function(e) FALSE)
#   if (!has_conda) {
#     message("Conda not found. Installing Miniconda (one time)…")
#     reticulate::install_miniconda()
#   }
#
#   # 2) Create env if missing
#   envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
#   if (!(env_name %in% envs)) {
#     message("Creating conda env '", env_name, "' with Python ", python_version, " …")
#     reticulate::conda_create(envname = env_name, packages = paste0("python=", python_version))
#   } else {
#     message("Using existing conda env '", env_name, "'.")
#   }
#
#   # 3) Activate env for this R session
#   reticulate::use_condaenv(env_name, required = TRUE)
#
#   # 4) Load our Python module properly and call setup_deps()
#   mod <- get_dl_module()
#   info <- mod$setup_deps(
#     prefer_gpu = isTRUE(prefer_gpu),
#     cuda = cuda,
#     torch_version = torch_version,
#     numpy_spec = "numpy>=1.24",
#     extra_packages = c("torchvision", "torchaudio"),
#     index_url = index_url
#   )
#
#   cat(sprintf(
#     "\n✔ Environment ready\n  - Python: %s\n  - torch: %s (compiled CUDA: %s)\n  - CUDA available: %s | GPU: %s | MPS: %s\n",
#     info$python, info$pytorch_version,
#     ifelse(is.null(info$cuda_compiled), "NA", info$cuda_compiled),
#     info$cuda_available, ifelse(is.null(info$gpu_name), "None", info$gpu_name),
#     info$mps_available
#   ))
#   invisible(info)
# }
