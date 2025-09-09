# R/ba_wrapper.R
# Minimal Bayes Alphabet wrapper (NumPyro/PyMC/RKHS), no DL refs.

ba_env <- new.env(parent = emptyenv())   # cache for bayes_alphabet_bglr_parity.py

#---- 1) Pin Python env + set JAX/BLAS flags ------------------------------------

# Call early (e.g., in .onLoad or user script) to select a Python env.
# ba_use_env <- function(env = "~/.virtualenvs/r-reticulate", conda = FALSE) {
#   if (conda) reticulate::use_condaenv(env, required = TRUE)
#   else       reticulate::use_virtualenv(env, required = TRUE)
#   invisible(TRUE)
# }

# Set JAX / XLA / BLAS thread knobs BEFORE importing Python.
ba_set_env <- function(
    jax_x64 = TRUE, preallocate = FALSE, mem_fraction = 0.85,
    omp_threads = 1L, mkl_threads = 1L, numexpr_threads = 1L
) {
  Sys.setenv(
    XLA_PYTHON_CLIENT_PREALLOCATE = if (isTRUE(preallocate)) "true" else "false",
    XLA_PYTHON_CLIENT_MEM_FRACTION = sprintf("%.2f", mem_fraction),
    XLA_PYTHON_CLIENT_ALLOCATOR = "platform",
    OMP_NUM_THREADS = as.integer(omp_threads),
    MKL_NUM_THREADS = as.integer(mkl_threads),
    NUMEXPR_NUM_THREADS = as.integer(numexpr_threads),
    JAX_ENABLE_X64 = if (isTRUE(jax_x64)) "1" else "0"
  )
  invisible(TRUE)
}

#---- 2) Import helpers ----------------------------------------------------------

ba_find_pyfile <- function(fname) {
  pkg <- tryCatch(utils::packageName(), error = function(e) "")
  candidates <- c(
    file.path("inst", "python", fname),
    file.path("python", fname),
    file.path("..", "inst", "python", fname)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(normalizePath(hit[1], winslash = "/"))
  if (nzchar(pkg)) {
    sys <- system.file(file.path("python", fname), package = pkg)
    if (nzchar(sys)) return(sys)
  }
  ""
}

get_bayes_module <- function(convert = TRUE) {
  if (!is.null(ba_env$mod) && isTRUE(ba_env$convert) == convert) return(ba_env$mod)
  if (exists("ba_set_env", mode = "function")) ba_set_env()

  pyfile <- ba_find_pyfile("bayes_alphabet_bglr_parity.py")
  if (!nzchar(pyfile) || !file.exists(pyfile))
    stop("bayes_alphabet_bglr_parity.py not found under inst/python/.", call. = FALSE)

  mod <- reticulate::import_from_path(
    "bayes_alphabet_bglr_parity",
    path = dirname(pyfile),
    delay_load = FALSE,
    convert = convert
  )
  ba_env$mod <- mod
  ba_env$convert <- convert
  mod
}

ba_reload_bayes_module <- function(convert = TRUE) {
  reticulate::py_run_string("import sys; sys.modules.pop('bayes_alphabet_bglr_parity', None)")
  ba_env$mod <- NULL; ba_env$convert <- NULL
  get_bayes_module(convert = convert)
}

get_perf_harness <- function(name = "stress_harness_bayes_alphabet.py", convert = TRUE) {
  pyfile <- ba_find_pyfile(name)
  if (!nzchar(pyfile)) return(NULL)
  reticulate::import_from_path(sub("\\.py$", "", name), path = dirname(pyfile),
                               delay_load = FALSE, convert = convert)
}

#---- 3) Environment bootstrap helpers ------------------------------------------

ba_init_python <- function(
    env = c("auto", "virtualenv", "conda", "system"),
    name = "r-bayes",
    allow_create = FALSE,
    allow_install = FALSE,
    prefer_gpu = TRUE,
    jax_x64 = TRUE, preallocate = FALSE, mem_fraction = 0.85,
    omp_threads = 1L, mkl_threads = 1L, numexpr_threads = 1L
) {
  env <- match.arg(env)
  ba_set_env(jax_x64, preallocate, mem_fraction, omp_threads, mkl_threads, numexpr_threads)

  rp <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(rp) && file.exists(rp)) return(invisible(reticulate::py_config()))

  try_use_conda_env <- function(name) {
    cl <- try(reticulate::conda_list(), silent = TRUE)
    if (inherits(cl, "try-error") || !is.data.frame(cl)) return(FALSE)
    nmcol <- intersect(c("name", "envname"), names(cl))
    has_env <- length(nmcol) > 0 && any(cl[[nmcol[1]]] == name)
    if (has_env) {
      ok <- !inherits(try(reticulate::use_condaenv(name, required = TRUE), silent = TRUE), "try-error")
      if (ok && isTRUE(allow_install)) ba_install_pydeps(prefer_gpu = prefer_gpu)
      return(ok)
    }
    if (isTRUE(allow_create) && env %in% c("auto","conda")) {
      ok <- !inherits(try(reticulate::conda_create(name, packages = "python=3.10"), silent = TRUE), "try-error")
      if (ok) {
        reticulate::use_condaenv(name, required = TRUE)
        if (isTRUE(allow_install)) ba_install_pydeps(prefer_gpu = prefer_gpu)
      }
      return(ok)
    }
    FALSE
  }

  try_use_virtualenv <- function(name) {
    ok <- !inherits(try(reticulate::use_virtualenv(name, required = TRUE), silent = TRUE), "try-error")
    if (ok) {
      if (isTRUE(allow_install)) ba_install_pydeps(prefer_gpu = prefer_gpu)
      return(TRUE)
    }
    if (isTRUE(allow_create) && env %in% c("auto","virtualenv")) {
      py <- Sys.which(if (.Platform$OS.type == "windows") "python" else "python3")
      if (!nzchar(py)) py <- "python"
      ok <- !inherits(try(reticulate::virtualenv_create(envname = name, python = py), silent = TRUE), "try-error")
      if (ok) {
        reticulate::use_virtualenv(name, required = TRUE)
        if (isTRUE(allow_install)) ba_install_pydeps(prefer_gpu = prefer_gpu)
      }
      return(ok)
    }
    FALSE
  }

  ok <- FALSE
  if (env == "conda")         ok <- try_use_conda_env(name)
  else if (env == "virtualenv") ok <- try_use_virtualenv(name)
  else if (env == "auto")     ok <- try_use_conda_env(name) || try_use_virtualenv(name)

  if (!ok && env != "system")
    warning("Falling back to system Python (no '", name, "' env found/created).")

  invisible(reticulate::py_config())
}

ba_install_pydeps <- function(prefer_gpu = TRUE) {
  pkgs <- c(
    "numpy>=1.23,<2.0", "scipy>=1.10", "pandas>=1.5",
    "jax>=0.4.20", "numpyro>=0.13.2", "pymc>=5.10", "arviz>=0.16",
    "xarray>=2024.1", "xarray-einstats>=0.6"
  )
  # Optional GPU hints (Linux CUDA / macOS Metal) can be appended here.
  reticulate::py_install(packages = pkgs, pip = TRUE)
  invisible(TRUE)
}

setup_bayes_env <- function(
    env_name = "r-bayes",
    python_version = "3.10",
    prefer_gpu = TRUE,
    cuda = c("auto", "cpu")
) {
  cuda <- match.arg(cuda)
  has_conda <- tryCatch({ reticulate::conda_binary(); TRUE }, error = function(e) FALSE)
  if (!has_conda) {
    message("Conda not found. Installing Miniconda (one time)…")
    reticulate::install_miniconda()
  }
  envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
  if (!(env_name %in% envs)) {
    message("Creating conda env '", env_name, "' with Python ", python_version, " …")
    reticulate::conda_create(envname = env_name, packages = paste0("python=", python_version))
  } else {
    message("Using existing conda env '", env_name, "'.")
  }
  reticulate::use_condaenv(env_name, required = TRUE)
  ba_set_env()

  mod <- get_bayes_module()
  info <- mod$setup_deps(
    prefer_gpu = isTRUE(prefer_gpu),
    try_numpyro = TRUE,
    try_pymc = TRUE,
    allow_install = TRUE,
    cuda = cuda
  )

  `%||%` <- function(a, b) if (is.null(a)) b else a
  cat(sprintf(
    "\n✔ Bayes env ready\n  - Python: %s\n  - jax: %s | jaxlib: %s | backend: %s | devices: %s\n  - numpyro: %s | pymc: %s\n",
    info$python,
    info$jax_version, info$jaxlib_version, info$default_backend %||% "NA",
    if (length(info$devices)) paste(info$devices, collapse = ",") else "none",
    info$numpyro_version %||% "NA", info$pymc_version %||% "NA"
  ))
  if (!is.null(info$note)) message(info$note)
  invisible(info)
}

# Ensure deps without mutating env unless allow_install=TRUE.
ensure_pydeps_bayes <- function(
    prefer_gpu = TRUE,
    cuda = c("auto", "cpu"),
    try_numpyro = TRUE,
    try_pymc = TRUE,
    allow_install = FALSE,
    use_python_setup = FALSE
) {
  cuda <- match.arg(cuda)
  cfg <- reticulate::py_config()
  ba_set_env()

  preflight <- c("arviz", "jax", "jaxlib")
  if (isTRUE(try_numpyro)) preflight <- c(preflight, "numpyro")
  if (isTRUE(try_pymc))    preflight <- c(preflight, "pymc")

  missing <- preflight[!vapply(preflight, reticulate::py_module_available, logical(1))]
  if (length(missing)) {
    if (isTRUE(allow_install)) {
      reticulate::py_install(missing, pip = TRUE)
    } else {
      stop(
        "Missing Python modules in the active env (", cfg$python, "): ",
        paste(missing, collapse = ", "),
        "\nInstall with reticulate::py_install(...), or run allow_install=TRUE once."
      )
    }
  }

  if (!isTRUE(use_python_setup)) {
    jax <- reticulate::import("jax", convert = TRUE)
    jaxlib <- reticulate::import("jaxlib", convert = TRUE)
    return(invisible(list(
      active_python = cfg$python, ok = TRUE, setup = "skipped",
      jax = jax$`__version__`, jaxlib = jaxlib$`__version__`
    )))
  }

  mod <- get_bayes_module(convert = TRUE)
  if (!reticulate::py_has_attr(mod, "setup_deps"))
    return(invisible(list(active_python = cfg$python, ok = TRUE, setup = "no-op")))

  attempts <- list(
    list(prefer_gpu = isTRUE(prefer_gpu), try_numpyro = isTRUE(try_numpyro),
         try_pymc = isTRUE(try_pymc), allow_install = isTRUE(allow_install), cuda = cuda),
    list(prefer_gpu = isTRUE(prefer_gpu), cuda = cuda),
    list(prefer_gpu = isTRUE(prefer_gpu)), list()
  )
  last_err <- NULL
  for (args in attempts) {
    res <- try(do.call(mod$setup_deps, args), silent = TRUE)
    if (!inherits(res, "try-error")) return(invisible(c(list(active_python = cfg$python), res)))
    last_err <- res
  }
  stop("setup_deps() failed: ", as.character(last_err))
}

#---- 4) Fit wrapper -------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

ba_fit <- function(ETA, y,
                   backend = c("numpyro","pymc"),
                   vi = FALSE,
                   n_iter = 2000L, tune = 1000L, chains = 2L, seed = 123L,
                   target_accept = 0.9, ...) {
  stopifnot(is.list(ETA), is.numeric(y))
  ba_set_env()
  mod <- get_bayes_module(convert = TRUE)
  backend <- match.arg(backend)

  res <- mod$bayesian_alphabet(
    ETA = ETA,
    y = as.numeric(y),
    backend = backend,
    vi = isTRUE(vi),
    n_iter = as.integer(n_iter),
    tune = as.integer(tune),
    chains = as.integer(chains),
    random_seed = as.integer(seed),
    target_accept = as.numeric(target_accept),
    allow_install = FALSE,
    disable_auto_install = TRUE,
    ...
  )
  class(res) <- c("ba_fit", class(res))
  res
}

#---- 5) S3 helpers --------------------------------------------------------------

#' @export
print.ba_fit <- function(x, ...) {
  cat("Bayes Alphabet fit (backend:", x$extra$backend, "| inference:", x$extra$inference, ")\n")
  cat(sprintf("mu=%.4f  varE=%.4f  varG=%.4f  h2=%.4f\n", x$mu, x$varE, x$varG, x$h2))
  if (!is.null(x$ETA) && length(x$ETA)) {
    for (i in seq_along(x$ETA)) {
      eff <- x$ETA[[i]]
      line <- sprintf("  ETA[%d] %-12s | p=%5d  varU=%.4f", i-1L, eff$model, eff$p, eff$varU %||% NA_real_)
      if (!is.null(eff$pi)) line <- paste0(line, sprintf("  pi=%.3f", eff$pi))
      cat(line, "\n")
    }
  }
  invisible(x)
}

#' @export
summary.ba_fit <- function(object, ...) {
  out <- list(
    backend = object$extra$backend,
    inference = object$extra$inference,
    mu = object$mu,
    varE = object$varE,
    varG = object$varG,
    h2 = object$h2,
    n_eta = length(object$ETA),
    diagnostics = object$extra$diagnostics
  )
  class(out) <- "summary_ba_fit"
  out
}

#' @export
print.summary_ba_fit <- function(x, ...) {
  cat("Bayes Alphabet summary\n")
  cat(" backend   :", x$backend, "\n")
  cat(" inference :", x$inference, "\n")
  cat(sprintf(" mu=%.4f  varE=%.4f  varG=%.4f  h2=%.4f\n", x$mu, x$varE, x$varG, x$h2))
  if (!is.null(x$diagnostics)) {
    cat(" diagnostics keys:", paste(names(x$diagnostics), collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
coef.ba_fit <- function(object, effect = 0L, ...) {
  eff <- object$ETA[[effect + 1L]]
  if (is.null(eff) || is.null(eff$b)) stop("Effect not found or has no coefficients.")
  eff$b
}

#' @export
fitted.ba_fit <- function(object, ...) object$yHat

#' @export
predict.ba_fit <- function(object, se.fit = FALSE, ...) {
  if (isTRUE(se.fit)) return(list(fit = object$yHat, se.fit = object$`SD.yHat`))
  object$yHat
}

#---- 6) ETA builders ------------------------------------------------------------

ba_eta_brr         <- function(X) list(X = as.matrix(X), method = "bayesian_ridge")
ba_eta_bayesa      <- function(X, df = 5, varB = NULL) { eta <- list(X = as.matrix(X), method = "bayes_a", df = as.numeric(df)); if (!is.null(varB)) eta$varB <- varB; eta }
ba_eta_bayesb      <- function(X, d = NULL) { eta <- list(X = as.matrix(X), method = "bayes_b"); if (!is.null(d)) eta$d <- as.numeric(d); eta }
ba_eta_bayesc      <- function(X, d = NULL) { eta <- list(X = as.matrix(X), method = "bayes_c"); if (!is.null(d)) eta$d <- as.numeric(d); eta }
ba_eta_bayescpi    <- function(X, fixed_pi = NULL, d = NULL) { eta <- list(X = as.matrix(X), method = "bayes_c_pi"); if (!is.null(fixed_pi)) eta$fixed_pi <- as.numeric(fixed_pi); if (!is.null(d)) eta$d <- as.numeric(d); eta }
ba_eta_lasso       <- function(X) list(X = as.matrix(X), method = "lasso")
ba_eta_elastic_net <- function(X, alpha = 0.5) list(X = as.matrix(X), method = "elastic_net", alpha = as.numeric(alpha))
ba_eta_rkhs_K      <- function(K) list(method = "rkhs", K = as.matrix(K))
ba_eta_rkhs_X      <- function(X) list(method = "rkhs", X = as.matrix(X))

ba_eta_stack <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1]]) && !is.matrix(args[[1]])) return(args[[1]])
  args
}

#---- 7) Convenience top-levels --------------------------------------------------

ba_fit_brr_bayesb <- function(X1, X2, y, backend = "numpyro", vi = FALSE, ...) {
  ETA <- ba_eta_stack(ba_eta_brr(X1), ba_eta_bayesb(X2))
  ba_fit(ETA, y, backend = backend, vi = vi, ...)
}
ba_fit_bayesa <- function(X, y, backend = "numpyro", vi = FALSE, ...) {
  ba_fit(ba_eta_stack(ba_eta_bayesa(X)), y, backend = backend, vi = vi, ...)
}
ba_fit_bayesa_cpi <- function(X1, X2, y, fixed_pi = 0.1, d = NULL,
                              backend = "numpyro", vi = TRUE, ...) {
  ETA <- ba_eta_stack(ba_eta_bayesa(X1), ba_eta_bayescpi(X2, fixed_pi = fixed_pi, d = d))
  ba_fit(ETA, y, backend = backend, vi = vi, ...)
}
ba_fit_rkhs_K <- function(K, y, backend = "numpyro", vi = FALSE, ...) {
  ba_fit(ba_eta_stack(ba_eta_rkhs_K(K)), y, backend = backend, vi = vi, ...)
}
ba_fit_rkhs_X <- function(X, y, backend = "numpyro", vi = FALSE, ...) {
  ba_fit(ba_eta_stack(ba_eta_rkhs_X(X)), y, backend = backend, vi = vi, ...)
}

#---- 8) Device info & quick JAX microbench -------------------------------------

ba_device_info <- function() {
  ba_set_env()
  jax <- tryCatch({
    mod <- get_bayes_module(convert = TRUE)
    if (reticulate::py_has_attr(mod, "jax")) mod$jax else reticulate::import("jax", convert = TRUE)
  }, error = function(e) NULL)
  if (is.null(jax)) return(list(jax_version = NA_character_, devices = character(), default_backend = NA_character_))
  devs <- try(jax$devices(), silent = TRUE)
  ds <- if (inherits(devs, "try-error")) character() else vapply(devs, function(d) reticulate::py_get_attr(d, "platform"), "")
  list(
    jax_version = tryCatch(jax$`__version__`, error = function(e) NA_character_),
    default_backend = tryCatch(jax$default_backend(), error = function(e) NA_character_),
    devices = ds
  )
}

ba_gpu_microbench <- function(n = 2048L, dtype = c("float32","float64")) {
  dtype <- match.arg(dtype)
  ba_set_env()
  jax <- reticulate::import("jax", convert = TRUE)
  jnp <- reticulate::import("jax.numpy", convert = TRUE)

  key <- jax$random$PRNGKey(0L)
  a <- jax$random$normal(key, c(as.integer(n), as.integer(n)), dtype = jnp[[dtype]])
  b <- jax$random$normal(key, c(as.integer(n), as.integer(n)), dtype = jnp[[dtype]])
  mm <- reticulate::py_eval("(__import__('jax').jit(lambda x, y: x @ y))")

  t0 <- proc.time()[["elapsed"]]; out <- mm(a, b)$block_until_ready()
  t1 <- proc.time()[["elapsed"]]; out <- mm(a, b)$block_until_ready()
  t2 <- proc.time()[["elapsed"]]

  list(
    n = n, dtype = dtype,
    compile_s = t1 - t0,
    exec_s    = t2 - t1,
    gflops    = (2.0 * n * n * n) / max(1e-9, (t2 - t1)) / 1e9
  )
}

#---- 9) Optional harness runner -------------------------------------------------

run_perf_harness <- function(module = c("stress_harness_bayes_alphabet.py", "rkhs_harness.py"),
                             args = character()) {
  module <- match.arg(module)
  ba_set_env()
  h <- get_perf_harness(name = module, convert = TRUE)
  if (is.null(h)) stop(module, " not found under inst/python/. Add a module with main(argv).", call. = FALSE)
  if (!reticulate::py_has_attr(h, "main"))
    stop("Python harness has no main(argv) function.")
  invisible(h$main(args))
}
