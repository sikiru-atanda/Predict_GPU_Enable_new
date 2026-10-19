# Python packages for PredictProR's Python backends, installed on first use.
#
# inst/python/predictpror_deps.py knows the packages and installs them (with a
# lock on the environment, so parallel workers never run pip at once). R
# decides which dependency groups a run needs, installs them once in the main
# session before any parallel worker starts (model_execute() calls
# gp_python_prepare_runtimes()), and hands the workers the interpreters and a
# list of what was verified (PREDICTPRO_PYTHON_DEPS_OK, propagated by
# gp_parallel_runtime_env()). The bridge runners check again, which is free
# for anything on that list. An interpreter that cannot take installs (not
# writable, or an OS-managed system Python) is replaced by a PredictProR
# virtualenv built from it.
#
# PREDICTPRO_AUTO_INSTALL=false turns installing off; missing packages then
# stop with the list of what is needed.

gp_python_deps_script_path <- function() {
  pyfile <- system.file("python/predictpror_deps.py", package = "PredictProR")
  if (!nzchar(pyfile)) {
    candidates <- c("inst/python/predictpror_deps.py", "python/predictpror_deps.py",
                    "../inst/python/predictpror_deps.py")
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) pyfile <- normalizePath(hit[1L], winslash = "/", mustWork = TRUE)
  }
  if (!nzchar(pyfile)) {
    stop("predictpror_deps.py not found under inst/python/.", call. = FALSE)
  }
  pyfile
}

gp_python_auto_install_enabled <- function() {
  value <- tolower(trimws(Sys.getenv("PREDICTPRO_AUTO_INSTALL", unset = "")))
  !value %in% c("0", "false", "no", "off")
}

gp_python_purpose_label <- function(purpose) {
  switch(purpose, ml = "machine-learning", dl = "deep-learning", gp = "Gaussian-process", purpose)
}

# Groups every bridge of a purpose needs; model-specific packages (xgboost,
# catboost, lightgbm, gpytorch) are added by gp_python_dependency_groups().
gp_python_base_groups <- function(purpose) {
  switch(purpose,
         ml = c("tabular", "sklearn"),
         dl = c("tabular", "torch", "sklearn"),
         gp = c("tabular", "torch", "gp"),
         stop("Unknown Python purpose: ", purpose, call. = FALSE))
}

# Python ML models (canonical names) -> extra group; Lasso is fitted in R.
gp_python_ml_model_groups <- function() {
  list(
    Xgboost = "xgboost", RandomForest = character(), CatBoost = "catboost",
    LightGBM = "lightgbm", PartialLeastSquare = character(),
    SupportVectorMachine = character(), `K-NearestNeighbors` = character(),
    Ridge_Regression = character()
  )
}

gp_python_gp_models <- function() {
  c("GP", "KRR", "GP_FA", "LowRankGP")
}

# Dependency groups per purpose for the requested (canonical) model names.
gp_python_dependency_groups <- function(models, gp_backend = "auto", hybrid_gp = FALSE) {
  models <- unique(as.character(models %||% character()))
  out <- list()
  ml_map <- gp_python_ml_model_groups()
  ml_hits <- intersect(models, names(ml_map))
  if (length(ml_hits)) {
    out$ml <- unique(c(gp_python_base_groups("ml"), unlist(ml_map[ml_hits], use.names = FALSE)))
  }
  dl_hits <- intersect(models, gp_deep_learning_supported_models())
  if (length(dl_hits)) {
    out$dl <- unique(c(gp_python_base_groups("dl"), if ("gp_dkl" %in% dl_hits) "gpytorch"))
  }
  gp_hits <- intersect(models, gp_python_gp_models())
  # hybrid_gp with gp_backend "auto"/"r" falls back to the R solver when the
  # Python backend is unavailable (documented), so it installs nothing
  if (length(gp_hits) && !(isTRUE(hybrid_gp) && !gp_hybrid_gp_require_python_backend(gp_backend))) {
    out$gp <- gp_python_base_groups("gp")
  }
  out
}

gp_python_detect_for <- function(purpose) {
  switch(purpose,
         ml = gp_detect_ml_python(),
         dl = gp_detect_dl_python(),
         gp = gp_detect_gp_python())
}

gp_python_no_runtime_message <- function(purpose) {
  switch(purpose,
         ml = paste0("No Python interpreter was found for Python ML models. Install Python (>= 3.10) ",
                     "or set PREDICTPRO_ML_PYTHON / PREDICTPRO_PYTHON to one; PredictProR installs ",
                     "numpy, pandas, scikit-learn and the model packages into it."),
         dl = paste0("No configured Python runtime was found for deep-learning models.\n",
                     "Set PREDICTPRO_DL_PYTHON to a Python interpreter (>= 3.10); PredictProR installs ",
                     "PyTorch and the other DL packages into it.\n",
                     "You can also run PredictProR::setup_predictdl_env(prefer_gpu = TRUE, cuda = 'auto') ",
                     "as a one-time provisioning helper."),
         gp = paste0("No configured Python runtime was found for Gaussian-process models.\n",
                     "Run setup_predictgp_env() or set PREDICTPRO_GP_PYTHON to a Python interpreter ",
                     "(>= 3.10); PredictProR installs PyTorch and GPyTorch into it."))
}

gp_python_set_purpose_python <- function(purpose, python_bin) {
  switch(purpose,
         ml = Sys.setenv(PREDICTPRO_ML_PYTHON = python_bin),
         dl = Sys.setenv(PREDICTPRO_DL_PYTHON = python_bin, PREDICTPRO_PYTHON_DL = python_bin),
         gp = Sys.setenv(PREDICTPRO_GP_PYTHON = python_bin))
  gp_python_forget_cached_probes()
  invisible(python_bin)
}

# Drop cached interpreter probes/choices so detection sees newly installed
# packages. For an interpreter just installed into, its probes on disk go too
# (its site-packages changed, so they would not be used again anyway).
gp_python_forget_cached_probes <- function(python_bin = NULL) {
  keys <- ls(envir = PredictProR_runtime_cache)
  drop <- grep("^(python_probe::|preferred_python::|ml_python::|conda_candidates::)", keys, value = TRUE)
  if (!is.null(python_bin)) {
    norm <- gp_python_path(python_bin, mustWork = FALSE)
    drop <- drop[!startsWith(drop, "python_probe::") | grepl(norm, drop, fixed = TRUE)]
    gp_python_discovery_cache_drop(paste0("python_probe::", norm, "::"))
  }
  if (length(drop)) rm(list = drop, envir = PredictProR_runtime_cache)
  invisible(NULL)
}

gp_python_deps_key <- function(python_bin, groups) {
  paste0(gp_python_path(python_bin, mustWork = FALSE), "::",
         paste(sort(unique(groups)), collapse = ","))
}

# Verified "interpreter::groups" entries; the environment variable carries
# them to parallel workers.
gp_python_verified_entries <- function() {
  ok <- strsplit(Sys.getenv("PREDICTPRO_PYTHON_DEPS_OK", unset = ""), ";", fixed = TRUE)[[1L]]
  ok[nzchar(ok)]
}

# TRUE when these groups (or a superset) were verified for this interpreter.
gp_python_deps_verified <- function(python_bin, groups) {
  path <- gp_python_path(python_bin, mustWork = FALSE)
  for (entry in gp_python_verified_entries()) {
    sep <- regexpr("::", entry, fixed = TRUE)
    if (sep < 1L || !identical(substr(entry, 1L, sep - 1L), path)) next
    have <- strsplit(substring(entry, sep + 2L), ",", fixed = TRUE)[[1L]]
    if (all(groups %in% have)) return(TRUE)
  }
  FALSE
}

gp_python_mark_verified <- function(python_bin, groups) {
  ok <- unique(c(gp_python_verified_entries(), gp_python_deps_key(python_bin, groups)))
  Sys.setenv(PREDICTPRO_PYTHON_DEPS_OK = paste(ok, collapse = ";"))
  invisible(TRUE)
}

# Run predictpror_deps.py; `stream` shows pip's progress as it happens.
# `timeout` (seconds, 0 = none) bounds a run without streaming.
gp_python_deps_run <- function(python_bin, command, groups, stream = FALSE, extra_args = character(),
                               timeout = 0) {
  json_out <- tempfile("predictpror_deps_", fileext = ".json")
  on.exit(unlink(json_out), add = TRUE)
  args <- c(gp_python_deps_script_path(), command, "--groups", paste(groups, collapse = ","),
            "--json-out", json_out, extra_args)
  output <- character()
  if (isTRUE(stream) && requireNamespace("processx", quietly = TRUE)) {
    show <- function(line, proc) message(line)
    res <- tryCatch(
      processx::run(python_bin, args, error_on_status = FALSE, windows_hide_window = TRUE,
                    stdout_line_callback = show, stderr_line_callback = show),
      error = function(e) list(status = 127L, stdout = "", stderr = conditionMessage(e))
    )
    status <- res$status
    if (!identical(as.integer(status), 0L)) output <- c(res$stdout, res$stderr)
  } else if (isTRUE(stream)) {
    status <- suppressWarnings(system2(python_bin, gp_quote_system_args(args), stdout = "", stderr = ""))
  } else {
    output <- suppressWarnings(system2(python_bin, gp_quote_system_args(args), stdout = TRUE, stderr = TRUE,
                                       timeout = timeout))
    status <- attr(output, "status") %||% 0L
  }
  result <- if (file.exists(json_out)) {
    tryCatch(jsonlite::fromJSON(json_out, simplifyVector = TRUE), error = function(e) NULL)
  }
  if (is.null(result)) {
    timed_out <- timeout > 0 && identical(as.integer(status), 124L)
    result <- list(
      status = "error", reason = if (timed_out) "timeout" else "no_result",
      message = if (timed_out) {
        sprintf("%s did not finish the PredictProR package %s within %d s.", python_bin, command, as.integer(timeout))
      } else {
        paste0("Could not run ", python_bin, " (exit status ", status, ").\n",
               paste(utils::tail(output, 40L), collapse = "\n"))
      }
    )
  }
  result
}

gp_python_env_home <- function() {
  home <- Sys.getenv("PREDICTPRO_PYTHON_ENV_HOME", unset = "")
  if (!nzchar(home)) home <- file.path(tools::R_user_dir("PredictProR", which = "cache"), "python-envs")
  normalizePath(home, winslash = "/", mustWork = FALSE)
}

gp_python_venv_python <- function(env_dir) {
  file.path(env_dir, if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
}

gp_python_runs <- function(python_bin) {
  if (!file.exists(python_bin)) return(FALSE)
  out <- tryCatch(suppressWarnings(system2(python_bin, "--version", stdout = TRUE, stderr = TRUE)),
                  error = function(e) structure("", status = 1L))
  identical(as.integer(attr(out, "status") %||% 0L), 0L)
}

# Hold a directory lock (atomic dir.create) while `code` runs; a lock older
# than `stale_minutes` is taken over (its owner died).
gp_python_with_dir_lock <- function(lock_dir, code, wait_minutes = 30, stale_minutes = 60) {
  deadline <- Sys.time() + wait_minutes * 60
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) break
    info <- file.info(lock_dir)
    if (!is.na(info$mtime) && difftime(Sys.time(), info$mtime, units = "mins") > stale_minutes) {
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }
    if (Sys.time() > deadline) {
      stop("Timed out waiting for another R process to finish setting up ", dirname(lock_dir),
           ". If none is running, delete ", lock_dir, ".", call. = FALSE)
    }
    Sys.sleep(2)
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  force(code)
}

# A PredictProR virtualenv built from `base_python` (which cannot take installs
# itself); it sees the base interpreter's packages and installs the rest.
gp_python_managed_env <- function(base_python, purpose, reason) {
  ver_script <- tempfile(fileext = ".py")
  on.exit(unlink(ver_script), add = TRUE)
  writeLines("import sys; print('%d.%d' % sys.version_info[:2])", ver_script)
  ver <- tryCatch(suppressWarnings(system2(base_python, gp_quote_system_args(ver_script),
                                           stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())
  ver <- utils::tail(ver[grepl("^[0-9]+\\.[0-9]+$", ver)], 1L)
  if (!length(ver)) {
    stop("Could not run ", base_python, " to create a PredictProR Python environment.", call. = FALSE)
  }
  env_dir <- file.path(gp_python_env_home(), paste0(purpose, "-py", ver))
  env_python <- gp_python_venv_python(env_dir)
  # written last: an environment without it was never finished (or a venv
  # call failed half-way) and is built again
  stamp <- file.path(env_dir, ".predictpror-env")
  dir.create(dirname(env_dir), recursive = TRUE, showWarnings = FALSE)
  gp_python_with_dir_lock(paste0(env_dir, ".lock"), {
    if (!file.exists(stamp) || !gp_python_runs(env_python)) {
      why <- switch(reason,
                    externally_managed = "is a system Python managed by the operating system",
                    not_writable = "cannot be written to",
                    reason)
      message(sprintf(
        "PredictProR: %s %s, so the %s packages go into a PredictProR environment at %s.",
        base_python, why, gp_python_purpose_label(purpose), env_dir
      ))
      made <- gp_python_make_venv(base_python, env_dir)
      if (!isTRUE(made$ok)) {
        stop("Could not create a Python environment at ", env_dir, " from ", base_python, ".\n",
             paste(utils::tail(made$output, 20L), collapse = "\n"),
             "\nSet PREDICTPRO_*_PYTHON to a Python you can install into.", call. = FALSE)
      }
      writeLines(base_python, stamp)
    }
  })
  gp_python_path(env_python, mustWork = TRUE)
}

# `python -m venv --system-site-packages env_dir`. Debian/Ubuntu leave
# ensurepip out of the system Python unless python3-venv is installed; venv
# then exits with an error but leaves a runnable environment without pip. In
# that case the environment is made with --without-pip and the installer adds
# pip (predictpror_deps.py: ensurepip, else pip's zipapp).
gp_python_make_venv <- function(base_python, env_dir) {
  env_python <- gp_python_venv_python(env_dir)
  attempt <- function(extra) {
    unlink(env_dir, recursive = TRUE, force = TRUE)
    out <- suppressWarnings(system2(
      base_python,
      gp_quote_system_args(c("-m", "venv", "--system-site-packages", extra, env_dir)),
      stdout = TRUE, stderr = TRUE
    ))
    status <- as.integer(attr(out, "status") %||% 0L)
    list(ok = identical(status, 0L) && gp_python_runs(env_python), output = out)
  }
  made <- attempt(character())
  if (!isTRUE(made$ok) && any(grepl("ensurepip", made$output, fixed = TRUE))) {
    retry <- attempt("--without-pip")
    retry$output <- c(made$output, retry$output)
    made <- retry
  }
  if (!isTRUE(made$ok)) {
    unlink(env_dir, recursive = TRUE, force = TRUE)
  }
  made
}

#' Install the Python packages PredictProR needs
#'
#' PredictProR installs the Python packages its models need on first use.
#' `model_execute()` does this automatically for the requested models; this
#' function does it explicitly, for example to prepare a server before running
#' jobs. It checks the interpreter PredictProR uses for `purpose`
#' (`PREDICTPRO_ML_PYTHON`/`PREDICTPRO_PYTHON`, `PREDICTPRO_DL_PYTHON`,
#' `PREDICTPRO_GP_PYTHON`, or the detected one) and installs what is missing.
#' PyTorch is installed from the PyTorch wheel index: the newest CUDA build the
#' NVIDIA driver supports, else the CPU build (`PYTORCH_INDEX_URL` forces an
#' index, `PREDICTPRO_TORCH_DEVICE = "cpu"` the CPU build). If the interpreter
#' cannot take installs (not writable, or an operating-system managed system
#' Python), a PredictProR virtualenv is created from it under
#' `tools::R_user_dir("PredictProR", "cache")` (or `PREDICTPRO_PYTHON_ENV_HOME`)
#' and used for the rest of the session. Set `PREDICTPRO_AUTO_INSTALL = "false"`
#' to turn installing off.
#'
#' @param purpose One or more of `"ml"`, `"dl"`, `"gp"`.
#' @param models Optional model names; adds the packages those models need
#'   (e.g. xgboost for `"Xgboost"`, gpytorch for `"GPNet"`).
#' @param python_bin Optional interpreter to use instead of the detected one
#'   (only with a single `purpose`).
#' @param quiet Logical; if `TRUE`, do not show pip's progress.
#'
#' @return Invisibly, a named character vector of the interpreter used per purpose.
#' @export
ensure_python_dependencies <- function(purpose = c("ml", "dl", "gp"), models = NULL,
                                       python_bin = NULL, quiet = FALSE) {
  purpose <- unique(match.arg(purpose, several.ok = TRUE))
  if (!is.null(python_bin) && length(purpose) != 1L) {
    stop("`python_bin` can only be given with a single `purpose`.", call. = FALSE)
  }
  plan <- if (length(models)) {
    gp_python_dependency_groups(gp_canonicalize_supported_model_names(models))
  } else list()
  used <- character()
  for (p in purpose) {
    groups <- unique(c(gp_python_base_groups(p), plan[[p]]))
    used[[p]] <- gp_python_ensure_dependencies(p, groups, python_bin = python_bin, quiet = quiet)
  }
  invisible(used)
}

# Make sure `python_bin` (default: the interpreter for `purpose`) has `groups`;
# returns the interpreter to use (a PredictProR virtualenv when the original
# cannot take installs).
gp_python_ensure_dependencies <- function(purpose, groups, python_bin = NULL, quiet = FALSE,
                                          managed = FALSE) {
  python_bin <- python_bin %||% gp_python_detect_for(purpose)
  if (is.null(python_bin) || !nzchar(python_bin) || !file.exists(python_bin)) {
    stop(gp_python_no_runtime_message(purpose), call. = FALSE)
  }
  python_bin <- gp_python_path(python_bin, mustWork = TRUE)
  groups <- sort(unique(groups))
  if (gp_python_deps_verified(python_bin, groups)) {
    return(invisible(python_bin))
  }
  # a check is find_spec calls plus one torch import; it never waits for pip
  status <- gp_python_deps_run(python_bin, "check", groups, timeout = 900)
  if (!identical(status$status, "ok")) {
    stop(status$message, call. = FALSE)
  }
  missing <- as.character(unlist(status$missing))
  if (length(missing)) {
    problem <- as.character(status$install_problem %||% "")[1L]
    if (!gp_python_auto_install_enabled()) {
      stop(sprintf(
        "The %s models need Python packages that are missing from %s: %s.\nPREDICTPRO_AUTO_INSTALL is false; install them, or unset it so PredictProR installs them.",
        gp_python_purpose_label(purpose), python_bin, paste(missing, collapse = ", ")
      ), call. = FALSE)
    }
    if (problem %in% c("not_writable", "externally_managed")) {
      # a PredictProR environment that cannot take installs either is not replaced again
      if (isTRUE(managed) || startsWith(python_bin, paste0(gp_python_env_home(), "/"))) {
        stop(sprintf(
          "The PredictProR Python environment %s cannot be written to (%s). Check the permissions of %s, or set PREDICTPRO_PYTHON_ENV_HOME to a writable directory.",
          python_bin, problem, gp_python_env_home()
        ), call. = FALSE)
      }
      env_python <- gp_python_managed_env(python_bin, purpose, problem)
      gp_python_set_purpose_python(purpose, env_python)
      return(gp_python_ensure_dependencies(purpose, groups, python_bin = env_python, quiet = quiet,
                                           managed = TRUE))
    }
    if (nzchar(problem) && !is.na(problem)) {
      res <- gp_python_deps_run(python_bin, "ensure", groups, stream = FALSE, extra_args = "--quiet")
      stop(res$message %||% paste("Cannot install into", python_bin, ":", problem), call. = FALSE)
    }
    message(sprintf(
      "PredictProR: installing Python packages for %s models into %s: %s.\nThis is done once per environment; PyTorch is a large download and can take several minutes.",
      gp_python_purpose_label(purpose), python_bin, paste(missing, collapse = ", ")
    ))
    res <- gp_python_deps_run(python_bin, "ensure", groups, stream = !isTRUE(quiet),
                              extra_args = if (isTRUE(quiet)) "--quiet" else character())
    if (!identical(res$status, "ok")) {
      stop(res$message, call. = FALSE)
    }
    gp_python_forget_cached_probes(python_bin)
    message(sprintf("PredictProR: Python packages for %s models are ready.", gp_python_purpose_label(purpose)))
  }
  gp_python_mark_verified(python_bin, groups)
  invisible(python_bin)
}

# Interpreter for a bridge call: verified once per session (and per worker,
# unless the main session already verified it), installing what is missing.
gp_python_ready <- function(purpose, python_bin = NULL) {
  python_bin <- python_bin %||% gp_python_detect_for(purpose)
  if (is.null(python_bin) || !nzchar(python_bin)) {
    return(python_bin)
  }
  gp_python_ensure_dependencies(purpose, gp_python_base_groups(purpose), python_bin = python_bin)
}

# model_execute(): install what the requested models need, in this session and
# before any parallel worker starts.
gp_python_prepare_runtimes <- function(models, gp_backend = "auto", hybrid_gp = FALSE, quiet = FALSE) {
  assign("python_runtimes::run", character(), envir = PredictProR_runtime_cache)
  plan <- gp_python_dependency_groups(models, gp_backend = gp_backend, hybrid_gp = hybrid_gp)
  used <- character()
  for (purpose in names(plan)) {
    used[[purpose]] <- gp_python_ensure_dependencies(purpose, plan[[purpose]], quiet = quiet)
  }
  assign("python_runtimes::run", used, envir = PredictProR_runtime_cache)
  invisible(used)
}

# Interpreters (by purpose) the current model_execute() run prepared; written
# to Run_metadata as python_ml / python_dl / python_gp.
gp_python_runtimes_used <- function() {
  get0("python_runtimes::run", envir = PredictProR_runtime_cache, inherits = FALSE) %||% character()
}
