gp_platform_family <- function(sys_name = Sys.info()[["sysname"]],
                               os_type = .Platform$OS.type) {
  sys_name <- tolower(trimws(as.character(sys_name %||% "")))
  os_type <- tolower(trimws(as.character(os_type %||% "")))

  if (identical(os_type, "windows") || identical(sys_name, "windows")) {
    return("windows")
  }
  if (identical(sys_name, "darwin")) {
    return("macos")
  }
  if (sys_name %in% c("linux", "unix") || identical(os_type, "unix")) {
    return("linux")
  }
  "other"
}

gp_python_executable_relpath <- function(os_type = .Platform$OS.type) {
  if (identical(tolower(os_type), "windows")) {
    return(file.path("Scripts", "python.exe"))
  }
  file.path("bin", "python")
}

# An interpreter path made absolute without resolving the interpreter file
# itself. A virtualenv's bin/python is a symlink to the base Python; the
# resolved path starts the base Python outside the virtualenv, so its
# packages are not seen and installs go into the base environment.
gp_python_path <- function(path, mustWork = FALSE, os_type = .Platform$OS.type) {
  if (!length(path)) {
    return(character())
  }
  if (identical(tolower(os_type), "windows")) {
    return(normalizePath(path, winslash = "/", mustWork = mustWork))
  }
  path <- path.expand(as.character(path))
  if (isTRUE(mustWork)) {
    normalizePath(path, mustWork = TRUE)   # the usual error for a missing file
  }
  out <- file.path(normalizePath(dirname(path), mustWork = FALSE), basename(path))
  out[!nzchar(path)] <- ""
  out
}

gp_system_python_commands <- function(os_type = .Platform$OS.type) {
  if (identical(tolower(os_type), "windows")) {
    c("python", "python3")
  } else {
    c("python3", "python")
  }
}

gp_system_python_candidates <- function(os_type = .Platform$OS.type) {
  candidates <- unname(Sys.which(gp_system_python_commands(os_type)))
  unique(candidates[nzchar(candidates)])
}

gp_quote_system_args <- function(args, os_type = .Platform$OS.type) {
  if (is.null(args) || !length(args)) {
    return(character())
  }
  args <- as.character(args)
  if (anyNA(args)) {
    stop("External-command arguments cannot contain NA values.", call. = FALSE)
  }
  quote_type <- if (identical(tolower(os_type), "windows")) "cmd" else "sh"
  vapply(args, shQuote, character(1L), type = quote_type, USE.NAMES = FALSE)
}

gp_detect_cores <- function(logical = TRUE, reserve = 0L, detected = NULL) {
  cores <- detected %||% parallel::detectCores(logical = logical)
  cores <- suppressWarnings(as.integer(cores)[1L])
  if (!is.finite(cores) || is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  reserve <- suppressWarnings(as.integer(reserve)[1L])
  if (!is.finite(reserve) || is.na(reserve) || reserve < 0L) {
    reserve <- 0L
  }
  max(1L, cores - reserve)
}

gp_parallel_safe_fork_preference <- function(prefer_fork = TRUE,
                                             needs_process_initializer = FALSE) {
  isTRUE(prefer_fork) && !isTRUE(needs_process_initializer) && isTRUE(gp_fork_safe_after_mirai())
}

# Forked workers crashed ("caught illegal operation") when this session had
# forked workers before, then ran mirai daemons, then forked again (Linux,
# mirai 2.7 / nanonext 1.9; measured 2026-10-10, reproducible). After
# PredictProR has started mirai daemons, or while daemons are set, workers are
# started as separate processes (PSOCK, future multisession) instead.
gp_fork_safe_after_mirai <- function() {
  if (exists("mirai_daemons_started", envir = PredictProR_runtime_cache, inherits = FALSE)) {
    return(FALSE)
  }
  !(isNamespaceLoaded("mirai") && isTRUE(tryCatch(mirai::daemons_set(), error = function(e) FALSE)))
}
