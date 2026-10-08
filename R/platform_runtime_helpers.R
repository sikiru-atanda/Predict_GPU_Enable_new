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
  isTRUE(prefer_fork) && !isTRUE(needs_process_initializer)
}
