PredictProR_py_payload_cache_env <- new.env(parent = emptyenv())

gp_bridge_csv_write <- function(x, path, quote = FALSE) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, quote = quote, na = "NA")
  } else {
    utils::write.csv(x, path, row.names = FALSE, quote = quote)
  }
  invisible(path)
}

gp_bridge_payload_cache_enabled <- function(prefix) {
  value <- tolower(Sys.getenv(
    paste0("PREDICTPRO_", toupper(prefix), "_PAYLOAD_CACHE"),
    unset = "true"
  ))
  !value %in% c("0", "false", "no", "off")
}

gp_bridge_payload_cache_int <- function(prefix, suffix, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(
    paste0("PREDICTPRO_", toupper(prefix), "_PAYLOAD_CACHE_", suffix),
    unset = as.character(default)
  )))
  if (!is.finite(value) || value < 0L) {
    return(as.integer(default))
  }
  as.integer(value)
}

gp_bridge_payload_cache_cells <- function(...) {
  objects <- list(...)
  sum(vapply(objects, function(x) {
    if (is.null(x)) {
      return(0)
    }
    if (is.matrix(x) || is.data.frame(x)) {
      return(length(x))
    }
    length(x)
  }, numeric(1L)), na.rm = TRUE)
}

gp_bridge_payload_cache_state <- function(prefix) {
  state_name <- paste0("payload_cache_", tolower(prefix))
  if (!exists(state_name, envir = PredictProR_py_payload_cache_env, inherits = FALSE)) {
    assign(state_name, new.env(parent = emptyenv()), envir = PredictProR_py_payload_cache_env)
  }
  get(state_name, envir = PredictProR_py_payload_cache_env, inherits = FALSE)
}

gp_bridge_payload_cache_root <- function(prefix, state) {
  if (!exists(".root", envir = state, inherits = FALSE)) {
    root <- file.path(
      tempdir(),
      paste0("predictpror_", tolower(prefix), "_payload_cache_", Sys.getpid())
    )
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    assign(".root", normalizePath(root, winslash = "/", mustWork = TRUE), envir = state)
  }
  get(".root", envir = state, inherits = FALSE)
}

gp_bridge_payload_cache_hash <- function(prefix, purpose, key_parts) {
  key_file <- tempfile(
    pattern = paste0("predictpror_", tolower(prefix), "_payload_key_"),
    fileext = ".rds"
  )
  on.exit(unlink(key_file, force = TRUE), add = TRUE)
  saveRDS(
    list(
      version = 1L,
      prefix = toupper(prefix),
      purpose = purpose,
      payload = key_parts
    ),
    key_file,
    version = 3L,
    compress = FALSE
  )
  paste0("entry_", unname(tools::md5sum(key_file)))
}

gp_bridge_payload_cache_prune <- function(prefix, state, max_entries) {
  keys <- ls(envir = state, all.names = FALSE)
  keys <- keys[startsWith(keys, "entry_")]
  if (length(keys) <= max_entries) {
    return(invisible(NULL))
  }

  last_used <- vapply(keys, function(key) {
    entry <- get(key, envir = state, inherits = FALSE)
    as.numeric(entry$last_used %||% entry$created %||% Sys.time())
  }, numeric(1L))
  drop_keys <- keys[order(last_used)][seq_len(length(keys) - max_entries)]
  root <- normalizePath(gp_bridge_payload_cache_root(prefix, state), winslash = "/", mustWork = TRUE)
  for (key in drop_keys) {
    entry <- get(key, envir = state, inherits = FALSE)
    path <- normalizePath(entry$path %||% "", winslash = "/", mustWork = FALSE)
    if (nzchar(path) && startsWith(path, root) && dir.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
    }
    rm(list = key, envir = state)
  }
  invisible(NULL)
}

gp_bridge_payload_cache_prepare <- function(prefix,
                                            purpose,
                                            key_parts,
                                            cells,
                                            build_fun,
                                            min_cells = 20000L,
                                            max_entries = 12L) {
  if (!gp_bridge_payload_cache_enabled(prefix)) {
    return(NULL)
  }
  min_cells <- gp_bridge_payload_cache_int(prefix, "MIN_CELLS", min_cells)
  max_entries <- gp_bridge_payload_cache_int(prefix, "MAX_ENTRIES", max_entries)
  if (max_entries <= 0L || cells < min_cells) {
    return(NULL)
  }

  key <- tryCatch(
    gp_bridge_payload_cache_hash(prefix, purpose, key_parts),
    error = function(e) NULL
  )
  if (is.null(key) || !nzchar(key)) {
    return(NULL)
  }

  state <- gp_bridge_payload_cache_state(prefix)
  if (exists(key, envir = state, inherits = FALSE)) {
    entry <- get(key, envir = state, inherits = FALSE)
    if (dir.exists(entry$path)) {
      entry$last_used <- Sys.time()
      assign(key, entry, envir = state)
      return(list(input_dir = entry$path, cache_key = key, cache_hit = TRUE))
    }
    rm(list = key, envir = state)
  }

  root <- gp_bridge_payload_cache_root(prefix, state)
  input_dir <- file.path(root, key)
  tmp_dir <- tempfile(pattern = paste0(key, "_tmp_"), tmpdir = root)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- FALSE
  on.exit({
    if (!ok && dir.exists(tmp_dir)) {
      unlink(tmp_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  build_fun(tmp_dir)
  if (!dir.exists(input_dir)) {
    renamed <- file.rename(tmp_dir, input_dir)
    if (!isTRUE(renamed)) {
      dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(tmp_dir, full.names = TRUE), input_dir, overwrite = TRUE)
      unlink(tmp_dir, recursive = TRUE, force = TRUE)
    }
  } else {
    unlink(tmp_dir, recursive = TRUE, force = TRUE)
  }
  ok <- TRUE

  entry <- list(
    path = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    created = Sys.time(),
    last_used = Sys.time()
  )
  assign(key, entry, envir = state)
  gp_bridge_payload_cache_prune(prefix, state, max_entries = max_entries)
  list(input_dir = entry$path, cache_key = key, cache_hit = FALSE)
}
