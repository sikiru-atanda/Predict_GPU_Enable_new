gp_mori_available <- function() {
  requireNamespace("mori", quietly = TRUE) &&
    utils::packageVersion("mori") >= package_version("0.2.0") &&
    exists("share", envir = asNamespace("mori"), mode = "function", inherits = FALSE)
}

gp_mori_enabled <- function() {
  flag <- tolower(Sys.getenv("GP_USE_MORI", "auto"))
  !(flag %in% c("false", "0", "no", "off"))
}

gp_mori_score_weights <- function() {
  list(
    size = suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_SIZE", "1.0"))),
    fanout = suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_FANOUT", "0.75"))),
    reuse = suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_REUSE", "0.5"))),
    copy_volume = suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_COPYVOL", "0.35"))),
    overhead = suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_OVERHEAD", "1.0")))
  )
}

gp_mori_backend_eligible <- function(backend = NULL, plan = NULL) {
  backend <- as.character(backend %||% "")
  if (identical(backend, "mirai")) {
    return(TRUE)
  }
  if (identical(backend, "future")) {
    if (is.null(plan)) {
      return(TRUE)
    }
    if (identical(plan, future::multisession)) {
      return(TRUE)
    }
    return(FALSE)
  }
  identical(backend, "base_parallel") || identical(backend, "foreach")
}

gp_mori_object_eligible <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  if (is.atomic(x) || is.matrix(x)) {
    return(TRUE)
  }
  if (is.data.frame(x)) {
    cols_ok <- vapply(x, function(col) is.atomic(col) || is.matrix(col), logical(1))
    return(all(cols_ok))
  }
  if (is.list(x)) {
    return(all(vapply(x, function(el) {
      is.null(el) || is.atomic(el) || is.matrix(el) || is.data.frame(el)
    }, logical(1))))
  }
  FALSE
}

gp_mori_min_mb <- function() {
  suppressWarnings(as.numeric(Sys.getenv("GP_MORI_MIN_MB", "128")))
}

gp_mori_score_threshold <- function() {
  suppressWarnings(as.numeric(Sys.getenv("GP_MORI_SCORE_THRESHOLD", "1.0")))
}

gp_mori_object_size_mb <- function(x) {
  as.numeric(utils::object.size(x)) / (1024^2)
}

gp_mori_decide_one <- function(x,
                               backend = NULL,
                               plan = NULL,
                               workers = 1L,
                               n_tasks = 1L,
                               name = NULL) {
  size_mb <- gp_mori_object_size_mb(x)
  min_mb <- gp_mori_min_mb()
  threshold <- gp_mori_score_threshold()
  workers <- max(1L, as.integer(workers %||% 1L))
  n_tasks <- max(1L, as.integer(n_tasks %||% 1L))
  fanout <- max(1L, min(workers, n_tasks))
  reuse_ratio <- n_tasks / fanout
  copy_volume_mb <- size_mb * fanout
  weights <- gp_mori_score_weights()

  gates <- list(
    enabled = gp_mori_enabled(),
    available = gp_mori_available(),
    backend = gp_mori_backend_eligible(backend = backend, plan = plan),
    object = gp_mori_object_eligible(x),
    size = is.finite(size_mb) && is.finite(min_mb) && size_mb >= min_mb
  )

  if (!all(unlist(gates))) {
    reason <- names(gates)[!unlist(gates)][[1]]
    return(list(
      use_mori = FALSE,
      score = -Inf,
      threshold = threshold,
      size_mb = size_mb,
      min_mb = min_mb,
      workers = workers,
      n_tasks = n_tasks,
      fanout = fanout,
      reuse_ratio = reuse_ratio,
      copy_volume_mb = copy_volume_mb,
      backend = as.character(backend %||% ""),
      name = name %||% "<unnamed>",
      reason = paste0("gate_", reason)
    ))
  }

  size_signal <- log1p(max(size_mb, 0) / max(min_mb, 1e-6))
  fanout_signal <- log1p(max(fanout - 1L, 0L))
  reuse_signal <- log1p(max(reuse_ratio - 1, 0))
  copy_signal <- log1p(max(copy_volume_mb, 0) / max(min_mb, 1e-6))
  overhead_penalty <- 1 / log1p(max(size_mb, 1))

  score <-
    weights$size * size_signal +
    weights$fanout * fanout_signal +
    weights$reuse * reuse_signal +
    weights$copy_volume * copy_signal -
    weights$overhead * overhead_penalty

  list(
    use_mori = is.finite(score) && score >= threshold,
    score = score,
    threshold = threshold,
    size_mb = size_mb,
    min_mb = min_mb,
    workers = workers,
    n_tasks = n_tasks,
    fanout = fanout,
    reuse_ratio = reuse_ratio,
    copy_volume_mb = copy_volume_mb,
    backend = as.character(backend %||% ""),
    name = name %||% "<unnamed>",
    reason = if (is.finite(score) && score >= threshold) "score_share" else "score_skip"
  )
}

gp_mori_should_share <- function(x,
                                 backend = NULL,
                                 plan = NULL,
                                 workers = 1L,
                                 n_tasks = 1L,
                                 name = NULL) {
  isTRUE(gp_mori_decide_one(
    x = x,
    backend = backend,
    plan = plan,
    workers = workers,
    n_tasks = n_tasks,
    name = name
  )$use_mori)
}

gp_mori_share_one <- function(x,
                              backend = NULL,
                              plan = NULL,
                              workers = 1L,
                              n_tasks = 1L,
                              name = NULL) {
  decision <- gp_mori_decide_one(
    x = x,
    backend = backend,
    plan = plan,
    workers = workers,
    n_tasks = n_tasks,
    name = name
  )

  if (!isTRUE(decision$use_mori)) {
    return(list(
      object = x,
      shared = FALSE,
      size_mb = decision$size_mb,
      name = decision$name,
      reason = decision$reason,
      score = decision$score,
      threshold = decision$threshold,
      fanout = decision$fanout,
      reuse_ratio = decision$reuse_ratio,
      copy_volume_mb = decision$copy_volume_mb
    ))
  }

  out <- tryCatch(
    mori::share(x),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(list(
      object = x,
      shared = FALSE,
      size_mb = decision$size_mb,
      name = decision$name,
      reason = paste0("share_failed: ", conditionMessage(out)),
      score = decision$score,
      threshold = decision$threshold,
      fanout = decision$fanout,
      reuse_ratio = decision$reuse_ratio,
      copy_volume_mb = decision$copy_volume_mb
    ))
  }

  list(
    object = out,
    shared = TRUE,
    size_mb = decision$size_mb,
    name = decision$name,
    reason = "shared",
    score = decision$score,
    threshold = decision$threshold,
    fanout = decision$fanout,
    reuse_ratio = decision$reuse_ratio,
    copy_volume_mb = decision$copy_volume_mb
  )
}

gp_mori_share_globals <- function(globals,
                                  backend = NULL,
                                  plan = NULL,
                                  workers = 1L,
                                  n_tasks = 1L) {
  if (!is.list(globals) || !length(globals)) {
    return(list(globals = globals, summary = data.frame()))
  }

  entries <- lapply(names(globals), function(nm) {
    gp_mori_share_one(
      globals[[nm]],
      backend = backend,
      plan = plan,
      workers = workers,
      n_tasks = n_tasks,
      name = nm
    )
  })
  names(entries) <- names(globals)

  out_globals <- globals
  for (nm in names(entries)) {
    out_globals[[nm]] <- entries[[nm]]$object
  }

  summary_df <- do.call(rbind, lapply(entries, function(ent) {
    data.frame(
      name = ent$name,
      shared = ent$shared,
      size_mb = ent$size_mb,
      reason = ent$reason,
      score = ent$score,
      threshold = ent$threshold,
      fanout = ent$fanout,
      reuse_ratio = ent$reuse_ratio,
      copy_volume_mb = ent$copy_volume_mb,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL

  list(globals = out_globals, summary = summary_df)
}

gp_mori_log_summary <- function(summary_df) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) {
    return(invisible(NULL))
  }
  n_shared <- sum(summary_df$shared, na.rm = TRUE)
  if (n_shared > 0L) {
    logger::log_info(
      "mori shared {n_shared}/{nrow(summary_df)} global objects: {paste(summary_df$name[summary_df$shared], collapse = ', ')}"
    )
  } else {
    logger::log_debug("mori active but no globals met sharing policy")
  }
  invisible(summary_df)
}
