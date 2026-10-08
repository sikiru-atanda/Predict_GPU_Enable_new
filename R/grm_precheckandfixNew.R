# A helper function for consistent messages with optional color
inform_user <- function(msg, color = NULL) {
  if (!is.null(color) && requireNamespace("insight", quietly = TRUE)) {
    message(insight::print_color(msg, color))
  } else {
    message(msg)
  }
}

gp_kernel_qc_int <- function(value, default, min_value = 1L) {
  value <- suppressWarnings(as.integer(value %||% default))
  if (length(value) != 1L || is.na(value) || value < min_value) {
    value <- as.integer(default)
  }
  value
}

gp_kernel_qc_num <- function(value, default, min_value = 0) {
  value <- suppressWarnings(as.numeric(value %||% default))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value < min_value) {
    value <- as.numeric(default)
  }
  value
}

gp_kernel_empty_duplicate_pairs <- function() {
  data.frame(
    Row = integer(0),
    Col = integer(0),
    Corr = numeric(0),
    RowName = character(0),
    ColName = character(0),
    stringsAsFactors = FALSE
  )
}

gp_kernel_cpp_duplicate_scan_available <- function() {
  !is.null(tryCatch(
    getNativeSymbolInfo("predictpror_kernel_duplicate_pairs", PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

gp_kernel_duplicate_selected_n <- function(n, row_index = NULL) {
  if (is.null(row_index)) {
    return(as.integer(n))
  }
  row_index <- sort(unique(as.integer(row_index)))
  row_index <- row_index[!is.na(row_index) & row_index >= 1L & row_index <= n]
  as.integer(length(row_index))
}

gp_kernel_use_cpp_duplicate_scan <- function(n = NULL, selected_n = NULL) {
  value <- tolower(trimws(Sys.getenv("PREDICTPRO_KERNEL_DUP_CPP", "auto")))
  if (value %in% c("0", "false", "no", "off", "r")) {
    return(FALSE)
  }
  if (!gp_kernel_cpp_duplicate_scan_available()) {
    return(FALSE)
  }
  if (value %in% c("1", "true", "yes", "on", "cpp")) {
    return(TRUE)
  }
  selected_n <- suppressWarnings(as.integer(selected_n %||% n %||% 0L))
  min_n <- gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_CPP_MIN_N"), 3000L)
  length(selected_n) == 1L && !is.na(selected_n) && selected_n >= min_n
}

gp_kernel_block_duplicate_pairs_cpp <- function(grm,
                                                threshold = 0.95,
                                                diag_epsilon = 1e-15,
                                                block_size = 1024L,
                                                max_pairs = Inf,
                                                row_index = NULL) {
  if (!is.matrix(grm)) {
    grm <- as.matrix(grm)
  }
  storage.mode(grm) <- "double"
  raw <- .Call(
    "predictpror_kernel_duplicate_pairs",
    grm,
    as.numeric(threshold),
    as.numeric(diag_epsilon),
    as.integer(gp_kernel_qc_int(block_size, 1024L)),
    as.numeric(max_pairs),
    if (is.null(row_index)) NULL else as.integer(row_index),
    PACKAGE = "PredictProR"
  )
  if (!length(raw$Row)) {
    return(gp_kernel_empty_duplicate_pairs())
  }
  out <- data.frame(
    Row = as.integer(raw$Row),
    Col = as.integer(raw$Col),
    Corr = as.numeric(raw$Corr),
    RowName = rownames(grm)[as.integer(raw$Row)],
    ColName = colnames(grm)[as.integer(raw$Col)],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

gp_kernel_block_duplicate_pairs <- function(grm,
                                            threshold = 0.95,
                                            diag_epsilon = 1e-15,
                                            block_size = 1024L,
                                            max_pairs = Inf,
                                            row_index = NULL) {
  stopifnot(is.matrix(grm))
  n <- nrow(grm)
  if (n < 2L || threshold >= 1) {
    return(gp_kernel_empty_duplicate_pairs())
  }

  if (is.null(row_index)) {
    row_index <- seq_len(n)
  } else {
    row_index <- sort(unique(as.integer(row_index)))
    row_index <- row_index[!is.na(row_index) & row_index >= 1L & row_index <= n]
  }
  if (length(row_index) < 2L) {
    return(gp_kernel_empty_duplicate_pairs())
  }

  block_size <- gp_kernel_qc_int(block_size, 1024L)
  max_pairs <- gp_kernel_qc_num(max_pairs, Inf)
  diag_vals <- diag(grm)
  diag_scale <- sqrt(pmax(diag_vals, 0) + diag_epsilon)
  valid_scale <- is.finite(diag_scale) & diag_scale > 0
  if (!any(valid_scale[row_index])) {
    return(gp_kernel_empty_duplicate_pairs())
  }

  chunks <- vector("list", 64L)
  chunk_count <- 0L
  pair_count <- 0L
  selected_n <- length(row_index)

  for (row_start in seq(1L, selected_n - 1L, by = block_size)) {
    row_end <- min(selected_n, row_start + block_size - 1L)
    rows <- row_index[row_start:row_end]
    rows <- rows[valid_scale[rows]]
    if (!length(rows)) {
      next
    }

    for (col_start in seq(row_start, selected_n, by = block_size)) {
      col_end <- min(selected_n, col_start + block_size - 1L)
      cols <- row_index[col_start:col_end]
      cols <- cols[valid_scale[cols]]
      if (!length(cols)) {
        next
      }

      corr_block <- grm[rows, cols, drop = FALSE] /
        outer(diag_scale[rows], diag_scale[cols])
      keep <- is.finite(corr_block) & corr_block > threshold
      if (row_start == col_start) {
        keep[lower.tri(keep, diag = TRUE)] <- FALSE
      } else {
        keep[outer(rows, cols, FUN = ">=")] <- FALSE
      }

      hits <- which(keep, arr.ind = TRUE)
      if (!nrow(hits)) {
        next
      }

      remaining <- max_pairs - pair_count
      if (is.finite(remaining) && remaining <= 0) {
        break
      }
      if (is.finite(remaining) && nrow(hits) > remaining) {
        hits <- hits[seq_len(remaining), , drop = FALSE]
      }

      chunk_count <- chunk_count + 1L
      if (chunk_count > length(chunks)) {
        length(chunks) <- length(chunks) * 2L
      }
      chunks[[chunk_count]] <- data.frame(
        Row = rows[hits[, 1]],
        Col = cols[hits[, 2]],
        Corr = corr_block[hits],
        RowName = rownames(grm)[rows[hits[, 1]]],
        ColName = colnames(grm)[cols[hits[, 2]]],
        stringsAsFactors = FALSE
      )
      pair_count <- pair_count + nrow(hits)
      if (is.finite(max_pairs) && pair_count >= max_pairs) {
        break
      }
    }
    if (is.finite(max_pairs) && pair_count >= max_pairs) {
      break
    }
  }

  if (!chunk_count) {
    return(gp_kernel_empty_duplicate_pairs())
  }
  out <- do.call(rbind, chunks[seq_len(chunk_count)])
  rownames(out) <- NULL
  out
}

#' Find pairs (i, j) in a GRM where correlation > 'threshold'
#' without computing the entire correlation matrix.
#'
#' @param grm A symmetric positive semi-definite matrix, e.g., a GRM.
#' @param threshold A numeric scalar. We only return (i, j) with corr(i, j) > threshold.
#' @param diag_epsilon If your diagonal has tiny floating point values (unlikely for typical GRM),
#'   you might add a small epsilon to avoid dividing by zero.
#' @return A data.frame with columns Row, Col, and Corr for pairs > threshold.
compute_correlations_above_threshold <- function(grm,
                                                 threshold = 0.95,
                                                 diag_epsilon = 1e-15,
                                                 block_size = NULL,
                                                 max_pairs = Inf,
                                                 row_index = NULL) {
  stopifnot(is.matrix(grm), isSymmetric(grm))
  selected_n <- gp_kernel_duplicate_selected_n(nrow(grm), row_index)
  duplicate_fun <- if (gp_kernel_use_cpp_duplicate_scan(nrow(grm), selected_n)) {
    gp_kernel_block_duplicate_pairs_cpp
  } else {
    gp_kernel_block_duplicate_pairs
  }
  duplicate_fun(
    grm = grm,
    threshold = threshold,
    diag_epsilon = diag_epsilon,
    block_size = block_size %||% gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_BLOCK_SIZE"), 1024L),
    max_pairs = max_pairs,
    row_index = row_index
  )
}

gp_kernel_duplicate_removals <- function(potential_duplicates, diag_vals, ids) {
  if (!nrow(potential_duplicates)) {
    return(list(remove = character(0), clusters = data.frame()))
  }

  pair_a <- match(potential_duplicates$Indiv_A, ids)
  pair_b <- match(potential_duplicates$Indiv_B, ids)
  keep_pair <- !is.na(pair_a) & !is.na(pair_b) & pair_a != pair_b
  pair_a <- pair_a[keep_pair]
  pair_b <- pair_b[keep_pair]
  if (!length(pair_a)) {
    return(list(remove = character(0), clusters = data.frame()))
  }

  parent <- seq_along(ids)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union_pair <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) {
      parent[rb] <<- ra
    }
  }
  for (k in seq_along(pair_a)) {
    union_pair(pair_a[k], pair_b[k])
  }

  roots <- vapply(seq_along(ids), find_root, integer(1))
  groups <- split(seq_along(ids), roots)
  groups <- groups[vapply(groups, length, integer(1)) > 1L]
  if (!length(groups)) {
    return(list(remove = character(0), clusters = data.frame()))
  }

  diag_target <- stats::median(diag_vals[is.finite(diag_vals)], na.rm = TRUE)
  if (!is.finite(diag_target)) {
    diag_target <- 1
  }

  remove_idx <- integer(0)
  cluster_rows <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    members <- groups[[i]]
    member_diag <- diag_vals[members]
    ord <- order(abs(member_diag - diag_target), ids[members], na.last = TRUE)
    representative <- members[ord[1]]
    remove_idx <- c(remove_idx, setdiff(members, representative))
    cluster_rows[[i]] <- data.frame(
      Cluster = i,
      Representative = ids[representative],
      Members = paste(ids[members], collapse = ","),
      Removed = paste(setdiff(ids[members], ids[representative]), collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  list(
    remove = ids[sort(unique(remove_idx))],
    clusters = do.call(rbind, cluster_rows)
  )
}

gp_kernel_identity_target <- function(grm_kernel_data) {
  diag_vals <- diag(grm_kernel_data)
  target <- stats::median(diag_vals[is.finite(diag_vals) & diag_vals > 0], na.rm = TRUE)
  if (!is.finite(target) || target <= 0) {
    target <- 1
  }
  diag(target, nrow(grm_kernel_data), ncol(grm_kernel_data))
}

gp_kernel_eigen_floor_value <- function(grm_kernel_data, bend_value = 0.01) {
  diag_vals <- diag(grm_kernel_data)
  target <- stats::median(diag_vals[is.finite(diag_vals) & diag_vals > 0], na.rm = TRUE)
  if (!is.finite(target) || target <= 0) {
    target <- 1
  }
  max(.Machine$double.eps, gp_kernel_qc_num(bend_value, 0.01) * target)
}

gp_kernel_blend_identity <- function(grm_kernel_data, blending_value = 0.02) {
  alpha <- gp_kernel_qc_num(blending_value, 0.02)
  alpha <- min(max(alpha, 0), 1)
  (1 - alpha) * grm_kernel_data + alpha * gp_kernel_identity_target(grm_kernel_data)
}

gp_kernel_pd_status <- function(grm_kernel_data,
                                check = c("exact", "sample", "skip"),
                                sample_size = 500L) {
  check <- match.arg(check)
  if (identical(check, "skip")) {
    return(list(is_positive_definite = NA, check = "skip", checked_n = 0L))
  }

  check_matrix <- grm_kernel_data
  checked_n <- nrow(check_matrix)
  if (identical(check, "sample")) {
    sample_n <- min(checked_n, gp_kernel_qc_int(sample_size, 500L))
    sample_idx <- unique(pmax(1L, pmin(checked_n, round(seq(1, checked_n, length.out = sample_n)))))
    check_matrix <- grm_kernel_data[sample_idx, sample_idx, drop = FALSE]
    checked_n <- nrow(check_matrix)
  }

  is_pd <- tryCatch({
    chol(check_matrix)
    TRUE
  }, error = function(e) FALSE)

  list(is_positive_definite = is_pd, check = check, checked_n = checked_n)
}

gp_kernel_is_positive_definite <- function(grm_kernel_data) {
  gp_kernel_pd_status(grm_kernel_data, check = "exact")$is_positive_definite
}

gp_kernel_cpp_repair_available <- function() {
  !is.null(tryCatch(
    getNativeSymbolInfo("predictpror_kernel_spd_repair", PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

gp_kernel_spd_repair_cpp <- function(grm_kernel_data,
                                     min_eigen = NULL,
                                     keep_diag = TRUE,
                                     bend_value = 0.01) {
  if (!gp_kernel_cpp_repair_available()) {
    stop("PredictProR compiled kernel repair routine is not available.", call. = FALSE)
  }
  if (!is.matrix(grm_kernel_data)) {
    grm_kernel_data <- as.matrix(grm_kernel_data)
  }
  storage.mode(grm_kernel_data) <- "double"
  dn <- dimnames(grm_kernel_data)
  if (is.null(min_eigen)) {
    min_eigen <- gp_kernel_eigen_floor_value(grm_kernel_data, bend_value)
  }
  out <- .Call(
    "predictpror_kernel_spd_repair",
    grm_kernel_data,
    as.numeric(min_eigen),
    isTRUE(keep_diag),
    PACKAGE = "PredictProR"
  )
  dimnames(out$matrix) <- dn
  out
}

gp_kernel_stabilize_pd <- function(grm_kernel_data,
                                   method = c("auto", "ridge", "nearPD_cpp", "nearPD", "none"),
                                   repair_priority = c("speed", "structure"),
                                   bend_value = 0.01,
                                   nearpd_size_limit = 2500L,
                                   cpp_size_limit = 3000L,
                                   cpp_keep_diag = TRUE,
                                   pd_check = c("exact", "sample"),
                                   pd_sample_size = 500L) {
  method <- match.arg(method)
  repair_priority <- match.arg(repair_priority)
  pd_check <- match.arg(pd_check)
  n <- nrow(grm_kernel_data)
  pd_status <- gp_kernel_pd_status(grm_kernel_data, check = pd_check, sample_size = pd_sample_size)
  if (isTRUE(pd_status$is_positive_definite)) {
    return(list(matrix = grm_kernel_data, method = "none", alpha = 0, pd_check = pd_check, repair_priority = repair_priority))
  }
  if (identical(method, "none")) {
    return(list(matrix = grm_kernel_data, method = "none_failed", alpha = NA_real_, pd_check = pd_check, repair_priority = repair_priority))
  }

  ridge_grid <- unique(c(
    gp_kernel_qc_num(bend_value, 0.01),
    0.02, 0.05, 0.10, 0.20
  ))
  ridge_grid <- ridge_grid[ridge_grid > 0 & ridge_grid <= 1]
  if (method %in% c("auto", "ridge")) {
    for (alpha in ridge_grid) {
      candidate <- gp_kernel_blend_identity(grm_kernel_data, alpha)
      candidate_status <- gp_kernel_pd_status(candidate, check = pd_check, sample_size = pd_sample_size)
      if (isTRUE(candidate_status$is_positive_definite)) {
        return(list(matrix = candidate, method = "ridge", alpha = alpha, pd_check = pd_check, repair_priority = repair_priority))
      }
    }
    if (identical(method, "ridge")) {
      return(list(matrix = candidate, method = "ridge_failed", alpha = tail(ridge_grid, 1), pd_check = pd_check, repair_priority = repair_priority))
    }
  }

  try_cpp_repair <- function() {
    if (n > cpp_size_limit || !gp_kernel_cpp_repair_available()) {
      return(NULL)
    }
    cpp_fix <- tryCatch(
      gp_kernel_spd_repair_cpp(
        grm_kernel_data,
        min_eigen = gp_kernel_eigen_floor_value(grm_kernel_data, bend_value),
        keep_diag = cpp_keep_diag,
        bend_value = bend_value
      ),
      error = function(e) {
        attr(e, "predictpror_cpp_repair_failed") <- TRUE
        e
      }
    )
    if (!inherits(cpp_fix, "error")) {
      candidate <- cpp_fix$matrix
      candidate_status <- gp_kernel_pd_status(candidate, check = pd_check, sample_size = pd_sample_size)
      if (isTRUE(candidate_status$is_positive_definite)) {
        return(list(
          matrix = candidate,
          method = "nearPD_cpp",
          alpha = NA_real_,
          pd_check = pd_check,
          repair_priority = repair_priority,
          cpp_info = cpp_fix[setdiff(names(cpp_fix), "matrix")]
        ))
      }
    }
    NULL
  }

  try_matrix_nearpd <- function() {
    if (n > nearpd_size_limit) {
      return(NULL)
    }
    candidate <- tryCatch(
      as.matrix(Matrix::nearPD(
        grm_kernel_data,
        posd.tol = gp_kernel_qc_num(bend_value, 0.01),
        trace = FALSE
      )$mat),
      error = function(e) NULL
    )
    if (is.null(candidate)) {
      return(NULL)
    }
    candidate_status <- gp_kernel_pd_status(candidate, check = pd_check, sample_size = pd_sample_size)
    if (isTRUE(candidate_status$is_positive_definite)) {
      return(list(matrix = candidate, method = "nearPD", alpha = NA_real_, pd_check = pd_check, repair_priority = repair_priority))
    }
    NULL
  }

  if (identical(method, "nearPD_cpp")) {
    cpp_res <- try_cpp_repair()
    if (!is.null(cpp_res)) {
      return(cpp_res)
    }
    return(list(
      matrix = grm_kernel_data,
      method = "nearPD_cpp_failed",
      alpha = NA_real_,
      pd_check = pd_check,
      repair_priority = repair_priority
    ))
  }

  if (identical(method, "nearPD")) {
    nearpd_res <- try_matrix_nearpd()
    if (!is.null(nearpd_res)) {
      return(nearpd_res)
    }
    return(list(
      matrix = grm_kernel_data,
      method = "nearPD_failed",
      alpha = NA_real_,
      pd_check = pd_check,
      repair_priority = repair_priority
    ))
  }

  # Ridge handles the common cheap case. If ridge fails, prefer the native
  # diagonal-preserving spectral repair within its configured size limit; keep
  # Matrix::nearPD as the fallback for unavailable, failed, or oversized C++ repair.
  repair_order <- c("nearPD_cpp", "nearPD")
  for (repair_method in repair_order) {
    res <- switch(
      repair_method,
      nearPD_cpp = try_cpp_repair(),
      nearPD = try_matrix_nearpd()
    )
    if (!is.null(res)) {
      return(res)
    }
  }

  list(matrix = grm_kernel_data, method = "unfixed_large_nearPD_skipped", alpha = NA_real_, pd_check = pd_check, repair_priority = repair_priority)
}

# ----------------------------------------------------------------
#                 grm_kernel_diagnostic_fix
# ----------------------------------------------------------------
grm_kernel_diagnostic_fix <- function(grm_kernel_data,
                                      high_diag_cut_off = 1.2,
                                      low_diag_cut_off = 0.8,
                                      duplicate_cut_off = 0.95,
                                      optimize_diagonal = FALSE,
                                      optimize_duplicate = FALSE,
                                      duplicate_scan = c("auto", "full", "sample", "none"),
                                      kernel_large_n_threshold = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_FULL_DUP_N"), 5000L),
                                      duplicate_sample_size = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_SAMPLE_SIZE"), 2000L),
                                      duplicate_block_size = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_BLOCK_SIZE"), 1024L),
                                      duplicate_max_pairs = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_MAX_PAIRS"), 10000L)) {
  # A prefix to include in messages or errors
  msg <- ""
  duplicate_scan <- match.arg(duplicate_scan)

  # Check input arguments
  if (duplicate_cut_off < 0 || duplicate_cut_off > 1) {
    stop(paste(msg, "Duplicate threshold value must be between 0 and 1."), call. = FALSE)
  }
  if (high_diag_cut_off < low_diag_cut_off) {
    stop(paste(
      msg,
      "Cut-off for large diagonal(s) must be >= the cut-off for small diagonal(s)."
    ), call. = FALSE)
  }
  if (high_diag_cut_off < 0 || low_diag_cut_off < 0) {
    stop(paste(
      msg,
      "Cut-off values for large and small diagonal(s) must be positive."
    ), call. = FALSE)
  }

  if (!is.matrix(grm_kernel_data)) {
    grm_kernel_data <- as.matrix(grm_kernel_data)
  }
  if (!isSymmetric(grm_kernel_data)) {
    stop(paste(msg, "Duplicate diagnostics require a symmetric matrix."), call. = FALSE)
  }

  n <- nrow(grm_kernel_data)
  diag_vals <- diag(grm_kernel_data)
  scan_used <- duplicate_scan
  scan_complete <- FALSE
  row_index <- NULL
  if (identical(scan_used, "auto")) {
    scan_used <- if (n <= kernel_large_n_threshold) "full" else "sample"
  }
  if (identical(scan_used, "full")) {
    scan_complete <- TRUE
  }
  if (identical(scan_used, "sample")) {
    sample_n <- min(n, duplicate_sample_size)
    row_index <- if (sample_n < n) {
      # Deterministic spread across the matrix, avoiding random-state side effects.
      unique(pmax(1L, pmin(n, round(seq(1, n, length.out = sample_n)))))
    } else {
      seq_len(n)
    }
    scan_complete <- length(row_index) == n
  }

  potential_duplicate_df <- if (identical(scan_used, "none")) {
    gp_kernel_empty_duplicate_pairs()
  } else {
    compute_correlations_above_threshold(
      grm_kernel_data,
      threshold = duplicate_cut_off,
      block_size = duplicate_block_size,
      max_pairs = duplicate_max_pairs,
      row_index = row_index
    )
  }

  if (nrow(potential_duplicate_df) > 0) {
    potential_duplicates <- data.frame(
      Indiv_A = potential_duplicate_df$RowName,
      Indiv_B = potential_duplicate_df$ColName,
      Corr    = potential_duplicate_df$Corr
    )
  } else{
    potential_duplicates <- data.frame(
      Indiv_A = character(),
      Indiv_B = character(),
      Corr = numeric(),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(potential_duplicates) > 0) {
    rownames(potential_duplicates) <- NULL
    potential_duplicates <- potential_duplicates[order(potential_duplicates$Corr, decreasing = TRUE), ]
  }

  # Identify diagonal outliers (too large or too small)
  diag_outliers <- diag_vals[diag_vals > high_diag_cut_off | diag_vals < low_diag_cut_off]
  diag_outliers <- sort(diag_outliers, decreasing = TRUE)
  diag_outliers_df <- data.frame(value = diag_outliers)

  remove_ids <- character(0)

  duplicate_clusters <- NULL
  if (isTRUE(optimize_duplicate) && nrow(potential_duplicates) > 0) {
    duplicate_plan <- gp_kernel_duplicate_removals(
      potential_duplicates = potential_duplicates,
      diag_vals = diag_vals,
      ids = rownames(grm_kernel_data)
    )
    remove_ids <- c(remove_ids, duplicate_plan$remove)
    duplicate_clusters <- duplicate_plan$clusters
  }

  if (isTRUE(optimize_diagonal) && length(diag_outliers) > 0) {
    remove_ids <- c(remove_ids, names(diag_outliers))
  }

  remove_ids <- unique(remove_ids[nzchar(remove_ids)])
  keep_idx <- !(rownames(grm_kernel_data) %in% remove_ids)
  clean_matrix <- grm_kernel_data[keep_idx, keep_idx, drop = FALSE]

  res <- list(
    clean_matrix = clean_matrix
  )

  if (nrow(potential_duplicates) > 0) {
    res[["potential_duplicates"]] <- potential_duplicates
  }
  if (length(diag_outliers) > 0) {
    res[["potential_diag_outliers"]] <- diag_outliers_df
  }
  if (!is.null(duplicate_clusters) && nrow(duplicate_clusters) > 0) {
    res[["duplicate_clusters"]] <- duplicate_clusters
  }
  if (length(remove_ids) > 0) {
    res[["removed_ids"]] <- remove_ids
  }
  res[["kernel_qc"]] <- list(
    n = n,
    duplicate_scan = scan_used,
    duplicate_scan_complete = scan_complete,
    duplicate_pairs_reported = nrow(potential_duplicates),
    duplicate_pair_limit_reached = nrow(potential_duplicates) >= duplicate_max_pairs,
    diagonal_outliers = length(diag_outliers),
    removed_ids = length(remove_ids)
  )

  return(res)
}

# ----------------------------------------------------------------
#                 grm_kernel_precheck
# ----------------------------------------------------------------
grm_kernel_precheck <- function(grm_kernel_data = NULL,
                                pedigree_matrix  = NULL,
                                bending          = TRUE,
                                bend_value       = 0.01,
                                blending         = FALSE,
                                blending_value   = 0.02,
                                high_diag_cut_off  = 1.2,
                                low_diag_cut_off   = 0.8,
                                duplicate_cut_off  = 0.95,
                                rcn_cutoff         = 1e-12,
                                optimize_diagonal  = FALSE,
                                optimize_duplicate = FALSE,
                                show_message       = TRUE,
                                kernel_check_level = c("auto", "light", "full"),
                                kernel_large_n_threshold = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_FULL_DUP_N"), 5000L),
                                duplicate_scan = c("auto", "full", "sample", "none"),
                                duplicate_sample_size = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_SAMPLE_SIZE"), 2000L),
                                duplicate_block_size = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_BLOCK_SIZE"), 1024L),
                                duplicate_max_pairs = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_DUP_MAX_PAIRS"), 10000L),
                                kernel_fix_method = c("auto", "ridge", "nearPD_cpp", "nearPD", "none"),
                                kernel_repair_priority = c("speed", "structure"),
                                kernel_rcn_check = c("auto", "exact", "skip"),
                                kernel_nearpd_size_limit = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_NEARPD_MAX_N"), 2500L),
                                kernel_cpp_repair_size_limit = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_CPP_REPAIR_MAX_N"), 3000L),
                                kernel_cpp_keep_diag = TRUE,
                                kernel_pd_check = c("auto", "exact", "sample", "skip"),
                                kernel_pd_sample_size = gp_kernel_qc_int(Sys.getenv("PREDICTPRO_KERNEL_PD_SAMPLE_SIZE"), 500L),
                                kernel_sanitize = c("auto", "always", "diagnostic", "none"),
                                kernel_sanitize_value = NULL,
                                ...) {
  msg <- ""
  dots <- list(...)
  if ("message" %in% names(dots)) {
    show_message <- isTRUE(dots$message)
  }
  kernel_check_level <- match.arg(kernel_check_level)
  duplicate_scan <- match.arg(duplicate_scan)
  kernel_fix_method <- match.arg(kernel_fix_method)
  kernel_repair_priority <- match.arg(kernel_repair_priority)
  kernel_rcn_check <- match.arg(kernel_rcn_check)
  kernel_pd_check <- match.arg(kernel_pd_check)
  kernel_sanitize <- match.arg(kernel_sanitize)
  kernel_sanitize_value <- gp_kernel_qc_num(
    kernel_sanitize_value %||% Sys.getenv("PREDICTPRO_KERNEL_SANITIZE_VALUE"),
    blending_value
  )

  # Helper for messages
  local_inform <- function(m) {
    if (isTRUE(show_message)) {
      inform_user(m, color = "blue")
    }
  }

  # ---------------------------
  # 1) Basic checks
  # ---------------------------
  if (!is.null(grm_kernel_data)) {
    if (!is.matrix(grm_kernel_data)) {
      grm_kernel_data <- as.matrix(grm_kernel_data)
    }
    # Check for NAs
    if (anyNA(grm_kernel_data)) {
      stop(paste(msg, "NA is not allowed in the GRM or kernel matrix."), call. = FALSE)
    }
    # Check row and column names
    if (is.null(rownames(grm_kernel_data))) {
      stop(paste(msg, "Rownames (individuals) are missing in the matrix."), call. = FALSE)
    }
    if (is.null(colnames(grm_kernel_data))) {
      stop(paste(msg, "Colnames (individuals) are missing in the matrix."), call. = FALSE)
    }
    # Check for "NA" as literal in row/col names
    if (any(colnames(grm_kernel_data) %in% c("NA", "Na", "na")) ||
        any(rownames(grm_kernel_data) %in% c("NA", "Na", "na"))) {
      stop(paste(msg, "Column or Row names cannot contain the string 'NA'."), call. = FALSE)
    }

    # Check for duplicated colnames
    dup_cols <- colnames(grm_kernel_data)[duplicated(colnames(grm_kernel_data))]
    if (length(dup_cols) > 0) {
      stop(paste(msg, "The relationship matrix contains duplicate column names (SNPs)."), call. = FALSE)
    }
    # Check for duplicated rownames
    dup_rows <- rownames(grm_kernel_data)[duplicated(rownames(grm_kernel_data))]
    if (length(dup_rows) > 0) {
      stop(paste(msg, "The relationship matrix contains duplicate row names (genotypes)."), call. = FALSE)
    }

    # Check for numeric only
    if (!is.numeric(grm_kernel_data)) {
      stop(paste(msg, "The relationship matrix contains non-numeric values."), call. = FALSE)
    }

    # Row/col names must match for a symmetric matrix
    if (!identical(colnames(grm_kernel_data), rownames(grm_kernel_data))) {
      stop(paste(msg, "Column names do not match row names."), call. = FALSE)
    }

    # ---------------------------
    # 2) Force symmetry if needed
    # ---------------------------
    if (!isSymmetric.matrix(grm_kernel_data)) {
      local_inform(paste(msg, "The relationship matrix is not symmetric. We will symmetrize it."))
      grm_kernel_data <- Matrix::forceSymmetric(grm_kernel_data)
      grm_kernel_data <- as.matrix(grm_kernel_data)
    }

    n_kernel <- nrow(grm_kernel_data)
    effective_check_level <- kernel_check_level
    if (identical(effective_check_level, "auto")) {
      effective_check_level <- if (n_kernel <= kernel_large_n_threshold) "full" else "light"
    }
    effective_duplicate_scan <- duplicate_scan
    if (identical(effective_duplicate_scan, "auto")) {
      effective_duplicate_scan <- if (identical(effective_check_level, "full")) "full" else "sample"
    }
    effective_pd_check <- kernel_pd_check
    if (identical(effective_pd_check, "auto")) {
      effective_pd_check <- if (identical(effective_check_level, "full")) "exact" else "sample"
    }

    # ---------------------------
    # 3) Stabilize positive definiteness
    # ---------------------------
    pd_fix <- list(matrix = grm_kernel_data, method = "none", alpha = 0, pd_check = effective_pd_check)
    pd_status <- gp_kernel_pd_status(
      grm_kernel_data,
      check = effective_pd_check,
      sample_size = kernel_pd_sample_size
    )
    if (isFALSE(pd_status$is_positive_definite)) {
      if (isTRUE(bending) && !is.null(bend_value)) {
        local_inform(paste(msg, "Matrix is not positive definite. We will stabilize it by adaptive ridge blending."))
        pd_fix <- gp_kernel_stabilize_pd(
          grm_kernel_data,
          method = kernel_fix_method,
          repair_priority = kernel_repair_priority,
          bend_value = bend_value,
          nearpd_size_limit = kernel_nearpd_size_limit,
          cpp_size_limit = kernel_cpp_repair_size_limit,
          cpp_keep_diag = kernel_cpp_keep_diag,
          pd_check = if (identical(effective_pd_check, "skip")) "sample" else effective_pd_check,
          pd_sample_size = kernel_pd_sample_size
        )
        grm_kernel_data <- pd_fix$matrix
      } else {
        local_inform(paste(msg, "Matrix is not positive definite. Set bending = TRUE to stabilize it."))
      }
    }

    # ---------------------------
    # 4) Blending / Duplicates / Diagonal Outliers
    # ---------------------------
    sanitizer_applied <- FALSE
    sanitizer_reason <- "none"
    sanitizer_alpha <- NA_real_
    if (!blending) {
      # 4a) Call the diagnostic fix function
      fix_res <- grm_kernel_diagnostic_fix(
        grm_kernel_data     = grm_kernel_data,
        high_diag_cut_off   = high_diag_cut_off,
        low_diag_cut_off    = low_diag_cut_off,
        duplicate_cut_off   = duplicate_cut_off,
        optimize_diagonal   = optimize_diagonal,
        optimize_duplicate  = optimize_duplicate,
        duplicate_scan = effective_duplicate_scan,
        kernel_large_n_threshold = kernel_large_n_threshold,
        duplicate_sample_size = duplicate_sample_size,
        duplicate_block_size = duplicate_block_size,
        duplicate_max_pairs = duplicate_max_pairs
      )

      # Evaluate the condition number of the "clean matrix"
      working_matrix <- fix_res[["clean_matrix"]]
      do_rcn <- identical(kernel_rcn_check, "exact") ||
        (identical(kernel_rcn_check, "auto") && nrow(working_matrix) <= kernel_large_n_threshold)
      rcn <- if (isTRUE(do_rcn)) rcond(working_matrix) else NA_real_

      # If we still see duplicates or the RCN is too low => blend
      if ("potential_duplicates" %in% names(fix_res) ||
          (!is.na(rcn) && rcn < rcn_cutoff)) {
        local_inform(paste(msg,
                           "Matrix still has duplicates or is ill-conditioned. We will blend with identity."
        ))
        blended_matrix <- gp_kernel_blend_identity(working_matrix, blending_value)
        sanitizer_applied <- TRUE
        sanitizer_reason <- "duplicates_or_low_rcond"
        sanitizer_alpha <- blending_value

        # Re-check after blending
        fix_res2 <- grm_kernel_diagnostic_fix(
          grm_kernel_data     = blended_matrix,
          high_diag_cut_off   = high_diag_cut_off,
          low_diag_cut_off    = low_diag_cut_off,
          duplicate_cut_off   = duplicate_cut_off,
          optimize_diagonal   = optimize_diagonal,
          optimize_duplicate  = optimize_duplicate,
          duplicate_scan = if (identical(effective_check_level, "full")) effective_duplicate_scan else "sample",
          kernel_large_n_threshold = kernel_large_n_threshold,
          duplicate_sample_size = duplicate_sample_size,
          duplicate_block_size = duplicate_block_size,
          duplicate_max_pairs = duplicate_max_pairs
        )

        do_rcn2 <- identical(kernel_rcn_check, "exact") ||
          (identical(kernel_rcn_check, "auto") && nrow(fix_res2[["clean_matrix"]]) <= kernel_large_n_threshold)
        rcn2 <- if (isTRUE(do_rcn2)) rcond(fix_res2[["clean_matrix"]]) else NA_real_

        if ("potential_duplicates" %in% names(fix_res2) ||
            (!is.na(rcn2) && rcn2 < rcn_cutoff)) {
          local_inform(paste(
            "WARNINGS\n",
            msg,
            "The matrix still contains duplicates or is ill-conditioned.\n",
            "Consider increasing blending_value (e.g., 0.05 or higher)."
          ))
          # Attempt second blend step
          grm_kernel_data <- gp_kernel_blend_identity(fix_res2[["clean_matrix"]], blending_value)
        } else {
          # It worked well
          if (isTRUE(show_message)) {
            local_inform(paste(
              "WARNINGS\n",
              msg,
              "Matrix had duplicates. We fixed it by blending with identity."
            ))
          }
          grm_kernel_data <- fix_res2[["clean_matrix"]]
        }
      } else {
        # No duplicates or stable enough => just use what fix_res gave us
        grm_kernel_data <- fix_res[["clean_matrix"]]
      }

    } else {
      # 4b) If blending = TRUE from the start, apply user-defined blending_stat
      if (!is.null(blending_value)) {
        grm_kernel_data <- blending_stat(
          grm_kernel_data  = grm_kernel_data,
          pedigree_matrix  = pedigree_matrix,
          blending         = TRUE,
          blending_value   = blending_value
        )
        sanitizer_applied <- TRUE
        sanitizer_reason <- "explicit_blending"
        sanitizer_alpha <- blending_value
        # Re-check duplicates
        fix_res3 <- grm_kernel_diagnostic_fix(
          grm_kernel_data     = grm_kernel_data,
          high_diag_cut_off   = high_diag_cut_off,
          low_diag_cut_off    = low_diag_cut_off,
          duplicate_cut_off   = duplicate_cut_off,
          optimize_diagonal   = optimize_diagonal,
          optimize_duplicate  = optimize_duplicate,
          duplicate_scan = effective_duplicate_scan,
          kernel_large_n_threshold = kernel_large_n_threshold,
          duplicate_sample_size = duplicate_sample_size,
          duplicate_block_size = duplicate_block_size,
          duplicate_max_pairs = duplicate_max_pairs
        )
        if ("potential_duplicates" %in% names(fix_res3)) {
          local_inform(paste(
            "WARNINGS\n",
            msg,
            "The relationship matrix contains duplicates.\n",
            "Try increasing blending_value (e.g., 0.05)."
          ))
        }
        # Update matrix after fix
        grm_kernel_data <- fix_res3[["clean_matrix"]]
      }
    }

    exact_rcn_skipped <- exists("do_rcn") && !isTRUE(do_rcn)
    light_diagnostics <- identical(effective_check_level, "light") ||
      identical(effective_pd_check, "sample") ||
      identical(effective_pd_check, "skip") ||
      isTRUE(exact_rcn_skipped)
    apply_prefit_sanitizer <- !isTRUE(sanitizer_applied) &&
      (identical(kernel_sanitize, "always") ||
         (identical(kernel_sanitize, "auto") && isTRUE(light_diagnostics)))
    if (isTRUE(apply_prefit_sanitizer)) {
      grm_kernel_data <- gp_kernel_blend_identity(grm_kernel_data, kernel_sanitize_value)
      sanitizer_applied <- TRUE
      sanitizer_reason <- if (identical(kernel_sanitize, "always")) {
        "always"
      } else {
        "large_kernel_light_diagnostics"
      }
      sanitizer_alpha <- kernel_sanitize_value
    }

    attr(grm_kernel_data, "kernel_qc") <- list(
      check_level = effective_check_level,
      duplicate_scan = effective_duplicate_scan,
      pd_check = effective_pd_check,
      pd_checked_n = pd_status$checked_n,
      pd_fix_method = pd_fix$method,
      pd_fix_alpha = pd_fix$alpha,
      pd_fix_repair_priority = pd_fix$repair_priority %||% kernel_repair_priority,
      pd_fix_cpp = pd_fix$cpp_info %||% NULL,
      rcond = if (exists("rcn")) rcn else NA_real_,
      full_duplicate_threshold_n = kernel_large_n_threshold,
      sanitizer_applied = sanitizer_applied,
      sanitizer_reason = sanitizer_reason,
      sanitizer_alpha = sanitizer_alpha
    )
  }

  # Mark final matrix
  attr(grm_kernel_data, "cleared") <- "for_model_fit"
  return(grm_kernel_data)
}
