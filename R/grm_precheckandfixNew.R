# A helper function for consistent messages with optional color
inform_user <- function(msg, color = NULL) {
  if (!is.null(color) && requireNamespace("insight", quietly = TRUE)) {
    message(insight::print_color(msg, color))
  } else {
    message(msg)
  }
}

#' Find pairs (i, j) in a GRM where correlation > 'threshold'
#' without computing the entire correlation matrix.
#'
#' @param grm A symmetric positive semi-definite matrix, e.g., a GRM.
#' @param threshold A numeric scalar. We only return (i, j) with corr(i, j) > threshold.
#' @param diag_epsilon If your diagonal has tiny floating point values (unlikely for typical GRM),
#'   you might add a small epsilon to avoid dividing by zero.
#' @return A data.frame with columns Row, Col, and Corr for pairs > threshold.
#'
compute_correlations_above_threshold <- function(grm, threshold = 0.95, diag_epsilon = 1e-15) {
  stopifnot(is.matrix(grm), isSymmetric(grm))
  n <- nrow(grm)

  # Precompute diagonal sqrt (the "std dev" for each row in the GRM)
  # Add a small epsilon in case of nearly zero diagonals.
  diag_vals <- sqrt(diag(grm) + diag_epsilon)

  # We'll store results in a list for efficiency, then convert to data.frame at the end
  result_list <- vector("list", 1000)
  count <- 0

  # For i < j, compute correlation on the fly:
  #    corr(i, j) = grm[i, j] / (diag_vals[i] * diag_vals[j])
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      corr_ij <- grm[i, j] / (diag_vals[i] * diag_vals[j])
      if (corr_ij > threshold) {
        count <- count + 1
        # Expand our list if needed
        if (count > length(result_list)) {
          length(result_list) <- length(result_list) * 2
        }
        # Save the triple
        result_list[[count]] <- c(Row = i, Col = j, Corr = corr_ij)
      }
    }
  }

  # Drop any unused slots
  result_list <- result_list[seq_len(count)]

  # Convert to data.frame
  if (count > 0) {
    df <- do.call(rbind, result_list)
    df <- as.data.frame(df)
    # Name the rows/cols properly
    rownames(df) <- NULL
    # Optionally, attach row/colnames from the GRM
    df$RowName <- rownames(grm)[df$Row]
    df$ColName <- colnames(grm)[df$Col]
    return(df)
  } else {
    # Return empty if no pairs above threshold
    return(data.frame(Row = integer(0), Col = integer(0), Corr = numeric(0),
                      RowName = character(0), ColName = character(0)))
  }
}

# ----------------------------------------------------------------
#                 grm_kernel_diagnostic_fix
# ----------------------------------------------------------------
grm_kernel_diagnostic_fix <- function(grm_kernel_data,
                                      high_diag_cut_off = 1.2,
                                      low_diag_cut_off = 0.8,
                                      duplicate_cut_off = 0.95,
                                      optimize_diagonal = FALSE,
                                      optimize_duplicate = FALSE) {
  # A prefix to include in messages or errors
  msg <- "\n==================================================\n"

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

  # Extract diagonal, correlation, and a "sparse" version
  diag_vals <- diag(grm_kernel_data)
  #corr_grmkernel <- stats::cov2cor(grm_kernel_data)
  sparse_grm <- sparse_matrix(grm_kernel_data)
  sparse_grm <- as.data.frame(sparse_grm)
  #sparse_corr <- sparse_matrix(corr_grmkernel)

  # Combine to get correlations for all off-diagonal pairs
  # Assuming sparse_matrix returns a data.frame with columns Row, Col, and Value
  #sparse_grm <- data.frame(sparse_grm, Corr = sparse_corr[, 3])
  off_diag <- sparse_grm[sparse_grm$Row != sparse_grm$Col, ]

  # Identify potential duplicates (off-diagonal correlation above threshold)
  #potential_duplicates <- off_diag[off_diag$Corr > duplicate_cut_off, ]

  # On-the-fly approach
  potential_duplicate_df <- compute_correlations_above_threshold(grm_kernel_data,
                                                                 threshold = duplicate_cut_off)

  if (nrow(potential_duplicate_df) > 0) {
    potential_duplicates <- data.frame(
      Indiv_A = potential_duplicate_df$RowName,
      Indiv_B = potential_duplicate_df$ColName,
      Corr    = potential_duplicate_df$Corr
    )
    # sort, etc.
  } else{
    potential_duplicates <- data.frame()
  }

  # Convert row/col indices to actual individual names, if any duplicates found
  if (nrow(potential_duplicates) > 0) {
    potential_duplicates <- data.frame(
      Indiv_A = rownames(grm_kernel_data)[potential_duplicates$Row],
      Indiv_B = colnames(grm_kernel_data)[potential_duplicates$Col],
      Corr    = potential_duplicates$Corr
    )
    rownames(potential_duplicates) <- NULL
    # Sort by correlation descending
    potential_duplicates <- potential_duplicates[order(potential_duplicates$Corr, decreasing = TRUE), ]
  }

  # Identify diagonal outliers (too large or too small)
  diag_outliers <- diag_vals[diag_vals > high_diag_cut_off | diag_vals < low_diag_cut_off]
  diag_outliers <- sort(diag_outliers, decreasing = TRUE)
  diag_outliers_df <- data.frame(value = diag_outliers)

  # Make a copy for potential removal steps
  grm_kernel_data_opti <- NULL

  # 1) Remove duplicates if requested
  if (isTRUE(optimize_duplicate) && nrow(potential_duplicates) > 0) {
    duplicates_remove <- unique(c(
      potential_duplicates$Indiv_A,
      potential_duplicates$Indiv_B
    ))
    # Subset matrix to remove these individuals
    keep_idx <- !(rownames(grm_kernel_data) %in% duplicates_remove)
    grm_kernel_data_opti <- grm_kernel_data[keep_idx, keep_idx, drop = FALSE]
  }

  # 2) Remove diagonal outliers if requested
  if (isTRUE(optimize_diagonal) && length(diag_outliers) > 0) {
    outliers_remove <- names(diag_outliers)
    # If we've already removed duplicates, we remove from the new matrix
    if (!is.null(grm_kernel_data_opti)) {
      keep_idx <- !(rownames(grm_kernel_data_opti) %in% outliers_remove)
      grm_kernel_data_opti <- grm_kernel_data_opti[keep_idx, keep_idx, drop = FALSE]
    } else {
      # Otherwise, remove from the original matrix
      keep_idx <- !(rownames(grm_kernel_data) %in% outliers_remove)
      grm_kernel_data_opti <- grm_kernel_data[keep_idx, keep_idx, drop = FALSE]
    }
  }

  # Prepare the output list
  # Always return something called "clean_matrix"
  # plus info on potential duplicates or diagonal outliers
  if (!is.null(grm_kernel_data_opti)) {
    clean_matrix <- grm_kernel_data_opti
  } else {
    clean_matrix <- grm_kernel_data
  }

  # Potentially return lists of who was identified as duplicates/diagonal outliers
  res <- list(
    clean_matrix = clean_matrix
  )

  if (nrow(potential_duplicates) > 0) {
    res[["potential_duplicates"]] <- potential_duplicates
  }
  if (length(diag_outliers) > 0) {
    res[["potential_diag_outliers"]] <- diag_outliers_df
  }

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
                                ...) {
  msg <- "\n==================================================\n"

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
    all_numeric <- all(apply(grm_kernel_data, c(1, 2), is.numeric))
    if (!all_numeric) {
      stop(paste(msg, "The relationship matrix contains non-numeric values."), call. = FALSE)
    }

    # Ensure it's a matrix
    if (!is.matrix(grm_kernel_data)) {
      grm_kernel_data <- as.matrix(grm_kernel_data)
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

    # ---------------------------
    # 3) Bending if requested
    # ---------------------------
    if (isTRUE(bending) && !is.null(bend_value)) {
      # Check if matrix is positive definite
      if (!matrix_diagonistic_check(grm_kernel_data, "is_positive_definite")) {
        local_inform(paste(msg, "Matrix is not positive definite. We will bend it via nearPD."))
        grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data,
                                                    posd.tol = bend_value,
                                                    trace = FALSE)$mat)
      }
    } else {
      # If bending = FALSE, check anyway and warn/fix if not positive definite
      if (!matrix_diagonistic_check(grm_kernel_data, "is_positive_definite")) {
        local_inform(paste(msg,
                           "Matrix is not positive definite. We fix it by bending (nearPD)."
        ))
        grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data,
                                                    posd.tol = bend_value,
                                                    trace = FALSE)$mat)
      }
    }

    # ---------------------------
    # 4) Blending / Duplicates / Diagonal Outliers
    # ---------------------------
    if (!blending) {
      # 4a) Call the diagnostic fix function
      fix_res <- grm_kernel_diagnostic_fix(
        grm_kernel_data     = grm_kernel_data,
        high_diag_cut_off   = high_diag_cut_off,
        low_diag_cut_off    = low_diag_cut_off,
        duplicate_cut_off   = duplicate_cut_off,
        optimize_diagonal   = optimize_diagonal,
        optimize_duplicate  = optimize_duplicate
      )

      # Evaluate the condition number of the "clean matrix"
      working_matrix <- fix_res[["clean_matrix"]]
      rcn <- rcond(working_matrix)

      # If we still see duplicates or the RCN is too low => blend
      if ("potential_duplicates" %in% names(fix_res) || rcn < rcn_cutoff) {
        local_inform(paste(msg,
                           "Matrix still has duplicates or is ill-conditioned. We will blend with identity."
        ))
        n_row <- nrow(working_matrix)
        # Blend
        blended_matrix <- (1 - blending_value) * working_matrix +
          blending_value * diag(x = 1, nrow = n_row, ncol = n_row)

        # Re-check after blending
        fix_res2 <- grm_kernel_diagnostic_fix(
          grm_kernel_data     = blended_matrix,
          high_diag_cut_off   = high_diag_cut_off,
          low_diag_cut_off    = low_diag_cut_off,
          duplicate_cut_off   = duplicate_cut_off,
          optimize_diagonal   = optimize_diagonal,
          optimize_duplicate  = optimize_duplicate
        )

        rcn2 <- rcond(fix_res2[["clean_matrix"]])

        if ("potential_duplicates" %in% names(fix_res2) || rcn2 < rcn_cutoff) {
          local_inform(paste(
            "WARNINGS\n",
            msg,
            "The matrix still contains duplicates or is ill-conditioned.\n",
            "Consider increasing blending_value (e.g., 0.05 or higher)."
          ))
          # Attempt second blend step
          n_row <- nrow(fix_res2[["clean_matrix"]])
          grm_kernel_data <- (1 - blending_value) * fix_res2[["clean_matrix"]] +
            blending_value * diag(x = 1, nrow = n_row, ncol = n_row)
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
        # Re-check duplicates
        fix_res3 <- grm_kernel_diagnostic_fix(
          grm_kernel_data     = grm_kernel_data,
          high_diag_cut_off   = high_diag_cut_off,
          low_diag_cut_off    = low_diag_cut_off,
          duplicate_cut_off   = duplicate_cut_off,
          optimize_diagonal   = optimize_diagonal,
          optimize_duplicate  = optimize_duplicate
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
  }

  # Mark final matrix
  attr(grm_kernel_data, "cleared") <- "for_model_fit"
  return(grm_kernel_data)
}
