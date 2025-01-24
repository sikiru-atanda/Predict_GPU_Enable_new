
evaluate_genetic_space <- function(
    M,
    train_ids,
    test_ids,
    local_coverage_threshold = NULL,
    n_pcs = 5,
    centroid_distance_factor = 2,
    distribution_overlap_threshold = 0.7,
    mahalanobis_prob = 0.95,
    do_plot = TRUE,
    coverage_prob_for_ellipse = 0.95,
    scale_marker_data = TRUE
) {

  msg <- "\n==================================================\n"
  # ----------------------------
  # 1) Preliminary checks
  # ----------------------------
  all_ids <- unique(c(train_ids, test_ids))

  is_symmetric <- isTRUE(all.equal(M, t(M), tolerance = 1e-8))

  if (!all(all_ids %in% rownames(M))) {
    stop(paste(msg, "Some train/test IDs are not found in rownames(M)."), call. = FALSE)
  }
  if (is_symmetric) {
    if (!all(all_ids %in% colnames(M))) {
      stop(paste(msg, "Some train/test IDs not found in colnames(M) for the GRM."), call. = FALSE)
    }
  }

  if (is_symmetric) {
    M_sub <- M[all_ids, all_ids, drop = FALSE]
  } else {
    M_sub <- M[all_ids, , drop = FALSE]
  }

  train_idx <- which(all_ids %in% train_ids)
  test_idx  <- which(all_ids %in% test_ids)

  # ----------------------------
  # 2) Compute a dynamic local coverage threshold (if needed)
  # ----------------------------
  if (is.null(local_coverage_threshold)) {
    if (is_symmetric) {
      S_tt <- M_sub[train_ids, train_ids, drop = FALSE]
      S_vec <- S_tt[upper.tri(S_tt)]
      # Suppose we pick the 80th percentile similarity
      local_coverage_threshold_auto <- quantile(S_vec, probs = 0.80, na.rm = TRUE)
      local_coverage_threshold <- local_coverage_threshold_auto
      message(paste(msg, sprintf("Auto-chosen GRM coverage threshold = %.3f (80th percentile).",
                      local_coverage_threshold_auto)))

    } else {
      # Marker matrix approach
      M_train_only <- M[train_ids, , drop = FALSE]
      pca_train <- prcomp(M_train_only, scale. = scale_marker_data)

      # Calculate the cumulative variance explained by the PCs
      cum_variance <- cumsum(pca_train$sdev^2) / sum(pca_train$sdev^2)

      # Determine the number of PCs that explain at least 90% of the variation
      n_pcs <- min(which(cum_variance >= 0.90))

      k_used <- min(n_pcs, ncol(pca_train$x))
      train_coords_k <- pca_train$x[, seq_len(k_used), drop = FALSE]

      # Compute nearest neighbor distance for each training genotype (within train)
      train_dists <- sapply(seq_len(nrow(train_coords_k)), function(i) {
        this_pt <- train_coords_k[i, ]
        # Distances to all other train points
        dists <- sqrt(rowSums((train_coords_k - this_pt)^2))
        # Exclude self-distance
        dists[i] <- Inf
        min(dists)
      })

      # Suppose we pick the 90th percentile as coverage threshold
      local_coverage_threshold_auto <- quantile(train_dists, probs = 0.90, na.rm = TRUE)

      local_coverage_threshold <- local_coverage_threshold_auto
      message(paste(msg, sprintf("Auto-chosen marker coverage threshold = %.3f (90th percentile).",
                      local_coverage_threshold_auto)))
    }
  }

  # ----------------------------
  # 3) Local Coverage
  # ----------------------------
  local_coverage_detail <- NULL  # store per-test info
  local_coverage_value  <- NA    # single summary in [0,1]
  pca_marker <- NULL             # store if marker-based

  if (is_symmetric) {
    # (A) GRM-based coverage
    # For each test genotype, find maximum similarity to the training set
    coverage_scores <- sapply(test_ids, function(tid) {
      max(M_sub[tid, train_ids])
    })
    # Fraction of test genotypes >= threshold
    local_coverage_value  <- mean(coverage_scores >= local_coverage_threshold)
    local_coverage_detail <- coverage_scores

  } else {
    # (B) Marker-based coverage: nearest neighbor distance in PCA space
    # We'll do PCA on M_sub (train+test)
    pca_marker <- prcomp(M_sub, scale. = scale_marker_data)

    cum_variance <- cumsum(pca_marker$sdev^2) / sum(pca_marker$sdev^2)

    # Determine the number of PCs that explain at least 90% of the variation
    #n_pcs <- min(which(cum_variance >= 0.90))
    n_pcs <- which(cum_variance >= 0.90)[1]
    # Number of PCs to keep
    k_used <- min(n_pcs, ncol(pca_marker$x))

    # Reorder rows to match all_ids if needed
    pca_scores_full <- pca_marker$x
    if (!identical(rownames(pca_scores_full), all_ids)) {
      pca_scores_full <- pca_scores_full[match(all_ids, rownames(pca_scores_full)), ]
    }

    # Subset the top k PCs
    pc_scores <- pca_scores_full[, seq_len(k_used), drop = FALSE]

    train_coords <- pc_scores[train_idx, , drop = FALSE]
    test_coords  <- pc_scores[test_idx, , drop = FALSE]

    # For each test genotype, find min distance to training
    coverage_scores <- apply(test_coords, 1, function(test_pt) {
      dists <- sqrt(rowSums((train_coords - test_pt)^2))
      min(dists)
    })

    # Local coverage fraction: how many test genotypes have distance <= threshold
    local_coverage_value  <- mean(coverage_scores <= local_coverage_threshold)
    local_coverage_detail <- coverage_scores
  }

  # ----------------------------
  # 4) Global Distribution (Centroid, Overlap, Mahalanobis)
  # ----------------------------
  # We'll produce consistent PCA coordinates for the entire subset (train+test).
  # If M is symmetric, do SVD, else re-use the prcomp from above.

  if (is_symmetric) {
    svd_res <- svd(M_sub)
    rownames_U <- rownames(M_sub)

    # variance explained
    d2 <- svd_res$d^2
    var_explained <- d2 / sum(d2)

    cum_variance <- cumsum(svd_res$d^2) / sum(svd_res$d^2)

    # Determine the number of PCs that explain at least 90% of the variation
    n_pcs <-  which(cum_variance >= 0.90)[1]

    # keep min(n_pcs, length(svd_res$d)) components
    k_used <- min(n_pcs, length(svd_res$d))

    # build PC scores akin to PC = U * D
    pc_scores_full <- svd_res$u[, seq_len(k_used), drop = FALSE]
    for (i in seq_len(k_used)) {
      pc_scores_full[, i] <- pc_scores_full[, i] * svd_res$d[i]
    }

    rownames(pc_scores_full) <- rownames_U

    # Reorder if needed
    if (!identical(rownames(pc_scores_full), all_ids)) {
      pc_scores_full <- pc_scores_full[match(all_ids, rownames(pc_scores_full)), , drop = FALSE]
    }
    # variance explained might be truncated if k_used < length(d2)
    var_explained <- var_explained[seq_len(k_used)]

  } else {
    # We have pca_marker from above
    if (is.null(pca_marker)) {
      pca_marker <- prcomp(M_sub, scale. = scale_marker_data)
    }
    d2 <- pca_marker$sdev^2
    var_explained <- d2 / sum(d2)

    # Calculate the cumulative variance explained by the PCs
    cum_variance <- cumsum(pca_marker$sdev^2) / sum(pca_marker$sdev^2)

    # Determine the number of PCs that explain at least 90% of the variation
    n_pcs <- which(cum_variance >= 0.90)[1]

    k_used <- min(n_pcs, ncol(pca_marker$x))
    pc_scores_full <- pca_marker$x

    if (!identical(rownames(pc_scores_full), all_ids)) {
      pc_scores_full <- pc_scores_full[match(all_ids, rownames(pc_scores_full)), , drop=FALSE]
    }
    pc_scores_full <- pc_scores_full[, seq_len(k_used), drop = FALSE]
    var_explained <- var_explained[seq_len(k_used)]
  }

  # 4A) Centroid distance in first 2 PCs
  n_dims_2 <- min(2, k_used)

  train_pc2 <- pc_scores_full[train_idx, seq_len(n_dims_2), drop = FALSE]
  test_pc2  <- pc_scores_full[test_idx,  seq_len(n_dims_2), drop = FALSE]

  train_centroid <- colMeans(train_pc2)
  test_centroid  <- colMeans(test_pc2)

  centroid_distance <- sqrt(sum((train_centroid - test_centroid)^2))
  # distribution overlap (simple exponential measure in 2D)
  distribution_overlap <- exp(-0.125 * centroid_distance^2)

  # 4B) Mahalanobis coverage
  train_coords_k <- pc_scores_full[train_idx, , drop = FALSE]
  test_coords_k  <- pc_scores_full[test_idx,  , drop = FALSE]

  train_cov <- cov(train_coords_k)
  inv_cov <- tryCatch(
    solve(train_cov),
    error = function(e) {
      message(paste(msg, "Cov matrix not invertible, using pseudo-inverse."))
      MASS::ginv(train_cov)
    }
  )
  mahal_d2 <- apply(test_coords_k, 1, function(x) {
    delta <- x - colMeans(train_coords_k)
    as.numeric(t(delta) %*% inv_cov %*% delta)
  })
  chi_sq_cutoff <- stats::qchisq(mahalanobis_prob, df = k_used)
  mahalanobis_coverage <- mean(mahal_d2 <= chi_sq_cutoff)

  # ----------------------------
  # 5) Composite Score & Flags
  # ----------------------------
  # Basic flags
  local_flag <- (local_coverage_value >= 0.5)

  # measure average training scale in 2D
  train_cov_2D <- cov(train_pc2)
  avg_train_scale <- sqrt(mean(diag(train_cov_2D)))  # typical stdev in 2D

  centroid_flag <- (centroid_distance <= avg_train_scale * centroid_distance_factor)
  overlap_flag  <- (distribution_overlap >= distribution_overlap_threshold)
  mahal_flag    <- (mahalanobis_coverage >= 0.8)

  # Summaries
  flags_summary_low <- c(
    "Local coverage is low.",
    "Centroid distance is large.",
    "Distribution overlap is low.",
    "Mahalanobis coverage is low."
  )
  flags_summary_high <- c(
    "Local coverage is high.",
    "Centroid distance is low.",
    "Distribution overlap is high.",
    "Mahalanobis coverage is high."
  )

  flags_summary <- c(
    if (!local_flag) flags_summary_low[1],
    if (!centroid_flag) flags_summary_low[2],
    if (!overlap_flag) flags_summary_low[3],
    if (!mahal_flag) flags_summary_low[4]
  )
  flag_index <- which(!flags_summary_low %in% flags_summary)

  # Recommendation
  recommendation <- ""
  if (length(flags_summary) == 0) {
    recommendation <- "All coverage metrics are good. The training set is likely adequate."
  } else if (length(flags_summary) == 4) {
    recommendation <- "All metrics indicate low coverage. Consider expanding the training set."
  } else {
    recommendation <- paste(
      paste(flags_summary, collapse = " "),
      "However, the following metrics suggest optimal coverage:",
      paste(flags_summary_high[flag_index], collapse = " ")
    )
  }

  # Composite score (example approach)
  # local_coverage_value in [0,1]
  # distribution_score = average of (normed centroid distance and distribution overlap)
  norm_centroid_term <- 1 - min(1, centroid_distance / (avg_train_scale * 3))
  distribution_score <- mean(c(norm_centroid_term, distribution_overlap))
  composite_score <- mean(c(local_coverage_value, distribution_score, mahalanobis_coverage))

  # ----------------------------
  # 6) Optional Plot (PC1 vs PC2)
  # ----------------------------
  the_plot <- NULL
  if (do_plot && k_used >= 2) {
    pc_data <- data.frame(
      PC1 = pc_scores_full[, 1],
      PC2 = pc_scores_full[, 2],
      Type = factor(ifelse(all_ids %in% train_ids, "Train", "Test"),
                    levels = c("Train", "Test"))
    )
    pct1 <- var_explained[1] * 100
    pct2 <- var_explained[2] * 100

    the_plot <- ggplot2::ggplot(pc_data, ggplot2::aes(x = PC1, y = PC2, color = Type)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::stat_ellipse(level = coverage_prob_for_ellipse, linetype = 2) +
      ggplot2::labs(
        title    = "Genetic Space: Local & Global Coverage",
        subtitle = sprintf("Composite Score = %.2f", composite_score),
        x = sprintf("PC1 (%.1f%%)", pct1),
        y = sprintf("PC2 (%.1f%%)", pct2)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
  }

  # ----------------------------
  # 7) Return results
  # ----------------------------
  return(list(
    is_symmetric         = is_symmetric,
    local_coverage       = local_coverage_value,
    local_coverage_detail= local_coverage_detail,
    local_threshold_used = local_coverage_threshold,
    centroid_distance    = centroid_distance,
    distribution_overlap = distribution_overlap,
    mahalanobis_coverage = mahalanobis_coverage,
    composite_score      = composite_score,
    flags = list(
      local_flag     = local_flag,
      centroid_flag  = centroid_flag,
      overlap_flag   = overlap_flag,
      mahal_flag     = mahal_flag
    ),
    recommendation = recommendation,
    plot           = the_plot,
    var_explained  = var_explained
  ))
}


# evaluate_genetic_space <- function(
#     M,
#     train_ids,
#     test_ids,
#
#     # Local coverage thresholds
#     # If M is a GRM, we interpret this as a similarity threshold.
#     # If M is a marker matrix, it's a maximum allowed PCA distance to be "covered".
#     local_coverage_threshold = 0.5,
#
#     # Number of PCs used for distribution checks
#     n_pcs = 5,
#
#     # For global distribution checks
#     centroid_distance_factor = 2,       # how many "avg stdev" away is too far
#     distribution_overlap_threshold = 0.7,
#
#     # For Mahalanobis coverage
#     # fraction of test we want inside the training distribution
#     mahalanobis_prob = 0.95,
#     # we use chi-square quantile with df = number of PCs
#
#     # Plotting settings
#     do_plot = TRUE,
#     # ellipse coverage in the 2D plot
#     coverage_prob_for_ellipse = 0.95,
#     # if M is marker data, should we scale columns for PCA or not
#     scale_marker_data = TRUE
# ) {
#   # ----------------------------
#   # 1) Preliminary checks
#   # ----------------------------
#   all_ids <- unique(c(train_ids, test_ids))
#
#   # 1A) Check if M is symmetric (GRM) or not (marker matrix)
#   is_symmetric <- isTRUE(all.equal(M, t(M), tolerance = 1e-8))
#
#   # 1B) Ensure rownames exist
#   if (!all(all_ids %in% rownames(M))) {
#     stop("Some train/test IDs are not found in rownames(M).")
#   }
#   if (is_symmetric) {
#     # Check colnames for a GRM
#     if (!all(all_ids %in% colnames(M))) {
#       stop("Some train/test IDs not found in colnames(M) for the GRM.")
#     }
#   }
#
#   # 1C) Subset M to these genotypes
#   if (is_symmetric) {
#     M_sub <- M[all_ids, all_ids, drop = FALSE]
#   } else {
#     M_sub <- M[all_ids, , drop = FALSE]  # marker matrix: subset rows only
#   }
#
#   # ----------------------------
#   # 2) Local Coverage
#   # ----------------------------
#   local_coverage_detail <- NULL  # store per-test info
#   local_coverage_value  <- NA    # single summary in [0,1]
#   pca_marker <- NULL
#
#   if (is_symmetric) {
#     # (A) GRM-based coverage
#     # For each test genotype, find maximum similarity to the training set
#     coverage_scores <- sapply(test_ids, function(tid) {
#       max(M_sub[tid, train_ids])
#     })
#     # Fraction of test genotypes above local_coverage_threshold
#     local_coverage_value  <- mean(coverage_scores >= local_coverage_threshold)
#     local_coverage_detail <- coverage_scores
#
#   } else {
#     # (B) Marker-based coverage: nearest neighbor distance in PCA space
#     # We'll do PCA on M_sub
#     pca_marker <- prcomp(M_sub, scale. = scale_marker_data)
#     # Number of PCs to keep
#     k_used <- min(n_pcs, ncol(pca_marker$x))
#
#     # Reorder rows of pca_marker$x to match all_ids if needed
#     # Typically prcomp() keeps the same row order as input, but let's be safe:
#     pca_scores_full <- pca_marker$x
#     if (!identical(rownames(pca_scores_full), all_ids)) {
#       # reorder to match all_ids
#       pca_scores_full <- pca_scores_full[match(all_ids, rownames(pca_scores_full)), ]
#     }
#
#     # Subset the top k PCs
#     pc_scores <- pca_scores_full[, seq_len(k_used), drop = FALSE]
#
#     # Identify train/test rows
#     train_idx <- which(all_ids %in% train_ids)
#     test_idx  <- which(all_ids %in% test_ids)
#
#     train_coords <- pc_scores[train_idx, , drop = FALSE]
#     test_coords  <- pc_scores[test_idx, , drop = FALSE]
#
#     # For each test genotype, find min distance to a training genotype
#     coverage_scores <- apply(test_coords, 1, function(test_pt) {
#       dists <- sqrt(rowSums((train_coords - test_pt)^2))
#       min(dists)
#     })
#
#     # Local coverage fraction: how many test genotypes are within "threshold" distance
#     # Because threshold was originally for similarity in GRM, you might choose a different
#     # numeric if these distances are typically around 2–3 for PC space, etc.
#     local_coverage_value  <- mean(coverage_scores <= local_coverage_threshold)
#     local_coverage_detail <- coverage_scores
#   }
#
#   # ----------------------------
#   # 3) Global Distribution via PCA or SVD
#   # ----------------------------
#   # We'll produce a consistent PCA coordinate set for the entire subset (train+test).
#   # Then compute centroid distance, distribution overlap, and Mahalanobis coverage.
#
#   # If M is symmetric, do SVD
#   # If M is marker-based, reuse the prcomp from above (or do it if we haven't yet).
#
#   if (is_symmetric) {
#     svd_res <- svd(M_sub)
#     # rownames for the SVD result
#     rownames_U <- rownames(M_sub)
#
#     # variance explained
#     d2 <- svd_res$d^2
#     var_explained <- d2 / sum(d2)
#
#     # We'll keep min(n_pcs, length(svd_res$d)) components
#     k_used <- min(n_pcs, length(svd_res$d))
#
#     # Build PC scores akin to PC = U * D
#     pc_scores_full <- svd_res$u[, seq_len(k_used), drop = FALSE]
#     for (i in seq_len(k_used)) {
#       pc_scores_full[, i] <- pc_scores_full[, i] * svd_res$d[i]
#     }
#
#     rownames(pc_scores_full) <- rownames_U
#
#     # Reorder if needed
#     if (!identical(rownames(pc_scores_full), all_ids)) {
#       pc_scores_full <- pc_scores_full[match(all_ids, rownames(pc_scores_full)), , drop=FALSE]
#     }
#
#   } else {
#     # We already have pca_marker from above if M is marker-based
#     # but let's ensure we only do it once
#     if (!is.null(pca_marker)) {
#       pca_marker <- prcomp(M_sub, scale. = scale_marker_data)
#     }
#     # variance explained
#     d2 <- pca_marker$sdev^2
#     var_explained <- d2 / sum(d2)
#
#     k_used <- min(n_pcs, ncol(pca_marker$x))
#     pc_scores_full <- pca_marker$x
#
#     # Reorder rows to match all_ids if needed
#     if (!identical(rownames(pc_scores_full), all_ids)) {
#       pc_scores_full <- pc_scores_full[match(all_ids, rownames(pc_scores_full)), , drop=FALSE]
#     }
#     pc_scores_full <- pc_scores_full[, seq_len(k_used), drop = FALSE]
#   }
#
#   # 3A) Centroid distance in first 2 PCs
#   n_dims_2 <- min(2, k_used)
#   train_idx <- which(all_ids %in% train_ids)
#   test_idx  <- which(all_ids %in% test_ids)
#
#   train_pc2 <- pc_scores_full[train_idx, seq_len(n_dims_2), drop = FALSE]
#   test_pc2  <- pc_scores_full[test_idx,  seq_len(n_dims_2), drop = FALSE]
#
#   train_centroid <- colMeans(train_pc2)
#   test_centroid  <- colMeans(test_pc2)
#
#   centroid_distance <- sqrt(sum((train_centroid - test_centroid)^2))
#
#   # distribution overlap (simple exponential measure)
#   distribution_overlap <- exp(-0.125 * centroid_distance^2)
#
#   # 3B) Mahalanobis coverage
#   # We'll use the top k PCs (pc_scores_full).
#   train_coords_k <- pc_scores_full[train_idx, , drop = FALSE]
#   test_coords_k  <- pc_scores_full[test_idx,  , drop = FALSE]
#
#   # Compute train covariance, invert it
#   train_cov <- cov(train_coords_k)
#   inv_cov <- tryCatch(
#     solve(train_cov),
#     error = function(e) {
#       message("Cov matrix not invertible, using pseudo-inverse.")
#       MASS::ginv(train_cov)
#     }
#   )
#
#   # Mahalanobis distances for each test genotype
#   mahal_d2 <- apply(test_coords_k, 1, function(x) {
#     delta <- x - colMeans(train_coords_k)
#     as.numeric(t(delta) %*% inv_cov %*% delta)
#   })
#
#   # cutoff for chi-sq distribution with df = k_used
#   chi_sq_cutoff <- stats::qchisq(mahalanobis_prob, df = k_used)
#   mahalanobis_coverage <- mean(mahal_d2 <= chi_sq_cutoff)
#
#   # ----------------------------
#   # 4) Composite Score & Flags
#   # ----------------------------
#   # We'll define some naive flags:
#   # - local coverage >= 0.8 is good
#   # - centroid_distance <= (avg training scale * centroid_distance_factor)
#   # - distribution_overlap >= distribution_overlap_threshold
#   # - mahalanobis_coverage >= 0.8 (or you can define your own threshold)
#
#   local_flag <- (local_coverage_value >= 0.5)
#
#   # measure average training scale in 2D
#   train_cov_2D <- cov(train_pc2)
#   avg_train_scale <- sqrt(mean(diag(train_cov_2D)))  # typical stdev in 2D
#
#   centroid_flag <- (centroid_distance <= avg_train_scale * centroid_distance_factor)
#   overlap_flag  <- (distribution_overlap >= distribution_overlap_threshold)
#   mahal_flag    <- (mahalanobis_coverage >= 0.8)
#
#   flags_summary_low <- c(
#     "Local coverage is low.",
#     "Centroid distance is large.",
#     "Distribution overlap is low.",
#     "Mahalanobis coverage is low."
#   )
#
#   flags_summary_high <- c(
#     "Local coverage is high.",
#     "Centroid distance is low.",
#     "Distribution overlap is high.",
#     "Mahalanobis coverage is high."
#   )
#
#   # Generate combined messages based on flag conditions
#   flags_summary <- c(
#     if (!local_flag) flags_summary_low[1],
#     if (!centroid_flag) flags_summary_low[2],
#     if (!overlap_flag) flags_summary_low[3],
#     if (!mahal_flag) flags_summary_low[4]
#   )
#
#   # Identify metrics that remain optimal
#   flag_index <- which(!flags_summary_low %in% flags_summary)
#
#   # Generate recommendation based on flag states
#   recommendation <- ""
#   if (length(flags_summary) == 0) {
#     recommendation <- "All coverage metrics are good. The training set is likely adequate."
#   } else if (length(flags_summary) == 4) {
#     recommendation <- "All metrics indicate low coverage. Consider expanding the training set."
#   } else {
#     recommendation <- paste(
#       paste(flags_summary, collapse = " "),
#       "However, the following metrics suggest optimal coverage:",
#       paste(flags_summary_high[flag_index], collapse = " ")
#     )
#   }
#
#
#   # Optionally combine into a single numeric "composite_score"
#   # Example approach:
#   #   local_coverage_value in [0,1]
#   #   distribution_score = average of (normed centroid and overlap)
#   norm_centroid_term <- 1 - min(1, centroid_distance / (avg_train_scale * 3))
#   distribution_score <- mean(c(norm_centroid_term, distribution_overlap))
#   #   mahalanobis_coverage in [0,1]
#   composite_score <- mean(c(local_coverage_value, distribution_score, mahalanobis_coverage))
#
#   # ----------------------------
#   # 5) Optional Plot (PC1 vs PC2)
#   # ----------------------------
#   the_plot <- NULL
#   if (do_plot && k_used >= 2) {
#     # Prepare a data.frame for ggplot
#     pc_data <- data.frame(
#       PC1 = pc_scores_full[, 1],
#       PC2 = pc_scores_full[, 2],
#       Type = factor(ifelse(all_ids %in% train_ids, "Train", "Test"),
#                     levels = c("Train", "Test"))
#     )
#     # Variance explained
#     pct1 <- var_explained[1] * 100
#     pct2 <- var_explained[2] * 100
#
#     the_plot <- ggplot2::ggplot(pc_data, ggplot2::aes(x = PC1, y = PC2, color = Type)) +
#       ggplot2::geom_point(alpha = 0.6) +
#       ggplot2::stat_ellipse(level = coverage_prob_for_ellipse, linetype = 2) +
#       ggplot2::labs(
#         title    = "Genetic Space: Local & Global & Mahalanobis Coverage",
#         subtitle = sprintf("Composite Score = %.2f", composite_score),
#         x = sprintf("PC1 (%.1f%%)", pct1),
#         y = sprintf("PC2 (%.1f%%)", pct2)
#       ) +
#       ggplot2::theme_minimal() +
#       ggplot2::theme(legend.position = "top")
#   }
#
#   # ----------------------------
#   # 6) Return
#   # ----------------------------
#   return(list(
#     is_symmetric         = is_symmetric,
#     local_coverage       = local_coverage_value,
#     local_coverage_detail= local_coverage_detail,
#     centroid_distance    = centroid_distance,
#     distribution_overlap = distribution_overlap,
#     mahalanobis_coverage = mahalanobis_coverage,
#     composite_score      = composite_score,
#     flags = list(
#       local_flag     = local_flag,
#       centroid_flag  = centroid_flag,
#       overlap_flag   = overlap_flag,
#       mahal_flag     = mahal_flag
#     ),
#     recommendation       = recommendation,
#     plot                 = the_plot,
#     var_explained        = var_explained[seq_len(k_used)]
#   ))
# }
#



# sik <- evaluate_genetic_space(M = Geno.data,
#                               train_ids = setdiff(rownames(Geno.data), tst_GID),
#                               test_ids = tst_GID)
# sik$recommendation
#
# sik$plot
#
# sik$mahalanobis_coverage
# sik$distribution_overlap


