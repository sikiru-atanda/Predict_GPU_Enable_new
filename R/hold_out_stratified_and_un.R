#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param response
#' @param test_size
#' @param random_state
#' @param replication
#' @param sampling_method
#' @param message
#'
#' @return
#' @export
#'
#' @examples


hold_out_stratified_and_un <- function (
    pheno_data,
    gen_name,
    response,
    test_size = 0.2,
    random_state = 123,
    replication = 1,
    sampling_method = c("stratified", "unstratified"),
    message = TRUE,
    ...
) {
  msg_bar <- "\n==================================================\n"
  stopifnot(is.data.frame(pheno_data))
  if (!response %in% names(pheno_data)) stop("Column `response` not found.")
  if (!gen_name %in% names(pheno_data)) stop("Column `gen_name` not found.")
  if (!is.numeric(test_size) || test_size <= 0 || test_size >= 1)
    stop("`test_size` must be in (0,1).")
  if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state)

  y_raw <- pheno_data[[response]]
  g_ids <- pheno_data[[gen_name]]

  # Basic NA handling
  keep <- !(is.na(y_raw) | is.na(g_ids))
  if (any(!keep)) {
    if (isTRUE(message)) warning(paste0(msg_bar, "Removed ", sum(!keep), " rows with NA in response or gen_name."), call. = FALSE)
    y_raw <- y_raw[keep]
    g_ids <- g_ids[keep]
  }

  # Friendly hints (do not gate any logic on these)
  if (isTRUE(message)) {
    # Duplicates per genotype ≠ “multiple environments”, but it does signal repeated measures
    if (length(y_raw) > length(unique(g_ids))) {
      warning(paste0(msg_bar,
                     "You have repeated measurements per genotype (e.g., multi-env or replicates).\n",
                     "This hold-out splits rows, not genotypes. Ensure this matches your intent."), call. = FALSE)
    }
    if (replication > 1) {
      warning(paste0(msg_bar, "Multiple replications requested; this will resample repeatedly."), call. = FALSE)
    }
  }

  method <- match.arg(sampling_method)
  n <- length(y_raw)
  n_test_target <- round(n * test_size)

  # ----- Build strata once (for efficiency) -----
  # Rules:
  # - If factor or very few unique values: use those levels
  # - Else: quantile bins; ensure unique cut points
  grp_sample <- if (is.factor(y_raw) || length(unique(y_raw)) == 2) 2 else 5

  if (method == "stratified") {
    if (is.factor(y_raw) || length(unique(y_raw)) <= grp_sample) {
      groups <- as.integer(factor(y_raw, exclude = NULL))
    } else {
      # numeric stratification via quantiles with dedup guard
      num_bins <- grp_sample
      cut_points <- quantile(y_raw, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
      while (any(duplicated(cut_points)) && num_bins > 1) {
        num_bins <- num_bins - 1
        cut_points <- quantile(y_raw, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
      }
      if (num_bins < 2 || any(duplicated(cut_points))) {
        if (isTRUE(message)) warning(paste0(msg_bar, "Could not form unique quantile bins; falling back to unstratified."), call. = FALSE)
        method <- "unstratified"
        groups <- rep(1L, n)
      } else {
        # Use findInterval (faster than cut)
        groups <- findInterval(y_raw, vec = cut_points, all.inside = TRUE)
      }
    }
  } else {
    groups <- rep(1L, n)
  }

  # Precompute indices by group
  split_idx <- split(seq_len(n), groups)
  G <- length(split_idx)

  # Helper to compute per-group take counts that (1) sum to n_test_target, (2) never take all from a group
  compute_take_counts <- function() {
    sizes <- vapply(split_idx, length, integer(1))
    ideal <- sizes * test_size
    base_take <- pmin(pmax(floor(ideal), 0L), pmax(sizes - 1L, 0L))  # at least 1 left in train when size>0
    sum_base <- sum(base_take)
    resid <- n_test_target - sum_base
    if (resid > 0) {
      # Largest remainders method with random tie-breaking
      rema <- ideal - base_take
      ord <- order(rema, runif(G), decreasing = TRUE)
      for (j in ord) {
        if (resid <= 0) break
        if (base_take[j] < max(sizes[j] - 1L, 0L)) {
          base_take[j] <- base_take[j] + 1L
          resid <- resid - 1L
        }
      }
    } else if (resid < 0) {
      # Too many allocated: remove from smallest remainders first
      rema <- ideal - base_take
      ord <- order(rema, runif(G), decreasing = FALSE)
      for (j in ord) {
        if (resid >= 0) break
        if (base_take[j] > 0L) {
          base_take[j] <- base_take[j] - 1L
          resid <- resid + 1L
        }
      }
    }
    base_take
  }

  out <- vector("list", replication)
  names(out) <- sprintf("Rep_test_Sample%03d", seq_len(replication))

  for (r in seq_len(replication)) {
    # Make each replication reproducible but distinct
    if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state + r - 1L)

    take_counts <- compute_take_counts()

    picks <- integer(0)
    for (g in seq_len(G)) {
      idx <- split_idx[[g]]
      k <- take_counts[g]
      if (k > 0L && length(idx) > 0L) {
        # sample.int is faster; sample() would also work
        sel <- idx[sample.int(length(idx), k, replace = FALSE)]
        picks <- c(picks, sel)
      }
    }
    out[[r]] <- sort(picks)
  }

  if (method == "unstratified" && isTRUE(message)) {
    warning(paste0(msg_bar, "Stratified is usually a better sampling option for imbalanced/heterogeneous responses."), call. = FALSE)
  }

  out
}

# hold_out_stratified_and_un <- function (
#     pheno_data,
#     gen_name,
#     response,
#     test_size = 0.2,
#     random_state = 123,
#     replication = 1,
#     sampling_method = c("stratified", "unstratified"),
#     message = TRUE
# ) {
#
#   msg <- "\n==================================================\n"
#
#   if (!is.null(random_state) && is.numeric(random_state)) {
#     set.seed(random_state)
#   }
#
#   selected_sampling_method <- match.arg(sampling_method)
#
#   y <- as.double(pheno_data[[response]])
#   ID_GIDs <- as.character(unique(pheno_data[[gen_name]]))
#
#   if (length(y) > length(ID_GIDs) && isTRUE(message)) {
#     warning(paste(msg, "You have more than one environment. This cross-validation method works best with one environment."), call. = FALSE)
#   }
#
#   if (replication > 1 && isTRUE(message)) {
#     warning(paste(msg, "You are using more than one replication. This might take time."), call. = FALSE)
#   }
#
#   grp_sample <- if (length(unique(y)) == 2) 2 else 5
#
#   test_res <- vector(mode = "list", length = replication)
#   names(test_res) <- sprintf("Rep_test_Sample%03d", seq_len(replication))
#
#   for (r in seq_len(replication)) {
#     if (selected_sampling_method == "stratified") {
#       if (is.factor(y) || length(unique(y)) <= grp_sample) {
#         strata <- split(seq_len(length(y)), y)
#       } else {
#         num_bins <- grp_sample
#         cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
#
#         # Ensure cut_points are unique
#         while (any(duplicated(cut_points))) {
#           num_bins <- num_bins - 1
#           if (num_bins < 2) {
#             warning(paste(msg, "Unable to create unique bins for stratified sampling, switching to unstratified sampling"), call. = FALSE)
#             tst_sample <- sample(seq_len(length(y)), floor(length(y) * test_size), replace = FALSE)
#             test_res[[r]] <- sort(tst_sample)
#             next
#           }
#           cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
#         }
#
#         strata <- split(seq_len(length(y)), cut(y, cut_points, include.lowest = TRUE))
#       }
#       tst_sample <- unlist(lapply(strata, function(x) sample(x, ceiling(length(x) * test_size), replace = FALSE)))
#     } else {
#       tst_sample <- sample(seq_len(length(y)), floor(length(y) * test_size), replace = FALSE)
#     }
#
#     test_res[[r]] <- sort(tst_sample)
#   }
#
#   if (selected_sampling_method == "unstratified" && isTRUE(message)) {
#     warning(paste(msg, "Stratified is a better sampling option, consider using it."), call. = FALSE)
#   }
#
#   return(test_res)
# }
#
#
