#' Title
#'
#' @param pheno_data
#' @param response
#' @param gen_name
#' @param nfolds
#' @param random_state
#' @param sampling_method
#' @param replication
#'
#' @return
#' @export
#'
#' @examples

kfolds_stratified_un <- function(
    pheno_data,
    response,
    gen_name,
    nfolds = 5,
    random_state = NULL,
    sampling_method = c("stratified", "unstratified"),
    replication = 1,
    message = TRUE,
    ...
) {
  msg_bar <- "\n==================================================\n"
  stopifnot(is.data.frame(pheno_data))
  if (!response %in% names(pheno_data)) stop("Column `response` not found.")
  if (!gen_name %in% names(pheno_data)) stop("Column `gen_name` not found.")
  if (!is.numeric(nfolds) || nfolds < 2) stop("`nfolds` must be >= 2.")
  if (!is.numeric(replication) || replication < 1) stop("`replication` must be >= 1.")

  y_raw <- pheno_data[[response]]
  g_ids <- pheno_data[[gen_name]]

  # Drop rows with NA in key fields
  keep <- !(is.na(y_raw) | is.na(g_ids))
  if (any(!keep)) {
    if (isTRUE(message)) warning(paste0(msg_bar, "Removed ", sum(!keep),
                                        " rows with NA in `response` or `gen_name`."),
                                 call. = FALSE)
    y_raw <- y_raw[keep]
    g_ids <- g_ids[keep]
  }
  n <- length(y_raw)
  if (nfolds > n) stop("`nfolds` cannot exceed the number of rows after NA filtering.")

  # Gentle hint: repeated measures per genotype (not necessarily ‘multi-environment’)
  if (isTRUE(message) && length(unique(g_ids)) < length(y_raw)) {
    warning(paste0(msg_bar,
                   "You appear to have repeated measurements per genotype (e.g., replicates or multi-env).\n",
                   "This K-fold splits rows, not genotypes. If you want genotype-level splits, say so."),
            call. = FALSE)
  }

  if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state)
  method <- match.arg(sampling_method)

  # ---------- Build strata ONCE ----------
  # Strategy:
  #  - If factor or binary/low-cardinality: use factor levels.
  #  - Else numeric: use quantile bins (default 5), dedup if ties in cutpoints.
  make_groups <- function(y) {
    # treat integer/binary as categories if <=2 unique values
    if (is.factor(y) || length(unique(y)) <= 2) {
      as.integer(factor(y, exclude = NULL))
    } else if (method == "unstratified") {
      rep(1L, length(y))
    } else {
      num_bins <- 5L
      cps <- stats::quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
      while (any(duplicated(cps)) && num_bins > 1L) {
        num_bins <- num_bins - 1L
        cps <- stats::quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
      }
      if (num_bins < 2L || any(duplicated(cps))) {
        if (isTRUE(message))
          warning(paste0(msg_bar, "Could not form unique quantile bins; falling back to unstratified."),
                  call. = FALSE)
        rep(1L, length(y))
      } else {
        # faster than cut()
        findInterval(y, vec = cps, all.inside = TRUE)
      }
    }
  }

  y <- if (is.numeric(y_raw)) as.double(y_raw) else y_raw
  groups <- make_groups(y)
  split_idx <- split(seq_len(n), groups)
  G <- length(split_idx)

  # Balanced per-stratum round-robin assignment
  assign_folds_balanced <- function(idx, k) {
    if (length(idx) == 0L) return(integer(0))
    # Shuffle within stratum, then assign 1..k round-robin to minimize imbalance
    idx_shuf <- idx[sample.int(length(idx))]
    fold_seq <- rep(seq_len(k), length.out = length(idx_shuf))
    out <- integer(length(idx_shuf))
    out[] <- fold_seq
    # Return a named integer vector: names = row indices, values = fold IDs
    stats::setNames(out, idx_shuf)
  }

  out <- vector("list", replication)
  names(out) <- sprintf("Rep%d_%dFold_CV", seq_len(replication), nfolds)

  for (r in seq_len(replication)) {
    # Make each replication reproducible but distinct
    if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state + r - 1L)

    # Build per-stratum assignments, then merge
    fold_assign_named <- unlist(
      lapply(split_idx, assign_folds_balanced, k = nfolds),
      use.names = TRUE
    )
    folds <- integer(n)
    folds[as.integer(names(fold_assign_named))] <- as.integer(fold_assign_named)

    # Safety: if any unassigned (shouldn’t happen), fill unstratified
    if (any(folds == 0L)) {
      left <- which(folds == 0L)
      folds[left] <- sample(rep(seq_len(nfolds), length.out = length(left)))
    }

    out[[r]] <- folds
  }

  out
}

# kfolds_stratified_un <- function(
#     pheno_data,
#     response,
#     gen_name,
#     nfolds = 5,
#     random_state = NULL,
#     sampling_method = c("stratified", "unstratified"),
#     replication = 1,
#     ...
# ) {
#
#   msg <- "\n==================================================\n"
#
#   if (!is.null(random_state) && is.numeric(random_state)) {
#     set.seed(random_state)
#   }
#
#   sampling_method <- match.arg(sampling_method)
#   y <- as.double(pheno_data[[response]])
#
#   if (nfolds < 2) stop(paste(msg, "Number of nfolds should be greater than one for k-fold Cross_validation"), call. = FALSE)
#   if (nfolds > length(y)) stop(paste(msg, "Y variable must be greater than number of nfolds"), call. = FALSE)
#
#   if (length(unique(pheno_data[[gen_name]])) < length(y)) {
#     warning(paste(msg, "You have more than one environment. This Cross_validation method works best with one environment."))
#   }
#
#   Rep_FoldCV <- vector("list", replication)
#
#   for (r in 1:replication) {
#     folds <- vector("integer", length(y))
#
#     if (sampling_method == "stratified") {
#       if (!is.numeric(y) || length(unique(y)) <= 2) {
#         # Handle binary and categorical variables
#         levels_y <- levels(factor(y))
#         for (lvl in levels_y) {
#           idx <- which(y == lvl)
#           fold_assignments <- sample(rep(1:nfolds, length.out = length(idx)))
#           folds[idx] <- fold_assignments
#         }
#       } else {
#         # Handle continuous variables using quantile-based stratification
#         num_bins <- min(max(floor(length(y) / nfolds), 2), 5)
#         cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
#
#         # Ensure cut_points are unique
#         while (any(duplicated(cut_points))) {
#           num_bins <- num_bins - 1
#           if (num_bins < 2) {
#             warning(paste(msg, "Unable to create unique bins for stratified sampling, switching to unstratified sampling."), call. = FALSE)
#             folds <- sample(rep(1:nfolds, length.out = length(y)))
#             Rep_FoldCV[[r]] <- folds
#             names(Rep_FoldCV)[r] <- sprintf("Rep%d_%dFold_CV", r, nfolds)
#             next
#           }
#           cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
#         }
#
#         y_binned <- cut(y, breaks = cut_points, include.lowest = TRUE)
#         levels_y <- levels(y_binned)
#         for (lvl in levels_y) {
#           idx <- which(y_binned == lvl)
#           fold_assignments <- sample(rep(1:nfolds, length.out = length(idx)))
#           folds[idx] <- fold_assignments
#         }
#       }
#     } else {
#       # Unstratified sampling
#       folds <- sample(rep(1:nfolds, length.out = length(y)))
#     }
#
#     Rep_FoldCV[[r]] <- folds
#     names(Rep_FoldCV)[r] <- sprintf("Rep%d_%dFold_CV", r, nfolds)
#   }
#
#   return(Rep_FoldCV)
# }



# test_folds = nfolds_stratified_un(pheno_data = pheno,
#                                   response = "Yield",
#                                   gen_name = "GID",
#                                   sampling_method = "stratified",
#                                   replication = 1,
#                                   nFolds = 5)
#
# nfolds <- as.double(strsplit(strsplit(names(test_folds)[1], c("_"))[[1]][2],"")[[1]][1])
