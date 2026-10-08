
#' Title
#'
#' @param i_vec
#' @param j_vec
#' @param ploidy One positive integer or `"auto"`. Above diploidy, `R2` is
#'   squared Pearson correlation between allele dosages and `Dprime` is not
#'   reported because unphased dosage does not identify haplotype D-prime.
#'
#' @return
#' @export
#'
#' @examples
ld_pair_vec <- function(i_vec, j_vec, ploidy = "auto") {
  keep <- !(is.na(i_vec) | is.na(j_vec))
  g1 <- i_vec[keep]; g2 <- j_vec[keep]
  n  <- length(g1)
  if (n == 0) return(c(NA_real_, NA_real_))
  pair <- cbind(g1, g2)
  resolved_ploidy <- gp_resolve_matrix_ploidy(pair, ploidy)
  gp_validate_alt_dosage(pair, resolved_ploidy, hard_calls = FALSE)
  if (resolved_ploidy != 2L) {
    dosage_r <- suppressWarnings(stats::cor(g1, g2))
    return(c(if (is.finite(dosage_r)) dosage_r^2 else NA_real_, NA_real_))
  }
  pA <- mean(g1) / resolved_ploidy
  pB <- mean(g2) / resolved_ploidy
  pAB <- mean((g1 == 2 & g2 == 2) +
                0.5 * (g1 == 2 & g2 == 1) +
                0.5 * (g1 == 1 & g2 == 2) +
                0.25 * (g1 == 1 & g2 == 1))
  D  <- pAB - pA * pB
  denom <- if (D >= 0)
    min(pA * (1 - pB), (1 - pA) * pB)
  else
    max(-pA * pB, -(1 - pA) * (1 - pB))
  Dprime <- if (denom == 0) NA_real_ else D / denom
  R2     <- if (pA * (1 - pA) * pB * (1 - pB) == 0)
    NA_real_
  else
    D * D / (pA * (1 - pA) * pB * (1 - pB))
  c(R2, Dprime)
}

#' Title
#'
#' @param G
#' @param snp_ids
#' @param window
#' @param progress
#' @param as.data.frame
#' @param ploidy One positive integer or `"auto"`; see [ld_pair_vec()].
#'
#' @return
#' @export
#'
#' @examples
ld_window <- function(G, snp_ids = colnames(G), window = 100L,
                      progress = TRUE, as.data.frame = TRUE,
                      ploidy = "auto") {
  resolved_ploidy <- gp_resolve_matrix_ploidy(G, ploidy)
  gp_validate_alt_dosage(G, resolved_ploidy, hard_calls = FALSE)
  M <- length(snp_ids)
  if (progress) pb <- txtProgressBar(min = 0, max = M, style = 3)
  out <- vector("list", M)
  for (i in seq_len(M - 1L)) {
    if (progress) setTxtProgressBar(pb, i)
    j_end <- min(M, i + window)
    if (j_end <= i) next
    js <- seq.int(i + 1L, j_end)
    res <- vapply(
      js,
      function(j) ld_pair_vec(G[, i], G[, j], ploidy = resolved_ploidy),
      numeric(2L)
    )
    out[[i]] <- data.table::data.table(
      SNP1   = rep.int(snp_ids[i], length(js)),
      SNP2   = snp_ids[js],
      R2     = res[1L, ],
      Dprime = res[2L, ]
    )
  }
  if (progress) close(pb)
  res <- data.table::rbindlist(out, use.names = TRUE)
  if (as.data.frame) res <- as.data.frame(res)
  res
}

ld_prune_cpp_available <- function() {
  !is.null(tryCatch(
    getNativeSymbolInfo("predictpror_ld_prune_graph", PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

ld_prune_graph_cpp <- function(geno_data, window = 100L, R2_threshold = 0.9,
                               maf_thresh = 0.01, ploidy = "auto") {
  if (!is.matrix(geno_data)) {
    geno_data <- as.matrix(geno_data)
  }
  if (!is.numeric(geno_data)) {
    storage.mode(geno_data) <- "double"
  }
  resolved_ploidy <- gp_resolve_matrix_ploidy(geno_data, ploidy)
  gp_validate_alt_dosage(geno_data, resolved_ploidy, hard_calls = FALSE)
  M <- ncol(geno_data)
  snp_ids <- colnames(geno_data)
  if (is.null(snp_ids)) {
    snp_ids <- seq_len(M)
    colnames(geno_data) <- snp_ids
  }
  if (M <= 1L) {
    return(snp_ids)
  }
  keep_idx <- .Call(
    "predictpror_ld_prune_graph",
    geno_data,
    as.integer(window),
    as.numeric(R2_threshold),
    as.numeric(maf_thresh),
    as.integer(resolved_ploidy),
    PACKAGE = "PredictProR"
  )
  snp_ids[as.integer(keep_idx)]
}

#' R fallback LD-pruning of marker data (graph / window approach)
#'
#' Pure-R fallback for the C++ LD-pruning routine: drops SNPs whose pairwise
#' `r^2` with a previously kept SNP exceeds `R2_threshold` within a sliding
#' window of `window` markers, after removing SNPs below `maf_thresh`.
#'
#' @param geno_data Numeric matrix of genomic markers (rows = individuals,
#'   columns = SNPs, coded as ALT dosage from 0 through `ploidy`).
#' @param window Integer; window size (number of markers) used for the
#'   pairwise `r^2` comparison.
#' @param R2_threshold Numeric in `[0, 1]`; SNPs with pairwise `r^2` above this
#'   are pruned.
#' @param maf_thresh Numeric in `[0, 0.5]`; SNPs below this minor-allele
#'   frequency are dropped before pruning.
#' @param ploidy One positive integer or `"auto"`. Polyploid LD pruning uses
#'   squared dosage correlation and ploidy-aware MAF.
#'
#' @return Integer vector of kept SNP column indices (or a pruned marker
#'   matrix, depending on the call site).
#' @export
#'
#' @examples
ld_prune_graph_r <- function(geno_data, window = 100L, R2_threshold = 0.9,
                             maf_thresh = 0.01, ploidy = "auto") {
  if (!is.matrix(geno_data)) geno_data <- as.matrix(geno_data)
  resolved_ploidy <- gp_resolve_matrix_ploidy(geno_data, ploidy)
  gp_validate_alt_dosage(geno_data, resolved_ploidy, hard_calls = FALSE)

  M <- ncol(geno_data)
  snp_ids <- colnames(geno_data)
  if (is.null(snp_ids)) {
    snp_ids <- seq_len(M)
    colnames(geno_data) <- snp_ids
  }
  if (M <= 1L) {
    return(snp_ids)
  }

  #message("Calculating LD...")
  message("Initialize geno optimization...")
  ld_df <- ld_window(
    as.matrix(geno_data), snp_ids = snp_ids, window = window,
    progress = TRUE, as.data.frame = TRUE, ploidy = resolved_ploidy
  )
  ld_df <- ld_df[ld_df$R2 > R2_threshold & !is.na(ld_df$R2), ]

  # compute MAF
  maf_vec <- colMeans(geno_data, na.rm = TRUE) / resolved_ploidy
  maf_vec <- pmin(maf_vec, 1 - maf_vec)

  # Remove SNPs with too low MAF
  keep_maf <- which(maf_vec >= maf_thresh)
  geno_data <- geno_data[, keep_maf, drop = FALSE]
  snp_ids <- colnames(geno_data)
  if (!length(snp_ids)) {
    return(character(0))
  }
  maf_vec <- maf_vec[snp_ids]

  if (!nrow(ld_df) || !all(c("SNP1", "SNP2") %in% names(ld_df))) {
    message(sprintf("Trimmed to %d SNPs from %d", length(snp_ids), M))
    return(snp_ids)
  }

  #message("Building LD graph...")
  g <- igraph::graph_from_data_frame(ld_df[, c("SNP1", "SNP2"), drop = FALSE], directed = FALSE)

  # Subset graph to SNPs that passed MAF threshold
  g <- igraph::induced_subgraph(g, vids = intersect(snp_ids, igraph::V(g)$name))

  #message("Pruning LD graph...")
  comps <- igraph::components(g)

  # One SNP per component
  keep <- c()
  for (i in seq_len(comps$no)) {
    snps <- names(comps$membership[comps$membership == i])
    if (length(snps) == 1) {
      keep <- c(keep, snps)
    } else {
      maf_subset <- maf_vec[snps]
      chosen <- snps[which.max(maf_subset)]
      keep <- c(keep, chosen)
    }
  }

  # Add singleton SNPs that weren't in any LD pair
  all_ld_snps <- unique(c(ld_df$SNP1, ld_df$SNP2))
  singleton_snps <- setdiff(snp_ids, all_ld_snps)
  keep <- unique(c(keep, singleton_snps))
  keep <- snp_ids[snp_ids %in% keep]

  #message(sprintf("Pruned to %d SNPs from %d", length(keep), ncol(geno_data)))
  message(sprintf("Trimmed to %d SNPs from %d", length(keep), ncol(geno_data)))
  return(keep)
}

#' Prune markers by linkage disequilibrium graph
#'
#' Selects a reduced marker set by LD threshold using the native backend when
#' available and falling back to the R implementation otherwise.
#'
#' @param geno_data Genotype marker matrix with markers in columns.
#' @param window Integer marker window size used for LD comparisons.
#' @param R2_threshold Maximum allowed pairwise LD R-squared value.
#' @param maf_thresh Minor allele frequency threshold used while pruning.
#' @param ploidy One positive integer or `"auto"`. Polyploid pruning uses
#'   ploidy-aware MAF and squared allele-dosage correlation.
#'
#' @return Character vector of retained marker names.
#' @export
ld_prune_graph <- function(geno_data, window = 100L, R2_threshold = 0.9,
                           maf_thresh = 0.01, ploidy = "auto") {
  resolved_ploidy <- gp_resolve_matrix_ploidy(geno_data, ploidy)
  gp_validate_alt_dosage(geno_data, resolved_ploidy, hard_calls = FALSE)
  backend <- tolower(trimws(Sys.getenv("PREDICTPRO_LD_PRUNE_BACKEND", "auto")))
  use_cpp <- backend %in% c("auto", "cpp", "c++", "native", "1", "true", "yes", "on")
  force_cpp <- backend %in% c("cpp", "c++", "native")
  if (isTRUE(use_cpp) && ld_prune_cpp_available()) {
    message("Initialize geno optimization with native LD pruning...")
    out <- tryCatch(
      ld_prune_graph_cpp(
        geno_data = geno_data,
        window = window,
        R2_threshold = R2_threshold,
        maf_thresh = maf_thresh,
        ploidy = resolved_ploidy
      ),
      error = function(e) {
        if (isTRUE(force_cpp)) {
          stop(e)
        }
        warning(
          "Native LD pruning failed; falling back to R implementation: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(out)) {
      total <- if (is.matrix(geno_data)) ncol(geno_data) else ncol(as.matrix(geno_data))
      message(sprintf("Trimmed to %d SNPs from %d", length(out), total))
      return(out)
    }
  } else if (isTRUE(force_cpp)) {
    stop("Native LD pruning is not available in this PredictProR build.", call. = FALSE)
  }
  ld_prune_graph_r(
    geno_data = geno_data,
    window = window,
    R2_threshold = R2_threshold,
    maf_thresh = maf_thresh,
    ploidy = resolved_ploidy
  )
}

# rm(list = ls()); ls()
# gc()
# cat('\014')
# graphics.off()
#
# load("barley_data.Rdata")
# keep_prunned_snp <- ld_prune_graph(geno_data)
#
# geno_data <- as.matrix(geno_data)
# geno_data <- geno_data[, keep_prunned_snp, drop=FALSE]
#
# dim(geno_data)
