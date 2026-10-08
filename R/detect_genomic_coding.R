# detect_genomic_coding <- function(object_geno = NULL) {
#   unique_values <- unique(object_geno)
#
#   if (inherits(unique_values, "character")) {
#     unique_values <- as.double(unique_values)
#   }
#
#   count_values <- table(unique_values)
#
#   coding_schemes <- list(
#     "Presence/Absence (0, 1)" = c(0, 1),
#     "Presence/Absence (0, 2)" = c(0, 2),
#     "SNP (-1, 0, 1)" = c(-1, 0, 1),
#     "SNP (0, 1, 2, -1)" = c(0, 1, 2, -1),
#     "SNP (0, 1, 2)" = c(0, 1, 2),
#     "SNP (0, 0.5, 1)" = c(0, 0.5, 1)
#   )
#
#
#   for (coding_scheme in names(coding_schemes)) {
#     if (length(count_values) == length(coding_schemes[[coding_scheme]]) &&
#         all(names(count_values) %in% coding_schemes[[coding_scheme]])) {
#       return(coding_scheme)
#     }
#   }
#
#   return("Unknown coding scheme")
# }
#' Detect genomic marker coding
#'
#' Scans marker values in chunks and returns the matching common coding scheme.
#'
#' @param object_geno Genomic marker matrix or matrix-like object.
#' @param chunk_size Number of rows to inspect per chunk.
#' @param ploidy One positive integer or `"auto"`. When supplied, the input is
#'   validated as ALT dosage in `[0, ploidy]` and reported without applying
#'   legacy diploid recoding heuristics.
#'
#' @return A character label for the detected coding scheme, or
#'   \code{"Unknown coding scheme"}.
#' @export
detect_genomic_coding <- function(object_geno = NULL, chunk_size = 1000,
                                  ploidy = "auto") {
  if (is.null(object_geno)) {
    return("Unknown coding scheme")
  }

  if (is.null(dim(object_geno))) {
    object_geno <- matrix(object_geno, ncol = 1L)
  } else if (!is.matrix(object_geno)) {
    object_geno <- as.matrix(object_geno)
  }

  requested_ploidy <- gp_validate_ploidy(ploidy)
  attributed_ploidy <- attr(object_geno, "ploidy", exact = TRUE)
  if (!identical(requested_ploidy, "auto") || !is.null(attributed_ploidy)) {
    resolved_ploidy <- gp_resolve_matrix_ploidy(object_geno, requested_ploidy)
    gp_validate_alt_dosage(object_geno, resolved_ploidy, hard_calls = FALSE)
    return(paste0("ALT dosage (0..", resolved_ploidy, "), ploidy ", resolved_ploidy))
  }

  coding_schemes <- list(
    "Presence/Absence (0, 1)" = c(0, 1),
    "Presence/Absence (0, 2)" = c(0, 2),
    "SNP (-1, 0, 1)" = c(-1, 0, 1),
    "SNP (0, 1, 2, -1)" = c(0, 1, 2, -1),
    "SNP (0, 1, 2)" = c(0, 1, 2),
    "SNP (0, 0.5, 1)" = c(0, 0.5, 1)
  )

  # Initialize a numeric vector to store unique values
  unique_values <- numeric()

  n <- nrow(object_geno)
  chunk_indices <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  for (indices in chunk_indices) {
    chunk <- object_geno[indices, ]
    chunk_unique_values <- unique(as.vector(chunk))
    chunk_unique_values <- chunk_unique_values[!is.na(chunk_unique_values)]
    unique_values <- unique(c(unique_values, chunk_unique_values))
  }

  unique_values <- sort(unique_values)

  if (length(unique_values) && is.numeric(unique_values) &&
      all(is.finite(unique_values)) && min(unique_values) >= 0 &&
      max(unique_values) > 2) {
    return("Potential polyploid ALT dosage; supply ploidy")
  }

  for (coding_scheme in names(coding_schemes)) {
    expected_values <- coding_schemes[[coding_scheme]]
    if (length(unique_values) == length(expected_values) &&
        all(unique_values == sort(expected_values))) {
      return(coding_scheme)
    }
  }

  return("Unknown coding scheme")
}
