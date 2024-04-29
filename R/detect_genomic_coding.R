#' Detect Genomic Coding Scheme
#'
#' This function detects the coding scheme used in a genomic dataset based on the unique values present.
#' It supports detection of common coding schemes for presence/absence and SNP data.
#'
#' @param object_geno Numeric or character vector containing genomic data from which to detect the coding scheme.
#'
#' @return A character string describing the detected coding scheme. Common schemes include:
#'   - "Presence/Absence (0, 1)"
#'   - "Presence/Absence (0, 2)"
#'   - "SNP (-1, 0, 1)"
#'   - "SNP (0, 1, 2, -1)"
#'   - "SNP (0, 1, 2)"
#' If the function cannot match the data to a known scheme, it returns "Unknown coding scheme".
#'
#' @examples
#' # Assuming object_geno is a numeric vector with SNP data
#' coding_scheme <- detect_genomic_coding(object_geno = c(0, 1, 2, 0, 1, 2))
#'
#' # Assuming object_geno is a character vector with presence/absence data
#' coding_scheme <- detect_genomic_coding(object_geno = c("0", "1", "0", "1"))
#'
#' @export

detect_genomic_coding <- function(object_geno = NULL) {
  unique_values <- unique(object_geno)

  if (inherits(unique_values, "character")) {
    unique_values <- as.double(unique_values)
  }

  count_values <- table(unique_values)

  coding_schemes <- list(
    "Presence/Absence (0, 1)" = c(0, 1),
    "Presence/Absence (0, 2)" = c(0, 2),
    "SNP (-1, 0, 1)" = c(-1, 0, 1),
    "SNP (0, 1, 2, -1)" = c(0, 1, 2, -1),
    "SNP (0, 1, 2)" = c(0, 1, 2),
    "SNP (0, 0.5, 1)" = c(0, 0.5, 1)
  )


  for (coding_scheme in names(coding_schemes)) {
    if (length(count_values) == length(coding_schemes[[coding_scheme]]) &&
        all(names(count_values) %in% coding_schemes[[coding_scheme]])) {
      return(coding_scheme)
    }
  }

  return("Unknown coding scheme")
}
