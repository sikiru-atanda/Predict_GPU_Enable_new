#' Title
#'
#' @param object_geno
#'
#' @return
#' @export
#'
#' @examples
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
    "SNP (0, 1, 2)" = c(0, 1, 2)
  )


  for (coding_scheme in names(coding_schemes)) {
    if (length(count_values) == length(coding_schemes[[coding_scheme]]) &&
        all(names(count_values) %in% coding_schemes[[coding_scheme]])) {
      return(coding_scheme)
    }
  }

  return("Unknown coding scheme")
}
