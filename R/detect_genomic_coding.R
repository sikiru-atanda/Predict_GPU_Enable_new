#' Title
#' This function check genomic data for different coding
#' schemes: presence/absence or SNP data (coded as 0, 1, 2)
#'
#' @param object_geno
#'
#' @return
#' @export
#'
#' @examples
#'
detect_genomic_coding <- function(object_geno = NULL) {
  #unique_values <- unique(unlist(object_geno)) ## For dataframe
  unique_values <- unique(object_geno)

  if(inherits(unique_values, "character")){
    unique_values = as.double(unique_values)
  }

  if (all(unique_values %in% c(0, 1))) {
    return("Presence/Absence (0, 1)")

  } else if (all(unique_values %in% c(0, 2))) {
    return("Presence/Absence (0, 2)")

  } else if (all(unique_values %in% c(-1, 0, 1))) {
    return("SNP (-1, 0, 1)")

  } else {
  if (all(unique_values %in% c(0, 1, 2))) {
    return("SNP (0, 1, 2)")
  }
}

}
# Example usage:


#coding_schemes <- detect_genomic_coding(genomic_data)

