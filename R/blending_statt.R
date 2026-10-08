#' Blend Genomic Relationship Matrix (GRM) with Identity or Pedigree Matrix
#'
#' This function blends a genomic relationship matrix (GRM) with either an identity matrix or a pedigree matrix,
#' based on the blending proportion specified by the user. This is often used to adjust the GRM for better
#' performance in genetic analyses.
#'
#' @param grm_kernel_data A square numeric matrix representing the genomic relationship matrix (GRM) of individuals.
#' @param pedigree_matrix An optional square numeric matrix representing the pedigree relationship matrix of individuals.
#'                        If provided, this matrix is used for blending with the GRM. It must have the same dimensions
#'                        and the row and column names as `grm_kernel_data`.
#' @param blending Logical, if `TRUE`, the blending operation is performed. If `FALSE`, the function returns the original
#'                 `grm_kernel_data` without any modification.
#' @param blending_value A numeric value between 0 and 1 specifying the proportion of the identity matrix or
#'                       pedigree matrix to blend with the GRM. The default value is 0.02.
#'
#' @return Returns a matrix which is the result of blending the input GRM with either an identity matrix or a pedigree
#'         matrix, depending on the input parameters.
#'
#' @examples
#' # Generate a mock GRM matrix
#' grm <- matrix(runif(100), ncol=10)
#' rownames(grm) <- colnames(grm) <- paste0("Ind", 1:10)
#'
#' # Blend the GRM with an identity matrix
#' blended_grm <- blending_stat(grm_kernel_data = grm, blending = TRUE, blending_value = 0.05)
#'
#' # Generate a mock pedigree matrix
#' pedigree <- matrix(runif(100), ncol=10)
#' rownames(pedigree) <- colnames(pedigree) <- paste0("Ind", 1:10)
#'
#' # Blend the GRM with the pedigree matrix
#' blended_grm_pedigree <- blending_stat(grm_kernel_data = grm, pedigree_matrix = pedigree,
#'                                       blending = TRUE, blending_value = 0.05)
#'
#' @export
#' @importFrom stats diag
#'
blending_stat <- function(grm_kernel_data = NULL,
                          pedigree_matrix = NULL,
                          blending = TRUE,
                          blending_value = 0.02
) {
  msg <- ""

  if (blending_value >= 1) {
    stop(msg, 'Consider lower value for blending')
  }

  ncol_nrow <- ncol(grm_kernel_data)

  if (isTRUE(blending) && is.null(pedigree_matrix)) {
    # Blend with I if requested
    grm_kernel_data <- (1 - blending_value) * grm_kernel_data + blending_value * diag(x = 1, nrow = ncol_nrow, ncol = ncol_nrow)

    if (isTRUE(message)) {
      message(paste(msg, "Matrix was BLENDED using an identity matrix."))
    }
  }

  if (blending && !is.null(pedigree_matrix)) {
    # Blend with pedigree_matrix if provided
    if (!identical(colnames(grm_kernel_data), colnames(pedigree_matrix)) ||
        !identical(rownames(grm_kernel_data), rownames(pedigree_matrix))) {
      stop(msg, 'colnames and rownames did not match for pedigree_matrix and genomic/omic relationship matrix')
    }

    # Ensure the order of genotypes/individuals in both matrices is equivalent
    pedigree_matrix <- pedigree_matrix[rownames(grm_kernel_data), colnames(grm_kernel_data)]

    grm_kernel_data <- (1 - blending_value) * grm_kernel_data + blending_value * pedigree_matrix

    if (isTRUE(message)) {
      message(paste(msg, "Matrix was BLENDED using pedigree_matrix."))
    }
  }

  return(grm_kernel_data)
}
