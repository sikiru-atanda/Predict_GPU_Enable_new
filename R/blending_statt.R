#' Title
#'
#' @param grm_kernel_data
#' @param pedigree_matrix
#' @param blending
#' @param blending_value
#'
#' @return
#' @export
#'
#' @examples
blending_stat <- function(grm_kernel_data = NULL,
                          pedigree_matrix = NULL,
                          blending = TRUE,
                          blending_value = 0.02
) {
  msg <- sprintf("==================================================\n")

  if (blending_value >= 1) {
    stop(print(paste(msg, 'Consider lower value for blending')), call. = FALSE)
  }

  ncol_nrow <- ncol(grm_kernel_data)

  if (blending && is.null(pedigree_matrix)) {
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
      stop(print(paste(msg, 'colnames and rownames did not match for pedigree_matrix and genomic/omic relationship matrix')), call. = FALSE)
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
