

#' Title
#'
#' Inference from ASRgenomics and modified to the suitability of objective in this package
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
blending_statOLD <- function(grm_kernel_data = NULL,
         pedigree_matrix = NULL,
         blending = TRUE,
         blending_value = 0.02
         ){

  msg <- "\n==================================================\n"

  if(blending_value>=1){
    stop(print(paste(msg,'Consider lower value for blending')), call. = FALSE)
    #message(paste(msg,"Consider lower value for blending."))
  }
ncol_nrow = ncol(grm_kernel_data)
if (isTRUE(blending)){
  # Performing blend with I if requested
  if (is.null(pedigree_matrix)) {
    grm_kernel_data <- (1-blending_value)*grm_kernel_data + blending_value*diag(x=1, nrow=ncol_nrow , ncol=ncol_nrow )
    if (isTRUE(message)) {
      message(paste(msg,"Matrix was BLENDED using an identity matrix."))

    }
  }
  # Performing blend with pedigree_matrix if provided by the user
  if (!is.null(pedigree_matrix)) {
    if (!identical(colnames(grm_kernel_data), colnames(pedigree_matrix)) & !identical(rownames(grm_kernel_data), rownames(pedigree_matrix))) {
      stop(print(paste(msg,'colnames and rownames did not match for pedigree_matrix ang genomic/omic relationsip matrix')), call. = FALSE)}
    ## To ensure the order of the genotype/individual in both matrix are equivalent
    pedigree_matrix = pedigree_matrix[rownames(grm_kernel_data), colnames(grm_kernel_data)]

    grm_kernel_data <- (1-blending_value)*grm_kernel_data + blending_value*pedigree_matrix
      if (isTRUE(message)) {
        message(paste(msg,"Matrix was BLENDED using pedigree_matrix."))

      }

  }
}

return(grm_kernel_data)

}
