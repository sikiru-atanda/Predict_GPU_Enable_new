

#' Title
#' This process was infer from ASRgenomics. It was modified and improved to suite the
#' objectives in this package
#'
#' @param grm_kernel_data
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param message
#'
#' @return
#' @export
#'
#' @examples
grm_kernel_diagnostic_fix <- function(grm_kernel_data = NULL,
                              high_diag_cut_off = 1.2,
                              low_diag_cut_off = 0.8,
                              duplicate_cut_off = 0.95,
                              optimize_diagonal = FALSE,
                              optimize_duplicate = FALSE
                              ){


  msg <- sprintf("==================================================\n")
  # Check input value
  if (duplicate_cut_off < 0 | duplicate_cut_off > 1) {
    stop(print(paste(msg,'Duplicate threshould value must be between 0 and 1.')), call. = FALSE)

  }
  if (high_diag_cut_off < low_diag_cut_off) {
    stop(print(paste(msg,'Cut-off value for large diagonal(s)  must be equal to or greater than the cut-off value for small diagonal(s)')), call. = FALSE)

  }
  if (high_diag_cut_off < 0 | low_diag_cut_off < 0) {
    stop(print(paste(msg,'Cut-off value for large and small diagonal(s)  must be positive.')), call. = FALSE)

  }

  ### Start Process to remove duplicates.
  #n <- nrow(grm_kernel_data)
  #indNames <- rownames(K)
  # Generate vector of of the diagonal element of the matrix
  diag_grmkernel <- diag(grm_kernel_data)
  corr_grmkernel <- stats::cov2cor(grm_kernel_data)
  sparse_grmkernel <- sparse_matrix(grm_kernel_data)
  sparse_corr <- sparse_matrix(corr_grmkernel)
  sparse_grmkernel <- data.frame(sparse_grmkernel, Corr=sparse_corr[,3])
  off_diag <- sparse_grmkernel[sparse_grmkernel$Row != sparse_grmkernel$Col,]
  rm(sparse_grmkernel, corr_grmkernel, sparse_corr)


  # Generating potential duplicates to remove
  potential_duplicate <- off_diag[off_diag$Corr > duplicate_cut_off,]
  if(nrow(potential_duplicate)>0){
  potential_duplicate <- data.frame(Indiv_A=rownames(grm_kernel_data)[potential_duplicate$Row],
                                    Indiv_B=colnames(grm_kernel_data)[potential_duplicate$Col],
                                    potential_duplicate[, -c(1:2)]
                                    )
  rownames(potential_duplicate) <- NULL
  potential_duplicate <- potential_duplicate[order(potential_duplicate$Corr, decreasing=TRUE),]

}

# Generating potential diagonal elements to remove

diag_element_remove <- data.frame(
  value = sort(diag_grmkernel[diag_grmkernel > high_diag_cut_off
                     | diag_grmkernel < low_diag_cut_off], decreasing=TRUE))


## Remove the duplicate element based on the threshold defined by the user
if(isTRUE(optimize_duplicate) &  nrow(potential_duplicate)>0){
offdiag_element_remove <- unique(c(potential_duplicate$Indiv_A, potential_duplicate$Indiv_B))


grm_kernel_data_opti <- grm_kernel_data[-which(rownames(grm_kernel_data) %in% offdiag_element_remove),
                                                -which(rownames(grm_kernel_data) %in% offdiag_element_remove)]

}


## Remove the diagona; element based on the threshold defined by the user
if (isTRUE(optimize_diagonal) & nrow(diag_element_remove) > 0){
  #if (nrow(diag_element_remove) > 0){
    if(exists('grm_kernel_data_opti')){
      grm_kernel_data_opti <- grm_kernel_data_opti[-which(rownames(grm_kernel_data_opti) %in% row.names(diag_element_remove)),
                -which(rownames(grm_kernel_data_opti) %in% row.names(grm_kernel_data_opti))]

    } else {

      grm_kernel_data_opti <- grm_kernel_data[-which(rownames(grm_kernel_data) %in% row.names(diag_element_remove)),
                                                   -which(rownames(grm_kernel_data) %in% row.names(grm_kernel_data))]

    }


  #}

}


  if((exists("potential_duplicate") & !exists("diag_element_remove")) &  nrow(potential_duplicate)>0){
  res = list(grm_kernel_data_opti,
             potential_duplicate)

  names(res) = c("clean_matrix",
                 "potential_off_diag_with_duplicate")

  } else if((exists("diag_element_remove") & !exists("potential_duplicate")) &  nrow(diag_element_remove)>0){

    res = list(grm_kernel_data_opti,
               diag_element_remove)

    names(res) =  c("clean_matrix",
                    "potential_diag_to_remove")

  } else if((exists("diag_element_remove") & exists("potential_duplicate")) &  ((nrow(diag_element_remove)>0) & (nrow(potential_duplicate)>0))){

    res = list(grm_kernel_data_opti,
               potential_duplicate,
               diag_element_remove)

    names(res) =  c("clean_matrix",
                    "potential_off_diag_with_duplicate",
                    "potential_diag_to_remove")
  } else {

    res = list(grm_kernel_data)

    names(res) =  c("clean_matrix")

  }

  return(res)





}
