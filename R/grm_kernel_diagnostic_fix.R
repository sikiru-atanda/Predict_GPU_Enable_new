

#' Title
#' Diagnose and Optimize Genomic Relationship Matrix (GRM) Kernel
#'
#' This function performs diagnostic checks on a genomic relationship matrix (GRM) kernel and
#' applies optimizations based on specified cutoff values for diagonal elements and correlation
#' thresholds for potential duplicates. It allows for the exclusion of outliers and potential
#' duplicate individuals to ensure the quality of the genomic relationship matrix.
#'
#' @param grm_kernel_data Numeric matrix, the genomic relationship matrix (GRM) kernel to be diagnosed.
#' @param high_diag_cut_off Numeric, the upper threshold for the diagonal elements of the GRM kernel.
#'        Diagonal elements above this threshold may indicate outliers and can be excluded.
#' @param low_diag_cut_off Numeric, the lower threshold for the diagonal elements of the GRM kernel.
#'        Diagonal elements below this threshold may indicate outliers and can be excluded.
#' @param duplicate_cut_off Numeric, the correlation threshold for identifying potential duplicate
#'        individuals within the GRM kernel. Pairs with a correlation above this threshold are
#'        considered potential duplicates.
#' @param optimize_diagonal Logical, if TRUE, the function will exclude diagonal elements outside
#'        the specified upper and lower cutoffs.
#' @param optimize_duplicate Logical, if TRUE, the function will exclude individuals identified
#'        as potential duplicates based on the correlation threshold.
#'
#' @return A list containing the optimized or cleaned GRM kernel matrix and optional diagnostic
#'         information about potential duplicates and diagonal elements to be removed. The list
#'         contains the following elements:
#'         - `clean_matrix`: The optimized GRM kernel matrix.
#'         - `potential_off_diag_with_duplicate`: A data frame of potential duplicates, if any, and their correlations.
#'         - `potential_diag_to_remove`: A data frame of diagonal elements suggested for removal, if any.
#'
#' @examples
#' # Example GRM kernel matrix
#' grm <- matrix(rnorm(100), ncol=10)
#' diag(grm) <- runif(10, 0.8, 1.2) # Example diagonal elements
#' # Diagnose and optimize the GRM kernel
#' result <- grm_kernel_diagnostic_fix(grm_kernel_data = grm,
#'                                     high_diag_cut_off = 1.2,
#'                                     low_diag_cut_off = 0.8,
#'                                     duplicate_cut_off = 0.95,
#'                                     optimize_diagonal = TRUE,
#'                                     optimize_duplicate = TRUE)
#' print(result$clean_matrix)
#'
#'This process was infer from ASRgenomics.
#'It was modified and improved to suite the objective in this package
#' @export
grm_kernel_diagnostic_fix <- function(grm_kernel_data = NULL,
                              high_diag_cut_off = 1.2,
                              low_diag_cut_off = 0.8,
                              duplicate_cut_off = 0.95,
                              optimize_diagonal = FALSE,
                              optimize_duplicate = FALSE
                              ){


  msg <- sprintf("==================================================\n")

  grm_kernel_data_opti <-  NULL
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

  potential_duplicate <- NULL
  diag_element_remove <- NULL
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

## Remove the diagonal; element based on the threshold defined by the user
if (isTRUE(optimize_diagonal) & nrow(diag_element_remove) > 0){
  #if (nrow(diag_element_remove) > 0){
    if(!is.null(grm_kernel_data_opti)){
      grm_kernel_data_opti <- grm_kernel_data_opti[-which(rownames(grm_kernel_data_opti) %in% row.names(diag_element_remove)),
                -which(rownames(grm_kernel_data_opti) %in% row.names(grm_kernel_data_opti))]

    } else {
      grm_kernel_data_opti <- grm_kernel_data[-which(rownames(grm_kernel_data) %in% row.names(diag_element_remove)),
                                                   -which(rownames(grm_kernel_data) %in% row.names(grm_kernel_data))]
    }

}


  if((!is.null(potential_duplicate) & is.null(diag_element_remove)) &  nrow(potential_duplicate)>0){
  res <-  list(clean_matrix = grm_kernel_data_opti,
               potential_off_diag_with_duplicate = potential_duplicate)


  } else if((!is.null(diag_element_remove) & is.null(potential_duplicate)) &  nrow(diag_element_remove)>0){

    res <-  list(lean_matrix = grm_kernel_data_opti,
               potential_diag_to_remove = diag_element_remove)

  } else if((!is.null(diag_element_remove) & !is.null(potential_duplicate)) &  ((nrow(diag_element_remove)>0) & (nrow(potential_duplicate)>0))){

    res <-  list(clean_matrix = grm_kernel_data_opti,
               potential_off_diag_with_duplicate = potential_duplicate,
               potential_diag_to_remove = diag_element_remove)

  } else {

    res <-  list(clean_matrix = grm_kernel_data)
  }

  return(res)

}
