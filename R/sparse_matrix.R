

#' Title
#' Convert dense matrix to sparse form
#'
#' @param grm_kernel_data
#' @param drop_zero
#'
#' @return
#' @export
#'
#' @examples
sparse_matrix <-  function(
    grm_kernel_data = NULL,
    inverse = NULL,
    drop_zero= TRUE){

if(isTRUE(drop_zero)) {
  which <- (grm_kernel_data != 0 & lower.tri(grm_kernel_data, diag = TRUE))
} else {
  which <- lower.tri(grm_kernel_data, diag = TRUE)
}

res_sparse <- data.frame(Row = t(row(grm_kernel_data))[t(which)],
                           Col = t(col(grm_kernel_data))[t(which)],
                           Value = t(grm_kernel_data)[t(which)])

res_sparse <- as.matrix(res_sparse)

# # Add attributes.
## This attribute will not work if grm_kernel_data did not have
## rowNames and ColNames as attributes
attr(res_sparse, "rowNames") <- rownames(grm_kernel_data)
attr(res_sparse, "colNames") <- colnames(grm_kernel_data)

return(res_sparse)

}
