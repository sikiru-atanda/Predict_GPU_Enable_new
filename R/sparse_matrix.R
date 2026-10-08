

#' Convert a dense (kernel) matrix to ASReml-style sparse triplet form
#'
#' Returns the lower triangle (including the diagonal) of a dense GRM / kernel
#' matrix in the three-column `(Row, Column, value)` sparse representation that
#' ASReml's `vm()` / `ginverse` machinery consumes. The output carries an
#' `INVERSE` attribute when the input represents an inverse relationship.
#'
#' @param grm_kernel_data Square numeric (G)RM or kernel matrix.
#' @param inverse Logical-like flag stamped onto the result as
#'   `attr(, "INVERSE")`; set to `TRUE` when the input is the inverse
#'   relationship matrix.
#' @param drop_zero Logical; when `TRUE` (default), zero entries are dropped
#'   from the sparse representation.
#'
#' @return A three-column matrix / data frame (`Row`, `Column`, `value`)
#'   representing the lower triangle of `grm_kernel_data`, with the `INVERSE`
#'   attribute set from `inverse`.
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
