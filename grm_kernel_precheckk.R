#' Title
#'
#' @param grm_kernel_data
#' @param pedigree_matrix
#' @param bending
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param rcn_cutoff
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
grm_kernel_precheckOLD <- function(
                                grm_kernel_data = NULL,
                                pedigree_matrix = NULL,
                                bending = TRUE,
                                bend_value = 0.01,
                                blending = FALSE,
                                blending_value = 0.02,
                                high_diag_cut_off = 1.2,
                                low_diag_cut_off = 0.8,
                                duplicate_cut_off = 0.95,
                                rcn_cutoff = 1e-12,
                                optimize_diagonal = FALSE,
                                optimize_duplicate = FALSE,
                                message = TRUE,
                                ...
                            ) {
  msg <- sprintf("==================================================\n")

  check_symmetry <- function(matrix) {
    if (!isSymmetric.matrix(matrix)) {
      message(insight::print_color(paste(msg, paste("Relationsip Matrix is not symmetric'. We fix it.")), "blue"))
      matrix <- Matrix::forceSymmetric(matrix)
    }
    return(matrix)
  }

  if (!is.null(grm_kernel_data)) {
    stopifnot(!anyNA(grm_kernel_data), !is.null(rownames(grm_kernel_data)), !is.null(colnames(grm_kernel_data)))
    grm_kernel_data <- as.matrix(grm_kernel_data)
    stopifnot(identical(colnames(grm_kernel_data), rownames(grm_kernel_data)))

    grm_kernel_data <- check_symmetry(grm_kernel_data)

    if (bending && !is.positive.definite(matrixcalc::is.positive.definite(grm_kernel_data))) {
      message(insight::print_color(paste(msg, paste("Relationsip Matrix is not positive definite. We fix it.")), "blue"))
      grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data, posd.tol = bend_value, trace = FALSE)$mat)
    }
  }

  ### Either blending is TRUE or FALSE, this step is important
  if (blending && !isFALSE(blending)) {
    grm_kernel_data <- blending_stat(grm_kernel_data = grm_kernel_data,
                                     pedigree_matrix = pedigree_matrix,
                                     blending = blending,
                                     blending_value = blending_value)
  }

  res <- grm_kernel_diagnostic_fix(
    grm_kernel_data = grm_kernel_data,
    high_diag_cut_off = high_diag_cut_off,
    low_diag_cut_off = low_diag_cut_off,
    duplicate_cut_off = duplicate_cut_off,
    optimize_diagonal = optimize_diagonal,
    optimize_duplicate = optimize_duplicate
  )

  rcn <- rcond(res$clean_matrix)

  if ("potential_off_diag_with_duplicate" %in% names(res) || rcn < rcn_cutoff) {
    grm_kernel_data <- res$clean_matrix
    ncol_nrow <- ncol(grm_kernel_data)
    grm_kernel_data <- (1 - blending_value) * grm_kernel_data + blending_value * diag(x = 1, nrow = ncol_nrow, ncol = ncol_nrow)

    res <- grm_kernel_diagnostic_fix(
      grm_kernel_data = grm_kernel_data,
      high_diag_cut_off = high_diag_cut_off,
      low_diag_cut_off = low_diag_cut_off,
      duplicate_cut_off = duplicate_cut_off,
      optimize_diagonal = optimize_diagonal,
      optimize_duplicate = optimize_duplicate
    )

    rcn <- rcond(grm_kernel_data)

    if ("potential_off_diag_with_duplicate" %in% names(res) || rcn < rcn_cutoff) {
      message(paste(insight::print_color("WARNINGS\n", "blue"),
                    insight::print_color(paste(msg, paste("Matrix contain duplicate(s) or still ill-conditioned which might be potential problem.\n \t Change the blending value eg. 0.05  etc.")), "blue")))
      ncol_nrow <- ncol(grm_kernel_data)
      grm_kernel_data <- (1 - blending_value) * grm_kernel_data + blending_value * diag(x = 1, nrow = ncol_nrow, ncol = ncol_nrow)
    } else {
      if (isTRUE(message)) {
        message(paste(insight::print_color("WARNINGS\n", "blue"),
                      insight::print_color(paste(msg, paste("Matrix contain duplicate(s) which might be potential problem.\n \t We fix it by blending using an identity matrix.")), "blue")))
      }
    }
  }

  if (exists("res")) {
    rm(res)
  }

  class(grm_kernel_data) <- c("matrix", "array", "krm_data")
  attr(grm_kernel_data, "cleared") <- "for_model_fit"
  return(grm_kernel_data)
}
