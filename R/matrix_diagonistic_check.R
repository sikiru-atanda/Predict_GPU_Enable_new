matrix_diagonistic_check <- function(x, matrix_diagnostic = c(
                                        "is_square",
                                        "is_symmetric",
                                        "is_singular",
                                        "is_positive_definite",
                                        "is_positive_semi_definite"),
                                        tol = 1e-8) {

  msg <- ""

  # Check if the input is a matrix
  if (!is.matrix(x)) {
    stop(paste(msg, "Input must be a matrix."))
  }

  # Check if the matrix is numeric
  if (!is.numeric(x)) {
    stop(paste(msg, "Matrix must be numeric."))
  }

  # Check if the matrix is square for the relevant checks
  if (any(matrix_diagnostic %in% c("is_symmetric", "is_singular", "is_positive_definite", "is_positive_semi_definite"))) {
    if (nrow(x) != ncol(x)) {
      stop(paste(msg, "Matrix is not square, which is required for the requested diagnostics."))
    }
  }

  # Check if the matrix is square
  if ("is_square" %in% matrix_diagnostic) {
    if (nrow(x) != ncol(x)) {
      stop(paste(msg, "Matrix is not square."))
    } else {
      return(TRUE)
    }
  }

  # Check if the matrix is symmetric
  if ("is_symmetric" %in% matrix_diagnostic) {
    if (!identical(x, t(x))) {
      stop(paste(msg, "Matrix is not symmetric."))
    } else {
      return(TRUE)
    }
  }

  # Check if the matrix is singular
  if ("is_singular" %in% matrix_diagnostic) {
    det_val <- det(x)
    #if (abs(det_val) < .Machine$double.eps) {
    if (abs(det_val) < tol) {
      stop(paste(msg, "Matrix is singular."))
    } else {
      return(FALSE)  # Non-singular matrix
    }
  }

  # Check if the matrix is positive definite
  if ("is_positive_definite" %in% matrix_diagnostic) {
    pos_def_check <- tryCatch({
      chol(x)
      TRUE
    }, error = function(e) {
      FALSE
    })

    if (!pos_def_check) {
      return(FALSE)
      #stop(paste(msg, "Matrix is not positive definite."))
    } else {
      return(TRUE)
    }
  }

  # Check if the matrix is positive semi-definite
  if ("is_positive_semi_definite" %in% matrix_diagnostic) {
    eigenvalues <- eigen(x, only.values = TRUE)$values
    if (all(eigenvalues >= 0)) {
      return(TRUE)
    } else {
      return(FALSE)
      #stop(paste(msg, "Matrix is not positive semi-definite."))
    }
  }

  return(NULL)
}
