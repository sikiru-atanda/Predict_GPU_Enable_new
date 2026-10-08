#' Diagnose a relationship or kernel matrix
#'
#' Checks the structural conditions required when a matrix is used as a
#' genomic relationship or kernel covariance: finite numeric entries, square
#' dimensions, symmetry, matching unique sample names, and positive
#' semidefiniteness.
#'
#' @param x Numeric relationship or kernel matrix.
#' @param tolerance Relative numerical tolerance used for symmetry and
#'   eigenvalue checks.
#' @param require_names Logical; require matching row and column names.
#'
#' @return A list with a `valid` flag, numerical diagnostics, and any issues.
#' @export
relationship_matrix_diagnostics <- function(x,
                                            tolerance = sqrt(.Machine$double.eps),
                                            require_names = TRUE) {
  issues <- character()
  if (!is.matrix(x)) {
    x <- tryCatch(as.matrix(x), error = function(e) NULL)
  }
  if (is.null(x) || !is.numeric(x)) {
    return(list(
      valid = FALSE,
      dimension = c(NA_integer_, NA_integer_),
      symmetric_error = NA_real_,
      minimum_eigenvalue = NA_real_,
      maximum_eigenvalue = NA_real_,
      effective_rank = NA_integer_,
      condition_number = NA_real_,
      issues = "matrix must be coercible to a numeric matrix"
    ))
  }

  nr <- nrow(x)
  nc <- ncol(x)
  if (nr != nc) {
    issues <- c(issues, "matrix must be square")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    issues <- c(issues, "matrix entries must all be finite")
  }

  rn <- rownames(x)
  cn <- colnames(x)
  if (isTRUE(require_names)) {
    if (is.null(rn) || is.null(cn)) {
      issues <- c(issues, "matching row and column names are required")
    } else {
      if (!identical(rn, cn)) {
        issues <- c(issues, "row and column names must match in the same order")
      }
      if (anyDuplicated(rn) || anyDuplicated(cn)) {
        issues <- c(issues, "sample names must be unique")
      }
      if (anyNA(rn) || anyNA(cn) || any(!nzchar(rn)) || any(!nzchar(cn))) {
        issues <- c(issues, "sample names must be non-missing and non-empty")
      }
    }
  }

  can_decompose <- nr == nc && !anyNA(x) && all(is.finite(x))
  symmetric_error <- if (can_decompose) max(abs(x - t(x))) else NA_real_
  matrix_scale <- if (can_decompose) max(1, max(abs(x))) else 1
  symmetry_tol <- abs(as.numeric(tolerance)[1L]) * matrix_scale
  if (can_decompose && symmetric_error > symmetry_tol) {
    issues <- c(issues, "matrix must be symmetric within tolerance")
  }

  eigenvalues <- if (can_decompose && symmetric_error <= symmetry_tol) {
    tryCatch(
      eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) numeric()
    )
  } else {
    numeric()
  }
  minimum_eigenvalue <- if (length(eigenvalues)) min(eigenvalues) else NA_real_
  maximum_eigenvalue <- if (length(eigenvalues)) max(eigenvalues) else NA_real_
  eigen_tol <- abs(as.numeric(tolerance)[1L]) * max(1, abs(maximum_eigenvalue))
  if (length(eigenvalues) && minimum_eigenvalue < -eigen_tol) {
    issues <- c(issues, "matrix is not positive semidefinite within tolerance")
  }
  positive_eigenvalues <- eigenvalues[eigenvalues > eigen_tol]
  effective_rank <- if (length(eigenvalues)) length(positive_eigenvalues) else NA_integer_
  condition_number <- if (length(positive_eigenvalues)) {
    max(positive_eigenvalues) / min(positive_eigenvalues)
  } else {
    NA_real_
  }

  list(
    valid = length(issues) == 0L,
    dimension = c(nr, nc),
    symmetric_error = symmetric_error,
    minimum_eigenvalue = minimum_eigenvalue,
    maximum_eigenvalue = maximum_eigenvalue,
    effective_rank = as.integer(effective_rank),
    condition_number = condition_number,
    issues = unique(issues)
  )
}

#' Validate a relationship or kernel matrix
#'
#' @inheritParams relationship_matrix_diagnostics
#' @param positive_definite Logical; require strictly positive eigenvalues
#'   rather than allowing a positive-semidefinite matrix.
#'
#' @return Invisibly returns the diagnostics when validation succeeds.
#' @rdname relationship_matrix_diagnostics
#' @export
validate_relationship_matrix <- function(x,
                                         tolerance = sqrt(.Machine$double.eps),
                                         require_names = TRUE,
                                         positive_definite = FALSE) {
  diagnostics <- relationship_matrix_diagnostics(
    x = x,
    tolerance = tolerance,
    require_names = require_names
  )
  issues <- diagnostics$issues
  if (isTRUE(positive_definite) &&
      is.finite(diagnostics$minimum_eigenvalue) &&
      diagnostics$minimum_eigenvalue <= abs(tolerance)) {
    issues <- c(issues, "matrix must be positive definite")
  }
  if (length(issues)) {
    stop(
      paste("Invalid relationship/kernel matrix:", paste(unique(issues), collapse = "; ")),
      call. = FALSE
    )
  }
  invisible(diagnostics)
}

gp_factor_analytic_covariance <- function(loadings,
                                          specific_variances = NULL,
                                          tolerance = sqrt(.Machine$double.eps)) {
  loadings <- as.matrix(loadings)
  storage.mode(loadings) <- "double"
  if (!nrow(loadings) || !ncol(loadings) || anyNA(loadings) || any(!is.finite(loadings))) {
    stop("Factor loadings must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  if (is.null(specific_variances)) {
    specific_variances <- rep(0, nrow(loadings))
  }
  specific_variances <- suppressWarnings(as.numeric(specific_variances))
  if (length(specific_variances) != nrow(loadings) ||
      anyNA(specific_variances) || any(!is.finite(specific_variances))) {
    stop("Provide one finite specific variance per loading row.", call. = FALSE)
  }
  tol <- abs(as.numeric(tolerance)[1L])
  if (any(specific_variances < -tol)) {
    stop("Factor-analytic specific variances cannot be negative.", call. = FALSE)
  }
  specific_variances[specific_variances < 0] <- 0
  covariance <- tcrossprod(loadings) + diag(specific_variances, nrow = nrow(loadings))
  dimnames(covariance) <- list(rownames(loadings), rownames(loadings))
  covariance
}

gp_correlation_matrix_from_parameters <- function(parameters,
                                                  n,
                                                  compound_symmetry = FALSE,
                                                  labels = NULL,
                                                  tolerance = sqrt(.Machine$double.eps)) {
  parameters <- suppressWarnings(as.numeric(parameters))
  n <- as.integer(n)[1L]
  if (!is.finite(n) || n < 1L || anyNA(parameters) || any(!is.finite(parameters))) {
    stop("Correlation parameters and dimension must be finite.", call. = FALSE)
  }
  correlation <- diag(n)
  if (isTRUE(compound_symmetry)) {
    if (length(parameters) != 1L) {
      stop("Compound-symmetry covariance requires exactly one correlation parameter.", call. = FALSE)
    }
    correlation[lower.tri(correlation)] <- parameters[[1L]]
    correlation[upper.tri(correlation)] <- parameters[[1L]]
  } else {
    expected <- n * (n - 1L) / 2L
    if (length(parameters) != expected) {
      stop(
        sprintf("Expected %d correlation parameters for %d environments, found %d.",
                expected, n, length(parameters)),
        call. = FALSE
      )
    }
    cursor <- 1L
    for (r in seq_len(n)) {
      if (r <= 1L) next
      for (c in seq_len(r - 1L)) {
        correlation[r, c] <- parameters[[cursor]]
        correlation[c, r] <- parameters[[cursor]]
        cursor <- cursor + 1L
      }
    }
  }
  if (any(abs(correlation) > 1 + tolerance)) {
    stop("Correlation parameters must be in [-1, 1].", call. = FALSE)
  }
  if (!is.null(labels)) dimnames(correlation) <- list(labels, labels)
  diagnostics <- relationship_matrix_diagnostics(
    correlation,
    tolerance = tolerance,
    require_names = !is.null(labels)
  )
  if (!isTRUE(diagnostics$valid)) {
    stop(
      paste("Reconstructed correlation matrix is invalid:", paste(diagnostics$issues, collapse = "; ")),
      call. = FALSE
    )
  }
  correlation
}

gp_covariance_from_correlation <- function(variances,
                                           correlation,
                                           labels = NULL,
                                           tolerance = sqrt(.Machine$double.eps)) {
  variances <- suppressWarnings(as.numeric(variances))
  correlation <- as.matrix(correlation)
  n <- nrow(correlation)
  if (ncol(correlation) != n || length(variances) != n ||
      anyNA(variances) || any(!is.finite(variances))) {
    stop("Provide one finite variance per correlation-matrix row.", call. = FALSE)
  }
  if (any(variances < -abs(tolerance))) {
    stop("Variance parameters cannot be negative.", call. = FALSE)
  }
  variances[variances < 0] <- 0
  if (!is.null(labels)) dimnames(correlation) <- list(labels, labels)
  diagnostics <- relationship_matrix_diagnostics(
    correlation,
    tolerance = tolerance,
    require_names = !is.null(labels)
  )
  if (!isTRUE(diagnostics$valid) || any(abs(diag(correlation) - 1) > tolerance)) {
    stop("A valid correlation matrix with a unit diagonal is required.", call. = FALSE)
  }
  scale_matrix <- diag(sqrt(variances), nrow = n)
  covariance <- scale_matrix %*% correlation %*% scale_matrix
  if (!is.null(labels)) dimnames(covariance) <- list(labels, labels)
  covariance
}

gp_unstructured_covariance_from_parameters <- function(parameters,
                                                       n,
                                                       labels = NULL,
                                                       tolerance = sqrt(.Machine$double.eps)) {
  parameters <- suppressWarnings(as.numeric(parameters))
  expected <- n * (n + 1L) / 2L
  if (length(parameters) != expected || anyNA(parameters) || any(!is.finite(parameters))) {
    stop(
      sprintf("Expected %d finite unstructured covariance parameters for %d environments, found %d.",
              expected, n, length(parameters)),
      call. = FALSE
    )
  }
  covariance <- matrix(0, nrow = n, ncol = n)
  cursor <- 1L
  for (r in seq_len(n)) {
    for (c in seq_len(r)) {
      covariance[r, c] <- parameters[[cursor]]
      covariance[c, r] <- parameters[[cursor]]
      cursor <- cursor + 1L
    }
  }
  if (!is.null(labels)) dimnames(covariance) <- list(labels, labels)
  diagnostics <- relationship_matrix_diagnostics(
    covariance,
    tolerance = tolerance,
    require_names = !is.null(labels)
  )
  if (!isTRUE(diagnostics$valid)) {
    stop(
      paste("Reconstructed unstructured covariance matrix is invalid:",
            paste(diagnostics$issues, collapse = "; ")),
      call. = FALSE
    )
  }
  covariance
}
