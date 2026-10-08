\name{relationship_matrix_diagnostics}
\alias{relationship_matrix_diagnostics}
\alias{validate_relationship_matrix}
\title{Diagnose and validate relationship or kernel matrices}
\usage{
relationship_matrix_diagnostics(
  x,
  tolerance = sqrt(.Machine$double.eps),
  require_names = TRUE
)

validate_relationship_matrix(
  x,
  tolerance = sqrt(.Machine$double.eps),
  require_names = TRUE,
  positive_definite = FALSE
)
}
\arguments{
\item{x}{Numeric relationship or kernel matrix.}
\item{tolerance}{Relative tolerance for symmetry and eigenvalue checks.}
\item{require_names}{Require unique, matching row and column sample names.}
\item{positive_definite}{Require strict positive definiteness instead of
positive semidefiniteness.}
}
\value{
The diagnostics function returns validity, dimension, symmetry error,
eigenvalue range, effective rank, condition number, and issues. The validator
returns the same object invisibly or stops on a violation.
}
\description{
Checks the covariance-matrix invariants required by genomic relationship and
multi-omic kernel models.
}
