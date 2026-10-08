% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{grm_kernel_diagnostic_fixOLD}
\alias{grm_kernel_diagnostic_fixOLD}
\title{Legacy GRM/kernel diagnostic repair helper}
\usage{
grm_kernel_diagnostic_fixOLD(grm_kernel_data = NULL, high_diag_cut_off = 1.2, low_diag_cut_off = 0.8,
    duplicate_cut_off = 0.95, optimize_diagonal = FALSE, optimize_duplicate = FALSE)
}
\arguments{
\item{grm_kernel_data}{Genomic relationship or kernel matrix to diagnose and repair.}

\item{high_diag_cut_off}{Upper diagonal threshold used by the legacy diagnostic.}

\item{low_diag_cut_off}{Lower diagonal threshold used by the legacy diagnostic.}

\item{duplicate_cut_off}{Pairwise duplicate threshold used by the legacy diagnostic.}

\item{optimize_diagonal}{Logical; request diagonal optimization in the legacy repair path.}

\item{optimize_duplicate}{Logical; request duplicate optimization in the legacy repair path.}
}
\description{
Legacy exported wrapper for GRM/kernel diagnostic repair.
}
