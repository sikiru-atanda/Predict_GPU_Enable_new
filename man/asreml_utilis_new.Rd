% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{asreml_utilis_new}
\alias{asreml_utilis_new}
\title{asreml_utilis_new}
\usage{
asreml_utilis_new(fixed = NULL, random = NULL, cova = NULL, GS_model = NULL, 
    response = NULL, pheno_data = NULL, gmatrix = NULL, gkernel = NULL, 
    omic1_kernel = NULL, omic2_kernel = NULL, omic3_kernel = NULL, 
    kernel_list = NULL, inverse = NULL, epsilon = TRUE, gen_name = NULL, heter_groups = NULL, 
    heter_resid = FALSE, var_cov_str = NULL, weights = NULL, 
    workspace = 1e+08, engine = NULL, pworkspace = 1e+06, maxit = 50, 
    cross_validation = FALSE, ...) 
}
\arguments{
\item{fixed}{Fixed-effect formula or terms.}

\item{random}{Random-effect formula or terms.}

\item{cova}{Parameter passed to the cova workflow option.}

\item{GS_model}{Genomic selection model identifier.}

\item{response}{Name of the response column.}

\item{pheno_data}{Phenotype data frame.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{gkernel}{Parameter passed to the gkernel workflow option.}

\item{omic1_kernel}{Parameter passed to the omic1_kernel workflow option.}

\item{omic2_kernel}{Parameter passed to the omic2_kernel workflow option.}

\item{omic3_kernel}{Parameter passed to the omic3_kernel workflow option.}

\item{kernel_list}{Optional named list of additional relationship or kernel matrices.}

\item{inverse}{Parameter passed to the inverse workflow option.}

\item{epsilon}{Parameter passed to the epsilon workflow option.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{heter_resid}{Parameter passed to the heter_resid workflow option.}

\item{var_cov_str}{Parameter passed to the var_cov_str workflow option.}

\item{weights}{Optional positive Stage 2 observation precisions, supplied as
a numeric vector, one-column table, or column name in \code{pheno_data}.
ASReml uses an internal data column and
\code{asr_gaussian(dispersion = 1)}, so the residual variance for row
\code{i} is exactly \code{1 / weights[i]}.}

\item{workspace}{Parameter passed to the workspace workflow option.}

\item{engine}{Parameter passed to the engine workflow option.}

\item{pworkspace}{Parameter passed to the pworkspace workflow option.}

\item{maxit}{Parameter passed to the maxit workflow option.}

\item{cross_validation}{Parameter passed to the cross_validation workflow option.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{asreml_utilis_new}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
