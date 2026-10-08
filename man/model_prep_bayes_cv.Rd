% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{model_prep_bayes_cv}
\alias{model_prep_bayes_cv}
\title{model_prep_bayes_cv}
\usage{
model_prep_bayes_cv(fixed = NULL, random = NULL, GS_model_cv = NULL, response = NULL, 
    gen_name = NULL, pheno_data = NULL, weights = NULL, fixed_term_model_bayesian = NULL, 
    rand_term_model_bayesian = NULL, test_set = NULL, nIter = NULL, 
    burnIn = NULL, thin = NULL, geno_data = NULL, omic1_data = NULL, 
    omic2_data = NULL, omic3_data = NULL, omics_data_label = NULL, 
    gmatrix = NULL, omic1_kernel = NULL, omic2_kernel = NULL, 
    omic3_kernel = NULL, kernel_list = NULL, heter_groups = NULL, 
    heter_resid = FALSE, bayes_kernel_heter_resid = NULL,
    omics_kernel_label = NULL, cross_validation = TRUE,
    response_family = "gaussian") 
}
\arguments{
\item{fixed}{Fixed-effect formula or terms.}

\item{random}{Random-effect formula or terms.}

\item{GS_model_cv}{Parameter passed to the GS_model_cv workflow option.}

\item{response}{Name of the response column.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{pheno_data}{Phenotype data frame.}

\item{weights}{Optional observation weights.}

\item{fixed_term_model_bayesian}{Parameter passed to the fixed_term_model_bayesian workflow option.}

\item{rand_term_model_bayesian}{Parameter passed to the rand_term_model_bayesian workflow option.}

\item{test_set}{Test set identifiers or indices.}

\item{nIter}{Number of MCMC iterations.}

\item{burnIn}{Number of burn-in iterations.}

\item{thin}{MCMC thinning interval.}

\item{geno_data}{Genotype data frame or matrix.}

\item{omic1_data}{Parameter passed to the omic1_data workflow option.}

\item{omic2_data}{Parameter passed to the omic2_data workflow option.}

\item{omic3_data}{Parameter passed to the omic3_data workflow option.}

\item{omics_data_label}{Parameter passed to the omics_data_label workflow option.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{omic1_kernel}{Parameter passed to the omic1_kernel workflow option.}

\item{omic2_kernel}{Parameter passed to the omic2_kernel workflow option.}

\item{omic3_kernel}{Parameter passed to the omic3_kernel workflow option.}

\item{kernel_list}{Optional named list of additional relationship or kernel matrices.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{heter_resid}{Logical indicating whether Bayesian Gaussian fits should estimate residual variance by \code{heter_groups} when supported.}

\item{bayes_kernel_heter_resid}{Optional logical override for heterogeneous
residual routing in kernel-Bayesian models. \code{NULL} inherits
\code{heter_resid}.}

\item{omics_kernel_label}{Parameter passed to the omics_kernel_label workflow option.}

\item{cross_validation}{Parameter passed to the cross_validation workflow option.}

\item{response_family}{Response family used by the model workflow.}

}
\description{
Documentation for \code{model_prep_bayes_cv}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
