% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{bayes_finalize_cross_val}
\alias{bayes_finalize_cross_val}
\title{bayes_finalize_cross_val}
\usage{
bayes_finalize_cross_val(fixed = NULL, random = NULL, GS_model = NULL, response = NULL, 
    weights = NULL, fixed_term_model_bayesian = NULL, rand_term_model_bayesian = NULL, 
    pheno_data = NULL, geno_data = NULL, omic1_data = NULL, omic2_data = NULL, 
    omic3_data = NULL, gen_name = NULL, nIter = NULL, burnIn = NULL, 
    thin = NULL, omics_data_label = NULL, ...) 
}
\arguments{
\item{fixed}{Fixed-effect formula or terms.}

\item{random}{Random-effect formula or terms.}

\item{GS_model}{Genomic selection model identifier.}

\item{response}{Name of the response column.}

\item{weights}{Optional observation weights.}

\item{fixed_term_model_bayesian}{Parameter passed to the fixed_term_model_bayesian workflow option.}

\item{rand_term_model_bayesian}{Parameter passed to the rand_term_model_bayesian workflow option.}

\item{pheno_data}{Phenotype data frame.}

\item{geno_data}{Genotype data frame or matrix.}

\item{omic1_data}{Parameter passed to the omic1_data workflow option.}

\item{omic2_data}{Parameter passed to the omic2_data workflow option.}

\item{omic3_data}{Parameter passed to the omic3_data workflow option.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{nIter}{Number of MCMC iterations.}

\item{burnIn}{Number of burn-in iterations.}

\item{thin}{MCMC thinning interval.}

\item{omics_data_label}{Parameter passed to the omics_data_label workflow option.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{bayes_finalize_cross_val}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
