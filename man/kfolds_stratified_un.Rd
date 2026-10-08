% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{kfolds_stratified_un}
\alias{kfolds_stratified_un}
\title{kfolds_stratified_un}
\usage{
kfolds_stratified_un(pheno_data, response, gen_name, nfolds = 5, random_state = NULL, 
    sampling_method = c("stratified", "unstratified"), replication = 1, 
    message = TRUE, ...) 
}
\arguments{
\item{pheno_data}{Phenotype data frame.}

\item{response}{Name of the response column.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{nfolds}{Parameter passed to the nfolds workflow option.}

\item{random_state}{Random seed used for reproducible splits.}

\item{sampling_method}{Parameter passed to the sampling_method workflow option.}

\item{replication}{Number of cross-validation replications.}

\item{message}{Logical indicating whether progress messages should be printed.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{kfolds_stratified_un}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
