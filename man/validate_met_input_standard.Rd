% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{validate_met_input_standard}
\alias{validate_met_input_standard}
\title{validate_met_input_standard}
\usage{
validate_met_input_standard(pheno_data = NULL, pheno_data_train = NULL, pheno_data_test = NULL, 
    geno_data = NULL, omic1_data = NULL, omic2_data = NULL, omic3_data = NULL, 
    gmatrix = NULL, omic1_kernel = NULL, omic2_kernel = NULL, omic3_kernel = NULL, 
    kernel_list = NULL, response = NULL, response_family = "gaussian",
    gen_name = NULL, heter_groups = NULL,
    met_ml_dl = FALSE, GS_model = NULL, GS_model_cv = NULL, 
    cross_validation = FALSE, cross_validation_meth = NULL) 
}
\arguments{
\item{pheno_data}{Phenotype data frame.}

\item{pheno_data_train}{Parameter passed to the pheno_data_train workflow option.}

\item{pheno_data_test}{Parameter passed to the pheno_data_test workflow option.}

\item{geno_data}{Genotype data frame or matrix.}

\item{omic1_data}{Parameter passed to the omic1_data workflow option.}

\item{omic2_data}{Parameter passed to the omic2_data workflow option.}

\item{omic3_data}{Parameter passed to the omic3_data workflow option.}

\item{gmatrix}{Genomic relationship or similarity kernel matrix keyed by \code{gen_name}.}

\item{omic1_kernel}{Optional omic1 relationship or similarity kernel matrix keyed by \code{gen_name}.}

\item{omic2_kernel}{Optional omic2 relationship or similarity kernel matrix keyed by \code{gen_name}.}

\item{omic3_kernel}{Optional omic3 relationship or similarity kernel matrix keyed by \code{gen_name}.}

\item{kernel_list}{Optional named list of relationship or similarity kernel matrices keyed by \code{gen_name}.}

\item{response}{Name of the response column.}

\item{response_family}{Response family for the MET request: \code{"gaussian"},
\code{"binary"}, \code{"multiclass"}, or \code{"ordinal"}. Model
applicability is checked against this family before fitting.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{met_ml_dl}{Parameter passed to the met_ml_dl workflow option.}

\item{GS_model}{Genomic selection model identifier.}

\item{GS_model_cv}{Parameter passed to the GS_model_cv workflow option.}

\item{cross_validation}{Parameter passed to the cross_validation workflow option.}

\item{cross_validation_meth}{Parameter passed to the cross_validation_meth workflow option.}

}
\description{
Validates a single-trait multi-environment request. Every requested model must
have an explicit MET implementation; mixed supported and unsupported model
lists are rejected before fitting rather than partially executed.
}
