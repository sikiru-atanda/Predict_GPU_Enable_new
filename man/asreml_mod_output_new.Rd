% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{asreml_mod_output_new}
\alias{asreml_mod_output_new}
\title{asreml_mod_output_new}
\usage{
asreml_mod_output_new(mod_asreml = NULL, pheno_data = NULL, response = NULL, 
    gmatrix = NULL, omic1_kernel = NULL, omic2_kernel = NULL, 
    omic3_kernel = NULL, kernel_list = NULL, omics_kernel_label = list(omic1_kernel = NULL, 
        omic2_kernel = NULL, omic3_kernel = NULL), heter_groups = NULL, 
    gen_name = NULL, var_cov_str = NULL, heter_resid = NULL, 
    workspace = NULL, pworkspace = 1e+15, maxit = 50, ...) 
}
\arguments{
\item{mod_asreml}{Parameter passed to the mod_asreml workflow option.}

\item{pheno_data}{Phenotype data frame.}

\item{response}{Name of the response column.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{omic1_kernel}{Parameter passed to the omic1_kernel workflow option.}

\item{omic2_kernel}{Parameter passed to the omic2_kernel workflow option.}

\item{omic3_kernel}{Parameter passed to the omic3_kernel workflow option.}

\item{kernel_list}{Optional named list of additional relationship or kernel matrices.}

\item{omics_kernel_label}{Parameter passed to the omics_kernel_label workflow option.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{var_cov_str}{Parameter passed to the var_cov_str workflow option.}

\item{heter_resid}{Parameter passed to the heter_resid workflow option.}

\item{workspace}{Optional ASReml workspace setting.}

\item{pworkspace}{Parameter passed to the pworkspace workflow option.}

\item{maxit}{Parameter passed to the maxit workflow option.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{asreml_mod_output_new}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
