% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{validate_multi_trait_input_standard}
\alias{validate_multi_trait_input_standard}
\title{validate_multi_trait_input_standard}
\usage{
validate_multi_trait_input_standard(pheno_data = NULL, pheno_data_train = NULL, pheno_data_test = NULL, 
    geno_data = NULL, omic1_data = NULL, omic2_data = NULL, omic3_data = NULL, 
    gmatrix = NULL, gkernel = NULL, omic1_kernel = NULL, omic2_kernel = NULL, 
    omic3_kernel = NULL, kernel_list = NULL, response = NULL, gen_name = NULL, 
    response_family = "auto", multi_trait_asreml = FALSE, multi_trait_ml = FALSE, 
    multi_trait_dl = FALSE, multi_trait_gp = FALSE, multi_trait_bayes = FALSE, 
    heter_groups = NULL, 
    GS_model = NULL, GS_model_cv = NULL, cross_validation = FALSE,
    cross_validation_meth = NULL) 
}
\arguments{
\item{pheno_data}{Phenotype data frame.}

\item{pheno_data_train}{Parameter passed to the pheno_data_train workflow option.}

\item{pheno_data_test}{Parameter passed to the pheno_data_test workflow option.}

\item{geno_data}{Genotype data frame or matrix.}

\item{omic1_data}{Parameter passed to the omic1_data workflow option.}

\item{omic2_data}{Parameter passed to the omic2_data workflow option.}

\item{omic3_data}{Parameter passed to the omic3_data workflow option.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{gkernel}{Optional genomic kernel matrix.}

\item{omic1_kernel}{Optional first omic kernel matrix.}

\item{omic2_kernel}{Optional second omic kernel matrix.}

\item{omic3_kernel}{Optional third omic kernel matrix.}

\item{kernel_list}{Optional named list of user-supplied kernels.}

\item{response}{Name of the response column.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{response_family}{Response family used by the model workflow.}

\item{multi_trait_asreml}{Parameter passed to the multi_trait_asreml workflow option.}

\item{multi_trait_ml}{Parameter passed to the multi_trait_ml workflow option.}

\item{multi_trait_dl}{Parameter passed to the multi_trait_dl workflow option.}

\item{multi_trait_gp}{Logical. If \code{TRUE}, validate the dedicated joint
Gaussian-process multi-trait workflow.}

\item{multi_trait_bayes}{Logical. If \code{TRUE}, validate the dedicated joint
Bayesian BGLR multi-trait workflow.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{GS_model}{Optional selected true-prediction model.}

\item{GS_model_cv}{Optional selected cross-validation model.}

\item{cross_validation}{Logical indicating whether cross-validation is active.}

\item{cross_validation_meth}{Optional cross-validation method name. MT-MET
workflows require CV0, CV1, CV2, or a repeated variant.}

}
\description{
Documentation for \code{validate_multi_trait_input_standard}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
