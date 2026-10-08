% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{validate_hybrid_input_standard}
\alias{validate_hybrid_input_standard}
\title{validate_hybrid_input_standard}
\usage{
validate_hybrid_input_standard(pheno_data = NULL, pheno_data_train = NULL, pheno_data_test = NULL, 
    geno_data = NULL, gmatrix = NULL, female_gmatrix = NULL, 
    male_gmatrix = NULL, female_geno_data = NULL, male_geno_data = NULL, 
    response = NULL, gen_name = NULL, female_parent = NULL, male_parent = NULL, 
    heter_groups = NULL, hybrid_asreml = FALSE, hybrid_bayes = FALSE, 
    hybrid_gp = FALSE, hybrid_ml = FALSE, hybrid_dl = FALSE, 
    GS_model = NULL, GS_model_cv = NULL, cross_validation = FALSE, 
    cross_validation_meth = NULL) 
}
\arguments{
\item{pheno_data}{Phenotype data frame.}

\item{pheno_data_train}{Parameter passed to the pheno_data_train workflow option.}

\item{pheno_data_test}{Parameter passed to the pheno_data_test workflow option.}

\item{geno_data}{Genotype data frame or matrix.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{female_gmatrix}{Parameter passed to the female_gmatrix workflow option.}

\item{male_gmatrix}{Parameter passed to the male_gmatrix workflow option.}

\item{female_geno_data}{Parameter passed to the female_geno_data workflow option.}

\item{male_geno_data}{Parameter passed to the male_geno_data workflow option.}

\item{response}{Name of the response column.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{female_parent}{Parameter passed to the female_parent workflow option.}

\item{male_parent}{Parameter passed to the male_parent workflow option.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{hybrid_asreml}{Parameter passed to the hybrid_asreml workflow option.}

\item{hybrid_bayes}{Parameter passed to the hybrid_bayes workflow option.}

\item{hybrid_gp}{Parameter passed to the hybrid_gp workflow option.}

\item{hybrid_ml}{Parameter passed to the hybrid_ml workflow option.}

\item{hybrid_dl}{Parameter passed to the hybrid_dl workflow option.}

\item{GS_model}{Optional selected true-prediction model.}

\item{GS_model_cv}{Optional selected cross-validation model.}

\item{cross_validation}{Logical indicating whether cross-validation is active.}

\item{cross_validation_meth}{Optional cross-validation method. Hybrid MET workflows require CV0, CV1, CV2, or a repeated variant.}

}
\description{
Documentation for \code{validate_hybrid_input_standard}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
