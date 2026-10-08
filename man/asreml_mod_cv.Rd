% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{asreml_mod_cv}
\alias{asreml_mod_cv}
\title{asreml_mod_cv}
\usage{
asreml_mod_cv(pheno_data, gen_name, heter_groups, asreml_models_prep_cv, 
    response, tst, weights = NULL, workspace = NULL, pworkspace = NULL, maxit = 50)
}
\arguments{
\item{pheno_data}{Phenotype data frame.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{asreml_models_prep_cv}{Parameter passed to the asreml_models_prep_cv workflow option.}

\item{response}{Name of the response column.}

\item{tst}{Parameter passed to the tst workflow option.}

\item{weights}{Optional positive Stage 2 observation precisions. ASReml uses
them as a data-column weight with Gaussian dispersion fixed at one.}

\item{workspace}{Optional ASReml workspace setting.}

\item{pworkspace}{Optional ASReml prediction workspace setting.}

\item{maxit}{Maximum number of ASReml iterations.}

}
\description{
Documentation for \code{asreml_mod_cv}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
