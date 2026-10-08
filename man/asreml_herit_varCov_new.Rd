% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{asreml_herit_varCov_new}
\alias{asreml_herit_varCov_new}
\title{asreml_herit_varCov_new}
\usage{
asreml_herit_varCov_new(model = NULL, heter_groups = NULL, var_cov_str = NULL, 
    heter_resid = FALSE, names_in_inv_list = NULL, inter_gen_pos = NULL, 
    gen_pos = NULL, response = NULL, gen_name = NULL, ...)
}
\arguments{
\item{model}{Model object or model identifier.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{var_cov_str}{ASReml variance-covariance structure. Supported forms include
\code{"us"}, \code{"fa"}, \code{"rr"}, \code{"corgh"}, \code{"corh"},
and \code{"corv"}.}

\item{heter_resid}{Parameter passed to the heter_resid workflow option.}

\item{names_in_inv_list}{Parameter passed to the names_in_inv_list workflow option.}

\item{inter_gen_pos}{Parameter passed to the inter_gen_pos workflow option.}

\item{gen_pos}{Parameter passed to the gen_pos workflow option.}

\item{response}{Optional response-column name used to exclude prediction-only rows from weighted residual-variance summaries.}

\item{gen_name}{Optional genotype identifier column used in the observation-level residual-variance table.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{asreml_herit_varCov_new}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
