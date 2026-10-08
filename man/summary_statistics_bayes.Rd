\name{summary_statistics_bayes}
\alias{summary_statistics_bayes}
\title{Summarize Bayesian Prediction Outputs}
\usage{
summary_statistics_bayes(mod = NULL, eval_metrics = NULL, GS_model = NULL, model_result = NULL, 
    gen_name = NULL, response_family = "gaussian", CI_width_thresholds = c(0.33, 
        0.66), confidence_level = 0.95, high_reliability_thres = 0.7, 
    low_reliability_thres = 0.4, system_database = FALSE, heter_groups = NULL, 
    ...) 
}
\arguments{
\item{mod}{Fitted model object.}

\item{eval_metrics}{Evaluation metrics to compute.}

\item{GS_model}{Genomic selection model identifier.}

\item{model_result}{Parameter passed to the model_result workflow option.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{response_family}{Response family used by the model workflow.}

\item{CI_width_thresholds}{Parameter passed to the CI_width_thresholds workflow option.}

\item{confidence_level}{Parameter passed to the confidence_level workflow option.}

\item{high_reliability_thres}{Parameter passed to the high_reliability_thres workflow option.}

\item{low_reliability_thres}{Parameter passed to the low_reliability_thres workflow option.}

\item{system_database}{Logical indicating whether system database style outputs should be used.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Returns a compact summary table for Bayesian model fits, including phenotype
range, residual variance, train/test counts, and user-facing definitions for
target-specific uncertainty metrics when Bayesian prediction outputs are
available.
}
\details{
For Bayesian prediction tables in this package, \code{Genetic_variance},
\code{PEV}, and \code{Reliability} are reported for the exact prediction target
being shown. In single-environment analyses this is the individual genotype
target. In multi-environment analyses this can be either the
environment-specific genotype target or the across-environment genotype-average
target.

A near-zero reliability, including for \code{BayesB}, does not by itself imply
the predictions are wrong. It means that the posterior target-specific
prediction error variance is close to or larger than the corresponding
target-specific genetic variance under the fitted model and data.
}
\value{
A list containing a \code{summary_statistics} data frame.
}
