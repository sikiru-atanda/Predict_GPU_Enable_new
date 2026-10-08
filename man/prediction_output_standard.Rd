\name{prediction_output_standard}
\alias{prediction_output_standard}
\alias{validate_prediction_output}
\title{Inspect and validate PredictProR output contracts}
\usage{
prediction_output_standard(
  response_family = c("gaussian", "classification", "auto"),
  task = c("single_trait", "multi_trait", "multi_environment",
    "multi_trait_multi_environment", "hybrid", "hybrid_multi_trait",
    "hybrid_multi_environment", "hybrid_multi_trait_multi_environment"),
  gen_name = "GID",
  heter_groups = "Env"
)

validate_prediction_output(
  x,
  response_family = c("auto", "gaussian", "classification"),
  gen_name = NULL,
  heter_groups = NULL,
  tolerance = 1e-08,
  strict = TRUE
)
}
\arguments{
\item{response_family}{Gaussian, classification, or automatically detected
response contract.}
\item{task}{Prediction layout, including optional trait, environment, and
hybrid identifiers.}
\item{gen_name}{Genotype identifier column name.}
\item{heter_groups}{Environment identifier column name.}
\item{x}{A canonical model result or the indexed result returned by
\code{model_execute()}.}
\item{tolerance}{Numerical tolerance for uncertainty and bound checks.}
\item{strict}{Stop when a violation is found. When false, return all findings.}
}
\value{
The standard function returns required fields and ordered table columns. The
validator invisibly returns one audit row per model and trait.
}
\description{
Defines and checks the common model-result contract across Bayesian, ASReml,
Gaussian-process, classical machine-learning, deep-learning, multi-trait, MET,
and hybrid prediction paths.
}
