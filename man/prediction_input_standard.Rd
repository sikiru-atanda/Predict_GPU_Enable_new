\name{prediction_input_standard}
\alias{prediction_input_standard}
\alias{validate_prediction_input}
\title{Inspect and validate PredictProR input contracts}
\usage{
prediction_input_standard(
  task = c("all", "general", "hybrid", "multi_trait", "multi_environment")
)

validate_prediction_input(
  task = c("auto", "general", "hybrid", "multi_trait", "multi_environment"),
  ...
)
}
\arguments{
\item{task}{Workflow contract to inspect or validate. The automatic validator
infers hybrid, multi-trait, MET, or general mode from the named arguments.}
\item{...}{Named arguments accepted by the selected workflow-specific
validator.}
}
\value{
\code{prediction_input_standard()} returns one contract or a named list of all
contracts. \code{validate_prediction_input()} invisibly returns the validation
result, resolved task, and validator name.
}
\description{
Provides one public interface over the general, hybrid, multi-trait, and
multi-environment input standards used by \code{model_execute()}.
}
