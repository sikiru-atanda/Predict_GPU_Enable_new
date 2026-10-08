\name{feature_select_top_k}
\alias{feature_select_top_k}
\title{Select top-k predictors from feature score metadata}
\usage{
feature_select_top_k(predictor_data, feature_score_metadata, trait, k)
}
\arguments{
\item{predictor_data}{Matrix or data frame to subset.}
\item{feature_score_metadata}{Output from \code{feature_score_predictors()}.}
\item{trait}{Trait name to use.}
\item{k}{Number of predictors to keep, or \code{"all"}.}
}
\value{
Predictor data subset to the top-k ranked columns.
}
\description{
Applies a feature score metadata table to a predictor matrix without changing row order.
}
