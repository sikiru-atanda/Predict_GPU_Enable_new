\name{feature_score_predictors}
\alias{feature_score_predictors}
\title{Score and rank predictors for dynamic top-k model fitting}
\usage{
feature_score_predictors(
  predictor_data,
  pheno_data,
  response,
  gen_name,
  scoring_model = c("Ridge_Regression", "BayesB", "RandomForest"),
  seed = 123L,
  source_block = NULL,
  ridge_lambda = 1,
  bayes_nIter = 1500L,
  bayes_burnIn = 500L,
  bayes_thin = 5L,
  ntree = 500L,
  mtry = NULL,
  nodesize = NULL,
  rf_n_jobs = 1L,
  response_family = "auto",
  selection_mode = "fixed",
  replication = NA_integer_,
  fold = NA_integer_
)
}
\arguments{
\item{predictor_data}{Numeric predictor matrix or data frame with individuals in rows and predictors in columns.}
\item{pheno_data}{Phenotype data containing \code{gen_name} and response columns.}
\item{response}{Character vector of response variables to score.}
\item{gen_name}{Column name in \code{pheno_data} containing individual IDs.}
\item{scoring_model}{One of \code{"Ridge_Regression"}, \code{"BayesB"}, or \code{"RandomForest"}.}
\item{seed}{Integer seed used by stochastic scoring backends.}
\item{source_block}{Optional predictor source labels, either named by predictor
or a single recycled label.}
\item{ridge_lambda}{Positive ridge penalty for ridge scoring.}
\item{bayes_nIter}{Iterations for BayesB scoring.}
\item{bayes_burnIn}{Burn-in for BayesB scoring.}
\item{bayes_thin}{Thinning for BayesB scoring.}
\item{ntree}{Number of trees for Python RandomForest scoring.}
\item{mtry}{Number of predictors sampled at each split.}
\item{nodesize}{Minimum node size.}
\item{rf_n_jobs}{Internal Python RandomForest jobs. Keep at 1 when the caller is already parallelizing over outer tasks.}
\item{response_family}{Response family. \code{"auto"} resolves each trait.
Ridge and BayesB scoring require Gaussian responses; use RandomForest for
binary, ordinal, or multiclass responses.}
\item{selection_mode}{Provenance label, normally \code{"fixed"} or
\code{"fold_internal"}.}
\item{replication, fold}{Optional resampling identifiers recorded in metadata.}
}
\value{
A \code{predictpror_feature_scores} data frame with trait, source, score, rank,
response family, selection mode, resampling identifiers, training count, and a
deterministic training-ID hash. Scores are predictive rankings, not genetic
variance components or causal effects.
}
\description{
Creates a reusable trait-specific predictor metadata table containing ranks,
normalized scores, and training-set provenance.
}
