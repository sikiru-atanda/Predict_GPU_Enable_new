% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{AI_randomforest_cv}
\alias{AI_randomforest_cv}
\title{AI_randomforest_cv}
\usage{
AI_randomforest_cv(y, omics, tst, response_family = "gaussian", ntree = 500, 
    mtry = NULL, maxnodes = NULL, nodesize = NULL, rf_n_jobs = 1L, 
    scaling = FALSE, centering = TRUE, omic_count) 
}
\arguments{
\item{y}{Parameter passed to the y workflow option.}

\item{omics}{Parameter passed to the omics workflow option.}

\item{tst}{Parameter passed to the tst workflow option.}

\item{response_family}{Response family used by the model workflow.}

\item{ntree}{Parameter passed to the ntree workflow option.}

\item{mtry}{Parameter passed to the mtry workflow option.}

\item{maxnodes}{Parameter passed to the maxnodes workflow option.}

\item{nodesize}{Parameter passed to the nodesize workflow option.}

\item{rf_n_jobs}{Number of internal Python RandomForest jobs. Keep this at 1 when PredictProR is already running models in parallel.}

\item{scaling}{Logical indicating whether predictors should be scaled.}

\item{centering}{Logical indicating whether predictors should be centered.}

\item{omic_count}{Parameter passed to the omic_count workflow option.}

}
\description{
Documentation for \code{AI_randomforest_cv}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
