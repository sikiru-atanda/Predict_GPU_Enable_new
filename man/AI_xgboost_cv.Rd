% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{AI_xgboost_cv}
\alias{AI_xgboost_cv}
\title{AI_xgboost_cv}
\usage{
AI_xgboost_cv(y, omics, tst, response_family = "gaussian", xgb_booster = "gbtree", 
    xgb_rate_drop = 0.1, xgb_skip_drop = 0.5, xgb_objective = "reg:squarederror", 
    xgb_sample_type = "uniform", xgb_normalize_type = "tree", 
    eta = 0.01, nrounds = 100, max_depth = 6, scaling = FALSE, 
    centering = TRUE, omic_count, xgb_gamma = 0.01, min_child_weight = 1, 
    subsample = 0.7, colsample_bytree = 0.7, xgb_alpha = 0.001, 
    xgb_lambda = 1, early_stop_for_iteration_xgb = FALSE, xgb_nthread = 1L, 
    random_state = NULL) 
}
\arguments{
\item{y}{Parameter passed to the y workflow option.}

\item{omics}{Parameter passed to the omics workflow option.}

\item{tst}{Parameter passed to the tst workflow option.}

\item{response_family}{Response family used by the model workflow.}

\item{xgb_booster}{Parameter passed to the xgb_booster workflow option.}

\item{xgb_rate_drop}{Parameter passed to the xgb_rate_drop workflow option.}

\item{xgb_skip_drop}{Parameter passed to the xgb_skip_drop workflow option.}

\item{xgb_objective}{Parameter passed to the xgb_objective workflow option.}

\item{xgb_sample_type}{Parameter passed to the xgb_sample_type workflow option.}

\item{xgb_normalize_type}{Parameter passed to the xgb_normalize_type workflow option.}

\item{eta}{Parameter passed to the eta workflow option.}

\item{nrounds}{Parameter passed to the nrounds workflow option.}

\item{max_depth}{Parameter passed to the max_depth workflow option.}

\item{scaling}{Logical indicating whether predictors should be scaled.}

\item{centering}{Logical indicating whether predictors should be centered.}

\item{omic_count}{Parameter passed to the omic_count workflow option.}

\item{xgb_gamma}{Parameter passed to the xgb_gamma workflow option.}

\item{min_child_weight}{Parameter passed to the min_child_weight workflow option.}

\item{subsample}{Parameter passed to the subsample workflow option.}

\item{colsample_bytree}{Parameter passed to the colsample_bytree workflow option.}

\item{xgb_alpha}{Parameter passed to the xgb_alpha workflow option.}

\item{xgb_lambda}{Parameter passed to the xgb_lambda workflow option.}

\item{early_stop_for_iteration_xgb}{Parameter passed to the early_stop_for_iteration_xgb workflow option.}

\item{xgb_nthread}{Parameter passed to the xgb_nthread workflow option.}

\item{random_state}{Random seed used for reproducible splits.}

}
\description{
Documentation for \code{AI_xgboost_cv}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
