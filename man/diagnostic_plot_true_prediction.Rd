% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{diagnostic_plot_true_prediction}
\alias{diagnostic_plot_true_prediction}
\title{diagnostic_plot_true_prediction}
\usage{
diagnostic_plot_true_prediction(boot_results = NULL, GID_names = NULL, Predicted_value_for_CI = NULL, 
    mod = NULL, CI_width_thresholds = c(0.33, 0.66), predictions = NULL, 
    standard_errors = NULL, prediction_error_var = NULL, genetic_var = NULL, 
    reference_variance = NULL, lower_bound = NULL, upper_bound = NULL, 
    observed_y = NULL, train_test_label = NULL, confidence_level = 0.95, 
    risk_calibration = NULL, require_heldout_calibration = TRUE,
    model_for_CI_cal = "ML", response_family = "gaussian", composite_reliability_score = NULL, 
    composite_reliability = NULL, composite_reliability_percentage = NULL, 
    high_reliability_thres = 0.7, low_reliability_thres = 0.4, 
    system_database = FALSE) 
}
\arguments{
\item{boot_results}{Parameter passed to the boot_results workflow option.}

\item{GID_names}{Parameter passed to the GID_names workflow option.}

\item{Predicted_value_for_CI}{Parameter passed to the Predicted_value_for_CI workflow option.}

\item{mod}{Fitted model object.}

\item{CI_width_thresholds}{Parameter passed to the CI_width_thresholds workflow option.}

\item{predictions}{Parameter passed to the predictions workflow option.}

\item{standard_errors}{Parameter passed to the standard_errors workflow option.}

\item{prediction_error_var}{Parameter passed to the prediction_error_var workflow option.}

\item{genetic_var}{Parameter passed to the genetic_var workflow option.}

\item{reference_variance}{Parameter passed to the reference_variance workflow option.}

\item{lower_bound}{Parameter passed to the lower_bound workflow option.}

\item{upper_bound}{Parameter passed to the upper_bound workflow option.}

\item{observed_y}{Parameter passed to the observed_y workflow option.}

\item{train_test_label}{Parameter passed to the train_test_label workflow option.}

\item{confidence_level}{Parameter passed to the confidence_level workflow option.}

\item{risk_calibration}{Held-out or cross-fitted calibration information for Gaussian ML/DL prediction-error and interval diagnostics.}

\item{require_heldout_calibration}{Logical. If \code{TRUE}, Gaussian ML/DL diagnostics fail rather than substituting bootstrap instability or in-sample standard errors when held-out calibration is unavailable.}

\item{model_for_CI_cal}{Parameter passed to the model_for_CI_cal workflow option.}

\item{response_family}{Response family used by the model workflow.}

\item{composite_reliability_score}{Parameter passed to the composite_reliability_score workflow option.}

\item{composite_reliability}{Parameter passed to the composite_reliability workflow option.}

\item{composite_reliability_percentage}{Parameter passed to the composite_reliability_percentage workflow option.}

\item{high_reliability_thres}{Parameter passed to the high_reliability_thres workflow option.}

\item{low_reliability_thres}{Parameter passed to the low_reliability_thres workflow option.}

\item{system_database}{Logical indicating whether system database style outputs should be used.}

}
\description{
Documentation for \code{diagnostic_plot_true_prediction}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
