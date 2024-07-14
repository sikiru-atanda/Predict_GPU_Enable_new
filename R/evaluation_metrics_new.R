#' Evaluate Prediction Performance Metrics
#'
#' This function computes various evaluation metrics to assess the performance of prediction models.
#'
#' @param y_observed A numeric vector of observed values.
#' @param y_predicted A numeric vector of predicted values, corresponding to `y_observed`.
#' @param eval_metrics A character vector specifying which evaluation metrics to compute. Possible values include "Accuracy", "Mean_Squared_Error", "Bias", "Root_Mean_Squared_Error", "Relative_Squared_Error", "Mean_Absolute_Error", and "Mean_Absolute_Percent_Error".
#' @return The function returns a numeric value representing the computed evaluation metric specified by `eval_metrics`.
#' @examples
#' \dontrun{
#' observed <- c(1, 2, 3, 4, 5)
#' predicted <- c(1.1, 1.9, 3.1, 3.9, 4.8)
#' mse <- evaluation_metrics(y_observed = observed,
#'                           y_predicted = predicted,
#'                           eval_metrics = "Mean_Squared_Error")
#' print(mse)
#'
#' mae <- evaluation_metrics(y_observed = observed,
#'                           y_predicted = predicted,
#'                           eval_metrics = "Mean_Absolute_Error")
#' print(mae)
#'}
#' @export
#'
evaluation_metrics <- function(y_observed = NULL,
                               y_predicted = NULL,
                               eval_metrics = NULL) {

  msg <- "\n==================================================\n"

  if (is.null(y_observed) || is.null(y_predicted)) {
    stop("Both observed and predicted values must be provided.")
  }
  if (length(y_observed) != length(y_predicted)) {
    stop("Observed and predicted values must have the same length.")
  }

  eval_metrics_available <- c("accuracy", "mean_squared_error", "bias",
                              "root_mean_squared_error", "relative_squared_error",
                              "mean_absolute_error", "mean_absolute_percent_error", "kendalls_tau")

  metric_functions <- list(
    accuracy = function(y, y_hat) cor(y, y_hat),
    mean_squared_error = function(y, y_hat) mean((y - y_hat)^2),
    bias = function(y, y_hat) mean(y - y_hat),
    root_mean_squared_error = function(y, y_hat) sqrt(mean((y - y_hat)^2)),
    relative_squared_error = function(y, y_hat) sum((y - y_hat)^2) / sum((y - mean(y))^2),
    mean_absolute_error = function(y, y_hat) mean(abs(y - y_hat)),
    mean_absolute_percent_error = function(y, y_hat) mean(abs((y - y_hat) / y)) * 100,
    kendalls_tau = kendalls_tau
  )

  # Initialize an empty list to store results
  #results <- list()

  for (metric in eval_metrics) {
    if (!metric %in% names(metric_functions)) {
      stop(msg, "Unknown evaluation metric. Choose from: ", paste(eval_metrics_available, collapse = ", "), call. = FALSE)

    }
    func_metric <- metric_functions[[metric]]
    #results[[metric]] <- func(y_observed, y_predicted)
    results <- func_metric(y_observed, y_predicted)
  }

  return(results)
}

