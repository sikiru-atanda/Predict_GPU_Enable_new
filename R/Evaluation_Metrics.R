
#' Title
#'
#' @param y_observed
#' @param y_predicted
#' @param Metrics
#'
#' @return
#' @export
#'
#' @examples
 Evaluation_Metrics <- function(y_observed, y_predicted,
                               Metrics = c("Accuracy",
                                           "Bias",
                                           "Percent_Bias",
                                           "Sum_Square_Error",
                                           "Mean_Squared_Error",


                                           )) {

   Accuracy <-  function(y_observed, y_predicted) {

      return(stats::cor(y_observed,y_predicted))
    }

   ### https://rdrr.io/cran/Metrics/src/R/regression.R
   # Bias
   # computes the average amount by which fitted value is greater than y_predicted.
      Bias <- function(y_observed, y_predicted) {
        return(mean(y_observed - y_predicted))
      }


      Percent_Bias <- function(y_observed, y_predicted) {
        return(mean((y_observed - y_predicted) / abs(y_observed)))
      }

      ## Root Square
      ## computes the element wise squared difference between y_observed and y_predicted.
      rs <- function(y_observed, y_predicted) {
        return((y_observed - y_predicted) ^ 2)
      }

      ## Sum of Squared Errors
      ## computes the sum of the squared differences between y_observed and y_predicted..
      Sum_Square_Error <- function(y_observed, y_predicted) {
        return(sum(rs(y_observed, y_predicted)))
      }


      ## Mean Squared Error
      # computes the average squared difference between y_observed and y_predicted.

      Mean_Squared_Error <- function(y_observed, y_predicted) {
        return(mean(rs(y_observed, y_predicted)))
      }

      # ## Root Mean Squared Error
      # ## computes the root mean squared error between y_observed and y_predicted.
      #
      # rmse <- function(y_observed, y_predicted) {
      #   return(sqrt(mse(y_observed, y_predicted)))
      # }
      #
      # ### Relative Squared Error
      # rse <- function (actual, predicted) {
      #   return(sse(actual, predicted) / sse(actual, mean(actual)))
      # }
      #
      # ## Absolute Error
      #
      # ## computes the elementwise absolute difference between y_observed and y_predicted
      #
      # ae <- function(y_observed, y_predicted) {
      #   return(abs(y_observed - y_predicted))
      # }
      #
      # ## Mean Absolute Error
      #
      # mae <- function(y_observed, y_predicted) {
      #   return(mean(ae(y_observed, y_predicted)))
      # }
      # mdae <- function(y_observed, y_predicted) {
      #   return(stats::median(ae(y_observed, y_predicted)))
      # }
      #
      # ## Absolute Percent Error
      # ape <- function(y_observed, y_predicted) {
      #   return(ae(y_observed, y_predicted) / abs(y_observed))
      # }
      #
      # ## Mean Absolute Percent Error
      # mape <- function(y_observed, y_predicted) {
      #   return(mean(ape(y_observed, y_predicted)))
      # }

      if (Metrics=="Accuracy"){
      res <- Accuracy(y_observed, y_predicted)
      }

      if (Metrics=="Bias"){
        res <- Bias(y_observed, y_predicted)
      }

      if (Metrics=="Mean_Squared_Error"){
        res = Mean_Squared_Error(y_observed, y_predicted)
      }
      return(res)
  }
