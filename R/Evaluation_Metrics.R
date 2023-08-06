
#' Title
#'
#' @param y_observed
#' @param eval_metrics
#' @param ...
#' @param y_predicted
#'
#' @return
#' @export
#'
#' @examples
 evaluation_metrics <- function(y_observed=NULL,
                                y_predicted=NULL,
                                eval_metrics = c("Accuracy",
                                                 "Mean_Squared_Error",
                                                 "Bias",
                                                 "Root_Mean_Squared_Error",
                                                 "Relative_Squared_Error",
                                                 #"Absolute_Error",
                                                 "Mean_Absolute_Error",
                                                 #"Absolute_Percent_Error",
                                                 "Mean_Absolute_Percent_Error"),
                                                   ... ) {

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
      Root_Square <- function(y_observed, y_predicted) {
        return((y_observed - y_predicted) ^ 2)
      }

      ## Sum of Squared Errors
      ## computes the sum of the squared differences between y_observed and y_predicted..
      Sum_Square_Error <- function(y_observed, y_predicted) {
        return(sum(Root_Square(y_observed, y_predicted)))
      }


      ## Mean Squared Error
      # computes the average squared difference between y_observed and y_predicted.

      Mean_Squared_Error <- function(y_observed, y_predicted) {
        return(mean(Root_Square(y_observed, y_predicted)))
      }

      ## Root Mean Squared Error
      ## computes the root mean squared error between y_observed and y_predicted.

      Root_Mean_Squared_Error <- function(y_observed, y_predicted) {
        ## sqrt(mean((y_observed - y_predicted)^2))
        return(sqrt(Mean_Squared_Error(y_observed, y_predicted)))
      }

      ### Relative Squared Error
      Relative_Squared_Error <- function (actual, predicted) {
        return(Sum_Square_Error(actual, predicted) / Sum_Square_Error(actual, mean(actual)))
      }

      ## Absolute Error

      ## computes the elementwise absolute difference between y_observed and y_predicted

      Absolute_Error <- function(y_observed, y_predicted) {
        return(abs(y_observed - y_predicted))
      }

      ## Mean Absolute Error

      Mean_Absolute_Error <- function(y_observed, y_predicted) {
        return(mean(Absolute_Error(y_observed, y_predicted)))
      }
      Median_Absolute_Error <- function(y_observed, y_predicted) {
        return(stats::median(Absolute_Error(y_observed, y_predicted)))
      }

      ## Absolute Percent Error
      Absolute_Percent_Error <- function(y_observed, y_predicted) {
        return(Absolute_Error(y_observed, y_predicted) / abs(y_observed))
      }

      ## Mean Absolute Percent Error
      Mean_Absolute_Percent_Error <- function(y_observed, y_predicted) {
        return(mean(Absolute_Percent_Error(y_observed, y_predicted)))
      }

      if (sum(eval_metrics%in%"Accuracy")){
      res <- Accuracy(y_observed, y_predicted)
      }

      if (sum(eval_metrics%in%"Bias")){
        res <- Bias(y_observed, y_predicted)
      }

      if (sum(eval_metrics%in%"Mean_Squared_Error")){
        res <-  Mean_Squared_Error(y_observed, y_predicted)
      }

      if (sum(eval_metrics%in%"Root_Mean_Squared_Error")){
        res <-  Root_Mean_Squared_Error(y_observed, y_predicted)
      }

      if (sum(eval_metrics%in%"Relative_Squared_Error")){
        res <-  Relative_Squared_Error(y_observed, y_predicted)
      }

      # if (sum(eval_metrics%in%"Absolute_Error")){
      #   res = Absolute_Error(y_observed, y_predicted)
      # }

      if (sum(eval_metrics%in%"Mean_Absolute_Error")){
        res <-  Mean_Absolute_Error(y_observed, y_predicted)
      }

      if (sum(eval_metrics%in%"Median_Absolute_Error")){
        res <-  Median_Absolute_Error(y_observed, y_predicted)
      }


      # if (sum(eval_metrics%in%"Absolute_Percent_Error")){
      #   res <-  Absolute_Percent_Error(y_observed, y_predicted)
      # }


      if (sum(eval_metrics%in%"Mean_Absolute_Percent_Error")){
        res <-  Mean_Absolute_Percent_Error(y_observed, y_predicted)
      }




      return(res)
  }
