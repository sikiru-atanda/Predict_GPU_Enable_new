#' Title
#'
#' @param boot_results
#' @param GID_names
#' @param CI_width_thresholds
#' @param predictions
#' @param standard_errors
#' @param prediction_error_var
#' @param genetic_var
#' @param confidence_level
#' @param model_for_CI_cal
#' @param composite_reliability_score
#' @param composite_reliability
#' @param composite_reliability_percentage
#' @param high_reliability_thres
#' @param low_reliability_thres
#' @param system_database
#'
#' @return
#' @export
#'
#' @examples
diagnostic_plot_true_prediction <- function(boot_results = NULL,
                                           GID_names = NULL,
                                           CI_width_thresholds = c(0.33, 0.66),
                                           predictions = NULL,
                                           standard_errors = NULL,
                                           prediction_error_var = NULL,
                                           genetic_var = NULL,
                                           confidence_level = 0.95,
                                           model_for_CI_cal = "ML",
                                           composite_reliability_score = NULL,
                                           composite_reliability = NULL,
                                           composite_reliability_percentage = NULL,
                                           #threshold = NULL,
                                           high_reliability_thres = 0.7,
                                           low_reliability_thres = 0.4,
                                           system_database = FALSE){

  msg <- sprintf("==================================================\n")

  if (is.null(composite_reliability) || any(is.na(composite_reliability)) ||
      is.null(composite_reliability_score) || any(is.na(composite_reliability_score))) {

    composite_reliability <- NULL
    composite_reliability_score <- NULL
    p111 <- NULL
  }

  if(model_for_CI_cal == "Bayes" | model_for_CI_cal == "GBLUP") {
  if(is.null(prediction_error_var) & !is.null(standard_errors)){
  prediction_error_var <- (standard_errors)^2

  } else {
    if(is.null(prediction_error_var) & is.null(standard_errors)){

      stop(print(paste(msg,'Provide vector of prediction error or standard error of the predictions.')), call. = FALSE)


    }
  }
    if(is.null(genetic_var))  stop(print(paste(msg,'Provide genetic variance.')), call. = FALSE)
}




  if(!is.null(predictions) & model_for_CI_cal == "Bayes" | model_for_CI_cal == "GBLUP"){
    res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                  predictions = predictions,
                                                  standard_errors = standard_errors,
                                                  confidence_level = confidence_level,
                                                  model_for_CI_cal = model_for_CI_cal,
                                                  boot_results = boot_results)

    predictions <- res_CI$predictions
    genetic_var <- var(predictions)
    prediction_error_var <- res_CI$prediction_error_var
  } else {

    res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                  predictions = NULL,
                                                  standard_errors = NULL,
                                                  confidence_level = confidence_level,
                                                  model_for_CI_cal = model_for_CI_cal,
                                                  boot_results = boot_results)

    predictions <- res_CI$predictions
    genetic_var <- var(predictions)
    prediction_error_var <- res_CI$prediction_error_var
    #stop(print(paste(msg,'Provide vector of the predicted values.')), call. = FALSE)

  }

  res <-  reliability_thresholds(prediction_error_var = prediction_error_var,
                                 genetic_var =  genetic_var,
                                 high_reliability_thres = high_reliability_thres,
                                 low_reliability_thres = low_reliability_thres)



  # if (res$reliability_percentage > 80) {
  #   message(insight::print_color(paste(msg,paste("Decision: Most predictions are reliable.\n")), "blue"))
  #   message("Decision: Most predictions are reliable.\n")
  # } else if (res$reliability_percentage > 50) {
  #   message(insight::print_color(paste(msg,paste("Decision: More than half of the predictions are reliable, but there is some uncertainty.\n")), "blue"))
  #
  # } else {
  #   if (res$reliability_percentage < 50) {
  #   warning(insight::print_color(paste(msg,paste("Decision: Less than half of the predictions are reliable. Consider reviewing the model or data.\n")), "red"))
  #    }
  # }


  predicted_values = data.frame(Predicted_value = predictions,
                                Reliability = res$reliability,
                                Remarks = res$remarks)

  unique_remakers = c("Reliable", "Acceptable", "Unreliable")

  # Convert the 'Remarks' column to a factor with specified levels
  predicted_values$Remarks <- factor(predicted_values$Remarks, levels = c("Reliable", "Acceptable", "Unreliable"))

  # Classify predictions based on this threshold
  #predicted_values$Remarks <- predicted_values$Reliability > acceptable_se_threshold

  if(!is.null(composite_reliability) & !is.null(composite_reliability_score)){
    # Count occurrences the three classes
    count_high <- sum(composite_reliability == 'high')
    count_moderate <- sum(composite_reliability == 'moderate')
    count_low <- sum(composite_reliability == 'low')

    # Total number of elements
    total <- length(composite_reliability)

    # Calculate proportions
    proportion_high_reliability <- count_high / total
    proportion_moderate_reliability <- count_moderate / total
    proportion_low_reliability <- count_low / total

    data <- data.frame(Predicted_value = predictions,
                       Reliability = composite_reliability_score,
                       Remarks = composite_reliability)


    p111 <- ggplot2::ggplot(data, ggplot2::aes(x = Reliability, y =Predicted_value , color = Remarks)) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::labs(title = "Prediction Confidence Analysis",
                          subtitle = sprintf("Proportion of composite reliability: high, moderate, and low: %.2f%%, %.2f%%, and %.2f%%",
                                             proportion_high_reliability * 100, proportion_moderate_reliability * 100, proportion_low_reliability * 100),
                          # subtitle = sprintf("Proportion of high reliability: %.2f%%, Proportion of medium reliability: %.2f%%, Proportion of low reliability: %.2f%%",
                          #                    res$proportion_high_reliability*100, res$proportion_medium_reliability*100, res$proportion_low_reliability*100),
                          x = "Composite Reliability",
                          y = "Predicted Values") +
            ggplot2::scale_color_manual(values = setNames(
              c("green", "lightblue",  "red"),
              unique_remakers
            )) +
            #scale_color_manual(values = c("red", "green"), labels = c("Unreliable", "Reliable")) +
            ggplot2::theme_minimal()+
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                           plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"))



  }

### Plot of Reliability VS predicted value

  p11 <- ggplot2::ggplot(predicted_values, ggplot2::aes(x = Reliability, y =Predicted_value , color = Remarks)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::labs(title = "Prediction Confidence Analysis",
                      subtitle = sprintf("Proportion of reliability: high, medium, and low: %.2f%%, %.2f%%, and %.2f%%",
                                         res$proportion_high_reliability * 100, res$proportion_medium_reliability * 100, res$proportion_low_reliability * 100),
                      # subtitle = sprintf("Proportion of high reliability: %.2f%%, Proportion of medium reliability: %.2f%%, Proportion of low reliability: %.2f%%",
                      #                    res$proportion_high_reliability*100, res$proportion_medium_reliability*100, res$proportion_low_reliability*100),
                      x = "Reliability",
                      y = "Predicted Values") +
        ggplot2::scale_color_manual(values = setNames(
          c("green", "lightblue",  "red"),
          unique_remakers
        )) +
        #scale_color_manual(values = c("red", "green"), labels = c("Unreliable", "Reliable")) +
        ggplot2::theme_minimal()+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"))


  testData <- data.frame(predictions = predictions,
                         lower_bound = res_CI$lower_bound,
                         upper_bound = res_CI$upper_bound,
                         reliability_remarks = res_CI$reliability_remarks,
                         relative_interval_width = res_CI$relative_interval_width,
                         Standard_Errors = res_CI$standard_errors,
                         threshold = as.double(res_CI$high_threshold))

  unique_remakers_CI <-  c("high", "moderate", "low")

  # Remap the values in the column to match reliability instance
  testData <- testData |>
    dplyr::mutate(reliability_remarks = dplyr::case_when(
      reliability_remarks == "low" ~ "high",
      reliability_remarks == "high" ~ "low",
      TRUE ~ reliability_remarks
    ))


  testData$reliability_remarks <- factor(testData$reliability_remarks, levels = unique_remakers_CI)



  p12 <- ggplot2::ggplot(testData, ggplot2::aes(x = predictions, y = predictions,
                                               ymin = lower_bound, ymax = upper_bound)) +
        ggplot2::geom_errorbar() +
        #geom_point() +
        ggplot2::geom_point(ggplot2::aes(color = reliability_remarks)) +
        ggplot2::scale_color_manual(values = stats::setNames(
          c("green", "lightblue",  "red"),
          unique_remakers_CI
        )) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Prediction Intervals and Reliability",
                      x = "Predicted Values",
                      y = "Predicted Values",
                      color = "Reliability")+ ## This change the title of the legend
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
              plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue")
              )


  # Scatter plot of Predictions vs. RIW

  # p13 <- ggplot2::ggplot(testData, ggplot2::aes(x = predictions, y = relative_interval_width)) +
  #   ggplot2::geom_point(ggplot2::aes(color = reliability_remarks), size = 2) +
  #   ggplot2::scale_color_manual(values = stats::setNames(
  #     c("green", "lightblue",  "red"),
  #     unique_remakers_CI
  #   )) +
  #   ggplot2::geom_hline(yintercept = testData$threshold, linetype = "dashed", color = "red") +
  #   ggplot2::labs(title = "Predictions vs. Relative Interval Width", x = "Predictions", y = "Relative Interval Width") +
  #   ggplot2::theme_minimal()


  # # Scatter plot of Standard Errors vs. RIW
  # p14 <- ggplot2::ggplot(testData, ggplot2::aes(x = Standard_Errors, y = relative_interval_width)) +
  #   ggplot2::geom_point(ggplot2::aes(color = reliability_remarks), size = 2) +
  #   ggplot2::scale_color_manual(values = stats::setNames(
  #     c("green", "lightblue",  "red"),
  #     unique_remakers_CI
  #   )) +
  #   ggplot2::geom_hline(yintercept = testData$threshold, linetype = "dashed", color = "red") +
  #   ggplot2::labs(title = "Standard Errors vs. Relative Interval Width", x = "Standard Errors", y = "Relative Interval Width") +
  #   ggplot2::theme_minimal()


  if(isTRUE(system_database)){
  return(list(predicted_vs_reliability = p11,
              prediction_inter_vs_reliability = p12,
              #predicted_vs_RIW = p13,
              predicted_vs_composite_reliability = p111
              #Standard_Errors_vs_RIW = p14

  ))

  } else {
    if(is.null(p111)){

    combined_plot <- gridExtra::grid.arrange(p11, p12, ncol = 2)

    }else{

      #combined_plot <- gridExtra::grid.arrange(p11, p111, p12, ncol = 2)
      combined_plot <- gridExtra::grid.arrange(p11, p12, ncol = 2)

    }

    # name_plot <- paste(paste("Cross_validation_diagonistic_plots", model, sep = "_"), trait, sep = "_")
    #
    # plot_save <- ggplot2::ggsave(paste(name_plot,"pdf", sep = "."), plot = combined_plot,
    #                              width = 15, height = 8, units = "in", dpi = 300)
    return(combined_plot)
    # Save the combined plot
    # return(ggplot2::ggsave("True_prediction_diagonistic_plots.pdf", plot = combined_plot,
    #                        width = 15, height = 8, units = "in", dpi = 300))
  }


}
