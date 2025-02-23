
# Function to calculate RIW for predictions for Bayes model using using the standard error and predictions value
calculate_intervals_bayes_GBLUP <- function(predictions,
                                            standard_errors,
                                            confidence_level = 0.95
                                            #threshold = NULL

                                      ) {
  z_value <- qnorm(1 - (1 - confidence_level) / 2)
  lower_bound <- predictions - z_value * standard_errors
  upper_bound <- predictions + z_value * standard_errors
  interval_width <- upper_bound - lower_bound
  #relative_interval_width <- interval_width / predictions


  return(list(lower_bound = lower_bound,
              upper_bound = upper_bound,
              interval_width = interval_width))

}


# Reliability categorization function
reliability_thresholds_MPIW_from_CI <- function(boot_results = NULL,
                                                Predicted_value_for_CI = NULL,
                                                cv_results = NULL,
                                                replication = NULL,
                                                mod = NULL,
                                                CI_width_thresholds = c(0.33, 0.66),
                                                interval_width_high_threshold = NULL,
                                                interval_width_low_threshold = NULL,
                                                y_scaler = NULL,
                                                predictions = NULL,
                                                standard_errors = NULL,
                                                confidence_level = 0.95,
                                                model_for_CI_cal = "ML") {

  msg <- "\n==================================================\n"
if(model_for_CI_cal=="ML"){
  if(is.null(boot_results)) stop("Boostrapping results is required for machine learning models result diagonistic")
  #boot_results$t <- revert_scaling(boot_results$t, y_scaler)
  lower_bound <- apply(boot_results$t, 2, quantile, probs = 0.05)
  upper_bound <- apply(boot_results$t, 2, quantile, probs = 0.95)
  # lower_bound <- revert_scaling(lower_bound, y_scaler)
  # upper_bound <- revert_scaling(upper_bound, y_scaler)
  interval_width <- upper_bound - lower_bound
  predictions <- apply(boot_results$t, 2, mean)
  standard_errors <- apply(boot_results$t, 2, sd)
  prediction_error_var <- apply(boot_results$t, 2, var)


} else if (model_for_CI_cal=="Bayes"){

  Predicted_value_for_CI <- Predicted_value_for_CI + mod$model$mu
  lower_bound <- apply(Predicted_value_for_CI, 1, quantile, probs = 0.05)
  upper_bound <- apply(Predicted_value_for_CI, 1, quantile, probs = 0.95)
  interval_width <- upper_bound - lower_bound
  predictions <- apply(Predicted_value_for_CI, 1, mean)
  standard_errors <- apply(Predicted_value_for_CI, 1, sd)
  prediction_error_var <- apply(Predicted_value_for_CI, 1, var)

} else {

  if(model_for_CI_cal == "RKHS" | model_for_CI_cal == "GBLUP"){
    if(is.null(standard_errors)) stop("Provide the standard error of the prediction for results diagonistic of Bayes and GBLUP model.")
    if (is.null(predictions)) stop("Provide predicted value for results diagnoistic of Bayes and GBLUP model.")

    res <- calculate_intervals_bayes_GBLUP(predictions = predictions,
                                           standard_errors = standard_errors,
                                           confidence_level = confidence_level)

  lower_bound <-  res$lower_bound
  upper_bound <-  res$upper_bound
  interval_width <- upper_bound - lower_bound

  prediction_error_var = (standard_errors)^2

  }
###################################
  ## please holder
######################################
  if(!is.null(cv_results) & !is.null(replication)){

    res <- calculate_intervals_bayes_GBLUP(predictions = predictions,
                                           standard_errors = standard_errors,
                                           confidence_level = confidence_level)

    lower_bound <-  res$lower_bound
    upper_bound <-  res$upper_bound
    interval_width <- upper_bound - lower_bound

    prediction_error_var = (standard_errors)^2

  }

}

if((is.null(interval_width_high_threshold) & is.null(interval_width_low_threshold)) & !is.null(CI_width_thresholds)){
  high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1])
  low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2])
} else if ((is.null(interval_width_high_threshold) | is.null(interval_width_low_threshold)) & !is.null(CI_width_thresholds)){
  high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1])
  low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2])
} else{
  if ((is.null(interval_width_high_threshold) | is.null(interval_width_low_threshold)) & is.null(CI_width_thresholds)) {
    stop(print(paste(msg,'Both interval_width_high_threshold and interval_width_low_threshold must be either NULL or not NULL.')), call. = FALSE)

  }
}

  rel_remarks <- ifelse(interval_width < high_threshold, "low",
                        ifelse(interval_width < low_threshold, "moderate", "high"))

  noise_reliability <- ifelse(interval_width < high_threshold, 1,
                              ifelse(interval_width < low_threshold, 0.5, 0))

  avg_interval_width <- mean(interval_width)
  #proportion_high_reliability <- mean(interval_width < quantile(interval_width, thresholds[1]))

  # Calculate proportions for each reliability level
  proportion_low_reliability <- mean(interval_width >= high_threshold)
  proportion_medium_reliability <- mean(interval_width > high_threshold & interval_width < high_threshold)
  proportion_high_reliability <- mean(interval_width <= high_threshold)

  reliability_percentage <- (proportion_high_reliability + proportion_medium_reliability)*100

  relative_interval_width <- interval_width / predictions



  return(list(lower_bound = lower_bound,
              upper_bound = upper_bound,
              MPIW = avg_interval_width,
              Uncertainty = interval_width,
              relative_interval_width = relative_interval_width,
              reliability_remarks = rel_remarks,
              proportion_high_reliability = proportion_high_reliability,
              proportion_medium_reliability = proportion_medium_reliability,
              proportion_low_reliability = proportion_low_reliability,
              reliability_percentage = reliability_percentage,
              predictions = predictions,
              standard_errors = standard_errors,
              prediction_error_var = prediction_error_var,
              high_threshold = as.double(high_threshold)))

}


reliability_thresholds <- function(prediction_error_var = NULL,
                                   genetic_var = NULL,
                                   high_reliability_thres = 0.6,
                                   low_reliability_thres = 0.2) {

  msg <- "\n==================================================\n"
  high_reliability_thres = 0.6
  low_reliability_thres = 0.2
  if(is.null(genetic_var)){
    stop(print(paste(msg,'Provide genetic or additive variance to calculate reliability.')), call. = FALSE)
  }
  reliability <- 1-(prediction_error_var/(genetic_var))

  reliability <- pmax(0, pmin(1, reliability))
  # The default is to find the threshold for the top 10 percent and as such,
  # the 90th percentile was calculated.
  # which implies the value above which 10% of the reliability scores fall.
  # # This helps maintain a higher standard for the reliability of your predictions
  if(is.null(high_reliability_thres)){
    high_reliability_thres <- stats::quantile(reliability, 0.8)
  }

  if(is.null(low_reliability_thres)){
    low_reliability_thres <- stats::quantile(reliability, 0.1)
  }

  rel_remarks <- ifelse(reliability >= high_reliability_thres, "Reliable",
                        ifelse(reliability <= low_reliability_thres, "Unreliable", "Acceptable"))

  proportion_high_reliability <- mean(reliability >= high_reliability_thres)

  proportion_medium_reliability <- mean(reliability > low_reliability_thres & reliability <= high_reliability_thres)
  proportion_low_reliability <- mean(reliability <= low_reliability_thres)

  reliability_percentage <- (proportion_high_reliability + proportion_medium_reliability)*100

  return(list(reliability = reliability,
              remarks = rel_remarks,
              proportion_high_reliability = proportion_high_reliability,
              proportion_medium_reliability = proportion_medium_reliability,
              proportion_low_reliability = proportion_low_reliability,
              reliability_percentage = reliability_percentage))
}

########
# Composite Reliability: Use composite reliability metrics that combine multiple
# aspects (prediction interval, individual similarity, etc.)
# to get a comprehensive measure of prediction reliability.

composite_reliability_tst <- function(geno_trn = NULL,
                                      geno_tst = NULL,
                                      geno_tst_trn = NULL,
                                      names_tst = NULL,
                                      names_trn = NULL,
                                      n_components = 20,
                                      threshold = 100,
                                      target = "test_set",
                                      CI_width_thresholds = c(0.33, 0.66),
                                      interval_width_high_threshold = NULL,
                                      interval_width_low_threshold = NULL,
                                      interval_width = NULL,
                                      apply_pca = TRUE) {

  msg <- "\n==================================================\n"

  if (setequal(rownames(geno_trn), rownames(geno_tst))) {
    return(list(trustworthiness=NA,
                proportion_high_reliability = NA,
                proportion_medium_reliability = NA,
                proportion_low_reliability = NA,
                reliability_percentage = NA,
                reliability_score = NA))
  }

  if (is.null(geno_tst)) {
    return(list(trustworthiness=NA,
                proportion_high_reliability = NA,
                proportion_medium_reliability = NA,
                proportion_low_reliability = NA,
                reliability_percentage = NA,
                reliability_score = NA))
  }
  # Helper function to calculate thresholds based on IQR
  # calculate_thresholds <- function(q25, q75, iqr_multiplier) {
  #   iqr <- q75 - q25
  #   high_threshold <- max(0, q25 - iqr_multiplier * (iqr / 2))
  #   moderate_threshold <- max(0, q75 - iqr_multiplier * (iqr / 2))
  #   return(list(high = high_threshold, moderate = moderate_threshold))
  # }

  # Calculate Mahalanobis distances
  # mahalanobis_reliability <- mahalanobis_distances_testSet(geno_trn = geno_trn,
  #                                                         geno_tst = geno_tst,
  #                                                         geno_tst_trn = geno_tst_trn,
  #                                                         names_tst = names_tst,
  #                                                         names_trn = names_trn,
  #                                                         n_components = n_components,
  #                                                         threshold = threshold,
  #                                                         target = target,
  #                                                         apply_pca = apply_pca)
  #
  # if(isTRUE(is.na(mahalanobis_reliability))){
  #   return(list(trustworthiness=NA,
  #               proportion_high_reliability = NA,
  #               proportion_medium_reliability = NA,
  #               proportion_low_reliability = NA,
  #               reliability_percentage = NA,
  #               reliability_score = NA))
  # }

  # Calculate Mahalanobis distance thresholds
  # q75 <- stats::quantile(mahalanobis_dist, 0.75, na.rm = TRUE)
  # q25 <- stats::quantile(mahalanobis_dist, 0.25, na.rm = TRUE)
  # low_similarity_threshold <- q75 + iqr_multiplier * (q75 - q25)
  # high_similarity_threshold <- calculate_thresholds(q25, q75, iqr_multiplier)$high

  # Calculate interval width thresholds if necessary
  if((is.null(interval_width_high_threshold) & is.null(interval_width_low_threshold)) & !is.null(CI_width_thresholds)){
    high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1])
    low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2])
  } else if ((is.null(interval_width_high_threshold) | is.null(interval_width_low_threshold)) & !is.null(CI_width_thresholds)){
    high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1])
    low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2])
  } else{
    if ((is.null(interval_width_high_threshold) | is.null(interval_width_low_threshold)) & is.null(CI_width_thresholds)) {
      stop(print(paste(msg,'Both interval_width_high_threshold and interval_width_low_threshold must be either NULL or not NULL.')), call. = FALSE)

    }
  }
  # if (is.null(interval_width_moderate_threshold) & is.null(interval_width_high_threshold)) {
  #   interval_width_q75 <- stats::quantile(interval_width, 0.75, na.rm = TRUE)
  #   interval_width_q25 <- stats::quantile(interval_width, 0.25, na.rm = TRUE)
  #   thresholds <- calculate_thresholds(interval_width_q25, interval_width_q75, iqr_multiplier)
  #   interval_width_high_threshold <- thresholds$high
  #   interval_width_moderate_threshold <- thresholds$moderate
  # } else if (is.null(interval_width_moderate_threshold) | is.null(interval_width_high_threshold)) {
  #   stop("Both interval_width_moderate_threshold and interval_width_high_threshold must be NULL or not NULL.")
  # }

  # # Combine metrics to make a decision on trustworthiness
  # trustworthiness <- ifelse(interval_width <= interval_width_high_threshold & mahalanobis_dist <= high_similarity_threshold, "High",
  #                           ifelse(interval_width <= interval_width_moderate_threshold & mahalanobis_dist > high_similarity_threshold & mahalanobis_dist <= low_similarity_threshold, "Moderate", "Low"))
  #
  # noise_reliability <- ifelse(interval_width < interval_width_high_threshold, 1,
  #                             ifelse(interval_width > interval_width_moderate_threshold, 0, 0.5))

  noise_reliability <- ifelse(interval_width < interval_width_high_threshold, 1,
                              ifelse(interval_width < interval_width_low_threshold, 0.5, 0))
  # Aggregate reliability scores
  #reliability_score <- rowMeans(cbind(noise_reliability, mahalanobis_reliability))
  reliability_score <- noise_reliability
  # Calculate the proportion of high reliability predictions
  proportion_high_reliability <- mean(reliability_score == 1)

  # Decision rules
  # decision <- ifelse(reliability_score == 1, "High Reliability: Trust the predictions",
  #                    ifelse(reliability_score >= 0.5, "Medium Reliability: Use predictions with caution", "Low Reliability: Do not trust the predictions"))

  decision <- ifelse(reliability_score == 1, "high",
                     ifelse(reliability_score >= 0.5, "moderate", "low"))

  # Count occurrences the three classes
  count_high <- sum(decision == 'high')
  count_moderate <- sum(decision == 'moderate')
  count_low <- sum(decision == 'low')

  # Total number of elements
  total <- length(decision)

  # Calculate proportions
  proportion_high_reliability <- count_high / total
  proportion_moderate_reliability <- count_moderate / total
  proportion_low_reliability <- count_low / total

  #reliability_percentage <- (proportion_high_reliability + proportion_moderate_reliability)*100
  reliability_percentage <- (proportion_high_reliability)*100

  return(list(trustworthiness=decision,
              proportion_high_reliability = proportion_high_reliability,
              proportion_medium_reliability = proportion_moderate_reliability,
              proportion_low_reliability = proportion_low_reliability,
              reliability_percentage = reliability_percentage,
              reliability_score = reliability_score))


}



########
##
# Trade-off Between Precision and Coverage
#
# The goal is to find a balance between precision (narrow intervals) and
# coverage (intervals that contain the true value).
# This balance is typically managed by the
# confidence level of the intervals (e.g., 90%, 95%).

# lower MPIW Preferable if the intervals are well-calibrated
# and provide adequate coverage.
# Function to compute standardized metrics
# compute_MPIW_other_metrics <- function(predictions, boot_results,
#                             thresholds = c(0.33, 0.66)
#                             #true_values
#                             ) {
#   lower_bound <- apply(boot_results$t, 2, quantile, probs = 0.05)
#   upper_bound <- apply(boot_results$t, 2, quantile, probs = 0.95)
#   interval_width <- upper_bound - lower_bound
#
#   threshold1 <- quantile(interval_width, thresholds[1])
#   threshold2 <- quantile(interval_width, thresholds[2])
#   #coverage <- mean((true_values >= lower_bound) & (true_values <= upper_bound))
#
#   avg_interval_width <- mean(interval_width)
#   #proportion_high_reliability <- mean(interval_width < quantile(interval_width, thresholds[1]))
#
#   # Calculate proportions for each reliability level
#   proportion_low_reliability <- mean(interval_width >= threshold2)
#   proportion_medium_reliability <- mean(interval_width > threshold1 & interval_width < threshold2)
#   proportion_high_reliability <- mean(interval_width <= threshold1)
#
#  return(list(
#     #PICP = coverage,
#     MPIW = avg_interval_width,
#     proportion_high_reliability = proportion_high_reliability,
#     proportion_medium_reliability = proportion_medium_reliability,
#     proportion_low_reliability = proportion_low_reliability
#  ))
#
# }


# The Mean Prediction Interval Width (MPIW) is a metric used to assess
# the precision of prediction intervals generated by a model.
# The interpretation of MPIW depends on the context and
# the trade-off between interval width and coverage probability.
#
# Interpretation of MPIW
#
#     Lower MPIW:
#         Better Precision: A lower MPIW indicates that the prediction intervals
#                           are narrower, which suggests higher precision.
#         Risk of Undercoverage: While narrower intervals are more precise,
#                               they may not cover the true value as often,
#                               leading to undercoverage if not properly calibrated.
#
#     Higher MPIW:
#         Better Coverage: A higher MPIW indicates that the prediction intervals are wider,
#                          which suggests better coverage of the true values.
#         Lower Precision: Wider intervals are less precise and may provide
#                          less actionable information, but they are more likely
#                          to include the true value.
#
# Trade-off Between Precision and Coverage
#
# The goal is to find a balance between precision (narrow intervals)
# and coverage (intervals that contain the true value).
# This balance is typically managed by the confidence level of
# the intervals (e.g., 90%, 95%).
# Comparison Across Models
#
# When comparing MPIW across different models:
#
#     Lower MPIW: Preferable if the intervals are still well-calibrated and
#                provide adequate coverage.
#     Higher MPIW: Acceptable if it ensures better coverage,
#                   but it should be as low as possible while
#                  maintaining the desired coverage level.
