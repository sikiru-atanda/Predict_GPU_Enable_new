#' Title
#'
#' @param boot_results
#' @param GID_names
#' @param CI_width_thresholds
#' @param predictions
#' @param standard_errors
#' @param prediction_error_var
#' @param genetic_var
#' @param reference_variance
#' @param confidence_level
#' @param risk_calibration Held-out or cross-fitted calibration information for
#'   Gaussian ML/DL prediction error and interval diagnostics.
#' @param require_heldout_calibration If `TRUE`, the default, Gaussian ML/DL
#'   diagnostics fail rather than substituting bootstrap instability or
#'   in-sample standard errors when held-out calibration is unavailable.
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
                                           Predicted_value_for_CI = NULL,
                                           mod = NULL,
                                           CI_width_thresholds = c(0.33, 0.66),
                                           predictions = NULL,
                                           standard_errors = NULL,
                                           prediction_error_var = NULL,
                                           genetic_var = NULL,
                                           reference_variance = NULL,
                                           lower_bound = NULL,
                                           upper_bound = NULL,
                                           observed_y = NULL,
                                           train_test_label = NULL,
                                           confidence_level = 0.95,
                                           risk_calibration = NULL,
                                           require_heldout_calibration = TRUE,
                                           model_for_CI_cal = "ML",
                                           response_family = "gaussian",
                                           composite_reliability_score = NULL,
                                           composite_reliability = NULL,
                                           composite_reliability_percentage = NULL,
                                           #threshold = NULL,
                                           high_reliability_thres = 0.7,
                                           low_reliability_thres = 0.4,
                                           system_database = FALSE){

  msg <- ""
  fam <- gp_resolve_response_family(response_family, y = predictions)
  ml_gaussian <- identical(model_for_CI_cal, "ML") && identical(fam, "gaussian")
  heldout_summary <- if (isTRUE(ml_gaussian)) {
    gp_ml_heldout_residual_summary(
      calibration = risk_calibration,
      confidence_level = confidence_level
    )
  } else {
    NULL
  }
  if (isTRUE(ml_gaussian) && isTRUE(require_heldout_calibration) &&
      is.null(heldout_summary)) {
    stop(
      paste0(
        "Gaussian ML/DL uncertainty diagnostics require held-out or cross-fitted ",
        "calibration. Supply `risk_calibration`, or set ",
        "`require_heldout_calibration = FALSE` to explicitly request the legacy ",
        "bootstrap/standard-error diagnostic, whose nominal coverage is not guaranteed."
      ),
      call. = FALSE
    )
  }
  cal_res <- NULL

  if (is.null(composite_reliability) || length(composite_reliability) == 0 || any(is.na(composite_reliability)) ||
      is.null(composite_reliability_score) || length(composite_reliability_score) == 0 || any(is.na(composite_reliability_score)) ||
      (!is.null(predictions) && length(composite_reliability) != length(predictions)) ||
      (!is.null(predictions) && length(composite_reliability_score) != length(predictions))) {

    composite_reliability <- NULL
    composite_reliability_score <- NULL
    p111 <- NULL
  }

  if(model_for_CI_cal == "RKHS" | model_for_CI_cal == "GBLUP") {
  if(is.null(prediction_error_var) & !is.null(standard_errors)){
  prediction_error_var <- (standard_errors)^2

  } else {
    if(is.null(prediction_error_var) & is.null(standard_errors)){

      stop(print(paste(msg,'Provide vector of prediction error or standard error of the predictions.')), call. = FALSE)


    }
  }
    if(is.null(genetic_var))  stop(print(paste(msg,'Provide genetic variance.')), call. = FALSE)
}

  if(model_for_CI_cal == "Bayes") {

    prediction_error_var <- (standard_errors)^2
    #prediction_error_var <-  apply(Predicted_value_for_CI, 1, var)
    #standard_errors <- apply(Predicted_value_for_CI, 1, sd)

  }


  if(!is.null(predictions) & model_for_CI_cal == "RKHS" | model_for_CI_cal == "GBLUP"){
    res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                  predictions = predictions,
                                                  standard_errors = standard_errors,
                                                  confidence_level = confidence_level,
                                                  model_for_CI_cal = model_for_CI_cal,
                                                  boot_results = boot_results)

    predictions <- res_CI$predictions
    #genetic_var <- var(predictions)
    prediction_error_var <- res_CI$prediction_error_var

  } else if(model_for_CI_cal == "Bayes"){

    res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                  predictions = predictions,
                                                  Predicted_value_for_CI = Predicted_value_for_CI,
                                                  mod = mod,
                                                  standard_errors = standard_errors,
                                                  confidence_level = confidence_level,
                                                  model_for_CI_cal = model_for_CI_cal,
                                                  boot_results = boot_results)

  } else if (isTRUE(ml_gaussian) && !is.null(predictions) &&
             (!is.null(standard_errors) || !is.null(heldout_summary))) {

    if (is.null(prediction_error_var)) {
      prediction_error_var <- standard_errors^2
    }
    stability_prediction_error_var <- prediction_error_var
    if (is.null(lower_bound) || is.null(upper_bound)) {
      z <- stats::qnorm(1 - (1 - confidence_level) / 2)
      lower_bound <- predictions - z * standard_errors
      upper_bound <- predictions + z * standard_errors
    }

    cal_res <- gp_ml_calibrate_gaussian_uncertainty(
      predicted_value = predictions,
      pred_se = standard_errors,
      pred_variances = prediction_error_var,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      observed_y = observed_y,
      train_test_label = train_test_label,
      boot_results = boot_results,
      heldout_calibration = risk_calibration,
      confidence_level = confidence_level,
      require_heldout_calibration = require_heldout_calibration
    )

    interval_width <- cal_res$calibrated_upper_bound -
      cal_res$calibrated_lower_bound
    finite_width <- interval_width[is.finite(interval_width)]
    high_threshold <- if (length(finite_width)) {
      as.numeric(stats::quantile(
        finite_width,
        probs = CI_width_thresholds[[1L]],
        names = FALSE
      ))
    } else {
      NA_real_
    }
    low_threshold <- if (length(finite_width)) {
      as.numeric(stats::quantile(
        finite_width,
        probs = CI_width_thresholds[[2L]],
        names = FALSE
      ))
    } else {
      NA_real_
    }
    interval_remarks <- ifelse(
      interval_width < high_threshold,
      "low",
      ifelse(interval_width < low_threshold, "moderate", "high")
    )

    res_CI <- list(
      predictions = predictions,
      lower_bound = cal_res$calibrated_lower_bound,
      upper_bound = cal_res$calibrated_upper_bound,
      reliability_remarks = interval_remarks,
      relative_interval_width = interval_width / predictions,
      standard_errors = cal_res$calibrated_standard_error,
      high_threshold = high_threshold
    )
    prediction_error_var <- stability_prediction_error_var
  } else {

    res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                  predictions = NULL,
                                                  standard_errors = NULL,
                                                  confidence_level = confidence_level,
                                                  model_for_CI_cal = model_for_CI_cal,
                                                  boot_results = boot_results)

    predictions <- res_CI$predictions
    #genetic_var <- var(predictions)
    prediction_error_var <- res_CI$prediction_error_var
    #stop(print(paste(msg,'Provide vector of the predicted values.')), call. = FALSE)

  }

  if (identical(model_for_CI_cal, "ML") && fam != "gaussian") {
    conf <- pmax(as.numeric(predictions), 1 - as.numeric(predictions))
    conf <- pmax(0, pmin(1, conf))
    remarks <- dplyr::case_when(
      conf >= high_reliability_thres ~ "High Confidence",
      conf <= low_reliability_thres ~ "Low Confidence",
      TRUE ~ "Moderate Confidence"
    )
    predicted_values <- data.frame(
      Predicted_value = predictions,
      Metric = conf,
      Remarks = factor(remarks, levels = c("High Confidence", "Moderate Confidence", "Low Confidence"))
    )
    p11 <- ggplot2::ggplot(predicted_values, ggplot2::aes(x = Metric, y = Predicted_value, color = Remarks)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::ggtitle("Prediction Confidence Analysis") +
      ggplot2::labs(
        subtitle = sprintf(
          "Proportion of prediction confidence: high, moderate, and low: %.2f%%, %.2f%%, and %.2f%%",
          mean(conf >= high_reliability_thres, na.rm = TRUE) * 100,
          mean(conf < high_reliability_thres & conf > low_reliability_thres, na.rm = TRUE) * 100,
          mean(conf <= low_reliability_thres, na.rm = TRUE) * 100
        ),
        x = "Prediction Confidence",
        y = "Predicted Probability"
      ) +
      ggplot2::scale_color_manual(values = c("High Confidence" = "green", "Moderate Confidence" = "lightblue", "Low Confidence" = "red")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue")
      )

    testData <- data.frame(
      predictions = predictions,
      lower_bound = pmax(0, pmin(1, res_CI$lower_bound)),
      upper_bound = pmax(0, pmin(1, res_CI$upper_bound)),
      confidence_remarks = factor(
        dplyr::case_when(
          res_CI$reliability_remarks == "low" ~ "High Confidence",
          res_CI$reliability_remarks == "moderate" ~ "Moderate Confidence",
          TRUE ~ "Low Confidence"
        ),
        levels = c("High Confidence", "Moderate Confidence", "Low Confidence")
      )
    )

    p12 <- ggplot2::ggplot(testData, ggplot2::aes(x = predictions, y = predictions, ymin = lower_bound, ymax = upper_bound)) +
      ggplot2::geom_errorbar() +
      ggplot2::geom_point(ggplot2::aes(color = confidence_remarks)) +
      ggplot2::scale_color_manual(values = c("High Confidence" = "green", "Moderate Confidence" = "lightblue", "Low Confidence" = "red")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Prediction Intervals and Confidence",
        x = "Predicted Probability",
        y = "Predicted Probability",
        color = "Prediction Confidence"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue")
      )

    if (isTRUE(system_database)) {
      return(list(
        predicted_vs_reliability = p11,
        prediction_inter_vs_reliability = p12,
        predicted_vs_composite_reliability = NULL
      ))
    }

    return(gp_arrange_grob_safely(p11, p12, ncol = 2))
  }

  if (identical(model_for_CI_cal, "ML")) {
    ref_var <- reference_variance %||% genetic_var
    if (is.null(ref_var) || anyNA(ref_var)) {
      ref_var <- stats::var(predictions, na.rm = TRUE)
    }
    stability_res <- gp_ml_gaussian_stability(
      prediction_error_var = prediction_error_var,
      reference_variance = ref_var,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )
    res <- list(
      reliability = stability_res$stability,
      remarks = stability_res$remarks,
      reliability_percentage = stability_res$stability_percentage,
      proportion_high_reliability = stability_res$proportion_high_reliability,
      proportion_medium_reliability = stability_res$proportion_medium_reliability,
      proportion_low_reliability = stability_res$proportion_low_reliability
    )
  } else {
    res <-  reliability_thresholds(prediction_error_var = prediction_error_var,
                                   genetic_var =  genetic_var,
                                   high_reliability_thres = high_reliability_thres,
                                   low_reliability_thres = low_reliability_thres)
  }

  metric_label <- if (identical(model_for_CI_cal, "ML")) "Prediction Stability" else "Reliability"
  metric_title <- if (identical(model_for_CI_cal, "ML")) "Prediction Stability Analysis" else "Prediction Confidence Analysis"
  subtitle_label <- if (identical(model_for_CI_cal, "ML")) "prediction stability" else "reliability"
  risk_basis_note <- NULL
  if (isTRUE(ml_gaussian) && !is.null(cal_res)) {
    target_uncertainty_status <- if (any(is.finite(
      cal_res$calibrated_prediction_error_var
    ))) {
      "Target-specific PEV available"
    } else {
      "Target-specific SE/PEV unavailable (variation <= 1%)"
    }
    risk_detail <- if (grepl(
      "row-specific bootstrap prediction instability",
      cal_res$prediction_error_estimand,
      fixed = TRUE
    )) {
      "row-specific bootstrap prediction instability"
    } else if (grepl(
      "unavailable",
      cal_res$prediction_error_estimand,
      ignore.case = TRUE
    )) {
      "Target-specific SE/PEV unavailable"
    } else {
      "target-specific cross-fitted predictive uncertainty"
    }
    risk_basis_note <- paste0(
      risk_detail,
      "; marginal cross-fitted interval; calibration n = ",
      cal_res$prediction_interval_calibration_n
    )
  } else if (isTRUE(ml_gaussian) && !isTRUE(require_heldout_calibration)) {
    risk_basis_note <- paste(
      "Legacy bootstrap/standard-error diagnostic;",
      "nominal prediction-interval coverage is not guaranteed"
    )
  }


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

  res$reliability <- pmax(0, pmin(1, res$reliability))
  predicted_values = data.frame(Predicted_value = predictions,
                                Metric = res$reliability,
                                Remarks = res$remarks)

  unique_remakers = if (identical(model_for_CI_cal, "ML")) {
    c("Stable", "Moderately Stable", "Unstable")
  } else {
    c("Reliable", "Acceptable", "Unreliable")
  }

  # Convert the 'Remarks' column to a factor with specified levels
  predicted_values$Remarks <- factor(predicted_values$Remarks, levels = unique_remakers)

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
                       Metric = composite_reliability_score,
                       Remarks = composite_reliability)


    p111 <- ggplot2::ggplot(data, ggplot2::aes(x = Metric, y =Predicted_value , color = Remarks)) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::labs(title = metric_title,
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

  p11_subtitle <- sprintf("Proportion of %s: high, medium, and low: %.2f%%, %.2f%%, and %.2f%%",
                          subtitle_label,
                          res$proportion_high_reliability * 100,
                          res$proportion_medium_reliability * 100,
                          res$proportion_low_reliability * 100)
  if (!is.null(risk_basis_note)) {
    p11_subtitle <- paste(p11_subtitle, risk_basis_note, sep = "\n")
  }

  p11 <- ggplot2::ggplot(predicted_values, ggplot2::aes(x = Metric, y =Predicted_value , color = Remarks)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::ggtitle(metric_title) +
        ggplot2::labs(subtitle = p11_subtitle,
                      # subtitle = sprintf("Proportion of high reliability: %.2f%%, Proportion of medium reliability: %.2f%%, Proportion of low reliability: %.2f%%",
                      #                    res$proportion_high_reliability*100, res$proportion_medium_reliability*100, res$proportion_low_reliability*100),
                      x = metric_label,
                      y = "Predicted Values") +
        ggplot2::scale_color_manual(values = setNames(
          c("green", "lightblue",  "red"),
          unique_remakers
        )) +
        #scale_color_manual(values = c("red", "green"), labels = c("Unreliable", "Reliable")) +
        ggplot2::theme_minimal()+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"),
                         plot.margin = ggplot2::margin(t = 18, r = 10, b = 10, l = 10))


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



  interval_title <- if (identical(model_for_CI_cal, "ML")) {
    "Prediction Intervals and Stability"
  } else {
    "Prediction Intervals and Reliability"
  }
  p12_title <- if (
    identical(model_for_CI_cal, "ML") &&
      !is.null(risk_basis_note) &&
      grepl(
        "nominal prediction-interval coverage is not guaranteed",
        risk_basis_note,
        fixed = TRUE
      )
  ) {
    paste(interval_title, risk_basis_note, sep = " - ")
  } else {
    interval_title
  }
  p12_subtitle <- risk_basis_note

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
        ggplot2::labs(title = p12_title,
                      subtitle = p12_subtitle,
                      x = "Predicted Values",
                      y = "Predicted Values",
                      color = metric_label)+ ## This change the title of the legend
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"),
          plot.margin = ggplot2::margin(t = 18, r = 10, b = 10, l = 10)
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

    combined_plot <- gp_arrange_grob_safely(p11, p12, ncol = 2)

    }else{

      #combined_plot <- gridExtra::grid.arrange(p11, p111, p12, ncol = 2)
      combined_plot <- gp_arrange_grob_safely(p11, p12, ncol = 2)

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
