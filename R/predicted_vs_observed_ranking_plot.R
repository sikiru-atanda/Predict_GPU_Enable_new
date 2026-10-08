

##
# Kendall's Tau is a good metric when when interested in the ranking of the
# predicted values versus the observed values.
#Kendall's Tau will help you understand how well your model preserves the ranking of observations,
#which is particularly important in many genomic prediction applications.

##indicating how well the predicted rankings match the observed rankings.
## This ensures that the correlation measure accurately reflects the presence of ties in the data.
# Function to calculate Kendall's Tau
kendalls_tau <- function(x, y) {
  #browser()
  n <- length(x)
  if (n != length(y)) {
    stop("x and y must have the same length")
  }
  if (n < 2) {
    return(NA_real_)
  }

  concordant <- 0
  discordant <- 0
  ties_x <- 0
  ties_y <- 0

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sign_x <- sign(x[i] - x[j])
      sign_y <- sign(y[i] - y[j])
      if (sign_x == 0 && sign_y == 0) {
        # Both are tied
        ties_x <- ties_x + 1
        ties_y <- ties_y + 1
      } else if (sign_x == 0) {
        # x tied, but y not tied
        ties_x <- ties_x + 1
      } else if (sign_y == 0) {
        # y tied, but x not tied
        ties_y <- ties_y + 1
      } else if (sign_x == sign_y) {
        # Concordant pair
        concordant <- concordant + 1
      } else {
        # Discordant pair
        discordant <- discordant + 1
      }
    }
  }

  tau <- (concordant - discordant) / sqrt((concordant + discordant + ties_x) * (concordant + discordant + ties_y))
  return(tau)
}

# function to generate labels with "Top", "Middle", "Bottom"
generate_label <- function(segment, percent, status) {
  sprintf("%s %d%% - %s", segment, percent, status)
}


predicted_vs_observed_ranking_plot <- function(
                                               #GID_names = NULL,
                                               observed_value = NULL,
                                               predicted_value = NULL,
                                               prediction_error_var = NULL,
                                               standard_errors = NULL,
                                               percent_top = 20,
                                               percent_bottom =  20,
                                               acceptable_se_threshold = 0.5,
                                               high_reliability_thres = 0.7,
                                               low_reliability_thres = 0.4,
                                               confidence_level = 0.95,
                                               CI_width_thresholds = c(0.33, 0.66),
                                               abs_very_close_threshold = 0.01,
                                               abs_close_threshold = 0.05,
                                               cv_results = NULL,
                                               replication = NULL,
                                               path_plot = NULL,
system_database = FALSE){
#browser()

safe_mean_logical <- function(x) {
  if (!length(x) || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

df <- data.frame(
  #genotype = GID_names,
  observed = observed_value,
  predicted = predicted_value,
  stringsAsFactors = FALSE
)


# Ranking data
df <- df |>
  dplyr::mutate(
    obs_rank =  rank(-observed, ties.method = "first"),
    pred_rank =  rank(-predicted, ties.method = "first")
  )

# User-defined percentages
percent_middle <- 100 - (percent_top + percent_bottom)

# Calculate quantiles based on user-defined percentages
df <- df |>
  dplyr::mutate(
    total =  dplyr::n(),
    top_quantile = total * (percent_top / 100),
    bottom_quantile = total * (1 - percent_bottom / 100),

    # Calculate the category with case_when first
    category =  dplyr::case_when(
      obs_rank <= top_quantile & pred_rank <= top_quantile ~ generate_label("Top", percent_top, "Consistent"),
      obs_rank <= top_quantile & pred_rank > top_quantile ~ generate_label("Top", percent_top, "Changed"),
      obs_rank > top_quantile & obs_rank <= (total - bottom_quantile) ~ generate_label("Middle", percent_middle, "Middle"),
      obs_rank > bottom_quantile & pred_rank > bottom_quantile ~ generate_label("Bottom", percent_bottom, "Consistent"),
      obs_rank > bottom_quantile & pred_rank <= bottom_quantile ~ generate_label("Bottom", percent_bottom, "Changed"),
      TRUE ~ generate_label("Middle", percent_middle, "Middle")
    )
  )

# Ensure unique factor levels before converting 'category' to a factor

unique_levels <- c(generate_label("Top", percent_top, "Consistent"),
                  generate_label("Top", percent_top, "Changed"),
                  generate_label("Middle", percent_middle, "Middle"),
                  generate_label("Bottom", percent_bottom, "Consistent"),
                  generate_label("Bottom", percent_bottom, "Changed"))

# Add additional mutations and factor assignment
df <- df |>
  dplyr::mutate(
    category = factor(category, levels = unique_levels),
    # Determine consistency within specific quantiles
    top_20_obs = obs_rank <= top_quantile,
    top_20_pred = pred_rank <= top_quantile,
    bottom_20_obs = obs_rank > bottom_quantile,
    bottom_20_pred = pred_rank > bottom_quantile
  )

# Filter for top and bottom 20% subsets
top_20_subset <- df |> dplyr::filter(top_20_obs | top_20_pred)
bottom_20_subset <- df |> dplyr::filter(bottom_20_obs | bottom_20_pred)

# Calculate the consistency rates for the top and bottom 20% subsets
top_20_consistency_rate <- safe_mean_logical(top_20_subset$top_20_obs == top_20_subset$top_20_pred)
bottom_20_consistency_rate <- safe_mean_logical(bottom_20_subset$bottom_20_obs == bottom_20_subset$bottom_20_pred)

#mydf$task <- factor(mydf$task, levels = levels_index)
# Calculate the consistency rates
# top_20_consistency_rate <- mean(df$top_20_obs == df$top_20_pred) * 100
# bottom_20_consistency_rate <- mean(df$bottom_20_obs == df$bottom_20_pred) * 100
#

# Example usage
tau_value <- kendalls_tau(df$observed, df$predicted)
# cat("Kendall's Tau:", tau, "\n")

# Calculate Kendall's Tau
# <- kendalls_tau(df$obs_rank, df$pred_rank)
# cat("Top 20% Consistency Rate: ", top_20_consistency_rate * 100, "%\n")
# cat("Bottom 20% Consistency Rate: ", bottom_20_consistency_rate * 100, "%\n")
# cat("Overall Kendall's Tau: ", tau, "\n")


#unique_levels <- as.character(unique(df$category))
# Plot


p1 <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted, color = category)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = stats::setNames(
        c("blue", "lightblue", "grey", "red", "pink"),
        unique_levels
      )) +
      # scale_color_manual(values = c("Top 20% - Consistent" = "blue",
      #                               "Top 20% - Changed" = "lightblue",
      #                               "Middle 60%" = "grey",
      #                               "Bottom 20% - Consistent" = "red",
      #                               "Bottom 20% - Changed" = "pink")) +
      ggplot2::labs(title = "Comparing Predicted and Observed Rankings Across Validation Sets",
           subtitle = sprintf("Kendall's Tau: %.2f | Top %d%% Consistency: %.2f%% | Bottom %d%% Consistency: %.2f%%",
                              tau_value, percent_top, top_20_consistency_rate*100, percent_bottom, bottom_20_consistency_rate*100),
           # subtitle = sprintf("Kendall's Tau: %.2f, Top 20%% Consistency: %.1f%%, Bottom 20%% Consistency: %.1f%%",
           #                    tau_value, top_20_consistency_rate, bottom_20_consistency_rate),
           x = "Observed Values",
           y = "Predicted Values",
           color = "Category") +
     ggplot2::theme_minimal() +
     ggplot2::theme(
     plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"),
     #legend.position = "right",
     legend.title = ggplot2::element_text(size = 12),
     legend.text = ggplot2::element_text(size = 10)) +
     ggplot2::ylim(c(min(df$predicted), max(df$predicted)))  # Ensure y-axis properly shows the range


# svm_predict <- function(data, indices, kernel_params) {
#
#   samplegeno <-  data[indices, -1]
#   y <-  data[indices, 1]
#   svm_model <- svm(y ~ ., data=samplegeno)
#   #svm_model <- svm(y = y,  x=samplegeno)
#   return(predict(svm_model, data[, -1]))
# }


# Define a threshold for acceptable standard error
#acceptable_se_threshold <- 0.5

# res <- reliability_thresholds(prediction_error_var = prediction_error_var,
#                                genetic_var = var(df$predicted),
#                                high_reliability_thres = high_reliability_thres,
#                                low_reliability_thres = low_reliability_thres)
#
#
# predicted_values = data.frame(Predicted_value = df$predicted,
#                               Reliability = res$reliability,
#                               Remarks = res$remarks)
#
# unique_remakers = c("Reliable", "Acceptable", "Unreliable")
#
# # Convert the 'Remarks' column to a factor with specified levels
# predicted_values$Remarks <- factor(predicted_values$Remarks, levels = c("Reliable", "Acceptable", "Unreliable"))
#
# # Classify predictions based on this threshold
# #predicted_values$Remarks <- predicted_values$Reliability > acceptable_se_threshold
#
#
#
# p2 <- ggplot2::ggplot(predicted_values, ggplot2::aes(x = Reliability, y =Predicted_value , color = Remarks)) +
#       ggplot2::geom_point(alpha = 0.5) +
#       ggplot2::labs(title = "Prediction Confidence Analysis",
#                     subtitle = sprintf("Proportion of reliability: high, medium, and low: %.2f%%, %.2f%%, and %.2f%%",
#                                        res$proportion_high_reliability * 100, res$proportion_medium_reliability * 100, res$proportion_low_reliability * 100),
#                     # subtitle = sprintf("Proportion of high reliability: %.2f%%, Proportion of medium reliability: %.2f%%, Proportion of low reliability: %.2f%%",
#                     #                    res$proportion_high_reliability*100, res$proportion_medium_reliability*100, res$proportion_low_reliability*100),
#            x = "Reliability",
#            y = "Predicted Values") +
#       ggplot2::scale_color_manual(values = setNames(
#         c("green", "lightblue",  "red"),
#         unique_remakers
#       )) +
#       #scale_color_manual(values = c("red", "green"), labels = c("Unreliable", "Reliable")) +
#       ggplot2::theme_minimal()+
#       ggplot2::theme(
#         plot.title = ggplot2::element_text(hjust = 0.5)
#         #plot.subtitle = ggplot2::element_text(hjust = 0.5)
#       )


# res_CI <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
#                                               predictions = df$predicted,
#                                               standard_errors = standard_errors,
#                                               confidence_level = confidence_level,
#                                               model_for_CI_cal = "None",
#                                               boot_results = NULL,
#                                               cv_results = cv_results,
#                                               replication = replication)
#
# testData <- data.frame(predictions = df$predicted,
#                        lower_bound = res_CI$lower_bound,
#                        upper_bound = res_CI$upper_bound,
#                        reliability_remarks = res_CI$reliability_remarks)
#
# unique_remakers_CI <-  c("high", "moderate", "low")
#
# testData$reliability_remarks <- factor(testData$reliability_remarks, levels = unique_remakers_CI)
#

# p3 <- ggplot2::ggplot(testData, ggplot2::aes(x = predictions, y = predictions,
#                      ymin = lower_bound, ymax = upper_bound)) +
#       ggplot2::geom_errorbar() +
#       #geom_point() +
#       ggplot2::geom_point(ggplot2::aes(color = reliability_remarks)) +
#       ggplot2::scale_color_manual(values = stats::setNames(
#         c("green", "lightblue",  "red"),
#         unique_remakers_CI
#       )) +
#       ggplot2::theme_minimal() +
#       ggplot2::labs(title = "Prediction Intervals and Reliability",
#            x = "Predicted Values",
#            y = "Predicted Values",
#            color = "Reliability") ## This change the title of the legend


# Define thresholds
#abs_very_close_threshold <- 0.01
#abs_close_threshold <- 0.05

# Calculate the absolute differences
df$difference <- abs(df$observed - df$predicted)

# Categorize differences
df$category <- ifelse(df$difference <= abs_very_close_threshold, "Highly Accurate",
                      ifelse(df$difference <= abs_close_threshold, "Moderate Error", "Significant Error"))

# Calculate proportions
proportion_highly_accurate <- sum(df$category == "Highly Accurate") / nrow(df)
proportion_moderate_error <- sum(df$category == "Moderate Error") / nrow(df)
proportion_significant_error <- sum(df$category == "Significant Error") / nrow(df)

# Calculate correlation and R-squared
correlation <- suppressWarnings(cor(df$observed, df$predicted, use = "complete.obs"))
r_squared <- correlation^2

# Create subtitle using sprintf
subtitle <- sprintf("Highly Accurate: %.2f%% | Moderate Error: %.2f%% | Significant Error: %.2f%%",
                    proportion_highly_accurate * 100, proportion_moderate_error * 100, proportion_significant_error * 100)

# Create the plot
fit_df <- data.frame(
  observed = df$observed,
  predicted = df$predicted,
  fit_type = "Linear Fit",
  stringsAsFactors = FALSE
)

can_loess <- nrow(df) >= 4 && length(unique(df$observed)) >= 3 && length(unique(df$predicted)) >= 3
if (can_loess) {
  fit_df <- rbind(
    fit_df,
    data.frame(
      observed = df$observed,
      predicted = df$predicted,
      fit_type = "LOESS Fit",
      stringsAsFactors = FALSE
    )
  )
}

p4 <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted)) +
      ggplot2::geom_point(ggplot2::aes(color = difference, shape = category), size = 2) +
      ggplot2::geom_smooth(
        data = fit_df[fit_df$fit_type == "Linear Fit", , drop = FALSE],
        ggplot2::aes(linetype = fit_type),
        formula = y ~ x,
        method = "lm",
        se = FALSE,
        color = "red"
      ) +
      ggplot2::scale_color_gradient(low = "blue", high = "red", name = "Abs.diff") +
      ggplot2::scale_shape_manual(values = c("Highly Accurate" = 16, "Moderate Error" = 17, "Significant Error" = 18), name = "Category") +
      ggplot2::labs(
        title = if (can_loess) "Evaluation of Observed vs Predicted Values with Linear and LOESS Fits" else "Evaluation of Observed vs Predicted Values with Linear Fit",
        subtitle = subtitle,
        y = "Predicted Values",
        x = "Observed Values",
        color = "Absolute Difference",
        shape = "Accuracy Category",
        linetype = "Fit Type"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 10),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      ) +
      ggplot2::annotate(
        "text",
        x = min(df$observed),
        y = max(df$predicted),
        label = paste("Predictive ability:", round(correlation, 3)),
        hjust = 0,
        vjust = 1,
        size = 4,
        fontface = "bold",
        color = "darkblue"
      )

if (can_loess) {
  p4 <- p4 + ggplot2::geom_smooth(
    data = fit_df[fit_df$fit_type == "LOESS Fit", , drop = FALSE],
    ggplot2::aes(linetype = fit_type),
    formula = y ~ x,
    method = "loess",
    se = FALSE,
    color = "green",
    span = 1
  )
}

# Display the plot


# Plot
# Calculate correlation coefficient and R-squared value
# correlation <- cor(df$observed, df$predicted, use = "complete.obs")
# r_squared <- correlation^2
# # Create a new column for the absolute difference
# df$difference <- abs(df$observed - df$predicted)
#
# # Create the plot
# p4 <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted)) +
#   ggplot2::geom_point(ggplot2::aes(color = difference), size = 2) +  # Scatter plot points colored by the difference
#   ggplot2::geom_smooth(ggplot2::aes(linetype = "Linear Fit"), formula = y ~ x, method = "lm", se = FALSE, color = "red") +  # Linear fit line (red)
#   ggplot2::geom_smooth(ggplot2::aes(linetype = "LOESS Fit"), formula = y ~ x, method = "loess", se = FALSE, color = "green", span = 0.5) + # LOESS fit line (green)
#   ggplot2::scale_color_gradient(low = "blue", high = "red", name = "Abs.diff") + # Gradient color scale
#   ggplot2::labs(
#     title = "Evaluation of Observed vs Predicted Values with Linear and LOESS Fits",
#     y = "Predicted Values",
#     x = "Observed Values",
#     color = "Absolute Difference",
#     linetype = "Fit Type"
#   ) +
#   ggplot2::theme_minimal() +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
#     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "darkblue"),
#     legend.position = "right",
#     legend.title = ggplot2::element_text(size = 12),
#     legend.text = ggplot2::element_text(size = 10),
#     plot.margin = ggplot2::margin(10, 10, 10, 10)
#   ) +
#   ggplot2::annotate("text", x = min(df$observed), y = max(df$predicted),
#            label = paste("R-squared:", round(r_squared, 3), "\nPredictive ability:", round(correlation, 3)),
#            hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "darkblue")



  # p4 <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted)) +
  #   ggplot2::geom_point(ggplot2::aes(color = difference), size = 2) +  # Scatter plot points colored by the difference
  #   ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = FALSE, color = "red", linetype = "solid") +  # Linear fit line (red)
  #   ggplot2::geom_smooth(formula = y ~ x, method = "loess", se = FALSE, color = "green", linetype = "dashed", span = 0.5) + # LOESS fit line (green)
  #   ggplot2::scale_color_gradient(low = "blue", high = "red", name = "Abs.diff") + # Gradient color scale
  #   ggplot2::labs(
  #     title = "Observed vs Predicted Values with Linear and LOESS Fits",
  #     y = "Predicted Values",
  #     x = "Observed Values"
  #   ) +
  #   ggplot2::theme_minimal() +
  #   ggplot2::theme(
  #     plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
  #     legend.position = "right",
  #     legend.title = ggplot2::element_text(size = 12),
  #     legend.text = ggplot2::element_text(size = 10),
  #     plot.margin = ggplot2::margin(10, 10, 10, 10)
  #   ) +
  #   ggplot2::annotate("text", x = min(df$observed), y = max(df$predicted),
  #                     label = paste("R-squared:", round(r_squared, 3), "\nPredictive ability:", round(correlation, 3)),
  #                     hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "darkblue")
  #

# p4 <- ggplot2::ggplot(df, ggplot2::aes(y = predicted, x = observed)) +
#       ggplot2::geom_point(color = "black") +                                # Scatter plot points
#       #geom_smooth(method = "lm", se = FALSE, color = "red") +      # Linear fit line (red)
#       ggplot2::geom_smooth(method = "lm", se = FALSE,  ggplot2::aes(color = "Linear Fit")) +
#       ggplot2::geom_smooth(method = "loess", se = FALSE, ggplot2::aes(color = "LOESS Fit"), span = 0.5) +
#       #geom_smooth(method = "loess", se = FALSE, color = "green", span = 0.5) + # LOESS fit line (green)
#       #geom_abline(slope = 1, intercept = 0, linetype = "dashed") + # Ideal line (diagonal)
#       ggplot2::scale_color_manual(values = c("Linear Fit" = "red", "LOESS Fit" = "green")) + # Define colors
#       ggplot2::labs(
#         title = "JMP plot",
#         y = "Predicted Values",
#         x = "Observed Values",
#         color = "Fit Type"
#       ) +
#       ggplot2::theme_minimal()


#if(isFALSE(system_database)){
combined_plot <- gp_arrange_grob_safely(p1, p4, ncol = 2)

# name_plot <- paste(paste("Cross_validation_diagonistic_plots", model, sep = "_"), trait, sep = "_")
#
# plot_save <- ggplot2::ggsave(paste(name_plot,"pdf", sep = "."), plot = combined_plot,
#                              width = 15, height = 8, units = "in", dpi = 300, path = path_plot)
# # Save the combined plot
return(combined_plot)


# } else {
#
#  return(list(predicted_vs_observed_ranking = p1,
#        #predicted_vs_reliability = p2,
#        #prediction_inter_vs_reliability = p3,
#        predicted_vs_observed = p4
#
#   )
#  )
# }

}


# kendalls_tau <- function(x, y) {
#   n <- length(x)
#   concordant <- 0
#   discordant <- 0
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       concordant <- concordant + (sign(x[i] - x[j]) == sign(y[i] - y[j]) && x[i] != x[j] && y[i] != y[j])
#       discordant <- discordant + (sign(x[i] - x[j]) != sign(y[i] - y[j]) && x[i] != x[j] && y[i] != y[j])
#     }
#   }
#   tau <- (concordant - discordant) / (n * (n - 1) / 2)
#   return(tau)
# }
#
# kendalls_tau <- function(x, y) {
#   n <- length(x)
#   concordant <- 0
#   discordant <- 0
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       if (x[i] != x[j] && y[i] != y[j]) {
#         concordant <- concordant + (sign(x[i] - x[j]) == sign(y[i] - y[j]))
#         discordant <- discordant + (sign(x[i] - x[j]) != sign(y[i] - y[j]))
#       }
#     }
#   }
#   tau <- (concordant - discordant) / (n * (n - 1) / 2)
#   return(tau)
# }
#
# # Example calculation
# tau_value <- kendalls_tau(df$obs_rank, df$pred_rank)
# cat("Overall Kendall's Tau: ", tau_value, "\n")


# Example reliability categorization function
# reliability_thresholds <- function(interval_width, thresholds = c(0.33, 0.66)) {
#   threshold1 <- quantile(interval_width, thresholds[1])
#   threshold2 <- quantile(interval_width, thresholds[2])
#
#   reliability <- ifelse(interval_width < threshold1, "high",
#                         ifelse(interval_width < threshold2, "moderate", "low"))
#   return(reliability)
# }
#
# # Function to compute standardized metrics
# compute_metrics <- function(predictions, lower_bound, upper_bound, true_values) {
#   interval_width <- upper_bound - lower_bound
#   #coverage <- mean((true_values >= lower_bound) & (true_values <= upper_bound))
#
#   avg_interval_width <- mean(interval_width)
#   proportion_high_reliability <- mean(interval_width < quantile(interval_width, 0.33))
#
#
#   list(
#     #PICP = coverage,
#     MPIW = avg_interval_width,
#     proportion_high_reliability = proportion_high_reliability
#   )
# }

# True values (for validation)
#true_values <- rnorm(100)

# Compute metrics for asreml model
# metrics_asreml <- compute_metrics(testData2$predictions, testData2$lower_bound,
#                                   testData2$upper_bound
#                                   #true_values
# )
#
# rownames(testData) <- rownames(geno_data)[test_index]
# rownames(testData2) <- rownames(geno_data)[test_index]

# sik = reliability_thresholds_MPIW_from_CI(boot_results = boot_results)
