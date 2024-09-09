
# Function to calculate statistics
calculate_statistics <- function(pred_matrix) {
  mean_pred <- rowMeans(pred_matrix)
  # lower_bound <- apply(pred_matrix, 1, quantile, probs = 0.05)
  # upper_bound <- apply(pred_matrix, 1, quantile, probs = 0.95)
  # interval_width <- upper_bound - lower_bound
  # std_error <- apply(pred_matrix, 1, sd)
  # pred_error_var <- apply(pred_matrix, 1, var)

  return(list(
    pred_mean = mean_pred
    # lower = lower_bound,
    # upper = upper_bound,
    # interval_width = interval_width,
    # std_error = std_error,
    # pred_error_var = pred_error_var
  ))

}

single_predicted_vs_observed_result_plots_process <- function(results,
                                                              pheno_data,
                                                              abs_very_close_threshold ,
                                                              abs_close_threshold){
#browser()
  combined_results <- data.frame()

  # Loop through each element in the results list
  for (i in seq_along(results)) {
    # Extract the current element
    current_result <- results[[i]]

    # Extract the ypred_cv_Reps_all data frame
    ypred_data <- as.data.frame(current_result$ypred_cv_Reps_all)

    # Add the trait, rep, and model columns
    ypred_data$trait <- current_result$trait
    ypred_data$rep <- current_result$rep
    ypred_data$model <- current_result$model

    # Combine the current ypred_data with the combined_results data frame
    combined_results <- rbind(combined_results, ypred_data)
  }


  # Extract unique models, traits, and reps
  models <- unique(combined_results$model)
  traits <- unique(combined_results$trait)
  reps <- unique(combined_results$rep)

  # Initialize lists to store results
  mod_res_per_trait_per_model <- list()

  # Iterate through each model
  for (mod in models) {
    # Filter data for the current model
    dat <- combined_results[combined_results[["model"]] == mod, ]

    traits <- unique(dat$trait)
    # Initialize a list to store results per trait for the current model
    mod_res_per_trait <- list()

    # Iterate through each trait
    for (tt in traits) {
      # Filter data for the current trait and extract 'yhat'
      dattt <- dat[dat[["trait"]] == tt, ]
      reps <- unique(dattt$rep)
      datt <- as.numeric(dat[dat[["trait"]] == tt, "yhat"])

      # Convert 'yhat' to matrix format
      mat_res <- matrix(datt, ncol = length(reps))

      # Calculate statistics
      stats <- calculate_statistics(mat_res)

      # Store the statistics in the list for the current trait
      mod_res_per_trait[[tt]] <- stats
    }

    # Store the trait results list in the list for the current model
    mod_res_per_trait_per_model[[mod]] <- mod_res_per_trait
  }

  #return(mod_res_per_trait_per_model)
  plot_reps_list <- lapply(traits, function(trait) {

    plot_all_models <- list()
    for (mod in models) {

      sik <-  predicted_vs_observed_ranking_plot(
        observed_value = as.numeric(pheno_data[[trait]]),
        predicted_value = as.numeric(mod_res_per_trait_per_model[[mod]][[trait]]$pred_mean),
        abs_very_close_threshold = abs_very_close_threshold,
        abs_close_threshold =abs_close_threshold
        )


      plot_all_models[[mod]] <- sik
    }


    list(traits = trait,
         predicted_vs_observed_plots = plot_all_models
    )
  })

  names(plot_reps_list) <-  traits

  return(list(predicted_vs_observed_plots = plot_reps_list,
              mod_res_per_trait_per_model = mod_res_per_trait_per_model))

}




