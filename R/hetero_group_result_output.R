# # Initialize an empty list to store data frames for each trait
# trait_dfs <- list()
#
# # Loop through the results to aggregate by trait
# for (i in seq_along(TT)) {
#   result <- TT[[i]]
#   trait <- result$trait
#   # Check if the trait already has a data frame in the list
#   if (!trait %in% names(trait_dfs)) {
#     # If not, create a new data frame for this trait
#     trait_dfs[[trait]] <- result$eval_metrics_reps
#   } else {
#     # If yes, bind the new rows to the existing data frame
#     trait_dfs[[trait]] <- rbind(trait_dfs[[trait]], result$eval_metrics_reps)
#   }
# }
#
# # Now `trait_dfs` contains a data frame for each trait with all replications
#
# # Initialize an empty list to store summary data frames for each trait
# summary_dfs <- list()
#
# # Loop through the trait data frames to calculate means
# for (trait in names(trait_dfs)) {
#   df <- trait_dfs[[trait]]
#   df <- df[, !(names(df) %in% c("Rep"))]
#   # Calculate mean of each metric, grouping by Env if necessary
#   summary_df <- aggregate(. ~ Env, data = df, FUN = mean)
#   # Add the summary data frame to the list
#   summary_dfs[[trait]] <- summary_df
# }
#
# # `summary_dfs` now contains a summary data frame for each trait with mean values
