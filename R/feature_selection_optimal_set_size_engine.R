
set_parallel_plan_rf <- function(replication = 1,
                                 num_cores = NULL,
                                 sys_name) {
  # Define the plan based on the system
  plan_type <- ifelse(sys_name == "Windows", "multisession", "multicore")

  # Check if parallel execution is beneficial
  if (replication > 1) {
    if(is.null(num_cores)){

      num_cores <-  parallel::detectCores()
      num_cores <- num_cores*0.7
    }
    future::plan(plan_type, workers = num_cores)
  } else {
    future::plan("sequential")
  }
}
# Function to perform CV with replication, aggregate feature importance, and stability selection
cv_with_replication_stability <- function(X, y, n_folds = 5,
                                          replication = 10,
                                          n_tree = 500,
                                          threshold = 0.8,
                                          num_cores = NULL) {
  set.seed(123)  # For reproducibility in CV

  # Initialize matrix to store selection frequency for each marker
  selection_matrix <- matrix(0, nrow = ncol(X), ncol = n_folds * replication)
  rownames(selection_matrix) <- colnames(X)

  # Perform k-fold cross-validation
  cv_folds <- caret::createFolds(y, k = n_folds, returnTrain = TRUE)

  sys_name <- Sys.info()["sysname"]
  if (!is.null(num_cores) && num_cores > 1) {
    #sys_name <- Sys.info()["sysname"]
    set_parallel_plan_rf(replication = replication,num_cores = num_cores,
                         sys_name = sys_name)
  } else {
    # Automatically determine the number of cores and use half of them
    detected_cores <- parallel::detectCores(logical = TRUE)
    # For non-Windows systems, consider physical cores only
    num_cores <- round(detected_cores * 0.5)

    set_parallel_plan_rf(replication = replication,
                         num_cores = num_cores,
                         sys_name = sys_name)
  }


  fold_index <- 1
  for (fold in cv_folds) {
    train_X <- X[fold, ]
    train_y <- y[fold]

    # Replicate the Random Forest model within each fold in parallel
    results <- future.apply::future_lapply(seq_len(replication), function(rep) {
      #results <- future_map(1:replication, function(rep) {
      set.seed(rep)  # Change seed for each replication
      rf_model <- randomForest::randomForest(train_X, train_y, ntree = n_tree, importance = TRUE)
      importance_scores <- randomForest::importance(rf_model, type = 1)  # Mean Decrease Accuracy
      return(importance_scores)
    }, future.seed = TRUE)

    future::plan("sequential")

    # Aggregate selection frequencies
    for (rep in 1:replication) {
      selection_matrix[, fold_index] <- ifelse(results[[rep]] > 0, 1, 0)
      fold_index <- fold_index + 1
    }
  }

  # Calculate the selection frequency for each marker
  selection_frequency <- rowMeans(selection_matrix)

  # Rank features based on selection frequency
  feature_ranking <- data.frame(
    Marker = names(selection_frequency),
    SelectionFrequency = selection_frequency
  ) |>
    dplyr::arrange(dplyr::desc(SelectionFrequency))

  # Filter markers based on the threshold
  selected_markers <- feature_ranking |>
    dplyr::filter(SelectionFrequency >= threshold)

  return(list(
    SelectedMarkers = selected_markers,
    AllMarkers = feature_ranking
  ))
}

sequential_feature_selection <- function(sorted_markers, start_marker_count = 10,
                                         X, y, n_tree = 500, num_cores = 10){
  # Define the range of marker subsets to evaluate
  marker_subsets <- seq(start_marker_count, length(sorted_markers), by = start_marker_count)  # Adjust step size as needed


  sys_name <- Sys.info()["sysname"]
  if (!is.null(num_cores) && num_cores > 1) {
    #sys_name <- Sys.info()["sysname"]
    set_parallel_plan_rf(replication = length(marker_subsets),num_cores = num_cores,
                         sys_name = sys_name)
  } else {
    # Automatically determine the number of cores and use half of them
    detected_cores <- parallel::detectCores(logical = TRUE)
    # For non-Windows systems, consider physical cores only
    num_cores <- round(detected_cores * 0.5)

    set_parallel_plan_rf(replication = length(marker_subsets),
                         num_cores = num_cores,
                         sys_name = sys_name)
  }
  # Use future_lapply to parallelize the loop
  results <- future.apply::future_lapply(marker_subsets, function(n) {
    subset_X <- X[, sorted_markers[1:n]]
    model <- randomForest::randomForest(subset_X, y, ntree = n_tree)
    pred <- stats::predict(model, X)
    rmse <- sqrt(mean((y - pred)^2))
    return(data.frame(NumMarkers = n, RMSE = rmse))
  }, future.seed = TRUE)

  future::plan("sequential")
  # Combine the results into a single data frame
  results <- do.call(rbind, results)

  results <- results[order(results$RMSE, decreasing = FALSE), ]

  return(results)

}


feature_selection_rf <- function(X, y, n_folds = 5, replication = 10,
                                 n_tree = 500, threshold = 0.8,
                                 num_cores = NULL,
                                 all_marker_set = TRUE,
                                 start_marker_count = 10,
                                 optimal_marker_set = TRUE){

  cv_results_rf <- cv_with_replication_stability (X, y, n_folds = 5, replication = 10,
                                                  n_tree = 500, threshold = 0.8,
                                                  num_cores = 10)

  if(isTRUE(all_marker_set) && isFALSE(optimal_marker_set)){

    all_markers <- cv_results_rf$AllMarkers[, -2, drop = FALSE]
    rownames(all_markers) <- NULL
    return(all_markers)

  } else {
    if(isTRUE(optimal_marker_set)){
      sorted_markers <- cv_results_rf$AllMarkers$Marker

      optimal_markers_size <- sequential_feature_selection(sorted_markers = sorted_markers, start_marker_count= start_marker_count,
                                                           X = X, y = y, n_tree = n_tree, num_cores = num_cores)
    }

    return(optimal_markers_size)
  }

}
