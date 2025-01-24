rm(list = ls()); ls()
gc()
setwd("D:/FILES PILE UP/2017_2FoldsCV")

library(data.table)


#setwd("/home/asa259/2017_2018AVG")

#setwd("F:/FILES PILE UP/2017_2018_Data")

#setwd("D:/PredictProR")
ref = read.csv("bluesGY_GSOptimalCombineTrial.csv", header = T, as.is = T, sep = ",", stringsAsFactors = F)


ref = subset(ref, select = c("GID", "GY_BLUE", "PopClass"))

ref$GID = as.character(ref$GID)


table(ref$PopClass)

#setwd("D:/FILES PILE UP/2017_2018_Data")
mark = read.csv("2017Kenya_Dan0.05.csv", header = T, as.is = T,
                stringsAsFactors = F, row.names = 1)

mark[mark == 1] <- 2


VanRaden <- function(geno_clean){
  freq <- colMeans(geno_clean) / 2
  geno_clean <- scale(geno_clean, center=T, scale=F)
  return(tcrossprod(geno_clean) / sum(2*freq*(1-freq)))
}


pops <- unique(ref[["PopClass"]])
tst_pops <- c("pop1", "pop12", "pop8", "pop4")

trn_pops <- pops[!pops%in%tst_pops]


tst_geno_pop12 = mark[rownames(mark)%in%ref[ref$PopClass%in%tst_pops[2], "GID"], ]
tst_y = ref[ref$PopClass%in%tst_pops[2], "GY_BLUE"]
#GRM_pop2 <- VanRaden(tst_geno_pop12)
trn_geno = mark[!rownames(mark)%in%ref[ref$PopClass%in%tst_pops, "GID"], ]
#GRM_trn <- VanRaden(trn_geno)

trn_y = ref[!ref$PopClass%in%tst_pops, "GY_BLUE"]


####################################
library(GA)
library(Matrix)

# Define the fitness function
fitness_function <- function(indices, relationship_matrix, train_indices, test_indices, lambda = 0.5) {
  selected_train <- relationship_matrix[train_indices[indices], train_indices[indices], drop = FALSE]
  selected_test <- relationship_matrix[train_indices[indices], test_indices, drop = FALSE]

  # Calculate the mean relationship between the selected training set and the test set
  mean_relationship_to_test <- mean(selected_test)

  # Calculate the mean relationship within the selected training set
  mean_relationship_within_train <- mean(selected_train)

  # Objective: maximize mean relationship to test set and minimize within training set
  fitness_value <- lambda * mean_relationship_to_test - (1 - lambda) * mean_relationship_within_train
  return(fitness_value)
}

# Function to optimize the training set selection
optimize_training_set <- function(rel_matrix, trn_indices, tst_indices, pop_size = 50, max_iter = 100, lambda = 0.5) {
  n_train <- length(trn_indices)

  if (ncol(rel_matrix) != nrow(rel_matrix)) {
    stop("The relationship matrix must be square.")
  }

  if (any(trn_indices > nrow(rel_matrix)) || any(tst_indices > nrow(rel_matrix))) {
    stop("Train or test indices are out of bounds.")
  }

  # Run the genetic algorithm
  ga_result <- GA::ga(
    type = "binary",
    fitness = function(indices) {
      selected_indices <- which(indices == 1)
      if (length(selected_indices) == 0) return(-Inf)
      fitness_function(selected_indices, rel_matrix, trn_indices, tst_indices, lambda)
    },
    nBits = n_train,
    popSize = pop_size,
    maxiter = max_iter,
    run = max_iter / 10
    #elitism = 0.1 * pop_size # Preserve the top 10% individuals
  )

  if (is.null(ga_result@solution)) {
    stop("Genetic algorithm did not find a solution.")
  }

  # Extract the best solution
  best_solution <- trn_indices[which(ga_result@solution == 1)]
  best_solution <- best_solution[!is.na(best_solution)]
  return(best_solution)
}
#####

GRM_combined <- VanRaden(mark)
train_indices <- which(rownames(GRM_combined)%in%ref[ref$PopClass%in%trn_pops, "GID"])
test_indices <- which(rownames(GRM_combined)%in%ref[ref$PopClass%in%tst_pops[2], "GID"])

# # Combine the relationship matrices
# GRM_combined <- rbind(cbind(GRM_trn, matrix(0, nrow(GRM_trn), ncol(GRM_pop2))),
#                       cbind(matrix(0, nrow(GRM_pop2), ncol(GRM_trn)), GRM_pop2))
#
# # Define indices for training and test sets
# train_indices <- 1:nrow(GRM_trn)
# test_indices <- (nrow(GRM_trn) + 1):nrow(GRM_combined)

# Optimize the training set selection

best_training_set_indices <- optimize_training_set(rel_matrix = GRM_combined,
                                                   trn_indices = train_indices,
                                                   tst_indices = test_indices,
                                                   pop_size = 200,
                                                   max_iter = 100, lambda = 0.5)



# Inspect the results
print(best_training_set_indices)

optimized_trn = rownames(GRM_combined)[best_training_set_indices]

length(optimized_trn)
length(which(optimized_trn2%in%optimized_trn))

which(rownames(tst_geno_pop12)%in%optimized_trn)


optimized_trn = rownames(GRM_combined)[best_training_set_indices]


opt_trn_geno = mark[rownames(mark)%in%optimized_trn,  ]

opt_trn_y = ref[ref$GID%in%optimized_trn, "GY_BLUE"]


historical_data <- cbind(y=opt_trn_y, opt_trn_geno)

# Initialize training set with a small subset of the historical data
initial_indices <- sample(1:nrow(historical_data), size = 50)  # Start with 50 samples
train_indices <- initial_indices
remaining_indices <- setdiff(1:nrow(historical_data), train_indices)

# Set the expected reliability threshold and maximum number of iterations
expected_reliability_threshold <- 0.80  # Define your reliability threshold
max_iterations <- 20  # Maximum number of iterations
max_train_size <- 200  # Maximum size of the training set

train_predict_xgboost <- function(data_label_geno, indices, test_geno, params, nrounds) {
  train_data <- data_label_geno[, -1]
  y_train <- data_label_geno[, 1]
  # Subset the data
  dtrain_boot <- xgboost::xgb.DMatrix(data = as.matrix(train_data[indices,]), label = y_train[indices])

  # Train the model
  model <- xgboost::xgboost(data = dtrain_boot, params = params, nrounds = nrounds, verbose = 0)

  # Predict on the original data

  pred <- stats::predict(model, as.matrix(test_geno), reshape = TRUE)


  return(pred)
}

# Get the best-tuned parameters
params <- list(booster = "dart",
               objective = "reg:squarederror",
               sample_type = "uniform",
               normalize_type = "tree",
               rate_drop = 0.1,
               skip_drop = 0.5,
               eta = 0.1,
               max_depth = 6,
               min_child_weight = 1,
               subsample = 0.8,
               colsample_bytree = 0.8,
               lambda = 1,
               alpha = 0)
# Function to train models and estimate uncertainty
train_and_estimate_uncertainty <- function(train_data, test_data,
                                           params = NULL,
                                           nrounds = 100, verbose = 0) {
  # Train bagging model

  boot_results <- boot::boot(
    data = train_data,
    statistic = train_predict_xgboost,
    R = 100,  # Number of bootstrap samples
    #sim = "ordinary",
    #dtrain = geno_omic_object_train,
    test_geno = test_data,
    params = params,
    nrounds = nrounds
  )

  pred_variances <- apply(boot_results$t, 2, var)

  pred_SE <- apply(boot_results$t, 2, sd)

  predmean <- apply(boot_results$t, 2, mean)

  lower_bound <- apply(boot_results$t, 2, quantile, probs = 0.05)
  upper_bound <- apply(boot_results$t, 2, quantile, probs = 0.95)
  interval_width <- upper_bound - lower_bound

  # Calculate reliability (e.g., 1 - (mean interval width / mean prediction))
  #interval_width <- prediction_intervals$upper_bound - prediction_intervals$lower_bound
  reliability <- 1 - (pred_variances / var(predmean))


  return(list(
    prediction_intervals = interval_width,
    sd_predictions = pred_SE,
    reliability = mean(reliability),
    predmean = predmean,
    rel = reliability
  ))
}

# Function to evaluate the impact of adding an individual to the training set
evaluate_addition <- function(train_indices, candidate_index, historical_data, geno_tst, params, nrounds) {
  temp_train_indices <- c(train_indices, candidate_index)
  temp_train_data <- historical_data[temp_train_indices, ]

  temp_results <- train_and_estimate_uncertainty(train_data = temp_train_data,
                                                 test_data = geno_tst,
                                                 params  = params,
                                                 nrounds = nrounds,
                                                 verbose = 0)
  return(temp_results$reliability)
}

# Function to evaluate the impact of removing an individual from the training set
evaluate_removal <- function(train_indices, candidate_index, historical_data, geno_tst, params, nrounds) {
  temp_train_indices <- setdiff(train_indices, candidate_index)
  temp_train_data <- historical_data[temp_train_indices, ]
  temp_results <- train_and_estimate_uncertainty(train_data = temp_train_data,
                                                 test_data = geno_tst,
                                                 params  = params,
                                                 nrounds = nrounds,
                                                 verbose = 0)
  return(temp_results$reliability)
}

# Iterative process to dynamically update the training set
previous_uncertainty <- Inf

for (iteration in 1:max_iterations) {
  cat("Iteration:", iteration, "\n")

  train_data <- historical_data[train_indices, ]

  # Train models and estimate uncertainty
  results <- train_and_estimate_uncertainty(train_data = train_data,
                                            test_data = tst_geno_pop12,
                                            params  = params,
                                            nrounds = nrounds,
                                            verbose = 0)

  current_uncertainty <- mean(results$sd_predictions)

  if (results$reliability >= expected_reliability_threshold ||
      (previous_uncertainty - current_uncertainty) < improvement_threshold) {
    cat("Stopping criteria met. Stopping iterations.\n")
    break
  }

  previous_uncertainty <- current_uncertainty

  additional_uncertainty <- numeric(length(remaining_indices))

  # Evaluate impact of each remaining individual
  for (i in seq_along(remaining_indices)) {
    test_idx <- remaining_indices[i]
    additional_uncertainty[i] <- evaluate_addition(train_indices, test_idx, historical_data, geno_tst = tst_geno_pop12,
                                                   params, nrounds)
  }

  best_candidate_idx <- remaining_indices[which.min(additional_uncertainty)]
  train_indices <- c(train_indices, best_candidate_idx)
  remaining_indices <- setdiff(remaining_indices, best_candidate_idx)

  # Optionally, remove least informative data points
  if (length(train_indices) > length(initial_indices)) {
    removal_uncertainty <- sapply(train_indices, function(idx) evaluate_removal(train_indices, idx, historical_data, geno_tst = tst_geno_pop12,
                                                                                params, nrounds))
    least_informative_indices <- train_indices[order(removal_uncertainty, decreasing = TRUE)[1:10]]  # Remove top 10 least informative
    train_indices <- setdiff(train_indices, least_informative_indices)
  }

  if (length(remaining_indices) == 0) break
}

# Final training set
final_train_data <- historical_data[train_indices, ]

# Final model training on the selected optimal training set
final_results <- train_and_estimate_uncertainty(train_data = final_train_data,
                                                test_data = tst_geno_pop12,
                                                params  = params,
                                                nrounds = nrounds,
                                                verbose = 0)
cor(final_results$predictions, tst_y)
# Display final prediction intervals and reliability
print(final_results$prediction_intervals)
cat("Final reliability:", final_results$reliability, "\n")

if (final_results$reliability < expected_reliability_threshold) {
  cat("The historical data did not contain individuals with high reliability for the test set.\n")
  cat("The test set might need to be planted for better prediction accuracy.\n")
}
