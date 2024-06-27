# Load necessary libraries
# library(caret)
# library(e1071)
# library(transport)
# library(stats)

calculate_thresholds <- function(q25, q75, iqr_multiplier= 1.5) {
  iqr <- q75 - q25
  high_threshold <- max(0, q25 - iqr_multiplier * (iqr / 2))
  moderate_threshold <- max(0, q75 - iqr_multiplier * (iqr / 2))
  return(list(high = high_threshold, moderate = moderate_threshold))
}

perform_pca <- function(train_data, test_data, n_components = 10, title = "PCA Plot") {

               pca_train <- stats::prcomp(train_data, center = FALSE, scale. = TRUE)
               pca_test <- stats::prcomp(test_data, center = FALSE, scale. = TRUE)

              pca_train_df <- data.frame(pca_train$x[, 1:n_components], Set = 'Train')
              pca_test_df <- data.frame(pca_test$x[, 1:n_components], Set = 'Test')
              pca_combined <- rbind(pca_train_df, pca_test_df)

              ggplot2::ggplot(pca_combined, ggplot2::aes(x = PC1, y = PC2, color = Set)) +
                ggplot2::geom_point() +
                ggplot2::ggtitle(title) +
                ggplot2::theme_minimal()

}

#pca_plot_snp <- perform_pca(train_data = as.matrix(trn_geno), test_data = as.matrix(tst_geno))
# pca_plot_snp <- perform_pca(train_data, test_data, "PCA Plot for SNP Data")
#
# print(pca_plot_snp)
# Function to apply PCA if the number of features exceeds the threshold
apply_pca_if_needed <- function(data, n_components = 50, threshold = 100, apply_pca = FALSE) {

  if (ncol(data) > threshold | isTRUE(apply_pca)) {

    if(!is.null(n_components)){
      if(n_components>nrow(data)){
        n_components <- nrow(data)*0.4
      }
    } else{
      n_components <- nrow(data)*0.4
    }

    pca_model <- stats::prcomp(data, center = TRUE, scale. = FALSE)
    pca_data <- data.frame(pca_model$x[, 1:n_components])
    return(pca_data)
  } else {
    return(data)
  }
}

# Function to calculate Mahalanobis distance
calculate_mahalanobis <- function(x, mean, inv_cov_matrix) {
  x_minus_mean <- x - mean
  sqrt(t(x_minus_mean) %*% inv_cov_matrix %*% x_minus_mean)
}

mahalanobis_distances_testSet <- function(geno_trn = NULL,
                                          geno_tst = NULL,
                                          geno_tst_trn = NULL,
                                          names_tst = NULL,
                                          names_trn = NULL,
                                          n_components = 20,
                                          threshold = 100,
                                          target = "test_set",
                                          apply_pca = TRUE) {

  # Return NA if geno_tst_trn, geno_trn, and geno_tst are all NULL
  if (is.null(geno_tst_trn) && is.null(geno_trn) && is.null(geno_tst)) {
    return(NA)
  }


  if (!is.null(geno_trn) && !is.null(geno_tst)) {
    geno_trn <- apply_pca_if_needed(geno_trn, n_components, threshold, apply_pca)
    geno_tst <- apply_pca_if_needed(geno_tst, n_components, threshold, apply_pca)
  }


  if (!is.null(geno_tst_trn)) {
    geno_tst_trn <- apply_pca_if_needed(geno_tst_trn, n_components, threshold, apply_pca)

    if (!is.null(names_trn) && !is.null(names_tst)) {
      train_indices <- match(names_trn, rownames(geno_tst_trn))
      geno_trn <- geno_tst_trn[train_indices, ]
      geno_tst <- geno_tst_trn[-train_indices, ]
    } else if (!is.null(names_trn) && is.null(names_tst)) {
      train_indices <- match(names_trn, rownames(geno_tst_trn))
      geno_trn <- geno_tst_trn[train_indices, ]
      geno_tst <- geno_tst_trn[-train_indices, ]
    } else if (is.null(names_trn) && !is.null(names_tst)) {
      test_indices <- match(names_tst, rownames(geno_tst_trn))
      geno_trn <- geno_tst_trn[-test_indices, ]
      geno_tst <- geno_tst_trn[test_indices, ]
    }
  }



  calculate_thresholds <- function(q25, q75, iqr_multiplier) {
    iqr <- q75 - q25
    return(list(low = q25 - iqr_multiplier * iqr, high = q75 + iqr_multiplier * iqr))
  }

  if (target == "test_set" && !is.null(geno_tst)) {
    mean_vector <- colMeans(geno_trn)
    cov_matrix <- corpcor::cov.shrink(geno_trn, verbose = FALSE)
    #
    if (inherits(cov_matrix, "shrinkage")) {
      inv_cov_matrix <- solve(cov_matrix)
    } else {
      inv_cov_matrix <- MASS::ginv(cov_matrix)
    }

    #training_distances <- apply(geno_trn, 1, function(x) calculate_mahalanobis(x, mean_vector, inv_cov_matrix))
    testing_distances <- apply(geno_tst, 1, function(x) calculate_mahalanobis(x, mean_vector, inv_cov_matrix))

    q75 <- stats::quantile(testing_distances, 0.75, na.rm = TRUE)
    q25 <- stats::quantile(testing_distances, 0.25, na.rm = TRUE)
    thresholds <- calculate_thresholds(q25, q75, 1.5)

    testing_results <- ifelse(testing_distances < thresholds$high, 1,
                              ifelse(testing_distances == thresholds$high, 0.5, 0))
    return(testing_results)

  } else if (target == "training_set" && !is.null(geno_trn)) {
    mean_vector <- colMeans(geno_tst)
    cov_matrix <- corpcor::cov.shrink(geno_tst, verbose = FALSE)
    #
    if (inherits(cov_matrix, "shrinkage")) {
      inv_cov_matrix <- solve(cov_matrix)
    } else {
      inv_cov_matrix <- MASS::ginv(cov_matrix)
    }

    training_distances <- apply(geno_trn, 1, function(x) calculate_mahalanobis(x, mean_vector, inv_cov_matrix))
    #testing_distances <- apply(geno_tst, 1, function(x) calculate_mahalanobis(x, mean_vector, inv_cov_matrix))

    q75 <- stats::quantile(training_distances, 0.75, na.rm = TRUE)
    q25 <- stats::quantile(training_distances, 0.25, na.rm = TRUE)
    thresholds <- calculate_thresholds(q25, q75, 1.5)

    testing_results <- ifelse(training_distances < thresholds$high, 1,
                              ifelse(training_distances == thresholds$high, 0.5, 0))
    return(testing_results)
  }

  return(NA) # Return NA if conditions are not met
}


###Kolmogorov-Smirnov Test Similarity
# Function to calculate similarity score using KS test with PCA if needed
# ks_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data <- apply_pca_if_needed(test_data, n_components, threshold)
#
#   scores <- sapply(1:ncol(train_data), function(i) {
#     stats::ks.test(train_data[, i], test_data[, i])$statistic
#   })
#   return(mean(scores))
# }
#
# # Function to calculate similarity score for each individual in the test set
# individual_ks_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data <- apply_pca_if_needed(test_data, n_components, threshold)
#   apply(test_data, 1, function(test_individual) {
#     scores <- sapply(1:ncol(train_data), function(i) {
#       stats::ks.test(train_data[, i], rep(test_individual[i], nrow(train_data)))$statistic
#     })
#     return(mean(scores))
#   })
# }
#
# # ind_similarity_snp_ks <- individual_ks_similarity_score(train_data = trn_geno,
# #                                                  test_data = tst_geno)
# #
# # similarity_snp_ks <- ks_similarity_score(train_data = trn_geno,
# #                                              test_data = tst_geno)
# ###
# # Function to calculate similarity score using Wasserstein distance with PCA if needed
# wasserstein_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data_pca <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data_pca <- apply_pca_if_needed(test_data, n_components, threshold)
#
#   scores <- sapply(1:ncol(train_data_pca), function(i) {
#     transport::wasserstein1d(train_data_pca[, i], test_data_pca[, i])
#   })
#   return(mean(scores))
# }
#
# # Function to calculate individual Wasserstein similarity scores with PCA if needed
# individual_wasserstein_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data_pca <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data_pca <- apply_pca_if_needed(test_data, n_components, threshold)
#
#   similarity_scores <- apply(test_data_pca, 1, function(test_individual) {
#     scores <- sapply(1:ncol(train_data_pca), function(i) {
#       transport::wasserstein1d(train_data_pca[, i], rep(test_individual[i], nrow(train_data_pca)))
#     })
#     return(mean(scores))
#   })
#   return(similarity_scores)
# }
#
# # similarity_snp_ws <- wasserstein_similarity_score(train_data = trn_geno,
# #                                                       test_data = tst_geno)
# #
# # ind_similarity_snp_ws <- individual_wasserstein_similarity_score(train_data = trn_geno,
# #                                                                  test_data = tst_geno)
# ######
# # Function to calculate similarity using a classifier with PCA if needed
# classifier_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data_pca <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data_pca <- apply_pca_if_needed(test_data, n_components, threshold)
#
#   # Label the data
#   train_labels <- rep(0, nrow(train_data_pca))
#   test_labels <- rep(1, nrow(test_data_pca))
#
#   # Combine the data
#   combined_data <- rbind(train_data_pca, test_data_pca)
#   combined_labels <- c(train_labels, test_labels)
#
#   # Create a dataframe
#   combined_df <- data.frame(combined_data, Label = factor(combined_labels))
#
#   # Split into training and validation sets
#   set.seed(123)
#   trainIndex <- caret::createDataPartition(combined_df$Label, p = .7, list = FALSE)
#   train_df <- combined_df[trainIndex,]
#   test_df <- combined_df[-trainIndex,]
#
#   # Train a classifier (e.g., random forest)
#   # model <- caret::train(Label ~ ., data = train_df, method = "rf", trControl = trainControl(method = "none"))
#   #
#   # # Predict the probability of being in the test set
#   # prediction <- predict(model, newdata = test_df, type = "prob")[, 2]
#   #predictions <- stats::predict(model, newdata = test_df)
#   #similarity_score <- 1 -prediction
#   # Train a classifier (e.g., logistic regression)
#   model <- caret::train(Label ~ ., data = train_df, method = "glm", family = "binomial")
#
#   # Predict on the validation set
#   predictions <- stats::predict(model, newdata = test_df)
#
#   #Calculate accuracy
#   confusion <- caret::confusionMatrix(predictions, test_df$Label)
#   accuracy <- confusion$overall['Accuracy']
#
#   # Return the similarity score (1 - accuracy)
#   similarity_score <- 1 - as.numeric(accuracy)
#
#   return(similarity_score)
# }
#
# # Function to calculate individual classifier-based similarity scores with PCA if needed
# individual_classifier_similarity_score <- function(train_data, test_data, threshold = 100, n_components = 10) {
#   train_data_pca <- apply_pca_if_needed(train_data, n_components, threshold)
#   test_data_pca <- apply_pca_if_needed(test_data, n_components, threshold)
#
#   # Ensure column names are consistent
#   colnames(train_data_pca) <- paste0("PC", 1:ncol(train_data_pca))
#   colnames(test_data_pca) <- paste0("PC", 1:ncol(test_data_pca))
#
#   similarity_scores <- sapply(1:nrow(test_data_pca), function(i) {
#     test_individual <- test_data_pca[i, , drop = FALSE]
#
#     # Label the data
#     train_labels <- rep(0, nrow(train_data_pca))
#     test_labels <- 1
#
#     # Combine the data
#     combined_data <- rbind(train_data_pca, test_individual)
#     combined_labels <- c(train_labels, test_labels)
#
#     # Create a dataframe
#     combined_df <- data.frame(combined_data, Label = factor(combined_labels))
#
#     # Train a classifier (e.g., random forest)
#     model <- caret::train(Label ~ ., data = combined_df, method = "rf", trControl = trainControl(method = "none"))
#
#     # Predict the probability of being in the test set
#     prediction <- stats::predict(model, newdata = test_individual, type = "prob")[, 2]
#
#     # # Train a classifier (e.g., logistic regression)
#     # model <- train(Label ~ ., data = combined_df, method = "glm", family = "binomial", trControl = trainControl(method = "none"))
#     #
#     # # Predict the probability of being in the test set
#     # prediction <- predict(model, newdata = test_individual, type = "prob")[, 2]
#     #
#
#     # Similarity score is the inverse probability of being in the test set
#     similarity_score <- 1 - prediction
#     return(similarity_score)
#   })
#   return(similarity_scores)
# }
#

# similarity_snp_cl <- classifier_similarity_score(train_data = trn_geno,
#                                                      test_data = tst_geno)
#
# ind_similarity_snp_cl <- individual_classifier_similarity_score(train_data = trn_geno,
#                                                                 test_data = tst_geno)
# #
# composite_reliability <- similarity / rel_width_interval

