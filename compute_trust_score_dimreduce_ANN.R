compute_trust_score_dimreduce_ANN <- function(
    X_train, 
    y_train,
    y_pred_train,
    X_test,
    k = 5,
    correctness_tol = NULL,
    reduce_dim = TRUE,
    n_components = 10
) {
  library(RANN)
  
  # 1. "correctness"
  if (is.null(correctness_tol)) {
    correctness_tol <- 0.5 * sd(y_train, na.rm = TRUE)
  }
  errors <- abs(y_pred_train - y_train)
  correctness <- (errors <= correctness_tol)
  
  X_train_correct <- X_train[correctness, , drop = FALSE]
  if (nrow(X_train_correct) < k) {
    warning("Fewer than k correct training points. Adjusting k.")
    k <- nrow(X_train_correct)
  }
  
  # 2. Dimensionality Reduction
  if (reduce_dim) {
    pca_model <- prcomp(X_train_correct, center = TRUE, scale. = TRUE)
    n_components <- min(n_components, ncol(X_train_correct))
    
    X_train_correct_reduced <- pca_model$x[, 1:n_components, drop = FALSE]
    
    X_test_scaled <- scale(X_test, center = pca_model$center, scale = pca_model$scale)
    X_test_reduced <- as.matrix(X_test_scaled) %*% pca_model$rotation[, 1:n_components]
    
    # Approximate NN in reduced space
    ann_results <- nn2(
      data  = X_train_correct_reduced,
      query = X_test_reduced,
      k     = k
    )
    
  } else {
    # Approximate NN in original space
    ann_results <- nn2(data = X_train_correct, query = X_test, k = k)
  }
  
  # 3. Trust Scores
  trust_scores <- apply(ann_results$nn.dists, 1, function(dists) {
    mean_dist <- mean(dists)
    1 / (1 + mean_dist)
  })
  
  return(trust_scores)
}

trust_scores_ANN <- compute_trust_score_dimreduce_ANN(X_train = SNP_train, 
                                  y_train = y_train_geno,
                                  y_pred_train = geno_preds_train, 
                                  X_test = SNP_test)
