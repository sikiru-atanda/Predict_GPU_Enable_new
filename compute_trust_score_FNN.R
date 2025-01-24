compute_trust_score <- function(
    X_train,          # n x p matrix of training features
    y_train,          # vector of training labels
    y_pred_train,     # vector of training predictions
    X_test,           # m x p matrix of test features
    k = 5,            # number of nearest neighbors
    correctness_tol = NULL, 
    reduce_dim = TRUE,
    n_components = 10  # number of PCs to keep
) {
  library(FNN)  # for get.knnx
  
  # 1. Determine "correctness" in training data
  if (is.null(correctness_tol)) {
    correctness_tol <- 0.5 * sd(y_train, na.rm = TRUE)
  }
  errors <- abs(y_pred_train - y_train)
  correctness <- (errors <= correctness_tol)
  
  # 2. Subset to "correct" training points
  X_train_correct <- X_train[correctness, , drop = FALSE]
  if (nrow(X_train_correct) < k) {
    warning("Fewer than k correct training points. Adjusting k.")
    k <- nrow(X_train_correct)
  }
  
  # 3. OPTIONAL: Dimensionality Reduction
  if (reduce_dim) {
    # Fit PCA on training set
    pca_model <- prcomp(X_train_correct, center = TRUE, scale. = TRUE)
    
    # Keep only the first n_components principal components
    # (Ensure n_components <= ncol(X_train_correct))
    n_components <- min(n_components, ncol(X_train_correct))
    
    X_train_correct_reduced <- pca_model$x[, 1:n_components, drop = FALSE]
    
    # Project X_test onto same PC loadings
    # Must apply same centering/scaling used in PCA
    X_test_scaled <- scale(X_test, center = pca_model$center, scale = pca_model$scale)
    X_test_reduced <- as.matrix(X_test_scaled) %*% pca_model$rotation[, 1:n_components]
    
    # 4. Find k-NN in the reduced space
    nn_results <- get.knnx(X_train_correct_reduced, X_test_reduced, k = k)
  } else {
    # No dimensionality reduction
    nn_results <- get.knnx(X_train_correct, X_test, k = k)
  }
  
  # 5. Compute Trust Score
  trust_scores <- apply(nn_results$nn.dist, 1, function(dists) {
    mean_dist <- mean(dists)
    1 / (1 + mean_dist)
  })
  
  return(trust_scores)
}


trust_scores <- compute_trust_score(X_train = SNP_train, 
                                    y_train = y_train_geno,
                                    y_pred_train = geno_preds_train, 
                                    X_test = SNP_test, k = 5)
