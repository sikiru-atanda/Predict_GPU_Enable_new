compute_mahalanobis <- function(
    train_data,
    test_data,
    dist_threshold = NULL,
    p_threshold = 1000,    # dimension limit above which we do PCA
    n_components = 50      # how many principal components to keep if p > p_threshold
) {
  # train_data: matrix/data.frame of training features, size (n x p)
  # test_data:  matrix/data.frame of test features, size (m x p)
  # dist_threshold: numeric or NULL. If NULL, auto-set from training distribution
  # p_threshold:    integer, dimension cutoff for deciding if we do PCA
  # n_components:   how many PCs to keep if we do PCA (or SVD)
  
  # Convert to matrix if data.frames were passed
  train_data <- as.matrix(train_data)
  test_data  <- as.matrix(test_data)
  
  n_train <- nrow(train_data)
  p       <- ncol(train_data)
  
  #-------------------------------------------------------
  # 0. Optional check for dimension
  # If p > p_threshold, do PCA on train_data & test_data
  #-------------------------------------------------------
  if (p > p_threshold) {
    message(sprintf(
      "p = %d > %d, applying PCA to reduce dimension to %d PCs before Mahalanobis.",
      p, p_threshold, n_components
    ))
    
    # 0a. Center & scale train_data for PCA
    train_mean <- colMeans(train_data)
    train_sd   <- apply(train_data, 2, sd)
    
    # Avoid division by zero if any column is constant
    train_sd[train_sd < 1e-15] <- 1e-15
    
    train_data_scaled <- scale(train_data, center = train_mean, scale = train_sd)
    test_data_scaled  <- scale(test_data,  center = train_mean, scale = train_sd)
    
    # 0b. PCA on training data
    pca_result <- prcomp(train_data_scaled, center = FALSE, scale. = FALSE)
    
    # Keep the top n_components
    d_use <- min(n_components, ncol(train_data_scaled), nrow(train_data_scaled) - 1)
    
    # Project train and test onto these PCs
    train_data_reduced <- pca_result$x[, 1:d_use, drop = FALSE]  # n x d_use
    # For test data, multiply by the rotation matrix
    test_data_reduced <- test_data_scaled %*% pca_result$rotation[, 1:d_use, drop=FALSE]
    
  } else {
    # If p <= p_threshold, just proceed with original data
    train_data_reduced <- train_data
    test_data_reduced  <- test_data
  }
  
  #-------------------------------------------------------
  # 1. Center the *reduced* training data
  #-------------------------------------------------------
  mu <- colMeans(train_data_reduced)
  centered_train <- sweep(train_data_reduced, 2, mu, "-")
  
  #-------------------------------------------------------
  # 2. Compute covariance matrix in reduced space
  #-------------------------------------------------------
  cov_mat <- cov(centered_train)
  
  #-------------------------------------------------------
  # 3. Invert the covariance matrix
  #-------------------------------------------------------
  cov_inv <- solve(cov_mat)
  
  #-------------------------------------------------------
  # 4. Center the *reduced* test data
  #-------------------------------------------------------
  centered_test <- sweep(test_data_reduced, 2, mu, "-")
  
  #-------------------------------------------------------
  # 5. Compute Mahalanobis distances for each test instance
  #-------------------------------------------------------
  dists <- apply(centered_test, 1, function(x) {
    as.numeric(t(x) %*% cov_inv %*% x)
  })
  
  #-------------------------------------------------------
  # 6. Determine threshold if needed
  #-------------------------------------------------------
  if (is.null(dist_threshold)) {
    # For illustration, choose 95th percentile of training distribution distances
    train_dists <- apply(centered_train, 1, function(x) {
      as.numeric(t(x) %*% cov_inv %*% x)
    })
    dist_threshold <- quantile(train_dists, 0.95)
  }
  
  #-------------------------------------------------------
  # 7. Classify as in-distribution (0) or out-of-distribution (1)
  #-------------------------------------------------------
  ood_flags <- as.numeric(dists > dist_threshold)
  
  #-------------------------------------------------------
  # Return results
  #-------------------------------------------------------
  return(list(
    distances    = dists,
    threshold    = dist_threshold,
    ood_flags    = ood_flags,  # 1 = OOD, 0 = in-distribution
    # For reference, you could also return the rotation matrix or PCA object if needed
    pca_used     = (p > p_threshold)
  ))
}
