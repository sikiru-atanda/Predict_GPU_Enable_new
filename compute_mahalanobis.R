compute_mahalanobis <- function(train_data, test_data, dist_threshold = NULL) {
  # train_data: matrix of training features (or a GRM if that's how you store it)
  # test_data: matrix of test features
  # dist_threshold: distance above which we consider a test point out-of-distribution
  
  # 1. Center the training data
  mu <- colMeans(train_data)
  centered_train <- sweep(train_data, 2, mu, "-")
  
  # 2. Compute covariance matrix
  cov_mat <- cov(centered_train)
  
  # 3. Invert the covariance matrix
  cov_inv <- solve(cov_mat)
  
  # 4. Center the test data
  centered_test <- sweep(test_data, 2, mu, "-")
  
  # 5. Compute Mahalanobis distances for each test instance
  dists <- apply(centered_test, 1, function(x) {
    as.numeric(t(x) %*% cov_inv %*% x)
  })
  
  # If no threshold is given, you could pick one (e.g., based on a quantile of training distances)
  # Or directly return the distances for further use.
  if (is.null(dist_threshold)) {
    # For illustration, choose 95th percentile of training distribution distances as threshold
    train_dists <- apply(centered_train, 1, function(x) as.numeric(t(x) %*% cov_inv %*% x))
    dist_threshold <- quantile(train_dists, 0.95)
  }
  
  # Classify as in-distribution (0) or out-of-distribution (1)
  ood_flags <- as.numeric(dists > dist_threshold)
  
  return(list(
    distances    = dists,
    threshold    = dist_threshold,
    ood_flags    = ood_flags  # 1 = OOD, 0 = in-distribution
  ))
}

mahal_result <- compute_mahalanobis(X_train, X_test, dist_threshold = NULL)
### Example usage:
# mahal_result <- compute_mahalanobis(X_train, X_test, dist_threshold = NULL)
# mahal_distances  <- mahal_result$distances
# mahal_threshold  <- mahal_result$threshold
# mahal_ood_flags  <- mahal_result$ood_flags
