compute_mahalanobis_grm <- function(
    G_all,                # (N x N) GRM for all individuals (train + test)
    train_indices,        # vector of row/col indices for training individuals
    test_indices,         # vector of row/col indices for test individuals
    n_components = 10,    # number of eigen components to keep
    dist_threshold = NULL # optional user-specified threshold
) {
  # 1) Eigen-decomposition (or SVD) of the full GRM
  #    - G_all should be symmetric (N x N).
  eig <- eigen(G_all, symmetric = TRUE)
  eigvals  <- eig$values
  eigvecs  <- eig$vectors
  
  # Some GRMs may have small negative eigenvalues from rounding; set them to 0
  eigvals[eigvals < 1e-12] <- 0
  
  # Are they sorted in descending order? Usually 'eigen()' does that. We'll assume yes.
  # Otherwise, reorder them with order(..., decreasing=TRUE).
  
  # 2) Choose how many components to keep
  #    Only keep positive eigenvalues
  positive_eigs <- eigvals > 1e-12
  d_max <- sum(positive_eigs)    # number of positive eigenvalues
  d_use <- min(n_components, d_max)
  if (d_use < 1) stop("No positive eigenvalues found. Cannot embed.")
  
  # 3) Build embedding for all individuals
  #    X_all = U_d * sqrt(Lambda_d), an (N x d_use) matrix
  lambda_d <- diag(sqrt(eigvals[1:d_use]), nrow = d_use, ncol = d_use)
  U_d      <- eigvecs[, 1:d_use, drop = FALSE]
  
  X_all <- U_d %*% lambda_d  # dimension: (N x d_use)
  
  # 4) Separate training vs. test rows
  X_train <- X_all[train_indices, , drop = FALSE]  # (n x d_use)
  X_test  <- X_all[test_indices,  , drop = FALSE]  # (m x d_use)
  
  # 5) Mahalanobis distance calculation
  #    a) Covariance of X_train
  mu_train <- colMeans(X_train)
  X_train_centered <- sweep(X_train, 2, mu_train, "-")
  cov_train <- cov(X_train_centered)
  
  #    b) Invert the covariance matrix
  cov_inv <- solve(cov_train)
  
  #    c) Compute distances for X_test
  #       Center X_test by mu_train
  X_test_centered <- sweep(X_test, 2, mu_train, "-")
  
  # Vectorized approach: distance_i = (x_i)^T * cov_inv * (x_i)
  # We'll do it row by row
  mahal_dists <- apply(X_test_centered, 1, function(x) {
    as.numeric(t(x) %*% cov_inv %*% x)
  })
  
  # 6) OOD threshold
  if (is.null(dist_threshold)) {
    # For illustration, define it as the 95th percentile of training distances
    # i.e. we see how training individuals deviate from the training mean.
    # Then test them against that cutoff.
    
    # Compute "self-Mahalanobis" distances on training set to define baseline
    train_mahal <- apply(X_train_centered, 1, function(x) {
      as.numeric(t(x) %*% cov_inv %*% x)
    })
    dist_threshold <- quantile(train_mahal, 0.95)
  }
  
  # 7) Flag as 1 = OOD if distance > threshold
  ood_flags <- as.numeric(mahal_dists > dist_threshold)
  
  return(list(
    distances    = mahal_dists,
    threshold    = dist_threshold,
    ood_flags    = ood_flags
  ))
}
