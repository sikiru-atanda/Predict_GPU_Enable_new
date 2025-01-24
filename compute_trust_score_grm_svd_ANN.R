compute_trust_score_grm_svd_ANN <- function(
    G_train_train,      # (n x n) GRM for training
    G_train_test,       # (n x m) cross-GRM for train vs. test
    y_train,            # vector of training labels
    y_pred_train,       # vector of predictions on training set
    k = 5,              # number of nearest neighbors
    correctness_tol = NULL,
    n_components = 5    # number of eigen-components to keep
) {
  # We'll use RANN for approximate k-NN
  library(RANN)
  
  n <- nrow(G_train_train)  # number of training individuals
  m <- ncol(G_train_test)    # number of test individuals
  
  #--------------------------
  # 1. Define correctness
  #--------------------------
  if (is.null(correctness_tol)) {
    correctness_tol <- 0.5 * sd(y_train, na.rm = TRUE)
  }
  errors <- abs(y_pred_train - y_train)
  correctness <- (errors <= correctness_tol)
  
  #--------------------------
  # 2. Eigendecomposition (or SVD) of G_train_train
  #    We assume G_train_train is symmetric and (ideally) centered.
  #--------------------------
  eig <- eigen(G_train_train, symmetric = TRUE)
  eigvals  <- eig$values
  eigvecs  <- eig$vectors
  
  # Filter out tiny or negative eigenvalues (in case GRM isn't perfectly PSD)
  eigvals[eigvals < 1e-12] <- 0
  
  # Ensure sorted in descending order (usually 'eigen()' does this)
  # But we might explicitly do:
  # order_idx <- order(eigvals, decreasing = TRUE)
  # eigvals  <- eigvals[order_idx]
  # eigvecs  <- eigvecs[, order_idx]
  
  # Decide how many components to keep
  d_use <- min(n_components, sum(eigvals > 0))
  if (d_use == 0) {
    stop("No positive eigenvalues found. Check your GRM or consider a different approach.")
  }
  
  #--------------------------
  # 3. Embed training set
  #    X_train_embedded = U_d * sqrt(Lambda_d)
  #--------------------------
  lambda_d <- diag(sqrt(eigvals[1:d_use]), d_use, d_use)
  U_d      <- eigvecs[, 1:d_use, drop = FALSE]
  
  X_train_embedded <- U_d %*% lambda_d   # (n x d_use)
  
  #--------------------------
  # 4. Embed test set
  #    X_test_embedded = G_train_test^T * U_d * Lambda_d^-1
  #
  #    -> dimension: (m x d_use)
  #--------------------------
  X_test_embedded <- t(G_train_test) %*% U_d  # (m x d_use)
  
  # invert sqrt(Lambda_d)
  lambda_d_inv <- diag(1 / diag(lambda_d), d_use, d_use)
  X_test_embedded <- X_test_embedded %*% lambda_d_inv
  
  #--------------------------
  # 5. Subset to correct training individuals
  #--------------------------
  idx_correct <- which(correctness)
  X_train_correct_embedded <- X_train_embedded[idx_correct, , drop = FALSE]
  
  if (nrow(X_train_correct_embedded) < k) {
    warning("Fewer than k correct training points. Adjusting k.")
    k <- nrow(X_train_correct_embedded)
  }
  
  if (k == 0) {
    # Edge case: if no correct individuals exist
    warning("No 'correct' training individuals found. Returning NA trust scores.")
    return(rep(NA, m))
  }
  
  #--------------------------
  # 6. Approximate k-NN search in the embedded space
  #--------------------------
  ann_results <- nn2(
    data  = X_train_correct_embedded,  # reference set
    query = X_test_embedded,           # test set
    k     = k
  )
  
  # ann_results$nn.dists is an (m x k) matrix of distances from each test sample
  # to its k approximate nearest neighbors among the correct training individuals.
  
  #--------------------------
  # 7. Compute Trust Scores
  #    trust_score_j = 1 / (1 + mean(distance_j))
  #--------------------------
  trust_scores <- apply(ann_results$nn.dists, 1, function(dists) {
    mean_dist <- mean(dists)
    1 / (1 + mean_dist)
  })
  
  return(trust_scores)
}
