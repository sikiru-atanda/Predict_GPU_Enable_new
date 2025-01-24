compute_trust_score_grm_svd <- function(
    G_train_train,      # (n x n) GRM for training
    G_train_test,       # (n x m) cross-GRM for train vs. test
    y_train,            # training labels
    y_pred_train,       # predictions on training
    k = 5,
    correctness_tol = NULL,
    n_components = 5    # number of eigen components
) {
  library(FNN)  # for k-NN in Euclidean space
  
  n <- nrow(G_train_train)
  m <- ncol(G_train_test)
  
  # 1. Define correctness
  if (is.null(correctness_tol)) {
    correctness_tol <- 0.5 * sd(y_train, na.rm = TRUE)
  }
  errors <- abs(y_pred_train - y_train)
  correctness <- (errors <= correctness_tol)
  
  # 2. Eigen-decompose (or SVD) the training GRM
  #    - G_train_train should be symmetric (n x n).
  eig <- eigen(G_train_train, symmetric = TRUE)
  # eig$values, eig$vectors are the eigenvalues/vectors
  
  # 3. Keep top n_components
  #    Filter out non-positive eigenvalues if they exist
  eigvals  <- eig$values
  eigvecs  <- eig$vectors
  
  # quick check: if some eigenvalues are negative or zero (happens if not perfectly PSD),
  # we might set them to 0. This is simplistic:
  eigvals[eigvals < 1e-12] <- 0
  
  # Sort eigenvalues in descending order
  # (eigen() in R typically returns them in decreasing order, but let's ensure)
  # We'll assume they're already sorted. If not, you would reorder them.
  
  d_use <- min(n_components, sum(eigvals > 0))  # can't use more than # positive eigenvalues
  
  # 4. Build the embedding for training individuals
  #    Let’s define X_train = U_d * sqrt(Lambda_d)
  lambda_d <- diag(sqrt(eigvals[1:d_use]), d_use, d_use)
  U_d      <- eigvecs[, 1:d_use, drop = FALSE]
  
  X_train_embedded <- U_d %*% lambda_d   # (n x d_use)
  
  # 5. Embedding for test individuals
  #    Typically: X_test = G_train_test^T * U_d * Lambda_d^-1
  #    or a variant. But a common kernel PCA formula is:
  #      X_test_j = (G_train_test^T_j . U_d) * 1 / sqrt(Lambda_d)
  
  # Let's compute that:
  # G_train_test is (n x m), so G_train_test^T is (m x n).
  # We'll get an (m x d_use) matrix for X_test_embedded.
  
  X_test_embedded <- t(G_train_test) %*% U_d  # (m x d_use)
  # Now scale by 1 / sqrt(Lambda_d)
  # which is just solve(lambda_d)
  lambda_d_inv <- diag(1 / diag(lambda_d), d_use, d_use)
  X_test_embedded <- X_test_embedded %*% lambda_d_inv  # (m x d_use)
  
  # 6. Subset to correct training individuals for k-NN
  idx_correct <- which(correctness)
  X_train_correct_embedded <- X_train_embedded[idx_correct, , drop = FALSE]
  if (nrow(X_train_correct_embedded) < k) {
    warning("Fewer than k correct training points. Adjusting k.")
    k <- nrow(X_train_correct_embedded)
  }
  
  # 7. k-NN in the embedded space
  # Convert to data.frames or matrices if needed
  nn_results <- get.knnx(data = X_train_correct_embedded, query = X_test_embedded, k = k)
  
  # 8. Trust Score = 1 / (1 + average distance to k neighbors)
  trust_scores <- apply(nn_results$nn.dist, 1, function(dists) {
    mean_dist <- mean(dists)
    1 / (1 + mean_dist)
  })
  
  return(trust_scores)
}


trust_scores_grm <- compute_trust_score_grm_svd( 
    G_train_train = G_train_train,      # (n x n) GRM for training
    G_train_test = G_train_test,       # (n x m) cross-GRM for train vs. test
    y_train = y_train_geno,            # training labels
    y_pred_train = geno_preds_train,       # predictions on training
    k = 5
)

