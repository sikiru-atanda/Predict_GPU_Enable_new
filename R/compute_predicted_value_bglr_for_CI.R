# Function to compute predicted values for one omic source
compute_predicted_value <- function(geno_data, bMat_file = NULL,
                                    bMat,
                                    mu_values, chunk_size = 5000) {
  # Read marker effects for this omic source
  #bMat <- BGLR::readBinMat(bMat_file)  # M x p matrix

  # Get dimensions
  M <- nrow(bMat)  # Number of MCMC iterations
  p <- ncol(bMat)  # Number of markers
  N <- nrow(geno_data)  # Number of individuals

  # Initialize matrix for predictions
  Predicted_value <- matrix(0, nrow = N, ncol = M)

  if(ncol(geno_data)>=chunk_size & ncol(geno_data)!= nrow(geno_data)){
  # Process in chunks to avoid memory overload
  for (i in seq(1, p, by = chunk_size)) {
    #cat(sprintf("Processing markers %d to %d...\n", i, min(i + chunk_size - 1, p)))

    # Define range for slicing
    end_idx <- min(i + chunk_size - 1, p)

    # Extract a chunk of bMat (M x chunk_size)
    bMat_chunk <- bMat[, i:end_idx, drop = FALSE]  # (M x chunk_size)

    # Extract corresponding geno_data slice (N x chunk_size)
    geno_chunk <- geno_data[, i:end_idx, drop = FALSE]  # (N x chunk_size)

    # Compute partial predictions for this chunk
    Predicted_value <- Predicted_value + (geno_chunk %*% t(bMat_chunk))  # (N x M)
  }
  } else {

    Predicted_value <- geno_data %*% t(bMat)
}
  # Add intercept (mu) for this omic source
  #Predicted_value <- Predicted_value + matrix(mu_values, nrow = N, ncol = M, byrow = TRUE)

  return(Predicted_value)
}
