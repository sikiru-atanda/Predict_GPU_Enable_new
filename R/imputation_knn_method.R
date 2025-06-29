set_parallel_plan_impute <- function(num_chunks = 1,
                                     num_cores = NULL,
                                     sys_name) {
  # Define the plan based on the system
  plan_type <- ifelse(sys_name == "Windows", "multisession", "multicore")

  # Check if parallel execution is beneficial
  if (num_chunks > 1) {
    if(is.null(num_cores)){

      num_cores <-  parallel::detectCores()
      num_cores <- num_cores*0.5
    }
    future::plan(plan_type, workers = num_cores)
  } else {
    future::plan("sequential")
  }
}

# Function to perform KNN imputation on a chunk
# impute_chunk <- function(chunk_data, k) {
#   # Apply KNN imputation
#   imputed_chunk <- VIM::kNN(chunk_data, k = k, imp_var = FALSE)
#   imputed_chunk <- as.data.frame(imputed_chunk)
#   return(imputed_chunk)
# }
impute_chunk <- function(chunk_data, k = 5) {
  tryCatch({
    # Initial KNN imputation with default k
    imputed_chunk <- VIM::kNN(chunk_data, k = k, imp_var = FALSE)
    imputed_chunk <- as.data.frame(imputed_chunk)
    return(imputed_chunk)
  }, error = function(e) {
    # If an error occurs, adjust k dynamically based on non-NA values
    non_na_count <- rowSums(!is.na(chunk_data))
    adjusted_k <- pmin(k, max(non_na_count))
    adjusted_k <- pmax(adjusted_k, 2)  # Ensure adjusted_k is at least 2

    # Try KNN imputation again with the adjusted k
    tryCatch({
      imputed_chunk <- VIM::kNN(chunk_data, k = adjusted_k, imp_var = FALSE)
      imputed_chunk <- as.data.frame(imputed_chunk)
      return(imputed_chunk)
    }, error = function(e) {
      # If an error occurs again, force k = 2
      imputed_chunk <- VIM::kNN(chunk_data, k = 2, imp_var = FALSE)
      imputed_chunk <- as.data.frame(imputed_chunk)
      return(imputed_chunk)
    })
  })
}

# Function to handle large-scale KNN imputation with parallel processing
handle_large_scale_knn <- function(data, k = 5, chunk_size = 50, num_cores = NULL) {

  msg <- "\n==================================================\n"

  num_cols <- ncol(data)
  row_names <- rownames(data)
  num_chunks <- ceiling(num_cols / chunk_size)

  # Create list of chunk indices
  chunk_indices <- lapply(1:num_chunks, function(i) {
    start_col <- (i - 1) * chunk_size + 1
    end_col <- min(i * chunk_size, num_cols)
    list(start_col = start_col, end_col = end_col)
  })

  sys_name <- Sys.info()["sysname"]
  if (!is.null(num_cores) && num_cores > 1) {
    #sys_name <- Sys.info()["sysname"]
    set_parallel_plan_impute(num_chunks = num_chunks,num_cores = num_cores,
                             sys_name = sys_name)
  } else {
    # Automatically determine the number of cores and use half of them
    detected_cores <- parallel::detectCores(logical = TRUE)
    # For non-Windows systems, consider physical cores only
    num_cores <- round(detected_cores * 0.5)

    set_parallel_plan_impute(num_chunks = num_chunks,
                             num_cores = num_cores,
                             sys_name = sys_name)
  }

  # Perform imputation on each chunk in parallel
  imputed_chunks <- future.apply::future_lapply(chunk_indices, function(indices) {
    chunk_data <- data[, indices$start_col:indices$end_col, drop = FALSE]

    # Impute the chunk
    imputed_chunk <- impute_chunk(chunk_data, k)

    # Ensure imputed_chunk is a matrix to avoid subscript issues
    imputed_chunk <- as.matrix(imputed_chunk)

    # Check dimensions
    if (nrow(imputed_chunk) != nrow(data) || ncol(imputed_chunk) != (indices$end_col - indices$start_col + 1)) {
      stop(paste(msg,"Dimensions of imputed_chunk do not match the expected dimensions."), call. = FALSE)
    }

    return(imputed_chunk)
  }, future.seed = TRUE)

  future::plan("sequential")
  # Combine all imputed chunks into the final data frame
  imputed_data <- do.call(cbind, imputed_chunks)

  rownames(imputed_data) <- row_names
  return(as.data.frame(imputed_data))
}

# # Example usage
# set.seed(123)
# n_rows <- 100  # Number of samples
# n_cols <- 1000  # Number of SNP markers
# snp_data <- matrix(sample(c(0, 1, NA), n_rows * n_cols, replace = TRUE, prob = c(0.4, 0.4, 0.2)),
#                    nrow = n_rows, ncol = n_cols)
# snp_data <- as.data.frame(snp_data)  # Convert to data frame
# colnames(snp_data) <- paste("snp", 1:ncol(snp_data), sep = "_")
#
# # Perform KNN imputation on the large dataset in parallel
# imputed_snp_data <- handle_large_scale_knn(snp_data, k = 5, chunk_size = 50, num_cores=NULL)
# #
# # Display a summary of the imputed data
# print("Imputed SNP Data (Chunk-Based KNN Imputation):")
# print(summary(imputed_snp_data))

