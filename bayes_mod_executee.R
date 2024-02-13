# No, the "multisession" mode in the future package is designed to work on
# Windows systems and is not supported on Linux. On Linux,
# you typically use the "multicore" mode or "multiprocess" mode
# for parallel processing with the future package.
#
# If you want to use multiple sessions on Linux,
# you can use the "multicore" mode.

# Install and load required packages
# if (!requireNamespace("BGLR", quietly = TRUE)) {
#   install.packages("BGLR")
# }
# if (!requireNamespace("doParallel", quietly = TRUE)) {
#   install.packages("doParallel")
# }
# if (!requireNamespace("future", quietly = TRUE)) {
#   install.packages("future")
# }
#
# # Load packages
# library(BGLR)
# library(doParallel)
# library(future)

# Function for parameter check
parameter_check <- function(...) {
  # Check if required parameters are present
  stopifnot(!is.null(bayes_para), !is.null(bayes_para$nIter), !is.null(bayes_para$burnIn), !is.null(bayes_para$thin))
}

# Main function
bayes_mod_executee <- function(pheno_data = NULL,
                              response = NULL,
                              weights = NULL,
                              ETA = NULL,
                              bayes_para = NULL,
                              verbose = FALSE,
                              num_cores = NULL,
                              ...) {

  # Parameter check
  parameter_check(bayes_para = bayes_para)

  # Detect number of cores based on the operating system
  if(length(response) > 1){
  if(is.null(num_cores)){
  if (Sys.info()["sysname"] == "Windows") {
    num_cores <- parallel::detectCores()
    num_cores <- round(num_cores*0.5)
    # Use future for parallel processing
    future::plan("multisession", workers = num_cores)
  } else {
    num_cores <- parallel::detectCores(logical = FALSE)
    num_cores <- round(num_cores*0.5)
    future::plan("multicore", workers = num_cores)
  }

  } else {
    if (Sys.info()["sysname"] == "Windows") {
      # Use future for parallel processing
      future::plan("multisession", workers = num_cores)
    } else {
      future::plan("multicore", workers = num_cores)
    }

  }

  }


  # Common parameters
  common_params <- list(
    ETA = ETA,
    nIter = bayes_para$nIter,
    burnIn = bayes_para$burnIn,
    thin = bayes_para$thin,
    verbose = FALSE
  )

  # Function for individual trait analysis


  analyze_trait <- function(trait) {
    current_date_time <- as.character(Sys.time())
    files_key <- gsub(" ", "", current_date_time)
    files_key <- gsub("[-:]", "_", files_key)
    files_key <- paste0(files_key, trait)

    tryCatch(
      {
    # Run BGLR
    if (is.null(weights)) {
      fm <- BGLR::BGLR(y = pheno_data[, trait],
                       ETA=ETA$ETA,
                       #weights = weights,
                       nIter= bayes_para$nIter,
                       burnIn= bayes_para$burnIn,
                       thin = bayes_para$thin,
                       verbose = FALSE,
                       saveAt = files_key)
    } else {
      fm <- BGLR::BGLR(y = pheno_data[, trait], ETA=ETA,
                       weights = weights,
                       nIter= bayes_para$nIter,
                       burnIn= bayes_para$burnIn,
                       thin = bayes_para$thin,
                       verbose = FALSE,
                       saveAt = files_key)
    }

    # Get output files names
    output_files_names <- list.files(pattern = files_key)

    # Return output
    return(list(model = fm, output_files_names = output_files_names))

    #Res = list(model = fm, output_files_names = output_files_names)

    #names(Res) <- trait
    #return(Res)
    },
    error = function(e) {
      # Handle the error, you can print a message or take other actions
      cat("Error in analyze_trait for trait", trait, ":", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or an appropriate value to indicate the failure
    }
    )
  }

  if (length(response) > 1){
  # Apply analyze_trait function in parallel or sequentially
    options(future::future.rng.onMisuse = "ignore")
  Univariate <- future.apply::future_lapply(response, analyze_trait)

  # Close parallel cluster
  if (length(response) > 1) {
    future::plan("sequential")
  }

  # Combine results
  #output <- do.call(c, Univariate)

  output <- Univariate

  names(output) <- response


  } else {

    output <- analyze_trait(response)

  }

  return(output)
  # Clean up
  # if (length(response) > 1) {
  #   return(output)
  # } else {
  #   return(output$model)
  # }
}

# Example usage
# Set up your data, response, weights, ETA, and bayes_para
result <- bayes_mod_executee(pheno_data = pheno,
                            response = c("Yield", "deBLUP", "BLUE" ),
                            weights = NULL,
                            ETA = ETA,
                            bayes_para = bayes_para)


result$BLUE
pheno[, "BLUE"] = "NA"

str(pheno)
