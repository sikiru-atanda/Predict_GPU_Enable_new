ensure_environment <- function() {
  script_path <- system.file("setup_environment.R", package = "PredictProR")
  message("Please run the following command to set up the Python environment:")
  message(paste0("source('", script_path, "')"))
}

### Setting Up the Python Environment

# After installing the `PredictProR` package,
# you need to set up the Python environment. Run the following R command
# to initialize the environment:
##PredictProR::ensure_environment()
