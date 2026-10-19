ensure_environment <- function() {
  message("Set PREDICTPRO_DL_PYTHON to a server/runtime Python (>= 3.10); PredictProR installs the DL packages into it on first use.")
  message("Alternatively, run PredictProR::setup_predictdl_env(prefer_gpu = TRUE, cuda = 'auto') once to provision it.")
}

### Setting Up the Python Environment

# After installing the `PredictProR` package,
# you need to set up the Python environment. Run the following R command
# to initialize the environment:
##PredictProR::ensure_environment()
