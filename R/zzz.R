# .onAttach <- function(libname, pkgname) {
#   # Path to the logo file
#   logo_path <- system.file("NDSUlogo.png", package = pkgname)
#
#   # Display the logo and custom message
#   packageStartupMessage("This is a product of sik
# ----------------------------------------
# ")
#   if (file.exists(logo_path)) {
#     logo <- png::readPNG(logo_path)
#     grid::grid.raster(logo)
#   }
#
#   packageStartupMessage("
# Welcome to North Dakota State University R Package!
# ----------------------------------------
# ")
# }


.onAttach <- function(libname, pkgname) {
  # ANSI escape codes for colors
  green <- "\033[32m"
  yellow <- "\033[33m"
  reset <- "\033[0m"

  # ASCII art of the NDSU logo with colors
  logo_ascii <- paste0(
    green, " _   _  ____  ____  _    _  ", reset, "\n",
    green, "| \\ | ||  _ \\|  _ \\| |  | | ", reset, "\n",
    green, "|  \\| || | | | | | | |  | | ", reset, "\n",
    green, "| . ` || |_| | |_| | |__| | ", reset, "\n",
    green, "|_|\\_||____/|____/ \\____/  ", reset, "\n",
    yellow, "  North Dakota State University", reset, "\n"
  )

  # Display the ASCII logo and custom message
  packageStartupMessage("
PredictProR Product of NDSU!
----------------------------------------
")
  packageStartupMessage(logo_ascii)

}

ensure_keras_and_dependencies <- function() {
  #library(reticulate)

  # Set up virtual environment
  venv_name <- "r-reticulate"

  tryCatch({
    # if (!reticulate::virtualenv_exists(venv_name)) {
    #   reticulate::use_virtualenv(venv_name, required = TRUE)
    # }



    # Install numpy
    if (!reticulate::py_module_available("numpy")) {
      reticulate::py_install("numpy==1.24.2", envname = venv_name, pip = TRUE)
    }

    # Install TensorFlow and Keras using keras::install_keras
    if (!reticulate::py_module_available("tensorflow")) {
      keras::install_keras(method = "virtualenv", envname = venv_name)
    }

    keras::install_keras()
    # Install keras-tuner
    if (!reticulate::py_module_available("keras_tuner")) {
      reticulate::py_install("keras-tuner", envname = venv_name, pip = TRUE)
    }
  }, error = function(e) {
    message("Error installing Python packages: ", e$message)
    message("Please ensure you have Python installed and accessible from R.")
  })
}

# Run the function to ensure the environment and dependencies are set up
ensure_keras_and_dependencies()

# Now you can proceed with the deep learning model script

### using conda as alternative

# ensure_keras_and_dependencies <- function() {
#   library(reticulate)
#   library(keras)
#
#   # Set up conda environment
#   conda_env_name <- "r-reticulate"
#
#   tryCatch({
#     # Check if conda is installed
#     if (is.null(reticulate::conda_binary())) {
#       stop("Conda is not installed. Please install Conda from https://docs.conda.io/en/latest/miniconda.html.")
#     }
#
#     # Create the conda environment if it doesn't exist
#     if (!conda_env_name %in% reticulate::conda_list()$name) {
#       reticulate::conda_create(envname = conda_env_name, packages = "python=3.8")
#     }
#
#     reticulate::use_condaenv(conda_env_name, required = TRUE)
#
#     # Install numpy
#     if (!reticulate::py_module_available("numpy")) {
#       reticulate::conda_install(envname = conda_env_name, packages = "numpy==1.24.2")
#     }
#
#     # Install TensorFlow and Keras using keras::install_keras
#     if (!reticulate::py_module_available("tensorflow")) {
#       keras::install_keras(method = "conda", conda = "auto", envname = conda_env_name)
#     }
#
#     # Install keras-tuner
#     if (!reticulate::py_module_available("keras_tuner")) {
#       reticulate::conda_install(envname = conda_env_name, packages = "keras-tuner")
#     }
#   }, error = function(e) {
#     message("Error installing Python packages: ", e$message)
#     message("Please ensure you have Conda installed and accessible from R.")
#   })
# }
#
# # Run the function to ensure the environment and dependencies are set up
# ensure_keras_and_dependencies()
