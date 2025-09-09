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


.onLoad <- function(libname, pkgname) {
  # Harmless default; don't perform installs here
  options(timeout = max(getOption("timeout", 60), 300))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "To set up a Python env with PyTorch, run:\n",
    "  setup_predictdl_env(prefer_gpu = TRUE, cuda = 'auto')"
  )
}


# .onLoad <- function(libname, pkgname) {
#   if (!reticulate::py_module_available("tensorflow") || !reticulate::py_module_available("keras")) {
#     packageStartupMessage("TensorFlow and/or Keras not found. Please run 'setup_environment.R' to install the necessary dependencies.")
#   }
# }


# .onLoad <- function(libname, pkgname) {
#   # install_if_missing <- function(pkg) {
#   #   if (!requireNamespace(pkg, quietly = TRUE)) {
#   #     install.packages(pkg)
#   #   }
#   # }
#   #
#   # required_packages <- c("reticulate", "tensorflow", "keras")
#   # lapply(required_packages, install_if_missing)
#   #
#   # library(reticulate)
#
#   # Check if Miniconda is already installed
#   conda_installed <- tryCatch({
#     reticulate::conda_binary()
#     TRUE
#   }, error = function(e) {
#     FALSE
#   })
#
#   if (!conda_installed) {
#     message("Python is not installed. Installing Miniconda Python...")
#     reticulate::install_miniconda()
#   } else {
#     message("Miniconda is already installed.")
#   }
#
#   # Create or use an existing virtual environment
#   env_name <- "myenv"
#   if (!reticulate::virtualenv_exists(env_name)) {
#     message("Creating a virtual environment...")
#     reticulate::virtualenv_create(env_name, python_version = "3.8")  # Specify a compatible Python version
#   } else {
#     message("Using existing virtual environment: ", env_name)
#   }
#
#   # Use the virtual environment
#   reticulate::use_virtualenv(env_name, required = TRUE)
#
#   install_python_package <- function(package_name, version = NULL) {
#     tryCatch({
#       message(paste("Installing", package_name, "in the Python environment..."))
#       reticulate::py_install(package_name, envname = env_name, pip = TRUE, version = version)
#       message(paste(package_name, "installed successfully."))
#     }, error = function(e) {
#       message(paste("Failed to install", package_name, ":", e$message))
#     })
#   }
#
#   ensure_numpy <- function(required_version = "1.26.4") {
#     tryCatch({
#       result <- reticulate::py_run_string("import numpy; version = numpy.__version__")
#       current_version <- as.character(result$version)
#       message(paste("Current NumPy version:", current_version))
#       if (utils::compareVersion(current_version, required_version) != 0) {
#         message(paste("Installing/upgrading NumPy to version:", required_version))
#         install_python_package("numpy", version = required_version)
#       } else {
#         message("Required NumPy version is already installed.")
#       }
#     }, error = function(e) {
#       message("Failed to ensure NumPy version:", e$message)
#       install_python_package("numpy", version = required_version)
#     })
#   }
#
#   ensure_tf_keras <- function() {
#     tryCatch({
#       reticulate::py_run_string("import tensorflow")
#       reticulate::py_run_string("import keras")
#     }, error = function(e) {
#       message("TensorFlow and/or Keras not found in the Python environment. Installing...")
#       install_python_package("tensorflow", version = "2.10.0")  # Specify compatible versions
#       install_python_package("keras", version = "2.10.0")
#     })
#   }
#
#   # Run the pip update and ensure necessary packages
#   update_pip <- function() {
#     tryCatch({
#       reticulate::py_run_string("import pip")
#       reticulate::py_install("pip", envname = env_name, pip = TRUE)
#     }, error = function(e) {
#       message("pip is not available. Installing pip...")
#       install_python_package("pip")
#     })
#   }
#
#   update_pip()
#   ensure_numpy()
#   ensure_tf_keras()
# }




# .onLoad <- function(libname, pkgname) {
#   # Helper function to install R packages if they are not already installed
#   install_if_missing <- function(pkg) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       install.packages(pkg)
#     }
#   }
#
#   # List of required R packages
#   required_packages <- c("reticulate", "tensorflow", "keras")
#
#   # Install required R packages
#   lapply(required_packages, install_if_missing)
#
#   # Load reticulate
#   library(reticulate)
#
#   # Function to run a shell command and capture output
#   run_shell_command <- function(command) {
#     tryCatch({
#       result <- system(command, intern = TRUE)
#       message(paste("Command output:", paste(result, collapse = "\n")))
#       return(result)
#     }, error = function(e) {
#       message(paste("Failed to run command:", command, "Error:", e$message))
#     })
#   }
#
#   # Function to install a package using pip in the Python environment
#   install_python_package <- function(package_name) {
#     tryCatch({
#       message(paste("Installing", package_name, "in the Python environment..."))
#       run_shell_command(paste("python -m pip install --upgrade", package_name))
#       message(paste(package_name, "installed successfully."))
#     }, error = function(e) {
#       message(paste("Failed to install", package_name, ":", e$message))
#     })
#   }
#
#   # Ensure pip is installed and up-to-date
#   ensure_pip <- function() {
#     tryCatch({
#       reticulate::py_run_string("import pip")
#     }, error = function(e) {
#       message("pip not found. Installing pip...")
#       install_python_package("pip")
#     })
#   }
#
#   update_pip <- function() {
#     tryCatch({
#       run_shell_command("python -m pip install --upgrade pip")
#       message("pip has been updated.")
#     }, error = function(e) {
#       message("Failed to update pip: ", e$message)
#     })
#   }
#
#   # Ensure NumPy is installed with the correct version
#   ensure_numpy <- function(required_version = "1.26.4") {
#     tryCatch({
#       # Check current numpy version
#       current_version <- reticulate::py_run_string("import numpy; print(numpy.__version__)")
#       current_version <- as.character(current_version)
#       message(paste("Current NumPy version:", current_version))
#       if (current_version != required_version) {
#         message(paste("Installing/upgrading NumPy to version:", required_version))
#         install_python_package(paste("numpy==", required_version, sep = ""))
#       } else {
#         message("Required NumPy version is already installed.")
#       }
#     }, error = function(e) {
#       message("Failed to ensure NumPy version:", e$message)
#     })
#   }
#
#   # Ensure TensorFlow and Keras are installed
#   ensure_tf_keras <- function() {
#     tryCatch({
#       reticulate::py_run_string("import tensorflow")
#       reticulate::py_run_string("import keras")
#     }, error = function(e) {
#       message("TensorFlow and/or Keras not found in the Python environment.")
#       install_python_package("tensorflow")
#       install_python_package("keras")
#     })
#   }
#
#   # Run the pip update and ensure necessary packages
#   ensure_pip()
#   update_pip()
#   ensure_numpy()
#   ensure_tf_keras()
# }
#
#
#
