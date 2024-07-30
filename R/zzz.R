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

# .onLoad <- function(libname, pkgname) {
#   required_packages <- c("reticulate", "tensorflow", "keras")
#
#   # Check for required packages and install if missing
#   for (pkg in required_packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       install.packages(pkg)
#     }
#   }
#
#   # Load the packages
#   lapply(required_packages, library, character.only = TRUE)
#
#   # Function to check if Python is installed and meets the minimum version requirement
#   check_python <- function(min_version = "3.0.0") {
#     python <- reticulate::py_discover_config()
#     if (is.null(python$python)) {
#       message("Python not found.")
#       return(FALSE)
#     } else {
#       current_version <- reticulate::py_version()
#       if (is.null(current_version)) {
#         message("Unable to determine Python version.")
#         return(FALSE)
#       } else {
#         message(paste("Python version found:", current_version))
#         return(compareVersion(current_version, min_version) >= 0)
#       }
#     }
#   }
#
#   # Ensure Python is installed and has the required version
#   if (!check_python()) {
#     message("Python 3 is not installed or the version is too low. Installing Miniconda with Python 3.")
#     tryCatch({
#       if (!file.exists(reticulate::miniconda_path())) {
#         reticulate::install_miniconda()
#       }
#       if (!any(reticulate::conda_list()$name == "r-reticulate")) {
#         reticulate::conda_create("r-reticulate", packages = c("python=3.10", "numpy"))
#       }
#       reticulate::use_condaenv("r-reticulate", required = TRUE)
#     }, error = function(e) {
#       message("Failed to install or configure Miniconda: ", e$message)
#       stop("Miniconda installation or configuration failed.")
#     })
#   } else {
#     reticulate::use_python(reticulate::py_discover_config()$python)
#   }
#
#   # Ensure numpy is installed in the Python environment
#   tryCatch({
#     reticulate::py_install("numpy", envname = "r-reticulate")
#   }, error = function(e) {
#     message("Failed to install numpy: ", e$message)
#     stop("Numpy installation failed.")
#   })
#
#   # Install TensorFlow and Keras if not installed
#   tryCatch({
#     tf_config <- tensorflow::tf_config()
#     if (is.null(tf_config) || is.null(tf_config$installed) || !tf_config$installed) {
#       message("TensorFlow not found. Installing TensorFlow and Keras.")
#       tensorflow::install_tensorflow(extra_packages = "tensorflow-probability")
#       keras::install_keras()
#     }
#   }, error = function(e) {
#     message("Failed to install TensorFlow or Keras: ", e$message)
#     stop("TensorFlow or Keras installation failed.")
#   })
#
#   message("All dependencies are loaded and configured.")
# }

