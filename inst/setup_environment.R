# # setup_environment.R
#
# # Function to install a package if it's missing
# install_if_missing <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
# }
#
# # Check and install reticulate
# install_if_missing("reticulate")
#
# # Load the reticulate library
# library(reticulate)
#
# # Check if Miniconda is already installed
# conda_installed <- tryCatch({
#   reticulate::conda_binary()
#   TRUE
# }, error = function(e) {
#   FALSE
# })
#
# if (!conda_installed) {
#   message("Python is not installed. Installing Miniconda Python...")
#   reticulate::install_miniconda()
# } else {
#   message("Miniconda is already installed.")
# }
#
# # Create or use an existing virtual environment
# env_name <- "myenv"
# if (!reticulate::virtualenv_exists(env_name)) {
#   message("Creating a virtual environment...")
#   reticulate::virtualenv_create(env_name, python_version = "3.8")  # Specify a compatible Python version
# } else {
#   message("Using existing virtual environment: ", env_name)
# }
#
# # Use the virtual environment
# reticulate::use_virtualenv(env_name, required = TRUE)
#
# # Function to install Python packages
# install_python_package <- function(package_name, version = NULL) {
#   tryCatch({
#     message(paste("Installing", package_name, "in the Python environment..."))
#     reticulate::py_install(package_name, envname = env_name, pip = TRUE, version = version)
#     message(paste(package_name, "installed successfully."))
#   }, error = function(e) {
#     message(paste("Failed to install", package_name, ":", e$message))
#   })
# }
#
# # Ensure NumPy is installed
# ensure_numpy <- function(required_version = "1.26.4") {
#   tryCatch({
#     result <- reticulate::py_run_string("import numpy; version = numpy.__version__")
#     current_version <- as.character(result$version)
#     message(paste("Current NumPy version:", current_version))
#     if (utils::compareVersion(current_version, required_version) != 0) {
#       message(paste("Installing/upgrading NumPy to version:", required_version))
#       install_python_package("numpy", version = required_version)
#     } else {
#       message("Required NumPy version is already installed.")
#     }
#   }, error = function(e) {
#     message("Failed to ensure NumPy version:", e$message)
#     install_python_package("numpy", version = required_version)
#   })
# }
#
# # Ensure TensorFlow and Keras are installed
# ensure_tf_keras <- function() {
#   tryCatch({
#     reticulate::py_run_string("import tensorflow")
#     reticulate::py_run_string("import keras")
#   }, error = function(e) {
#     message("TensorFlow and/or Keras not found in the Python environment. Installing...")
#     install_python_package("tensorflow", version = "2.10.0")  # Specify compatible versions
#     install_python_package("keras", version = "2.10.0")
#   })
# }
#
# # Run the pip update and ensure necessary packages
# update_pip <- function() {
#   tryCatch({
#     reticulate::py_run_string("import pip")
#     reticulate::py_install("pip", envname = env_name, pip = TRUE)
#   }, error = function(e) {
#     message("pip is not available. Installing pip...")
#     install_python_package("pip")
#   })
# }
#
# # Run the setup
# update_pip()
# ensure_numpy()
# ensure_tf_keras()
#
# message("Setup complete. Your Python environment is ready.")
#
#
# # Run this command after installing the package
# # source(system.file("setup_environment.R", package = "PredictProR"))


# setup_environment.R

# # Increase the timeout for network connections to avoid timeouts during large file downloads
# options(timeout = 300)  # Increase timeout to 300 seconds
#
# # Function to install a package if it's missing
# install_if_missing <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
# }
#
# # Check and install reticulate
# install_if_missing("reticulate")
#
# # Load the reticulate library
# library(reticulate)
#
# # Function to detect the operating system
# detect_os <- function() {
#   os <- Sys.info()["sysname"]
#   if (os == "Linux") {
#     return("linux")
#   } else if (os == "Darwin") {  # macOS
#     return("mac")
#   } else if (os == "Windows") {
#     return("windows")
#   } else {
#     stop("Unsupported operating system")
#   }
# }
#
# # Function to install Miniconda based on OS
# install_miniconda_dynamic <- function() {
#   os <- detect_os()
#   conda_installed <- tryCatch({
#     reticulate::conda_binary()
#     TRUE
#   }, error = function(e) {
#     FALSE
#   })
#
#   if (!conda_installed) {
#     message("Miniconda is not installed. Installing Miniconda...")
#
#     if (os == "linux") {
#       reticulate::install_miniconda("https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh")
#     } else if (os == "mac") {
#       reticulate::install_miniconda("https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh")
#     } else if (os == "windows") {
#       reticulate::install_miniconda()  # Default URL works for Windows
#     }
#
#   } else {
#     message("Miniconda is already installed.")
#   }
# }
#
# # Detect the OS
# os <- detect_os()
#
# # Install Miniconda dynamically based on OS
# install_miniconda_dynamic()
#
# # Create or use an existing virtual environment
# env_name <- "myenv"
# if (!reticulate::virtualenv_exists(env_name)) {
#   message("Creating a virtual environment...")
#   reticulate::virtualenv_create(env_name, python_version = "3.8")  # Specify a compatible Python version
# } else {
#   message("Using existing virtual environment: ", env_name)
# }
#
# # Use the virtual environment
# reticulate::use_virtualenv(env_name, required = TRUE)
#
# # Function to install Python packages
# install_python_package <- function(package_name, version = NULL) {
#   tryCatch({
#     message(paste("Installing", package_name, "in the Python environment..."))
#     reticulate::py_install(package_name, envname = env_name, pip = TRUE, version = version)
#     message(paste(package_name, "installed successfully."))
#   }, error = function(e) {
#     message(paste("Failed to install", package_name, ":", e$message))
#   })
# }
#
# # Ensure NumPy is installed
# ensure_numpy <- function(required_version = "1.26.4") {
#   tryCatch({
#     result <- reticulate::py_run_string("import numpy; version = numpy.__version__")
#     current_version <- as.character(result$version)
#     message(paste("Current NumPy version:", current_version))
#     if (utils::compareVersion(current_version, required_version) != 0) {
#       message(paste("Installing/upgrading NumPy to version:", required_version))
#       install_python_package("numpy", version = required_version)
#     } else {
#       message("Required NumPy version is already installed.")
#     }
#   }, error = function(e) {
#     message("Failed to ensure NumPy version:", e$message)
#     install_python_package("numpy", version = required_version)
#   })
# }
#
# # Ensure TensorFlow and Keras are installed
# ensure_tf_keras <- function() {
#   tryCatch({
#     reticulate::py_run_string("import tensorflow")
#     reticulate::py_run_string("import keras")
#   }, error = function(e) {
#     message("TensorFlow and/or Keras not found in the Python environment. Installing...")
#     install_python_package("tensorflow", version = "2.10.0")  # Specify compatible versions
#     install_python_package("keras", version = "2.10.0")
#   })
# }
#
# # Run the pip update and ensure necessary packages
# update_pip <- function() {
#   tryCatch({
#     reticulate::py_run_string("import pip")
#     reticulate::py_install("pip", envname = env_name, pip = TRUE)
#   }, error = function(e) {
#     message("pip is not available. Installing pip...")
#     install_python_package("pip")
#   })
# }
#
# # Run the setup
# update_pip()
# ensure_numpy()
# ensure_tf_keras()
#
# message("Setup complete. Your Python environment is ready.")


# Increase the timeout for network connections to avoid timeouts during large file downloads
options(timeout = 300)  # Increase timeout to 300 seconds

# Function to install a package if it's missing
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Check and install reticulate
install_if_missing("reticulate")

# Load the reticulate library
library(reticulate)

# Function to detect the operating system
detect_os <- function() {
  os <- Sys.info()["sysname"]
  if (os == "Linux") {
    return("linux")
  } else if (os == "Darwin") {  # macOS
    return("mac")
  } else if (os == "Windows") {
    return("windows")
  } else {
    stop("Unsupported operating system")
  }
}

# Function to install Miniconda based on OS
install_miniconda_dynamic <- function() {
  os <- detect_os()
  conda_installed <- tryCatch({
    reticulate::conda_binary()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!conda_installed) {
    message("Miniconda is not installed. Installing Miniconda...")

    if (os == "linux") {
      reticulate::install_miniconda("https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh")
    } else if (os == "mac") {
      reticulate::install_miniconda("https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh")
    } else if (os == "windows") {
      reticulate::install_miniconda()  # Default URL works for Windows
    }

  } else {
    message("Miniconda is already installed.")
  }
}

# Function to check if a specific Python version is installed within Miniconda
check_python_version <- function(env_name, required_version = "3.8") {
  py_path <- tryCatch({
    reticulate::virtualenv_python(env_name)
  }, error = function(e) {
    return(NULL)
  })

  if (!is.null(py_path)) {
    py_version <- system2(py_path, "--version", stdout = TRUE, stderr = TRUE)
    py_version <- sub("Python ", "", py_version)
    if (startsWith(py_version, required_version)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Detect the OS
os <- detect_os()

# Install Miniconda dynamically based on OS
install_miniconda_dynamic()

# Create or use an existing virtual environment
env_name <- "myenv_specific_to_project"
python_version <- "3.8"

if (!reticulate::virtualenv_exists(env_name) || !check_python_version(env_name, python_version)) {
  message("Creating a virtual environment...")
  reticulate::virtualenv_create(env_name, python_version = python_version)  # Specify a compatible Python version
} else {
  message("Using existing virtual environment: ", env_name)
}

# Use the virtual environment
reticulate::use_virtualenv(env_name, required = TRUE)

# Function to install Python packages
install_python_package <- function(package_name, version = NULL) {
  tryCatch({
    message(paste("Installing", package_name, "in the Python environment..."))
    reticulate::py_install(package_name, envname = env_name, pip = TRUE, version = version)
    message(paste(package_name, "installed successfully."))
  }, error = function(e) {
    message(paste("Failed to install", package_name, ":", e$message))
  })
}

# Ensure NumPy is installed
ensure_numpy <- function(required_version = "1.26.4") {
  tryCatch({
    result <- reticulate::py_run_string("import numpy; version = numpy.__version__")
    current_version <- as.character(result$version)
    message(paste("Current NumPy version:", current_version))
    if (utils::compareVersion(current_version, required_version) != 0) {
      message(paste("Installing/upgrading NumPy to version:", required_version))
      install_python_package("numpy", version = required_version)
    } else {
      message("Required NumPy version is already installed.")
    }
  }, error = function(e) {
    message("Failed to ensure NumPy version:", e$message)
    install_python_package("numpy", version = required_version)
  })
}

# Ensure TensorFlow and Keras are installed
ensure_tf_keras <- function() {
  tryCatch({
    reticulate::py_run_string("import tensorflow")
    reticulate::py_run_string("import keras")
  }, error = function(e) {
    message("TensorFlow and/or Keras not found in the Python environment. Installing...")
    install_python_package("tensorflow", version = "2.10.0")  # Specify compatible versions
    install_python_package("keras", version = "2.10.0")
  })
}

# Run the pip update and ensure necessary packages
update_pip <- function() {
  tryCatch({
    reticulate::py_run_string("import pip")
    reticulate::py_install("pip", envname = env_name, pip = TRUE)
  }, error = function(e) {
    message("pip is not available. Installing pip...")
    install_python_package("pip")
  })
}

# Function to log actions and errors
log_message <- function(msg) {
  log_file <- file("setup_log.txt", open = "a")
  writeLines(paste(Sys.time(), "-", msg), con = log_file)
  close(log_file)
  message(msg)
}

# Run the setup
log_message("Starting Python environment setup...")
update_pip()
ensure_numpy()
ensure_tf_keras()
log_message("Setup complete. Your Python environment is ready.")
