

# setup_environment.R

# Increase the timeout for network connections to avoid timeouts during large file downloads
options(timeout = 300) # Increase timeout to 300 seconds

# Function to install a package if it's missing
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Check and install necessary packages
install_if_missing("reticulate")
install_if_missing("semver")
library(semver)
library(reticulate)

# Function to detect the operating system
detect_os <- function() {
  os <- Sys.info()["sysname"]
  if (os == "Linux") {
    return("linux")
  } else if (os == "Darwin") { # macOS
    return("mac")
  } else if (os == "Windows") {
    return("windows")
  } else {
    stop("Unsupported operating system")
  }
}

# Function to install Miniconda
install_miniconda_dynamic <- function() {
  conda_installed <- tryCatch({
    reticulate::conda_binary()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!conda_installed) {
    message("Miniconda is not installed. Installing Miniconda...")
    reticulate::install_miniconda(force = TRUE)
  } else {
    message("Miniconda is already installed.")
  }
}

# Function to check if a specific Python version is installed within Conda
check_python_version <- function(env_name, required_version = "3.8") {
  py_path <- tryCatch({
    reticulate::conda_python(env_name)
  }, error = function(e) {
    message("Error retrieving Python path: ", e$message)
    return(NULL)
  })

  if (is.null(py_path) || py_path == "") {
    stop("Failed to retrieve Python path.")
  }

  py_version <- system2(py_path, "--version", stdout = TRUE, stderr = TRUE)
  if (grepl("Python", py_version)) {
    py_version <- sub("Python ", "", py_version)
  } else {
    stop("Failed to retrieve Python version.")
  }

  parsed_required_version <- semver::parse_version(required_version)
  parsed_py_version <- semver::parse_version(py_version)

  if (parsed_py_version$major == parsed_required_version$major &&
      parsed_py_version$minor == parsed_required_version$minor) {
    return(TRUE)
  }
  return(FALSE)
}

# Function to check if a Conda environment exists
conda_env_exists <- function(env_name) {
  envs <- reticulate::conda_list()
  return(env_name %in% envs$name)
}

# Function to ensure the Conda environment has the correct Python version
ensure_python_version <- function(env_name, required_version = "3.8") {
  if (conda_env_exists(env_name)) {
    if (!check_python_version(env_name, required_version)) {
      message("The existing Conda environment has a different Python version.")
      message("Removing the existing Conda environment: ", env_name)
      reticulate::conda_remove(env_name)
      message("Conda environment removed. Creating a new environment with Python ", required_version, "...")
      reticulate::conda_create(env_name, python_version = required_version)
    } else {
      message("The existing Conda environment matches the required Python version.")
    }
  } else {
    message("No existing Conda environment found. Creating a new environment with Python ", required_version, "...")
    reticulate::conda_create(env_name, python_version = required_version)
  }
}

# Function to install TensorFlow, Keras, and NumPy using pip
install_tensorflow_keras_numpy <- function(env_name) {
  tryCatch({
    message("Installing TensorFlow, Keras, and NumPy in the Python environment using pip...")
    reticulate::py_install(
      packages = c("tensorflow==2.13.0", "keras==2.13.1", "numpy==1.24.3"),
      envname = env_name,
      pip = TRUE
    )
    message("All packages installed successfully using pip.")
  }, error = function(e) {
    message("Failed to install packages: ", e$message)
    stop("One or more Python packages failed to install.")
  })
}

# Function to verify that TensorFlow, Keras, and NumPy are installed
verify_installation <- function() {
  tryCatch({
    reticulate::py_run_string("import tensorflow; import keras; import numpy")
    message("TensorFlow, Keras, and NumPy are available in the environment.")
  }, error = function(e) {
    message("Failed to verify installation: ", e$message)
    stop("TensorFlow, Keras, and/or NumPy are not correctly installed.")
  })
}

# Function to log actions and errors
log_message <- function(msg) {
  log_file <- file("setup_log.txt", open = "a")
  writeLines(paste(Sys.time(), "-", msg), con = log_file)
  close(log_file)
  message(msg)
}

# Main setup process

# Detect the OS
os <- detect_os()

# Define environment name and Python version
systime <- format(Sys.time(), "%m-%d-%Y_%I-%M%p")
env_name <- paste("myenv", systime, sep = "_")
python_version <- "3.8"

# Ensure the Conda environment has the correct Python version
ensure_python_version(env_name, python_version)

# Install Miniconda dynamically
install_miniconda_dynamic()

# Use the Conda environment
reticulate::use_condaenv(env_name, required = TRUE)

# Install required packages using pip
log_message("Starting Python environment setup...")
install_tensorflow_keras_numpy(env_name)
log_message("Package installation complete. Verifying packages...")

# Verify that the packages are installed
verify_installation()

log_message("Setup complete. Your Python environment is ready.")

# Additional improvements:

# Allow users to specify the desired Python version
# python_version <- readline("Enter the desired Python version (e.g., 3.8): ")
python_version <- "3.8"

# Use a requirements.txt file for dependency management
requirements_file <- "requirements.txt"
if (file.exists(requirements_file)) {
  message("Using requirements file:", requirements_file)
  reticulate::py_install(requirements = requirements_file, envname = env_name, pip = TRUE)
} else {
  message("No requirements file found.")
}

# Function to verify that TensorFlow, Keras, and NumPy are installed
verify_installation <- function() {
  tryCatch({
    reticulate::py_run_string("import tensorflow; import keras; import numpy")
    return(TRUE) # Return TRUE if successful
  }, error = function(e) {
    message("Failed to verify installation: ", e$message)
    return(FALSE) # Return FALSE if there's an error
  })
}

# ... (rest of the code)

# Added suggestion: If verification fails, retry with a different TensorFlow version
if (!verify_installation()) {
  message("Verification failed. Retrying with TensorFlow 2.12.0...")
  reticulate::py_install(packages = c("tensorflow==2.12.0", "keras==2.13.1", "numpy==1.24.3"), envname = env_name, pip = TRUE)
  verify_installation()
}


# Function to extract timestamp from environment name
extract_timestamp <- function(env_name) {
  # Assuming the env_name format is like "myenv_MM-DD-YYYY_HH-MMAM"
  timestamp_pattern <- "myenv_(\\d{2}-\\d{2}-\\d{4}_\\d{2}-\\d{2}[APM]{2})"
  matches <- regmatches(env_name, regexec(timestamp_pattern, env_name))
  if (length(matches[[1]]) > 1) {
    # Convert to POSIXct
    env_time <- as.POSIXct(matches[[1]][2], format = "%m-%d-%Y_%I-%M%p", tz = "UTC")
    return(env_time)
  } else {
    return(NA)
  }
}

# List conda environments
conda_envs <- reticulate::conda_list()

# Extract timestamps from environment names
conda_envs$created_at <- sapply(conda_envs$name, extract_timestamp)

# Get the current time
current_time <- Sys.time()

# Calculate the time difference in minutes
conda_envs$time_diff <- difftime(current_time, conda_envs$created_at, units = "mins")

# Filter environments older than 30 minutes
old_envs <- conda_envs[!is.na(conda_envs$created_at) & conda_envs$time_diff > 30, ]

# Remove old environments
if (nrow(old_envs) > 0) {
  for (env_name in old_envs$name) {
    message("Removing environment: ", env_name)
    reticulate::conda_remove(env_name)
  }
} else {
  message("No environments older than 30 minutes found.")
}





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
# options(timeout = 300) # Increase timeout to 300 seconds
#
# # Function to install a package if it's missing
# install_if_missing <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
# }
#
# # Check and install necessary packages
# install_if_missing("reticulate")
# install_if_missing("semver")
# library(semver)
# library(reticulate)
#
# # Function to detect the operating system
# detect_os <- function() {
#   os <- Sys.info()["sysname"]
#   if (os == "Linux") {
#     return("linux")
#   } else if (os == "Darwin") { # macOS
#     return("mac")
#   } else if (os == "Windows") {
#     return("windows")
#   } else {
#     stop("Unsupported operating system")
#   }
# }
#
# # Function to install Miniconda
# install_miniconda_dynamic <- function() {
#   conda_installed <- tryCatch({
#     reticulate::conda_binary()
#     TRUE
#   }, error = function(e) {
#     FALSE
#   })
#
#   if (!conda_installed) {
#     message("Miniconda is not installed. Installing Miniconda...")
#     reticulate::install_miniconda(force = TRUE)
#   } else {
#     message("Miniconda is already installed.")
#   }
# }
#
# # Function to check if a specific Python version is installed within Conda
# check_python_version <- function(env_name, required_version = "3.8") {
#   # Retrieve the Python path, with error handling
#   py_path <- tryCatch({
#     reticulate::conda_python(env_name)
#   }, error = function(e) {
#     message("Error retrieving Python path: ", e$message)
#     return(NULL)
#   })
#
#   # If Python path is NULL or empty, stop execution
#   if (is.null(py_path) || py_path == "") {
#     stop("Failed to retrieve Python path.")
#   }
#
#   # Retrieve the Python version and compare it with the required version
#   py_version <- system2(py_path, "--version", stdout = TRUE, stderr = TRUE)
#   if (grepl("Python", py_version)) {
#     py_version <- sub("Python ", "", py_version)
#   } else {
#     stop("Failed to retrieve Python version.")
#   }
#
#   parsed_required_version <- semver::parse_version(required_version)
#   parsed_py_version <- semver::parse_version(py_version)
#
#   if (parsed_py_version$major == parsed_required_version$major &&
#       parsed_py_version$minor == parsed_required_version$minor) {
#     return(TRUE)
#   }
#   return(FALSE)
# }
#
#
# # Function to check if a Conda environment exists
# conda_env_exists <- function(env_name) {
#   envs <- reticulate::conda_list()
#   return(env_name %in% envs$name)
# }
#
# # Function to ensure the Conda environment has the correct Python version
# ensure_python_version <- function(env_name, required_version = "3.8") {
#   if (conda_env_exists(env_name)) {
#     if (!check_python_version(env_name, required_version)) {
#       message("The existing Conda environment has a different Python version.")
#       message("Removing the existing Conda environment: ", env_name)
#       reticulate::conda_remove(env_name)
#       message("Conda environment removed. Creating a new environment with Python ", required_version, "...")
#       reticulate::conda_create(env_name, python_version = required_version)
#     } else {
#       message("The existing Conda environment matches the required Python version.")
#     }
#   } else {
#     message("No existing Conda environment found. Creating a new environment with Python ", required_version, "...")
#     reticulate::conda_create(env_name, python_version = required_version)
#   }
# }
#
# # Detect the OS
# os <- detect_os()
#
# systime <- format(Sys.time(), "%m-%d-%Y_%I-%M%p")
# # Define environment name and Python version
# env_name <- paste("myenv", systime, sep = "_")
# python_version <- "3.8"
#
# # Ensure the Conda environment has the correct Python version
# ensure_python_version(env_name, python_version)
#
# # Install Miniconda dynamically
# install_miniconda_dynamic()
#
# # Use the Conda environment
# reticulate::use_condaenv(env_name, required = TRUE)
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
#     if (semver::compare_versions(current_version, required_version) != 0) {
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
#     install_python_package("tensorflow", version = "2.10.0") # Specify compatible versions
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
# # Function to log actions and errors
# log_message <- function(msg) {
#   log_file <- file("setup_log.txt", open = "a")
#   writeLines(paste(Sys.time(), "-", msg), con = log_file)
#   close(log_file)
#   message(msg)
# }
#
# # Run the setup
# log_message("Starting Python environment setup...")
# update_pip()
# ensure_numpy()
# ensure_tf_keras()
# log_message("Setup complete. Your Python environment is ready.")
#
# # Additional improvements:
#
# # Allow users to specify the desired Python version
# # python_version <- readline("Enter the desired Python version (e.g., 3.8): ")
# python_version <- "3.8"
#
# # Use a requirements.txt file for dependency management
# requirements_file <- "requirements.txt"
# if (file.exists(requirements_file)) {
#   message("Using requirements file:", requirements_file)
#   reticulate::py_install(requirements = requirements_file, envname = env_name, pip = TRUE)
# } else {
#   message("No requirements file found.")
# }
