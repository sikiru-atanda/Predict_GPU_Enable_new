# Function to check if PLINK is installed and meets the required version
check_and_install_plink_gud_but_nt_use <- function(min_version = "1.90",
                                                   plink_path = NULL) {

  # Automatically correct the path
  plink_path <- "D:/PredictProR/plink"
  # Step 1: Find all instances of plink
  # Determine the operating system
  os_type <- .Platform$OS.type
  system_command <- ifelse(.Platform$OS.type == "windows", "where", "which -a")
  if(!is.null(plink_path)){
    plink_paths <- plink_path
    if(os_type == "windows"){
      plink_paths <- paste(plink_paths, "plink.exe", sep = "/")
    }

  } else{
  plink_paths <- system(paste(system_command, "plink"), intern = TRUE)
}
  if (length(plink_paths) != 0) {
    versions <- sapply(plink_paths, function(path) {
      version_info <- system(paste(shQuote(path), "--version"), intern = TRUE, ignore.stderr = TRUE)
      match <- character(0)  # Default empty character vector if no match is found

      if(length(grep("v", version_info)) >= 1) {
        match <- regmatches(version_info, regexpr("v([0-9]+\\.[0-9]+b?[0-9]*\\.?[0-9]*)", version_info))
      } else {
        match <- regmatches(version_info, regexpr("Release\\s+([0-9]+\\.[0-9]+)", version_info))
      }
    })

    if (length(versions) > 0) {
      # Convert to numeric version for comparison
      numeric_versions <- as.numeric(sapply(versions, function(version) {
        version <- version[version != ""]
        as.numeric(regmatches(version, regexpr("[0-9]+\\.[0-9]+", version)))
      }))
      valid_indices <- which(numeric_versions >= as.numeric(min_version))

      if (length(valid_indices) > 0) {
        # Extract the numeric versions that are valid
        valid_versions <- numeric_versions[valid_indices]

        # Find the index of the highest version among the valid ones
        max_version_index <- valid_indices[which.max(valid_versions)]
        # Identify the path with the max version
        max_version_path <- plink_paths[max_version_index]

        #max_version_path <- plink_paths[which.max(numeric_versions)]
        # Step 2: Check versions and select
        if (length(max_version_path) == 0) {
          # Define PLINK URLs for different platforms
          if (os_type == "unix" || tolower(Sys.info()["sysname"]) == "linux") {
            #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210328.zip"
            plink_url <- "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
          } else if (os_type == "windows") {
            #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_win64_20210328.zip"
            plink_url <-  "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip"
          }

          if (is.null(plink_url)) {
            stop("Unsupported operating system.")
          }

          # Specify the directory where PLINK will be installed (e.g., within the package directory)
          #package_dir <- system.file(package = "YourPackageName")
          #package_dir <- system.file(package = "PredictProR")
          package_dir <- "D:/PredictProR"
          plink_dir <- file.path(package_dir, "plink")
          if (!dir.exists(plink_dir)) {
            dir.create(plink_dir)
          }

          # Define the temporary file path and download PLINK
          temp_file <- tempfile()
          download.file(plink_url, temp_file, mode = "wb")

          # Unzip PLINK
          unzip(temp_file, exdir = plink_dir)

          # Set PLINK binary to be executable (primarily for Unix-like systems)
          if (os_type == "unix") {
            system(paste("chmod +x", file.path(plink_dir, "plink")))
          }

          # Cleanup the temporary file
          unlink(temp_file)

          message("PLINK has been successfully installed at ", plink_dir)

          return(plink_dir)
        } else {
          plink_dir <- max_version_path
          message(sprintf("Using PLINK at '%s', version %s", valid_indices[which.max(numeric_versions)]))
          return(plink_dir)
        }
      } else{
        # Define PLINK URLs for different platforms
        if (os_type == "unix" || tolower(Sys.info()["sysname"]) == "linux") {
          #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210328.zip"
          plink_url <- "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
        } else if (os_type == "windows") {
          #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_win64_20210328.zip"
          plink_url <-  "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip"
        }

        if (is.null(plink_url)) {
          stop("Unsupported operating system.")
        }

        # Specify the directory where PLINK will be installed (e.g., within the package directory)
        #package_dir <- system.file(package = "YourPackageName")
        #package_dir <- system.file(package = "PredictProR")
        package_dir <- "D:/PredictProR"
        plink_dir <- file.path(package_dir, "plink")
        if (!dir.exists(plink_dir)) {
          dir.create(plink_dir)
        }

        # Define the temporary file path and download PLINK
        temp_file <- tempfile()
        download.file(plink_url, temp_file, mode = "wb")

        # Unzip PLINK
        unzip(temp_file, exdir = plink_dir)

        # Set PLINK binary to be executable (primarily for Unix-like systems)
        if (os_type == "unix") {
          system(paste("chmod +x", file.path(plink_dir, "plink")))
        }

        # Cleanup the temporary file
        unlink(temp_file)

        message("PLINK has been successfully installed at ", plink_dir)

        return(plink_dir)
      }

    }

  } else{
    # Define PLINK URLs for different platforms
    if (os_type == "unix" || tolower(Sys.info()["sysname"]) == "linux") {
      #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210328.zip"
      plink_url <- "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
    } else if (os_type == "windows") {
      #plink_url <- "http://s3.amazonaws.com/plink1-assets/plink_win64_20210328.zip"
      plink_url <-  "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip"
    }

    if (is.null(plink_url)) {
      stop("Unsupported operating system.")
    }

    # Specify the directory where PLINK will be installed (e.g., within the package directory)
    #package_dir <- system.file(package = "YourPackageName")
    #package_dir <- system.file(package = "PredictProR")
    package_dir <- "D:/PredictProR"
    plink_dir <- file.path(package_dir, "plink")
    if (!dir.exists(plink_dir)) {
      dir.create(plink_dir)
    }

    # Define the temporary file path and download PLINK
    temp_file <- tempfile()
    download.file(plink_url, temp_file, mode = "wb")

    # Unzip PLINK
    unzip(temp_file, exdir = plink_dir)

    # Set PLINK binary to be executable (primarily for Unix-like systems)
    if (os_type == "unix") {
      system(paste("chmod +x", file.path(plink_dir, "plink")))
    }

    # Cleanup the temporary file
    unlink(temp_file)

    message("PLINK has been successfully installed at ", plink_dir)

    return(plink_dir)
  }


}
