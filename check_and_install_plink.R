#
# # Function to remove existing PLINK versions
# remove_existing_plink <- function(plink_paths) {
#   for (path in plink_paths) {
#     tryCatch({
#       file.remove(path)
#       message(sprintf("Removed PLINK at '%s'", path))
#     }, error = function(e) {
#       message(sprintf("Failed to remove PLINK at '%s'", path))
#     })
#   }
# }
#
# # Helper function to get PLINK download URL based on OS
# get_plink_download_url <- function() {
#   msg <- "\n==================================================\n"
#   if (tolower(Sys.info()["sysname"]) %in% c("linux", "unix")) {
#     return("https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip")
#   } else if (.Platform$OS.type == "windows") {
#     return("https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip")
#   } else {
#
#     stop(paste(msg, "Unsupported operating system."), call. = FALSE)
#   }
# }
#
# install_plink <- function(vcf_file_path) {
#   plink_url <- get_plink_download_url()
#   #package_dir <- "D:/PredictProR" # Consider making this an argument or dynamically determining the path
#   package_dir <- vcf_file_path
#   plink_dir <- file.path(package_dir, "plink")
#
#   if (!dir.exists(plink_dir)) {
#     dir.create(plink_dir, recursive = TRUE)
#   }
#
#   temp_file <- tempfile()
#   download.file(plink_url, temp_file, mode = "wb")
#   unzip(temp_file, exdir = plink_dir)
#   unlink(temp_file)
#
#   if (.Platform$OS.type != "windows") {
#     system(paste("chmod +x", file.path(plink_dir, "plink")))
#   }
#
#   message("PLINK has been successfully installed at ", plink_dir)
#   return(plink_dir)
# }
#
#
# check_and_install_plink <- function(min_version = "1.90", vcf_file_path) {
#
#   # Define the command to find PLINK based on OS
#   system_command <- ifelse(.Platform$OS.type == "windows", "where", "which -a")
#   plink_paths <- tryCatch(system(paste(system_command, "plink"), intern = TRUE), error = function(e) character(0))
#
#   if (length(plink_paths) == 0) {
#     plink_dir <- install_plink(vcf_file_path)
#     return(plink_dir)
#   }
#
#   # Extract and check versions
#   versions <- lapply(plink_paths, function(path) {
#     version_info <- tryCatch(system(paste(shQuote(path), "--version"), intern = TRUE, ignore.stderr = TRUE), error = function(e) character(0))
#     if (length(version_info) > 0) {
#       regmatches(version_info, regexpr("([0-9]+\\.[0-9]+)", version_info))[1]
#     } else {
#       NA
#     }
#   })
#
#   valid_versions <- sapply(versions, function(v) { if (!is.na(v)) as.numeric(v) else NA })
#   valid_indices <- which(!is.na(valid_versions) & valid_versions >= as.numeric(min_version))
#
#   if (length(valid_indices) > 0) {
#     # Use the highest valid version
#     highest_version_index <- valid_indices[which.max(valid_versions[valid_indices])]
#     ### This part is added to ensure that only the verison that the command line is
#     ## written is used. In the future we can accomodate version plink2 if stable
#     ## version is available
#     # Remove lower versions
#     lower_versions <- setdiff(plink_paths, plink_paths[highest_version_index])
#     remove_existing_plink(lower_versions)
#
#     plink_dir <- install_plink(vcf_file_path)
#    # plink_dir <- dirname(plink_paths[highest_version_index])
#     #message(sprintf("Using PLINK at '%s'", plink_paths[highest_version_index]))
#
#   } else {
#     remove_existing_plink(plink_paths)
#     plink_dir <- install_plink(vcf_file_path)
#   }
#
#   return(plink_dir)
# }
#
