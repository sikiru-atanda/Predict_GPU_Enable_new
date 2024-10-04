check_phasing <- function(vcf_file) {
  # Read the VCF file
  vcf_lines <- readLines(vcf_file)

  # Filter out comment lines (those starting with #)
  vcf_data <- vcf_lines[!grepl("^#", vcf_lines)]

  # Check a few rows for phasing (first 100 rows in this example)
  phased <- FALSE
  for (line in head(vcf_data, 100)) {
    # Split the genotype columns
    genotype_data <- unlist(strsplit(line, "\t"))[-(1:9)]  # First 9 columns are non-genotype fields
    if (any(grepl("\\|", genotype_data))) {
      phased <- TRUE
      break
    }
  }

  return(phased)
  # if (phased) {
  #   return(TRUE)
  #   #cat("Data appears to be phased.\n")
  # } else {
  #   return(FALSE)
  #   #cat("Data appears to be unphased.\n")
  # }
}

# Example usage



# Function to remove existing PLINK versions
remove_existing_plink <- function(plink_paths) {

  msg <- "\n==================================================\n"

  for (path in plink_paths) {
    tryCatch({
      if (file.exists(path)) {
        file.remove(path)
        message(paste(msg, sprintf("Removed PLINK at '%s'", path)))
      }
    }, error = function(e) {
      message(paste(msg, sprintf("Failed to remove PLINK at '%s' due to: %s", path, e$message)))
    })
  }
}

# Helper function to get PLINK 2 download URL based on OS
get_plink_download_url <- function() {
  msg <- "\n==================================================\n"
  if (tolower(Sys.info()["sysname"]) %in% c("linux", "unix")) {
    return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_amd_avx2_20240625.zip")
  } else if (.Platform$OS.type == "windows") {
    return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_win64_20240625.zip")
  } else if (tolower(Sys.info()["sysname"]) == "darwin") {
    return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_mac_avx2_20240625.zip")

  } else {
    stop(paste(msg, "Unsupported operating system.", call. = FALSE), call. = FALSE)
  }
}

## installing plink if not found
install_plink <- function(vcf_file_path) {
  msg <- "\n==================================================\n"
  plink_url <- get_plink_download_url()
  plink_dir <- file.path(vcf_file_path, "plink2")

  if (!dir.exists(plink_dir)) {
    #message(paste(msg, "Plink2 was not found, thus we initiate installation."))
    dir.create(plink_dir, recursive = TRUE)
  }

  temp_file <- tempfile()
  download.file(plink_url, temp_file, mode = "wb")

  unzip(temp_file, exdir = plink_dir)
  unlink(temp_file)

  if (.Platform$OS.type != "windows") {
    system(paste("chmod +x", file.path(plink_dir, "plink2")))
  }

  if (!file.exists(file.path(plink_dir, "plink2.exe")) && !file.exists(file.path(plink_dir, "plink2"))) {
    stop(paste(msg, "PLINK2 executable not found in the unzipped directory."), call. = FALSE)
  }

  message(paste(msg, "PLINK2 has been successfully installed at ", plink_dir))
  return(plink_dir)
}

check_and_install_plink <- function(min_version = "2.0", vcf_file_path) {

  msg <- "\n==================================================\n"
  # Manually set the plink_dir
  plink_dir <- file.path(vcf_file_path, "plink2")
  plink_path <- if (.Platform$OS.type == "windows") {
    file.path(plink_dir, "plink2.exe")
  } else {
    file.path(plink_dir, "plink2")
  }

  if (file.exists(plink_path)) {
    return(plink_dir)
  }

  message(paste(msg, "Plink2 was not found, we initiate installation."))
  plink_dir <- install_plink(vcf_file_path)

  return(plink_dir)
}

run_plink_qc <- function(input_file,
                         output_name,
                         output_format = "vcf",
                         vcf_file_path = NULL,
                         min_version = "2.0",
                         remove_monomorphic = TRUE,
                         maf_threshold = 0.05,
                         heterozygosity = NULL,
                         snp_call_rate = 0.1,
                         allow_extra_chr = FALSE,
                         individual_call_rate = NULL,
                         ld_pruning = FALSE,         # LD pruning option
                         ld_pruning_method = "indep-pairwise", # LD pruning method
                         window_size = 50,           # Window size for LD pruning
                         step_size = 5,              # Step size for LD pruning
                         r2_threshold = 0.2,         # r² threshold for LD pruning
                         use_kb_window = FALSE,      # Use kb for window size in LD pruning
                         phased = FALSE,             # Option for phased LD pruning
                         use_founders = FALSE        # Option to only consider founders
) {

  msg <- "\n==================================================\n"

  tryCatch({
    # Ensure vcf_file_path is provided
    if (is.null(vcf_file_path)) {
      stop(paste(msg, "vcf_file_path must be provided."), call. = FALSE)

    }

    setwd(vcf_file_path)
    # Adjust snp_call_rate if necessary
    if(!is.null(snp_call_rate)){
      if (snp_call_rate > 0.5) {
        snp_call_rate <- 1 - snp_call_rate
      }
    }

    # Check if PLINK2 is installed and accessible
    plink_dir <- check_and_install_plink(min_version, vcf_file_path)

    # Determine the PLINK executable path based on the operating system
    plink_path <- if (.Platform$OS.type == "windows") {
      file.path(plink_dir, "plink2.exe")
    } else {
      file.path(plink_dir, "plink2")
    }

    # Ensure full path to input VCF file
    input_file_full_path <- file.path(getwd(), input_file)

    # Construct the base PLINK2 command using --vcf (since you're working with a VCF file)
    plink_cmd <- sprintf('%s --vcf "%s"', plink_path, input_file_full_path)

    # Add the --allow-extra-chr flag if specified
    if (isTRUE(allow_extra_chr)) {
      plink_cmd <- sprintf('%s --allow-extra-chr', plink_cmd)
    }

    # Add the --min-alleles and --max-alleles flags to remove monomorphic markers
    if (isTRUE(remove_monomorphic)) {
      plink_cmd <- sprintf('%s --min-alleles 2 --max-alleles 2', plink_cmd)
    }

    # Add the --geno flag for SNP call rate threshold
    if (!is.null(snp_call_rate)) {
      plink_cmd <- sprintf('%s --geno %f', plink_cmd, snp_call_rate)
    }

    # Add the --maf flag for MAF threshold
    if (!is.null(maf_threshold)) {
      plink_cmd <- sprintf('%s --maf %f', plink_cmd, maf_threshold)
    }

    # Add individual call rate threshold if specified
    if (!is.null(individual_call_rate)) {
      plink_cmd <- sprintf('%s --mind %f', plink_cmd, individual_call_rate)
    }

    # If founders should be used, add the --make-founders flag
    if (isTRUE(use_founders)) {
      plink_cmd <- sprintf('%s --make-founders', plink_cmd)
    }

    # Handle LD pruning directly within the same process
    if (isTRUE(ld_pruning)) {
      message(paste(msg, "Performing LD pruning..."))

      # Adjust step size to 1 if using kb window (required by PLINK2)
      if (use_kb_window) {
        step_size <- 1
        window_param <- paste0(window_size, "kb")
      } else {
        window_param <- window_size
      }

      # Add the LD pruning options based on the selected method
      if (ld_pruning_method == "indep-pairwise") {
        plink_cmd <- sprintf('%s --indep-pairwise %s %d %f', plink_cmd, window_param, step_size, r2_threshold)
      } else if (ld_pruning_method == "indep-pairphase") {
        phase_check <- check_phasing(input_file)
        if(phase_check==FALSE){
          stop("snps need to be phased to use indep-pairphase. Use indep-pairwise instead or phase using BEAGLE")
        }
        plink_cmd <- sprintf('%s --indep-pairphase %s %d %f', plink_cmd, window_param, step_size, r2_threshold)
      } else {
        stop(paste(msg, "Unsupported LD pruning method: ", ld_pruning_method), call. = FALSE)
      }

      plink_cmd <- sprintf('%s --out "%s_pruned"', plink_cmd, output_name)

    } else {
      # Add output format if not pruning
      output_cmd <- switch(output_format,
                           "vcf.gz" = "--export vcf bgz",
                           "vcf" = "--export vcf",
                           "bed" = "--make-bed",
                           stop(paste(msg, "Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."), call. = FALSE))
      plink_cmd <- sprintf('%s %s --out "%s"', plink_cmd, output_cmd, output_name)
    }

    # Execute the combined QC and optional LD pruning command
    message("Executing command: ", plink_cmd)
    output <- system(plink_cmd, intern = TRUE)

    # Check for errors in the output
    if (length(grep("Error", output, ignore.case = TRUE)) > 0) {
      message(paste(msg, "PLINK2 encountered an error during QC or LD pruning."))
      print(output)
      stop(paste(msg, "PLINK2 QC/LD pruning failed."), call. = FALSE)
    } else {
      message(paste(msg, "PLINK2 QC and optional LD pruning completed successfully."))
    }

    return("Success")

  }, error = function(e) {
    message("An error occurred: ", e$message)
    return("Failed")
  })
}

# Example usage
# run_plink_qc(input_file = "2021_NDSU_AYT.vcf",
#              output_name = "pea_whole_genome_indep",
#              output_format = "vcf",
#              vcf_file_path = "D:/PredictProR_no_use_file",
#              remove_monomorphic = TRUE,
#              maf_threshold = 0.05,
#              heterozygosity = NULL,
#              snp_call_rate = 0.1,
#              allow_extra_chr = TRUE,
#              individual_call_rate = NULL,
#              ld_pruning = TRUE,         # LD pruning option
#              ld_pruning_method = "indep-pairphase", # LD pruning method
#              window_size = 50,           # Window size for LD pruning
#              step_size = 5,              # Step size for LD pruning
#              r2_threshold = 0.2,         # r² threshold for LD pruning
#              use_kb_window = TRUE,      # Use kb for window size in LD pruning
#              phased = TRUE,             # Option for phased LD pruning
#              use_founders = FALSE )
#


# # Function to remove existing PLINK versions
# remove_existing_plink <- function(plink_paths) {
#
#   msg <- "\n==================================================\n"
#
#   for (path in plink_paths) {
#     tryCatch({
#       if (file.exists(path)) {
#         file.remove(path)
#         message(paste(msg, sprintf("Removed PLINK at '%s'", path)))
#       }
#     }, error = function(e) {
#       message(paste(msg, sprintf("Failed to remove PLINK at '%s' due to: %s", path, e$message)))
#     })
#   }
# }
#
# # Helper function to get PLINK 2 download URL based on OS
# get_plink_download_url <- function() {
#   msg <- "\n==================================================\n"
#   if (tolower(Sys.info()["sysname"]) %in% c("linux", "unix")) {
#     return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_amd_avx2_20240625.zip")
#   } else if (.Platform$OS.type == "windows") {
#     return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_win64_20240625.zip")
#   } else if (tolower(Sys.info()["sysname"]) == "darwin") {
#     return("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_mac_avx2_20240625.zip")
#
#   } else {
#     stop(paste(msg, "Unsupported operating system.", call. = FALSE), call. = FALSE)
#   }
# }
#
# ## installing plink if not found
# install_plink <- function(vcf_file_path) {
#   msg <- "\n==================================================\n"
#   plink_url <- get_plink_download_url()
#   plink_dir <- file.path(vcf_file_path, "plink2")
#
#   if (!dir.exists(plink_dir)) {
#     #message(paste(msg, "Plink2 was not found, thus we initiate installation."))
#     dir.create(plink_dir, recursive = TRUE)
#   }
#
#   temp_file <- tempfile()
#   download.file(plink_url, temp_file, mode = "wb")
#
#   unzip(temp_file, exdir = plink_dir)
#   unlink(temp_file)
#
#   if (.Platform$OS.type != "windows") {
#     system(paste("chmod +x", file.path(plink_dir, "plink2")))
#   }
#
#   if (!file.exists(file.path(plink_dir, "plink2.exe")) && !file.exists(file.path(plink_dir, "plink2"))) {
#     stop(paste(msg, "PLINK2 executable not found in the unzipped directory."), call. = FALSE)
#   }
#
#   message(paste(msg, "PLINK2 has been successfully installed at ", plink_dir))
#   return(plink_dir)
# }
#
# check_and_install_plink <- function(min_version = "2.0", vcf_file_path) {
#
#   msg <- "\n==================================================\n"
#   # Manually set the plink_dir
#   plink_dir <- file.path(vcf_file_path, "plink2")
#   plink_path <- if (.Platform$OS.type == "windows") {
#     file.path(plink_dir, "plink2.exe")
#   } else {
#     file.path(plink_dir, "plink2")
#   }
#
#   if (file.exists(plink_path)) {
#     return(plink_dir)
#   }
#
#   message(paste(msg, "Plink2 was not found, we initiate installation."))
#   plink_dir <- install_plink(vcf_file_path)
#
#   return(plink_dir)
# }
#
# run_plink_qc <- function(input_file,
#                          output_name,
#                          output_format = "vcf",
#                          vcf_file_path = NULL,
#                          min_version = "2.0",
#                          remove_monomorphic = TRUE,
#                          maf_threshold = 0.05,
#                          heterozygosity = NULL,
#                          snp_call_rate = 0.1,
#                          allow_extra_chr = FALSE,
#                          individual_call_rate = NULL) {
#
#   msg <- "\n==================================================\n"
#
#   tryCatch({
#   # Ensure vcf_file_path is provided
#   if (is.null(vcf_file_path)) {
#     stop(paste(msg, "vcf_file_path must be provided."), call. = FALSE)
#   }
#
#   # Adjust snp_call_rate if necessary
#   if(!is.null(snp_call_rate)){
#   if (snp_call_rate > 0.5) {
#     snp_call_rate <- 1 - snp_call_rate
#   }
#
#   }
#
#   # Check if PLINK2 is installed and accessible
#   plink_dir <- check_and_install_plink(min_version, vcf_file_path)
#
#   # Determine the PLINK executable path based on the operating system
#   plink_path <- if (.Platform$OS.type == "windows") {
#     file.path(plink_dir, "plink2.exe")
#   } else {
#     file.path(plink_dir, "plink2")
#   }
#
#   # Ensure full path to input VCF file
#   input_file_full_path <- file.path(getwd(), input_file)
#
#   # Construct the base PLINK2 command according to the file extension
#   plink_cmd <- sprintf('%s --vcf "%s"', plink_path, input_file_full_path)
#
#   # Add the --allow-extra-chr flag if specified
#   if (isTRUE(allow_extra_chr)) {
#     plink_cmd <- sprintf('%s --allow-extra-chr', plink_cmd)
#   }
#
#   # Add the --min-alleles and --max-alleles flags to remove monomorphic markers
#   if (isTRUE(remove_monomorphic)) {
#     plink_cmd <- sprintf('%s --min-alleles 2 --max-alleles 2', plink_cmd)
#   }
#
#   # Add the --geno flag for SNP call rate threshold
#   if (!is.null(snp_call_rate)) {
#     plink_cmd <- sprintf('%s --geno %f', plink_cmd, snp_call_rate)
#   }
#
#   # Add the --maf flag for MAF threshold
#   if (!is.null(maf_threshold)) {
#     plink_cmd <- sprintf('%s --maf %f', plink_cmd, maf_threshold)
#   }
#
#   # Add individual call rate threshold if specified
#   if (!is.null(individual_call_rate)) {
#     plink_cmd <- sprintf('%s --mind %f', plink_cmd, individual_call_rate)
#   }
#
#   # Add output format
#   output_cmd <- switch(output_format,
#                        "vcf.gz" = "--export vcf bgz",
#                        "vcf" = "--export vcf",
#                        "bed" = "--make-bed",
#                        stop(paste(msg, "Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."), call. = FALSE))
#   plink_cmd <- sprintf('%s %s --out "%s"', plink_cmd, output_cmd, output_name)
#
#   # Execute the command
#   message("Executing command: ", plink_cmd)
#   output <- system(plink_cmd, intern = TRUE)
#
#   # Check if the output contains any errors
#   if (length(grep("Error", output, ignore.case = TRUE)) > 0) {
#     message(paste(msg, "PLINK2 encountered an error during QC."))
#     print(output)
#     stop(paste(msg, "PLINK2 QC failed."), call. = FALSE)
#   } else {
#     message(paste(msg,"PLINK2 QC completed successfully."))
#   }
#
#   # Extract and log the PLINK2 version
#   version_line <- grep("PLINK v2", output, value = TRUE)
#   if (length(version_line) > 0) {
#     version <- regmatches(version_line, regexpr("v[0-9]+\\.[0-9]+[a-z]*", version_line))[1]
#     message(paste(msg, "PLINK2 version: ", version))
#   }
#
#   message(paste(msg, "PLINK2 QC and recoding complete. Output files are prefixed with '", output_name, "'"))
#
#   return("Success")
#
# }, error = function(e) {
#   message("An error occurred: ", e$message)
#   return("Failed")
# })
#
# }

# Example usage
# run_plink_qc(input_file = "initial_qc_data.vcf.gz",
#              output_name = "pea_whole_genome",
#              output_format = "vcf",
#              vcf_file_path = "D:/RICA",
#              remove_monomorphic = TRUE,
#              maf_threshold = 0.05,
#              heterozygosity = NULL,
#              snp_call_rate = 0.1,
#              allow_extra_chr = TRUE,
#              individual_call_rate = NULL)
