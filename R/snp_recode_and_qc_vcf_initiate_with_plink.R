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
                         output_format = "vcf.gz",
                         vcf_file_path = NULL,
                         min_version = "2.0",
                         remove_monomorphic = TRUE,
                         maf_threshold = 0.05,
                         heterozygosity = NULL,
                         snp_call_rate = 0.1,
                         allow_extra_chr = FALSE,
                         individual_call_rate = NULL) {

  msg <- "\n==================================================\n"
  # Ensure vcf_file_path is provided
  if (is.null(vcf_file_path)) {
    stop(paste(msg, "vcf_file_path must be provided."), call. = FALSE)
  }

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

  # Construct the base PLINK2 command according to the file extension
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

  # Add output format
  output_cmd <- switch(output_format,
                       "vcf.gz" = "--export vcf bgz",
                       "vcf" = "--export vcf",
                       "bed" = "--make-bed",
                       stop(paste(msg, "Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."), call. = FALSE))
  plink_cmd <- sprintf('%s %s --out "%s"', plink_cmd, output_cmd, output_name)

  # Execute the command
  message("Executing command: ", plink_cmd)
  output <- system(plink_cmd, intern = TRUE)

  # Check if the output contains any errors
  if (length(grep("Error", output, ignore.case = TRUE)) > 0) {
    message(paste(msg, "PLINK2 encountered an error during QC."))
    print(output)
    stop(paste(msg, "PLINK2 QC failed."), call. = FALSE)
  } else {
    message(paste(msg,"PLINK2 QC completed successfully."))
  }

  # Extract and log the PLINK2 version
  version_line <- grep("PLINK v2", output, value = TRUE)
  if (length(version_line) > 0) {
    version <- regmatches(version_line, regexpr("v[0-9]+\\.[0-9]+[a-z]*", version_line))[1]
    message(paste(msg, "PLINK2 version: ", version))
  }

  message(paste(msg, "PLINK2 QC and recoding complete. Output files are prefixed with '", output_name, "'"))
}

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
