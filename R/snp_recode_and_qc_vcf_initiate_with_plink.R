check_phasing <- function(vcf_file) {
  # Read the VCF file
  vcf_lines <- vcf_read_all_lines(vcf_file)

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

  msg <- ""

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
  msg <- ""
  override <- Sys.getenv("PREDICTPRO_PLINK2_URL", unset = "")
  if (nzchar(override)) {
    return(override)
  }

  os_name <- tolower(Sys.info()[["sysname"]] %||% "")
  machine <- tolower(Sys.info()[["machine"]] %||% "")
  if (os_name %in% c("linux", "unix")) {
    return("https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260504.zip")
  } else if (.Platform$OS.type == "windows") {
    return("https://s3.amazonaws.com/plink2-assets/plink2_win64_20260504.zip")
  } else if (os_name == "darwin" && machine %in% c("arm64", "aarch64")) {
    return("https://s3.amazonaws.com/plink2-assets/plink2_mac_arm64_20260504.zip")
  } else if (os_name == "darwin") {
    return("https://s3.amazonaws.com/plink2-assets/plink2_mac_20260504.zip")

  } else {
    stop(paste(msg, "Unsupported operating system.", call. = FALSE), call. = FALSE)
  }
}

predictpror_plink2_cache_dir <- function() {
  root <- Sys.getenv("PREDICTPRO_PLINK2_DIR", unset = "")
  if (!nzchar(root)) {
    root <- file.path(tools::R_user_dir("PredictProR", which = "cache"), "plink2")
  }
  root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) {
    root <- file.path(tempdir(), "PredictProR", "plink2")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

plink2_executable_name <- function() {
  if (.Platform$OS.type == "windows") "plink2.exe" else "plink2"
}

## installing plink if not found
install_plink <- function(vcf_file_path = predictpror_plink2_cache_dir()) {
  msg <- ""
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
    plink_unix <- file.path(plink_dir, "plink2")
    if (file.exists(plink_unix)) {
      Sys.chmod(plink_unix, mode = "0755")
    }
  }

  if (!file.exists(file.path(plink_dir, "plink2.exe")) && !file.exists(file.path(plink_dir, "plink2"))) {
    stop(paste(msg, "PLINK2 executable not found in the unzipped directory."), call. = FALSE)
  }

  message(paste(msg, "PLINK2 has been successfully installed at ", plink_dir))
  return(plink_dir)
}

resolve_plink2_path <- function(vcf_file_path = NULL, install_if_missing = TRUE) {
  exe <- plink2_executable_name()
  env_path <- Sys.getenv("PREDICTPRO_PLINK2_PATH", unset = "")
  if (nzchar(env_path) && !file.exists(env_path)) {
    env_path <- unname(Sys.which(env_path))
  }
  candidates <- c(
    env_path,
    unname(Sys.which("plink2")),
    if (!is.null(vcf_file_path)) file.path(vcf_file_path, "plink2", exe) else character(),
    file.path(predictpror_plink2_cache_dir(), "plink2", exe)
  )
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    return(normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE))
  }
  if (!isTRUE(install_if_missing)) {
    return(NULL)
  }
  plink_dir <- install_plink(predictpror_plink2_cache_dir())
  plink_path <- file.path(plink_dir, exe)
  if (!file.exists(plink_path)) {
    stop("PLINK2 installation completed but the executable was not found.", call. = FALSE)
  }
  normalizePath(plink_path, winslash = "/", mustWork = TRUE)
}

check_and_install_plink <- function(min_version = "2.0", vcf_file_path) {

  dirname(resolve_plink2_path(vcf_file_path = vcf_file_path, install_if_missing = TRUE))
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
                         r2_threshold = 0.2,         # r^2 threshold for LD pruning
                         use_kb_window = FALSE,      # Use kb for window size in LD pruning
                         phased = FALSE,             # Option for phased LD pruning
                         use_founders = FALSE        # Option to only consider founders
) {

  msg <- ""

  tryCatch({
    # Ensure vcf_file_path is provided
    if (is.null(vcf_file_path)) {
      stop(paste(msg, "vcf_file_path must be provided."), call. = FALSE)

    }

    vcf_file_path <- normalizePath(vcf_file_path, winslash = "/", mustWork = TRUE)
    # Adjust snp_call_rate if necessary
    if(!is.null(snp_call_rate)){
      if (snp_call_rate > 0.5) {
        snp_call_rate <- 1 - snp_call_rate
      }
    }

    # Check if PLINK2 is installed and accessible
    plink_path <- resolve_plink2_path(vcf_file_path = vcf_file_path, install_if_missing = TRUE)

    is_absolute_path <- function(path) {
      grepl("^([A-Za-z]:|/|\\\\\\\\|~)", path)
    }

    input_file_path <- if (is_absolute_path(input_file)) {
      path.expand(input_file)
    } else {
      file.path(vcf_file_path, input_file)
    }
    input_file_full_path <- normalizePath(input_file_path, winslash = "/", mustWork = TRUE)
    sanitized_vcf <- sanitize_vcf_for_external_tools(input_file_full_path, mode = "plink")
    if (!identical(sanitized_vcf$path, input_file_full_path)) {
      on.exit(unlink(sanitized_vcf$path), add = TRUE)
      input_file_full_path <- sanitized_vcf$path
    }
    if (isTRUE(sanitized_vcf$changed)) {
      message(
        "VCF preflight for PLINK kept ",
        sanitized_vcf$metrics$output_variants,
        " of ",
        sanitized_vcf$metrics$input_variants,
        " variants."
      )
    }

    output_prefix <- if (is_absolute_path(output_name)) {
      path.expand(output_name)
    } else {
      file.path(vcf_file_path, output_name)
    }

    plink_args <- c("--vcf", input_file_full_path)

    # Add the --allow-extra-chr flag if specified
    if (isTRUE(allow_extra_chr)) {
      plink_args <- c(plink_args, "--allow-extra-chr")
    }

    # Add the --min-alleles and --max-alleles flags to remove monomorphic markers
    if (isTRUE(remove_monomorphic)) {
      plink_args <- c(plink_args, "--min-alleles", "2", "--max-alleles", "2")
    }

    # Add the --geno flag for SNP call rate threshold
    if (!is.null(snp_call_rate)) {
      plink_args <- c(plink_args, "--geno", as.character(snp_call_rate))
    }

    # Add the --maf flag for MAF threshold
    if (!is.null(maf_threshold)) {
      plink_args <- c(plink_args, "--maf", as.character(maf_threshold))
    }

    # Add individual call rate threshold if specified
    if (!is.null(individual_call_rate)) {
      plink_args <- c(plink_args, "--mind", as.character(individual_call_rate))
    }

    # If founders should be used, add the --make-founders flag
    if (isTRUE(use_founders)) {
      plink_args <- c(plink_args, "--make-founders")
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
        plink_args <- c(plink_args, "--indep-pairwise", as.character(window_param), as.character(step_size), as.character(r2_threshold))
      } else if (ld_pruning_method == "indep-pairphase") {
        phase_check <- check_phasing(input_file_full_path)
        if(phase_check==FALSE){
          stop("snps need to be phased to use indep-pairphase. Use indep-pairwise instead or phase using BEAGLE")
        }
        plink_args <- c(plink_args, "--indep-pairphase", as.character(window_param), as.character(step_size), as.character(r2_threshold))
      } else {
        stop(paste(msg, "Unsupported LD pruning method: ", ld_pruning_method), call. = FALSE)
      }

      plink_args <- c(plink_args, "--out", paste0(output_prefix, "_pruned"))

    } else {
      # Add output format if not pruning
      output_args <- switch(output_format,
                            "vcf.gz" = c("--export", "vcf", "bgz"),
                            "vcf" = c("--export", "vcf"),
                            "bed" = c("--make-bed"),
                            stop(paste(msg, "Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."), call. = FALSE))
      plink_args <- c(plink_args, output_args, "--out", output_prefix)
    }

    # Execute the combined QC and optional LD pruning command
    message("Executing PLINK2: ", paste(c(shQuote(plink_path), shQuote(plink_args)), collapse = " "))
    output <- system2(
      plink_path,
      args = gp_quote_system_args(plink_args),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }

    # Check for errors in the output
    if (!identical(as.integer(status), 0L) || length(grep("Error", output, ignore.case = TRUE)) > 0) {
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

