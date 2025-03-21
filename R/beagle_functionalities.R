#' @import data.table
#' @importFrom data.table ., .SD, :=

# Function to dynamically install Java
check_and_install_java <- function() {
  java_check <- system("java -version", intern = TRUE, ignore.stderr = TRUE)

  if (length(java_check) == 0) {
    message("Java is not installed. Installing Java...")

    install_java <- function() {
      os_type <- Sys.info()["sysname"]

      if (os_type == "Linux") {
        message("Installing Java on Linux...")
        system("sudo apt update", wait = TRUE)
        system("sudo apt install -y default-jre", wait = TRUE)
      } else if (os_type == "Darwin") {
        message("Installing Java on macOS using Homebrew...")
        system("brew install openjdk@11", wait = TRUE)
      } else if (os_type == "Windows") {
        message("Installing Java on Windows using Chocolatey...")
        system("choco install jdk11", wait = TRUE)
      } else {
        stop("Unsupported operating system. Please install Java manually.")
      }

      system("java -version", wait = TRUE)
    }

    install_java()
  } else {
    message("Java is already installed.")
  }
}

# Function to download Beagle
download_beagle <- function(output_dir = getwd(), beagle_version = "5.4") {
  beagle_url <- "https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar"
  beagle_jar_path <- file.path(output_dir, paste0("beagle.", beagle_version, ".jar"))

  if (!file.exists(beagle_jar_path)) {
    tryCatch({
      download.file(beagle_url, destfile = beagle_jar_path, mode = "wb")
      message("Beagle version 06Aug24.a91 downloaded successfully.")
    }, error = function(e) {
      message("Failed to download Beagle: ", e$message)
    })
  } else {
    message("Beagle is already downloaded.")
  }

  return(beagle_jar_path)
}


# Function to generate a new filename by inserting 'res' before the file extension
create_res_output <- function(output_vcf) {
  # Separate the base name and the extension (assumes .vcf is the extension)
  file_parts <- strsplit(output_vcf, split = "\\.vcf")[[1]]

  # Combine the base name with '_res' and add the extension back
  res_output_vcf <- paste0(file_parts[1], "_res.vcf")

  return(res_output_vcf)
}

# Function to filter and dynamically remove SNPs with identical REF and ALT or missing values in REF/ALT
filter_vcf <- function(input_vcf, output_vcf) {

  # Read the VCF metadata (lines starting with '##')
  # vcf_metadata <- data.table::fread(input_vcf, header = FALSE, sep = "\n", quote = "",
  #                                   nrows = grep("#CHROM", readLines(input_vcf)) - 1)

  vcf_metadata <- data.table::fread(input_vcf, skip = "#CHROM", header = TRUE)[, 1:9]


  # Read the VCF data (starting from the #CHROM line onwards)
  vcf_data <- data.table::fread(input_vcf, skip = "#CHROM", header = TRUE)

  # Define missing values (e.g., NA, ".", or any other missing indication)
  missing_values <- c(NA, ".", "")

  # Identify SNPs where either REF == ALT or where REF or ALT is missing
  #problematic_snps <- vcf_data[REF == ALT | REF %in% missing_values | ALT %in% missing_values]
  problematic_snps <- vcf_data[
    vcf_data[["REF"]] == vcf_data[["ALT"]] |
      vcf_data[["REF"]] %in% missing_values |
      vcf_data[["ALT"]] %in% missing_values
  ]
  # If there are problematic SNPs, filter and write the new VCF
  if (nrow(problematic_snps) > 0) {
    # Filter out SNPs where REF == ALT or where REF or ALT is missing
    filtered_vcf <- vcf_data[REF != ALT & !(REF %in% missing_values) & !(ALT %in% missing_values)]

    # Write the filtered VCF back to file
    data.table::fwrite(vcf_metadata, output_vcf, col.names = FALSE, quote = FALSE)  # Write metadata
    data.table::fwrite(filtered_vcf, output_vcf, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)  # Write data

    # Print result
    cat("Filtered VCF file saved as:", output_vcf, "\n")
    cat("Number of SNPs removed:", nrow(vcf_data) - nrow(filtered_vcf), "\n")

    return(TRUE)
  } else {
    # No problematic SNPs, so no file is written
    cat("No problematic SNPs found. No file was written.\n")
    return(FALSE)
  }
}

# Function to extrapolate a genetic map if not provided
generate_genetic_map <- function(vcf_file) {
  message("Generating a simple extrapolated genetic map...")

  # Read the VCF file to extract CHROM and POS columns
  vcf_data <- data.table::fread(vcf_file, select = c("#CHROM", "POS"))

  # Generate a placeholder variant ID (e.g., chr:pos)
  vcf_data[, variant_id := paste0(`#CHROM`, ":", POS)]

  # Generate genetic positions assuming 1 cM per 1 Mb
  vcf_data <- vcf_data[, .(chromosome = `#CHROM`, variant_id, cM = POS / 1e6, position = POS)]

  # Write the extrapolated genetic map to a file without column headers
  map_file <- file.path(getwd(), "extrapolated_genetic_map.map")

  # Write the data without headers (Beagle does not expect headers)
  data.table::fwrite(vcf_data, map_file, sep = "\t", col.names = FALSE)

  return(map_file)
}

# Function to check and clean CHROM and POS columns in a VCF file
clean_vcf_chrom_pos <- function(input_vcf, output_vcf) {

  # Read the first 100 lines (or fewer if the file is smaller)
  vcf_sample <- data.table::fread(input_vcf, skip = "#CHROM", header = TRUE, nrows = 100)


  # Extract CHROM and POS columns from the sample
  chrom_sample <- as.character(vcf_sample$`#CHROM`)
  pos_sample <- vcf_sample$POS

  # Function to extract numeric part from CHROM
  clean_chrom <- function(chrom_value) {
    # Look for 'CHR' and extract the number that follows
    matches <- regmatches(chrom_value, regexpr("CHR[0-9]+", chrom_value, ignore.case = TRUE))

    if (length(matches) > 0) {
      # Extract numeric portion after 'CHR'
      return(as.numeric(gsub("[^0-9]", "", matches)))
    } else {
      # If no 'CHR' is found, return the numeric portion from the original value
      return(as.numeric(gsub("[^0-9]", "", chrom_value)))
    }
  }

  # Apply the clean_chrom function to the CHROM sample column
  cleaned_chrom_sample <- sapply(chrom_sample, clean_chrom)

  # Check if there are any non-numeric values in the CHROM or POS columns in the sample
  if(length(cleaned_chrom_sample)>0){
    chrom_issue <- any(is.na(cleaned_chrom_sample))
    pos_issue <- any(is.na(pos_sample) | !is.numeric(pos_sample))

    if (!chrom_issue && !pos_issue) {
      #cat("No issues found in the first 100 lines. No file will be written.\n")

    } else {
      stop("Provide numeric CHROM and POS.", call. = FALSE)
    }
  } else {
    return(FALSE)
  }


  # Read the full VCF metadata (lines starting with '##')
  # vcf_metadata <- data.table::fread(input_vcf, header = FALSE, sep = "\n", quote = "",
  #                                   nrows = grep("#CHROM", readLines(input_vcf)) - 1)

  vcf_metadata <- data.table::fread(input_vcf, skip = "#CHROM", header = TRUE)[, 1:9]

  # Read the full VCF data (starting from the #CHROM line onwards)
  vcf_data <- data.table::fread(input_vcf, skip = "#CHROM", header = TRUE)

  # Extract CHROM and POS columns
  chrom <- as.character(vcf_data$`#CHROM`)
  pos <- vcf_data$POS

  # Apply the clean_chrom function to the full CHROM column
  cleaned_chrom <- sapply(chrom, clean_chrom)

  # Check if CHROM and POS are numeric
  if (any(is.na(cleaned_chrom))) {
    stop("Error: CHROM contains non-numeric values. Please provide a valid numeric format for CHROM.")
  }

  if (!is.numeric(pos) || any(is.na(pos))) {
    stop("Error: POS contains non-numeric values. Please provide a valid numeric format for POS.")
  }

  # Replace the original CHROM column with the cleaned numeric version
  vcf_data$`#CHROM` <- cleaned_chrom

  # Write the cleaned VCF back to file only if there were issues
  data.table::fwrite(vcf_metadata, output_vcf, col.names = FALSE, quote = FALSE)  # Write metadata
  data.table::fwrite(vcf_data, output_vcf, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)  # Write data

  # Print result
  cat("Cleaned VCF file saved as:", output_vcf, "\n")

  return(TRUE)
}

# Function to validate reference panel file (VCF or bref3 format)
validate_ref_file <- function(ref_file) {
  if (!file.exists(ref_file)) {
    stop("Reference panel file does not exist.")
  }
  ext <- tools::file_ext(ref_file)
  if (!ext %in% c("vcf", "vcf.gz", "bref3")) {
    stop("Reference panel must be in VCF, gzipped VCF, or bref3 format.")
  }
  message("Reference file validated successfully.")
}

# Function to validate markers file (one marker per line)
validate_markers_file <- function(markers_file) {
  if (!file.exists(markers_file)) {
    stop("Markers file does not exist.")
  }
  markers_data <- readLines(markers_file)
  if (length(markers_data) == 0 || !all(nzchar(markers_data))) {
    stop("Markers file is empty or contains invalid lines.")
  }
  message("Markers file validated successfully.")
}

# Function to validate PLINK PED file format
validate_ped_file <- function(ped_file) {
  if (!file.exists(ped_file)) {
    stop("PED file does not exist.")
  }
  ped_data <- data.table::fread(ped_file, header = FALSE)
  if (ncol(ped_data) < 6) {
    stop("PED file must have at least 6 columns: Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype.")
  }
  message("PED file validated successfully.")
}

# Function to validate the map file format
validate_map_file <- function(map_file) {
  # Read the first few rows of the map file
  map_data <- data.table::fread(map_file, header = FALSE, sep = "\t", nrows = 10)

  # Check if the file has exactly 4 columns (chromosome, variant ID, cM, position)
  if (ncol(map_data) != 4) {
    stop("Map file must have exactly 4 columns: chromosome, variant ID, cM, and position (bp).")
  }

  # Check if the first row contains any non-numeric values (e.g., headers)
  if (!is.numeric(as.numeric(map_data[[1]]))) {
    stop("Map file contains a header. Please remove the header and try again.")
  }

  # Validate that cM column is numeric
  if (!all(sapply(map_data[[3]], is.numeric))) {
    stop("The cM column in the map file must contain numeric values.")
  }

  message("Map file validated successfully.")
}

# Function to run Beagle for imputation or phasing
run_beagle <- function(vcf_file, output_prefix, beagle_jar,
                       map_file = NULL, ref_file = NULL, generate_map = TRUE,
                       markers_file = NULL, ped_file = NULL,
                       ibd = FALSE, ne = 10000, nthreads = 4,
                       window = 40.0, overlap = 2.0, burnin = 3,
                       iterations = 12, phase_states = 280,
                       impute = TRUE, imp_states = 1600, imp_segment = 6.0,
                       gp = TRUE, gprobs = FALSE, gt = TRUE
                       ) {

  # Validate inputs
  if (!file.exists(vcf_file)) stop("VCF file does not exist.")

  vcf_file_use <- create_res_output(vcf_file)
  res_filter <- filter_vcf(vcf_file, vcf_file_use)
  if(res_filter==TRUE) vcf_file <- vcf_file_use

  res_chrom_pos <- clean_vcf_chrom_pos(vcf_file, vcf_file_use)
  if(res_chrom_pos==TRUE) vcf_file <- vcf_file_use

  # Check if map_file is provided and validate the format
  if (!is.null(map_file)) {
    if (!file.exists(map_file)) {
      stop("Map file does not exist. Please provide a valid map file.")
    } else {
      validate_map_file(map_file)
    }
  } else if (isTRUE(generate_map)) {
    map_file <- generate_genetic_map(vcf_file)
  }

  # Validate reference panel file if provided
  if (!is.null(ref_file)) {
    validate_ref_file(ref_file)
  }

  # Validate markers file if provided
  if (!is.null(markers_file)) {
    validate_markers_file(markers_file)
  }

  # Validate PED file if provided
  if (!is.null(ped_file)) {
    validate_ped_file(ped_file)
  }


  # Construct the Beagle command
  beagle_cmd <- sprintf("java -Xmx4g -jar %s gt=%s out=%s ne=%d nthreads=%d window=%.1f overlap=%.1f burnin=%d iterations=%d phase-states=%d",
                        beagle_jar, vcf_file, output_prefix, ne, nthreads, window, overlap, burnin, iterations, phase_states)

  # Add optional files
  if (!is.null(map_file)) beagle_cmd <- sprintf("%s map=%s", beagle_cmd, map_file)
  if (!is.null(ref_file)) beagle_cmd <- sprintf("%s ref=%s", beagle_cmd, ref_file)
  if (!is.null(markers_file)) beagle_cmd <- sprintf("%s markers=%s", beagle_cmd, markers_file)
  if (!is.null(ped_file)) beagle_cmd <- sprintf("%s ped=%s", beagle_cmd, ped_file)

  # Add imputation-specific flags
  if (impute) {
    beagle_cmd <- sprintf("%s impute=true imp-states=%d imp-segment=%.1f", beagle_cmd, imp_states, imp_segment)
  } else {
    beagle_cmd <- sprintf("%s impute=false", beagle_cmd)
  }

  # Add IBD detection, GP, and GT options
  #if (ibd) beagle_cmd <- sprintf("%s ibd=true", beagle_cmd)
  if (gp) beagle_cmd <- sprintf("%s gp=true", beagle_cmd)
  if (gprobs) beagle_cmd <- sprintf("%s gprobs=true", beagle_cmd)
  if (!gt) beagle_cmd <- sprintf("%s gt=false", beagle_cmd)



  # Log Beagle execution and capture any potential errors
  log_file <- paste0(output_prefix, "_log.txt")
  message("Running Beagle...")
  result <- system(beagle_cmd, intern = TRUE)

  # Write Beagle output log to a file for troubleshooting
  writeLines(result, log_file)

  # Check if the VCF output file was created
  phased_vcf <- paste0(output_prefix, ".vcf.gz")
  if (file.exists(phased_vcf)) {
    message("Phased VCF file created: ", phased_vcf)
  } else {
    message("Beagle run completed, but no phased VCF file was found.")
    message("Check the log file for more details: ", log_file)
  }

  if (file.exists(vcf_file_use)) {
  unlink(vcf_file_use)

  }
}

# Function to check Java, download Beagle, and run the appropriate Beagle process
impute_vcf_with_beagle <- function(vcf_file, output_prefix = "imputed",
                                   map_file = NULL, ref_file = NULL,
                                   markers_file = NULL, ped_file = NULL,
                                   generate_map = TRUE,
                                   ibd = FALSE, ne = 10000, nthreads = 4,
                                   window = 40.0, overlap = 2.0, burnin = 3,
                                   iterations = 12, phase_states = 280,
                                   impute = TRUE, imp_states = 1600, imp_segment = 6.0,
                                   gp = TRUE, gprobs = FALSE, gt = TRUE,
                                   java_installed = FALSE,
                                   beagle_version = "06Aug24.a91") {

  # Step 1: Check and install Java if needed
  if (!java_installed) {
    check_and_install_java()
  }

  # Step 2: Download Beagle
  beagle_jar <- download_beagle(beagle_version = beagle_version)

  # Step 3: Run Beagle with customizable parameters
  run_beagle(vcf_file = vcf_file, output_prefix = output_prefix,
             generate_map = generate_map,
             beagle_jar = beagle_jar, map_file = map_file,
             ref_file = ref_file, markers_file = markers_file,
             ped_file = ped_file, ibd = ibd, ne = ne,
             nthreads = nthreads, window = window,
             overlap = overlap, burnin = burnin,
             iterations = iterations, phase_states = phase_states,
             impute = impute, imp_states = imp_states,
             imp_segment = imp_segment, gp = gp, gprobs = gprobs, gt = gt)
}

# # Example usage: Phasing only (no imputation)
# impute_vcf_with_beagle("output_vcf_file.vcf",
#                        output_prefix = "impute6_output",
#                        map_file = NULL,
#                        generate_map = FALSE,
#                        impute = TRUE,
#                        #imp_states = 1600, imp_segment = 6.0,
#                        gp= TRUE)
#
#
# impute_vcf_with_beagle(vcf_file = "2021_NDSU_AYT.vcf",
#            output_prefix = "ibd_output",
#            #beagle_jar = "path/to/beagle.jar",
#            #ibd = TRUE,
#            impute = TRUE,
#            gp= TRUE,
#            gt = TRUE)
