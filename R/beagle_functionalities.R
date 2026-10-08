#' @import data.table
#' @importFrom data.table ., .SD, :=

beagle_system2 <- function(command, args, stdout = TRUE, stderr = TRUE) {
  system2(command, args = args, stdout = stdout, stderr = stderr)
}

# Function to check Java without mutating the user's operating system
check_and_install_java <- function() {
  java_bin <- Sys.getenv("PREDICTPRO_JAVA", unset = "")
  if (nzchar(java_bin) && !file.exists(java_bin)) {
    java_bin <- unname(Sys.which(java_bin))
  }
  if (!nzchar(java_bin)) {
    java_bin <- unname(Sys.which("java"))
  }
  if (!nzchar(java_bin) || !file.exists(java_bin)) {
    stop(
      "Java was not found. Install Java 8 or newer and make `java` available on PATH, ",
      "or set PREDICTPRO_JAVA to the Java executable path.",
      call. = FALSE
    )
  }

  java_check <- beagle_system2(
    java_bin,
    args = gp_quote_system_args("-version"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(java_check, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop("Java was found but `java -version` failed:\n", paste(java_check, collapse = "\n"), call. = FALSE)
  }
  invisible(normalizePath(java_bin, winslash = "/", mustWork = TRUE))
}

predictpror_beagle_cache_dir <- function() {
  root <- Sys.getenv("PREDICTPRO_BEAGLE_DIR", unset = "")
  if (!nzchar(root)) {
    root <- file.path(tools::R_user_dir("PredictProR", which = "cache"), "beagle")
  }
  root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) {
    root <- file.path(tempdir(), "PredictProR", "beagle")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

validate_beagle_jar <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  info <- file.info(path)
  if (is.na(info$size) || info$size < 4L) {
    stop("The Beagle jar is empty or truncated: ", path, call. = FALSE)
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  signature <- readBin(con, what = "raw", n = 2L)
  if (!identical(as.integer(signature), c(0x50L, 0x4bL))) {
    stop("The Beagle jar does not have a valid ZIP/JAR signature: ", path, call. = FALSE)
  }
  path
}

# Function to download Beagle
download_beagle <- function(output_dir = predictpror_beagle_cache_dir(), beagle_version = "29Oct24.c8e") {
  beagle_jar_env <- Sys.getenv("PREDICTPRO_BEAGLE_JAR", unset = "")
  if (nzchar(beagle_jar_env) && file.exists(beagle_jar_env)) {
    return(validate_beagle_jar(beagle_jar_env))
  }

  beagle_url <- Sys.getenv("PREDICTPRO_BEAGLE_URL", unset = "")
  if (!nzchar(beagle_url)) {
    beagle_url <- paste0("https://faculty.washington.edu/browning/beagle/beagle.", beagle_version, ".jar")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  beagle_jar_path <- file.path(output_dir, paste0("beagle.", beagle_version, ".jar"))

  if (!file.exists(beagle_jar_path)) {
    temp_jar <- tempfile("beagle-download-", tmpdir = output_dir, fileext = ".jar")
    on.exit(unlink(temp_jar, force = TRUE), add = TRUE)
    tryCatch({
      utils::download.file(beagle_url, destfile = temp_jar, mode = "wb", quiet = TRUE)
      validate_beagle_jar(temp_jar)
      if (!file.rename(temp_jar, beagle_jar_path)) {
        if (!file.copy(temp_jar, beagle_jar_path, overwrite = FALSE)) {
          stop("could not move the verified download into the Beagle cache")
        }
      }
      message("Beagle downloaded successfully: ", beagle_jar_path)
    }, error = function(e) {
      stop("Failed to download Beagle from ", beagle_url, ": ", e$message, call. = FALSE)
    })
  } else {
    message("Beagle is already downloaded.")
  }

  validate_beagle_jar(beagle_jar_path)
}


# Function to generate a new filename by inserting 'res' before the file extension
create_res_output <- function(output_vcf) {
  output_vcf <- sub("\\.gz$", "", output_vcf, ignore.case = TRUE)
  output_vcf <- sub("\\.vcf$", "", output_vcf, ignore.case = TRUE)
  paste0(output_vcf, "_res.vcf")
}

# Function to filter and dynamically remove SNPs with identical REF and ALT or missing values in REF/ALT
filter_vcf <- function(input_vcf, output_vcf) {
  res <- sanitize_vcf_for_external_tools(
    input_vcf = input_vcf,
    output_vcf = output_vcf,
    mode = "beagle"
  )
  if (isTRUE(res$changed)) {
    removed <- res$metrics$input_variants - res$metrics$output_variants
    message("Filtered VCF file saved as: ", res$path)
    message("Number of variants removed or normalized: ", removed + res$metrics$normalized_fields)
  } else {
    message("No problematic VCF rows found. No file was written.")
  }
  isTRUE(res$changed)
}

# Function to extrapolate a genetic map if not provided
generate_genetic_map <- function(vcf_file, output_dir = dirname(normalizePath(vcf_file, winslash = "/", mustWork = TRUE))) {
  message("Generating a simple extrapolated genetic map...")

  # Read the VCF file to extract CHROM and POS columns
  vcf_data <- vcf_read_table(vcf_file)
  chrom_col <- vcf_find_column(vcf_data, c("#CHROM", "CHROM"))
  pos_col <- vcf_find_column(vcf_data, "POS")
  if (is.null(chrom_col) || is.null(pos_col)) {
    stop("VCF must contain CHROM and POS columns to generate a genetic map.", call. = FALSE)
  }
  pos <- suppressWarnings(as.numeric(vcf_data[[pos_col]]))

  # Generate a placeholder variant ID (e.g., chr:pos)
  variant_id <- paste0(vcf_data[[chrom_col]], ":", vcf_data[[pos_col]])

  # Generate genetic positions assuming 1 cM per 1 Mb
  vcf_data <- data.table::data.table(
    chromosome = vcf_data[[chrom_col]],
    variant_id = variant_id,
    cM = pos / 1e6,
    position = pos
  )

  # Write the extrapolated genetic map to a file without column headers
  map_file <- file.path(output_dir, "extrapolated_genetic_map.map")

  # Write the data without headers (Beagle does not expect headers)
  data.table::fwrite(vcf_data, map_file, sep = "\t", col.names = FALSE)

  return(map_file)
}

# Function to check and clean CHROM and POS columns in a VCF file
clean_vcf_chrom_pos <- function(input_vcf, output_vcf) {
  res <- sanitize_vcf_for_external_tools(
    input_vcf = input_vcf,
    output_vcf = output_vcf,
    mode = "beagle"
  )
  if (isTRUE(res$changed)) {
    message("Cleaned VCF file saved as: ", res$path)
  }
  isTRUE(res$changed)
}

# Function to validate reference panel file (VCF or bref3 format)
validate_ref_file <- function(ref_file) {
  if (!file.exists(ref_file)) {
    stop("Reference panel file does not exist.")
  }
  ref_lower <- tolower(ref_file)
  if (!grepl("\\.(vcf|vcf\\.gz|bref3)$", ref_lower)) {
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

  # cM and bp must be numeric; a header row or text there fails. (Chromosome
  # names may legitimately be text, e.g. "chr1A", so column 1 is not checked.)
  numeric_ok <- function(x) all(!is.na(suppressWarnings(as.numeric(as.character(x)))))
  if (!numeric_ok(map_data[[3]]) || !numeric_ok(map_data[[4]])) {
    stop(
      "Map file columns 3 (cM) and 4 (bp position) must be numeric. ",
      "Remove any header row and check the column order.",
      call. = FALSE
    )
  }

  message("Map file validated successfully.")
}

# Beagle needs positions in ascending order within a chromosome; otherwise
# it fails with an obscure "Window has only one position" error.
gp_beagle_check_sorted_positions <- function(vcf_file) {
  pos <- data.table::fread(
    vcf_file, skip = "#CHROM", header = TRUE, sep = "\t", select = 1:2,
    colClasses = "character", showProgress = FALSE
  )
  if (!nrow(pos)) return(invisible(TRUE))
  chrom <- as.character(pos[[1L]])
  bp <- suppressWarnings(as.numeric(pos[[2L]]))
  same_chrom <- chrom[-1L] == chrom[-length(chrom)]
  bad <- which(same_chrom & diff(bp) < 0)
  if (length(bad)) {
    stop(
      "VCF positions are not sorted within chromosome '", chrom[bad[1L]], "' (POS ",
      format(bp[bad[1L]], scientific = FALSE), " is followed by ",
      format(bp[bad[1L] + 1L], scientific = FALSE), "). ",
      "Sort the VCF by CHROM and POS (e.g. `bcftools sort`) before Beagle imputation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Function to run Beagle for imputation or phasing
run_beagle <- function(vcf_file, output_prefix, beagle_jar,
                       map_file = NULL, ref_file = NULL, generate_map = FALSE,
                       markers_file = NULL, ped_file = NULL,
                       ibd = FALSE, ne = 10000, nthreads = 4,
                       window = 40.0, overlap = 2.0, burnin = 3,
                       iterations = 12, phase_states = 280,
                       impute = TRUE, imp_states = 1600, imp_segment = 6.0,
                       gp = TRUE, gprobs = FALSE, gt = TRUE,
                       java_bin = NULL, java_memory = "4g",
                       seed = -99999L, err = NULL, em = TRUE,
                       chrom = NULL, excludesamples = NULL,
                       excludemarkers = NULL, imp_step = NULL,
                       imp_nsteps = NULL, cluster = NULL, ap = FALSE
                       ) {

  # Validate inputs
  if (!file.exists(vcf_file)) stop("VCF file does not exist.")
  if (!file.exists(beagle_jar)) stop("Beagle jar file does not exist.", call. = FALSE)
  if (!is.null(markers_file) || !is.null(ped_file) || isTRUE(ibd) || isTRUE(gprobs) || !isTRUE(gt)) {
    stop(
      "markers_file, ped_file, ibd, gprobs, and gt=FALSE are not Beagle 5.4 arguments. ",
      "Use excludemarkers/excludesamples or a supported Beagle 5.4 option instead.",
      call. = FALSE
    )
  }
  java_bin <- java_bin %||% check_and_install_java()
  vcf_file <- normalizePath(vcf_file, winslash = "/", mustWork = TRUE)
  beagle_jar <- validate_beagle_jar(beagle_jar)
  output_dir <- normalizePath(dirname(output_prefix), winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  output_prefix <- output_dir |>
    file.path(basename(output_prefix))

  sanitized_vcf <- sanitize_vcf_for_external_tools(vcf_file, mode = "beagle")
  if (!identical(sanitized_vcf$path, vcf_file)) {
    on.exit(unlink(sanitized_vcf$path), add = TRUE)
    vcf_file <- sanitized_vcf$path
  }
  gp_beagle_check_sorted_positions(vcf_file)
  if (isTRUE(sanitized_vcf$changed)) {
    message(
      "VCF preflight for Beagle kept ",
      sanitized_vcf$metrics$output_variants,
      " of ",
      sanitized_vcf$metrics$input_variants,
      " variants."
    )
  }

  # Check if map_file is provided and validate the format
  if (!is.null(map_file)) {
    if (!file.exists(map_file)) {
      stop("Map file does not exist. Please provide a valid map file.")
    } else {
      validate_map_file(map_file)
    }
  } else if (isTRUE(generate_map)) {
    map_file <- generate_genetic_map(vcf_file, output_dir = dirname(output_prefix))
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

  beagle_args <- c(
    paste0("-Xmx", java_memory),
    "-jar", beagle_jar,
    paste0("gt=", vcf_file),
    paste0("out=", output_prefix),
    paste0("ne=", as.integer(ne)),
    paste0("nthreads=", as.integer(nthreads)),
    paste0("window=", format(as.numeric(window), scientific = FALSE)),
    paste0("overlap=", format(as.numeric(overlap), scientific = FALSE)),
    paste0("burnin=", as.integer(burnin)),
    paste0("iterations=", as.integer(iterations)),
    paste0("phase-states=", as.integer(phase_states)),
    paste0("seed=", as.integer(seed)),
    paste0("em=", tolower(as.character(isTRUE(em))))
  )

  if (!is.null(err)) beagle_args <- c(beagle_args, paste0("err=", format(as.numeric(err), scientific = TRUE)))
  if (!is.null(chrom)) beagle_args <- c(beagle_args, paste0("chrom=", as.character(chrom)))
  if (!is.null(excludesamples)) {
    validate_markers_file(excludesamples)
    beagle_args <- c(beagle_args, paste0("excludesamples=", normalizePath(excludesamples, winslash = "/", mustWork = TRUE)))
  }
  if (!is.null(excludemarkers)) {
    validate_markers_file(excludemarkers)
    beagle_args <- c(beagle_args, paste0("excludemarkers=", normalizePath(excludemarkers, winslash = "/", mustWork = TRUE)))
  }

  # Add optional files
  if (!is.null(map_file)) beagle_args <- c(beagle_args, paste0("map=", normalizePath(map_file, winslash = "/", mustWork = TRUE)))
  if (!is.null(ref_file)) beagle_args <- c(beagle_args, paste0("ref=", normalizePath(ref_file, winslash = "/", mustWork = TRUE)))
  if (!is.null(markers_file)) beagle_args <- c(beagle_args, paste0("markers=", normalizePath(markers_file, winslash = "/", mustWork = TRUE)))
  if (!is.null(ped_file)) beagle_args <- c(beagle_args, paste0("ped=", normalizePath(ped_file, winslash = "/", mustWork = TRUE)))

  # Add imputation-specific flags
  if (impute) {
    beagle_args <- c(
      beagle_args,
      "impute=true",
      paste0("imp-states=", as.integer(imp_states)),
      paste0("imp-segment=", format(as.numeric(imp_segment), scientific = FALSE))
    )
    if (!is.null(imp_step)) beagle_args <- c(beagle_args, paste0("imp-step=", format(as.numeric(imp_step), scientific = FALSE)))
    if (!is.null(imp_nsteps)) beagle_args <- c(beagle_args, paste0("imp-nsteps=", as.integer(imp_nsteps)))
    if (!is.null(cluster)) beagle_args <- c(beagle_args, paste0("cluster=", format(as.numeric(cluster), scientific = FALSE)))
  } else {
    beagle_args <- c(beagle_args, "impute=false")
  }

  # Add IBD detection, GP, and GT options
  #if (ibd) beagle_args <- c(beagle_args, "ibd=true")
  if (gp) beagle_args <- c(beagle_args, "gp=true")
  if (ap) beagle_args <- c(beagle_args, "ap=true")



  # Log Beagle execution and capture any potential errors
  log_file <- paste0(output_prefix, "_log.txt")
  message("Running Beagle: ", paste(c(shQuote(java_bin), shQuote(beagle_args)), collapse = " "))
  result <- beagle_system2(
    java_bin,
    args = gp_quote_system_args(beagle_args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (is.null(status)) {
    status <- 0L
  }

  # Write Beagle output log to a file for troubleshooting
  writeLines(result, log_file)

  if (!identical(as.integer(status), 0L)) {
    log_tail <- utils::tail(result[nzchar(trimws(result))], 6L)
    stop(
      "Beagle failed. Last log lines:\n  ", paste(log_tail, collapse = "\n  "),
      "\nFull log: ", log_file, call. = FALSE
    )
  }

  # Check if the VCF output file was created
  phased_vcf <- paste0(output_prefix, ".vcf.gz")
  if (file.exists(phased_vcf)) {
    message("Phased VCF file created: ", phased_vcf)
  } else {
    stop(
      "Beagle returned status 0 but did not create the expected VCF: ",
      phased_vcf, ". Check ", log_file, ".",
      call. = FALSE
    )
  }

  invisible(list(
    output_vcf = normalizePath(phased_vcf, winslash = "/", mustWork = TRUE),
    log_file = normalizePath(log_file, winslash = "/", mustWork = TRUE),
    beagle_log_file = if (file.exists(paste0(output_prefix, ".log"))) {
      normalizePath(paste0(output_prefix, ".log"), winslash = "/", mustWork = TRUE)
    } else NULL,
    status = as.integer(status),
    command_args = beagle_args,
    input_preflight = sanitized_vcf$metrics,
    beagle_jar = beagle_jar,
    java = normalizePath(java_bin, winslash = "/", mustWork = TRUE)
  ))
}

# Function to check Java, download Beagle, and run the appropriate Beagle process
impute_vcf_with_beagle <- function(vcf_file, output_prefix = "imputed",
                                   map_file = NULL, ref_file = NULL,
                                   markers_file = NULL, ped_file = NULL,
                                   generate_map = FALSE,
                                   ibd = FALSE, ne = 10000, nthreads = 4,
                                   window = 40.0, overlap = 2.0, burnin = 3,
                                   iterations = 12, phase_states = 280,
                                   impute = TRUE, imp_states = 1600, imp_segment = 6.0,
                                   gp = TRUE, gprobs = FALSE, gt = TRUE,
                                   java_installed = FALSE,
                                   beagle_version = "29Oct24.c8e",
                                   beagle_jar = NULL,
                                   beagle_dir = NULL,
                                   java_bin = NULL,
                                   java_memory = "4g",
                                   seed = -99999L, err = NULL, em = TRUE,
                                   chrom = NULL, excludesamples = NULL,
                                   excludemarkers = NULL, imp_step = NULL,
                                   imp_nsteps = NULL, cluster = NULL,
                                   ap = FALSE) {

  # Step 1: Check Java if needed
  if (!isTRUE(java_installed) || is.null(java_bin)) {
    java_bin <- check_and_install_java()
  }

  # Step 2: Download Beagle
  if (is.null(beagle_jar)) {
    beagle_jar <- download_beagle(
      output_dir = beagle_dir %||% predictpror_beagle_cache_dir(),
      beagle_version = beagle_version
    )
  }

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
             imp_segment = imp_segment, gp = gp, gprobs = gprobs, gt = gt,
             java_bin = java_bin, java_memory = java_memory,
             seed = seed, err = err, em = em, chrom = chrom,
             excludesamples = excludesamples, excludemarkers = excludemarkers,
             imp_step = imp_step, imp_nsteps = imp_nsteps, cluster = cluster,
             ap = ap)
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
