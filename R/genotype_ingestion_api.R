read_hapmap_file <- function(filepath) {
  if (length(filepath) != 1L || !is.character(filepath) || !nzchar(filepath)) {
    stop("filepath must be one non-empty file path.", call. = FALSE)
  }
  filepath <- normalizePath(filepath, winslash = "/", mustWork = TRUE)
  lower <- tolower(filepath)
  read_one <- function(path) {
    data.table::fread(
      path,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      data.table = TRUE,
      na.strings = "NA",
      showProgress = FALSE
    )
  }

  read_gzip <- function(path) {
    expanded <- tempfile("predictpror-hapmap-gunzip-", fileext = ".txt")
    on.exit(unlink(expanded, force = TRUE), add = TRUE)
    input_con <- gzfile(path, open = "rb")
    output_con <- file(expanded, open = "wb")
    input_open <- TRUE
    output_open <- TRUE
    on.exit(if (input_open) close(input_con), add = TRUE)
    on.exit(if (output_open) close(output_con), add = TRUE)
    repeat {
      chunk <- readBin(input_con, what = "raw", n = 1024L * 1024L)
      if (!length(chunk)) break
      writeBin(chunk, output_con)
    }
    close(input_con)
    input_open <- FALSE
    close(output_con)
    output_open <- FALSE
    read_one(expanded)
  }

  if (grepl("\\.zip$", lower)) {
    listing <- utils::unzip(filepath, list = TRUE)
    files <- listing$Name[!grepl("[/\\\\]$", listing$Name)]
    supported <- files[grepl("(?i)\\.(hmp(\\.txt)?|hapmap|txt|tsv)(\\.gz)?$", files, perl = TRUE)]
    if (length(supported) != 1L) {
      stop("A HapMap zip archive must contain exactly one HapMap text file.", call. = FALSE)
    }
    temp_dir <- tempfile("predictpror-hapmap-")
    dir.create(temp_dir, recursive = TRUE)
    on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    utils::unzip(filepath, files = supported, exdir = temp_dir)
    extracted <- normalizePath(file.path(temp_dir, supported), winslash = "/", mustWork = TRUE)
    root <- paste0(normalizePath(temp_dir, winslash = "/", mustWork = TRUE), "/")
    if (!startsWith(extracted, root)) {
      stop("Unsafe path found inside the HapMap zip archive.", call. = FALSE)
    }
    hapmap <- if (grepl("(?i)\\.gz$", extracted, perl = TRUE)) read_gzip(extracted) else read_one(extracted)
  } else if (grepl("(?i)\\.(hmp(\\.txt)?|hapmap|txt|tsv)(\\.gz)?$", lower, perl = TRUE)) {
    hapmap <- if (grepl("(?i)\\.gz$", lower, perl = TRUE)) read_gzip(filepath) else read_one(filepath)
  } else {
    stop(
      "Unsupported HapMap filename. Use .hmp, .hmp.txt, .hapmap, .txt, or .tsv, optionally gzip-compressed or in a .zip archive.",
      call. = FALSE
    )
  }
  gp_normalize_hapmap_columns(hapmap)
}

#' Quality control and alternate-allele recoding for HapMap genotypes
#'
#' Reads a HapMap file (or accepts a HapMap table), validates its 11 metadata
#' columns, decodes diploid or explicit polyploid allele-copy calls, applies
#' marker and sample quality control, and returns a sample-by-marker matrix.
#'
#' The coding is deterministic across files: values from zero through `ploidy`
#' count the second allele in the HapMap `alleles` field. Centered dosage
#' subtracts `ploidy / 2`. Polyploid heterozygotes must contain all allele
#' copies (for example `AAAG`, `AAGG`, or `AGGG` at ploidy four); an IUPAC
#' heterozygote code cannot identify polyploid dosage and is rejected.
#'
#' @param hapmap_file_name Optional HapMap filename.
#' @param hapmap_file_path Directory containing `hapmap_file_name`.
#' @param hapmap Optional in-memory HapMap data frame or data table. Supply
#'   either this argument or a filename, not both.
#' @param maf_threshold Minimum minor-allele frequency in `[0, 1]`, or `NULL`.
#' @param het_threshold Maximum marker heterozygosity in `[0, 1]`, or `NULL`.
#' @param ind_call_rate_threshold Minimum sample call rate in `[0, 1]`, or
#'   `NULL`.
#' @param snp_call_rate_threshold Minimum marker call rate in `[0, 1]`, or
#'   `NULL`.
#' @param impute Logical; impute remaining missing alternate-allele dosages.
#' @param imputation_method Native dosage imputation method: `"mean"`,
#'   `"median"`, `"mode"`, or `"knn"`.
#' @param impute_knn_k Number of neighbours for native KNN imputation.
#' @param ploidy One positive integer or `"auto"`. Automatic HapMap inference
#'   requires at least one complete allele-copy call; otherwise diploid is the
#'   backward-compatible default.
#' @param recode_format `"alt_dosage"`, `"centered_dosage"`,
#'   `"allele_frequency"`, or a legacy diploid alias.
#' @param out_put_map Logical; include the retained HapMap metadata table.
#' @param message Logical; print a concise processing summary.
#' @param invalid_call Whether invalid or allele-inconsistent genotype calls
#'   cause an error or are set to missing.
#' @param non_biallelic Whether non-biallelic marker rows are dropped or cause
#'   an error.
#' @param ... Reserved for compatibility.
#'
#' @return A list with `snps_matrix`, optional `snp_map`,
#'   `qc_metrics_and_summary_stat`, `coding`, and `format_report`.
#' @export
hmp_qc_recode <- function(hapmap_file_name = NULL,
                          hapmap_file_path = NULL,
                          hapmap = NULL,
                          maf_threshold = 0.01,
                          het_threshold = 0.1,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          impute = TRUE,
                          recode_format = "0,1,2",
                          ploidy = "auto",
                          imputation_method = c("mean", "median", "mode", "knn"),
                          impute_knn_k = 5L,
                          out_put_map = TRUE,
                          message = TRUE,
                          invalid_call = c("error", "missing"),
                          non_biallelic = c("drop", "error"),
                          ...) {
  invalid_call <- match.arg(invalid_call)
  non_biallelic <- match.arg(non_biallelic)
  has_file <- !is.null(hapmap_file_name)
  if (has_file && !is.null(hapmap)) {
    stop("Supply either a HapMap file or an in-memory hapmap table, not both.", call. = FALSE)
  }
  if (has_file) {
    hapmap_file_path <- hapmap_file_path %||% getwd()
    hapmap <- read_hapmap_file(file.path(hapmap_file_path, hapmap_file_name))
  } else if (is.null(hapmap)) {
    stop("HapMap data are missing. Supply hapmap or hapmap_file_name.", call. = FALSE)
  }

  parsed <- gp_hapmap_to_alt_dosage(
    hapmap,
    invalid_call = invalid_call,
    non_biallelic = non_biallelic,
    ploidy = ploidy
  )
  result <- gp_qc_alt_dosage(
    dosage = parsed$dosage,
    map = parsed$map,
    maf_threshold = maf_threshold,
    het_threshold = het_threshold,
    ind_call_rate_threshold = ind_call_rate_threshold,
    snp_call_rate_threshold = snp_call_rate_threshold,
    impute = impute,
    recode_format = recode_format,
    ploidy = parsed$report$ploidy,
    imputation_method = imputation_method,
    impute_knn_k = impute_knn_k
  )
  result$format_report <- parsed$report
  if (!isTRUE(out_put_map)) {
    result$snp_map <- NULL
  }
  if (isTRUE(message)) {
    base::message(
      "HapMap processing retained ", ncol(result$snps_matrix), " markers and ",
      nrow(result$snps_matrix), " samples; coding: ", result$coding, "."
    )
  }
  result
}

#' Quality control and alternate-allele recoding for VCF genotypes
#'
#' Supports plain or gzip-compressed VCF with uniform-ploidy biallelic GT calls. GT
#' may appear alone or within fields such as `GT:DS`, `GT:GP`, or `AD:GT:DP`.
#' Phased and unphased calls are accepted. Multiallelic, symbolic, mixed-ploidy,
#' and allele-index-greater-than-one rows are rejected or excluded explicitly.
#'
#' Values from zero through `ploidy` count the VCF ALT allele.
#'
#' @param vcf_file_name VCF or VCF.GZ filename.
#' @param vcf_file_path Directory containing the VCF file.
#' @inheritParams hmp_qc_recode
#' @param batch_size,num_cores Compatibility parameters used by the legacy
#'   PLINK engine.
#' @param ld_pruning,ld_pruning_method,window_size,step_size,r2_threshold,use_kb_window,phased,use_founders
#'   Compatibility parameters used only when `qc_engine = "plink"`.
#' @param qc_engine `"r"` for the deterministic built-in reader (default), or
#'   `"plink"` for the legacy external-PLINK workflow.
#'
#' @return A list with `snps_matrix`, optional `snp_map`,
#'   `qc_metrics_and_summary_stat`, `coding`, and `format_report`.
#' @export
vcf_qc_recode <- function(vcf_file_name = NULL,
                          vcf_file_path = NULL,
                          maf_threshold = 0.05,
                          het_threshold = 0.2,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          impute = FALSE,
                          recode_format = "0,1,2",
                          ploidy = "auto",
                          imputation_method = c("mean", "median", "mode", "knn"),
                          impute_knn_k = 5L,
                          out_put_map = TRUE,
                          batch_size = 2000,
                          num_cores = NULL,
                          ld_pruning = FALSE,
                          ld_pruning_method = "indep-pairwise",
                          window_size = 50,
                          step_size = 5,
                          r2_threshold = 0.2,
                          use_kb_window = TRUE,
                          phased = TRUE,
                          use_founders = FALSE,
                          message = TRUE,
                          qc_engine = c("r", "plink")) {
  qc_engine <- match.arg(qc_engine)
  if (identical(qc_engine, "plink")) {
    resolved_ploidy <- gp_validate_ploidy(ploidy)
    if (!identical(resolved_ploidy, "auto") && resolved_ploidy != 2L) {
      stop("The legacy PLINK QC engine is restricted to diploid data; use qc_engine = 'r' for polyploid dosage.", call. = FALSE)
    }
    return(vcf_qc_recode_legacy(
      vcf_file_name = vcf_file_name,
      vcf_file_path = vcf_file_path,
      maf_threshold = maf_threshold,
      het_threshold = het_threshold,
      ind_call_rate_threshold = ind_call_rate_threshold,
      snp_call_rate_threshold = snp_call_rate_threshold,
      impute = impute,
      recode_format = recode_format,
      out_put_map = out_put_map,
      batch_size = batch_size,
      num_cores = num_cores,
      ld_pruning = ld_pruning,
      ld_pruning_method = ld_pruning_method,
      window_size = window_size,
      step_size = step_size,
      r2_threshold = r2_threshold,
      use_kb_window = use_kb_window,
      phased = phased,
      use_founders = use_founders,
      message = message
    ))
  }
  if (isTRUE(ld_pruning)) {
    stop(
      "ld_pruning = TRUE requires qc_engine = 'plink'; the built-in R engine does not silently substitute an LD algorithm.",
      call. = FALSE
    )
  }
  if (is.null(vcf_file_name)) {
    stop("vcf_file_name is required.", call. = FALSE)
  }
  vcf_file_path <- vcf_file_path %||% getwd()
  input_vcf <- normalizePath(
    file.path(vcf_file_path, vcf_file_name), winslash = "/", mustWork = TRUE
  )
  sanitized <- sanitize_vcf_for_external_tools(input_vcf, mode = "recode", ploidy = ploidy)
  if (!identical(sanitized$path, input_vcf)) {
    on.exit(unlink(sanitized$path, force = TRUE), add = TRUE)
  }
  vcf_data <- vcf_read_table(sanitized$path)
  id_col <- vcf_find_column(vcf_data, "ID")
  chrom_col <- vcf_find_column(vcf_data, c("#CHROM", "CHROM"))
  pos_col <- vcf_find_column(vcf_data, "POS")
  ref_col <- vcf_find_column(vcf_data, "REF")
  alt_col <- vcf_find_column(vcf_data, "ALT")
  marker_ids <- as.character(vcf_data[[id_col]])
  replace_id <- is.na(marker_ids) | !nzchar(marker_ids) | marker_ids == "." |
    duplicated(marker_ids) | duplicated(marker_ids, fromLast = TRUE)
  if (any(replace_id)) {
    coordinate_ids <- paste0(
      vcf_data[[chrom_col]], ":", vcf_data[[pos_col]], ":",
      vcf_data[[ref_col]], ">", vcf_data[[alt_col]]
    )
    marker_ids[replace_id] <- coordinate_ids[replace_id]
    vcf_data[[id_col]] <- marker_ids
  }
  if (anyDuplicated(marker_ids)) {
    stop("VCF variants do not have unique IDs or unique CHROM:POS:REF:ALT coordinates.", call. = FALSE)
  }
  dosage <- gp_parse_vcf_alt_dosage(vcf_data, ploidy = sanitized$metrics$ploidy)
  result <- gp_qc_alt_dosage(
    dosage = dosage,
    map = vcf_data[, 1:9],
    marker_id_col = 3L,
    maf_threshold = maf_threshold,
    het_threshold = het_threshold,
    ind_call_rate_threshold = ind_call_rate_threshold,
    snp_call_rate_threshold = snp_call_rate_threshold,
    impute = impute,
    recode_format = recode_format,
    ploidy = sanitized$metrics$ploidy,
    imputation_method = imputation_method,
    impute_knn_k = impute_knn_k
  )
  result$format_report <- sanitized$metrics
  result$format_report$generated_marker_ids <- as.integer(sum(replace_id))
  if (!isTRUE(out_put_map)) {
    result$snp_map <- NULL
  }
  if (isTRUE(message)) {
    base::message(
      "VCF processing retained ", ncol(result$snps_matrix), " markers and ",
      nrow(result$snps_matrix), " samples; coding: ", result$coding, "."
    )
  }
  result
}

# Call-rate QC must use the OBSERVED calls: after Beagle every call rate is
# 1, so filtering later removes nothing. Markers are dropped first, then
# samples are judged on the kept markers. Thresholds follow vcf_qc_recode():
# a value >= 0.5 is the minimum call rate, a value < 0.5 the maximum missing
# rate. Beagle identifies markers by CHROM:POS, which survives ID sanitising.
gp_beagle_call_rate_exclusions <- function(vcf_path,
                                           snp_call_rate_threshold = NULL,
                                           ind_call_rate_threshold = NULL,
                                           output_prefix) {
  min_call <- function(x) {
    if (is.null(x) || !length(x) || is.na(x[1L])) return(NULL)
    x <- as.numeric(x[1L])
    if (x < 0.5) 1 - x else x
  }
  snp_min <- min_call(snp_call_rate_threshold)
  ind_min <- min_call(ind_call_rate_threshold)
  out <- list(excludemarkers = NULL, excludesamples = NULL,
              markers_removed = 0L, samples_removed = 0L)
  if (is.null(snp_min) && is.null(ind_min)) {
    return(out)
  }
  vcf <- vcf_read_table(vcf_path)
  cols <- vcf_required_columns(vcf)
  samples <- names(vcf)[10:ncol(vcf)]
  formats <- as.character(vcf[[cols$format]])
  formats[is.na(formats) | !nzchar(formats)] <- "GT"
  missing <- matrix(FALSE, nrow(vcf), length(samples))
  for (fmt in unique(formats)) {
    rows <- which(formats == fmt)
    gt_index <- match("GT", strsplit(fmt, ":", fixed = TRUE)[[1L]])
    for (j in seq_along(samples)) {
      gt <- vcf_extract_gt_field(vcf[[samples[j]]][rows], gt_index)
      missing[rows, j] <- is.na(gt) | !nzchar(gt) | grepl(".", gt, fixed = TRUE)
    }
  }
  marker_id <- paste(vcf[[1L]], vcf[[2L]], sep = ":")
  keep_marker <- rep(TRUE, nrow(vcf))
  if (!is.null(snp_min)) {
    keep_marker <- (1 - rowMeans(missing)) >= snp_min
  }
  keep_sample <- rep(TRUE, length(samples))
  if (!is.null(ind_min) && any(keep_marker)) {
    keep_sample <- (1 - colMeans(missing[keep_marker, , drop = FALSE])) >= ind_min
  }
  if (!any(keep_marker) || !any(keep_sample)) {
    stop(
      "Call-rate QC before Beagle would remove every ",
      if (!any(keep_marker)) "marker" else "sample",
      " (snp_call_rate_threshold = ", format(snp_call_rate_threshold %||% "NULL"),
      ", ind_call_rate_threshold = ", format(ind_call_rate_threshold %||% "NULL"), ").",
      call. = FALSE
    )
  }
  if (any(!keep_marker)) {
    out$excludemarkers <- paste0(output_prefix, "_callrate_excluded_markers.txt")
    writeLines(marker_id[!keep_marker], out$excludemarkers)
    out$markers_removed <- sum(!keep_marker)
  }
  if (any(!keep_sample)) {
    out$excludesamples <- paste0(output_prefix, "_callrate_excluded_samples.txt")
    writeLines(samples[!keep_sample], out$excludesamples)
    out$samples_removed <- sum(!keep_sample)
  }
  out
}


#' Impute a raw VCF or HapMap file with Beagle 5.4
#'
#' This is the public end-to-end entry point for Beagle preprocessing. HapMap
#' input is first converted to VCF 4.3 with stable declared-allele orientation;
#' VCF input is validated and sanitized for diploid biallelic GT processing.
#' The result records the converted input, Beagle logs, command arguments,
#' format-cleanup counts, and (optionally) a recoded genotype matrix.
#'
#' Beagle imputes sporadically missing genotypes during phasing. Imputation of
#' markers absent from the target requires a phased reference panel supplied
#' through `ref_file`.
#'
#' @param input A raw VCF/VCF.GZ filename, HapMap filename, CSV/TXT genotype
#'   table filename (see [convert_csv_to_vcf()]), or in-memory HapMap table.
#' @param input_format One of `"auto"`, `"vcf"`, `"hapmap"`, or `"csv"`.
#'   `"auto"` uses the file name (`.vcf`/`.vcf.gz`) and, for other text files,
#'   the header: VCF marker columns (`CHROM`, `POS`, `REF`, `ALT`, ...) mean a
#'   CSV/TXT genotype table, otherwise HapMap. CSV/TXT tables and HapMap are
#'   converted to VCF before Beagle runs.
#' @param csv_input_coding Numeric coding of a CSV/TXT genotype table:
#'   `"alt_dosage"` (0/1/2) or `"centered_dosage"` (-1/0/1).
#' @param output_prefix Output prefix for Beagle. Beagle creates `.vcf.gz` and
#'   `.log` files from this prefix (plus `_log.txt`, and `_input.vcf` for
#'   HapMap input). `NULL` (default) uses `<input name>_beagle` in the current
#'   working directory; pass e.g. `file.path(tempdir(), "beagle")` to keep the
#'   files out of it.
#' @param beagle_jar Optional path to a Beagle 5.4 JAR. If omitted, the package
#'   uses `PREDICTPRO_BEAGLE_JAR` or downloads the pinned official version to
#'   the user cache.
#' @param return_genotypes Logical; also run `vcf_qc_recode()` on the Beagle
#'   output.
#' @param recode_format Numeric coding requested when `return_genotypes` is
#'   true.
#' @param qc_args Named list of additional `vcf_qc_recode()` arguments.
#' @param invalid_call,non_biallelic HapMap conversion policies; see
#'   [convert_hapmap_to_vcf()].
#' @param ploidy Must resolve to two. Beagle 5.4 is a diploid backend; use
#'   `hmp_qc_recode()` or `vcf_qc_recode()` with native imputation for
#'   polyploid dosage.
#' @param snp_call_rate_threshold,ind_call_rate_threshold Optional call-rate
#'   filters applied to the OBSERVED calls before imputation (after Beagle every
#'   call rate is 1). Same convention as [vcf_qc_recode()]: a value >= 0.5 is the
#'   minimum call rate, a value < 0.5 the maximum missing rate. Failing markers
#'   (by `CHROM:POS`) and samples are passed to Beagle as
#'   `excludemarkers`/`excludesamples`; markers are filtered first.
#' @param ... Beagle 5.4 options passed to [impute_vcf_with_beagle()].
#'
#' @return An object of class `predictpror_beagle_result` containing input and
#'   output paths, conversion and Beagle provenance, and optional genotypes.
#' @export
impute_genotypes_with_beagle <- function(input,
                                         input_format = c("auto", "vcf", "hapmap", "csv"),
                                         output_prefix = NULL,
                                         beagle_jar = NULL,
                                         return_genotypes = TRUE,
                                         recode_format = "0,1,2",
                                         qc_args = list(),
                                         invalid_call = c("error", "missing"),
                                         non_biallelic = c("drop", "error"),
                                         ploidy = 2L,
                                         snp_call_rate_threshold = NULL,
                                         ind_call_rate_threshold = NULL,
                                         csv_input_coding = c("alt_dosage", "centered_dosage"),
                                         ...) {
  input_format <- match.arg(input_format)
  csv_input_coding <- match.arg(csv_input_coding)
  invalid_call <- match.arg(invalid_call)
  non_biallelic <- match.arg(non_biallelic)
  ploidy <- gp_validate_ploidy(ploidy)
  is_table <- is.data.frame(input)
  if (identical(input_format, "auto")) {
    if (is_table) {
      input_format <- "hapmap"
    } else if (is.character(input) && length(input) == 1L &&
               grepl("(?i)\\.vcf(\\.gz)?$", input, perl = TRUE)) {
      input_format <- "vcf"
    } else if (is.character(input) && length(input) == 1L && gp_is_genotype_table_file(input)) {
      input_format <- "csv"
    } else {
      input_format <- "hapmap"
    }
  }
  if (is_table && !identical(input_format, "hapmap")) {
    stop("In-memory input is supported only for input_format = 'hapmap'.", call. = FALSE)
  }
  if (!is_table) {
    if (!is.character(input) || length(input) != 1L || !nzchar(input)) {
      stop("input must be one genotype filename or an in-memory HapMap table.", call. = FALSE)
    }
    input <- normalizePath(input, winslash = "/", mustWork = TRUE)
  }

  if (identical(input_format, "vcf") && identical(ploidy, "auto")) {
    raw_vcf <- vcf_read_table(input)
    raw_cols <- vcf_required_columns(raw_vcf)
    raw_samples <- names(raw_vcf)[10:ncol(raw_vcf)]
    raw_formats <- as.character(raw_vcf[[raw_cols$format]])
    raw_formats[is.na(raw_formats) | !nzchar(raw_formats)] <- "GT"
    raw_gt <- list()
    for (fmt in unique(raw_formats)) {
      rows <- which(raw_formats == fmt)
      gt_index <- match("GT", strsplit(fmt, ":", fixed = TRUE)[[1L]])
      if (is.na(gt_index)) stop("Beagle VCF input contains a row without GT.", call. = FALSE)
      for (sample in raw_samples) {
        raw_gt[[length(raw_gt) + 1L]] <- unique(vcf_extract_gt_field(raw_vcf[[sample]][rows], gt_index))
      }
    }
    # Ploidy only depends on the distinct calls; collect them without
    # growing one long vector cell by cell.
    ploidy <- gp_resolve_gt_ploidy(unique(unlist(raw_gt, use.names = FALSE)), "auto", context = "Beagle VCF")
  }
  if (identical(ploidy, "auto") && identical(input_format, "hapmap")) {
    ploidy <- gp_beagle_hapmap_call_ploidy(input) %||% "auto"
  }
  if (identical(ploidy, "auto") && identical(input_format, "csv")) {
    gp_beagle_check_csv_dosage(input, csv_input_coding)
  }
  if (!identical(ploidy, "auto") && ploidy != 2L) {
    gp_beagle_polyploid_stop(ploidy)
  }

  if (is.null(output_prefix)) {
    stem <- if (is_table) "hapmap" else tools::file_path_sans_ext(
      tools::file_path_sans_ext(basename(input))
    )
    output_prefix <- file.path(getwd(), paste0(stem, "_beagle"))
  }
  output_dir <- normalizePath(dirname(output_prefix), winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  output_prefix <- file.path(output_dir, basename(output_prefix))

  conversion_report <- NULL
  beagle_input <- input
  if (identical(input_format, "csv")) {
    beagle_input <- paste0(output_prefix, "_input.vcf")
    convert_csv_to_vcf(input, beagle_input, ploidy = if (identical(ploidy, "auto")) 2L else ploidy,
                       input_coding = csv_input_coding)
    beagle_input <- normalizePath(beagle_input, winslash = "/", mustWork = TRUE)
  } else if (identical(input_format, "hapmap")) {
    beagle_input <- paste0(output_prefix, "_input.vcf")
    converted <- convert_hapmap_to_vcf(
      input_hapmap = input,
      output_vcf = beagle_input,
      invalid_call = invalid_call,
      non_biallelic = non_biallelic,
      ploidy = ploidy
    )
    conversion_report <- attr(converted, "conversion_report", exact = TRUE)
    beagle_input <- unclass(converted)
  }

  beagle_args <- list(...)
  call_rate_qc <- gp_beagle_call_rate_exclusions(
    vcf_path = beagle_input,
    snp_call_rate_threshold = snp_call_rate_threshold,
    ind_call_rate_threshold = ind_call_rate_threshold,
    output_prefix = output_prefix
  )
  for (key in c("excludemarkers", "excludesamples")) {
    if (is.null(call_rate_qc[[key]])) next
    if (!is.null(beagle_args[[key]])) {
      # keep the user's own exclusion list and add the call-rate failures
      writeLines(unique(c(readLines(beagle_args[[key]], warn = FALSE), readLines(call_rate_qc[[key]]))),
                 call_rate_qc[[key]])
    }
    beagle_args[[key]] <- call_rate_qc[[key]]
  }
  if (call_rate_qc$markers_removed || call_rate_qc$samples_removed) {
    message(
      "Call-rate QC on observed calls before Beagle: removed ",
      call_rate_qc$markers_removed, " markers and ", call_rate_qc$samples_removed, " samples."
    )
  }

  beagle_result <- do.call(impute_vcf_with_beagle, c(list(
    vcf_file = beagle_input,
    output_prefix = output_prefix,
    beagle_jar = beagle_jar
  ), beagle_args))
  genotypes <- NULL
  if (isTRUE(return_genotypes)) {
    defaults <- list(
      vcf_file_name = basename(beagle_result$output_vcf),
      vcf_file_path = dirname(beagle_result$output_vcf),
      maf_threshold = NULL,
      het_threshold = NULL,
      ind_call_rate_threshold = 0,
      snp_call_rate_threshold = 0,
      impute = FALSE,
      recode_format = recode_format,
      ploidy = 2L,
      out_put_map = TRUE,
      message = FALSE,
      qc_engine = "r"
    )
    if (!is.list(qc_args) || (length(qc_args) && is.null(names(qc_args)))) {
      stop("qc_args must be a named list.", call. = FALSE)
    }
    genotypes <- do.call(
      vcf_qc_recode,
      utils::modifyList(defaults, qc_args, keep.null = TRUE)
    )
  }

  result <- list(
    input_format = input_format,
    input = if (is_table) "<in-memory HapMap>" else input,
    beagle_input_vcf = normalizePath(beagle_input, winslash = "/", mustWork = TRUE),
    conversion_report = conversion_report,
    beagle = beagle_result,
    output_vcf = beagle_result$output_vcf,
    call_rate_prefilter = call_rate_qc[c("markers_removed", "samples_removed")],
    genotypes = genotypes
  )
  class(result) <- c("predictpror_beagle_result", "list")
  result
}


# ---- Beagle is diploid-only -------------------------------------------------
gp_beagle_polyploid_stop <- function(ploidy = NULL) {
  what <- if (is.null(ploidy)) "polyploid (dosage above 2)" else paste0("ploidy ", ploidy)
  range <- if (is.null(ploidy)) "0..ploidy" else paste0("0..", ploidy)
  stop(
    "Beagle imputation works on diploid genotypes only (it phases two haplotypes per line); ",
    "these genotypes are ", what, ". Use imputation_method = \"knn\" (ploidy-aware KNN ",
    "imputation of dosages ", range, ") or \"mean\"",
    if (is.null(ploidy)) ", and give the ploidy, e.g. ploidy = 4L" else "",
    ".",
    call. = FALSE
  )
}

# Ploidy implied by polyploid HapMap calls (AAGG = 4 copies); NULL for diploid
# or unreadable input. Only the first rows are inspected.
gp_beagle_hapmap_call_ploidy <- function(input, n_rows = 200L) {
  tab <- if (is.data.frame(input)) utils::head(input, n_rows) else tryCatch(
    as.data.frame(data.table::fread(input, nrows = n_rows, colClasses = "character",
                                    check.names = FALSE, header = TRUE)),
    error = function(e) NULL
  )
  if (is.null(tab) || ncol(tab) < 12L) return(NULL)
  calls <- toupper(trimws(unlist(tab[, 12:ncol(tab)], use.names = FALSE)))
  calls <- calls[!is.na(calls) & nzchar(calls) & !grepl("^N+$|^-+$", calls)]
  calls <- gsub("[/|]", "", calls)
  if (!length(calls)) return(NULL)
  width <- max(nchar(calls))
  if (width > 2L) as.integer(width) else NULL
}

# A CSV/TXT dosage table with values above 2 (or beyond +-1 when centred) is
# polyploid; stop before it is converted as diploid.
gp_beagle_check_csv_dosage <- function(input, input_coding = "alt_dosage", n_rows = 200L) {
  tab <- tryCatch(
    as.data.frame(data.table::fread(input, nrows = n_rows, colClasses = "character", check.names = FALSE,
                                    header = TRUE, na.strings = c("", "NA", "N", "."),
                                    skip = gp_genotype_table_header_line(input) - 1L)),
    error = function(e) NULL
  )
  if (is.null(tab)) return(invisible(NULL))
  names(tab) <- sub("^#", "", names(tab))
  samples <- setdiff(names(tab), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
  if (!length(samples)) return(invisible(NULL))
  values <- trimws(unlist(tab[, samples, drop = FALSE], use.names = FALSE))
  gt <- values[!is.na(values) & grepl("^[01.]([/|][01.])+$", values)]
  if (length(gt) && max(lengths(strsplit(gt, "[/|]"))) > 2L) {
    gp_beagle_polyploid_stop(max(lengths(strsplit(gt, "[/|]"))))
  }
  d <- suppressWarnings(as.numeric(values))
  d <- d[is.finite(d)]
  limit <- if (identical(input_coding, "centered_dosage")) 1 else 2
  if (length(d) && max(abs(d)) > limit) gp_beagle_polyploid_stop(NULL)
  invisible(NULL)
}