#' Convert a CSV or TXT genotype table to VCF
#'
#' Writes a VCF 4.3 file (GT field) from a delimited text table, so numeric
#' genotype data can be imputed with Beagle ([impute_genotypes_with_beagle()])
#' or read with [vcf_qc_recode()]. The separator (comma, tab, ...) is detected
#' automatically. [model_execute()] calls this itself for `csv_file_name`, and
#' [impute_genotypes_with_beagle()] for CSV/TXT input.
#'
#' The table needs the nine VCF marker columns `CHROM`, `POS`, `ID`, `REF`,
#' `ALT`, `QUAL`, `FILTER`, `INFO`, `FORMAT` (a leading `#` on `CHROM` is
#' accepted), followed by one column per sample. Genotypes may be numeric
#' dosages of the ALT allele (`0`, `1`, `2` for diploids; `0..ploidy` in
#' general), centred dosages (`input_coding = "centered_dosage"`, e.g. -1/0/1),
#' or VCF GT calls such as `0/1`. `NA`, `"N"` and `"."` are written as missing.
#'
#' @param input_csv Path to the CSV/TXT genotype table.
#' @param output_vcf Path of the VCF to write. `NULL` returns the VCF lines.
#' @param ploidy Ploidy of the genotype calls (default 2).
#' @param input_coding `"alt_dosage"` (0..ploidy ALT copies) or
#'   `"centered_dosage"` (dosage minus ploidy / 2).
#' @return The VCF lines, or (invisibly) `output_vcf` when it is given.
#' @seealso [convert_hapmap_to_vcf()], [impute_genotypes_with_beagle()]
#' @examples
#' tab <- data.frame(CHROM = 1, POS = c(100, 200), ID = c("m1", "m2"), REF = "A", ALT = "G",
#'                   QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
#'                   L1 = c(0, 2), L2 = c(1, NA))
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(tab, f, row.names = FALSE, na = ".")
#' convert_csv_to_vcf(f)
#' @export
convert_csv_to_vcf <- function(input_csv, output_vcf = NULL, ploidy = 2L,
                               input_coding = c("alt_dosage", "centered_dosage")) {

  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  input_coding <- match.arg(input_coding)
  required_headers <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

  snp_data <- as.data.frame(data.table::fread(
    input_csv, header = TRUE, colClasses = "character", check.names = FALSE,
    na.strings = c("", "NA", "N", "."), skip = gp_genotype_table_header_line(input_csv) - 1L
  ))
  names(snp_data) <- sub("^#", "", names(snp_data))
  missing_headers <- setdiff(required_headers, names(snp_data))
  if (length(missing_headers)) {
    stop("Missing required VCF marker columns: ", paste(missing_headers, collapse = ", "),
         ". The table needs CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT followed by one column per sample.",
         call. = FALSE)
  }
  sample_names <- setdiff(names(snp_data), required_headers)
  if (!length(sample_names)) stop("The genotype table has no sample columns.", call. = FALSE)
  meta <- snp_data[, required_headers, drop = FALSE]
  meta[] <- lapply(meta, function(x) { x[is.na(x) | !nzchar(x)] <- "."; x })
  meta$FORMAT <- "GT"

  # All genotypes at once (the former per-marker loop took minutes on real panels).
  geno <- as.matrix(snp_data[, sample_names, drop = FALSE])
  values <- trimws(as.vector(geno))
  missing_gt <- paste(rep(".", ploidy), collapse = "/")
  gt_pattern <- if (ploidy == 1L) "^(0|1)$" else paste0("^(0|1)", strrep("[/|](0|1)", ploidy - 1L), "$")
  is_missing <- is.na(values)
  is_gt <- !is_missing & grepl(gt_pattern, values)
  is_dosage <- !is_missing & !is_gt
  dosage <- suppressWarnings(as.numeric(values[is_dosage]))
  if (identical(input_coding, "centered_dosage")) dosage <- dosage + ploidy / 2
  bad <- !is.finite(dosage) | dosage < 0 | dosage > ploidy | dosage != round(dosage)
  if (any(bad)) {
    stop("Unexpected genotype value for ", input_coding, " at ploidy ", ploidy, ": ",
         values[is_dosage][which(bad)[1L]], call. = FALSE)
  }
  dosage_to_gt <- vapply(0:ploidy, function(d) {
    paste(c(rep("0", ploidy - d), rep("1", d)), collapse = "/")
  }, character(1L))
  out <- rep(missing_gt, length(values))
  out[is_gt] <- values[is_gt]
  out[is_dosage] <- dosage_to_gt[as.integer(dosage) + 1L]
  out <- matrix(out, nrow = nrow(geno), ncol = ncol(geno))

  body <- do.call(paste, c(unname(as.list(meta)), lapply(seq_len(ncol(out)), function(j) out[, j]), sep = "\t"))
  vcf_content <- c(
    "##fileformat=VCFv4.3",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", required_headers[-1L], sample_names), collapse = "\t"),
    body
  )

  if (!is.null(output_vcf)) {
    writeLines(vcf_content, output_vcf)
    message("VCF file written: ", output_vcf)
    return(invisible(output_vcf))
  }
  vcf_content
}

# Line number of the column-header row (the first line not starting with "##").
gp_genotype_table_header_line <- function(path) {
  head_lines <- readLines(path, n = 200L, warn = FALSE)
  i <- which(!startsWith(head_lines, "##"))[1L]
  if (is.na(i)) 1L else i
}

# TRUE when a text file is a genotype table for convert_csv_to_vcf(): its
# header has the VCF marker columns CHROM, POS, REF and ALT (HapMap headers use
# lower-case chrom/pos and have no REF/ALT).
gp_is_genotype_table_file <- function(path) {
  if (!file.exists(path) || grepl("(?i)\\.vcf(\\.gz)?$", path, perl = TRUE)) return(FALSE)
  head_lines <- tryCatch(readLines(path, n = 200L, warn = FALSE), error = function(e) character())
  header <- head_lines[!startsWith(head_lines, "##")][1L]
  if (is.na(header)) return(FALSE)
  fields <- sub("^#", "", trimws(gsub("\"", "", strsplit(header, "[\t,;]")[[1L]])))
  all(c("CHROM", "POS", "REF", "ALT") %in% fields)
}
