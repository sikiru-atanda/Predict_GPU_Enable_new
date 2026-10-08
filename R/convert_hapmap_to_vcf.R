#' Convert HapMap genotypes to a Beagle-compatible VCF
#'
#' Converts the standard 11-column HapMap representation to VCF 4.3 without
#' changing allele orientation. The first declared HapMap allele becomes REF,
#' the second becomes ALT, and each GT call is encoded relative to that pair.
#' This makes conversion and PredictProR numeric recoding use the same
#' alternate-allele dosage definition.
#'
#' @param input_hapmap A HapMap filename or an in-memory HapMap data frame.
#' @param output_vcf Optional output filename. It must end in `.vcf` or
#'   `.vcf.gz`. If omitted, the VCF lines are returned.
#' @param invalid_call Whether invalid or allele-inconsistent calls cause an
#'   error or are written as missing.
#' @param non_biallelic Whether markers that do not declare exactly two A/C/G/T
#'   alleles are dropped or cause an error.
#' @param ploidy One positive integer or `"auto"`. Polyploid calls must provide
#'   unambiguous allele-copy dosage, such as `AAGG` for dosage two at ploidy
#'   four. IUPAC heterozygote codes are dosage-ambiguous above ploidy two.
#'
#' @return If `output_vcf` is `NULL`, a character vector of VCF lines with a
#'   `conversion_report` attribute. Otherwise, the normalized output path is
#'   returned invisibly with the same report attribute.
#' @export
convert_hapmap_to_vcf <- function(input_hapmap,
                                  output_vcf = NULL,
                                  invalid_call = c("error", "missing"),
                                  non_biallelic = c("drop", "error"),
                                  ploidy = "auto") {
  invalid_call <- match.arg(invalid_call)
  non_biallelic <- match.arg(non_biallelic)
  hapmap <- if (is.character(input_hapmap) && length(input_hapmap) == 1L) {
    read_hapmap_file(input_hapmap)
  } else if (is.data.frame(input_hapmap)) {
    data.table::as.data.table(input_hapmap)
  } else {
    stop("input_hapmap must be a HapMap filename or data frame.", call. = FALSE)
  }

  parsed <- gp_hapmap_to_alt_dosage(
    hapmap,
    invalid_call = invalid_call,
    non_biallelic = non_biallelic,
    ploidy = ploidy
  )
  map <- parsed$map
  dosage <- parsed$dosage
  ploidy <- parsed$report$ploidy
  sample_names <- colnames(dosage)

  missing_gt <- paste(rep(".", ploidy), collapse = "/")
  gt <- matrix(missing_gt, nrow = nrow(dosage), ncol = ncol(dosage))
  gt_by_dosage <- vapply(0:ploidy, function(alt_count) {
    paste(c(rep("0", ploidy - alt_count), rep("1", alt_count)), collapse = "/")
  }, character(1L))
  for (alt_count in 0:ploidy) {
    gt[!is.na(dosage) & dosage == alt_count] <- gt_by_dosage[[alt_count + 1L]]
  }

  allele_pairs <- parsed$allele_pairs
  chrom <- as.character(map[[3L]])
  pos <- suppressWarnings(as.integer(as.character(map[[4L]])))
  marker_id <- as.character(map[[1L]])
  invalid_coord <- is.na(chrom) | !nzchar(trimws(chrom)) | grepl("\\s", chrom) |
    is.na(pos) | pos < 1L
  if (any(invalid_coord)) {
    stop(
      "HapMap chrom and pos must define non-empty chromosomes and positive integer positions. Invalid row(s): ",
      paste(utils::head(parsed$retained_rows[invalid_coord], 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(is.na(marker_id) | !nzchar(marker_id)) || anyDuplicated(marker_id)) {
    stop("HapMap marker IDs must be non-empty and unique for VCF conversion.", call. = FALSE)
  }

  vcf <- data.table::data.table(
    `#CHROM` = chrom,
    POS = as.character(pos),
    ID = marker_id,
    REF = vapply(allele_pairs, `[[`, character(1L), 1L),
    ALT = vapply(allele_pairs, `[[`, character(1L), 2L),
    QUAL = ".",
    FILTER = "PASS",
    INFO = ".",
    FORMAT = "GT"
  )
  for (j in seq_along(sample_names)) {
    vcf[[sample_names[[j]]]] <- gt[, j]
  }
  meta <- c(
    "##fileformat=VCFv4.3",
    "##source=PredictProR_HapMap_to_VCF",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  )
  header <- paste(names(vcf), collapse = "\t")
  body <- do.call(paste, c(as.data.frame(vcf), sep = "\t"))
  lines <- c(meta, header, body)
  report <- parsed$report
  report$output_format <- "VCFv4.3"
  report$ref_allele <- "first declared HapMap allele"
  report$alt_allele <- "second declared HapMap allele"
  report$ploidy <- ploidy
  attr(lines, "conversion_report") <- report

  if (is.null(output_vcf)) {
    return(lines)
  }
  if (!grepl("(?i)\\.vcf(\\.gz)?$", output_vcf, perl = TRUE)) {
    stop("output_vcf must end in .vcf or .vcf.gz.", call. = FALSE)
  }
  output_vcf <- normalizePath(output_vcf, winslash = "/", mustWork = FALSE)
  dir.create(dirname(output_vcf), recursive = TRUE, showWarnings = FALSE)
  con <- if (grepl("(?i)\\.gz$", output_vcf, perl = TRUE)) {
    gzfile(output_vcf, open = "wt")
  } else {
    file(output_vcf, open = "wt")
  }
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, useBytes = TRUE)
  close(con)
  on.exit(NULL, add = FALSE)
  output_vcf <- normalizePath(output_vcf, winslash = "/", mustWork = TRUE)
  attr(output_vcf, "conversion_report") <- report
  invisible(output_vcf)
}
