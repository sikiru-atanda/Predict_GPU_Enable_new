vcf_is_gz <- function(path) {
  grepl("\\.gz$", path, ignore.case = TRUE)
}

vcf_connection <- function(path, open = "rt") {
  if (vcf_is_gz(path)) {
    gzfile(path, open = open)
  } else {
    file(path, open = open)
  }
}

vcf_read_header <- function(path) {
  con <- vcf_connection(path, open = "rt")
  on.exit(close(con), add = TRUE)

  meta_lines <- character()
  header_line <- NULL
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      break
    }
    if (startsWith(line, "##")) {
      meta_lines <- c(meta_lines, line)
    } else if (startsWith(line, "#CHROM") || startsWith(line, "#chrom")) {
      header_line <- line
      break
    }
  }

  if (is.null(header_line)) {
    stop("VCF header line beginning with #CHROM was not found.", call. = FALSE)
  }

  version_lines <- grep("^##fileformat=VCFv", meta_lines, value = TRUE)
  if (length(version_lines) != 1L ||
      !grepl("^##fileformat=VCFv4\\.[1-5]$", version_lines[[1L]])) {
    stop(
      "VCF input must declare exactly one supported fileformat line (VCFv4.1 through VCFv4.5).",
      call. = FALSE
    )
  }

  list(
    meta_lines = meta_lines,
    header_line = header_line,
    col_names = strsplit(header_line, "\t", fixed = TRUE)[[1L]],
    fileformat = sub("^##fileformat=", "", version_lines[[1L]])
  )
}

vcf_read_all_lines <- function(path) {
  con <- vcf_connection(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

vcf_empty_table <- function(col_names) {
  data.table::as.data.table(stats::setNames(rep(list(character()), length(col_names)), col_names))
}

vcf_read_table <- function(path) {
  header <- vcf_read_header(path)
  table <- tryCatch(
    data.table::fread(
      path,
      skip = "#CHROM",
      header = TRUE,
      sep = "\t",
      quote = "",
      fill = TRUE,
      check.names = FALSE,
      data.table = TRUE,
      colClasses = "character",
      na.strings = "NA",
      showProgress = FALSE
    ),
    error = function(e) {
      lines <- vcf_read_all_lines(path)
      header_idx <- grep("^#CHROM\t|^#CHROM$", lines)
      if (!length(header_idx)) {
        header_idx <- grep("^#chrom\t|^#chrom$", lines, ignore.case = TRUE)
      }
      if (!length(header_idx)) {
        stop("VCF header line beginning with #CHROM was not found.", call. = FALSE)
      }
      body_lines <- if (header_idx[[1L]] < length(lines)) {
        lines[(header_idx[[1L]] + 1L):length(lines)]
      } else {
        character()
      }
      body_lines <- body_lines[nzchar(body_lines)]
      if (!length(body_lines)) {
        return(vcf_empty_table(header$col_names))
      }
      data.table::fread(
        text = paste(c(lines[[header_idx[[1L]]]], body_lines), collapse = "\n"),
        header = TRUE,
        sep = "\t",
        quote = "",
        fill = TRUE,
        check.names = FALSE,
        data.table = TRUE,
        colClasses = "character",
        na.strings = "NA",
        showProgress = FALSE
      )
    }
  )

  if (!ncol(table)) {
    table <- vcf_empty_table(header$col_names)
  }
  table
}

vcf_write_table <- function(vcf_data, output_vcf, meta_lines = character()) {
  dir.create(dirname(normalizePath(output_vcf, winslash = "/", mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
  if (length(meta_lines)) {
    writeLines(meta_lines, output_vcf, useBytes = TRUE)
    append <- TRUE
  } else {
    append <- FALSE
  }
  data.table::fwrite(
    vcf_data,
    file = output_vcf,
    sep = "\t",
    quote = FALSE,
    na = ".",
    append = append,
    col.names = TRUE
  )
  normalizePath(output_vcf, winslash = "/", mustWork = TRUE)
}

vcf_find_column <- function(vcf_data, candidates) {
  nms <- names(vcf_data)
  idx <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) {
    return(NULL)
  }
  nms[[idx[[1L]]]]
}

vcf_required_columns <- function(vcf_data) {
  cols <- list(
    chrom = vcf_find_column(vcf_data, c("#CHROM", "CHROM")),
    pos = vcf_find_column(vcf_data, "POS"),
    ref = vcf_find_column(vcf_data, "REF"),
    alt = vcf_find_column(vcf_data, "ALT"),
    format = vcf_find_column(vcf_data, "FORMAT")
  )
  missing_cols <- names(cols)[vapply(cols, is.null, logical(1))]
  if (length(missing_cols)) {
    stop("VCF is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (ncol(vcf_data) < 10L) {
    stop("VCF must contain at least one sample genotype column for PLINK/Beagle processing.", call. = FALSE)
  }
  cols
}

vcf_valid_pos <- function(pos) {
  pos_chr <- trimws(as.character(pos))
  !is.na(pos_chr) & grepl("^[0-9]+$", pos_chr) & suppressWarnings(as.numeric(pos_chr) > 0)
}

vcf_valid_allele <- function(x) {
  x <- toupper(trimws(as.character(x)))
  !is.na(x) & nzchar(x) & x != "." & grepl("^[ACGT]+$", x)
}

vcf_extract_gt_field <- function(values, gt_index) {
  values <- as.character(values)
  if (identical(gt_index, 1L)) {
    return(sub(":.*$", "", values))
  }
  vapply(strsplit(values, ":", fixed = TRUE), function(parts) {
    if (length(parts) >= gt_index) {
      parts[[gt_index]]
    } else {
      NA_character_
    }
  }, character(1))
}

vcf_gt_supported <- function(gt, ploidy = 2L) {
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  gt <- trimws(as.character(gt))
  gt[is.na(gt)] <- "."
  pattern <- if (ploidy == 1L) {
    "^(0|1|\\.)$"
  } else {
    paste0("^(0|1|\\.)", paste(rep("[/|](0|1|\\.)", ploidy - 1L), collapse = ""), "$")
  }
  gt == "." | grepl(pattern, gt)
}

vcf_collect_gt <- function(vcf_data, format_col, sample_cols) {
  formats <- as.character(vcf_data[[format_col]])
  formats[is.na(formats) | !nzchar(trimws(formats))] <- "GT"
  out <- matrix(NA_character_, nrow = nrow(vcf_data), ncol = length(sample_cols))
  for (fmt in unique(formats)) {
    rows <- which(formats == fmt)
    fmt_parts <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
    gt_index <- match("GT", fmt_parts)
    if (is.na(gt_index)) next
    for (j in seq_along(sample_cols)) {
      out[rows, j] <- vcf_extract_gt_field(vcf_data[[sample_cols[[j]]]][rows], gt_index)
    }
  }
  out
}

vcf_unsupported_gt_rows <- function(vcf_data, format_col, sample_cols, ploidy = 2L) {
  n <- nrow(vcf_data)
  invalid <- rep(FALSE, n)
  formats <- as.character(vcf_data[[format_col]])
  formats[is.na(formats) | !nzchar(trimws(formats))] <- "GT"

  for (fmt in unique(formats)) {
    rows <- which(formats == fmt)
    fmt_parts <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
    gt_index <- match("GT", fmt_parts)
    if (is.na(gt_index)) {
      invalid[rows] <- TRUE
      next
    }
    for (sample_col in sample_cols) {
      gt <- vcf_extract_gt_field(vcf_data[[sample_col]][rows], gt_index)
      invalid[rows] <- invalid[rows] | !vcf_gt_supported(gt, ploidy = ploidy)
    }
  }
  invalid
}

vcf_output_path <- function(input_vcf, output_vcf = NULL, suffix = "predictpro_sanitized") {
  if (!is.null(output_vcf) && nzchar(output_vcf)) {
    return(normalizePath(output_vcf, winslash = "/", mustWork = FALSE))
  }
  tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(sub("\\.gz$", "", input_vcf, ignore.case = TRUE))), "_", suffix, "_"),
    tmpdir = tempdir(),
    fileext = ".vcf"
  )
}

sanitize_vcf_for_external_tools <- function(input_vcf,
                                            output_vcf = NULL,
                                            mode = c("plink", "beagle", "recode"),
                                            write_if_changed = TRUE,
                                            strict = FALSE,
                                            ploidy = "auto") {
  mode <- match.arg(mode)
  input_vcf <- normalizePath(input_vcf, winslash = "/", mustWork = TRUE)
  header <- vcf_read_header(input_vcf)
  vcf_data <- vcf_read_table(input_vcf)
  cols <- vcf_required_columns(vcf_data)

  metrics <- list(
    mode = mode,
    input_variants = nrow(vcf_data),
    output_variants = nrow(vcf_data),
    normalized_fields = 0L,
    removed_bad_chromosome = 0L,
    removed_bad_position = 0L,
    removed_missing_allele = 0L,
    removed_ref_alt_equal = 0L,
    removed_non_acgt_or_multiallelic = 0L,
    removed_unsupported_genotype = 0L,
    removed_conflicting_alleles = 0L,
    removed_duplicate_marker = 0L
  )

  if (!nrow(vcf_data)) {
    stop("VCF contains no variants.", call. = FALSE)
  }

  sample_cols <- setdiff(names(vcf_data)[10:ncol(vcf_data)], character())
  requested_ploidy <- gp_validate_ploidy(ploidy)
  if (mode %in% c("plink", "beagle")) {
    if (!identical(requested_ploidy, "auto") && requested_ploidy != 2L) {
      stop(mode, " processing is restricted to diploid GT calls.", call. = FALSE)
    }
    resolved_ploidy <- 2L
  } else {
    gt_matrix <- vcf_collect_gt(vcf_data, cols$format, sample_cols)
    resolved_ploidy <- gp_resolve_gt_ploidy(
      as.vector(gt_matrix), requested_ploidy, context = "VCF recoding"
    )
  }
  metrics$ploidy <- resolved_ploidy
  keep <- rep(TRUE, nrow(vcf_data))

  drop_rows <- function(mask, metric) {
    mask <- mask & keep
    metrics[[metric]] <<- as.integer(sum(mask, na.rm = TRUE))
    keep[mask] <<- FALSE
  }

  chrom_raw <- as.character(vcf_data[[cols$chrom]])
  chrom_clean <- trimws(chrom_raw)
  bad_chrom <- is.na(chrom_clean) | !nzchar(chrom_clean) | chrom_clean == "." | grepl("\\s", chrom_clean)
  drop_rows(bad_chrom, "removed_bad_chromosome")

  pos_raw <- as.character(vcf_data[[cols$pos]])
  pos_clean <- trimws(pos_raw)
  bad_pos <- !vcf_valid_pos(pos_clean)
  drop_rows(bad_pos, "removed_bad_position")

  ref_raw <- as.character(vcf_data[[cols$ref]])
  alt_raw <- as.character(vcf_data[[cols$alt]])
  ref_clean <- toupper(trimws(ref_raw))
  alt_clean <- toupper(trimws(alt_raw))

  missing_allele <- is.na(ref_clean) | is.na(alt_clean) |
    !nzchar(ref_clean) | !nzchar(alt_clean) |
    ref_clean == "." | alt_clean == "."
  drop_rows(missing_allele, "removed_missing_allele")

  ref_alt_equal <- ref_clean == alt_clean
  ref_alt_equal[is.na(ref_alt_equal)] <- FALSE
  drop_rows(ref_alt_equal, "removed_ref_alt_equal")

  invalid_allele <- !vcf_valid_allele(ref_clean) | !vcf_valid_allele(alt_clean) |
    grepl(",", alt_clean, fixed = TRUE) |
    grepl("[<>\\[\\]*]", alt_clean)
  drop_rows(invalid_allele, "removed_non_acgt_or_multiallelic")

  gt_bad <- vcf_unsupported_gt_rows(
    vcf_data, cols$format, sample_cols, ploidy = resolved_ploidy
  )
  drop_rows(gt_bad, "removed_unsupported_genotype")

  coord_key <- paste(chrom_clean, pos_clean, sep = "\r")
  allele_key <- paste(ref_clean, alt_clean, sep = "\r")
  kept_idx <- which(keep)
  if (length(kept_idx)) {
    coord_has_conflict <- vapply(split(allele_key[kept_idx], coord_key[kept_idx]), function(x) {
      length(unique(x)) > 1L
    }, logical(1))
    conflict_coords <- names(coord_has_conflict)[coord_has_conflict]
    if (length(conflict_coords)) {
      drop_rows(coord_key %in% conflict_coords, "removed_conflicting_alleles")
    }
  }

  kept_idx <- which(keep)
  if (length(kept_idx)) {
    marker_key <- paste(chrom_clean, pos_clean, ref_clean, alt_clean, sep = "\r")
    duplicate_rows <- rep(FALSE, nrow(vcf_data))
    duplicate_rows[kept_idx] <- duplicated(marker_key[kept_idx])
    drop_rows(duplicate_rows, "removed_duplicate_marker")
  }

  normalized <- (chrom_clean != chrom_raw) |
    (pos_clean != pos_raw) |
    (ref_clean != ref_raw) |
    (alt_clean != alt_raw)
  normalized[is.na(normalized)] <- FALSE
  metrics$normalized_fields <- as.integer(sum(normalized & keep))

  changed <- any(!keep) || metrics$normalized_fields > 0L
  if (isTRUE(strict) && changed) {
    stop("VCF preflight found rows or fields that require cleanup before ", mode, " processing.", call. = FALSE)
  }

  if (!any(keep)) {
    stop("VCF preflight removed all variants; check chromosome, position, allele, and genotype coding.", call. = FALSE)
  }

  if (!changed || !isTRUE(write_if_changed)) {
    metrics$output_variants <- as.integer(sum(keep))
    return(list(path = input_vcf, changed = changed, metrics = metrics))
  }

  cleaned <- data.table::as.data.table(vcf_data[keep, , drop = FALSE])
  cleaned[[cols$chrom]] <- chrom_clean[keep]
  cleaned[[cols$pos]] <- pos_clean[keep]
  cleaned[[cols$ref]] <- ref_clean[keep]
  cleaned[[cols$alt]] <- alt_clean[keep]
  for (sample_col in sample_cols) {
    x <- as.character(cleaned[[sample_col]])
    x[is.na(x) | !nzchar(trimws(x))] <- "."
    cleaned[[sample_col]] <- x
  }

  output_path <- vcf_output_path(input_vcf, output_vcf)
  output_path <- vcf_write_table(cleaned, output_path, meta_lines = header$meta_lines)
  metrics$output_variants <- nrow(cleaned)
  list(path = output_path, changed = TRUE, metrics = metrics)
}
