gp_validate_unit_interval <- function(x, name, allow_null = TRUE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) || x < 0 || x > 1) {
    stop(name, " must be a single number between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

gp_normalize_hapmap_columns <- function(hapmap) {
  existing_aliases <- attr(hapmap, "hapmap_header_aliases_normalized", exact = TRUE)
  hapmap <- data.table::copy(data.table::as.data.table(hapmap))
  if (ncol(hapmap) < 12L) {
    stop("HapMap input must contain 11 metadata columns followed by at least one sample column.", call. = FALSE)
  }
  aliases <- list(
    c("rs#", "rs", "snp", "marker", "markerid", "marker_id"),
    "alleles",
    c("chrom", "chr", "chromosome"),
    c("pos", "position", "bp"),
    "strand",
    c("assembly#", "assembly"),
    "center",
    c("protLSID", "prot_lsid"),
    c("assayLSID", "assay_lsid"),
    c("panelLSID", "panel_lsid", "panel"),
    c("QCcode", "qc_code", "qc")
  )
  canonical <- c(
    "rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
    "protLSID", "assayLSID", "panelLSID", "QCcode"
  )
  actual <- names(hapmap)[1:11]
  normalize_name <- function(x) gsub("[^a-z0-9#]", "", tolower(x))
  valid <- vapply(seq_along(aliases), function(i) {
    normalize_name(actual[[i]]) %in% normalize_name(aliases[[i]])
  }, logical(1L))
  if (!all(valid)) {
    stop(
      "The first 11 HapMap columns are not recognized in the required order. Expected: ",
      paste(canonical, collapse = ", "), ". Found: ",
      paste(actual, collapse = ", "), ".",
      call. = FALSE
    )
  }
  renamed <- actual != canonical
  data.table::setnames(hapmap, old = actual, new = canonical)
  current_aliases <- stats::setNames(actual[renamed], canonical[renamed])
  attr(hapmap, "hapmap_header_aliases_normalized") <- c(
    existing_aliases[setdiff(names(existing_aliases), names(current_aliases))],
    current_aliases
  )
  hapmap
}

gp_hapmap_missing_codes <- function() {
  c("", "NA", "N", "NN", "B", "D", "H", "V", ".", "-", "?", "0")
}

gp_hapmap_missing_call <- function(x) {
  x <- toupper(trimws(as.character(x)))
  is.na(x) | x %in% gp_hapmap_missing_codes() |
    grepl("^[N.?-]+$", gsub("[/|:,_[:space:]]", "", x))
}

gp_hapmap_iupac_alleles <- function() {
  list(
    A = "A", C = "C", G = "G", T = "T", U = "T",
    R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"),
    W = c("A", "T"), K = c("G", "T"), M = c("A", "C")
  )
}

gp_parse_hapmap_alleles <- function(value) {
  value <- toupper(trimws(as.character(value)))
  value <- gsub("U", "T", value, fixed = TRUE)
  if (is.na(value) || !nzchar(value)) {
    return(character())
  }
  parts <- unlist(strsplit(value, "[/|:,_[:space:]]+"), use.names = FALSE)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 1L && nchar(parts) == 2L && grepl("^[ACGT]{2}$", parts)) {
    parts <- strsplit(parts, "", fixed = TRUE)[[1L]]
  }
  if (!length(parts) || any(!grepl("^[ACGT]+$", parts))) {
    return(character())
  }
  unique(parts)
}

gp_decode_hapmap_call <- function(value, declared, ploidy = 2L) {
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  value <- toupper(trimws(as.character(value)))
  value <- gsub("U", "T", value, fixed = TRUE)
  if (gp_hapmap_missing_call(value)) {
    return(list(dosage = NA_real_, missing = TRUE, valid = TRUE, ambiguous = FALSE))
  }

  compact <- gsub("[/|:,_[:space:]]", "", value)
  if (nchar(compact) == 1L && compact %in% declared) {
    return(list(
      dosage = if (identical(compact, declared[[2L]])) ploidy else 0,
      missing = FALSE, valid = TRUE, ambiguous = FALSE
    ))
  }
  if (ploidy == 2L && nchar(compact) == 1L) {
    decoded <- gp_hapmap_iupac_alleles()[[compact]]
    if (is.null(decoded)) {
      return(list(dosage = NA_real_, missing = FALSE, valid = FALSE, ambiguous = FALSE))
    }
    if (length(decoded) == 1L) {
      decoded <- rep(decoded, 2L)
    }
    valid <- all(decoded %in% declared)
    return(list(
      dosage = if (valid) sum(decoded == declared[[2L]]) else NA_real_,
      missing = FALSE, valid = valid, ambiguous = FALSE
    ))
  }
  alleles <- if (grepl("^[ACGT]+$", compact)) {
    strsplit(compact, "", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  if (length(alleles) && length(unique(alleles)) == 1L && alleles[[1L]] %in% declared) {
    return(list(
      dosage = if (identical(alleles[[1L]], declared[[2L]])) ploidy else 0,
      missing = FALSE, valid = TRUE, ambiguous = FALSE
    ))
  }
  if (length(alleles) == ploidy && all(alleles %in% declared)) {
    return(list(
      dosage = sum(alleles == declared[[2L]]),
      missing = FALSE, valid = TRUE, ambiguous = FALSE
    ))
  }
  ambiguous <- ploidy > 2L && nchar(compact) == 1L &&
    compact %in% names(gp_hapmap_iupac_alleles()) && !compact %in% declared
  list(dosage = NA_real_, missing = FALSE, valid = FALSE, ambiguous = ambiguous)
}

gp_resolve_hapmap_ploidy <- function(genotype_calls, ploidy = "auto") {
  requested <- gp_validate_ploidy(ploidy)
  if (!identical(requested, "auto")) return(requested)
  calls <- toupper(trimws(as.character(genotype_calls)))
  calls <- gsub("U", "T", calls, fixed = TRUE)
  compact <- gsub("[/|:,_[:space:]]", "", calls)
  missing <- gp_hapmap_missing_call(calls)
  candidate <- nchar(compact)
  candidate <- candidate[!missing & candidate > 2L & grepl("^[ACGT]+$", compact)]
  candidate <- sort(unique(candidate))
  if (!length(candidate)) return(2L)
  if (length(candidate) != 1L) {
    stop(
      "HapMap calls imply multiple ploidies (", paste(candidate, collapse = ", "),
      "). Supply one explicit crop ploidy and use complete allele-copy calls.",
      call. = FALSE
    )
  }
  as.integer(candidate[[1L]])
}

gp_hapmap_to_alt_dosage <- function(hapmap,
                                     invalid_call = c("error", "missing"),
                                     non_biallelic = c("drop", "error"),
                                     ploidy = "auto") {
  invalid_call <- match.arg(invalid_call)
  non_biallelic <- match.arg(non_biallelic)
  hapmap <- gp_normalize_hapmap_columns(hapmap)
  normalized_headers <- attr(hapmap, "hapmap_header_aliases_normalized", exact = TRUE)
  hapmap <- as.data.frame(hapmap, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(hapmap) < 12L) {
    stop("HapMap input must contain at least one sample column after the 11 metadata columns.", call. = FALSE)
  }

  sample_names <- names(hapmap)[12:ncol(hapmap)]
  if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop("HapMap sample names must be non-empty and unique.", call. = FALSE)
  }

  allele_pairs <- lapply(hapmap[[2L]], gp_parse_hapmap_alleles)
  biallelic <- vapply(allele_pairs, function(x) {
    length(x) == 2L && all(nchar(x) == 1L)
  }, logical(1L))
  if (any(!biallelic) && identical(non_biallelic, "error")) {
    bad <- which(!biallelic)
    stop(
      "HapMap markers must declare exactly two A/C/G/T alleles. Invalid row(s): ",
      paste(utils::head(bad, 10L), collapse = ", "),
      if (length(bad) > 10L) " ..." else "",
      call. = FALSE
    )
  }

  keep <- which(biallelic)
  if (!length(keep)) {
    stop("No biallelic HapMap markers remain for processing.", call. = FALSE)
  }
  hapmap_kept <- hapmap[keep, , drop = FALSE]
  allele_pairs <- allele_pairs[keep]
  dosage <- matrix(
    NA_real_,
    nrow = length(keep),
    ncol = length(sample_names),
    dimnames = list(as.character(hapmap_kept[[1L]]), sample_names)
  )
  genotype_calls <- as.matrix(hapmap_kept[, 12:ncol(hapmap_kept), drop = FALSE])
  ploidy <- gp_resolve_hapmap_ploidy(genotype_calls, ploidy)

  invalid <- list()
  invalid_count <- 0L
  ambiguous_count <- 0L
  pair_to_iupac <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")
  for (i in seq_len(nrow(hapmap_kept))) {
    declared <- allele_pairs[[i]]
    ref <- declared[[1L]]
    alt <- declared[[2L]]
    calls <- toupper(trimws(as.character(genotype_calls[i, ])))
    calls <- gsub("U", "T", calls, fixed = TRUE)
    compact <- gsub("[/|:,_[:space:]]", "", calls)
    missing <- gp_hapmap_missing_call(calls)
    is_ref <- grepl(paste0("^", ref, "+$"), compact)
    is_alt <- grepl(paste0("^", alt, "+$"), compact)
    pair_key <- paste0(sort(declared), collapse = "")
    full_call <- nchar(compact) == ploidy & grepl("^[ACGT]+$", compact)
    declared_only <- grepl(paste0("^[", ref, alt, "]+$"), compact)
    explicit_polyploid <- full_call & declared_only
    iupac_het <- compact %in% unname(pair_to_iupac[[pair_key]])
    diploid_pair <- ploidy == 2L & compact %in% c(paste0(ref, alt), paste0(alt, ref))
    is_het <- (ploidy == 2L & iupac_het) | diploid_pair
    ambiguous <- !missing & ploidy > 2L & iupac_het
    ambiguous_count <- ambiguous_count + sum(ambiguous)
    valid <- missing | is_ref | is_alt | is_het | explicit_polyploid
    if (any(!valid)) {
      bad <- which(!valid)
      invalid_count <- invalid_count + length(bad)
      room <- max(0L, 10L - length(invalid))
      if (room > 0L) {
        for (j in utils::head(bad, room)) {
          invalid[[length(invalid) + 1L]] <- c(row = keep[[i]], sample = sample_names[[j]])
        }
      }
    }
    dosage[i, is_ref & !missing] <- 0
    dosage[i, is_het & !missing] <- 1
    dosage[i, is_alt & !missing] <- 2
    if (ploidy != 2L) {
      dosage[i, is_alt & !missing] <- ploidy
      poly_idx <- which(explicit_polyploid & !is_ref & !is_alt & !missing)
      if (length(poly_idx)) {
        dosage[i, poly_idx] <- nchar(gsub(alt, "", compact[poly_idx], fixed = TRUE))
        dosage[i, poly_idx] <- ploidy - dosage[i, poly_idx]
      }
    }
  }

  if (invalid_count > 0L && identical(invalid_call, "error")) {
    labels <- vapply(invalid, function(x) {
      paste0("row ", x[["row"]], " sample ", x[["sample"]])
    }, character(1L))
    stop(
      "HapMap genotype calls contain alleles that are invalid, dosage-ambiguous for ploidy ",
      ploidy, ", or inconsistent with the declared allele pair: ",
      paste(labels, collapse = "; "),
      if (invalid_count > 10L) "; ..." else "",
      call. = FALSE
    )
  }

  list(
    dosage = dosage,
    map = data.table::as.data.table(hapmap_kept[, 1:11, drop = FALSE]),
    allele_pairs = allele_pairs,
    retained_rows = keep,
    report = list(
      input_markers = nrow(hapmap),
      retained_markers = nrow(hapmap_kept),
      dropped_non_biallelic = sum(!biallelic),
      invalid_calls_set_missing = if (identical(invalid_call, "missing")) invalid_count else 0L,
      ambiguous_polyploid_calls = ambiguous_count,
      samples = length(sample_names),
      ploidy = ploidy,
      normalized_header_aliases = normalized_headers,
      dosage_allele = "second declared HapMap allele"
    )
  )
}

gp_parse_vcf_alt_dosage <- function(vcf_data, ploidy = "auto") {
  cols <- vcf_required_columns(vcf_data)
  sample_names <- names(vcf_data)[10:ncol(vcf_data)]
  if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop("VCF sample names must be non-empty and unique.", call. = FALSE)
  }
  dosage <- matrix(
    NA_real_, nrow = nrow(vcf_data), ncol = length(sample_names),
    dimnames = list(as.character(vcf_data[[3L]]), sample_names)
  )
  formats <- as.character(vcf_data[[cols$format]])
  formats[is.na(formats) | !nzchar(trimws(formats))] <- "GT"
  sample_block <- if (inherits(vcf_data, "data.table")) {
    vcf_data[, sample_names, with = FALSE]
  } else {
    vcf_data[, sample_names, drop = FALSE]
  }
  sample_matrix <- as.matrix(sample_block)
  gt_matrix <- matrix(
    NA_character_, nrow = nrow(vcf_data), ncol = length(sample_names),
    dimnames = list(NULL, sample_names)
  )
  for (format_value in unique(formats)) {
    rows <- which(formats == format_value)
    fmt <- strsplit(format_value, ":", fixed = TRUE)[[1L]]
    gt_index <- match("GT", fmt)
    if (is.na(gt_index)) {
      stop("VCF row ", rows[[1L]], " does not contain a GT FORMAT field.", call. = FALSE)
    }
    values <- as.character(sample_matrix[rows, , drop = FALSE])
    gt <- if (identical(gt_index, 1L)) {
      sub(":.*$", "", values)
    } else {
      vapply(strsplit(values, ":", fixed = TRUE), function(parts) {
        if (length(parts) >= gt_index) parts[[gt_index]] else NA_character_
      }, character(1L))
    }
    gt_matrix[rows, ] <- matrix(gt, nrow = length(rows), ncol = length(sample_names))
  }
  ploidy <- gp_resolve_gt_ploidy(as.vector(gt_matrix), ploidy, context = "VCF")
  pattern <- if (ploidy == 1L) {
    "^(0|1|\\.)$"
  } else {
    paste0("^(0|1|\\.)", paste(rep("[/|](0|1|\\.)", ploidy - 1L), collapse = ""), "$")
  }
  for (i in seq_len(nrow(gt_matrix))) {
    gt <- gt_matrix[i, ]
    missing <- is.na(gt) | gt == "." | grepl("\\.", gt)
    invalid <- !missing & !grepl(pattern, gt)
    if (any(invalid)) {
      first <- which(invalid)[[1L]]
      stop(
        "VCF row ", i, " sample '", sample_names[[first]],
        "' is not a ploidy-", ploidy, " biallelic GT call: ", gt[[first]],
        call. = FALSE
      )
    }
    values_dosage <- rep(NA_real_, length(gt))
    called <- which(!missing)
    if (length(called)) {
      values_dosage[called] <- vapply(
        strsplit(gt[called], "[/|]", perl = TRUE),
        function(alleles) sum(alleles == "1"), numeric(1L)
      )
    }
    dosage[i, ] <- values_dosage
  }
  gp_attach_ploidy(dosage, ploidy)
}

gp_qc_alt_dosage <- function(dosage,
                              map,
                              marker_id_col = 1L,
                              maf_threshold = 0.01,
                              het_threshold = 0.1,
                              ind_call_rate_threshold = 0.9,
                              snp_call_rate_threshold = 0.9,
                              impute = TRUE,
                              recode_format = "0,1,2",
                              ploidy = "auto",
                              imputation_method = c("mean", "median", "mode", "knn"),
                              impute_knn_k = 5L) {
  maf_threshold <- gp_validate_unit_interval(maf_threshold, "maf_threshold")
  het_threshold <- gp_validate_unit_interval(het_threshold, "het_threshold")
  ind_call_rate_threshold <- gp_validate_unit_interval(
    ind_call_rate_threshold, "ind_call_rate_threshold"
  )
  snp_call_rate_threshold <- gp_validate_unit_interval(
    snp_call_rate_threshold, "snp_call_rate_threshold"
  )
  dosage <- as.matrix(dosage)
  storage.mode(dosage) <- "double"
  ploidy <- gp_resolve_matrix_ploidy(dosage, ploidy)
  gp_validate_alt_dosage(dosage, ploidy, hard_calls = FALSE)
  recode_contract <- gp_normalize_recode_format(recode_format)
  imputation_method <- match.arg(tolower(imputation_method), c("mean", "median", "mode", "knn"))
  if (nrow(dosage) != nrow(map)) {
    stop("The dosage rows and marker-map rows are not aligned.", call. = FALSE)
  }
  input_markers <- nrow(dosage)
  input_samples <- ncol(dosage)

  observed_levels <- apply(dosage, 1L, function(x) length(unique(x[!is.na(x)])))
  monomorphic <- observed_levels <= 1L
  monomorphic[is.na(monomorphic)] <- TRUE
  total_monomorphic <- sum(monomorphic)
  if (any(monomorphic)) {
    dosage <- dosage[!monomorphic, , drop = FALSE]
    map <- map[!monomorphic, , drop = FALSE]
  }

  marker_call_rate <- rowMeans(!is.na(dosage))
  low_marker <- if (is.null(snp_call_rate_threshold)) rep(FALSE, nrow(dosage)) else
    marker_call_rate < snp_call_rate_threshold
  total_low_marker <- sum(low_marker)
  if (any(low_marker)) {
    dosage <- dosage[!low_marker, , drop = FALSE]
    map <- map[!low_marker, , drop = FALSE]
  }
  if (!nrow(dosage)) {
    stop("No markers remain after monomorphic and SNP call-rate filtering.", call. = FALSE)
  }

  individual_call_rate <- colMeans(!is.na(dosage))
  low_individual <- if (is.null(ind_call_rate_threshold)) rep(FALSE, ncol(dosage)) else
    individual_call_rate < ind_call_rate_threshold
  removed_individuals <- colnames(dosage)[low_individual]
  if (any(low_individual)) {
    dosage <- dosage[, !low_individual, drop = FALSE]
  }
  if (!ncol(dosage)) {
    stop("No samples remain after individual call-rate filtering.", call. = FALSE)
  }

  heterozygosity <- apply(dosage, 1L, function(x) {
    called <- !is.na(x)
    if (!any(called)) NA_real_ else mean(x[called] > 0 & x[called] < ploidy)
  })
  high_heterozygosity <- if (is.null(het_threshold)) rep(FALSE, nrow(dosage)) else
    !is.na(heterozygosity) & heterozygosity > het_threshold
  total_high_heterozygosity <- sum(high_heterozygosity)
  if (any(high_heterozygosity)) {
    dosage <- dosage[!high_heterozygosity, , drop = FALSE]
    map <- map[!high_heterozygosity, , drop = FALSE]
  }

  alt_frequency <- rowMeans(dosage, na.rm = TRUE) / ploidy
  minor_allele_frequency <- pmin(alt_frequency, 1 - alt_frequency)
  low_maf <- if (is.null(maf_threshold)) rep(FALSE, nrow(dosage)) else
    is.na(minor_allele_frequency) | minor_allele_frequency < maf_threshold
  total_low_maf <- sum(low_maf)
  if (any(low_maf)) {
    dosage <- dosage[!low_maf, , drop = FALSE]
    map <- map[!low_maf, , drop = FALSE]
  }
  if (!nrow(dosage)) {
    stop("No markers remain after heterozygosity and MAF filtering.", call. = FALSE)
  }

  imputed_cells <- 0L
  if (isTRUE(impute) && anyNA(dosage)) {
    imputed <- gp_impute_alt_dosage(
      dosage, ploidy = ploidy, method = imputation_method, k = impute_knn_k
    )
    dosage <- imputed$dosage
    imputed_cells <- imputed$imputed_cells
  }

  marker_names <- as.character(map[[marker_id_col]])
  if (any(!nzchar(marker_names)) || anyDuplicated(marker_names)) {
    stop("Retained marker IDs must be non-empty and unique.", call. = FALSE)
  }
  result <- gp_recode_alt_dosage(t(dosage), ploidy, recode_contract)
  colnames(result) <- marker_names
  class(result) <- c("matrix", "array", "genotype")

  metrics <- data.frame(
    metric = c(
      "snp_call_rate_threshold", "total_snp_removed",
      "ind_call_rate_threshold", "total_genotypes_removed",
      "het_threshold", "total_het_snps_removed", "maf_threshold",
      "total_maf_snps_remove", "total_monomorphic_snps_removed",
      "input_markers", "retained_markers", "input_samples",
      "retained_samples", "imputed_cells", "ploidy"
    ),
    stat = c(
      snp_call_rate_threshold %||% NA_real_, total_low_marker,
      ind_call_rate_threshold %||% NA_real_, length(removed_individuals),
      het_threshold %||% NA_real_, total_high_heterozygosity,
      maf_threshold %||% NA_real_,
      total_low_maf, total_monomorphic, input_markers, nrow(dosage),
      input_samples, ncol(dosage), imputed_cells, ploidy
    ),
    stringsAsFactors = FALSE
  )

  list(
    snps_matrix = result,
    snp_map = as.data.frame(map, check.names = FALSE),
    qc_metrics_and_summary_stat = metrics,
    removed_samples = removed_individuals,
    coding = attr(result, "genotype_coding"),
    coding_contract = attr(result, "genotype_coding_contract"),
    ploidy = ploidy,
    imputation_method = if (isTRUE(impute)) imputation_method else NULL
  )
}
