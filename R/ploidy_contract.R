gp_validate_ploidy <- function(ploidy = 2L, allow_auto = TRUE) {
  if (isTRUE(allow_auto) && is.character(ploidy) &&
      length(ploidy) == 1L && identical(tolower(trimws(ploidy)), "auto")) {
    return("auto")
  }
  if (length(ploidy) != 1L || !is.numeric(ploidy) || is.na(ploidy) ||
      !is.finite(ploidy) || ploidy < 1 || ploidy != floor(ploidy)) {
    stop("ploidy must be 'auto' or one positive integer.", call. = FALSE)
  }
  as.integer(ploidy)
}

gp_resolve_matrix_ploidy <- function(x, ploidy = "auto") {
  ploidy <- gp_validate_ploidy(ploidy)
  attributed <- attr(x, "ploidy", exact = TRUE)
  if (!identical(ploidy, "auto")) {
    if (!is.null(attributed) && !identical(as.integer(attributed), ploidy)) {
      stop(
        "Explicit ploidy (", ploidy, ") conflicts with genotype metadata (",
        attributed, ").", call. = FALSE
      )
    }
    return(ploidy)
  }
  if (!is.null(attributed)) {
    return(gp_validate_ploidy(attributed, allow_auto = FALSE))
  }
  observed <- suppressWarnings(as.numeric(x))
  observed <- observed[is.finite(observed)]
  if (length(observed) && max(observed) > 2) {
    stop(
      "A numeric genotype matrix contains dosage above 2 but has no ploidy metadata. ",
      "Supply ploidy explicitly; the maximum observed dosage is not a reliable ploidy estimator.",
      call. = FALSE
    )
  }
  2L
}

gp_attach_ploidy <- function(x, ploidy) {
  attr(x, "ploidy") <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  x
}

gp_validate_alt_dosage <- function(x, ploidy, hard_calls = FALSE, name = "genotype dosage") {
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  values <- suppressWarnings(as.numeric(x))
  observed <- !is.na(values)
  invalid <- observed & (!is.finite(values) | values < 0 | values > ploidy)
  if (isTRUE(hard_calls)) {
    invalid <- invalid | (observed & abs(values - round(values)) > sqrt(.Machine$double.eps))
  }
  if (any(invalid)) {
    example <- values[which(invalid)[[1L]]]
    stop(
      name, " must be ", if (isTRUE(hard_calls)) "integer " else "numeric ",
      "ALT-allele dosage in [0, ", ploidy, "]; found ", example, ".",
      call. = FALSE
    )
  }
  invisible(ploidy)
}

gp_normalize_recode_format <- function(recode_format) {
  if (is.null(recode_format)) return("alt_dosage")
  value <- tolower(trimws(as.character(recode_format[[1L]])))
  aliases <- c(
    "0,1,2" = "alt_dosage",
    "dosage" = "alt_dosage",
    "alt_dosage" = "alt_dosage",
    "-1,0,1" = "centered_dosage",
    "centered" = "centered_dosage",
    "centered_dosage" = "centered_dosage",
    "frequency" = "allele_frequency",
    "allele_frequency" = "allele_frequency"
  )
  normalized <- unname(aliases[[value]])
  if (is.null(normalized)) {
    stop(
      "Invalid recode format. Use 'alt_dosage', 'centered_dosage', ",
      "'allele_frequency', or the legacy diploid aliases '0,1,2' and '-1,0,1'.",
      call. = FALSE
    )
  }
  normalized
}

gp_recode_alt_dosage <- function(x, ploidy, recode_format = "alt_dosage") {
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  coding <- gp_normalize_recode_format(recode_format)
  out <- switch(
    coding,
    alt_dosage = x,
    centered_dosage = x - ploidy / 2,
    allele_frequency = x / ploidy
  )
  attr(out, "ploidy") <- ploidy
  attr(out, "genotype_coding") <- switch(
    coding,
    alt_dosage = "alternate_allele_dosage",
    centered_dosage = "alternate_allele_dosage_centered_on_half_ploidy",
    allele_frequency = "alternate_allele_frequency"
  )
  attr(out, "genotype_coding_contract") <- coding
  out
}

gp_impute_alt_dosage <- function(dosage,
                                  ploidy,
                                  method = c("mean", "median", "mode", "knn"),
                                  k = 5L) {
  method <- match.arg(tolower(method), c("mean", "median", "mode", "knn"))
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  dosage <- as.matrix(dosage)
  storage.mode(dosage) <- "double"
  gp_validate_alt_dosage(dosage, ploidy, hard_calls = FALSE)
  missing_before <- is.na(dosage)
  if (!any(missing_before)) {
    return(list(dosage = gp_attach_ploidy(dosage, ploidy), imputed_cells = 0L, method = method))
  }
  all_missing <- which(rowSums(!is.na(dosage)) == 0L)
  if (length(all_missing)) {
    stop(
      "Cannot impute marker rows with no observed dosage. Row(s): ",
      paste(utils::head(all_missing, 10L), collapse = ", "),
      if (length(all_missing) > 10L) " ..." else "", call. = FALSE
    )
  }

  if (identical(method, "knn")) {
    samples_by_markers <- t(dosage)
    if (geno_impute_knn_cpp_available()) {
      samples_by_markers <- geno_impute_knn_cpp(samples_by_markers, k = k)
    } else {
      samples_by_markers <- handle_large_scale_knn(
        data = samples_by_markers, k = k, chunk_size = 50, num_cores = NULL
      )
    }
    dosage <- t(samples_by_markers)
  } else {
    replacement <- apply(dosage, 1L, function(x) {
      observed <- x[!is.na(x)]
      switch(
        method,
        mean = mean(observed),
        median = stats::median(observed),
        mode = {
          counts <- table(observed)
          as.numeric(names(counts)[which.max(counts)])
        }
      )
    })
    for (i in which(rowSums(missing_before) > 0L)) {
      dosage[i, missing_before[i, ]] <- replacement[[i]]
    }
  }
  dosage[dosage < 0] <- 0
  dosage[dosage > ploidy] <- ploidy
  if (anyNA(dosage)) {
    stop("Native dosage imputation left missing values.", call. = FALSE)
  }
  gp_validate_alt_dosage(dosage, ploidy, hard_calls = FALSE)
  list(
    dosage = gp_attach_ploidy(dosage, ploidy),
    imputed_cells = as.integer(sum(missing_before)),
    method = method
  )
}

#' Native ploidy-aware genotype dosage imputation
#'
#' Imputes missing values in an individual-by-marker ALT-dosage matrix without
#' an external imputation program. This is a dosage-level imputer; it does not
#' phase haplotypes or add variants from a reference panel.
#'
#' @param geno_data Numeric matrix or data frame with individuals in rows and
#'   markers in columns. Values are ALT-allele dosage in `[0, ploidy]`.
#' @param ploidy One positive integer or `"auto"`. Unannotated polyploid
#'   matrices require an explicit value.
#' @param method One of `"mean"`, `"median"`, `"mode"`, or `"knn"`.
#' @param k Number of neighbours for `method = "knn"`.
#' @param recode_format Output coding: ALT dosage, dosage centered on
#'   `ploidy / 2`, or allele frequency. Legacy diploid aliases are accepted.
#'
#' @return A list containing `snps_matrix`, `imputed_cells`, `method`,
#'   `ploidy`, and `coding_contract`.
#' @export
impute_genotypes_native <- function(geno_data,
                                    ploidy = "auto",
                                    method = c("mean", "median", "mode", "knn"),
                                    k = 5L,
                                    recode_format = "alt_dosage") {
  if (is.data.frame(geno_data)) geno_data <- as.matrix(geno_data)
  if (!is.matrix(geno_data) || !is.numeric(geno_data)) {
    stop("geno_data must be a numeric individual-by-marker matrix or data frame.", call. = FALSE)
  }
  resolved_ploidy <- gp_resolve_matrix_ploidy(geno_data, ploidy)
  gp_validate_alt_dosage(geno_data, resolved_ploidy, hard_calls = FALSE)
  original_dimnames <- dimnames(geno_data)
  imputed <- gp_impute_alt_dosage(
    dosage = t(geno_data), ploidy = resolved_ploidy,
    method = method, k = k
  )
  out <- t(imputed$dosage)
  dimnames(out) <- original_dimnames
  out <- gp_recode_alt_dosage(out, resolved_ploidy, recode_format)
  result <- list(
    snps_matrix = out,
    imputed_cells = imputed$imputed_cells,
    method = imputed$method,
    ploidy = resolved_ploidy,
    coding_contract = gp_normalize_recode_format(recode_format)
  )
  class(result) <- c("predictpror_genotype_imputation", "list")
  result
}

gp_extract_gt_ploidy <- function(gt) {
  gt <- trimws(as.character(gt))
  gt[is.na(gt) | !nzchar(gt)] <- "."
  # Vectorised over the distinct calls: allele count = separators + 1.
  # (A per-cell strsplit took minutes on large VCFs.)
  values <- unique(gt)
  counts <- nchar(gsub("[^/|]", "", values)) + 1L
  counts[values == "."] <- NA_integer_
  as.integer(counts[match(gt, values)])
}

gp_resolve_gt_ploidy <- function(gt, ploidy = "auto", context = "VCF") {
  requested <- gp_validate_ploidy(ploidy)
  observed <- gp_extract_gt_ploidy(gt)
  observed_unique <- sort(unique(observed[!is.na(observed)]))
  if (identical(requested, "auto")) {
    if (!length(observed_unique)) return(2L)
    if (length(observed_unique) != 1L) {
      stop(
        context, " contains mixed GT ploidies (",
        paste(observed_unique, collapse = ", "),
        "). PredictProR currently requires one uniform crop ploidy per analysis.",
        call. = FALSE
      )
    }
    return(as.integer(observed_unique[[1L]]))
  }
  incompatible <- !is.na(observed) & observed != requested
  if (any(incompatible)) {
    stop(
      context, " GT calls do not match ploidy = ", requested,
      "; first incompatible call has ", observed[which(incompatible)[[1L]]],
      " allele fields.", call. = FALSE
    )
  }
  requested
}
