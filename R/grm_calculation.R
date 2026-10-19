gp_grm_dense_cpp_available <- function() {
  !is.null(tryCatch(
    getNativeSymbolInfo("predictpror_grm_dense", PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

gp_grm_backend <- function(backend = NULL, method = NULL) {
  if (is.null(backend)) {
    backend <- Sys.getenv("PREDICTPRO_GRM_BACKEND", Sys.getenv("PREDICTPRO_KERNEL_BACKEND", "auto"))
  }
  backend <- as.character(backend)
  if (!identical(length(backend), 1L) || !(backend %in% c("auto", "cpp", "r"))) {
    stop("backend must be one of: auto, cpp, r.", call. = FALSE)
  }
  if (identical(backend, "auto") && gp_grm_dense_cpp_available()) {
    method <- if (is.null(method)) "" else as.character(method)
    # The native backend builds the covariates in one pass and forms the
    # product with BLAS dsyrk: as fast as R for VanRaden, 0.3-1.1 s faster for
    # the others at 1,000 x 10,000, without R's temporaries (0.20.211).
    native_auto_methods <- c(
      "VanRaden",
      "Weighted_VanRaden",
      "Epistasis",
      "Yang",
      "Dominance",
      "Dominance_Vitezica",
      "Dominance_Su",
      "Dominance_Heterozygosity"
    )
    if (method %in% native_auto_methods) {
      return("cpp")
    }
  }
  if (identical(backend, "auto")) {
    return("r")
  }
  if (identical(backend, "cpp") && !gp_grm_dense_cpp_available()) {
    stop("PredictProR native dense GRM backend is not available.", call. = FALSE)
  }
  backend
}

#' Calculate Genomic Relationship Matrix (GRM)
#'
#' This function computes dense genomic relationship matrices from hard-call
#' genotype data using additive, dominance, Yang, or epistasis methods.
#' Additive VanRaden, weighted VanRaden, and additive-by-additive epistasis are
#' available for uniform diploid or autopolyploid dosage. Yang and the current
#' dominance parameterizations remain explicitly diploid.
#' By default, PredictProR selects the fastest dense route for the selected
#' method unless a backend is forced.
#'
#' @param geno_clean A matrix of genotype data where rows represent individuals
#'   and columns represent SNPs. Genotype data must be ALT dosage from zero
#'   through `ploidy`.
#' @param ploidy One positive integer or `"auto"`. Automatic resolution uses
#'   matrix metadata; an unannotated matrix defaults to diploid only when all
#'   observed dosages are at most two.
#' @param weight Optional SNP weights for `Weighted_VanRaden`. A one-column
#'   marker-named matrix or a square diagonal matrix is accepted. Weights must
#'   be finite and non-negative; the implemented scaling is
#'   `Z \%*\% diag(weight) \%*\% t(Z) / sum(ploidy * p * (1 - p))`.
#' @param method A character scalar specifying the GRM method. Options are
#'   `"VanRaden"`, `"Weighted_VanRaden"`, `"Yang"`, `"Epistasis"`,
#'   `"Dominance_Vitezica"`, and `"Dominance_Su"`. `"Dominance"` is an alias
#'   for `"Dominance_Vitezica"`; `"Dominance_Heterozygosity"` is an alias for
#'   `"Dominance_Su"`.
#' @param backend Backend selector. One of `"auto"`, `"cpp"`, or `"r"`. `NULL`
#'   uses `PREDICTPRO_GRM_BACKEND`, `PREDICTPRO_KERNEL_BACKEND`, or `"auto"`.
#'   The automatic route uses the native backend, which forms each product
#'   with BLAS dsyrk for every ploidy the method allows.
#'
#' @return A matrix representing the genomic relationship matrix calculated
#'   based on the specified method.
#'
#' @examples
#' geno_clean <- matrix(c(0, 1, 2, 1, 0, 1, 2, 2, 1), nrow = 3, byrow = TRUE)
#' colnames(geno_clean) <- c("SNP1", "SNP2", "SNP3")
#' rownames(geno_clean) <- c("Ind1", "Ind2", "Ind3")
#' grm_vanraden <- PredictProR:::grm_calculation(geno_clean, method = "VanRaden")
grm_calculation <- function(
    geno_clean = NULL,
    weight = NULL,
    method = NULL,
    backend = NULL,
    ploidy = "auto"
) {

  msg <- ""

  if (!is.matrix(geno_clean)) {
    stop(paste(msg, "snp/marker data must be a matrix."), call. = FALSE)
  }
  if (!is.numeric(geno_clean)) {
    stop(paste(msg, "SNP data must contain numeric ALT-allele dosage."), call. = FALSE)
  }
  ploidy <- gp_resolve_matrix_ploidy(geno_clean, ploidy)
  gp_validate_alt_dosage(geno_clean, ploidy, hard_calls = FALSE, name = "GRM genotype dosage")
  if (!is.null(weight) && !is.matrix(weight)) {
    stop(paste(msg, "weight must be a matrix if provided."), call. = FALSE)
  }

  gmatrix_method_available <- c(
    "VanRaden",
    "Weighted_VanRaden",
    "Yang",
    "Epistasis",
    "Dominance",
    "Dominance_Vitezica",
    "Dominance_Su",
    "Dominance_Heterozygosity"
  )
  method <- if (is.null(method)) NULL else as.character(method)
  if (is.null(method) || !length(method) || anyNA(method) || !all(method %in% gmatrix_method_available)) {
    stop(
      paste(msg, "Invalid genomic relationship method. Choose from: ",
            paste(gmatrix_method_available, collapse = ", ")),
      call. = FALSE
    )
  }
  method <- unique(vapply(method, function(x) switch(
    x,
    Dominance = "Dominance_Vitezica",
    Dominance_Heterozygosity = "Dominance_Su",
    x
  ), character(1)))

  diploid_only <- c("Yang", "Dominance_Vitezica", "Dominance_Su")
  incompatible <- intersect(method, diploid_only)
  if (ploidy != 2L && length(incompatible)) {
    stop(
      paste(msg, paste(incompatible, collapse = ", "),
            "is currently a diploid-only GRM parameterization. For ploidy", ploidy,
            "use VanRaden, Weighted_VanRaden, Epistasis, or a generic marker kernel."),
      call. = FALSE
    )
  }


  storage.mode(geno_clean) <- "double"
  pre_marker_names <- colnames(geno_clean)
  marker_has_missing <- colSums(is.na(geno_clean)) > 0
  # a marker varies when any line differs from the first line (vectorised;
  # replaces a per-marker unique() scan, same result)
  marker_varies <- colSums(geno_clean != rep(geno_clean[1L, ], each = nrow(geno_clean)), na.rm = TRUE) > 0
  keep_markers <- !marker_has_missing & marker_varies
  if (!all(keep_markers)) {
    geno_clean <- geno_clean[, keep_markers, drop = FALSE]
  }
  if (!ncol(geno_clean)) {
    stop(paste(msg, "No informative SNP markers remain after genotype cleanup."), call. = FALSE)
  }
  if (anyNA(geno_clean)) {
    stop(paste(msg, "Missing SNP data remain after genotype cleanup."), call. = FALSE)
  }
  if (!is.null(weight) && is.null(colnames(geno_clean))) {
    stop(paste(msg, "SNP marker names are required when weight is provided."), call. = FALSE)
  }

  gp_grm_weight_vector <- function(weight, marker_names) {
    if (is.null(weight)) {
      return(NULL)
    }
    weight <- as.matrix(weight)
    if (nrow(weight) == length(marker_names) && ncol(weight) == length(marker_names)) {
      has_row_names <- !is.null(rownames(weight))
      has_col_names <- !is.null(colnames(weight))
      if (has_row_names && has_col_names) {
        if (!all(marker_names %in% rownames(weight)) || !all(marker_names %in% colnames(weight))) {
          stop(paste(msg, "square weight matrix names must match SNP marker names."), call. = FALSE)
        }
        weight <- weight[marker_names, marker_names, drop = FALSE]
      } else if (has_row_names && !identical(rownames(weight), marker_names)) {
        stop(paste(msg, "square weight matrix col names are required when row order differs from SNP marker order."), call. = FALSE)
      } else if (has_col_names && !identical(colnames(weight), marker_names)) {
        stop(paste(msg, "square weight matrix row names are required when column order differs from SNP marker order."), call. = FALSE)
      }
      off_diag <- weight
      diag(off_diag) <- 0
      if (anyNA(off_diag) || any(abs(off_diag) > sqrt(.Machine$double.eps))) {
        stop(paste(msg, "square weight matrix must be diagonal for Weighted_VanRaden."), call. = FALSE)
      }
      return(as.numeric(diag(weight)))
    }
    if (is.null(rownames(weight))) {
      stop(paste(msg, "weight row names must match SNP marker names."), call. = FALSE)
    }
    weight <- weight[match(marker_names, rownames(weight)), , drop = FALSE]
    if (!identical(marker_names, rownames(weight))) {
      stop(paste(msg, "SNP order in weight does not match geno_clean after marker cleanup."), call. = FALSE)
    }
    as.numeric(weight[, 1])
  }

  weights_vec <- gp_grm_weight_vector(weight, colnames(geno_clean))
  if (any(method %in% "Weighted_VanRaden") && is.null(weights_vec)) {
    stop(paste(msg, "weight is required for Weighted_VanRaden."), call. = FALSE)
  }
  if (!is.null(weights_vec) && (length(weights_vec) != ncol(geno_clean) || anyNA(weights_vec) || any(!is.finite(weights_vec)))) {
    stop(paste(msg, "weight must provide one finite value per retained SNP marker."), call. = FALSE)
  }
  if (!is.null(weights_vec) && (any(weights_vec < 0) || !any(weights_vec > 0))) {
    stop(paste(msg, "Weighted_VanRaden weights must be non-negative with at least one positive value."), call. = FALSE)
  }

  freq <- colMeans(geno_clean) / ploidy
  denom <- sum(ploidy * freq * (1 - freq))
  if (!is.finite(denom) || denom <= 0) {
    stop(paste(msg, "GRM denominator must be positive after genotype cleanup."), call. = FALSE)
  }

  # Computed only for the methods that use them: the additive product was
  # formed for every method (dominance, weighted and Yang then formed their
  # own; the native backend used neither): 8 s wasted at 1,000 x 10,000.
  centered_cache <- NULL
  base_cache <- NULL
  centered_markers <- function() {
    if (is.null(centered_cache)) centered_cache <<- scale(geno_clean, center = TRUE, scale = FALSE)
    centered_cache
  }
  base_grm <- function() {
    if (is.null(base_cache)) base_cache <<- tcrossprod(centered_markers()) / denom
    base_cache
  }
  q <- 1 - freq

  dominance_vitezica_grm <- function() {
    dominance_covariates <- matrix(0, nrow = nrow(geno_clean), ncol = ncol(geno_clean))
    for (k in seq_len(ncol(geno_clean))) {
      marker <- geno_clean[, k]
      dominance_covariates[marker == 0, k] <- -2 * freq[[k]]^2
      dominance_covariates[marker == 1, k] <- 2 * freq[[k]] * q[[k]]
      dominance_covariates[marker == 2, k] <- -2 * q[[k]]^2
    }
    dominance_denom <- sum((2 * freq * q)^2)
    if (!is.finite(dominance_denom) || dominance_denom <= 0) {
      stop(paste(msg, "Dominance_Vitezica denominator must be positive after genotype cleanup."), call. = FALSE)
    }
    tcrossprod(dominance_covariates) / dominance_denom
  }

  dominance_su_grm <- function() {
    heterozygosity_covariates <- matrix(0, nrow = nrow(geno_clean), ncol = ncol(geno_clean))
    for (k in seq_len(ncol(geno_clean))) {
      marker <- geno_clean[, k]
      pq2 <- 2 * freq[[k]] * q[[k]]
      heterozygosity_covariates[marker == 0, k] <- -pq2
      heterozygosity_covariates[marker == 1, k] <- 1 - pq2
      heterozygosity_covariates[marker == 2, k] <- -pq2
    }
    dominance_denom <- sum(2 * freq * q * (1 - 2 * freq * q))
    if (!is.finite(dominance_denom) || dominance_denom <= 0) {
      stop(paste(msg, "Dominance_Su denominator must be positive after genotype cleanup."), call. = FALSE)
    }
    tcrossprod(heterozygosity_covariates) / dominance_denom
  }

  calculate_one_grm <- function(method_one) {
    selected_backend <- gp_grm_backend(backend, method_one)
    if (identical(selected_backend, "cpp")) {
      Ga <- .Call(
        "predictpror_grm_dense",
        geno_clean,
        as.character(method_one),
        weights_vec,
        list(ploidy = ploidy),
        PACKAGE = "PredictProR"
      )
      rownames(Ga) <- rownames(geno_clean)
      colnames(Ga) <- rownames(geno_clean)
      attr(Ga, "ploidy") <- ploidy
      return(Ga)
    }

    # one-argument tcrossprod() is a BLAS dsyrk; with two arguments it is a
    # dgemm over both triangles (12.4 s vs 5-8 s at 1,000 x 10,000)
    Ga <- switch(
      method_one,
      VanRaden = base_grm(),
      Weighted_VanRaden = {
        tcrossprod(sweep(centered_markers(), 2, sqrt(weights_vec), `*`)) / denom
      },
      Epistasis = base_grm() * base_grm(),
      Yang = {
        n_marker <- ncol(geno_clean)
        n_individuals <- nrow(geno_clean)
        inv_var <- 1 / (2 * freq * (1 - freq))
        locus <- (1 / n_marker) * tcrossprod(sweep(centered_markers(), 2, sqrt(inv_var), `*`))
        locus[lower.tri(locus, diag = TRUE)] <- 0
        locus <- locus + t(locus)
        multiplier <- geno_clean^2 -
          t(t(geno_clean) * (1 + 2 * freq)) +
          matrix(rep(2 * freq^2, each = n_individuals), ncol = n_marker)
        diag(locus) <- 1 + (1 / n_marker) * colSums(t(multiplier) * inv_var)
        locus
      },
      Dominance_Vitezica = dominance_vitezica_grm(),
      Dominance_Su = dominance_su_grm(),
      stop(paste(msg, "Select method to calculate genomic relationship matrix"), call. = FALSE)
    )

    Ga <- as.matrix(Ga)
    rownames(Ga) <- rownames(geno_clean)
    colnames(Ga) <- rownames(geno_clean)
    attr(Ga, "ploidy") <- ploidy
    Ga
  }

  if (length(method) > 1L) {
    return(gp_named_method_list(lapply(method, calculate_one_grm), method))
  }

  calculate_one_grm(method)
}
