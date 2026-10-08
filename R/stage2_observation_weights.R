# Stage 2 observation-weight contract ---------------------------------------

gp_stage2_weight_column <- function() {
  ".PredictProR_stage2_precision_weight"
}

gp_stage2_weight_values <- function(weights,
                                    pheno_data = NULL,
                                    context = "Stage 2 observation weights") {
  if (is.null(weights)) {
    return(NULL)
  }

  source <- "numeric_vector"
  value <- weights
  if (is.character(weights) && length(weights) == 1L) {
    if (is.null(pheno_data) || !weights %in% names(pheno_data)) {
      stop(
        context, ": `weights = \"", weights,
        "\"` does not identify a column in `pheno_data`.",
        call. = FALSE
      )
    }
    source <- paste0("pheno_data$", weights)
    value <- pheno_data[[weights]]
  } else if (is.data.frame(weights) || is.matrix(weights)) {
    if (ncol(weights) != 1L) {
      stop(context, ": a data-frame or matrix weight input must have exactly one column.", call. = FALSE)
    }
    source_name <- colnames(weights)[1L]
    if (is.null(source_name) || is.na(source_name) || !nzchar(source_name)) {
      source_name <- "single_column_input"
    }
    source <- source_name
    value <- weights[, 1L]
  }

  if (is.factor(value)) {
    value <- as.character(value)
  }
  if (!is.atomic(value) || is.list(value)) {
    stop(context, ": weights must be a numeric vector, a phenotype-column name, or a one-column table.", call. = FALSE)
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) != length(value) ||
      any(is.na(numeric_value) & !is.na(value))) {
    stop(context, ": weights must be numeric or losslessly coercible to numeric.", call. = FALSE)
  }
  if (!is.null(pheno_data) && length(numeric_value) != nrow(pheno_data)) {
    stop(
      context, ": weight length (", length(numeric_value),
      ") does not match `nrow(pheno_data)` (", nrow(pheno_data), ").",
      call. = FALSE
    )
  }

  list(values = numeric_value, source = source)
}

gp_stage2_test_mask <- function(pheno_data,
                                response = NULL,
                                gen_name = NULL,
                                test_set = NULL,
                                test_mask = NULL) {
  n <- nrow(pheno_data)
  if (!is.null(test_mask)) {
    test_mask <- as.logical(test_mask)
    if (length(test_mask) != n) {
      stop("Stage 2 weight test mask length does not match `pheno_data`.", call. = FALSE)
    }
    test_mask[is.na(test_mask)] <- FALSE
    return(test_mask)
  }

  out <- rep(FALSE, n)
  if (!is.null(test_set) && !is.null(gen_name) && gen_name %in% names(pheno_data)) {
    if (is.data.frame(test_set) || is.matrix(test_set)) {
      test_set <- test_set[, 1L]
    }
    out <- as.character(pheno_data[[gen_name]]) %in% as.character(test_set)
  }
  if (!is.null(response)) {
    response <- intersect(as.character(response), names(pheno_data))
    if (length(response)) {
      response_missing <- Reduce(
        `&`,
        lapply(response, function(one) {
          value <- pheno_data[[one]]
          is.na(value) | (is.numeric(value) & !is.finite(value))
        })
      )
      out <- out | response_missing
    }
  }
  out
}

gp_resolve_stage2_precision_weights <- function(weights,
                                                pheno_data,
                                                response = NULL,
                                                gen_name = NULL,
                                                test_set = NULL,
                                                test_mask = NULL,
                                                context = "Stage 2 observation weights") {
  if (is.null(weights)) {
    return(list(
      supplied = FALSE,
      precision = NULL,
      source = NULL,
      test_mask = rep(FALSE, nrow(pheno_data))
    ))
  }
  resolved <- gp_stage2_weight_values(weights, pheno_data = pheno_data, context = context)
  prediction_rows <- gp_stage2_test_mask(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    test_set = test_set,
    test_mask = test_mask
  )
  invalid <- !is.finite(resolved$values) | resolved$values <= 0
  invalid_training <- invalid & !prediction_rows
  if (any(invalid_training)) {
    bad <- which(invalid_training)
    stop(
      context, ": every analyzed/training row must have a finite weight greater than zero. Invalid row(s): ",
      paste(utils::head(bad, 10L), collapse = ", "),
      if (length(bad) > 10L) " ..." else "",
      ".",
      call. = FALSE
    )
  }

  # Prediction-only rows do not enter the likelihood. A neutral positive value
  # keeps ASReml/BGLR/GP model frames valid without inventing training weights.
  precision <- resolved$values
  precision[invalid & prediction_rows] <- 1
  list(
    supplied = TRUE,
    precision = precision,
    source = resolved$source,
    test_mask = prediction_rows
  )
}

gp_attach_stage2_precision_weights <- function(pheno_data,
                                               weights,
                                               response = NULL,
                                               gen_name = NULL,
                                               test_set = NULL,
                                               context = "model_execute Stage 2 observation weights") {
  if (is.null(weights)) {
    return(list(pheno_data = pheno_data, weights = NULL, source = NULL))
  }
  resolved <- gp_resolve_stage2_precision_weights(
    weights = weights,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    test_set = test_set,
    context = context
  )
  weight_col <- gp_stage2_weight_column()
  if (weight_col %in% names(pheno_data)) {
    stop(
      "`pheno_data` contains the reserved internal column `", weight_col,
      "`; rename that column before fitting.",
      call. = FALSE
    )
  }
  pheno_data[[weight_col]] <- resolved$precision
  list(pheno_data = pheno_data, weights = weight_col, source = resolved$source)
}

gp_bglr_weights_from_precision <- function(weights,
                                           pheno_data = NULL,
                                           response = NULL,
                                           y = NULL,
                                           gen_name = NULL,
                                           test_set = NULL,
                                           test_mask = NULL,
                                           context = "BGLR Stage 2 observation weights") {
  if (is.null(weights)) {
    return(NULL)
  }
  if (is.null(pheno_data)) {
    if (is.null(y)) {
      stop(context, ": `pheno_data` or `y` is required to align weights.", call. = FALSE)
    }
    pheno_data <- data.frame(.stage2_response = y)
    response <- ".stage2_response"
  } else if (!is.null(y)) {
    if (length(y) != nrow(pheno_data)) {
      stop(context, ": response length does not match `pheno_data`.", call. = FALSE)
    }
    response_col <- ".PredictProR_stage2_response_for_weight_alignment"
    pheno_data[[response_col]] <- y
    response <- response_col
  }
  resolved <- gp_resolve_stage2_precision_weights(
    weights = weights,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    test_set = test_set,
    test_mask = test_mask,
    context = context
  )
  # BGLR defines Var(e_i) proportional to 1 / bglr_weight_i^2. Public
  # PredictProR weights are precisions, so sqrt(precision) gives 1/precision.
  sqrt(resolved$precision)
}

gp_asreml_weighted_data_clause <- function(data_expr,
                                           weight_column = NULL,
                                           response_family = "gaussian") {
  if (is.null(weight_column)) {
    return(paste0('na.action=list(x="include",y="include"),data=', data_expr, ')'))
  }
  if (!identical(tolower(as.character(response_family)[1L]), "gaussian")) {
    stop("Stage 2 precision weights are currently supported for Gaussian ASReml fits only.", call. = FALSE)
  }
  paste0(
    'na.action=list(x="include",y="include"), weights = ', weight_column,
    ', family = asr_gaussian(dispersion = 1), data=', data_expr, ')'
  )
}

gp_validate_stage2_weight_route <- function(ctx) {
  if (is.null(ctx$weights)) {
    return(invisible(TRUE))
  }
  requested <- unique(stats::na.omit(as.character(c(ctx$GS_model, ctx$GS_model_cv))))
  rkhs_heterogeneous <- if (!is.null(ctx$bayes_kernel_heter_resid)) {
    isTRUE(ctx$bayes_kernel_heter_resid)
  } else {
    isTRUE(ctx$heter_resid)
  }
  if ("RKHS" %in% requested && isTRUE(rkhs_heterogeneous)) {
    stop(
      paste(
        "RKHS cannot combine Stage 2 observation weights with heterogeneous-environment mode",
        "because BGLR::Multitrait has no observation-weight argument.",
        "Remove `weights` to estimate environment-specific genetic and residual variances,",
        "or set `bayes_kernel_heter_resid = FALSE` for weighted environment-aware RKHS",
        "with common genetic and residual variance parameters."
      ),
      call. = FALSE
    )
  }
  unsupported_flags <- c(
    multi_trait_gp = isTRUE(ctx$multi_trait_gp),
    multi_trait_bayes = isTRUE(ctx$multi_trait_bayes),
    multi_trait_asreml = isTRUE(ctx$multi_trait_asreml),
    multi_trait_ml = isTRUE(ctx$multi_trait_ml),
    multi_trait_dl = isTRUE(ctx$multi_trait_dl),
    met_ml_dl = isTRUE(ctx$met_ml_dl),
    hybrid_gp = isTRUE(ctx$hybrid_gp),
    hybrid_bayes = isTRUE(ctx$hybrid_bayes),
    hybrid_asreml = isTRUE(ctx$hybrid_asreml),
    hybrid_ml = isTRUE(ctx$hybrid_ml),
    hybrid_dl = isTRUE(ctx$hybrid_dl)
  )
  if (any(unsupported_flags)) {
    stop(
      "Stage 2 observation weights are not implemented for route(s): ",
      paste(names(unsupported_flags)[unsupported_flags], collapse = ", "),
      ". Use a supported single-trait GP, ASReml GBLUP, or univariate BGLR route; weights will not be silently ignored.",
      call. = FALSE
    )
  }

  allowed <- unique(stats::na.omit(as.character(c(
    ctx$gp_valid_models,
    ctx$bayes_valid_models,
    ctx$bayes_gblup_valid_models,
    ctx$asreml_model
  ))))
  unsupported_models <- setdiff(requested, allowed)
  if (length(unsupported_models)) {
    stop(
      "Stage 2 observation weights are not supported by model(s): ",
      paste(unsupported_models, collapse = ", "),
      ". Supported weighted families are GP, ASReml GBLUP, and univariate BGLR; weights will not be silently ignored.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
