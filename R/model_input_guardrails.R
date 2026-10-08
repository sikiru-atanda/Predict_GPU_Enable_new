gp_model_input_guardrail_stop <- function(label, detail) {
  msg <- ""
  stop(
    paste(
      msg,
      "PredictProR input guardrail failed for",
      label,
      "\n",
      detail
    ),
    call. = FALSE
  )
}

gp_kernel_container <- function(x) {
  is.list(x) && !is.data.frame(x) && !is.matrix(x)
}

gp_coerce_named_numeric_kernel <- function(x, label = "kernel") {
  if (is.null(x)) {
    return(NULL)
  }

  source_kernel <- attr(x, "predictpror_source_kernel", exact = TRUE)
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    if (length(dim(x)) == 2L) {
      x <- as.matrix(x)
    } else {
      gp_model_input_guardrail_stop(
        label,
        "Provide a square matrix-like object with row and column names."
      )
    }
  }

  if (length(dim(x)) != 2L || nrow(x) != ncol(x)) {
    gp_model_input_guardrail_stop(
      label,
      "Kernels and relationship matrices must be square."
    )
  }

  rn <- rownames(x)
  cn <- colnames(x)
  if (is.null(rn) || anyNA(rn) || any(!nzchar(as.character(rn)))) {
    gp_model_input_guardrail_stop(
      label,
      "Row names must contain the individual IDs used by gen_name."
    )
  }
  rn <- as.character(rn)
  if (is.null(cn)) {
    cn <- rn
    colnames(x) <- cn
  }
  cn <- as.character(cn)
  if (anyNA(cn) || any(!nzchar(cn))) {
    gp_model_input_guardrail_stop(
      label,
      "Column names must contain the individual IDs used by gen_name."
    )
  }
  if (anyDuplicated(rn) || anyDuplicated(cn)) {
    gp_model_input_guardrail_stop(
      label,
      "Row and column names must be unique individual IDs."
    )
  }
  if (!setequal(rn, cn)) {
    gp_model_input_guardrail_stop(
      label,
      "Row and column names must describe the same set of individuals."
    )
  }
  if (!identical(rn, cn)) {
    x <- x[, rn, drop = FALSE]
    cn <- colnames(x)
  }

  na_before <- is.na(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  if (any(is.na(x) & !na_before)) {
    gp_model_input_guardrail_stop(
      label,
      "All kernel or relationship matrix entries must be numeric."
    )
  }
  if (anyNA(x)) {
    gp_model_input_guardrail_stop(
      label,
      "Kernels and relationship matrices cannot contain missing values."
    )
  }
  if (any(!is.finite(x))) {
    gp_model_input_guardrail_stop(
      label,
      "Kernels and relationship matrices must contain finite numeric values."
    )
  }

  rownames(x) <- rn
  colnames(x) <- rn
  if (!is.null(source_kernel)) {
    attr(x, "predictpror_source_kernel") <- source_kernel
  }
  x
}

gp_normalize_kernel_guardrail <- function(x, label = "kernel") {
  if (is.null(x)) {
    return(NULL)
  }
  if (gp_kernel_container(x)) {
    x_names <- names(x)
    out <- lapply(seq_along(x), function(i) {
      nm <- if (!is.null(x_names) && length(x_names) >= i) x_names[[i]] else ""
      if (is.null(nm) || !nzchar(nm)) {
        nm <- paste0("kernel_", i)
      }
      gp_normalize_kernel_guardrail(x[[i]], paste(label, nm, sep = "$"))
    })
    names(out) <- x_names
    return(out)
  }
  gp_coerce_named_numeric_kernel(x, label = label)
}

gp_normalize_env_similarity_guardrail <- function(x,
                                                  env_ids = NULL,
                                                  label = "env_similarity") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    if (length(dim(x)) == 2L) {
      x <- as.matrix(x)
    } else {
      gp_model_input_guardrail_stop(
        label,
        "Provide a square matrix-like object for the environment similarity kernel."
      )
    }
  }
  if (length(dim(x)) != 2L || nrow(x) != ncol(x)) {
    gp_model_input_guardrail_stop(
      label,
      "Environment similarity kernels must be square."
    )
  }

  ids <- as.character(env_ids %||% character())
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (is.null(rownames(x)) && length(ids) == nrow(x)) {
    rownames(x) <- ids
  }
  if (is.null(colnames(x)) && length(ids) == ncol(x)) {
    colnames(x) <- ids
  }

  gp_coerce_named_numeric_kernel(x, label = label)
}

gp_coerce_feature_matrix_guardrail <- function(x,
                                               label = "feature matrix",
                                               id_col = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    id_col <- as.character(id_col %||% "")
    if (nzchar(id_col) && id_col %in% names(x)) {
      ids <- as.character(x[[id_col]])
      if (anyNA(ids) || any(!nzchar(ids))) {
        gp_model_input_guardrail_stop(
          label,
          paste0("The ", id_col, " column must contain non-missing individual IDs.")
        )
      }
      if (anyDuplicated(ids)) {
        gp_model_input_guardrail_stop(
          label,
          paste0("The ", id_col, " column must contain unique individual IDs.")
        )
      }
      existing_rn <- rownames(x)
      default_rn <- is.null(existing_rn) ||
        identical(as.character(existing_rn), as.character(seq_len(nrow(x))))
      if (!isTRUE(default_rn) && !identical(as.character(existing_rn), ids)) {
        gp_model_input_guardrail_stop(
          label,
          paste0("Row names and the ", id_col, " column describe different individual IDs.")
        )
      }
      x[[id_col]] <- NULL
      rownames(x) <- ids
    }
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    if (length(dim(x)) == 2L) {
      x <- as.matrix(x)
    } else {
      gp_model_input_guardrail_stop(
        label,
        "Provide a matrix-like feature table with row names keyed by gen_name."
      )
    }
  }

  rn <- rownames(x)
  if (is.null(rn) || anyNA(rn) || any(!nzchar(as.character(rn)))) {
    gp_model_input_guardrail_stop(
      label,
      "Feature matrices must have row names keyed by gen_name."
    )
  }
  if (anyDuplicated(as.character(rn))) {
    gp_model_input_guardrail_stop(
      label,
      "Feature matrix row names must be unique individual IDs."
    )
  }
  cn <- colnames(x)
  if (is.null(cn) || anyNA(cn) || any(!nzchar(as.character(cn)))) {
    gp_model_input_guardrail_stop(
      label,
      "Feature matrices must have non-missing predictor names."
    )
  }
  if (anyDuplicated(as.character(cn))) {
    gp_model_input_guardrail_stop(
      label,
      "Feature matrices must not contain duplicated predictor names."
    )
  }

  na_before <- is.na(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  if (any(is.na(x) & !na_before)) {
    gp_model_input_guardrail_stop(
      label,
      "Feature matrix entries must be numeric."
    )
  }
  rownames(x) <- as.character(rn)
  colnames(x) <- as.character(cn)
  x
}

gp_normalize_feature_bank_guardrail <- function(x,
                                                label = "feature_bank",
                                                id_col = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  if (gp_kernel_container(x)) {
    x_names <- names(x)
    out <- lapply(seq_along(x), function(i) {
      nm <- if (!is.null(x_names) && length(x_names) >= i) x_names[[i]] else ""
      if (is.null(nm) || !nzchar(nm)) {
        nm <- paste0("feature_", i)
      }
      gp_normalize_feature_bank_guardrail(
        x[[i]],
        paste(label, nm, sep = "$"),
        id_col = id_col
      )
    })
    names(out) <- x_names
    return(out)
  }
  gp_coerce_feature_matrix_guardrail(x, label = label, id_col = id_col)
}

gp_asreml_order_required <- function(ctx) {
  engine <- tolower(as.character(ctx$engine %||% ""))
  shared_requires_asreml <- ctx$predictpror_shared_requires_asreml %||%
    ctx$.predictpror_shared_requires_asreml %||%
    FALSE
  isTRUE(ctx$multi_trait_asreml) ||
    isTRUE(ctx$hybrid_asreml) ||
    isTRUE(shared_requires_asreml) ||
    identical(engine, "asreml")
}

gp_reject_obsolete_asreml_structure <- function(var_cov_str) {
  if (is.null(var_cov_str)) {
    return(invisible(var_cov_str))
  }
  values <- tolower(trimws(as.character(var_cov_str)))
  if (any(!is.na(values) & values == "corgv")) {
    stop(
      paste0(
        "`var_cov_str = \"corgv\"` is not supported by current ASReml-R. ",
        "Use `var_cov_str = \"corgh\"`."
      ),
      call. = FALSE
    )
  }
  invisible(var_cov_str)
}

gp_order_pheno_for_model_guardrail <- function(pheno_data,
                                               gen_name,
                                               heter_groups = NULL,
                                               order_for_asreml = FALSE) {
  if (is.null(pheno_data) || !is.data.frame(pheno_data)) {
    return(pheno_data)
  }
  if (is.null(gen_name) || !gen_name %in% names(pheno_data)) {
    return(pheno_data)
  }
  if (!isTRUE(order_for_asreml)) {
    return(pheno_data)
  }

  ph <- pheno_data
  cleared <- attr(ph, "cleared", exact = TRUE)
  order_cols <- gen_name
  heter_valid <- !is.null(heter_groups) &&
    length(heter_groups) == 1L &&
    !is.na(heter_groups) &&
    nzchar(heter_groups) &&
    heter_groups %in% names(ph)
  if (isTRUE(heter_valid)) {
    order_cols <- c(heter_groups, gen_name)
  }

  ord <- do.call(
    order,
    c(lapply(order_cols, function(nm) as.character(ph[[nm]])), list(na.last = TRUE))
  )
  ph <- ph[ord, , drop = FALSE]
  rownames(ph) <- NULL
  if (!is.null(cleared)) {
    attr(ph, "cleared") <- cleared
  }
  ph
}

#' Warn when a single-kernel GP model is used with a multi-environment panel
#'
#' Internal guardrail. KRR and GP use a direct single-environment covariance
#' structure, although each can consume an additive bank of genomic/omics
#' kernels. The legacy LowRankGP single-trait route is a KRR alias and is
#' rejected at task selection. On a multi-environment panel these direct models
#' fit a near-zero main genetic variance and push the genetic signal into the
#' residual, severely under-estimating per-environment heritability (validated:
#' per-env h2 ~0.08-0.16 vs ~0.5-0.7 for GP_FA / ASReml corgh). This emits a
#' steering warning toward the factor-analytic MET model. GP_FA is exempt.
#'
#' @param GS_model,GS_model_cv Requested model name(s) (canonical) for true
#'   prediction and cross-validation.
#' @param cross_validation Whether the run is cross-validation (selects the
#'   active model name).
#' @param is_met Whether the phenotype panel is multi-environment.
#' @return `NULL`, invisibly (called for its warning side effect).
#' @keywords internal
#' @noRd
gp_warn_single_kernel_gp_met <- function(GS_model, GS_model_cv, cross_validation, is_met) {
  if (!isTRUE(is_met)) return(invisible(NULL))
  single_kernel_gp <- c("KRR", "GP", "LowRankGP")
  models <- if (isTRUE(cross_validation)) GS_model_cv else GS_model
  hit <- intersect(as.character(models), single_kernel_gp)
  if (length(hit)) {
    warning(sprintf(
      paste0("Gaussian-process model(s) %s can combine multiple genomic/omics ",
             "kernels, but use a single-environment covariance structure; with a ",
             "multi-environment panel they tend to collapse the genetic signal ",
             "into the residual and under-estimate per-environment heritability. ",
             "For multi-environment (MET) genomic prediction use GS_model = \"GP_FA\" ",
             "(FA-GBLUP), or engine = \"asreml\" GBLUP with a variance-covariance ",
             "structure (e.g. var_cov_str = \"corgh\")."),
      paste(hit, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

gp_validate_single_trait_met_model_scope <- function(ctx, is_met = FALSE) {
  if (!is.list(ctx) || !isTRUE(is_met)) {
    return(invisible(TRUE))
  }
  specialized <- any(c(
    isTRUE(ctx[["multi_trait_gp"]]),
    isTRUE(ctx[["multi_trait_bayes"]]),
    isTRUE(ctx[["multi_trait_asreml"]]),
    isTRUE(ctx[["multi_trait_ml"]]),
    isTRUE(ctx[["multi_trait_dl"]]),
    isTRUE(ctx[["hybrid_asreml"]]),
    isTRUE(ctx[["hybrid_bayes"]]),
    isTRUE(ctx[["hybrid_gp"]]),
    isTRUE(ctx[["hybrid_ml"]]),
    isTRUE(ctx[["hybrid_dl"]])
  ))
  if (specialized) {
    return(invisible(TRUE))
  }

  requested <- unique(stats::na.omit(c(
    as.character(ctx[["GS_model"]] %||% character()),
    as.character(ctx[["GS_model_cv"]] %||% character())
  )))
  requested <- tryCatch(
    gp_canonicalize_supported_model_names(requested),
    error = function(e) requested
  )
  if (!length(requested)) {
    return(invisible(TRUE))
  }

  response_values <- NULL
  response_name <- as.character(ctx[["response"]] %||% character())
  pheno <- ctx[["pheno_data"]]
  if (length(response_name) == 1L && is.data.frame(pheno) &&
      response_name %in% names(pheno)) {
    response_values <- pheno[[response_name]]
  }
  response_family <- tryCatch(
    gp_resolve_response_family(ctx[["response_family"]] %||% "auto", y = response_values),
    error = function(e) "gaussian"
  )
  supported <- if (identical(response_family, "gaussian")) {
    gp_multi_environment_supported_models()
  } else {
    gp_met_classification_supported_models(response_family)
  }
  unsupported <- setdiff(requested, supported)
  if (length(unsupported)) {
    marker_bayes <- intersect(unsupported, c("BRR", "BayesA", "BayesB", "BayesC", "BL"))
    unsupported_display <- gp_display_supported_model_names(unsupported)
    supported_display <- gp_display_supported_model_names(supported)
    if (length(marker_bayes)) {
      stop(
        paste(
          "Bayesian marker-regression models are currently single-environment only.",
          paste0("Unsupported single-trait multi-environment model(s) for response_family='", response_family, "':"),
          paste(unsupported_display, collapse = ", "),
          "The request was stopped before fitting; no model was silently dropped.",
          "Use one of:", paste(supported_display, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    stop(
      paste(
        paste0("Unsupported single-trait multi-environment model(s) for response_family='", response_family, "':"),
        paste(unsupported_display, collapse = ", "),
        "The request was stopped before fitting; no model was silently dropped.",
        "Supported models are:", paste(supported_display, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  kernel_bayes <- intersect(requested, c("GBLUP_BRR", "RKHS"))
  if (!identical(response_family, "gaussian") && length(kernel_bayes)) {
    formula_terms <- function(x) {
      if (is.null(x)) return(character())
      entries <- if (is.list(x) && !inherits(x, "formula")) x else list(x)
      unique(unlist(lapply(entries, function(entry) {
        tryCatch(attr(stats::terms(entry), "term.labels"), error = function(e) character())
      }), use.names = FALSE))
    }
    fixed_terms <- formula_terms(ctx[["fixed"]])
    random_terms <- formula_terms(ctx[["random"]])
    gen_name <- as.character(ctx[["gen_name"]] %||% "")[[1L]]
    heter_groups <- as.character(ctx[["heter_groups"]] %||% "")[[1L]]
    mentions <- function(term, value) {
      identical(term, value) || value %in% tryCatch(
        all.vars(stats::as.formula(paste("~", term))),
        error = function(e) character()
      )
    }
    fixed_has_env <- nzchar(heter_groups) && any(vapply(fixed_terms, mentions, logical(1L), value = heter_groups))
    random_has_gen <- nzchar(gen_name) && gen_name %in% random_terms
    random_has_gxe <- nzchar(gen_name) && nzchar(heter_groups) && any(vapply(
      random_terms,
      function(term) {
        pieces <- trimws(strsplit(term, ":", fixed = TRUE)[[1L]])
        length(pieces) >= 2L && all(c(gen_name, heter_groups) %in% pieces)
      },
      logical(1L)
    ))
    if (!(fixed_has_env && random_has_gen && random_has_gxe)) {
      stop(
        paste0(
          "Single-trait multi-environment response_family='", response_family,
          "' kernel-Bayesian models [", paste(kernel_bayes, collapse = ", "),
          "] require an explicit environment main effect and GxE structure. ",
          "Use fixed = ~ ", heter_groups, " and random = ~ ", gen_name,
          " + ", gen_name, ":", heter_groups, ". The request was stopped before fitting."
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

gp_apply_model_input_guardrails <- function(ctx) {
  if (!is.list(ctx)) {
    return(ctx)
  }
  if (exists("gp_is_multi_environment_trait_panel", mode = "function")) {
    is_met <- tryCatch(
      isTRUE(gp_is_multi_environment_trait_panel(
        pheno_data = ctx[["pheno_data"]], gen_name = ctx[["gen_name"]],
        heter_groups = ctx[["heter_groups"]], response = ctx[["response"]],
        response_family = ctx[["response_family"]])),
      error = function(e) FALSE)
    gp_validate_single_trait_met_model_scope(ctx, is_met = is_met)
    specialized_met_route <- any(vapply(
      c(
        "multi_trait_gp", "multi_trait_bayes", "multi_trait_asreml",
        "multi_trait_ml", "multi_trait_dl", "hybrid_asreml",
        "hybrid_bayes", "hybrid_gp", "hybrid_ml", "hybrid_dl"
      ),
      function(flag) isTRUE(ctx[[flag]]),
      logical(1L)
    ))
    if (!isTRUE(specialized_met_route)) {
      gp_warn_single_kernel_gp_met(ctx[["GS_model"]], ctx[["GS_model_cv"]],
                                   ctx[["cross_validation"]], is_met)
    }
  }

  # Phase 3.16: auto-promote met_ml_dl=TRUE when the user listed any
  # MET-capable ML/DL model on a MET panel (heter_groups set, repeated
  # genotype rows). Run BEFORE the kernel-build auto-promote so that
  # downstream gets the flipped flag.
  ctx <- gp_auto_promote_met_ml_dl_flag(ctx)

  # When met_ml_dl=TRUE and the user supplied raw geno_data (and/or omic
  # data) but no kernel and no gmatrix_method, auto-promote the GRM build
  # so the MET ML/DL kernel-feature pipeline has something to work with.
  # Previously this silently bailed per-fold with "kernel_feature_table
  # must contain at least one genotype" because process_geno_data() only
  # builds a GRM when gmatrix_method is set.
  ctx <- gp_auto_promote_met_ml_dl_kernel_build(ctx)

  kernel_fields <- c(
    "gmatrix",
    "gkernel",
    "female_gmatrix",
    "male_gmatrix",
    "pedigree_matrix",
    "omic1_kernel",
    "omic2_kernel",
    "omic3_kernel",
    "kernel_list"
  )
  for (field in kernel_fields) {
    if (!is.null(ctx[[field]])) {
      ctx[[field]] <- gp_normalize_kernel_guardrail(ctx[[field]], field)
    }
  }
  if (!is.null(ctx[["env_similarity"]])) {
    ctx[["env_similarity"]] <- gp_normalize_env_similarity_guardrail(
      ctx[["env_similarity"]],
      env_ids = ctx[["env_ids"]],
      label = "env_similarity"
    )
  }

  if (!is.null(ctx$pheno_clean) && is.list(ctx$pheno_clean) &&
      !is.null(ctx$pheno_clean[["pheno_clean_data"]])) {
    ctx$pheno_clean[["pheno_clean_data"]] <- gp_order_pheno_for_model_guardrail(
      pheno_data = ctx$pheno_clean[["pheno_clean_data"]],
      gen_name = ctx$gen_name,
      heter_groups = ctx$heter_groups,
      order_for_asreml = gp_asreml_order_required(ctx)
    )
  }

  ctx
}

# Auto-promote a GRM build when MET ML/DL is requested with raw geno/omic data
# but no kernel/no gmatrix_method. Returns ctx unchanged if any of:
#   * met_ml_dl is not TRUE,
#   * the user already supplied a kernel (gmatrix/gkernel/omic*_kernel/
#     kernel_list with at least one non-NULL),
#   * the user already set gmatrix_method explicitly,
#   * geno_data (and train/test variants) are all NULL (nothing to build from
#     -- the contract will stop with a clear message in that case).
# When triggered, sets ctx$gmatrix_method to a sensible default ("Yang") so
# process_geno_data() builds the kernel downstream, and records a one-line
# note on ctx so the run metadata reflects the auto-promotion. The default
# is documented and overridable via ctx$gmatrix_method upfront.
# TRUE when the analysis ploidy (argument or the genotype matrix's "ploidy"
# attribute) is above two.
gp_ctx_is_polyploid <- function(ctx) {
  p <- suppressWarnings(as.integer(ctx[["ploidy"]])[1L])
  if (is.na(p)) {
    for (nm in c("geno_data", "train_geno_data", "test_geno_data")) {
      a <- attr(ctx[[nm]], "ploidy")
      if (!is.null(a)) { p <- suppressWarnings(as.integer(a)[1L]); break }
    }
  }
  isTRUE(!is.na(p) && p > 2L)
}

gp_auto_promote_met_ml_dl_kernel_build <- function(ctx) {
  if (!is.list(ctx) || !isTRUE(ctx[["met_ml_dl"]])) {
    return(ctx)
  }
  if (!is.null(ctx[["gmatrix_method"]]) && nzchar(as.character(ctx[["gmatrix_method"]]))) {
    return(ctx)
  }
  has_kernel <- function(x) !is.null(x) && (
    (is.list(x) && length(x) > 0L) ||
    (is.matrix(x) || inherits(x, "Matrix") || is.array(x))
  )
  any_kernel <- has_kernel(ctx[["gmatrix"]]) || has_kernel(ctx[["gkernel"]]) ||
    has_kernel(ctx[["omic1_kernel"]]) || has_kernel(ctx[["omic2_kernel"]]) ||
    has_kernel(ctx[["omic3_kernel"]]) || has_kernel(ctx[["kernel_list"]])
  if (any_kernel) {
    return(ctx)
  }
  has_geno <- !is.null(ctx[["geno_data"]]) ||
    !is.null(ctx[["train_geno_data"]]) ||
    !is.null(ctx[["test_geno_data"]])
  if (!has_geno) {
    # Nothing to build from. Let the input contract stop with its standard
    # MET-ML/DL message; do not silently no-op-promote.
    return(ctx)
  }
  # Yang is diploid-only; polyploid dosage uses the ploidy-aware VanRaden GRM
  auto_method <- if (gp_ctx_is_polyploid(ctx)) "VanRaden" else "Yang"
  ctx[["gmatrix_method"]] <- auto_method
  notes <- ctx[["model_input_guardrail_notes"]] %||% character()
  ctx[["model_input_guardrail_notes"]] <- c(
    notes,
    paste0("met_ml_dl: auto-set gmatrix_method = '", auto_method, "' because no kernel was supplied; the MET ML/DL kernel-feature pipeline requires a kernel and one is now built from geno_data.")
  )
  message(
    "PredictProR MET ML/DL: no kernel supplied. Auto-building a ", auto_method, " GRM from geno_data ",
    "so the kernel-feature pipeline has a kernel to PCA over. ",
    "Set gmatrix_method explicitly (e.g. 'VanRaden') or pass a precomputed gmatrix/kernel_list to override."
  )
  ctx
}

# Phase 3.16: auto-promote met_ml_dl = TRUE when the user lists any
# MET-capable ML/DL model on a MET panel. Before this, users had to know to
# set met_ml_dl=TRUE explicitly; otherwise ML/DL models silently fell through
# to the single-env predict_with_model() path and produced all-NA predictions
# on MET data (CV: silent drop in gp_filter_crossval_results;
# true-prediction: error "only works for single location"). This sister
  # function complements gp_auto_promote_met_ml_dl_kernel_build by flipping
  # the gate flag, not the kernel builder. Unsupported mixed requests are
  # rejected by gp_validate_single_trait_met_model_scope before this helper.
gp_auto_promote_met_ml_dl_flag <- function(ctx) {
  if (!is.list(ctx)) return(ctx)
  # Already on -> no-op
  if (isTRUE(ctx[["met_ml_dl"]])) return(ctx)

  # MET context probe: heter_groups column present and has >= 2 unique levels
  # in the pheno frame; pheno also has repeated genotype rows.
  heter_groups <- ctx[["heter_groups"]]
  if (is.null(heter_groups) || !nzchar(as.character(heter_groups)[1L])) return(ctx)
  pheno <- ctx[["pheno_data"]]
  if (is.null(pheno) || !is.data.frame(pheno)) return(ctx)
  hg_col <- as.character(heter_groups)[1L]
  if (!(hg_col %in% names(pheno))) return(ctx)
  hg_levels <- unique(as.character(pheno[[hg_col]]))
  hg_levels <- hg_levels[!is.na(hg_levels) & nzchar(hg_levels)]
  if (length(hg_levels) < 2L) return(ctx)
  gen_name <- ctx[["gen_name"]]
  if (is.null(gen_name) || !nzchar(as.character(gen_name)[1L])) return(ctx)
  gid_col <- as.character(gen_name)[1L]
  if (!(gid_col %in% names(pheno))) return(ctx)
  has_repeated_geno <- anyDuplicated(as.character(pheno[[gid_col]])) > 0L
  if (!has_repeated_geno) return(ctx)

  gp_validate_single_trait_met_model_scope(ctx, is_met = TRUE)

  # MET-capable ML/DL models (canonical names after canonicalization).
  met_capable <- tryCatch(gp_met_supported_models(),
                          error = function(e) c("CatBoost", "LightGBM",
                                                "Xgboost", "RandomForest",
                                                "mlp", "ft_transformer",
                                                "saint", "tabnet", "moe"))
  user_models <- c(as.character(ctx[["GS_model"]] %||% character()),
                   as.character(ctx[["GS_model_cv"]] %||% character()))
  user_models_canon <- tryCatch(
    gp_canonicalize_supported_model_names(user_models),
    error = function(e) user_models)
  met_in_list <- intersect(user_models_canon, met_capable)
  if (length(met_in_list) == 0L) return(ctx)

  # Phase 3.16 / Phase 3.20 follow-up: PredictProR has four canonical
  # fit-context cases:
  #   1. single-env, single-trait
  #   2. single-env, multi-trait   (user flips multi_trait_*)
  #   3. MET, single-trait         (auto-promote met_ml_dl)
  #   4. MET, multi-trait          (per-trait iteration by default; joint
  #      via multi_trait_gp / multi_trait_bayes / multi_trait_asreml)
  # We auto-promote met_ml_dl=TRUE for cases 3 AND 4 (per-trait). The CV
  # dispatcher in main_crossvalidation_execution_logic.R iterates per-trait
  # per fold and passes `response = trait` (single string) into
  # gp_met_python_cv_predict / gp_met_dl_cv_predict, so the MET ML/DL
  # pipeline only ever sees a single response per call -- multi-trait MET
  # is just a sequence of single-trait MET fits, one per trait.
  # When a joint flag is set, the joint engine owns dispatch and we DEFER.
  is_specialized_path <- any(c(
    isTRUE(ctx[["multi_trait_gp"]]),
    isTRUE(ctx[["multi_trait_bayes"]]),
    isTRUE(ctx[["multi_trait_asreml"]]),
    isTRUE(ctx[["multi_trait_ml"]]),
    isTRUE(ctx[["multi_trait_dl"]]),
    isTRUE(ctx[["hybrid_asreml"]]),
    isTRUE(ctx[["hybrid_bayes"]]),
    isTRUE(ctx[["hybrid_gp"]]),
    isTRUE(ctx[["hybrid_ml"]]),
    isTRUE(ctx[["hybrid_dl"]])
  ))
  if (is_specialized_path) {
    # A dedicated multi-trait or hybrid engine owns dispatch. Defer generic
    # MET ML/DL auto-promotion even when the specialized response is scalar.
    notes <- ctx[["model_input_guardrail_notes"]] %||% character()
    ctx[["model_input_guardrail_notes"]] <- c(
      notes,
      "met_ml_dl auto-promote DEFERRED to a dedicated multi-trait or hybrid workflow; that engine owns MET dispatch."
    )
    return(ctx)
  }

  # Case 3 (single-trait MET) and Case 4 (per-trait multi-trait MET): safe to flip.
  ctx[["met_ml_dl"]] <- TRUE
  notes <- ctx[["model_input_guardrail_notes"]] %||% character()
  ctx[["model_input_guardrail_notes"]] <- c(
    notes,
    sprintf(
      "met_ml_dl: auto-set to TRUE because heter_groups='%s' (%d levels), repeated genotype rows present, and MET-capable ML/DL model(s) were listed: %s",
      hg_col, length(hg_levels), paste(met_in_list, collapse = ", ")
    )
  )
  message(sprintf(
    "PredictProR MET ML/DL: auto-set `met_ml_dl = TRUE` because you listed MET-capable ML/DL model(s) [%s] alongside heter_groups='%s' with repeated genotype rows. Override by setting `met_ml_dl = FALSE` explicitly.",
    paste(met_in_list, collapse = ", "), hg_col
  ))

  ctx
}
