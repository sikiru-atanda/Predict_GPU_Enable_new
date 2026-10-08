#' Score and rank predictors for dynamic top-k model fitting
#'
#' @param predictor_data Numeric predictor matrix or data frame with individuals
#'   in rows and predictors in columns.
#' @param pheno_data Phenotype data containing `gen_name` and response columns.
#' @param response Character vector of response variables to score.
#' @param gen_name Column name in `pheno_data` containing individual IDs.
#' @param scoring_model One of `"Ridge_Regression"`, `"BayesB"`, or
#'   `"RandomForest"`.
#' @param seed Integer seed used by stochastic scoring backends.
#' @param source_block Optional predictor source labels. Either a named vector
#'   keyed by predictor name, or a single label recycled to all predictors.
#' @param ridge_lambda Positive ridge penalty for ridge scoring.
#' @param bayes_nIter Iterations for BayesB scoring.
#' @param bayes_burnIn Burn-in for BayesB scoring.
#' @param bayes_thin Thinning for BayesB scoring.
#' @param ntree Number of trees for Python RandomForest scoring.
#' @param mtry Number of predictors sampled at each split.
#' @param nodesize Minimum node size.
#' @param rf_n_jobs Internal Python RandomForest jobs. Keep this at 1 when the
#'   caller is already parallelizing over traits, models, folds, or k values.
#' @param response_family Response family used for scoring. `"auto"` infers it
#'   separately for every trait. Ridge and BayesB scoring are restricted to
#'   Gaussian responses; use RandomForest for binary, ordinal, or multiclass
#'   responses.
#' @param selection_mode Provenance label for the ranking, normally `"fixed"`
#'   or `"fold_internal"`.
#' @param replication,fold Optional resampling identifiers recorded in the
#'   returned metadata.
#'
#' @return A `predictpror_feature_scores` data frame with one row per
#'   trait/predictor and columns for source, score, rank, response family,
#'   selection mode, replication/fold, training-sample count, and a deterministic
#'   training-ID hash. These are predictive feature scores, not genetic variance
#'   components or causal effects.
#' @export
feature_score_predictors <- function(predictor_data,
                                     pheno_data,
                                     response,
                                     gen_name,
                                     scoring_model = c("Ridge_Regression", "BayesB", "RandomForest"),
                                     seed = 123L,
                                     source_block = NULL,
                                     ridge_lambda = 1,
                                     bayes_nIter = 1500L,
                                     bayes_burnIn = 500L,
                                     bayes_thin = 5L,
                                     ntree = 500L,
                                     mtry = NULL,
                                     nodesize = NULL,
                                     rf_n_jobs = 1L,
                                     response_family = "auto",
                                     selection_mode = "fixed",
                                     replication = NA_integer_,
                                     fold = NA_integer_) {
  scoring_model <- match.arg(scoring_model)
  selection_mode <- gp_feature_normalize_cv_policy(selection_mode, allow_both = FALSE)
  if (is.null(predictor_data) || !NROW(predictor_data) || !NCOL(predictor_data)) {
    stop("predictor_data must contain at least one row and one predictor column.", call. = FALSE)
  }
  if (is.null(rownames(predictor_data))) {
    stop("predictor_data must have row names containing individual IDs.", call. = FALSE)
  }
  if (is.null(pheno_data) || !gen_name %in% names(pheno_data)) {
    stop("pheno_data must contain the gen_name column.", call. = FALSE)
  }
  response <- as.character(response %||% character())
  if (!length(response) || !all(response %in% names(pheno_data))) {
    stop("All response variables must exist in pheno_data.", call. = FALSE)
  }

  x <- as.matrix(predictor_data)
  storage.mode(x) <- "double"
  predictor_names <- colnames(x)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("feature_", seq_len(ncol(x)))
    colnames(x) <- predictor_names
  }
  source_block <- gp_feature_source_block(predictor_names, source_block)
  pheno_ids <- as.character(pheno_data[[gen_name]])
  common_ids <- intersect(rownames(x), pheno_ids)
  if (!length(common_ids)) {
    stop("No common individual IDs were found between predictor_data and pheno_data.", call. = FALSE)
  }

  out <- lapply(response, function(trait) {
    y_all <- pheno_data[[trait]]
    trait_family <- gp_resolve_response_family(response_family, y = y_all)
    if (!identical(trait_family, "gaussian") && !identical(scoring_model, "RandomForest")) {
      stop(
        scoring_model,
        " feature scoring is only defined for gaussian responses. Use scoring_model = 'RandomForest' for ",
        trait_family,
        " trait '", trait, "'.",
        call. = FALSE
      )
    }
    valid <- if (is.numeric(y_all) || is.integer(y_all)) {
      is.finite(as.double(y_all))
    } else {
      !is.na(y_all)
    }
    valid <- valid & pheno_ids %in% rownames(x)
    ids <- pheno_ids[valid]
    if (!length(ids)) {
      stop("Trait ", trait, " has no observed training values for feature scoring.", call. = FALSE)
    }
    x_train <- x[ids, , drop = FALSE]
    y_train <- y_all[valid]
    scores <- switch(
      scoring_model,
      Ridge_Regression = gp_feature_score_ridge(x_train, y_train, lambda = ridge_lambda),
      BayesB = gp_feature_score_bayesb(
        x_train,
        y_train,
        seed = seed,
        nIter = bayes_nIter,
        burnIn = bayes_burnIn,
        thin = bayes_thin
      ),
      RandomForest = gp_feature_score_python_randomforest(
        x_train,
        y_train,
        response_family = trait_family,
        seed = seed,
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        n_jobs = rf_n_jobs
      )
    )
    gp_feature_metadata_table(
      trait = trait,
      predictors = predictor_names,
      source_block = source_block,
      raw_score = scores,
      scoring_model = scoring_model,
      seed = seed,
      n_train = length(ids),
      response_family = trait_family,
      selection_mode = selection_mode,
      replication = replication,
      fold = fold,
      training_id_hash = gp_feature_training_id_hash(ids)
    )
  })

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  class(ans) <- c("predictpror_feature_scores", class(ans))
  ans
}

#' Select top-k predictors from feature score metadata
#'
#' @param predictor_data Matrix or data frame to subset.
#' @param feature_score_metadata Output from `feature_score_predictors()`.
#' @param trait Trait name to use.
#' @param k Number of predictors to keep, or `"all"`.
#'
#' @return Predictor data subset to the top-k ranked columns.
#' @export
feature_select_top_k <- function(predictor_data,
                                 feature_score_metadata,
                                 trait,
                                 k) {
  if (is.null(feature_score_metadata) || !NROW(feature_score_metadata)) {
    return(predictor_data)
  }
  trait <- as.character(trait[[1L]])
  meta <- feature_score_metadata[as.character(feature_score_metadata$trait) == trait, , drop = FALSE]
  if (!nrow(meta)) {
    return(predictor_data)
  }
  meta <- meta[order(meta$rank, meta$predictor), , drop = FALSE]
  k_val <- gp_feature_normalize_k(k, nrow(meta))
  keep <- as.character(meta$predictor[seq_len(k_val)])
  keep <- keep[keep %in% colnames(predictor_data)]
  predictor_data[, keep, drop = FALSE]
}

gp_feature_source_aliases <- function(source = NULL) {
  aliases <- list(
    geno_data = c(
      "geno_data", "train_geno_data", "test_geno_data", "geno_model_ready",
      "female_geno_data", "male_geno_data",
      "geno", "genotype", "genotypes", "marker", "markers", "snp", "snps", "m"
    ),
    omic1_data = c(
      "omic1_data", "train_omic1_data", "test_omic1_data", "omic1_model_ready",
      "omic1", "omics1", "omic", "omics", "omic_data", "omics_data"
    ),
    omic2_data = c(
      "omic2_data", "train_omic2_data", "test_omic2_data", "omic2_model_ready",
      "omic2", "omics2"
    ),
    omic3_data = c(
      "omic3_data", "train_omic3_data", "test_omic3_data", "omic3_model_ready",
      "omic3", "omics3"
    )
  )
  if (is.null(source)) {
    return(unique(unlist(aliases, use.names = FALSE)))
  }
  source <- tolower(as.character(source[[1L]]))
  canonical <- names(aliases)[vapply(aliases, function(x) source %in% tolower(x), logical(1L))]
  if (!length(canonical)) {
    return(unique(c(source)))
  }
  unique(c(source, aliases[[canonical[[1L]]]]))
}

gp_feature_match_named_entry <- function(nms, keys) {
  if (is.null(nms) || !length(nms) || !length(keys)) {
    return(NULL)
  }
  hit <- match(tolower(keys), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    return(NULL)
  }
  nms[[hit[[1L]]]]
}

gp_feature_list_has_source_keys <- function(nms) {
  if (is.null(nms) || !length(nms)) {
    return(FALSE)
  }
  any(tolower(nms) %in% tolower(gp_feature_source_aliases()))
}

gp_feature_has_source_scope <- function(feature_selected) {
  if (is.data.frame(feature_selected)) {
    return(any(c(
      "source", "source_block", "block", "input", "data_source", "feature_source"
    ) %in% names(feature_selected)))
  }
  if (!is.list(feature_selected)) return(FALSE)
  if (gp_feature_list_has_source_keys(names(feature_selected))) return(TRUE)
  any(vapply(feature_selected, gp_feature_has_source_scope, logical(1L)))
}

gp_feature_selected_vector <- function(feature_selected, trait = NULL, source = NULL) {
  if (is.null(feature_selected)) {
    return(NULL)
  }
  if (is.data.frame(feature_selected)) {
    dat <- feature_selected
    if (!is.null(trait) && "trait" %in% names(dat)) {
      dat <- dat[as.character(dat[["trait"]]) == as.character(trait), , drop = FALSE]
    }
    if (!is.null(source)) {
      source_col <- intersect(
        c("source", "source_block", "block", "input", "data_source", "feature_source"),
        names(dat)
      )
      if (length(source_col)) {
        source_aliases <- tolower(gp_feature_source_aliases(source))
        dat <- dat[tolower(as.character(dat[[source_col[[1L]]]])) %in% source_aliases, , drop = FALSE]
      }
    }
    predictor_col <- intersect(
      c("predictor", "feature", "marker", "Marker", "Feature"),
      names(dat)
    )
    selected_flag_col <- intersect(c("selected", "Selected"), names(dat))
    if (length(predictor_col) && length(selected_flag_col) &&
        is.logical(dat[[selected_flag_col[[1L]]]])) {
      flag <- dat[[selected_flag_col[[1L]]]]
      dat <- dat[!is.na(flag) & flag, , drop = FALSE]
    }
    if (!nrow(dat)) {
      return(NULL)
    }
    value_col <- intersect(
      c("predictor", "feature", "marker", "Marker", "Feature", "selected", "Selected"),
      names(dat)
    )
    if (length(value_col)) {
      return(dat[[value_col[[1L]]]])
    }
    if (ncol(dat)) {
      return(dat[[1L]])
    }
    return(NULL)
  }
  if (is.list(feature_selected) && !is.data.frame(feature_selected)) {
    nms <- names(feature_selected)
    if (!is.null(source) && length(nms)) {
      source_hit <- gp_feature_match_named_entry(nms, gp_feature_source_aliases(source))
      if (!is.null(source_hit)) {
        return(gp_feature_selected_vector(feature_selected[[source_hit]], trait = trait, source = NULL))
      }
      if (gp_feature_list_has_source_keys(nms)) {
        return(NULL)
      }
    }
    if (!is.null(trait) && length(nms) && as.character(trait) %in% nms) {
      return(gp_feature_selected_vector(feature_selected[[as.character(trait)]], trait = NULL, source = source))
    }
    if (!is.null(source) && any(vapply(feature_selected, gp_feature_has_source_scope, logical(1L)))) {
      nested <- lapply(feature_selected, function(value) {
        gp_feature_selected_vector(value, trait = NULL, source = source)
      })
      nested <- Filter(length, nested)
      if (!length(nested)) return(NULL)
      return(unlist(nested, use.names = FALSE))
    }
    common_keys <- c("predictor", "predictors", "feature", "features", "marker", "markers", "selected")
    hit <- intersect(common_keys, nms %||% character())
    if (length(hit)) {
      return(gp_feature_selected_vector(feature_selected[[hit[[1L]]]], trait = NULL, source = source))
    }
    return(unlist(feature_selected, use.names = FALSE))
  }
  feature_selected
}

gp_feature_resolve_selected_names <- function(x,
                                              feature_selected,
                                              trait = NULL,
                                              source = NULL,
                                              id_col = NULL) {
  if (is.null(x) || is.null(feature_selected) || is.null(colnames(x))) {
    return(character())
  }
  selected <- gp_feature_selected_vector(
    feature_selected,
    trait = trait,
    source = source
  )
  if (is.null(selected) || !length(selected)) {
    return(character())
  }
  predictor_names <- colnames(x)
  id_col <- as.character(id_col %||% character())
  id_col <- id_col[id_col %in% predictor_names]
  feature_names <- setdiff(predictor_names, id_col)
  if (is.logical(selected)) {
    if (length(selected) != length(feature_names)) {
      stop("feature_selected logical vector must have one value per predictor column.", call. = FALSE)
    }
    selected[is.na(selected)] <- FALSE
    keep <- feature_names[selected]
  } else if (is.numeric(selected) || is.integer(selected)) {
    idx <- as.integer(selected)
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(feature_names)]
    keep <- feature_names[idx]
  } else {
    selected <- as.character(selected)
    selected <- selected[!is.na(selected) & nzchar(selected)]
    keep <- selected[selected %in% feature_names]
  }
  unique(keep)
}

gp_feature_select_explicit <- function(x,
                                       feature_selected,
                                       trait = NULL,
                                       context = "predictor_data",
                                       source = NULL,
                                       id_col = NULL) {
  if (is.null(x) || is.null(feature_selected) || is.null(colnames(x))) {
    return(x)
  }
  predictor_names <- colnames(x)
  id_col <- as.character(id_col %||% character())
  id_col <- id_col[id_col %in% predictor_names]
  keep <- gp_feature_resolve_selected_names(
    x = x,
    feature_selected = feature_selected,
    trait = trait,
    source = source,
    id_col = id_col
  )
  declared <- gp_feature_selected_vector(feature_selected, trait = trait, source = source)
  if (is.numeric(declared) || is.integer(declared)) {
    idx <- suppressWarnings(as.integer(declared))
    invalid <- !is.finite(idx) | idx < 1L | idx > length(setdiff(predictor_names, id_col))
    if (any(invalid)) {
      stop("feature_selected contains predictor indices outside the available columns in ", context, ".", call. = FALSE)
    }
  } else if (is.character(declared)) {
    declared <- unique(declared[!is.na(declared) & nzchar(declared)])
    missing <- setdiff(declared, setdiff(predictor_names, id_col))
    if (length(missing)) {
      stop(
        "feature_selected columns are missing from ", context, ": ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!length(keep)) {
    stop("feature_selected did not match any predictor columns in ", context, ".", call. = FALSE)
  }
  out <- x[, c(id_col, keep), drop = FALSE]
  attr(out, "predictpror_feature_selected") <- keep
  out
}

gp_feature_select_raw_input <- function(x,
                                        feature_selected,
                                        source,
                                        gen_name = NULL) {
  if (is.null(x) || is.null(feature_selected)) {
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    if ("snps_matrix" %in% names(x)) {
      out <- x
      selected_matrix <- gp_feature_select_raw_input(
        x = x[["snps_matrix"]],
        feature_selected = feature_selected,
        source = source,
        gen_name = gen_name
      )
      if (is.null(selected_matrix)) return(NULL)
      out[["snps_matrix"]] <- selected_matrix
      return(out)
    }
    return(x)
  }
  if (!is.data.frame(x) && !is.matrix(x)) {
    return(x)
  }
  keep <- gp_feature_resolve_selected_names(
    x = x,
    feature_selected = feature_selected,
    trait = NULL,
    source = source,
    id_col = gen_name
  )
  if (!length(keep)) return(NULL)
  id_col <- as.character(gen_name %||% character())
  id_col <- id_col[id_col %in% colnames(x)]
  out <- x[, c(id_col, keep), drop = FALSE]
  attr(out, "predictpror_feature_selected") <- keep
  out
}

gp_feature_apply_raw_selection_to_inputs <- function(ctx) {
  if (!is.list(ctx) || is.null(ctx$feature_selected)) {
    return(ctx)
  }
  source_fields <- list(
    geno_data = c("geno_data", "train_geno_data", "test_geno_data", "female_geno_data", "male_geno_data"),
    omic1_data = c("omic1_data", "train_omic1_data", "test_omic1_data"),
    omic2_data = c("omic2_data", "train_omic2_data", "test_omic2_data"),
    omic3_data = c("omic3_data", "train_omic3_data", "test_omic3_data")
  )
  all_fields <- unique(unlist(source_fields, use.names = FALSE))
  had_raw_predictors <- any(vapply(all_fields, function(field) !is.null(ctx[[field]]), logical(1L)))
  source_columns <- lapply(source_fields, function(fields) {
    unique(unlist(lapply(fields, function(field) {
      x <- ctx[[field]]
      if (is.list(x) && !is.data.frame(x) && "snps_matrix" %in% names(x)) {
        x <- x[["snps_matrix"]]
      }
      cols <- colnames(x) %||% character()
      setdiff(cols, as.character(ctx$gen_name %||% character()))
    }), use.names = FALSE))
  })
  if (gp_feature_has_source_scope(ctx$feature_selected)) {
    for (source in names(source_columns)) {
      declared <- gp_feature_selected_vector(ctx$feature_selected, trait = NULL, source = source)
      if (is.character(declared)) {
        missing <- setdiff(unique(declared[!is.na(declared) & nzchar(declared)]), source_columns[[source]])
        if (length(missing)) {
          stop(
            "feature_selected columns are missing from ", source, ": ",
            paste(missing, collapse = ", "),
            call. = FALSE
          )
        }
      }
    }
  } else {
    declared <- gp_feature_selected_vector(ctx$feature_selected, trait = NULL, source = NULL)
    if (is.character(declared)) {
      missing <- setdiff(
        unique(declared[!is.na(declared) & nzchar(declared)]),
        unique(unlist(source_columns, use.names = FALSE))
      )
      if (length(missing)) {
        stop(
          "feature_selected columns are missing from the supplied marker/omics inputs: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }
  for (source in names(source_fields)) {
    for (field in source_fields[[source]]) {
      if (!is.null(ctx[[field]])) {
        ctx[[field]] <- gp_feature_select_raw_input(
          x = ctx[[field]],
          feature_selected = ctx$feature_selected,
          source = source,
          gen_name = ctx$gen_name
        )
      }
    }
  }
  has_selected_predictors <- any(vapply(all_fields, function(field) !is.null(ctx[[field]]), logical(1L)))
  if (isTRUE(had_raw_predictors) && !isTRUE(has_selected_predictors)) {
    stop(
      "feature_selected did not match any named columns in the supplied marker/omics inputs.",
      call. = FALSE
    )
  }
  ctx
}

gp_feature_apply_explicit_to_ml_dat_res <- function(ml_dat_res,
                                                    feature_selected,
                                                    trait = NULL) {
  if (is.null(ml_dat_res) || is.null(feature_selected)) {
    return(ml_dat_res)
  }
  out <- ml_dat_res
  train_x <- out[["merged_data"]][["merge_data"]]
  if (is.null(train_x)) {
    return(out)
  }
  selected_train <- gp_feature_select_explicit(
    train_x,
    feature_selected = feature_selected,
    trait = trait,
    context = "merged ML predictor data"
  )
  keep <- colnames(selected_train)
  out[["merged_data"]][["merge_data"]] <- selected_train
  if (!is.null(out[["merged_data_test"]])) {
    missing_test <- setdiff(keep, colnames(out[["merged_data_test"]]))
    if (length(missing_test)) {
      stop("feature_selected columns are missing from the ML test predictor block: ",
           paste(missing_test, collapse = ", "), call. = FALSE)
    }
    out[["merged_data_test"]] <- out[["merged_data_test"]][, keep, drop = FALSE]
  }
  out[["feature_selected"]] <- keep
  out[["feature_k"]] <- length(keep)
  out
}

gp_feature_source_block <- function(predictors, source_block = NULL) {
  if (is.null(source_block)) {
    return(stats::setNames(rep("merged", length(predictors)), predictors))
  }
  if (length(source_block) == 1L && is.null(names(source_block))) {
    return(stats::setNames(rep(as.character(source_block), length(predictors)), predictors))
  }
  out <- as.character(source_block[predictors])
  out[is.na(out) | !nzchar(out)] <- "merged"
  stats::setNames(out, predictors)
}

gp_feature_normalize_cv_policy <- function(policy = "fixed", allow_both = TRUE) {
  policy <- tolower(trimws(as.character(policy %||% "fixed")[[1L]]))
  aliases <- c(
    global = "fixed",
    full = "fixed",
    fixed = "fixed",
    fold = "fold_internal",
    internal = "fold_internal",
    nested = "fold_internal",
    fold_internal = "fold_internal",
    final_refit = "final_refit",
    both = "both"
  )
  normalized <- unname(aliases[policy])
  valid <- c("fixed", "fold_internal", "final_refit", if (isTRUE(allow_both)) "both")
  if (is.na(normalized) || !normalized %in% valid) {
    stop(
      "feature_scoring_cv must be one of: ", paste(valid, collapse = ", "),
      call. = FALSE
    )
  }
  normalized
}

gp_feature_cv_modes <- function(policy = "fixed") {
  policy <- gp_feature_normalize_cv_policy(policy, allow_both = TRUE)
  if (identical(policy, "final_refit")) {
    return("fixed")
  }
  if (identical(policy, "both")) c("fixed", "fold_internal") else policy
}

gp_feature_specialized_cv_config <- function(feature_score_metadata = NULL,
                                             feature_k = NULL,
                                             feature_k_grid = NULL,
                                             feature_scoring_cv = "fixed",
                                             context = "specialized model") {
  if (is.null(feature_score_metadata)) {
    return(list(active = FALSE, k = NULL, mode = "fixed"))
  }
  k_values <- feature_k %||% feature_k_grid
  k_values <- unique(stats::na.omit(as.integer(k_values)))
  if (length(k_values) != 1L) {
    stop(
      context,
      " requires one pre-specified `feature_k` (or a length-one `feature_k_grid`) per run. Multiple k values require nested model tuning and are not silently averaged.",
      call. = FALSE
    )
  }
  mode <- gp_feature_normalize_cv_policy(feature_scoring_cv, allow_both = TRUE)
  if (identical(mode, "both")) {
    stop(
      context,
      " requires choosing either `feature_scoring_cv = 'fixed'` or `'fold_internal'`; run both modes separately so their estimates are not pooled.",
      call. = FALSE
    )
  }
  list(active = TRUE, k = k_values[[1L]], mode = mode)
}

gp_feature_training_id_hash <- function(ids) {
  ids <- sort(unique(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  bytes <- utf8ToInt(paste(ids, collapse = "\r"))
  hash <- 5381
  if (length(bytes)) {
    for (value in bytes) {
      hash <- (hash * 33 + value) %% 2147483647
    }
  }
  sprintf("n%s-h%08x", length(ids), as.integer(hash))
}

gp_feature_source_map_from_inputs <- function(predictor_data,
                                              source_matrices = list()) {
  predictors <- colnames(predictor_data) %||% character()
  out <- stats::setNames(rep("merged", length(predictors)), predictors)
  if (!length(predictors) || !length(source_matrices)) {
    return(out)
  }
  membership <- lapply(source_matrices, function(x) {
    if (is.null(x) || is.null(colnames(x))) character() else colnames(x)
  })
  for (predictor in predictors) {
    hits <- names(membership)[vapply(membership, function(x) predictor %in% x, logical(1L))]
    if (length(hits) == 1L) {
      out[[predictor]] <- hits[[1L]]
    } else if (length(hits) > 1L) {
      stop(
        "Predictor '", predictor,
        "' occurs in multiple feature sources. Make predictor names unique before feature selection.",
        call. = FALSE
      )
    }
  }
  out
}

gp_feature_score_training_ids <- function(predictor_data,
                                          pheno_data,
                                          response,
                                          gen_name,
                                          training_ids,
                                          scoring_model,
                                          seed,
                                          source_block = NULL,
                                          response_family = "auto",
                                          selection_mode = "fold_internal",
                                          replication = NA_integer_,
                                          fold = NA_integer_,
                                          ridge_lambda = 1,
                                          bayes_nIter = 1500L,
                                          bayes_burnIn = 500L,
                                          bayes_thin = 5L,
                                          ntree = 500L,
                                          mtry = NULL,
                                          nodesize = NULL,
                                          rf_n_jobs = 1L) {
  training_ids <- unique(as.character(training_ids))
  ph <- pheno_data[as.character(pheno_data[[gen_name]]) %in% training_ids, , drop = FALSE]
  if (!nrow(ph)) {
    stop("Fold-internal feature scoring received no training phenotypes.", call. = FALSE)
  }
  feature_score_predictors(
    predictor_data = predictor_data,
    pheno_data = ph,
    response = response,
    gen_name = gen_name,
    scoring_model = scoring_model,
    seed = seed,
    source_block = source_block,
    ridge_lambda = ridge_lambda,
    bayes_nIter = bayes_nIter,
    bayes_burnIn = bayes_burnIn,
    bayes_thin = bayes_thin,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    rf_n_jobs = rf_n_jobs,
    response_family = response_family,
    selection_mode = selection_mode,
    replication = replication,
    fold = fold
  )
}

gp_feature_selected_metadata <- function(feature_score_metadata,
                                         trait,
                                         k) {
  if (is.null(feature_score_metadata) || !nrow(feature_score_metadata)) {
    if (is.null(feature_score_metadata)) return(data.frame())
    return(feature_score_metadata[FALSE, , drop = FALSE])
  }
  meta <- feature_score_metadata[
    as.character(feature_score_metadata[["trait"]]) == as.character(trait),
    ,
    drop = FALSE
  ]
  meta <- meta[order(meta[["rank"]], meta[["predictor"]]), , drop = FALSE]
  if (!nrow(meta)) {
    return(meta)
  }
  meta[seq_len(gp_feature_normalize_k(k, nrow(meta))), , drop = FALSE]
}

gp_feature_selected_by_source <- function(feature_score_metadata,
                                          trait,
                                          k) {
  selected <- gp_feature_selected_metadata(feature_score_metadata, trait = trait, k = k)
  if (is.null(selected) || !nrow(selected)) {
    return(list())
  }
  sources <- as.character(selected[["source_block"]] %||% rep("merged", nrow(selected)))
  sources[is.na(sources) | !nzchar(sources)] <- "merged"
  split(as.character(selected[["predictor"]]), sources)
}

gp_feature_explicit_selected_by_source <- function(feature_selected,
                                                    trait = NULL,
                                                    source_matrices = list(),
                                                    context = "feature_selected") {
  canonical_sources <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")
  source_matrices <- source_matrices %||% list()
  available <- canonical_sources[vapply(
    canonical_sources,
    function(source) !is.null(source_matrices[[source]]),
    logical(1L)
  )]
  selected_by_source <- list()
  for (source in available) {
    keep <- gp_feature_resolve_selected_names(
      x = source_matrices[[source]],
      feature_selected = feature_selected,
      trait = trait,
      source = source
    )
    if (length(keep)) selected_by_source[[source]] <- keep
  }
  if (length(available) && !length(selected_by_source)) {
    stop(
      context,
      " did not match any named columns for trait '",
      as.character(trait %||% "joint"),
      "' in the supplied marker/omics inputs.",
      call. = FALSE
    )
  }
  selected_by_source
}

gp_feature_declared_by_source <- function(feature_selected, trait = NULL) {
  if (is.null(feature_selected)) return(list())
  if (gp_feature_has_source_scope(feature_selected)) {
    out <- list()
    for (source in c("geno_data", "omic1_data", "omic2_data", "omic3_data")) {
      declared <- gp_feature_selected_vector(feature_selected, trait = trait, source = source)
      if (is.character(declared)) {
        declared <- unique(declared[!is.na(declared) & nzchar(declared)])
        if (length(declared)) out[[source]] <- declared
      }
    }
    return(out)
  }
  declared <- gp_feature_selected_vector(feature_selected, trait = trait)
  if (!is.character(declared)) return(list())
  declared <- unique(declared[!is.na(declared) & nzchar(declared)])
  if (length(declared)) list(merged = declared) else list()
}

gp_feature_explicit_metadata <- function(selected_by_source,
                                         trait,
                                         replication = NA_integer_,
                                         fold = NA_integer_,
                                         prediction_model = NA_character_,
                                         n_train = NA_integer_) {
  selected_by_source <- Filter(length, selected_by_source %||% list())
  if (!length(selected_by_source)) return(NULL)
  source_block <- rep(names(selected_by_source), lengths(selected_by_source))
  predictors <- unlist(selected_by_source, use.names = FALSE)
  total <- length(predictors)
  data.frame(
    trait = rep(as.character(trait), total),
    predictor = as.character(predictors),
    source_block = as.character(source_block),
    rank = seq_len(total),
    raw_score = rep(NA_real_, total),
    normalized_score = rep(NA_real_, total),
    scoring_model = rep("user_supplied", total),
    response_family = rep(NA_character_, total),
    selection_mode = rep("explicit", total),
    replication = rep(as.integer(replication), total),
    fold = rep(as.integer(fold), total),
    training_id_hash = rep(NA_character_, total),
    seed = rep(NA_integer_, total),
    n_train = rep(as.integer(n_train), total),
    selected = rep(TRUE, total),
    feature_k = rep(as.integer(total), total),
    prediction_model = rep(as.character(prediction_model), total),
    stringsAsFactors = FALSE
  )
}

gp_feature_model_needs_kernel <- function(model, met_ml_dl = FALSE) {
  if (is.null(model) || !length(model) || is.na(model[[1L]])) return(FALSE)
  model <- as.character(model[[1L]])
  gp_models <- if (exists("gp_lowrank_supported_models", mode = "function")) {
    gp_lowrank_supported_models()
  } else {
    character()
  }
  met_ai_models <- c(
    "CatBoost", "LightGBM", "Xgboost", "RandomForest",
    "mlp", "ft_transformer", "saint", "tabnet", "moe"
  )
  model %in% c("GBLUP_BRR", "RKHS", "GBLUP", gp_models) ||
    (isTRUE(met_ml_dl) && model %in% met_ai_models)
}

gp_feature_fold_view <- function(predictor_data,
                                 pheno_data,
                                 response,
                                 gen_name,
                                 test_rows,
                                 k,
                                 scoring_model,
                                 seed,
                                 source_block = NULL,
                                 response_family = "auto",
                                 replication = NA_integer_,
                                 fold = NA_integer_,
                                 ridge_lambda = 1,
                                 bayes_nIter = 1500L,
                                 bayes_burnIn = 500L,
                                 bayes_thin = 5L,
                                 ntree = 500L,
                                 mtry = NULL,
                                 nodesize = NULL,
                                 rf_n_jobs = 1L) {
  if (is.null(predictor_data) || !ncol(predictor_data)) {
    stop("Fold-internal feature selection requires a non-empty predictor matrix.", call. = FALSE)
  }
  test_rows <- as.integer(test_rows)
  test_rows <- test_rows[is.finite(test_rows) & test_rows >= 1L & test_rows <= nrow(pheno_data)]
  train_rows <- setdiff(seq_len(nrow(pheno_data)), test_rows)
  y <- pheno_data[[response]]
  observed <- if (is.numeric(y) || is.integer(y)) is.finite(as.double(y)) else !is.na(y)
  train_rows <- train_rows[observed[train_rows]]
  training_ids <- unique(as.character(pheno_data[[gen_name]][train_rows]))
  metadata <- gp_feature_score_training_ids(
    predictor_data = predictor_data,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    training_ids = training_ids,
    scoring_model = scoring_model,
    seed = seed,
    source_block = source_block,
    response_family = response_family,
    selection_mode = "fold_internal",
    replication = replication,
    fold = fold,
    ridge_lambda = ridge_lambda,
    bayes_nIter = bayes_nIter,
    bayes_burnIn = bayes_burnIn,
    bayes_thin = bayes_thin,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    rf_n_jobs = rf_n_jobs
  )
  selected <- gp_feature_selected_metadata(metadata, trait = response, k = k)
  keep <- as.character(selected[["predictor"]])
  missing <- setdiff(keep, colnames(predictor_data))
  if (length(missing)) {
    stop(
      "Fold-internal selected predictors are missing from the model matrix: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  metadata[["selected"]] <- metadata[["predictor"]] %in% keep
  metadata[["feature_k"]] <- as.integer(length(keep))
  list(
    predictor_data = predictor_data[, keep, drop = FALSE],
    metadata = metadata,
    selected_predictors = keep,
    selected_by_source = gp_feature_selected_by_source(metadata, response, k),
    feature_k = as.integer(length(keep))
  )
}

gp_feature_multitrait_view <- function(predictor_data,
                                       pheno_data,
                                       response,
                                       gen_name,
                                       k,
                                       selection_mode = "fixed",
                                       feature_score_metadata = NULL,
                                       test_rows = integer(),
                                       scoring_model = "Ridge_Regression",
                                       seed = 123L,
                                       source_block = NULL,
                                       response_family = "gaussian",
                                       replication = NA_integer_,
                                       fold = NA_integer_,
                                       ridge_lambda = 1,
                                       bayes_nIter = 1500L,
                                       bayes_burnIn = 500L,
                                       bayes_thin = 5L,
                                       ntree = 500L,
                                       mtry = NULL,
                                       nodesize = NULL,
                                       rf_n_jobs = 1L) {
  selection_mode <- gp_feature_normalize_cv_policy(selection_mode, allow_both = FALSE)
  if (!selection_mode %in% c("fixed", "fold_internal")) {
    stop("Multi-trait feature selection requires fixed or fold_internal mode.", call. = FALSE)
  }
  response <- as.character(response)
  views <- lapply(seq_along(response), function(i) {
    trait <- response[[i]]
    if (identical(selection_mode, "fixed")) {
      selected <- gp_feature_selected_metadata(feature_score_metadata, trait = trait, k = k)
      if (!nrow(selected)) {
        stop("No fixed feature scores were found for trait '", trait, "'.", call. = FALSE)
      }
      meta <- feature_score_metadata[
        as.character(feature_score_metadata[["trait"]]) == trait,
        ,
        drop = FALSE
      ]
      keep <- as.character(selected[["predictor"]])
      meta[["selected"]] <- meta[["predictor"]] %in% keep
      meta[["feature_k"]] <- as.integer(length(keep))
      meta[["replication"]] <- replication
      meta[["fold"]] <- fold
      meta[["selection_mode"]] <- "fixed"
      list(metadata = meta, selected_predictors = keep)
    } else {
      gp_feature_fold_view(
        predictor_data = predictor_data,
        pheno_data = pheno_data,
        response = trait,
        gen_name = gen_name,
        test_rows = test_rows,
        k = k,
        scoring_model = scoring_model,
        seed = as.integer(seed + i - 1L),
        source_block = source_block,
        response_family = response_family,
        replication = replication,
        fold = fold,
        ridge_lambda = ridge_lambda,
        bayes_nIter = bayes_nIter,
        bayes_burnIn = bayes_burnIn,
        bayes_thin = bayes_thin,
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        rf_n_jobs = rf_n_jobs
      )
    }
  })
  keep <- unique(unlist(lapply(views, `[[`, "selected_predictors"), use.names = FALSE))
  missing <- setdiff(keep, colnames(predictor_data))
  if (length(missing)) {
    stop("Multi-trait selected predictors are missing from the model matrix: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metadata <- gp_feature_bind_metadata(lapply(views, `[[`, "metadata"))
  metadata[["joint_union_selected"]] <- metadata[["predictor"]] %in% keep
  metadata[["joint_union_k"]] <- as.integer(length(keep))
  list(
    predictor_data = predictor_data[, keep, drop = FALSE],
    metadata = metadata,
    selected_predictors = keep,
    selected_by_trait = stats::setNames(
      lapply(views, `[[`, "selected_predictors"),
      response
    ),
    joint_union_k = as.integer(length(keep))
  )
}

gp_feature_apply_multitrait_to_ml_dat_res <- function(ml_dat_res,
                                                      pheno_data,
                                                      response,
                                                      gen_name,
                                                      feature_score_metadata,
                                                      k) {
  if (is.null(ml_dat_res) || is.null(feature_score_metadata) || is.null(k)) {
    return(ml_dat_res)
  }
  train_x <- ml_dat_res[["merged_data"]][["merge_data"]] %||% NULL
  if (is.null(train_x)) return(ml_dat_res)
  view <- gp_feature_multitrait_view(
    predictor_data = train_x,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    k = k,
    selection_mode = "fixed",
    feature_score_metadata = feature_score_metadata
  )
  out <- ml_dat_res
  keep <- view$selected_predictors
  out[["merged_data"]][["merge_data"]] <- view$predictor_data
  if (!is.null(out[["merged_data_test"]])) {
    missing <- setdiff(keep, colnames(out[["merged_data_test"]]))
    if (length(missing)) {
      stop("Selected multi-trait predictors are missing from the test feature block: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    out[["merged_data_test"]] <- out[["merged_data_test"]][, keep, drop = FALSE]
  }
  out[["feature_selected_by_trait"]] <- view$selected_by_trait
  out[["feature_selected"]] <- keep
  out[["feature_k"]] <- as.integer(length(keep))
  out[["feature_selection_metadata"]] <- view$metadata
  out
}

gp_feature_genomic_kernel_from_view <- function(predictor_data,
                                                source_map = NULL,
                                                gmatrix_method = NULL,
                                                ploidy = "auto",
                                                context = "feature-selected genomic model") {
  if (is.null(predictor_data) || !ncol(predictor_data)) {
    stop(context, " has no selected marker columns.", call. = FALSE)
  }
  source_map <- gp_feature_source_block(colnames(predictor_data), source_map)
  sources <- unique(as.character(source_map[colnames(predictor_data)]))
  if (!all(tolower(sources) %in% tolower(gp_feature_source_aliases("geno_data")))) {
    stop(
      context,
      " currently accepts selected genotype markers only; the selected set also contains: ",
      paste(setdiff(sources, "geno_data"), collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(gmatrix_method) || length(gmatrix_method) != 1L) {
    stop(context, " requires exactly one `gmatrix_method`.", call. = FALSE)
  }
  grm_calculation(
    geno_clean = as.matrix(predictor_data),
    method = gmatrix_method,
    ploidy = ploidy
  )
}

gp_feature_joint_kernel_bank <- function(predictor_data,
                                         pheno_data,
                                         response,
                                         gen_name,
                                         feature_score_metadata,
                                         k,
                                         source_map,
                                         source_matrices,
                                         gmatrix_method = NULL,
                                         kernel_method = NULL,
                                         ploidy = "auto",
                                         scaling = TRUE,
                                         centering = FALSE,
                                         context = "joint model") {
  view <- gp_feature_multitrait_view(
    predictor_data = predictor_data,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    k = k,
    selection_mode = "fixed",
    feature_score_metadata = feature_score_metadata
  )
  source_map <- gp_feature_source_block(colnames(predictor_data), source_map)
  selected_sources <- as.character(source_map[view$selected_predictors])
  selected_sources[is.na(selected_sources) | !nzchar(selected_sources)] <- "merged"
  selected_by_source <- split(view$selected_predictors, selected_sources)
  subset_sources <- gp_feature_subset_source_matrices(
    source_matrices = source_matrices,
    selected_by_source = selected_by_source,
    context = context
  )
  bank <- gp_feature_build_kernel_bank(
    selected_sources = subset_sources,
    gmatrix_method = gmatrix_method,
    kernel_method = kernel_method,
    ploidy = ploidy,
    scaling = scaling,
    centering = centering,
    context = context
  )
  list(view = view, kernel_bank = bank)
}

gp_feature_bind_metadata <- function(records) {
  records <- Filter(function(x) is.data.frame(x) && nrow(x), records %||% list())
  if (!length(records)) return(NULL)
  all_names <- unique(unlist(lapply(records, names), use.names = FALSE))
  records <- lapply(records, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

gp_feature_subset_source_matrices <- function(source_matrices,
                                              selected_by_source,
                                              context = "feature selection") {
  source_matrices <- source_matrices %||% list()
  selected_by_source <- selected_by_source %||% list()
  canonical <- c(
    geno_data = "geno_data",
    omic1_data = "omic1_data",
    omic2_data = "omic2_data",
    omic3_data = "omic3_data"
  )
  out <- stats::setNames(vector("list", length(canonical)), names(canonical))
  for (source in names(canonical)) {
    x <- source_matrices[[source]] %||% NULL
    keep <- gp_feature_selected_vector(selected_by_source, source = source)
    if (is.null(x) || !length(keep)) {
      out[[source]] <- NULL
      next
    }
    missing <- setdiff(keep, colnames(x))
    if (length(missing)) {
      stop(
        context, " selected predictors missing from ", source, ": ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    source_ploidy <- attr(x, "ploidy", exact = TRUE)
    out[[source]] <- as.matrix(x[, keep, drop = FALSE])
    storage.mode(out[[source]]) <- "double"
    if (identical(source, "geno_data") && !is.null(source_ploidy)) {
      attr(out[[source]], "ploidy") <- source_ploidy
    }
  }
  out
}

gp_feature_build_kernel_bank <- function(selected_sources,
                                         gmatrix_method = NULL,
                                         kernel_method = NULL,
                                         ploidy = "auto",
                                         scaling = TRUE,
                                         centering = FALSE,
                                         context = "feature selection") {
  selected_sources <- selected_sources %||% list()
  out <- list()
  geno <- selected_sources[["geno_data"]] %||% NULL
  if (!is.null(geno)) {
    if (is.null(gmatrix_method) || !length(gmatrix_method)) {
      stop(
        context,
        " must rebuild the genomic relationship matrix from the selected markers; provide `gmatrix_method`.",
        call. = FALSE
      )
    }
    if (any(as.character(gmatrix_method) == "Weighted_VanRaden")) {
      stop(
        context,
        " does not rebuild Weighted_VanRaden without marker-specific weights. Use a different `gmatrix_method` or supply an external selected feature set and matching kernel.",
        call. = FALSE
      )
    }
    out[["gmatrix_model_ready"]] <- grm_calculation(
      geno_clean = geno,
      method = gmatrix_method,
      ploidy = ploidy
    )
  }
  for (source in c("omic1_data", "omic2_data", "omic3_data")) {
    x <- selected_sources[[source]] %||% NULL
    if (is.null(x)) next
    if (is.null(kernel_method) || !length(kernel_method)) {
      stop(
        context,
        " must rebuild the omics kernel from selected features; provide `kernel_method`.",
        call. = FALSE
      )
    }
    key <- sub("_data$", "_kernel_model_ready", source)
    out[[key]] <- kernel_calculation(
      M_matrix_clean = x,
      scaling = scaling,
      centering = centering,
      method = kernel_method,
      message = FALSE
    )
  }
  out
}

gp_feature_prepare_cv_model_inputs <- function(model,
                                               trait,
                                               pheno_data,
                                               gen_name,
                                               selected_by_source,
                                               source_matrices,
                                               params) {
  model <- as.character(model[[1L]])
  marker_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  kernel_bayes_models <- c("GBLUP_BRR", "RKHS")
  needs_marker_design <- model %in% marker_models
  needs_kernel <- gp_feature_model_needs_kernel(model, met_ml_dl = params$met_ml_dl)
  if (needs_kernel && isTRUE(params$feature_precomputed_kernel_supplied)) {
    stop(
      "Automatic feature selection cannot subset a user-supplied precomputed kernel for ", model,
      ". Provide the corresponding raw marker/omics matrix and omit the precomputed kernel, or provide an externally rebuilt kernel matching `feature_selected`.",
      call. = FALSE
    )
  }
  selected_sources <- gp_feature_subset_source_matrices(
    source_matrices = source_matrices,
    selected_by_source = selected_by_source,
    context = paste0("Trait '", trait, "' ", model)
  )
  if ((needs_marker_design || needs_kernel) && !length(Filter(Negate(is.null), selected_sources))) {
    stop(
      "Automatic feature selection for ", model, " requires raw marker/omics matrices with named columns. Precomputed kernels alone cannot be feature-selected.",
      call. = FALSE
    )
  }
  kernel_bank <- if (needs_kernel) {
    gp_feature_build_kernel_bank(
      selected_sources = selected_sources,
      gmatrix_method = params$gmatrix_method,
      kernel_method = params$kernel_method,
      ploidy = params$ploidy %||% "auto",
      scaling = params$scaling %||% TRUE,
      centering = params$centering %||% FALSE,
      context = paste0("Trait '", trait, "' ", model)
    )
  } else {
    list()
  }
  model_prep <- NULL
  if (model %in% c(marker_models, kernel_bayes_models)) {
    model_prep <- model_prep_bayes_cv(
      fixed = params$fixed,
      random = params$random,
      GS_model_cv = model,
      response = trait,
      gen_name = gen_name,
      pheno_data = pheno_data,
      test_set = NULL,
      weights = params$weights,
      fixed_term_model_bayesian = params$fixed_term_model_bayesian,
      rand_term_model_bayesian = params$rand_term_model_bayesian,
      nIter = params$nIter,
      burnIn = params$burnIn,
      thin = params$thin,
      geno_data = selected_sources[["geno_data"]],
      omic1_data = selected_sources[["omic1_data"]],
      omic2_data = selected_sources[["omic2_data"]],
      omic3_data = selected_sources[["omic3_data"]],
      omics_data_label = params$omics_data_label,
      gmatrix = kernel_bank[["gmatrix_model_ready"]],
      omic1_kernel = kernel_bank[["omic1_kernel_model_ready"]],
      omic2_kernel = kernel_bank[["omic2_kernel_model_ready"]],
      omic3_kernel = kernel_bank[["omic3_kernel_model_ready"]],
      kernel_list = kernel_bank,
      heter_groups = params$heter_groups,
      heter_resid = params$heter_resid %||% FALSE,
      bayes_kernel_heter_resid = params$bayes_kernel_heter_resid,
      omics_kernel_label = params$omics_kernel_label,
      cross_validation = TRUE,
      response_family = params$response_family %||% "gaussian"
    )
  }
  asreml_prep <- NULL
  if (identical(model, "GBLUP")) {
    asreml_prep <- asreml_utilis_new(
      fixed = params$fixed,
      random = params$random,
      engine = params$engine,
      cova = params$cova,
      GS_model = "GBLUP",
      response = trait,
      pheno_data = pheno_data,
      gmatrix = kernel_bank[["gmatrix_model_ready"]],
      omic1_kernel = kernel_bank[["omic1_kernel_model_ready"]],
      omic2_kernel = kernel_bank[["omic2_kernel_model_ready"]],
      omic3_kernel = kernel_bank[["omic3_kernel_model_ready"]],
      kernel_list = kernel_bank,
      inverse = params$inverse,
      epsilon = params$epsilon %||% 1e-6,
      gen_name = gen_name,
      heter_groups = params$heter_groups,
      heter_resid = params$heter_resid %||% FALSE,
      var_cov_str = params$var_cov_str,
      weights = params$weights,
      workspace = params$workspace %||% 1e8,
      pworkspace = params$pworkspace %||% 1e6,
      maxit = params$maxit %||% 50L,
      cross_validation = TRUE
    )
  }
  list(
    selected_sources = selected_sources,
    kernel_bank = kernel_bank,
    model_prep_all_bayes_cv = model_prep,
    asreml_models_prep_cv = asreml_prep
  )
}

gp_feature_rebuild_true_kernel_bank <- function(model,
                                                trait,
                                                selected_by_source,
                                                source_matrices,
                                                gmatrix_method = NULL,
                                                kernel_method = NULL,
                                                ploidy = "auto",
                                                scaling = TRUE,
                                                centering = FALSE,
                                                met_ml_dl = FALSE,
                                                precomputed_kernel_supplied = FALSE) {
  model <- as.character(model[[1L]])
  needs_kernel <- gp_feature_model_needs_kernel(model, met_ml_dl = met_ml_dl)
  if (!needs_kernel) return(NULL)
  if (isTRUE(precomputed_kernel_supplied)) {
    stop(
      "Automatic feature selection cannot subset a user-supplied precomputed kernel for ", model,
      ". Provide raw marker/omics inputs without that kernel, or provide an externally rebuilt matching kernel.",
      call. = FALSE
    )
  }
  selected_sources <- gp_feature_subset_source_matrices(
    source_matrices = source_matrices,
    selected_by_source = selected_by_source,
    context = paste0("Trait '", trait, "' ", model)
  )
  if (!length(Filter(Negate(is.null), selected_sources))) {
    stop(
      "Automatic feature selection for ", model, " requires raw marker/omics matrices with named columns. Precomputed kernels alone cannot be feature-selected.",
      call. = FALSE
    )
  }
  gp_feature_build_kernel_bank(
    selected_sources = selected_sources,
    gmatrix_method = gmatrix_method,
    kernel_method = kernel_method,
    ploidy = ploidy,
    scaling = scaling,
    centering = centering,
    context = paste0("Trait '", trait, "' ", model)
  )
}

gp_feature_scale_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  means <- colMeans(x, na.rm = TRUE)
  for (j in seq_len(ncol(x))) {
    miss <- !is.finite(x[, j])
    if (any(miss)) {
      x[miss, j] <- means[[j]]
    }
  }
  sds <- apply(x, 2, stats::sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  x <- sweep(x, 2, means, "-")
  sweep(x, 2, sds, "/")
}

gp_feature_score_ridge <- function(x, y, lambda = 1) {
  y <- as.double(y)
  if (!all(is.finite(y))) {
    keep <- is.finite(y)
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
  }
  x_scaled <- gp_feature_scale_matrix(x)
  y_scaled <- as.double(scale(y))
  y_scaled[!is.finite(y_scaled)] <- 0
  lambda <- suppressWarnings(as.numeric(lambda %||% 1))
  if (!is.finite(lambda) || lambda <= 0) {
    lambda <- 1
  }
  n <- nrow(x_scaled)
  p <- ncol(x_scaled)
  beta <- if (p > n) {
    gram <- tcrossprod(x_scaled)
    as.vector(crossprod(x_scaled, qr.solve(gram + diag(lambda, n), y_scaled)))
  } else {
    xtx <- crossprod(x_scaled)
    as.vector(qr.solve(xtx + diag(lambda, p), crossprod(x_scaled, y_scaled)))
  }
  names(beta) <- colnames(x_scaled)
  abs(beta)
}

gp_feature_score_bayesb <- function(x,
                                    y,
                                    seed = 123L,
                                    nIter = 1500L,
                                    burnIn = 500L,
                                    thin = 5L) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BayesB feature scoring requires the BGLR package.", call. = FALSE)
  }
  y <- as.double(y)
  keep <- is.finite(y)
  x <- gp_feature_scale_matrix(x[keep, , drop = FALSE])
  y <- y[keep]
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  restore_seed <- seed_exists && is.integer(old_seed) && length(old_seed) > 1L
  # BGLR chain files go to the session temp dir, not the user's working dir
  # (prefix drawn before seeding so the seeded fit is unchanged)
  save_prefix <- gp_bglr_save_prefix(model_name = "feature_BayesB", response = "trait")
  on.exit(unlink(Sys.glob(paste0(save_prefix, "*"))), add = TRUE)
  gp_set_seed(as.integer(seed %||% 123L))
  fit <- BGLR::BGLR(
    y = y,
    ETA = list(markers = list(X = x, model = "BayesB")),
    nIter = as.integer(nIter),
    burnIn = as.integer(burnIn),
    thin = as.integer(thin),
    verbose = FALSE,
    saveAt = save_prefix
  )
  if (restore_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  effects <- fit$ETA$markers$b %||% rep(0, ncol(x))
  names(effects) <- colnames(x)
  abs(as.numeric(effects))
}

gp_feature_score_python_randomforest <- function(x,
                                                 y,
                                                 response_family = "auto",
                                                 seed = 123L,
                                                 ntree = 500L,
                                                 mtry = NULL,
                                                 nodesize = NULL,
                                                 n_jobs = 1L) {
  x <- gp_feature_scale_matrix(x)
  params <- list(
    ntree = as.integer(ntree %||% 500L),
    mtry = mtry,
    nodesize = nodesize,
    n_jobs = as.integer(n_jobs %||% 1L),
    random_state = as.integer(seed %||% 123L)
  )
  fam <- gp_resolve_response_family(response_family, y = y)
  class_levels <- if (identical(fam, "gaussian")) NULL else {
    if (is.factor(y) || is.ordered(y)) as.character(levels(y)) else sort(unique(as.character(y)))
  }
  imp <- gp_py_ml_feature_importance(
    model_type = "randomforest",
    X_train = x,
    y_train = y,
    response_family = fam,
    model_params = params,
    class_levels = class_levels
  )
  imp <- as.numeric(imp)
  imp[!is.finite(imp)] <- 0
  abs(imp)
}

gp_feature_metadata_table <- function(trait,
                                      predictors,
                                      source_block,
                                      raw_score,
                                      scoring_model,
                                      seed,
                                      n_train,
                                      response_family = "gaussian",
                                      selection_mode = "fixed",
                                      replication = NA_integer_,
                                      fold = NA_integer_,
                                      training_id_hash = NA_character_) {
  raw_score <- as.numeric(raw_score)
  raw_score[!is.finite(raw_score)] <- 0
  max_score <- max(raw_score, na.rm = TRUE)
  normalized <- if (is.finite(max_score) && max_score > 0) raw_score / max_score else rep(0, length(raw_score))
  ord <- order(-normalized, predictors)
  ranks <- integer(length(predictors))
  ranks[ord] <- seq_along(ord)
  data.frame(
    trait = trait,
    predictor = predictors,
    source_block = unname(source_block[predictors]),
    rank = ranks,
    raw_score = raw_score,
    normalized_score = normalized,
    scoring_model = scoring_model,
    response_family = response_family,
    selection_mode = selection_mode,
    replication = as.integer(replication),
    fold = as.integer(fold),
    training_id_hash = as.character(training_id_hash),
    seed = as.integer(seed %||% NA_integer_),
    n_train = as.integer(n_train),
    stringsAsFactors = FALSE
  )[ord, , drop = FALSE]
}

gp_feature_normalize_k <- function(k, n_predictors) {
  if (is.null(k) || identical(k, "all")) {
    return(as.integer(n_predictors))
  }
  k_val <- suppressWarnings(as.integer(k[[1L]]))
  if (!is.finite(k_val) || k_val < 1L) {
    stop("feature k must be a positive integer or 'all'.", call. = FALSE)
  }
  as.integer(min(k_val, n_predictors))
}

gp_feature_k_grid <- function(k_grid, n_predictors, include_all = TRUE) {
  if (is.null(k_grid)) {
    vals <- unique(c(10L, 50L, 100L, 250L, 500L, as.integer(n_predictors)))
  } else {
    vals <- unlist(k_grid, use.names = FALSE)
    vals <- ifelse(as.character(vals) == "all", n_predictors, vals)
    vals <- suppressWarnings(as.integer(vals))
    # Always compare against the full predictor set, so CV can conclude that
    # selection does not help (e.g. BayesB already shrinks most markers).
    # Specialized hybrid/joint multi-trait routes need exactly one k per run.
    if (isTRUE(include_all)) {
      vals <- c(vals, as.integer(n_predictors))
    }
  }
  vals <- vals[is.finite(vals) & vals > 0L]
  vals <- unique(as.integer(pmin(vals, n_predictors)))
  vals[order(vals)]
}

gp_feature_apply_to_ml_dat_res <- function(ml_dat_res,
                                           feature_score_metadata,
                                           trait,
                                           k) {
  if (is.null(ml_dat_res) || is.null(feature_score_metadata) || is.null(k)) {
    return(ml_dat_res)
  }
  out <- ml_dat_res
  train_x <- out[["merged_data"]][["merge_data"]]
  selected_train <- feature_select_top_k(train_x, feature_score_metadata, trait = trait, k = k)
  keep <- colnames(selected_train)
  out[["merged_data"]][["merge_data"]] <- selected_train
  if (!is.null(out[["merged_data_test"]])) {
    missing <- setdiff(keep, colnames(out[["merged_data_test"]]))
    if (length(missing)) {
      stop(
        "Selected predictors are missing from the test feature block: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    out[["merged_data_test"]] <- out[["merged_data_test"]][, keep, drop = FALSE]
  }
  trait_meta <- feature_score_metadata[
    as.character(feature_score_metadata[["trait"]]) == as.character(trait),
    ,
    drop = FALSE
  ]
  trait_meta[["selected"]] <- trait_meta[["predictor"]] %in% keep
  trait_meta[["feature_k"]] <- as.integer(length(keep))
  trait_meta[["selection_mode"]] <- "fixed"
  out[["feature_selected"]] <- keep
  out[["feature_k"]] <- as.integer(length(keep))
  out[["feature_selection_metadata"]] <- trait_meta
  out
}

gp_feature_select_matrix_by_metadata <- function(x,
                                                 feature_score_metadata,
                                                 trait,
                                                 k,
                                                 source = NULL) {
  if (is.null(x) || is.null(feature_score_metadata) || is.null(k) || is.null(colnames(x))) {
    return(x)
  }
  meta <- feature_score_metadata[as.character(feature_score_metadata$trait) == as.character(trait), , drop = FALSE]
  if (!nrow(meta)) {
    return(x)
  }
  meta <- meta[order(meta$rank, meta$predictor), , drop = FALSE]
  keep <- as.character(meta$predictor[seq_len(gp_feature_normalize_k(k, nrow(meta)))])
  if (!is.null(source) && "source_block" %in% names(meta)) {
    source_aliases <- tolower(gp_feature_source_aliases(source))
    selected_rows <- meta[["predictor"]] %in% keep &
      tolower(as.character(meta[["source_block"]])) %in% source_aliases
    keep <- as.character(meta[["predictor"]][selected_rows])
    if (!length(keep)) return(NULL)
  }
  keep <- keep[keep %in% colnames(x)]
  if (!length(keep)) {
    return(x)
  }
  x[, keep, drop = FALSE]
}
