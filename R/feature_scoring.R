#' Score and rank predictors for dynamic top-k model fitting
#'
#' @param predictor_data Numeric predictor matrix or data frame with individual
#'   IDs as row names and unique predictor names as columns.
#' @param pheno_data Phenotype data containing `gen_name` and response columns.
#' @param response Character vector of response variables to score.
#' @param gen_name Column name containing individual IDs.
#' @param scoring_model One of `"Ridge_Regression"`, `"BayesB"` or
#'   `"RandomForest"`. With `joint = TRUE` these use BGLR Multitrait BRR,
#'   BGLR Multitrait SpikeSlab, or Python multi-output RandomForest AIR.
#' @param seed Nonnegative integer seed for stochastic scoring.
#' @param source_block Optional source label or named vector keyed by predictor.
#' @param ridge_lambda Positive penalty for single-trait closed-form ridge.
#'   Ignored by joint BRR, which estimates marker covariance and shrinkage.
#' @param bayes_nIter,bayes_burnIn,bayes_thin MCMC controls for BayesB and
#'   joint BRR/SpikeSlab scoring. Defaults are 1500, 500 and 5; they do not
#'   establish convergence or stable cut-off ranks. Compare seeds and diagnose
#'   mixing; longer chains may help. Thinning does not guarantee independence.
#' @param ntree Number of trees for RandomForest importance scoring.
#' @param mtry Number of predictors sampled at each split.
#' @param nodesize Minimum node size.
#' @param rf_n_jobs Internal Python RandomForest jobs; use 1 when the caller
#'   parallelizes over traits, models, folds or k values.
#' @param response_family Response family. `"auto"` infers it per trait.
#'   Ridge and BayesB require Gaussian responses. Separate-trait RandomForest
#'   also supports binary, ordinal and multiclass responses. Joint scoring
#'   requires at least two distinct numeric Gaussian traits.
#' @param selection_mode Ranking provenance, normally `"fixed"` or
#'   `"fold_internal"`.
#' @param replication,fold Optional resampling identifiers.
#' @param fixed Optional design formula for separate-trait scoring, such as
#'   `~ Env + Rep`. Missing covariates are dropped and counted. Terms base R
#'   cannot evaluate warn and are omitted; matrix-valued terms are supported.
#'   Joint scoring accepts only NULL or `~1` and sets its own trait intercepts
#'   and environment effects from `heter_groups`.
#' @param heter_groups Optional environment column used to adjust environment
#'   means among the visible records.
#' @param joint Logical, default FALSE. TRUE estimates one multi-trait ranking
#'   for shared-set selection; FALSE preserves separate trait rankings.
#'
#' @return A `predictpror_feature_scores` data frame with one row per predictor
#'   for joint scoring, or per trait/predictor otherwise. `scoring_scope` is
#'   `"joint"` or `"trait"`; joint trait labels join names with `"+"`, which is
#'   reserved in response names. Ranks follow `raw_score` with locale-independent
#'   name-based ties. RandomForest ranks zero-variation predictors last;
#'   `normalized_score` is raw/max(abs(raw)). Metadata includes source, engine,
#'   seed, record counts, adjusted terms, ranking source and training-record hash.
#'   Joint metadata adds `min_pair_coobserved`, chain controls, `effect_sd`,
#'   `scoring_status` and covariance-prior regularization fields. The hash
#'   identifies retained records, not the visible trait-cell pattern. Scores
#'   are effect-size or predictive rankings, not variance components or causal
#'   effects. A top-k cut-off tie is recorded and warns once per run.
#' @details
#' Separate-trait ridge adjusts y and X by the Frisch-Waugh-Lovell theorem,
#' giving ridge with unpenalized design effects. Bayesian FIXED terms instead
#' have diffuse finite-variance Gaussian priors. Forests residualize Gaussian
#' responses, or use design covariates as unranked classification predictors.
#'
#' Joint Bayesian scoring standardizes each trait by its visible, design-adjusted
#' phenotypic SD and each marker by its visible adjusted SD. Its raw score is
#' sqrt(sum_t posterior_mean(beta_jt)^2), the norm of the posterior means,
#' rather than the posterior mean of a norm. It is invariant to changing trait
#' units, but not to general linear reparameterization or duplicating correlated
#' traits. Higher-heritability traits may contribute more under this scaling;
#' this is not guaranteed marker by marker. Under LD, genetic variance contains
#' cross-marker covariance terms: squared scores do not add to heritability.
#' SpikeSlab uses trait-specific inclusion indicators with correlated Gaussian
#' slabs (probIn = 0.5, counts = 2); it does not enforce group inclusion.
#' This is a different prior/model from the separate-trait BGLR BayesB scorer;
#' their common public label does not make the methods matched controls.
#'
#' Marker and residual covariances are unstructured with explicit inverse-Wishart
#' priors calibrated from visible, scaled, design-adjusted responses. A pairwise
#' calibration correlation with a minimum eigenvalue below 1e-6 is shrunk toward
#' identity to that floor, with a warning and metadata. The marker prior uses
#' df0 = T + 1 and half the calibrated covariance divided by the marker MSx,
#' times df0 + T + 1; residual df0 = max(5,T) uses half the covariance times
#' df0 + T + 1. This calibrates inverse-Wishart modes, not prior means
#' (which need not exist) or a fixed posterior heritability. Perfectly or
#' numerically singular mean-imputed response startup covariance stops.
#' A zero standardized marker column has centered prior mean zero: it receives
#' score zero, `effect_sd = NA` and `scoring_status = "no_predictor_variation"`.
#' All-missing predictor columns and panels with no informative predictors stop.
#'
#' Bayesian `effect_sd` is sqrt(sum_t posterior_variance(beta_jt)), equivalently
#' the square root of the trace of coefficient posterior covariance. It is not
#' the posterior SD or a credible interval of the nonlinear score, a Monte Carlo
#' standard error, or a convergence diagnostic. Those require additional draws
#' and dependence information. RandomForest `effect_sd` is NA.
#'
#' BGLR retains records with any visible response and samples missing traits.
#' Ignoring missingness requires MAR given the modeled observed information,
#' distinct missingness parameters and appropriate prior factorization, as well
#' as a suitable outcome model; MAR alone does not guarantee unbiased estimates.
#' With no co-observed trait pair there is no direct residual cross-covariance
#' likelihood contribution, although coupled covariance priors may affect the
#' posterior. Records are independent conditional on markers and design effects;
#' scoring has no additional line effect or GxE term. Repeated records do not
#' create independent genetic information. Joint RandomForest keeps complete
#' response records only, warns about incomplete records and targets that retained
#' subset unless its selection is representative.
#'
#' RandomForest AIR subtracts the impurity importance of a shadow predictor,
#' with one shared row permutation for all shadows (Nembrini et al. 2018).
#' Joint scoring fits the scaled, design-residualized response matrix. sklearn's
#' squared-error criterion averages across outputs; switching to their sum is a
#' common factor and leaves split choices and ranking unchanged. AIR can be
#' negative. It reduces impurity bias but does not remove LD proxy competition,
#' establish causal attribution or provide a valid p-value. The scoring forest
#' is fitted separately from any downstream prediction forest. Zero-variation
#' predictors retain raw score zero and status `no_predictor_variation`, but rank
#' after all informative predictors, including those with negative AIR. Cutoff
#' ties compare both this eligibility and the raw score. Supplied tied ranks are
#' preserved across locales; newly generated ties use radix predictor-name order.
#'
#' Selecting with the same phenotypes used to fit a model can inflate variance
#' and heritability estimates and make uncertainty optimistic; `model_execute()`
#' labels selection-conditional outputs. Fold-internal scoring must receive
#' only visible outcome cells. See `system.file("feature-selection-contract.md",
#' package = "PredictProR")` for route-specific masks, shared-set rules and limits.
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
                                     fold = NA_integer_,
                                     fixed = NULL,
                                     heter_groups = NULL,
                                     joint = FALSE) {
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

  if (!is.logical(joint) || length(joint) != 1L || is.na(joint)) {
    stop("joint must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(joint)) {
    return(gp_feature_score_joint(
      x, pheno_data, response, gen_name, scoring_model, seed, source_block,
      bayes_nIter, bayes_burnIn, bayes_thin, ntree, mtry, nodesize, rf_n_jobs,
      response_family, selection_mode, replication, fold, fixed, heter_groups
    ))
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
    rows <- which(valid)
    if (!length(rows)) {
      stop("Trait ", trait, " has no observed training values for feature scoring.", call. = FALSE)
    }
    design <- gp_feature_covariate_design(pheno_data[rows, , drop = FALSE],
                                          fixed = fixed, heter_groups = heter_groups)
    rows <- rows[design$keep]
    if (!length(rows)) {
      stop("Trait ", trait, " has no records with complete covariates for feature scoring.",
           call. = FALSE)
    }
    covariates <- if (isTRUE(design$trivial)) NULL else design$C
    ids <- pheno_ids[rows]
    x_train <- x[ids, , drop = FALSE]
    y_train <- y_all[rows]
    scores <- switch(
      scoring_model,
      Ridge_Regression = gp_feature_score_ridge(x_train, y_train, lambda = ridge_lambda,
                                                covariates = covariates),
      BayesB = gp_feature_score_bayesb(
        x_train,
        y_train,
        seed = seed,
        nIter = bayes_nIter,
        burnIn = bayes_burnIn,
        thin = bayes_thin,
        covariates = covariates
      ),
      RandomForest = gp_feature_score_python_randomforest(
        x_train,
        y_train,
        response_family = trait_family,
        seed = seed,
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        n_jobs = rf_n_jobs,
        covariates = covariates
      )
    )
    record_keys <- if (!is.null(heter_groups) && length(heter_groups) == 1L &&
                       heter_groups %in% names(pheno_data)) {
      paste(ids, as.character(pheno_data[[heter_groups]][rows]), sep = "\r")
    } else {
      ids
    }
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
      training_id_hash = gp_feature_training_id_hash(record_keys),
      scoring_engine = gp_feature_scoring_engine(scoring_model),
      covariates_adjusted = paste(design$terms, collapse = "+"),
      n_records_dropped = design$n_dropped,
      extra = if (identical(scoring_model, "RandomForest")) {
        list(scoring_status = attr(scores, "scoring_status"))
      } else NULL
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
  top <- gp_feature_top_k_rows(meta, k)
  keep <- as.character(top$predictor)
  keep <- keep[keep %in% colnames(predictor_data)]
  out <- predictor_data[, keep, drop = FALSE]
  attr(out, "tied_at_cutoff") <- isTRUE(attr(top, "tied_at_cutoff"))
  out
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
  # Matrix subsetting drops custom metadata; observed maxima cannot recover
  # crop ploidy, and the dosage coding contract is independent of marker count.
  for (name in c("ploidy", "genotype_coding", "genotype_coding_contract")) {
    value <- attr(x, name, exact = TRUE)
    if (!is.null(value)) attr(out, name) <- value
  }
  attr(out, "predictpror_feature_selected") <- keep
  out
}

gp_feature_apply_raw_selection_to_inputs <- function(ctx) {
  if (!is.list(ctx) || is.null(ctx$feature_selected)) {
    return(ctx)
  }
  joint_route <- any(vapply(ctx[c("multi_trait_bayes", "multi_trait_gp",
    "multi_trait_asreml", "multi_trait_ml", "multi_trait_dl")], isTRUE, logical(1L)))
  source_fields <- list(
    geno_data = c("geno_data", "train_geno_data", "test_geno_data", "female_geno_data", "male_geno_data"),
    omic1_data = c("omic1_data", "train_omic1_data", "test_omic1_data"),
    omic2_data = c("omic2_data", "train_omic2_data", "test_omic2_data"),
    omic3_data = c("omic3_data", "train_omic3_data", "test_omic3_data")
  )
  all_fields <- unique(unlist(source_fields, use.names = FALSE))
  had_raw_predictors <- any(vapply(all_fields, function(field) !is.null(ctx[[field]]), logical(1L)))
  if (!had_raw_predictors) {
    if (joint_route) {
      declarations <- if (gp_feature_has_source_scope(ctx$feature_selected)) {
        stats::setNames(lapply(names(source_fields), function(source) {
          gp_feature_selected_vector(ctx$feature_selected, source = source)
        }), names(source_fields))
      } else list(merged = gp_feature_selected_vector(ctx$feature_selected))
      if (any(vapply(declarations, function(z) is.numeric(z) || is.logical(z), logical(1L)))) {
        stop("Positional joint feature_selected requires raw predictor columns; ",
             "a precomputed kernel can only carry declared predictor names.", call. = FALSE)
      }
      declarations <- lapply(declarations, function(z) unique(as.character(z)))
      declared <- unlist(declarations, use.names = FALSE)
      if (anyNA(declared) || any(!nzchar(declared)) || anyDuplicated(declared)) {
        stop("Joint selection requires unique, nonempty predictor names across declared sources.", call. = FALSE)
      }
      ctx$feature_selected_source_map <- stats::setNames(
        rep(names(declarations), lengths(declarations)), declared)
    }
    # These source labels are also user declarations: kernel rows/columns cannot
    # verify or resolve the originating selected predictors.
    return(ctx)
  }
  source_columns <- lapply(source_fields, function(fields) {
    unique(unlist(lapply(fields, function(field) {
      x <- ctx[[field]]
      if (is.list(x) && !is.data.frame(x) && "snps_matrix" %in% names(x)) {
        x <- x[["snps_matrix"]]
      }
      cols <- colnames(x) %||% character()
      cols <- cols[!cols %in% as.character(ctx$gen_name %||% character())]
      if (joint_route && (anyDuplicated(cols) || anyNA(cols) || any(!nzchar(cols)))) {
        stop("Joint selection requires unique, nonempty predictor names.", call. = FALSE)
      }
      cols
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
  if (joint_route) {
    # Resolve source-specific positions before any raw columns disappear.
    resolved <- lapply(names(source_columns), function(source) {
      declared <- gp_feature_selected_vector(ctx$feature_selected, source = source)
      cols <- source_columns[[source]]
      if (is.null(declared) || !length(declared)) return(character())
      if (!length(cols) && !gp_feature_has_source_scope(ctx$feature_selected)) return(character())
      if (is.numeric(declared)) {
        if (any(!is.finite(declared)) || any(declared != floor(declared))) {
          stop("Joint feature_selected indices must be finite whole numbers.", call. = FALSE)
        }
        if (any(declared < 1 | declared > length(cols))) {
          stop("feature_selected contains predictor indices outside the available columns in ",
               source, ".", call. = FALSE)
        }
        return(unique(cols[as.integer(declared)]))
      }
      if (is.logical(declared)) {
        if (length(declared) != length(cols)) {
          stop("feature_selected logical vector must have one value per predictor column in ",
               source, ".", call. = FALSE)
        }
        declared[is.na(declared)] <- FALSE
        return(cols[declared])
      }
      unique(as.character(declared)[!is.na(declared) & as.character(declared) %in% cols])
    })
    names(resolved) <- names(source_columns)
    keep <- unlist(resolved, use.names = FALSE)
    if (anyDuplicated(keep)) {
      stop("Joint selection requires unique predictor names across retained sources.", call. = FALSE)
    }
    ctx$feature_selected_resolved <- keep
    ctx$feature_selected_source_map <- stats::setNames(
      rep(names(resolved), lengths(resolved)), keep)
    ctx$feature_selected <- resolved[lengths(resolved) > 0L]
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

# feature_scoring_cv = "auto" (spec B6): fold-internal when the package scores the features (no
# leakage of held-out phenotypes into selection), fixed when the user supplies an external ranking
# the package cannot reproduce per fold.
gp_feature_resolve_cv_policy <- function(policy = "auto", ranking_supplied = FALSE) {
  key <- tolower(trimws(as.character(policy %||% "auto")[[1L]]))
  if (identical(key, "auto")) {
    return(if (isTRUE(ranking_supplied)) "fixed" else "fold_internal")
  }
  gp_feature_normalize_cv_policy(key, allow_both = TRUE)
}

gp_feature_normalize_cv_policy <- function(policy = "fixed", allow_both = TRUE) {
  policy <- tolower(trimws(as.character(policy %||% "fixed")[[1L]]))
  if (identical(policy, "auto")) {
    stop("feature_scoring_cv = 'auto' must be resolved with gp_feature_resolve_cv_policy() before use.",
         call. = FALSE)
  }
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
  if (length(k_values) == 1L && is.character(k_values) &&
      identical(tolower(trimws(k_values)), "all")) {
    k_values <- "all"
  } else {
    valid <- is.numeric(k_values) && length(k_values) == 1L &&
      !is.na(k_values) && is.finite(k_values) && k_values >= 1 &&
      k_values == floor(k_values) && k_values <= .Machine$integer.max
    if (!valid) {
      stop(context, " requires one pre-specified positive whole feature_k or 'all'; multiple or invalid k values require a separate run.",
           call. = FALSE)
    }
    k_values <- as.integer(k_values)
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

# Scores exactly the given phenotype records (a fold's training part). Records, not line IDs: a
# line with held-out records in the same fold (multi-environment CV0/CV2, repeated records) must
# not have them scored.
gp_feature_score_training_records <- function(predictor_data,
                                              pheno_rows,
                                              response,
                                              gen_name,
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
                                              rf_n_jobs = 1L,
                                              fixed = NULL,
                                              heter_groups = NULL) {
  if (!NROW(pheno_rows)) {
    stop("Fold-internal feature scoring received no training phenotypes.", call. = FALSE)
  }
  feature_score_predictors(
    predictor_data = predictor_data,
    pheno_data = pheno_rows,
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
    fold = fold,
    fixed = fixed,
    heter_groups = heter_groups
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
  gp_feature_top_k_rows(meta, k)
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
                                 rf_n_jobs = 1L,
                                 fixed = NULL,
                                 heter_groups = NULL,
                                 train_rows = NULL,
                                 predictor_hash = NULL) {
  # predictor_hash: a content hash of predictor_data, or of a matrix predictor_data is a column
  # subset of by name (the cache key adds predictor_data's dimnames). Given, the fold's ranking is
  # memoised in this process, so the CV tasks (models x k) sharing a fold score it once.
  if (is.null(predictor_data) || !ncol(predictor_data)) {
    stop("Fold-internal feature selection requires a non-empty predictor matrix.", call. = FALSE)
  }
  test_rows <- as.integer(test_rows)
  test_rows <- test_rows[is.finite(test_rows) & test_rows >= 1L & test_rows <= nrow(pheno_data)]
  # Spec B1: score exactly the fold's visible training records. Without train_rows those are all
  # rows but the test rows; a scenario that keeps further rows out of training (hybrid
  # Both_New_Parents: the test parents' half-sib hybrids) passes its own training rows.
  train_rows <- if (is.null(train_rows)) {
    setdiff(seq_len(nrow(pheno_data)), test_rows)
  } else {
    train_rows <- as.integer(train_rows)
    sort(unique(train_rows[is.finite(train_rows) & train_rows >= 1L & train_rows <= nrow(pheno_data)]))
  }
  y <- pheno_data[[response]]
  observed <- if (is.numeric(y) || is.integer(y)) is.finite(as.double(y)) else !is.na(y)
  train_rows <- train_rows[observed[train_rows]]
  pheno_rows <- pheno_data[train_rows, , drop = FALSE]
  settings <- list(
    response = response, gen_name = gen_name, scoring_model = scoring_model, seed = seed,
    source_block = source_block, response_family = response_family, replication = replication,
    fold = fold, ridge_lambda = ridge_lambda, bayes_nIter = bayes_nIter,
    bayes_burnIn = bayes_burnIn, bayes_thin = bayes_thin, ntree = ntree, mtry = mtry,
    nodesize = nodesize, rf_n_jobs = rf_n_jobs, heter_groups = heter_groups
  )
  score <- function() {
    gp_feature_score_training_records(
      predictor_data = predictor_data,
      pheno_rows = pheno_rows,
      response = response,
      gen_name = gen_name,
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
      rf_n_jobs = rf_n_jobs,
      fixed = fixed,
      heter_groups = heter_groups
    )
  }
  metadata <- if (is.null(predictor_hash)) {
    score()
  } else {
    key <- gp_feature_fold_cache_key(predictor_hash, predictor_data, pheno_rows, settings, fixed)
    cached <- gp_feature_fold_cache_get(key)
    if (is.null(cached)) {
      cached <- score()
      gp_feature_fold_cache_set(key, cached)
    }
    cached
  }
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
  metadata[["tied_at_cutoff"]] <- isTRUE(attr(selected, "tied_at_cutoff"))
  list(
    predictor_data = predictor_data[, keep, drop = FALSE],
    metadata = metadata,
    selected_predictors = keep,
    selected_by_source = gp_feature_selected_by_source(metadata, response, k),
    feature_k = as.integer(length(keep))
  )
}

gp_feature_trait_view <- function(predictor_data,
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
                                       rf_n_jobs = 1L,
                                       fixed = NULL,
                                       heter_groups = NULL,
                                       train_rows = NULL) {
  selection_mode <- gp_feature_normalize_cv_policy(selection_mode, allow_both = FALSE)
  response <- as.character(response)
  if (length(response) != 1L) {
    stop("The trait view requires one response; joint models require the joint view.", call. = FALSE)
  }
  if (!selection_mode %in% c("fixed", "fold_internal")) {
    stop("Trait feature selection requires fixed or fold_internal mode.", call. = FALSE)
  }
  trait <- response[[1L]]
  if (identical(selection_mode, "fold_internal")) {
    return(gp_feature_fold_view(
      predictor_data = predictor_data, pheno_data = pheno_data,
      response = trait, gen_name = gen_name, test_rows = test_rows, k = k,
      scoring_model = scoring_model, seed = seed, source_block = source_block,
      response_family = response_family, replication = replication, fold = fold,
      ridge_lambda = ridge_lambda, bayes_nIter = bayes_nIter,
      bayes_burnIn = bayes_burnIn, bayes_thin = bayes_thin,
      ntree = ntree, mtry = mtry, nodesize = nodesize, rf_n_jobs = rf_n_jobs,
      fixed = fixed, heter_groups = heter_groups, train_rows = train_rows
    ))
  }
  selected <- gp_feature_selected_metadata(feature_score_metadata, trait = trait, k = k)
  if (!nrow(selected)) {
    stop("No fixed feature scores were found for trait '", trait, "'.", call. = FALSE)
  }
  meta <- feature_score_metadata[as.character(feature_score_metadata$trait) == trait, , drop = FALSE]
  keep <- as.character(selected$predictor)
  missing <- setdiff(keep, colnames(predictor_data))
  if (length(missing)) {
    stop("Selected predictors missing from the model matrix: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  meta$selected <- meta$predictor %in% keep
  meta$feature_k <- as.integer(length(keep))
  meta$tied_at_cutoff <- isTRUE(attr(selected, "tied_at_cutoff"))
  meta$replication <- replication
  meta$fold <- fold
  meta$selection_mode <- "fixed"
  list(predictor_data = predictor_data[, keep, drop = FALSE], metadata = meta,
       selected_predictors = keep, feature_k = as.integer(length(keep)))
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
  view <- gp_feature_joint_view(
    predictor_data = train_x,
    pheno_visible = pheno_data,
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
  out[["feature_selected_by_trait"]] <- NULL
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
  view <- gp_feature_joint_view(
    predictor_data = predictor_data,
    pheno_visible = pheno_data,
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

# Covariate design C for scoring (spec B2): model.matrix(fixed) plus environment indicators (when
# heter_groups is given and not already in fixed), full column rank, over the scored records.
# Records with a missing covariate are dropped (keep = FALSE). trivial = TRUE when C is the
# intercept only, in which case scoring is exactly as without adjustment. A term base R cannot
# evaluate (an ASReml special such as at(), or poly() of a column with missing values) is left out
# of C with a warning. terms lists the terms C keeps after the rank reduction.
gp_feature_covariate_design <- function(pheno_rows, fixed = NULL, heter_groups = NULL) {
  n <- nrow(pheno_rows)
  rhs <- character()
  vars <- character()
  mf <- NULL
  term_vars <- NULL
  fenv <- if (inherits(fixed, "formula")) environment(fixed) %||% globalenv() else globalenv()
  if (inherits(fixed, "formula")) {
    tt <- stats::delete.response(stats::terms(fixed))
    rhs <- attr(tt, "term.labels")
    missing_cols <- setdiff(all.vars(tt), names(pheno_rows))
    if (length(missing_cols)) {
      stop("Feature scoring cannot adjust for `fixed`: column(s) not in the phenotype data: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    rhs <- gp_feature_evaluable_terms(rhs, pheno_rows, fenv)
    if (length(rhs)) {
      tt <- stats::terms(stats::reformulate(rhs, env = fenv))
      vars <- all.vars(tt)
      # Evaluated variables (so factor(EnvCode), log(x) and the like are checked as they enter C)
      mf <- stats::model.frame(tt, data = pheno_rows, na.action = stats::na.pass)
      term_vars <- attr(tt, "factors")
    }
  }
  env <- NULL
  if (!is.null(heter_groups) && length(heter_groups) == 1L && nzchar(heter_groups) &&
      heter_groups %in% names(pheno_rows) && !heter_groups %in% vars) {
    env <- heter_groups
  }
  keep <- rep(TRUE, n)
  if (!is.null(mf)) {
    for (v in names(mf)) {
      col <- mf[[v]]
      complete <- if (is.numeric(col)) is.finite(col) else !is.na(col)
      # a matrix-valued term (poly(), ns()) is complete when its whole row is
      if (is.matrix(complete)) complete <- rowSums(!complete) == 0L
      keep <- keep & complete
    }
  }
  if (!is.null(env)) keep <- keep & !is.na(pheno_rows[[env]])
  dat <- pheno_rows[keep, , drop = FALSE]
  if (!is.null(mf)) {
    mf <- mf[keep, , drop = FALSE]
    single <- names(mf)[vapply(mf, function(col) {
      (is.factor(col) || is.character(col) || is.logical(col)) &&
        length(unique(as.character(col))) < 2L
    }, logical(1))]
    if (length(single)) {
      drop_term <- colSums(term_vars[single, , drop = FALSE]) > 0
      rhs <- rhs[!drop_term[rhs]]
    }
    # a term can evaluate on all scored records and still fail on the complete ones (poly() of a
    # covariate left with too few distinct values)
    rhs <- gp_feature_evaluable_terms(rhs, dat, fenv, build = TRUE)
  }
  C <- matrix(1, nrow(dat), 1L, dimnames = list(NULL, "(Intercept)"))
  col_terms <- "(Intercept)"
  if (length(rhs)) {
    mm_formula <- stats::reformulate(rhs, env = fenv)
    C <- stats::model.matrix(mm_formula, data = dat)
    col_terms <- c("(Intercept)", attr(stats::terms(mm_formula), "term.labels"))[attr(C, "assign") + 1L]
  }
  if (!is.null(env)) {
    e <- droplevels(factor(as.character(dat[[env]])))
    if (nlevels(e) >= 2L) {
      env_cols <- stats::model.matrix(~ e)[, -1L, drop = FALSE]
      C <- cbind(C, env_cols)
      col_terms <- c(col_terms, rep(env, ncol(env_cols)))
    }
  }
  qr_c <- qr(C)
  kept_cols <- qr_c$pivot[seq_len(qr_c$rank)]
  C <- C[, kept_cols, drop = FALSE]
  # the terms with a column left in C (a constant or aliased term adjusts for nothing)
  terms_used <- intersect(c(rhs, env), col_terms[kept_cols])
  list(C = C, keep = keep, terms = terms_used, n_dropped = sum(!keep), trivial = ncol(C) <= 1L)
}

# Terms of `fixed` base R can evaluate on the scored records. The prediction routes accept terms
# such as ASReml's at() that only their engine understands; scoring leaves those out of C, with
# one warning per term and run, rather than failing. build = TRUE checks the design columns
# (model.matrix) on records that are already complete.
gp_feature_evaluable_terms <- function(rhs, pheno_rows, fenv = globalenv(), build = FALSE) {
  if (!length(rhs)) return(rhs)
  ok <- vapply(rhs, function(lbl) {
    tryCatch({
      f <- stats::reformulate(lbl, env = fenv)
      if (isTRUE(build)) {
        stats::model.matrix(f, data = pheno_rows)
      } else {
        stats::model.frame(f, data = pheno_rows, na.action = stats::na.pass)
      }
      TRUE
    }, error = function(e) FALSE)
  }, logical(1L))
  bad <- rhs[!ok]
  new_bad <- setdiff(bad, gp_feature_state$fixed_terms_warned)
  if (length(new_bad)) {
    gp_feature_state$fixed_terms_warned <- c(gp_feature_state$fixed_terms_warned, new_bad)
    warning(sprintf(paste(
      "Feature scoring cannot evaluate the fixed term(s) %s in base R on the scored records (for",
      "example an ASReml special such as at()); scoring does not adjust for them. The prediction",
      "model still uses them."), paste(new_bad, collapse = ", ")), call. = FALSE)
  }
  rhs[ok]
}

gp_feature_residualize <- function(m, covariates) {
  q <- qr(covariates)
  m - qr.fitted(q, m)
}

# No information left after removing C: no variance to begin with (or not finite), or residual SD
# at rounding size relative to the raw SD. Relative, not absolute: a QTL measured in units of 1e-9
# keeps its information.
gp_feature_flat_after_adjustment <- function(sd_res, sd_raw) {
  !is.finite(sd_raw) | sd_raw == 0 | !is.finite(sd_res) |
    sd_res <= sqrt(.Machine$double.eps) * sd_raw
}

gp_feature_scoring_engine <- function(scoring_model) {
  switch(as.character(scoring_model),
         Ridge_Regression = "ridge (closed form)",
         BayesB = "BGLR BayesB",
         RandomForest = "sklearn RandomForest AIR",
         NA_character_)
}

gp_feature_scale_matrix <- function(x, covariates = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  means <- colMeans(x, na.rm = TRUE)
  for (j in seq_len(ncol(x))) {
    miss <- !is.finite(x[, j])
    if (any(miss)) {
      x[miss, j] <- means[[j]]
    }
  }
  if (!is.null(covariates) && ncol(covariates) > 1L) {
    # Frisch-Waugh-Lovell: residualizing on C makes ridge with unpenalized C (spec B2)
    sd_raw <- apply(x, 2, stats::sd)
    x <- gp_feature_residualize(x, covariates)
    sds <- apply(x, 2, stats::sd)
    # A column inside span(C) has residuals of rounding size, not exactly 0: a tolerance relative
    # to the column's own scale (omics units can be tiny) sets it to exactly 0, so it scores 0 as
    # in the unadjusted path. A column with no variance at all is flat whatever its residuals.
    flat <- gp_feature_flat_after_adjustment(sds, sd_raw)
    x[, flat] <- 0
    sds[flat] <- 1
    return(sweep(x, 2, sds, "/"))
  }
  sds <- apply(x, 2, stats::sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  x <- sweep(x, 2, means, "-")
  sweep(x, 2, sds, "/")
}

gp_feature_score_ridge <- function(x, y, lambda = 1, covariates = NULL) {
  y <- as.double(y)
  if (!all(is.finite(y))) {
    keep <- is.finite(y)
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(covariates)) covariates <- covariates[keep, , drop = FALSE]
  }
  if (!is.null(covariates) && ncol(covariates) > 1L) {
    x_scaled <- gp_feature_scale_matrix(x, covariates = covariates)
    y_res <- gp_feature_residualize(y, covariates)
    y_sd <- stats::sd(y_res)
    y_flat <- gp_feature_flat_after_adjustment(y_sd, stats::sd(y))
    y_scaled <- if (y_flat) rep(0, length(y_res)) else y_res / y_sd
  } else {
    x_scaled <- gp_feature_scale_matrix(x)
    y_scaled <- as.double(scale(y))
  }
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
                                    thin = 5L,
                                    covariates = NULL) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BayesB feature scoring requires the BGLR package.", call. = FALSE)
  }
  y <- as.double(y)
  keep <- is.finite(y)
  if (!is.null(covariates)) covariates <- covariates[keep, , drop = FALSE]
  adjusted <- !is.null(covariates) && ncol(covariates) > 1L
  # Under adjustment the markers are standardized by their SD within C, as for ridge, so the
  # BayesB prior acts on the same marker scale; C stays in the model as a FIXED term
  x <- if (adjusted) {
    gp_feature_scale_matrix(x[keep, , drop = FALSE], covariates = covariates)
  } else {
    gp_feature_scale_matrix(x[keep, , drop = FALSE])
  }
  y <- y[keep]
  eta <- list(markers = list(X = x, model = "BayesB"))
  if (adjusted) {
    eta$covariates <- list(
      X = covariates[, colnames(covariates) != "(Intercept)", drop = FALSE],
      model = "FIXED"
    )
  }
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
    ETA = eta,
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
                                                 n_jobs = 1L,
                                                 covariates = NULL) {
  x <- gp_feature_scale_matrix(x)
  fam <- gp_resolve_response_family(response_family, y = y)
  class_levels <- if (identical(fam, "gaussian")) NULL else {
    if (is.factor(y) || is.ordered(y)) as.character(levels(y)) else sort(unique(as.character(y)))
  }
  unranked <- NULL
  if (!is.null(covariates) && ncol(covariates) > 1L) {
    if (identical(fam, "gaussian")) {
      y <- gp_feature_residualize(as.double(y), covariates)   # spec B2: residualized trait
    } else {
      unranked <- covariates[, colnames(covariates) != "(Intercept)", drop = FALSE]
    }
  }
  params <- list(
    ntree = as.integer(ntree %||% 500L),
    mtry = mtry,
    nodesize = nodesize,
    n_jobs = as.integer(n_jobs %||% 1L),
    random_state = as.integer(seed %||% 123L),
    importance = "air",
    air_unranked = if (is.null(unranked)) 0L else as.integer(ncol(unranked))
  )
  imp <- gp_py_ml_feature_importance(
    model_type = "randomforest",
    X_train = if (is.null(unranked)) x else cbind(x, unranked),
    y_train = y,
    response_family = fam,
    model_params = params,
    class_levels = class_levels
  )
  imp <- as.numeric(imp)
  imp[!is.finite(imp)] <- 0
  active <- colSums(x^2) > 0
  imp[!active] <- 0
  attr(imp, "scoring_status") <- ifelse(active, "estimated", "no_predictor_variation")
  imp   # AIR keeps its sign (spec B3): no abs()
}

# Eligibility is separate from the signed AIR estimate. Constant predictors
# carry finite zero scores, but follow every informative predictor in a ranking.
gp_feature_inactive_rows <- function(meta) {
  if (!all(c("scoring_model", "scoring_status") %in% names(meta))) {
    return(rep(FALSE, nrow(meta)))
  }
  meta$scoring_model %in% "RandomForest" &
    meta$scoring_status %in% "no_predictor_variation"
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
                                      training_id_hash = NA_character_,
                                      scoring_engine = NA_character_,
                                      covariates_adjusted = "",
                                      n_records_dropped = 0L,
                                      ranking_source = "package",
                                      scoring_scope = "trait",
                                      extra = NULL) {
  raw_score <- as.numeric(raw_score)
  raw_score[!is.finite(raw_score)] <- 0
  # Ranks follow the raw score (AIR scores can all be negative; ranking by a normalized score
  # that is 0 for every predictor would order them by name). normalized keeps the sign.
  max_abs <- max(abs(raw_score), na.rm = TRUE)
  normalized <- if (is.finite(max_abs) && max_abs > 0) raw_score / max_abs else rep(0, length(raw_score))
  ranks <- integer(length(predictors))
  meta <- data.frame(
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
    scoring_engine = as.character(scoring_engine),
    covariates_adjusted = as.character(covariates_adjusted),
    n_records_dropped = as.integer(n_records_dropped),
    ranking_source = as.character(ranking_source),
    stringsAsFactors = FALSE
  )
  meta[["scoring_scope"]] <- scoring_scope
  if (!is.null(extra)) {
    if (!is.list(extra) || is.null(names(extra)) || anyNA(names(extra)) ||
        any(!nzchar(names(extra))) || anyDuplicated(names(extra)) ||
        any(names(extra) %in% names(meta))) {
      stop("Metadata extras must have unique new column names.", call. = FALSE)
    }
    for (nm in names(extra)) {
      if (!length(extra[[nm]]) %in% c(1L, length(predictors))) {
        stop("Metadata extra ", nm, " has invalid length.", call. = FALSE)
      }
      meta[[nm]] <- extra[[nm]]
    }
  }
  ord <- order(gp_feature_inactive_rows(meta), -raw_score, predictors, method = "radix")
  meta$rank[ord] <- seq_along(ord)
  meta[ord, , drop = FALSE]
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

gp_feature_state <- new.env(parent = emptyenv())

# Fold-internal rankings memoised within a process (spec B6 cost). Every CV task (model x k) that
# shares a fold's training records used to rescore them; now the first scores and the others reuse
# the ranking. Keyed on everything that determines the ranking (gp_feature_fold_cache_key()), so a
# hit returns exactly what scoring would. A parallel worker keeps its own.
gp_feature_fold_cache <- new.env(parent = emptyenv())

gp_feature_fold_cache_clear <- function() {
  rm(list = ls(gp_feature_fold_cache, all.names = TRUE), envir = gp_feature_fold_cache)
  invisible(NULL)
}

gp_feature_fold_cache_get <- function(key) {
  get0(key, envir = gp_feature_fold_cache, inherits = FALSE)
}

gp_feature_fold_cache_set <- function(key, value) {
  assign(key, value, envir = gp_feature_fold_cache)
  invisible(value)
}

# Everything that determines a fold's ranking: the predictor matrix (content hash, dimensions and
# names), the scored records in order (row keys, response, environment and covariate values), the
# fixed formula (as text: its environment is not part of the ranking), and every scorer setting.
gp_feature_fold_cache_key <- function(predictor_hash, predictor_data, pheno_rows, settings,
                                      fixed = NULL) {
  record_cols <- unique(c(settings$gen_name, settings$response, settings$heter_groups,
                          if (inherits(fixed, "formula")) all.vars(fixed)))
  record_cols <- intersect(record_cols, names(pheno_rows))
  gp_bridge_payload_cache_hash("feature_fold", "fold_ranking", list(
    predictors = list(hash = predictor_hash, dim = dim(predictor_data),
                      dimnames = dimnames(predictor_data)),
    records = pheno_rows[, record_cols, drop = FALSE],
    fixed = if (inherits(fixed, "formula")) paste(deparse(fixed), collapse = " ") else NULL,
    settings = settings
  ))
}

# Content hash of a predictor matrix for gp_feature_fold_view(predictor_hash =).
gp_feature_predictor_hash <- function(x) {
  gp_bridge_payload_cache_hash("feature_fold", "predictors", list(x))
}

# Ties at the k cut-off found on parallel workers warned there, not here: once the CV's results
# are back in the main process, warn once from the returned records if this run has not warned.
gp_feature_warn_returned_ties <- function(results) {
  if (isTRUE(gp_feature_state$tie_warned)) return(invisible(FALSE))
  for (res in results) {
    meta <- if (is.list(res)) res[["feature_selection_metadata"]] else NULL
    if (!is.data.frame(meta) || !nrow(meta) || !"tied_at_cutoff" %in% names(meta)) next
    first <- which(meta[["tied_at_cutoff"]] %in% TRUE)
    if (!length(first)) next
    r <- first[[1L]]
    same <- rep(TRUE, nrow(meta))
    for (col in intersect(c("trait", "replication", "fold", "selection_mode"), names(meta))) {
      same <- same & meta[[col]] %in% meta[[col]][[r]]
    }
    grp <- meta[same, , drop = FALSE]
    grp <- grp[order(grp[["rank"]], grp[["predictor"]], method = "radix"), , drop = FALSE]
    k <- suppressWarnings(as.integer(grp[["feature_k"]][[1L]]))
    s <- as.numeric(grp[["raw_score"]])
    inactive <- gp_feature_inactive_rows(grp)
    n_tied <- if (is.finite(k) && k >= 1L && k <= length(s)) {
      sum(s == s[[k]] & inactive == inactive[[k]], na.rm = TRUE)
    } else NA_integer_
    gp_feature_warn_tie(n_tied, k, as.character(grp[["trait"]][[1L]]))
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

# Top-k rows preserve recorded ranks. New rankings break ties by radix predictor name.
# Equal raw scores within the same eligibility group make the cut arbitrary:
# attr "tied_at_cutoff" is TRUE and
# the first such tie in a run warns (later ones are only recorded).
gp_feature_top_k_rows <- function(meta, k) {
  meta <- meta[order(meta[["rank"]], meta[["predictor"]], method = "radix"), , drop = FALSE]
  if (!nrow(meta)) {
    attr(meta, "tied_at_cutoff") <- FALSE
    return(meta)
  }
  k_val <- gp_feature_normalize_k(k, nrow(meta))
  top <- meta[seq_len(k_val), , drop = FALSE]
  tied <- FALSE
  if (k_val < nrow(meta) && "raw_score" %in% names(meta)) {
    s <- as.numeric(meta[["raw_score"]])
    inactive <- gp_feature_inactive_rows(meta)
    tied <- isTRUE(s[[k_val]] == s[[k_val + 1L]]) &&
      identical(inactive[[k_val]], inactive[[k_val + 1L]])
    if (tied) {
      gp_feature_warn_tie(sum(s == s[[k_val]] & inactive == inactive[[k_val]], na.rm = TRUE), k_val, as.character(meta[["trait"]][[1L]]))
    }
  }
  attr(top, "tied_at_cutoff") <- tied
  top
}

gp_feature_warn_tie <- function(n_tied, k, trait) {
  if (isTRUE(gp_feature_state$tie_warned)) return(invisible(FALSE))
  gp_feature_state$tie_warned <- TRUE
  warning(sprintf(paste(
    "Feature selection for '%s': %d predictors share the score at the top-%d cut-off, so which",
    "of them are kept follows the recorded ranking (new ties use predictor name; recorded as tied_at_cutoff). Later ties in",
    "this run are recorded without a warning."), trait, as.integer(n_tied), as.integer(k)),
    call. = FALSE)
  invisible(TRUE)
}

gp_feature_reset_tie_warning <- function() {
  gp_feature_state$tie_warned <- FALSE
  invisible(NULL)
}

# model_execute() calls this when it starts and when it ends, so a run's once-per-run warnings and
# memoised fold rankings neither carry into the next run nor stay in the session.
gp_feature_reset_run_state <- function() {
  gp_feature_reset_tie_warning()
  gp_feature_state$fixed_terms_warned <- character()
  gp_feature_fold_cache_clear()
  invisible(NULL)
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
  trait_meta[["tied_at_cutoff"]] <- isTRUE(attr(selected_train, "tied_at_cutoff"))
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
  # rank order, so the source filter below keeps the refit's columns in rank order (BGLR/RF
  # results depend on the column order) whatever order the metadata rows came in
  meta <- meta[order(meta[["rank"]], meta[["predictor"]], method = "radix"), , drop = FALSE]
  keep <- as.character(gp_feature_top_k_rows(meta, k)$predictor)
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

# Whether a fit used a selected subset (spec B5): "response" when a package-scored ranking cut the
# predictors to k < all, "supplied_ranking" when the ranking came from the user
# (ranking_source = "user": the package cannot know whether it used these phenotypes), "explicit"
# for a user-supplied feature_selected, else "none".
gp_feature_selection_kind <- function(feature_score_metadata = NULL,
                                      feature_k = NULL,
                                      feature_selected = NULL,
                                      trait = NULL) {
  if (!is.null(feature_score_metadata) && NROW(feature_score_metadata) && !is.null(feature_k) &&
      !identical(as.character(feature_k)[[1L]], "all")) {
    meta <- feature_score_metadata
    if (!is.null(trait) && "trait" %in% names(meta)) {
      meta <- meta[as.character(meta[["trait"]]) %in% as.character(trait), , drop = FALSE]
    }
    n_pred <- if (nrow(meta)) max(table(as.character(meta[["trait"]]))) else 0L
    k_val <- suppressWarnings(as.integer(feature_k[[1L]]))
    if (n_pred > 0L && is.finite(k_val) && k_val < n_pred) {
      source <- as.character(meta[["ranking_source"]] %||% character())
      return(if (length(source) && all(source %in% "user")) "supplied_ranking" else "response")
    }
  }
  if (!is.null(feature_selected)) return("explicit")
  "none"
}

gp_feature_selection_remark_text <- function(kind) {
  switch(
    kind,
    response = paste(
      "Estimated from markers selected with these phenotypes (response-based feature selection):",
      "genetic variance and heritability are overestimated, genetic correlations reflect the",
      "selected markers, and prediction-error variances, standard errors and reliabilities",
      "derived from them are optimistic. These are not population estimates."
    ),
    supplied_ranking = paste(
      "Estimated from markers selected with a supplied ranking (feature_score_metadata). If that",
      "ranking was computed using these phenotypes, genetic variance and heritability are",
      "overestimated and the prediction-error variances, standard errors and reliabilities",
      "derived from them are optimistic."
    ),
    explicit = paste(
      "Estimated from a user-supplied marker set (feature_selected). If the set was chosen using",
      "these phenotypes, genetic variance and heritability are overestimated and the",
      "prediction-error variances, standard errors and reliabilities derived from them are",
      "optimistic."
    ),
    NULL
  )
}

# Tables the remark goes on: variance components and predictions, and every covariance or
# correlation table of the result contracts (genetic, GxE, environment, residual, prediction and
# prediction-error, their SEs, ASReml's Covariance/Correlation lists): all are estimated from the
# fit on the selected markers.
gp_feature_remark_target <- function(nm) {
  nm %in% c("variance_components", "Variance_components", "predicted_values", "Predicted_value",
            "covariance_components", "correlation_components", "Covariance", "Correlation") ||
    grepl("^(Genetic|Total_genetic|Residual|Prediction|Prediction_error)_(covariance|correlation)_", nm) ||
    grepl("^(genetic|gxe|environment|residual)_(covariance|correlation)$", nm)
}

# A data frame gets a Feature_selection_remark column; a matrix (a correlation or covariance
# matrix) keeps its shape and values and carries the remark as an attribute, which print() shows
# and gp_feature_remark_export_table() turns into a column of its CSV.
gp_feature_mark_table <- function(el, text) {
  if (is.data.frame(el)) {
    if (nrow(el) && !"Feature_selection_remark" %in% names(el)) el[["Feature_selection_remark"]] <- text
  } else if (is.matrix(el) && is.null(attr(el, "Feature_selection_remark", exact = TRUE))) {
    attr(el, "Feature_selection_remark") <- text
  }
  el
}

gp_feature_remark_export_table <- function(x) {
  remark <- attr(x, "Feature_selection_remark", exact = TRUE)
  if (is.null(remark) || !is.matrix(x)) return(x)
  attr(x, "Feature_selection_remark") <- NULL
  out <- as.data.frame(x, stringsAsFactors = FALSE)
  out[["Feature_selection_remark"]] <- rep(remark, nrow(out))
  out
}

# Recursion is limited to plain lists and "predictpror_model_result": other classed objects
# (plots, fitted models) can have length()/[[ methods that disagree with names(), which made
# an unguarded walk fail with "subscript out of bounds". Unnamed lists are tolerated.
gp_feature_add_selection_remark <- function(x, kind, depth = 0L, inside_target = FALSE) {
  text <- gp_feature_selection_remark_text(kind)
  if (is.null(text) || !is.list(x) || is.data.frame(x) || depth > 6L) return(x)
  nms <- names(x)
  if (length(nms) != length(x)) nms <- rep("", length(x))
  nms[is.na(nms)] <- ""
  for (i in seq_along(x)) {
    nm <- nms[[i]]
    if (nzchar(nm) && startsWith(nm, "cv_")) next
    el <- x[[i]]
    # inside a target list (ASReml's Covariance list of matrices, a per-environment list of
    # predictions) every table is marked, whatever its own name
    target <- inside_target || gp_feature_remark_target(nm)
    if (is.data.frame(el) || is.matrix(el)) {
      if (target) x[[i]] <- gp_feature_mark_table(el, text)
    } else if (is.list(el) && (is.null(oldClass(el)) || any(oldClass(el) %in% c("list", "predictpror_model_result")))) {
      x[[i]] <- gp_feature_add_selection_remark(el, kind, depth + 1L, inside_target = target)
    }
  }
  x
}

# The selection kind of each trait's refit, from the k its best model used (best_models$feature_k,
# else the run's feature_k). One computation for the exported files (results_handling()) and the
# returned tables (gp_feature_apply_final_selection_remarks()), so both carry the same remark.
gp_feature_trait_selection_kinds <- function(traits,
                                             best_models = NULL,
                                             feature_score_metadata = NULL,
                                             feature_k = NULL,
                                             feature_selected = NULL) {
  traits <- unique(as.character(traits))
  kind_for <- function(tr) {
    k_tr <- feature_k
    if (is.data.frame(best_models) && "feature_k" %in% names(best_models) && "trait" %in% names(best_models)) {
      k_rows <- best_models[["feature_k"]][as.character(best_models[["trait"]]) %in% tr]
      k_rows <- k_rows[!is.na(k_rows)]
      if (length(k_rows)) k_tr <- k_rows[[1L]]
    }
    gp_feature_selection_kind(feature_score_metadata, k_tr, feature_selected, tr)
  }
  stats::setNames(vapply(traits, kind_for, character(1)), traits)
}

# Final stage (spec B5): the public formatters rebuild tables from column whitelists, so the
# remark is added to the finished tables, per trait, using the k each trait's best model used.
gp_feature_apply_final_selection_remarks <- function(final_results,
                                                     best_models = NULL,
                                                     feature_score_metadata = NULL,
                                                     feature_k = NULL,
                                                     feature_selected = NULL) {
  if (!is.list(final_results)) return(final_results)
  by_trait <- final_results[["model_results_by_trait"]]
  traits <- names(by_trait)
  if (!length(traits)) return(final_results)
  kinds <- gp_feature_trait_selection_kinds(traits, best_models, feature_score_metadata,
                                            feature_k, feature_selected)
  for (tr in traits) {
    final_results[["model_results_by_trait"]][[tr]] <-
      gp_feature_add_selection_remark(final_results[["model_results_by_trait"]][[tr]], kinds[[tr]])
    by_model <- final_results[["model_results_by_model"]]
    for (m in names(by_model)) {
      if (is.list(by_model[[m]]) && tr %in% names(by_model[[m]])) {
        final_results[["model_results_by_model"]][[m]][[tr]] <-
          gp_feature_add_selection_remark(final_results[["model_results_by_model"]][[m]][[tr]], kinds[[tr]])
      }
    }
  }
  # With several traits model_results duplicates one per-trait entry, which is already marked.
  if (length(traits) == 1L && !is.null(final_results[["model_results"]])) {
    final_results[["model_results"]] <-
      gp_feature_add_selection_remark(final_results[["model_results"]], kinds[[1L]])
  }
  final_results
}
