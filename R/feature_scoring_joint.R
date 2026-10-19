# Joint Gaussian feature scoring (Part A). No exported helpers.
gp_feature_joint_label <- function(response) paste(response, collapse = "+")

gp_feature_joint_inputs <- function(x, pheno_data, response, gen_name,
                                    heter_groups = NULL, complete_only = FALSE,
                                    residualize_markers = TRUE) {
  if (any(is.infinite(x))) {
    stop("Infinite predictor values are not allowed in joint scoring.", call. = FALSE)
  }
  y0 <- as.matrix(pheno_data[, response, drop = FALSE])
  storage.mode(y0) <- "double"
  if (any(is.infinite(y0))) {
    bad <- colnames(y0)[colSums(is.infinite(y0)) > 0L]
    stop("Infinite response values in trait(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  matched <- as.character(pheno_data[[gen_name]]) %in% rownames(x)
  observed <- !is.na(y0)
  for (tr in response) {
    if (!any(matched & observed[, tr])) {
      stop("Trait '", tr, "' has no visible values for joint feature scoring.",
           call. = FALSE)
    }
  }
  any_rows <- matched & rowSums(observed) > 0L
  incomplete <- any_rows & rowSums(observed) < length(response)
  n_incomplete <- if (isTRUE(complete_only)) sum(incomplete) else 0L
  rows <- which(if (isTRUE(complete_only)) {
    any_rows & !incomplete
  } else any_rows)
  if (!length(rows)) {
    stop("Joint RandomForest feature scoring has no complete records.",
         call. = FALSE)
  }
  design <- gp_feature_covariate_design(
    pheno_data[rows, , drop = FALSE], fixed = NULL,
    heter_groups = heter_groups
  )
  rows <- rows[design$keep]
  if (!length(rows)) {
    stop("Joint feature scoring has no records with complete environment values.",
         call. = FALSE)
  }
  C <- design$C
  y <- y0[rows, , drop = FALSE]
  sd_t <- stats::setNames(numeric(length(response)), response)
  for (tr in response) {
    vis <- is.finite(y[, tr])
    if (!any(vis)) {
      stop("Trait '", tr, "' has no visible values after environment adjustment.",
           call. = FALSE)
    }
    # Subsetting can make environment columns zero/aliased. lm.fit is rank-aware.
    fit <- stats::lm.fit(x = C[vis, , drop = FALSE], y = y[vis, tr])
    s <- stats::sd(fit$residuals)
    if (gp_feature_flat_after_adjustment(s, stats::sd(y[vis, tr]))) {
      stop("Trait '", tr, "' is constant or flat after environment adjustment.",
           call. = FALSE)
    }
    sd_t[[tr]] <- s
    y[, tr] <- y[, tr] / s
  }
  ids <- as.character(pheno_data[[gen_name]][rows])
  xr <- x[ids, , drop = FALSE]
  no_values <- colSums(is.finite(xr)) == 0L
  if (any(no_values)) {
    stop("Predictors have no finite values on scored records: ",
         paste(colnames(xr)[no_values], collapse = ", "), call. = FALSE)
  }
  xr <- gp_feature_scale_matrix(
    xr, covariates = if (isTRUE(residualize_markers) &&
                         !isTRUE(design$trivial)) C else NULL
  )
  pairs <- utils::combn(seq_along(response), 2L)
  pair_n <- apply(pairs, 2L, function(z) {
    sum(is.finite(y[, z[1L]]) & is.finite(y[, z[2L]]))
  })
  keys <- if (!is.null(heter_groups)) {
    paste(ids, as.character(pheno_data[[heter_groups]][rows]), sep = "\r")
  } else ids
  list(x = xr, y = y, C = C, rows = rows, trait_sd = sd_t,
       min_pair_coobserved = as.integer(min(pair_n)),
       record_keys = keys, n_train = length(rows),
       n_incomplete = as.integer(n_incomplete),
       n_unmatched = as.integer(sum(!matched)),
       n_all_missing = as.integer(sum(matched & rowSums(observed) == 0L)),
       n_dropped = as.integer(n_incomplete + design$n_dropped),
       covariates_adjusted = paste(design$terms, collapse = "+"))
}

gp_feature_joint_prior_setup <- function(x, y, C, R2 = 0.5, eigen_floor = 1e-6) {
  T <- ncol(y)
  # Prior calibration only: observed responses passed to BGLR are unchanged,
  # and C remains in the model. This is not plug-in removal of nuisance effects.
  adjusted <- y
  for (j in seq_len(T)) {
    observed <- is.finite(y[, j])
    adjusted[observed, j] <- stats::lm.fit(
      C[observed, , drop = FALSE], y[observed, j]
    )$residuals
  }
  S <- stats::cov(adjusted, use = "pairwise.complete.obs")
  S[is.na(S)] <- 0
  if (any(!is.finite(S)) || any(diag(S) <= 0)) {
    stop("Cannot calibrate joint covariance priors from visible adjusted traits.",
         call. = FALSE)
  }
  D <- diag(sqrt(diag(S)), T)
  correlation <- S / outer(diag(D), diag(D))
  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1
  minimum <- min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values)
  alpha <- if (minimum < eigen_floor) {
    min(1, (eigen_floor - minimum) / (1 - minimum))
  } else 0
  if (alpha > 0) {
    warning(sprintf("Joint scoring regularized its prior covariance scale (shrinkage %.6g).", alpha),
            call. = FALSE)
    correlation <- (1 - alpha) * correlation + alpha * diag(T)
  }
  S_pd <- D %*% correlation %*% D
  dimnames(S_pd) <- dimnames(S)
  # Match BGLR's default calibration for an intercept-only, regular panel.
  # One random marker term gets R2; FIXED terms do not divide its allocation.
  MSx <- sum(colSums(x^2)) / nrow(x) - sum(colMeans(x)^2)
  if (!is.finite(MSx) || MSx <= 0) {
    stop("Joint scoring has no informative standardized marker columns.",
         call. = FALSE)
  }
  marker_df <- T + 1L
  residual_df <- max(5L, T) # a proper residual IW shape for any supported T
  list(marker = list(type = "UN", df0 = marker_df,
                     S0 = R2 * S_pd / MSx * (marker_df + T + 1L)),
       residual = list(type = "UN", df0 = residual_df,
                       S0 = (1 - R2) * S_pd * (residual_df + T + 1L)),
       shrinkage = alpha, minimum_eigenvalue = minimum,
       covariance = S_pd)
}

gp_feature_score_joint_bglr <- function(x, y, C, model = c("BRR", "SpikeSlab"),
                                      seed = 123L, nIter = 1500L,
                                      burnIn = 500L, thin = 5L) {
  model <- match.arg(model)
  controls <- list(nIter = nIter, burnIn = burnIn, thin = thin)
  for (nm in names(controls)) {
    v <- controls[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v != floor(v) || v < (if (nm == "burnIn") 0 else 1) ||
        v > .Machine$integer.max) {
      stop(nm, " must be a finite integer in its valid range.", call. = FALSE)
    }
  }
  if (burnIn >= nIter) stop("burnIn must be smaller than nIter.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed >= .Machine$integer.max) {
    stop("seed must be a nonnegative integer below the R integer limit.", call. = FALSE)
  }
  retained <- floor(nIter / thin) - floor(burnIn / thin)
  if (retained < 2L) {
    stop("Joint scoring requires at least two retained posterior draws.",
         call. = FALSE)
  }
  # Constant columns contain no likelihood information; integrate their
  # zero-mean effects analytically instead of ranking finite-chain prior noise.
  active <- colSums(x^2) > 0
  if (!any(active)) stop("No informative predictors for joint scoring.", call. = FALSE)
  xa <- x[, active, drop = FALSE]
  # BGLR 1.1.4 sets initial R from the mean-imputed response covariance,
  # even if explicit S0 is supplied. A valid prior cannot fix singular startup.
  startup_y <- y
  for (j in seq_len(ncol(y))) {
    startup_y[!is.finite(startup_y[, j]), j] <- mean(y[, j], na.rm = TRUE)
  }
  startup_cov <- stats::cov(startup_y)
  vals <- eigen(startup_cov, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(vals)) || min(vals) <= 1e-12 * max(vals)) {
    stop("Joint BGLR cannot initialize this singular or numerically singular trait panel; ",
         "check redundant traits and visible sample size.", call. = FALSE)
  }
  prior <- gp_feature_joint_prior_setup(xa, y, C)
  eta <- list(markers = list(X = xa, model = model, Cov = prior$marker))
  if (identical(model, "SpikeSlab")) {
    eta$markers$inclusionProb <- list(
      probIn = rep(0.5, ncol(y)), counts = rep(2, ncol(y))
    )
  }
  non_intercept <- colnames(C) != "(Intercept)"
  if (any(non_intercept)) {
    eta$covariates <- list(X = C[, non_intercept, drop = FALSE], model = "FIXED")
  }
  prefix <- tempfile("feature_joint_")
  on.exit({
    # showConnections() forces GC and can emit unrelated unused-connection
    # warnings. Inspect live connections directly and normalize Windows separators.
    for (id in getAllConnections()) {
      con <- getConnection(id)
      description <- summary(con)$description
      if (startsWith(chartr("\\", "/", description), chartr("\\", "/", prefix))) {
        try(close(con), silent = TRUE)
      }
    }
    unlink(Sys.glob(paste0(prefix, "*")))
  }, add = TRUE)
  invisible(utils::capture.output(
    suppressMessages(fit <- gp_with_pinned_seed(seed, BGLR::Multitrait(
      y = y, ETA = eta, resCov = prior$residual, nIter = as.integer(nIter),
      burnIn = as.integer(burnIn), thin = as.integer(thin),
      verbose = FALSE, saveAt = prefix
    )))
  ))
  beta <- as.matrix(fit$ETA$markers$beta)
  sd_beta <- as.matrix(fit$ETA$markers$SD.beta)
  target_dim <- c(ncol(xa), ncol(y))
  if (!identical(dim(beta), target_dim) ||
      !identical(dim(sd_beta), target_dim) ||
      any(!is.finite(beta)) || any(!is.finite(sd_beta))) {
    stop("BGLR joint feature scoring returned invalid effect matrices.", call. = FALSE)
  }
  score <- rep(0, ncol(x))
  effect_sd <- rep(NA_real_, ncol(x))
  score[active] <- sqrt(rowSums(beta^2))
  effect_sd[active] <- sqrt(rowSums(sd_beta^2))
  list(score = score, effect_sd = effect_sd,
       scoring_status = ifelse(active, "estimated", "no_predictor_variation"),
       prior_covariance_shrinkage = prior$shrinkage,
       prior_minimum_eigenvalue = prior$minimum_eigenvalue,
       prior_residual_df0 = prior$residual$df0,
       retained_draws = as.integer(retained))
}

# Internal Gaussian multi-response scoring; all outcome-dependent calculations
# use the phenotype cells supplied by the caller (masked in fold-internal CV).
gp_feature_score_joint <- function(x, pheno_data, response, gen_name,
                                   scoring_model, seed, source_block,
                                   bayes_nIter, bayes_burnIn, bayes_thin,
                                   ntree, mtry, nodesize, rf_n_jobs,
                                   response_family, selection_mode,
                                   replication, fold, fixed, heter_groups) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed >= .Machine$integer.max) {
    stop("seed must be a nonnegative integer below the R integer limit.", call. = FALSE)
  }
  if (length(response) < 2L || anyDuplicated(response)) {
    stop("Joint feature scoring requires at least two distinct responses.", call. = FALSE)
  }
  if (any(grepl("+", response, fixed = TRUE))) {
    stop("Joint feature labels reserve '+'; rename response columns containing '+'.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x)) || anyNA(rownames(x)) ||
      any(!nzchar(rownames(x))) || anyDuplicated(colnames(x)) ||
      anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("Joint predictors require unique, nonempty individual and predictor names.", call. = FALSE)
  }
  for (tr in response) {
    fam <- gp_resolve_response_family(response_family, y = pheno_data[[tr]])
    if (!identical(fam, "gaussian") || !is.numeric(pheno_data[[tr]])) {
      stop("Joint feature scoring requires numeric Gaussian traits; check '", tr, "'.",
           call. = FALSE)
    }
  }
  if (!is.null(fixed) &&
      (!inherits(fixed, "formula") ||
       length(attr(stats::terms(fixed), "term.labels")) ||
       attr(stats::terms(fixed), "intercept") != 1L)) {
    stop("Joint feature scoring sets its own intercepts and environment effects; ",
         "`fixed` must be NULL or an intercept-only formula.", call. = FALSE)
  }
  if (!is.null(heter_groups) &&
      (length(heter_groups) != 1L || is.na(heter_groups) ||
       !nzchar(heter_groups) || !heter_groups %in% names(pheno_data))) {
    stop("Joint scoring requires an existing, single environment column.", call. = FALSE)
  }
  is_rf <- identical(scoring_model, "RandomForest")
  input <- gp_feature_joint_inputs(
    x, pheno_data, response, gen_name, heter_groups,
    complete_only = is_rf, residualize_markers = !is_rf
  )
  scored <- if (is_rf) {
    if (input$n_incomplete > 0L) {
      warning(sprintf("Joint RandomForest scoring dropped %d incomplete-response records.",
                      input$n_incomplete), call. = FALSE)
    }
    gp_feature_score_joint_randomforest(input$x, input$y, input$C,
                                       seed, ntree, mtry, nodesize, rf_n_jobs)
  } else {
    gp_feature_score_joint_bglr(
      input$x, input$y, input$C,
      model = if (identical(scoring_model, "BayesB")) "SpikeSlab" else "BRR",
      seed = seed, nIter = bayes_nIter, burnIn = bayes_burnIn, thin = bayes_thin
    )
  }
  engine <- switch(scoring_model,
    Ridge_Regression = "BGLR Multitrait BRR",
    BayesB = "BGLR Multitrait SpikeSlab",
    RandomForest = "sklearn RandomForest AIR (multi-output)"
  )
  extra <- list(
    min_pair_coobserved = input$min_pair_coobserved,
    effect_sd = scored$effect_sd,
    chain_nIter = if (is_rf) NA_integer_ else as.integer(bayes_nIter),
    chain_burnIn = if (is_rf) NA_integer_ else as.integer(bayes_burnIn),
    chain_thin = if (is_rf) NA_integer_ else as.integer(bayes_thin),
    scoring_status = scored$scoring_status %||% rep("estimated", ncol(x)),
    prior_covariance_shrinkage = scored$prior_covariance_shrinkage %||% NA_real_,
    prior_minimum_eigenvalue = scored$prior_minimum_eigenvalue %||% NA_real_,
    prior_residual_df0 = scored$prior_residual_df0 %||% NA_integer_,
    retained_draws = scored$retained_draws %||% NA_integer_,
    n_incomplete_records = input$n_incomplete,
    n_unmatched_records = input$n_unmatched,
    n_all_missing_records = input$n_all_missing
  )
  meta <- gp_feature_metadata_table(
    trait = gp_feature_joint_label(response), predictors = colnames(x),
    source_block = source_block, raw_score = scored$score,
    scoring_model = scoring_model, seed = seed, n_train = input$n_train,
    response_family = "gaussian", selection_mode = selection_mode,
    replication = replication, fold = fold,
    training_id_hash = gp_feature_training_id_hash(input$record_keys),
    scoring_engine = engine, covariates_adjusted = input$covariates_adjusted,
    n_records_dropped = input$n_dropped, scoring_scope = "joint", extra = extra
  )
  rownames(meta) <- NULL
  class(meta) <- c("predictpror_feature_scores", class(meta))
  meta
}


gp_feature_score_joint_randomforest <- function(x, y, C, seed = 123L,
                                               ntree = 500L, mtry = NULL,
                                               nodesize = NULL, n_jobs = 1L) {
  active <- colSums(x^2) > 0
  if (!any(active)) stop("No informative predictors for joint scoring.", call. = FALSE)
  y_res <- vapply(seq_len(ncol(y)), function(j) {
    stats::lm.fit(C, y[, j])$residuals
  }, numeric(nrow(y)))
  colnames(y_res) <- colnames(y)
  imp <- gp_py_ml_feature_importance(
    model_type = "randomforest", X_train = x, y_train = y_res,
    response_family = "gaussian",
    model_params = list(
      ntree = as.integer(ntree), mtry = mtry, nodesize = nodesize,
      n_jobs = as.integer(n_jobs), random_state = as.integer(seed),
      importance = "air", air_unranked = 0L
    ),
    class_levels = NULL
  )
  imp <- as.numeric(imp)
  if (length(imp) != ncol(x) || any(!is.finite(imp))) {
    stop("Joint RandomForest returned invalid AIR importances.", call. = FALSE)
  }
  imp[!active] <- 0
  list(score = imp, effect_sd = rep(NA_real_, ncol(x)),
       scoring_status = ifelse(active, "estimated", "no_predictor_variation"))
}


gp_feature_joint_rows <- function(meta, response) {
  advice <- "Pass a joint ranking from feature_score_predictors(joint = TRUE)."
  required <- c("trait", "scoring_scope", "predictor", "rank", "raw_score")
  if (!is.data.frame(meta) || !nrow(meta) || !all(required %in% names(meta))) {
    stop("Joint feature selection requires a complete joint score table. ", advice, call. = FALSE)
  }
  if (anyNA(meta$scoring_scope) || any(as.character(meta$scoring_scope) != "joint")) {
    stop("Per-trait rankings cannot be used by a joint model. ", advice, call. = FALSE)
  }
  labels <- unique(as.character(meta$trait))
  response <- as.character(response)
  label_traits <- if (length(labels) == 1L && !is.na(labels)) {
    strsplit(labels, "+", fixed = TRUE)[[1L]]
  } else character()
  if (length(labels) != 1L || is.na(labels) || length(response) < 2L ||
      anyDuplicated(response) || anyDuplicated(label_traits) ||
      length(label_traits) != length(response) || !setequal(label_traits, response)) {
    stop("The joint ranking's traits must match the run's responses.", call. = FALSE)
  }
  predictors <- as.character(meta$predictor)
  if (anyNA(predictors) || any(!nzchar(predictors)) || anyDuplicated(predictors)) {
    stop("Joint ranking predictor names must be unique and nonempty.", call. = FALSE)
  }
  if (!is.numeric(meta$rank) || any(!is.finite(meta$rank)) ||
      any(meta$rank < 1 | meta$rank != floor(meta$rank)) || anyDuplicated(meta$rank)) {
    stop("Joint ranking ranks must be unique positive finite integers.", call. = FALSE)
  }
  if (!is.numeric(meta$raw_score) || any(!is.finite(meta$raw_score))) {
    stop("Joint ranking raw scores must be finite numeric values.", call. = FALSE)
  }
  ranked <- order(meta$rank, method = "radix")
  inactive <- gp_feature_inactive_rows(meta)
  # Supplied ranks resolve ties already: accept historical/cross-locale tie orders.
  # RF constants follow active columns even when the latter have negative AIR.
  if (any(diff(as.integer(inactive[ranked])) < 0) ||
      any(meta$raw_score[inactive] != 0) ||
      any(diff(meta$raw_score[ranked[!inactive[ranked]]]) > 0)) {
    stop("Joint ranking ranks must follow descending raw scores among informative predictors, with zero-variation RandomForest predictors last.", call. = FALSE)
  }
  meta
}

gp_feature_check_joint_ranking <- function(meta, response) {
  gp_feature_joint_rows(meta, response)
  invisible(meta)
}

gp_feature_joint_view <- function(predictor_data, pheno_visible, response,
                                  gen_name, k, selection_mode = "fixed",
                                  feature_score_metadata = NULL,
                                  scoring_model = "Ridge_Regression",
                                  seed = 123L, source_block = NULL,
                                  response_family = "gaussian",
                                  replication = NA_integer_, fold = NA_integer_,
                                  ridge_lambda = 1, bayes_nIter = 1500L,
                                  bayes_burnIn = 500L, bayes_thin = 5L,
                                  ntree = 500L, mtry = NULL, nodesize = NULL,
                                  rf_n_jobs = 1L, heter_groups = NULL,
                                  held_out_trait = NA_character_) {
  if (!identical(k, "all") &&
      (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
       k < 1 || k != floor(k))) {
    stop("Joint feature k must be a positive integer or all.", call. = FALSE)
  }
  selection_mode <- gp_feature_normalize_cv_policy(selection_mode,
                                                    allow_both = FALSE)
  if (!selection_mode %in% c("fixed", "fold_internal")) {
    stop("Joint feature selection requires fixed or fold_internal mode.",
         call. = FALSE)
  }
  meta <- if (identical(selection_mode, "fixed")) {
    gp_feature_joint_rows(feature_score_metadata, response)
  } else {
    feature_score_predictors(
      predictor_data, pheno_visible, response, gen_name,
      scoring_model = scoring_model, seed = seed, source_block = source_block,
      ridge_lambda = ridge_lambda, bayes_nIter = bayes_nIter,
      bayes_burnIn = bayes_burnIn, bayes_thin = bayes_thin,
      ntree = ntree, mtry = mtry, nodesize = nodesize, rf_n_jobs = rf_n_jobs,
      response_family = response_family, selection_mode = "fold_internal",
      replication = replication, fold = fold, fixed = NULL,
      heter_groups = heter_groups, joint = TRUE
    )
  }
  missing <- setdiff(as.character(meta$predictor), colnames(predictor_data))
  unranked <- setdiff(colnames(predictor_data), as.character(meta$predictor))
  if (length(missing) || length(unranked)) {
    example_names <- head(c(missing, unranked), 10L)
    stop("Joint ranking must contain every predictor in the model matrix: ",
         length(missing), " predictors in the ranking are absent from the model matrix; ",
         length(unranked), " predictors in the model matrix are absent from the ranking. ",
         "Examples (at most 10): ", paste(example_names, collapse = ", "), ". ",
         "Score the post-QC matrix, or use qc_filtering = FALSE when intentionally using the raw matrix.",
         call. = FALSE)
  }
  selected <- gp_feature_top_k_rows(meta, k)
  keep <- as.character(selected$predictor)
  meta$selected <- meta$predictor %in% keep
  meta$feature_k <- as.integer(length(keep))
  meta$tied_at_cutoff <- isTRUE(attr(selected, "tied_at_cutoff"))
  meta$replication <- as.integer(replication)
  meta$fold <- as.integer(fold)
  meta$selection_mode <- selection_mode
  meta$held_out_trait <- as.character(held_out_trait)
  list(predictor_data = predictor_data[, keep, drop = FALSE],
       metadata = meta, selected_predictors = keep,
       feature_k = as.integer(length(keep)))
}

gp_feature_joint_leaf_signature <- function(value) {
  if (is.factor(value)) value <- as.character(value)
  if (is.logical(value)) {
    if (anyNA(value)) stop("Missing logical values cannot define a shared predictor set.", call. = FALSE)
    return(list(type = "mask", size = length(value), value = which(value)))
  }
  if (is.numeric(value)) {
    if (any(!is.finite(value)) || any(value < 1 | value != floor(value))) {
      stop("Shared predictor indices must be positive finite integers.", call. = FALSE)
    }
    return(list(type = "indices", value = sort(unique(as.double(value)))))
  }
  value <- as.character(value %||% character())
  if (anyNA(value) || any(!nzchar(value))) {
    stop("Shared predictor names must be nonempty and nonmissing.", call. = FALSE)
  }
  list(type = "names", value = sort(unique(value)))
}

gp_feature_joint_set_signature <- function(value) {
  if (gp_feature_has_source_scope(value)) {
    sources <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")
    allowed <- tolower(gp_feature_source_aliases())
    if (is.data.frame(value)) {
      source_column <- intersect(c("source", "source_block", "block", "input",
                                    "data_source", "feature_source"), names(value))[[1L]]
      flags <- intersect(c("selected", "Selected"), names(value))
      keep <- rep(TRUE, nrow(value))
      if (length(flags) && is.logical(value[[flags[[1L]]]])) {
        keep <- !is.na(value[[flags[[1L]]]]) & value[[flags[[1L]]]]
      }
      declared <- as.character(value[[source_column]][keep])
    } else declared <- names(value)
    if (anyNA(declared) || any(!tolower(declared) %in% allowed)) {
      stop("Unrecognized source key in a per-trait predictor set.", call. = FALSE)
    }
    if (is.list(value) && !is.data.frame(value)) {
      nms <- names(value)
      canonical <- vapply(nms, function(nm) {
        hit <- sources[vapply(sources, function(s) {
          tolower(nm) %in% tolower(gp_feature_source_aliases(s))
        }, logical(1L))]
        if (length(hit)) hit[[1L]] else tolower(nm)
      }, character(1L))
      if (anyDuplicated(canonical)) {
        stop("Repeated source aliases make a per-trait predictor set ambiguous.", call. = FALSE)
      }
    }
    entries <- lapply(sources, function(s) {
      gp_feature_joint_leaf_signature(gp_feature_selected_vector(value, source = s))
    })
    names(entries) <- sources
    entries <- entries[vapply(entries, function(x) length(x$value) > 0L, logical(1L))]
    return(list(scope = "source", value = entries))
  }
  gp_feature_joint_leaf_signature(gp_feature_selected_vector(value))
}

gp_feature_joint_selected <- function(feature_selected, response) {
  if (is.null(feature_selected)) return(NULL)
  response <- as.character(response)
  message <- paste(
    "A joint model uses one predictor set for all traits;",
    "supply one shared set or identical sets for every response."
  )
  pieces <- NULL
  result <- feature_selected
  if (is.data.frame(feature_selected) && "trait" %in% names(feature_selected)) {
    tr <- as.character(feature_selected$trait)
    labels <- unique(tr)
    joint_traits <- if (length(labels) == 1L && !is.na(labels)) {
      strsplit(labels, "+", fixed = TRUE)[[1L]]
    } else character()
    if (length(response) > 1L && !anyDuplicated(joint_traits) &&
        length(joint_traits) == length(response) && setequal(joint_traits, response)) {
      feature_selected$trait <- NULL
      return(feature_selected)
    }
    missing <- setdiff(response, tr)
    if (length(missing)) {
      stop(message, " Missing trait entry: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    if (anyNA(tr) || any(!tr %in% response)) {
      stop(message, " Trait entries must match the requested responses.", call. = FALSE)
    }
    pieces <- lapply(response, function(r) {
      x <- feature_selected[tr == r, , drop = FALSE]
      x$trait <- NULL
      x
    })
    result <- pieces[[1L]]
  } else if (is.list(feature_selected) && !is.data.frame(feature_selected)) {
    nms <- names(feature_selected)
    if (is.null(nms) && length(feature_selected) == 1L) {
      return(gp_feature_joint_selected(feature_selected[[1L]], response))
    }
    if (any(nms %in% response)) {
      if (any(tolower(nms[nms %in% response]) %in% tolower(gp_feature_source_aliases()))) {
        stop("Ambiguous trait/source keys in feature_selected; supply one shared vector.", call. = FALSE)
      }
      missing <- setdiff(response, nms)
      if (length(missing)) {
        stop(message, " Missing trait entry: ", paste(missing, collapse = ", "), call. = FALSE)
      }
      if (anyDuplicated(nms) || any(!nms %in% response)) {
        stop(message, " Trait entries must match the requested responses.", call. = FALSE)
      }
      pieces <- feature_selected[response]
      result <- pieces[[1L]]
    } else {
      common <- c(gp_feature_source_aliases(), "merged", "predictor", "predictors",
                  "feature", "features", "selected")
      if (is.null(nms) || anyNA(nms) || any(!tolower(nms) %in% tolower(common))) {
        stop(message, " Use response-keyed entries or recognized source keys.", call. = FALSE)
      }
    }
  }
  if (!is.null(pieces)) {
    signatures <- lapply(pieces, gp_feature_joint_set_signature)
    if (any(vapply(signatures, function(x) !length(x$value), logical(1L)))) {
      stop(message, " Per-trait sets cannot be empty.", call. = FALSE)
    }
    if (!all(vapply(signatures, function(x) identical(x, signatures[[1L]]), logical(1L)))) {
      stop(message, call. = FALSE)
    }
  }
  result
}

# Called only by the joint true-prediction branches.
gp_feature_attach_joint_output <- function(result, response, selection = NULL,
                                           ml_dat_res = NULL, feature_selected = NULL,
                                           source_map = NULL) {
  view <- if (is.null(selection)) NULL else selection$view
  meta <- if (is.null(view)) ml_dat_res[["feature_selection_metadata"]] else view$metadata
  keep <- if (is.null(view)) ml_dat_res[["feature_selected"]] else view$selected_predictors
  if (is.null(keep) && !is.null(feature_selected)) {
    keep <- unique(as.character(gp_feature_selected_vector(feature_selected)))
  }
  if (is.null(keep) || !length(keep)) return(result)
  keep <- as.character(keep)
  if (is.null(meta)) {
    sources <- gp_feature_source_block(keep, source_map)
    meta <- data.frame(
      trait = gp_feature_joint_label(response), predictor = keep,
      source_block = unname(sources[keep]), rank = NA_integer_,
      raw_score = NA_real_, scoring_model = "user_supplied",
      scoring_scope = "joint", ranking_source = "user",
      selection_mode = "explicit", selected = TRUE,
      feature_k = as.integer(length(keep)), tied_at_cutoff = FALSE,
      stringsAsFactors = FALSE
    )
  }
  result[["feature_selection_metadata"]] <- meta
  result[["feature_selected"]] <- keep
  result[["feature_k"]] <- as.integer(length(keep))
  result
}
