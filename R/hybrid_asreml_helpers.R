gp_hybrid_asreml_supported_models <- function() c("GBLUP")

gp_hybrid_cv_supported_methods <- function() {
  c(
    "Hybrid_Known_Parents",
    "Hybrid_One_New_Parent",
    "Hybrid_Both_New_Parents"
  )
}

gp_hybrid_asreml_summary_statistics <- function(pred_df,
                                                response,
                                                gen_name,
                                                female_parent,
                                                male_parent,
                                                model_type) {
  data.frame(
    stat = c(
      "mode",
      "response_family",
      "model_type",
      "scope",
      "response",
      "observed_hybrids",
      "predicted_hybrids",
      "n_female_parents",
      "n_male_parents",
      "female_parent_column",
      "male_parent_column",
      "hybrid_components"
    ),
    summary = c(
      "hybrid_asreml",
      "gaussian",
      model_type,
      "hybrid Gaussian true prediction via GCA_f + GCA_m + SCA + ASReml-R",
      response,
      sum(pred_df$Train_Test_Label == "Train", na.rm = TRUE),
      sum(pred_df$Train_Test_Label == "Test", na.rm = TRUE),
      length(unique(as.character(pred_df[[female_parent]]))),
      length(unique(as.character(pred_df[[male_parent]]))),
      female_parent,
      male_parent,
      "female_gca;male_gca;sca"
    ),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_asreml_diagnostic_plot <- function(pred_df) {
  obs <- pred_df[!is.na(pred_df$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Hybrid Gaussian ASReml-R observed vs predicted",
      subtitle = "Hybrid predictions include female GCA, male GCA, and SCA components",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_hybrid_asreml_cv_plot <- function(pred_df) {
  obs <- pred_df[!is.na(pred_df$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = cv_scenario)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~ cv_scenario, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Hybrid ASReml-R cross-validation observed vs predicted",
      subtitle = "Hybrid scenarios: known parents, one new parent, both parents new",
      x = "Observed value",
      y = "Predicted value",
      color = "Scenario"
    )
}

gp_hybrid_cv_normalize_method <- function(x) {
  if (is.null(x) || !length(x)) {
    return("hybrid_known_parents")
  }
  x <- tolower(as.character(x)[1])
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

gp_hybrid_observed_rows <- function(pheno_object, response) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  keep <- !is.na(ph[[response]])
  ph[keep, , drop = FALSE]
}

gp_hybrid_row_key <- function(dat, gen_name, heter_groups = NULL) {
  if (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(dat)) {
    return(as.character(dat[[gen_name]]))
  }
  paste(as.character(dat[[gen_name]]), as.character(dat[[heter_groups]]), sep = "::__env__::")
}

# Phenotypes for one CV fold refit: every row outside the fold's training set
# is masked, not only the held-out rows. For Hybrid_Both_New_Parents the
# training set excludes all hybrids sharing either parent with the held-out
# cross; masking only the held-out cross let kernel models train on those
# parents, so they were not new. Other scenarios train on every non-test row,
# so for them this equals masking the test rows.
gp_hybrid_cv_mask_responses <- function(ph, responses, train_rows) {
  not_train <- setdiff(seq_len(nrow(ph)), train_rows)
  for (resp in responses) {
    ph[[resp]][not_train] <- NA_real_
  }
  ph
}

# Predictions for one CV fold: only the fold's held-out rows. A fold refit
# labels every NA-response row "Test", which also includes hybrids that were
# already NA in the input (true-prediction targets); those have no observed
# value and must not be counted or scored as held-out predictions.
gp_hybrid_fold_test_predictions <- function(pred_all,
                                            ph,
                                            test_rows,
                                            gen_name,
                                            heter_groups = NULL) {
  if (is.null(pred_all) || !nrow(pred_all)) {
    return(pred_all)
  }
  test_keys <- gp_hybrid_row_key(ph[test_rows, , drop = FALSE],
                                 gen_name = gen_name, heter_groups = heter_groups)
  pred_keys <- gp_hybrid_row_key(pred_all, gen_name = gen_name, heter_groups = heter_groups)
  pred_all[pred_keys %in% test_keys, , drop = FALSE]
}

gp_hybrid_match_observed_values <- function(pred_df,
                                            ph,
                                            response,
                                            gen_name,
                                            heter_groups = NULL) {
  if (is.null(pred_df) || !nrow(pred_df)) {
    return(pred_df)
  }
  pred_key <- gp_hybrid_row_key(pred_df, gen_name = gen_name, heter_groups = heter_groups)
  ph_key <- gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups)
  pred_df$Observed_value <- ph[[response]][match(pred_key, ph_key)]
  pred_df
}

gp_hybrid_cv_folds_known_parents <- function(ph,
                                             gen_name,
                                             female_parent,
                                             male_parent,
                                             nfolds = 5L,
                                             random_state = 123L,
                                             replication = 1L) {
  n <- nrow(ph)
  if (n < 2L) {
    return(list())
  }
  gp_set_seed(as.integer(random_state) + as.integer(replication) * 1000L)
  fold_ids <- sample(rep(seq_len(min(nfolds, n)), length.out = n))
  out <- list()
  for (fold in sort(unique(fold_ids))) {
    test_idx <- which(fold_ids == fold)
    train_idx <- setdiff(seq_len(n), test_idx)
    train_f <- unique(as.character(ph[[female_parent]][train_idx]))
    train_m <- unique(as.character(ph[[male_parent]][train_idx]))
    keep_test <- test_idx[
      as.character(ph[[female_parent]][test_idx]) %in% train_f &
        as.character(ph[[male_parent]][test_idx]) %in% train_m
    ]
    if (!length(keep_test)) {
      next
    }
    out[[length(out) + 1L]] <- list(
      scenario = "Known_Parents_New_Hybrid",
      fold = paste0("known_parents_fold_", fold),
      test_ids = as.character(ph[[gen_name]][keep_test])
    )
  }
  out
}

gp_hybrid_cv_folds_one_new_parent <- function(ph,
                                              gen_name,
                                              female_parent,
                                              male_parent) {
  out <- list()

  female_levels <- unique(as.character(ph[[female_parent]]))
  for (f in female_levels) {
    test_idx <- which(as.character(ph[[female_parent]]) == f)
    train_idx <- setdiff(seq_len(nrow(ph)), test_idx)
    if (!length(test_idx) || !length(train_idx)) next
    train_f <- unique(as.character(ph[[female_parent]][train_idx]))
    train_m <- unique(as.character(ph[[male_parent]][train_idx]))
    male_test <- unique(as.character(ph[[male_parent]][test_idx]))
    if (f %in% train_f) next
    if (!all(male_test %in% train_m)) next
    out[[length(out) + 1L]] <- list(
      scenario = "One_New_Parent",
      fold = paste0("new_female_", f),
      heldout_parent = f,
      heldout_side = "female",
      test_ids = as.character(ph[[gen_name]][test_idx])
    )
  }

  male_levels <- unique(as.character(ph[[male_parent]]))
  for (m in male_levels) {
    test_idx <- which(as.character(ph[[male_parent]]) == m)
    train_idx <- setdiff(seq_len(nrow(ph)), test_idx)
    if (!length(test_idx) || !length(train_idx)) next
    train_f <- unique(as.character(ph[[female_parent]][train_idx]))
    train_m <- unique(as.character(ph[[male_parent]][train_idx]))
    female_test <- unique(as.character(ph[[female_parent]][test_idx]))
    if (m %in% train_m) next
    if (!all(female_test %in% train_f)) next
    out[[length(out) + 1L]] <- list(
      scenario = "One_New_Parent",
      fold = paste0("new_male_", m),
      heldout_parent = m,
      heldout_side = "male",
      test_ids = as.character(ph[[gen_name]][test_idx])
    )
  }

  out
}

gp_hybrid_cv_folds_both_new_parents <- function(ph,
                                                gen_name,
                                                female_parent,
                                                male_parent) {
  out <- list()
  pairs <- unique(ph[, c(female_parent, male_parent), drop = FALSE])
  for (i in seq_len(nrow(pairs))) {
    f <- as.character(pairs[[female_parent]][i])
    m <- as.character(pairs[[male_parent]][i])
    test_idx <- which(as.character(ph[[female_parent]]) == f & as.character(ph[[male_parent]]) == m)
    train_idx <- which(as.character(ph[[female_parent]]) != f & as.character(ph[[male_parent]]) != m)
    if (!length(test_idx) || !length(train_idx)) next
    train_f <- unique(as.character(ph[[female_parent]][train_idx]))
    train_m <- unique(as.character(ph[[male_parent]][train_idx]))
    if (f %in% train_f || m %in% train_m) next
    out[[length(out) + 1L]] <- list(
      scenario = "Both_New_Parents",
      fold = paste0("new_pair_", f, "_", m),
      heldout_female = f,
      heldout_male = m,
      test_ids = as.character(ph[[gen_name]][test_idx]),
      train_ids = as.character(ph[[gen_name]][train_idx])
    )
  }
  out
}

gp_hybrid_build_cv_scenarios <- function(pheno_object,
                                         response,
                                         gen_name,
                                         female_parent,
                                         male_parent,
                                         heter_groups = NULL,
                                         cross_validation_meth,
                                         nfolds = 5L,
                                         random_state = 123L,
                                         replication = 1L) {
  ph_all <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  observed_rows <- which(!is.na(ph_all[[response]]))
  ph <- ph_all[observed_rows, , drop = FALSE]
  method <- gp_hybrid_cv_normalize_method(cross_validation_meth)
  if (!nrow(ph)) {
    return(list())
  }
  met_methods <- c(
    "cv0", "cv1", "cv2",
    "repeated_cv0", "repeated_cv1", "repeated_cv2"
  )
  is_met <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    !is.na(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph) &&
    length(unique(stats::na.omit(as.character(ph[[heter_groups]])))) > 1L
  if (isTRUE(is_met)) {
    if (!method %in% met_methods) {
      stop(
        paste(
          "Hybrid multi-environment cross-validation requires CV0, CV1, CV2,",
          "Repeated_CV0, Repeated_CV1, or Repeated_CV2. Parent-novelty",
          "scenarios remain available for single-environment hybrid analysis."
        ),
        call. = FALSE
      )
    }
    assignments <- gp_met_cv_assignments(
      pheno_data = ph,
      gen_name = gen_name,
      heter_groups = heter_groups,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      random_state = as.integer(random_state %||% 123L) + as.integer(replication - 1L),
      replication = 1L,
      message = FALSE
    )[[1L]]
    fold_values <- sort(unique(assignments[is.finite(assignments) & assignments > 0L]))
    scenario_name <- toupper(sub("^repeated_", "", method))
    return(lapply(fold_values, function(fold_value) {
      test_local <- which(assignments == fold_value)
      train_local <- which(!is.na(assignments) & assignments != fold_value)
      test_rows <- observed_rows[test_local]
      train_rows <- observed_rows[train_local]
      list(
        scenario = scenario_name,
        fold = paste0(tolower(scenario_name), "_fold_", fold_value),
        test_rows = test_rows,
        train_rows = train_rows,
        test_ids = unique(as.character(ph[[gen_name]][test_local])),
        train_ids = unique(as.character(ph[[gen_name]][train_local])),
        test_envs = unique(as.character(ph[[heter_groups]][test_local]))
      )
    }))
  }

  ph <- unique(ph[, c(gen_name, female_parent, male_parent), drop = FALSE])
  if (identical(method, "hybrid_known_parents")) {
    return(gp_hybrid_cv_folds_known_parents(
      ph = ph,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      nfolds = nfolds,
      random_state = random_state,
      replication = replication
    ))
  }
  if (identical(method, "hybrid_one_new_parent")) {
    return(gp_hybrid_cv_folds_one_new_parent(
      ph = ph,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent
    ))
  }
  if (identical(method, "hybrid_both_new_parents")) {
    return(gp_hybrid_cv_folds_both_new_parents(
      ph = ph,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent
    ))
  }
  stop(
    paste(
      "Unsupported hybrid cross-validation method.",
      "Choose from:",
      paste(gp_hybrid_cv_supported_methods(), collapse = ", ")
    ),
    call. = FALSE
  )
}

gp_hybrid_cv_scenario_rows <- function(sc, ph, gen_name) {
  n <- nrow(ph)
  test_rows <- suppressWarnings(as.integer(sc$test_rows %||% integer()))
  test_rows <- unique(test_rows[is.finite(test_rows) & test_rows >= 1L & test_rows <= n])
  if (!length(test_rows)) {
    test_ids <- unique(as.character(sc$test_ids %||% character()))
    test_rows <- which(as.character(ph[[gen_name]]) %in% test_ids)
  }
  train_rows <- suppressWarnings(as.integer(sc$train_rows %||% integer()))
  train_rows <- unique(train_rows[is.finite(train_rows) & train_rows >= 1L & train_rows <= n])
  if (!length(train_rows)) {
    if (!is.null(sc$train_ids)) {
      train_ids <- unique(as.character(sc$train_ids))
      train_rows <- which(as.character(ph[[gen_name]]) %in% train_ids)
    } else {
      train_rows <- setdiff(seq_len(n), test_rows)
    }
  }
  train_rows <- setdiff(train_rows, test_rows)
  list(train = train_rows, test = test_rows)
}

gp_hybrid_asreml_cv_metric_table <- function(pred_df,
                                             eval_metrics,
                                             response_family = "gaussian") {
  if (!nrow(pred_df)) {
    return(data.frame())
  }
  split_keys <- unique(pred_df[, c("cv_scenario", "fold", "rep"), drop = FALSE])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    idx <- pred_df$cv_scenario == key$cv_scenario &
      pred_df$fold == key$fold &
      pred_df$rep == key$rep
    sub <- pred_df[idx, , drop = FALSE]
    row <- data.frame(
      cv_scenario = key$cv_scenario,
      fold = key$fold,
      rep = key$rep,
      stringsAsFactors = FALSE
    )
    for (m in eval_metrics) {
      row[[m]] <- safe_metric_value(
        y_true = sub$Observed_value,
        y_pred = sub$Predicted_value,
        metric = m,
        response_family = response_family
      )
    }
    row
  })
  do.call(rbind, out)
}

gp_hybrid_asreml_cv_process <- function(pred_df,
                                        eval_df,
                                        response,
                                        model_type = "GBLUP") {
  agg <- if (is.null(eval_df) || !nrow(eval_df)) {
    data.frame()
  } else {
    numeric_cols <- names(eval_df)[vapply(eval_df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, "rep")
    stats::aggregate(
      eval_df[, numeric_cols, drop = FALSE],
      by = list(cv_scenario = eval_df$cv_scenario),
      FUN = mean,
      na.rm = TRUE
    )
  }
  if (nrow(agg)) {
    agg$model <- as.character(model_type)[1L]
    agg <- agg[, c("cv_scenario", "model", setdiff(names(agg), c("cv_scenario", "model"))), drop = FALSE]
  }

  counts <- if (!nrow(pred_df)) {
    data.frame()
  } else {
    stats::aggregate(
      pred_df$HybridID,
      by = list(cv_scenario = pred_df$cv_scenario),
      FUN = length
    )
  }
  if (nrow(counts)) {
    names(counts)[names(counts) == "x"] <- "n_test_hybrids"
  }

  list(
    hybrid_metric_summary = agg,
    hybrid_prediction_counts = counts,
    hybrid_cv_predictions = pred_df,
    hybrid_cv_plot = gp_hybrid_asreml_cv_plot(pred_df)
  )
}

gp_hybrid_prepare_parent_gmatrix <- function(parent_ids,
                                             gmatrix = NULL,
                                             geno_data = NULL,
                                             gmatrix_method = NULL,
                                             ploidy = "auto",
                                             label = "parent") {
  parent_ids <- unique(as.character(parent_ids))
  if (!length(parent_ids)) {
    stop(paste("No", label, "parent IDs were supplied."), call. = FALSE)
  }

  if (!is.null(gmatrix)) {
    g_use <- as.matrix(gmatrix)
    if (is.null(rownames(g_use)) || is.null(colnames(g_use))) {
      stop(paste("The", label, "gmatrix must have row and column names."), call. = FALSE)
    }
    missing_ids <- setdiff(parent_ids, rownames(g_use))
    if (length(missing_ids)) {
      stop(
        paste(
          "The following", label, "parents are missing from the supplied gmatrix:",
          paste(utils::head(missing_ids, 10), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(g_use[parent_ids, parent_ids, drop = FALSE])
  }

  if (is.null(geno_data)) {
    stop(
      paste(
        "Provide either", paste0(label, "_gmatrix"),
        "or a parent genotype matrix that can be converted to a gmatrix."
      ),
      call. = FALSE
    )
  }

  geno_use <- as.matrix(geno_data)
  if (is.null(rownames(geno_use))) {
    stop(paste("The", label, "genotype matrix must have row names equal to parent IDs."), call. = FALSE)
  }
  missing_ids <- setdiff(parent_ids, rownames(geno_use))
  if (length(missing_ids)) {
    stop(
      paste(
        "The following", label, "parents are missing from the supplied genotype matrix:",
        paste(utils::head(missing_ids, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (is.null(gmatrix_method) || !nzchar(gmatrix_method)) {
    stop("Provide gmatrix_method when using genotype data for hybrid ASReml-R.", call. = FALSE)
  }

  geno_use <- geno_use[parent_ids, , drop = FALSE]
  g_out <- grm_calculation(
    geno_clean = geno_use, method = gmatrix_method, ploidy = ploidy
  )
  rownames(g_out) <- parent_ids
  colnames(g_out) <- parent_ids
  g_out
}

gp_hybrid_sca_kernel <- function(female_levels,
                                 male_levels,
                                 female_gmatrix,
                                 male_gmatrix,
                                 cross_levels) {
  female_levels <- as.character(female_levels)
  male_levels <- as.character(male_levels)
  cross_levels <- as.character(cross_levels)

  cross_df <- do.call(
    rbind,
    strsplit(cross_levels, "__x__", fixed = TRUE)
  )
  if (ncol(cross_df) != 2L) {
    stop("Hybrid cross labels could not be split into female and male parent IDs.", call. = FALSE)
  }
  female_idx <- match(cross_df[, 1], female_levels)
  male_idx <- match(cross_df[, 2], male_levels)
  if (any(is.na(female_idx)) || any(is.na(male_idx))) {
    stop("Hybrid cross labels did not map cleanly to female and male parent levels.", call. = FALSE)
  }

  kf <- female_gmatrix[female_idx, female_idx, drop = FALSE]
  km <- male_gmatrix[male_idx, male_idx, drop = FALSE]
  k_sca <- kf * km
  rownames(k_sca) <- cross_levels
  colnames(k_sca) <- cross_levels
  k_sca
}

# Map ASReml fixed-effect solutions onto model.matrix columns. model.matrix
# names a factor level as `<term><level>` (e.g. "as.factor(Env)North") while
# ASReml names it `<term>_<level>` (e.g. "as.factor(Env)_North"), so a plain
# name match kept only the intercept and dropped every environment effect.
gp_hybrid_asreml_match_fixed_effects <- function(mm_fixed, fixed_terms, fixed_vec) {
  norm <- function(x) gsub("[`[:space:]]", "", x)
  beta <- stats::setNames(rep(0, ncol(mm_fixed)), colnames(mm_fixed))
  if (!length(fixed_vec)) {
    return(beta)
  }
  fixed_vec <- fixed_vec[is.finite(fixed_vec)]
  asreml_names <- norm(names(fixed_vec))
  term_labels <- norm(attr(fixed_terms, "term.labels"))
  assign_idx <- attr(mm_fixed, "assign")
  for (j in seq_along(beta)) {
    col <- norm(colnames(mm_fixed)[j])
    candidates <- col
    if (!is.null(assign_idx) && assign_idx[j] > 0L) {
      term <- term_labels[assign_idx[j]]
      if (startsWith(col, term)) {
        candidates <- c(candidates, paste0(term, "_", substring(col, nchar(term) + 1L)))
      }
    }
    hit <- match(candidates, asreml_names)
    hit <- hit[!is.na(hit)]
    if (length(hit)) {
      beta[j] <- fixed_vec[[hit[1L]]]
    }
  }
  beta
}

gp_hybrid_asreml_random_coefs <- function(model) {
  # Small testable seam over summary(asreml, coef=TRUE)$coef.random so the
  # warning paths in gp_hybrid_asreml_extract_random_effect() can be
  # exercised without depending on the asreml package being loaded.
  summary(model, coef = TRUE)$coef.random
}

gp_hybrid_asreml_extract_random_effect <- function(model,
                                                   factor_name,
                                                   inv_name,
                                                   id_col_name,
                                                   effect_col_name) {
  # Defensive fallback chain. The hybrid prediction path doesn't go through
  # asreml::predict() at all (so it's immune to the workspace/classify
  # failures the single-trait predict-or-extract helper guards against). It
  # instead reads random effects directly out of mod$coefficients via
  # summary(mod, coef=TRUE)$coef.random and matches rows by vm() label.
  # Each "return(data.frame(...))" below is a graceful degradation point:
  # the caller sees an empty result, sets the corresponding component
  # (Female_GCA / Male_GCA / SCA_effect) to NA, and the hybrid
  # Predicted_value falls back to the fixed-effect intercept alone.
  # The warnings below make that silent degradation visible so users know
  # which component disappeared and why.
  error_reported <- FALSE
  random_coefs <- tryCatch(gp_hybrid_asreml_random_coefs(model),
                           error = function(e) {
                             warning(
                               "Hybrid ASReml extraction for ", effect_col_name,
                               ": summary(model, coef=TRUE) failed (",
                               conditionMessage(e),
                               "); this component will be set to NA in the prediction.",
                               call. = FALSE)
                             error_reported <<- TRUE
                             NULL
                           })
  if (is.null(random_coefs)) {
    # Avoid a second warning when the tryCatch already reported one above.
    if (!error_reported) {
      warning(
        "Hybrid ASReml extraction for ", effect_col_name,
        ": summary(model, coef=TRUE)$coef.random returned NULL; ",
        "this component will be set to NA.",
        call. = FALSE)
    }
    return(data.frame(stringsAsFactors = FALSE))
  }
  if (!nrow(random_coefs)) {
    warning(
      "Hybrid ASReml extraction for ", effect_col_name,
      ": coef.random has no rows; this component will be set to NA.",
      call. = FALSE)
    return(data.frame(stringsAsFactors = FALSE))
  }

  coef_rows <- rownames(random_coefs)
  solution_col <- c("solution", "effect")[c("solution", "effect") %in% colnames(random_coefs)][1]
  if (is.na(solution_col) || is.null(solution_col)) {
    warning(
      "Hybrid ASReml extraction for ", effect_col_name,
      ": neither 'solution' nor 'effect' column found in coef.random; ",
      "this component will be set to NA.",
      call. = FALSE)
    return(data.frame(stringsAsFactors = FALSE))
  }
  se_col <- c("std.error", "Standard_error")[c("std.error", "Standard_error") %in% colnames(random_coefs)][1]
  prefix <- paste0("vm\\(", pp_escape_regex(factor_name), "\\s*,\\s*", pp_escape_regex(inv_name), "\\)_")
  matched_rows <- coef_rows[grepl(prefix, coef_rows)]
  if (!length(matched_rows)) {
    warning(
      "Hybrid ASReml extraction for ", effect_col_name,
      ": no coef.random rows matched vm(", factor_name, ", ", inv_name,
      "); this component will be set to NA.",
      call. = FALSE)
    return(data.frame(stringsAsFactors = FALSE))
  }

  ids <- sub(prefix, "", matched_rows)
  out <- data.frame(
    stringsAsFactors = FALSE,
    id_value = ids,
    effect_value = as.numeric(random_coefs[matched_rows, solution_col])
  )
  names(out) <- c(id_col_name, effect_col_name)
  out$Standard_error <- if (!is.na(se_col) && !is.null(se_col)) as.numeric(random_coefs[matched_rows, se_col]) else NA_real_
  out$Prediction_error_variance <- out$Standard_error^2
  out
}

gp_hybrid_asreml_gaussian_model <- function(pheno_object,
                                            response,
                                            gen_name,
                                            female_parent,
                                            male_parent,
                                            heter_groups = NULL,
                                            gmatrix = NULL,
                                            female_gmatrix = NULL,
                                            male_gmatrix = NULL,
                                            geno_data = NULL,
                                            female_geno_data = NULL,
                                            male_geno_data = NULL,
                                            gmatrix_method = NULL,
                                            ploidy = "auto",
                                            include_sca = TRUE,
                                            inverse = TRUE,
                                            epsilon = 1e-6,
                                            engine = "asreml",
                                            workspace = 1e08,
                                            maxit = 50) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (!requireNamespace(engine, quietly = TRUE)) {
    stop("You need to install asreml-R to use hybrid ASReml-R.", call. = FALSE)
  }

  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  needed_cols <- c(gen_name, female_parent, male_parent, response)
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    needed_cols <- c(needed_cols, heter_groups)
  }
  missing_cols <- setdiff(needed_cols, names(ph))
  if (length(missing_cols)) {
    stop(paste("Hybrid phenotype data is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (is.null(heter_groups) || !nzchar(heter_groups)) {
    if (anyDuplicated(ph[[gen_name]])) {
      stop("Hybrid ASReml-R currently requires one phenotype row per hybrid ID unless heter_groups is supplied for multi-environment data.", call. = FALSE)
    }
  } else if (anyDuplicated(gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups))) {
    stop("Hybrid ASReml-R currently requires one phenotype row per hybrid-by-environment combination.", call. = FALSE)
  }

  ph[[female_parent]] <- as.character(ph[[female_parent]])
  ph[[male_parent]] <- as.character(ph[[male_parent]])
  ph[[gen_name]] <- as.character(ph[[gen_name]])
  ph$HybridCross <- paste(ph[[female_parent]], ph[[male_parent]], sep = "__x__")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    ph[[heter_groups]] <- as.factor(ph[[heter_groups]])
  }

  female_levels <- unique(ph[[female_parent]])
  male_levels <- unique(ph[[male_parent]])
  cross_levels <- unique(ph$HybridCross)

  female_g <- gp_hybrid_prepare_parent_gmatrix(
    parent_ids = female_levels,
    gmatrix = female_gmatrix %||% gmatrix,
    geno_data = if (!is.null(female_geno_data)) {
      female_geno_data
    } else if (is.null(female_gmatrix) && is.null(gmatrix)) {
      geno_data
    } else {
      NULL
    },
    gmatrix_method = gmatrix_method,
    ploidy = ploidy,
    label = "female"
  )
  male_g <- gp_hybrid_prepare_parent_gmatrix(
    parent_ids = male_levels,
    gmatrix = male_gmatrix %||% gmatrix,
    geno_data = if (!is.null(male_geno_data)) {
      male_geno_data
    } else if (is.null(male_gmatrix) && is.null(gmatrix)) {
      geno_data
    } else {
      NULL
    },
    gmatrix_method = gmatrix_method,
    ploidy = ploidy,
    label = "male"
  )

  sca_g <- NULL
  if (isTRUE(include_sca)) {
    sca_g <- gp_hybrid_sca_kernel(
      female_levels = female_levels,
      male_levels = male_levels,
      female_gmatrix = female_g,
      male_gmatrix = male_g,
      cross_levels = cross_levels
    )
  }

  ph[[female_parent]] <- factor(ph[[female_parent]], levels = female_levels)
  ph[[male_parent]] <- factor(ph[[male_parent]], levels = male_levels)
  ph$HybridCross <- factor(ph$HybridCross, levels = cross_levels)

  female_inv_name <- paste0("female_gca_inv_", as.integer(Sys.getpid()))
  male_inv_name <- paste0("male_gca_inv_", as.integer(Sys.getpid()))
  assign(female_inv_name, compute_inverse_and_sparse(kernel = female_g, epsilon = epsilon, inverse = inverse), envir = .GlobalEnv)
  assign(male_inv_name, compute_inverse_and_sparse(kernel = male_g, epsilon = epsilon, inverse = inverse), envir = .GlobalEnv)
  cleanup_names <- c(female_inv_name, male_inv_name)
  if (!is.null(sca_g)) {
    sca_inv_name <- paste0("hybrid_sca_inv_", as.integer(Sys.getpid()))
    assign(sca_inv_name, compute_inverse_and_sparse(kernel = sca_g, epsilon = epsilon, inverse = inverse), envir = .GlobalEnv)
    cleanup_names <- c(cleanup_names, sca_inv_name)
  } else {
    sca_inv_name <- NULL
  }
  on.exit(remove_from_global(cleanup_names), add = TRUE)

  asreml_fn <- getFromNamespace("asreml", engine)
  update_fn <- getFromNamespace("update.asreml", engine)
  asreml_options_fn <- tryCatch(getFromNamespace("asreml.options", engine), error = function(e) NULL)
  if (!is.null(asreml_options_fn)) {
    old_ai_sing <- tryCatch(asreml_options_fn()$ai.sing, error = function(e) NULL)
    try(asreml_options_fn(ai.sing = TRUE), silent = TRUE)
    if (!is.null(old_ai_sing)) {
      on.exit(try(asreml_options_fn(ai.sing = old_ai_sing), silent = TRUE), add = TRUE)
    }
  }

  random_terms <- c(
    paste0("vm(", female_parent, ",", female_inv_name, ")"),
    paste0("vm(", male_parent, ",", male_inv_name, ")")
  )
  if (!is.null(sca_inv_name)) {
    random_terms <- c(random_terms, paste0("vm(HybridCross,", sca_inv_name, ")"))
  }

  fixed_formula <- if (!is.null(heter_groups) && nzchar(heter_groups)) {
    stats::as.formula(
      paste(
        response,
        "~",
        paste0("as.factor(`", heter_groups, "`)")
      )
    )
  } else {
    stats::as.formula(paste(response, "~ 1"))
  }
  fit_args <- list(
    fixed = fixed_formula,
    random = stats::as.formula(paste("~", paste(random_terms, collapse = " + "))),
    data = ph,
    na.action = list(x = "include", y = "include"),
    workspace = workspace,
    maxit = maxit,
    ai.sing = TRUE
  )

  mod <- tryCatch(
    do.call(asreml_fn, fit_args),
    error = function(e) {
      stop(paste("Hybrid ASReml-R fit failed:", conditionMessage(e)), call. = FALSE)
    }
  )

  for (i in seq_len(3L)) {
    if (!isTRUE(mod$converge)) {
      mod <- update_fn(mod, ai.sing = TRUE)
    } else {
      break
    }
  }
  if (!isTRUE(mod$converge)) {
    warning(
      "Hybrid ASReml-R model did not converge after 3 update() calls; ",
      "predictions and variance components may be unreliable. Consider a larger maxit.",
      call. = FALSE
    )
  }

  female_eff <- gp_hybrid_asreml_extract_random_effect(
    model = mod,
    factor_name = female_parent,
    inv_name = female_inv_name,
    id_col_name = female_parent,
    effect_col_name = "Female_GCA"
  )
  male_eff <- gp_hybrid_asreml_extract_random_effect(
    model = mod,
    factor_name = male_parent,
    inv_name = male_inv_name,
    id_col_name = male_parent,
    effect_col_name = "Male_GCA"
  )
  sca_eff <- if (!is.null(sca_inv_name)) {
    gp_hybrid_asreml_extract_random_effect(
      model = mod,
      factor_name = "HybridCross",
      inv_name = sca_inv_name,
      id_col_name = "HybridCross",
      effect_col_name = "SCA_effect"
    )
  } else {
    data.frame(stringsAsFactors = FALSE)
  }

    fixed_coefs <- tryCatch(coef(mod)$fixed, error = function(e) NULL)
    fixed_vec <- numeric(0)
    if (!is.null(fixed_coefs)) {
      if (is.matrix(fixed_coefs) || is.data.frame(fixed_coefs)) {
        fixed_vec <- as.numeric(fixed_coefs[, 1])
        names(fixed_vec) <- rownames(fixed_coefs)
      } else {
        fixed_vec <- as.numeric(fixed_coefs)
        names(fixed_vec) <- names(fixed_coefs)
      }
    }
    fixed_terms <- stats::delete.response(stats::terms(fixed_formula))
    mm_fixed <- stats::model.matrix(fixed_terms, data = ph)
    beta_fixed <- gp_hybrid_asreml_match_fixed_effects(mm_fixed, fixed_terms, fixed_vec)
    fixed_part <- as.numeric(mm_fixed %*% beta_fixed)

  out_cols <- c(gen_name, female_parent, male_parent, response, "HybridCross")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    out_cols <- c(out_cols, heter_groups)
  }
  out <- ph[, out_cols, drop = FALSE]
  names(out)[names(out) == response] <- "Observed_value"
  out$Train_Test_Label <- ifelse(is.na(out$Observed_value), "Test", "Train")
  if (nrow(female_eff)) {
    out$Female_GCA <- female_eff$Female_GCA[match(out[[female_parent]], female_eff[[female_parent]])]
  } else {
    out$Female_GCA <- NA_real_
  }
  if (nrow(male_eff)) {
    out$Male_GCA <- male_eff$Male_GCA[match(out[[male_parent]], male_eff[[male_parent]])]
  } else {
    out$Male_GCA <- NA_real_
  }
  if (nrow(sca_eff)) {
    out$SCA_effect <- sca_eff$SCA_effect[match(out$HybridCross, sca_eff$HybridCross)]
  } else {
    out$SCA_effect <- NA_real_
  }
  out$Prediction_intercept <- fixed_part
  out$Predicted_value <- fixed_part +
      ifelse(is.na(out$Female_GCA), 0, out$Female_GCA) +
      ifelse(is.na(out$Male_GCA), 0, out$Male_GCA) +
      ifelse(is.na(out$SCA_effect), 0, out$SCA_effect)
  female_se <- if (nrow(female_eff) && "Standard_error" %in% names(female_eff)) female_eff$Standard_error[match(out[[female_parent]], female_eff[[female_parent]])] else rep(NA_real_, nrow(out))
  male_se <- if (nrow(male_eff) && "Standard_error" %in% names(male_eff)) male_eff$Standard_error[match(out[[male_parent]], male_eff[[male_parent]])] else rep(NA_real_, nrow(out))
  sca_se <- if (nrow(sca_eff) && "Standard_error" %in% names(sca_eff)) sca_eff$Standard_error[match(out$HybridCross, sca_eff$HybridCross)] else rep(NA_real_, nrow(out))
  # NOTE (rigor): the hybrid prediction is f_hat + m_hat + s_hat (female GCA +
  # male GCA + SCA). Its exact prediction error variance is
  #   Var(f+m+s) = Var(f) + Var(m) + Var(s)
  #                + 2 Cov(f,m) + 2 Cov(f,s) + 2 Cov(m,s),
  # where the covariances come from the off-diagonal blocks of the mixed-model
  # coefficient (C-inverse) matrix. The expression below keeps only the diagonal
  # Var terms and DROPS the cross-covariances, i.e. it treats the three BLUP
  # components as independent.
  #
  # This was investigated end-to-end against live ASReml fits (2026-05-27):
  #   * predict(mod, classify = "HybridCross", present = c(female, male,
  #     "HybridCross"))$pvals reconstructs f_hat+m_hat+s_hat EXACTLY and its
  #     std.error is the covariance-correct joint PEV -- but ONLY for OBSERVED
  #     crosses. Crosses with no data record are returned "Aliased" (NA), i.e.
  #     exactly the unobserved hybrids that true prediction targets.
  #   * asreml-R does not expose the full coefficient C-inverse: mod$vcoeff gives
  #     only the diagonal coefficient variances (== std.error^2 already used), and
  #     mod$aom gives residual/BLUP diagnostics, not the cross-term blocks. So the
  #     +2 Cov terms are NOT obtainable for unobserved crosses via the public API.
  # Empirically the correction was modest (~2%) and sign-varied (covariances were
  # net negative in the test, so independence slightly OVER-stated the SE), so the
  # documented independence approximation is retained rather than a partial fix
  # that would be exact for observed crosses but not for the unobserved ones.
  out$Standard_error <- sqrt(
    rowSums(cbind(
      ifelse(is.na(female_se), 0, female_se^2),
      ifelse(is.na(male_se), 0, male_se^2),
      ifelse(is.na(sca_se), 0, sca_se^2)
    ))
  )
  out$Standard_error[is.na(female_se) & is.na(male_se) & is.na(sca_se)] <- NA_real_
  out$Prediction_error_variance <- out$Standard_error^2
  rownames(out) <- NULL

  vc <- tryCatch(asreml_varcomp_table(summary(mod)), error = function(e) NULL)

  list(
    Asreml_model = mod,
    Predicted_value = out,
    predicted_values = out,
    female_gca_effects = female_eff,
    male_gca_effects = male_eff,
    sca_effects = sca_eff,
    variance_components = vc,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    diagnostic_plots = gp_hybrid_asreml_diagnostic_plot(out)
  )
}

gp_hybrid_asreml_gaussian_cv <- function(pheno_object,
                                         response,
                                         gen_name,
                                         female_parent,
                                         male_parent,
                                         heter_groups = NULL,
                                         cross_validation_meth,
                                         eval_metrics,
                                         gmatrix = NULL,
                                         female_gmatrix = NULL,
                                         male_gmatrix = NULL,
                                         geno_data = NULL,
                                         female_geno_data = NULL,
                                         male_geno_data = NULL,
                                         gmatrix_method = NULL,
                                         ploidy = "auto",
                                         include_sca = TRUE,
                                         inverse = TRUE,
                                         epsilon = 1e-6,
                                         engine = "asreml",
                                         workspace = 1e08,
                                         maxit = 50,
                                         nfolds = 5L,
                                         random_state = 123L,
                                         replication = 1L) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  all_preds <- list()
  all_eval <- list()
  all_raw <- list()

  for (rep_i in seq_len(replication)) {
    scenarios <- gp_hybrid_build_cv_scenarios(
      pheno_object = ph,
      response = response,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      random_state = random_state,
      replication = rep_i
    )
    if (!length(scenarios)) {
      next
    }

    for (sc in scenarios) {
      scenario_rows <- gp_hybrid_cv_scenario_rows(sc, ph, gen_name)
      tst <- scenario_rows$test
      if (!length(tst) || !length(scenario_rows$train)) next
      ph_cv <- gp_hybrid_cv_mask_responses(ph, response, scenario_rows$train)

      fit <- gp_hybrid_asreml_gaussian_model(
        pheno_object = ph_cv,
        response = response,
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        gmatrix = gmatrix,
        female_gmatrix = female_gmatrix,
        male_gmatrix = male_gmatrix,
        geno_data = geno_data,
        female_geno_data = female_geno_data,
        male_geno_data = male_geno_data,
        gmatrix_method = gmatrix_method,
        ploidy = ploidy,
        include_sca = include_sca,
        inverse = inverse,
        epsilon = epsilon,
        engine = engine,
        workspace = workspace,
        maxit = maxit
      )

      pred_test <- gp_hybrid_fold_test_predictions(
        fit$predicted_values, ph, tst, gen_name = gen_name, heter_groups = heter_groups
      )
      pred_test <- gp_hybrid_match_observed_values(
        pred_df = pred_test,
        ph = ph,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      pred_test$cv_scenario <- sc$scenario
      pred_test$fold <- sc$fold
      pred_test$rep <- rep_i

      eval_df <- gp_hybrid_asreml_cv_metric_table(
        pred_df = pred_test,
        eval_metrics = eval_metrics,
        response_family = "gaussian"
      )

      all_preds[[length(all_preds) + 1L]] <- pred_test
      all_eval[[length(all_eval) + 1L]] <- eval_df
      all_raw[[length(all_raw) + 1L]] <- list(
        trait = response,
        rep = rep_i,
        model = "GBLUP",
        eval_metrics_reps = eval_df,
        ypred_cv_Reps_all = pred_test,
        cv_info = list(
          method = cross_validation_meth,
          scenario = sc$scenario,
          fold = sc$fold,
          n_test = nrow(pred_test),
          n_train = sum(!is.na(ph_cv[[response]])),
          n_test_rows = nrow(pred_test)
        )
      )
    }
  }

  pred_df <- if (length(all_preds)) do.call(rbind, all_preds) else data.frame()
  eval_df <- if (length(all_eval)) do.call(rbind, all_eval) else data.frame()
  processed <- gp_hybrid_asreml_cv_process(
    pred_df = pred_df,
    eval_df = eval_df,
    response = response,
    model_type = "GBLUP"
  )

  list(
    predicted_values = pred_df,
    hybrid_cv_metrics = eval_df,
    cv_results_processed = processed,
    cv_results_raw = all_raw,
    diagnostic_plots = processed$hybrid_cv_plot
  )
}
