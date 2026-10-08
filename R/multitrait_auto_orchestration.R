gp_multitrait_auto_flag_names <- function() {
  c(
    "multi_trait_asreml",
    "multi_trait_gp",
    "multi_trait_bayes",
    "multi_trait_ml",
    "multi_trait_dl"
  )
}

gp_multitrait_auto_cv_method <- function(cross_validation_meth,
                                         is_met = FALSE) {
  method <- as.character(cross_validation_meth %||% character())
  if (length(method) != 1L || is.na(method) || !nzchar(trimws(method))) {
    choices <- if (isTRUE(is_met)) {
      "CV0, CV1, CV2, Repeated_CV0, Repeated_CV1, or Repeated_CV2"
    } else {
      "K-Folds or Repeated_K-Folds"
    }
    stop(
      "Automatic multi-trait cross-validation requires exactly one CV method: ",
      choices, ".",
      call. = FALSE
    )
  }

  token <- normalize_cv_token(method)
  canonical <- if (isTRUE(is_met)) {
    switch(
      token,
      cv0 = "CV0",
      cv1 = "CV1",
      cv2 = "CV2",
      repeated_cv0 = "Repeated_CV0",
      repeated_cv1 = "Repeated_CV1",
      repeated_cv2 = "Repeated_CV2",
      NULL
    )
  } else {
    switch(
      token,
      k_folds = "K-Folds",
      repeated_k_folds = "Repeated_K-Folds",
      NULL
    )
  }
  if (is.null(canonical)) {
    choices <- if (isTRUE(is_met)) {
      "CV0, CV1, CV2, Repeated_CV0, Repeated_CV1, and Repeated_CV2"
    } else {
      "K-Folds and Repeated_K-Folds"
    }
    stop(
      "Automatic ", if (isTRUE(is_met)) "MT-MET" else "multi-trait",
      " cross-validation supports ", choices, " only; received `",
      method, "` (normalized token `",
      token, "`).",
      call. = FALSE
    )
  }
  canonical
}

gp_multitrait_auto_inventory <- function(cross_validation = FALSE,
                                         is_met = FALSE) {
  if (isTRUE(is_met)) {
    if (isTRUE(cross_validation)) {
      return(list(
        multi_trait_asreml = gp_multitrait_asreml_supported_models(),
        multi_trait_gp = gp_multitrait_gp_supported_models(),
        multi_trait_bayes = gp_multitrait_bayes_supported_models()
      ))
    }
    return(list(
      multi_trait_asreml = gp_multitrait_asreml_supported_models(),
      multi_trait_gp = gp_multitrait_gp_supported_models(),
      multi_trait_bayes = gp_multitrait_bayes_supported_models()
    ))
  }

  if (isTRUE(cross_validation)) {
    return(list(
      multi_trait_asreml = gp_multitrait_asreml_supported_models(),
      multi_trait_gp = gp_multitrait_gp_supported_models(),
      multi_trait_bayes = gp_multitrait_bayes_supported_models(),
      multi_trait_ml = gp_multitrait_ml_supported_models(),
      multi_trait_dl = gp_multitrait_dl_supported_models()
    ))
  }

  list(
    multi_trait_asreml = gp_multitrait_asreml_supported_models(),
    multi_trait_gp = gp_multitrait_gp_supported_models(),
    multi_trait_bayes = gp_multitrait_bayes_supported_models(),
    multi_trait_ml = gp_multitrait_ml_supported_models(),
    multi_trait_dl = gp_multitrait_dl_supported_models()
  )
}

gp_multitrait_auto_plan <- function(models,
                                    cross_validation = FALSE,
                                    is_met = FALSE) {
  requested <- as.character(models %||% character())
  requested <- requested[!is.na(requested) & nzchar(trimws(requested))]
  if (!length(requested)) {
    stage <- if (isTRUE(cross_validation)) "GS_model_cv" else "GS_model"
    stop(
      "multi_trait = TRUE requires at least one model in `", stage, "`.",
      call. = FALSE
    )
  }

  canonical <- gp_canonicalize_supported_model_names(requested)
  keep <- !duplicated(canonical)
  requested <- requested[keep]
  canonical <- canonical[keep]
  inventory <- gp_multitrait_auto_inventory(
    cross_validation = cross_validation,
    is_met = is_met
  )
  allowed <- unique(unlist(inventory, use.names = FALSE))

  if (!length(allowed)) {
    stop(
      "No automatic joint multi-trait multi-environment model supports the ",
      "requested cross-validation configuration.",
      call. = FALSE
    )
  }

  unsupported <- !canonical %in% allowed
  if (any(unsupported)) {
    scope <- if (isTRUE(is_met)) {
      "joint multi-trait multi-environment prediction"
    } else if (isTRUE(cross_validation)) {
      "joint multi-trait single-environment cross-validation"
    } else {
      "joint multi-trait single-environment prediction"
    }
    stop(
      scope,
      " does not support: ",
      paste(requested[unsupported], collapse = ", "),
      ". Supported models are: ",
      paste(gp_display_supported_model_names(allowed), collapse = ", "),
      ". No requested model was fitted.",
      call. = FALSE
    )
  }

  route <- vapply(canonical, function(model) {
    hits <- names(inventory)[vapply(
      inventory,
      function(supported) model %in% supported,
      logical(1L)
    )]
    if (length(hits) != 1L) {
      stop(
        "PredictProR could not resolve one multi-trait route for model `",
        gp_public_model_label(model),
        "`.",
        call. = FALSE
      )
    }
    hits[[1L]]
  }, character(1L))

  data.frame(
    requested_model = requested,
    model = vapply(canonical, gp_public_model_label, character(1L)),
    canonical_model = canonical,
    route = route,
    stringsAsFactors = FALSE
  )
}

gp_multitrait_auto_is_met <- function(args) {
  pheno <- args$pheno_data %||% args$pheno_data_train %||% args$pheno_data_test
  if (is.null(pheno) || is.null(args$gen_name) || !length(args$response)) {
    return(FALSE)
  }
  tryCatch(
    gp_is_multi_environment_trait_panel(
      pheno_data = pheno,
      gen_name = args$gen_name,
      heter_groups = args$heter_groups,
      response = args$response,
      response_family = args$response_family
    ),
    error = function(e) FALSE
  )
}

gp_multitrait_auto_child_args <- function(args,
                                          plan_row,
                                          cross_validation = FALSE,
                                          cv_evaluation_only = TRUE,
                                          final_prediction = FALSE) {
  child <- args
  child$multi_trait <- FALSE
  for (flag in gp_multitrait_auto_flag_names()) {
    child[[flag]] <- identical(as.character(plan_row$route[[1L]]), flag)
  }
  child$met_ml_dl <- FALSE
  child$cross_validation <- isTRUE(cross_validation)
  child$cv_evaluation_only <- isTRUE(cv_evaluation_only)
  if (isTRUE(cross_validation)) {
    child$GS_model <- NULL
    child$GS_model_cv <- as.character(plan_row$model[[1L]])
  } else {
    child$GS_model <- as.character(plan_row$model[[1L]])
    child$GS_model_cv <- NULL
  }

  if (identical(as.character(plan_row$route[[1L]]), "multi_trait_asreml")) {
    child$engine <- "asreml"
  }

  kernel_route <- as.character(plan_row$route[[1L]]) %in%
    c("multi_trait_asreml", "multi_trait_gp", "multi_trait_bayes")
  supplied_kernel <- any(vapply(
    list(
      child$gmatrix, child$gkernel, child$kernel_list,
      child$omic1_kernel, child$omic2_kernel, child$omic3_kernel
    ),
    Negate(is.null),
    logical(1L)
  ))
  if (isTRUE(kernel_route) && !isTRUE(supplied_kernel) &&
      !is.null(child$geno_data) && is.null(child$gmatrix_method)) {
    child$gmatrix_method <- "Yang"
  }

  suffix <- if (isTRUE(final_prediction)) "final" else if (isTRUE(cross_validation)) "cv" else "prediction"
  base_name <- as.character(child$plot_filename %||% "multi_trait")[[1L]]
  # results_handling() appends the exact public model label uniformly for all
  # workflows. Keep this value limited to the user's context and run stage so
  # automatic child routes do not duplicate or alter the model spelling.
  child$plot_filename <- paste(base_name, suffix, sep = "_")
  child
}

gp_multitrait_auto_bind_rows <- function(parts) {
  parts <- Filter(function(x) is.data.frame(x) && nrow(x), parts)
  if (!length(parts)) {
    return(data.frame())
  }
  columns <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(columns, names(x))
    for (column in missing) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

gp_multitrait_auto_metric_summary <- function(runs, plan) {
  parts <- lapply(seq_along(runs), function(i) {
    processed <- runs[[i]][["cv_results_processed"]] %||% list()
    tab <- processed[["multitrait_metric_summary"]]
    if (!is.data.frame(tab) || !nrow(tab)) return(NULL)
    tab$model <- as.character(plan$model[[i]])
    tab
  })
  gp_multitrait_auto_bind_rows(parts)
}

gp_multitrait_auto_ranking <- function(metric_summary,
                                       metric_for_ranking = "auto",
                                       ranking_tie_breakers = NULL,
                                       response_family = "gaussian") {
  empty <- list(
    primary_metric = NA_character_,
    overall_model_ranking = data.frame(),
    best_models_by_trait = data.frame(),
    selected_model = NULL
  )
  if (!is.data.frame(metric_summary) || !nrow(metric_summary) ||
      !all(c("trait", "model", "metric", "value") %in% names(metric_summary))) {
    return(empty)
  }

  family <- gp_resolve_response_family(response_family)
  primary <- gp_resolve_metric_for_ranking(metric_for_ranking, family)
  tie_breakers <- gp_resolve_ranking_tie_breakers(
    ranking_tie_breakers = ranking_tie_breakers,
    response_family = family,
    metric_for_ranking = primary
  )
  metrics <- unique(c(primary, tie_breakers))
  long <- metric_summary[metric_summary$metric %in% metrics, , drop = FALSE]
  long$value <- suppressWarnings(as.numeric(long$value))
  long <- long[is.finite(long$value), , drop = FALSE]
  if (!nrow(long) || !primary %in% long$metric) {
    empty$primary_metric <- primary
    return(empty)
  }

  trait_long <- stats::aggregate(
    value ~ trait + model + metric,
    data = long,
    FUN = mean,
    na.rm = TRUE
  )
  trait_wide <- stats::reshape(
    trait_long,
    idvar = c("trait", "model"),
    timevar = "metric",
    direction = "wide"
  )
  names(trait_wide) <- sub("^value\\.", "", names(trait_wide))
  best_by_trait <- tryCatch(
    rank_models(trait_wide, primary, tie_breakers = tie_breakers),
    error = function(e) data.frame()
  )

  overall_long <- stats::aggregate(
    value ~ model + metric,
    data = trait_long,
    FUN = mean,
    na.rm = TRUE
  )
  overall <- stats::reshape(
    overall_long,
    idvar = "model",
    timevar = "metric",
    direction = "wide"
  )
  names(overall) <- sub("^value\\.", "", names(overall))
  usable_ties <- tie_breakers[tie_breakers %in% names(overall)]
  order_columns <- c(primary, usable_ties)
  order_keys <- lapply(order_columns, function(metric) {
    values <- suppressWarnings(as.numeric(overall[[metric]]))
    if (identical(gp_eval_metric_direction(metric), "maximize")) -values else values
  })
  order_keys[[length(order_keys) + 1L]] <- tolower(as.character(overall$model))
  ord <- do.call(order, c(order_keys, list(na.last = TRUE)))
  overall <- overall[ord, , drop = FALSE]
  overall$rank <- seq_len(nrow(overall))
  overall$selection_scope <- "mean_primary_metric_across_traits"
  rownames(overall) <- NULL

  list(
    primary_metric = primary,
    overall_model_ranking = overall,
    best_models_by_trait = best_by_trait,
    selected_model = if (nrow(overall)) as.character(overall$model[[1L]]) else NULL
  )
}

gp_multitrait_auto_result <- function(runs,
                                      plan,
                                      args,
                                      final_prediction = NULL) {
  model_names <- as.character(plan$model)
  names(runs) <- model_names
  is_cv <- isTRUE(args$cross_validation)
  metric_summary <- if (is_cv) gp_multitrait_auto_metric_summary(runs, plan) else data.frame()
  ranking <- if (is_cv) {
    gp_multitrait_auto_ranking(
      metric_summary = metric_summary,
      metric_for_ranking = args$metric_for_ranking,
      ranking_tie_breakers = args$ranking_tie_breakers,
      response_family = args$response_family
    )
  } else {
    list(
      primary_metric = NA_character_,
      overall_model_ranking = data.frame(),
      best_models_by_trait = data.frame(),
      selected_model = NULL
    )
  }

  model_results <- stats::setNames(lapply(runs, function(x) {
    x[["model_results"]] %||% x
  }), model_names)
  cv_by_model <- stats::setNames(lapply(runs, `[[`, "cv_results_processed"), model_names)
  cv_raw <- stats::setNames(lapply(runs, `[[`, "cv_results_raw"), model_names)
  processed <- if (is_cv) {
    list(
      multitrait_metric_summary = metric_summary,
      overall_model_ranking = ranking$overall_model_ranking,
      best_models_by_trait = ranking$best_models_by_trait,
      primary_metric = ranking$primary_metric,
      by_model = cv_by_model
    )
  } else {
    NULL
  }

  preprocess_summary <- gp_shared_preprocess_summary(
    args$.predictpror_preprocess_cache %||% NULL
  )
  # model_results is the final prediction when the selected model was refitted
  # (as in every other workflow: res$model_results$predicted_values); the
  # per-model runs stay in model_results_by_model.
  final_model_results <- if (is.list(final_prediction)) final_prediction[["model_results"]] else NULL
  out <- list(
    model_results = final_model_results %||% model_results,
    model_results_by_model = model_results,
    cv_results_raw = if (is_cv) cv_raw else NULL,
    cv_results_processed = processed,
    multi_trait_runs = runs,
    orchestration_plan = plan,
    selected_model = ranking$selected_model,
    final_prediction = final_prediction,
    run_metadata = list(
      task_unit = if (is_cv) "multi_trait_auto_cv" else "multi_trait_auto_prediction",
      model_count = nrow(plan),
      models = model_names,
      automatic_route_detection = TRUE,
      cross_validation = is_cv,
      cv_evaluation_only = isTRUE(args$cv_evaluation_only),
      final_prediction_performed = !is.null(final_prediction),
      selection_scope = if (is_cv) "mean_primary_metric_across_traits" else NA_character_,
      shared_preprocessing = isTRUE(preprocess_summary$enabled),
      preprocessing_stages_built = sum(preprocess_summary$stages_built),
      preprocessing_stage_reuse_hits = sum(preprocess_summary$stage_reuse_hits),
      preprocessing_detail = preprocess_summary
    )
  )
  class(out) <- c("PredictProR_multi_trait_result", "list")
  out
}

gp_model_execute_multitrait_orchestrate <- function(args,
                                                    execute_fn = model_execute) {
  explicit_flags <- vapply(
    gp_multitrait_auto_flag_names(),
    function(flag) isTRUE(args[[flag]]),
    logical(1L)
  )
  if (any(explicit_flags)) {
    stop(
      "When `multi_trait = TRUE`, leave multi_trait_asreml, multi_trait_gp, ",
      "multi_trait_bayes, multi_trait_ml, and multi_trait_dl set to FALSE. ",
      "PredictProR selects the protected route separately for each model.",
      call. = FALSE
    )
  }
  hybrid_flags <- c(
    "hybrid_asreml", "hybrid_bayes", "hybrid_gp", "hybrid_ml", "hybrid_dl"
  )
  if (isTRUE(args$met_ml_dl) || any(vapply(
    hybrid_flags,
    function(flag) isTRUE(args[[flag]]),
    logical(1L)
  ))) {
    stop(
      "`multi_trait = TRUE` is a separate high-level workflow; do not combine ",
      "it with met_ml_dl or hybrid workflow flags.",
      call. = FALSE
    )
  }
  response_input <- as.character(args$response %||% character())
  response_input <- response_input[!is.na(response_input) & nzchar(response_input)]
  response <- unique(response_input)
  if (length(response) < 2L) {
    stop(
      "multi_trait = TRUE requires at least two response column names.",
      call. = FALSE
    )
  }
  if (length(response) < length(response_input) && isTRUE(args$message)) {
    message(
      "`response` contained duplicate trait names; automatic multi-trait ",
      "execution will use: ", paste(response, collapse = ", "), "."
    )
  }
  # Write the sanitized response vector back before constructing any child
  # route.  This prevents duplicated traits from being fitted repeatedly and
  # from appearing repeatedly in missing-response diagnostics.
  args$response <- response
  family <- gp_resolve_response_family_for_responses(
    pheno_data = args$pheno_data,
    response = response,
    response_family = args$response_family
  )
  if (!identical(family, "gaussian")) {
    stop(
      "Automatic multi-trait orchestration currently supports gaussian traits only.",
      call. = FALSE
    )
  }
  args$response_family <- family

  is_cv <- isTRUE(args$cross_validation)
  is_met <- gp_multitrait_auto_is_met(args)
  if (is_cv) {
    # Canonicalize the public method once at the shared boundary. Single-
    # environment joint models use genotype K-folds; MT-MET models use the
    # shared CV0/CV1/CV2 scenario assignments.
    args$cross_validation_meth <- gp_multitrait_auto_cv_method(
      args$cross_validation_meth,
      is_met = is_met
    )
    sampling_token <- normalize_cv_token(args$sampling_method %||% "")
    if (identical(sampling_token, "stratified") && isTRUE(args$message)) {
      message(
        if (isTRUE(is_met)) {
          "Joint MT-MET CV uses CV0/CV1/CV2 scenario assignments; "
        } else {
          "Joint multi-trait CV uses unstratified genotype-level folds; "
        },
        "`sampling_method = \"stratified\"` has been replaced by ",
        "`sampling_method = \"unstratified\"`."
      )
    }
    args$sampling_method <- "unstratified"
  }
  models <- if (is_cv) args$GS_model_cv else args$GS_model
  plan <- gp_multitrait_auto_plan(
    models = models,
    cross_validation = is_cv,
    is_met = is_met
  )

  # Every protected child receives the same process-local cache.  The first
  # child prepares an immutable payload for the union of requested model
  # families; later CV children and the selected final-prediction child reuse
  # it.  Internal fold/model parallelism therefore starts from the same single
  # prepared payload already used by ordinary model_execute() workflows.
  shared_models <- unique(as.character(plan$canonical_model))
  shared_requires_asreml <- any(plan$route %in% "multi_trait_asreml")
  shared_cache <- gp_shared_preprocess_cache_new(
    models = shared_models,
    requires_asreml = shared_requires_asreml,
    requires_final_prediction = is_cv && !isTRUE(args$cv_evaluation_only)
  )
  args$.predictpror_preprocess_cache <- shared_cache
  args$.predictpror_shared_models <- shared_models
  args$.predictpror_shared_requires_asreml <- shared_requires_asreml
  args$.predictpror_shared_final_prediction <- is_cv &&
    !isTRUE(args$cv_evaluation_only)

  shared_kernel_route <- any(plan$route %in% c(
    "multi_trait_asreml", "multi_trait_gp", "multi_trait_bayes"
  ))
  supplied_kernel <- any(vapply(
    list(
      args$gmatrix, args$gkernel, args$kernel_list,
      args$omic1_kernel, args$omic2_kernel, args$omic3_kernel
    ),
    Negate(is.null),
    logical(1L)
  ))
  raw_genotype_supplied <- !is.null(args$geno_data) ||
    !is.null(args$train_geno_data)
  if (isTRUE(shared_kernel_route) && !isTRUE(supplied_kernel) &&
      isTRUE(raw_genotype_supplied) && is.null(args$gmatrix_method)) {
    args$gmatrix_method <- "Yang"
  }
  if (nrow(plan) > 1L && isTRUE(args$message)) {
    message(
      "Preparing phenotype, genomic/omic, and kernel inputs once for ",
      nrow(plan),
      " automatic multi-trait model routes; child models reuse the shared payload."
    )
  }

  runs <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    child_args <- gp_multitrait_auto_child_args(
      args = args,
      plan_row = plan[i, , drop = FALSE],
      cross_validation = is_cv,
      cv_evaluation_only = if (is_cv) TRUE else FALSE
    )
    runs[[i]] <- tryCatch(
      do.call(execute_fn, child_args),
      error = function(e) {
        stop(
          "Automatic multi-trait execution failed for model `",
          plan$model[[i]],
          "` through route `",
          plan$route[[i]],
          "`: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  preview <- gp_multitrait_auto_result(runs, plan, args)
  final_prediction <- NULL
  if (is_cv && !isTRUE(args$cv_evaluation_only)) {
    selected <- preview$selected_model
    if (is.null(selected) && nrow(plan) == 1L) {
      selected <- as.character(plan$model[[1L]])
    }
    if (is.null(selected) || !length(selected)) {
      stop(
        "Multi-trait CV completed, but no model could be selected for final prediction. ",
        "Inspect `cv_results_processed$multitrait_metric_summary` or set ",
        "cv_evaluation_only = TRUE.",
        call. = FALSE
      )
    }
    selected_row <- plan[match(selected, plan$model), , drop = FALSE]
    final_args <- gp_multitrait_auto_child_args(
      args = args,
      plan_row = selected_row,
      cross_validation = FALSE,
      cv_evaluation_only = FALSE,
      final_prediction = TRUE
    )
    final_prediction <- tryCatch(
      do.call(execute_fn, final_args),
      error = function(e) {
        stop(
          "Multi-trait CV selected `",
          selected,
          "`, but its final prediction fit failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  gp_multitrait_auto_result(
    runs = runs,
    plan = plan,
    args = args,
    final_prediction = final_prediction
  )
}
