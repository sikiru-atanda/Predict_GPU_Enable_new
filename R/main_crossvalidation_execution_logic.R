
# Main cross-validation execution logic

#' Execute cross-validation for multiple models
#' @param pheno_data Data frame with phenotypic data
#' @param params Configuration list for model hyperparameters and settings
#' @param verbose Logical; if TRUE, enable verbose logging
#' @return List of cross-validation results
models_execute_crossval <- function(pheno_data = NULL,
                                    params = list(),
                                    verbose = FALSE) {
  log_threshold <- logger::log_threshold
  log_info <- logger::log_info
  log_debug <- logger::log_debug
  log_warn <- logger::log_warn
  log_error <- logger::log_error
  INFO <- logger::INFO
  WARN <- logger::WARN

  log_threshold(if (verbose) INFO else WARN)
  log_info("Starting cross-validation")

  # Folds are seeded individually (set.seed per fold); restore the caller's
  # RNG kind and stream afterwards so CV does not move the user's session.
  old_rng_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_rng_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  on.exit({
    log_info("Cleaning up resources")
    stop_daemons()
  }, add = TRUE)

  `%||%` <- function(a, b) if (is.null(a)) b else a
  gp_py_bin <- if (exists("gp_detect_python", mode = "function")) gp_detect_python() else NULL
  init_py <- function() if (exists("gp_init_python_once", mode = "function")) gp_init_python_once(gp_py_bin)
  metric_needs_probabilities <- function(metric, fam) {
    metric <- tolower(metric)
    if (identical(fam, "binary")) {
      return(metric %in% c("log_loss", "brier_score", "ece"))
    }
    if (fam %in% c("ordinal", "multiclass")) {
      return(metric %in% c("log_loss", "ece"))
    }
    FALSE
  }
  class_levels_for_y <- function(y, fam) {
    if (identical(fam, "gaussian")) {
      return(NULL)
    }
    if (is.factor(y) || is.ordered(y)) {
      return(as.character(levels(y)))
    }
    sort(unique(as.character(stats::na.omit(y))))
  }

  # Extract parameters
  test_set <- params$test_set %||% NULL
  response <- params$response %||% NULL
  gen_name <- params$gen_name %||% NULL
  test_size <- params$test_size %||% NULL
  random_state <- params$random_state %||% NULL
  replication <- params$replication %||% 1L
  weights <- params$weights %||% NULL
  selected_raw <- params$selected_raw %||% NULL
  feature_score_metadata <- params$feature_score_metadata %||% NULL
  feature_k_grid <- params$feature_k_grid %||% NULL
  feature_scoring_model <- params$feature_scoring_model %||% NA_character_
  feature_scoring_cv <- params$feature_scoring_cv %||% "fixed"
  feature_scoring_seed <- as.integer(params$feature_scoring_seed %||% random_state %||% 123L)
  feature_ridge_lambda <- params$feature_ridge_lambda %||% 1
  feature_bayes_nIter <- params$feature_bayes_nIter %||% 1500L
  feature_bayes_burnIn <- params$feature_bayes_burnIn %||% 500L
  feature_bayes_thin <- params$feature_bayes_thin %||% 5L
  feature_source_map <- params$feature_source_map %||% NULL
  feature_source_matrices <- params$feature_source_matrices %||% list()
  cross_validation_meth <- params$cross_validation_meth %||% NULL
  sampling_method <- params$sampling_method %||% NULL
  eval_metrics <- params$eval_metrics %||% NULL
  response_family <- params$response_family %||% "gaussian"
  positive_class <- params$positive_class %||% NULL
  GS_model_cv <- params$GS_model_cv %||% NULL
  GS_model_cv_display <- params$GS_model_cv_display %||% NULL
  model_prep_all_bayes_cv <- params$model_prep_all_bayes_cv %||% NULL
  asreml_models_prep_cv <- params$asreml_models_prep_cv %||% NULL
  ml_dat_res <- params$ml_dat_res %||% NULL
  met_ml_dl <- params$met_ml_dl %||% FALSE
  gmatrix_model_ready <- params$gmatrix_model_ready %||% NULL
  omic1_kernel_model_ready <- params$omic1_kernel_model_ready %||% NULL
  omic2_kernel_model_ready <- params$omic2_kernel_model_ready %||% NULL
  omic3_kernel_model_ready <- params$omic3_kernel_model_ready %||% NULL
  kernel_list <- params$kernel_list %||% NULL
  omics_kernel_label <- params$omics_kernel_label %||% list(
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL
  )
  met_kernel_var_explained <- params$met_kernel_var_explained %||% 0.95
  met_kernel_min_ev <- params$met_kernel_min_ev %||% 1e-8
  met_kernel_max_pcs <- params$met_kernel_max_pcs %||% NULL
  lowrank_eps_trace <- params$lowrank_eps_trace %||% 1e-6
  lowrank_max_rank <- params$lowrank_max_rank %||% NULL
  lowrank_jitter <- params$lowrank_jitter %||% NULL
  lowrank_noise_grid <- params$lowrank_noise_grid %||% NULL
  lowrank_kernel_weights <- params$lowrank_kernel_weights %||% NULL
  gp_backend <- params$gp_backend %||% "auto"
  gp_output_level <- params$gp_output_level %||% "predict_only"
  gp_return_se <- params$gp_return_se %||% FALSE
  gp_full_vc <- params$gp_full_vc %||% FALSE
  gp_varcomp_mode <- params$gp_varcomp_mode %||% "mom"
  gp_fa_rank <- params$gp_fa_rank %||% 1L
  gp_prediction_output <- params$gp_prediction_output %||% "test_only"
  gp_factor_cache <- params$gp_factor_cache %||% NULL
  gp_iters <- params$gp_iters %||% NULL
  gp_lr <- params$gp_lr %||% NULL
  gp_exact_fast_cv <- params$gp_exact_fast_cv %||% NULL
  env_similarity <- params$env_similarity %||% NULL
  env_ids <- params$env_ids %||% NULL
  env_covariates <- params$env_covariates %||% NULL
  reaction_norm_feature_qc <- params$reaction_norm_feature_qc %||% TRUE
  kenv_kernel <- params$kenv_kernel %||% "matern32"
  kenv_bandwidth <- params$kenv_bandwidth %||% 1.0
  kenv_kernel_kwargs <- params$kenv_kernel_kwargs %||% NULL
  heter_groups <- params$heter_groups %||% NULL
  num_cores <- params$num_cores %||% NULL
  nfolds <- params$nfolds %||% 5L
  parallel_mode <- params$parallel_mode %||% "auto"
  parallel_backend_prefer_fork <- params$parallel_backend_prefer_fork %||% TRUE
  sequential_models <- params$sequential_models %||% NULL
  docker_nd_usage <- params$docker_nd_usage %||% FALSE
  params$pheno_data <- pheno_data
  gp_output_level <- "predict_only"
  gp_return_se <- FALSE
  gp_full_vc <- FALSE
  gp_prediction_output <- "test_only"
  params$gp_return_trait_correlations <- FALSE
  params$return_trait_correlations <- FALSE
  gp_validate_eval_metrics(eval_metrics = eval_metrics, response_family = response_family)

  # Validate cross-validation method
  msg <- ""
  if (!is.null(cross_validation_meth) && length(cross_validation_meth) > 1) {
    log_error("Only one cross-validation method allowed at a time")
    stop(msg, "Use only one cross_validation method at a time.", call. = FALSE)
  }

  sampling_method <- gp_resolve_cv_sampling_method(cross_validation_meth, sampling_method)

  # Process omics data
  if (!is.null(ml_dat_res)) {
    omics_data <- ml_dat_res[["merged_data"]][["merge_data"]]
    if (!"merged_data_test" %in% names(ml_dat_res) && !is.null(test_set)) {
      omics_data <- omics_data[!rownames(omics_data) %in% test_set, ]
    }
    omic_count <- if ("omic_count" %in% names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL
  }

  # Model name mappings
  AI_valid_models <- gp_classical_ml_supported_models()
  bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
  gp_valid_models <- gp_lowrank_supported_models()
  canonical_names <- gp_deep_learning_supported_models()
  friendly_names <- gp_deep_learning_friendly_models()

  GS_model_cv_can <- gp_canonicalize_supported_model_names(GS_model_cv)
  dp_models <- intersect(GS_model_cv_can, canonical_names)
  if (length(dp_models) > 0L) {
    AI_valid_models <- c(AI_valid_models, dp_models)
  }
  GS_model_cv_labels <- GS_model_cv_display
  if (is.null(GS_model_cv_labels) || length(GS_model_cv_labels) != length(GS_model_cv_can)) {
    GS_model_cv_labels <- gp_display_supported_model_names(GS_model_cv)
  }
  if (!is.null(sequential_models)) {
    sequential_models <- gp_canonicalize_supported_model_names(sequential_models)
  }
  AI_valid_models <- unique(c(AI_valid_models, canonical_names))

  # Validate CV method
  cv_token <- normalize_cv_token(cross_validation_meth)
  holds_tokens <- c("hold_out", "stratified_hold_out", "repeated_hold_out", "repeated_stratified_hold_out")
  kfold_tokens <- c("k_folds", "stratified_k_folds", "repeated_k_folds", "repeated_stratified_k_folds")
  multi_tokens <- c("cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2")
  loo_tokens <- c("leave_one_out")

  is_holdout <- cv_token %in% holds_tokens
  is_kfold <- cv_token %in% kfold_tokens
  is_multi <- cv_token %in% multi_tokens
  is_loo <- cv_token %in% loo_tokens
  if (!(is_holdout || is_kfold || is_multi || is_loo)) {
    log_error("Unsupported cross-validation method: {cross_validation_meth}")
    stop(paste0("Unsupported cross-validation method: ", cross_validation_meth,
                "\nAllowed: Hold-Out, K-Folds, Leave-One-Out, or CV0/CV1/CV2 variants."), call. = FALSE)
  }

  sys_name <- Sys.info()[["sysname"]]
  if (docker_nd_usage) sys_name <- "Windows"

  if (any(GS_model_cv_can %in% gp_valid_models) &&
      !is.null(gen_name) &&
      !anyDuplicated(as.character(pheno_data[[gen_name]]))) {
    params$gp_backend_cache <- tryCatch(
      gp_backend_prepare_kernel_cache(
        ids = as.character(pheno_data[[gen_name]]),
        gmatrix = gmatrix_model_ready,
        omic1_kernel = omic1_kernel_model_ready,
        omic2_kernel = omic2_kernel_model_ready,
        omic3_kernel = omic3_kernel_model_ready,
        kernel_list = kernel_list,
        kernel_weights = lowrank_kernel_weights
      ),
      error = function(e) NULL
    )
  }
  params$gp_backend <- gp_backend
  params$gp_output_level <- gp_output_level
  params$gp_return_se <- gp_return_se
  params$gp_full_vc <- gp_full_vc
  params$gp_varcomp_mode <- gp_varcomp_mode
  params$gp_fa_rank <- gp_fa_rank
  params$gp_prediction_output <- gp_prediction_output
  params$gp_factor_cache <- gp_factor_cache
  params$gp_iters <- gp_iters
  params$gp_lr <- gp_lr
  params$gp_exact_fast_cv <- gp_exact_fast_cv
  params$kernel_list <- kernel_list
  params$env_similarity <- env_similarity
  params$env_ids <- env_ids
  params$env_covariates <- env_covariates
  params$reaction_norm_feature_qc <- reaction_norm_feature_qc
  params$kenv_kernel <- kenv_kernel
  params$kenv_bandwidth <- kenv_bandwidth
  params$kenv_kernel_kwargs <- kenv_kernel_kwargs

  # Create task grid. Canonical model IDs and display labels must stay paired;
  # crossing them creates duplicated model tasks and repeats expensive CV fits.
  model_idx <- seq_along(GS_model_cv_can)
  tasks_full <- expand.grid(response = response,
                            replication = seq_len(replication),
                            model_idx = model_idx,
                            feature_selection_mode = if (!is.null(feature_score_metadata) && !is.null(feature_k_grid)) {
                              gp_feature_cv_modes(feature_scoring_cv)
                            } else {
                              NA_character_
                            },
                            feature_k = if (!is.null(feature_score_metadata) && !is.null(feature_k_grid)) {
                              as.integer(feature_k_grid)
                            } else {
                              NA_integer_
                            },
                            stringsAsFactors = FALSE)
  tasks_full$modell <- as.character(GS_model_cv_can[tasks_full$model_idx])
  tasks_full$modell_label <- as.character(GS_model_cv_labels[tasks_full$model_idx])
  tasks_full$model_idx <- NULL
  tasks_full <- unique(tasks_full)
  tasks_full$fold_only <- NA_integer_
  tasks_full$fold_n_split <- NA_integer_
  tasks_unsplit <- tasks_full

  # Split each trait x replication x model x k task into one task per fold so
  # parallel workers are not bounded by the slowest model's sequential folds.
  # Not split: hold-out and LOO, GP models (already batched across folds in
  # one bridge call), MET ML/DL batch models, a single core, or when the user
  # asked for sequential execution or disabled it (cv_fold_tasks = FALSE).
  met_batch_models <- c("CatBoost", "LightGBM", "Xgboost", "RandomForest",
                        "mlp", "ft_transformer", "saint", "tabnet", "moe")
  fold_split_enabled <- (is_kfold || is_multi) && !is_loo &&
    isTRUE(as.integer(nfolds) > 1L) &&
    !identical(parallel_mode, "sequential") &&
    !isFALSE(params$cv_fold_tasks) &&
    !tolower(Sys.getenv("GP_CV_FOLD_TASKS", "true")) %in% c("false", "0", "no") &&
    (if (is.null(num_cores)) gp_detect_cores(logical = TRUE, reserve = 1L) else as.integer(num_cores)) > 1L
  if (isTRUE(fold_split_enabled)) {
    splittable <- !(tasks_full$modell %in% gp_valid_models) &
      !(isTRUE(met_ml_dl) & tasks_full$modell %in% met_batch_models)
    if (any(splittable)) {
      n_split <- as.integer(nfolds)
      expanded <- tasks_full[rep(which(splittable), each = n_split), , drop = FALSE]
      expanded$fold_only <- rep(seq_len(n_split), times = sum(splittable))
      expanded$fold_n_split <- n_split
      tasks_full <- rbind(tasks_full[!splittable, , drop = FALSE], expanded)
      rownames(tasks_full) <- NULL
    }
  }
  # Inputs for the parallel policy's cost model (see sp_decide_policy).
  cv_task_cost <- local({
    n_rows <- NROW(pheno_data)
    k <- if (is_kfold || is_multi) max(as.integer(nfolds), 2L) else 5L
    bayes_nIter <- tryCatch({
      vals <- vapply(model_prep_all_bayes_cv %||% list(), function(x) {
        as.numeric(x[["bayes_para"]][["nIter"]] %||% NA_real_)
      }, numeric(1))
      vals <- vals[is.finite(vals)]
      if (length(vals)) max(vals) else NA_real_
    }, error = function(e) NA_real_)
    n_markers <- tryCatch({
      if (exists("omics_data") && !is.null(omics_data)) NCOL(omics_data) else
        if (exists("gmatrix_model_ready") && !is.null(gmatrix_model_ready)) NCOL(gmatrix_model_ready) else 1000
    }, error = function(e) 1000)
    list(
      n_train = n_rows * (1 - 1 / k),
      n_markers = n_markers,
      nIter = if (is.finite(bayes_nIter)) bayes_nIter else (params$nIter %||% 6000),
      ntree = params$ntree %||% 500,
      iterations = params$iteration %||% params$catboost_iterations %||% 300,
      folds_per_task = if (any(!is.na(tasks_full$fold_only))) 1 else if (is_kfold || is_multi) k else 1
    )
  })
  task_group_key <- function(tbl) {
    paste(tbl$response, tbl$replication, tbl$modell, tbl$modell_label,
          tbl$feature_selection_mode, tbl$feature_k, sep = "\r")
  }
  gp_cv_direct_bridge_only <- all(tasks_full$modell %in% gp_valid_models) &&
    exists("gp_bridge_direct_execution_enabled", mode = "function") &&
    isTRUE(gp_bridge_direct_execution_enabled())
  # Models that never need an embedded (reticulate) Python: CLI-backed ML/GP
  # and the pure-R BGLR/ASReml fits. Initializing reticulate for a Bayes task
  # fails on machines whose python3 cannot be embedded (no shared libpython).
  direct_python_cv_models <- unique(c(AI_valid_models, gp_valid_models, gp_parallel_r_only_models()))
  cv_scheduler_needs_python_init <- any(!tasks_full$modell %in% direct_python_cv_models)

  policy_model_params <- lapply(tasks_full$modell, function(mod) {
    gp_parallel_extract_model_policy_params(mod, params)
  })

  # Single task execution
  # Fold-level scheduling: a task row with `fold_only` fits only the folds at
  # positions p with ((p - 1) %% fold_n_split) + 1 == fold_only and returns
  # its out-of-fold predictions (a "partial"). The partials of one
  # trait x replication x model x k are merged and passed back as `merged`,
  # which skips all fitting and computes the metrics exactly as an unsplit
  # task would.
  run_one_task <- function(task_row, merged = NULL) {
    log_debug("Running task for {task_row$response}, rep {task_row$replication}, model {task_row$modell}")
    # Fold splits and fold seeds come from set.seed(); parallel workers run
    # L'Ecuyer-CMRG streams, under which the same seed gives different folds.
    # Pin the generator inside the task so folds (and results) are identical
    # on every backend.
    old_rng_kind_task <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    on.exit(do.call(RNGkind, as.list(old_rng_kind_task)), add = TRUE)
    fold_only <- suppressWarnings(as.integer(task_row[["fold_only"]] %||% NA_integer_))
    if (!length(fold_only)) fold_only <- NA_integer_
    fold_n_split <- suppressWarnings(as.integer(task_row[["fold_n_split"]] %||% NA_integer_))
    fit_fold <- function(pos) {
      is.null(merged) && (is.na(fold_only) || ((pos - 1L) %% fold_n_split) + 1L == fold_only)
    }
    if (is.null(merged) &&
        !as.character(task_row$modell) %in% direct_python_cv_models && is.function(init_py)) init_py()
    trait <- as.character(task_row$response)
    rep_i <- as.integer(task_row$replication)
    model_can <- as.character(task_row$modell)
    model_label <- as.character(task_row$modell_label)
    model_prep_all_bayes_cv_task <- model_prep_all_bayes_cv
    asreml_models_prep_cv_task <- asreml_models_prep_cv
    gmatrix_model_ready_task <- gmatrix_model_ready
    omic1_kernel_model_ready_task <- omic1_kernel_model_ready
    omic2_kernel_model_ready_task <- omic2_kernel_model_ready
    omic3_kernel_model_ready_task <- omic3_kernel_model_ready
    kernel_list_task <- kernel_list
    gmatrix_model_ready <- gmatrix_model_ready_task
    omic1_kernel_model_ready <- omic1_kernel_model_ready_task
    omic2_kernel_model_ready <- omic2_kernel_model_ready_task
    omic3_kernel_model_ready <- omic3_kernel_model_ready_task
    kernel_list <- kernel_list_task
    feature_k <- suppressWarnings(as.integer(task_row$feature_k %||% NA_integer_))
    feature_selection_mode <- as.character(task_row$feature_selection_mode %||% NA_character_)
    if (is.na(feature_selection_mode) || !nzchar(feature_selection_mode)) {
      feature_selection_mode <- "fixed"
    }
    omics_data <- if (exists("omics_data")) omics_data else NULL
    if (!is.null(omics_data) && !is.null(selected_raw)) {
      omics_data <- gp_feature_select_explicit(
        x = omics_data,
        feature_selected = selected_raw,
        trait = trait,
        context = "cross-validation predictor data"
      )
    }
    if (!is.null(omics_data) && !is.null(feature_score_metadata) && is.finite(feature_k) &&
        identical(feature_selection_mode, "fixed")) {
      omics_data <- feature_select_top_k(
        predictor_data = omics_data,
        feature_score_metadata = feature_score_metadata,
        trait = trait,
        k = feature_k
      )
    }
    omics_data_base <- omics_data
    feature_selection_records <- list()
    feature_model_adapter_required <- model_can %in% c(
      bayes_valid_models,
      bayes_gblup_valid_models,
      "GBLUP",
      gp_valid_models
    ) || isTRUE(met_ml_dl)
    apply_selected_model_inputs <- function(selected_by_source) {
      if (!isTRUE(feature_model_adapter_required)) return(invisible(NULL))
      artifacts <- gp_feature_prepare_cv_model_inputs(
        model = model_can,
        trait = trait,
        pheno_data = pheno_data,
        gen_name = gen_name,
        selected_by_source = selected_by_source,
        source_matrices = feature_source_matrices,
        params = params
      )
      if (!is.null(artifacts[["model_prep_all_bayes_cv"]])) {
        model_prep_all_bayes_cv_task <<- artifacts[["model_prep_all_bayes_cv"]]
      }
      if (!is.null(artifacts[["asreml_models_prep_cv"]])) {
        asreml_models_prep_cv_task <<- artifacts[["asreml_models_prep_cv"]]
      }
      bank <- artifacts[["kernel_bank"]] %||% list()
      if (length(bank)) {
        gmatrix_model_ready_task <<- bank[["gmatrix_model_ready"]] %||% NULL
        omic1_kernel_model_ready_task <<- bank[["omic1_kernel_model_ready"]] %||% NULL
        omic2_kernel_model_ready_task <<- bank[["omic2_kernel_model_ready"]] %||% NULL
        omic3_kernel_model_ready_task <<- bank[["omic3_kernel_model_ready"]] %||% NULL
        kernel_list_task <<- bank
        gmatrix_model_ready <<- gmatrix_model_ready_task
        omic1_kernel_model_ready <<- omic1_kernel_model_ready_task
        omic2_kernel_model_ready <<- omic2_kernel_model_ready_task
        omic3_kernel_model_ready <<- omic3_kernel_model_ready_task
        kernel_list <<- kernel_list_task
        params$gmatrix_model_ready <<- gmatrix_model_ready
        params$omic1_kernel_model_ready <<- omic1_kernel_model_ready
        params$omic2_kernel_model_ready <<- omic2_kernel_model_ready
        params$omic3_kernel_model_ready <<- omic3_kernel_model_ready
        params$kernel_list <<- kernel_list
        params$gp_backend_cache <<- NULL
      }
      params$asreml_models_prep_cv <<- asreml_models_prep_cv_task
      invisible(artifacts)
    }
    if (!is.null(selected_raw)) {
      explicit_selected_by_source <- gp_feature_explicit_selected_by_source(
        feature_selected = selected_raw,
        trait = trait,
        source_matrices = feature_source_matrices,
        context = "feature_selected"
      )
      precomputed_kernel_supplied <- isTRUE(params$feature_precomputed_kernel_supplied)
      if (isTRUE(feature_model_adapter_required) && length(explicit_selected_by_source) &&
          !(gp_feature_model_needs_kernel(model_can, met_ml_dl = met_ml_dl) &&
            precomputed_kernel_supplied)) {
        apply_selected_model_inputs(explicit_selected_by_source)
      }
      metadata_sources <- explicit_selected_by_source
      if (!length(metadata_sources) && !is.null(omics_data)) {
        metadata_sources <- list(merged = colnames(omics_data))
      }
      if (!length(metadata_sources)) {
        metadata_sources <- gp_feature_declared_by_source(selected_raw, trait = trait)
      }
      explicit_meta <- gp_feature_explicit_metadata(
        selected_by_source = metadata_sources,
        trait = trait,
        replication = rep_i,
        prediction_model = model_label,
        n_train = sum(!is.na(pheno_data[[trait]]))
      )
      if (!is.null(explicit_meta)) {
        feature_selection_records[[length(feature_selection_records) + 1L]] <- explicit_meta
      }
      if (is.null(feature_score_metadata)) {
        feature_k <- as.integer(sum(lengths(metadata_sources)))
        feature_selection_mode <- "explicit"
        feature_scoring_model <- "user_supplied"
      }
    }
    if (is.finite(feature_k) && identical(feature_selection_mode, "fixed")) {
      fixed_meta <- feature_score_metadata[
        as.character(feature_score_metadata[["trait"]]) == trait,
        ,
        drop = FALSE
      ]
      if (nrow(fixed_meta)) {
        fixed_selected <- gp_feature_selected_metadata(fixed_meta, trait = trait, k = feature_k)
        fixed_meta[["selected"]] <- fixed_meta[["predictor"]] %in% fixed_selected[["predictor"]]
        fixed_meta[["feature_k"]] <- as.integer(sum(fixed_meta[["selected"]]))
        fixed_meta[["replication"]] <- rep_i
        fixed_meta[["fold"]] <- NA_integer_
        fixed_meta[["selection_mode"]] <- "fixed"
        fixed_meta[["prediction_model"]] <- model_label
        feature_selection_records[[1L]] <- fixed_meta
      }
    }
    apply_fold_internal_features <- function(tst, fold_id) {
      if (!identical(feature_selection_mode, "fold_internal") || !is.finite(feature_k)) {
        omics_data <<- omics_data_base
        return(invisible(NULL))
      }
      # k covering every predictor is the "all predictors" candidate: fold
      # scoring cannot change the set, and the model-ready kernels were
      # already built from all predictors, so skip rescoring and rebuilding.
      if (!is.null(omics_data_base) && feature_k >= NCOL(omics_data_base)) {
        omics_data <<- omics_data_base
        return(invisible(NULL))
      }
      fold_number <- suppressWarnings(as.integer(fold_id))
      if (!length(fold_number) || is.na(fold_number)) {
        fold_number <- 1L
      }
      view <- gp_feature_fold_view(
        predictor_data = omics_data_base,
        pheno_data = pheno_data,
        response = trait,
        gen_name = gen_name,
        test_rows = tst,
        k = feature_k,
        scoring_model = feature_scoring_model,
        seed = feature_scoring_seed + rep_i * 10000L + fold_number,
        source_block = feature_source_map,
        response_family = response_family,
        replication = rep_i,
        fold = fold_number,
        ridge_lambda = feature_ridge_lambda,
        bayes_nIter = feature_bayes_nIter,
        bayes_burnIn = feature_bayes_burnIn,
        bayes_thin = feature_bayes_thin,
        ntree = params$ntree %||% 500L,
        mtry = params$mtry,
        nodesize = params$nodesize,
        rf_n_jobs = params$rf_n_jobs %||% 1L
      )
      omics_data <<- view$predictor_data
      apply_selected_model_inputs(view$selected_by_source)
      meta <- view$metadata
      meta[["prediction_model"]] <- model_label
      feature_selection_records[[length(feature_selection_records) + 1L]] <<- meta
      invisible(view)
    }

    repp <- rep_i
    base_seed <- if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1)
      as.integer(random_state) else 123L
    new_seed <- (base_seed + rep_i * 10000L) %% .Machine$integer.max
    # Each fold gets its own seed so its fit does not depend on which other
    # folds ran before it in the same process (identical results whether
    # folds run in one task or as separate parallel tasks).
    seed_fold <- function(pos) {
      set.seed((new_seed + as.integer(pos) * 7919L) %% .Machine$integer.max)
    }
    params$gid_vector <- as.character(pheno_data[[gen_name]])
    params$gmatrix_model_ready <- gmatrix_model_ready
    params$omic1_kernel_model_ready <- omic1_kernel_model_ready
    params$omic2_kernel_model_ready <- omic2_kernel_model_ready
    params$omic3_kernel_model_ready <- omic3_kernel_model_ready
    params$kernel_list <- kernel_list
    params$lowrank_eps_trace <- lowrank_eps_trace
    params$lowrank_max_rank <- lowrank_max_rank
    params$lowrank_jitter <- lowrank_jitter
    params$lowrank_noise_grid <- lowrank_noise_grid
    params$lowrank_kernel_weights <- lowrank_kernel_weights
    params$response <- trait
    if (is.null(merged) && is.finite(feature_k) && identical(feature_selection_mode, "fixed")) {
      apply_selected_model_inputs(
        gp_feature_selected_by_source(feature_score_metadata, trait = trait, k = feature_k)
      )
    }

    trait_values <- pheno_data[[trait]]
    valid_idx <- if (identical(response_family, "gaussian")) {
      which(is.finite(as.double(trait_values)))
    } else {
      which(!is.na(trait_values))
    }
    if (length(valid_idx) == 0L) {
      log_warn("Trait {trait} has no finite observations; skipping")
      return(NULL)
    }
    if (length(valid_idx) < nrow(pheno_data) && isTRUE(verbose)) {
      log_warn(
        "Trait {trait}: cross-validation will use {length(valid_idx)} of {nrow(pheno_data)} rows with non-missing observations"
      )
    }
    pheno_data_split <- pheno_data[valid_idx, , drop = FALSE]

    # Generate CV splits
    effective_nfolds <- nfolds
    if (is_holdout) {
      test_set_val <- hold_out_stratified_and_un(
        pheno_data = pheno_data_split, gen_name = gen_name, response = trait,
        test_size = test_size, random_state = new_seed, replication = repp,
        sampling_method = sampling_method
      )
    } else if (is_loo) {
      effective_nfolds <- nrow(pheno_data_split)
      test_set_val <- stats::setNames(
        replicate(repp, seq_len(effective_nfolds), simplify = FALSE),
        sprintf("Rep%d_%dFold_LOO", seq_len(repp), effective_nfolds)
      )
    } else if (is_kfold) {
      test_set_val <- kfolds_stratified_un(
        pheno_data = pheno_data_split, gen_name = gen_name, response = trait,
        nfolds = nfolds, random_state = new_seed, replication = repp,
        sampling_method = sampling_method
      )
    } else {
      CVi <- as.integer(sub("^.*?([0-2])$", "\\1", tolower(cross_validation_meth)))
      test_set_val <- CV0_CV1_CV2_for_multi_environment(
        pheno_data = pheno_data_split, gen_name = gen_name,
        heter_groups = heter_groups,
        CV = CVi, nfolds = nfolds,
        random_state = new_seed, replication = repp, message = verbose
      )
      if (isTRUE(verbose)) {
        n_ids <- length(unique(as.character(pheno_data_split[[gen_name]])))
        log_info("Trait {trait}: multi-environment CV is being carried out on {n_ids} unique individual(s)")
      }
    }

    len_y <- nrow(pheno_data)
    y <- if (identical(response_family, "gaussian")) as.double(pheno_data[[trait]]) else pheno_data[[trait]]
    class_levels <- class_levels_for_y(pheno_data[[trait]], response_family)
    handle_error <- FALSE
    # Phase 3.17: warn-once dictionary for Bayes CV tryCatch handlers --
    # silent log_error calls used to swallow real engine errors for
    # RKHS / GBLUP_BRR / BRR; now they emit one warning() per (trait,
    # model) so users see the cause.
    bayes_warned <- list()
    yhat_init <- if (response_family %in% c("gaussian", "binary")) {
      rep(NA_real_, len_y)
    } else {
      rep(NA_character_, len_y)
    }
    ypred_cv <- data.frame(
      row_id = seq_len(len_y),
      GID = if (!is.null(gen_name) && gen_name %in% names(pheno_data)) {
        as.character(pheno_data[[gen_name]])
      } else {
        as.character(seq_len(len_y))
      },
      y = y,
      yhat = yhat_init,
      cv_role = "train",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    need_prob_storage <- !identical(response_family, "gaussian") ||
      any(vapply(eval_metrics, metric_needs_probabilities, logical(1), fam = response_family))
    yprob_cv <- NULL
    results_eval_metrics_reps <- NULL
    gp_batch_cv_info <- NULL

    if (is_multi) {
      if (is.null(heter_groups)) {
        log_error("heter_groups required for CV1/CV2")
        stop("For CV1/CV2 you must supply `heter_groups`.", call. = FALSE)
      }
      ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
    }

    assign_probabilities <- function(tst, pred_prob) {
      if (is.null(pred_prob)) return()
      pred_prob <- as.matrix(pred_prob)
      if (!nrow(pred_prob)) return()
      prob_levels <- colnames(pred_prob)
      if (is.null(prob_levels)) {
        prob_levels <- if (ncol(pred_prob) == 1L && identical(response_family, "binary")) {
          tail(class_levels %||% "1", 1L)
        } else {
          class_levels %||% as.character(seq_len(ncol(pred_prob)))
        }
      }
      colnames(pred_prob) <- prob_levels
      if (is.null(yprob_cv)) {
        yprob_cv <<- matrix(NA_real_, nrow = len_y, ncol = ncol(pred_prob))
        colnames(yprob_cv) <<- prob_levels
      }
      missing_cols <- setdiff(prob_levels, colnames(yprob_cv))
      if (length(missing_cols)) {
        expanded <- matrix(NA_real_, nrow = nrow(yprob_cv), ncol = ncol(yprob_cv) + length(missing_cols))
        colnames(expanded) <- c(colnames(yprob_cv), missing_cols)
        expanded[, colnames(yprob_cv)] <- yprob_cv
        yprob_cv <<- expanded
      }
      yprob_cv[tst, prob_levels] <<- pred_prob
    }

    metric_pred_input <- function(idx, metric) {
      # A classification backend can return both a convenience class label and
      # a probability matrix. Use the probability payload as the single source
      # for every classification metric so accuracy/F1 and reliability/ECE
      # cannot disagree about the predicted class.
      if (!identical(response_family, "gaussian") && !is.null(yprob_cv)) {
        prob <- yprob_cv[idx, , drop = FALSE]
        if (identical(response_family, "binary") && ncol(prob) == 1L) {
          return(as.numeric(prob[, 1L, drop = TRUE]))
        }
        return(prob)
      }
      if (!metric_needs_probabilities(metric, response_family)) {
        return(ypred_cv$yhat[idx])
      }
      if (!is.null(yprob_cv)) {
        prob <- yprob_cv[idx, , drop = FALSE]
        if (identical(response_family, "binary") && ncol(prob) == 1L) {
          return(as.numeric(prob[, 1, drop = TRUE]))
        }
        return(prob)
      }
      ypred_cv$yhat[idx]
    }

    assign_preds <- function(tst, preds) {
      if (is.null(preds)) return()
      pred_prob <- attr(preds, "probabilities", exact = TRUE)
      if (is.matrix(preds) || is.data.frame(preds)) {
        pred_prob <- as.matrix(preds)
        preds <- if (identical(response_family, "binary")) {
          if (ncol(pred_prob) == 1L) {
            as.numeric(pred_prob[, 1L, drop = TRUE])
          } else {
            positive_col <- match(tail(class_levels, 1L), colnames(pred_prob))
            if (is.na(positive_col)) positive_col <- ncol(pred_prob)
            as.numeric(pred_prob[, positive_col, drop = TRUE])
          }
        } else {
          gp_pred_to_class_labels(
            pred_prob,
            positive_class = positive_class,
            class_levels = class_levels
          )
        }
      } else if (need_prob_storage && !is.null(pred_prob)) {
        assign_probabilities(tst, pred_prob)
      }
      if (!is.null(pred_prob) && (is.matrix(pred_prob) || is.data.frame(pred_prob))) {
        assign_probabilities(tst, pred_prob)
      }
      cast_preds <- if (identical(response_family, "gaussian") ||
                        (identical(response_family, "binary") && is.numeric(preds))) {
        as.double
      } else {
        as.character
      }
      if (length(preds) == length(tst)) {
        ypred_cv$yhat[tst] <<- cast_preds(preds)
        return()
      }
      idx <- NULL
      if (!is.null(names(preds))) {
        nm <- names(preds)
        if (!is.null(rownames(pheno_data))) {
          idx <- match(nm, rownames(pheno_data))
        }
      }
      if (!is.null(idx) && any(!is.na(idx))) {
        ok_idx <- which(!is.na(idx))
        ypred_cv$yhat[idx[ok_idx]] <<- cast_preds(preds[ok_idx])
        return()
      }
      m <- min(length(preds), length(tst))
      if (m > 0) {
        ypred_cv$yhat[tst[seq_len(m)]] <<- cast_preds(preds[seq_len(m)])
      }
      if (isTRUE(verbose)) log_warn("Preds length {length(preds)} != tst length {length(tst)}; assigned {m}")
    }

    assign_gp_batch_preds <- function(fold_list, context = "CV") {
      if (!is.null(merged) ||
          !(model_can %in% gp_valid_models) ||
          identical(feature_selection_mode, "fold_internal") ||
          !length(fold_list) ||
          !exists("gp_backend_cv_predict_batch", mode = "function") ||
          !isTRUE(gp_backend_cv_batch_enabled())) {
        return(FALSE)
      }
      fold_list <- lapply(fold_list, function(idx) as.integer(idx[!is.na(idx)]))
      fold_list <- fold_list[lengths(fold_list) > 0L]
      if (!length(fold_list)) {
        return(FALSE)
      }
      tryCatch({
        batch_preds <- gp_backend_cv_predict_batch(
          model_name = model_can,
          y = y,
          folds = fold_list,
          additional_params = params
        )
        if (!is.list(batch_preds) || length(batch_preds) != length(fold_list)) {
          stop("Batched GP CV returned an unexpected fold count.", call. = FALSE)
        }
        info <- attr(batch_preds, "gp_info", exact = TRUE)
        if (is.list(info)) {
          keep <- c(
            "package_cv_batch", "bridge_command", "cv_fast_path",
            "cv_fast_path_skip_reason", "cv_fast_path_approximate",
            "cv_fast_path_exact_refit", "cv_fast_path_solver",
            "cv_fast_path_shared_hyperparameters",
            "cv_fast_path_likelihood_noise_variance",
            "cv_fast_path_likelihood_noise_std", "cv_fast_path_self_checked",
            "cv_fast_path_self_check_reason", "cv_preprocessed_context",
            "cv_preprocessed_context_rows",
            "cv_preprocessed_context_tensor_cache_entries", "gp_exact_fast_cv"
          )
          gp_batch_cv_info <<- info[intersect(keep, names(info))]
        }
        for (ii in seq_along(fold_list)) {
          tst_i <- fold_list[[ii]]
          ypred_cv$cv_role[tst_i] <<- "test"
          assign_preds(tst_i, batch_preds[[ii]])
        }
        TRUE
      }, error = function(e) {
        warning(
          sprintf(
            "Batched GP %s failed for %s; falling back to per-fold bridge: %s",
            context,
            trait,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        if (isTRUE(verbose)) {
          log_warn("Batched GP {context} failed for {trait}; falling back to per-fold bridge: {conditionMessage(e)}")
        }
        FALSE
      })
    }

    assign_met_ai_batch_preds <- function(fold_list, context = "CV") {
      met_models <- c("CatBoost", "LightGBM", "Xgboost", "RandomForest",
                      "mlp", "ft_transformer", "saint", "tabnet", "moe")
      met_dl_models <- c("mlp", "ft_transformer", "saint", "tabnet", "moe")
      if (!is.null(merged) ||
          !isTRUE(met_ml_dl) ||
          !(model_can %in% met_models) ||
          identical(feature_selection_mode, "fold_internal") ||
          !length(fold_list) ||
          is.null(heter_groups) ||
          !nzchar(as.character(heter_groups)[1L])) {
        return(FALSE)
      }
      fold_list <- lapply(fold_list, function(idx) as.integer(idx[!is.na(idx)]))
      fold_list <- fold_list[lengths(fold_list) > 0L]
      if (!length(fold_list)) {
        return(FALSE)
      }
      tryCatch({
        batch_preds <- if (model_can %in% met_dl_models) {
          gp_met_dl_cv_predict_batch(
            model_type = model_can,
            pheno_data = pheno_data,
            response = trait,
            gen_name = gen_name,
            heter_groups = heter_groups,
            response_family = response_family,
            gmatrix = gmatrix_model_ready,
            omic1_kernel = omic1_kernel_model_ready,
            omic2_kernel = omic2_kernel_model_ready,
            omic3_kernel = omic3_kernel_model_ready,
            kernel_list = kernel_list,
            omics_kernel_label = omics_kernel_label,
            folds = fold_list,
            model_params = gp_met_dl_model_params(model_can, params),
            met_kernel_var_explained = met_kernel_var_explained,
            met_kernel_min_ev = met_kernel_min_ev,
            met_kernel_max_pcs = met_kernel_max_pcs
          )
        } else {
          gp_met_python_cv_predict_batch(
            pheno_data = pheno_data,
            response = trait,
            gen_name = gen_name,
            heter_groups = heter_groups,
            response_family = response_family,
            model_type = model_can,
            gmatrix = gmatrix_model_ready,
            omic1_kernel = omic1_kernel_model_ready,
            omic2_kernel = omic2_kernel_model_ready,
            omic3_kernel = omic3_kernel_model_ready,
            kernel_list = kernel_list,
            omics_kernel_label = omics_kernel_label,
            folds = fold_list,
            model_params = gp_met_python_model_params(model_can, params),
            met_kernel_var_explained = met_kernel_var_explained,
            met_kernel_min_ev = met_kernel_min_ev,
            met_kernel_max_pcs = met_kernel_max_pcs
          )
        }
        if (!is.list(batch_preds) || length(batch_preds) != length(fold_list)) {
          stop("Batched MET ML/DL CV returned an unexpected fold count.", call. = FALSE)
        }
        for (ii in seq_along(fold_list)) {
          tst_i <- fold_list[[ii]]
          ypred_cv$cv_role[tst_i] <<- "test"
          assign_preds(tst_i, batch_preds[[ii]])
        }
        TRUE
      }, error = function(e) {
        warning(
          sprintf(
            "Batched MET ML/DL %s failed for %s; falling back to per-fold bridge: %s",
            context,
            trait,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        if (isTRUE(verbose)) {
          log_warn("Batched MET ML/DL {context} failed for {trait}; falling back to per-fold bridge: {conditionMessage(e)}")
        }
        FALSE
      })
    }

    # After the fold loop: a fold task returns its partial out-of-fold
    # predictions; a merge call adopts the merged predictions of all fold
    # tasks so the metric code below runs unchanged.
    finish_fold_phase <- function() {
      if (!is.null(merged)) {
        ypred_cv <<- merged$ypred_cv
        yprob_cv <<- merged$yprob_cv
        feature_selection_records <<- c(feature_selection_records, merged$feature_selection_records)
        gp_batch_cv_info <<- merged$gp_batch_cv_info %||% gp_batch_cv_info
        return(NULL)
      }
      if (is.na(fold_only)) {
        return(NULL)
      }
      fold_records <- lapply(feature_selection_records, function(r) {
        if (is.data.frame(r) && "fold" %in% names(r)) r[!is.na(r$fold), , drop = FALSE] else NULL
      })
      list(
        fold_partial = TRUE,
        ypred_cv = ypred_cv,
        yprob_cv = yprob_cv,
        feature_selection_records = Filter(Negate(is.null), fold_records),
        gp_batch_cv_info = gp_batch_cv_info
      )
    }

    if (is_holdout) {
      tst <- valid_idx[test_set_val[[repp]]]
      seed_fold(1L)
      apply_fold_internal_features(tst, 1L)
      yNA <- y; yNA[tst] <- NA
      ypred_cv$cv_role[tst] <- "test"
      gp_batched_done <- assign_gp_batch_preds(list(tst), context = "Hold_Out")
      met_ai_batched_done <- if (!isTRUE(gp_batched_done)) {
        assign_met_ai_batch_preds(list(tst), context = "Hold_Out")
      } else {
        FALSE
      }

      if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            params$bayes_model <- if (model_can == "GBLUP_BRR") "BRR" else model_can
            params$bayes_trait <- trait
            params$ETA <- model_prep_all_bayes_cv_task[[model_can]][["bayes_ETA"]][["ETA"]]
            params$bayes_para <- model_prep_all_bayes_cv_task[[model_can]][["bayes_para"]]
            params$bayes_groups <- model_prep_all_bayes_cv_task[[model_can]][["bayes_groups"]]
            params$met_kernel_cv_meta <- model_prep_all_bayes_cv_task[[model_can]][["met_kernel_cv_meta"]]
            assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = params))
          }, error = function(e) {
            log_error("Bayes {model_can} failed for {trait}: {conditionMessage(e)}")
            # Phase 3.17: make the Bayes silent-drop visible to the user.
            warning(sprintf("Bayes %s failed for trait '%s': %s",
                            model_can, trait, conditionMessage(e)),
                    call. = FALSE)
            handle_error <<- TRUE
        })
      }

      if (!is.null(params$engine) && model_can == "GBLUP" && params$engine == "asreml") {
        tryCatch({
          params$response <- trait
          assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = params))
        }, error = function(e) {
          log_error("GBLUP(asreml) failed for {trait}: {conditionMessage(e)}")
          handle_error <<- TRUE
        })
      }

      if (model_can %in% gp_valid_models && !isTRUE(gp_batched_done)) {
        tryCatch({
          assign_preds(tst, predict_with_model(model_can, y = yNA, tst = tst, additional_params = params))
        }, error = function(e) {
          log_error("{model_can} failed for {trait}: {conditionMessage(e)}")
          handle_error <<- TRUE
        })
      }

        if (model_can %in% AI_valid_models && !isTRUE(met_ai_batched_done)) {
          tryCatch({
            if (isTRUE(met_ml_dl) && model_can %in% c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "mlp", "ft_transformer", "saint", "tabnet", "moe") && is_multi) {
              pred_met <- if (model_can %in% c("mlp", "ft_transformer", "saint", "tabnet", "moe")) {
                gp_met_dl_cv_predict(
                  model_type = model_can,
                  pheno_data = pheno_data,
                  response = trait,
                  gen_name = gen_name,
                  heter_groups = heter_groups,
                  response_family = response_family,
                  gmatrix = gmatrix_model_ready,
                  omic1_kernel = omic1_kernel_model_ready,
                  omic2_kernel = omic2_kernel_model_ready,
                  omic3_kernel = omic3_kernel_model_ready,
                  kernel_list = kernel_list,
                  omics_kernel_label = omics_kernel_label,
                  tst = tst,
                  model_params = gp_met_dl_model_params(model_can, params),
                  met_kernel_var_explained = met_kernel_var_explained,
                  met_kernel_min_ev = met_kernel_min_ev,
                  met_kernel_max_pcs = met_kernel_max_pcs
                )
              } else {
                gp_met_python_cv_predict(
                  pheno_data = pheno_data,
                  response = trait,
                  gen_name = gen_name,
                  heter_groups = heter_groups,
                  response_family = response_family,
                  model_type = model_can,
                  gmatrix = gmatrix_model_ready,
                  omic1_kernel = omic1_kernel_model_ready,
                  omic2_kernel = omic2_kernel_model_ready,
                  omic3_kernel = omic3_kernel_model_ready,
                  kernel_list = kernel_list,
                  omics_kernel_label = omics_kernel_label,
                  tst = tst,
                  model_params = gp_met_python_model_params(model_can, params),
                  met_kernel_var_explained = met_kernel_var_explained,
                  met_kernel_min_ev = met_kernel_min_ev,
                  met_kernel_max_pcs = met_kernel_max_pcs
                )
              }
              assign_preds(tst, pred_met)
            } else {
              assign_preds(tst, predict_with_model(
                model = model_can, y = yNA,
                omics_data = if (exists("omics_data")) omics_data else NULL,
                tst = tst, additional_params = params))
            }
          }, error = function(e) {
            log_error("AI model {model_can} failed for {trait}: {conditionMessage(e)}")
          handle_error <<- TRUE
        })
      }

      test_idx <- which(ypred_cv$cv_role == "test")
      y_t <- ypred_cv$y[test_idx]
      y_p <- ypred_cv$yhat[test_idx]
      ok <- gp_metric_valid_rows(y_t, y_p, response_family = response_family)
      n_ok <- sum(ok)
      results_eval_metrics_reps <- data.frame(Rep = repp, check.names = FALSE)
      if (is.finite(feature_k)) {
        results_eval_metrics_reps$feature_k <- feature_k
        results_eval_metrics_reps$feature_scoring_model <- as.character(feature_scoring_model)
        results_eval_metrics_reps$feature_scoring_cv <- feature_selection_mode
      }
      degenerate_pred <- response_family == "gaussian" &&
        isTRUE(stats::sd(y_p[ok]) == 0)
      if (n_ok < 2 || degenerate_pred) {
        for (m in eval_metrics) {
          results_eval_metrics_reps[[m]] <- if (identical(response_family, "gaussian")) {
            metric_worst_value(m, y_true = y_t[ok])
          } else {
            NA_real_
          }
        }
      } else {
        for (m in eval_metrics) {
          results_eval_metrics_reps[[m]] <- safe_metric_value(
            y_true = y_t,
            y_pred = metric_pred_input(test_idx, m),
            metric = m,
            response_family = response_family,
            positive_class = positive_class
          )
        }
      }
    }

    if (is_kfold || is_loo) {
      group <- rep(NA_integer_, nrow(pheno_data))
      group[valid_idx] <- test_set_val[[repp]]
      gp_fold_list <- lapply(seq_len(effective_nfolds), function(j) which(group == j))
      gp_batched_done <- assign_gp_batch_preds(gp_fold_list, context = if (is_loo) "LOO" else "K-Folds")
      met_ai_batched_done <- if (!isTRUE(gp_batched_done)) {
        assign_met_ai_batch_preds(gp_fold_list, context = if (is_loo) "LOO" else "K-Folds")
      } else {
        FALSE
      }
      if (!isTRUE(gp_batched_done) && !isTRUE(met_ai_batched_done)) {
        for (j in seq_len(effective_nfolds)) {
          tst <- which(group == j)
          if (!length(tst)) next
          if (!fit_fold(j)) next
          seed_fold(j)
          apply_fold_internal_features(tst, j)
          yNA <- y; yNA[tst] <- NA
          ypred_cv$cv_role[tst] <- "test"

          if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
            tryCatch({
              params$bayes_model <- if (model_can == "GBLUP_BRR") "BRR" else model_can
              params$bayes_trait <- trait
              params$ETA <- model_prep_all_bayes_cv_task[[model_can]][["bayes_ETA"]][["ETA"]]
              params$bayes_para <- model_prep_all_bayes_cv_task[[model_can]][["bayes_para"]]
              params$bayes_groups <- model_prep_all_bayes_cv_task[[model_can]][["bayes_groups"]]
              params$met_kernel_cv_meta <- model_prep_all_bayes_cv_task[[model_can]][["met_kernel_cv_meta"]]
              assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = params))
            }, error = function(e) {
              log_error("Bayes {model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
              # Phase 3.17: warn-once per (trait, model) so users see the
              # cause without spamming once per fold.
              if (!isTRUE(bayes_warned[[paste(trait, model_can, sep="::")]])) {
                warning(sprintf("Bayes %s failed for trait '%s' (fold %s; same error in subsequent folds suppressed): %s",
                                model_can, trait, j, conditionMessage(e)),
                        call. = FALSE)
                bayes_warned[[paste(trait, model_can, sep="::")]] <<- TRUE
              }
              handle_error <<- TRUE
            })
          }

          if (!is.null(params$engine) && model_can == "GBLUP" && params$engine == "asreml") {
            tryCatch({
              params$response <- trait
              assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = params))
            }, error = function(e) {
              log_error("GBLUP(asreml) failed for {trait} fold {j}: {conditionMessage(e)}")
              handle_error <<- TRUE
            })
          }

          if (model_can %in% gp_valid_models) {
            tryCatch({
              assign_preds(tst, predict_with_model(model_can, y = yNA, tst = tst, additional_params = params))
            }, error = function(e) {
              log_error("{model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
              handle_error <<- TRUE
            })
          }

          if (model_can %in% AI_valid_models) {
            tryCatch({
              if (isTRUE(met_ml_dl) && model_can %in% c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "mlp", "ft_transformer", "saint", "tabnet", "moe") && is_multi) {
                pred_met <- if (model_can %in% c("mlp", "ft_transformer", "saint", "tabnet", "moe")) {
                  gp_met_dl_cv_predict(
                    model_type = model_can,
                    pheno_data = pheno_data,
                    response = trait,
                    gen_name = gen_name,
                    heter_groups = heter_groups,
                    response_family = response_family,
                    gmatrix = gmatrix_model_ready,
                    omic1_kernel = omic1_kernel_model_ready,
                    omic2_kernel = omic2_kernel_model_ready,
                    omic3_kernel = omic3_kernel_model_ready,
                    kernel_list = kernel_list,
                    omics_kernel_label = omics_kernel_label,
                    tst = tst,
                    model_params = gp_met_dl_model_params(model_can, params),
                    met_kernel_var_explained = met_kernel_var_explained,
                    met_kernel_min_ev = met_kernel_min_ev,
                    met_kernel_max_pcs = met_kernel_max_pcs
                  )
                } else {
                  gp_met_python_cv_predict(
                    pheno_data = pheno_data,
                    response = trait,
                    gen_name = gen_name,
                    heter_groups = heter_groups,
                    response_family = response_family,
                    model_type = model_can,
                    gmatrix = gmatrix_model_ready,
                    omic1_kernel = omic1_kernel_model_ready,
                    omic2_kernel = omic2_kernel_model_ready,
                    omic3_kernel = omic3_kernel_model_ready,
                    kernel_list = kernel_list,
                    omics_kernel_label = omics_kernel_label,
                    tst = tst,
                    model_params = gp_met_python_model_params(model_can, params),
                    met_kernel_var_explained = met_kernel_var_explained,
                    met_kernel_min_ev = met_kernel_min_ev,
                    met_kernel_max_pcs = met_kernel_max_pcs
                  )
                }
                assign_preds(tst, pred_met)
              } else {
                assign_preds(tst, predict_with_model(
                  model = model_can, y = yNA,
                  omics_data = if (exists("omics_data")) omics_data else NULL,
                  tst = tst, additional_params = params))
              }
            }, error = function(e) {
              log_error("AI model {model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
              handle_error <<- TRUE
            })
          }

          if (isTRUE(verbose)) {
            n_now <- sum(!is.na(ypred_cv$yhat[tst]))
            if (n_now != length(tst)) log_warn("Fold {j}/{effective_nfolds}: predicted {n_now} of {length(tst)}")
          }
        }
      }
      fold_partial <- finish_fold_phase()
      if (!is.null(fold_partial)) {
        return(fold_partial)
      }

      oof <- which(ypred_cv$cv_role == "test" & !is.na(ypred_cv$yhat))
      results_eval_metrics_reps <- data.frame(Rep = repp, check.names = FALSE)
      if (is.finite(feature_k)) {
        results_eval_metrics_reps$feature_k <- feature_k
        results_eval_metrics_reps$feature_scoring_model <- as.character(feature_scoring_model)
        results_eval_metrics_reps$feature_scoring_cv <- feature_selection_mode
      }
      for (m in eval_metrics) {
        results_eval_metrics_reps[[m]] <- tryCatch(
          safe_metric_value(
            y_true = ypred_cv$y[oof],
            y_pred = metric_pred_input(oof, m),
            metric = m,
            response_family = response_family,
            positive_class = positive_class
          ),
          error = function(e) NA_real_
        )
      }
      if (isTRUE(verbose) && length(oof) < sum(ypred_cv$cv_role == "test"))
        log_info("OOF coverage: {length(oof)}/{sum(ypred_cv$cv_role == 'test')} rows got predictions")
    }

    if (is_multi) {
      if (is.null(heter_groups)) {
        log_error("heter_groups required for CV0/CV1/CV2")
        stop("For CV0/CV1/CV2 you must supply `heter_groups`.", call. = FALSE)
      }
      if (!heter_groups %in% names(ypred_cv)) {
        ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
      }

      group <- rep(NA_integer_, nrow(pheno_data))
      group[valid_idx] <- test_set_val[[repp]]
      if (length(group) != nrow(pheno_data)) {
        log_error("Split length mismatch: got {length(group)}, expected {nrow(pheno_data)}")
        stop(sprintf("Split length mismatch: got %d, expected %d.", length(group), nrow(pheno_data)))
      }

      uf <- sort(unique(group))
      uf <- uf[uf > 0L]

      if (!length(uf) && isTRUE(verbose)) {
        log_info("No test folds (>0) produced by the ME scheme; all rows remain train")
      }
      gp_fold_list <- lapply(uf, function(j) which(group == j))
      gp_batched_done <- assign_gp_batch_preds(gp_fold_list, context = cross_validation_meth)
      met_ai_batched_done <- if (!isTRUE(gp_batched_done)) {
        assign_met_ai_batch_preds(gp_fold_list, context = cross_validation_meth)
      } else {
        FALSE
      }

      if (!isTRUE(gp_batched_done) && !isTRUE(met_ai_batched_done)) {
      for (j in uf) {
        tst <- which(group == j)
        if (!length(tst)) next
        if (!fit_fold(match(j, uf))) next
        seed_fold(match(j, uf))
        apply_fold_internal_features(tst, match(j, uf))
        yNA <- y
        yNA[tst] <- NA_real_
        ypred_cv$cv_role[tst] <- "test"

        if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            params$bayes_model <- if (model_can == "GBLUP_BRR") "BRR" else model_can
            params$bayes_trait <- trait
            params$ETA <- model_prep_all_bayes_cv_task[[model_can]][["bayes_ETA"]][["ETA"]]
            params$bayes_para <- model_prep_all_bayes_cv_task[[model_can]][["bayes_para"]]
            params$bayes_groups <- model_prep_all_bayes_cv_task[[model_can]][["bayes_groups"]]
            params$met_kernel_cv_meta <- model_prep_all_bayes_cv_task[[model_can]][["met_kernel_cv_meta"]]
            assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = params))
          }, error = function(e) {
            log_error("Bayes {model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
            # Phase 3.17: warn-once per (trait, model) so the silent
            # Bayes CV drop becomes visible.
            if (!isTRUE(bayes_warned[[paste(trait, model_can, sep="::")]])) {
              warning(sprintf("Bayes %s failed for trait '%s' (fold %s; same error in subsequent folds suppressed): %s",
                              model_can, trait, j, conditionMessage(e)),
                      call. = FALSE)
              bayes_warned[[paste(trait, model_can, sep="::")]] <<- TRUE
            }
            handle_error <<- TRUE
          })
        }

        if (!is.null(params$engine) && model_can == "GBLUP" && params$engine == "asreml") {
          tryCatch({
            params$response <- trait
            assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = params))
          }, error = function(e) {
            log_error("GBLUP(asreml) failed for {trait} fold {j}: {conditionMessage(e)}")
            handle_error <<- TRUE
          })
        }

        if (model_can %in% gp_valid_models && !isTRUE(gp_batched_done)) {
          tryCatch({
            assign_preds(tst, predict_with_model(model_can, y = yNA, tst = tst, additional_params = params))
          }, error = function(e) {
            log_error("{model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
            handle_error <<- TRUE
          })
        }

        if (model_can %in% AI_valid_models) {
          tryCatch({
            if (isTRUE(met_ml_dl) && model_can %in% c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "mlp", "ft_transformer", "saint", "tabnet", "moe")) {
              pred_met <- if (model_can %in% c("mlp", "ft_transformer", "saint", "tabnet", "moe")) {
                gp_met_dl_cv_predict(
                  model_type = model_can,
                  pheno_data = pheno_data,
                  response = trait,
                  gen_name = gen_name,
                  heter_groups = heter_groups,
                  response_family = response_family,
                  gmatrix = gmatrix_model_ready,
                  omic1_kernel = omic1_kernel_model_ready,
                  omic2_kernel = omic2_kernel_model_ready,
                  omic3_kernel = omic3_kernel_model_ready,
                  kernel_list = kernel_list,
                  omics_kernel_label = omics_kernel_label,
                  tst = tst,
                  model_params = gp_met_dl_model_params(model_can, params),
                  met_kernel_var_explained = met_kernel_var_explained,
                  met_kernel_min_ev = met_kernel_min_ev,
                  met_kernel_max_pcs = met_kernel_max_pcs
                )
              } else {
                gp_met_python_cv_predict(
                  pheno_data = pheno_data,
                  response = trait,
                  gen_name = gen_name,
                  heter_groups = heter_groups,
                  response_family = response_family,
                  model_type = model_can,
                  gmatrix = gmatrix_model_ready,
                  omic1_kernel = omic1_kernel_model_ready,
                  omic2_kernel = omic2_kernel_model_ready,
                  omic3_kernel = omic3_kernel_model_ready,
                  kernel_list = kernel_list,
                  omics_kernel_label = omics_kernel_label,
                  tst = tst,
                  model_params = gp_met_python_model_params(model_can, params),
                  met_kernel_var_explained = met_kernel_var_explained,
                  met_kernel_min_ev = met_kernel_min_ev,
                  met_kernel_max_pcs = met_kernel_max_pcs
                )
              }
              assign_preds(tst, pred_met)
            } else {
              assign_preds(tst, predict_with_model(
                model = model_can, y = yNA,
                omics_data = if (exists("omics_data")) omics_data else NULL,
                tst = tst, additional_params = params))
            }
          }, error = function(e) {
            log_error("AI model {model_can} failed for {trait} fold {j}: {conditionMessage(e)}")
            handle_error <<- TRUE
          })
        }

        if (isTRUE(verbose)) {
          n_now <- sum(!is.na(ypred_cv$yhat[tst]))
          if (n_now != length(tst))
            log_warn("Fold {j}: predicted {n_now} of {length(tst)}")
        }
      }
      }
      fold_partial <- finish_fold_phase()
      if (!is.null(fold_partial)) {
        return(fold_partial)
      }

      ypred_cv$yhat[ypred_cv$cv_role != "test"] <- NA_real_
      if (!any(ypred_cv$cv_role == "test")) {
        log_error("No 'test' rows marked in CV")
        stop("[cv] No 'test' rows were marked. Check that group_vec has valid fold indices.", call. = FALSE)
      }

      test_envs <- sort(unique(ypred_cv[[heter_groups]][ypred_cv$cv_role == "test"]))
      by_env <- do.call(
        rbind,
        lapply(test_envs, function(env_i) {
          idx <- which(ypred_cv$cv_role == "test" & ypred_cv[[heter_groups]] == env_i)
          row <- list(Rep = repp)
          row[[heter_groups]] <- env_i
          if (is.finite(feature_k)) {
            row[["feature_k"]] <- feature_k
            row[["feature_scoring_model"]] <- as.character(feature_scoring_model)
            row[["feature_scoring_cv"]] <- feature_selection_mode
          }
          for (m in eval_metrics) {
            ok <- gp_metric_valid_rows(ypred_cv$y[idx], ypred_cv$yhat[idx], response_family = response_family)
            if (any(ok)) {
              pred_input <- metric_pred_input(idx, m)
              if (is.matrix(pred_input) || is.data.frame(pred_input)) {
                pred_input <- pred_input[ok, , drop = FALSE]
              } else {
                pred_input <- pred_input[ok]
              }
              row[[m]] <- safe_metric_value(
                y_true = ypred_cv$y[idx][ok], y_pred = pred_input, metric = m,
                response_family = response_family, positive_class = positive_class
              )
            } else {
              row[[m]] <- NA_real_
            }
          }
          as.data.frame(row, check.names = FALSE)
        })
      )
      oof <- which(ypred_cv$cv_role == "test")
      ok <- gp_metric_valid_rows(ypred_cv$y[oof], ypred_cv$yhat[oof], response_family = response_family)
      overall <- data.frame(Rep = repp, check.names = FALSE)
      overall[[heter_groups]] <- "ALL"
      if (is.finite(feature_k)) {
        overall$feature_k <- feature_k
        overall$feature_scoring_model <- as.character(feature_scoring_model)
        overall$feature_scoring_cv <- feature_selection_mode
      }
      for (m in eval_metrics) {
        if (any(ok)) {
          pred_input <- metric_pred_input(oof, m)
          if (is.matrix(pred_input) || is.data.frame(pred_input)) {
            pred_input <- pred_input[ok, , drop = FALSE]
          } else {
            pred_input <- pred_input[ok]
          }
          overall[[m]] <- safe_metric_value(
            ypred_cv$y[oof][ok], pred_input, m,
            response_family = response_family, positive_class = positive_class
          )
        } else {
          overall[[m]] <- NA_real_
        }
      }
      results_eval_metrics_reps <- rbind(by_env, overall)
    }

    list(
      trait = trait,
      rep = rep_i,
      model = model_label,
      model_canonical = model_can,
      response_family = response_family,
      class_levels = class_levels,
      positive_class = positive_class,
      eval_metrics_reps = results_eval_metrics_reps,
      results_eval_metrics_reps = results_eval_metrics_reps,
      ypred_cv_Reps_all = ypred_cv,
      yprob_cv_Reps_all = if (!is.null(yprob_cv)) as.data.frame(yprob_cv, stringsAsFactors = FALSE) else NULL,
      feature_selection_metadata = gp_feature_bind_metadata(feature_selection_records),
      cv_info = c(
        list(
          method = cross_validation_meth,
          token = cv_token,
          nfolds = if (is_kfold || is_loo) effective_nfolds else NA_integer_,
          replication = repp,
          n_test = sum(ypred_cv$cv_role == "test"),
          n_train = sum(ypred_cv$cv_role == "train"),
          feature_scoring_model = as.character(feature_scoring_model),
          feature_k = if (is.finite(feature_k)) feature_k else NA_integer_,
          feature_scoring_cv = feature_selection_mode
        ),
        gp_batch_cv_info %||% list()
      )
    )
  }

  # Merge fold-task partials back into one result per unsplit task, in the
  # original task order, by computing metrics on the pooled out-of-fold
  # predictions (the same code path as an unsplit task).
  merge_fold_task_results <- function(results) {
    if (all(is.na(tasks_full$fold_only))) {
      return(results)
    }
    policy <- attr(results, "gp_execution_policy", exact = TRUE)
    keys <- task_group_key(tasks_full)
    out <- lapply(seq_len(nrow(tasks_unsplit)), function(u) {
      base_row <- tasks_unsplit[u, , drop = FALSE]
      idx <- which(keys == task_group_key(base_row))
      if (length(idx) == 1L && is.na(tasks_full$fold_only[idx])) {
        return(results[[idx]])
      }
      parts <- Filter(function(p) is.list(p) && isTRUE(p$fold_partial), results[idx])
      if (!length(parts)) {
        return(NULL)
      }
      ypred <- parts[[1L]]$ypred_cv
      prob_cols <- unique(unlist(lapply(parts, function(p) colnames(p$yprob_cv))))
      yprob <- if (length(prob_cols)) {
        matrix(NA_real_, nrow(ypred), length(prob_cols), dimnames = list(NULL, prob_cols))
      } else NULL
      for (p in parts) {
        test_rows <- which(p$ypred_cv$cv_role == "test")
        if (!length(test_rows)) next
        ypred$cv_role[test_rows] <- "test"
        ypred$yhat[test_rows] <- p$ypred_cv$yhat[test_rows]
        if (!is.null(yprob) && !is.null(p$yprob_cv)) {
          cols <- colnames(p$yprob_cv)
          yprob[test_rows, cols] <- as.matrix(p$yprob_cv)[test_rows, cols, drop = FALSE]
        }
      }
      run_one_task(base_row, merged = list(
        ypred_cv = ypred,
        yprob_cv = yprob,
        feature_selection_records = unlist(lapply(parts, `[[`, "feature_selection_records"), recursive = FALSE),
        gp_batch_cv_info = parts[[1L]]$gp_batch_cv_info
      ))
    })
    attr(out, "gp_execution_policy") <- policy
    out
  }

  if (isTRUE(gp_cv_direct_bridge_only) && nrow(tasks_full) == 1L) {
    results <- list(run_one_task(tasks_full[1L, , drop = FALSE]))
    attr(results, "gp_execution_policy") <- list(
      backend = "sequential",
      workers = 1L,
      plan = NULL,
      chunk_size = NULL,
      backend_score = Inf,
      decision_reason = "single_direct_gp_cv_task",
      num_gpus = NA_integer_,
      model_ids = tasks_full$modell,
      user_sequential_models = unique(sequential_models %||% character()),
      internal_flags = TRUE,
      internal_sequential_count = 1L,
      parallel_count = 0L
    )
  } else {
    results <- gp_execute_partitioned_tasks(
      model_ids = tasks_full$modell,
      run_task = function(i) run_one_task(tasks_full[i, , drop = FALSE]),
      sequential_models = sequential_models,
      decision_args = list(
        model_params_list = policy_model_params,
        # What workers actually receive: run_task serialises this whole frame.
        # The former hand-picked list sized 0.81 GB of a 1.79 GB payload (mice
        # CV), and, kept in the frame, was itself shipped (826 MB). Built
        # inline so it is not bound here; shares the frame's objects (no copy).
        globals = as.list(environment(), all.names = FALSE),
        user_mode = parallel_mode,
        num_cores = num_cores,
        globals_max_GB = params$globals_max_GB %||% as.numeric(Sys.getenv("GP_GLOBALS_MAX_GB", "4")),
        worker_memory_gb = params$worker_memory_gb,
        memory_budget_gb = params$memory_budget_gb,
        task_cost = cv_task_cost,
        verbose = verbose,
        sys_name = sys_name,
        prefer_fork = gp_parallel_safe_fork_preference(
          prefer_fork = parallel_backend_prefer_fork,
          needs_process_initializer = cv_scheduler_needs_python_init
        )
      ),
      packages = c("dplyr", "logger"),
      seed = TRUE,
      initializer = if (isTRUE(cv_scheduler_needs_python_init)) init_py else NULL,
      fail_fast = params$fail_fast %||% FALSE
    )
    results <- merge_fold_task_results(results)
  }

  filtered_results <- gp_filter_crossval_results(results)
  attr(filtered_results, "gp_execution_policy") <- attr(results, "gp_execution_policy")

  log_info("Cross-validation completed with {length(filtered_results)} valid results")
  return(filtered_results)
}
