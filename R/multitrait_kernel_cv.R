gp_multitrait_joint_cv_fold_plan <- function(pheno_data,
                                              response,
                                              gen_name,
                                              nfolds = 5L,
                                              replication = 1L,
                                              random_seed = 123L,
                                              min_train_per_trait = 2L,
                                              max_attempts = 200L) {
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  if (!gen_name %in% names(ph) || !all(response %in% names(ph))) {
    stop("Joint multi-trait CV could not find the genotype or response columns.", call. = FALSE)
  }
  ids <- as.character(ph[[gen_name]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(
      "Joint multi-trait single-environment CV requires one non-missing phenotype row per genotype.",
      call. = FALSE
    )
  }
  observed <- vapply(
    ph[, response, drop = FALSE],
    function(x) is.finite(suppressWarnings(as.numeric(x))),
    logical(nrow(ph))
  )
  if (is.null(dim(observed))) {
    observed <- matrix(observed, ncol = length(response))
  }
  colnames(observed) <- response
  eligible <- which(rowSums(observed) > 0L)
  if (length(eligible) < 2L) {
    stop("Joint multi-trait CV requires at least two genotypes with observed traits.", call. = FALSE)
  }

  nfolds <- suppressWarnings(as.integer(nfolds)[1L])
  replication <- suppressWarnings(as.integer(replication)[1L])
  min_train_per_trait <- suppressWarnings(as.integer(min_train_per_trait)[1L])
  max_attempts <- suppressWarnings(as.integer(max_attempts)[1L])
  if (!is.finite(nfolds) || nfolds < 2L || !is.finite(replication) || replication < 1L) {
    stop("Joint multi-trait CV requires nfolds >= 2 and replication >= 1.", call. = FALSE)
  }
  k <- min(nfolds, length(eligible))
  trait_counts <- colSums(observed[eligible, , drop = FALSE])
  if (any(trait_counts <= min_train_per_trait)) {
    bad <- response[trait_counts <= min_train_per_trait]
    stop(
      "Joint multi-trait CV cannot leave at least ", min_train_per_trait,
      " training observations for every trait. Too-sparse traits: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_existed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (isTRUE(seed_existed) && is.integer(old_seed) && length(old_seed) > 1L) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  plans <- vector("list", replication)
  manifest <- vector("list", replication)
  base_seed <- suppressWarnings(as.integer(random_seed %||% 123L)[1L])
  if (!is.finite(base_seed)) base_seed <- 123L
  for (rep_idx in seq_len(replication)) {
    accepted <- NULL
    for (attempt in seq_len(max_attempts)) {
      gp_set_seed(as.integer(base_seed + (rep_idx - 1L) * 1009L + attempt - 1L))
      shuffled <- sample(eligible, length(eligible), replace = FALSE)
      fold_id <- rep(seq_len(k), length.out = length(shuffled))
      candidate <- split(shuffled, fold_id)
      valid <- all(vapply(candidate, function(test_idx) {
        train_idx <- setdiff(eligible, test_idx)
        length(test_idx) > 0L &&
          all(colSums(observed[train_idx, , drop = FALSE]) >= min_train_per_trait)
      }, logical(1L)))
      if (isTRUE(valid)) {
        accepted <- lapply(candidate, sort)
        break
      }
    }
    if (is.null(accepted)) {
      counts <- paste(paste0(response, "=", trait_counts), collapse = ", ")
      stop(
        "Joint multi-trait CV could not create leakage-free genotype folds while retaining ",
        min_train_per_trait, " training observations per trait after ", max_attempts,
        " attempts. Observed trait counts: ", counts,
        ". Reduce nfolds or provide more overlapping trait observations.",
        call. = FALSE
      )
    }
    plans[[rep_idx]] <- accepted
    manifest[[rep_idx]] <- do.call(rbind, lapply(seq_along(accepted), function(fold_idx) {
      idx <- accepted[[fold_idx]]
      data.frame(
        rep = rep_idx,
        fold = fold_idx,
        GID = ids[idx],
        observed_trait_count = rowSums(observed[idx, , drop = FALSE]),
        stringsAsFactors = FALSE
      )
    }))
  }
  list(
    folds = plans,
    manifest = do.call(rbind, manifest),
    observed = observed,
    eligible_rows = eligible,
    effective_nfolds = k
  )
}

gp_multitrait_joint_met_cv_fold_plan <- function(pheno_data,
                                                  response,
                                                  gen_name,
                                                  heter_groups,
                                                  cross_validation_meth,
                                                  nfolds = 5L,
                                                  replication = 1L,
                                                  random_seed = 123L,
                                                  min_train_per_trait = 2L) {
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  required <- c(gen_name, heter_groups, response)
  missing <- setdiff(required, names(ph))
  if (length(missing)) {
    stop(
      "Joint MT-MET CV could not find required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  ids <- as.character(ph[[gen_name]])
  envs <- as.character(ph[[heter_groups]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyNA(envs) || any(!nzchar(envs))) {
    stop("Joint MT-MET CV requires non-missing genotype and environment identifiers.", call. = FALSE)
  }
  cell_key <- paste(ids, envs, sep = "\r")
  if (anyDuplicated(cell_key)) {
    stop(
      "Joint MT-MET CV requires one phenotype row per genotype-by-environment cell.",
      call. = FALSE
    )
  }
  if (length(unique(envs)) < 2L) {
    stop("Joint MT-MET CV requires at least two environments.", call. = FALSE)
  }

  observed <- vapply(
    ph[, response, drop = FALSE],
    function(x) is.finite(suppressWarnings(as.numeric(x))),
    logical(nrow(ph))
  )
  if (is.null(dim(observed))) {
    observed <- matrix(observed, ncol = length(response))
  }
  colnames(observed) <- response
  if (any(colSums(observed) <= min_train_per_trait)) {
    bad <- response[colSums(observed) <= min_train_per_trait]
    stop(
      "Joint MT-MET CV cannot retain at least ", min_train_per_trait,
      " training observations for traits: ", paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }

  assignments <- gp_met_cv_assignments(
    pheno_data = ph,
    gen_name = gen_name,
    heter_groups = heter_groups,
    cross_validation_meth = cross_validation_meth,
    nfolds = nfolds,
    random_state = random_seed,
    replication = replication,
    message = FALSE
  )
  folds <- vector("list", length(assignments))
  manifests <- vector("list", length(assignments))
  for (rep_idx in seq_along(assignments)) {
    assignment <- suppressWarnings(as.integer(assignments[[rep_idx]]))
    fold_ids <- sort(unique(assignment[is.finite(assignment) & assignment > 0L]))
    if (!length(fold_ids)) {
      stop("Joint MT-MET CV produced no test folds in replication ", rep_idx, ".", call. = FALSE)
    }
    folds[[rep_idx]] <- lapply(fold_ids, function(fold_id) which(assignment == fold_id))
    for (fold_idx in seq_along(folds[[rep_idx]])) {
      test_idx <- folds[[rep_idx]][[fold_idx]]
      train_idx <- setdiff(seq_len(nrow(ph)), test_idx)
      train_counts <- colSums(observed[train_idx, , drop = FALSE])
      if (any(train_counts < min_train_per_trait)) {
        bad <- response[train_counts < min_train_per_trait]
        stop(
          "Joint MT-MET ", cross_validation_meth, " cannot retain at least ",
          min_train_per_trait, " training observations for every trait in replication ",
          rep_idx, ", fold ", fold_idx, ". Too-sparse traits: ",
          paste(bad, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
    manifests[[rep_idx]] <- do.call(rbind, lapply(seq_along(folds[[rep_idx]]), function(fold_idx) {
      idx <- folds[[rep_idx]][[fold_idx]]
      data.frame(
        rep = rep_idx,
        fold = fold_idx,
        row_index = idx,
        GID = ids[idx],
        Env = envs[idx],
        observed_trait_count = rowSums(observed[idx, , drop = FALSE]),
        stringsAsFactors = FALSE
      )
    }))
  }
  list(
    folds = folds,
    manifest = do.call(rbind, manifests),
    observed = observed,
    effective_nfolds = max(vapply(folds, length, integer(1L))),
    cv_token = normalize_cv_token(cross_validation_meth)
  )
}

gp_multitrait_joint_cv_prediction_table <- function(predictions,
                                                     target,
                                                     model_type,
                                                     rep_idx,
                                                     fold_idx,
                                                     gen_name = "GID",
                                                     heter_groups = NULL) {
  pred <- as.data.frame(predictions, stringsAsFactors = FALSE)
  id_col <- unique(c("GID", "gid", gen_name, "Name", "name", "ID", "id"))
  id_col <- id_col[id_col %in% names(pred)][1L]
  trait_col <- c("Trait", "trait")
  trait_col <- trait_col[trait_col %in% names(pred)][1L]
  value_col <- c("Predicted_value", "Prediction", "prediction", "yhat")
  value_col <- value_col[value_col %in% names(pred)][1L]
  if (!length(id_col) || !length(trait_col) || !length(value_col) ||
      is.na(id_col) || is.na(trait_col) || is.na(value_col)) {
    stop(
      "Joint multi-trait CV fold predictions are missing genotype, trait, or predicted-value columns.",
      call. = FALSE
    )
  }

  include_env <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    !is.na(heter_groups) && nzchar(heter_groups)
  if (isTRUE(include_env)) {
    pred_env_col <- c(heter_groups, "Env", "env")
    pred_env_col <- pred_env_col[pred_env_col %in% names(pred)][1L]
    target_env_col <- c(heter_groups, "Env", "env")
    target_env_col <- target_env_col[target_env_col %in% names(target)][1L]
    if (!length(pred_env_col) || !length(target_env_col) ||
        is.na(pred_env_col) || is.na(target_env_col)) {
      stop("Joint MT-MET CV fold predictions are missing the environment column.", call. = FALSE)
    }
    pred_key <- paste(
      as.character(pred[[id_col]]), as.character(pred[[pred_env_col]]),
      as.character(pred[[trait_col]]), sep = "\r"
    )
    target_key <- paste(
      as.character(target$GID), as.character(target[[target_env_col]]),
      as.character(target$Trait), sep = "\r"
    )
  } else {
    pred_key <- paste(as.character(pred[[id_col]]), as.character(pred[[trait_col]]), sep = "\r")
    target_key <- paste(as.character(target$GID), as.character(target$Trait), sep = "\r")
  }
  relevant <- pred_key %in% target_key
  if (anyDuplicated(pred_key[relevant])) {
    stop(
      "Joint multi-trait CV produced duplicate genotype-trait predictions in rep ",
      rep_idx, ", fold ", fold_idx, ".",
      call. = FALSE
    )
  }
  match_idx <- match(target_key, pred_key)
  if (anyNA(match_idx)) {
    missing <- target_key[is.na(match_idx)]
    stop(
      "Joint multi-trait CV did not return every held-out observed genotype-trait prediction in rep ",
      rep_idx, ", fold ", fold_idx, ". Missing: ",
      paste(utils::head(gsub("\r", "/", missing, fixed = TRUE), 6L), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  predicted <- suppressWarnings(as.numeric(pred[[value_col]][match_idx]))
  if (any(!is.finite(predicted))) {
    stop(
      "Joint multi-trait CV returned non-finite held-out predictions in rep ",
      rep_idx, ", fold ", fold_idx, ".",
      call. = FALSE
    )
  }

  out <- data.frame(
    GID = as.character(target$GID),
    Trait = as.character(target$Trait),
    model = as.character(model_type),
    Observed_value = suppressWarnings(as.numeric(target$Observed_value)),
    Predicted_value = predicted,
    Train_Test_Label = "Test",
    fold = as.integer(fold_idx),
    rep = as.integer(rep_idx),
    cv_role = "test",
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_env)) {
    out[[as.character(heter_groups)]] <- as.character(target[[target_env_col]])
    leading <- c("GID", as.character(heter_groups), setdiff(names(out), c("GID", as.character(heter_groups))))
    out <- out[, leading, drop = FALSE]
  }
  optional <- c(
    "Standard_error", "PEV", "PEV_basis", "Prediction_uncertainty_source",
    "Prediction_interval_method", "Prediction_interval_nominal_coverage",
    "Prediction_interval_calibration_n", "lower_bound", "upper_bound",
    "Uncertainty", "Uncertainty_remarks", "Reliability",
    "Reliability_variance_input", "Reliability_reference_variance",
    "Reliability_basis", "Reliability_remarks"
  )
  for (column in optional[optional %in% names(pred)]) {
    out[[column]] <- pred[[column]][match_idx]
  }
  out
}

gp_multitrait_joint_gaussian_cv <- function(pheno_object,
                                             response,
                                             gen_name,
                                             model_type,
                                             fit_fold,
                                             heter_groups = NULL,
                                             cross_validation_meth = "K-Folds",
                                             nfolds = 5L,
                                             replication = 1L,
                                             eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                             random_seed = 123L,
                                             kernel_names = NULL,
                                             kernel_strategy = NA_character_,
                                             kernel_role = "independent_genetic_covariance_term") {
  cv_token <- normalize_cv_token(cross_validation_meth)
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  is_met <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    !is.na(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph) &&
    length(unique(stats::na.omit(as.character(ph[[heter_groups]])))) > 1L
  if (isTRUE(is_met)) {
    if (!cv_token %in% c("cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2")) {
      stop(
        "Joint MT-MET GP/Bayesian CV supports CV0, CV1, CV2, Repeated_CV0, Repeated_CV1, and Repeated_CV2 only.",
        call. = FALSE
      )
    }
    plan <- gp_multitrait_joint_met_cv_fold_plan(
      pheno_data = ph,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      random_seed = random_seed
    )
  } else {
    if (!cv_token %in% c("k_folds", "repeated_k_folds")) {
      stop(
        "Joint multi-trait single-environment GP/Bayesian CV supports K-Folds and Repeated_K-Folds only.",
        call. = FALSE
      )
    }
    plan <- gp_multitrait_joint_cv_fold_plan(
      pheno_data = ph,
      response = response,
      gen_name = gen_name,
      nfolds = nfolds,
      replication = replication,
      random_seed = random_seed
    )
  }
  ids <- as.character(ph[[gen_name]])
  cv_results <- vector("list", length(plan$folds))
  for (rep_idx in seq_along(plan$folds)) {
    pred_blocks <- vector("list", length(plan$folds[[rep_idx]]))
    for (fold_idx in seq_along(plan$folds[[rep_idx]])) {
      test_idx <- plan$folds[[rep_idx]][[fold_idx]]
      test_ids <- ids[test_idx]
      ph_masked <- ph
      ph_masked[test_idx, response] <- NA_real_
      fold_seed <- as.integer((random_seed %||% 123L) + rep_idx * 10000L + fold_idx)
      fold_predictions <- tryCatch(
        fit_fold(
          pheno_masked = ph_masked,
          test_ids = test_ids,
          fold_seed = fold_seed,
          rep_idx = rep_idx,
          fold_idx = fold_idx
        ),
        error = function(e) {
          stop(
            "Joint multi-trait CV fit failed for model `", model_type,
            "` in rep ", rep_idx, ", fold ", fold_idx, ": ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      target_parts <- lapply(response, function(trait) {
        keep <- is.finite(suppressWarnings(as.numeric(ph[[trait]][test_idx])))
        if (!any(keep)) return(NULL)
        data.frame(
          GID = test_ids[keep],
          Trait = trait,
          Observed_value = suppressWarnings(as.numeric(ph[[trait]][test_idx][keep])),
          stringsAsFactors = FALSE
        )
      })
      if (isTRUE(is_met)) {
        target_parts <- lapply(seq_along(response), function(trait_index) {
          trait <- response[[trait_index]]
          keep <- is.finite(suppressWarnings(as.numeric(ph[[trait]][test_idx])))
          if (!any(keep)) return(NULL)
          out <- data.frame(
            GID = test_ids[keep],
            Trait = trait,
            Observed_value = suppressWarnings(as.numeric(ph[[trait]][test_idx][keep])),
            stringsAsFactors = FALSE
          )
          out[[as.character(heter_groups)]] <- as.character(ph[[heter_groups]][test_idx][keep])
          out[, c("GID", as.character(heter_groups), "Trait", "Observed_value"), drop = FALSE]
        })
      }
      target_parts <- Filter(Negate(is.null), target_parts)
      target <- do.call(rbind, target_parts)
      pred_blocks[[fold_idx]] <- gp_multitrait_joint_cv_prediction_table(
        predictions = fold_predictions,
        target = target,
        model_type = model_type,
        rep_idx = rep_idx,
        fold_idx = fold_idx,
        gen_name = gen_name,
        heter_groups = if (isTRUE(is_met)) heter_groups else NULL
      )
    }
    pred_long <- do.call(rbind, pred_blocks)
    rownames(pred_long) <- NULL
    held_out_rows <- sort(unique(unlist(plan$folds[[rep_idx]], use.names = FALSE)))
    observed_expected <- sum(plan$observed[held_out_rows, , drop = FALSE])
    prediction_key <- if (isTRUE(is_met)) {
      paste(pred_long$GID, pred_long[[heter_groups]], pred_long$Trait, sep = "\r")
    } else {
      paste(pred_long$GID, pred_long$Trait, sep = "\r")
    }
    if (nrow(pred_long) != observed_expected || anyDuplicated(prediction_key)) {
      stop(
        "Joint multi-trait CV out-of-fold coverage is incomplete or duplicated in rep ",
        rep_idx, ".",
        call. = FALSE
      )
    }
    metric_table <- if (isTRUE(is_met)) {
      metric_predictions <- pred_long
      metric_predictions$Env <- as.character(metric_predictions[[heter_groups]])
      gp_multitrait_asreml_mtmet_metric_table(
        pred_long = metric_predictions,
        eval_metrics = eval_metrics,
        model_type = model_type,
        rep = rep_idx
      )
    } else {
      gp_multitrait_ml_metric_table(
        pred_long = pred_long,
        eval_metrics = eval_metrics,
        model_type = model_type,
        rep = rep_idx
      )
    }
    cv_results[[rep_idx]] <- list(
      trait = paste(response, collapse = ","),
      model = model_type,
      rep = rep_idx,
      eval_metrics_reps = metric_table,
      ypred_cv_Reps_all = pred_long,
      yprob_cv_Reps_all = NULL,
      fold_assignments = plan$manifest[plan$manifest$rep == rep_idx, , drop = FALSE]
    )
  }

  processed <- gp_multitrait_ml_cv_process(cv_results)
  processed[["multitrait_oof_predictions"]] <- do.call(
    rbind,
    lapply(cv_results, `[[`, "ypred_cv_Reps_all")
  )
  processed[["fold_assignments"]] <- plan$manifest
  kernel_names <- as.character(kernel_names %||% character())
  if (length(kernel_names)) {
    processed[["kernel_configuration"]] <-
      gp_multitrait_kernel_configuration(kernel_names, role = kernel_role)
    processed[["model_parameters"]] <- gp_multi_kernel_parameter_rows(
      kernel_names = kernel_names,
      strategy = kernel_strategy
    )
  }
  if (isTRUE(is_met)) {
    metrics_all <- do.call(rbind, lapply(cv_results, `[[`, "eval_metrics_reps"))
    processed[["mtmet_metric_summary_by_environment"]] <- if (is.data.frame(metrics_all) && nrow(metrics_all)) {
      stats::aggregate(
        value ~ trait + Env + model + metric,
        data = metrics_all,
        FUN = function(x) mean(x, na.rm = TRUE)
      )
    } else {
      data.frame()
    }
  }
  processed[["cv_scheme"]] <- list(
    unit = if (!isTRUE(is_met)) {
      "genotype"
    } else if (cv_token %in% c("cv0", "repeated_cv0")) {
      "environment"
    } else if (cv_token %in% c("cv1", "repeated_cv1")) {
      "genotype"
    } else {
      "genotype_environment_cell"
    },
    all_traits_masked_together = TRUE,
    method = cv_token,
    environment_column = if (isTRUE(is_met)) heter_groups else NULL,
    requested_nfolds = as.integer(nfolds),
    effective_nfolds = plan$effective_nfolds,
    replication = as.integer(replication),
    variance_source = "model_refit_per_fold; no predictive variance relabeled as genetic variance"
  )
  list(cv_results = cv_results, cv_results_processed = processed)
}

gp_multitrait_gp_gaussian_cv <- function(model_type,
                                         pheno_object,
                                         response,
                                         gen_name,
                                         gmatrix,
                                         kernel_list = NULL,
                                         kernel_weights = NULL,
                                         estimate_kernel_weights = FALSE,
                                         force_prediction_se = FALSE,
                                         heter_groups = NULL,
                                         cross_validation_meth = "K-Folds",
                                         nfolds = 5L,
                                         replication = 1L,
                                         eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                         random_seed = 123L,
                                         gp_factor_cache = NULL,
                                         gp_fa_rank = 1L,
                                         gp_return_se = TRUE,
                                         gp_full_vc = FALSE,
                                         gp_return_trait_correlations = FALSE,
                                         gp_iters = NULL,
                                         para_tunning = FALSE,
                                         env_similarity = NULL,
                                         env_ids = NULL,
                                         env_covariates = NULL,
                                         reaction_norm_feature_qc = TRUE,
                                         kenv_kernel = "matern32",
                                         kenv_bandwidth = 1,
                                         kenv_kernel_kwargs = NULL) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  is_met <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    !is.na(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph) &&
    length(unique(stats::na.omit(as.character(ph[[heter_groups]])))) > 1L
  route <- gp_multitrait_gp_route_spec(
    model_name = model_type,
    is_met = is_met,
    fa_rank = gp_fa_rank
  )
  model_public <- gp_public_model_label(model_type)[1L]
  max_iter <- suppressWarnings(as.integer(gp_iters %||% 100L)[1L])
  if (!is.finite(max_iter) || max_iter < 1L) max_iter <- 100L
  fit_fold <- function(pheno_masked, test_ids, test_rows = NULL, test_envs = NULL,
                       fold_seed, rep_idx, fold_idx) {
    fit <- if (isTRUE(is_met)) {
      gp_multi_trait_met_model(
        force_prediction_se = isTRUE(force_prediction_se),
        pheno_data = pheno_masked,
        gmatrix = gmatrix,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups,
        test_set = NULL,
        kernel_list = kernel_list,
        kernel_weights = kernel_weights,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        varcomp_mode = route$varcomp_mode,
        gp_factor_cache = gp_factor_cache,
        tune_gp = if (isTRUE(para_tunning)) "auto" else FALSE,
        return_se = isTRUE(gp_return_se) || isTRUE(gp_full_vc),
        prediction_output = "test_only",
        return_trait_correlations = isTRUE(gp_return_trait_correlations),
        trait_structure = route$trait_structure,
        trait_fa_rank = route$trait_fa_rank,
        gxe_trait_structure = route$gxe_trait_structure,
        gxe_trait_fa_rank = route$gxe_trait_fa_rank,
        seed = fold_seed
      )
    } else {
      gp_multi_trait_model(
        estimate_kernel_weights = isTRUE(estimate_kernel_weights),
        force_prediction_se = isTRUE(force_prediction_se),
        pheno_data = pheno_masked,
        gmatrix = gmatrix,
        response = response,
        gen_name = gen_name,
        test_set = test_ids,
        kernel_list = kernel_list,
        kernel_weights = kernel_weights,
        varcomp_mode = route$varcomp_mode,
        gp_factor_cache = gp_factor_cache,
        tune_gp = if (isTRUE(para_tunning)) "auto" else FALSE,
        return_se = isTRUE(gp_return_se) || isTRUE(gp_full_vc),
        prediction_output = "test_only",
        return_trait_correlations = isTRUE(gp_return_trait_correlations),
        trait_structure = route$trait_structure,
        trait_fa_rank = route$trait_fa_rank,
        max_iter = max_iter,
        seed = fold_seed
      )
    }
    gp_model_execute_multitrait_predicted_values(
      predictions = fit[["predictions"]],
      pheno_data = pheno_masked,
      response = response,
      gen_name = gen_name,
      heter_groups = if (isTRUE(is_met)) heter_groups else NULL
    )
  }
  gp_mixture_fitted <- isTRUE(estimate_kernel_weights) && !isTRUE(is_met) &&
    length(gp_collect_kernel_inputs(gmatrix = gmatrix, kernel_list = kernel_list)) > 1L
  if (isTRUE(estimate_kernel_weights) && isTRUE(is_met)) {
    warning(
      "estimate_kernel_weights = TRUE is not available for joint multi-trait MET GP; ",
      "the MET backend fits variance components by method of moments and carries the ",
      "kernel mixture as fixed weights. The MET folds used fixed weights.",
      call. = FALSE
    )
  }
  out <- gp_multitrait_joint_gaussian_cv(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    model_type = model_public,
    fit_fold = fit_fold,
    heter_groups = if (isTRUE(is_met)) heter_groups else NULL,
    cross_validation_meth = cross_validation_meth,
    nfolds = nfolds,
    replication = replication,
    eval_metrics = eval_metrics,
    random_seed = random_seed,
    kernel_names = names(gp_collect_kernel_inputs(
      gmatrix = gmatrix, kernel_list = kernel_list
    )),
    # With estimation off, mt_gp.py / mt_gp_met.py sum the aligned bank into a
    # single K_geno using fixed kernel_weights, so the per-kernel genetic
    # variances are not separately identified. With estimation on (single
    # environment only) the mixture is fitted by REML profile likelihood.
    kernel_strategy = if (isTRUE(gp_mixture_fitted)) {
      "reml_estimated_weight_combined_kernel"
    } else {
      "fixed_weight_combined_kernel"
    },
    kernel_role = if (isTRUE(gp_mixture_fitted)) {
      "estimated_weight_component_of_combined_kernel"
    } else {
      "fixed_weight_component_of_combined_kernel"
    }
  )
  out$cv_results_processed$gp_route <- route
  out
}

gp_multitrait_bayes_gaussian_cv <- function(model_type,
                                            pheno_object,
                                            response,
                                            gen_name,
                                            kernels,
                                            bayes_para,
                                            heter_groups = NULL,
                                            cross_validation_meth = "K-Folds",
                                            nfolds = 5L,
                                            replication = 1L,
                                            eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                            random_seed = 123L,
                                            heter_resid = FALSE,
                                            confidence_level = 0.95,
                                            CI_width_thresholds = c(0.33, 0.66),
                                            high_reliability_thres = 0.9,
                                            low_reliability_thres = 0.5) {
  model_type <- gp_canonicalize_supported_model_names(model_type)
  if (length(model_type) != 1L || !model_type %in% gp_multitrait_bayes_supported_models()) {
    stop(
      "Multi-trait Bayesian CV supports: ",
      paste(gp_multitrait_bayes_supported_models(), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  model_public <- gp_public_model_label(model_type)[1L]
  fit_fold <- function(pheno_masked, test_ids, test_rows = NULL, test_envs = NULL,
                       fold_seed, rep_idx, fold_idx) {
    seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (seed_existed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (isTRUE(seed_existed) && is.integer(old_seed) && length(old_seed) > 1L) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    gp_set_seed(fold_seed)
    fit <- bayes_multitrait_joint_fit(
      pheno_data = pheno_masked,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups,
      kernels = kernels,
      GS_model = model_type,
      bayes_para = bayes_para,
      heter_resid = heter_resid,
      confidence_level = confidence_level,
      CI_width_thresholds = CI_width_thresholds,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      verbose = FALSE
    )
    fit[["bayes_result"]][["predicted_values"]]
  }
  gp_multitrait_joint_gaussian_cv(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    model_type = model_public,
    fit_fold = fit_fold,
    heter_groups = heter_groups,
    cross_validation_meth = cross_validation_meth,
    nfolds = nfolds,
    replication = replication,
    eval_metrics = eval_metrics,
    random_seed = random_seed,
    kernel_names = names(kernels),
    # BGLR receives one ETA term per kernel, so each kernel's covariance
    # contribution is fitted independently.
    kernel_strategy = "independent_bayesian_eta_term_per_kernel",
    kernel_role = "independent_genetic_covariance_term"
  )
}
