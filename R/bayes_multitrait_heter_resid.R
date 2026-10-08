bayes_make_env_trait_matrix <- function(pheno_data,
                                        response,
                                        gen_name,
                                        heter_groups,
                                        target_ids = NULL) {
  if (is.null(target_ids)) {
    target_ids <- unique(as.character(pheno_data[[gen_name]]))
  }
  env_levels <- unique(as.character(pheno_data[[heter_groups]]))
  y_mat <- matrix(
    NA_real_,
    nrow = length(target_ids),
    ncol = length(env_levels),
    dimnames = list(as.character(target_ids), env_levels)
  )

  key <- paste(as.character(pheno_data[[gen_name]]), as.character(pheno_data[[heter_groups]]), sep = "\r")
  if (anyDuplicated(key)) {
    warning(
      "BGLR Multitrait heterogeneous residual path found replicated genotype-environment records; ",
      "phenotypes were averaged before fitting the environment-as-trait model.",
      call. = FALSE
    )
  }
  agg <- stats::aggregate(
    pheno_data[[response]],
    by = list(
      gid = as.character(pheno_data[[gen_name]]),
      env = as.character(pheno_data[[heter_groups]])
    ),
    FUN = function(x) mean(as.numeric(x), na.rm = TRUE)
  )
  names(agg)[3] <- "value"
  agg$value[is.nan(agg$value)] <- NA_real_
  keep <- agg$gid %in% rownames(y_mat) & agg$env %in% colnames(y_mat)
  if (any(keep)) {
    y_mat[cbind(agg$gid[keep], agg$env[keep])] <- agg$value[keep]
  }
  y_mat
}

bayes_rkhs_standardize_heterogeneous_variance_components <- function(x) {
  if (!is.data.frame(x) ||
      !all(c("Env", "Component", "Components", "Standard_error") %in% names(x))) {
    return(x)
  }
  env <- as.character(x[["Env"]])
  component <- as.character(x[["Component"]])
  if (!length(env)) {
    return(data.frame(
      Component = character(),
      Components = numeric(),
      Standard_error = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  labels <- ifelse(
    env == "Across_environment",
    paste0("across_environment_", component),
    paste0(component, "_", env)
  )
  out <- data.frame(
    Component = labels,
    Components = suppressWarnings(as.numeric(x[["Components"]])),
    Standard_error = suppressWarnings(as.numeric(x[["Standard_error"]])),
    stringsAsFactors = FALSE
  )
  rownames(out) <- make.unique(labels)
  out
}

bayes_multitrait_env_heter_fit <- function(pheno_data,
                                           response,
                                           gen_name,
                                           heter_groups,
                                           kernels,
                                           GS_model = "RKHS",
                                           bayes_para,
                                           omics_kernel_label = NULL,
                                           confidence_level = 0.95,
                                           CI_width_thresholds = c(0.33, 0.66),
                                           high_reliability_thres = 0.9,
                                           low_reliability_thres = 0.5) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR is required for Bayesian Multitrait heterogeneous residual fits.", call. = FALSE)
  }
  if (is.null(kernels) || !length(kernels)) {
    stop("BGLR Multitrait heterogeneous residual RKHS requires at least one kernel matrix.", call. = FALSE)
  }
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) {
    stop("BGLR Multitrait heterogeneous residual fits require heter_groups.", call. = FALSE)
  }
  if (!gp_bayes_has_multiple_residual_groups(pheno_data, heter_groups)) {
    stop("BGLR Multitrait heterogeneous residual fits require at least two residual groups.", call. = FALSE)
  }
  if (!GS_model %in% c("RKHS", "BRR")) {
    stop("BGLR Multitrait heterogeneous residual fallback supports GS_model 'RKHS' or 'BRR' (alias 'GBLUP_BRR').", call. = FALSE)
  }
  # Kernel-prior fit (model="RKHS" in BGLR::Multitrait) is mathematically the
  # multivariate generalisation of GBLUP_BRR/BRR with X = eigen-sqrt(K): both
  # induce the same N(0, sigma^2 * K) prior on genetic effects. Routing
  # GBLUP_BRR MET through this path yields identical variance components to
  # RKHS MET on the same kernel, while sidestepping BGLR's BRR-vs-RKHS prior-
  # scale calibration difference that depressed BRR-only MET fits.
  GS_model_label <- if (identical(GS_model, "BRR")) "GBLUP_BRR" else "RKHS"

  target_ids <- Reduce(intersect, lapply(kernels, rownames))
  if (!length(target_ids)) {
    stop("Kernel matrices do not share genotype row names for BGLR Multitrait fitting.", call. = FALSE)
  }
  pheno_ids <- unique(as.character(pheno_data[[gen_name]]))
  target_ids <- unique(c(pheno_ids[pheno_ids %in% target_ids], setdiff(target_ids, pheno_ids)))
  y_mat <- bayes_make_env_trait_matrix(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    target_ids = target_ids
  )

  # Phase 3.20 fix: BGLR::Multitrait initialises the residual covariance with
  # `cov(y, use="complete.obs")`, which throws "no complete element pairs" when
  # no row of y_mat has all columns non-NA. That happens routinely in MET CV2:
  # the scheduler masks at least one (gid, env) cell per genotype per fold, so
  # every partially-observed row has at least one NA, and the user's optional
  # by-GID masking adds fully-NA rows that don't help. We sidestep the failure
  # by appending k = ncol(y_mat) + 1 ballast genotypes whose phenotype rows are
  # drawn so that cov(ballast) approximates the pairwise-complete sample
  # covariance of the real data. Kernel entries are identity (uncorrelated
  # with real lines), so the ballast does not bias any real genotype's
  # posterior genetic effect. They are stripped from output after the fit.
  phantom_ids <- character(0L)
  if (sum(stats::complete.cases(y_mat)) == 0L) {
    n_env <- ncol(y_mat)
    n_phantom <- n_env + 1L
    col_means <- colMeans(y_mat, na.rm = TRUE)
    col_means[!is.finite(col_means)] <- 0
    # Use pairwise-complete-obs cov as the target for the phantom covariance
    # so BGLR's Sy_init is on the right scale. Anything pathological (NA,
    # non-finite, indefinite) falls back to a diagonal of the per-env SD^2.
    pair_cov <- suppressWarnings(stats::cov(y_mat, use = "pairwise.complete.obs"))
    col_vars <- diag(pair_cov)
    col_vars[!is.finite(col_vars) | col_vars <= 0] <- 1
    fallback_diag <- diag(col_vars, nrow = n_env)
    target_cov <- if (any(!is.finite(pair_cov))) fallback_diag else {
      # Project onto PSD: clip negative eigenvalues to a small positive floor.
      ev <- tryCatch(eigen(pair_cov, symmetric = TRUE), error = function(e) NULL)
      if (is.null(ev)) fallback_diag else {
        lam <- ev$values
        floor_val <- max(1e-8, max(lam) * 1e-6)
        lam[lam < floor_val] <- floor_val
        ev$vectors %*% diag(lam, n_env) %*% t(ev$vectors)
      }
    }
    L <- tryCatch(t(chol(target_cov)), error = function(e) NULL)
    if (is.null(L)) L <- diag(sqrt(col_vars), nrow = n_env)
    Z <- matrix(stats::rnorm(n_phantom * n_env), n_phantom, n_env)
    # Center Z so the ballast empirical mean is exactly col_means and
    # cov approximates target_cov (not exactly, but close for n_phantom = E+1).
    Z <- sweep(Z, 2L, colMeans(Z), `-`)
    ballast_mat <- Z %*% t(L) +
      matrix(col_means, nrow = n_phantom, ncol = n_env, byrow = TRUE)
    phantom_ids <- paste0("__BGLR_INIT_BALLAST_",
                          format(as.integer(Sys.time()), scientific = FALSE),
                          "_", sample.int(.Machine$integer.max, 1L), "_",
                          seq_len(n_phantom))
    rownames(ballast_mat) <- phantom_ids
    colnames(ballast_mat) <- colnames(y_mat)
    y_mat <- rbind(y_mat, ballast_mat)
    kernels <- lapply(kernels, function(K) {
      K <- as.matrix(K)
      n_old <- nrow(K)
      K_new <- matrix(0, nrow = n_old + n_phantom, ncol = n_old + n_phantom)
      K_new[seq_len(n_old), seq_len(n_old)] <- K
      diag(K_new)[(n_old + 1L):(n_old + n_phantom)] <- 1
      rownames(K_new) <- colnames(K_new) <- c(rownames(K), phantom_ids)
      K_new
    })
    target_ids <- c(target_ids, phantom_ids)
  }

  build_eta <- function(cov_type = "UN") {
    ETA <- lapply(names(kernels), function(nm) {
      K <- as.matrix(kernels[[nm]][target_ids, target_ids, drop = FALSE])
      rownames(K) <- colnames(K) <- target_ids
      list(K = K, model = "RKHS", Cov = list(type = cov_type))
    })
    names(ETA) <- names(kernels)
    ETA
  }
  fit_with_covariance <- function(cov_type) {
    ETA <- build_eta(cov_type)
    save_prefix <- gp_bglr_save_prefix(
      model_name = paste0("BGLR_Multitrait_", GS_model_label, "_", cov_type),
      response = response %||% "trait"
    )
    cleanup_prefix <- function() {
      con <- tryCatch(showConnections(all = TRUE), error = function(e) NULL)
      if (!is.null(con) && nrow(con)) {
        desc <- as.character(con[, "description"])
        hit <- grepl(basename(save_prefix), desc, fixed = TRUE)
        if (any(hit)) {
          for (id in as.integer(rownames(con)[hit])) {
            try(close(getConnection(id)), silent = TRUE)
          }
        }
      }
      unlink(
        list.files(
          path = dirname(save_prefix),
          pattern = paste0("^", basename(save_prefix)),
          full.names = TRUE
        )
      )
      invisible(NULL)
    }
    fit <- NULL
    tryCatch(
      invisible(utils::capture.output({
        fit <- BGLR::Multitrait(
          y = y_mat,
          ETA = ETA,
          resCov = list(type = "DIAG"),
          nIter = bayes_para[["nIter"]],
          burnIn = bayes_para[["burnIn"]],
          thin = bayes_para[["thin"]],
          verbose = FALSE,
          saveAt = save_prefix
        )
      })),
      error = function(e) {
        cleanup_prefix()
        stop(e)
      }
    )
    output_files_names <- list.files(
      path = dirname(save_prefix),
      pattern = paste0("^", basename(save_prefix)),
      full.names = TRUE
    )
    list(
      fit = fit,
      ETA = ETA,
      cov_type = cov_type,
      output_files_names = output_files_names
    )
  }

  fit_attempt <- tryCatch(
    fit_with_covariance("UN"),
    error = function(e) e
  )
  if (inherits(fit_attempt, "error")) {
    warning(
      "BGLR Multitrait RKHS heterogeneous residual fit could not initialize ",
      "the unstructured genetic covariance; retrying with diagonal genetic ",
      "covariance. Original error: ",
      conditionMessage(fit_attempt),
      call. = FALSE
    )
    fit_attempt <- fit_with_covariance("DIAG")
  }
  fit <- fit_attempt$fit
  ETA <- fit_attempt$ETA
  cov_type_used <- fit_attempt$cov_type
  output_files_names <- fit_attempt$output_files_names
  on.exit(unlink(output_files_names), add = TRUE)

  # Phase 3.20 fix: strip the BGLR-init ballast genotypes from y_mat and from
  # all BGLR matrices before downstream summaries so real predictions are
  # untouched. Mu (global intercept) is unaffected; SD.mu has length = ncol(y)
  # so it does not carry phantom rows.
  if (length(phantom_ids) > 0L) {
    real_ids <- setdiff(rownames(y_mat), phantom_ids)
    y_mat <- y_mat[real_ids, , drop = FALSE]
    fit$ETAHat <- fit$ETAHat[real_ids, , drop = FALSE]
    if (!is.null(fit$SD.ETAHat)) {
      fit$SD.ETAHat <- fit$SD.ETAHat[real_ids, , drop = FALSE]
    }
    if (!is.null(fit$y) && is.matrix(fit$y) && nrow(fit$y) == length(real_ids) + length(phantom_ids)) {
      fit$y <- fit$y[real_ids, , drop = FALSE]
    }
    target_ids <- real_ids
  }

  pred_mat <- sweep(as.matrix(fit$ETAHat), 2, as.double(fit$mu), FUN = "+")
  env_levels <- colnames(pred_mat)
  eta_for_covariance <- lapply(ETA, function(term) {
    if (!is.null(term$K)) {
      term$K <- as.matrix(term$K)[target_ids, target_ids, drop = FALSE]
    }
    term
  })
  covariance_bundle <- bayes_multitrait_covariance_bundle(
    fit = fit,
    ETA = eta_for_covariance,
    term_roles = stats::setNames(rep("genetic", length(ETA)), names(ETA)),
    labels = env_levels,
    draw_indices = gp_bayes_saved_draw_indices(
      nIter = bayes_para[["nIter"]],
      burnIn = bayes_para[["burnIn"]],
      thin = bayes_para[["thin"]]
    )
  )
  kernel_contract_matrices <- lapply(kernels, function(K) {
    as.matrix(K)[target_ids, target_ids, drop = FALSE]
  })
  kernel_contract <- bayes_multitrait_kernel_covariance_contract(
    covariance_bundle = covariance_bundle,
    term_kernel = stats::setNames(names(ETA), names(ETA)),
    kernel_matrices = kernel_contract_matrices,
    labels = env_levels,
    is_met = FALSE,
    axis_label = "Env"
  )
  genetic_covariance <- covariance_bundle$main$covariance
  residual_covariance <- covariance_bundle$residual$covariance
  sd_mu <- as.double(fit$SD.mu %||% rep(0, length(env_levels)))
  if (length(sd_mu) != length(env_levels)) {
    sd_mu <- rep(0, length(env_levels))
  }
  se_mat <- sqrt(as.matrix(fit$SD.ETAHat) ^ 2 +
                   matrix(sd_mu ^ 2, nrow = nrow(pred_mat), ncol = ncol(pred_mat), byrow = TRUE))
  pev_mat <- se_mat ^ 2
  alpha <- (1 - confidence_level) / 2
  z <- stats::qnorm(1 - alpha)

  grid <- expand.grid(
    gid = rownames(pred_mat),
    env = env_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx <- cbind(grid$gid, grid$env)
  observed <- y_mat[idx]
  pred <- pred_mat[idx]
  se <- se_mat[idx]
  pev <- pev_mat[idx]
  lower <- pred - z * se
  upper <- pred + z * se
  genetic_mean <- as.matrix(fit$ETAHat)[idx]
  env_genetic_var <- diag(genetic_covariance)
  names(env_genetic_var) <- env_levels
  genetic_var <- as.double(env_genetic_var[grid$env])
  genetic_var[!is.finite(genetic_var) | genetic_var <= 0] <- mean(genetic_var[is.finite(genetic_var) & genetic_var > 0], na.rm = TRUE)
  genetic_var[!is.finite(genetic_var)] <- 0
  rel <- reliability_thresholds(
    prediction_error_var = pev,
    genetic_var = genetic_var,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  int_stats <- bayes_classify_interval_widths(
    predictions = pred,
    lower_bound = lower,
    upper_bound = upper,
    CI_width_thresholds = CI_width_thresholds
  )

  predicted_value <- data.frame(
    name = grid$gid,
    Predicted_value = pred,
    Train_Test_Label = ifelse(is.na(observed), "Test", "Train"),
    Observed_value = observed,
    Env = grid$env,
    Standard_error = se,
    PEV = pev,
    lower_bound = lower,
    upper_bound = upper,
    Uncertainty = int_stats$Uncertainty,
    Uncertainty_remarks = int_stats$reliability_remarks,
    Genetic_variance = genetic_var,
    Reliability = rel$reliability,
    Reliability_remarks = rel$remarks,
    Reliability_percentage = rel$reliability_percentage,
    Prediction_target = "Environment_specific",
    stringsAsFactors = FALSE
  )
  names(predicted_value)[names(predicted_value) == "name"] <- gen_name
  names(predicted_value)[names(predicted_value) == "Env"] <- heter_groups

  residual_idx <- which(!is.na(observed))
  residual_value <- data.frame(
    name = grid$gid[residual_idx],
    Env = grid$env[residual_idx],
    Predicted_value = pred[residual_idx],
    Residual_value = observed[residual_idx] - pred[residual_idx],
    stringsAsFactors = FALSE
  )
  names(residual_value)[1:2] <- c(gen_name, heter_groups)

  n_env <- ncol(pred_mat)
  across_pred <- rowMeans(pred_mat, na.rm = TRUE)
  across_genetic <- rowMeans(as.matrix(fit$ETAHat), na.rm = TRUE)
  across_se <- sqrt(rowSums(se_mat ^ 2, na.rm = TRUE)) / n_env
  across_pev <- across_se ^ 2
  across_lower <- across_pred - z * across_se
  across_upper <- across_pred + z * across_se
  across_weights <- rep(1 / n_env, n_env)
  across_genetic_variance <- drop(crossprod(across_weights, genetic_covariance %*% across_weights))
  across_genetic_var <- rep(across_genetic_variance, length(across_pred))
  across_genetic_var[!is.finite(across_genetic_var)] <- 0
  across_rel <- reliability_thresholds(
    prediction_error_var = across_pev,
    genetic_var = across_genetic_var,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  across_int <- bayes_classify_interval_widths(
    predictions = across_pred,
    lower_bound = across_lower,
    upper_bound = across_upper,
    CI_width_thresholds = CI_width_thresholds
  )
  gid_test <- rownames(y_mat)[rowSums(is.na(y_mat)) > 0]
  across_observed <- rowMeans(y_mat, na.rm = TRUE)
  across_observed[is.nan(across_observed)] <- NA_real_
  across_env_predicted_value <- data.frame(
    name = rownames(pred_mat),
    Predicted_value = across_pred,
    Train_Test_Label = ifelse(rownames(pred_mat) %in% gid_test, "Test", "Train"),
    Observed_value = across_observed,
    Standard_error = across_se,
    PEV = across_pev,
    lower_bound = across_lower,
    upper_bound = across_upper,
    Uncertainty = across_int$Uncertainty,
    Uncertainty_remarks = across_int$reliability_remarks,
    Genetic_variance = across_genetic_var,
    Reliability = across_rel$reliability,
    Reliability_remarks = across_rel$remarks,
    Reliability_percentage = across_rel$reliability_percentage,
    Prediction_target = "Across_environment_average",
    stringsAsFactors = FALSE
  )
  names(across_env_predicted_value)[names(across_env_predicted_value) == "name"] <- gen_name

  ebv <- predicted_value[, c(gen_name, heter_groups), drop = FALSE]
  ebv$Estimated_breeding_value <- genetic_mean
  ebv$Standard_error <- predicted_value$Standard_error
  ebv$Prediction_error_variance <- predicted_value$PEV
  ebv$Genetic_variance <- predicted_value$Genetic_variance
  ebv$Reliability <- predicted_value$Reliability

  resid_var <- diag(residual_covariance)
  names(resid_var) <- env_levels
  per_env_genetic <- diag(genetic_covariance)
  names(per_env_genetic) <- env_levels
  genetic_component <- mean(per_env_genetic, na.rm = TRUE)
  residual_component <- mean(resid_var, na.rm = TRUE)
  genetic_draws <- bayes_multitrait_diagonal_draws(
    covariance_bundle$main$draws,
    length(env_levels)
  )
  residual_draws <- bayes_multitrait_diagonal_draws(
    covariance_bundle$residual$draws,
    length(env_levels)
  )
  variance_rows <- lapply(seq_along(env_levels), function(j) {
    h2_draws <- NULL
    if (!is.null(genetic_draws) && !is.null(residual_draws)) {
      n_draws <- min(nrow(genetic_draws), nrow(residual_draws))
      gd <- genetic_draws[seq_len(n_draws), j]
      rd <- residual_draws[seq_len(n_draws), j]
      h2_draws <- gd / (gd + rd)
    }
    h2 <- if (!is.null(h2_draws)) mean(h2_draws, na.rm = TRUE) else
      per_env_genetic[[j]] / (per_env_genetic[[j]] + resid_var[[j]])
    data.frame(
      Env = env_levels[[j]],
      Component = c("genetic_variance", "residual_variance", "heritability"),
      Components = c(per_env_genetic[[j]], resid_var[[j]], h2),
      Standard_error = c(
        if (is.null(genetic_draws)) NA_real_ else stats::sd(genetic_draws[, j], na.rm = TRUE),
        if (is.null(residual_draws)) NA_real_ else stats::sd(residual_draws[, j], na.rm = TRUE),
        if (is.null(h2_draws)) NA_real_ else stats::sd(h2_draws, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
  })
  across_genetic_draws <- bayes_multitrait_quadratic_draws(
    covariance_bundle$main$draws,
    across_weights
  )
  across_residual_draws <- bayes_multitrait_quadratic_draws(
    covariance_bundle$residual$draws,
    across_weights
  )
  across_residual_variance <- drop(crossprod(
    across_weights,
    residual_covariance %*% across_weights
  ))
  across_h2_draws <- NULL
  if (!is.null(across_genetic_draws) && !is.null(across_residual_draws)) {
    n_draws <- min(length(across_genetic_draws), length(across_residual_draws))
    across_h2_draws <- across_genetic_draws[seq_len(n_draws)] /
      (across_genetic_draws[seq_len(n_draws)] + across_residual_draws[seq_len(n_draws)])
  }
  variance_rows[[length(variance_rows) + 1L]] <- data.frame(
    Env = "Across_environment",
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(
      across_genetic_variance,
      across_residual_variance,
      if (is.null(across_h2_draws)) {
        across_genetic_variance / (across_genetic_variance + across_residual_variance)
      } else {
        mean(across_h2_draws, na.rm = TRUE)
      }
    ),
    Standard_error = c(
      if (is.null(across_genetic_draws)) NA_real_ else stats::sd(across_genetic_draws, na.rm = TRUE),
      if (is.null(across_residual_draws)) NA_real_ else stats::sd(across_residual_draws, na.rm = TRUE),
      if (is.null(across_h2_draws)) NA_real_ else stats::sd(across_h2_draws, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  variance_components <- do.call(rbind, variance_rows)
  rownames(variance_components) <- NULL
  residual_by_env <- data.frame(
    Environment = names(resid_var),
    Residual_variance = as.double(resid_var),
    Standard_error = if (is.null(residual_draws)) {
      NA_real_
    } else {
      apply(residual_draws, 2, stats::sd, na.rm = TRUE)
    },
    stringsAsFactors = FALSE
  )

  model_adapter <- list(
    model = list(
      y = as.vector(y_mat),
      yHat = as.vector(pred_mat),
      varE = residual_component,
      ETA = fit$ETA,
      weights = rep(1, length(y_mat)),
      multitrait_fit = fit
    ),
    output_files_names = output_files_names
  )

  result <- list(
    Coefficients = list(),
    Estimated_breeding_value = setNames(
      list(ebv),
      paste0("Multitrait_", GS_model_label, "_estimated_breeding_value")
    ),
    Total_estimated_breeding_value = ebv,
    Predicted_value = predicted_value,
    Total_Predicted_value = across_env_predicted_value,
    Residual_value = residual_value,
    Variance_components = variance_components,
    kernel_variance_components = kernel_contract$kernel_variance_components,
    kernel_configuration = kernel_contract$kernel_configuration,
    Genetic_covariance_by_kernel = kernel_contract$Genetic_covariance_by_kernel,
    genetic_covariance_by_kernel_long = kernel_contract$genetic_covariance_by_kernel_long,
    kernel_combination = kernel_contract$kernel_combination,
    Residual_variance_by_environment = residual_by_env,
    Genetic_covariance_environments = genetic_covariance,
    Genetic_correlation_environments = bayes_multitrait_safe_cov2cor(
      genetic_covariance,
      env_levels
    ),
    Genetic_covariance_environments_SE = bayes_multitrait_covariance_draw_sd(
      covariance_bundle$main$draws,
      env_levels
    ),
    Genetic_correlation_environments_SE = bayes_multitrait_correlation_draw_sd(
      covariance_bundle$main$draws,
      env_levels
    ),
    Residual_covariance_environments = residual_covariance,
    Residual_correlation_environments = bayes_multitrait_safe_cov2cor(
      residual_covariance,
      env_levels
    ),
    Residual_covariance_environments_SE = bayes_multitrait_covariance_draw_sd(
      covariance_bundle$residual$draws,
      env_levels
    ),
    Residual_correlation_environments_SE = bayes_multitrait_correlation_draw_sd(
      covariance_bundle$residual$draws,
      env_levels
    ),
    M_matrix_model_ready = kernels,
    var_u_total = genetic_component,
    diagnostic_tst_plot = NULL,
    model_notes = data.frame(
      note = c(
        paste0(
          "BGLR Multitrait (", GS_model_label,
          ") used environments as traits with ETA genetic covariance type '",
          cov_type_used,
          "' and resCov=list(type='DIAG')."
        ),
        "Genetic and residual covariance estimates are BGLR posterior means; variance and heritability SEs are posterior SDs from saved covariance draws.",
        "Environment-specific prediction SE combines BGLR marginal ETA and intercept posterior SDs in quadrature because their posterior covariance is not exposed.",
        "Across-environment SE is computed from marginal posterior SEs and assumes zero posterior covariance among environment predictions.",
        paste0(
          "MET heter_resid GBLUP_BRR and RKHS share BGLR::Multitrait's ",
          "kernel-prior sampler (model='RKHS' internally; ",
          "sigma^2*X*X' = sigma^2*K under X = eigen-sqrt(K)), so they are ",
          "numerically identical given the same RNG seed. The GS_model ",
          "label is preserved here only for traceability. Single-env fits ",
          "still use the model-specific BRR/RKHS branches in BGLR::BGLR ",
          "and diverge as expected."
        )
      ),
      stringsAsFactors = FALSE
    ),
    prediction_uncertainty_source = paste(
      "Quadrature of BGLR marginal ETA and intercept posterior SDs",
      "(their covariance is unavailable); across-environment aggregation",
      "also assumes zero posterior covariance across environments"
    )
  )
  if (identical(GS_model_label, "RKHS")) {
    result$Variance_components <-
      bayes_rkhs_standardize_heterogeneous_variance_components(
        result$Variance_components
      )
    result$model_parameters <- data.frame(
      stat = c(
        "bayesian_kernel_model",
        "stage2_observation_weighting",
        "stage2_weight_supplied",
        "stage2_weight_source",
        "stage2_weight_count",
        "stage2_weight_nonunit_count",
        "stage2_precision_min",
        "stage2_precision_max",
        "bglr_native_weight_transform",
        "rkhs_environment_mode",
        "rkhs_heterogeneous_environment_variances",
        "rkhs_genetic_environment_structure",
        "rkhs_residual_environment_structure",
        "bayes_kernel_heter_resid_effective",
        "rkhs_stage2_weight_compatibility"
      ),
      summary = c(
        "RKHS",
        "none_equal_precision",
        "FALSE",
        "none",
        as.character(nrow(pheno_data)),
        "0",
        "1",
        "1",
        "none",
        "heterogeneous_environment",
        "TRUE",
        paste0("BGLR_Multitrait_", cov_type_used),
        "BGLR_Multitrait_DIAG",
        "TRUE",
        "weights_not_supported_in_heterogeneous_environment_mode"
      ),
      stringsAsFactors = FALSE
    )
  }
  if (is.null(result$model_parameters)) {
    result$model_parameters <- data.frame(
      stat = c("bayesian_kernel_model", "rkhs_environment_mode"),
      summary = c(GS_model_label, "heterogeneous_environment"),
      stringsAsFactors = FALSE
    )
  }
  result$model_parameters <- rbind(
    result$model_parameters,
    gp_multi_kernel_parameter_rows(
      kernel_names = names(kernels),
      strategy = paste0(
        "independent_BGLR_Multitrait_", cov_type_used,
        "_kernel_ETA_covariances"
      )
    )
  )

  list(bayes_result = result, bayes_model = model_adapter, bglr_multitrait_model = fit)
}

# Decide whether a kernel-Bayesian CV fold should be served by the multitrait
# kernel-prior sampler (the path the true-prediction code uses). True when:
#   * the user asked for heter_resid Gaussian fits,
#   * the model is RKHS or BRR (the alias for GBLUP_BRR after collapse),
#   * heter_groups labels >= 2 environments in the phenotype panel, and
#   * the same genotype appears in more than one row (i.e. an actual MET).
# Returning TRUE means the CV loop must skip the univariate BGLR::BGLR fold
# path and call bayes_mod_cv_met_kernel_predict() instead, with the meta
# carried in `met_kernel_cv_meta`.
gp_bayes_should_route_met_kernel_cv <- function(GS_model,
                                                heter_resid,
                                                pheno_data,
                                                gen_name,
                                                heter_groups,
                                                response_family = "gaussian") {
  if (!isTRUE(heter_resid)) {
    return(FALSE)
  }
  if (!identical(gp_normalize_response_family(response_family), "gaussian")) {
    return(FALSE)
  }
  model <- as.character(GS_model)
  if (length(model) != 1L || !(model %in% c("RKHS", "BRR", "GBLUP_BRR"))) {
    return(FALSE)
  }
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) {
    return(FALSE)
  }
  if (!gp_bayes_has_multiple_residual_groups(pheno_data, heter_groups)) {
    return(FALSE)
  }
  if (is.null(gen_name) || !gen_name %in% names(pheno_data)) {
    return(FALSE)
  }
  anyDuplicated(as.character(pheno_data[[gen_name]])) > 0
}

# Per-fold prediction for kernel-Bayesian MET heter_resid CV. Wraps
# bayes_multitrait_env_heter_fit so each fold uses the same kernel-prior
# Multitrait sampler the true-prediction path uses (commit 35a3bcb fix
# applied uniformly to CV). `y` is the response vector with test rows
# already masked to NA (this is what the CV scheduler passes today).
# Returns predictions aligned to y[tst].
bayes_mod_cv_met_kernel_predict <- function(y, tst, met_kernel_cv_meta) {
  if (is.null(met_kernel_cv_meta)) {
    stop("bayes_mod_cv_met_kernel_predict requires met_kernel_cv_meta.", call. = FALSE)
  }
  required <- c("pheno_data", "response", "gen_name", "heter_groups",
                "kernels", "GS_model", "bayes_para")
  missing_keys <- setdiff(required, names(met_kernel_cv_meta))
  if (length(missing_keys)) {
    stop("met_kernel_cv_meta is missing entries: ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  ph <- as.data.frame(met_kernel_cv_meta$pheno_data, stringsAsFactors = FALSE)
  if (length(y) != nrow(ph)) {
    stop("bayes_mod_cv_met_kernel_predict: length(y) does not match nrow(pheno_data).",
         call. = FALSE)
  }
  ph[[met_kernel_cv_meta$response]] <- as.numeric(y)

  out <- bayes_multitrait_env_heter_fit(
    pheno_data = ph,
    response = met_kernel_cv_meta$response,
    gen_name = met_kernel_cv_meta$gen_name,
    heter_groups = met_kernel_cv_meta$heter_groups,
    kernels = met_kernel_cv_meta$kernels,
    GS_model = met_kernel_cv_meta$GS_model,
    bayes_para = met_kernel_cv_meta$bayes_para,
    omics_kernel_label = met_kernel_cv_meta$omics_kernel_label
  )
  pv <- out$bayes_result$Predicted_value
  gen_col <- met_kernel_cv_meta$gen_name
  env_col <- met_kernel_cv_meta$heter_groups
  key_ph   <- paste(as.character(ph[[gen_col]]),  as.character(ph[[env_col]]),  sep = "\r")
  key_pred <- paste(as.character(pv[[gen_col]]),  as.character(pv[[env_col]]),  sep = "\r")
  idx_in_pred <- match(key_ph[tst], key_pred)
  as.numeric(pv$Predicted_value[idx_in_pred])
}
