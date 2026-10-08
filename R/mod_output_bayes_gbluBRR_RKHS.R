#' Combine per-omics Bayesian EBV components into one per-genotype table
#'
#' Aggregates a list of per-omics estimated breeding-value (EBV) data frames
#' from BGLR fits into a single per-genotype EBV table (joined on the
#' genotype-ID column, summed across omics layers) and computes the total
#' genetic variance.
#'
#' @param ebv_list A list of per-omics EBV data frames produced by the
#'   Bayesian-finalize routines (each carrying a genotype-ID column and an
#'   `Estimated_breeding_value`-style column).
#' @param gen_name Name of the genotype-ID column shared by the EBV tables.
#' @param heter_groups Optional name of the heterogeneous-group / environment
#'   column to preserve in the join (for MET runs).
#' @param var_u_total Optional pre-computed total genetic variance to attach
#'   to the result; if `NULL`, it is summed from the omics components.
#'
#' @return A single combined EBV data frame with the per-genotype (and
#'   per-env, if `heter_groups` is supplied) summed breeding values.
#' @export
#'
#' @examples
combine_bayesian_ebv_components <- function(ebv_list, gen_name, heter_groups = NULL, var_u_total = NULL) {
  if (length(ebv_list) == 0) {
    return(NULL)
  }

  id_cols <- c(gen_name, if (!is.null(heter_groups)) heter_groups)
  base_idx <- which(vapply(ebv_list, function(x) all(id_cols %in% names(x)), logical(1)))[1]
  if (is.na(base_idx)) {
    base_idx <- 1L
  }
  base_ebv <- ebv_list[[base_idx]]
  present_id_cols <- intersect(id_cols, names(base_ebv))
  out <- base_ebv[, present_id_cols, drop = FALSE]
  out$Estimated_breeding_value <- Reduce(`+`, lapply(ebv_list, function(x) x$Estimated_breeding_value))

  if (all(vapply(ebv_list, function(x) "Prediction_error_variance" %in% names(x), logical(1)))) {
    out$Prediction_error_variance <- Reduce(`+`, lapply(ebv_list, function(x) x$Prediction_error_variance))
    out$Standard_error <- sqrt(out$Prediction_error_variance)
    out$Reliability <- NA_real_
  } else {
    out$Standard_error <- NA_real_
    out$Prediction_error_variance <- NA_real_
    out$Reliability <- NA_real_
  }

  out
}

bayes_safe_inverse <- function(mat) {
  if (length(mat) == 0 || any(dim(mat) == 0)) {
    return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
  }
  tryCatch(
    solve(mat),
    error = function(e) {
      qr.solve(mat, diag(nrow(mat)))
    }
  )
}

bayes_expected_target_variance_from_draws <- function(draws) {
  if (is.null(draws) || nrow(draws) <= 1 || ncol(draws) == 0) {
    return(0)
  }
  mean(apply(draws, 2, stats::var, na.rm = TRUE), na.rm = TRUE)
}

bayes_expected_target_variance_from_cov <- function(mean_vec, cov_mat) {
  n_target <- length(mean_vec)
  if (n_target <= 1) {
    return(0)
  }

  centering <- diag(n_target) - matrix(1 / n_target, n_target, n_target)
  as.numeric((drop(t(mean_vec) %*% centering %*% mean_vec) + sum(centering * cov_mat)) / n_target)
}

bayes_classify_interval_widths <- function(predictions,
                                           lower_bound,
                                           upper_bound,
                                           CI_width_thresholds = c(0.33, 0.66),
                                           interval_width_high_threshold = NULL,
                                           interval_width_low_threshold = NULL) {
  interval_width <- upper_bound - lower_bound

  if ((is.null(interval_width_high_threshold) && is.null(interval_width_low_threshold)) &&
      !is.null(CI_width_thresholds)) {
    high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1], na.rm = TRUE)
    low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2], na.rm = TRUE)
  } else if ((is.null(interval_width_high_threshold) || is.null(interval_width_low_threshold)) &&
             !is.null(CI_width_thresholds)) {
    high_threshold <- stats::quantile(interval_width, CI_width_thresholds[1], na.rm = TRUE)
    low_threshold <- stats::quantile(interval_width, CI_width_thresholds[2], na.rm = TRUE)
  } else {
    high_threshold <- interval_width_high_threshold
    low_threshold <- interval_width_low_threshold
  }

  rel_remarks <- ifelse(interval_width < high_threshold, "low",
                        ifelse(interval_width < low_threshold, "moderate", "high"))

  list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    Uncertainty = interval_width,
    relative_interval_width = interval_width / predictions,
    reliability_remarks = rel_remarks
  )
}

bayes_assign_target_genetic_variance <- function(groups, draw_mat = NULL, mean_vec = NULL, cov_mat = NULL) {
  groups <- as.character(groups)
  out <- rep(NA_real_, length(groups))

  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    grp_var <- if (!is.null(draw_mat)) {
      bayes_expected_target_variance_from_draws(draw_mat[idx, , drop = FALSE])
    } else {
      bayes_expected_target_variance_from_cov(
        mean_vec = mean_vec[idx],
        cov_mat = cov_mat[idx, idx, drop = FALSE]
      )
    }
    out[idx] <- grp_var
  }

  out
}

bayes_build_target_summary <- function(predicted_mean,
                                       genetic_mean,
                                       variance_groups,
                                       draw_mat = NULL,
                                       cov_mat = NULL,
                                       prior_genetic_variance = NULL,
                                       confidence_level = 0.95,
                                       CI_width_thresholds = c(0.33, 0.66),
                                       interval_width_high_threshold = NULL,
                                       interval_width_low_threshold = NULL,
                                       high_reliability_thres = 0.9,
                                       low_reliability_thres = 0.5) {
  if (!is.null(draw_mat)) {
    standard_error <- apply(draw_mat, 1, stats::sd, na.rm = TRUE)
    prediction_error_var <- apply(draw_mat, 1, stats::var, na.rm = TRUE)
    alpha <- (1 - confidence_level) / 2
    lower_bound <- apply(draw_mat, 1, stats::quantile, probs = alpha, na.rm = TRUE)
    upper_bound <- apply(draw_mat, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    interval_stats <- bayes_classify_interval_widths(
      predictions = predicted_mean,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      CI_width_thresholds = CI_width_thresholds,
      interval_width_high_threshold = interval_width_high_threshold,
      interval_width_low_threshold = interval_width_low_threshold
    )
    if (!is.null(prior_genetic_variance)) {
      genetic_var <- as.double(prior_genetic_variance)
    } else {
      genetic_var <- bayes_assign_target_genetic_variance(groups = variance_groups, draw_mat = draw_mat)
    }
  } else {
    prediction_error_var <- pmax(0, diag(cov_mat))
    standard_error <- sqrt(prediction_error_var)
    interval_stats <- reliability_thresholds_MPIW_from_CI(
      predictions = predicted_mean,
      standard_errors = standard_error,
      confidence_level = confidence_level,
      CI_width_thresholds = CI_width_thresholds,
      interval_width_high_threshold = interval_width_high_threshold,
      interval_width_low_threshold = interval_width_low_threshold,
      model_for_CI_cal = "RKHS",
      boot_results = NULL
    )
    if (!is.null(prior_genetic_variance)) {
      genetic_var <- as.double(prior_genetic_variance)
    } else {
      genetic_var <- bayes_assign_target_genetic_variance(
        groups = variance_groups,
        mean_vec = genetic_mean,
        cov_mat = cov_mat
      )
    }
  }

  rel <- reliability_thresholds(
    prediction_error_var = prediction_error_var,
    genetic_var = genetic_var,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )

  list(
    predicted_mean = predicted_mean,
    genetic_mean = genetic_mean,
    Standard_error = standard_error,
    PEV = prediction_error_var,
    lower_bound = interval_stats$lower_bound,
    upper_bound = interval_stats$upper_bound,
    Uncertainty = interval_stats$Uncertainty,
    Uncertainty_remarks = interval_stats$reliability_remarks,
    Genetic_variance = genetic_var,
    Reliability = rel$reliability,
    Reliability_remarks = rel$remarks,
    Reliability_percentage = rel$reliability_percentage
  )
}

bayes_build_aggregation_matrix <- function(ids) {
  ids <- as.character(ids)
  unique_ids <- unique(ids)
  A <- matrix(0, nrow = length(unique_ids), ncol = length(ids),
              dimnames = list(unique_ids, NULL))
  for (i in seq_along(unique_ids)) {
    idx <- which(ids == unique_ids[i])
    A[i, idx] <- 1 / length(idx)
  }
  A
}

bayes_environment_covariance_basis <- function(basis, ids, environments) {
  basis <- as.matrix(basis)
  ids <- as.character(ids)
  environments <- as.character(environments)
  env_levels <- unique(environments)
  ids_by_env <- lapply(env_levels, function(env) unique(ids[environments == env]))
  common_ids <- Reduce(intersect, ids_by_env)
  if (!length(common_ids)) {
    return(NULL)
  }
  out <- matrix(NA_real_, length(env_levels), length(env_levels),
                dimnames = list(env_levels, env_levels))
  for (i in seq_along(env_levels)) {
    for (j in seq_len(i)) {
      values <- vapply(common_ids, function(gid) {
        row_i <- which(ids == gid & environments == env_levels[[i]])
        row_j <- which(ids == gid & environments == env_levels[[j]])
        mean(basis[row_i, row_j, drop = FALSE], na.rm = TRUE)
      }, numeric(1L))
      out[i, j] <- out[j, i] <- mean(values, na.rm = TRUE)
    }
  }
  out
}

bayes_eta_variance_draws <- function(output_files_names,
                                     eta_index,
                                     eta_model,
                                     draw_indices = NULL) {
  gp_bayes_read_eta_variance_draws(
    output_files_names = output_files_names,
    eta_index = eta_index,
    eta_model = eta_model,
    draw_indices = draw_indices
  )
}

bayes_environment_covariance_summary <- function(mod,
                                                  ETA,
                                                  pheno_data,
                                                  gen_name,
                                                  heter_groups,
                                                  draw_indices = NULL) {
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) {
    return(NULL)
  }
  eta_input <- ETA[["ETA"]]
  eta_fit <- mod$model$ETA
  eta_models <- vapply(eta_input, function(x) x$model, character(1L))
  random_idx <- which(eta_models != "FIXED")
  covariance <- NULL
  term_draws <- list()
  term_basis <- list()
  for (eta_index in random_idx) {
    term <- eta_input[[eta_index]]
    basis <- if (!is.null(term$K)) {
      as.matrix(term$K)
    } else if (!is.null(term$X)) {
      tcrossprod(as.matrix(term$X))
    } else {
      NULL
    }
    if (is.null(basis)) {
      next
    }
    env_basis <- bayes_environment_covariance_basis(
      basis,
      ids = pheno_data[[gen_name]],
      environments = pheno_data[[heter_groups]]
    )
    if (is.null(env_basis)) {
      next
    }
    variance <- eta_fit[[eta_index]]$varU %||% eta_fit[[eta_index]]$varB %||% NA_real_
    if (is.finite(variance)) {
      covariance <- if (is.null(covariance)) {
        variance * env_basis
      } else {
        covariance + variance * env_basis
      }
    }
    draws <- bayes_eta_variance_draws(
      output_files_names = mod[["output_files_names"]],
      eta_index = eta_index,
      eta_model = eta_models[[eta_index]],
      draw_indices = draw_indices
    )
    if (!is.null(draws)) {
      term_draws[[as.character(eta_index)]] <- draws
      term_basis[[as.character(eta_index)]] <- env_basis
    }
  }
  if (is.null(covariance)) {
    return(NULL)
  }
  covariance <- (covariance + t(covariance)) / 2
  labels <- rownames(covariance)
  lower_draws <- NULL
  if (length(term_draws)) {
    n_draws <- min(vapply(term_draws, length, integer(1L)))
    lower_index <- bayes_multitrait_lower_triangle_index(length(labels))
    lower_draws <- matrix(0, n_draws, nrow(lower_index))
    for (nm in names(term_draws)) {
      basis_values <- term_basis[[nm]][lower_index]
      lower_draws <- lower_draws + outer(
        term_draws[[nm]][seq_len(n_draws)],
        basis_values
      )
    }
  }
  list(
    covariance = covariance,
    correlation = bayes_multitrait_safe_cov2cor(covariance, labels),
    covariance_se = bayes_multitrait_covariance_draw_sd(lower_draws, labels),
    correlation_se = bayes_multitrait_correlation_draw_sd(lower_draws, labels)
  )
}

bayes_compute_prediction_targets <- function(mod,
                                             ETA,
                                             GS_model,
                                             pheno_data,
                                             gen_name,
                                             heter_groups = NULL,
                                             output_files_names = NULL,
                                             draw_indices = NULL,
                                             confidence_level = 0.95,
                                             CI_width_thresholds = c(0.33, 0.66),
                                             interval_width_high_threshold = NULL,
                                             interval_width_low_threshold = NULL,
                                             high_reliability_thres = 0.9,
                                             low_reliability_thres = 0.5) {
  env_col <- heter_groups
  if (is.null(env_col) && "Env" %in% names(pheno_data)) {
    env_col <- "Env"
  }
  eta_input <- ETA[["ETA"]]
  eta_model <- mod$model$ETA
  eta_types <- vapply(eta_input, function(x) x$model, character(1))
  random_idx <- which(eta_types != "FIXED")
  fixed_idx <- which(eta_types == "FIXED")
  random_types <- eta_types[random_idx]
  all_random_brr <- length(random_idx) > 0 && all(random_types == "BRR")

  if (length(random_idx) == 0) {
    return(NULL)
  }

  n_obs <- nrow(pheno_data)
  observed_trait_rows <- !is.na(mod$model$y)
  if (length(observed_trait_rows) != n_obs) {
    observed_trait_rows <- rep(TRUE, n_obs)
  }
  has_met <- anyDuplicated(as.character(pheno_data[[gen_name]][observed_trait_rows])) > 0
  prediction_mean_obs <- as.double(mod$model$yHat)
  obs_groups <- if (is.null(env_col)) rep("Overall", n_obs) else as.character(pheno_data[[env_col]])

  if (all_random_brr) {
    random_draws <- vector("list", length(random_idx))
    genetic_mean_obs <- rep(0, n_obs)
    prior_cov_obs <- matrix(0, n_obs, n_obs)

    for (k in seq_along(random_idx)) {
      eta_idx <- random_idx[k]
      b_file <- output_files_names[grepl(paste0("ETA_", eta_idx, "_b.bin"), output_files_names)]
      if (length(b_file) != 1) {
        stop(paste("Missing BRR posterior draw file for ETA index", eta_idx), call. = FALSE)
      }
      beta_draws <- BGLR::readBinMat(b_file)
      if (length(draw_indices) && max(draw_indices) <= nrow(beta_draws)) {
        beta_draws <- beta_draws[draw_indices, , drop = FALSE]
      }
      X_eta <- as.matrix(eta_input[[eta_idx]]$X)
      eta_draws <- X_eta %*% t(beta_draws)
      random_draws[[k]] <- eta_draws
      genetic_mean_obs <- genetic_mean_obs + drop(X_eta %*% eta_model[[eta_idx]]$b)
      sigma_b <- as.double(eta_model[[eta_idx]]$varB)
      if (!is.finite(sigma_b)) {
        sigma_b <- 0
      }
      prior_cov_obs <- prior_cov_obs + sigma_b * tcrossprod(X_eta)
    }

    total_draws <- Reduce(`+`, random_draws)
    prediction_draws <- total_draws
    if (length(fixed_idx)) {
      for (eta_idx in fixed_idx) {
        b_file <- output_files_names[grepl(
          paste0("ETA_", eta_idx, "_b\\.bin$"),
          basename(output_files_names)
        )]
        if (length(b_file) == 1L) {
          beta_draws <- as.matrix(BGLR::readBinMat(b_file))
          if (length(draw_indices) && max(draw_indices) <= nrow(beta_draws)) {
            beta_draws <- beta_draws[draw_indices, , drop = FALSE]
          }
          fixed_draws <- as.matrix(eta_input[[eta_idx]]$X) %*% t(beta_draws)
          n_draws <- min(ncol(prediction_draws), ncol(fixed_draws))
          prediction_draws <- prediction_draws[, seq_len(n_draws), drop = FALSE] +
            fixed_draws[, seq_len(n_draws), drop = FALSE]
        }
      }
    }
    mu_draws <- gp_bayes_read_mu_draws(
      output_files_names,
      draw_indices = draw_indices,
      n_draws = ncol(prediction_draws)
    )
    if (is.null(mu_draws)) {
      mu_draws <- rep(as.double(mod$model$mu), ncol(prediction_draws))
    }
    prediction_draws <- sweep(prediction_draws, 2, mu_draws, FUN = "+")
    # Keep the public point estimate, SE/PEV, and interval on one posterior
    # sample contract. BGLR's yHat is a running mean over all post-burn-in
    # iterations, whereas saveEffects records the thinned draws used below.
    # With short or strongly thinned chains those Monte Carlo summaries can
    # differ enough that yHat falls outside the interval computed from the
    # saved draws. The mean of the same retained target draws is the
    # reproducible point estimate corresponding to the reported uncertainty.
    prediction_mean_obs <- rowMeans(prediction_draws, na.rm = TRUE)
    obs_summary <- bayes_build_target_summary(
      predicted_mean = prediction_mean_obs,
      genetic_mean = genetic_mean_obs,
      variance_groups = obs_groups,
      draw_mat = prediction_draws,
      prior_genetic_variance = pmax(0, diag(prior_cov_obs)),
      confidence_level = confidence_level,
      CI_width_thresholds = CI_width_thresholds,
      interval_width_high_threshold = interval_width_high_threshold,
      interval_width_low_threshold = interval_width_low_threshold,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )

    if (has_met) {
      A <- bayes_build_aggregation_matrix(pheno_data[[gen_name]])
      across_draws <- A %*% prediction_draws
      prior_cov_across <- A %*% prior_cov_obs %*% t(A)
      prior_cov_across <- (prior_cov_across + t(prior_cov_across)) / 2
      across_summary <- bayes_build_target_summary(
        predicted_mean = drop(A %*% prediction_mean_obs),
        genetic_mean = drop(A %*% genetic_mean_obs),
        variance_groups = rep("Across_environment", nrow(A)),
        draw_mat = across_draws,
        prior_genetic_variance = pmax(0, diag(prior_cov_across)),
        confidence_level = confidence_level,
        CI_width_thresholds = CI_width_thresholds,
        interval_width_high_threshold = interval_width_high_threshold,
        interval_width_low_threshold = interval_width_low_threshold,
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres
      )
      across_summary[[gen_name]] <- rownames(A)
    } else {
      across_summary <- NULL
    }
  } else {
    total_G <- matrix(0, n_obs, n_obs)
    genetic_mean_obs <- rep(0, n_obs)

    for (eta_idx in random_idx) {
      eta_term <- eta_input[[eta_idx]]
      eta_fit <- eta_model[[eta_idx]]
      eta_type <- eta_types[[eta_idx]]

      if (!is.null(eta_term$K) && !is.null(eta_fit$varU) && !is.null(eta_fit$u)) {
        K_eta <- as.matrix(eta_term$K)
        sigma_u <- as.double(eta_fit$varU)
        if (!is.finite(sigma_u)) {
          sigma_u <- 0
        }
        total_G <- total_G + sigma_u * K_eta
        genetic_mean_obs <- genetic_mean_obs + as.double(eta_fit$u)
      } else if (identical(eta_type, "BRR") && !is.null(eta_term$X) && !is.null(eta_fit$varB) && !is.null(eta_fit$b)) {
        X_eta <- as.matrix(eta_term$X)
        sigma_b <- as.double(eta_fit$varB)
        if (!is.finite(sigma_b)) {
          sigma_b <- 0
        }
        total_G <- total_G + sigma_b * tcrossprod(X_eta)
        genetic_mean_obs <- genetic_mean_obs + drop(X_eta %*% eta_fit$b)
      } else {
        stop("Unsupported Bayesian ETA term encountered while computing gaussian target metrics.", call. = FALSE)
      }
    }

    weights <- mod$model$weights
    if (is.null(weights)) {
      weights <- rep(1, n_obs)
    }
    residual_diag <- gp_bayes_observation_residual_variance(
      varE = mod$model$varE,
      groups = if (is.null(env_col)) NULL else pheno_data[[env_col]],
      weights = weights,
      n = n_obs
    )
    V <- total_G + diag(residual_diag, nrow = n_obs)
    V_inv <- bayes_safe_inverse(V)

    X_fixed <- matrix(1, nrow = n_obs, ncol = 1L)
    colnames(X_fixed) <- "(Intercept)"
    if (length(fixed_idx) > 0) {
      fixed_eta <- do.call(cbind, lapply(fixed_idx, function(idx) as.matrix(eta_input[[idx]]$X)))
      if (!is.null(fixed_eta) && ncol(fixed_eta) > 0) {
        X_fixed <- cbind(X_fixed, fixed_eta)
      }
    }

    posterior_cov_obs <- total_G - total_G %*% V_inv %*% total_G
    if (!is.null(X_fixed)) {
      XtVinvX <- crossprod(X_fixed, V_inv %*% X_fixed)
      XtVinvX_inv <- bayes_safe_inverse(XtVinvX)
      posterior_cov_obs <- posterior_cov_obs +
        total_G %*% V_inv %*% X_fixed %*% XtVinvX_inv %*% crossprod(X_fixed, V_inv %*% total_G)
    }
    posterior_cov_obs <- (posterior_cov_obs + t(posterior_cov_obs)) / 2

    obs_summary <- bayes_build_target_summary(
      predicted_mean = prediction_mean_obs,
      genetic_mean = genetic_mean_obs,
      variance_groups = obs_groups,
      cov_mat = posterior_cov_obs,
      prior_genetic_variance = pmax(0, diag(total_G)),
      confidence_level = confidence_level,
      CI_width_thresholds = CI_width_thresholds,
      interval_width_high_threshold = interval_width_high_threshold,
      interval_width_low_threshold = interval_width_low_threshold,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )

    if (has_met) {
      A <- bayes_build_aggregation_matrix(pheno_data[[gen_name]])
      across_cov <- A %*% posterior_cov_obs %*% t(A)
      across_cov <- (across_cov + t(across_cov)) / 2
      prior_cov_across <- A %*% total_G %*% t(A)
      prior_cov_across <- (prior_cov_across + t(prior_cov_across)) / 2
      across_summary <- bayes_build_target_summary(
        predicted_mean = drop(A %*% prediction_mean_obs),
        genetic_mean = drop(A %*% genetic_mean_obs),
        variance_groups = rep("Across_environment", nrow(A)),
        cov_mat = across_cov,
        prior_genetic_variance = pmax(0, diag(prior_cov_across)),
        confidence_level = confidence_level,
        CI_width_thresholds = CI_width_thresholds,
        interval_width_high_threshold = interval_width_high_threshold,
        interval_width_low_threshold = interval_width_low_threshold,
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres
      )
      across_summary[[gen_name]] <- rownames(A)
    } else {
      across_summary <- NULL
    }
  }

  list(observation = obs_summary, across = across_summary)
}

#' Format Bayesian GBLUP, BRR, and RKHS model output
#'
#' Builds prediction, uncertainty, reliability, and diagnostic summaries from
#' Bayesian relationship/kernel model fits.
#'
#' @param mod Fitted Bayesian model object.
#' @param ETA Bayesian ETA list used for the model.
#' @param gen_name Genotype identifier column name.
#' @param GS_model Model name.
#' @param gkernel Optional genomic kernel matrix.
#' @param gmatrix Optional genomic relationship matrix.
#' @param omic1_kernel Optional first omic kernel matrix.
#' @param omic2_kernel Optional second omic kernel matrix.
#' @param omic3_kernel Optional third omic kernel matrix.
#' @param kernel_list Optional list of kernels.
#' @param omics_kernel_label Optional omic kernel labels.
#' @param pheno_data Phenotype data frame used by the model.
#' @param response Response column name used to preserve original ordered class
#'   levels for Bayesian classification output.
#' @param heter_groups Optional environment/grouping column name.
#' @param bayes_para Bayesian parameter list.
#' @param CI_width_thresholds Quantile thresholds for confidence interval width
#'   categories.
#' @param confidence_level Confidence level for prediction intervals.
#' @param high_reliability_thres High reliability threshold.
#' @param low_reliability_thres Low reliability threshold.
#' @param n_components Number of components used in diagnostic summaries.
#' @param threshold Numeric threshold used by diagnostic helpers.
#' @param target Prediction target label.
#' @param iqr_multiplier IQR multiplier for outlier diagnostics.
#' @param interval_width_high_threshold Optional absolute high interval-width
#'   threshold.
#' @param interval_width_low_threshold Optional absolute low interval-width
#'   threshold.
#' @param interval_width_moderate_threshold Optional absolute moderate
#'   interval-width threshold.
#' @param system_database Logical indicating whether database-mode output is
#'   active.
#' @param response_family Response family for the model output.
#' @param ... Additional arguments reserved for compatibility.
#'
#' @return A list of formatted Bayesian model output tables and diagnostics.
#' @export
mod_output_bayes_gbluBRR_RKHS <- function(mod=NULL,
                                         ETA=NULL,
                                         gen_name=NULL,
                                         GS_model = NULL,
                                         gkernel=NULL,
                                         gmatrix = NULL,
                                         omic1_kernel=NULL,
                                         omic2_kernel=NULL,
                                         omic3_kernel=NULL,
                                         kernel_list = NULL,
                                         omics_kernel_label = list(omic1_kernel = NULL,
                                                                   omic2_kernel = NULL,
                                                                   omic3_kernel = NULL),
                                          pheno_data = NULL,
                                          heter_groups = NULL,
                                         bayes_para = NULL,
                                         CI_width_thresholds = c(0.33, 0.66),
                                         confidence_level = 0.95,
                                         high_reliability_thres = 0.9,
                                         low_reliability_thres = 0.5,
                                         n_components = 20,
                                         threshold = 100,
                                         target = "test_set",
                                         iqr_multiplier = 1.5,
                                         interval_width_high_threshold = NULL,
                                         interval_width_low_threshold = NULL,
                                         interval_width_moderate_threshold = NULL,
                                         system_database = FALSE,
                                          response_family = "gaussian",
                                          response = NULL,
                                          ...){

  fam <- gp_resolve_response_family(response_family, y = mod$model$y)
  if (!identical(fam, "gaussian")) {
    res <- gp_bayes_classification_output(
      mod = mod,
      pheno_data = pheno_data,
      gen_name = gen_name,
      response = response,
      response_family = fam,
      heter_groups = heter_groups,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      ETA = ETA,
      bayes_para = bayes_para,
      GS_model = GS_model,
      confidence_level = confidence_level
    )
    unlink(mod[["output_files_names"]])
    return(res)
  }

  diagnostic_plots <- NULL
  tst <- which(is.na(mod$model$y))

  if (is.null(heter_groups) && "Env" %in% names(pheno_data)) {
    heter_groups <- "Env"
  }

  if(length(tst)==0) tst <- NULL
  train_test_label <-  ifelse(is.na(mod$model$y), "Test", "Train")
  GID_unique_all <- unique(as.character(pheno_data[, gen_name, drop = TRUE]))
  GID_tst <- unique(as.character(pheno_data[tst, gen_name, drop = TRUE]))
  train_test_label_across_env <- ifelse(GID_unique_all %in% GID_tst, "Test", "Train")

  # if(length(tst)>1){
  # Standard_error <- mod$model$SD.yHat
  # PEV <- (mod$model$SD.yHat)^2
  # Reliability <- 1 - (PEV/var(mod$model$yHat))
  #
  # } else {
  #   Standard_error <- mod$model$SD.yHat
  #   PEV <- (mod$model$SD.yHat)^2
  #   Reliability <- 1 - (PEV/var(mod$model$yHat))
  #
  # }

  msg <- ""
  Zg <- NULL
  ### Check if the user provide lable/name for the omics data

  if(inherits(omics_kernel_label,'list')){
    if(!all(sapply(omics_kernel_label, function(x){ is.null(x)}))!=FALSE){
      label <-  which(sapply(omics_kernel_label, function(x) !is.null(x)))
      print_lable <-  omics_kernel_label[label]
    } else {
      print_lable <-  NULL
    }

  }else {
    if(inherits(omics_kernel_label,"character")){
      print_lable <-  omics_kernel_label
    }
    if(is.null(omics_kernel_label)){
      print_lable <-  NULL
    }
  }


  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){

    ### Residual value is only estimable for response value without NA

    ### incidence matrix for main eff. of the genotypes

    Zg <- TRUE  # flag: record-level (genotype x environment) output rows

    ### Extract all environments in MET
    all_envs_for_met <-  as.character(pheno_data[,heter_groups])

  } else{
    Zg <- NULL
    all_envs_for_met <- NULL
  }

  datasets <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    gkernel = gkernel,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  dataset_names <- names(datasets)
  eta_models <- vapply(ETA[["ETA"]], function(x) x$model, character(1))
  random_eta_idx <- which(eta_models != "FIXED")
  if (length(random_eta_idx) < length(datasets)) {
    stop(paste(msg, "Unable to align Bayesian ETA terms with the supplied kernels."), call. = FALSE)
  }
  kernel_eta_idx <- random_eta_idx[seq_along(datasets)]

  train_test_label <-  ifelse(is.na(mod$model$y), "Test", "Train")
  tst <- which(is.na(mod$model$y))
  tst_GID <- unique(as.character(pheno_data[tst, gen_name]))

  # Subset datasets if tst has more than 1 element
  data_trn <- NULL
  data_tst <- NULL

  if (length(tst) > 1) {
    if(!is.null(heter_groups)){
      data_trn <- lapply(datasets, function(dataset) dataset[!rownames(dataset)%in%tst_GID, !colnames(dataset)%in%tst_GID])
      data_tst <- lapply(datasets, function(dataset) dataset[rownames(dataset)%in%tst_GID, colnames(dataset)%in%tst_GID])
    }else{
    data_trn <- lapply(datasets, function(dataset) dataset[-tst, -tst])
    data_tst <- lapply(datasets, function(dataset) dataset[tst, tst])

    }
  }

  if(length(datasets)>1 & length(tst)>1){
    data_trn <- do.call(cbind, data_trn)
    data_trn <- scale(data_trn)
    #
    data_tst <- do.call(cbind, data_tst)
    data_tst <- scale(data_tst)

    dataset <- do.call(cbind, datasets)
  } else{
    dataset <- do.call(cbind, datasets)

  }

  ### Check bayes_parameter_check function in bayesians_preprocess for details
  nIter <- bayes_para[["nIter"]]
  burnIn <- bayes_para[["burnIn"]]

  posindex <- gp_bayes_saved_draw_indices(nIter = nIter, burnIn = burnIn, thin = bayes_para[["thin"]])

  #######################################
  eta_models_present <- vapply(ETA[["ETA"]], function(x) x$model, character(1))
  BIN <- mod[["output_files_names"]][grepl("ETA_[0-9]+_b\\.bin$", basename(mod[["output_files_names"]]))]
  ### Extract Error variance
  var_residual <- gp_bayes_read_varE_draws(mod[["output_files_names"]], draw_indices = posindex)
  #var_residual <- var_residual[posindex]
  var_residual <- var_residual

  se_var_residual <- standard_deviation(var_residual)
  #########
  varB_files <- mod[["output_files_names"]][grepl("ETA_[0-9]+_varB.dat$", basename(mod[["output_files_names"]]))]
  varU_files <- mod$output_files_names[grepl("_varU.dat", mod[["output_files_names"]])]

  #########

  #Var_U_omics_list <- vector(mode = "list", length = length(BIN))
  var_u_mean_omics_list <- list()
  se_var_u_omics_list <-  list()
  var_u_total <- 0
  genomic_h2 <- 0
  se_genomic_h2 <- 0
  Predicted_value_for_CI_total <- 0
  datasets <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    gkernel = gkernel,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  dataset_names <- names(datasets)
  posterior_list <-   list()
  res_coeff_ebv_pev_rel_se_list <-  list()
  coefficients_list <-   list()
  estimated_breeding_value_list <-  list()
  m_matrix_model_ready_list <-   list()
  sum_posterior <-  0
  sum_estimated_breeding_value <-  0

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    eta_idx <- kernel_eta_idx[i]
    term_model <- ETA[["ETA"]][[eta_idx]]$model
    coef_file <- if (term_model == "BRR") {
      BIN[grepl(paste0("ETA_", eta_idx, "_b\\.bin$"), basename(BIN))]
    } else {
      NULL
    }
    varb_file <- if (term_model == "RKHS") {
      varU_files[which(kernel_eta_idx == eta_idx)[1]]
    } else {
      NULL
    }

    if (!is.null(dataset)) {
      if (term_model == "BRR" && length(coef_file) != 1) {
        stop(paste(msg, "Unable to find the BRR posterior coefficient file for", dataset_names[i]), call. = FALSE)
      }
      if (term_model == "RKHS" && length(varb_file) != 1) {
        stop(paste(msg, "Unable to find the RKHS variance file for", dataset_names[i]), call. = FALSE)
      }
        #Var_U_omics_list[[dataset_names[i]]] <- process_varU(varB_files[i], posindex, GS_model)
        #gid_name <- rownames(geno_data)
        #var_u_omics <- process_var_u(varB_files[i], posindex, GS_model)
        if(term_model %in% c("BRR")){
          eta_matrix <- as.matrix(ETA[["ETA"]][[eta_idx]]$X)
          beta_draws <- BGLR::readBinMat(coef_file)
          if (length(posindex) && max(posindex) <= nrow(beta_draws)) {
            beta_draws <- beta_draws[posindex, , drop = FALSE]
          }
          Predicted_value_for_CI <- eta_matrix %*% t(beta_draws)
          train_idx <- if (length(tst) > 0) setdiff(seq_len(nrow(eta_matrix)), tst) else seq_len(nrow(eta_matrix))
          train_u <- Predicted_value_for_CI[train_idx, , drop = FALSE]
          train_y <- mod$model$y[train_idx]

          var_u_omics <- apply(train_u, 2, stats::var, na.rm = TRUE)
          # Phase 3.24: do NOT overwrite var_residual here. The BGLR varE
          # chain read at line ~630 is the correct residual variance.
          # var(u_main - y) is not the residual -- it lumps GxE + intercept
          # + true residual together and inflates Ve in the final
          # variance_components table.

          coeff_mean <- colMeans(beta_draws)
          coeff_labels <- colnames(eta_matrix)
          if (is.null(coeff_labels) || length(coeff_labels) != length(coeff_mean)) {
            coeff_labels <- seq_along(coeff_mean)
          }
          coefficient_df <- data.frame(
            X_variables = coeff_labels,
            Coeff = coeff_mean,
            stringsAsFactors = FALSE
          )

          ebv_mean <- drop(eta_matrix %*% mod$model$ETA[[eta_idx]]$b)
          ebv_se <- apply(Predicted_value_for_CI, 1, stats::sd, na.rm = TRUE)
          ebv_pev <- apply(Predicted_value_for_CI, 1, stats::var, na.rm = TRUE)
          ebv_rel <- pmax(0, pmin(1, 1 - (ebv_pev / mean(var_u_omics, na.rm = TRUE))))

          if (is.null(Zg)) {
            ebv_df <- data.frame(
              name = rownames(dataset),
              Estimated_breeding_value = ebv_mean,
              Standard_error = ebv_se,
              Prediction_error_variance = ebv_pev,
              Reliability = ebv_rel,
              stringsAsFactors = FALSE
            )
            names(ebv_df)[1] <- gen_name
          } else {
            ebv_df <- data.frame(
              name = pheno_data[, gen_name],
              Env = all_envs_for_met,
              Estimated_breeding_value = ebv_mean,
              Standard_error = ebv_se,
              Prediction_error_variance = ebv_pev,
              Reliability = ebv_rel,
              stringsAsFactors = FALSE
            )
            names(ebv_df)[1:2] <- c(gen_name, heter_groups)
          }

          res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]] <- list(
            Posterior = Predicted_value_for_CI,
            Coefficient = coefficient_df,
            Estimated_breeding_value = ebv_df,
            PEV = ebv_pev,
            Reliability = ebv_rel,
            Standard_error = ebv_se,
            Predicted_value_for_CI = Predicted_value_for_CI
          )
        } else{
          if(term_model %in% c("RKHS")){
          var_u_omics <- process_var_u(varb_file, posindex, term_model)
          eta_kernel <- as.matrix(ETA[["ETA"]][[eta_idx]]$K)
          train_idx <- if (length(tst) > 0) setdiff(seq_len(nrow(eta_kernel)), tst) else seq_len(nrow(eta_kernel))
          kernel_scale <- mean(diag(eta_kernel)[train_idx], na.rm = TRUE)
          if (is.finite(kernel_scale) && kernel_scale > 0) {
            var_u_omics <- var_u_omics * kernel_scale
          }
          }
        }
        #posterior_list[[dataset_names[i]]] <- BGLR::readBinMat(BIN[aa])
        if(term_model != "BRR") {
          if(term_model=="RKHS"){
        res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]] <-   cal_coeff_ebv_pev_rel_RHKS_glub(
                                                        mod =  mod,
                                                        gmatrix = dataset,  # not used by the helper (reads mod$model$ETA)
                                                        gen_name = gen_name,
                                                        var_u = mean(var_u_omics),
                                                        var_E = mean(var_residual),
                                                        hetero = all_envs_for_met,
                                                        heter_groups = heter_groups,
                                                        gid_name = if (!is.null(Zg)) as.character(pheno_data[, gen_name]) else rownames(dataset),
                                                        eta_index = eta_idx,
                                                        )


          }
        }
        var_u_mean_omics_list[[dataset_names[i]]] <- mean(var_u_omics)
        se_var_u_omics_list[[dataset_names[i]]] <- standard_deviation(var_u_omics)
        coefficients_list[[paste("coefficient",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Coefficient"]]
        estimated_breeding_value_list[[paste("estimated_breeding_value",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]]
        if(term_model=="BRR"){
        sum_posterior <- sum_posterior + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Posterior"]]
        }
        sum_estimated_breeding_value <- sum_estimated_breeding_value + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]][,"Estimated_breeding_value"]
        if(dataset_names[i]=="gmatrix"){
        m_matrix_model_ready_list[[paste(gsub("gmatrix", "geno", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
        } else{
          m_matrix_model_ready_list[[paste(gsub("_kernel", "", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
        }
        var_u_total <- var_u_total + var_u_omics
        if(term_model=="BRR"){
        Predicted_value_for_CI_total <- Predicted_value_for_CI_total  + Predicted_value_for_CI
}
        if(length(datasets)==1){
          # predicted_value[, 1] <- rownames(dataset)
          # #PEV <- apply(g_ebv, 1, var)
          # predicted_value <- predicted_value |>
          #   dplyr::mutate(
          #                 #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
          #                 Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
          #                 Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
          #                 Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])

          if(is.null(Zg)){
          sum_ebv <- data.frame(name = rownames(dataset),
                                Estimated_breeding_value = sum_estimated_breeding_value,
                                stringsAsFactors = FALSE)

          colnames(sum_ebv)[1] <- gen_name

          } else{
            sum_ebv <- data.frame(name = as.character(pheno_data[, gen_name]),  # one row per record
                                  Env = all_envs_for_met,
                                  Estimated_breeding_value = sum_estimated_breeding_value,
                                  stringsAsFactors = FALSE)

            colnames(sum_ebv)[1:2] <- c(gen_name, heter_groups)
          }
          sum_ebv <- sum_ebv |>
            dplyr::mutate(
                          #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
                          Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
                          Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
                          Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])


          if(length(tst)>1){
            #residual_value[, 1] <- rownames(dataset)[tst]
            sum_ebv <- sum_ebv[tst, ]


          }
          # else {
          #
          #   residual_value[, 1] <- rownames(dataset)
          # }
          #### MET
          # if(!is.null(Zg)){
          #   sep_pev_rel <- sep_pev_rel_gblup(geno_object = dataset,
          #                                    va = mean(var_u_omics),
          #                                    ve = mean(var_residual))
          #
          #   across_env_predicted_value <- predicted_value |>
          #     dplyr::ungroup() |>
          #     dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
          #     dplyr::summarise(
          #       Predicted_value = mean(Predicted_value, na.rm = TRUE),
          #       Standard_error  = mean(Standard_error,  na.rm = TRUE),
          #       PEV             = mean(PEV,             na.rm = TRUE),
          #       Reliability     = mean(Reliability,     na.rm = TRUE),
          #       .groups         = "drop"
          #     )
          #
          #   across_env_predicted_value <- na.omit(across_env_predicted_value)
          #
          #
          # }
          ### MET End
        } else{

          if(length(datasets)>1){
            if(i==1) gid_name <- rownames(dataset)
            #predicted_value[, 1] <- gid_name
            #residual_value[, 1] <- gid_name
            ##### Treat sum_EBV
            if(i==length(datasets)){
            if(GS_model=="BRR"){
            pev <- apply(sum_posterior, 1, var)
            rel <- 1 - (pev / mean(var_u_total))
            #rel <- ifelse(rel<0, NA, rel)

            } else {
              pev <- NA
              rel <- NA


            }

            if(is.null(Zg)){
              sum_ebv <- data.frame(name = gid_name,
                                    Estimated_breeding_value = sum_estimated_breeding_value,
                                    stringsAsFactors = FALSE)

              colnames(sum_ebv)[1] <- gen_name

            } else{
              sum_ebv <- data.frame(name = as.character(pheno_data[, gen_name]),  # one row per record
                                    Env = all_envs_for_met,
                                    Estimated_breeding_value = sum_estimated_breeding_value,
                                    stringsAsFactors = FALSE)

              colnames(sum_ebv)[1:2] <- c(gen_name, heter_groups)
            }

              #if(length(tst)>1) sum_ebv <- sum_ebv[tst, ]
            sum_ebv <- sum_ebv |>
              dplyr::mutate(
                           #Standard_error = ifelse(!is.na(pev), sqrt(pev), NA),
                            Standard_error = NA,
                            Prediction_error_variance =NA,
                            Reliability = NA)

            # predicted_value <- predicted_value |>
            #   dplyr::mutate(
            #                 #Standard_error = ifelse(!is.na(pev), sqrt(pev), NA),
            #                 Standard_error = Standard_error,
            #                 Prediction_error_variance = PEV,
            #                 Reliability = Reliability)

            ## MET
            # if(!is.null(Zg) & length(datasets)>1){
            #
            #
            #   across_env_predicted_value <- predicted_value[tst, ] |>
            #     dplyr::ungroup() |>
            #     dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
            #     dplyr::summarise(
            #       Predicted_value = mean(Predicted_value, na.rm = TRUE),
            #       Standard_error  = mean(Standard_error,  na.rm = TRUE),
            #       PEV             = mean(PEV,             na.rm = TRUE),
            #       Reliability     = mean(Reliability,     na.rm = TRUE),
            #       .groups         = "drop"
            #     )
            #
            #   across_env_predicted_value <- na.omit(across_env_predicted_value)
            #
            #
            # }        ## End MET

            }
          } ##
        }

    } ## End

  }

  # Phase 3.24: include variance contributions from ALL random ETA terms,
  # not just the first N corresponding to user-supplied kernels. For a
  # formula like ~GID + GID:Loc the ETA has 2 random terms (GID main and
  # GID:Loc Hadamard) but only 1 user kernel (gmatrix). The main loop
  # iterated ETA[[1]] only and missed ETA[[2]]'s Vg, deflating var_u_total
  # and collapsing h^2. Walk the remaining random ETA indices here and
  # add each one's variance contribution. Each branch is defensive --
  # silent NULL on any file-read error so the main result still returns.
  missed_eta_idx <- setdiff(random_eta_idx, kernel_eta_idx)
  # Phase 3.24 follow-up: harmonise chain lengths. The main loop's BRR
  # branch builds var_u_omics from beta_draws (post-burn-in length =
  # nIter - burnIn). BGLR's varU.dat files for RKHS terms contain the
  # FULL chain (length nIter), so process_var_u() returns a longer
  # vector. Mixing them in `var_u_total + extra_vu` triggers
  # "longer object length is not a multiple of shorter object length".
  # Subset any oversized chain by posindex so all terms align on the
  # post-burn-in window.
  align_to_posindex <- function(x) {
    if (is.null(x) || !length(x)) return(x)
    if (length(x) > length(posindex) &&
        all(is.finite(posindex)) &&
        max(posindex) <= length(x)) {
      x[posindex]
    } else {
      x
    }
  }
  # Also align var_u_total + var_residual themselves if they were left at
  # the full-chain length by the main loop's RKHS branch.
  var_u_total <- align_to_posindex(var_u_total)
  var_residual <- align_to_posindex(var_residual)
  se_var_residual <- standard_deviation(var_residual)

  for (extra_idx in missed_eta_idx) {
    term_model_extra <- ETA[["ETA"]][[extra_idx]]$model
    extra_vu <- NULL
    if (identical(term_model_extra, "RKHS")) {
      extra_varU <- varU_files[grepl(paste0("ETA_", extra_idx, "_varU\\.dat$"),
                                     basename(varU_files))]
      if (length(extra_varU) == 1L && file.exists(extra_varU)) {
        extra_vu <- tryCatch(
          process_var_u(extra_varU, posindex, term_model_extra),
          error = function(e) NULL
        )
        eta_kernel_extra <- as.matrix(ETA[["ETA"]][[extra_idx]]$K)
        train_idx_extra <- if (length(tst) > 0) {
          setdiff(seq_len(nrow(eta_kernel_extra)), tst)
        } else {
          seq_len(nrow(eta_kernel_extra))
        }
        kernel_scale_extra <- mean(diag(eta_kernel_extra)[train_idx_extra], na.rm = TRUE)
        if (!is.null(extra_vu) && is.finite(kernel_scale_extra) && kernel_scale_extra > 0) {
          extra_vu <- extra_vu * kernel_scale_extra
        }
      }
    } else if (identical(term_model_extra, "BRR")) {
      extra_coef <- BIN[grepl(paste0("ETA_", extra_idx, "_b\\.bin$"),
                              basename(BIN))]
      if (length(extra_coef) == 1L && file.exists(extra_coef)) {
        extra_vu <- tryCatch({
          beta_draws_extra <- BGLR::readBinMat(extra_coef)
          eta_X_extra <- as.matrix(ETA[["ETA"]][[extra_idx]]$X)
          pred_extra <- eta_X_extra %*% t(beta_draws_extra)
          train_idx_extra <- if (length(tst) > 0) setdiff(seq_len(nrow(eta_X_extra)), tst) else seq_len(nrow(eta_X_extra))
          apply(pred_extra[train_idx_extra, , drop = FALSE], 2, stats::var, na.rm = TRUE)
        }, error = function(e) NULL)
      }
    }
    extra_vu <- align_to_posindex(extra_vu)
    if (!is.null(extra_vu) && length(extra_vu) > 0L) {
      extra_label <- paste0("genetic_variance_eta_", extra_idx)
      var_u_mean_omics_list[[extra_label]] <- mean(extra_vu, na.rm = TRUE)
      se_var_u_omics_list[[extra_label]] <- standard_deviation(extra_vu)
      var_u_total <- var_u_total + extra_vu
    }
  }

  sum_ebv <- combine_bayesian_ebv_components(
    ebv_list = estimated_breeding_value_list,
    gen_name = gen_name,
    heter_groups = heter_groups,
    var_u_total = var_u_total
  )


  observed_trait_rows <- !is.na(mod$model$y)
  if (length(observed_trait_rows) != nrow(pheno_data)) {
    observed_trait_rows <- rep(TRUE, nrow(pheno_data))
  }
  has_met_input <- anyDuplicated(as.character(pheno_data[[gen_name]][observed_trait_rows])) > 0

  predicted_value <- data.frame(
    name = pheno_data[, gen_name],
    Predicted_value = mod$model$yHat,
    Train_Test_Label = train_test_label,
    stringsAsFactors = FALSE
  )
  names(predicted_value)[names(predicted_value) == "name"] <- gen_name
  if (isTRUE(has_met_input)) {
    predicted_value[[heter_groups]] <- as.character(pheno_data[, heter_groups])
  }

  observed_idx <- which(!is.na(mod$model$y))
  residual_value <- data.frame(
    name = pheno_data[observed_idx, gen_name],
    Predicted_value = mod$model$yHat[observed_idx],
    Residual_value = mod$model$y[observed_idx] - mod$model$yHat[observed_idx],
    stringsAsFactors = FALSE
  )
  names(residual_value)[names(residual_value) == "name"] <- gen_name
  if (isTRUE(has_met_input)) {
    residual_value[[heter_groups]] <- as.character(pheno_data[observed_idx, heter_groups])
    residual_value <- residual_value[, c(gen_name, heter_groups, "Predicted_value", "Residual_value"), drop = FALSE]
  }

  random_models <- eta_models_present[eta_models_present != "FIXED"]
  if (length(random_models) && all(random_models == "BRR") &&
      is.matrix(sum_posterior) && ncol(sum_posterior) > 1L) {
    variance_rows <- if (length(tst)) {
      setdiff(seq_len(nrow(sum_posterior)), tst)
    } else {
      seq_len(nrow(sum_posterior))
    }
    var_u_total <- apply(
      sum_posterior[variance_rows, , drop = FALSE],
      2,
      stats::var,
      na.rm = TRUE
    )
  }

  variance_components <- bayes_variance_componentsnew(var_u_mean_omics_list,
                                                      se_var_u_omics_list,
                                                      var_u_total,
                                                      var_residual,
                                                      se_var_residual)

  has_met <- anyDuplicated(as.character(pheno_data[[gen_name]][observed_trait_rows])) > 0
  env_col <- heter_groups
  if (is.null(env_col) && "Env" %in% names(pheno_data)) {
    env_col <- "Env"
  }
  bayes_target_metrics <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = ETA,
    GS_model = GS_model,
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    output_files_names = mod[["output_files_names"]],
    draw_indices = posindex,
    confidence_level = confidence_level,
    CI_width_thresholds = CI_width_thresholds,
    interval_width_high_threshold = interval_width_high_threshold,
    interval_width_low_threshold = interval_width_low_threshold,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  has_met <- has_met || (!is.null(bayes_target_metrics) && !is.null(bayes_target_metrics$across))
  met_covariance_summary <- if (isTRUE(has_met) && !is.null(env_col)) {
    bayes_environment_covariance_summary(
      mod = mod,
      ETA = ETA,
      pheno_data = pheno_data,
      gen_name = gen_name,
      heter_groups = env_col,
      draw_indices = posindex
    )
  } else {
    NULL
  }

  # Phase 3.29: unify the kernel-MET variance_components with the ASReml / GP
  # engines by appending per-environment genetic_variance_<env> /
  # residual_variance_<env> / heritability_<env> rows (see
  # bayes_variance_components_add_per_env). Without this the table only exposed
  # per-ETA-term genetic_variance_<i> rows + one pooled heritability, which does
  # not line up with the environments a user sees in the ASReml/GP outputs.
  if (isTRUE(has_met) && !is.null(env_col) && env_col %in% names(pheno_data) &&
      !is.null(bayes_target_metrics)) {
    variance_components <- bayes_variance_components_add_per_env(
      variance_components = variance_components,
      env_labels = pheno_data[[env_col]],
      marginal_genetic_variance_per_obs = bayes_target_metrics$observation$Genetic_variance,
      varE = mod$model$varE,
      weights = mod$model$weights
    )
  }

  if (!is.null(bayes_target_metrics)) {
    predicted_value$Predicted_value <- bayes_target_metrics$observation$predicted_mean
    if (has_met && !is.null(env_col) && !env_col %in% names(predicted_value)) {
      predicted_value[[env_col]] <- pheno_data[[env_col]]
    }
    predicted_value$Standard_error <- bayes_target_metrics$observation$Standard_error
    predicted_value$PEV <- bayes_target_metrics$observation$PEV
    predicted_value$lower_bound <- bayes_target_metrics$observation$lower_bound
    predicted_value$upper_bound <- bayes_target_metrics$observation$upper_bound
    predicted_value$Uncertainty <- bayes_target_metrics$observation$Uncertainty
    predicted_value$Uncertainty_remarks <- bayes_target_metrics$observation$Uncertainty_remarks
    predicted_value$Genetic_variance <- bayes_target_metrics$observation$Genetic_variance
    predicted_value$Reliability <- bayes_target_metrics$observation$Reliability
    predicted_value$Reliability_remarks <- bayes_target_metrics$observation$Reliability_remarks
    predicted_value$Reliability_percentage <- bayes_target_metrics$observation$Reliability_percentage
    predicted_value <- gp_bayes_attach_interval_provenance(
      prediction_table = predicted_value,
      model_type = GS_model,
      confidence_level = confidence_level,
      posterior_draw_count = length(posindex)
    )
    if (has_met && !is.null(env_col) && env_col %in% names(predicted_value)) {
      reordered_cols <- c(
        gen_name,
        "Predicted_value",
        "Train_Test_Label",
        env_col,
        setdiff(names(predicted_value), c(gen_name, "Predicted_value", "Train_Test_Label", env_col))
      )
      predicted_value <- predicted_value[, unique(reordered_cols), drop = FALSE]
    }

    sum_ebv$Estimated_breeding_value <- bayes_target_metrics$observation$genetic_mean
    sum_ebv$Standard_error <- bayes_target_metrics$observation$Standard_error
    sum_ebv$Prediction_error_variance <- bayes_target_metrics$observation$PEV
    sum_ebv$Genetic_variance <- bayes_target_metrics$observation$Genetic_variance
    sum_ebv$Reliability <- bayes_target_metrics$observation$Reliability

    diagnostic_plots <- diagnostic_plot_true_prediction(
      boot_results = NULL,
      GID_names = pheno_data[, gen_name],
      CI_width_thresholds = CI_width_thresholds,
      predictions = predicted_value$Predicted_value,
      standard_errors = predicted_value$Standard_error,
      prediction_error_var = predicted_value$PEV,
      genetic_var = predicted_value$Genetic_variance,
      confidence_level = confidence_level,
      model_for_CI_cal = "RKHS",
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      system_database = system_database
    )
  }

  if (!is.null(bayes_target_metrics) && has_met && !is.null(bayes_target_metrics$across)) {
    across_env_predicted_value <- data.frame(
      name = bayes_target_metrics$across[[gen_name]],
      Predicted_value = bayes_target_metrics$across$predicted_mean,
      Train_Test_Label = train_test_label_across_env,
      Standard_error = bayes_target_metrics$across$Standard_error,
      PEV = bayes_target_metrics$across$PEV,
      lower_bound = bayes_target_metrics$across$lower_bound,
      upper_bound = bayes_target_metrics$across$upper_bound,
      Uncertainty = bayes_target_metrics$across$Uncertainty,
      Uncertainty_remarks = bayes_target_metrics$across$Uncertainty_remarks,
      Genetic_variance = bayes_target_metrics$across$Genetic_variance,
      Reliability = bayes_target_metrics$across$Reliability,
      Reliability_remarks = bayes_target_metrics$across$Reliability_remarks,
      Reliability_percentage = bayes_target_metrics$across$Reliability_percentage,
      stringsAsFactors = FALSE
    )
    names(across_env_predicted_value)[names(across_env_predicted_value) == "name"] <- gen_name
    across_env_predicted_value <- gp_bayes_attach_interval_provenance(
      prediction_table = across_env_predicted_value,
      model_type = GS_model,
      confidence_level = confidence_level,
      posterior_draw_count = length(posindex)
    )
  }

  if(!"geno_model_ready" %in%names(m_matrix_model_ready_list)){
    if(length(m_matrix_model_ready_list)>1){
      names(m_matrix_model_ready_list) <- paste(paste("Omics", seq_along(m_matrix_model_ready_list), sep = ""), "model_ready", sep = "_")
      names(coefficients_list) <- paste(paste("Omics", seq_along(coefficients_list), sep = ""), "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- paste(paste("Omics", seq_along(estimated_breeding_value_list), sep = ""), "estimated_breeding_value", sep = "_")
    } else {
      names(m_matrix_model_ready_list) <- paste("Omics", "model_ready", sep = "_")
      names(coefficients_list) <- "Coefficient" # paste("Omics", "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- "Estimated_breeding_value"  # paste("Omics", "estimated_breeding_value", sep = "_")

    }
  } else if(length(grep("omic", names(m_matrix_model_ready_list)))==0 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
    geno_index <- grep("geno", names(m_matrix_model_ready_list))
    names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
    names(coefficients_list)[geno_index] <- "Coefficient" # paste("Geno", "coefficient", sep = "_")
    names(estimated_breeding_value_list)[geno_index] <- "Estimated_breeding_value" # paste("Geno", "estimated_breeding_value", sep = "_")

  } else{
    if(length(grep("omic", names(m_matrix_model_ready_list)))>=1 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
      omic_index <- grep("omic", names(m_matrix_model_ready_list))
      geno_index <- grep("geno", names(m_matrix_model_ready_list))
      if(length(omic_index)>1){
        names(m_matrix_model_ready_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "estimated_breeding_value", sep = "_")
      } else {
        names(m_matrix_model_ready_list)[omic_index] <- paste("Omics", "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste("Omics", "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste("Omics", "estimated_breeding_value", sep = "_")

      }
      ####
      names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
      names(coefficients_list)[geno_index] <- paste("Geno", "coefficient", sep = "_")
      names(estimated_breeding_value_list)[geno_index] <- paste("Geno", "estimated_breeding_value", sep = "_")

    }

  }


  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)

    # The genomic layer here is a kernel (gmatrix / gkernel), not geno_data:
    # this function has no geno_data, and the former check stopped every
    # labelled multi-omics kernel fit ("object 'geno_data' not found").
    if (!is.null(gmatrix) || !is.null(gkernel)) {
      print_lable <- c("Genomic", print_lable)
    }
    if(length(print_lable)>length(m_matrix_model_ready_list) | length(print_lable)<length(m_matrix_model_ready_list)){
      message("The number of omics labels does not match the number of data layers; default layer names were used.")



    } else {
      if(length(print_lable)==length(m_matrix_model_ready_list)){
        rownames(variance_components)[1:length(print_lable)] <- paste(print_lable, "variance", sep="_")
        ######
        names(coefficients_list) <-  paste(print_lable, "coefficient", sep="_")
        names(estimated_breeding_value_list) <- paste(print_lable, "estimated_breeding_value", sep="_")
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(m_matrix_model_ready_list) <- paste(print_lable, "model_ready", sep="_")
        ####
      }
    }

  } else {

    if (length(m_matrix_model_ready_list) > 1L) message("No omics labels were provided; default layer names were used.")


  }


  if(!has_met){

    #if(!is.null(diagnostic_plots)){
    res <- list(Coefficients = coefficients_list,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Residual_value = residual_value,
                diagnostic_tst_plot = diagnostic_plots,
                Variance_components = variance_components,
                var_u_total = var_u_total,
                M_matrix_model_ready =  m_matrix_model_ready_list)
    # } else {
    #   res <- list(Coefficients = coefficients_list,
    #               Estimated_breeding_value = estimated_breeding_value_list,
    #               Total_estimated_breeding_value = sum_ebv,
    #               Predicted_value =  predicted_value,
    #               Residual_value = residual_value,
    #               #diagnostic_plots = diagnostic_plots,
    #               Variance_components = variance_components,
    #               M_matrix_model_ready =  m_matrix_model_ready_list)
    # }

  } else {

    # across_env_predicted_value <- predicted_value |>
    #   dplyr::ungroup() |>
    #   dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
    #   dplyr::summarise(
    #     Predicted_value = mean(Predicted_value, na.rm = TRUE),
    #     Standard_error  = mean(Standard_error,  na.rm = TRUE),
    #     PEV             = mean(PEV,             na.rm = TRUE),
    #     Reliability     = mean(Reliability,     na.rm = TRUE),
    #     .groups         = "drop"
    #   )
    #
    # across_env_predicted_value <- na.omit(across_env_predicted_value)
    across_env_predicted_value <- as.data.frame(across_env_predicted_value)

    res <- list(Coefficients = coefficients_list,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Total_Predicted_value = across_env_predicted_value,
                Residual_value = residual_value,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list,
                var_u_total = var_u_total,
                diagnostic_tst_plot = diagnostic_plots
    )

  }
  if (!is.null(met_covariance_summary)) {
    env_levels <- rownames(met_covariance_summary$covariance)
    residual_variance_by_env <- tapply(
      gp_bayes_observation_residual_variance(
        varE = mod$model$varE,
        groups = pheno_data[[env_col]],
        weights = mod$model$weights,
        n = nrow(pheno_data)
      ),
      as.character(pheno_data[[env_col]]),
      mean,
      na.rm = TRUE
    )[env_levels]
    residual_covariance <- diag(
      x = as.numeric(residual_variance_by_env),
      nrow = length(env_levels),
      ncol = length(env_levels)
    )
    dimnames(residual_covariance) <- list(env_levels, env_levels)
    residual_covariance_se <- diag(
      x = rep(se_var_residual, length(env_levels)),
      nrow = length(env_levels),
      ncol = length(env_levels)
    )
    dimnames(residual_covariance_se) <- list(env_levels, env_levels)
    res$Genetic_covariance_environments <- met_covariance_summary$covariance
    res$Genetic_correlation_environments <- met_covariance_summary$correlation
    res$Genetic_covariance_environments_SE <- met_covariance_summary$covariance_se
    res$Genetic_correlation_environments_SE <- met_covariance_summary$correlation_se
    res$Residual_covariance_environments <- residual_covariance
    res$Residual_correlation_environments <- diag(1, nrow = length(env_levels))
    dimnames(res$Residual_correlation_environments) <- list(env_levels, env_levels)
    res$Residual_covariance_environments_SE <- residual_covariance_se
    res$Residual_correlation_environments_SE <- matrix(
      0,
      nrow = length(env_levels),
      ncol = length(env_levels),
      dimnames = list(env_levels, env_levels)
    )
    res$prediction_uncertainty_source <-
      "BGLR posterior draws for BRR targets; analytic mixed-model covariance for RKHS targets"
  }
  ### remove the generated output files from the working directory
  unlink(mod[["output_files_names"]])

  return(res)

}
