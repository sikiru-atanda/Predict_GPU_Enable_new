gp_hybrid_bayes_supported_models <- function() c("GBLUP_BRR", "RKHS")

gp_hybrid_bayes_summary_statistics <- function(pred_df,
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
      "hybrid_bayes",
      "gaussian",
      model_type,
      "hybrid Gaussian Bayesian kernel prediction via GCA_f + GCA_m + SCA",
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

gp_hybrid_bayes_diagnostic_plot <- function(pred_df, model_label) {
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
      title = paste("Hybrid Gaussian", model_label, "observed vs predicted"),
      subtitle = "Bayesian hybrid kernels separate female GCA, male GCA, and SCA",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_hybrid_bayes_cv_plot <- function(pred_df, model_label) {
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
      title = paste("Hybrid", model_label, "cross-validation observed vs predicted"),
      subtitle = "Hybrid scenarios: known parents, one new parent, both parents new",
      x = "Observed value",
      y = "Predicted value",
      color = "Scenario"
    )
}

gp_hybrid_bayes_obs_kernel <- function(levels_vec, gmat) {
  idx <- match(as.character(levels_vec), rownames(gmat))
  if (any(is.na(idx))) {
    stop("Hybrid parent levels did not align with the supplied parent relationship matrix.", call. = FALSE)
  }
  out <- gmat[idx, idx, drop = FALSE]
  rownames(out) <- as.character(levels_vec)
  colnames(out) <- as.character(levels_vec)
  out
}

gp_hybrid_bayes_obs_brr_matrix <- function(levels_vec, gmat, prefix) {
  idx <- match(as.character(levels_vec), rownames(gmat))
  if (any(is.na(idx))) {
    stop("Hybrid parent levels did not align with the supplied parent relationship matrix.", call. = FALSE)
  }
  out <- gmat[idx, , drop = FALSE]
  colnames(out) <- paste0(prefix, "__", colnames(out))
  out
}

gp_hybrid_bayes_eta <- function(ph,
                                female_parent,
                                male_parent,
                                female_gmatrix,
                                male_gmatrix,
                                sca_gmatrix = NULL,
                                heter_groups = NULL,
                                model_type = "GBLUP_BRR",
                                include_sca = TRUE) {
  model_type <- as.character(model_type)[1]
  eta <- list()
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    env_mm <- stats::model.matrix(
      stats::reformulate(heter_groups),
      data = ph
    )
    eta[[length(eta) + 1L]] <- list(
      X = env_mm,
      model = "FIXED",
      saveEffects = TRUE
    )
  }
  if (identical(model_type, "GBLUP_BRR")) {
    eta <- c(eta, list(
      list(
        X = gp_hybrid_bayes_obs_brr_matrix(ph[[female_parent]], female_gmatrix, "female_gca"),
        model = "BRR",
        saveEffects = TRUE
      ),
      list(
        X = gp_hybrid_bayes_obs_brr_matrix(ph[[male_parent]], male_gmatrix, "male_gca"),
        model = "BRR",
        saveEffects = TRUE
      )
    ))
    if (isTRUE(include_sca) && !is.null(sca_gmatrix)) {
      eta[[length(eta) + 1L]] <- list(
        X = gp_hybrid_bayes_obs_brr_matrix(ph$HybridCross, sca_gmatrix, "sca"),
        model = "BRR",
        saveEffects = TRUE
      )
    }
    return(eta)
  }

  eta <- c(eta, list(
    list(
      K = gp_hybrid_bayes_obs_kernel(ph[[female_parent]], female_gmatrix),
      model = "RKHS",
      saveEffects = TRUE
    ),
    list(
      K = gp_hybrid_bayes_obs_kernel(ph[[male_parent]], male_gmatrix),
      model = "RKHS",
      saveEffects = TRUE
    )
  ))
  if (isTRUE(include_sca) && !is.null(sca_gmatrix)) {
    eta[[length(eta) + 1L]] <- list(
      K = gp_hybrid_bayes_obs_kernel(ph$HybridCross, sca_gmatrix),
      model = "RKHS",
      saveEffects = TRUE
    )
  }
  eta
}

gp_hybrid_bayes_component_means <- function(model_fit, eta, model_type) {
  out <- list()
  offset <- 0L
  if (length(eta) && identical(eta[[1]]$model %||% "", "FIXED")) {
    offset <- 1L
  }
  if (identical(model_type, "GBLUP_BRR")) {
    out$female <- drop(eta[[offset + 1L]]$X %*% model_fit$ETA[[offset + 1L]]$b)
    out$male <- drop(eta[[offset + 2L]]$X %*% model_fit$ETA[[offset + 2L]]$b)
    out$sca <- if (length(eta) >= (offset + 3L)) drop(eta[[offset + 3L]]$X %*% model_fit$ETA[[offset + 3L]]$b) else rep(0, nrow(eta[[offset + 1L]]$X))
  } else {
    out$female <- as.numeric(model_fit$ETA[[offset + 1L]]$u)
    out$male <- as.numeric(model_fit$ETA[[offset + 2L]]$u)
    out$sca <- if (length(eta) >= (offset + 3L)) as.numeric(model_fit$ETA[[offset + 3L]]$u) else rep(0, length(out$female))
  }
  out
}

gp_hybrid_bayes_fixed_component <- function(ph, model_fit, eta, mu = 0) {
  if (!length(eta) || !identical(eta[[1]]$model %||% "", "FIXED")) {
    return(rep(as.numeric(mu %||% 0), nrow(ph)))
  }
  beta <- model_fit$ETA[[1]]$b
  if (is.null(beta)) {
    return(rep(as.numeric(mu %||% 0), nrow(ph)))
  }
  as.numeric(mu %||% 0) + drop(eta[[1]]$X %*% beta)
}

gp_hybrid_bayes_brr_draws <- function(mod, eta) {
  bins <- mod$output_files_names[grepl("ETA_[0-9]+_b\\.bin$", basename(mod$output_files_names))]
  if (!length(bins)) {
    return(NULL)
  }
  draw_total <- NULL
  for (eta_idx in seq_along(eta)) {
    coef_file <- bins[grepl(paste0("ETA_", eta_idx, "_b\\.bin$"), basename(bins))]
    if (length(coef_file) != 1L) {
      return(NULL)
    }
    beta_draws <- BGLR::readBinMat(coef_file)
    comp_draws <- eta[[eta_idx]]$X %*% t(beta_draws)
    if (is.null(draw_total)) {
      draw_total <- comp_draws
    } else {
      draw_total <- draw_total + comp_draws
    }
  }
  draw_total
}

gp_hybrid_bayes_variance_components <- function(mod,
                                                 eta,
                                                 model_type,
                                                 draw_indices = NULL) {
  eta_models <- vapply(eta, function(x) as.character(x$model), character(1L))
  random_idx <- which(eta_models != "FIXED")
  if (!length(random_idx)) {
    stop("Hybrid Bayesian variance extraction found no random ETA terms.", call. = FALSE)
  }

  component_labels <- c("female_gca_variance", "male_gca_variance", "sca_variance")
  if (length(random_idx) > length(component_labels)) {
    component_labels <- c(
      component_labels,
      paste0("hybrid_random_component_", seq_len(length(random_idx) - length(component_labels)))
    )
  }
  component_labels <- component_labels[seq_along(random_idx)]
  train_idx <- which(!is.na(mod$model$y))
  if (length(train_idx) < 2L) {
    stop("Hybrid Bayesian variance extraction requires at least two observed rows.", call. = FALSE)
  }

  component_draws <- vector("list", length(random_idx))
  total_effect_draws <- NULL
  for (j in seq_along(random_idx)) {
    eta_idx <- random_idx[[j]]
    term <- eta[[eta_idx]]
    if (identical(model_type, "GBLUP_BRR")) {
      coefficient_file <- mod$output_files_names[grepl(
        paste0("ETA_", eta_idx, "_b\\.bin$"),
        basename(mod$output_files_names)
      )]
      if (length(coefficient_file) != 1L) {
        stop(paste("Unable to locate hybrid BRR coefficient draws for ETA", eta_idx), call. = FALSE)
      }
      coefficient_draws <- as.matrix(BGLR::readBinMat(coefficient_file))
      if (length(draw_indices) && max(draw_indices) <= nrow(coefficient_draws)) {
        coefficient_draws <- coefficient_draws[draw_indices, , drop = FALSE]
      }
      effect_draws <- as.matrix(term$X) %*% t(coefficient_draws)
      component_draws[[j]] <- apply(
        effect_draws[train_idx, , drop = FALSE],
        2L,
        stats::var,
        na.rm = TRUE
      )
      total_effect_draws <- if (is.null(total_effect_draws)) {
        effect_draws
      } else {
        total_effect_draws + effect_draws
      }
    } else {
      variance_draws <- bayes_eta_variance_draws(
        output_files_names = mod$output_files_names,
        eta_index = eta_idx,
        eta_model = "RKHS",
        draw_indices = draw_indices
      )
      if (is.null(variance_draws) || !length(variance_draws)) {
        stop(paste("Unable to locate hybrid RKHS variance draws for ETA", eta_idx), call. = FALSE)
      }
      kernel_scale <- mean(diag(as.matrix(term$K))[train_idx], na.rm = TRUE)
      component_draws[[j]] <- as.numeric(variance_draws) * kernel_scale
    }
  }

  total_genetic_draws <- if (identical(model_type, "GBLUP_BRR")) {
    apply(total_effect_draws[train_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
  } else {
    Reduce(`+`, component_draws)
  }
  residual_draws <- gp_bayes_read_varE_draws(
    mod$output_files_names,
    draw_indices = draw_indices
  )
  common_n <- min(length(total_genetic_draws), length(residual_draws))
  heritability_draws <- total_genetic_draws[seq_len(common_n)] /
    (total_genetic_draws[seq_len(common_n)] + residual_draws[seq_len(common_n)])

  all_draws <- c(
    component_draws,
    list(total_genetic_draws, residual_draws, heritability_draws)
  )
  labels <- c(
    component_labels,
    "total_genetic_variance",
    "residual_variance",
    "heritability"
  )
  data.frame(
    Component = labels,
    Components = vapply(all_draws, function(x) mean(x, na.rm = TRUE), numeric(1L)),
    Standard_error = vapply(all_draws, function(x) stats::sd(x, na.rm = TRUE), numeric(1L)),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_bayes_effect_tables <- function(pred_df,
                                          female_parent,
                                          male_parent) {
  female_eff <- stats::aggregate(
    pred_df$Female_GCA,
    by = list(parent = pred_df[[female_parent]]),
    FUN = mean,
    na.rm = TRUE
  )
  names(female_eff) <- c(female_parent, "Female_GCA")

  male_eff <- stats::aggregate(
    pred_df$Male_GCA,
    by = list(parent = pred_df[[male_parent]]),
    FUN = mean,
    na.rm = TRUE
  )
  names(male_eff) <- c(male_parent, "Male_GCA")

  sca_eff <- stats::aggregate(
    pred_df$SCA_effect,
    by = list(HybridCross = pred_df$HybridCross),
    FUN = mean,
    na.rm = TRUE
  )
  names(sca_eff) <- c("HybridCross", "SCA_effect")

  list(
    female_gca_effects = female_eff,
    male_gca_effects = male_eff,
    sca_effects = sca_eff
  )
}

gp_hybrid_bayes_prediction_frame <- function(ph,
                                             response,
                                             gen_name,
                                             female_parent,
                                             male_parent,
                                             heter_groups = NULL,
                                             fixed_component,
                                             female_comp,
                                             male_comp,
                                             sca_comp,
                                             model_fit,
                                             draw_total = NULL) {
  out_cols <- c(gen_name, female_parent, male_parent, response, "HybridCross")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    out_cols <- c(out_cols, heter_groups)
  }
  out <- ph[, out_cols, drop = FALSE]
  names(out)[names(out) == response] <- "Observed_value"
  out$Train_Test_Label <- ifelse(is.na(out$Observed_value), "Test", "Train")
  out$Prediction_intercept <- as.numeric(fixed_component)
  out$Female_GCA <- as.numeric(female_comp)
  out$Male_GCA <- as.numeric(male_comp)
  out$SCA_effect <- as.numeric(sca_comp)
  out$Female_additive_contribution <- out$Female_GCA
  out$Male_additive_contribution <- out$Male_GCA
  out$Hybrid_interaction_contribution <- out$SCA_effect
  out$Female_GCA_like_contribution <- out$Female_GCA
  out$Male_GCA_like_contribution <- out$Male_GCA
  out$SCA_like_contribution <- out$SCA_effect
  out$Predictive_decomposition_basis <- "hybrid_bayesian_kernel_decomposition"
  out$Predictive_decomposition_remarks <- paste(
    "Bayesian hybrid decomposition:",
    "intercept + female GCA kernel + male GCA kernel + SCA kernel = predicted value"
  )
  out$Predicted_value <- as.numeric(model_fit$yHat)
  if (!is.null(draw_total)) {
    out$Standard_error <- apply(draw_total, 1, stats::sd, na.rm = TRUE)
    out$Prediction_error_variance <- apply(draw_total, 1, stats::var, na.rm = TRUE)
  } else {
    out$Standard_error <- NA_real_
    out$Prediction_error_variance <- NA_real_
  }
  out
}

gp_hybrid_bayes_gaussian_model <- function(pheno_object,
                                           response,
                                           gen_name,
                                           female_parent,
                                           male_parent,
                                           heter_groups = NULL,
                                           model_type = "GBLUP_BRR",
                                           gmatrix = NULL,
                                           female_gmatrix = NULL,
                                           male_gmatrix = NULL,
                                           geno_data = NULL,
                                           female_geno_data = NULL,
                                           male_geno_data = NULL,
                                           gmatrix_method = NULL,
                                           include_sca = TRUE,
                                           inverse = TRUE,
                                           epsilon = 1e-6,
                                           nIter = NULL,
                                           burnIn = NULL,
                                           thin = NULL,
                                           confidence_level = 0.95,
                                           CI_width_thresholds = c(0.33, 0.66),
                                           interval_width_high_threshold = NULL,
                                           interval_width_low_threshold = NULL,
                                           high_reliability_thres = 0.9,
                                           low_reliability_thres = 0.5,
                                           random_state = NULL) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("Install the BGLR package to use hybrid Bayesian models.", call. = FALSE)
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
      stop("Hybrid Bayesian models currently require one phenotype row per hybrid ID unless heter_groups is supplied for multi-environment data.", call. = FALSE)
    }
  } else if (anyDuplicated(gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups))) {
    stop("Hybrid Bayesian models currently require one phenotype row per hybrid-by-environment combination.", call. = FALSE)
  }

  ph[[female_parent]] <- as.character(ph[[female_parent]])
  ph[[male_parent]] <- as.character(ph[[male_parent]])
  ph[[gen_name]] <- as.character(ph[[gen_name]])
  ph$HybridCross <- paste(ph[[female_parent]], ph[[male_parent]], sep = "__x__")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    ph[[heter_groups]] <- as.character(ph[[heter_groups]])
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

  eta <- gp_hybrid_bayes_eta(
    ph = ph,
    female_parent = female_parent,
    male_parent = male_parent,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    heter_groups = heter_groups,
    model_type = model_type,
    include_sca = include_sca
  )
  bayes_para <- bayes_parameter_check(
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    message = FALSE
  )
  # The MCMC chain is seeded from random_state so repeat runs are identical
  # (it was unseeded: two identical runs differed by up to 0.16).
  mod <- gp_with_pinned_seed(random_state, bayes_mod_execute(
    pheno_data = ph,
    response = response,
    ETA = eta,
    bayes_para = bayes_para,
    GS_model = if (identical(model_type, "GBLUP_BRR")) "BRR" else "RKHS",
    response_family = "gaussian"
  ))
  on.exit(unlink(mod$output_files_names, force = TRUE), add = TRUE)
  draw_indices <- gp_bayes_saved_draw_indices(
    nIter = bayes_para[["nIter"]],
    burnIn = bayes_para[["burnIn"]],
    thin = bayes_para[["thin"]]
  )

  components <- gp_hybrid_bayes_component_means(
    model_fit = mod$model,
    eta = eta,
    model_type = model_type
  )
  draw_total <- if (identical(model_type, "GBLUP_BRR")) {
    gp_hybrid_bayes_brr_draws(mod = mod, eta = eta)
  } else {
    NULL
  }
  fixed_component <- gp_hybrid_bayes_fixed_component(
    ph = ph,
    model_fit = mod$model,
    eta = eta,
    mu = mod$model$mu %||% 0
  )
  pred_df <- gp_hybrid_bayes_prediction_frame(
    ph = ph,
    response = response,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    fixed_component = fixed_component,
    female_comp = components$female,
    male_comp = components$male,
    sca_comp = components$sca,
    model_fit = mod$model,
    draw_total = draw_total
  )
  target_metrics <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = list(ETA = eta),
    GS_model = if (identical(model_type, "GBLUP_BRR")) "BRR" else "RKHS",
    pheno_data = ph,
    gen_name = gen_name,
    heter_groups = heter_groups,
    output_files_names = mod$output_files_names,
    draw_indices = draw_indices,
    confidence_level = confidence_level,
    CI_width_thresholds = CI_width_thresholds,
    interval_width_high_threshold = interval_width_high_threshold,
    interval_width_low_threshold = interval_width_low_threshold,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )[["observation"]]
  pred_df$Predicted_value <- target_metrics$predicted_mean
  pred_df$Standard_error <- target_metrics$Standard_error
  pred_df$PEV <- target_metrics$PEV
  pred_df$Prediction_error_variance <- target_metrics$PEV
  pred_df$lower_bound <- target_metrics$lower_bound
  pred_df$upper_bound <- target_metrics$upper_bound
  pred_df$Uncertainty <- target_metrics$Uncertainty
  pred_df$Uncertainty_remarks <- target_metrics$Uncertainty_remarks
  pred_df$Reliability <- target_metrics$Reliability
  pred_df$Reliability_remarks <- target_metrics$Reliability_remarks
  pred_df$Reliability_variance_input <- target_metrics$PEV
  pred_df$Reliability_reference_variance <- target_metrics$Genetic_variance
  pred_df$PEV_basis <- if (identical(model_type, "GBLUP_BRR")) {
    "BGLR posterior target-draw variance"
  } else {
    "BGLR RKHS analytic posterior prediction variance"
  }
  pred_df$Prediction_uncertainty_source <- pred_df$PEV_basis
  pred_df$Prediction_interval_method <- if (identical(model_type, "GBLUP_BRR")) {
    "BGLR posterior target-draw quantiles"
  } else {
    "Gaussian approximation from BGLR RKHS posterior covariance"
  }
  pred_df$Prediction_interval_nominal_coverage <- confidence_level
  pred_df$Prediction_interval_calibration_n <- if (identical(model_type, "GBLUP_BRR")) {
    length(draw_indices)
  } else {
    NA_integer_
  }
  pred_df$Reliability_basis <- paste0(
    "BGLR model-implied hybrid reliability = max(0, min(1, 1 - PEV/genetic_variance)); ",
    "PEV and genetic_variance are target-specific posterior quantities; ",
    "if that form degenerates to a flat zero while PEV varies, the formatter uses ",
    "genetic_variance/(genetic_variance + PEV) and marks the affected rows in Reliability_remarks"
  )
  variance_components <- gp_hybrid_bayes_variance_components(
    mod = mod,
    eta = eta,
    model_type = model_type,
    draw_indices = draw_indices
  )
  eff_tables <- gp_hybrid_bayes_effect_tables(
    pred_df = pred_df,
    female_parent = female_parent,
    male_parent = male_parent
  )

  list(
    bayes_model = mod$model,
    predicted_values = pred_df,
    Predicted_value = pred_df,
    female_gca_effects = eff_tables$female_gca_effects,
    male_gca_effects = eff_tables$male_gca_effects,
    sca_effects = eff_tables$sca_effects,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    variance_components = variance_components,
    bayes_parameter = bayes_para,
    diagnostic_plots = gp_hybrid_bayes_diagnostic_plot(pred_df, model_label = model_type)
  )
}

gp_hybrid_bayes_cv_metric_table <- function(pred_df,
                                            eval_metrics,
                                            response_family = "gaussian") {
  gp_hybrid_asreml_cv_metric_table(
    pred_df = pred_df,
    eval_metrics = eval_metrics,
    response_family = response_family
  )
}

gp_hybrid_bayes_cv_process <- function(pred_df,
                                       eval_df,
                                       response,
                                       model_type) {
  out <- gp_hybrid_asreml_cv_process(
    pred_df = pred_df,
    eval_df = eval_df,
    response = response,
    model_type = model_type
  )
  out$hybrid_cv_plot <- gp_hybrid_bayes_cv_plot(pred_df, model_label = model_type)
  out
}

gp_hybrid_bayes_gaussian_cv <- function(pheno_object,
                                        response,
                                        gen_name,
                                        female_parent,
                                        male_parent,
                                        heter_groups = NULL,
                                        cross_validation_meth,
                                        eval_metrics,
                                        model_type = "GBLUP_BRR",
                                        gmatrix = NULL,
                                        female_gmatrix = NULL,
                                        male_gmatrix = NULL,
                                        geno_data = NULL,
                                        female_geno_data = NULL,
                                        male_geno_data = NULL,
                                        gmatrix_method = NULL,
                                        include_sca = TRUE,
                                        inverse = TRUE,
                                        epsilon = 1e-6,
                                        nIter = NULL,
                                        burnIn = NULL,
                                        thin = NULL,
                                        confidence_level = 0.95,
                                        CI_width_thresholds = c(0.33, 0.66),
                                        interval_width_high_threshold = NULL,
                                        interval_width_low_threshold = NULL,
                                        high_reliability_thres = 0.9,
                                        low_reliability_thres = 0.5,
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

      fit <- gp_hybrid_bayes_gaussian_model(
        pheno_object = ph_cv,
        response = response,
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        model_type = model_type,
        gmatrix = gmatrix,
        female_gmatrix = female_gmatrix,
        male_gmatrix = male_gmatrix,
        geno_data = geno_data,
        female_geno_data = female_geno_data,
        male_geno_data = male_geno_data,
        gmatrix_method = gmatrix_method,
        include_sca = include_sca,
        inverse = inverse,
        epsilon = epsilon,
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        confidence_level = confidence_level,
        CI_width_thresholds = CI_width_thresholds,
        interval_width_high_threshold = interval_width_high_threshold,
        interval_width_low_threshold = interval_width_low_threshold,
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres,
        random_state = (as.numeric(random_state %||% 123L) + rep_i * 10000 +
                          match(list(sc), scenarios)) %% .Machine$integer.max
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

      eval_df <- gp_hybrid_bayes_cv_metric_table(
        pred_df = pred_test,
        eval_metrics = eval_metrics,
        response_family = "gaussian"
      )

      all_preds[[length(all_preds) + 1L]] <- pred_test
      all_eval[[length(all_eval) + 1L]] <- eval_df
      all_raw[[length(all_raw) + 1L]] <- list(
        trait = response,
        rep = rep_i,
        model = model_type,
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
  processed <- gp_hybrid_bayes_cv_process(
    pred_df = pred_df,
    eval_df = eval_df,
    response = response,
    model_type = model_type
  )

  list(
    predicted_values = pred_df,
    hybrid_cv_metrics = eval_df,
    cv_results_processed = processed,
    cv_results_raw = all_raw,
    diagnostic_plots = processed$hybrid_cv_plot
  )
}
