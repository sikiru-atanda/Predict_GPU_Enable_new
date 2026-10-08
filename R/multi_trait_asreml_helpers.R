gp_multitrait_asreml_supported_models <- function() c("GBLUP")

gp_multitrait_asreml_summary_statistics <- function(pred_long, response, model_type) {
  counts <- do.call(rbind, lapply(response, function(tr) {
    dat <- pred_long[pred_long$Trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      observed_count = sum(dat$Train_Test_Label == "Train", na.rm = TRUE),
      predicted_count = sum(dat$Train_Test_Label == "Test", na.rm = TRUE),
      total_count = nrow(dat),
      stringsAsFactors = FALSE
    )
  }))
  data.frame(
    stat = c(
      "mode",
      "response_family",
      "model_type",
      "n_traits",
      "traits",
      "scope",
      paste0("observed_count_", counts$trait),
      paste0("predicted_count_", counts$trait)
    ),
    summary = c(
      "multi_trait_asreml",
      "gaussian",
      model_type,
      length(response),
      paste(response, collapse = ", "),
      "unbalanced gaussian true prediction via multivariate GBLUP + ASReml-R",
      counts$observed_count,
      counts$predicted_count
    ),
    stringsAsFactors = FALSE
  )
}

gp_multitrait_asreml_diagnostic_plot <- function(pred_long) {
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~Trait, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Multi-trait Gaussian ASReml-R observed vs predicted",
      subtitle = "Observed cells contribute to the multivariate fit; missing cells are predicted from the shared multivariate GBLUP model",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_multitrait_asreml_long_data <- function(pheno_object, response, gen_name) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  if (!gen_name %in% names(ph)) {
    stop("The genotype ID column was not found in the phenotype data.", call. = FALSE)
  }
  if (anyDuplicated(ph[[gen_name]])) {
    stop("Multi-trait ASReml-R currently requires one phenotype row per genotype.", call. = FALSE)
  }

  keep_cols <- c(gen_name, response)
  ph <- ph[, keep_cols, drop = FALSE]
  long <- stats::reshape(
    ph,
    varying = response,
    v.names = "Observed_value",
    timevar = "Trait",
    times = response,
    idvar = gen_name,
    direction = "long"
  )
  rownames(long) <- NULL
  names(long)[names(long) == gen_name] <- gen_name
  genotype_levels <- as.character(ph[[gen_name]])
  long[[gen_name]] <- factor(as.character(long[[gen_name]]), levels = genotype_levels)
  long$Trait <- factor(as.character(long$Trait), levels = response)

  # ASReml's stacked multivariate residual expects every unit's traits to be
  # adjacent. stats::reshape() emits trait-major blocks, so establish the
  # genotype-major, trait-minor order before the data reaches asreml().
  long <- long[order(long[[gen_name]], long$Trait), , drop = FALSE]
  long[[gen_name]] <- as.character(long[[gen_name]])
  long$Trait <- as.character(long$Trait)
  rownames(long) <- NULL
  long
}

gp_multitrait_asreml_prediction_wide <- function(pred_long, gen_name) {
  out <- stats::reshape(
    pred_long[, c(gen_name, "Trait", "Predicted_value"), drop = FALSE],
    idvar = gen_name,
    timevar = "Trait",
    direction = "wide"
  )
  names(out) <- sub("^Predicted_value\\.", "", names(out))
  out
}

gp_multitrait_asreml_trait_counts <- function(pred_long, response) {
  do.call(rbind, lapply(response, function(tr) {
    dat <- pred_long[pred_long$Trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      observed_count = sum(dat$Train_Test_Label == "Train", na.rm = TRUE),
      predicted_count = sum(dat$Train_Test_Label == "Test", na.rm = TRUE),
      total_count = nrow(dat),
      stringsAsFactors = FALSE
    )
  }))
}

gp_multitrait_asreml_regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

gp_multitrait_asreml_component_trait <- function(component, response) {
  component <- as.character(component)
  vapply(component, function(label) {
    hits <- response[vapply(response, function(tr) {
      pattern <- paste0("(^|[^[:alnum:]])", gp_multitrait_asreml_regex_escape(tr), "([^[:alnum:]]|$)")
      grepl(pattern, label)
    }, logical(1L))]
    if (length(hits) == 1L) hits[[1L]] else NA_character_
  }, character(1L))
}

gp_multitrait_asreml_component_trait_pair <- function(component, response) {
  component <- as.character(component)
  lapply(component, function(label) {
    hits <- do.call(rbind, lapply(response, function(tr) {
      pattern <- paste0("(^|[^[:alnum:]])", gp_multitrait_asreml_regex_escape(tr), "([^[:alnum:]]|$)")
      loc <- gregexpr(pattern, label, perl = TRUE)[[1L]]
      loc <- loc[loc > 0L]
      if (!length(loc)) {
        return(NULL)
      }
      data.frame(trait = tr, position = loc, stringsAsFactors = FALSE)
    }))
    if (!is.data.frame(hits) || !nrow(hits)) {
      return(character())
    }
    hits <- hits[order(hits$position), , drop = FALSE]
    unique_hits <- unique(as.character(hits$trait))
    if (length(unique_hits) >= 2L) {
      unique_hits[seq_len(2L)]
    } else {
      rep(unique_hits[[1L]], 2L)
    }
  })
}

gp_multitrait_asreml_covariance_outputs <- function(model,
                                                     response,
                                                     kernel_specs = NULL,
                                                     trait_covariance = NULL) {
  raw <- tryCatch(asreml_varcomp_table(summary(model)), error = function(e) NULL)
  vc <- gp_format_variance_components(raw, default_estimation_method = "ASReml_REML")
  if (!is.data.frame(vc) || !nrow(vc)) {
    return(list())
  }

  if (!is.null(kernel_specs) && nrow(kernel_specs)) {
    trait_covariance <- tolower(as.character(trait_covariance %||% "us")[[1L]])
    genetic_by_kernel <- stats::setNames(vector("list", nrow(kernel_specs)), kernel_specs$Kernel)
    for (i in seq_len(nrow(kernel_specs))) {
      rows <- which(grepl(kernel_specs$inverse_name[[i]], rownames(raw), fixed = TRUE))
      if (!length(rows)) {
        stop(
          "No ASReml variance components matched kernel `", kernel_specs$Kernel[[i]],
          "` (", kernel_specs$inverse_name[[i]], ").", call. = FALSE
        )
      }
      genetic_by_kernel[[i]] <- gp_multitrait_asreml_covariance_from_rows(
        raw = raw,
        response = response,
        rows = rows,
        structure = trait_covariance
      )
    }
    genetic_cov <- Reduce(`+`, genetic_by_kernel)
    residual_rows <- which(
      grepl("units|!r|resid|residual", rownames(raw), ignore.case = TRUE) &
        !grepl("^units:trait!r$", gsub("\\s+", "", rownames(raw)), ignore.case = TRUE)
    )
    residual_cov <- gp_multitrait_asreml_covariance_from_rows(
      raw = raw,
      response = response,
      rows = residual_rows,
      structure = trait_covariance
    )
    build_cor <- function(mat) {
      if (all(is.finite(diag(mat))) && all(diag(mat) > 0)) stats::cov2cor(mat) else NULL
    }
    out <- list(
      Genetic_covariance_traits = genetic_cov,
      Genetic_covariance_by_kernel = genetic_by_kernel,
      Residual_covariance_traits = residual_cov
    )
    genetic_cor <- build_cor(genetic_cov)
    residual_cor <- build_cor(residual_cov)
    if (!is.null(genetic_cor)) out[["Genetic_correlation_traits"]] <- genetic_cor
    if (!is.null(residual_cor)) out[["Residual_correlation_traits"]] <- residual_cor
    return(out)
  }

  component_raw <- as.character(vc[["Component"]])
  component_label <- tolower(component_raw)
  component_pairs <- gp_multitrait_asreml_component_trait_pair(component_raw, response = response)
  values <- suppressWarnings(as.numeric(vc[["Components"]]))
  residual <- grepl("units|!r|resid|residual", component_label)
  genetic <- !residual & grepl("vm\\(|gmatrix|genetic|gid|geno", component_label)

  build_cov <- function(mask) {
    mat <- matrix(NA_real_, nrow = length(response), ncol = length(response))
    dimnames(mat) <- list(response, response)
    for (i in which(mask & is.finite(values))) {
      pair <- component_pairs[[i]]
      if (length(pair) != 2L || !all(pair %in% response)) {
        next
      }
      mat[pair[[1L]], pair[[2L]]] <- values[[i]]
      mat[pair[[2L]], pair[[1L]]] <- values[[i]]
    }
    if (!any(is.finite(mat))) {
      return(NULL)
    }
    mat
  }

  build_cor <- function(mat) {
    if (is.null(mat)) {
      return(NULL)
    }
    diag_vals <- diag(mat)
    if (!all(is.finite(diag_vals) & diag_vals > 0)) {
      return(NULL)
    }
    out <- tryCatch(stats::cov2cor(mat), error = function(e) NULL)
    if (is.null(out) || !any(is.finite(out))) {
      return(NULL)
    }
    out
  }

  genetic_cov <- build_cov(genetic)
  residual_cov <- build_cov(residual)
  out <- list()
  if (!is.null(genetic_cov)) {
    out[["Genetic_covariance_traits"]] <- genetic_cov
    genetic_cor <- build_cor(genetic_cov)
    if (!is.null(genetic_cor)) {
      out[["Genetic_correlation_traits"]] <- genetic_cor
    }
  }
  if (!is.null(residual_cov)) {
    out[["Residual_covariance_traits"]] <- residual_cov
    residual_cor <- build_cor(residual_cov)
    if (!is.null(residual_cor)) {
      out[["Residual_correlation_traits"]] <- residual_cor
    }
  }
  out
}

gp_multitrait_asreml_variance_components <- function(model,
                                                       response,
                                                       kernel_specs = NULL) {
  raw <- tryCatch(asreml_varcomp_table(summary(model)), error = function(e) NULL)
  vc <- gp_format_variance_components(raw, default_estimation_method = "ASReml_REML")
  if (!is.data.frame(vc) || !nrow(vc)) {
    return(vc)
  }

  component_raw <- as.character(vc[["Component"]])
  # In sigma-parameterized multivariate fits ASReml reports the fixed scale
  # row `units:Trait!R = 1` alongside the estimable US residual covariance
  # parameters.  It is a bookkeeping/dispersion parameter (bound `F`), not a
  # residual variance or covariance.  Exclude it from the public biological
  # variance-component table instead of relabelling it residual_covariance.
  fixed_residual_scale <- grepl(
    "^units:trait!r$",
    gsub("\\s+", "", component_raw),
    ignore.case = TRUE
  )
  if (any(fixed_residual_scale)) {
    vc <- vc[!fixed_residual_scale, , drop = FALSE]
    component_raw <- as.character(vc[["Component"]])
  }
  trait <- gp_multitrait_asreml_component_trait(component_raw, response = response)
  label <- tolower(component_raw)
  residual <- grepl("units|!r|resid|residual", label)
  genetic <- !residual & grepl("vm\\(|gmatrix|genetic|gid|geno", label)

  component_public <- component_raw
  component_public[genetic & !is.na(trait)] <- "genetic_variance"
  component_public[residual & !is.na(trait)] <- "residual_variance"
  component_public[genetic & is.na(trait)] <- "genetic_covariance"
  component_public[residual & is.na(trait)] <- "residual_covariance"

  vc[["Trait"]] <- trait
  vc[["Component"]] <- component_public
  vc[["Estimation_method"]] <- "ASReml_REML"
  if (!is.null(kernel_specs) && nrow(kernel_specs)) {
    kernel <- rep(NA_character_, length(component_raw))
    for (i in seq_len(nrow(kernel_specs))) {
      kernel[grepl(kernel_specs$inverse_name[[i]], component_raw, fixed = TRUE)] <-
        kernel_specs$Kernel[[i]]
    }
    vc[["Kernel"]] <- kernel
    vc <- vc[, c("Trait", "Kernel", setdiff(names(vc), c("Trait", "Kernel"))), drop = FALSE]
  } else {
    vc <- vc[, c("Trait", setdiff(names(vc), "Trait")), drop = FALSE]
  }
  rownames(vc) <- NULL
  vc
}

gp_multitrait_asreml_prediction_from_pvals <- function(pred_raw, response, gen_name) {
  pred_raw <- as.data.frame(pred_raw, stringsAsFactors = FALSE)
  pred_cols <- names(pred_raw)
  trait_col <- pred_cols[tolower(pred_cols) %in% "trait"][1]
  gid_col <- pred_cols[tolower(pred_cols) %in% tolower(gen_name)][1]
  value_col <- pred_cols[tolower(pred_cols) %in% c("predicted.value", "predicted_value")][1]
  se_col <- pred_cols[tolower(pred_cols) %in% c("std.error", "standard_error")][1]
  if (is.na(trait_col) || is.na(gid_col) || is.na(value_col)) {
    stop("Could not identify Trait/GID/predicted.value columns in the ASReml multivariate prediction table.", call. = FALSE)
  }

  pred_df <- data.frame(
    id_value = as.character(pred_raw[[gid_col]]),
    Trait = as.character(pred_raw[[trait_col]]),
    Predicted_value = as.numeric(pred_raw[[value_col]]),
    Standard_error = if (!is.na(se_col)) as.numeric(pred_raw[[se_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  names(pred_df)[1] <- gen_name
  pred_df <- pred_df[pred_df$Trait %in% response, , drop = FALSE]
  pred_df$Trait <- factor(pred_df$Trait, levels = response)
  pred_df <- pred_df[order(pred_df$Trait, pred_df[[gen_name]]), , drop = FALSE]
  pred_df$Trait <- as.character(pred_df$Trait)
  pred_df$Prediction_error_variance <- pred_df$Standard_error^2
  rownames(pred_df) <- NULL
  pred_df
}

gp_multitrait_asreml_strip_trait_from_genotype <- function(x, trait) {
  x <- as.character(x)
  trait <- gp_multitrait_asreml_regex_escape(trait)
  x <- sub(paste0("^Trait[_:.]", trait, "[_:.]"), "", x)
  x <- sub(paste0("^", trait, "[_:.]"), "", x)
  x <- sub(paste0("[_:.]Trait[_:.]", trait, "$"), "", x)
  x <- sub(paste0("[_:.]", trait, "$"), "", x)
  x
}

gp_multitrait_asreml_genotype_from_component <- function(component, trait, gen_name) {
  component <- as.character(component)
  vm_rx <- paste0(".*vm\\(", gp_multitrait_asreml_regex_escape(gen_name), "\\s*,\\s*[^)]*\\)_")
  suffix <- sub(vm_rx, "", component)
  suffix <- sub(":.*$", "", suffix)
  mapply(
    gp_multitrait_asreml_strip_trait_from_genotype,
    x = suffix,
    trait = trait,
    USE.NAMES = FALSE
  )
}

gp_multitrait_asreml_fixed_trait_effects <- function(fixed_coefs, response) {
  effects <- stats::setNames(rep(0, length(response)), response)
  if (is.null(fixed_coefs) || !nrow(fixed_coefs)) {
    return(effects)
  }

  fixed_effect_col <- pp_asreml_coef_col(fixed_coefs, c("effect", "solution"), "fixed-effect solution")
  intercept <- if ("(Intercept)" %in% rownames(fixed_coefs)) {
    as.numeric(fixed_coefs["(Intercept)", fixed_effect_col])
  } else {
    0
  }
  effects[] <- intercept

  fixed_rows <- setdiff(rownames(fixed_coefs), "(Intercept)")
  if (!length(fixed_rows)) {
    return(effects)
  }

  trait <- gp_multitrait_asreml_component_trait(fixed_rows, response = response)
  keep <- !is.na(trait)
  if (any(keep)) {
    effects[trait[keep]] <- intercept + as.numeric(fixed_coefs[fixed_rows[keep], fixed_effect_col])
  }
  effects
}

gp_multitrait_asreml_aggregate_random <- function(random_df) {
  if (!nrow(random_df)) {
    return(random_df)
  }
  stats::aggregate(
    random_df[, c("Effect", "PEV_component"), drop = FALSE],
    random_df[, c("Trait", "Genotype"), drop = FALSE],
    function(x) sum(x, na.rm = TRUE)
  )
}

# Shared by true prediction and CV so both accept the same (case-insensitive)
# single-environment trait covariance structures.
gp_multitrait_asreml_trait_covariance <- function(var_cov_str) {
  trait_covariance <- if (is.null(var_cov_str)) "us" else tolower(as.character(var_cov_str)[[1L]])
  if (length(var_cov_str %||% trait_covariance) != 1L ||
      !trait_covariance %in% c("us", "corgh", "diag")) {
    stop("Single-environment multi-trait ASReml-R var_cov_str must be one of: us, corgh, diag.", call. = FALSE)
  }
  trait_covariance
}

gp_multitrait_asreml_workspace_error <- function(error) {
  grepl("insufficient workspace|pworkspace", as.character(error), ignore.case = TRUE)
}

gp_multitrait_asreml_retry_pworkspace <- function(pworkspace) {
  value <- as.character(pworkspace %||% "")[[1L]]
  if (!nzchar(value)) return("1gb")
  token <- tolower(gsub("\\s+", "", value))
  # ASReml-R reads a bare number as double-precision words (8 bytes each).
  bytes <- suppressWarnings(as.numeric(token)) * 8
  if (!is.finite(bytes)) {
    unit <- sub("^[0-9.]+", "", token)
    amount <- suppressWarnings(as.numeric(sub("(kb|mb|gb)$", "", token)))
    multiplier <- switch(unit, kb = 1e3, mb = 1e6, gb = 1e9, NA_real_)
    bytes <- amount * multiplier
  }
  if (!is.finite(bytes) || bytes >= 1e9) return(NULL)
  "1gb"
}

gp_multitrait_asreml_predict_pvals <- function(mod, predict_fn, classify,
                                               pworkspace = NULL, levels = NULL) {
  call_predict <- function(memory) {
    if (!is.null(levels) && asreml_has_arg_value(memory)) {
      return(predict_fn(mod, classify = classify, sed = FALSE,
                        levels = levels, pworkspace = memory[[1L]])$pvals)
    }
    if (!is.null(levels)) {
      return(predict_fn(mod, classify = classify, sed = FALSE,
                        levels = levels)$pvals)
    }
    if (asreml_has_arg_value(memory)) {
      return(predict_fn(mod, classify = classify, sed = FALSE,
                        pworkspace = memory[[1L]])$pvals)
    }
    predict_fn(mod, classify = classify, sed = FALSE)$pvals
  }
  first <- tryCatch(call_predict(pworkspace), error = identity)
  if (!inherits(first, "error")) return(first)
  if (!gp_multitrait_asreml_workspace_error(conditionMessage(first))) stop(first)
  if (!isTRUE(mod$converge)) {
    stop("ASReml multi-trait prediction workspace retry is unsafe because the fitted model did not converge.",
         call. = FALSE)
  }
  asreml_validate_coefficient_fallback_model(
    mod, context = "ASReml multi-trait prediction workspace retry"
  )
  retry_memory <- gp_multitrait_asreml_retry_pworkspace(pworkspace)
  if (is.null(retry_memory)) stop(first)
  second <- tryCatch(call_predict(retry_memory), error = identity)
  if (inherits(second, "error")) {
    stop("ASReml prediction failed at the supplied pworkspace and at ",
         retry_memory, ": ", conditionMessage(second), call. = FALSE)
  }
  second
}

extract_asreml_multitrait_prediction <- function(model,
                                                 response,
                                                 gen_name = "GID",
                                                 kernel_ids = NULL) {
  asreml_validate_coefficient_fallback_model(
    model = model,
    context = "ASReml multi-trait coefficient fallback"
  )
  random_coefs <- summary(model, coef = TRUE)$coef.random
  fixed_coefs <- coef(model)$fixed

  if (is.null(random_coefs) || !nrow(random_coefs)) {
    stop("No ASReml random coefficients are available for multi-trait fallback prediction extraction.", call. = FALSE)
  }

  coef_rows <- rownames(random_coefs)
  solution_col <- pp_asreml_coef_col(random_coefs, c("solution", "effect"), "solution")
  se_col <- if (any(c("std.error", "Standard_error") %in% colnames(random_coefs))) {
    pp_asreml_coef_col(random_coefs, c("std.error", "Standard_error"), "standard error")
  } else {
    NULL
  }

  geno_rx <- paste0("vm\\(", gp_multitrait_asreml_regex_escape(gen_name), "\\s*,")
  vm_rows <- coef_rows[grepl(geno_rx, coef_rows)]
  if (!is.null(kernel_ids) && length(kernel_ids)) {
    kernel_rx <- paste(gp_multitrait_asreml_regex_escape(as.character(kernel_ids)), collapse = "|")
    vm_rows <- vm_rows[grepl(
      paste0("vm\\(", gp_multitrait_asreml_regex_escape(gen_name), "\\s*,\\s*(", kernel_rx, ")\\)"),
      vm_rows
    )]
  }
  if (!length(vm_rows)) {
    stop("No ASReml random coefficient rows matched the multi-trait genotype term.", call. = FALSE)
  }

  trait <- gp_multitrait_asreml_component_trait(vm_rows, response = response)
  keep <- !is.na(trait)
  if (!any(keep)) {
    stop("No ASReml random coefficient rows could be assigned to the requested traits.", call. = FALSE)
  }

  vm_rows <- vm_rows[keep]
  trait <- trait[keep]
  genotype <- gp_multitrait_asreml_genotype_from_component(
    component = vm_rows,
    trait = trait,
    gen_name = gen_name
  )
  random_df <- data.frame(
    Trait = trait,
    Genotype = genotype,
    Effect = as.numeric(random_coefs[vm_rows, solution_col]),
    PEV_component = if (!is.null(se_col)) as.numeric(random_coefs[vm_rows, se_col])^2 else NA_real_,
    stringsAsFactors = FALSE
  )
  random_df <- random_df[nzchar(random_df$Genotype), , drop = FALSE]
  if (!nrow(random_df)) {
    stop("No genotype labels could be parsed from the ASReml multi-trait coefficient rows.", call. = FALSE)
  }
  if (!is.null(se_col) && any(!is.finite(random_df$PEV_component))) {
    stop(
      paste(
        "ASReml multi-trait coefficient fallback is unsafe because at least one",
        "matched random-effect coefficient is aliased or non-estimable.",
        "Use direct ASReml predict() for this model."
      ),
      call. = FALSE
    )
  }

  genetic <- gp_multitrait_asreml_aggregate_random(random_df)
  names(genetic)[names(genetic) == "Genotype"] <- gen_name
  fixed_effects <- gp_multitrait_asreml_fixed_trait_effects(fixed_coefs, response = response)
  genetic$Predicted_value <- fixed_effects[genetic$Trait] + genetic$Effect
  genetic$Prediction_error_variance <- ifelse(
    is.na(genetic$PEV_component),
    NA_real_,
    pmax(0, genetic$PEV_component)
  )
  genetic$Standard_error <- sqrt(genetic$Prediction_error_variance)
  genetic$Trait <- factor(genetic$Trait, levels = response)
  genetic <- genetic[order(genetic$Trait, genetic[[gen_name]]), , drop = FALSE]
  genetic$Trait <- as.character(genetic$Trait)
  out <- genetic[, c(gen_name, "Trait", "Predicted_value", "Standard_error", "Prediction_error_variance"), drop = FALSE]
  rownames(out) <- NULL
  out
}

gp_multitrait_asreml_predict_or_extract <- function(mod,
                                                    predict_fn,
                                                    classify,
                                                    response,
                                                    gen_name,
                                                    kernel_ids = NULL,
                                                    pworkspace = NULL) {
  direct_prediction <- tryCatch({
    pred_raw <- gp_multitrait_asreml_predict_pvals(
      mod = mod, predict_fn = predict_fn, classify = classify,
      pworkspace = pworkspace
    )
    list(
      mode = "predict",
      prediction = gp_multitrait_asreml_prediction_from_pvals(
        pred_raw = pred_raw,
        response = response,
        gen_name = gen_name
      ),
      error = NULL
    )
  }, error = function(e) {
    list(mode = "predict_error", prediction = NULL, error = conditionMessage(e))
  })

  if (!is.null(direct_prediction$prediction)) {
    return(direct_prediction)
  }

  extracted_error <- NULL
  extracted_prediction <- tryCatch({
    extract_asreml_multitrait_prediction(
      model = mod,
      response = response,
      gen_name = gen_name,
      kernel_ids = kernel_ids
    )
  }, error = function(e) {
    extracted_error <<- conditionMessage(e)
    NULL
  })

  if (!is.null(extracted_prediction)) {
    extracted_prediction <- asreml_mark_prediction_uncertainty_unavailable(extracted_prediction)
    return(list(
      mode = "extract",
      prediction = extracted_prediction,
      error = direct_prediction$error
    ))
  }

  list(
    mode = "predict_error",
    prediction = NULL,
    error = paste(
      c(
        direct_prediction$error,
        extracted_error
      ),
      collapse = " | "
    )
  )
}

gp_multitrait_asreml_kernel_bank <- function(gmatrix = NULL,
                                              kernel_list = NULL,
                                              phenotype_ids,
                                              context = "Multi-trait ASReml") {
  if (gp_is_kernel_list(gmatrix)) {
    kernel_list <- c(gmatrix, kernel_list %||% list())
    gmatrix <- NULL
  }
  kernels <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    kernel_list = kernel_list
  )
  # A model-ready bank may contain only canonical keys.  gp_collect_kernel_inputs
  # normally receives those through its named arguments; accept the complete
  # bank here as well so internal callers cannot silently lose its first kernel.
  if (!length(kernels) && gp_is_kernel_list(kernel_list)) {
    kernels <- kernel_list
    names(kernels) <- vapply(seq_along(kernels), function(i) {
      gp_kernel_source_name(
        kernels[[i]],
        (names(kernel_list) %||% paste0("kernel", seq_along(kernel_list)))[[i]]
      )
    }, character(1L))
    names(kernels) <- make.unique(gp_sanitize_kernel_name(names(kernels)), sep = "_")
  }
  if (!length(kernels)) {
    stop(context, " requires at least one named relationship/kernel matrix.", call. = FALSE)
  }

  phenotype_ids <- unique(as.character(phenotype_ids))
  keep_ids <- phenotype_ids
  for (i in seq_along(kernels)) {
    K <- kernels[[i]]
    if (!gp_is_kernel_matrix(K) || nrow(K) != ncol(K)) {
      stop(context, " kernel `", names(kernels)[[i]], "` must be a named numeric square matrix.", call. = FALSE)
    }
    if (!setequal(rownames(K), colnames(K))) {
      stop(context, " kernel `", names(kernels)[[i]], "` has different row and column genotype IDs.", call. = FALSE)
    }
    keep_ids <- intersect(keep_ids, intersect(rownames(K), colnames(K)))
  }
  if (!length(keep_ids)) {
    stop(context, " found no genotype IDs shared by the phenotype and every kernel.", call. = FALSE)
  }
  dropped <- setdiff(phenotype_ids, keep_ids)
  if (length(dropped)) {
    warning(
      context, " retained ", length(keep_ids), " genotype(s) shared by every kernel and dropped ",
      length(dropped), " phenotype genotype(s) absent from at least one kernel.",
      call. = FALSE
    )
  }
  kernels <- lapply(kernels, function(K) {
    K <- as.matrix(K[keep_ids, keep_ids, drop = FALSE])
    storage.mode(K) <- "double"
    K
  })
  names(kernels) <- make.unique(gp_sanitize_kernel_name(names(kernels)), sep = "_")
  list(kernels = kernels, keep_ids = keep_ids)
}

# Describes the kernels the CV folds actually fit. Feature selection rebuilds a
# single genomic kernel per fold and drops any extra kernels, so the report must
# collapse to that one kernel rather than echo the requested bank.
gp_multitrait_asreml_cv_kernel_report <- function(kernel_names,
                                                  trait_covariance = "us",
                                                  feature_active = FALSE) {
  if (isTRUE(feature_active)) {
    fitted_names <- "feature_selected_genomic_kernel"
    role <- "rebuilt_per_fold_from_selected_features"
    strategy <- paste0(
      "feature_selection_active_single_rebuilt_kernel; ",
      length(kernel_names),
      " supplied kernel(s) replaced by one fold-rebuilt genomic kernel"
    )
  } else {
    fitted_names <- as.character(kernel_names)
    role <- "independent_genetic_covariance_term"
    strategy <- paste0(
      "independent_", trait_covariance, "_trait_covariance_per_kernel_summed"
    )
  }
  list(
    configuration = gp_multitrait_kernel_configuration(fitted_names, role = role),
    parameters = gp_multi_kernel_parameter_rows(
      kernel_names = fitted_names,
      strategy = strategy
    )
  )
}

gp_multitrait_asreml_assign_kernel_inverses <- function(kernels,
                                                        prefix,
                                                        inverse,
                                                        epsilon) {
  inverse_names <- paste0(
    gp_sanitize_kernel_name(prefix), "_", seq_along(kernels), "_",
    as.integer(Sys.getpid())
  )
  for (i in seq_along(kernels)) {
    assign(
      inverse_names[[i]],
      compute_inverse_and_sparse(
        kernel = kernels[[i]],
        epsilon = epsilon,
        inverse = inverse
      ),
      envir = .GlobalEnv
    )
  }
  data.frame(
    Kernel = names(kernels),
    inverse_name = inverse_names,
    stringsAsFactors = FALSE
  )
}

gp_multitrait_asreml_random_formula <- function(index_factor,
                                                 covariance,
                                                 gen_name,
                                                 kernel_specs) {
  covariance_term <- if (grepl("\\(", covariance, fixed = FALSE)) {
    covariance
  } else {
    paste0(covariance, "(", index_factor, ")")
  }
  terms <- paste0(
    covariance_term, ":vm(", gen_name, ",",
    kernel_specs$inverse_name, ")"
  )
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

gp_multitrait_asreml_label_kernel_components <- function(vc, kernel_specs) {
  if (!is.data.frame(vc) || !nrow(vc) || is.null(kernel_specs) || !nrow(kernel_specs)) {
    return(vc)
  }
  raw_component <- as.character(vc[["Raw_component"]] %||% vc[["Component"]])
  kernel <- rep(NA_character_, length(raw_component))
  for (i in seq_len(nrow(kernel_specs))) {
    hit <- grepl(kernel_specs$inverse_name[[i]], raw_component, fixed = TRUE)
    kernel[hit] <- kernel_specs$Kernel[[i]]
  }
  vc[["Kernel"]] <- kernel
  vc
}

gp_multitrait_asreml_covariance_from_rows <- function(raw,
                                                       response,
                                                       rows,
                                                       structure = c("us", "corgh", "diag")) {
  structure <- match.arg(structure)
  component_names <- rownames(raw)
  values <- suppressWarnings(as.numeric(raw[["component"]]))
  pairs <- gp_multitrait_asreml_component_trait_pair(component_names, response)
  find_one <- function(trait_1, trait_2, correlation = FALSE) {
    candidates <- rows[vapply(rows, function(index) {
      pair <- pairs[[index]]
      same_pair <- length(pair) == 2L &&
        identical(sort(pair), sort(c(trait_1, trait_2)))
      is_cor <- grepl("\\.cor($|[^[:alnum:]])", component_names[[index]], ignore.case = TRUE)
      same_pair && identical(is_cor, correlation)
    }, logical(1L))]
    if (length(candidates) != 1L) {
      stop(
        "Expected exactly one ASReml ", if (correlation) "correlation" else "variance/covariance",
        " component for `", trait_1, "` and `", trait_2, "`; found ",
        length(candidates), ".", call. = FALSE
      )
    }
    candidates[[1L]]
  }
  out <- matrix(0, length(response), length(response), dimnames = list(response, response))
  diagonal_indices <- setNames(integer(length(response)), response)
  for (trait in response) {
    diagonal_indices[[trait]] <- find_one(trait, trait, correlation = FALSE)
    out[trait, trait] <- values[[diagonal_indices[[trait]]]]
  }
  if (!identical(structure, "diag") && length(response) > 1L) {
    for (i in seq_len(length(response) - 1L)) {
      for (j in (i + 1L):length(response)) {
        tr1 <- response[[i]]
        tr2 <- response[[j]]
        if (identical(structure, "corgh")) {
          index <- find_one(tr1, tr2, correlation = TRUE)
          value <- values[[index]] * sqrt(pmax(out[tr1, tr1], 0) * pmax(out[tr2, tr2], 0))
        } else {
          index <- find_one(tr1, tr2, correlation = FALSE)
          value <- values[[index]]
        }
        out[tr1, tr2] <- out[tr2, tr1] <- value
      }
    }
  }
  out
}

gp_multitrait_asreml_gaussian_model <- function(pheno_object,
                                                response,
                                                gen_name,
                                                gmatrix = NULL,
                                                kernel_list = NULL,
                                                inverse = TRUE,
                                                epsilon = 1e-6,
                                                engine = "asreml",
                                                workspace = 1e08,
                                                pworkspace = NULL,
                                                maxit = 50,
                                                ai_sing = TRUE,
                                                trait_covariance = c("us", "corgh", "diag"),
                                                extract_variance = TRUE,
                                                max_stability_updates = 4L,
                                                require_variance_stability = extract_variance) {
  if (!engine %in% rownames(installed.packages())) {
    stop("You need to install asreml-R to use multi-trait ASReml-R.", call. = FALSE)
  }
  trait_covariance <- match.arg(trait_covariance)

  pheno_object <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  long_pheno <- gp_multitrait_asreml_long_data(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name
  )

  gids <- as.character(unique(pheno_object[[gen_name]]))
  kernel_bundle <- gp_multitrait_asreml_kernel_bank(
    gmatrix = gmatrix,
    kernel_list = kernel_list,
    phenotype_ids = gids,
    context = "Multi-trait ASReml"
  )
  keep_ids <- kernel_bundle$keep_ids

  pheno_object <- pheno_object[match(keep_ids, as.character(pheno_object[[gen_name]])), , drop = FALSE]
  long_pheno <- long_pheno[long_pheno[[gen_name]] %in% keep_ids, , drop = FALSE]
  long_pheno[[gen_name]] <- factor(long_pheno[[gen_name]], levels = keep_ids)
  long_pheno$Trait <- factor(long_pheno$Trait, levels = response)

  kernel_specs <- gp_multitrait_asreml_assign_kernel_inverses(
    kernels = kernel_bundle$kernels,
    prefix = "pp_mt_inv",
    inverse = inverse,
    epsilon = epsilon
  )
  on.exit({
    for (inverse_name in kernel_specs$inverse_name) remove_from_global(inverse_name)
  }, add = TRUE)

  asreml_fn <- getFromNamespace("asreml", engine)
  predict_fn <- getOption("PredictProR.asreml_predict_impl", getFromNamespace("predict.asreml", engine))
  update_fn <- getFromNamespace("update.asreml", engine)
  asreml_options_fn <- getFromNamespace("asreml.options", engine)

  if (isTRUE(ai_sing)) {
    old_options <- tryCatch(asreml_options_fn(), error = function(e) NULL)
    old_ai_sing <- if (is.list(old_options) && "ai.sing" %in% names(old_options)) old_options[["ai.sing"]] else NULL
    try(asreml_options_fn(ai.sing = TRUE), silent = TRUE)
    on.exit({
      if (!is.null(old_ai_sing)) {
        try(asreml_options_fn(ai.sing = old_ai_sing), silent = TRUE)
      }
    }, add = TRUE)
  }

  fit_args <- list(
    fixed = stats::as.formula("Observed_value ~ Trait"),
    random = gp_multitrait_asreml_random_formula(
      index_factor = "Trait",
      covariance = trait_covariance,
      gen_name = gen_name,
      kernel_specs = kernel_specs
    ),
    residual = stats::as.formula(paste0("~ units:", trait_covariance, "(Trait)")),
    data = long_pheno,
    asmv = "Trait",
    na.action = list(x = "include", y = "include"),
    workspace = workspace,
    maxit = maxit,
    ai.sing = isTRUE(ai_sing)
  )

  mod <- tryCatch(
    do.call(asreml_fn, fit_args),
    error = function(e) {
      stop(paste("Multi-trait ASReml-R fit failed:", conditionMessage(e)), call. = FALSE)
    }
  )

  convergence_trace <- NULL
  if (nrow(kernel_specs) > 1L) {
    stabilized <- gp_multitrait_asreml_fit_until_stable(
      model = mod,
      update_fn = update_fn,
      max_updates = max_stability_updates,
      max_percent_change = if (isTRUE(require_variance_stability)) 1 else Inf
    )
    mod <- stabilized$model
    convergence_trace <- stabilized$trace
    if (!isTRUE(tail(stabilized$trace$converged, 1L))) {
      stop(
        "Multi-trait multi-kernel ASReml GBLUP did not converge; predictions and variance output were not returned.",
        call. = FALSE
      )
    }
    if (isTRUE(require_variance_stability) && !isTRUE(stabilized$stable)) {
      final <- tail(stabilized$trace, 1L)
      stop(
        "Multi-trait multi-kernel ASReml GBLUP did not reach a variance-stable endpoint after ",
        max_stability_updates, " continuation fits (maximum component change = ",
        signif(final$maximum_percent_change, 6), "%). Variance output was not promoted.",
        call. = FALSE
      )
    }
  } else {
    for (i in seq_len(3L)) {
      if (!isTRUE(mod$converge)) mod <- update_fn(mod) else break
    }
    gp_asreml_warn_if_not_converged(mod, context = "multi-trait GBLUP")
  }

  prediction_result <- gp_multitrait_asreml_predict_or_extract(
    mod = mod,
    predict_fn = predict_fn,
    classify = paste("Trait", gen_name, sep = ":"),
    response = response,
    gen_name = gen_name,
    kernel_ids = kernel_specs$inverse_name,
    pworkspace = pworkspace
  )
  if (is.null(prediction_result$prediction)) {
    stop(paste("Multi-trait ASReml-R prediction failed:", prediction_result$error), call. = FALSE)
  }
  pred_df <- prediction_result$prediction

  obs_long <- gp_multitrait_asreml_long_data(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name
  )
  pred_long <- merge(
    obs_long[, c(gen_name, "Trait", "Observed_value"), drop = FALSE],
    pred_df,
    by = c(gen_name, "Trait"),
    all.x = TRUE,
    sort = FALSE
  )
  pred_long$Train_Test_Label <- ifelse(is.na(pred_long$Observed_value), "Test", "Train")
  pred_long <- pred_long[order(match(pred_long[[gen_name]], keep_ids), match(pred_long$Trait, response)), , drop = FALSE]
  rownames(pred_long) <- NULL
  vc <- if (isTRUE(extract_variance)) {
    gp_multitrait_asreml_variance_components(
      mod,
      response = response,
      kernel_specs = kernel_specs
    )
  } else {
    data.frame()
  }
  cov_outputs <- if (isTRUE(extract_variance)) {
    gp_multitrait_asreml_covariance_outputs(
      mod,
      response = response,
      kernel_specs = kernel_specs,
      trait_covariance = trait_covariance
    )
  } else {
    list()
  }

  c(list(
    Asreml_model = mod,
    Predicted_value = pred_long,
    predicted_values = pred_long,
    multitrait_prediction_wide = gp_multitrait_asreml_prediction_wide(pred_long, gen_name = gen_name),
    multitrait_trait_counts = gp_multitrait_asreml_trait_counts(pred_long, response),
    variance_components = vc,
    diagnostic_plots = gp_multitrait_asreml_diagnostic_plot(pred_long),
    asreml_trait_covariance = trait_covariance,
    asreml_ai_sing = isTRUE(ai_sing),
    asreml_prediction_mode = prediction_result$mode,
    asreml_convergence_trace = convergence_trace,
    kernel_configuration = gp_multitrait_kernel_configuration(kernel_bundle$kernels),
    model_parameters = gp_multi_kernel_parameter_rows(
      kernel_names = kernel_specs$Kernel,
      strategy = paste0("independent_", trait_covariance, "_trait_covariance_per_kernel_summed")
    )
  ), cov_outputs)
}

# Multi-trait ASReml cross-validation with the same predict-or-extract
# fallback the true-prediction path uses. Per fold: mask y[tst, ] = NA for
# all trait columns jointly, refit the multi-trait model, and let
# gp_multitrait_asreml_predict_or_extract() handle the predict step (it
# falls through to extract_asreml_multitrait_prediction() when
# asreml::predict() fails -- e.g. workspace overflow, classify mismatch).
# Returns the same cv_results / cv_results_processed shape used by
# gp_multitrait_ml_gaussian_cv() so the route layer can reuse the
# existing aggregator and downstream summary plumbing.
#
# Scope note: gp_multitrait_asreml_gaussian_model() currently requires one
# row per genotype (single-environment multi-trait), so this CV is the
# single-env case. Multi-trait MET on the ASReml side would need the model
# layer to lift that constraint first; this CV would then extend naturally.
gp_multitrait_asreml_gaussian_cv <- function(pheno_object,
                                             response,
                                             gen_name,
                                             gmatrix = NULL,
                                             kernel_list = NULL,
                                             cross_validation_meth = "K-Folds",
                                             nfolds = 5L,
                                             replication = 1L,
                                             eval_metrics = c("root_mean_squared_error",
                                                              "mean_absolute_error",
                                                              "kendalls_tau"),
                                             trait_covariance = c("us", "corgh", "diag"),
                                             inverse = TRUE,
                                             epsilon = 1e-6,
                                             engine = "asreml",
                                             workspace = 1e08,
                                             pworkspace = NULL,
                                             maxit = 50,
                                             ai_sing = TRUE,
                                             random_seed = 123L,
                                             feature_predictor_data = NULL,
                                             feature_score_metadata = NULL,
                                             feature_k = NULL,
                                             feature_scoring_cv = "fixed",
                                             feature_scoring_model = "Ridge_Regression",
                                             feature_source_map = NULL,
                                             gmatrix_method = NULL,
                                             ploidy = "auto",
                                             feature_ridge_lambda = 1,
                                             feature_bayes_nIter = 1500L,
                                             feature_bayes_burnIn = 500L,
                                             feature_bayes_thin = 5L,
                                             ntree = 500L,
                                             mtry = NULL,
                                             nodesize = NULL,
                                             rf_n_jobs = 1L) {
  if (!engine %in% rownames(utils::installed.packages())) {
    stop("You need to install asreml-R to use multi-trait ASReml-R CV.", call. = FALSE)
  }
  trait_covariance <- match.arg(trait_covariance)
  cv_token <- normalize_cv_token(cross_validation_meth)
  if (!cv_token %in% c("k_folds", "repeated_k_folds")) {
    stop("Multi-trait ASReml CV currently supports K-Folds only.", call. = FALSE)
  }

  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  ids <- as.character(ph[[gen_name]])
  if (anyDuplicated(ids) > 0L) {
    stop("Multi-trait ASReml CV currently requires one phenotype row per genotype.", call. = FALSE)
  }
  kernel_bundle <- gp_multitrait_asreml_kernel_bank(
    gmatrix = gmatrix,
    kernel_list = kernel_list,
    phenotype_ids = ids,
    context = "Multi-trait ASReml CV"
  )
  keep_ids <- kernel_bundle$keep_ids
  ph <- ph[match(keep_ids, ids), c(gen_name, response), drop = FALSE]
  rownames(ph) <- NULL
  ids <- as.character(ph[[gen_name]])
  n <- nrow(ph)
  feature_active <- !is.null(feature_score_metadata) && !is.null(feature_k)
  feature_scoring_cv <- gp_feature_normalize_cv_policy(feature_scoring_cv, allow_both = FALSE)
  if (isTRUE(feature_active)) {
    if (is.null(feature_predictor_data)) {
      stop("Multi-trait ASReml feature selection requires raw marker predictors.", call. = FALSE)
    }
    missing_feature_ids <- setdiff(ids, rownames(feature_predictor_data))
    if (length(missing_feature_ids)) {
      stop("Multi-trait ASReml marker predictors are missing phenotype IDs.", call. = FALSE)
    }
    feature_predictor_data <- feature_predictor_data[ids, , drop = FALSE]
  }

  cv_results <- vector("list", length = replication)
  for (rep_idx in seq_len(replication)) {
    pred_blocks <- list()
    feature_records <- list()
    block_idx <- 1L
    folds_rel <- gp_tuning_folds(
      y = seq_len(n),
      response_family = "gaussian",
      nfolds = min(nfolds, n),
      random_state = as.integer((random_seed %||% 123L) + rep_idx - 1L)
    )
    for (fold_idx in seq_along(folds_rel)) {
      test_idx <- sort(unique(as.integer(folds_rel[[fold_idx]])))
      train_idx <- setdiff(seq_len(n), test_idx)
      if (!length(train_idx) || !length(test_idx)) next
      gmatrix_fold <- kernel_bundle$kernels[[1L]]
      kernel_list_fold <- if (length(kernel_bundle$kernels) > 1L) {
        kernel_bundle$kernels[-1L]
      } else {
        NULL
      }
      if (isTRUE(feature_active)) {
        feature_view <- gp_feature_multitrait_view(
          predictor_data = feature_predictor_data,
          pheno_data = ph,
          response = response,
          gen_name = gen_name,
          k = feature_k,
          selection_mode = feature_scoring_cv,
          feature_score_metadata = feature_score_metadata,
          test_rows = test_idx,
          scoring_model = feature_scoring_model,
          seed = as.integer((random_seed %||% 123L) + rep_idx * 10000L + fold_idx),
          source_block = feature_source_map,
          response_family = "gaussian",
          replication = rep_idx,
          fold = fold_idx,
          ridge_lambda = feature_ridge_lambda,
          bayes_nIter = feature_bayes_nIter,
          bayes_burnIn = feature_bayes_burnIn,
          bayes_thin = feature_bayes_thin,
          ntree = ntree,
          mtry = mtry,
          nodesize = nodesize,
          rf_n_jobs = rf_n_jobs
        )
        gmatrix_fold <- gp_feature_genomic_kernel_from_view(
          predictor_data = feature_view$predictor_data,
          source_map = feature_source_map,
          gmatrix_method = gmatrix_method,
          ploidy = ploidy,
          context = "Multi-trait ASReml CV"
        )
        kernel_list_fold <- NULL
        feature_records[[length(feature_records) + 1L]] <- feature_view$metadata
      }
      ph_masked <- ph
      ph_masked[test_idx, response] <- NA_real_

      fold_fit <- tryCatch(
        gp_multitrait_asreml_gaussian_model(
          pheno_object = ph_masked,
          response = response,
          gen_name = gen_name,
          gmatrix = gmatrix_fold,
          kernel_list = kernel_list_fold,
          inverse = inverse,
          epsilon = epsilon,
          engine = engine,
          workspace = workspace,
          pworkspace = pworkspace,
          maxit = maxit,
          ai_sing = ai_sing,
          trait_covariance = trait_covariance,
          extract_variance = FALSE
        ),
        error = function(e) {
          warning(
            "Multi-trait ASReml CV fold ", fold_idx, " (rep ", rep_idx,
            ") fit failed; predictions for this fold left as NA: ",
            conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )

      if (is.null(fold_fit)) {
        for (tr in response) {
          pred_blocks[[block_idx]] <- data.frame(
            stringsAsFactors = FALSE,
            id_value = ids[test_idx],
            Trait = tr,
            model = "Multi_trait_ASReml",
            Observed_value = ph[[tr]][test_idx],
            Predicted_value = NA_real_,
            Train_Test_Label = "Test",
            fold = fold_idx,
            rep = rep_idx,
            cv_role = "test"
          )
          names(pred_blocks[[block_idx]])[1] <- gen_name
          if (isTRUE(feature_active)) {
            pred_blocks[[block_idx]][["feature_k"]] <- as.integer(feature_k)
            pred_blocks[[block_idx]][["feature_scoring_model"]] <- feature_scoring_model
            pred_blocks[[block_idx]][["feature_scoring_cv"]] <- feature_scoring_cv
          }
          block_idx <- block_idx + 1L
        }
        next
      }

      pred_long_fold <- fold_fit$predicted_values
      pred_long_fold <- pred_long_fold[
        as.character(pred_long_fold[[gen_name]]) %in% ids[test_idx], ,
        drop = FALSE]
      for (tr in response) {
        rows <- pred_long_fold[as.character(pred_long_fold$Trait) == tr, , drop = FALSE]
        match_idx <- match(ids[test_idx], as.character(rows[[gen_name]]))
        pred_vals <- as.numeric(rows$Predicted_value[match_idx])
        pred_blocks[[block_idx]] <- data.frame(
          stringsAsFactors = FALSE,
          id_value = ids[test_idx],
          Trait = tr,
          model = "Multi_trait_ASReml",
          Observed_value = ph[[tr]][test_idx],
          Predicted_value = pred_vals,
          Train_Test_Label = "Test",
          fold = fold_idx,
          rep = rep_idx,
          cv_role = "test"
        )
        names(pred_blocks[[block_idx]])[1] <- gen_name
        if (isTRUE(feature_active)) {
          pred_blocks[[block_idx]][["feature_k"]] <- as.integer(feature_k)
          pred_blocks[[block_idx]][["feature_scoring_model"]] <- feature_scoring_model
          pred_blocks[[block_idx]][["feature_scoring_cv"]] <- feature_scoring_cv
        }
        block_idx <- block_idx + 1L
      }
    }
    pred_long <- if (length(pred_blocks)) do.call(rbind, pred_blocks) else data.frame()
    metric_table <- gp_multitrait_ml_metric_table(
      pred_long = pred_long,
      eval_metrics = eval_metrics,
      model_type = "Multi_trait_ASReml",
      rep = rep_idx
    )
    if (isTRUE(feature_active) && nrow(metric_table)) {
      metric_table[["feature_k"]] <- as.integer(feature_k)
      metric_table[["feature_scoring_model"]] <- feature_scoring_model
      metric_table[["feature_scoring_cv"]] <- feature_scoring_cv
    }
    cv_results[[rep_idx]] <- list(
      trait = paste(response, collapse = ","),
      model = "Multi_trait_ASReml",
      rep = rep_idx,
      eval_metrics_reps = metric_table,
      ypred_cv_Reps_all = pred_long,
      yprob_cv_Reps_all = NULL,
      feature_selection_metadata = gp_feature_bind_metadata(feature_records)
    )
  }

  processed <- gp_multitrait_ml_cv_process(cv_results)
  processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(cv_results, `[[`, "feature_selection_metadata")
  )
  kernel_report <- gp_multitrait_asreml_cv_kernel_report(
    kernel_names = names(kernel_bundle$kernels),
    trait_covariance = trait_covariance,
    feature_active = feature_active
  )
  processed$kernel_configuration <- kernel_report$configuration
  processed$model_parameters <- kernel_report$parameters
  list(cv_results = cv_results, cv_results_processed = processed)
}

# Grouped multi-trait, multi-environment ASReml GBLUP -----------------------

gp_multitrait_asreml_mtmet_structure <- function(var_cov_str = NULL,
                                                  n_groups) {
  n_groups <- as.integer(n_groups)[1L]
  if (!is.finite(n_groups) || n_groups < 2L) {
    stop("Grouped MT-MET ASReml requires at least two trait-environment groups.", call. = FALSE)
  }
  token <- tolower(gsub("\\s+", "", as.character(var_cov_str %||% "")[[1L]]))
  if (!nzchar(token)) {
    token <- paste0("fa", min(3L, n_groups - 1L))
  }
  rank_text <- if (grepl("^fa[1-9][0-9]*$", token)) {
    sub("^fa", "", token)
  } else if (grepl("^fa\\([1-9][0-9]*\\)$", token)) {
    sub("^fa\\(([1-9][0-9]*)\\)$", "\\1", token)
  } else {
    stop(
      "Grouped MT-MET ASReml currently uses the validated factor-analytic ",
      "genetic structure. Set var_cov_str to fa1, fa2, fa3, ...; ",
      "for example var_cov_str = 'fa3'.",
      call. = FALSE
    )
  }
  rank <- as.integer(rank_text)
  if (!is.finite(rank) || rank < 1L || rank >= n_groups) {
    stop(
      "The MT-MET ASReml FA rank must be at least 1 and smaller than the ",
      "number of trait-environment groups (", n_groups, ").",
      call. = FALSE
    )
  }
  list(token = paste0("fa", rank), rank = rank)
}

gp_multitrait_asreml_mtmet_long_data <- function(pheno_object,
                                                  response,
                                                  gen_name,
                                                  heter_groups) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(gen_name, heter_groups, response)
  missing <- setdiff(required, names(ph))
  if (length(missing)) {
    stop(
      "Grouped MT-MET ASReml data are missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  response <- unique(as.character(response))
  if (length(response) < 2L) {
    stop("Grouped MT-MET ASReml requires at least two response columns.", call. = FALSE)
  }
  keep <- !is.na(ph[[gen_name]]) & nzchar(as.character(ph[[gen_name]])) &
    !is.na(ph[[heter_groups]]) & nzchar(as.character(ph[[heter_groups]]))
  ph <- ph[keep, required, drop = FALSE]
  if (!nrow(ph)) {
    stop("No rows have both genotype and environment identifiers.", call. = FALSE)
  }
  env_levels <- unique(as.character(ph[[heter_groups]]))
  genotype_levels <- unique(as.character(ph[[gen_name]]))
  if (length(env_levels) < 2L) {
    stop("Grouped MT-MET ASReml requires at least two environments.", call. = FALSE)
  }

  ph[[gen_name]] <- as.character(ph[[gen_name]])
  ph[[heter_groups]] <- as.character(ph[[heter_groups]])
  ph$PPUnit <- sprintf("U%0*d", nchar(nrow(ph)), seq_len(nrow(ph)))
  names(ph)[names(ph) == heter_groups] <- "PPEnv"

  long <- stats::reshape(
    ph,
    varying = response,
    v.names = "Observed_value",
    timevar = "PPTrait",
    times = response,
    idvar = c(gen_name, "PPEnv", "PPUnit"),
    direction = "long"
  )
  rownames(long) <- NULL
  long[[gen_name]] <- factor(as.character(long[[gen_name]]), levels = genotype_levels)
  long$PPEnv <- factor(as.character(long$PPEnv), levels = env_levels)
  long$PPTrait <- factor(as.character(long$PPTrait), levels = response)
  long$PPUnit <- factor(long$PPUnit)

  group_lookup <- expand.grid(
    Trait = response,
    Env = env_levels,
    stringsAsFactors = FALSE
  )
  width <- max(3L, nchar(nrow(group_lookup)))
  group_lookup$PPGroup <- sprintf(paste0("G%0", width, "d"), seq_len(nrow(group_lookup)))
  lookup_key <- paste(group_lookup$Trait, group_lookup$Env, sep = "\r")
  row_key <- paste(as.character(long$PPTrait), as.character(long$PPEnv), sep = "\r")
  long$PPGroup <- factor(
    group_lookup$PPGroup[match(row_key, lookup_key)],
    levels = group_lookup$PPGroup
  )
  if (anyNA(long$PPGroup)) {
    stop("Could not map every long phenotype row to one trait-environment group.", call. = FALSE)
  }
  long <- long[order(long$PPEnv, long$PPUnit, long$PPTrait), , drop = FALSE]
  rownames(long) <- NULL

  list(
    data = long,
    group_lookup = group_lookup[, c("PPGroup", "Trait", "Env"), drop = FALSE],
    genotype_levels = genotype_levels,
    environment_levels = env_levels,
    response = response
  )
}

gp_multitrait_asreml_endpoint_diagnostics <- function(model,
                                                      summary_fn = summary) {
  vc <- summary_fn(model)$varcomp
  change <- suppressWarnings(as.numeric(vc[["%ch"]]))
  change <- change[is.finite(change)]
  data.frame(
    converged = isTRUE(model$converge),
    log_likelihood = suppressWarnings(as.numeric(model$loglik)[1L]),
    maximum_percent_change = if (length(change)) max(change) else NA_real_,
    components_above_one_percent = if (length(change)) sum(change > 1) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

gp_multitrait_asreml_fit_until_stable <- function(model,
                                                  update_fn,
                                                  max_updates = 4L,
                                                  max_percent_change = 1) {
  trace <- list()
  for (i in 0:as.integer(max_updates)) {
    diagnostic <- gp_multitrait_asreml_endpoint_diagnostics(model)
    diagnostic$update_number <- i
    trace[[length(trace) + 1L]] <- diagnostic
    stable <- isTRUE(diagnostic$converged) &&
      is.finite(diagnostic$maximum_percent_change) &&
      diagnostic$maximum_percent_change <= max_percent_change
    if (stable || i == as.integer(max_updates)) {
      break
    }
    model <- update_fn(model)
  }
  trace <- do.call(rbind, trace)
  rownames(trace) <- NULL
  list(
    model = model,
    trace = trace,
    stable = isTRUE(tail(trace$converged, 1L)) &&
      is.finite(tail(trace$maximum_percent_change, 1L)) &&
      tail(trace$maximum_percent_change, 1L) <= max_percent_change
  )
}

gp_multitrait_asreml_vpredict <- function(model, rhs, vpredict_fn) {
  value <- vpredict_fn(
    model,
    stats::as.formula(paste("derived ~", rhs))
  )
  value <- as.data.frame(value, stringsAsFactors = FALSE, check.names = FALSE)
  estimate_col <- gp_contract_first_name(names(value), c("Estimate", "estimate"))
  se_col <- gp_contract_first_name(names(value), c("SE", "se", "std.error", "Standard_error"))
  if (is.na(estimate_col) || !nrow(value)) {
    stop("ASReml vpredict did not return a derived estimate.", call. = FALSE)
  }
  c(
    estimate = as.numeric(value[[estimate_col]][1L]),
    standard_error = if (!is.na(se_col)) as.numeric(value[[se_col]][1L]) else NA_real_
  )
}

gp_multitrait_asreml_mtmet_variance <- function(model,
                                                 group_lookup,
                                                 fa_rank,
                                                 vpredict_fn = NULL,
                                                 summary_fn = summary,
                                                 tolerance = 1e-8) {
  if (is.null(vpredict_fn)) {
    vpredict_fn <- getFromNamespace("vpredict", "asreml")
  }
  raw <- summary_fn(model)$varcomp
  component_names <- rownames(raw)
  if (is.null(component_names) || !length(component_names)) {
    stop("ASReml did not return named variance components.", call. = FALSE)
  }
  require_one <- function(index, description) {
    if (length(index) != 1L) {
      stop(
        "Expected exactly one ", description, "; found ", length(index), ".",
        call. = FALSE
      )
    }
    index
  }
  genetic_rows <- which(grepl("fa\\(PPGroup", component_names) & grepl("vm\\(", component_names))
  parameter_types <- c("var", paste0("fa", seq_len(fa_rank)))
  parameter_map <- do.call(rbind, lapply(group_lookup$PPGroup, function(group_value) {
    do.call(rbind, lapply(parameter_types, function(parameter_type) {
      index <- require_one(
        genetic_rows[endsWith(component_names[genetic_rows], paste0("!", group_value, "!", parameter_type))],
        paste0("genetic ", parameter_type, " for group `", group_value, "`")
      )
      data.frame(
        PPGroup = group_value,
        parameter = parameter_type,
        Component = component_names[index],
        component_index = index,
        Components = as.numeric(raw[index, "component"]),
        Standard_error = as.numeric(raw[index, "std.error"]),
        Bound = as.character(raw[index, "bound"]),
        Percent_change = as.numeric(raw[index, "%ch"]),
        stringsAsFactors = FALSE
      )
    }))
  }))

  groups <- as.character(group_lookup$PPGroup)
  indices <- matrix(
    NA_integer_,
    nrow = length(groups),
    ncol = length(parameter_types),
    dimnames = list(groups, parameter_types)
  )
  estimates <- indices
  storage.mode(estimates) <- "double"
  for (group_value in groups) {
    map <- parameter_map[parameter_map$PPGroup == group_value, , drop = FALSE]
    indices[group_value, map$parameter] <- map$component_index
    estimates[group_value, map$parameter] <- map$Components
  }
  if (any(estimates[, "var"] < -sqrt(.Machine$double.eps))) {
    stop("ASReml returned a negative FA specific variance.", call. = FALSE)
  }

  loadings <- estimates[, paste0("fa", seq_len(fa_rank)), drop = FALSE]
  genetic_covariance <- tcrossprod(loadings) + diag(estimates[, "var"], nrow = length(groups))
  dimnames(genetic_covariance) <- list(groups, groups)
  if (any(!is.finite(diag(genetic_covariance))) ||
      any(diag(genetic_covariance) <= tolerance)) {
    stop(
      "ASReml did not return strictly positive finite genetic variances for every trait-environment group; variance output was not promoted.",
      call. = FALSE
    )
  }
  genetic_covariance_se <- genetic_covariance
  genetic_correlation <- stats::cov2cor(genetic_covariance)
  genetic_correlation_se <- genetic_covariance
  genetic_long <- list()
  long_index <- 1L
  genetic_variance_rhs <- setNames(character(length(groups)), groups)
  genetic_covariance_rhs <- matrix(
    NA_character_, length(groups), length(groups), dimnames = list(groups, groups)
  )
  for (group_value in groups) {
    genetic_variance_rhs[group_value] <- paste(
      c(
        sprintf("V%d", indices[group_value, "var"]),
        sprintf("V%d^2", indices[group_value, paste0("fa", seq_len(fa_rank))])
      ),
      collapse = " + "
    )
  }

  for (i in seq_along(groups)) {
    for (j in i:length(groups)) {
      group_1 <- groups[i]
      group_2 <- groups[j]
      if (i == j) {
        terms <- c(
          sprintf("V%d", indices[group_1, "var"]),
          sprintf("V%d^2", indices[group_1, paste0("fa", seq_len(fa_rank))])
        )
        rhs <- paste(terms, collapse = " + ")
      } else {
        terms <- sprintf(
          "V%d*V%d",
          indices[group_1, paste0("fa", seq_len(fa_rank))],
          indices[group_2, paste0("fa", seq_len(fa_rank))]
        )
        rhs <- paste(terms, collapse = " + ")
      }
      genetic_covariance_rhs[group_1, group_2] <- rhs
      genetic_covariance_rhs[group_2, group_1] <- rhs
      native <- gp_multitrait_asreml_vpredict(model, rhs, vpredict_fn)
      manual <- genetic_covariance[group_1, group_2]
      if (!is.finite(native[["estimate"]]) ||
          abs(native[["estimate"]] - manual) > tolerance * (1 + abs(manual))) {
        stop(
          "Named FA covariance reconstruction disagrees with ASReml vpredict for ",
          group_1, " and ", group_2, ".",
          call. = FALSE
        )
      }
      genetic_covariance_se[group_1, group_2] <- native[["standard_error"]]
      genetic_covariance_se[group_2, group_1] <- native[["standard_error"]]
      correlation_native <- if (i == j) {
        c(estimate = 1, standard_error = 0)
      } else {
        gp_multitrait_asreml_vpredict(
          model,
          paste0(
            "(", rhs, ")/sqrt((", genetic_variance_rhs[group_1], ")*(",
            genetic_variance_rhs[group_2], "))"
          ),
          vpredict_fn
        )
      }
      genetic_correlation_se[group_1, group_2] <- correlation_native[["standard_error"]]
      genetic_correlation_se[group_2, group_1] <- correlation_native[["standard_error"]]
      row_1 <- group_lookup[match(group_1, group_lookup$PPGroup), , drop = FALSE]
      row_2 <- group_lookup[match(group_2, group_lookup$PPGroup), , drop = FALSE]
      genetic_long[[long_index]] <- data.frame(
        Trait_1 = row_1$Trait,
        Env_1 = row_1$Env,
        Trait_2 = row_2$Trait,
        Env_2 = row_2$Env,
        Covariance = native[["estimate"]],
        Standard_error = native[["standard_error"]],
        Correlation = correlation_native[["estimate"]],
        Correlation_standard_error = correlation_native[["standard_error"]],
        stringsAsFactors = FALSE
      )
      long_index <- long_index + 1L
    }
  }
  genetic_eigen <- eigen(genetic_covariance, symmetric = TRUE, only.values = TRUE)$values
  if (min(genetic_eigen) < -tolerance) {
    stop("The reconstructed ASReml genetic covariance matrix is not positive semidefinite.", call. = FALSE)
  }

  residual_covariance <- list()
  residual_covariance_se <- list()
  residual_correlation <- list()
  residual_correlation_se <- list()
  residual_long <- list()
  residual_index <- 1L
  public_rows <- list()
  public_index <- 1L
  residual_variance_rhs <- setNames(character(length(groups)), groups)
  traits <- unique(as.character(group_lookup$Trait))
  environments <- unique(as.character(group_lookup$Env))

  for (environment in environments) {
    variance_index <- setNames(integer(length(traits)), traits)
    variance_value <- setNames(numeric(length(traits)), traits)
    for (trait in traits) {
      component_name <- paste0("PPEnv_", environment, "!PPTrait_", trait)
      index <- require_one(which(component_names == component_name), paste0("residual variance `", component_name, "`"))
      variance_index[trait] <- index
      variance_value[trait] <- as.numeric(raw[index, "component"])
    }
    if (any(variance_value < -sqrt(.Machine$double.eps))) {
      stop("ASReml returned a negative residual variance in environment `", environment, "`.", call. = FALSE)
    }
    correlation_matrix <- diag(1, nrow = length(traits))
    dimnames(correlation_matrix) <- list(traits, traits)
    correlation_index <- matrix(NA_integer_, length(traits), length(traits), dimnames = list(traits, traits))
    if (length(traits) > 1L) {
      for (i in seq_len(length(traits) - 1L)) {
        for (j in (i + 1L):length(traits)) {
          candidates <- c(
            paste0("PPEnv_", environment, "!PPTrait!", traits[j], ":!PPTrait!", traits[i], ".cor"),
            paste0("PPEnv_", environment, "!PPTrait!", traits[i], ":!PPTrait!", traits[j], ".cor")
          )
          found <- candidates[candidates %in% component_names]
          index <- require_one(match(found, component_names, nomatch = 0L), paste0("residual correlation for `", traits[i], "` and `", traits[j], "` in `", environment, "`"))
          correlation_matrix[i, j] <- correlation_matrix[j, i] <- as.numeric(raw[index, "component"])
          correlation_index[i, j] <- correlation_index[j, i] <- index
        }
      }
    }
    residual_sd <- sqrt(pmax(variance_value, 0))
    covariance_matrix <- diag(residual_sd) %*% correlation_matrix %*% diag(residual_sd)
    dimnames(covariance_matrix) <- list(traits, traits)
    covariance_se_matrix <- covariance_matrix
    correlation_se_matrix <- covariance_matrix
    diag(correlation_se_matrix) <- 0

    for (i in seq_along(traits)) {
      for (j in i:length(traits)) {
        if (i == j) {
          cov_rhs <- sprintf("V%d", variance_index[traits[i]])
          cor_native <- c(estimate = 1, standard_error = 0)
        } else {
          cov_rhs <- sprintf(
            "V%d*sqrt(V%d*V%d)",
            correlation_index[i, j], variance_index[traits[i]], variance_index[traits[j]]
          )
          cor_native <- gp_multitrait_asreml_vpredict(
            model,
            sprintf("V%d", correlation_index[i, j]),
            vpredict_fn
          )
        }
        cov_native <- gp_multitrait_asreml_vpredict(model, cov_rhs, vpredict_fn)
        manual <- covariance_matrix[i, j]
        if (!is.finite(cov_native[["estimate"]]) ||
            abs(cov_native[["estimate"]] - manual) > tolerance * (1 + abs(manual))) {
          stop("Named residual covariance reconstruction disagrees with ASReml vpredict.", call. = FALSE)
        }
        covariance_se_matrix[i, j] <- covariance_se_matrix[j, i] <- cov_native[["standard_error"]]
        correlation_se_matrix[i, j] <- correlation_se_matrix[j, i] <- cor_native[["standard_error"]]
        residual_long[[residual_index]] <- data.frame(
          Env = environment,
          Trait_1 = traits[i],
          Trait_2 = traits[j],
          Covariance = cov_native[["estimate"]],
          Standard_error = cov_native[["standard_error"]],
          Correlation = cor_native[["estimate"]],
          Correlation_standard_error = cor_native[["standard_error"]],
          stringsAsFactors = FALSE
        )
        residual_index <- residual_index + 1L
      }
    }
    if (min(eigen(covariance_matrix, symmetric = TRUE, only.values = TRUE)$values) < -tolerance) {
      stop("The reconstructed residual covariance matrix is not positive semidefinite in `", environment, "`.", call. = FALSE)
    }
    if (any(!is.finite(diag(covariance_matrix))) ||
        any(diag(covariance_matrix) <= tolerance)) {
      stop(
        "ASReml did not return strictly positive finite residual variances in `",
        environment,
        "`; variance output was not promoted.",
        call. = FALSE
      )
    }
    residual_covariance[[environment]] <- covariance_matrix
    residual_covariance_se[[environment]] <- covariance_se_matrix
    residual_correlation[[environment]] <- correlation_matrix
    residual_correlation_se[[environment]] <- correlation_se_matrix
  }

  endpoint <- gp_multitrait_asreml_endpoint_diagnostics(model, summary_fn = summary_fn)
  for (i in seq_len(nrow(group_lookup))) {
    group_value <- group_lookup$PPGroup[i]
    trait <- group_lookup$Trait[i]
    environment <- group_lookup$Env[i]
    residual_index_value <- which(component_names == paste0("PPEnv_", environment, "!PPTrait_", trait))
    vg_rhs <- genetic_variance_rhs[group_value]
    ve_rhs <- sprintf("V%d", residual_index_value)
    residual_variance_rhs[group_value] <- ve_rhs
    vg <- gp_multitrait_asreml_vpredict(model, vg_rhs, vpredict_fn)
    ve <- gp_multitrait_asreml_vpredict(model, ve_rhs, vpredict_fn)
    h2 <- gp_multitrait_asreml_vpredict(
      model,
      paste0("(", vg_rhs, ")/((", vg_rhs, ")+", ve_rhs, ")"),
      vpredict_fn
    )
    public_rows[[public_index]] <- data.frame(
      Trait = trait,
      Env = environment,
      Component = c("genetic_variance", "residual_variance", "heritability"),
      Components = c(vg[["estimate"]], ve[["estimate"]], h2[["estimate"]]),
      Standard_error = c(vg[["standard_error"]], ve[["standard_error"]], h2[["standard_error"]]),
      Estimation_method = "ASReml_REML_vpredict",
      Variance_estimand = c(
        "additive_genetic_variance_for_trait_environment_group",
        "observation_level_residual_variance_within_environment",
        "observation_level_narrow_sense_h2=Vg/(Vg+Ve)"
      ),
      Fit_converged = endpoint$converged,
      Endpoint_max_percent_change = endpoint$maximum_percent_change,
      stringsAsFactors = FALSE
    )
    public_index <- public_index + 1L
  }

  public_variance <- do.call(rbind, public_rows)
  if (any(!is.finite(public_variance$Components)) ||
      any(!is.finite(public_variance$Standard_error))) {
    stop(
      "ASReml vpredict returned a non-finite variance, heritability, or standard error; variance output was not promoted.",
      call. = FALSE
    )
  }

  list(
    variance_components = public_variance,
    fa_parameters = merge(parameter_map, group_lookup, by = "PPGroup", all.x = TRUE, sort = FALSE),
    Genetic_covariance_trait_environment = genetic_covariance,
    Genetic_covariance_trait_environment_SE = genetic_covariance_se,
    Genetic_correlation_trait_environment = genetic_correlation,
    Genetic_correlation_trait_environment_SE = genetic_correlation_se,
    genetic_covariance_long = do.call(rbind, genetic_long),
    Residual_covariance_by_environment = residual_covariance,
    Residual_covariance_by_environment_SE = residual_covariance_se,
    Residual_correlation_by_environment = residual_correlation,
    Residual_correlation_by_environment_SE = residual_correlation_se,
    residual_covariance_long = do.call(rbind, residual_long),
    endpoint_diagnostics = endpoint,
    genetic_variance_rhs = genetic_variance_rhs,
    genetic_covariance_rhs = genetic_covariance_rhs,
    residual_variance_rhs = residual_variance_rhs
  )
}

gp_multitrait_asreml_translate_v_rhs <- function(rhs, index_map) {
  out <- as.character(rhs)
  for (index in rev(seq_along(index_map))) {
    out <- gsub(
      paste0("\\bV", index, "\\b"),
      paste0("V", index_map[[index]]),
      out,
      perl = TRUE
    )
  }
  out
}

gp_multitrait_asreml_mtmet_variance_multikernel <- function(model,
                                                             group_lookup,
                                                             fa_rank,
                                                             kernel_specs,
                                                             vpredict_fn = NULL,
                                                             summary_fn = summary,
                                                             tolerance = 1e-8) {
  if (is.null(vpredict_fn)) vpredict_fn <- getFromNamespace("vpredict", "asreml")
  raw_full <- summary_fn(model)$varcomp
  component_names <- rownames(raw_full)
  if (is.null(component_names) || !length(component_names)) {
    stop("ASReml did not return named multi-kernel variance components.", call. = FALSE)
  }
  if (is.null(kernel_specs) || !nrow(kernel_specs)) {
    stop("Multi-kernel MT-MET extraction requires kernel names and inverse-object names.", call. = FALSE)
  }

  all_genetic <- which(grepl("fa\\(PPGroup", component_names) & grepl("vm\\(", component_names))
  residual_rows <- setdiff(seq_len(nrow(raw_full)), all_genetic)
  per_kernel <- stats::setNames(vector("list", nrow(kernel_specs)), kernel_specs$Kernel)
  translated_variance_rhs <- stats::setNames(vector("list", nrow(kernel_specs)), kernel_specs$Kernel)
  translated_covariance_rhs <- stats::setNames(vector("list", nrow(kernel_specs)), kernel_specs$Kernel)

  for (kernel_index in seq_len(nrow(kernel_specs))) {
    inverse_name <- kernel_specs$inverse_name[[kernel_index]]
    genetic_rows <- all_genetic[grepl(inverse_name, component_names[all_genetic], fixed = TRUE)]
    if (!length(genetic_rows)) {
      stop(
        "No grouped MT-MET ASReml components matched kernel `",
        kernel_specs$Kernel[[kernel_index]], "` (", inverse_name, ").",
        call. = FALSE
      )
    }
    keep <- sort(unique(c(genetic_rows, residual_rows)))
    summary_subset <- function(object, ...) list(varcomp = raw_full[keep, , drop = FALSE])
    vpredict_subset <- function(object, xform, ...) {
      rhs <- paste(deparse(xform[[3L]], width.cutoff = 500L), collapse = "")
      rhs <- gp_multitrait_asreml_translate_v_rhs(rhs, keep)
      vpredict_fn(object, stats::as.formula(paste("derived ~", rhs)), ...)
    }
    extracted <- gp_multitrait_asreml_mtmet_variance(
      model = model,
      group_lookup = group_lookup,
      fa_rank = fa_rank,
      vpredict_fn = vpredict_subset,
      summary_fn = summary_subset,
      tolerance = tolerance
    )
    extracted$fa_parameters$Kernel <- kernel_specs$Kernel[[kernel_index]]
    per_kernel[[kernel_index]] <- extracted
    translated_variance_rhs[[kernel_index]] <- vapply(
      extracted$genetic_variance_rhs,
      gp_multitrait_asreml_translate_v_rhs,
      character(1L),
      index_map = keep
    )
    translated_covariance_rhs[[kernel_index]] <- matrix(
      vapply(
        as.vector(extracted$genetic_covariance_rhs),
        gp_multitrait_asreml_translate_v_rhs,
        character(1L),
        index_map = keep
      ),
      nrow = nrow(extracted$genetic_covariance_rhs),
      ncol = ncol(extracted$genetic_covariance_rhs),
      dimnames = dimnames(extracted$genetic_covariance_rhs)
    )
  }

  groups <- as.character(group_lookup$PPGroup)
  genetic_by_kernel <- lapply(per_kernel, `[[`, "Genetic_covariance_trait_environment")
  genetic_covariance <- Reduce(`+`, genetic_by_kernel)
  genetic_covariance_se <- genetic_covariance
  genetic_correlation <- stats::cov2cor(genetic_covariance)
  genetic_correlation_se <- genetic_covariance
  total_variance_rhs <- stats::setNames(vapply(groups, function(group_value) {
    paste(
      vapply(translated_variance_rhs, `[[`, character(1L), group_value),
      collapse = " + "
    )
  }, character(1L)), groups)
  total_covariance_rhs <- matrix(NA_character_, length(groups), length(groups), dimnames = list(groups, groups))
  genetic_long <- list()
  long_index <- 1L
  for (i in seq_along(groups)) {
    for (j in i:length(groups)) {
      group_1 <- groups[[i]]
      group_2 <- groups[[j]]
      rhs <- paste(
        vapply(translated_covariance_rhs, function(x) x[group_1, group_2], character(1L)),
        collapse = " + "
      )
      total_covariance_rhs[group_1, group_2] <- rhs
      total_covariance_rhs[group_2, group_1] <- rhs
      native <- gp_multitrait_asreml_vpredict(model, rhs, vpredict_fn)
      manual <- genetic_covariance[group_1, group_2]
      if (!is.finite(native[["estimate"]]) ||
          abs(native[["estimate"]] - manual) > tolerance * (1 + abs(manual))) {
        stop(
          "Summed multi-kernel FA covariance disagrees with ASReml vpredict for ",
          group_1, " and ", group_2, ".", call. = FALSE
        )
      }
      genetic_covariance[group_1, group_2] <- genetic_covariance[group_2, group_1] <- native[["estimate"]]
      genetic_covariance_se[group_1, group_2] <- genetic_covariance_se[group_2, group_1] <- native[["standard_error"]]
      correlation_native <- if (i == j) {
        c(estimate = 1, standard_error = 0)
      } else {
        gp_multitrait_asreml_vpredict(
          model,
          paste0(
            "(", rhs, ")/sqrt((", total_variance_rhs[[group_1]], ")*(",
            total_variance_rhs[[group_2]], "))"
          ),
          vpredict_fn
        )
      }
      genetic_correlation[group_1, group_2] <- genetic_correlation[group_2, group_1] <- correlation_native[["estimate"]]
      genetic_correlation_se[group_1, group_2] <- genetic_correlation_se[group_2, group_1] <- correlation_native[["standard_error"]]
      row_1 <- group_lookup[match(group_1, group_lookup$PPGroup), , drop = FALSE]
      row_2 <- group_lookup[match(group_2, group_lookup$PPGroup), , drop = FALSE]
      genetic_long[[long_index]] <- data.frame(
        Trait_1 = row_1$Trait, Env_1 = row_1$Env,
        Trait_2 = row_2$Trait, Env_2 = row_2$Env,
        Covariance = native[["estimate"]],
        Standard_error = native[["standard_error"]],
        Correlation = correlation_native[["estimate"]],
        Correlation_standard_error = correlation_native[["standard_error"]],
        stringsAsFactors = FALSE
      )
      long_index <- long_index + 1L
    }
  }
  if (min(eigen(genetic_covariance, symmetric = TRUE, only.values = TRUE)$values) < -tolerance) {
    stop("The summed multi-kernel ASReml genetic covariance is not positive semidefinite.", call. = FALSE)
  }

  first <- per_kernel[[1L]]
  residual_rhs <- vapply(
    first$residual_variance_rhs,
    gp_multitrait_asreml_translate_v_rhs,
    character(1L),
    index_map = sort(unique(c(
      all_genetic[grepl(kernel_specs$inverse_name[[1L]], component_names[all_genetic], fixed = TRUE)],
      residual_rows
    )))
  )
  endpoint <- gp_multitrait_asreml_endpoint_diagnostics(model, summary_fn = summary_fn)
  public_rows <- vector("list", nrow(group_lookup))
  kernel_variance_rows <- vector("list", nrow(group_lookup) * nrow(kernel_specs))
  kernel_row_index <- 1L
  for (i in seq_len(nrow(group_lookup))) {
    group_value <- as.character(group_lookup$PPGroup[[i]])
    vg_rhs <- total_variance_rhs[[group_value]]
    ve_rhs <- residual_rhs[[group_value]]
    vg <- gp_multitrait_asreml_vpredict(model, vg_rhs, vpredict_fn)
    ve <- gp_multitrait_asreml_vpredict(model, ve_rhs, vpredict_fn)
    h2 <- gp_multitrait_asreml_vpredict(
      model,
      paste0("(", vg_rhs, ")/((", vg_rhs, ")+", ve_rhs, ")"),
      vpredict_fn
    )
    public_rows[[i]] <- data.frame(
      Trait = group_lookup$Trait[[i]],
      Env = group_lookup$Env[[i]],
      Component = c("genetic_variance", "residual_variance", "heritability"),
      Components = c(vg[["estimate"]], ve[["estimate"]], h2[["estimate"]]),
      Standard_error = c(vg[["standard_error"]], ve[["standard_error"]], h2[["standard_error"]]),
      Estimation_method = "ASReml_REML_vpredict",
      Variance_estimand = c(
        "sum_of_independent_kernel_genetic_variances_for_trait_environment_group",
        "observation_level_residual_variance_within_environment",
        "observation_level_narrow_sense_h2=sum(Vg_kernel)/(sum(Vg_kernel)+Ve)"
      ),
      Fit_converged = endpoint$converged,
      Endpoint_max_percent_change = endpoint$maximum_percent_change,
      stringsAsFactors = FALSE
    )
    for (kernel_index in seq_len(nrow(kernel_specs))) {
      rhs <- translated_variance_rhs[[kernel_index]][[group_value]]
      value <- gp_multitrait_asreml_vpredict(model, rhs, vpredict_fn)
      kernel_variance_rows[[kernel_row_index]] <- data.frame(
        Kernel = kernel_specs$Kernel[[kernel_index]],
        Trait = group_lookup$Trait[[i]],
        Env = group_lookup$Env[[i]],
        Component = "kernel_genetic_variance",
        Components = value[["estimate"]],
        Standard_error = value[["standard_error"]],
        Estimation_method = "ASReml_REML_vpredict",
        stringsAsFactors = FALSE
      )
      kernel_row_index <- kernel_row_index + 1L
    }
  }
  public_variance <- do.call(rbind, public_rows)
  if (any(!is.finite(public_variance$Components)) ||
      any(!is.finite(public_variance$Standard_error))) {
    stop(
      "ASReml multi-kernel vpredict returned a non-finite total variance, heritability, or standard error; variance output was not promoted.",
      call. = FALSE
    )
  }

  kernel_long <- do.call(rbind, lapply(seq_along(per_kernel), function(i) {
    out <- per_kernel[[i]]$genetic_covariance_long
    out$Kernel <- kernel_specs$Kernel[[i]]
    out[, c("Kernel", setdiff(names(out), "Kernel")), drop = FALSE]
  }))
  fa_parameters <- do.call(rbind, lapply(per_kernel, `[[`, "fa_parameters"))
  rownames(fa_parameters) <- NULL
  list(
    variance_components = public_variance,
    kernel_variance_components = do.call(rbind, kernel_variance_rows),
    fa_parameters = fa_parameters,
    Genetic_covariance_trait_environment = genetic_covariance,
    Genetic_covariance_trait_environment_SE = genetic_covariance_se,
    Genetic_correlation_trait_environment = genetic_correlation,
    Genetic_correlation_trait_environment_SE = genetic_correlation_se,
    Genetic_covariance_by_kernel = genetic_by_kernel,
    genetic_covariance_long = do.call(rbind, genetic_long),
    genetic_covariance_by_kernel_long = kernel_long,
    Residual_covariance_by_environment = first$Residual_covariance_by_environment,
    Residual_covariance_by_environment_SE = first$Residual_covariance_by_environment_SE,
    Residual_correlation_by_environment = first$Residual_correlation_by_environment,
    Residual_correlation_by_environment_SE = first$Residual_correlation_by_environment_SE,
    residual_covariance_long = first$residual_covariance_long,
    endpoint_diagnostics = endpoint,
    genetic_variance_rhs = total_variance_rhs,
    genetic_covariance_rhs = total_covariance_rhs,
    residual_variance_rhs = residual_rhs
  )
}

gp_multitrait_asreml_mtmet_observed_grid <- function(long_bundle, gen_name) {
  long <- long_bundle$data
  key <- paste(
    as.character(long[[gen_name]]),
    as.character(long$PPGroup),
    sep = "\r"
  )
  rows <- split(seq_len(nrow(long)), key, drop = TRUE)
  out <- lapply(rows, function(index) {
    values <- suppressWarnings(as.numeric(long$Observed_value[index]))
    values <- values[is.finite(values)]
    data.frame(
      id_value = as.character(long[[gen_name]][index[1L]]),
      PPGroup = as.character(long$PPGroup[index[1L]]),
      Observed_value = if (length(values)) mean(values) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  names(out)[names(out) == "id_value"] <- gen_name
  rownames(out) <- NULL
  out
}

gp_multitrait_asreml_mtmet_prediction <- function(pred_raw,
                                                   long_bundle,
                                                   gen_name) {
  pred_raw <- as.data.frame(pred_raw, stringsAsFactors = FALSE, check.names = FALSE)
  gid_col <- gp_contract_first_name(names(pred_raw), c(gen_name, "GID", "gen"))
  group_col <- gp_contract_first_name(names(pred_raw), c("PPGroup", "group"))
  value_col <- gp_contract_first_name(names(pred_raw), c("predicted.value", "Predicted_value"))
  se_col <- gp_contract_first_name(names(pred_raw), c("std.error", "Standard_error"))
  if (anyNA(c(gid_col, group_col, value_col))) {
    stop(
      "Could not identify genotype, PPGroup, and predicted.value columns in ",
      "the grouped MT-MET ASReml prediction table.",
      call. = FALSE
    )
  }
  pred <- data.frame(
    id_value = as.character(pred_raw[[gid_col]]),
    PPGroup = as.character(pred_raw[[group_col]]),
    Predicted_value = as.numeric(pred_raw[[value_col]]),
    Standard_error = if (!is.na(se_col)) as.numeric(pred_raw[[se_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  names(pred)[1L] <- gen_name
  if (anyDuplicated(pred[c(gen_name, "PPGroup")])) {
    stop("ASReml returned duplicate genotype-by-group predictions.", call. = FALSE)
  }
  pred <- merge(
    pred,
    long_bundle$group_lookup,
    by = "PPGroup",
    all.x = TRUE,
    sort = FALSE
  )
  observed <- gp_multitrait_asreml_mtmet_observed_grid(long_bundle, gen_name)
  pred <- merge(
    pred,
    observed,
    by = c(gen_name, "PPGroup"),
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(pred$Trait) || anyNA(pred$Env)) {
    stop("An ASReml prediction group was absent from the trait-environment lookup.", call. = FALSE)
  }
  pred$Prediction_error_variance <- pred$Standard_error^2
  pred$Train_Test_Label <- ifelse(is.na(pred$Observed_value), "Test", "Train")
  pred$PPGroup <- factor(pred$PPGroup, levels = long_bundle$group_lookup$PPGroup)
  pred[[gen_name]] <- factor(pred[[gen_name]], levels = long_bundle$genotype_levels)
  pred <- pred[order(pred$PPGroup, pred[[gen_name]]), , drop = FALSE]
  pred$PPGroup <- as.character(pred$PPGroup)
  pred[[gen_name]] <- as.character(pred[[gen_name]])
  pred <- pred[, c(
    gen_name, "Trait", "Env", "Observed_value", "Predicted_value",
    "Standard_error", "Prediction_error_variance", "Train_Test_Label", "PPGroup"
  ), drop = FALSE]
  rownames(pred) <- NULL
  pred
}

gp_multitrait_asreml_mtmet_predict_or_chunk <- function(mod, predict_fn,
                                                        long_bundle, gen_name,
                                                        pworkspace = NULL) {
  classify <- paste(gen_name, "PPGroup", sep = ":")
  full <- tryCatch(
    gp_multitrait_asreml_predict_pvals(
      mod = mod, predict_fn = predict_fn, classify = classify,
      pworkspace = pworkspace
    ),
    error = identity
  )
  if (!inherits(full, "error")) {
    return(list(pvals = full, mode = "predict"))
  }
  if (!gp_multitrait_asreml_workspace_error(conditionMessage(full))) {
    stop("Grouped MT-MET ASReml prediction failed: ", conditionMessage(full), call. = FALSE)
  }
  asreml_validate_coefficient_fallback_model(
    mod, context = "Grouped MT-MET ASReml prediction fallback"
  )
  groups <- as.character(long_bundle$group_lookup$PPGroup)
  chunks <- lapply(groups, function(group) {
    tryCatch(
      gp_multitrait_asreml_predict_pvals(
        mod = mod, predict_fn = predict_fn, classify = classify,
        pworkspace = pworkspace, levels = list(PPGroup = group)
      ),
      error = function(e) {
        stop("Grouped MT-MET ASReml prediction failed for PPGroup `", group,
             "`: ", conditionMessage(e), call. = FALSE)
      }
    )
  })
  combined <- do.call(rbind, chunks)
  gid_col <- gp_contract_first_name(names(combined), c(gen_name, "GID", "gen"))
  group_col <- gp_contract_first_name(names(combined), c("PPGroup", "group"))
  if (anyNA(c(gid_col, group_col))) {
    stop("Grouped MT-MET ASReml chunk predictions lacked genotype or PPGroup labels.",
         call. = FALSE)
  }
  actual <- paste(as.character(combined[[gid_col]]),
                  as.character(combined[[group_col]]), sep = "\r")
  expected <- as.vector(outer(
    as.character(long_bundle$genotype_levels), groups,
    paste, sep = "\r"
  ))
  if (anyDuplicated(actual) || length(actual) != length(expected) ||
      !setequal(actual, expected)) {
    stop("Grouped MT-MET ASReml chunk predictions did not cover each genotype-by-group cell exactly once.",
         call. = FALSE)
  }
  rownames(combined) <- NULL
  list(pvals = combined, mode = "predict_group_chunks")
}

gp_multitrait_asreml_mtmet_gaussian_model <- function(pheno_object,
                                                       response,
                                                       gen_name,
                                                       heter_groups,
                                                       gmatrix = NULL,
                                                       kernel_list = NULL,
                                                       heter_resid = TRUE,
                                                       var_cov_str = NULL,
                                                       inverse = TRUE,
                                                       epsilon = 1e-6,
                                                       engine = "asreml",
                                                       workspace = 1e08,
                                                       pworkspace = 1e06,
                                                       maxit = 50L,
                                                       ai_sing = TRUE,
                                                       max_stability_updates = 4L,
                                                       extract_variance = TRUE,
                                                       require_variance_stability = extract_variance) {
  if (!engine %in% rownames(utils::installed.packages())) {
    stop("You need to install asreml-R to use grouped MT-MET GBLUP.", call. = FALSE)
  }
  if (!isTRUE(heter_resid)) {
    stop(
      "Grouped MT-MET ASReml currently requires heter_resid = TRUE so each ",
      "environment receives its own cross-trait residual covariance.",
      call. = FALSE
    )
  }
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE, check.names = FALSE)
  ids <- unique(as.character(ph[[gen_name]]))
  kernel_bundle <- gp_multitrait_asreml_kernel_bank(
    gmatrix = gmatrix,
    kernel_list = kernel_list,
    phenotype_ids = ids,
    context = "Grouped MT-MET ASReml"
  )
  keep_ids <- kernel_bundle$keep_ids
  ph <- ph[as.character(ph[[gen_name]]) %in% keep_ids, , drop = FALSE]
  long_bundle <- gp_multitrait_asreml_mtmet_long_data(
    pheno_object = ph,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  structure_spec <- gp_multitrait_asreml_mtmet_structure(
    var_cov_str = var_cov_str,
    n_groups = nrow(long_bundle$group_lookup)
  )
  long_bundle$genotype_levels <- keep_ids
  long_bundle$data[[gen_name]] <- factor(
    as.character(long_bundle$data[[gen_name]]),
    levels = keep_ids
  )

  kernel_specs <- gp_multitrait_asreml_assign_kernel_inverses(
    kernels = kernel_bundle$kernels,
    prefix = "pp_mtmet_inv",
    inverse = inverse,
    epsilon = epsilon
  )
  on.exit({
    for (inverse_name in kernel_specs$inverse_name) remove_from_global(inverse_name)
  }, add = TRUE)

  asreml_fn <- getFromNamespace("asreml", engine)
  predict_fn <- getOption("PredictProR.asreml_predict_impl", getFromNamespace("predict.asreml", engine))
  update_fn <- getFromNamespace("update.asreml", engine)
  asreml_options_fn <- getFromNamespace("asreml.options", engine)
  if (isTRUE(ai_sing)) {
    old_options <- tryCatch(asreml_options_fn(), error = function(e) NULL)
    old_ai_sing <- if (is.list(old_options) && "ai.sing" %in% names(old_options)) old_options[["ai.sing"]] else NULL
    try(asreml_options_fn(ai.sing = TRUE), silent = TRUE)
    on.exit({
      if (!is.null(old_ai_sing)) try(asreml_options_fn(ai.sing = old_ai_sing), silent = TRUE)
    }, add = TRUE)
  }

  long_data <- long_bundle$data
  fit_args <- list(
    fixed = stats::as.formula("Observed_value ~ PPGroup"),
    random = gp_multitrait_asreml_random_formula(
      index_factor = "PPGroup",
      covariance = paste0("fa(PPGroup,", structure_spec$rank, ")"),
      gen_name = gen_name,
      kernel_specs = kernel_specs
    ),
    residual = stats::as.formula("~ dsum(~ PPUnit:corgh(PPTrait) | PPEnv)"),
    data = long_data,
    na.action = list(x = "include", y = "include"),
    workspace = workspace,
    maxit = as.integer(maxit),
    ai.sing = isTRUE(ai_sing)
  )
  mod <- tryCatch(
    do.call(asreml_fn, fit_args),
    error = function(e) {
      stop("Grouped MT-MET ASReml GBLUP fit failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  stabilized <- gp_multitrait_asreml_fit_until_stable(
    model = mod,
    update_fn = update_fn,
    max_updates = max_stability_updates,
    max_percent_change = if (isTRUE(require_variance_stability)) 1 else Inf
  )
  mod <- stabilized$model
  if (!isTRUE(tail(stabilized$trace$converged, 1L))) {
    stop(
      "Grouped MT-MET ASReml GBLUP did not converge; predictions and variance ",
      "output were not returned.",
      call. = FALSE
    )
  }
  if (isTRUE(require_variance_stability) && !isTRUE(stabilized$stable)) {
    final <- tail(stabilized$trace, 1L)
    stop(
      "Grouped MT-MET ASReml GBLUP did not reach a variance-stable endpoint ",
      "after ", max_stability_updates, " continuation fits (converged = ",
      final$converged, ", maximum component change = ",
      signif(final$maximum_percent_change, 6), "%). Variance output was not promoted.",
      call. = FALSE
    )
  }

  prediction_result <- gp_multitrait_asreml_mtmet_predict_or_chunk(
    mod = mod, predict_fn = predict_fn, long_bundle = long_bundle,
    gen_name = gen_name, pworkspace = pworkspace
  )
  pred <- gp_multitrait_asreml_mtmet_prediction(
    pred_raw = prediction_result$pvals,
    long_bundle = long_bundle,
    gen_name = gen_name
  )
  variance <- if (isTRUE(extract_variance)) {
    if (nrow(kernel_specs) > 1L) {
      gp_multitrait_asreml_mtmet_variance_multikernel(
        model = mod,
        group_lookup = long_bundle$group_lookup,
        fa_rank = structure_spec$rank,
        kernel_specs = kernel_specs
      )
    } else {
      gp_multitrait_asreml_mtmet_variance(
        model = mod,
        group_lookup = long_bundle$group_lookup,
        fa_rank = structure_spec$rank
      )
    }
  } else {
    list(
      variance_components = data.frame(
        Trait = character(), Env = character(), Component = character(),
        Components = numeric(), Standard_error = numeric(),
        Estimation_method = character(), stringsAsFactors = FALSE
      )
    )
  }
  total_pred <- gp_multitrait_across_environment_prediction(
    pred,
    gen_name = gen_name,
    heter_groups = "Env"
  )

  c(list(
    Asreml_model = mod,
    Predicted_value = pred,
    predicted_values = pred,
    Total_Predicted_value = total_pred,
    total_predicted_values = total_pred,
    Variance_components = variance$variance_components,
    variance_components = variance$variance_components,
    group_lookup = long_bundle$group_lookup,
    asreml_mtmet_structure = structure_spec$token,
    asreml_prediction_mode = prediction_result$mode,
    asreml_convergence_trace = stabilized$trace,
    kernel_configuration = gp_multitrait_kernel_configuration(kernel_bundle$kernels),
    model_parameters = rbind(data.frame(
      stat = c("mode", "model", "genetic_structure", "residual_structure",
               "variance_endpoint", "asreml_prediction_mode"),
      summary = c(
        "multi_trait_multi_environment_asreml",
        "GBLUP",
        paste0(
          nrow(kernel_specs), " independent fa(PPGroup,", structure_spec$rank,
          "):vm(", gen_name, ",Kinv) terms; covariance summed"
        ),
        "dsum(PPUnit:corgh(PPTrait)|PPEnv)",
        if (isTRUE(extract_variance) && isTRUE(require_variance_stability)) {
          "converged_and_maximum_component_change_le_1_percent"
        } else {
          "converged_prediction_fit_variance_not_extracted"
        },
        prediction_result$mode
      ),
      stringsAsFactors = FALSE
    ), gp_multi_kernel_parameter_rows(
      kernel_names = kernel_specs$Kernel,
      strategy = "independent_FA_covariance_per_kernel_summed"
    )),
    prediction_uncertainty_source = "ASReml prediction standard errors and PEV",
    model_notes = data.frame(
      note = c(
        "PPGroup is an internal collision-safe key mapped explicitly to Trait and Env; component extraction never splits labels on underscores.",
        "Genetic variances and covariances are derived from named FA specific variances and loadings and verified against ASReml vpredict.",
        "With multiple kernels, each kernel has an independent FA genetic term; the reported genetic covariance and heritability use the vpredict-verified sum across kernels.",
        "Residual section scale rows fixed at one are ASReml sigma-parameterization bookkeeping and are excluded.",
        "Reported heritability is observation-level Vg/(Vg+Ve), not entry-mean heritability or prediction reliability."
      ),
      stringsAsFactors = FALSE
    )
  ), variance[setdiff(names(variance), "variance_components")])
}

gp_multitrait_asreml_mtmet_metric_table <- function(pred_long,
                                                     eval_metrics,
                                                     model_type = "GBLUP",
                                                     rep = 1L) {
  rows <- list()
  index <- 1L
  groups <- unique(pred_long[c("Trait", "Env")])
  for (group_index in seq_len(nrow(groups))) {
    trait <- as.character(groups$Trait[group_index])
    environment <- as.character(groups$Env[group_index])
    dat <- pred_long[
      pred_long$Trait == trait & pred_long$Env == environment &
        is.finite(pred_long$Observed_value) & is.finite(pred_long$Predicted_value),
      ,
      drop = FALSE
    ]
    for (fold in sort(unique(dat$fold))) {
      fold_data <- dat[dat$fold == fold, , drop = FALSE]
      for (metric in eval_metrics) {
        value <- tryCatch(
          evaluation_metrics(
            y_observed = fold_data$Observed_value,
            y_predicted = fold_data$Predicted_value,
            eval_metrics = metric,
            response_family = "gaussian"
          ),
          error = function(e) NA_real_
        )
        rows[[index]] <- data.frame(
          trait = trait,
          Env = environment,
          model = model_type,
          rep = rep,
          fold = fold,
          metric = metric,
          value = as.numeric(value),
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      trait = character(), Env = character(), model = character(),
      rep = integer(), fold = integer(), metric = character(), value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

gp_multitrait_asreml_mtmet_gaussian_cv <- function(pheno_object,
                                                    response,
                                                    gen_name,
                                                    heter_groups,
                                                    gmatrix = NULL,
                                                    kernel_list = NULL,
                                                    cross_validation_meth = "CV1",
                                                    nfolds = 5L,
                                                    replication = 1L,
                                                    eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                                    heter_resid = TRUE,
                                                    var_cov_str = NULL,
                                                    inverse = TRUE,
                                                    epsilon = 1e-6,
                                                    engine = "asreml",
                                                    workspace = 1e08,
                                                    pworkspace = 1e06,
                                                    maxit = 50L,
                                                    ai_sing = TRUE,
                                                    random_seed = 123L) {
  cv_token <- normalize_cv_token(cross_validation_meth)
  if (!cv_token %in% c("cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2")) {
    stop(
      "Grouped MT-MET ASReml CV supports CV0, CV1, CV2, Repeated_CV0, Repeated_CV1, and Repeated_CV2 only.",
      call. = FALSE
    )
  }
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE, check.names = FALSE)
  kernel_bundle <- gp_multitrait_asreml_kernel_bank(
    gmatrix = gmatrix,
    kernel_list = kernel_list,
    phenotype_ids = unique(as.character(ph[[gen_name]])),
    context = "Grouped MT-MET ASReml CV"
  )
  ids <- kernel_bundle$keep_ids
  if (length(ids) < 2L) {
    stop("Grouped MT-MET ASReml CV requires at least two genotypes shared with gmatrix.", call. = FALSE)
  }
  ph <- ph[as.character(ph[[gen_name]]) %in% ids, , drop = FALSE]
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
  original_bundle <- gp_multitrait_asreml_mtmet_long_data(
    ph,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  original_observed <- gp_multitrait_asreml_mtmet_observed_grid(original_bundle, gen_name)
  original_observed <- merge(
    original_observed,
    original_bundle$group_lookup,
    by = "PPGroup",
    all.x = TRUE,
    sort = FALSE
  )

  cv_results <- vector("list", length(plan$folds))
  for (rep_index in seq_along(plan$folds)) {
    prediction_blocks <- list()
    block_index <- 1L
    for (fold_index in seq_along(plan$folds[[rep_index]])) {
      test_idx <- sort(unique(as.integer(plan$folds[[rep_index]][[fold_index]])))
      ph_masked <- ph
      ph_masked[test_idx, response] <- NA_real_
      fold_fit <- gp_multitrait_asreml_mtmet_gaussian_model(
        pheno_object = ph_masked,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups,
        gmatrix = kernel_bundle$kernels[[1L]],
        kernel_list = if (length(kernel_bundle$kernels) > 1L) kernel_bundle$kernels[-1L] else NULL,
        heter_resid = heter_resid,
        var_cov_str = var_cov_str,
        inverse = inverse,
        epsilon = epsilon,
        engine = engine,
        workspace = workspace,
        pworkspace = pworkspace,
        maxit = maxit,
        ai_sing = ai_sing,
        extract_variance = FALSE,
        require_variance_stability = FALSE
      )
      pred <- fold_fit$predicted_values
      test_cell_key <- paste(
        as.character(ph[[gen_name]][test_idx]),
        as.character(ph[[heter_groups]][test_idx]),
        sep = "\r"
      )
      prediction_cell_key <- paste(
        as.character(pred[[gen_name]]),
        as.character(pred$Env),
        sep = "\r"
      )
      pred <- pred[prediction_cell_key %in% test_cell_key, , drop = FALSE]
      observed_key <- paste(
        as.character(original_observed[[gen_name]]),
        as.character(original_observed$PPGroup),
        sep = "\r"
      )
      prediction_key <- paste(as.character(pred[[gen_name]]), as.character(pred$PPGroup), sep = "\r")
      pred$Observed_value <- original_observed$Observed_value[match(prediction_key, observed_key)]
      expected <- is.finite(pred$Observed_value)
      if (any(expected & !is.finite(pred$Predicted_value))) {
        stop(
          "Grouped MT-MET ASReml CV produced missing predictions for observed ",
          "held-out cells in replication ", rep_index, ", fold ", fold_index, ".",
          call. = FALSE
        )
      }
      pred$model <- "GBLUP"
      pred$fold <- fold_index
      pred$rep <- rep_index
      pred$cv_role <- "test"
      pred$cv_scenario <- toupper(sub("^repeated_", "", cv_token))
      pred$Train_Test_Label <- "Test"
      prediction_blocks[[block_index]] <- pred
      block_index <- block_index + 1L
    }
    pred_long <- do.call(rbind, prediction_blocks)
    rownames(pred_long) <- NULL
    metrics <- gp_multitrait_asreml_mtmet_metric_table(
      pred_long,
      eval_metrics = eval_metrics,
      model_type = "GBLUP",
      rep = rep_index
    )
    cv_results[[rep_index]] <- list(
      trait = paste(response, collapse = ","),
      model = "GBLUP",
      rep = rep_index,
      eval_metrics_reps = metrics,
      ypred_cv_Reps_all = pred_long,
      yprob_cv_Reps_all = NULL
    )
  }
  processed <- gp_multitrait_ml_cv_process(cv_results)
  processed$multitrait_oof_predictions <- do.call(
    rbind,
    lapply(cv_results, `[[`, "ypred_cv_Reps_all")
  )
  processed$fold_assignments <- plan$manifest
  processed$kernel_configuration <-
    gp_multitrait_kernel_configuration(kernel_bundle$kernels)
  processed$model_parameters <- gp_multi_kernel_parameter_rows(
    kernel_names = names(kernel_bundle$kernels),
    strategy = "independent_FA_covariance_per_kernel_summed"
  )
  processed$cv_scheme <- list(
    unit = if (cv_token %in% c("cv0", "repeated_cv0")) {
      "environment"
    } else if (cv_token %in% c("cv1", "repeated_cv1")) {
      "genotype"
    } else {
      "genotype_environment_cell"
    },
    all_traits_masked_together = TRUE,
    method = cv_token,
    environment_column = heter_groups,
    requested_nfolds = as.integer(nfolds),
    effective_nfolds = plan$effective_nfolds,
    replication = as.integer(replication),
    kernel_count = length(kernel_bundle$kernels),
    kernel_names = paste(names(kernel_bundle$kernels), collapse = ";"),
    variance_source = "model_refit_per_fold; biological variance components come from the direct/final fit"
  )
  processed$mtmet_metric_summary_by_environment <- if (length(cv_results)) {
    metrics_all <- do.call(rbind, lapply(cv_results, `[[`, "eval_metrics_reps"))
    stats::aggregate(
      value ~ trait + Env + model + metric,
      data = metrics_all,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
  } else {
    data.frame()
  }
  list(cv_results = cv_results, cv_results_processed = processed)
}
