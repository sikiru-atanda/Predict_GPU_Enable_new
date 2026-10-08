pp_escape_regex <- function(x) {
  gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x)
}

pp_asreml_coef_col <- function(x, candidates, what) {
  hit <- candidates[candidates %in% colnames(x)][1]
  if (is.na(hit) || is.null(hit)) {
    stop("Unable to locate ", what, " column in ASReml coefficient table.")
  }
  hit
}

pp_asreml_aggregate_effects <- function(df, group_cols, effect_col, se_col = NULL) {
  if (!nrow(df)) {
    return(df)
  }

  if (!is.null(se_col) && se_col %in% names(df)) {
    stats::aggregate(
      df[, c(effect_col, se_col), drop = FALSE],
      df[, group_cols, drop = FALSE],
      function(x) sum(x, na.rm = TRUE)
    )
  } else {
    out <- stats::aggregate(
      df[, effect_col, drop = FALSE],
      df[, group_cols, drop = FALSE],
      function(x) sum(x, na.rm = TRUE)
    )
    out[["SE"]] <- NA_real_
    out
  }
}

extract_asreml_prediction <- function(model,
                                      genotype_term = "GID",
                                      environment_levels = NULL,
                                      heter_grp = "Env",
                                      kernel_ids = NULL) {
  asreml_validate_coefficient_fallback_model(
    model = model,
    context = "ASReml single-trait coefficient fallback"
  )
  random_coefs <- summary(model, coef = TRUE)$coef.random
  fixed_coefs <- coef(model)$fixed

  if (is.null(random_coefs) || !nrow(random_coefs)) {
    stop("No ASReml random coefficients are available for fallback prediction extraction.")
  }

  coef_rows <- rownames(random_coefs)
  solution_col <- pp_asreml_coef_col(random_coefs, c("solution", "effect"), "solution")
  se_col <- if (any(c("std.error", "Standard_error") %in% colnames(random_coefs))) {
    pp_asreml_coef_col(random_coefs, c("std.error", "Standard_error"), "standard error")
  } else {
    NULL
  }

  fixed_effect_col <- if (!is.null(fixed_coefs) && nrow(fixed_coefs)) {
    pp_asreml_coef_col(fixed_coefs, c("effect", "solution"), "fixed-effect solution")
  } else {
    NULL
  }

  intercept <- if (!is.null(fixed_coefs) &&
                   "(Intercept)" %in% rownames(fixed_coefs) &&
                   !is.null(fixed_effect_col)) {
    fixed_coefs["(Intercept)", fixed_effect_col]
  } else {
    0
  }

  geno_rx <- paste0("vm\\(", pp_escape_regex(genotype_term), "\\s*,")
  vm_rows <- coef_rows[grepl(geno_rx, coef_rows)]

  if (!is.null(kernel_ids) && length(kernel_ids)) {
    kernel_rx <- paste(pp_escape_regex(as.character(kernel_ids)), collapse = "|")
    vm_rows <- vm_rows[grepl(
      paste0("vm\\(", pp_escape_regex(genotype_term), "\\s*,\\s*(", kernel_rx, ")\\)"),
      vm_rows
    )]
  }

  if (!length(vm_rows)) {
    stop("No ASReml random coefficient rows matched the genotype term.")
  }

  environment_levels <- unique(trimws(as.character(environment_levels)))
  multi_env <- !is.null(environment_levels) &&
    length(environment_levels) > 1 &&
    !is.null(heter_grp) &&
    length(heter_grp) == 1L &&
    nzchar(heter_grp)

  main_rows <- vm_rows[
    grepl(paste0("^", geno_rx), vm_rows) &
      !grepl(":fa\\(", vm_rows)
  ]

  if (length(main_rows)) {
    g_df <- data.frame(
      Term = main_rows,
      Genotype = sub("^.*\\)_", "", main_rows),
      G = random_coefs[main_rows, solution_col],
      stringsAsFactors = FALSE
    )
    if (!is.null(se_col)) {
      g_df$SE <- random_coefs[main_rows, se_col]
      g_df <- pp_asreml_aggregate_effects(g_df, "Genotype", "G", "SE")
    } else {
      g_df <- stats::aggregate(G ~ Genotype, g_df, sum)
      g_df$SE <- NA_real_
    }
  } else {
    g_df <- data.frame(Genotype = character(), G = numeric(), SE = numeric(), stringsAsFactors = FALSE)
  }

  if (!multi_env) {
    if (!nrow(g_df)) {
      stop("No genotype main-effect rows were available for single-environment fallback prediction extraction.")
    }
    g_df$PEV <- ifelse(is.na(g_df$SE), NA_real_, g_df$SE^2)
    g_df$Predicted <- intercept + g_df$G
    return(g_df[, c("Genotype", "G", "SE", "PEV", "Predicted")])
  }

  heter_rx <- pp_escape_regex(heter_grp)
  env_vm_rows <- vm_rows[grepl(paste0("^", heter_rx, "_[^:]+:vm\\("), vm_rows)]
  fa_rows <- vm_rows[grepl(":fa\\(", vm_rows)]
  # Reduced-rank terms are fitted as `fa(Env, k):vm(GID, K)` / `rr(Env, k):vm(...)`.
  # ASReml names the per-environment genetic effects
  # `fa(Env, k)_<env>:vm(GID, K)_<gid>` (loadings x scores + specific effect);
  # the `_Comp<j>` rows are the latent scores and are not added again.
  rr_rx <- paste0("^(fa|rr)\\(\\s*", heter_rx, "\\s*,\\s*[0-9]+\\s*\\)_([^:]+):vm\\(")
  rr_rows <- vm_rows[grepl(rr_rx, vm_rows)]
  rr_env <- sub(":.*$", "", sub(rr_rx, "\\2:", rr_rows))
  rr_rows <- rr_rows[!grepl("^Comp[0-9]+$", rr_env)]
  rr_env <- rr_env[!grepl("^Comp[0-9]+$", rr_env)]

  gxe_parts <- list()

  if (length(rr_rows)) {
    rr_df <- data.frame(
      Term = rr_rows,
      Genotype = sub("^.*\\)_", "", rr_rows),
      Environment = rr_env,
      GxE_Effect = random_coefs[rr_rows, solution_col],
      stringsAsFactors = FALSE
    )
    if (!is.null(se_col)) {
      rr_df$SE <- random_coefs[rr_rows, se_col]
    }
    gxe_parts[[length(gxe_parts) + 1L]] <- rr_df
  }

  if (length(env_vm_rows)) {
    env_vm_df <- data.frame(
      Term = env_vm_rows,
      Genotype = sub("^.*\\)_", "", env_vm_rows),
      Environment = sub(paste0("^", heter_rx, "_"), "", sub(":.*$", "", env_vm_rows)),
      GxE_Effect = random_coefs[env_vm_rows, solution_col],
      stringsAsFactors = FALSE
    )
    if (!is.null(se_col)) {
      env_vm_df$SE <- random_coefs[env_vm_rows, se_col]
    }
    gxe_parts[[length(gxe_parts) + 1L]] <- env_vm_df
  }

  if (length(fa_rows)) {
    fa_df <- data.frame(
      Term = fa_rows,
      Genotype = sub(":fa\\(.*$", "", sub("^.*\\)_", "", fa_rows)),
      Environment = sub("^.*:fa\\([^)]*\\)_", "", fa_rows),
      GxE_Effect = random_coefs[fa_rows, solution_col],
      stringsAsFactors = FALSE
    )
    if (!is.null(se_col)) {
      fa_df$SE <- random_coefs[fa_rows, se_col]
    }
    gxe_parts[[length(gxe_parts) + 1L]] <- fa_df
  }

  if (!length(gxe_parts)) {
    stop("No multi-environment ASReml random coefficient rows were matched for fallback prediction extraction.")
  }

  gxe_df <- do.call(rbind, gxe_parts)
  gxe_df$Environment <- trimws(as.character(gxe_df$Environment))
  gxe_df <- gxe_df[gxe_df$Environment %in% environment_levels, , drop = FALSE]

  if (!nrow(gxe_df)) {
    stop("No valid multi-environment ASReml rows remained after environment filtering.")
  }

  if (!is.null(se_col) && "SE" %in% names(gxe_df)) {
    gxe_df <- pp_asreml_aggregate_effects(gxe_df, c("Genotype", "Environment"), "GxE_Effect", "SE")
  } else {
    gxe_df <- stats::aggregate(GxE_Effect ~ Genotype + Environment, gxe_df, sum)
    gxe_df$SE <- NA_real_
  }

  if (!nrow(g_df)) {
    g_df <- data.frame(Genotype = unique(gxe_df$Genotype), G = 0, SE = NA_real_, stringsAsFactors = FALSE)
  }

  gxe_df <- merge(gxe_df, g_df[, c("Genotype", "G"), drop = FALSE], by = "Genotype", all.x = TRUE)

  env_effects <- stats::setNames(rep(intercept, length(environment_levels)), environment_levels)

  if (!is.null(fixed_coefs) && nrow(fixed_coefs) && !is.null(fixed_effect_col)) {
    env_fixed_rows <- grep(paste0("^", heter_rx, "_"), rownames(fixed_coefs), value = TRUE)
    if (length(env_fixed_rows)) {
      fixed_env_map <- stats::setNames(
        fixed_coefs[env_fixed_rows, fixed_effect_col],
        sub(paste0("^", heter_rx, "_"), "", env_fixed_rows)
      )
      missing_envs <- setdiff(environment_levels, names(fixed_env_map))
      if (length(missing_envs) == 1L) {
        env_effects[missing_envs] <- intercept
      }
      env_effects[names(fixed_env_map)] <- intercept + fixed_env_map
    }
  } else {
    env_random_rows <- coef_rows[
      grepl(paste0("^", heter_rx, "_"), coef_rows) &
        !grepl(":vm\\(", coef_rows)
    ]
    if (length(env_random_rows)) {
      env_random_map <- stats::setNames(
        random_coefs[env_random_rows, solution_col],
        sub(paste0("^", heter_rx, "_"), "", env_random_rows)
      )
      env_effects[names(env_random_map)] <- intercept + env_random_map
    }
  }

  gxe_df$Env_Effect <- env_effects[gxe_df$Environment]
  gxe_df$PEV <- ifelse(is.na(gxe_df$SE), NA_real_, gxe_df$SE^2)
  gxe_df$Predicted <- gxe_df$Env_Effect + gxe_df$G + gxe_df$GxE_Effect

  gxe_df[, c("Genotype", "Environment", "G", "GxE_Effect", "SE", "PEV", "Env_Effect", "Predicted")]
}
