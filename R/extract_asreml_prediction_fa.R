extract_asreml_prediction <- function(model,
                                      genotype_term = "GID",
                                      environment_levels,
                                      kernel_ids = NULL) {
  #random_coefs <- coef(model)$random
  random_coefs <- summary(model, coef=TRUE)$coef.random

  fixed_coefs  <- coef(model)$fixed

  # 1. Intercept
  intercept <- if ("(Intercept)" %in% rownames(fixed_coefs)) fixed_coefs["(Intercept)", "effect"] else 0

  # 2. Extract all G×E terms from FA components across kernels
  gxe_rows <- grep(paste0("^vm\\(", genotype_term, ","), rownames(random_coefs), value = TRUE)
  gxe_rows <- gxe_rows[grepl(":fa\\(", gxe_rows)]  # only G×E FA terms

  if (length(gxe_rows) == 0) stop("No FA-based G×E terms found.")

  has_se <- "std.error" %in% colnames(random_coefs)

  gxe_df <- data.frame(
    Term = gxe_rows,
    GxE_Effect = random_coefs[gxe_rows, "solution"],
    SE = if (has_se) random_coefs[gxe_rows, "std.error"] else NA_real_
  )

  # Extract Genotype and Environment
  gxe_df$Genotype <- sub(":fa\\(.*", "", sub(paste0("vm\\(", genotype_term, ",.*?\\)_"), "", gxe_df$Term))
  gxe_df$Environment <- sub(".*:fa\\(.*?\\)_", "", gxe_df$Term)

  # Clean up environment matching
  environment_levels <- trimws(as.character(environment_levels))
  gxe_df$Environment <- trimws(as.character(gxe_df$Environment))
  gxe_df <- gxe_df[gxe_df$Environment %in% environment_levels, ]

  if (nrow(gxe_df) == 0) {
    stop("No valid GxE rows remain after filtering. Check environment_levels for mismatches.")
  }
  # Sum G×E over multiple kernels
  if ("SE" %in% colnames(gxe_df) && any(!is.na(gxe_df$SE))) {
    gxe_df <- aggregate(cbind(GxE_Effect, SE) ~ Genotype + Environment, data = gxe_df, FUN = function(x) sum(x, na.rm = TRUE))
  } else {
    gxe_df <- aggregate(GxE_Effect ~ Genotype + Environment, data = gxe_df, FUN = function(x) sum(x, na.rm = TRUE))
    gxe_df$SE <- NA_real_
  }

  gxe_df$PEV <- if (!is.na(gxe_df$SE[1])) gxe_df$SE^2 else NA_real_

  # 3. Extract G main effects from multiple kernels
  g_rows <- grep(paste0("^vm\\(", genotype_term, ","), rownames(random_coefs), value = TRUE)
  g_rows <- g_rows[!grepl(":fa\\(", g_rows)]  # exclude G×E terms

  has_G <- length(g_rows) > 0
  if (has_G) {
    g_df <- data.frame(
      Term = g_rows,
      G_value = random_coefs[g_rows, "solution"]
    )
    g_df$Genotype <- sub(".*\\)_", "", g_df$Term)

    # Sum across kernels
    g_df <- aggregate(G_value ~ Genotype, g_df, sum)
    colnames(g_df)[2] <- "G"
  } else {
    g_df <- data.frame(Genotype = unique(gxe_df$Genotype), G = 0)
  }

  gxe_df <- merge(gxe_df, g_df, by = "Genotype", all.x = TRUE)

  # 4. Environment effect (fixed or random)
  env_in_fixed <- any(grepl("^Env", rownames(fixed_coefs)))
  env_in_random <- any(grepl("^Env_", rownames(random_coefs)))

  if (env_in_fixed) {
    env_blue_terms <- grep("^Env_", rownames(fixed_coefs), value = TRUE)
    env_blue_vals <- fixed_coefs[env_blue_terms, "effect"]
    env_names <- sub("^Env_", "", env_blue_terms)
    fixed_env_map <- setNames(env_blue_vals, env_names)

    all_envs <- unique(environment_levels)
    missing_envs <- setdiff(all_envs, names(fixed_env_map))
    if (length(missing_envs) == 1) {
      baseline_env <- missing_envs
    } else {
      baseline_env <- names(fixed_env_map)[which(abs(fixed_env_map) < 1e-8)]
      if (length(baseline_env) != 1) stop("Unable to confidently identify baseline environment.")
      fixed_env_map <- fixed_env_map[setdiff(names(fixed_env_map), baseline_env)]
    }

    env_effects <- setNames(numeric(length(all_envs)), all_envs)
    env_effects[baseline_env] <- intercept
    env_effects[names(fixed_env_map)] <- intercept + fixed_env_map
    gxe_df$Env_Effect <- env_effects[gxe_df$Environment]

  } else if (env_in_random) {
    env_blups <- random_coefs[grep("^Env_", rownames(random_coefs)), "effect"]
    env_names <- gsub("^Env_", "", rownames(random_coefs)[grep("^Env_", rownames(random_coefs))])
    env_effects <- setNames(env_blups, env_names)
    gxe_df$Env_Effect <- env_effects[gxe_df$Environment] + intercept
  } else {
    gxe_df$Env_Effect <- intercept
  }

  # 5. Final predicted value
  gxe_df$Predicted <- gxe_df$Env_Effect + gxe_df$G + gxe_df$GxE_Effect

  # Return
  gxe_df <- gxe_df[, c("Genotype", "Environment", "G", "GxE_Effect", "SE", "PEV", "Env_Effect", "Predicted")]
  return(gxe_df)
}







# # ASReml model with multiple kernels + FA or CORGH structure
# resultt <- extract_asreml_prediction(
#   model = res_env,
#   genotype_term = "GID",  # or "Genotype" if that’s your term
#   environment_levels = unique(pheno2$Env),
#   kernel_ids = c("grm")  # names in Gu
# )
#
# head(result)
