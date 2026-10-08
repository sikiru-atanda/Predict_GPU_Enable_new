#' Remove objects from the global environment
#'
#' @param var_names Character vector of object names to remove.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
remove_from_global <- function(var_names) {
  for (var_name in var_names) {
    if (exists(var_name, envir = .GlobalEnv)) {
      rm(list = var_name, envir = .GlobalEnv)
    }
  }
}

#' Calculate ASReml coefficients from EBV output
#'
#' @param gmatrix Relationship or design matrix used for coefficient recovery.
#' @param ebv Estimated breeding values.
#' @param heter_groups Optional environment/grouping column name.
#' @param heter_grp Optional environment/grouping values.
#' @param gid_name Genotype identifier column name.
#'
#' @return A data frame of recovered coefficients.
#' @export
cal_coeff_asreml <- function(gmatrix,
                             ebv,
                             heter_groups,
                             heter_grp,
                             gid_name) {
  gmatrix <- as.matrix(gmatrix)
  storage.mode(gmatrix) <- "double"
  ebv <- as.numeric(ebv)
  if (!nrow(gmatrix) || !ncol(gmatrix) || length(ebv) != nrow(gmatrix)) {
    stop(
      "ASReml coefficient recovery requires one finite EBV per kernel row.",
      call. = FALSE
    )
  }
  if (any(!is.finite(gmatrix)) || any(!is.finite(ebv))) {
    stop("ASReml coefficient recovery received non-finite inputs.", call. = FALSE)
  }

  # Recover b from K b = u. Kernel matrices are commonly positive
  # semidefinite, so an ordinary inverse (and especially elementwise `*` in
  # place of matrix multiplication) is not a valid or stable solution.
  coeff <- as.numeric(MASS::ginv(gmatrix) %*% ebv)

  if (!is.null(heter_groups)) {
    row_aligned <- length(gid_name) == length(coeff) && length(heter_grp) == length(coeff)
    coeff <- data.frame(
      # Row-aligned labels (one genotype/environment per kernel row) when the
      # caller supplies them; otherwise the complete-grid labelling.
      x_variables = if (row_aligned) gid_name else rep(gid_name, length(unique(heter_grp))),
      Env = if (row_aligned) heter_grp else rep(unique(heter_grp), each = length(gid_name)),
      coeff = coeff,
      stringsAsFactors = FALSE
    )
    names(coeff)[2] <- heter_groups
  } else {
    coeff <- data.frame(
      x_variables = rownames(gmatrix),
      coeff = coeff,
      stringsAsFactors = FALSE
    )
  }

  coeff
}

asreml_kernel_output_names <- function(gmatrix = NULL,
                                       omic1_kernel = NULL,
                                       omic2_kernel = NULL,
                                       omic3_kernel = NULL,
                                       kernel_list = NULL) {
  output_names <- character()
  append_role <- function(value, role) {
    if (is.null(value)) return()
    count <- if (gp_is_kernel_list(value)) length(value) else 1L
    value_names <- if (gp_is_kernel_list(value)) names(value) else NULL
    if (count == 1L) {
      output_names <<- c(output_names, role)
    } else {
      if (is.null(value_names) || any(!nzchar(value_names))) {
        value_names <- paste0("kernel", seq_len(count))
      }
      output_names <<- c(
        output_names,
        paste(role, gp_sanitize_kernel_name(value_names), sep = "_")
      )
    }
  }

  append_role(gmatrix, "gmatrix")
  append_role(omic1_kernel, "omic1_kernel")
  append_role(omic2_kernel, "omic2_kernel")
  append_role(omic3_kernel, "omic3_kernel")

  if (gp_is_kernel_list(kernel_list)) {
    canonical <- c(
      "gmatrix_model_ready", "omic1_kernel_model_ready",
      "omic2_kernel_model_ready", "omic3_kernel_model_ready"
    )
    extras <- kernel_list[setdiff(names(kernel_list), canonical)]
    extras <- extras[vapply(extras, gp_is_kernel_matrix, logical(1))]
    if (length(extras)) {
      extra_names <- names(extras)
      if (is.null(extra_names) || any(!nzchar(extra_names))) {
        extra_names <- paste0("kernel", seq_along(extras))
      }
      output_names <- c(output_names, gp_sanitize_kernel_name(extra_names))
    }
  }

  make.unique(output_names, sep = "_")
}

asreml_single_kernel_varcomp_rows <- function(row_names,
                                              gen_name,
                                              kernel_name) {
  compact <- gsub("\\s+", "", as.character(row_names))
  target <- paste0("vm(", gen_name, ",", kernel_name, ")")
  which(compact == target)
}

asreml_build_options <- function(workspace = NULL,
                                 pworkspace = NULL,
                                 maxit = 50,
                                 trace = FALSE) {
  opts <- list(trace = isTRUE(trace))
  if (!is.null(maxit) && length(maxit) > 0L && !is.na(maxit[1L])) {
    opts$maxit <- maxit[1L]
  }
  if (!is.null(workspace) && length(workspace) > 0L && !is.na(workspace[1L])) {
    opts$workspace <- workspace[1L]
  }
  if (!is.null(pworkspace) && length(pworkspace) > 0L && !is.na(pworkspace[1L])) {
    opts$pworkspace <- pworkspace[1L]
  }
  opts
}

# ASReml workspace in 8-byte words. Accepts numbers (words) and strings such
# as "800mb" / "4gb"; NULL means ASReml's default of 1e8 words.
gp_asreml_workspace_words <- function(workspace) {
  if (is.null(workspace) || !length(workspace) || is.na(workspace[1L])) return(1e8)
  w <- workspace[1L]
  if (is.numeric(w)) return(as.numeric(w))
  m <- regmatches(tolower(w), regexec("^\\s*([0-9.]+)\\s*(kb|mb|gb)?\\s*$", tolower(w)))[[1L]]
  if (length(m) < 2L) return(1e8)
  bytes <- as.numeric(m[2L]) * switch(m[3L], kb = 1024, mb = 1024^2, gb = 1024^3, 1)
  bytes / 8
}

# Next workspace for a retry: 4x larger, capped at 60% of available memory;
# NULL when no larger workspace fits.
gp_asreml_next_workspace_words <- function(current_words) {
  avail <- tryCatch(as.numeric(ps::ps_system_memory()[["avail"]]), error = function(e) NA_real_)
  cap <- if (is.finite(avail)) 0.6 * avail / 8 else 4 * current_words
  nxt <- min(4 * current_words, cap)
  if (!is.finite(nxt) || nxt <= current_words * 1.05) return(NULL)
  nxt
}

asreml_has_arg_value <- function(x) {
  !is.null(x) && length(x) > 0L && !is.na(x[1L])
}

asreml_validate_coefficient_fallback_model <- function(model,
                                                       context = "ASReml coefficient fallback") {
  if (!is.null(model$converge) && length(model$converge) > 0L && !isTRUE(model$converge[[1L]])) {
    stop(context, " is unsafe because the fitted ASReml model did not converge.", call. = FALSE)
  }
  if (!is.null(model$ifault) && length(model$ifault) > 0L) {
    ifault <- suppressWarnings(as.integer(model$ifault[[1L]]))
    if (is.finite(ifault) && ifault != 0L) {
      stop(context, " is unsafe because the fitted ASReml model has ifault = ", ifault, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# asreml.options() is session-global. The previous values of the options set
# here are restored when `restore_frame` (by default the caller) exits, so a
# PredictProR fit does not change the user's own ASReml settings.
asreml_apply_options <- function(workspace = NULL,
                                 pworkspace = NULL,
                                 maxit = 50,
                                 trace = FALSE,
                                 engine = "asreml",
                                 restore_frame = parent.frame()) {
  opts <- asreml_build_options(
    workspace = workspace,
    pworkspace = pworkspace,
    maxit = maxit,
    trace = trace
  )
  if (!requireNamespace(engine, quietly = TRUE)) {
    stop("You need to install asreml-R to use ASReml-R options.", call. = FALSE)
  }
  options_fn <- getFromNamespace("asreml.options", engine)
  previous <- tryCatch(options_fn()[names(opts)], error = function(e) NULL)
  if (is.list(previous) && length(previous) && !is.null(restore_frame)) {
    previous <- previous[!vapply(previous, is.null, logical(1))]
    restore <- function() try(do.call(options_fn, previous), silent = TRUE)
    # after = FALSE: when one caller applies options more than once (e.g. a
    # workspace retry), the first-registered restore -- the user's own
    # values -- runs last.
    do.call(base::on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE), envir = restore_frame)
  }
  do.call(options_fn, opts)
  invisible(opts)
}

asreml_predict_pvals <- function(mod,
                                 classify,
                                 workspace = NULL,
                                 pworkspace = NULL) {
  predict_impl <- getOption("PredictProR.asreml_predict_impl", asreml::predict.asreml)
  has_workspace <- asreml_has_arg_value(workspace)
  has_pworkspace <- asreml_has_arg_value(pworkspace)

  pred <- if (has_workspace && has_pworkspace) {
    predict_impl(
      mod,
      classify = classify,
      sed = FALSE,
      workspace = workspace[1L],
      pworkspace = pworkspace[1L]
    )
  } else if (has_workspace) {
    predict_impl(
      mod,
      classify = classify,
      sed = FALSE,
      workspace = workspace[1L]
    )
  } else if (has_pworkspace) {
    predict_impl(
      mod,
      classify = classify,
      sed = FALSE,
      pworkspace = pworkspace[1L]
    )
  } else {
    predict_impl(mod, classify = classify, sed = FALSE)
  }
  pred$pvals
}

# Environment level of each ASReml random-coefficient row, parsed from its
# name: `Env_E1:vm(GID, K)_g1`, `vm(GID, K)_g1:Env_E1`, `fa(Env, 1)_E1:vm(...)`
# or `corgh(Env)_E1:vm(...)`. Genotype main-effect rows (`vm(GID, K)_g1`)
# return NA. ASReml orders coefficients by factor level, not by the order the
# environments appear in the data, so labels must not be assigned by position.
asreml_coef_row_env <- function(row_names, heter_groups) {
  if (is.null(row_names) || !length(row_names)) {
    return(character())
  }
  hx <- gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", heter_groups)
  plain_rx <- paste0("^", hx, "_")
  struct_rx <- paste0("^[A-Za-z]+\\(\\s*", hx, "[^)]*\\)_")
  vapply(strsplit(as.character(row_names), ":", fixed = TRUE), function(parts) {
    parts <- parts[!grepl("^vm\\(", parts)]
    for (q in parts) {
      if (grepl(plain_rx, q)) return(sub(plain_rx, "", q))
      if (grepl(struct_rx, q)) return(sub("^.*\\)_", "", q))
    }
    NA_character_
  }, character(1))
}

# Genotype-by-environment genetic values for the legacy BLUP-table fallback:
# per-environment rows are summed across kernels and each genotype's main
# effect (Env = NA) is added to every environment.
asreml_legacy_met_genetic_values <- function(df, gen_name, heter_groups, env_levels) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  gen <- as.character(df[[gen_name]])
  env <- as.character(df[[heter_groups]])
  blup <- as.numeric(df[["BLUP"]])
  main <- tapply(blup[is.na(env)], gen[is.na(env)], sum)
  gxe_rows <- !is.na(env)
  gxe <- if (any(gxe_rows)) {
    stats::aggregate(list(BLUP = blup[gxe_rows]), list(g = gen[gxe_rows], e = env[gxe_rows]), sum)
  } else {
    data.frame(g = character(), e = character(), BLUP = numeric())
  }
  grid <- expand.grid(g = unique(gen), e = as.character(env_levels), stringsAsFactors = FALSE)
  key <- paste(grid$g, grid$e, sep = "\r")
  gxe_val <- gxe$BLUP[match(key, paste(gxe$g, gxe$e, sep = "\r"))]
  main_val <- as.numeric(main[grid$g])
  has_any <- !is.na(gxe_val) | !is.na(main_val)
  value <- ifelse(is.na(gxe_val), 0, gxe_val) + ifelse(is.na(main_val), 0, main_val)
  out <- data.frame(grid$g, grid$e, value, stringsAsFactors = FALSE)[has_any, , drop = FALSE]
  names(out) <- c(gen_name, heter_groups, "BLUP")
  rownames(out) <- NULL
  out
}

asreml_reformat_met_ebv <- function(df, gen_name, heter_groups, heter_grp) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!gen_name %in% names(df)) {
    stop("Expected genotype column was not found in the ASReml coefficient table.", call. = FALSE)
  }

  se_col <- if ("Standard_error" %in% names(df)) "Standard_error" else NULL
  ignore_cols <- c(gen_name, heter_groups, se_col, "Prediction_error_variance", "Reliability", "z.ratio")
  value_candidates <- setdiff(names(df), ignore_cols)
  value_candidates <- value_candidates[vapply(df[value_candidates], is.numeric, logical(1))]
  if (!length(value_candidates)) {
    stop("Could not identify the ASReml BLUP column for MET output.", call. = FALSE)
  }

  blup_col <- if ("solution" %in% value_candidates) "solution" else value_candidates[[1]]
  if (is.null(heter_grp) || !length(heter_grp)) {
    env_values <- rep(NA_character_, nrow(df))
  } else {
    env_values <- asreml_coef_row_env(rownames(df), heter_groups)
    env_values[!env_values %in% as.character(heter_grp)] <- NA_character_
  }

  out <- data.frame(
    stringsAsFactors = FALSE,
    genotype = as.character(df[[gen_name]]),
    env = env_values,
    BLUP = as.numeric(df[[blup_col]]),
    Standard_error = if (!is.null(se_col)) as.numeric(df[[se_col]]) else rep(NA_real_, nrow(df))
  )
  names(out)[1:2] <- c(gen_name, heter_groups)
  out[["Prediction_error_variance"]] <- out[["Standard_error"]]^2
  out
}

# Helper: pull the global intercept and per-environment fixed offsets from
# an asreml model so a BLUP-based fallback prediction lands on the original
# phenotype scale (intercept + env_fixed + BLUP), not the deviation scale.
#
# Returns either:
#   * a named numeric vector keyed by environment level (when env_levels is
#     supplied and heter_groups is non-empty); each entry is the asreml
#     intercept + fixed env effect for that level; or
#   * a scalar (when env_levels is NULL) representing the ACROSS-env mean
#     offset to apply to a single-loc gid prediction. This is the across-env
#     mean of the per-env offsets (so a BLUP averaged across envs lands at
#     the overall phenotype scale, not the asreml reference-level intercept).
#
# Falls back to mean(pheno_data[[response]], na.rm = TRUE) when the model
# is unavailable or no fixed-coefficient table can be retrieved.
asreml_extract_fixed_offsets <- function(model = NULL,
                                         heter_groups = NULL,
                                         env_levels = NULL,
                                         pheno_data = NULL,
                                         response = NULL) {
  proxy <- 0
  if (!is.null(pheno_data) && !is.null(response) && response %in% names(pheno_data)) {
    p <- suppressWarnings(mean(as.numeric(pheno_data[[response]]), na.rm = TRUE))
    if (is.finite(p)) proxy <- p
  }
  fixed_coefs <- if (is.null(model)) NULL else
    tryCatch(stats::coef(model)$fixed, error = function(e) NULL)
  fixed_col <- if (!is.null(fixed_coefs))
    intersect(c("effect", "solution"), colnames(fixed_coefs))[1L] else NA
  if (is.null(fixed_coefs) || !nrow(fixed_coefs) ||
      is.na(fixed_col) || is.null(fixed_col)) {
    if (is.null(env_levels)) return(proxy)
    return(stats::setNames(rep(proxy, length(env_levels)), env_levels))
  }

  intercept <- if ("(Intercept)" %in% rownames(fixed_coefs)) {
    as.numeric(fixed_coefs["(Intercept)", fixed_col])
  } else 0
  if (!is.finite(intercept)) intercept <- 0

  # Gather per-env offsets when heter_groups is given and a fixed term keyed
  # by that column exists. Treatment contrasts in asreml put the reference
  # level's mean into (Intercept) and the other levels' deltas into the
  # "{heter}_{level}" rows -- this captures both.
  per_env <- NULL
  heter_chr <- as.character(heter_groups %||% "")
  if (nzchar(heter_chr)) {
    heter_rx <- gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", heter_chr)
    fixed_env_rows <- grep(paste0("^", heter_rx, "_"), rownames(fixed_coefs), value = TRUE)
    if (length(fixed_env_rows)) {
      fixed_env_map <- stats::setNames(
        as.numeric(fixed_coefs[fixed_env_rows, fixed_col]),
        sub(paste0("^", heter_rx, "_"), "", fixed_env_rows)
      )
      # If we know the env_levels, build the full per-env vector.
      # If env_levels is NULL but pheno_data has the column, derive it.
      derived_levels <- if (!is.null(env_levels)) as.character(env_levels) else if (
        !is.null(pheno_data) && heter_chr %in% names(pheno_data)
      ) unique(as.character(pheno_data[[heter_chr]])) else names(fixed_env_map)
      per_env <- stats::setNames(rep(intercept, length(derived_levels)), derived_levels)
      common <- intersect(names(fixed_env_map), derived_levels)
      if (length(common)) per_env[common] <- intercept + fixed_env_map[common]
    }
  }

  if (!is.null(env_levels)) {
    if (!is.null(per_env)) {
      out <- stats::setNames(rep(intercept, length(env_levels)), as.character(env_levels))
      common <- intersect(names(per_env), as.character(env_levels))
      if (length(common)) out[common] <- per_env[common]
      return(out)
    }
    return(stats::setNames(rep(intercept, length(env_levels)), as.character(env_levels)))
  }

  # env_levels NULL -> caller wants a SCALAR across-env offset. If we have
  # per-env offsets, average them so a BLUP averaged across envs lands at
  # the across-env phenotype mean. Otherwise return the intercept alone.
  if (!is.null(per_env) && length(per_env)) {
    return(mean(per_env, na.rm = TRUE))
  }
  intercept
}

asreml_prediction_from_ebv_tables <- function(ebv_tables,
                                              gen_name,
                                              pheno_data,
                                              response,
                                              heter_groups = NULL,
                                              model = NULL,
                                              model_heter_groups = heter_groups) {
  if (is.null(ebv_tables) || !length(ebv_tables)) {
    return(NULL)
  }

  combined_df <- dplyr::bind_rows(ebv_tables)
  if (!nrow(combined_df)) {
    return(NULL)
  }

  # Reproject BLUPs onto the original phenotype scale by adding the global
  # intercept (and per-env fixed offsets when present). Without this, the
  # legacy BLUP-mean fallback returns values centered around zero instead of
  # around the response mean, which collapses cor()/MSE against observed.
  # `model_heter_groups` lets a single-loc summary call (output heter=NULL)
  # still consult the MET env column on the model's fixed side so the
  # offset is the across-env mean, not the asreml reference-level intercept.
  # Genetic values are the SUM of the kernels' BLUPs (and, for MET, of the
  # genotype main effect and its environment-specific effect); averaging them
  # shrank multi-kernel and main + GxE predictions.
  env_col <- if (!is.null(heter_groups) && heter_groups %in% names(combined_df)) {
    heter_groups
  } else if (!is.null(model_heter_groups) && model_heter_groups %in% names(combined_df)) {
    model_heter_groups
  } else {
    NULL
  }
  if (!is.null(env_col)) {
    env_levels <- if (!is.null(pheno_data) && env_col %in% names(pheno_data)) {
      unique(as.character(pheno_data[[env_col]]))
    } else {
      unique(stats::na.omit(as.character(combined_df[[env_col]])))
    }
    cell <- asreml_legacy_met_genetic_values(combined_df, gen_name, env_col, env_levels)
    offsets <- asreml_extract_fixed_offsets(
      model = model, heter_groups = model_heter_groups %||% env_col,
      env_levels = env_levels,
      pheno_data = pheno_data, response = response
    )
    cell$Predicted_value <- cell$BLUP + as.numeric(offsets[as.character(cell[[env_col]])])
  }
  if (!is.null(heter_groups) && !is.null(env_col)) {
    pred <- cell[, c(gen_name, env_col, "Predicted_value"), drop = FALSE]
    names(pred)[2] <- heter_groups
  } else if (!is.null(env_col)) {
    # Across-environment summary of a MET model: mean genetic value over envs.
    pred <- stats::aggregate(
      list(Predicted_value = cell$Predicted_value),
      stats::setNames(list(cell[[gen_name]]), gen_name),
      mean
    )
  } else {
    pred <- stats::aggregate(
      list(BLUP = as.numeric(combined_df$BLUP)),
      stats::setNames(list(as.character(combined_df[[gen_name]])), gen_name),
      sum
    )
    intercept <- asreml_extract_fixed_offsets(
      model = model,
      heter_groups = model_heter_groups,
      env_levels = NULL,
      pheno_data = pheno_data, response = response
    )
    pred$Predicted_value <- pred$BLUP + as.numeric(intercept)
    pred$BLUP <- NULL
  }
  pred$Standard_error <- NA_real_
  pred$Prediction_error_variance <- NA_real_

  pred <- as.data.frame(pred, stringsAsFactors = FALSE)
  pred <- asreml_add_train_test_labels(
    predicted_df = pred,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  asreml_mark_prediction_uncertainty_unavailable(
    pred,
    reason = paste(
      "Unavailable: ASReml predict failed;",
      "BLUP-table fallback returns point predictions only"
    )
  )
}

asreml_extract_ebv_tables <- function(blup,
                                      names_in_inv_list,
                                      gen_name,
                                      heter_groups = NULL,
                                      heter_grp = NULL,
                                      var_cov_str = NULL,
                                      inter_gen_pos = NULL) {
  ebv_tables <- list()

  for (b in seq_along(names_in_inv_list)) {
    ebv_tables[[names_in_inv_list[b]]] <- blup[grep(paste(names_in_inv_list[b], "\\)", sep = ""), rownames(blup)), ]
  }

  is_met <- !is.null(inter_gen_pos)

  for (bb in seq_along(names_in_inv_list)) {
    ebv_df <- ebv_tables[[names_in_inv_list[bb]]]
    if (isTRUE(grepl("fa", var_cov_str)) || isTRUE(grepl("rr", var_cov_str))) {
      ebv_df <- ebv_df[
        !rownames(ebv_df) %in% rownames(ebv_df[grep("Comp", rownames(ebv_df)), , drop = FALSE]),
        ,
        drop = FALSE
      ]
      ebv_df <- as.data.frame(ebv_df)
      ebv_df[, gen_name] <- as.character(stringr::str_split_fixed(rownames(ebv_df), "\\)_", 3)[, 3])
    } else {
      ebv_df <- as.data.frame(ebv_df)
      ebv_df[, gen_name] <- as.character(stringr::str_split_fixed(rownames(ebv_df), "\\)_", 3)[, 2])
    }

    if (is_met) {
      ebv_df <- asreml_reformat_met_ebv(
        df = ebv_df,
        gen_name = gen_name,
        heter_groups = heter_groups,
        heter_grp = heter_grp
      )
    } else {
      rownames(ebv_df) <- NULL
      ebv_df <- ebv_df[, c(4, 1:2)]
      colnames(ebv_df)[1:2] <- c(gen_name, "BLUP")
    }

    rownames(ebv_df) <- NULL
    ebv_tables[[names_in_inv_list[bb]]] <- ebv_df
  }

  ebv_tables
}

asreml_apply_met_reliability <- function(ebv_df, heter_groups, heter_grp, variance_by_env) {
  n <- nrow(ebv_df)
  pev <- if ("PEV" %in% names(ebv_df)) {
    suppressWarnings(as.numeric(ebv_df[["PEV"]]))
  } else if ("Prediction_error_variance" %in% names(ebv_df)) {
    suppressWarnings(as.numeric(ebv_df[["Prediction_error_variance"]]))
  } else if ("Standard_error" %in% names(ebv_df)) {
    suppressWarnings(as.numeric(ebv_df[["Standard_error"]]))^2
  } else {
    rep(NA_real_, n)
  }
  if (length(pev) != n) {
    pev <- rep(NA_real_, n)
  }

  ebv_df[["Reliability"]] <- rep(NA_real_, n)
  ebv_df[["Reliability_variance_input"]] <- pev
  ebv_df[["Reliability_reference_variance"]] <- rep(NA_real_, n)
  ebv_df[["Reliability_basis"]] <- rep(
    "unavailable: complete positive ASReml environment genetic variance is required",
    n
  )
  ebv_df[["Reliability_remarks"]] <- rep(
    "Reliability unavailable because environment genetic variance is incomplete or non-positive",
    n
  )

  if (!heter_groups %in% names(ebv_df) ||
      !asreml_met_reliability_variance_valid(variance_by_env, heter_grp)) {
    return(ebv_df)
  }

  variance <- suppressWarnings(as.numeric(variance_by_env))
  names(variance) <- as.character(heter_grp)
  reference <- unname(variance[match(
    as.character(ebv_df[[heter_groups]]),
    as.character(heter_grp)
  )])
  valid <- is.finite(pev) & pev >= 0 & is.finite(reference) & reference > 0
  reliability <- rep(NA_real_, n)
  reliability[valid] <- pmax(0, pmin(1, 1 - pev[valid] / reference[valid]))

  ebv_df[["Reliability"]] <- reliability
  ebv_df[["Reliability_reference_variance"]] <- reference
  ebv_df[["Reliability_basis"]] <- rep(
    paste(
      "model-specific genetic variance;",
      "Reliability = max(0, min(1, 1 - PEV / environment genetic variance))"
    ),
    n
  )
  ebv_df[["Reliability_remarks"]] <- ifelse(
    !is.finite(reliability),
    "Reliability unavailable because PEV is not finite",
    ifelse(reliability >= 0.6, "Reliable",
           ifelse(reliability <= 0.2, "Unreliable", "Acceptable"))
  )

  ebv_df
}

asreml_met_reliability_variance_valid <- function(variance_by_env, environment_levels) {
  value <- suppressWarnings(as.numeric(variance_by_env))
  length(value) == length(environment_levels) &&
    length(value) > 0L &&
    all(is.finite(value)) &&
    all(value > 0)
}

asreml_add_uncertainty_fields <- function(predicted_df,
                                          pheno_data,
                                          response,
                                          high_reliability_thres = 0.6,
                                          low_reliability_thres = 0.2) {
  if (is.null(predicted_df) || !is.data.frame(predicted_df) || !nrow(predicted_df)) {
    return(predicted_df)
  }

  if (!"Prediction_error_variance" %in% names(predicted_df) &&
      "Standard_error" %in% names(predicted_df)) {
    predicted_df[["Prediction_error_variance"]] <- suppressWarnings(
      as.numeric(predicted_df[["Standard_error"]])
    )^2
  }
  if (!"PEV" %in% names(predicted_df) &&
      "Prediction_error_variance" %in% names(predicted_df)) {
    predicted_df[["PEV"]] <- suppressWarnings(as.numeric(predicted_df[["Prediction_error_variance"]]))
  }
  if (!"Prediction_error_variance" %in% names(predicted_df) &&
      "PEV" %in% names(predicted_df)) {
    predicted_df[["Prediction_error_variance"]] <- suppressWarnings(as.numeric(predicted_df[["PEV"]]))
  }

  y <- suppressWarnings(as.numeric(pheno_data[[response]]))
  reference_variance <- suppressWarnings(stats::var(y[is.finite(y)], na.rm = TRUE))
  if (!is.finite(reference_variance) || reference_variance <= 0) {
    reference_variance <- NA_real_
  }
  if (!"Reliability_reference_variance" %in% names(predicted_df)) {
    predicted_df[["Reliability_reference_variance"]] <- rep(reference_variance, nrow(predicted_df))
  }

  pev <- if ("PEV" %in% names(predicted_df)) {
    suppressWarnings(as.numeric(predicted_df[["PEV"]]))
  } else {
    rep(NA_real_, nrow(predicted_df))
  }
  ref <- suppressWarnings(as.numeric(predicted_df[["Reliability_reference_variance"]]))
  reliability <- rep(NA_real_, nrow(predicted_df))
  valid <- is.finite(pev) & is.finite(ref) & ref > 0
  reliability[valid] <- pmax(0, pmin(1, 1 - pev[valid] / ref[valid]))
  if (!"Reliability" %in% names(predicted_df)) {
    predicted_df[["Reliability"]] <- reliability
  } else {
    old <- suppressWarnings(as.numeric(predicted_df[["Reliability"]]))
    replace_idx <- !is.finite(old) & is.finite(reliability)
    old[replace_idx] <- reliability[replace_idx]
    predicted_df[["Reliability"]] <- old
  }
  if (!"Reliability_remarks" %in% names(predicted_df)) {
    predicted_df[["Reliability_remarks"]] <- ifelse(
      is.na(predicted_df[["Reliability"]]),
      NA_character_,
      ifelse(
        predicted_df[["Reliability"]] >= high_reliability_thres,
        "Reliable",
        ifelse(predicted_df[["Reliability"]] <= low_reliability_thres, "Unreliable", "Acceptable")
      )
    )
  }
  if (!"Reliability_percentage" %in% names(predicted_df)) {
    rel <- suppressWarnings(as.numeric(predicted_df[["Reliability"]]))
    predicted_df[["Reliability_percentage"]] <- if (any(is.finite(rel))) {
      rep(mean(rel[is.finite(rel)] > low_reliability_thres) * 100, nrow(predicted_df))
    } else {
      rep(NA_real_, nrow(predicted_df))
    }
  }
  if (!"Reliability_basis" %in% names(predicted_df)) {
    predicted_df[["Reliability_basis"]] <- rep("observed_response_variance", nrow(predicted_df))
  }
  predicted_df
}

asreml_mark_prediction_uncertainty_unavailable <- function(predicted_df,
                                                           reason = paste(
                                                             "Unavailable:",
                                                             "ASReml predict failed;",
                                                             "coefficient fallback returns point predictions only"
                                                           )) {
  if (is.null(predicted_df) || !is.data.frame(predicted_df) || !nrow(predicted_df)) {
    return(predicted_df)
  }

  numeric_cols <- c(
    "Standard_error",
    "Prediction_error_variance",
    "PEV",
    "Reliability",
    "Reliability_reference_variance",
    "Reliability_percentage",
    "lower_bound",
    "upper_bound",
    "Uncertainty"
  )
  required_numeric_cols <- c("Standard_error", "Prediction_error_variance", "PEV", "Reliability")
  for (col in unique(c(required_numeric_cols, intersect(numeric_cols, names(predicted_df))))) {
    predicted_df[[col]] <- rep(NA_real_, nrow(predicted_df))
  }

  predicted_df[["Reliability_remarks"]] <- rep(reason, nrow(predicted_df))
  predicted_df[["Uncertainty_remarks"]] <- rep(reason, nrow(predicted_df))
  predicted_df[["Reliability_basis"]] <- rep("unavailable_asreml_predict_failed", nrow(predicted_df))
  predicted_df[["Prediction_uncertainty_source"]] <- rep("unavailable_asreml_predict_failed", nrow(predicted_df))
  predicted_df
}

asreml_prediction_value_column <- function(predicted_df) {
  if (!is.null(predicted_df) && "Predicted_value" %in% names(predicted_df)) {
    return("Predicted_value")
  }
  if (!is.null(predicted_df) && "BLUP" %in% names(predicted_df)) {
    return("BLUP")
  }
  NULL
}

asreml_status_from_test_flags <- function(x) {
  if (length(x) == 0) {
    return(NA_character_)
  }
  if (all(x)) {
    return("Test")
  }
  if (any(x)) {
    return("Mixed")
  }
  "Train"
}

asreml_add_train_test_labels <- function(predicted_df,
                                         pheno_data,
                                         response,
                                         gen_name,
                                         heter_groups = NULL) {
  if (is.null(predicted_df) || !nrow(predicted_df) || !gen_name %in% names(predicted_df)) {
    return(predicted_df)
  }

  test_flags <- is.na(pheno_data[[response]])
  if (!any(test_flags)) {
    predicted_df[, "Train_Test_Label"] <- "Train"
    return(asreml_add_uncertainty_fields(
      predicted_df = predicted_df,
      pheno_data = pheno_data,
      response = response
    ))
  }

  if (!is.null(heter_groups) &&
      heter_groups %in% names(predicted_df) &&
      heter_groups %in% names(pheno_data)) {
    obs_keys <- data.frame(
      key = paste(pheno_data[[gen_name]], pheno_data[[heter_groups]], sep = "\r"),
      is_test = test_flags,
      stringsAsFactors = FALSE
    )
    split_flags <- split(obs_keys$is_test, obs_keys$key)
    status_lookup <- vapply(split_flags, asreml_status_from_test_flags, character(1))
    pred_keys <- paste(predicted_df[[gen_name]], predicted_df[[heter_groups]], sep = "\r")
    predicted_df[, "Train_Test_Label"] <- unname(status_lookup[pred_keys])
  } else {
    obs_keys <- data.frame(
      key = pheno_data[[gen_name]],
      is_test = test_flags,
      stringsAsFactors = FALSE
    )
    split_flags <- split(obs_keys$is_test, obs_keys$key)
    status_lookup <- vapply(split_flags, asreml_status_from_test_flags, character(1))
    predicted_df[, "Train_Test_Label"] <- unname(status_lookup[as.character(predicted_df[[gen_name]])])
  }

  predicted_df[, "Train_Test_Label"][is.na(predicted_df[, "Train_Test_Label"])] <- "Train"
  asreml_add_uncertainty_fields(
    predicted_df = predicted_df,
    pheno_data = pheno_data,
    response = response
  )
}

asreml_build_residuals <- function(predicted_df,
                                   pheno_data,
                                   response,
                                   gen_name,
                                   heter_groups = NULL) {
  pred_col <- asreml_prediction_value_column(predicted_df)
  if (is.null(predicted_df) || !nrow(predicted_df) || is.null(pred_col) || !gen_name %in% names(predicted_df)) {
    return(NULL)
  }

  observed <- pheno_data[!is.na(pheno_data[[response]]), , drop = FALSE]
  if (!nrow(observed)) {
    return(NULL)
  }

  predicted_df[[".pred_order"]] <- seq_len(nrow(predicted_df))

  if (!is.null(heter_groups) &&
      heter_groups %in% names(predicted_df) &&
      heter_groups %in% names(observed)) {
    observed_small <- observed[, c(gen_name, heter_groups, response), drop = FALSE]
    merged <- merge(
      predicted_df,
      observed_small,
      by = c(gen_name, heter_groups),
      all = FALSE,
      sort = FALSE
    )
    if (!nrow(merged)) {
      return(NULL)
    }
    merged <- merged[order(merged[[".pred_order"]]), , drop = FALSE]
    out <- data.frame(
      stringsAsFactors = FALSE,
      merged[, c(gen_name, heter_groups), drop = FALSE],
      Predicted_value = merged[[pred_col]],
      Residual_value = merged[[response]] - merged[[pred_col]]
    )
    return(out)
  }

  observed_small <- stats::aggregate(
    observed[[response]],
    by = setNames(list(observed[[gen_name]]), gen_name),
    FUN = mean,
    na.rm = TRUE
  )
  names(observed_small)[2] <- response
  merged <- merge(predicted_df, observed_small, by = gen_name, all = FALSE, sort = FALSE)
  if (!nrow(merged)) {
    return(NULL)
  }
  merged <- merged[order(merged[[".pred_order"]]), , drop = FALSE]
  data.frame(
    stringsAsFactors = FALSE,
    merged[, gen_name, drop = FALSE],
    Predicted_value = merged[[pred_col]],
    Residual_value = merged[[response]] - merged[[pred_col]]
  )
}

asreml_standardize_extracted_predictions <- function(extracted_df,
                                                     pheno_data,
                                                     response,
                                                     gen_name,
                                                     heter_groups = NULL) {
  if (is.null(extracted_df) || !nrow(extracted_df)) {
    return(list(single_prediction = NULL, across_env_prediction = NULL))
  }

  if (!is.null(heter_groups) && "Environment" %in% names(extracted_df)) {
    across_env_predicted_value <- as.data.frame(extracted_df, stringsAsFactors = FALSE)
    names(across_env_predicted_value)[
      names(across_env_predicted_value) %in% c("Genotype", "Environment", "SE", "PEV", "Predicted")
    ] <- c(gen_name, heter_groups, "Standard_error", "Prediction_error_variance", "Predicted_value")

    across_env_predicted_value <- asreml_add_train_test_labels(
      predicted_df = across_env_predicted_value,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
    across_env_predicted_value <- asreml_mark_prediction_uncertainty_unavailable(across_env_predicted_value)

    single_prediction <- across_env_predicted_value |>
      dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
      dplyr::summarise(
        Predicted_value = mean(Predicted_value, na.rm = TRUE),
        Standard_error = NA_real_,
        Prediction_error_variance = NA_real_,
        .groups = "drop"
      )

    single_prediction <- asreml_add_train_test_labels(
      predicted_df = single_prediction,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
    single_prediction <- asreml_mark_prediction_uncertainty_unavailable(single_prediction)

    return(list(
      single_prediction = single_prediction,
      across_env_prediction = across_env_predicted_value
    ))
  }

  single_prediction <- as.data.frame(extracted_df, stringsAsFactors = FALSE)
  names(single_prediction)[
    names(single_prediction) %in% c("Genotype", "SE", "PEV", "Predicted")
  ] <- c(gen_name, "Standard_error", "Prediction_error_variance", "Predicted_value")

  single_prediction <- asreml_add_train_test_labels(
    predicted_df = single_prediction,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  single_prediction <- asreml_mark_prediction_uncertainty_unavailable(single_prediction)

  list(single_prediction = single_prediction, across_env_prediction = NULL)
}

asreml_run_with_model_data_context <- function(mod, pheno_data, expr) {
  data_sym <- NULL
  if (!is.null(mod) && !is.null(mod$call) && !is.null(mod$call$data)) {
    data_expr <- mod$call$data
    if (is.symbol(data_expr)) {
      data_sym <- as.character(data_expr)
    }
  }
  if (is.null(data_sym) || !nzchar(data_sym)) {
    data_sym <- "pheno_data"
  }

  had_existing <- exists(data_sym, envir = .GlobalEnv, inherits = FALSE)
  old_value <- if (had_existing) get(data_sym, envir = .GlobalEnv, inherits = FALSE) else NULL

  assign(data_sym, pheno_data, envir = .GlobalEnv)
  on.exit({
    if (had_existing) {
      assign(data_sym, old_value, envir = .GlobalEnv)
    } else if (exists(data_sym, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = data_sym, envir = .GlobalEnv)
    }
  }, add = TRUE)

  force(expr)
}

asreml_predict_or_extract <- function(mod,
                                      classify,
                                      pheno_data,
                                      response,
                                      gen_name,
                                      heter_groups = NULL,
                                      names_in_inv_list = NULL,
                                      workspace = NULL,
                                      pworkspace = NULL) {
  direct_prediction <- tryCatch({
    predicted_df <- asreml_run_with_model_data_context(
      mod = mod,
      pheno_data = pheno_data,
      expr = asreml_predict_pvals(
        mod = mod,
        classify = classify,
        workspace = workspace,
        pworkspace = pworkspace
      )
    )
    predicted_df <- predicted_df[, -ncol(predicted_df), drop = FALSE]
    colnames(predicted_df)[colnames(predicted_df) %in% c("predicted.value", "std.error")] <- c("Predicted_value", "Standard_error")
    predicted_df[, "Prediction_error_variance"] <- predicted_df[, "Standard_error"]^2
    predicted_df <- asreml_add_train_test_labels(
      predicted_df = predicted_df,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
    if (!is.null(heter_groups) && heter_groups %in% names(predicted_df)) {
      single_prediction <- predicted_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
        dplyr::summarise(
          Predicted_value = mean(Predicted_value, na.rm = TRUE),
          Standard_error = mean(Standard_error, na.rm = TRUE),
          Prediction_error_variance = mean(Prediction_error_variance, na.rm = TRUE),
          .groups = "drop"
        )
      single_prediction <- asreml_add_train_test_labels(
        predicted_df = as.data.frame(single_prediction, stringsAsFactors = FALSE),
        pheno_data = pheno_data,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      list(
        mode = "predict",
        prediction = single_prediction,
        across_env_prediction = predicted_df,
        error = NULL
      )
    } else {
      list(mode = "predict", prediction = predicted_df, across_env_prediction = NULL, error = NULL)
    }
  }, error = function(e) {
    list(mode = "predict_error", prediction = NULL, across_env_prediction = NULL, error = conditionMessage(e))
  })

  if (!is.null(direct_prediction$prediction)) {
    return(direct_prediction)
  }

  extracted_error <- NULL
  extracted_prediction <- tryCatch({
    extract_asreml_prediction(
      model = mod,
      genotype_term = gen_name,
      environment_levels = if (!is.null(heter_groups)) unique(pheno_data[[heter_groups]]) else NULL,
      heter_grp = heter_groups,
      kernel_ids = names_in_inv_list
    )
  }, error = function(e) {
    extracted_error <<- conditionMessage(e)
    NULL
  })

  if (!is.null(extracted_prediction)) {
    standardized <- asreml_standardize_extracted_predictions(
      extracted_df = extracted_prediction,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
    return(list(
      mode = "extract",
      prediction = standardized$single_prediction,
      across_env_prediction = standardized$across_env_prediction,
      error = direct_prediction$error
    ))
  }

  list(
    mode = "legacy_fallback",
    prediction = NULL,
    across_env_prediction = NULL,
    error = paste(
      c(
        direct_prediction$error,
        extracted_error
      ),
      collapse = " | "
    )
  )
}
