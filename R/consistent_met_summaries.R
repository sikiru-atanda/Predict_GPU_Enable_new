## Phase 3.16: post-processor that harmonises the MET output structure
## across all model families.
##
## For every model that produced per-env predictions, attach:
##   - environment_correlation : k x k empirical correlation matrix of
##     per-env predictions (works for any model -- this is the key
##     consistency win the user asked for).
##   - environment_covariance  : populated when the engine natively
##     reports a (genetic) covariance matrix across envs (GP / Bayesian /
##     ASReml); otherwise NA. ML/DL models cannot give a covariance
##     matrix on the original-trait scale, so we set NA.
##
## Attached at two locations so users get a consistent path no matter
## which engine they ran:
##   - sik$model_results$environment_correlation  (true prediction)
##   - sik$model_results$environment_covariance
##   - sik$cv_results_processed$environment_correlation  (CV)
##   - sik$cv_results_processed$environment_covariance
## In both cases the value is a NAMED LIST keyed by model name.

#' Reshape long predictions to a wide (gid x env) numeric matrix.
#'
#' Internal helper. Used by both env correlation and env covariance.
#'
#' @keywords internal
gp_long_predictions_to_wide_matrix <- function(preds, gid_col, env_col, value_col) {
  if (is.null(preds) || !is.data.frame(preds) || nrow(preds) == 0L) {
    return(NULL)
  }
  if (!all(c(gid_col, env_col, value_col) %in% names(preds))) {
    return(NULL)
  }
  df <- preds[, c(gid_col, env_col, value_col), drop = FALSE]
  names(df) <- c("gid", "env", "val")
  df$gid <- as.character(df$gid)
  df$env <- as.character(df$env)
  df$val <- suppressWarnings(as.numeric(df$val))
  df <- df[is.finite(df$val), , drop = FALSE]
  if (nrow(df) < 2L) return(NULL)
  envs <- sort(unique(df$env))
  if (length(envs) < 2L) return(NULL)
  agg <- stats::aggregate(val ~ gid + env, data = df, FUN = mean,
                          na.rm = TRUE, na.action = stats::na.pass)
  wide <- stats::reshape(agg, idvar = "gid", timevar = "env",
                         direction = "wide")
  mat <- as.matrix(wide[, setdiff(names(wide), "gid"), drop = FALSE])
  colnames(mat) <- sub("^val\\.", "", colnames(mat))
  mat
}

#' Compute empirical environment correlation matrix from per-env predictions.
#'
#' @keywords internal
gp_empirical_env_correlation <- function(preds, gid_col, env_col, value_col) {
  mat <- gp_long_predictions_to_wide_matrix(preds, gid_col, env_col, value_col)
  if (is.null(mat)) return(NULL)
  cor_mat <- suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs"))
  if (!is.matrix(cor_mat)) return(NULL)
  cor_mat
}

#' Compute empirical environment covariance matrix from per-env predictions.
#'
#' Sister of gp_empirical_env_correlation. Used for engines that natively
#' report a genetic-covariance-across-envs (GP / ASReml GBLUP corgh /
#' Bayesian kernel RKHS / GBLUP_BRR / BRR) but whose native covariance is
#' not surfaced at the top level of the assembled result. The empirical
#' covariance of per-env predictions captures the same env-by-env signal
#' the engine measured, on the original prediction scale.
#'
#' @keywords internal
gp_empirical_env_covariance <- function(preds, gid_col, env_col, value_col) {
  mat <- gp_long_predictions_to_wide_matrix(preds, gid_col, env_col, value_col)
  if (is.null(mat)) return(NULL)
  cov_mat <- suppressWarnings(stats::cov(mat, use = "pairwise.complete.obs"))
  if (!is.matrix(cov_mat)) return(NULL)
  cov_mat
}

#' Engines whose per-model output natively reports an environment
#' covariance matrix.
#'
#' GP family always exposes interaction_cov natively.
#' GBLUP (asreml) always exposes a genetic-cov-across-envs structure.
#' Kernel-Bayes (RKHS / GBLUP_BRR / BRR) only exposes a per-env covariance
#' when run through the BGLR::Multitrait env-as-trait path -- in the
#' BGLR::BGLR univariate (homogeneous, single sigma_g) path there is no
#' env-by-env genetic covariance to report. Callers that need the
#' homogeneous-mode exclusion must pass `bayes_kernel_heter` = FALSE.
#'
#' @keywords internal
gp_models_with_native_env_covariance <- function(bayes_kernel_heter = TRUE) {
  base <- c(
    # GP family: gp_exact / gp_icm_fa / krr_exact populate interaction_cov
    "Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP",
    "LowRankGP", "KRR", "GP", "GP_FA",
    # ASReml-engine GBLUP under engine="asreml" with var_cov_str =
    # "corgh" / "us" / "fa1" emits genetic_covariance via summary(fit)
    "GBLUP"
  )
  if (isTRUE(bayes_kernel_heter)) {
    # Bayesian kernel-prior MET fits via BGLR::Multitrait emit per-env
    # covariance through bayes_multitrait_env_heter_fit.
    base <- c(base, "RKHS", "GBLUP_BRR", "BRR")
  }
  base
}

#' Post-process the assembled `model_execute` result to add a consistent
#' environment_correlation / environment_covariance pair for every model.
#'
#' Idempotent: if the fields are already present we leave them.
#' Defensive: any error inside is caught and the function returns the
#' input unchanged. Never blocks the main result return.
#'
#' @param sik the final result returned by `gp_finalize_model_execute_results`.
#' @param ctx the model_execute environment (or list) with `gen_name`,
#'   `heter_groups`, `cross_validation`, etc.
#' @keywords internal
gp_add_consistent_met_summaries <- function(sik, ctx = NULL) {
  if (!is.list(sik)) return(sik)
  heter_groups <- if (is.list(ctx)) ctx[["heter_groups"]] else NULL
  gen_name     <- if (is.list(ctx)) ctx[["gen_name"]]     else NULL
  pheno_data   <- if (is.list(ctx)) ctx[["pheno_data"]]   else NULL

  # Phase 3.23 contract:
  # 1. environment_correlation / environment_covariance are TRUE-PREDICTION
  #    properties of the fitted model, NOT a CV-output. Strip them from the
  #    CV side unconditionally.
  # 2. On the true-prediction side, populate ONLY for engines that natively
  #    expose an env-by-env genetic covariance (GP family, GBLUP via
  #    asreml, and kernel-Bayes when run in heterogeneous BGLR::Multitrait
  #    mode). ML / DL / Lasso / Ridge / PLS / kernel-Bayes-homogeneous have
  #    no model-level env covariance and are intentionally omitted.
  if (is.list(sik[["cv_results_processed"]])) {
    sik[["cv_results_processed"]][["environment_correlation"]] <- NULL
    sik[["cv_results_processed"]][["environment_covariance"]] <- NULL
  }

  if (is.null(heter_groups) || !nzchar(as.character(heter_groups)[1L])) {
    # Single-env run -- per-env covariance isn't meaningful.
    return(sik)
  }
  env_col <- as.character(heter_groups)[1L]
  gid_col <- as.character(gen_name %||% "GID")[1L]

  # Resolve whether kernel-Bayes was run in the heterogeneous Multitrait
  # path. `bayes_kernel_heter_resid` (Phase 3.21) overrides `heter_resid`
  # for the kernel-Bayes branch only; default is to follow heter_resid.
  bayes_kernel_heter <- if (is.list(ctx)) {
    bkhr <- ctx[["bayes_kernel_heter_resid"]]
    if (is.null(bkhr)) isTRUE(ctx[["heter_resid"]]) else isTRUE(bkhr)
  } else TRUE
  models_with_cov <- gp_models_with_native_env_covariance(
    bayes_kernel_heter = bayes_kernel_heter
  )

  per_model_corr <- list()
  per_model_cov  <- list()

  # ---- CV side: build empirical correlation per model from cv_results_raw
  # which has one entry per (model, trait, rep) with a $ypred_cv_Reps_all
  # data.frame (columns row_id, y, yhat, cv_role, Loc). We aggregate yhat
  # across all (trait, rep) for the same model into a single (GID, env)
  # long frame, then call the empirical correlation/covariance helpers.
  #
  # row_id is the index into the original pheno_data, so we use it to
  # look up GID. If pheno_data isn't available we fall back to using
  # row_id itself as a pseudo-GID (still gives per-env signal if every
  # (row_id, env) pair is unique).
  cv_raw <- sik[["cv_results_raw"]] %||% sik[["cv_results"]]
  if (is.list(cv_raw) && length(cv_raw) > 0L) {
    long_by_model <- list()
    for (entry in cv_raw) {
      if (!is.list(entry)) next
      # Prefer the user-supplied friendly name so the keys match what's in
      # the accuracy table (e.g. "FA-GBLUP" not "GP_FA"). Fall back to
      # canonical if friendly is absent.
      model_label <- entry[["model"]] %||% entry[["model_canonical"]]
      if (is.null(model_label) || !nzchar(model_label)) next
      model_label <- as.character(model_label)[1L]
      yp <- entry[["ypred_cv_Reps_all"]]
      if (!is.data.frame(yp) || nrow(yp) == 0L) next
      env_col_in <- if (env_col %in% names(yp)) env_col
                    else if ("Env" %in% names(yp)) "Env"
                    else if ("Loc" %in% names(yp)) "Loc"
                    else NULL
      if (is.null(env_col_in) || !("yhat" %in% names(yp))) next
      keep <- if ("cv_role" %in% names(yp)) yp[["cv_role"]] == "test" else rep(TRUE, nrow(yp))
      sub <- yp[keep, , drop = FALSE]
      if (nrow(sub) == 0L) next
      # Resolve GID for each row_id from pheno_data, else fall back to row_id.
      gid_vec <- if ("row_id" %in% names(sub) && !is.null(pheno_data) &&
                     is.data.frame(pheno_data) && gid_col %in% names(pheno_data)) {
        rid <- as.integer(sub[["row_id"]])
        rid[rid < 1L | rid > nrow(pheno_data)] <- NA_integer_
        as.character(pheno_data[rid, gid_col])
      } else if (gid_col %in% names(sub)) {
        as.character(sub[[gid_col]])
      } else {
        as.character(sub[["row_id"]] %||% seq_len(nrow(sub)))
      }
      df_long <- data.frame(
        gid = gid_vec,
        env = as.character(sub[[env_col_in]]),
        yhat = suppressWarnings(as.numeric(sub[["yhat"]])),
        stringsAsFactors = FALSE
      )
      df_long <- df_long[is.finite(df_long$yhat) & !is.na(df_long$gid), , drop = FALSE]
      if (nrow(df_long) == 0L) next
      long_by_model[[model_label]] <- rbind(long_by_model[[model_label]], df_long)
    }
    for (m in names(long_by_model)) {
      # Phase 3.23: skip models with no native env covariance entirely --
      # we don't emit ML/DL/Lasso/Ridge/PLS or kernel-Bayes-homogeneous.
      if (!(m %in% models_with_cov)) next
      sub <- long_by_model[[m]]
      names(sub) <- c(gid_col, env_col, "Predicted_value")
      per_model_corr[[m]] <- gp_empirical_env_correlation(
        sub, gid_col = gid_col, env_col = env_col, value_col = "Predicted_value")
      per_model_cov[[m]] <- gp_empirical_env_covariance(
        sub, gid_col = gid_col, env_col = env_col, value_col = "Predicted_value")
    }
  }

  # ---- True-prediction side: build empirical correlation/covariance from
  # sik$model_results$predictions OR per-model_results_by_model entries.
  mr <- sik[["model_results"]]
  if (is.list(mr)) {
    preds_tp <- tryCatch(mr[["predictions"]], error = function(e) NULL)
    if (is.data.frame(preds_tp) && nrow(preds_tp) > 0L) {
      val_col <- if ("Predicted_value" %in% names(preds_tp)) "Predicted_value"
                 else if ("predicted_value" %in% names(preds_tp)) "predicted_value"
                 else NULL
      mod_col <- if ("model" %in% names(preds_tp)) "model"
                 else if ("Model" %in% names(preds_tp)) "Model"
                 else NULL
      if (!is.null(val_col) && !is.null(mod_col) &&
          gid_col %in% names(preds_tp) && env_col %in% names(preds_tp)) {
        for (m in unique(as.character(preds_tp[[mod_col]]))) {
          # Same gate: only models with native env covariance.
          if (!(m %in% models_with_cov)) next
          if (!(m %in% names(per_model_corr))) {
            sub <- preds_tp[as.character(preds_tp[[mod_col]]) == m, , drop = FALSE]
            per_model_corr[[m]] <- gp_empirical_env_correlation(
              sub, gid_col = gid_col, env_col = env_col, value_col = val_col)
            per_model_cov[[m]] <- gp_empirical_env_covariance(
              sub, gid_col = gid_col, env_col = env_col, value_col = val_col)
          }
        }
      }
    }
  }

  # ---- Attach to TRUE PREDICTION side only (not CV).
  # Don't overwrite if engine-native value already exists.
  if (length(per_model_corr) > 0L || length(per_model_cov) > 0L) {
    mr <- sik[["model_results"]]
    if (!is.list(mr)) mr <- list()
    if (is.null(mr[["environment_correlation"]])) {
      mr[["environment_correlation"]] <- per_model_corr
    }
    if (is.null(mr[["environment_covariance"]])) {
      mr[["environment_covariance"]] <- per_model_cov
    }
    sik[["model_results"]] <- mr
  }
  sik
}
