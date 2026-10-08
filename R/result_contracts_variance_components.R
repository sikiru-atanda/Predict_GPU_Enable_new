gp_variance_component_columns <- function(include_trait = FALSE, include_env = FALSE,
                                          include_boundary = FALSE) {
  c(
    if (isTRUE(include_trait)) "Trait",
    if (isTRUE(include_env)) "Env",
    "Component",
    "Components",
    "Standard_error",
    # Optional boundary flag, appended only when present (see
    # gp_format_gp_multi_env_variance_components). Mirrors the on-demand
    # Trait/Env columns so boundary-free tables keep the canonical 3 columns.
    if (isTRUE(include_boundary)) "Boundary"
  )
}

gp_covariance_component_columns <- function() {
  c("Matrix", "Component", "Row", "Column", "Value", "Estimation_method")
}

gp_empty_variance_components <- function() {
  data.frame(
    Component = character(),
    Components = numeric(),
    Standard_error = numeric(),
    stringsAsFactors = FALSE
  )
}

gp_vc_finite_var <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) > 1L) {
    stats::var(x, na.rm = TRUE)
  } else {
    NA_real_
  }
}

gp_vc_finite_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) {
    mean(x, na.rm = TRUE)
  } else {
    NA_real_
  }
}

gp_vc_finite_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) > 1L) {
    stats::sd(x, na.rm = TRUE)
  } else {
    NA_real_
  }
}

gp_vc_simple_table <- function(genetic_variance,
                               residual_variance,
                               heritability,
                               genetic_se = NA_real_,
                               residual_se = NA_real_,
                               heritability_se = NA_real_,
                               estimation_method = "model_reported") {
  component <- c("genetic_variance", "residual_variance", "heritability")
  out <- data.frame(
    Component = component,
    Components = c(genetic_variance, residual_variance, heritability),
    Standard_error = c(genetic_se, residual_se, heritability_se),
    stringsAsFactors = FALSE
  )
  rownames(out) <- component
  out
}

#' Append per-environment genetic variance + heritability to a MET varcomp table
#'
#' Internal helper for cross-engine MET output consistency. Given an engine's
#' variance-components table and per-environment genetic and residual variances
#' (named numeric vectors aligned by environment), this appends
#' `genetic_variance_<env>`, `residual_variance_<env>` and `heritability_<env>`
#' rows (h2 = Vg / (Vg + Ve)) so GP, Bayesian and ASReml MET outputs all carry
#' per-env Vg + h2. Existing rows and the table's column layout are preserved.
#'
#' @param vc A variance-components data frame (typically with `Components` and
#'   `Standard_error` columns, optionally a `Component` label column).
#' @param env_genetic,env_residual Named numeric vectors (names = environment
#'   labels, same order) of per-env genetic and residual variances.
#' @return `vc` with per-env rows appended.
#' @keywords internal
#' @noRd
gp_variance_components_add_per_env <- function(vc, env_genetic, env_residual,
                                                env_genetic_se = NULL,
                                                env_residual_se = NULL,
                                                env_h2_se = NULL) {
  if (!is.data.frame(vc) || is.null(env_genetic) || !length(env_genetic)) return(vc)
  envs <- names(env_genetic)
  if (is.null(envs) || !all(envs == names(env_residual))) {
    # fall back to positional alignment if names are absent/mismatched
    envs <- if (!is.null(names(env_residual))) names(env_residual) else paste0("Env", seq_along(env_genetic))
  }
  eg <- as.double(env_genetic); er <- as.double(env_residual)
  h2 <- eg / (eg + er)

  # Per-env SEs. When the caller supplies them we use them directly; otherwise
  # leave the corresponding row's Standard_error at NA. h^2 SE: if a precomputed
  # per-env h2 SE vector is supplied (e.g. from the FA delta-method) we use it;
  # otherwise compute the 2-term delta-method SE from genetic + residual SEs
  # assuming Cov(sigma2_g, sigma2_e) = 0 (standard approximation, see
  # gp_format_gp_single_env_variance_components for the analogous single-env
  # formula).
  align <- function(vec, names_target) {
    if (is.null(vec)) return(rep(NA_real_, length(names_target)))
    out <- suppressWarnings(as.numeric(vec))
    if (!is.null(names(vec))) {
      out <- out[match(names_target, names(vec))]
    } else if (length(out) != length(names_target)) {
      out <- rep(NA_real_, length(names_target))
    }
    out
  }
  eg_se <- align(env_genetic_se,  envs)
  er_se <- align(env_residual_se, envs)
  h2_se_in <- align(env_h2_se,    envs)

  total <- eg + er
  h2_se_delta <- rep(NA_real_, length(envs))
  ok <- is.finite(eg) & is.finite(er) & is.finite(eg_se) & is.finite(er_se) & total > 0
  if (any(ok)) {
    d_g <- er[ok] / total[ok]^2
    d_r <- -eg[ok] / total[ok]^2
    h2_se_delta[ok] <- sqrt((d_g^2 * eg_se[ok]^2) + (d_r^2 * er_se[ok]^2))
  }
  h2_se_out <- ifelse(is.finite(h2_se_in), h2_se_in, h2_se_delta)

  rn <- c(paste0("genetic_variance_", envs),
          paste0("residual_variance_", envs),
          paste0("heritability_", envs))
  vals    <- c(eg,    er,    h2)
  se_vals <- c(eg_se, er_se, h2_se_out)
  cols <- names(vc)
  add <- as.data.frame(matrix(NA, nrow = length(rn), ncol = length(cols)), stringsAsFactors = FALSE)
  names(add) <- cols
  if ("Components" %in% cols)     add[["Components"]]     <- vals
  if ("Standard_error" %in% cols) add[["Standard_error"]] <- se_vals
  if ("Component" %in% cols)      add[["Component"]]      <- rn
  rownames(add) <- rn
  rbind(vc, add)
}

#' Reconstruct per-environment Vg + h2 from a GP MET varcomp table
#'
#' Internal helper covering the two MET variance-structure layouts the GP backend
#' produces, both with per-env residual rows `<env>!R`:
#'
#' * Factor-analytic (GP_FA): per-env loadings `...!<env>!fa<k>` and specific
#'   variances `...!<env>!var`. Per-env genetic variance is the diagonal of
#'   `Lambda Lambda' + Psi`: `Vg_env = sum_k loading_k_env^2 + specific_env`.
#' * Compound-symmetry / homogeneous (GP, KRR, LowRankGP): a single main genetic
#'   variance `vm(GID,...)!var` (+ optional GxE `Env!var`), homogeneous across
#'   environments: `Vg_env = sum of the non-residual `!var` rows` (boundary NA
#'   components contribute 0).
#'
#' In both cases `h2_env = Vg_env / (Vg_env + Ve_env)`, appended via
#' [gp_variance_components_add_per_env]. Tables without per-env residual rows are
#' returned unchanged.
#'
#' @param vc A GP MET variance-components data frame.
#' @return `vc` with per-env genetic/residual/heritability rows appended, else `vc`.
#' @keywords internal
#' @noRd
gp_met_varcomp_augment_per_env <- function(vc, fa_derived = NULL) {
  if (!is.data.frame(vc) || !nrow(vc)) return(vc)
  rn <- rownames(vc)
  if (is.null(rn)) return(vc)
  resid_idx <- grep("!R$", rn)
  if (!length(resid_idx)) return(vc)                    # no per-env residual structure
  valcol <- if ("Components" %in% names(vc)) "Components" else
    names(vc)[vapply(vc, is.numeric, logical(1))][1]
  if (is.na(valcol)) return(vc)
  val <- suppressWarnings(as.numeric(vc[[valcol]]))
  se  <- if ("Standard_error" %in% names(vc)) suppressWarnings(as.numeric(vc[["Standard_error"]])) else rep(NA_real_, length(val))
  esc <- function(e) gsub("([.^$*+?()\\[\\]{}|\\\\])", "\\\\\\1", e, perl = TRUE)
  resid_env <- sub("^Env_", "", sub("!R$", "", rn[resid_idx]))
  env_residual    <- stats::setNames(val[resid_idx], resid_env)
  env_residual_se <- stats::setNames(se[resid_idx],  resid_env)

  env_h2_se <- NULL  # set only by the FA branch (delta-method via fa_derived)

  fa_idx <- grep("!fa[0-9]+$", rn)
  if (length(fa_idx)) {
    # Factor-analytic: per-env Vg = sum_k loading_k^2 + specific_var.
    # Estimates come from varcomp (well-defined). SEs are NOT a simple
    # function of the loading/specific SEs -- they need the delta-method
    # propagation through the AI inverse, which the Python backend already
    # computes and surfaces via `fa_derived` (see
    # varcomp_asreml.compute_fa_derived_summary, plumbed in 0.20.13). When
    # fa_derived is available we use its per-env genetic SE and h^2 SE.
    envs <- unique(sub("^.*!([^!]+)!(fa[0-9]+|var)$", "\\1",
                       rn[grep("!([^!]+)!(fa[0-9]+|var)$", rn)]))
    envs <- envs[nzchar(envs)]
    if (!length(envs)) return(vc)
    env_genetic    <- stats::setNames(rep(NA_real_, length(envs)), envs)
    env_genetic_se <- stats::setNames(rep(NA_real_, length(envs)), envs)
    for (e in envs) {
      loads <- val[grep(paste0("!", esc(e), "!fa[0-9]+$"), rn)]
      spec  <- val[grep(paste0("!", esc(e), "!var$"), rn)]
      spec  <- if (length(spec) && is.finite(spec[1])) spec[1] else 0
      env_genetic[e] <- sum(loads^2, na.rm = TRUE) + spec
    }
    env_residual    <- env_residual[names(env_genetic)]
    env_residual_se <- env_residual_se[names(env_genetic)]

    if (is.list(fa_derived) && length(fa_derived[["env_labels"]])) {
      fa_labels <- as.character(fa_derived[["env_labels"]])
      fa_g      <- as.numeric(fa_derived[["sigma2_g_per_env"]]    %||% rep(NA_real_, length(fa_labels)))
      fa_g_se   <- as.numeric(fa_derived[["sigma2_g_se_per_env"]] %||% rep(NA_real_, length(fa_labels)))
      fa_h2_se  <- as.numeric(fa_derived[["h2_se_per_env"]]       %||% rep(NA_real_, length(fa_labels)))
      idx <- match(envs, fa_labels)
      matched <- !is.na(idx)
      derived_g <- rep(NA_real_, length(envs))
      derived_g[matched] <- fa_g[idx[matched]]
      use_derived_g <- is.finite(derived_g)
      env_genetic[use_derived_g] <- derived_g[use_derived_g]
      env_genetic_se[!is.na(idx)] <- fa_g_se[idx[!is.na(idx)]]
      env_h2_se <- stats::setNames(rep(NA_real_, length(envs)), envs)
      env_h2_se[!is.na(idx)] <- fa_h2_se[idx[!is.na(idx)]]
    }
  } else {
    # Compound-symmetry / homogeneous: genetic = sum of non-residual !var
    # rows, the same in every environment. The SE on the shared sigma2_g
    # row carries through to every env (since per-env sigma2_g IS sigma2_g
    # in this parameterisation; no derived quantity, no delta-method).
    gen_idx <- setdiff(grep("!var$", rn), resid_idx)
    if (!length(gen_idx)) return(vc)
    # Phase 3.28: aggregate only the kernels with a finite, well-estimated
    # weight (finite value AND finite SE). A boundary-clipped kernel
    # (value finite, SE = NA) contributes ~0 to total Vg by construction;
    # adding it to the sum makes essentially no numerical difference but
    # propagating its missing SE would NA out vg_se, hiding the well-
    # estimated kernel's SE from users.
    g_val <- val[gen_idx]; g_se <- se[gen_idx]
    keep_g <- is.finite(g_val) & is.finite(g_se)
    if (!any(keep_g)) {
      # All kernels are boundary or NA; preserve the legacy behavior:
      # sum what we can (na.rm), report NA SE.
      vg     <- sum(g_val, na.rm = TRUE)
      vg_se  <- NA_real_
    } else {
      vg     <- sum(g_val[keep_g])
      vg_se  <- sqrt(sum(g_se[keep_g] ^ 2))
    }
    env_genetic    <- stats::setNames(rep(vg,    length(resid_env)), resid_env)
    env_genetic_se <- stats::setNames(rep(vg_se, length(resid_env)), resid_env)
  }
  # An environment whose residual is not estimable (estimated at zero) has no
  # separable genetic variance either: its genetic variance and h2 are NA too.
  boundary_flag <- if ("Boundary" %in% names(vc)) as.character(vc[["Boundary"]]) else rep("", nrow(vc))
  not_estimable_env <- resid_env[boundary_flag[resid_idx] %in% "not_estimable"]
  if (length(not_estimable_env)) {
    drop <- names(env_genetic) %in% not_estimable_env
    env_genetic[drop] <- NA_real_
    env_genetic_se[names(env_genetic_se) %in% not_estimable_env] <- NA_real_
    if (!is.null(env_h2_se)) env_h2_se[names(env_h2_se) %in% not_estimable_env] <- NA_real_
  }
  if ((all(is.na(env_residual)) || all(is.na(env_genetic))) && !length(not_estimable_env)) return(vc)
  out <- gp_variance_components_add_per_env(
    vc,
    env_genetic     = env_genetic,
    env_residual    = env_residual,
    env_genetic_se  = env_genetic_se,
    env_residual_se = env_residual_se,
    env_h2_se       = env_h2_se
  )
  if (length(not_estimable_env)) {
    if (!"Boundary" %in% names(out)) out[["Boundary"]] <- ""
    out[["Boundary"]][is.na(out[["Boundary"]])] <- ""
    flagged <- rownames(out) %in% c(paste0("genetic_variance_", not_estimable_env),
                                    paste0("residual_variance_", not_estimable_env),
                                    paste0("heritability_", not_estimable_env))
    out[["Boundary"]][flagged] <- "not_estimable"
  }
  out
}

#' Relabel internal "Env" tokens in a GP variance-component table
#'
#' The GP backend standardises the environment factor to the internal name
#' "Env", so variance-component rows come back labelled `Env!var` (GxE / env-main
#' variance) and `Env_<level>!R` (per-env residual) regardless of the column the
#' user supplied. This rewrites those labels to the user's `heter_groups` name
#' (e.g. "YYY") in both the `Component` column and the rownames, so every model's
#' output uses the user-defined variable name. Per-env summary rows
#' (`genetic_variance_<level>` etc.) carry the env LEVEL, not the token "Env", so
#' they are untouched. No-op when `heter_groups` is NULL/empty or already "Env".
#'
#' MUST run AFTER [gp_met_varcomp_augment_per_env], which keys per-env residual
#' detection off the `Env_` prefix.
#'
#' @keywords internal
#' @noRd
gp_vc_relabel_env_to_heter_groups <- function(vc, heter_groups) {
  if (!is.data.frame(vc) || is.null(heter_groups) ||
      !nzchar(heter_groups[1L]) || identical(heter_groups[1L], "Env")) {
    return(vc)
  }
  hg <- heter_groups[1L]
  relabel <- function(x) {
    x <- sub("^Env!", paste0(hg, "!"), x)            # Env!var      -> <hg>!var
    x <- sub("^Env_", paste0(hg, "_"), x)            # Env_<lvl>!R  -> <hg>_<lvl>!R
    x <- sub("^Env:", paste0(hg, ":"), x)            # Env:vm(...)  -> <hg>:vm(...)
    # Variance-structure wrappers around the env factor, e.g.
    # fa(Env,1):vm(GID,G)!.., idv(Env):.., us(Env):.., corgh(Env):..
    gsub("\\(Env([,):])", paste0("(", hg, "\\1"), x)
  }
  if ("Component" %in% names(vc)) {
    vc[["Component"]] <- relabel(as.character(vc[["Component"]]))
  }
  rn <- rownames(vc)
  if (!is.null(rn)) {
    rownames(vc) <- relabel(rn)
  }
  vc
}

gp_unavailable_variance_components <- function(reason = "variance_components_not_estimated",
                                               estimation_method = "not_estimated") {
  out <- data.frame(
    Component = reason,
    Components = NA_real_,
    Standard_error = NA_real_,
    stringsAsFactors = FALSE
  )
  rownames(out) <- reason
  out
}

gp_vc_has_nonempty_labels <- function(x) {
  any(!is.na(x) & nzchar(as.character(x)))
}

gp_format_variance_component_columns <- function(out,
                                                 include_trait = NULL,
                                                 include_env = NULL,
                                                 include_boundary = NULL) {
  if (!is.data.frame(out)) {
    return(gp_empty_variance_components())
  }
  include_trait <- if (is.null(include_trait)) {
    "Trait" %in% names(out) && gp_vc_has_nonempty_labels(out[["Trait"]])
  } else {
    isTRUE(include_trait)
  }
  include_env <- if (is.null(include_env)) {
    "Env" %in% names(out) && gp_vc_has_nonempty_labels(out[["Env"]])
  } else {
    isTRUE(include_env)
  }
  # Boundary is a pass-through optional column: keep it only when the input
  # already carries non-empty boundary flags; never fabricate it.
  include_boundary <- if (is.null(include_boundary)) {
    "Boundary" %in% names(out) && gp_vc_has_nonempty_labels(out[["Boundary"]])
  } else {
    isTRUE(include_boundary)
  }
  if (isTRUE(include_trait) && !"Trait" %in% names(out)) {
    out[["Trait"]] <- NA_character_
  }
  if (isTRUE(include_env) && !"Env" %in% names(out)) {
    out[["Env"]] <- NA_character_
  }
  required <- c("Component", "Components", "Standard_error")
  for (nm in required) {
    if (!nm %in% names(out)) {
      out[[nm]] <- if (nm %in% c("Components", "Standard_error")) NA_real_ else NA_character_
    }
  }
  out <- out[, gp_variance_component_columns(include_trait = include_trait,
                                             include_env = include_env,
                                             include_boundary = include_boundary), drop = FALSE]
  rownames(out) <- NULL
  out
}

gp_ml_gaussian_variance_components <- function(boot_matrix = NULL,
                                               train_test_label = NULL,
                                               observed_y = NULL,
                                               predicted_value = NULL,
                                               reference_variance = NULL,
                                               prediction_error_var = NULL,
                                               marginal_prediction_mse = NULL) {
  boot <- NULL
  if (!is.null(boot_matrix)) {
    boot <- tryCatch(as.matrix(boot_matrix), error = function(e) NULL)
    if (!is.null(boot)) {
      storage.mode(boot) <- "double"
      if (is.null(dim(boot))) {
        boot <- matrix(boot, nrow = 1L)
      }
    }
  }

  n_pred <- length(predicted_value %||% numeric())
  if (!n_pred && !is.null(boot)) {
    n_pred <- ncol(boot)
  }
  if (!n_pred && !is.null(train_test_label)) {
    n_pred <- length(train_test_label)
  }
  if (!n_pred) {
    return(gp_empty_variance_components())
  }

  train_idx <- seq_len(n_pred)
  if (!is.null(train_test_label) && length(train_test_label) == n_pred) {
    labels <- as.character(train_test_label)
    train_idx <- which(is.na(labels) | labels != "Test")
    if (!length(train_idx)) {
      train_idx <- seq_len(n_pred)
    }
  }

  observed_full <- rep(NA_real_, n_pred)
  if (!is.null(observed_y)) {
    observed_y <- suppressWarnings(as.numeric(observed_y))
    if (length(observed_y) == n_pred) {
      observed_full <- observed_y
    } else if (length(observed_y) == length(train_idx)) {
      observed_full[train_idx] <- observed_y
    }
  }

  boot_prediction_variance <- numeric()
  if (!is.null(boot) && ncol(boot) >= max(train_idx)) {
    boot_prediction_variance <- apply(
      boot[, train_idx, drop = FALSE],
      1L,
      gp_vc_finite_var
    )
  }

  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  phenotype_variance <- suppressWarnings(as.numeric(reference_variance))[1L]
  if (!is.finite(phenotype_variance) || phenotype_variance < 0) {
    phenotype_variance <- gp_vc_finite_var(observed_full[train_idx])
  }

  prediction_variance <- gp_vc_finite_var(predicted_value[train_idx])
  if (!is.finite(prediction_variance) || prediction_variance < 0) {
    prediction_variance <- gp_vc_finite_mean(boot_prediction_variance)
  }

  prediction_error <- gp_vc_finite_mean(prediction_error_var[train_idx])
  if (!is.finite(prediction_error)) {
    marginal_prediction_mse <- suppressWarnings(
      as.numeric(marginal_prediction_mse)[1L]
    )
    if (is.finite(marginal_prediction_mse) && marginal_prediction_mse >= 0) {
      prediction_error <- marginal_prediction_mse
    }
  }

  component <- c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable",
    "phenotypic_reference_variance",
    "prediction_variance",
    "estimated_prediction_mse"
  )
  out <- data.frame(
    Component = component,
    Components = c(
      NA_real_, NA_real_, NA_real_, phenotype_variance,
      prediction_variance, prediction_error
    ),
    Standard_error = c(
      NA_real_, NA_real_, NA_real_, NA_real_,
      gp_vc_finite_sd(boot_prediction_variance),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
  rownames(out) <- component
  gp_format_variance_component_columns(out)
}

gp_prediction_table_variance_components <- function(pred) {
  if (!is.data.frame(pred) ||
      !"Predicted_value" %in% names(pred) ||
      !"Train_Test_Label" %in% names(pred)) {
    return(gp_empty_variance_components())
  }

  pred <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)
  group_cols <- character()
  if ("Trait" %in% names(pred) &&
      any(!is.na(pred[["Trait"]]) & nzchar(as.character(pred[["Trait"]]))) &&
      gp_contract_multiple_nonempty_values(pred[["Trait"]])) {
    group_cols <- c(group_cols, "Trait")
  }
  if ("Env" %in% names(pred) &&
      any(!is.na(pred[["Env"]]) & nzchar(as.character(pred[["Env"]]))) &&
      gp_contract_multiple_nonempty_values(pred[["Env"]])) {
    group_cols <- c(group_cols, "Env")
  }
  pred_groups <- if (length(group_cols)) {
    split(pred, interaction(pred[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE), drop = TRUE)
  } else {
    list(pred)
  }
  pieces <- lapply(pred_groups, function(dat) {
    train <- as.character(dat[["Train_Test_Label"]]) != "Test"
    if (!any(train, na.rm = TRUE)) {
      train <- rep(TRUE, nrow(dat))
    }
    pred_train <- suppressWarnings(as.numeric(dat[["Predicted_value"]][train]))
    obs_train <- if ("Observed_value" %in% names(dat)) {
      suppressWarnings(as.numeric(dat[["Observed_value"]][train]))
    } else {
      rep(NA_real_, length(pred_train))
    }
    prediction_variance <- gp_vc_finite_var(pred_train)
    phenotype_variance <- gp_vc_finite_var(obs_train)
    prediction_error <- if ("PEV" %in% names(dat)) {
      gp_vc_finite_mean(dat[["PEV"]][train])
    } else {
      NA_real_
    }
    if ((!is.finite(prediction_error) || prediction_error < 0) &&
        "Prediction_error_variance" %in% names(dat)) {
      prediction_error <- gp_vc_finite_mean(
        dat[["Prediction_error_variance"]][train]
      )
    }
    component <- c(
      "genetic_variance_not_identifiable",
      "residual_variance_not_identifiable",
      "heritability_not_identifiable",
      "phenotypic_reference_variance",
      "prediction_variance",
      "estimated_prediction_mse"
    )
    out <- data.frame(
      Component = component,
      Components = c(
        NA_real_, NA_real_, NA_real_, phenotype_variance,
        prediction_variance, prediction_error
      ),
      Standard_error = NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(group_cols)) {
      labels <- stats::setNames(
        lapply(group_cols, function(col) as.character(dat[[col]][[1L]])),
        group_cols
      )
      for (nm in names(labels)) {
        out[[nm]] <- labels[[nm]]
      }
      rownames(out) <- NULL
    }
    out <- gp_format_variance_component_columns(
      out,
      include_trait = "Trait" %in% group_cols,
      include_env = "Env" %in% group_cols
    )
    out
  })
  pieces <- pieces[vapply(pieces, is.data.frame, logical(1L))]
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (!length(pieces)) {
    return(gp_empty_variance_components())
  }
  out <- do.call(rbind, pieces)
  out <- gp_format_variance_component_columns(
    out,
    include_trait = "Trait" %in% group_cols,
    include_env = "Env" %in% group_cols
  )
  if (!length(group_cols)) {
    rownames(out) <- make.unique(out[["Component"]])
  }
  out
}

gp_ml_classification_variance_components <- function(pred) {
  if (!is.data.frame(pred)) {
    return(gp_empty_variance_components())
  }
  pred <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(pred)) {
    return(gp_empty_variance_components())
  }

  group_cols <- character()
  for (nm in c("Trait", "Env")) {
    if (nm %in% names(pred) &&
        any(!is.na(pred[[nm]]) & nzchar(as.character(pred[[nm]]))) &&
        gp_contract_multiple_nonempty_values(pred[[nm]])) {
      group_cols <- c(group_cols, nm)
    }
  }
  pred_groups <- if (length(group_cols)) {
    split(
      pred,
      interaction(pred[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE),
      drop = TRUE
    )
  } else {
    list(pred)
  }

  pieces <- lapply(pred_groups, function(dat) {
    train <- rep(TRUE, nrow(dat))
    if ("Train_Test_Label" %in% names(dat)) {
      labels <- tolower(trimws(as.character(dat[["Train_Test_Label"]])))
      train <- is.na(labels) | labels != "test"
      if (!any(train, na.rm = TRUE)) {
        train <- rep(TRUE, nrow(dat))
      }
    }

    confidence <- rep(NA_real_, nrow(dat))
    if ("Prediction_confidence" %in% names(dat)) {
      confidence <- suppressWarnings(as.numeric(dat[["Prediction_confidence"]]))
    }
    probability_cols <- grep("^Probability_", names(dat), value = TRUE)
    if (length(probability_cols)) {
      probability_matrix <- suppressWarnings(as.matrix(data.frame(
        lapply(dat[, probability_cols, drop = FALSE], as.numeric),
        check.names = FALSE
      )))
      storage.mode(probability_matrix) <- "double"
      probability_confidence <- apply(probability_matrix, 1L, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) max(x) else NA_real_
      })
      fill <- !is.finite(confidence) & is.finite(probability_confidence)
      confidence[fill] <- probability_confidence[fill]
    }

    uncertainty <- rep(NA_real_, nrow(dat))
    if ("Classification_uncertainty" %in% names(dat)) {
      uncertainty <- suppressWarnings(as.numeric(dat[["Classification_uncertainty"]]))
    }
    fill_uncertainty <- !is.finite(uncertainty) & is.finite(confidence)
    uncertainty[fill_uncertainty] <- 1 - confidence[fill_uncertainty]

    confidence_train <- confidence[train]
    uncertainty_train <- uncertainty[train]
    confidence_finite <- confidence_train[is.finite(confidence_train)]
    uncertainty_finite <- uncertainty_train[is.finite(uncertainty_train)]
    confidence_mean_se <- if (length(confidence_finite) > 1L) {
      stats::sd(confidence_finite) / sqrt(length(confidence_finite))
    } else {
      NA_real_
    }
    uncertainty_mean_se <- if (length(uncertainty_finite) > 1L) {
      stats::sd(uncertainty_finite) / sqrt(length(uncertainty_finite))
    } else {
      NA_real_
    }

    component <- c(
      "genetic_variance_not_identifiable",
      "residual_variance_not_identifiable",
      "heritability_not_identifiable",
      "mean_prediction_confidence",
      "prediction_confidence_variance",
      "mean_classification_uncertainty"
    )
    out <- data.frame(
      Component = component,
      Components = c(
        NA_real_, NA_real_, NA_real_,
        gp_vc_finite_mean(confidence_train),
        gp_vc_finite_var(confidence_train),
        gp_vc_finite_mean(uncertainty_train)
      ),
      Standard_error = c(
        NA_real_, NA_real_, NA_real_, confidence_mean_se, NA_real_, uncertainty_mean_se
      ),
      stringsAsFactors = FALSE
    )
    if (length(group_cols)) {
      for (nm in group_cols) {
        out[[nm]] <- as.character(dat[[nm]][[1L]])
      }
    }
    gp_format_variance_component_columns(
      out,
      include_trait = "Trait" %in% group_cols,
      include_env = "Env" %in% group_cols
    )
  })
  pieces <- pieces[vapply(pieces, is.data.frame, logical(1L))]
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (!length(pieces)) {
    return(gp_empty_variance_components())
  }
  out <- do.call(rbind, pieces)
  out <- gp_format_variance_component_columns(
    out,
    include_trait = "Trait" %in% group_cols,
    include_env = "Env" %in% group_cols
  )
  if (!length(group_cols)) {
    rownames(out) <- make.unique(out[["Component"]])
  }
  out
}

gp_predictive_model_variance_components <- function(pred,
                                                      response_family = "auto") {
  fam <- tryCatch(
    gp_normalize_response_family(response_family),
    error = function(e) "auto"
  )
  if (identical(fam, "auto")) {
    is_classification <- (exists("gp_is_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(pred))) ||
      (is.data.frame(pred) && any(grepl("^Probability_", names(pred))))
    fam <- if (isTRUE(is_classification)) "multiclass" else "gaussian"
  }
  if (identical(fam, "gaussian")) {
    return(gp_prediction_table_variance_components(pred))
  }
  gp_ml_classification_variance_components(pred)
}

gp_vc_model_from_parameters <- function(model_parameters) {
  if (!is.data.frame(model_parameters) ||
      !all(c("stat", "summary") %in% names(model_parameters))) {
    return(NULL)
  }
  stat <- tolower(trimws(as.character(model_parameters[["stat"]])))
  summary <- as.character(model_parameters[["summary"]])
  for (candidate in c(
    "requested_model", "gp_model_canonical", "gp_model", "model_type", "model"
  )) {
    idx <- which(stat == candidate & !is.na(summary) & nzchar(trimws(summary)))
    if (length(idx)) {
      return(summary[[idx[[1L]]]])
    }
  }
  NULL
}

gp_vc_resolve_model <- function(x, model = NULL) {
  explicit <- as.character(model %||% character())
  explicit <- explicit[!is.na(explicit) & nzchar(trimws(explicit))]
  if (length(explicit)) {
    return(explicit[[1L]])
  }
  if (!is.list(x)) {
    return(NULL)
  }
  for (nm in c("model_name", "GS_model", "model_type")) {
    candidate <- x[[nm]]
    if (is.atomic(candidate) && length(candidate) == 1L &&
        !is.na(candidate) && nzchar(trimws(as.character(candidate)))) {
      return(as.character(candidate))
    }
  }
  candidate <- gp_vc_model_from_parameters(x[["model_parameters"]])
  if (!is.null(candidate)) {
    return(candidate)
  }
  for (nm in c("bayes_result", "asreml_result", "gp_result", "result", "model_result")) {
    nested <- x[[nm]]
    if (is.list(nested)) {
      candidate <- gp_vc_resolve_model(nested)
      if (!is.null(candidate)) {
        return(candidate)
      }
    }
  }
  NULL
}

gp_vc_resolve_response_family <- function(x, response_family = NULL) {
  fam <- tryCatch(
    gp_normalize_response_family(response_family),
    error = function(e) "auto"
  )
  if (!identical(fam, "auto")) {
    return(fam)
  }
  if (is.list(x)) {
    params <- x[["model_parameters"]]
    if (is.data.frame(params) && all(c("stat", "summary") %in% names(params))) {
      stat <- tolower(trimws(as.character(params[["stat"]])))
      idx <- which(stat == "response_family")
      if (length(idx)) {
        candidate <- tryCatch(
          gp_normalize_response_family(params[["summary"]][[idx[[1L]]]]),
          error = function(e) "auto"
        )
        if (!identical(candidate, "auto")) {
          return(candidate)
        }
      }
    }
  }
  pred <- if (exists("gp_vc_extract_prediction_table", mode = "function")) {
    gp_vc_extract_prediction_table(x)
  } else if (is.list(x)) {
    x[["predicted_values"]] %||% x[["Predicted_value"]]
  } else {
    NULL
  }
  is_classification <- (exists("gp_is_classification_prediction_table", mode = "function") &&
    isTRUE(gp_is_classification_prediction_table(pred))) ||
    (is.data.frame(pred) && any(grepl("^Probability_", names(pred))))
  if (isTRUE(is_classification)) "multiclass" else "gaussian"
}

gp_ml_long_matrix <- function(pred,
                              id_col,
                              dimension_col,
                              value_col,
                              unit_cols = id_col) {
  if (!is.data.frame(pred) ||
      !all(c(unit_cols, dimension_col, value_col) %in% names(pred))) {
    return(NULL)
  }
  dat <- pred[, unique(c(unit_cols, dimension_col, value_col)), drop = FALSE]
  dat[[dimension_col]] <- as.character(dat[[dimension_col]])
  dat[[value_col]] <- suppressWarnings(as.numeric(dat[[value_col]]))
  keep <- is.finite(dat[[value_col]]) & !is.na(dat[[dimension_col]]) &
    nzchar(dat[[dimension_col]])
  dat <- dat[keep, , drop = FALSE]
  if (!nrow(dat) || length(unique(dat[[dimension_col]])) < 2L) {
    return(NULL)
  }
  dat[[".unit"]] <- do.call(
    paste,
    c(lapply(dat[, unit_cols, drop = FALSE], as.character), sep = "\r")
  )
  agg <- stats::aggregate(
    dat[[value_col]],
    by = list(unit = dat[[".unit"]], dimension = dat[[dimension_col]]),
    FUN = mean,
    na.rm = TRUE
  )
  names(agg)[3L] <- "value"
  wide <- stats::reshape(
    agg,
    idvar = "unit",
    timevar = "dimension",
    direction = "wide"
  )
  value_names <- paste0("value.", sort(unique(agg[["dimension"]])))
  missing_names <- setdiff(value_names, names(wide))
  for (nm in missing_names) {
    wide[[nm]] <- NA_real_
  }
  mat <- as.matrix(wide[, value_names, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- sub("^value\\.", "", value_names)
  rownames(mat) <- as.character(wide[["unit"]])
  mat
}

gp_ml_covariance_pair <- function(mat) {
  mat <- tryCatch(as.matrix(mat), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) < 2L || ncol(mat) < 2L) {
    return(NULL)
  }
  covariance <- suppressWarnings(stats::cov(mat, use = "pairwise.complete.obs"))
  if (!is.matrix(covariance)) {
    return(NULL)
  }
  variance <- diag(covariance)
  denom <- sqrt(outer(variance, variance))
  correlation <- covariance / denom
  diag(correlation)[is.finite(variance) & variance > 0] <- 1
  correlation[!is.finite(correlation)] <- NA_real_
  list(covariance = covariance, correlation = correlation)
}

gp_ml_covariance_bootstrap_se <- function(mat,
                                          n_bootstrap = 200L,
                                          random_seed = 123L) {
  pair <- gp_ml_covariance_pair(mat)
  if (is.null(pair)) {
    return(list(covariance_se = NULL, correlation_se = NULL))
  }
  n_bootstrap <- suppressWarnings(as.integer(n_bootstrap[[1L]]))
  if (!is.finite(n_bootstrap) || n_bootstrap < 2L || nrow(mat) < 3L) {
    na_mat <- pair$covariance
    na_mat[] <- NA_real_
    return(list(covariance_se = na_mat, correlation_se = na_mat))
  }
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(random_seed[[1L]]))

  p <- ncol(mat)
  cov_draws <- matrix(NA_real_, nrow = n_bootstrap, ncol = p * p)
  cor_draws <- matrix(NA_real_, nrow = n_bootstrap, ncol = p * p)
  for (b in seq_len(n_bootstrap)) {
    idx <- sample.int(nrow(mat), nrow(mat), replace = TRUE)
    draw <- gp_ml_covariance_pair(mat[idx, , drop = FALSE])
    if (!is.null(draw)) {
      cov_draws[b, ] <- as.vector(draw$covariance)
      cor_draws[b, ] <- as.vector(draw$correlation)
    }
  }
  finite_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1L) stats::sd(x) else NA_real_
  }
  covariance_se <- matrix(
    apply(cov_draws, 2L, finite_sd),
    nrow = p,
    dimnames = dimnames(pair$covariance)
  )
  correlation_se <- matrix(
    apply(cor_draws, 2L, finite_sd),
    nrow = p,
    dimnames = dimnames(pair$correlation)
  )
  list(covariance_se = covariance_se, correlation_se = correlation_se)
}

gp_ml_predictive_matrix_bundle <- function(pred,
                                           id_col,
                                           dimension_col,
                                           error_pred = NULL,
                                           n_bootstrap = 200L,
                                           random_seed = 123L) {
  suffix <- if (grepl("trait", dimension_col, ignore.case = TRUE)) {
    "traits"
  } else {
    "environments"
  }
  pred_mat <- gp_ml_long_matrix(
    pred = pred,
    id_col = id_col,
    dimension_col = dimension_col,
    value_col = "Predicted_value"
  )
  pred_pair <- gp_ml_covariance_pair(pred_mat)
  if (is.null(pred_pair)) {
    return(list())
  }
  pred_se <- gp_ml_covariance_bootstrap_se(
    pred_mat,
    n_bootstrap = n_bootstrap,
    random_seed = random_seed
  )
  out <- stats::setNames(
    list(
      pred_pair$covariance,
      pred_pair$correlation,
      pred_se$covariance_se,
      pred_se$correlation_se
    ),
    c(
      paste0("Prediction_covariance_", suffix),
      paste0("Prediction_correlation_", suffix),
      paste0("Prediction_covariance_", suffix, "_SE"),
      paste0("Prediction_correlation_", suffix, "_SE")
    )
  )

  error_pred <- error_pred %||% NULL
  if (is.data.frame(error_pred) &&
      all(c("Observed_value", "Predicted_value") %in% names(error_pred))) {
    err <- error_pred
    err[["Prediction_error"]] <- suppressWarnings(
      as.numeric(err[["Observed_value"]]) - as.numeric(err[["Predicted_value"]])
    )
    unit_cols <- id_col
    if ("rep" %in% names(err)) {
      unit_cols <- c(unit_cols, "rep")
    }
    err_mat <- gp_ml_long_matrix(
      pred = err,
      id_col = id_col,
      dimension_col = dimension_col,
      value_col = "Prediction_error",
      unit_cols = unit_cols
    )
    err_pair <- gp_ml_covariance_pair(err_mat)
    if (!is.null(err_pair)) {
      err_se <- gp_ml_covariance_bootstrap_se(
        err_mat,
        n_bootstrap = n_bootstrap,
        random_seed = as.integer(random_seed) + 1L
      )
      out[[paste0("Prediction_error_covariance_", suffix)]] <- err_pair$covariance
      out[[paste0("Prediction_error_correlation_", suffix)]] <- err_pair$correlation
      out[[paste0("Prediction_error_covariance_", suffix, "_SE")]] <- err_se$covariance_se
      out[[paste0("Prediction_error_correlation_", suffix, "_SE")]] <- err_se$correlation_se
    }
  }
  out
}

gp_ml_apply_cv_uncertainty <- function(pred,
                                       cv_pred,
                                       id_col,
                                       group_col = "Trait",
                                       confidence_level = 0.95,
                                       high_reliability_thres = 0.9,
                                       low_reliability_thres = 0.5,
                                       target_pred = NULL) {
  required <- c(id_col, group_col, "Observed_value", "Predicted_value")
  if (!is.data.frame(pred) || !nrow(pred) || !is.data.frame(cv_pred) ||
      !nrow(cv_pred) || !all(required %in% names(pred)) ||
      !all(required %in% names(cv_pred))) {
    return(pred)
  }
  out <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)
  cv <- as.data.frame(cv_pred, stringsAsFactors = FALSE, check.names = FALSE)
  cv[["Observed_value"]] <- suppressWarnings(as.numeric(cv[["Observed_value"]]))
  cv[["Predicted_value"]] <- suppressWarnings(as.numeric(cv[["Predicted_value"]]))
  cv <- cv[
    is.finite(cv[["Observed_value"]]) & is.finite(cv[["Predicted_value"]]),
    ,
    drop = FALSE
  ]
  if (!nrow(cv)) {
    return(out)
  }
  cv[[".key"]] <- paste(
    as.character(cv[[id_col]]),
    as.character(cv[[group_col]]),
    sep = "\r"
  )
  out_key <- paste(
    as.character(out[[id_col]]),
    as.character(out[[group_col]]),
    sep = "\r"
  )

  key_split <- split(seq_len(nrow(cv)), cv[[".key"]], drop = TRUE)
  key_summary <- lapply(key_split, function(idx) {
    values <- cv[["Predicted_value"]][idx]
    observed <- cv[["Observed_value"]][idx]
    mean_prediction <- gp_vc_finite_mean(values)
    observed_value <- gp_vc_finite_mean(observed)
    data.frame(
      key = cv[[".key"]][idx[[1L]]],
      group = as.character(cv[[group_col]][idx[[1L]]]),
      repeated_prediction_variance = gp_vc_finite_var(values),
      mean_prediction = mean_prediction,
      observed_value = observed_value,
      prediction_error_mse = (observed_value - mean_prediction)^2,
      absolute_prediction_error = abs(observed_value - mean_prediction),
      stringsAsFactors = FALSE
    )
  })
  key_summary <- do.call(rbind, key_summary)
  rownames(key_summary) <- NULL

  # Cross-fitted evaluation rows calibrate marginal prediction error, while
  # repeated all-target predictions from the same fold fits provide the
  # target-specific resampling variance.  Keeping these roles separate lets
  # genuinely unobserved Test GIDs receive SE/PEV without using their outcomes.
  target_summary <- NULL
  if (is.data.frame(target_pred) && nrow(target_pred) &&
      all(c(id_col, group_col, "Predicted_value") %in% names(target_pred))) {
    target <- as.data.frame(target_pred, stringsAsFactors = FALSE, check.names = FALSE)
    target[["Predicted_value"]] <- suppressWarnings(as.numeric(target[["Predicted_value"]]))
    target <- target[is.finite(target[["Predicted_value"]]), , drop = FALSE]
    if (nrow(target)) {
      target[[".key"]] <- paste(
        as.character(target[[id_col]]),
        as.character(target[[group_col]]),
        sep = "\r"
      )
      target_split <- split(seq_len(nrow(target)), target[[".key"]], drop = TRUE)
      target_summary <- do.call(rbind, lapply(target_split, function(idx) {
        data.frame(
          key = target[[".key"]][idx[[1L]]],
          repeated_prediction_variance = gp_vc_finite_var(target[["Predicted_value"]][idx]),
          stringsAsFactors = FALSE
        )
      }))
      rownames(target_summary) <- NULL
    }
  }

  groups <- unique(as.character(cv[[group_col]]))
  out_group <- as.character(out[[group_col]])
  out_prediction <- suppressWarnings(as.numeric(out[["Predicted_value"]]))
  group_summary <- lapply(groups, function(group) {
    idx <- which(key_summary[["group"]] == group)
    group_mse <- gp_vc_finite_mean(key_summary[["prediction_error_mse"]][idx])
    reference_variance <- gp_vc_finite_var(key_summary[["observed_value"]][idx])
    out_idx <- which(out_group == group & is.finite(out_prediction))
    fitted_reference_variance <- gp_vc_finite_var(out_prediction[out_idx])
    fitted_reference_source <- "all_final_prediction_variance_within_group"
    if (!is.finite(fitted_reference_variance) || fitted_reference_variance <= 0) {
      fitted_reference_variance <- reference_variance
      fitted_reference_source <- "training_phenotype_variance_fallback"
    }
    if (!is.finite(fitted_reference_variance) || fitted_reference_variance <= 0) {
      fitted_reference_variance <- 1
      fitted_reference_source <- "unit_variance_fallback"
    }
    resampling_variance <- gp_vc_finite_mean(
      key_summary[["repeated_prediction_variance"]][idx]
    )
    interval_radius <- gp_ml_finite_sample_quantile(
      key_summary[["absolute_prediction_error"]][idx],
      probability = confidence_level
    )
    data.frame(
      group = group,
      prediction_error_mse = group_mse,
      reference_variance = reference_variance,
      fitted_reference_variance = fitted_reference_variance,
      fitted_reference_source = fitted_reference_source,
      resampling_variance = resampling_variance,
      interval_radius = interval_radius,
      calibration_n = sum(is.finite(key_summary[["absolute_prediction_error"]][idx])),
      stringsAsFactors = FALSE
    )
  })
  group_summary <- do.call(rbind, group_summary)
  rownames(group_summary) <- NULL

  key_match <- match(out_key, key_summary[["key"]])
  group_match <- match(as.character(out[[group_col]]), group_summary[["group"]])
  raw_var <- if (is.data.frame(target_summary) && nrow(target_summary)) {
    target_summary[["repeated_prediction_variance"]][match(out_key, target_summary[["key"]])]
  } else {
    key_summary[["repeated_prediction_variance"]][key_match]
  }
  group_mse <- group_summary[["prediction_error_mse"]][group_match]
  pev <- rep(NA_real_, nrow(out))
  uncertainty_source <- rep(
    "unavailable_no_target_specific_predictive_uncertainty",
    nrow(out)
  )
  pev_basis <- rep(gp_ml_unavailable_target_uncertainty_estimand(), nrow(out))

  # A group-level held-out MSE is a valid marginal calibration diagnostic, but
  # repeating it for every target does not make it a target-specific PEV.  Use
  # one groupwise multiplicative calibration factor so that mean target PEV
  # matches held-out MSE while preserving the relative target-to-target
  # resampling uncertainty.  An additive common residual term would flatten
  # the SEs and recreate the misleading near-constant uncertainty pattern.
  for (group in unique(out_group)) {
    idx <- which(out_group == group)
    row_var <- raw_var[idx]
    valid <- is.finite(row_var) & row_var >= 0
    mse <- gp_vc_finite_mean(group_mse[idx])
    mean_row_var <- gp_vc_finite_mean(row_var[valid])
    if (!any(valid) || !is.finite(mse) || mse < 0 ||
        !is.finite(mean_row_var) || mean_row_var <= 0) {
      next
    }

    variance_scale <- mse / mean_row_var
    candidate_pev <- row_var[valid] * variance_scale
    candidate_pev <- pmax(candidate_pev, 0)
    candidate_se <- sqrt(candidate_pev)
    if (!gp_ml_has_substantive_target_variation(candidate_se)) {
      next
    }

    target_idx <- idx[valid]
    pev[target_idx] <- candidate_pev
    uncertainty_source[target_idx] <-
      "heldout_cross_fitted_target_resampling_calibrated"
    pev_basis[target_idx] <- paste(
      "target-specific repeated-prediction variance multiplicatively",
      "calibrated so its group mean equals held-out predictive MSE;",
      "not mixed-model PEV"
    )
  }
  se <- sqrt(pev)
  reference_variance <- group_summary[["reference_variance"]][group_match]
  reliability_reference_variance <-
    group_summary[["fitted_reference_variance"]][group_match]
  reliability_reference_source <-
    group_summary[["fitted_reference_source"]][group_match]
  reliability_variance_input <- raw_var
  group_resampling_variance <- group_summary[["resampling_variance"]][group_match]
  reliability_source <- ifelse(
    is.finite(reliability_variance_input) & reliability_variance_input >= 0,
    "repeated_prediction_se_squared_marker_adjustment",
    "group_mean_repeated_prediction_se_squared_marker_adjustment"
  )
  missing_reliability_variance <-
    !is.finite(reliability_variance_input) | reliability_variance_input < 0
  reliability_variance_input[missing_reliability_variance] <-
    group_resampling_variance[missing_reliability_variance]
  missing_reliability_variance <-
    !is.finite(reliability_variance_input) | reliability_variance_input < 0
  min_reliability_variance <- pmax(
    reliability_reference_variance * 1e-8,
    .Machine$double.eps
  )
  reliability_variance_input[missing_reliability_variance] <-
    min_reliability_variance[missing_reliability_variance]
  reliability_source[missing_reliability_variance] <-
    "minimal_resampling_variance_fallback_marker_adjustment"
  interval_radius <- group_summary[["interval_radius"]][group_match]
  calibration_n <- group_summary[["calibration_n"]][group_match]

  stability <- gp_ml_gaussian_stability(
    prediction_error_var = reliability_variance_input,
    reference_variance = reliability_reference_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  predictive_precision <- gp_ml_predictive_reliability(
    prediction_error_var = pev,
    reference_variance = reference_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  pred_value <- suppressWarnings(as.numeric(out[["Predicted_value"]]))
  out[["Standard_error"]] <- se
  out[["PEV"]] <- pev
  out[["Prediction_error_variance"]] <- pev
  out[["lower_bound"]] <- pred_value - interval_radius
  out[["upper_bound"]] <- pred_value + interval_radius
  out[["Uncertainty"]] <- 2 * interval_radius
  out[["Uncertainty_remarks"]] <-
    "groupwise_heldout_cross_fitted_marginal_interval"
  out[["Prediction_stability"]] <- stability$stability
  out[["Prediction_stability_remarks"]] <- stability$remarks
  out[["Prediction_stability_reference_variance"]] <- reliability_reference_variance
  out[["Reliability_variance_input"]] <- predictive_precision$variance_input
  out[["Reliability"]] <- predictive_precision$reliability
  out[["Reliability_remarks"]] <- predictive_precision$remarks
  out[["Reliability_reference_variance"]] <- reference_variance
  out[["Reliability_basis"]] <- predictive_precision$basis
  out[["Prediction_uncertainty_source"]] <- uncertainty_source
  out[["Repeated_prediction_variance"]] <- raw_var
  out[["PEV_basis"]] <- pev_basis
  out[["Prediction_interval_method"]] <-
    "groupwise cross-fitted absolute-residual interval with finite-sample quantile correction"
  out[["Prediction_interval_nominal_coverage"]] <- confidence_level
  out[["Prediction_interval_calibration_n"]] <- calibration_n
  out
}

gp_vc_find_named_data_frame <- function(x, names_to_find) {
  if (!is.list(x)) {
    return(NULL)
  }
  for (nm in names_to_find) {
    obj <- x[[nm]]
    if (is.data.frame(obj) && nrow(obj) > 0L) {
      return(as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE))
    }
  }
  for (obj in x) {
    if (is.list(obj)) {
      hit <- gp_vc_find_named_data_frame(obj, names_to_find)
      if (is.data.frame(hit) && nrow(hit) > 0L) {
        return(hit)
      }
    }
  }
  NULL
}

gp_vc_is_gp_result <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  }
  if (any(c("gp_result", "gp_info", "raw_python_result") %in% names(x))) {
    return(TRUE)
  }
  diagnostics <- x[["diagnostics"]]
  if (is.list(diagnostics)) {
    method <- tolower(as.character(diagnostics[["method"]] %||% diagnostics[["model"]] %||% ""))
    if (length(method) && any(grepl("gp|gaussian|kernel|lowrank", method))) {
      return(TRUE)
    }
  }
  for (nm in c("result", "model_result")) {
    nested <- x[[nm]]
    if (is.list(nested) && isTRUE(gp_vc_is_gp_result(nested))) {
      return(TRUE)
    }
  }
  FALSE
}

gp_vc_extract_prediction_table <- function(x) {
  if (!is.list(x)) {
    return(NULL)
  }
  for (nm in c("predicted_values", "Predicted_value", "predictions")) {
    obj <- x[[nm]]
    if (is.data.frame(obj)) {
      return(obj)
    }
  }
  for (nm in c("bayes_result", "asreml_result", "gp_result", "result")) {
    nested <- x[[nm]]
    if (is.list(nested)) {
      hit <- gp_vc_extract_prediction_table(nested)
      if (is.data.frame(hit)) {
        return(hit)
      }
    }
  }
  NULL
}

gp_vc_is_multi_environment <- function(x, heter_groups = NULL) {
  pred <- gp_vc_extract_prediction_table(x)
  if (!is.data.frame(pred)) {
    return(FALSE)
  }
  # Phase 3.27: respect user-supplied heter_groups name. Phase 3.22 renamed
  # the public-output env col from "Env" to the user's heter_groups (e.g.
  # "YYY"); a hardcoded candidate list silently misdetects MET as single-env.
  candidates <- unique(stats::na.omit(c(
    as.character(heter_groups),
    "Env", "env", "Environment", "environment", "Location", "location"
  )))
  candidates <- candidates[nzchar(candidates)]
  env_col <- gp_contract_first_name(names(pred), candidates)
  !is.na(env_col) && gp_contract_multiple_nonempty_values(pred[[env_col]])
}

gp_vc_numeric_column <- function(df, candidates) {
  hit <- gp_contract_first_name(names(df), candidates)
  if (is.na(hit)) {
    return(rep(NA_real_, nrow(df)))
  }
  gp_contract_as_numeric(df[[hit]])
}

gp_vc_component_labels <- function(df) {
  comp_col <- gp_contract_first_name(names(df), c("Component", "component", "term", "source"))
  if (is.na(comp_col)) {
    rep("", nrow(df))
  } else {
    tolower(as.character(df[[comp_col]]))
  }
}

gp_vc_component_types <- function(df) {
  type_col <- gp_contract_first_name(names(df), c("component_type", "Component_type", "type"))
  if (is.na(type_col)) {
    rep("", nrow(df))
  } else {
    tolower(as.character(df[[type_col]]))
  }
}

gp_vc_component_kernel <- function(df) {
  kernel_col <- gp_contract_first_name(names(df), c("kernel", "Kernel"))
  if (is.na(kernel_col)) {
    rep("", nrow(df))
  } else {
    tolower(as.character(df[[kernel_col]]))
  }
}

#' Locate the FA delta-method summary dict on a GP result list.
#'
#' Returns the dict written by `varcomp_asreml.compute_fa_derived_summary` and
#' surfaced via `summary[["fa_derived"]]` by `build_fa_icm_ai_varcomp`. Returns
#' NULL when not present (i.e. the model is not an FA fit).
#'
#' @keywords internal
#' @noRd
gp_vc_extract_fa_derived <- function(x) {
  is_match <- function(obj) {
    is.list(obj) && "sigma2_g_per_env" %in% names(obj)
  }
  paths <- list(
    c("fa_derived"),
    c("summary", "fa_derived"),
    c("gp_result", "summary", "fa_derived"),
    c("result", "summary", "fa_derived"),
    c("gp_info", "summary", "fa_derived")
  )
  for (path in paths) {
    obj <- x; ok <- TRUE
    for (key in path) {
      if (!is.list(obj) || is.null(obj[[key]])) { ok <- FALSE; break }
      obj <- obj[[key]]
    }
    if (isTRUE(ok) && isTRUE(is_match(obj))) {
      return(obj)
    }
  }
  NULL
}

#' Find the AI-derived Cov(sigma2_g, sigma2_e) dict from a GP result
#'
#' Returns a list with named numeric scalars `var_g`, `var_e`, `cov_g_e`
#' (matching the Python schema in `varcomp_asreml.compute_genetic_residual_covariance`),
#' or NULL when no such structure is available. The Python backend attaches
#' this dict under `summary[["genetic_residual_covariance"]]` on the AI-REML
#' output of `build_gp_exact_ai_varcomp` (and the generic
#' `build_varcomp_outputs`) for non-FA models that carry a `var_kernel`
#' (genetic) + `resid_env` (residual) parameterisation.
#'
#' @keywords internal
#' @noRd
gp_vc_extract_genetic_residual_covariance <- function(x) {
  required <- c("var_g", "var_e", "cov_g_e")
  is_match <- function(obj) {
    is.list(obj) && all(required %in% names(obj))
  }
  paths <- list(
    c("summary", "genetic_residual_covariance"),
    c("gp_result", "summary", "genetic_residual_covariance"),
    c("result", "summary", "genetic_residual_covariance"),
    c("gp_info", "summary", "genetic_residual_covariance"),
    c("genetic_residual_covariance")
  )
  for (path in paths) {
    obj <- x
    ok <- TRUE
    for (key in path) {
      if (!is.list(obj) || is.null(obj[[key]])) { ok <- FALSE; break }
      obj <- obj[[key]]
    }
    if (isTRUE(ok) && isTRUE(is_match(obj))) {
      return(obj)
    }
  }
  NULL
}

gp_vc_find_scalar <- function(x, candidates) {
  if (!is.list(x)) {
    return(NULL)
  }
  for (nm in candidates) {
    if (!is.null(x[[nm]]) && length(x[[nm]]) > 0L && !is.list(x[[nm]])) {
      val <- as.character(x[[nm]][[1L]])
      if (!is.na(val) && nzchar(val)) {
        return(val)
      }
    }
  }
  for (obj in x) {
    if (is.list(obj)) {
      hit <- gp_vc_find_scalar(obj, candidates)
      if (!is.null(hit)) {
        return(hit)
      }
    }
  }
  NULL
}

gp_vc_gp_estimation_method <- function(x) {
  mode <- gp_vc_find_scalar(
    x,
    c("gp_varcomp_mode", "varcomp_mode", "variance_component_method", "Estimation_method", "estimation_method")
  )
  mode <- tolower(as.character(mode %||% ""))
  if (grepl("reml", mode)) {
    return("REML")
  }
  if (grepl("^mom$|method_of_moments|moment", mode)) {
    return("method_of_moments")
  }
  if (grepl("krr_gcv_variance_ratio", mode, fixed = TRUE)) {
    return("KRR_GCV_variance_ratio")
  }
  "backend_reported"
}

gp_vc_apply_estimation_method <- function(df, estimation_method) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(df)
  }
  if (!"Estimation_method" %in% names(df)) {
    df[["Estimation_method"]] <- estimation_method
  } else {
    existing <- as.character(df[["Estimation_method"]])
    missing <- is.na(existing) | !nzchar(existing) |
      existing %in% c("covariance_diagonal", "derived_from_covariance_diagonals")
    df[["Estimation_method"]][missing] <- estimation_method
  }
  df
}

gp_vc_kind_mask <- function(df, kind = c("genetic", "residual")) {
  kind <- match.arg(kind)
  labels <- gp_vc_component_labels(df)
  types <- gp_vc_component_types(df)

  residual <- grepl("sigma2_noise|!r$|(^|[_:])r$|resid|residual|noise", labels) |
    grepl("resid|residual|noise", types)
  gxe <- grepl("w_ge|g:env|gxe|interaction", labels) |
    grepl("interaction", types)
  env <- !residual & (grepl("^w_e$|environment|^env$|env_main", labels) |
                        grepl("environment", types))
  genetic <- !residual & !gxe & !env &
    (grepl("^g$|w_g|genetic|vm\\(|gid|geno|kernel|gmatrix", labels) |
       grepl("genetic_main", types))

  if (identical(kind, "genetic")) genetic else residual
}

gp_vc_boundary_mask <- function(df, zero_tol = sqrt(.Machine$double.eps)) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(logical(0L))
  }
  bound_col <- gp_contract_first_name(
    names(df),
    c("bound", "Bound", "boundary", "Boundary", "constraint", "Constraint")
  )
  bound <- if (!is.na(bound_col)) {
    toupper(trimws(as.character(df[[bound_col]])))
  } else {
    rep("", nrow(df))
  }
  estimate <- gp_vc_numeric_column(
    df,
    c("estimate", "Estimate", "Components", "component", "value", "Variance")
  )
  se <- gp_vc_numeric_column(
    df,
    c("std.error", "std_error", "Std.Error", "SE", "se", "standard_error", "Standard_error")
  )
  boundary_flag <- bound %in% c("B", "F", "?", "S")
  zero_without_se <- is.finite(estimate) &
    abs(estimate) <= as.numeric(zero_tol)[1L] &
    !is.finite(se)
  boundary_flag | zero_without_se
}

gp_vc_boundary_status <- function(df, zero_tol = sqrt(.Machine$double.eps)) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(character())
  }
  bound_col <- gp_contract_first_name(
    names(df),
    c("bound", "Bound", "boundary", "Boundary", "constraint", "Constraint")
  )
  bound <- if (!is.na(bound_col)) {
    toupper(trimws(as.character(df[[bound_col]])))
  } else {
    rep("", nrow(df))
  }
  status <- rep("", nrow(df))
  status[bound == "B"] <- "bound_at_zero"
  status[bound == "F"] <- "fixed"
  status[bound == "C"] <- "constrained"
  status[bound == "?"] <- "unstable"
  status[bound == "S"] <- "singular"

  estimate <- gp_vc_numeric_column(
    df,
    c("estimate", "Estimate", "Components", "component", "value", "Variance")
  )
  se <- gp_vc_numeric_column(
    df,
    c("std.error", "std_error", "Std.Error", "SE", "se", "standard_error", "Standard_error")
  )
  inferred <- !nzchar(status) & is.finite(estimate) &
    abs(estimate) <= as.numeric(zero_tol)[1L] & !is.finite(se)
  status[inferred] <- "bound_at_zero"
  status
}

# Residual rows estimated AT the zero boundary. Genetic and residual variance
# are then not separable (the genetic part absorbs the residual), so genetic
# variance, residual variance, heritability and reliability are reported as
# not estimable for that trait / environment.
gp_vc_residual_at_zero <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(logical(0L))
  }
  gp_vc_kind_mask(df, "residual") & gp_vc_boundary_status(df) == "bound_at_zero"
}

gp_vc_component_boundary <- function(df, kind = c("genetic", "residual")) {
  kind <- match.arg(kind)
  if (!is.data.frame(df) || !nrow(df)) {
    return(FALSE)
  }
  keep <- gp_vc_kind_mask(df, kind)
  any(keep & gp_vc_boundary_mask(df), na.rm = TRUE)
}

gp_vc_all_components_boundary <- function(df, kind = c("genetic", "residual")) {
  kind <- match.arg(kind)
  if (!is.data.frame(df) || !nrow(df)) {
    return(FALSE)
  }
  keep <- gp_vc_kind_mask(df, kind)
  if (!any(keep)) {
    return(FALSE)
  }
  all(gp_vc_boundary_mask(df)[keep])
}

gp_vc_boundary_component_names <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(character())
  }
  comp_col <- gp_contract_first_name(names(df), c("Component", "component", "term", "source"))
  comp <- if (!is.na(comp_col)) {
    as.character(df[[comp_col]])
  } else {
    rownames(df)
  }
  comp[gp_vc_boundary_mask(df) & !is.na(comp) & nzchar(comp)]
}

gp_vc_collapse_component <- function(df, kind = c("genetic", "residual")) {
  kind <- match.arg(kind)
  if (!is.data.frame(df) || !nrow(df)) {
    return(c(estimate = NA_real_, se = NA_real_))
  }
  estimate <- gp_vc_numeric_column(df, c("estimate", "Estimate", "Components", "component", "value", "Variance"))
  se <- gp_vc_numeric_column(df, c("std.error", "std_error", "Std.Error", "SE", "se", "standard_error", "Standard_error"))
  keep <- gp_vc_kind_mask(df, kind)

  if (identical(kind, "genetic") && any(keep)) {
    kernel <- gp_vc_component_kernel(df)
    total_keep <- keep & kernel == "total"
    if (any(total_keep)) {
      keep <- total_keep
    }
  }

  vals <- estimate[keep & is.finite(estimate)]
  ses <- se[keep & is.finite(se)]
  c(
    estimate = if (length(vals)) sum(vals) else NA_real_,
    se = if (length(ses)) sqrt(sum(ses^2)) else NA_real_
  )
}

gp_vc_has_component_estimate <- function(df, kind = c("genetic", "residual")) {
  kind <- match.arg(kind)
  if (!is.data.frame(df) || !nrow(df)) {
    return(FALSE)
  }
  estimate <- gp_vc_numeric_column(df, c("estimate", "Estimate", "Components", "component", "value", "Variance"))
  keep <- gp_vc_kind_mask(df, kind)
  any(keep & is.finite(estimate))
}

gp_vc_reml_incomplete_table <- function(reason = "reml_variance_components_incomplete") {
  gp_unavailable_variance_components(
    reason = reason,
    estimation_method = "REML"
  )
}

gp_format_gp_multi_env_variance_components <- function(x, heter_groups = NULL) {
  if (!is.list(x) || !isTRUE(gp_vc_is_gp_result(x)) ||
      !isTRUE(gp_vc_is_multi_environment(x, heter_groups = heter_groups))) {
    return(gp_empty_variance_components())
  }

  estimation_method <- gp_vc_gp_estimation_method(x)
  strict_reml <- identical(estimation_method, "REML")
  varcomp <- gp_vc_find_named_data_frame(x, c("varcomp", "var_components_ai"))
  if (!is.data.frame(varcomp) || !nrow(varcomp)) {
    if (isTRUE(strict_reml)) {
      return(gp_vc_reml_incomplete_table("reml_variance_components_not_returned"))
    }
    return(gp_empty_variance_components())
  }
  if (isTRUE(strict_reml) &&
      (!gp_vc_has_component_estimate(varcomp, "genetic") ||
       !gp_vc_has_component_estimate(varcomp, "residual"))) {
    return(gp_vc_reml_incomplete_table())
  }

  out <- gp_format_variance_components(varcomp, default_estimation_method = estimation_method)
  if (!is.data.frame(out) || !nrow(out)) {
    if (isTRUE(strict_reml)) {
      return(gp_vc_reml_incomplete_table())
    }
    return(gp_empty_variance_components())
  }
  raw_comp_col <- gp_contract_first_name(names(varcomp), c("Component", "component", "term", "source"))
  raw_component <- if (!is.na(raw_comp_col)) {
    as.character(varcomp[[raw_comp_col]])
  } else {
    rownames(varcomp)
  }
  raw_status <- gp_vc_boundary_status(varcomp)
  row_status <- if (length(raw_component) && length(raw_status)) {
    s <- raw_status[match(as.character(out[["Component"]]), raw_component)]
    s[is.na(s)] <- ""
    s
  } else {
    rep("", nrow(out))
  }
  residual_not_estimable <- rep(FALSE, nrow(out))
  if (isTRUE(strict_reml)) {
    boundary_components <- gp_vc_boundary_component_names(varcomp)
    if (length(boundary_components)) {
      boundary_rows <- out[["Component"]] %in% boundary_components
      # A variance estimated at the zero boundary is a valid REML estimate
      # (e.g. no GxE, or an FA environment-specific variance of 0): keep the
      # 0, flag it "bound_at_zero", and withhold its degenerate SE. The
      # exception is a RESIDUAL at zero: genetic and residual variance are
      # then not separable, so the residual (and, per environment, genetic
      # variance and heritability) is reported as not estimable. Unstable,
      # singular or fixed parameters stay unavailable.
      residual_rows <- grepl("!R$|resid|residual", as.character(out[["Component"]]),
                             ignore.case = TRUE)
      residual_not_estimable <- boundary_rows & residual_rows & row_status == "bound_at_zero"
      genetic_kernel_rows <- boundary_rows &
        grepl("^vm\\(.*\\)!var$", as.character(out[["Component"]]))
      keep_value <- boundary_rows & !residual_not_estimable &
        (row_status == "bound_at_zero" | genetic_kernel_rows)
      out[["Standard_error"]][boundary_rows] <- NA_real_
      out[["Components"]][boundary_rows & !keep_value] <- NA_real_
    }
  }

  # Phase 3.29: surface an explicit boundary flag, independent of the estimation
  # method. A boundary-clipped genetic kernel keeps its value (~0) but NA SE
  # (REML path) or is returned as 0 / NA by the backend (KRR / non-REML path);
  # without a flag a user reading `vm(GID,omic1_kernel)!var = 0` cannot tell
  # whether the kernel was estimated AT the variance boundary or silently
  # dropped. The `Boundary` column ("bound_at_zero" vs "") makes that explicit.
  # It is a CONDITIONAL column -- added only when a boundary component exists --
  # mirroring the on-demand Trait/Env columns, so the pinned 3-column contract
  # for boundary-free tables is unchanged.
  # Preserve the raw optimizer status by component name. Falling back to the
  # formatted zero-without-SE heuristic covers non-REML backends that do not
  # return an explicit status column. Fixed, constrained, unstable and singular
  # parameters must not be mislabeled as estimates bound at zero.
  boundary_status <- row_status
  inferred_boundary <- !nzchar(boundary_status) & gp_vc_boundary_mask(out)
  boundary_status[inferred_boundary] <- "bound_at_zero"
  boundary_status[residual_not_estimable] <- "not_estimable"
  if (any(nzchar(boundary_status))) {
    out[["Boundary"]] <- boundary_status
  }
  out
}

gp_format_gp_single_env_variance_components <- function(x, heter_groups = NULL) {
  if (!is.list(x) || !isTRUE(gp_vc_is_gp_result(x)) ||
      isTRUE(gp_vc_is_multi_environment(x, heter_groups = heter_groups))) {
    return(gp_empty_variance_components())
  }

  estimation_method <- gp_vc_gp_estimation_method(x)
  strict_reml <- identical(estimation_method, "REML")
  varcomp <- gp_vc_find_named_data_frame(x, c("varcomp", "var_components_ai"))
  component_summary <- gp_vc_find_named_data_frame(x, c("var_components_summary", "var_components"))
  if (isTRUE(strict_reml) && (!is.data.frame(varcomp) || !nrow(varcomp))) {
    return(gp_vc_reml_incomplete_table("reml_variance_components_not_returned"))
  }

  genetic <- gp_vc_collapse_component(varcomp, "genetic")
  residual <- gp_vc_collapse_component(varcomp, "residual")
  genetic_summary <- gp_vc_collapse_component(component_summary, "genetic")

  # FA delta-method path: sigma2_g is NOT a direct REML parameter in factor-
  # analytic models, so varcomp carries only fa_psi / fa_lambda rows. The
  # principled sigma2_g + SE come from var_components_summary, where the
  # Python backend (gp_framework._attach_fa_derived_se_to_summary) writes a
  # std.error column derived by propagating uncertainty through (Lambda, Psi)
  # via J^T Sigma_theta J. Prefer those when available -- they bypass the
  # varcomp aggregation mis-match (loadings would otherwise be summed into
  # the variance) and the active-set z-floor SE suppression on fa_psi rows.
  genetic_from_summary <- is.finite(genetic_summary[["estimate"]]) &&
    is.finite(genetic_summary[["se"]])
  if (isTRUE(genetic_from_summary)) {
    genetic[["estimate"]] <- genetic_summary[["estimate"]]
    genetic[["se"]] <- genetic_summary[["se"]]
  } else if (!is.finite(genetic[["estimate"]]) && !isTRUE(strict_reml)) {
    genetic[["estimate"]] <- genetic_summary[["estimate"]]
    genetic[["se"]] <- genetic_summary[["se"]]
  }
  if (!is.finite(residual[["estimate"]]) && !isTRUE(strict_reml)) {
    residual_summary <- gp_vc_collapse_component(component_summary, "residual")
    residual[["estimate"]] <- residual_summary[["estimate"]]
    residual[["se"]] <- residual_summary[["se"]]
  }

  # Boundary detection still applies to varcomp rows, but only for components
  # we did not source from the (post-processed) delta-method summary. A
  # collapsed multi-kernel genetic total remains estimable when one kernel is
  # at zero and at least one other genetic kernel is in the interior. Only an
  # all-boundary genetic block makes the collapsed total unavailable.
  # Behaviour: the boundary component's estimate AND SE are NA'd (matches
  # the existing test contract and avoids letting a clipped 0 propagate
  # into a misleading h^2). The PARTNER component, when well-estimated, is
  # preserved -- previously the code wiped both, turning a single boundary
  # hit on residual into an all-NA table even when sigma2_g was fine.
  # Heritability is NA whenever either underlying component is at boundary,
  # so we don't claim h^2 = 0 or h^2 = 1 from a clipped value.
  genetic_boundary <- isTRUE(strict_reml) && !isTRUE(genetic_from_summary) &&
    gp_vc_all_components_boundary(varcomp, "genetic")
  residual_boundary <- isTRUE(strict_reml) && gp_vc_component_boundary(varcomp, "residual")
  # Residual estimated at zero: genetic and residual variance are not
  # separable, so genetic variance, residual variance and heritability are all
  # reported as not estimable (the genetic estimate would include the residual).
  if (isTRUE(strict_reml) && any(gp_vc_residual_at_zero(varcomp))) {
    out <- gp_vc_simple_table(
      genetic_variance = NA_real_,
      residual_variance = NA_real_,
      heritability = NA_real_,
      estimation_method = estimation_method
    )
    out[["Boundary"]] <- "not_estimable"
    return(out)
  }
  # Genetic variance estimated at zero (residual estimable): a valid REML
  # estimate. Report 0 (h2 = 0), flagged, with no SE.
  genetic_at_zero <- isTRUE(genetic_boundary) && {
    gmask <- gp_vc_kind_mask(varcomp, "genetic")
    any(gmask) && all(gp_vc_boundary_status(varcomp)[gmask] == "bound_at_zero")
  }
  if (isTRUE(genetic_at_zero)) {
    genetic[["estimate"]] <- 0
    genetic[["se"]] <- NA_real_
    genetic_boundary <- FALSE
  }
  if (isTRUE(genetic_boundary)) {
    genetic[["estimate"]] <- NA_real_
    genetic[["se"]] <- NA_real_
  }
  if (isTRUE(residual_boundary)) {
    residual[["estimate"]] <- NA_real_
    residual[["se"]] <- NA_real_
  }

  if (isTRUE(strict_reml) &&
      isTRUE(genetic_boundary) && isTRUE(residual_boundary) &&
      !is.finite(genetic[["estimate"]]) && !is.finite(residual[["estimate"]])) {
    return(gp_vc_simple_table(
      genetic_variance = NA_real_,
      residual_variance = NA_real_,
      heritability = NA_real_,
      genetic_se = NA_real_,
      residual_se = NA_real_,
      heritability_se = NA_real_,
      estimation_method = "REML_boundary"
    ))
  }

  if (isTRUE(strict_reml) &&
      (!is.finite(genetic[["estimate"]]) || !is.finite(residual[["estimate"]])) &&
      !isTRUE(genetic_boundary) &&
      !isTRUE(residual_boundary)) {
    return(gp_vc_reml_incomplete_table())
  }

  if (!is.finite(genetic[["estimate"]]) && !is.finite(residual[["estimate"]])) {
    return(gp_empty_variance_components())
  }

  total <- genetic[["estimate"]] + residual[["estimate"]]
  h2 <- if (is.finite(total) && total > 0) genetic[["estimate"]] / total else NA_real_
  h2_se <- NA_real_
  if (is.finite(h2) &&
      is.finite(genetic[["se"]]) &&
      is.finite(residual[["se"]])) {
    d_g <- residual[["estimate"]] / total^2
    d_r <- -genetic[["estimate"]] / total^2
    h2_se <- sqrt((d_g^2 * genetic[["se"]]^2) + (d_r^2 * residual[["se"]]^2))
    # Upgrade to the full 3-term delta-method if Cov(sigma2_g, sigma2_e) was
    # supplied by the Python backend. The 2-term formula above assumes
    # Cov = 0; including the off-diagonal entry of the active-set AI inverse
    # is the right thing for REML. Falls back to the 2-term value silently
    # when the cov is absent or NA, so non-AI-SE result paths are unaffected.
    g_r_cov <- gp_vc_extract_genetic_residual_covariance(x)
    if (!is.null(g_r_cov)) {
      var_g <- suppressWarnings(as.numeric(g_r_cov[["var_g"]]))
      var_e <- suppressWarnings(as.numeric(g_r_cov[["var_e"]]))
      cov_g_e <- suppressWarnings(as.numeric(g_r_cov[["cov_g_e"]]))
      # Phase 3.24 defensive guard: any of these may be `numeric(0)` when the
      # GP backend didn't return the off-diagonal cov entry (single-kernel GP
      # models like Kernel-GBLUP / KRR don't always populate cov_g_e).
      # is.finite(numeric(0)) is logical(0), which trips `&&` with
      # "missing value where TRUE/FALSE needed". Length-check first.
      if (length(var_g) == 1L && length(var_e) == 1L && length(cov_g_e) == 1L &&
          is.finite(var_g) && is.finite(var_e) && is.finite(cov_g_e)) {
        full_var <- (d_g^2 * var_g) + (d_r^2 * var_e) + (2 * d_g * d_r * cov_g_e)
        if (is.finite(full_var) && full_var > 0) {
          h2_se <- sqrt(full_var)
        }
      }
    }
  }

  out <- gp_vc_simple_table(
    genetic_variance = genetic[["estimate"]],
    residual_variance = residual[["estimate"]],
    heritability = h2,
    genetic_se = genetic[["se"]],
    residual_se = residual[["se"]],
    heritability_se = h2_se,
    estimation_method = estimation_method
  )
  if (isTRUE(genetic_at_zero)) {
    out[["Boundary"]] <- c("bound_at_zero", "", "")
  }
  out
}

gp_format_variance_components <- function(x, joint = FALSE, default_estimation_method = "model_reported") {
  if (is.null(x)) {
    return(gp_empty_variance_components())
  }
  vc <- if (is.data.frame(x)) {
    x
  } else if (is.list(x)) {
    gp_find_first_data_frame(
      x,
      c(
        "Variance_components", "variance_components", "varcomp",
        "var_components_ai", "var_components_summary", "var_components",
        "Residual_variance_by_environment",
        "bayes_result", "asreml_result", "gp_result"
      )
    )
  } else {
    NULL
  }
  if (!is.data.frame(vc)) {
    return(gp_empty_variance_components())
  }
  out <- as.data.frame(vc, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(out)) {
    return(gp_empty_variance_components())
  }

  comp_col <- gp_contract_first_name(names(out), c("Component", "component", "term", "source"))
  rn <- rownames(out)
  has_informative_rownames <- !is.null(rn) && length(rn) == nrow(out) &&
    !all(is.na(rn) | !nzchar(rn) | rn == as.character(seq_len(nrow(out))))
  component_col_is_estimate <- !is.na(comp_col) &&
    identical(comp_col, "component") &&
    any(is.finite(gp_contract_as_numeric(out[[comp_col]]))) &&
    isTRUE(has_informative_rownames)
  component <- if (!is.na(comp_col) && !isTRUE(component_col_is_estimate)) {
    as.character(out[[comp_col]])
  } else {
    if (isTRUE(has_informative_rownames)) {
      as.character(rn)
    } else if (nrow(out) == 3L) {
      c("genetic_variance", "residual_variance", "heritability")
    } else {
      paste0("component_", seq_len(nrow(out)))
    }
  }
  component[is.na(component) | !nzchar(component)] <- paste0("component_", which(is.na(component) | !nzchar(component)))

  if (!"Components" %in% names(out)) {
    estimate_col <- gp_contract_first_name(names(out), c("estimate", "Estimate", "value", "Variance", "variance", "varcomp"))
    if (!is.na(estimate_col)) {
      out[["Components"]] <- gp_contract_as_numeric(out[[estimate_col]])
    } else if ("component" %in% names(out) &&
               any(grepl("^[0-9.+-]", as.character(out[["component"]]))) &&
               (!identical(comp_col, "component") || isTRUE(component_col_is_estimate))) {
      out[["Components"]] <- gp_contract_as_numeric(out[["component"]])
    } else {
      out[["Components"]] <- NA_real_
    }
  } else {
    out[["Components"]] <- gp_contract_as_numeric(out[["Components"]])
  }

  if (!"Standard_error" %in% names(out)) {
    se_col <- gp_contract_first_name(names(out), c("std.error", "std_error", "Std.Error", "SE", "se", "standard_error"))
    out[["Standard_error"]] <- if (!is.na(se_col)) gp_contract_as_numeric(out[[se_col]]) else NA_real_
  } else {
    out[["Standard_error"]] <- gp_contract_as_numeric(out[["Standard_error"]])
  }

  has_trait <- "Trait" %in% names(out) && gp_vc_has_nonempty_labels(out[["Trait"]])
  has_env <- "Env" %in% names(out) && gp_vc_has_nonempty_labels(out[["Env"]])
  out[["Component"]] <- component

  method_col <- gp_contract_first_name(
    names(out),
    c("Estimation_method", "estimation_method", "Variance_component_method", "variance_component_method", "Method", "method")
  )
  out[["Estimation_method"]] <- if (!is.na(method_col)) {
    as.character(out[[method_col]])
  } else {
    rep(default_estimation_method, nrow(out))
  }
  out[["Estimation_method"]][is.na(out[["Estimation_method"]]) | !nzchar(out[["Estimation_method"]])] <-
    default_estimation_method

  out <- gp_format_variance_component_columns(out)
  if (!has_trait && !has_env) {
    rownames(out) <- make.unique(out[["Component"]])
  } else {
    rownames(out) <- NULL
  }
  out
}

gp_covariance_diag_components <- function(mat, component, label_name = "Trait", labels = NULL) {
  if (is.null(mat)) {
    return(NULL)
  }
  mat <- tryCatch(as.matrix(mat), error = function(e) NULL)
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) {
    return(NULL)
  }
  storage.mode(mat) <- "double"
  labs <- if (!is.null(labels)) as.character(unlist(labels, use.names = FALSE)) else NULL
  if (is.null(labs) || length(labs) != min(dim(mat)) || any(!nzchar(labs))) {
    labs <- colnames(mat)
  }
  if (is.null(labs) || length(labs) != min(dim(mat)) || any(!nzchar(labs))) {
    labs <- rownames(mat)
  }
  if (is.null(labs) || length(labs) != min(dim(mat)) || any(!nzchar(labs))) {
    labs <- paste0(tolower(label_name), seq_len(min(dim(mat))))
  }
  n <- min(length(labs), nrow(mat), ncol(mat))
  out <- data.frame(
    Component = rep(component, n),
    Components = diag(mat)[seq_len(n)],
    Standard_error = rep(NA_real_, n),
    Estimation_method = rep("covariance_diagonal", n),
    stringsAsFactors = FALSE
  )
  out[[label_name]] <- labs[seq_len(n)]
  gp_format_variance_component_columns(
    out,
    include_trait = identical(label_name, "Trait"),
    include_env = identical(label_name, "Env")
  )
}

gp_variance_components_from_covariances <- function(x) {
  if (!is.list(x)) {
    return(gp_empty_variance_components())
  }
  # The subprocess bridge keeps fitted covariance matrices under its private
  # `fit` payload.  Search recursively and prefer response-scale matrices so
  # the public trait variances are on the same scale as the predictions.
  find_covariance <- function(aliases) {
    if (exists("gp_find_named_object", mode = "function")) {
      return(gp_find_named_object(x, aliases))
    }
    hit <- aliases[aliases %in% names(x)][1L]
    if (!is.na(hit)) x[[hit]] else NULL
  }
  trait_labels <- find_covariance(c("trait_levels", "traits", "response"))
  genetic <- find_covariance(c(
    "Genetic_covariance_traits", "genetic_covariance", "trait_covariance_matrix",
    "Sigma_G_response_scale", "Sigma_G"
  ))
  gxe <- find_covariance(c(
    "Genetic_covariance_trait_environment", "gxe_covariance", "gxe_trait_covariance",
    "Sigma_GE_response_scale", "Sigma_GE"
  ))
  residual <- find_covariance(c(
    "Residual_covariance_traits", "residual_covariance", "residual_trait_covariance",
    "Sigma_eps_response_scale", "Sigma_eps"
  ))
  env <- find_covariance(c(
    "Genetic_covariance_environments", "environment_covariance", "env_covariance"
  ))

  genetic_piece <- gp_covariance_diag_components(
    genetic, "genetic_variance", "Trait", labels = trait_labels
  )
  gxe_piece <- gp_covariance_diag_components(
    gxe, "gxe_variance", "Trait", labels = trait_labels
  )
  residual_piece <- gp_covariance_diag_components(
    residual, "residual_variance", "Trait", labels = trait_labels
  )
  env_piece <- gp_covariance_diag_components(env, "environment_variance", "Env")

  heritability_piece <- NULL
  if (is.data.frame(genetic_piece) && nrow(genetic_piece) &&
      is.data.frame(residual_piece) && nrow(residual_piece) &&
      all(c("Trait", "Components") %in% names(genetic_piece)) &&
      all(c("Trait", "Components") %in% names(residual_piece))) {
    labels <- as.character(genetic_piece[["Trait"]])
    residual_match <- match(labels, as.character(residual_piece[["Trait"]]))
    vg <- suppressWarnings(as.numeric(genetic_piece[["Components"]]))
    # A joint MT-MET cell contains both the main genomic effect and its G x E
    # effect, so the heritability numerator must include both exported
    # covariance diagonals. Reliability already uses this same Vg + Vgxe
    # reference variance.
    vge <- rep(0, length(vg))
    has_gxe <- is.data.frame(gxe_piece) && nrow(gxe_piece) &&
      all(c("Trait", "Components") %in% names(gxe_piece))
    if (has_gxe) {
      gxe_match <- match(labels, as.character(gxe_piece[["Trait"]]))
      matched_vge <- suppressWarnings(as.numeric(gxe_piece[["Components"]][gxe_match]))
      valid_vge <- is.finite(matched_vge) & matched_vge >= 0
      vge[valid_vge] <- matched_vge[valid_vge]
    }
    ve <- suppressWarnings(as.numeric(residual_piece[["Components"]][residual_match]))
    total_genetic <- vg + vge
    denom <- total_genetic + ve
    h2 <- rep(NA_real_, length(vg))
    estimable <- is.finite(total_genetic) & total_genetic >= 0 &
      is.finite(ve) & ve >= 0 &
      is.finite(denom) & denom > 0
    h2[estimable] <- total_genetic[estimable] / denom[estimable]
    heritability_piece <- data.frame(
      Trait = labels,
      Component = "heritability",
      Components = h2,
      Standard_error = NA_real_,
      Estimation_method = if (has_gxe) {
        "derived_from_genetic_gxe_residual_covariance_diagonals"
      } else {
        "derived_from_covariance_diagonals"
      },
      stringsAsFactors = FALSE
    )
    heritability_piece <- gp_format_variance_component_columns(
      heritability_piece,
      include_trait = TRUE
    )
  }

  pieces <- list(
    genetic_piece,
    gxe_piece,
    residual_piece,
    heritability_piece,
    env_piece
  )
  pieces <- pieces[vapply(pieces, is.data.frame, logical(1L))]
  if (!length(pieces)) {
    return(gp_empty_variance_components())
  }
  include_trait <- any(vapply(pieces, function(df) "Trait" %in% names(df) && gp_vc_has_nonempty_labels(df[["Trait"]]), logical(1L)))
  include_env <- any(vapply(pieces, function(df) "Env" %in% names(df) && gp_vc_has_nonempty_labels(df[["Env"]]), logical(1L)))
  pieces <- lapply(pieces, gp_format_variance_component_columns, include_trait = include_trait, include_env = include_env)
  out <- do.call(rbind, pieces)
  gp_format_variance_component_columns(out, include_trait = include_trait, include_env = include_env)
}

gp_empty_covariance_components <- function() {
  data.frame(
    Matrix = character(),
    Component = character(),
    Row = character(),
    Column = character(),
    Value = numeric(),
    Estimation_method = character(),
    stringsAsFactors = FALSE
  )
}

gp_as_numeric_square_matrix <- function(x) {
  if (is.data.frame(x)) {
    numeric_df <- all(vapply(x, is.numeric, logical(1L)))
    if (numeric_df && nrow(x) == ncol(x)) {
      x <- as.matrix(x)
    }
  }
  if (!is.matrix(x) || !nrow(x) || !ncol(x)) {
    return(NULL)
  }
  out <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(out)) {
    return(NULL)
  }
  storage.mode(out) <- "double"
  out
}

gp_find_named_object <- function(x, aliases, depth = 0L, max_depth = 4L) {
  if (!is.list(x) || depth > max_depth) {
    return(NULL)
  }
  hit <- aliases[aliases %in% names(x)][1L]
  if (!is.na(hit)) {
    return(x[[hit]])
  }
  for (obj in x) {
    if (is.list(obj)) {
      found <- gp_find_named_object(obj, aliases, depth = depth + 1L, max_depth = max_depth)
      if (!is.null(found)) {
        return(found)
      }
    }
  }
  NULL
}

gp_collect_named_matrices <- function(x, default_name = "matrix", depth = 0L, max_depth = 3L) {
  mat <- gp_as_numeric_square_matrix(x)
  if (!is.null(mat)) {
    return(stats::setNames(list(mat), default_name))
  }
  if (!is.list(x) || depth > max_depth) {
    return(list())
  }
  out <- list()
  x_names <- names(x)
  for (i in seq_along(x)) {
    nm <- if (!is.null(x_names) && nzchar(x_names[[i]] %||% "")) {
      x_names[[i]]
    } else {
      paste0(default_name, "_", i)
    }
    mat <- gp_as_numeric_square_matrix(x[[i]])
    if (!is.null(mat)) {
      out[[nm]] <- mat
    } else if (is.list(x[[i]])) {
      nested <- gp_collect_named_matrices(
        x[[i]],
        default_name = nm,
        depth = depth + 1L,
        max_depth = max_depth
      )
      if (length(nested)) {
        for (nested_nm in names(nested)) {
          out[[nested_nm]] <- nested[[nested_nm]]
        }
      }
    }
  }
  out
}

gp_matrix_aliases <- function() {
  # Phase 3.27: include the env-by-env matrix names actually emitted by
  # gp_framework.py (env_covariance_fa from the FA loading factorisation
  # Sigma_g = Lambda Lambda' + Psi, env_covariance_est from the explicit
  # env-covariance module). Both surface as named square matrices on the
  # Python `result` dict and were previously invisible to the alias scan,
  # so covariance_components came back empty for every GP MET model.
  list(
    Genetic_covariance_traits = c("Genetic_covariance_traits", "genetic_covariance", "trait_covariance_matrix"),
    Genetic_correlation_traits = c("Genetic_correlation_traits", "genetic_correlation", "trait_correlation"),
    Genetic_covariance_environments = c("Genetic_covariance_environments", "environment_covariance", "env_covariance",
                                         "env_covariance_fa", "env_covariance_est"),
    Genetic_correlation_environments = c("Genetic_correlation_environments", "environment_correlation", "env_correlation",
                                           "env_correlation_fa"),
    Genetic_covariance_trait_environment = c("Genetic_covariance_trait_environment", "gxe_trait_covariance", "gxe_covariance"),
    Genetic_correlation_trait_environment = c("Genetic_correlation_trait_environment", "gxe_trait_correlation", "gxe_correlation"),
    Total_genetic_covariance_traits = c("Total_genetic_covariance_traits", "total_genetic_trait_covariance"),
    Total_genetic_correlation_traits = c("Total_genetic_correlation_traits", "total_genetic_trait_correlation"),
    Residual_covariance_traits = c("Residual_covariance_traits", "residual_covariance", "residual_trait_covariance"),
    Residual_correlation_traits = c("Residual_correlation_traits", "residual_correlation", "residual_trait_correlation"),
    Residual_covariance_environments = c("Residual_covariance_environments", "residual_environment_covariance"),
    Residual_correlation_environments = c("Residual_correlation_environments", "residual_environment_correlation"),
    Prediction_covariance_traits = c("Prediction_covariance_traits"),
    Prediction_correlation_traits = c("Prediction_correlation_traits"),
    Prediction_error_covariance_traits = c("Prediction_error_covariance_traits"),
    Prediction_error_correlation_traits = c("Prediction_error_correlation_traits"),
    Prediction_covariance_environments = c("Prediction_covariance_environments"),
    Prediction_correlation_environments = c("Prediction_correlation_environments"),
    Prediction_error_covariance_environments = c("Prediction_error_covariance_environments"),
    Prediction_error_correlation_environments = c("Prediction_error_correlation_environments"),
    Genetic_covariance = c("Genetic_covariance"),
    Genetic_correlation = c("Genetic_correlation")
  )
}

gp_matrix_component_name <- function(public_name, value_type = c("covariance", "correlation")) {
  value_type <- match.arg(value_type)
  nm <- tolower(as.character(public_name %||% ""))
  root <- if (grepl("prediction_error", nm)) {
    "prediction_error"
  } else if (grepl("prediction", nm)) {
    "prediction"
  } else if (grepl("trait_environment|gxe|g:env|interaction", nm)) {
    "gxe"
  } else if (grepl("residual|resid|noise", nm)) {
    "residual"
  } else if (grepl("environment|\\benv\\b", nm)) {
    "environment"
  } else if (grepl("genetic|trait", nm)) {
    "genetic"
  } else {
    "model"
  }
  paste(root, value_type, sep = "_")
}

gp_matrix_to_component_table <- function(mat,
                                         matrix_name,
                                         component,
                                         estimation_method = "model_reported",
                                         env_labels = NULL) {
  mat <- gp_as_numeric_square_matrix(mat)
  if (is.null(mat) || nrow(mat) < 2L) {
    return(gp_empty_covariance_components())
  }
  # Phase 3.27: apply caller-supplied env labels when present and matching
  # the matrix size (rescues numpy-array Cov/Cor that came from Python
  # without dimnames -- env_covariance_fa is the canonical case).
  if (!is.null(env_labels) && length(env_labels) == nrow(mat) &&
      length(env_labels) == ncol(mat)) {
    dimnames(mat) <- list(as.character(env_labels), as.character(env_labels))
  }
  row_labs <- rownames(mat)
  col_labs <- colnames(mat)
  if (is.null(col_labs) || length(col_labs) != ncol(mat) || any(!nzchar(col_labs))) {
    col_labs <- paste0("col", seq_len(ncol(mat)))
  }
  # Phase 3.27: when reticulate converts a pandas DataFrame, it often
  # preserves the column index (env labels) but resets the row index to
  # an integer 0/1/.../n-1 range. The resulting rownames look like "1",
  # "2", ... -- not informative. For square matrices (typical of env-by-env
  # Cov/Cor), inherit colnames onto rows in that case.
  cols_are_real <- !is.null(col_labs) && length(col_labs) == ncol(mat) &&
    all(nzchar(col_labs)) && !grepl("^col[0-9]+$", col_labs[[1L]])
  rows_are_integer_strings <- !is.null(row_labs) && length(row_labs) == nrow(mat) &&
    all(nzchar(row_labs)) && all(grepl("^[0-9]+$", row_labs))
  if (is.null(row_labs) || length(row_labs) != nrow(mat) || any(!nzchar(row_labs)) ||
      (rows_are_integer_strings && cols_are_real && nrow(mat) == ncol(mat))) {
    row_labs <- if (nrow(mat) == ncol(mat) && cols_are_real) {
      col_labs
    } else {
      paste0("row", seq_len(nrow(mat)))
    }
  }
  row_idx <- rep(seq_len(nrow(mat)), times = ncol(mat))
  col_idx <- rep(seq_len(ncol(mat)), each = nrow(mat))
  data.frame(
    Matrix = rep(as.character(matrix_name)[1L], length(row_idx)),
    Component = rep(as.character(component)[1L], length(row_idx)),
    Row = row_labs[row_idx],
    Column = col_labs[col_idx],
    Value = as.numeric(mat[cbind(row_idx, col_idx)]),
    Estimation_method = rep(as.character(estimation_method)[1L], length(row_idx)),
    stringsAsFactors = FALSE
  )
}

gp_format_matrix_component_table <- function(x,
                                             value_type = c("covariance", "correlation"),
                                             default_estimation_method = "model_reported") {
  value_type <- match.arg(value_type)
  if (!is.data.frame(x) || !nrow(x)) {
    return(gp_empty_covariance_components())
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Matrix" %in% names(out)) {
    out[["Matrix"]] <- value_type
  }
  if (!"Component" %in% names(out)) {
    out[["Component"]] <- gp_matrix_component_name(out[["Matrix"]][[1L]], value_type = value_type)
  }
  if (!"Row" %in% names(out)) {
    out[["Row"]] <- NA_character_
  }
  if (!"Column" %in% names(out)) {
    out[["Column"]] <- NA_character_
  }
  if (!"Value" %in% names(out)) {
    value_col <- gp_contract_first_name(names(out), c("value", "estimate", "Estimate", "Components"))
    out[["Value"]] <- if (!is.na(value_col)) gp_contract_as_numeric(out[[value_col]]) else NA_real_
  } else {
    out[["Value"]] <- gp_contract_as_numeric(out[["Value"]])
  }
  method_col <- gp_contract_first_name(
    names(out),
    c("Estimation_method", "estimation_method", "Method", "method")
  )
  out[["Estimation_method"]] <- if (!is.na(method_col)) {
    as.character(out[[method_col]])
  } else {
    default_estimation_method
  }
  out[["Estimation_method"]][is.na(out[["Estimation_method"]]) | !nzchar(out[["Estimation_method"]])] <-
    default_estimation_method
  out <- out[, gp_covariance_component_columns(), drop = FALSE]
  rownames(out) <- NULL
  out
}

gp_extract_matrix_component_table <- function(x,
                                              value_type = c("covariance", "correlation"),
                                              default_estimation_method = "model_reported",
                                              env_labels = NULL) {
  value_type <- match.arg(value_type)
  if (!is.list(x)) {
    return(gp_empty_covariance_components())
  }

  existing_name <- if (identical(value_type, "covariance")) {
    c("covariance_components", "Covariance_components")
  } else {
    c("correlation_components", "Correlation_components")
  }
  existing <- gp_vc_find_named_data_frame(x, existing_name)
  pieces <- list()
  if (is.data.frame(existing) && nrow(existing)) {
    pieces[[length(pieces) + 1L]] <- gp_format_matrix_component_table(
      existing,
      value_type = value_type,
      default_estimation_method = default_estimation_method
    )
  }

  aliases <- gp_matrix_aliases()
  for (public_name in names(aliases)) {
    if (!grepl(value_type, public_name, ignore.case = TRUE)) {
      next
    }
    obj <- gp_find_named_object(x, aliases[[public_name]])
    mats <- gp_collect_named_matrices(obj, default_name = public_name)
    mats <- mats[vapply(mats, function(mat) nrow(mat) >= 2L, logical(1L))]
    if (!length(mats)) {
      next
    }
    component <- gp_matrix_component_name(public_name, value_type = value_type)
    for (matrix_name in names(mats)) {
      pieces[[length(pieces) + 1L]] <- gp_matrix_to_component_table(
        mats[[matrix_name]],
        matrix_name = matrix_name,
        component = component,
        estimation_method = default_estimation_method,
        env_labels = env_labels
      )
    }
  }

  generic_name <- if (identical(value_type, "covariance")) {
    c("Covariance", "covariance")
  } else {
    c("Correlation", "correlation")
  }
  generic <- gp_find_named_object(x, generic_name)
  mats <- gp_collect_named_matrices(generic, default_name = generic_name[[1L]])
  if (length(mats)) {
    for (matrix_name in names(mats)) {
      pieces[[length(pieces) + 1L]] <- gp_matrix_to_component_table(
        mats[[matrix_name]],
        matrix_name = matrix_name,
        component = gp_matrix_component_name(matrix_name, value_type = value_type),
        estimation_method = default_estimation_method,
        env_labels = env_labels
      )
    }
  }

  # Phase 3.27: drill into gp_framework's `interaction_correlation_summary`
  # dict. Per-term it carries `cov_total` + `cor_total` (env x env genetic
  # covariance / correlation) plus `cov_by_kernel` + `cor_by_kernel` (the
  # per-kernel decomposition). Pull out the value_type-appropriate subset
  # so MET GP cov/corr tables are no longer empty.
  ics <- gp_find_named_object(x, c("interaction_correlation_summary"))
  if (is.list(ics) && length(ics)) {
    sub_keys <- if (identical(value_type, "covariance")) {
      c("cov_total", "cov_by_kernel")
    } else {
      c("cor_total", "cor_by_kernel")
    }
    for (term_name in names(ics)) {
      term <- ics[[term_name]]
      if (!is.list(term)) next
      for (sub_key in sub_keys) {
        sub <- term[[sub_key]]
        if (is.null(sub)) next
        if (is.data.frame(sub) || (is.numeric(sub) && is.matrix(sub))) {
          mat_label <- paste(term_name, sub_key, sep = ":")
          pieces[[length(pieces) + 1L]] <- gp_matrix_to_component_table(
            sub,
            matrix_name = mat_label,
            component = gp_matrix_component_name(
              "Genetic_covariance_environments",
              value_type = value_type
            ),
            estimation_method = default_estimation_method,
            env_labels = env_labels
          )
        } else if (is.list(sub)) {
          for (k in names(sub)) {
            kmat <- sub[[k]]
            if (is.null(kmat)) next
            pieces[[length(pieces) + 1L]] <- gp_matrix_to_component_table(
              kmat,
              matrix_name = paste(term_name, sub_key, k, sep = ":"),
              component = gp_matrix_component_name(
                "Genetic_covariance_environments",
                value_type = value_type
              ),
              estimation_method = default_estimation_method,
              env_labels = env_labels
            )
          }
        }
      }
    }
  }

  pieces <- pieces[vapply(pieces, is.data.frame, logical(1L))]
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (!length(pieces)) {
    return(gp_empty_covariance_components())
  }
  out <- do.call(rbind, pieces)
  out <- out[, gp_covariance_component_columns(), drop = FALSE]
  rownames(out) <- NULL
  unique(out)
}

gp_extract_correlation_matrices <- function(x) {
  if (!is.list(x)) {
    return(list())
  }
  out <- list()
  aliases <- gp_matrix_aliases()
  for (public_name in names(aliases)) {
    obj <- gp_find_named_object(x, aliases[[public_name]])
    mats <- gp_collect_named_matrices(obj, default_name = public_name)
    mats <- mats[vapply(mats, function(mat) nrow(mat) >= 2L, logical(1L))]
    if (!length(mats)) {
      next
    }
    out[[public_name]] <- mats[[1L]]
  }
  out
}
