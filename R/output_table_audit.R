# Content checks for every output table of a model_execute() result.
#
# The scenario suite used to check only that prediction tables existed and
# that NA variances were flagged. On 2026-10-19 that suite showed 317/317 PASS
# while joint Bayesian multi-trait predictions sat at twice the phenotype
# level and MET GBLUP_BRR predictions left out the environment effects.
# gp_audit_output_tables() checks what the tables say: prediction identities
# and level, variance-component ranges and identities, and the covariance /
# correlation matrices. inst/examples/scenarios/00_setup.R fails a scenario on
# any violation.

gp_audit_numeric <- function(x) suppressWarnings(as.numeric(x))

gp_audit_close <- function(a, b, tol = 1e-6) {
  abs(a - b) <= tol * pmax(1, abs(a), abs(b))
}

# Model results inside a model_execute() return value: every list holding a
# predicted_values data frame, with the path to it.
gp_audit_collect_results <- function(x, path = "", depth = 0L, acc = list()) {
  # plain lists only (plots and fitted model objects are not walked)
  walkable <- is.list(x) && !is.data.frame(x) &&
    (identical(class(x), "list") || inherits(x, "predictpror_model_result"))
  if (depth > 6L || !walkable) return(acc)
  if (is.data.frame(x[["predicted_values"]])) {
    acc[[if (nzchar(path)) path else "/"]] <- x
  }
  skip <- c("cv_results_raw", "cv_results_predicted_vs_observed", "raw_python_result",
            "gp_result", "gp_info", "bglr_multitrait_model", "bayes_model")
  nms <- names(x) %||% rep("", length(x))
  for (i in seq_along(x)) {
    nm <- if (nzchar(nms[i])) nms[i] else as.character(i)
    if (nm %in% skip) next
    el <- tryCatch(x[[i]], error = function(e) NULL)
    acc <- gp_audit_collect_results(el, paste0(path, "/", nm), depth + 1L, acc)
  }
  acc
}

gp_audit_is_bayesian <- function(result) {
  mp <- result[["model_parameters"]]
  is.data.frame(mp) && any(grepl("^bayes|BGLR", c(as.character(mp[[1L]]), as.character(mp[[2L]] %||% "")),
                                 ignore.case = TRUE))
}

#' Check the contents of a model_execute() result's output tables
#'
#' Internal. Returns one row per check: `path` (where the model result sits),
#' `table`, `check`, `ok` (`TRUE`, `FALSE`, or `NA` when informational) and
#' `detail`.
#'
#' @param x A `model_execute()` result (or one model result).
#' @param level_tol_overall Largest allowed |mean(observed - predicted)| on
#'   training rows, per trait, in phenotypic SD units.
#' @param level_tol_group The same per trait x environment.
#' @return A data frame of checks.
#' @keywords internal
#' @noRd
gp_audit_output_tables <- function(x, level_tol_overall = 0.75, level_tol_group = 1.0) {
  rows <- list()
  add <- function(path, table, check, ok, detail = "") {
    ok <- if (length(ok)) as.logical(ok[[1L]]) else NA
    detail <- if (length(detail)) substr(paste(detail, collapse = " "), 1L, 300L) else ""
    rows[[length(rows) + 1L]] <<- data.frame(path = path, table = table, check = check, ok = ok,
                                             detail = detail, stringsAsFactors = FALSE)
  }
  results <- gp_audit_collect_results(x)
  seen <- list()
  for (path in names(results)) {
    res <- results[[path]]
    # by_trait / by_model hold the same result as model_results: check once
    if (any(vapply(seen, identical, logical(1), res))) next
    seen[[length(seen) + 1L]] <- res
    bayes <- gp_audit_is_bayesian(res)
    gp_audit_predictions(res, path, add, level_tol_overall, level_tol_group)
    vc <- res[["variance_components"]]
    if (is.data.frame(vc) && nrow(vc)) gp_audit_variance_components(vc, path, add, bayes)
    cov_m <- gp_audit_matrix_tables(res[["covariance_components"]], path, "covariance_components", add)
    cor_m <- gp_audit_matrix_tables(res[["correlation_components"]], path, "correlation_components", add)
    gp_audit_cov_cor_consistency(cov_m, cor_m, path, add)
    if (is.data.frame(vc) && nrow(vc)) gp_audit_cov_vs_vc(cov_m, vc, path, add)
    gp_audit_other_tables(res, path, add)
    gp_audit_met_environments(res, path, add)
  }
  if (!length(rows)) {
    return(data.frame(path = character(), table = character(), check = character(),
                      ok = logical(), detail = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

gp_audit_predictions <- function(res, path, add, level_tol_overall, level_tol_group) {
  pv <- res[["predicted_values"]]
  t <- "predicted_values"
  report <- tryCatch(
    gp_output_contract_validate_one(res, path = path),
    error = function(e) data.frame(valid = FALSE, issues = conditionMessage(e))
  )
  add(path, t, "prediction output contract", isTRUE(all(report$valid)), report$issues)
  if (!nrow(pv)) return(invisible())
  has <- function(cols) all(cols %in% names(pv))
  if ("Predicted_value" %in% names(pv) && is.numeric(pv$Predicted_value)) {
    add(path, t, "predicted values finite", all(is.finite(pv$Predicted_value)),
        sprintf("%d non-finite of %d", sum(!is.finite(pv$Predicted_value)), nrow(pv)))
  }
  if ("Predicted_class" %in% names(pv)) {
    pr <- as.matrix(pv[, grep("^Probability_", names(pv)), drop = FALSE])
    storage.mode(pr) <- "double"
    if (ncol(pr)) {
      s <- rowSums(pr)
      add(path, t, "class probabilities in [0, 1] summing to 1",
          !anyNA(pr) && all(pr >= -1e-10 & pr <= 1 + 1e-10) && all(abs(s - 1) < 1e-6),
          sprintf("max |sum - 1| %.2g", max(abs(s - 1), na.rm = TRUE)))
      cls <- sub("^Probability_", "", colnames(pr))
      top <- apply(pr, 1L, max)
      chosen <- pr[cbind(seq_len(nrow(pr)), match(as.character(pv$Predicted_class), cls))]
      add(path, t, "predicted class has the highest probability",
          all(!is.na(chosen) & chosen >= top - 1e-8), sprintf("%d rows", sum(is.na(chosen) | chosen < top - 1e-8)))
    }
  }
  if (has(c("PEV", "Standard_error"))) {
    pev <- gp_audit_numeric(pv$PEV); se <- gp_audit_numeric(pv$Standard_error)
    k <- is.finite(pev) & is.finite(se)
    add(path, t, "PEV = SE^2", all(gp_audit_close(pev[k], se[k]^2)), sprintf("%d rows", sum(k)))
    add(path, t, "PEV and SE missing together", identical(is.finite(pev), is.finite(se)), "")
    add(path, t, "SE >= 0", all(se[k] >= 0), "")
  }
  if (has(c("lower_bound", "upper_bound", "Predicted_value"))) {
    lo <- gp_audit_numeric(pv$lower_bound); hi <- gp_audit_numeric(pv$upper_bound)
    pred <- gp_audit_numeric(pv$Predicted_value)
    k <- is.finite(lo) & is.finite(hi)
    add(path, t, "lower bound <= upper bound", all(lo[k] <= hi[k] + 1e-10), "")
    method <- as.character(pv$Prediction_interval_method %||% rep(NA_character_, nrow(pv)))
    se <- gp_audit_numeric(pv$Standard_error %||% rep(NA_real_, nrow(pv)))
    g <- k & grepl("^Gaussian approximation", method) & is.finite(se)
    if (any(g)) {
      cov <- gp_audit_numeric(pv$Prediction_interval_nominal_coverage %||% 0.95)[g]
      cov[!is.finite(cov)] <- 0.95
      z <- stats::qnorm(1 - (1 - cov) / 2)
      ok <- gp_audit_close(lo[g], pred[g] - z * se[g]) & gp_audit_close(hi[g], pred[g] + z * se[g])
      add(path, t, "Gaussian bounds = prediction -/+ z SE", all(ok), sprintf("%d of %d rows", sum(!ok), sum(g)))
    }
    q <- k & !grepl("^Gaussian approximation", method)
    if (any(q)) {
      inside <- lo[q] <= pred[q] + 1e-8 & pred[q] <= hi[q] + 1e-8
      add(path, t, "interval contains the prediction", all(inside), sprintf("%d of %d rows outside", sum(!inside), sum(q)))
    }
  }
  if ("Reliability" %in% names(pv)) {
    r <- gp_audit_numeric(pv$Reliability); kk <- is.finite(r)
    add(path, t, "reliability in [0, 1]", all(r[kk] >= -1e-10 & r[kk] <= 1 + 1e-10), "")
    if (has(c("Reliability_basis", "Reliability_variance_input", "Reliability_reference_variance"))) {
      basis <- as.character(pv$Reliability_basis)
      gen <- kk & (grepl("^genetic_(plus_gxe_)?prior_variance(:|$)", basis))
      if (any(gen)) {
        e <- pmin(1, pmax(0, 1 - gp_audit_numeric(pv$Reliability_variance_input[gen]) /
                            gp_audit_numeric(pv$Reliability_reference_variance[gen])))
        ok <- abs(r[gen] - e) < 1e-8
        add(path, t, "genetic reliability = 1 - genetic PEV / prior", all(ok), sprintf("%d of %d rows", sum(!ok), sum(gen)))
      }
    }
    vc <- res[["variance_components"]]
    if (is.data.frame(vc) && all(c("Component", "Boundary") %in% names(vc))) {
      ne <- which(vc$Boundary %in% "not_estimable" & grepl("^residual_variance", vc$Component))
      for (i in ne) {
        env <- sub("^residual_variance_?", "", vc$Component[i])
        s <- if (nzchar(env) && "Env" %in% names(pv)) pv$Env == env else rep(TRUE, nrow(pv))
        if ("Trait" %in% names(vc) && "Trait" %in% names(pv)) s <- s & pv$Trait == vc$Trait[i]
        add(path, t, "no reliability where the residual is at zero", !any(kk[s]),
            sprintf("%s: %d of %d rows report reliability", if (nzchar(env)) env else "all", sum(kk[s]), sum(s)))
      }
    }
  }
  # model-based (mixed-model / Bayesian / GP) prediction SEs stay on the
  # phenotype's scale; held-out ML/DL risk estimates are not held to this
  if (has(c("Standard_error", "Observed_value"))) {
    source <- as.character(pv$Prediction_uncertainty_source %||% rep("", nrow(pv)))
    model_based <- !grepl("heldout|cross_fitted|local|unavailable", source)
    trait <- if ("Trait" %in% names(pv)) as.character(pv$Trait) else rep("", nrow(pv))
    se <- gp_audit_numeric(pv$Standard_error); obs <- gp_audit_numeric(pv$Observed_value)
    for (tr in unique(trait)) {
      s <- trait == tr
      sdv <- stats::sd(obs[s], na.rm = TRUE)
      k <- s & model_based & is.finite(se)
      if (!is.finite(sdv) || sdv <= 0 || !any(k)) next
      add(path, t, "model-based SE below twice the phenotypic SD", max(se[k]) <= 2 * sdv,
          sprintf("%s largest SE = %.2f phenotypic SD", if (nzchar(tr)) tr else "trait", max(se[k]) / sdv))
    }
  }
  # level: predictions on the phenotype scale (training rows, per trait / env)
  if (has(c("Observed_value", "Predicted_value")) && is.numeric(pv$Predicted_value)) {
    obs <- gp_audit_numeric(pv$Observed_value); pred <- gp_audit_numeric(pv$Predicted_value)
    train <- if ("Train_Test_Label" %in% names(pv)) pv$Train_Test_Label == "Train" else rep(TRUE, nrow(pv))
    train <- train & is.finite(obs) & is.finite(pred)
    trait <- if ("Trait" %in% names(pv)) as.character(pv$Trait) else rep("", nrow(pv))
    env <- if ("Env" %in% names(pv)) as.character(pv$Env) else NULL
    for (tr in unique(trait[train])) {
      s <- train & trait == tr
      if (sum(s) < 5L) next
      sdv <- stats::sd(obs[s])
      if (!is.finite(sdv) || sdv <= 0) next
      bias <- mean(obs[s] - pred[s]) / sdv
      add(path, t, "prediction level matches observed (training rows)", abs(bias) <= level_tol_overall,
          sprintf("%s mean(observed - predicted) = %.2f phenotypic SD", if (nzchar(tr)) tr else "trait", bias))
      if (!is.null(env)) {
        by_env <- tapply(obs[s] - pred[s], env[s], mean) / sdv
        n_env <- tapply(obs[s], env[s], length)
        by_env <- by_env[n_env >= 5L]
        if (length(by_env)) {
          worst <- names(by_env)[which.max(abs(by_env))]
          add(path, t, "prediction level matches observed in each environment",
              all(abs(by_env) <= level_tol_group),
              sprintf("%s %s: %.2f phenotypic SD", if (nzchar(tr)) tr else "trait", worst, by_env[[worst]]))
        }
      }
    }
  }
  invisible()
}

gp_audit_variance_components <- function(vc, path, add, bayes = FALSE) {
  t <- "variance_components"
  if (!all(c("Component", "Components") %in% names(vc))) {
    add(path, t, "variance table has Component and Components columns", FALSE, paste(names(vc), collapse = ","))
    return(invisible())
  }
  comp <- as.character(vc$Component); v <- gp_audit_numeric(vc$Components)
  se <- gp_audit_numeric(vc$Standard_error %||% rep(NA_real_, nrow(vc)))
  nid <- grepl("not_identifiable", comp)
  if (any(nid)) add(path, t, "not-identifiable rows are NA", all(is.na(v[nid])), paste(comp[nid & !is.na(v)], collapse = ","))
  variance <- grepl("variance|!var$|!R$", comp) &
    !grepl("covariance|not_identifiable|lambda|prediction_confidence|classification_uncertainty", comp)
  bad <- variance & is.finite(v) & v < -1e-10
  add(path, t, "variances >= 0", !any(bad), paste(sprintf("%s=%.3g", comp[bad], v[bad]), collapse = ", "))
  h <- grepl("heritab", comp) & !nid
  badh <- h & is.finite(v) & (v < -1e-10 | v > 1 + 1e-10)
  if (any(h)) add(path, t, "heritability in [0, 1]", !any(badh), paste(sprintf("%s=%.3g", comp[badh], v[badh]), collapse = ", "))
  add(path, t, "standard errors >= 0", !any(is.finite(se) & se < 0), "")
  group_cols <- intersect(c("Trait", "Env", "Hybrid_component", "Kernel"), names(vc))
  key <- if (length(group_cols)) do.call(paste, c(vc[group_cols], sep = "|")) else rep("", nrow(vc))
  dup <- duplicated(paste(key, comp))
  add(path, t, "no duplicated component", !any(dup), paste(unique(paste(key, comp)[dup]), collapse = ","))
  # h2 identity (Bayesian h2 is the posterior mean of the per-draw ratio)
  if (!bayes) {
    for (g in unique(key)) {
      s <- key == g; cc <- comp[s]; vv <- stats::setNames(v[s], cc)
      for (hn in cc[grepl("^heritability", cc) & !grepl("not_identifiable|from_mean", cc)]) {
        suf <- sub("^heritability", "", hn)
        gn <- paste0("genetic_variance", suf); en <- paste0("residual_variance", suf)
        if (!(gn %in% cc && en %in% cc)) {
          if (identical(suf, "") && all(c("total_genetic_variance", "residual_variance") %in% cc)) {
            gn <- "total_genetic_variance"
          } else {
            next
          }
        }
        gx <- paste0("gxe_variance", suf)
        gg <- vv[[gn]] + if (gx %in% cc && is.finite(vv[[gx]])) vv[[gx]] else 0
        if (!all(is.finite(c(gg, vv[[en]], vv[[hn]])))) next
        e <- gg / (gg + vv[[en]])
        add(path, t, "heritability = genetic / (genetic + residual)", gp_audit_close(vv[[hn]], e, 1e-4),
            sprintf("%s %s = %.6g, ratio %.6g", g, hn, vv[[hn]], e))
      }
    }
  }
  invisible()
}

gp_audit_matrix_tables <- function(tbl, path, table, add) {
  out <- list()
  if (!is.data.frame(tbl) || !nrow(tbl) || !all(c("Row", "Column", "Value") %in% names(tbl))) return(out)
  kind <- if (grepl("correlation", table)) "correlation" else "covariance"
  key_cols <- intersect(c("Matrix", "Component"), names(tbl))
  key <- if (length(key_cols)) do.call(paste, c(tbl[key_cols], sep = "|")) else rep("matrix", nrow(tbl))
  for (g in unique(key)) {
    d <- tbl[key == g, , drop = FALSE]
    r <- as.character(d$Row); cl <- as.character(d$Column)
    labels <- unique(c(r, cl))
    add(path, table, "real row/column labels", !any(labels %in% c("not_identifiable", "", NA)), g)
    add(path, table, "no duplicated cells", !anyDuplicated(paste(r, cl)), g)
    M <- matrix(NA_real_, length(labels), length(labels), dimnames = list(labels, labels))
    M[cbind(r, cl)] <- gp_audit_numeric(d$Value)
    # a correlation is undefined (NA) where a variance is zero
    add(path, table, "matrix complete", !anyNA(M) || kind == "correlation",
        sprintf("%s: %d missing of %d", g, sum(is.na(M)), length(M)))
    if (anyNA(M)) {
      # a correlation left NA by a zero variance still pairs with its covariance
      if (kind == "correlation") out[[g]] <- M
      next
    }
    add(path, table, "symmetric", max(abs(M - t(M))) <= 1e-8 * max(1, max(abs(M))), g)
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    add(path, table, "positive semi-definite", min(ev) >= -1e-8 * max(1, max(abs(ev))),
        sprintf("%s: smallest eigenvalue %.3g", g, min(ev)))
    if (kind == "covariance") {
      add(path, table, "variances on the diagonal >= 0", all(diag(M) >= -1e-10), g)
    } else {
      add(path, table, "correlation diagonal = 1", all(abs(diag(M) - 1) < 1e-8), g)
      add(path, table, "correlations in [-1, 1]", all(abs(M) <= 1 + 1e-8), g)
    }
    out[[g]] <- M
  }
  out
}

gp_audit_cov_cor_consistency <- function(cov_m, cor_m, path, add) {
  # every covariance matrix has its correlation matrix
  cor_for <- function(g) gsub(":cov_", ":cor_", gsub("covariance", "correlation", g))
  for (g in names(cov_m)) {
    add(path, "correlation_components", "every covariance matrix has a correlation matrix",
        cor_for(g) %in% names(cor_m), g)
  }
  for (g in names(cor_m)) {
    gc <- gsub(":cor_", ":cov_", gsub("correlation", "covariance", g))
    C <- cov_m[[gc]]
    if (is.null(C)) next
    R <- cor_m[[g]]
    if (!setequal(rownames(C), rownames(R))) {
      add(path, "correlation_components", "same labels as the covariance", FALSE, g)
      next
    }
    R <- R[rownames(C), colnames(C)]
    d <- sqrt(pmax(diag(C), 0)); zero <- d == 0
    E <- C / outer(d, d)
    k <- !outer(zero, zero, "|")
    k <- k & is.finite(R)
    add(path, "correlation_components", "correlation = cov2cor(covariance)",
        !any(k) || max(abs(R[k] - E[k])) < 1e-6, g)
  }
}

gp_audit_cov_vs_vc <- function(cov_m, vc, path, add) {
  if (!all(c("Component", "Components") %in% names(vc))) return(invisible())
  comp <- as.character(vc$Component); v <- gp_audit_numeric(vc$Components)
  lookup <- function(level, type, by) {
    nm <- paste0(type, "_variance")
    if (by %in% names(vc)) {
      s <- as.character(vc[[by]]) == level & comp == nm
      if (sum(s) == 1L) return(v[s])
    }
    s <- comp == paste0(nm, "_", level)
    if (sum(s) == 1L) v[s] else NA_real_
  }
  for (g in names(cov_m)) {
    mtx <- sub("\\|.*", "", g)
    type <- if (grepl("^Genetic_covariance_(traits|environments)$", mtx)) "genetic" else
      if (grepl("^Residual_covariance_(traits|environments)$", mtx)) "residual" else next
    by <- if (grepl("traits", mtx)) "Trait" else "Env"
    M <- cov_m[[g]]
    ref <- vapply(rownames(M), lookup, numeric(1), type = type, by = by)
    k <- is.finite(ref)
    if (!any(k)) next
    add(path, "covariance_components", "diagonal = variance-table variances",
        all(gp_audit_close(diag(M)[k], ref[k], 1e-4)), g)
  }
  invisible()
}

gp_audit_other_tables <- function(res, path, add) {
  kvc <- res[["kernel_variance_components"]]
  if (is.data.frame(kvc) && all(c("Component", "Components") %in% names(kvc))) {
    comp <- as.character(kvc$Component); v <- gp_audit_numeric(kvc$Components)
    bad <- grepl("variance", comp) & !grepl("covariance", comp) & is.finite(v) & v < -1e-10
    add(path, "kernel_variance_components", "variances >= 0", !any(bad), paste(comp[bad], collapse = ","))
    if ("Kernel_weight" %in% names(kvc)) {
      w <- gp_audit_numeric(kvc$Kernel_weight)
      add(path, "kernel_variance_components", "kernel weights in [0, 1]",
          all(is.na(w) | (w >= -1e-10 & w <= 1 + 1e-10)), "")
    }
  }
  cp <- res[["correlation_parameters"]]
  if (is.data.frame(cp) && "estimate" %in% names(cp)) {
    e <- gp_audit_numeric(cp$estimate)
    add(path, "correlation_parameters", "correlations in [-1, 1]", all(is.na(e) | abs(e) <= 1 + 1e-8), "")
  }
  vci <- res[["variance_component_intervals"]]
  if (is.data.frame(vci) && all(c("Lower_credible_limit", "Upper_credible_limit") %in% names(vci))) {
    lo <- gp_audit_numeric(vci$Lower_credible_limit); hi <- gp_audit_numeric(vci$Upper_credible_limit)
    k <- is.finite(lo) & is.finite(hi)
    add(path, "variance_component_intervals", "lower <= upper", all(lo[k] <= hi[k]), "")
    h <- grepl("heritab", vci$Component)
    add(path, "variance_component_intervals", "heritability limits in [0, 1]",
        all(!h | !k | (lo >= -1e-10 & hi <= 1 + 1e-10)), "")
  }
  rv <- res[["Residual_variance_by_environment"]]
  if (is.data.frame(rv) && "Residual_variance" %in% names(rv)) {
    v <- gp_audit_numeric(rv$Residual_variance)
    add(path, "Residual_variance_by_environment", "residual variances >= 0", all(is.na(v) | v >= 0), "")
  }
  invisible()
}

# MET results: per-environment variance rows name real environments, and a
# mixed model with per-environment genetic variances reports the environment
# genetic correlation.
gp_audit_met_environments <- function(res, path, add) {
  pv <- res[["predicted_values"]]; vc <- res[["variance_components"]]
  if (!is.data.frame(pv) || !"Env" %in% names(pv) || !is.data.frame(vc) || !"Component" %in% names(vc)) {
    return(invisible())
  }
  envs <- unique(as.character(pv$Env))
  if (length(envs) < 2L || "Env" %in% names(vc)) return(invisible())
  comp <- as.character(vc$Component)
  per_env <- grepl("^(residual_variance|heritability)_.+$", comp) &
    !grepl("not_identifiable|from_mean_variances|across_environments", comp)
  if (!any(per_env)) return(invisible())
  suffix <- sub("^(residual_variance|heritability)_", "", comp[per_env])
  add(path, "variance_components", "per-environment rows name prediction environments",
      all(suffix %in% envs), paste(setdiff(unique(suffix), envs), collapse = ","))
  genetic_env <- sum(paste0("genetic_variance_", envs) %in% comp)
  if (genetic_env >= 2L) {
    cor_tbl <- res[["correlation_components"]]
    mats <- if (is.data.frame(cor_tbl) && "Matrix" %in% names(cor_tbl)) unique(as.character(cor_tbl$Matrix)) else character()
    add(path, "correlation_components", "environment genetic correlation reported",
        any(grepl("correlation_environments|G:Env:cor", mats) & !grepl("^Prediction|^Residual", mats)),
        paste(mats, collapse = ","))
  }
  invisible()
}
