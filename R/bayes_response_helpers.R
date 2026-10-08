gp_bayes_attach_interval_provenance <- function(prediction_table,
                                                model_type,
                                                confidence_level = 0.95,
                                                posterior_draw_count = NA_integer_) {
  if (is.null(prediction_table) || !is.data.frame(prediction_table) ||
      !nrow(prediction_table)) {
    return(prediction_table)
  }

  model_type <- toupper(trimws(as.character(model_type %||% ""))[1L])
  draw_based <- !identical(model_type, "RKHS")
  if (draw_based) {
    prediction_table$PEV_basis <- "BGLR posterior target-draw variance"
    prediction_table$Prediction_uncertainty_source <-
      "BGLR posterior target draws"
    prediction_table$Prediction_interval_method <-
      "BGLR posterior target-draw quantiles"
    prediction_table$Prediction_interval_calibration_n <-
      rep(as.integer(posterior_draw_count)[1L], nrow(prediction_table))
  } else {
    prediction_table$PEV_basis <-
      "BGLR RKHS analytic posterior prediction variance"
    prediction_table$Prediction_uncertainty_source <-
      "BGLR RKHS analytic posterior prediction covariance"
    prediction_table$Prediction_interval_method <-
      "Gaussian approximation from BGLR RKHS posterior covariance"
    prediction_table$Prediction_interval_calibration_n <-
      rep(NA_integer_, nrow(prediction_table))
  }
  prediction_table$Prediction_interval_nominal_coverage <-
    rep(as.numeric(confidence_level)[1L], nrow(prediction_table))
  prediction_table
}

gp_bayes_response_type <- function(response_family = "gaussian", y = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  switch(
    fam,
    gaussian = "gaussian",
    binary = "ordinal",
    ordinal = "ordinal",
    stop(
      "BGLR in PredictProR currently supports gaussian, binary, and ordinal response families only. ",
      "Unordered multiclass/nominal traits are not supported by BGLR.",
      call. = FALSE
    )
  )
}

# BGLR::BGLR() runs as.vector(y) and then factor(y, ordered = TRUE), so class
# labels are ordered ALPHABETICALLY regardless of their factor levels (e.g.
# low/medium/high becomes high < low < medium). Encode classes as zero-padded
# codes in the declared order, which sort correctly as strings, and keep the
# labels in the "class_levels" attribute for gp_bayes_bglr_fit() to decode.
gp_bayes_prepare_response <- function(y, response_family = "gaussian",
                                      class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  if (identical(fam, "gaussian")) {
    return(as.numeric(y))
  }
  class_levels <- as.character(class_levels %||% gp_response_class_levels(y, fam))
  labels <- as.character(y)
  unknown <- setdiff(stats::na.omit(labels), class_levels)
  if (length(unknown)) {
    stop("Bayesian classification response contains labels outside class_levels: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  index <- match(labels, class_levels)
  codes <- ifelse(
    is.na(index), NA_character_,
    formatC(index, width = nchar(length(class_levels)), flag = "0")
  )
  attr(codes, "class_levels") <- class_levels
  codes
}

gp_bayes_decode_class_fit <- function(fit, class_levels) {
  decode <- function(v) {
    if (is.null(v)) return(v)
    idx <- suppressWarnings(as.integer(as.character(v)))
    out <- class_levels[idx]
    out[is.na(v)] <- NA_character_
    out
  }
  if (!is.null(fit$probs)) colnames(fit$probs) <- decode(colnames(fit$probs))
  if (!is.null(fit$SD.probs)) colnames(fit$SD.probs) <- decode(colnames(fit$SD.probs))
  if (!is.null(fit$levels)) fit$levels <- decode(fit$levels)
  if (!is.null(fit$y)) fit$y <- decode(fit$y)
  fit
}

gp_bayes_bglr_fit <- function(fam = "gaussian", ...) {
  args <- list(...)
  # A precomputed RKHS eigen-decomposition (V, d) is of the unweighted kernel;
  # BGLR weights only a kernel it decomposes itself, so drop V/d then.
  w <- args[["weights"]]
  if (!is.null(w) && any(as.double(w) != 1)) {
    args$ETA <- gp_bayes_drop_precomputed_eigen(args$ETA)
  }
  class_levels <- attr(args[["y"]], "class_levels", exact = TRUE)
  run_bglr <- function(bglr_args) {
    withCallingHandlers(
      do.call(BGLR::BGLR, bglr_args),
      warning = function(w) {
        # BGLR ordinal/binary fits can emit this from internal threshold
        # probability calculations on small or sparse class splits while still
        # returning a valid posterior probability matrix.
        if (!identical(fam, "gaussian") && identical(conditionMessage(w), "NaNs produced")) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  fit <- if (!identical(fam, "gaussian") && anyNA(args[["y"]])) {
    gp_bayes_bglr_fit_observed_rows(args, run_bglr)
  } else {
    run_bglr(args)
  }
  if (!identical(fam, "gaussian") && length(class_levels)) {
    fit <- gp_bayes_decode_class_fit(fit, class_levels)
  }
  fit
}

gp_bayes_drop_precomputed_eigen <- function(eta) {
  lapply(eta, function(term) {
    if (is.list(term) && !is.null(term[["K"]])) {
      term[["V"]] <- NULL
      term[["d"]] <- NULL
    }
    term
  })
}

# BGLR (1.1.0 through at least 1.1.4) initialises ordinal/binary thresholds as
# qnorm(cumsum(table(z)) / n) where n also counts NA rows. The top class then
# keeps a finite upper bound qnorm(n_observed / n) that is never updated, the
# latent scale is compressed and NA rows inflate the top class. Classification
# fits therefore run BGLR on the observed rows only and predict the NA rows here:
#   X terms (markers, GBLUP_BRR eigen designs): X_new %*% b over saved effect
#     draws (posterior mean b for FIXED terms, which BGLR does not save);
#   K terms (RKHS): kernel conditional mean K[new, obs] K[obs, obs]^+ u_hat,
#     with conditional variance varU * (K[new, new] - K[new, obs] K^+ K[obs, new])
#     added to the unit probit residual.
# Class probabilities average the probit class masses over the saved draws of
# the thresholds and variances, matching how BGLR forms its own `probs`.
gp_bayes_bglr_fit_observed_rows <- function(args, run_bglr) {
  y <- args[["y"]]
  obs <- !is.na(y)
  new <- which(!obs)
  eta <- args[["ETA"]]
  obs_args <- args
  obs_args$y <- as.vector(y)[obs]
  obs_args$ETA <- lapply(eta, function(term) {
    if (!is.null(term[["K"]])) {
      term[["K"]] <- as.matrix(term[["K"]])[obs, obs, drop = FALSE]
      term[["V"]] <- NULL  # eigenpairs of the full kernel do not fit the subset
      term[["d"]] <- NULL
    } else if (!is.null(term[["X"]])) {
      term[["X"]] <- as.matrix(term[["X"]])[obs, , drop = FALSE]
      if (!identical(term[["model"]], "FIXED")) term[["saveEffects"]] <- TRUE
    }
    term
  })
  fit <- run_bglr(obs_args)

  prefix <- obs_args$saveAt %||% ""
  prefix_base <- sub("x$", "", basename(paste0(prefix, "x")))
  files <- list.files(dirname(paste0(prefix, "x")), full.names = TRUE)
  files <- files[startsWith(basename(files), prefix_base)]
  draw_index <- gp_bayes_saved_draw_indices(
    obs_args$nIter %||% 1500L, obs_args$burnIn %||% 500L, obs_args$thin %||% 5L
  )
  thr_path <- paste0(prefix, "thresholds.dat")
  thresholds <- as.matrix(utils::read.table(thr_path, header = FALSE))
  if (length(draw_index) && max(draw_index) <= nrow(thresholds)) {
    thresholds <- thresholds[draw_index, , drop = FALSE]
  }

  n_new <- length(new)
  mean_parts <- list()
  var_parts <- list()
  for (j in seq_along(eta)) {
    term <- eta[[j]]
    if (!is.null(term[["K"]])) {
      K <- as.matrix(term[["K"]])
      K_oo <- K[obs, obs, drop = FALSE]
      K_no <- K[new, obs, drop = FALSE]
      e <- eigen(K_oo, symmetric = TRUE)
      keep <- e$values > max(e$values) * 1e-10
      K_pinv <- e$vectors[, keep, drop = FALSE] %*%
        (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
      A <- K_no %*% K_pinv
      u_hat <- fit$ETA[[j]]$u %||% rep(0, sum(obs))
      mean_parts[[length(mean_parts) + 1L]] <- matrix(as.vector(A %*% u_hat), ncol = 1L)
      schur <- pmax(diag(K)[new] - rowSums(A * K_no), 0)
      var_u <- gp_bayes_read_eta_variance_draws(files, j, "RKHS", draw_index)
      if (!length(var_u)) var_u <- fit$ETA[[j]]$varU %||% 0
      var_parts[[length(var_parts) + 1L]] <- list(schur = schur, var = as.double(var_u))
    } else if (!is.null(term[["X"]])) {
      X_new <- as.matrix(term[["X"]])[new, , drop = FALSE]
      draws <- if (!identical(term[["model"]], "FIXED")) {
        gp_bayes_read_eta_effect_draws(files, j, X_new, draw_index)
      }
      if (is.null(draws)) {
        b_hat <- fit$ETA[[j]]$b %||% rep(0, ncol(X_new))
        draws <- matrix(as.vector(X_new %*% b_hat), ncol = 1L)
      }
      mean_parts[[length(mean_parts) + 1L]] <- draws
    }
  }

  n_draw <- nrow(thresholds)
  for (part in mean_parts) if (ncol(part) > 1L) n_draw <- min(n_draw, ncol(part))
  for (part in var_parts) if (length(part$var) > 1L) n_draw <- min(n_draw, length(part$var))
  tail_cols <- function(m) if (ncol(m) > 1L) m[, utils::tail(seq_len(ncol(m)), n_draw), drop = FALSE] else m
  mean_parts <- lapply(mean_parts, tail_cols)
  thresholds <- thresholds[utils::tail(seq_len(nrow(thresholds)), n_draw), , drop = FALSE]
  intercept <- fit$mu %||% 0

  n_class <- ncol(thresholds) + 1L
  prob_new <- matrix(0, n_new, n_class)
  latent_sum <- numeric(n_new)
  for (s in seq_len(n_draw)) {
    mu_s <- rep(intercept, n_new)
    for (part in mean_parts) mu_s <- mu_s + part[, min(s, ncol(part))]
    var_s <- rep(1, n_new)
    for (part in var_parts) {
      v <- if (length(part$var) > 1L) utils::tail(part$var, n_draw)[s] else part$var
      var_s <- var_s + v * part$schur
    }
    cuts <- c(-Inf, thresholds[s, ], Inf)
    cdf <- vapply(cuts, function(t) stats::pnorm((t - mu_s) / sqrt(var_s)), numeric(n_new))
    cdf <- matrix(cdf, nrow = n_new)
    prob_new <- prob_new + (cdf[, -1L, drop = FALSE] - cdf[, -ncol(cdf), drop = FALSE])
    latent_sum <- latent_sum + mu_s
  }
  prob_new <- prob_new / n_draw

  n <- length(y)
  probs <- matrix(NA_real_, n, n_class, dimnames = list(NULL, colnames(fit$probs)))
  probs[obs, ] <- fit$probs
  probs[new, ] <- prob_new
  sd_probs <- probs
  sd_probs[] <- NA_real_
  if (!is.null(fit$SD.probs)) sd_probs[obs, ] <- fit$SD.probs
  y_hat <- rep(NA_real_, n)
  y_hat[obs] <- fit$yHat
  y_hat[new] <- latent_sum / n_draw

  fit$probs <- probs
  fit$SD.probs <- sd_probs
  fit$yHat <- y_hat
  fit$y <- as.vector(y)
  fit$whichNa <- new
  fit$heldout_prediction <- list(
    method = "observed_rows_fit_with_conditional_prediction",
    reason = "BGLR ordinal thresholds are biased when the response contains NA rows",
    n_observed = sum(obs),
    n_predicted = n_new,
    n_draws = n_draw
  )
  fit
}

gp_bayes_eta_supports_groups <- function(ETA) {
  if (is.null(ETA) || !length(ETA)) {
    return(FALSE)
  }
  eta_models <- vapply(ETA, function(x) as.character(x[["model"]] %||% ""), character(1))
  # Phase 3.19: include "RKHS" so kernel-prior fits (used by GBLUP_BRR and
  # RKHS on MET) can also produce per-env residual groups via
  # BGLR::BGLR(groups=...). Without this, RKHS was silently sent down the
  # NULL-bayes_groups path which made the per-fold CV produce all-NA
  # predictions in the multi-trait case.
  all(eta_models %in% c("FIXED", "BRR", "BayesB", "BayesC", "BRR_sparse", "RKHS"))
}

gp_bayes_residual_group_count <- function(pheno_data, heter_groups = NULL) {
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) {
    return(0L)
  }
  groups <- as.character(pheno_data[[heter_groups]])
  groups <- groups[!is.na(groups)]
  length(unique(groups))
}

gp_bayes_has_multiple_residual_groups <- function(pheno_data, heter_groups = NULL) {
  gp_bayes_residual_group_count(pheno_data, heter_groups) > 1L
}

gp_bayes_residual_groups <- function(pheno_data,
                                     heter_groups = NULL,
                                     heter_resid = FALSE,
                                     ETA = NULL,
                                     response_family = "gaussian") {
  if (!isTRUE(heter_resid)) {
    return(NULL)
  }
  fam <- gp_resolve_response_family(response_family, y = NULL)
  if (!identical(fam, "gaussian")) {
    stop("BGLR heterogeneous residual groups are currently supported for gaussian responses only.", call. = FALSE)
  }
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) {
    stop("heter_resid=TRUE for Bayesian BGLR requires a valid heter_groups column.", call. = FALSE)
  }
  if (!gp_bayes_has_multiple_residual_groups(pheno_data, heter_groups)) {
    return(NULL)
  }
  if (!gp_bayes_eta_supports_groups(ETA)) {
    return(NULL)
  }
  as.factor(pheno_data[[heter_groups]])
}

gp_bayes_observation_residual_variance <- function(varE,
                                                   groups = NULL,
                                                   weights = NULL,
                                                   n = NULL) {
  varE <- as.double(varE)
  varE_names <- names(varE)
  if (is.null(n)) {
    n <- max(length(groups %||% numeric()), length(weights %||% numeric()), length(varE))
  }
  if (!length(varE)) {
    stop("BGLR residual variance is empty.", call. = FALSE)
  }

  if (length(varE) == 1L) {
    out <- rep(varE, n)
  } else {
    if (is.null(groups)) {
      stop("Grouped BGLR residual variances require observation groups for row-wise expansion.", call. = FALSE)
    }
    group_chr <- as.character(groups)
    if (length(group_chr) != n) {
      stop("BGLR residual group vector length does not match the number of observations.", call. = FALSE)
    }
    if (!is.null(varE_names) && length(varE_names) == length(varE) &&
        all(group_chr %in% varE_names)) {
      out <- varE[match(group_chr, varE_names)]
    } else {
      group_levels <- unique(group_chr)
      if (length(group_levels) != length(varE)) {
        stop(
          "Unable to align grouped BGLR residual variances to observation groups.",
          call. = FALSE
        )
      }
      out <- varE[match(group_chr, group_levels)]
    }
  }

  if (is.null(weights)) {
    weights <- rep(1, n)
  } else if (length(weights) == 1L) {
    weights <- rep(as.double(weights), n)
  } else if (length(weights) != n) {
    stop("BGLR weights length does not match the number of observations.", call. = FALSE)
  }
  weights <- as.double(weights)
  weights[!is.finite(weights) | weights == 0] <- 1
  out / (weights ^ 2)
}

gp_bayes_saved_draw_indices <- function(nIter, burnIn, thin) {
  nIter <- suppressWarnings(as.integer(nIter)[1L])
  burnIn <- suppressWarnings(as.integer(burnIn)[1L])
  thin <- suppressWarnings(as.integer(thin)[1L])
  if (!is.finite(nIter) || !is.finite(burnIn) || !is.finite(thin) ||
      nIter <= 0L || burnIn < 0L || thin <= 0L || burnIn >= nIter) {
    return(integer())
  }
  first <- floor(burnIn / thin) + 1L
  last <- floor(nIter / thin)
  if (first > last) integer() else seq.int(first, last)
}

gp_bayes_read_varE_draws <- function(output_files_names, draw_indices = NULL) {
  varE_file <- output_files_names[grepl("varE\\.dat$", basename(output_files_names))]
  if (length(varE_file) != 1L) {
    stop("Unable to locate the BGLR residual variance draw file.", call. = FALSE)
  }
  tbl <- tryCatch(
    utils::read.table(varE_file, header = FALSE),
    error = function(e) NULL
  )
  if (!is.null(tbl) && nrow(tbl) > 0 && ncol(tbl) > 0) {
    mat <- as.matrix(tbl)
    storage.mode(mat) <- "double"
    if (length(draw_indices) && all(is.finite(draw_indices)) &&
        min(draw_indices) >= 1L && max(draw_indices) <= nrow(mat)) {
      mat <- mat[draw_indices, , drop = FALSE]
    }
    if (ncol(mat) == 1L) {
      return(as.double(mat[, 1]))
    }
    out <- rowMeans(mat, na.rm = TRUE)
    attr(out, "group_draws") <- mat
    return(out)
  }
  out <- scan(varE_file, what = numeric(), sep = "\n", quiet = TRUE)
  if (length(draw_indices) && all(is.finite(draw_indices)) &&
      min(draw_indices) >= 1L && max(draw_indices) <= length(out)) {
    out <- out[draw_indices]
  }
  out
}

gp_bayes_read_mu_draws <- function(output_files_names,
                                   draw_indices = NULL,
                                   n_draws = NULL) {
  mu_file <- output_files_names[grepl("mu\\.dat$", basename(output_files_names))]
  if (length(mu_file) != 1L) {
    return(NULL)
  }
  tbl <- tryCatch(
    as.matrix(utils::read.table(mu_file, header = FALSE)),
    error = function(e) NULL
  )
  if (is.null(tbl) || !nrow(tbl) || !ncol(tbl)) {
    return(NULL)
  }
  storage.mode(tbl) <- "double"
  out <- as.double(tbl[, 1L])
  if (length(draw_indices) && max(draw_indices) <= length(out)) {
    out <- out[draw_indices]
  }
  if (!is.null(n_draws) && is.finite(n_draws) && n_draws > 0L) {
    n_draws <- as.integer(n_draws)
    if (length(out) < n_draws) {
      return(NULL)
    }
    out <- utils::tail(out, n_draws)
  }
  out
}

gp_bayes_read_eta_variance_draws <- function(output_files_names,
                                              eta_index,
                                              eta_model,
                                              draw_indices = NULL) {
  eta_model <- as.character(eta_model %||% "")[[1L]]
  suffix <- if (identical(eta_model, "RKHS")) "varU\\.dat$" else "varB\\.dat$"
  path <- output_files_names[grepl(
    paste0("ETA_", as.integer(eta_index), "_", suffix),
    basename(output_files_names)
  )]
  if (length(path) != 1L) {
    return(NULL)
  }
  values <- tryCatch(
    as.double(utils::read.table(path, header = FALSE)[[1L]]),
    error = function(e) NULL
  )
  if (is.null(values)) {
    return(NULL)
  }
  if (length(draw_indices) && all(is.finite(draw_indices)) &&
      min(draw_indices) >= 1L && max(draw_indices) <= length(values)) {
    values <- values[draw_indices]
  }
  values[is.finite(values)]
}

gp_bayes_read_eta_effect_draws <- function(output_files_names,
                                            eta_index,
                                            X,
                                            draw_indices = NULL) {
  path <- output_files_names[grepl(
    paste0("ETA_", as.integer(eta_index), "_b\\.bin$"),
    basename(output_files_names)
  )]
  if (length(path) != 1L || !requireNamespace("BGLR", quietly = TRUE)) {
    return(NULL)
  }
  beta <- tryCatch(as.matrix(BGLR::readBinMat(path[[1L]])), error = function(e) NULL)
  if (is.null(beta) || !nrow(beta) || !ncol(beta)) {
    return(NULL)
  }
  X <- as.matrix(X)
  if (ncol(beta) != ncol(X) && nrow(beta) == ncol(X)) {
    beta <- t(beta)
  }
  if (ncol(beta) != ncol(X)) {
    return(NULL)
  }
  if (length(draw_indices) && max(draw_indices) <= nrow(beta)) {
    beta <- beta[draw_indices, , drop = FALSE]
  }
  as.matrix(X %*% t(beta))
}

gp_bayes_eta_source_labels <- function(ETA, random_index = NULL) {
  eta_input <- ETA[["ETA"]] %||% ETA
  if (is.null(random_index)) {
    eta_models <- vapply(
      eta_input,
      function(term) as.character(term[["model"]] %||% ""),
      character(1L)
    )
    random_index <- which(!eta_models %in% c("", "FIXED"))
  }
  labels <- ETA[["ETA_element_name"]] %||% character()
  if (length(labels) != length(random_index)) {
    eta_names <- names(eta_input)[random_index]
    if (is.null(eta_names) || length(eta_names) != length(random_index) ||
        any(!nzchar(eta_names))) {
      eta_names <- paste0("eta_", random_index)
    }
    labels <- eta_names
  }
  make.unique(gp_sanitize_kernel_name(as.character(labels)), sep = "_")
}

gp_bayes_classification_eta_role <- function(basis, environments = NULL) {
  if (is.null(environments) || length(unique(as.character(environments))) <= 1L) {
    return("genetic")
  }
  basis <- as.matrix(basis)
  environments <- as.character(environments)
  if (nrow(basis) != length(environments) || ncol(basis) != length(environments)) {
    return("genetic")
  }
  cross_environment <- outer(environments, environments, FUN = "!=")
  within_environment <- !cross_environment
  diag(within_environment) <- FALSE
  cross_mean <- if (any(cross_environment)) {
    mean(abs(basis[cross_environment]), na.rm = TRUE)
  } else {
    0
  }
  within_mean <- if (any(within_environment)) {
    mean(abs(basis[within_environment]), na.rm = TRUE)
  } else {
    mean(abs(diag(basis)), na.rm = TRUE)
  }
  reference <- max(c(within_mean, mean(abs(diag(basis)), na.rm = TRUE), 1),
                   na.rm = TRUE)
  if (is.finite(cross_mean) && cross_mean <= 1e-8 * reference) {
    "gxe"
  } else {
    "genetic_main"
  }
}

gp_bayes_classification_liability_variance <- function(mod,
                                                        ETA,
                                                        pheno_data,
                                                        gen_name,
                                                        heter_groups = NULL,
                                                        bayes_para = NULL,
                                                        confidence_level = 0.95) {
  eta_input <- ETA[["ETA"]] %||% ETA
  if (!is.list(eta_input) || !length(eta_input)) {
    stop("Liability-scale variance extraction requires fitted Bayesian ETA terms.",
         call. = FALSE)
  }
  eta_models <- vapply(
    eta_input,
    function(term) as.character(term[["model"]] %||% ""),
    character(1L)
  )
  supported_models <- c("RKHS", "BRR", "BayesA", "BayesB", "BayesC", "BL")
  random_index <- which(eta_models %in% supported_models)
  if (!length(random_index)) {
    stop("No supported Bayesian genetic ETA term was available for liability-scale variance extraction.",
         call. = FALSE)
  }
  source_labels <- gp_bayes_eta_source_labels(ETA, random_index)
  n_iter <- bayes_para[["nIter"]] %||% NA_integer_
  burn_in <- bayes_para[["burnIn"]] %||% NA_integer_
  thin <- bayes_para[["thin"]] %||% NA_integer_
  draw_index <- gp_bayes_saved_draw_indices(n_iter, burn_in, thin)
  environments <- if (!is.null(heter_groups) &&
                      length(heter_groups) == 1L &&
                      heter_groups %in% names(pheno_data)) {
    as.character(pheno_data[[heter_groups]])
  } else {
    NULL
  }

  observed_index <- !is.na(mod[["model"]][["y"]])
  if (length(observed_index) != nrow(pheno_data)) {
    observed_index <- rep(TRUE, nrow(pheno_data))
  }
  pieces <- lapply(seq_along(random_index), function(piece_index) {
    i <- random_index[[piece_index]]
    term <- eta_input[[i]]
    marker_effect_model <- eta_models[[i]] %in% c("BayesA", "BayesB", "BayesC", "BL")
    effect_draws <- NULL
    if (isTRUE(marker_effect_model)) {
      if (is.null(term[["X"]])) {
        stop("Bayes Alphabet liability-scale variance requires an X matrix for ETA ",
             i, ".", call. = FALSE)
      }
      effect_draws <- gp_bayes_read_eta_effect_draws(
        output_files_names = mod[["output_files_names"]],
        eta_index = i,
        X = term[["X"]],
        draw_indices = draw_index
      )
      if (is.null(effect_draws) || ncol(effect_draws) < 2L) {
        stop("Unable to read posterior marker-effect draws for Bayesian ETA ", i, ".",
             call. = FALSE)
      }
      draws <- apply(
        effect_draws[observed_index, , drop = FALSE],
        2L,
        stats::var,
        na.rm = TRUE
      )
      diagonal <- rep(1, nrow(pheno_data))
      basis <- tcrossprod(as.matrix(term[["X"]]))
    } else {
      basis <- if (!is.null(term[["K"]])) {
        as.matrix(term[["K"]])
      } else if (!is.null(term[["X"]])) {
        tcrossprod(as.matrix(term[["X"]]))
      } else {
        NULL
      }
      if (is.null(basis) || nrow(basis) != nrow(pheno_data) ||
          ncol(basis) != nrow(pheno_data)) {
        stop("Unable to construct the liability-scale covariance basis for Bayesian ETA ",
             i, ".", call. = FALSE)
      }
      draws <- gp_bayes_read_eta_variance_draws(
        output_files_names = mod[["output_files_names"]],
        eta_index = i,
        eta_model = eta_models[[i]],
        draw_indices = draw_index
      )
      if (is.null(draws) || !length(draws)) {
        stop("Unable to read posterior variance draws for Bayesian ETA ", i, ".",
             call. = FALSE)
      }
      diagonal <- as.double(diag(basis))
      if (length(diagonal) != nrow(pheno_data) || any(!is.finite(diagonal))) {
        stop("The covariance diagonal for Bayesian ETA ", i,
             " is unavailable or non-finite.", call. = FALSE)
      }
    }
    list(
      eta_index = i,
      eta_model = eta_models[[i]],
      source_label = source_labels[[piece_index]],
      role = gp_bayes_classification_eta_role(basis, environments),
      draws = pmax(0, as.double(draws)),
      diagonal = pmax(0, diagonal),
      effect_draws = effect_draws
    )
  })
  n_draw <- min(vapply(pieces, function(x) length(x$draws), integer(1L)))
  if (!is.finite(n_draw) || n_draw < 2L) {
    stop("At least two retained posterior variance draws are required for liability-scale summaries.",
         call. = FALSE)
  }
  pieces <- lapply(pieces, function(x) {
    x$draws <- utils::tail(x$draws, n_draw)
    x
  })

  confidence_level <- suppressWarnings(as.double(confidence_level)[1L])
  if (!is.finite(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    confidence_level <- 0.95
  }
  alpha <- (1 - confidence_level) / 2
  vc_rows <- list()
  interval_rows <- list()
  add_summary <- function(component,
                          draws,
                          env = NA_character_,
                          method = "BGLR_posterior_variance_draws",
                          residual_estimated = NA) {
    draws <- as.double(draws)
    draws <- draws[is.finite(draws)]
    if (!length(draws)) {
      stop("No finite posterior draws were available for ", component, ".",
           call. = FALSE)
    }
    posterior_mean <- mean(draws)
    posterior_sd <- if (length(draws) > 1L) stats::sd(draws) else 0
    limits <- stats::quantile(
      draws, probs = c(alpha, 1 - alpha), names = FALSE, na.rm = TRUE,
      type = 8
    )
    vc_rows[[length(vc_rows) + 1L]] <<- data.frame(
      Env = env,
      Component = component,
      Components = posterior_mean,
      Standard_error = posterior_sd,
      stringsAsFactors = FALSE
    )
    interval_rows[[length(interval_rows) + 1L]] <<- data.frame(
      Env = env,
      Component = component,
      Posterior_mean = posterior_mean,
      Posterior_SD = posterior_sd,
      Lower_credible_limit = limits[[1L]],
      Upper_credible_limit = limits[[2L]],
      Credible_level = confidence_level,
      Scale = "latent_liability",
      Link_function = "probit",
      Estimation_method = method,
      Residual_variance_estimated = residual_estimated,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  scaled_draws <- function(piece, index = seq_len(nrow(pheno_data))) {
    scale <- mean(piece$diagonal[index], na.rm = TRUE)
    if (!is.finite(scale) || scale < 0) {
      stop("A non-finite liability-scale kernel diagonal was encountered.",
           call. = FALSE)
    }
    piece$draws * scale
  }
  overall_piece_draws <- lapply(pieces, scaled_draws)
  if (length(pieces) > 1L) {
    for (j in seq_along(pieces)) {
      role <- pieces[[j]]$role
      label <- paste0(
        "liability_scale_", role, "_variance_", pieces[[j]]$source_label
      )
      add_summary(label, overall_piece_draws[[j]])
    }
  }
  all_effect_draws <- all(vapply(
    pieces,
    function(piece) !is.null(piece$effect_draws),
    logical(1L)
  ))
  overall_genetic <- if (isTRUE(all_effect_draws)) {
    n_effect_draw <- min(vapply(pieces, function(piece) {
      ncol(piece$effect_draws)
    }, integer(1L)))
    total_effect <- Reduce(`+`, lapply(pieces, function(piece) {
      piece$effect_draws[, seq_len(n_effect_draw), drop = FALSE]
    }))
    apply(
      total_effect[observed_index, , drop = FALSE],
      2L,
      stats::var,
      na.rm = TRUE
    )
  } else {
    Reduce(`+`, overall_piece_draws)
  }
  overall_genetic_label <- if (length(pieces) == 1L) {
    "liability_scale_genetic_variance"
  } else {
    "liability_scale_total_genetic_variance"
  }
  add_summary(overall_genetic_label, overall_genetic)
  residual_draws <- rep(1, n_draw)
  add_summary(
    "liability_scale_residual_variance_fixed",
    residual_draws,
    method = "probit_identification_constraint",
    residual_estimated = FALSE
  )
  add_summary(
    "liability_scale_heritability",
    overall_genetic / (overall_genetic + residual_draws)
  )

  if (!is.null(environments) && length(unique(environments)) > 1L) {
    for (env in unique(environments)) {
      target_index <- which(environments == env)
      env_piece_draws <- lapply(pieces, scaled_draws, index = target_index)
      roles <- unique(vapply(pieces, `[[`, character(1L), "role"))
      for (role in roles) {
        role_index <- which(vapply(pieces, `[[`, character(1L), "role") == role)
        role_draws <- Reduce(`+`, env_piece_draws[role_index])
        add_summary(
          paste0("liability_scale_", role, "_variance"),
          role_draws,
          env = env
        )
      }
      env_genetic <- Reduce(`+`, env_piece_draws)
      add_summary("liability_scale_total_genetic_variance", env_genetic, env = env)
      add_summary(
        "liability_scale_residual_variance_fixed",
        residual_draws,
        env = env,
        method = "probit_identification_constraint",
        residual_estimated = FALSE
      )
      add_summary(
        "liability_scale_heritability",
        env_genetic / (env_genetic + residual_draws),
        env = env
      )
    }
  }

  variance_components <- do.call(rbind, vc_rows)
  variance_intervals <- do.call(rbind, interval_rows)
  if (all(is.na(variance_components$Env))) {
    variance_components$Env <- NULL
    variance_intervals$Env <- NULL
  }
  rownames(variance_components) <- NULL
  rownames(variance_intervals) <- NULL
  list(
    variance_components = variance_components,
    variance_component_intervals = variance_intervals
  )
}

gp_bayes_prob_matrix <- function(mod, class_levels = NULL) {
  prob <- mod[["model"]][["probs"]] %||% mod[["probs"]]
  if (is.null(prob)) {
    stop("BGLR classification fit did not return posterior class probabilities.", call. = FALSE)
  }
  prob <- as.matrix(prob)
  fit_levels <- mod[["model"]][["levels"]] %||% mod[["levels"]] %||% colnames(prob)
  if (!is.null(fit_levels) && length(fit_levels) == ncol(prob)) {
    colnames(prob) <- as.character(fit_levels)
  }
  if (!is.null(class_levels) && length(class_levels)) {
    class_levels <- as.character(class_levels)
    fit_levels <- colnames(prob)
    if (is.null(fit_levels) && length(class_levels) == ncol(prob)) {
      colnames(prob) <- class_levels
    } else if (!is.null(fit_levels) && setequal(fit_levels, class_levels)) {
      prob <- prob[, class_levels, drop = FALSE]
    } else {
      stop(
        "BGLR class-probability levels do not match the original response levels.",
        call. = FALSE
      )
    }
  }
  prob
}

gp_bayes_classification_output <- function(mod,
                                            pheno_data,
                                            gen_name,
                                            response = NULL,
                                            response_family = "binary",
                                            heter_groups = NULL,
                                            high_reliability_thres = 0.9,
                                            low_reliability_thres = 0.5,
                                            ETA = NULL,
                                            bayes_para = NULL,
                                            GS_model = NULL,
                                            confidence_level = 0.95) {
  fam <- gp_resolve_response_family(response_family, y = mod[["model"]][["y"]])
  original_class_levels <- NULL
  if (!is.null(response) && length(response) == 1L &&
      response %in% names(pheno_data)) {
    original_class_levels <- gp_response_class_levels(
      pheno_data[[response]],
      response_family = fam
    )
  }
  prob <- gp_bayes_prob_matrix(mod, class_levels = original_class_levels)
  class_levels <- colnames(prob) %||% as.character(seq_len(ncol(prob)))
  prob_cols <- gp_prob_column_names(class_levels)
  prob_df <- stats::setNames(as.data.frame(prob, stringsAsFactors = FALSE), prob_cols)
  pred_summary <- gp_classification_prediction_summary(
    prob = prob,
    class_levels = class_levels,
    high_confidence = high_reliability_thres,
    low_confidence = low_reliability_thres
  )

  train_test_label <- ifelse(is.na(mod[["model"]][["y"]]), "Test", "Train")
  predicted_value <- data.frame(
    name = as.character(pheno_data[[gen_name]]),
    Predicted_value = pred_summary$Predicted_class,
    Train_Test_Label = train_test_label,
    pred_summary,
    prob_df,
    stringsAsFactors = FALSE
  )
  predicted_value$Predicted_class <- NULL
  names(predicted_value)[names(predicted_value) == "name"] <- gen_name

  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    predicted_value[[heter_groups]] <- pheno_data[[heter_groups]]
    predicted_value <- predicted_value[, c(
      gen_name, "Predicted_value", "Train_Test_Label", heter_groups,
      setdiff(names(predicted_value), c(gen_name, "Predicted_value", "Train_Test_Label", heter_groups))
    ), drop = FALSE]
  }

  observed_value <- mod[["model"]][["y"]]
  residual_value <- data.frame(
    name = as.character(pheno_data[[gen_name]]),
    Observed_value = ifelse(is.na(observed_value), NA_character_, as.character(observed_value)),
    Predicted_value = predicted_value$Predicted_value,
    Train_Test_Label = train_test_label,
    stringsAsFactors = FALSE
  )
  names(residual_value)[names(residual_value) == "name"] <- gen_name
  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    residual_value[[heter_groups]] <- pheno_data[[heter_groups]]
  }

  total_predicted_value <- NULL
  observed_rows <- !is.na(observed_value)
  if (length(observed_rows) != nrow(pheno_data)) {
    observed_rows <- rep(TRUE, nrow(pheno_data))
  }
  has_met <- anyDuplicated(as.character(pheno_data[[gen_name]][observed_rows])) > 0
  if (has_met) {
    prob_by_gid <- stats::aggregate(prob_df, by = list(GID = as.character(pheno_data[[gen_name]])), FUN = mean, na.rm = TRUE)
    prob_by_gid_mat <- as.matrix(prob_by_gid[, -1, drop = FALSE])
    colnames(prob_by_gid_mat) <- class_levels
    across_summary <- gp_classification_prediction_summary(
      prob = prob_by_gid_mat,
      class_levels = class_levels,
      high_confidence = high_reliability_thres,
      low_confidence = low_reliability_thres
    )
    gid_test <- unique(as.character(pheno_data[[gen_name]][is.na(observed_value)]))
    total_predicted_value <- data.frame(
      name = prob_by_gid[[1]],
      Predicted_value = across_summary$Predicted_class,
      Train_Test_Label = ifelse(prob_by_gid[[1]] %in% gid_test, "Test", "Train"),
      across_summary,
      stats::setNames(as.data.frame(prob_by_gid_mat, stringsAsFactors = FALSE), prob_cols),
      stringsAsFactors = FALSE
    )
    total_predicted_value$Predicted_class <- NULL
    names(total_predicted_value)[names(total_predicted_value) == "name"] <- gen_name
  }

  liability_variance <- if (!is.null(ETA) &&
                            as.character(GS_model %||% "") %in%
                              c("RKHS", "BRR", "GBLUP_BRR", "BayesA", "BayesB", "BayesC", "BL")) {
    gp_bayes_classification_liability_variance(
      mod = mod,
      ETA = ETA,
      pheno_data = pheno_data,
      gen_name = gen_name,
      heter_groups = heter_groups,
      bayes_para = bayes_para,
      confidence_level = confidence_level
    )
  } else {
    list(variance_components = NULL, variance_component_intervals = NULL)
  }

  eta_input <- ETA[["ETA"]] %||% ETA
  eta_models <- if (is.list(eta_input) && length(eta_input)) {
    vapply(eta_input, function(term) as.character(term[["model"]] %||% ""), character(1L))
  } else {
    character()
  }
  random_index <- which(!eta_models %in% c("", "FIXED"))
  source_labels <- if (length(random_index)) {
    gp_bayes_eta_source_labels(ETA, random_index)
  } else {
    character()
  }
  model_ready <- if (length(random_index)) {
    out <- lapply(random_index, function(i) {
      eta_input[[i]][["X"]] %||% eta_input[[i]][["K"]]
    })
    names(out) <- paste0(source_labels, "_model_ready")
    out
  } else {
    list()
  }
  model_parameters <- data.frame(
    stat = c(
      "bayesian_model",
      "bayesian_feature_block_count",
      "bayesian_feature_block_names",
      "bayesian_feature_block_strategy",
      "classification_variance_scale"
    ),
    summary = c(
      as.character(GS_model %||% NA_character_),
      as.character(length(source_labels)),
      if (length(source_labels)) paste(source_labels, collapse = ";") else NA_character_,
      "separate_BGLR_ETA_components_summed_for_prediction",
      "latent_liability_probit_residual_variance_fixed_to_1"
    ),
    stringsAsFactors = FALSE
  )

  list(
    Coefficients = list(),
    Estimated_breeding_value = list(),
    Total_estimated_breeding_value = NULL,
    Predicted_value = predicted_value,
    Total_Predicted_value = total_predicted_value,
    Residual_value = residual_value,
    Variance_components = liability_variance$variance_components,
    variance_component_intervals = liability_variance$variance_component_intervals,
    M_matrix_model_ready = model_ready,
    model_parameters = model_parameters,
    diagnostic_tst_plot = NULL
  )
}
