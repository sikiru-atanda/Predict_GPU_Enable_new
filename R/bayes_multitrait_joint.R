bayes_multitrait_cov_prior <- function(type = "UN", n_traits) {
  type <- toupper(as.character(type %||% "UN")[1])
  if (!type %in% c("UN", "DIAG", "FA", "REC")) {
    type <- "UN"
  }
  out <- list(type = type)
  if (type %in% c("UN", "DIAG")) {
    out$S0 <- diag(n_traits)
    out$df0 <- n_traits + 1L
  }
  out
}

bayes_multitrait_kernel_features <- function(K, tol = 1e-8) {
  K <- as.matrix(K)
  K <- (K + t(K)) / 2
  ev <- eigen(K, symmetric = TRUE)
  keep <- which(is.finite(ev$values) & ev$values > tol)
  if (!length(keep)) {
    stop("BGLR multi-trait BRR kernel feature expansion found no positive eigenvalues.", call. = FALSE)
  }
  X <- ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), nrow = length(keep))
  rownames(X) <- rownames(K)
  colnames(X) <- paste0("eig", seq_len(ncol(X)))
  X
}

bayes_multitrait_prepare_kernels <- function(kernels) {
  if (is.null(kernels) || !length(kernels)) {
    stop("BGLR multi-trait Bayesian prediction requires at least one relationship/kernel matrix.", call. = FALSE)
  }
  kernels <- kernels[!vapply(kernels, is.null, logical(1))]
  kernels <- lapply(kernels, function(K) {
    K <- as.matrix(K)
    if (is.null(rownames(K)) || is.null(colnames(K))) {
      stop("All BGLR multi-trait Bayesian kernels must have row and column names keyed by genotype ID.", call. = FALSE)
    }
    ids <- intersect(rownames(K), colnames(K))
    if (!length(ids)) {
      stop("All BGLR multi-trait Bayesian kernels must share row and column IDs.", call. = FALSE)
    }
    K <- K[ids, ids, drop = FALSE]
    storage.mode(K) <- "numeric"
    (K + t(K)) / 2
  })
  common_ids <- Reduce(intersect, lapply(kernels, rownames))
  if (!length(common_ids)) {
    stop("BGLR multi-trait Bayesian kernels do not share genotype IDs.", call. = FALSE)
  }
  kernels <- lapply(kernels, function(K) K[common_ids, common_ids, drop = FALSE])
  names(kernels) <- names(kernels) %||% paste0("kernel", seq_along(kernels))
  names(kernels) <- make.unique(gp_sanitize_kernel_name(names(kernels)), sep = "_")
  kernels
}

bayes_multitrait_expand_pheno <- function(pheno_data,
                                          response,
                                          gen_name,
                                          heter_groups = NULL,
                                          kernel_ids) {
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  ph_ids <- unique(as.character(ph[[gen_name]]))
  extra_ids <- setdiff(kernel_ids, ph_ids)
  if (!length(extra_ids)) {
    return(ph)
  }
  has_env <- !is.null(heter_groups) && heter_groups %in% names(ph) &&
    length(unique(stats::na.omit(as.character(ph[[heter_groups]])))) > 1L
  extra <- if (isTRUE(has_env)) {
    env_levels <- unique(as.character(ph[[heter_groups]]))
    expand.grid(
      .gid = extra_ids,
      .env = env_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(.gid = extra_ids, stringsAsFactors = FALSE)
  }
  add <- as.data.frame(matrix(NA, nrow = nrow(extra), ncol = ncol(ph)), stringsAsFactors = FALSE)
  names(add) <- names(ph)
  add[[gen_name]] <- extra$.gid
  if (isTRUE(has_env)) {
    add[[heter_groups]] <- extra$.env
  }
  for (trait in response) {
    add[[trait]] <- NA_real_
  }
  rbind(ph, add)
}

bayes_multitrait_build_eta <- function(pheno_data,
                                       response,
                                       gen_name,
                                       heter_groups,
                                       kernels,
                                       GS_model,
                                       trait_cov_type = "UN") {
  n_traits <- length(response)
  trait_cov <- bayes_multitrait_cov_prior(trait_cov_type, n_traits)
  ids <- as.character(pheno_data[[gen_name]])
  is_met <- !is.null(heter_groups) && heter_groups %in% names(pheno_data) &&
    length(unique(stats::na.omit(as.character(pheno_data[[heter_groups]])))) > 1L
  ETA <- list()
  designs <- list()
  term_roles <- character()
  term_kernel <- character()

  if (isTRUE(is_met)) {
    ZE_fixed <- stats::model.matrix(
      ~ .env,
      data = data.frame(.env = factor(as.character(pheno_data[[heter_groups]])))
    )
    ZE_fixed <- ZE_fixed[, colnames(ZE_fixed) != "(Intercept)", drop = FALSE]
    ETA$ENV <- list(X = ZE_fixed, model = "FIXED")
    term_roles["ENV"] <- "fixed_environment"
  }

  for (nm in names(kernels)) {
    K <- kernels[[nm]]
    K_base <- K[ids, ids, drop = FALSE]
    rownames(K_base) <- paste0("row", seq_len(nrow(K_base)))
    colnames(K_base) <- rownames(K_base)
    if (identical(GS_model, "RKHS")) {
      ETA[[paste0(nm, "_G")]] <- list(K = K_base, model = "RKHS", Cov = trait_cov)
      term_roles[paste0(nm, "_G")] <- "genetic"
    } else {
      X <- bayes_multitrait_kernel_features(K_base)
      ETA[[paste0(nm, "_G")]] <- list(X = X, model = "BRR", Cov = trait_cov)
      designs[[paste0(nm, "_G")]] <- X
      term_roles[paste0(nm, "_G")] <- "genetic"
    }
    term_kernel[paste0(nm, "_G")] <- nm

    if (isTRUE(is_met)) {
      ZE <- stats::model.matrix(
        ~ .env - 1,
        data = data.frame(.env = factor(as.character(pheno_data[[heter_groups]])))
      )
      K_gxe <- K_base * tcrossprod(ZE)
      rownames(K_gxe) <- colnames(K_gxe) <- rownames(K_base)
      if (identical(GS_model, "RKHS")) {
        ETA[[paste0(nm, "_GxE")]] <- list(K = K_gxe, model = "RKHS", Cov = trait_cov)
      } else {
        X_gxe <- bayes_multitrait_kernel_features(K_gxe)
        ETA[[paste0(nm, "_GxE")]] <- list(X = X_gxe, model = "BRR", Cov = trait_cov)
        designs[[paste0(nm, "_GxE")]] <- X_gxe
      }
      term_roles[paste0(nm, "_GxE")] <- "gxe"
      term_kernel[paste0(nm, "_GxE")] <- nm
    }
  }

  list(
    ETA = ETA,
    designs = designs,
    term_roles = term_roles,
    term_kernel = term_kernel,
    is_met = is_met
  )
}

bayes_multitrait_random_component <- function(fit, ETA, designs, term_roles) {
  n <- nrow(fit$ETAHat)
  t <- ncol(fit$ETAHat)
  random_hat <- matrix(0, nrow = n, ncol = t)
  colnames(random_hat) <- colnames(fit$ETAHat)
  components <- list()
  for (nm in names(ETA)) {
    if (identical(term_roles[[nm]], "fixed_environment")) {
      next
    }
    term_fit <- fit$ETA[[nm]]
    comp <- NULL
    if (!is.null(term_fit$u)) {
      comp <- as.matrix(term_fit$u)
    } else if (!is.null(term_fit$beta) && !is.null(designs[[nm]])) {
      comp <- as.matrix(designs[[nm]]) %*% as.matrix(term_fit$beta)
    }
    if (!is.null(comp) && all(dim(comp) == c(n, t))) {
      random_hat <- random_hat + comp
      components[[nm]] <- comp
    }
  }
  list(random_hat = random_hat, components = components)
}

bayes_multitrait_safe_cov2cor <- function(x, labels = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x <- (x + t(x)) / 2
  if (!is.null(labels) && length(labels) == nrow(x)) {
    dimnames(x) <- list(as.character(labels), as.character(labels))
  }
  scale <- sqrt(pmax(0, diag(x)))
  denom <- outer(scale, scale)
  out <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  ok <- is.finite(denom) & denom > 0
  out[ok] <- x[ok] / denom[ok]
  diag(out)[is.finite(scale) & scale > 0] <- 1
  out[is.finite(out) & out > 1] <- 1
  out[is.finite(out) & out < -1] <- -1
  out
}

bayes_multitrait_eta_scale <- function(eta_term) {
  if (!is.null(eta_term$K)) {
    value <- mean(diag(as.matrix(eta_term$K)), na.rm = TRUE)
  } else if (!is.null(eta_term$X)) {
    value <- mean(rowSums(as.matrix(eta_term$X) ^ 2), na.rm = TRUE)
  } else {
    value <- NA_real_
  }
  if (is.finite(value) && value > 0) value else 1
}

bayes_multitrait_lower_triangle_index <- function(n_traits) {
  which(lower.tri(matrix(0, n_traits, n_traits), diag = TRUE), arr.ind = TRUE)
}

bayes_multitrait_covariance_draws <- function(covariance_object,
                                              n_traits,
                                              draw_indices = NULL) {
  file <- covariance_object$fName_Omega %||% covariance_object$fName_R %||% NULL
  if (is.null(file) || !file.exists(file)) {
    return(NULL)
  }
  out <- tryCatch(
    as.matrix(utils::read.table(file, header = FALSE)),
    error = function(e) NULL
  )
  expected <- n_traits * (n_traits + 1L) / 2L
  if (is.null(out) || !nrow(out) || ncol(out) < expected) {
    return(NULL)
  }
  storage.mode(out) <- "double"
  out <- out[, seq_len(expected), drop = FALSE]
  if (length(draw_indices) && max(draw_indices) <= nrow(out)) {
    out <- out[draw_indices, , drop = FALSE]
  }
  out
}

bayes_multitrait_sum_draws <- function(draws, scales) {
  keep <- !vapply(draws, is.null, logical(1L))
  draws <- draws[keep]
  scales <- scales[keep]
  if (!length(draws)) {
    return(NULL)
  }
  n_draws <- min(vapply(draws, nrow, integer(1L)))
  if (!is.finite(n_draws) || n_draws < 2L) {
    return(NULL)
  }
  Reduce(`+`, Map(function(x, scale) {
    x[seq_len(n_draws), , drop = FALSE] * scale
  }, draws, scales))
}

bayes_multitrait_covariance_bundle <- function(fit,
                                               ETA,
                                               term_roles,
                                               labels,
                                               draw_indices = NULL) {
  n_traits <- length(labels)
  covariance_for_role <- function(role) {
    term_names <- names(term_roles)[term_roles == role]
    matrices <- list()
    draws <- list()
    scales <- numeric()
    for (nm in term_names) {
      omega <- fit$ETA[[nm]]$Cov$Omega %||% NULL
      if (is.null(omega)) {
        next
      }
      omega <- as.matrix(omega)
      if (!all(dim(omega) == c(n_traits, n_traits))) {
        next
      }
      scale <- bayes_multitrait_eta_scale(ETA[[nm]])
      matrices[[nm]] <- scale * omega
      matrices[[nm]] <- (matrices[[nm]] + t(matrices[[nm]])) / 2
      dimnames(matrices[[nm]]) <- list(labels, labels)
      draws[[nm]] <- bayes_multitrait_covariance_draws(
        fit$ETA[[nm]]$Cov,
        n_traits,
        draw_indices = draw_indices
      )
      scales[[nm]] <- scale
    }
    covariance <- if (length(matrices)) Reduce(`+`, matrices) else NULL
    if (!is.null(covariance)) {
      covariance <- (covariance + t(covariance)) / 2
      dimnames(covariance) <- list(labels, labels)
    }
    draws_by_term <- lapply(names(draws), function(nm) {
      one <- draws[[nm]]
      if (is.null(one)) NULL else one * scales[[nm]]
    })
    names(draws_by_term) <- names(draws)
    list(
      covariance = covariance,
      draws = bayes_multitrait_sum_draws(draws, scales),
      covariance_by_term = matrices,
      draws_by_term = draws_by_term
    )
  }

  main <- covariance_for_role("genetic")
  gxe <- covariance_for_role("gxe")
  residual <- as.matrix(fit$resCov$R)
  residual <- (residual + t(residual)) / 2
  dimnames(residual) <- list(labels, labels)
  list(
    main = main,
    gxe = gxe,
    residual = list(
      covariance = residual,
      draws = bayes_multitrait_covariance_draws(
        fit$resCov,
        n_traits,
        draw_indices = draw_indices
      )
    )
  )
}

bayes_multitrait_kernel_covariance_contract <- function(covariance_bundle,
                                                        term_kernel,
                                                        kernel_matrices,
                                                        labels,
                                                        is_met = FALSE,
                                                        axis_label = "Trait") {
  kernel_names <- names(kernel_matrices)
  if (!length(kernel_names)) {
    return(NULL)
  }
  term_names <- names(term_kernel)
  term_kernel <- as.character(term_kernel)
  names(term_kernel) <- term_names
  axis_label <- match.arg(axis_label, c("Trait", "Env"))
  estimation_method <- "BGLR_Multitrait_independent_kernel_covariance_posterior"

  combine_role <- function(role_bundle, kernel_name) {
    term_names <- names(term_kernel)[term_kernel == kernel_name]
    term_names <- intersect(term_names, names(role_bundle$covariance_by_term))
    matrices <- role_bundle$covariance_by_term[term_names]
    draws <- role_bundle$draws_by_term[term_names]
    list(
      covariance = if (length(matrices)) Reduce(`+`, matrices) else NULL,
      draws = if (length(draws)) {
        bayes_multitrait_sum_draws(draws, rep(1, length(draws)))
      } else {
        NULL
      }
    )
  }
  posterior_diag_sd <- function(draws) {
    diagonal <- bayes_multitrait_diagonal_draws(draws, length(labels))
    if (is.null(diagonal) || nrow(diagonal) < 2L) {
      return(rep(NA_real_, length(labels)))
    }
    apply(diagonal, 2L, stats::sd, na.rm = TRUE)
  }
  long_matrix <- function(x, kernel_name, component, diagonal_mean) {
    grid <- expand.grid(
      axis_row = labels,
      axis_column = labels,
      stringsAsFactors = FALSE
    )
    names(grid)[1:2] <- paste0(axis_label, c("_row", "_column"))
    grid$Kernel <- kernel_name
    grid$Component <- component
    grid$Covariance <- as.numeric(x[cbind(
      match(grid[[paste0(axis_label, "_row")]], labels),
      match(grid[[paste0(axis_label, "_column")]], labels)
    )])
    grid$Kernel_weight <- NA_real_
    grid$Kernel_diagonal_mean <- diagonal_mean
    grid$Estimation_method <- estimation_method
    grid$Independent_kernel_estimate <- TRUE
    grid$Fit_converged <- NA
    grid[, c(
      "Kernel", "Component", paste0(axis_label, "_row"),
      paste0(axis_label, "_column"), "Covariance", "Kernel_weight",
      "Kernel_diagonal_mean", "Estimation_method",
      "Independent_kernel_estimate", "Fit_converged"
    ), drop = FALSE]
  }

  main_by_kernel <- list()
  gxe_by_kernel <- list()
  total_by_kernel <- list()
  variance_rows <- list()
  covariance_rows <- list()
  for (kernel_name in kernel_names) {
    main <- combine_role(covariance_bundle$main, kernel_name)
    if (is.null(main$covariance)) {
      next
    }
    gxe <- if (isTRUE(is_met)) {
      combine_role(covariance_bundle$gxe, kernel_name)
    } else {
      list(covariance = NULL, draws = NULL)
    }
    diagonal_mean <- mean(diag(as.matrix(kernel_matrices[[kernel_name]])), na.rm = TRUE)
    main_by_kernel[[kernel_name]] <- main$covariance
    if (!is.null(gxe$covariance)) {
      gxe_by_kernel[[kernel_name]] <- gxe$covariance
    }
    total_by_kernel[[kernel_name]] <- main$covariance +
      if (is.null(gxe$covariance)) 0 else gxe$covariance

    make_variance_rows <- function(piece, component, estimand) {
      axis <- data.frame(label = labels, stringsAsFactors = FALSE)
      names(axis) <- axis_label
      cbind(
        data.frame(Kernel = kernel_name, stringsAsFactors = FALSE),
        axis,
        data.frame(
          Component = component,
          Components = diag(piece$covariance),
          Standard_error = posterior_diag_sd(piece$draws),
          Estimation_method = estimation_method,
          Kernel_weight = NA_real_,
          Kernel_diagonal_mean = diagonal_mean,
          Independent_kernel_estimate = TRUE,
          Variance_estimand = estimand,
          Fit_converged = NA,
          stringsAsFactors = FALSE
        )
      )
    }
    variance_rows[[length(variance_rows) + 1L]] <- make_variance_rows(
      main,
      "kernel_genetic_variance",
      "independently_fitted_kernel_trait_covariance_diagonal"
    )
    covariance_rows[[length(covariance_rows) + 1L]] <- long_matrix(
      main$covariance, kernel_name, "genetic_covariance", diagonal_mean
    )
    if (!is.null(gxe$covariance)) {
      variance_rows[[length(variance_rows) + 1L]] <- make_variance_rows(
        gxe,
        "kernel_gxe_variance",
        "independently_fitted_kernel_gxe_covariance_diagonal"
      )
      covariance_rows[[length(covariance_rows) + 1L]] <- long_matrix(
        gxe$covariance, kernel_name, "gxe_covariance", diagonal_mean
      )
      covariance_rows[[length(covariance_rows) + 1L]] <- long_matrix(
        total_by_kernel[[kernel_name]], kernel_name,
        "total_genetic_covariance", diagonal_mean
      )
    }
  }
  if (!length(total_by_kernel)) {
    return(NULL)
  }
  list(
    kernel_variance_components = do.call(rbind, variance_rows),
    # Each kernel entered BGLR as its own ETA term, so its covariance
    # contribution is fitted independently rather than split from a shared one.
    kernel_configuration = gp_multitrait_kernel_configuration(
      names(total_by_kernel)
    ),
    Genetic_covariance_by_kernel = total_by_kernel,
    genetic_covariance_by_kernel_long = do.call(rbind, covariance_rows),
    kernel_combination = list(
      kernel_names = names(total_by_kernel),
      kernel_weights = stats::setNames(rep(NA_real_, length(total_by_kernel)), names(total_by_kernel)),
      kernel_diagonal_means = stats::setNames(
        vapply(kernel_matrices[names(total_by_kernel)], function(x) {
          mean(diag(as.matrix(x)), na.rm = TRUE)
        }, numeric(1L)),
        names(total_by_kernel)
      ),
      covariance_basis = "sum_of_independently_fitted_BGLR_kernel_covariances",
      independent_kernel_covariances_estimated = TRUE
    )
  )
}

bayes_multitrait_diagonal_draws <- function(draws, n_traits) {
  if (is.null(draws)) {
    return(NULL)
  }
  index <- bayes_multitrait_lower_triangle_index(n_traits)
  diag_cols <- which(index[, 1L] == index[, 2L])
  draws[, diag_cols, drop = FALSE]
}

bayes_multitrait_quadratic_draws <- function(draws, weights) {
  if (is.null(draws)) {
    return(NULL)
  }
  weights <- as.numeric(weights)
  index <- bayes_multitrait_lower_triangle_index(length(weights))
  coefficient <- weights[index[, 1L]] * weights[index[, 2L]]
  coefficient[index[, 1L] != index[, 2L]] <-
    2 * coefficient[index[, 1L] != index[, 2L]]
  drop(draws %*% coefficient)
}

bayes_multitrait_covariance_draw_sd <- function(draws, labels) {
  n_traits <- length(labels)
  out <- matrix(NA_real_, n_traits, n_traits, dimnames = list(labels, labels))
  if (is.null(draws) || nrow(draws) < 2L) {
    return(out)
  }
  index <- bayes_multitrait_lower_triangle_index(n_traits)
  value <- apply(draws, 2, stats::sd, na.rm = TRUE)
  out[index] <- value
  out[cbind(index[, 2L], index[, 1L])] <- value
  out
}

bayes_multitrait_correlation_draw_sd <- function(draws, labels) {
  n_traits <- length(labels)
  out <- matrix(NA_real_, n_traits, n_traits, dimnames = list(labels, labels))
  if (is.null(draws) || nrow(draws) < 2L) {
    return(out)
  }
  index <- bayes_multitrait_lower_triangle_index(n_traits)
  diag_cols <- which(index[, 1L] == index[, 2L])
  diag_by_trait <- match(seq_len(n_traits), index[diag_cols, 1L])
  diag_by_trait <- diag_cols[diag_by_trait]
  for (k in seq_len(nrow(index))) {
    i <- index[k, 1L]
    j <- index[k, 2L]
    if (i == j) {
      value <- 0
    } else {
      denom <- sqrt(draws[, diag_by_trait[[i]]] * draws[, diag_by_trait[[j]]])
      correlation_draw <- draws[, k] / denom
      correlation_draw[!is.finite(correlation_draw)] <- NA_real_
      value <- stats::sd(correlation_draw, na.rm = TRUE)
    }
    out[i, j] <- out[j, i] <- value
  }
  out
}

bayes_multitrait_long_prediction <- function(pred_mat,
                                             se_mat,
                                             random_hat,
                                             y_mat,
                                             pheno_data,
                                             response,
                                             gen_name,
                                             heter_groups,
                                             confidence_level,
                                             CI_width_thresholds,
                                             high_reliability_thres,
                                             low_reliability_thres,
                                             trait_genetic_variance = NULL) {
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  # Honor user-supplied gen_name / heter_groups as the meta column names.
  gid_col_name <- if (!is.null(gen_name) && nzchar(as.character(gen_name))) as.character(gen_name) else "GID"
  meta <- data.frame(
    .gid = as.character(pheno_data[[gen_name]]),
    stringsAsFactors = FALSE
  )
  names(meta)[names(meta) == ".gid"] <- gid_col_name
  include_env <- !is.null(heter_groups) && heter_groups %in% names(pheno_data)
  if (isTRUE(include_env)) {
    meta[[as.character(heter_groups)]] <- as.character(pheno_data[[heter_groups]])
  }
  row_grid <- expand.grid(
    row_index = seq_len(nrow(pred_mat)),
    Trait = response,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx <- cbind(row_grid$row_index, match(row_grid$Trait, response))
  pred <- as.numeric(pred_mat[idx])
  se <- as.numeric(se_mat[idx])
  pev <- se ^ 2
  lower <- pred - z * se
  upper <- pred + z * se
  observed <- as.numeric(y_mat[idx])
  trait_genetic_var <- as.numeric(trait_genetic_variance)
  if (length(trait_genetic_var) != length(response)) {
    trait_genetic_var <- apply(random_hat, 2, stats::var, na.rm = TRUE)
  }
  names(trait_genetic_var) <- response
  trait_genetic_var[!is.finite(trait_genetic_var) | trait_genetic_var <= 0] <- 0
  genetic_var <- as.numeric(trait_genetic_var[match(row_grid$Trait, response)])
  rel <- reliability_thresholds(
    prediction_error_var = pev,
    genetic_var = genetic_var,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  int_stats <- bayes_classify_interval_widths(
    predictions = pred,
    lower_bound = lower,
    upper_bound = upper,
    CI_width_thresholds = CI_width_thresholds
  )

  out <- cbind(meta[row_grid$row_index, , drop = FALSE], data.frame(
    Trait = row_grid$Trait,
    Predicted_value = pred,
    Train_Test_Label = ifelse(is.na(observed), "Test", "Train"),
    Observed_value = observed,
    Standard_error = se,
    PEV = pev,
    lower_bound = lower,
    upper_bound = upper,
    Uncertainty = int_stats$Uncertainty,
    Uncertainty_remarks = int_stats$reliability_remarks,
    Reliability = rel$reliability,
    Reliability_remarks = rel$remarks,
    stringsAsFactors = FALSE
  ))
  gp_format_gaussian_prediction_table(
    out,
    gen_name = gen_name,
    heter_groups = heter_groups,
    include_env = include_env,
    include_trait = TRUE
  )
}

bayes_multitrait_across_environment_prediction <- function(pred_long,
                                                           gen_name = "GID",
                                                           confidence_level = 0.95,
                                                           trait_genetic_variance = NULL,
                                                           high_reliability_thres = 0.9,
                                                           low_reliability_thres = 0.5) {
  if (!is.data.frame(pred_long) || !"Env" %in% names(pred_long)) {
    return(NULL)
  }
  split_key <- paste(pred_long$GID, pred_long$Trait, sep = "\r")
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  pieces <- lapply(split(pred_long, split_key), function(dat) {
    pred <- mean(dat$Predicted_value, na.rm = TRUE)
    n_target <- sum(is.finite(dat$Standard_error))
    se <- if (n_target) sqrt(sum(dat$Standard_error ^ 2, na.rm = TRUE)) / n_target else NA_real_
    obs <- mean(dat$Observed_value, na.rm = TRUE)
    if (!is.finite(obs)) obs <- NA_real_
    pev <- se ^ 2
    genetic_var <- as.numeric(trait_genetic_variance[dat$Trait[[1]]])
    rel <- reliability_thresholds(
      prediction_error_var = pev,
      genetic_var = genetic_var,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )
    data.frame(
      GID = dat$GID[[1]],
      Trait = dat$Trait[[1]],
      Predicted_value = pred,
      Train_Test_Label = ifelse(any(dat$Train_Test_Label == "Test"), "Test", "Train"),
      Observed_value = obs,
      Standard_error = se,
      PEV = pev,
      lower_bound = pred - z * se,
      upper_bound = pred + z * se,
      Uncertainty = 2 * z * se,
      Uncertainty_remarks = ifelse(
        is.finite(se),
        "posterior_marginal_independence_approximation",
        NA_character_
      ),
      Reliability_variance_input = pev,
      Reliability_reference_variance = genetic_var,
      Reliability = rel$reliability,
      Reliability_remarks = rel$remarks,
      Reliability_basis = paste(
        "1 - model-reported prediction-error variance / model-specific",
        "genetic variance; clipped to [0,1]"
      ),
      PEV_basis = paste(
        "model-reported prediction-error variance aggregated across",
        "environments under a posterior-marginal independence approximation"
      ),
      Prediction_uncertainty_source = "posterior_marginal_independence_approximation",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  gp_format_gaussian_prediction_table(out, include_trait = TRUE)
}

bayes_multitrait_variance_components <- function(covariance_bundle, labels, is_met = FALSE) {
  n_traits <- length(labels)
  main <- diag(covariance_bundle$main$covariance)
  residual <- diag(covariance_bundle$residual$covariance)
  gxe <- if (isTRUE(is_met) && !is.null(covariance_bundle$gxe$covariance)) {
    diag(covariance_bundle$gxe$covariance)
  } else {
    rep(0, n_traits)
  }
  main_draws <- bayes_multitrait_diagonal_draws(covariance_bundle$main$draws, n_traits)
  residual_draws <- bayes_multitrait_diagonal_draws(covariance_bundle$residual$draws, n_traits)
  gxe_draws <- bayes_multitrait_diagonal_draws(covariance_bundle$gxe$draws, n_traits)

  posterior_sd <- function(x, trait_index) {
    if (is.null(x) || nrow(x) < 2L) NA_real_ else stats::sd(x[, trait_index], na.rm = TRUE)
  }
  rows <- vector("list", n_traits)
  for (j in seq_len(n_traits)) {
    total_genetic <- main[[j]] + gxe[[j]]
    h2_draws <- NULL
    if (!is.null(main_draws) && !is.null(residual_draws)) {
      n_draws <- min(nrow(main_draws), nrow(residual_draws))
      gd <- main_draws[seq_len(n_draws), j]
      if (isTRUE(is_met) && !is.null(gxe_draws)) {
        n_draws <- min(n_draws, nrow(gxe_draws))
        gd <- main_draws[seq_len(n_draws), j] + gxe_draws[seq_len(n_draws), j]
      } else {
        gd <- gd[seq_len(n_draws)]
      }
      rd <- residual_draws[seq_len(n_draws), j]
      h2_draws <- gd / (gd + rd)
    }
    h2 <- if (!is.null(h2_draws) && any(is.finite(h2_draws))) {
      mean(h2_draws, na.rm = TRUE)
    } else {
      total_genetic / (total_genetic + residual[[j]])
    }
    component <- c("genetic_variance", "residual_variance", "heritability")
    estimate <- c(main[[j]], residual[[j]], h2)
    standard_error <- c(
      posterior_sd(main_draws, j),
      posterior_sd(residual_draws, j),
      if (is.null(h2_draws)) NA_real_ else stats::sd(h2_draws, na.rm = TRUE)
    )
    if (isTRUE(is_met)) {
      component <- c("genetic_variance", "gxe_variance", "residual_variance", "heritability")
      estimate <- c(main[[j]], gxe[[j]], residual[[j]], h2)
      standard_error <- c(
        posterior_sd(main_draws, j),
        posterior_sd(gxe_draws, j),
        posterior_sd(residual_draws, j),
        if (is.null(h2_draws)) NA_real_ else stats::sd(h2_draws, na.rm = TRUE)
      )
    }
    rows[[j]] <- data.frame(
      Trait = labels[[j]],
      Component = component,
      Components = estimate,
      Standard_error = standard_error,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

bayes_multitrait_joint_fit <- function(pheno_data,
                                       response,
                                       gen_name,
                                       heter_groups = NULL,
                                       kernels,
                                       GS_model = c("GBLUP_BRR", "RKHS"),
                                       bayes_para,
                                       trait_cov_type = "UN",
                                       residual_cov_type = "UN",
                                       heter_resid = FALSE,
                                       confidence_level = 0.95,
                                       CI_width_thresholds = c(0.33, 0.66),
                                       high_reliability_thres = 0.9,
                                       low_reliability_thres = 0.5,
                                       verbose = FALSE,
                                       random_state = NULL) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR is required for joint multi-trait Bayesian prediction.", call. = FALSE)
  }
  GS_model <- match.arg(GS_model)
  if (identical(GS_model, "GBLUP_BRR")) {
    GS_model <- "BRR"
  }
  kernels <- bayes_multitrait_prepare_kernels(kernels)
  kernel_ids <- Reduce(intersect, lapply(kernels, rownames))
  ph <- bayes_multitrait_expand_pheno(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_ids = kernel_ids
  )
  if (!all(response %in% names(ph))) {
    stop("All response columns must be present in pheno_data for joint multi-trait Bayesian prediction.", call. = FALSE)
  }
  ph <- ph[as.character(ph[[gen_name]]) %in% kernel_ids, , drop = FALSE]
  y_mat <- as.matrix(ph[, response, drop = FALSE])
  storage.mode(y_mat) <- "numeric"
  colnames(y_mat) <- response
  if (ncol(y_mat) < 2L) {
    stop("Joint multi-trait Bayesian prediction requires at least two response columns.", call. = FALSE)
  }
  eta_bundle <- bayes_multitrait_build_eta(
    pheno_data = ph,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernels = kernels,
    GS_model = GS_model,
    trait_cov_type = trait_cov_type
  )
  if (isTRUE(heter_resid) && isTRUE(eta_bundle$is_met)) {
    message(
      "BGLR joint row-wise multi-trait MET uses resCov for trait residual covariance; ",
      "environment-specific residual variance is not identifiable in this layout and was ignored."
    )
  }
  res_cov <- bayes_multitrait_cov_prior(residual_cov_type, length(response))
  save_prefix <- gp_bglr_save_prefix(
    model_name = paste0("BGLR_Multitrait_", GS_model),
    response = paste(response, collapse = "_")
  )
  cleanup_prefix <- function() {
    unlink(
      list.files(
        path = dirname(save_prefix),
        pattern = paste0("^", basename(save_prefix)),
        full.names = TRUE
      )
    )
    invisible(NULL)
  }
  fit <- NULL
  tryCatch(
    invisible(utils::capture.output({
      # A NULL random_state keeps the caller's stream (CV folds seed it).
      fit <- gp_with_pinned_seed(random_state, BGLR::Multitrait(
        y = y_mat,
        ETA = eta_bundle$ETA,
        resCov = res_cov,
        nIter = bayes_para[["nIter"]],
        burnIn = bayes_para[["burnIn"]],
        thin = bayes_para[["thin"]],
        verbose = isTRUE(verbose),
        saveAt = save_prefix
      ))
    })),
    error = function(e) {
      cleanup_prefix()
      stop(e)
    }
  )
  output_files_names <- list.files(
    path = dirname(save_prefix),
    pattern = paste0("^", basename(save_prefix)),
    full.names = TRUE
  )
  on.exit(unlink(output_files_names), add = TRUE)
  pred_mat <- sweep(as.matrix(fit$ETAHat), 2, as.double(fit$mu), FUN = "+")
  sd_mu <- as.double(fit$SD.mu %||% rep(0, ncol(pred_mat)))
  if (length(sd_mu) != ncol(pred_mat)) {
    sd_mu <- rep(0, ncol(pred_mat))
  }
  se_mat <- sqrt(as.matrix(fit$SD.ETAHat) ^ 2 +
                   matrix(sd_mu ^ 2, nrow = nrow(pred_mat), ncol = ncol(pred_mat), byrow = TRUE))
  random <- bayes_multitrait_random_component(
    fit = fit,
    ETA = eta_bundle$ETA,
    designs = eta_bundle$designs,
    term_roles = eta_bundle$term_roles
  )
  covariance_bundle <- bayes_multitrait_covariance_bundle(
    fit = fit,
    ETA = eta_bundle$ETA,
    term_roles = eta_bundle$term_roles,
    labels = response,
    draw_indices = gp_bayes_saved_draw_indices(
      nIter = bayes_para[["nIter"]],
      burnIn = bayes_para[["burnIn"]],
      thin = bayes_para[["thin"]]
    )
  )
  main_covariance <- covariance_bundle$main$covariance
  gxe_covariance <- covariance_bundle$gxe$covariance
  total_covariance <- main_covariance
  total_covariance_draws <- covariance_bundle$main$draws
  if (!is.null(gxe_covariance)) {
    total_covariance <- total_covariance + gxe_covariance
    total_covariance_draws <- bayes_multitrait_sum_draws(
      list(covariance_bundle$main$draws, covariance_bundle$gxe$draws),
      c(1, 1)
    )
  }
  trait_genetic_variance <- diag(total_covariance)
  names(trait_genetic_variance) <- response
  pred_long <- bayes_multitrait_long_prediction(
    pred_mat = pred_mat,
    se_mat = se_mat,
    random_hat = random$random_hat,
    y_mat = y_mat,
    pheno_data = ph,
    response = response,
    gen_name = gen_name,
    heter_groups = if (isTRUE(eta_bundle$is_met)) heter_groups else NULL,
    confidence_level = confidence_level,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    trait_genetic_variance = trait_genetic_variance
  )
  n_env <- if (isTRUE(eta_bundle$is_met)) {
    length(unique(stats::na.omit(as.character(ph[[heter_groups]]))))
  } else {
    1L
  }
  across_trait_genetic_variance <- diag(main_covariance)
  if (!is.null(gxe_covariance) && n_env > 0L) {
    across_trait_genetic_variance <- across_trait_genetic_variance + diag(gxe_covariance) / n_env
  }
  names(across_trait_genetic_variance) <- response
  total_pred <- bayes_multitrait_across_environment_prediction(
    pred_long = pred_long,
    gen_name = gen_name,
    confidence_level = confidence_level,
    trait_genetic_variance = across_trait_genetic_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  vc <- bayes_multitrait_variance_components(
    covariance_bundle = covariance_bundle,
    labels = response,
    is_met = eta_bundle$is_met
  )
  kernel_contract <- bayes_multitrait_kernel_covariance_contract(
    covariance_bundle = covariance_bundle,
    term_kernel = eta_bundle$term_kernel,
    kernel_matrices = kernels,
    labels = response,
    is_met = eta_bundle$is_met,
    axis_label = "Trait"
  )
  main_correlation <- bayes_multitrait_safe_cov2cor(main_covariance, response)
  residual_correlation <- bayes_multitrait_safe_cov2cor(
    covariance_bundle$residual$covariance,
    response
  )
  total_correlation <- bayes_multitrait_safe_cov2cor(total_covariance, response)
  public_model <- if (identical(GS_model, "BRR")) "GBLUP_BRR" else GS_model
  model_parameters <- data.frame(
    stat = c("mode", "model", "bglr_function", "eta_cov_argument", "residual_cov_type", "trait_cov_type"),
    summary = c(
      if (isTRUE(eta_bundle$is_met)) "multi_trait_multi_environment_bayes" else "multi_trait_bayes",
      public_model,
      "BGLR::Multitrait",
      "Cov",
      toupper(residual_cov_type),
      toupper(trait_cov_type)
    ),
    eta_cov_argument = "Cov",
    stringsAsFactors = FALSE
  )
  model_parameters <- rbind(
    model_parameters,
    transform(
      gp_multi_kernel_parameter_rows(
        kernel_names = names(kernels),
        strategy = "independent_BGLR_Multitrait_kernel_ETA_covariances"
      ),
      eta_cov_argument = NA_character_
    )
  )
  result <- list(
    Predicted_value = pred_long,
    predicted_values = pred_long,
    Variance_components = vc,
    variance_components = vc,
    model_parameters = model_parameters,
    kernel_variance_components = kernel_contract$kernel_variance_components,
    kernel_configuration = kernel_contract$kernel_configuration,
    Genetic_covariance_by_kernel = kernel_contract$Genetic_covariance_by_kernel,
    genetic_covariance_by_kernel_long = kernel_contract$genetic_covariance_by_kernel_long,
    kernel_combination = kernel_contract$kernel_combination,
    Genetic_covariance_traits = main_covariance,
    Genetic_correlation_traits = main_correlation,
    Genetic_covariance_traits_SE = bayes_multitrait_covariance_draw_sd(
      covariance_bundle$main$draws,
      response
    ),
    Genetic_correlation_traits_SE = bayes_multitrait_correlation_draw_sd(
      covariance_bundle$main$draws,
      response
    ),
    Residual_covariance_traits = covariance_bundle$residual$covariance,
    Residual_correlation_traits = residual_correlation,
    Residual_covariance_traits_SE = bayes_multitrait_covariance_draw_sd(
      covariance_bundle$residual$draws,
      response
    ),
    Residual_correlation_traits_SE = bayes_multitrait_correlation_draw_sd(
      covariance_bundle$residual$draws,
      response
    ),
    Total_genetic_covariance_traits = total_covariance,
    Total_genetic_correlation_traits = total_correlation,
    Total_genetic_covariance_traits_SE = bayes_multitrait_covariance_draw_sd(
      total_covariance_draws,
      response
    ),
    Total_genetic_correlation_traits_SE = bayes_multitrait_correlation_draw_sd(
      total_covariance_draws,
      response
    ),
    M_matrix_model_ready = kernels,
    model_notes = data.frame(
      note = c(
        "ETA covariance is passed with the documented named argument Cov=list(...).",
        "Genetic and residual variance-covariance estimates are posterior means of BGLR covariance parameters; variance-component SEs are posterior SDs from saved covariance draws.",
        "For row-wise joint MT-MET, resCov models residual covariance among traits, not heterogeneous residual variance by environment.",
        "Row-level prediction SE combines BGLR marginal ETA and intercept posterior SDs in quadrature because their posterior covariance is not exposed.",
        "Across-environment prediction SE uses marginal posterior variances and assumes zero posterior covariance between environment-specific predictions."
      ),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(total_pred)) {
    result$Total_Predicted_value <- total_pred
    result$total_predicted_values <- total_pred
  }
  if (!is.null(gxe_covariance)) {
    result$Genetic_covariance_trait_environment <- gxe_covariance
    result$Genetic_correlation_trait_environment <- bayes_multitrait_safe_cov2cor(
      gxe_covariance,
      response
    )
    result$Genetic_covariance_trait_environment_SE <- bayes_multitrait_covariance_draw_sd(
      covariance_bundle$gxe$draws,
      response
    )
    result$Genetic_correlation_trait_environment_SE <- bayes_multitrait_correlation_draw_sd(
      covariance_bundle$gxe$draws,
      response
    )
    result$prediction_uncertainty_source <-
      paste(
        "Quadrature of BGLR marginal ETA and intercept posterior SDs",
        "(their covariance is unavailable); across-environment aggregation",
        "also assumes zero posterior covariance across environment rows"
      )
  } else {
    result$prediction_uncertainty_source <- paste(
      "Quadrature of BGLR marginal ETA and intercept posterior SDs",
      "because their posterior covariance is unavailable"
    )
  }
  model_adapter <- list(
    model = list(
      y = y_mat,
      yHat = pred_mat,
      ETA = fit$ETA,
      multitrait_fit = fit
    ),
    output_files_names = output_files_names
  )
  list(bayes_result = result, bayes_model = model_adapter, bglr_multitrait_model = fit)
}
