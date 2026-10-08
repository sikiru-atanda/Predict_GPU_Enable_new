gp_hybrid_gp_supported_models <- function() "GP"

gp_hybrid_gp_model_label <- function(model_type) {
  model_type <- as.character(model_type %||% "GP")[1]
  if (identical(model_type, "LowRankGP")) {
    return("KRR")
  }
  model_type
}

gp_hybrid_gp_summary_statistics <- function(pred_df,
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
      "backend_runtime",
      "backend_device",
      "hybrid_components"
    ),
    summary = c(
      "hybrid_gp",
      "gaussian",
      model_type,
      "hybrid Gaussian-process prediction via additive female GCA, male GCA, and SCA kernels",
      response,
      sum(pred_df$Train_Test_Label == "Train", na.rm = TRUE),
      sum(pred_df$Train_Test_Label == "Test", na.rm = TRUE),
      length(unique(as.character(pred_df[[female_parent]]))),
      length(unique(as.character(pred_df[[male_parent]]))),
      female_parent,
      male_parent,
      paste(unique(as.character(pred_df$hybrid_gp_backend %||% "R")), collapse = ";"),
      paste(unique(as.character(pred_df$hybrid_gp_device %||% "cpu")), collapse = ";"),
      "female_gca;male_gca;sca"
    ),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_gp_diagnostic_plot <- function(pred_df, model_label) {
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
      title = paste("Hybrid Gaussian-process", model_label, "observed vs predicted"),
      subtitle = "Hybrid GP kernels separate female GCA, male GCA, and SCA",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_hybrid_gp_cv_plot <- function(pred_df, model_label) {
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
      title = paste("Hybrid GP", model_label, "cross-validation observed vs predicted"),
      subtitle = "Hybrid scenarios: known parents, one new parent, both parents new",
      x = "Observed value",
      y = "Predicted value",
      color = "Scenario"
    )
}

gp_hybrid_gp_component_weights <- function(weights = NULL,
                                           include_sca = TRUE,
                                           include_gxe = FALSE) {
  defaults <- c(
    female_gca = 1,
    male_gca = 1,
    sca = if (isTRUE(include_sca)) 1 else 0,
    gxe = if (isTRUE(include_gxe)) 1 else 0
  )
  if (is.null(weights)) {
    return(defaults)
  }
  if (is.null(names(weights))) {
    weights <- stats::setNames(as.numeric(weights), names(defaults)[seq_along(weights)])
  }
  out <- defaults
  overlap <- intersect(names(out), names(weights))
  out[overlap] <- suppressWarnings(as.numeric(weights[overlap]))
  out[!is.finite(out) | is.na(out) | out < 0] <- defaults[!is.finite(out) | is.na(out) | out < 0]
  if (!isTRUE(include_sca)) {
    out[["sca"]] <- 0
  }
  if (!isTRUE(include_gxe)) {
    out[["gxe"]] <- 0
  }
  out
}

gp_hybrid_gp_lambda_grid <- function(lambda_grid = NULL) {
  if (is.null(lambda_grid)) {
    lambda_grid <- c(0.01, 0.03, 0.1, 0.3, 1)
  }
  lambda_grid <- suppressWarnings(as.numeric(lambda_grid))
  lambda_grid <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  if (!length(lambda_grid)) {
    lambda_grid <- c(0.03, 0.1, 0.3)
  }
  sort(unique(lambda_grid))
}

gp_hybrid_gp_fixed_design <- function(ph, heter_groups = NULL, train_idx = NULL) {
  n <- nrow(ph)
  if (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph)) {
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "(Intercept)"
    return(X)
  }
  env <- as.character(ph[[heter_groups]])
  train_idx <- train_idx %||% seq_len(n)
  train_envs <- unique(env[train_idx])
  all_envs <- unique(env)
  if (!all(all_envs %in% train_envs)) {
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "(Intercept)"
    return(X)
  }
  env_factor <- factor(env, levels = all_envs)
  stats::model.matrix(~ env_factor)
}

gp_hybrid_gp_env_covariate_frame <- function(env_covariates,
                                             env_levels,
                                             heter_groups = NULL) {
  env_levels <- as.character(env_levels)
  dat <- as.data.frame(env_covariates, stringsAsFactors = FALSE)
  env_col <- NULL
  candidates <- unique(c(heter_groups, "Env", "env", "environment", "Environment"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  env_col <- candidates[candidates %in% names(dat)][1L] %||% NULL
  if (is.null(env_col)) {
    rn <- rownames(dat)
    if (is.null(rn) || !length(rn)) {
      stop("env_covariates must contain the heter_groups/Env column or row names keyed by environment IDs.", call. = FALSE)
    }
    dat$.__env_id__ <- rn
    env_col <- ".__env_id__"
  }
  dat[[env_col]] <- as.character(dat[[env_col]])
  missing_env <- setdiff(env_levels, dat[[env_col]])
  if (length(missing_env)) {
    stop(
      "env_covariates is missing environments used by the phenotype data: ",
      paste(missing_env, collapse = ", "),
      call. = FALSE
    )
  }
  dat <- dat[match(env_levels, dat[[env_col]]), , drop = FALSE]
  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, env_col)
  if (!length(numeric_cols)) {
    stop("env_covariates must contain at least one numeric covariate after excluding the environment ID column.", call. = FALSE)
  }
  X <- as.matrix(dat[, numeric_cols, drop = FALSE])
  storage.mode(X) <- "double"
  keep <- colSums(is.finite(X)) > 1L
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) {
    stop("No environment covariates passed finite-value QC.", call. = FALSE)
  }
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    fill <- mean(xj, na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    xj[!is.finite(xj)] <- fill
    X[, j] <- xj
  }
  vars <- apply(X, 2, stats::var)
  keep <- is.finite(vars) & vars > 0
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) {
    stop("No environment covariates remained after variance filtering.", call. = FALSE)
  }
  X <- scale(X)
  X[!is.finite(X)] <- 0
  list(
    X = as.matrix(X),
    env_ids = env_levels,
    n_features_input = length(numeric_cols),
    n_features_after_qc = ncol(X)
  )
}

gp_hybrid_gp_kernel_from_distance <- function(D,
                                              kernel = "matern32",
                                              bandwidth = 1.0,
                                              kernel_kwargs = NULL) {
  kernel <- tolower(as.character(kernel %||% "matern32")[1L])
  if (!kernel %in% c("linear", "rbf", "matern32", "matern52")) {
    stop("Unsupported environment kernel. Use linear, rbf, matern32, or matern52.", call. = FALSE)
  }
  bw <- suppressWarnings(as.numeric(kernel_kwargs$bandwidth %||% bandwidth)[1L])
  if (!is.finite(bw) || bw <= 0) bw <- 1
  if (identical(kernel, "rbf")) {
    return(exp(-(D^2) / (2 * bw^2)))
  }
  if (identical(kernel, "matern32")) {
    z <- sqrt(3) * D / bw
    return((1 + z) * exp(-z))
  }
  if (identical(kernel, "matern52")) {
    z <- sqrt(5) * D / bw
    return((1 + z + z^2 / 3) * exp(-z))
  }
  stop("Internal error: linear kernel is handled before distance conversion.", call. = FALSE)
}

gp_hybrid_gp_prepare_env_similarity <- function(ph,
                                                heter_groups = NULL,
                                                env_similarity = NULL,
                                                env_ids = NULL,
                                                env_covariates = NULL,
                                                reaction_norm_feature_qc = TRUE,
                                                kenv_kernel = "matern32",
                                                kenv_bandwidth = 1.0,
                                                kenv_kernel_kwargs = NULL) {
  if (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph)) {
    return(list(
      K_env = NULL,
      env_ids = character(),
      source = "none",
      detail = data.frame()
    ))
  }
  if (!is.null(env_similarity) && !is.null(env_covariates)) {
    stop("Pass exactly one of env_similarity or env_covariates for hybrid GP.", call. = FALSE)
  }
  env_levels <- unique(as.character(ph[[heter_groups]]))
  if (!length(env_levels)) {
    stop("No environment levels were found in heter_groups.", call. = FALSE)
  }
  source <- "identity_environment_kernel"
  detail <- list(
    source = source,
    n_env = length(env_levels),
    n_features_input = NA_integer_,
    n_features_after_qc = NA_integer_,
    kenv_kernel = NA_character_
  )
  if (!is.null(env_similarity)) {
    K <- as.matrix(env_similarity)
    storage.mode(K) <- "double"
    if (is.null(rownames(K)) || is.null(colnames(K))) {
      ids <- as.character(env_ids %||% rownames(K) %||% colnames(K) %||% env_levels)
      if (length(ids) != nrow(K)) {
        stop("env_similarity without row/column names must provide env_ids with one ID per row.", call. = FALSE)
      }
      rownames(K) <- colnames(K) <- ids
    }
    missing_env <- setdiff(env_levels, rownames(K))
    if (length(missing_env)) {
      stop("env_similarity is missing environments: ", paste(missing_env, collapse = ", "), call. = FALSE)
    }
    K <- K[env_levels, env_levels, drop = FALSE]
    source <- "user_env_similarity"
    detail$source <- source
  } else if (!is.null(env_covariates)) {
    cov <- gp_hybrid_gp_env_covariate_frame(
      env_covariates = env_covariates,
      env_levels = env_levels,
      heter_groups = heter_groups
    )
    if (identical(tolower(as.character(kenv_kernel)[1L]), "linear")) {
      K <- tcrossprod(cov$X) / max(ncol(cov$X), 1L)
      K <- K - min(K, na.rm = TRUE)
      if (max(K, na.rm = TRUE) > 0) K <- K / max(K, na.rm = TRUE)
      diag(K) <- 1
    } else {
      D <- as.matrix(stats::dist(cov$X))
      K <- gp_hybrid_gp_kernel_from_distance(
        D,
        kernel = kenv_kernel,
        bandwidth = kenv_bandwidth,
        kernel_kwargs = kenv_kernel_kwargs
      )
    }
    rownames(K) <- colnames(K) <- env_levels
    source <- "env_covariates_kernel"
    detail$source <- source
    detail$n_features_input <- cov$n_features_input
    detail$n_features_after_qc <- cov$n_features_after_qc
    detail$kenv_kernel <- as.character(kenv_kernel)[1L]
  } else {
    if (length(env_levels) > 5L) {
      warning(
        "Hybrid GP MET has more than 5 environments and no env_similarity/env_covariates. ",
        "Using a stable identity environment kernel; provide env_covariates for better borrowing across environments.",
        call. = FALSE
      )
    }
    K <- diag(length(env_levels))
    rownames(K) <- colnames(K) <- env_levels
  }
  K[!is.finite(K)] <- 0
  K <- (K + t(K)) / 2
  diag(K) <- pmax(diag(K), 1e-8)
  d <- sqrt(pmax(diag(K), 1e-8))
  K <- K / tcrossprod(d)
  diag(K) <- 1
  K <- gp_hybrid_gp_spd_correlation(K)
  dimnames(K) <- list(env_levels, env_levels)
  list(
    K_env = K,
    env_ids = env_levels,
    source = source,
    detail = data.frame(detail, stringsAsFactors = FALSE)
  )
}

gp_hybrid_gp_component_kernels <- function(ph,
                                           female_parent,
                                           male_parent,
                                           female_gmatrix,
                                           male_gmatrix,
                                           sca_gmatrix = NULL,
                                           include_sca = TRUE,
                                           component_weights = NULL,
                                           heter_groups = NULL,
                                           env_similarity = NULL) {
  include_gxe <- !is.null(env_similarity) &&
    !is.null(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph)
  weights <- gp_hybrid_gp_component_weights(
    component_weights,
    include_sca = include_sca,
    include_gxe = include_gxe
  )
  female_ids <- as.character(ph[[female_parent]])
  male_ids <- as.character(ph[[male_parent]])
  cross_ids <- as.character(ph$HybridCross)

  k_female <- weights[["female_gca"]] *
    female_gmatrix[female_ids, female_ids, drop = FALSE]
  k_male <- weights[["male_gca"]] *
    male_gmatrix[male_ids, male_ids, drop = FALSE]
  k_sca <- matrix(0, nrow = nrow(ph), ncol = nrow(ph))
  if (isTRUE(include_sca) && !is.null(sca_gmatrix) && weights[["sca"]] > 0) {
    k_sca <- weights[["sca"]] * sca_gmatrix[cross_ids, cross_ids, drop = FALSE]
  }
  k_base <- k_female + k_male + k_sca
  k_gxe <- matrix(0, nrow = nrow(ph), ncol = nrow(ph))
  if (isTRUE(include_gxe) && weights[["gxe"]] > 0) {
    env_ids <- as.character(ph[[heter_groups]])
    K_env <- env_similarity[env_ids, env_ids, drop = FALSE]
    k_gxe <- weights[["gxe"]] * k_base * K_env
  }
  rownames(k_female) <- colnames(k_female) <- seq_len(nrow(ph))
  rownames(k_male) <- colnames(k_male) <- seq_len(nrow(ph))
  rownames(k_sca) <- colnames(k_sca) <- seq_len(nrow(ph))
  rownames(k_gxe) <- colnames(k_gxe) <- seq_len(nrow(ph))
  list(
    female = as.matrix(k_female),
    male = as.matrix(k_male),
    sca = as.matrix(k_sca),
    gxe = as.matrix(k_gxe),
    base = as.matrix(k_base),
    total = as.matrix(k_base + k_gxe),
    weights = weights
  )
}

gp_hybrid_gp_solve_spd <- function(A, B, jitter = 1e-8, max_tries = 6L) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  diag_n <- nrow(A)
  for (i in seq_len(max_tries)) {
    add <- jitter * 10^(i - 1L)
    out <- tryCatch(
      solve(A + diag(add, diag_n), B),
      error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) {
      return(out)
    }
  }
  MASS::ginv(A) %*% B
}

gp_hybrid_gp_fit_core <- function(y,
                                  kernels,
                                  X,
                                  train_idx,
                                  lambda) {
  y <- as.numeric(y)
  train_idx <- as.integer(train_idx)
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda <= 0) {
    lambda <- 0.1
  }
  y_mean <- mean(y[train_idx], na.rm = TRUE)
  y_sd <- stats::sd(y[train_idx], na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) {
    y_sd <- 1
  }
  y_std <- (y - y_mean) / y_sd
  y_train <- matrix(y_std[train_idx], ncol = 1L)
  X <- as.matrix(X)
  X_train <- X[train_idx, , drop = FALSE]
  K_train <- kernels$total[train_idx, train_idx, drop = FALSE]
  V <- K_train + diag(lambda, length(train_idx))
  Vinv_y <- gp_hybrid_gp_solve_spd(V, y_train)
  Vinv_X <- gp_hybrid_gp_solve_spd(V, X_train)
  Xt_Vinv_X <- crossprod(X_train, Vinv_X)
  Xt_Vinv_y <- crossprod(X_train, Vinv_y)
  beta <- gp_hybrid_gp_solve_spd(
    Xt_Vinv_X + diag(1e-8, ncol(Xt_Vinv_X)),
    Xt_Vinv_y
  )
  alpha <- Vinv_y - Vinv_X %*% beta
  K_all_train <- kernels$total[, train_idx, drop = FALSE]
  fixed_std <- drop(X %*% beta)
  female_std <- drop(kernels$female[, train_idx, drop = FALSE] %*% alpha)
  male_std <- drop(kernels$male[, train_idx, drop = FALSE] %*% alpha)
  sca_std <- drop(kernels$sca[, train_idx, drop = FALSE] %*% alpha)
  gxe_std <- drop((kernels$gxe %||% matrix(0, nrow = length(y), ncol = length(y)))[, train_idx, drop = FALSE] %*% alpha)
  total_std <- fixed_std + female_std + male_std + sca_std + gxe_std

  Vinv_K <- gp_hybrid_gp_solve_spd(V, t(K_all_train))
  k_vinv_k <- rowSums(K_all_train * t(Vinv_K))
  base_var <- pmax(diag(kernels$total) - k_vinv_k, 0)
  M_inv <- gp_hybrid_gp_solve_spd(
    Xt_Vinv_X + diag(1e-8, ncol(Xt_Vinv_X)),
    diag(1, ncol(Xt_Vinv_X))
  )
  x_delta <- X - K_all_train %*% Vinv_X
  fixed_var <- rowSums((x_delta %*% M_inv) * x_delta)
  latent_var_std <- pmax(base_var + fixed_var, 0)
  observed_var_std <- pmax(latent_var_std + lambda, 0)

  list(
    lambda = lambda,
    y_mean = y_mean,
    y_sd = y_sd,
    beta = drop(beta),
    alpha = drop(alpha),
    fixed_component = y_mean + y_sd * fixed_std,
    female_component = y_sd * female_std,
    male_component = y_sd * male_std,
    sca_component = y_sd * sca_std,
    gxe_component = y_sd * gxe_std,
    prediction = y_mean + y_sd * total_std,
    latent_var = latent_var_std * y_sd^2,
    observed_var = observed_var_std * y_sd^2,
    train_idx = train_idx,
    backend_runtime = "R",
    backend_device = "cpu"
  )
}

gp_hybrid_gp_backend_mode <- function(gp_backend = "auto") {
  mode <- tolower(as.character(gp_backend %||% "auto")[1L])
  if (is.na(mode) || !nzchar(mode)) {
    mode <- "auto"
  }
  mode
}

gp_hybrid_gp_use_r_backend <- function(gp_backend = "auto") {
  gp_hybrid_gp_backend_mode(gp_backend) %in% c("r", "base", "base_r", "native", "native_r")
}

gp_hybrid_gp_require_python_backend <- function(gp_backend = "auto") {
  gp_hybrid_gp_backend_mode(gp_backend) %in% c("python", "py", "torch", "gpu", "cuda")
}

gp_hybrid_gp_policy <- function(ph,
                                train_idx,
                                female_parent = NULL,
                                male_parent = NULL,
                                heter_groups = NULL,
                                env_kernel_source = NULL,
                                gp_backend = "auto",
                                enabled = TRUE) {
  env <- if (!is.null(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph)) {
    as.character(ph[[heter_groups]])
  } else {
    rep("ENV1", nrow(ph))
  }
  test_mask <- !seq_len(nrow(ph)) %in% as.integer(train_idx)
  device <- if (exists("gp_policy_device_route", mode = "function")) {
    tryCatch(gp_policy_device_route(), error = function(e) list(route = "cpu", reason = "gp_device_policy_unavailable"))
  } else {
    list(route = "cpu", reason = "gp_device_policy_unavailable")
  }
  list(
    enabled = isTRUE(enabled),
    model_family = "hybrid_gp",
    backend_request = as.character(gp_backend %||% "auto")[1L],
    device = device,
    decisions = if (isTRUE(enabled)) {
      c(
        "device_route_selected",
        if (gp_hybrid_gp_use_r_backend(gp_backend)) "r_backend_requested",
        if (!is.null(env_kernel_source) && nzchar(env_kernel_source) && !identical(env_kernel_source, "none")) {
          paste0("env_kernel:", env_kernel_source)
        }
      )
    } else {
      "auto_policy_disabled"
    },
    profile = list(
      n_rows = as.integer(nrow(ph)),
      n_train = as.integer(length(train_idx)),
      n_test = as.integer(sum(test_mask)),
      n_female_parents = as.integer(if (!is.null(female_parent) && female_parent %in% names(ph)) {
        length(unique(as.character(ph[[female_parent]])))
      } else {
        NA_integer_
      }),
      n_male_parents = as.integer(if (!is.null(male_parent) && male_parent %in% names(ph)) {
        length(unique(as.character(ph[[male_parent]])))
      } else {
        NA_integer_
      }),
      n_env = as.integer(length(unique(env)))
    )
  )
}

gp_hybrid_gp_bridge_write_numeric_vector <- function(x, path) {
  x <- as.numeric(x)
  values <- ifelse(
    is.finite(x),
    format(x, digits = 17L, scientific = TRUE, trim = TRUE),
    "nan"
  )
  writeLines(values, con = path, useBytes = TRUE)
  invisible(path)
}

gp_hybrid_gp_bridge_write_integer_vector <- function(x, path) {
  writeLines(as.character(as.integer(x)), con = path, useBytes = TRUE)
  invisible(path)
}

gp_hybrid_gp_bridge_matrix_ref <- function(x, name, input_dir) {
  bin_path <- file.path(input_dir, paste0(name, ".bin"))
  meta_path <- file.path(input_dir, paste0(name, "_shape.txt"))
  gp_bridge_write_matrix_bin(x, bin_path, meta_path)
  list(
    bin = normalizePath(bin_path, winslash = "/", mustWork = TRUE),
    meta = normalizePath(meta_path, winslash = "/", mustWork = TRUE)
  )
}

gp_hybrid_gp_bridge_vector_ref <- function(x, name, input_dir, integer = FALSE) {
  path <- file.path(input_dir, paste0(name, ".txt"))
  if (isTRUE(integer)) {
    gp_hybrid_gp_bridge_write_integer_vector(x, path)
  } else {
    gp_hybrid_gp_bridge_write_numeric_vector(x, path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

gp_hybrid_gp_bridge_read_fit <- function(out_dir) {
  fit_path <- file.path(out_dir, "fit.json")
  if (!file.exists(fit_path)) {
    stop("Hybrid GP Python bridge did not write fit.json.", call. = FALSE)
  }
  fit <- gp_bridge_read_json_simplified(fit_path)
  if (!is.list(fit) || !length(fit)) {
    stop("Hybrid GP Python bridge wrote an invalid fit.json.", call. = FALSE)
  }
  fit
}

gp_hybrid_gp_bridge_run_spec <- function(spec,
                                         command,
                                         route,
                                         python_bin = NULL) {
  gp_bridge_write_json(spec$payload, spec$spec_json)
  gp_with_gp_device_route(
    route,
    gp_bridge_run_cli(command, c("--spec", spec$spec_json), python_bin = python_bin)
  )
  fit <- gp_hybrid_gp_bridge_read_fit(spec$out_dir)
  meta <- gp_bridge_read_json_simplified(file.path(spec$out_dir, "meta.json"))
  if (is.list(meta) && length(meta)) {
    fit$bridge_meta <- meta
  }
  fit$bridge_dir <- spec$bridge_dir
  fit
}

gp_hybrid_gp_prepare_additive_spec <- function(y,
                                               kernels,
                                               X,
                                               train_idx,
                                               lambda,
                                               route,
                                               gp_dtype = "float64",
                                               random_state = 123L,
                                               project_root = NULL,
                                               dir_path = tempfile("predictpror_hybrid_gp_")) {
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  kernel_refs <- list(
    female = gp_hybrid_gp_bridge_matrix_ref(kernels$female, "K_female", input_dir),
    male = gp_hybrid_gp_bridge_matrix_ref(kernels$male, "K_male", input_dir),
    sca = gp_hybrid_gp_bridge_matrix_ref(kernels$sca, "K_sca", input_dir)
  )
  if (!is.null(kernels$gxe)) {
    kernel_refs$gxe <- gp_hybrid_gp_bridge_matrix_ref(kernels$gxe, "K_gxe", input_dir)
  }
  spec_json <- file.path(input_dir, "fit_spec.json")
  route_value <- as.character(route %||% NA_character_)[1L]
  payload <- list(
    project_root = normalizePath(project_root %||% gp_project_root(), winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    inputs = list(
      y = gp_hybrid_gp_bridge_vector_ref(y, "y", input_dir),
      X = gp_hybrid_gp_bridge_matrix_ref(X, "X", input_dir),
      train_idx = gp_hybrid_gp_bridge_vector_ref(as.integer(train_idx) - 1L, "train_idx", input_dir, integer = TRUE),
      kernels = kernel_refs
    ),
    lambda_value = as.numeric(lambda)[1L],
    dtype = as.character(gp_dtype %||% "float64")[1L],
    seed = as.integer(random_state %||% 123L)
  )
  if (!is.na(route_value) && nzchar(route_value)) {
    payload$device <- route_value
  }
  list(
    payload = payload,
    spec_json = normalizePath(spec_json, winslash = "/", mustWork = FALSE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  )
}

gp_hybrid_gp_python_fit_core <- function(y,
                                         kernels,
                                         X,
                                         train_idx,
                                         lambda,
                                         policy,
                                         gp_backend = "auto",
                                         gp_dtype = "float64",
                                         random_state = 123L,
                                         python_bin = NULL,
                                         project_root = NULL) {
  route <- if (isTRUE(policy$enabled)) policy$device$route else NA_character_
  spec <- gp_hybrid_gp_prepare_additive_spec(
    y = y,
    kernels = kernels,
    X = X,
    train_idx = train_idx,
    lambda = lambda,
    route = route,
    gp_dtype = gp_dtype,
    random_state = random_state,
    project_root = project_root %||% gp_project_root()
  )
  fit <- gp_hybrid_gp_bridge_run_spec(
    spec = spec,
    command = "fit-hybrid-additive-spec",
    route = route,
    python_bin = python_bin
  )
  fit$lambda <- as.numeric(fit$lambda)[1L]
  fit$y_mean <- as.numeric(fit$y_mean)[1L]
  fit$y_sd <- as.numeric(fit$y_sd)[1L]
  for (nm in c(
    "beta", "alpha", "fixed_component", "female_component", "male_component",
    "sca_component", "gxe_component", "prediction", "latent_var", "observed_var",
    "train_idx"
  )) {
    fit[[nm]] <- as.numeric(fit[[nm]])
  }
  fit$train_idx <- as.integer(round(fit$train_idx))
  fit$backend_runtime <- "PythonTorch"
  fit$backend_request <- as.character(gp_backend %||% "auto")[1L]
  fit$backend_device <- as.character(fit$device_info$device %||% route %||% "auto")[1L]
  fit$gp_auto_policy <- policy
  fit
}

gp_hybrid_gp_fit_backend <- function(y,
                                     kernels,
                                     X,
                                     train_idx,
                                     lambda,
                                     ph = NULL,
                                     female_parent = NULL,
                                     male_parent = NULL,
                                     heter_groups = NULL,
                                     env_kernel_source = NULL,
                                     gp_backend = "auto",
                                     gp_auto_policy = TRUE,
                                     gp_dtype = "float64",
                                     random_state = 123L,
                                     python_bin = NULL,
                                     project_root = NULL) {
  policy <- gp_hybrid_gp_policy(
    ph = ph %||% data.frame(row_id = seq_along(y)),
    train_idx = train_idx,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    env_kernel_source = env_kernel_source,
    gp_backend = gp_backend,
    enabled = isTRUE(gp_auto_policy)
  )
  if (gp_hybrid_gp_use_r_backend(gp_backend)) {
    fit <- gp_hybrid_gp_fit_core(y, kernels, X, train_idx, lambda)
    fit$backend_request <- as.character(gp_backend %||% "r")[1L]
    fit$gp_auto_policy <- policy
    return(fit)
  }
  fit_py <- tryCatch(
    gp_hybrid_gp_python_fit_core(
      y = y,
      kernels = kernels,
      X = X,
      train_idx = train_idx,
      lambda = lambda,
      policy = policy,
      gp_backend = gp_backend,
      gp_dtype = gp_dtype,
      random_state = random_state,
      python_bin = python_bin,
      project_root = project_root
    ),
    error = function(e) e
  )
  if (!inherits(fit_py, "error")) {
    return(fit_py)
  }
  if (gp_hybrid_gp_require_python_backend(gp_backend)) {
    stop(
      "Hybrid GP Python/Torch backend failed: ",
      conditionMessage(fit_py),
      call. = FALSE
    )
  }
  fit <- gp_hybrid_gp_fit_core(y, kernels, X, train_idx, lambda)
  fit$backend_request <- as.character(gp_backend %||% "auto")[1L]
  fit$backend_fallback_reason <- conditionMessage(fit_py)
  fit$gp_auto_policy <- policy
  fit
}

gp_hybrid_gp_make_folds <- function(train_idx, nfolds = 3L, seed = 123L) {
  train_idx <- as.integer(train_idx)
  n <- length(train_idx)
  if (n < 6L) {
    return(list(train_idx))
  }
  k <- min(as.integer(nfolds), n)
  gp_set_seed(as.integer(seed))
  shuffled <- sample(train_idx)
  split(shuffled, rep(seq_len(k), length.out = n))
}

gp_hybrid_gp_select_lambda <- function(y,
                                       kernels,
                                       X,
                                       train_idx,
                                       lambda,
                                       lambda_grid = NULL,
                                       random_state = 123L,
                                       tuning_max_train = NULL) {
  if (!is.character(lambda) || !identical(tolower(lambda[1]), "auto")) {
    lam <- suppressWarnings(as.numeric(lambda)[1L])
    if (is.finite(lam) && lam > 0) {
      return(list(lambda = lam, status = "fixed", detail = data.frame()))
    }
  }
  grid <- gp_hybrid_gp_lambda_grid(lambda_grid)
  tune_train_idx <- train_idx
  max_train <- if (exists("gp_policy_tuning_limit", mode = "function")) {
    gp_policy_tuning_limit(
      tuning_max_train,
      "PREDICTPRO_HYBRID_GP_TUNING_MAX_ROWS",
      4000L
    )
  } else {
    suppressWarnings(as.integer(tuning_max_train %||% 4000L))
  }
  if (!is.finite(max_train) || is.na(max_train) || max_train < 1L) {
    max_train <- 4000L
  }
  bounded_tuning <- length(tune_train_idx) > max_train
  if (isTRUE(bounded_tuning)) {
  gp_set_seed(as.integer(random_state))
    tune_train_idx <- sort(sample(tune_train_idx, max_train))
  }
  folds <- gp_hybrid_gp_make_folds(tune_train_idx, nfolds = 3L, seed = random_state)
  if (length(folds) < 2L) {
    chosen <- grid[ceiling(length(grid) / 2L)]
    return(list(lambda = chosen, status = "auto_skipped_small_training_set", detail = data.frame()))
  }
  detail <- do.call(rbind, lapply(grid, function(lam) {
    fold_rmse <- vapply(folds, function(val_idx) {
      cal_idx <- setdiff(train_idx, val_idx)
      if (length(cal_idx) < 2L) {
        return(NA_real_)
      }
      fit <- gp_hybrid_gp_fit_core(
        y = y,
        kernels = kernels,
        X = X,
        train_idx = cal_idx,
        lambda = lam
      )
      err <- y[val_idx] - fit$prediction[val_idx]
      sqrt(mean(err^2, na.rm = TRUE))
    }, numeric(1))
    data.frame(
      lambda = lam,
      mean_rmse = mean(fold_rmse, na.rm = TRUE),
      n_folds = sum(is.finite(fold_rmse)),
      stringsAsFactors = FALSE
    )
  }))
  detail <- detail[is.finite(detail$mean_rmse), , drop = FALSE]
  chosen <- if (nrow(detail)) {
    detail$lambda[which.min(detail$mean_rmse)]
  } else {
    grid[ceiling(length(grid) / 2L)]
  }
  list(
    lambda = chosen,
    status = if (isTRUE(bounded_tuning)) "auto_cv_bounded_training_subset" else "auto_cv",
    detail = detail
  )
}

gp_hybrid_gp_effect_tables <- function(pred_df,
                                       female_parent,
                                       male_parent) {
  gp_hybrid_bayes_effect_tables(
    pred_df = pred_df,
    female_parent = female_parent,
    male_parent = male_parent
  )
}

gp_hybrid_gp_prediction_frame <- function(ph,
                                          response,
                                          gen_name,
                                          female_parent,
                                          male_parent,
                                          heter_groups = NULL,
                                          fit,
                                          model_type,
                                          decomposition_basis = "hybrid_gp_additive_kernel_decomposition") {
  out_cols <- c(gen_name, female_parent, male_parent, response, "HybridCross")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    out_cols <- c(out_cols, heter_groups)
  }
  out <- ph[, out_cols, drop = FALSE]
  names(out)[names(out) == response] <- "Observed_value"
  out$Train_Test_Label <- ifelse(is.na(out$Observed_value), "Test", "Train")
  out$Prediction_intercept <- as.numeric(fit$fixed_component)
  out$Female_GCA <- as.numeric(fit$female_component)
  out$Male_GCA <- as.numeric(fit$male_component)
  out$SCA_effect <- as.numeric(fit$sca_component)
  out$GxE_effect <- as.numeric(fit$gxe_component %||% rep(0, nrow(out)))
  out$Female_additive_contribution <- out$Female_GCA
  out$Male_additive_contribution <- out$Male_GCA
  out$Hybrid_interaction_contribution <- out$SCA_effect + out$GxE_effect
  out$Environment_interaction_contribution <- out$GxE_effect
  out$Female_GCA_like_contribution <- out$Female_GCA
  out$Male_GCA_like_contribution <- out$Male_GCA
  out$SCA_like_contribution <- out$SCA_effect
  out$GxE_like_contribution <- out$GxE_effect
  out$Predictive_decomposition_basis <- decomposition_basis
  out$Predictive_decomposition_remarks <- paste(
    "Hybrid GP decomposition:",
    "fixed mean/env effect + female GCA kernel + male GCA kernel + SCA kernel + optional hybrid-by-environment kernel = predicted value"
  )
  out$Predicted_value <- as.numeric(fit$prediction)
  out$Standard_error <- sqrt(pmax(as.numeric(fit$observed_var), 0))
  out$Prediction_error_variance <- pmax(as.numeric(fit$observed_var), 0)
  out$Prediction_SE_latent <- sqrt(pmax(as.numeric(fit$latent_var), 0))
  out$Prediction_Var_latent <- pmax(as.numeric(fit$latent_var), 0)
  out$hybrid_gp_model <- model_type
  out$hybrid_gp_lambda <- fit$lambda
  out$hybrid_gp_backend <- as.character(fit$backend_runtime %||% "R")[1L]
  out$hybrid_gp_device <- as.character(fit$backend_device %||% "cpu")[1L]
  out
}

gp_hybrid_gp_variance_components <- function(pred_df,
                                             kernels,
                                             fit,
                                             lambda_info) {
  component_var <- c(
    female_gca = stats::var(pred_df$Female_GCA, na.rm = TRUE),
    male_gca = stats::var(pred_df$Male_GCA, na.rm = TRUE),
    sca = stats::var(pred_df$SCA_effect, na.rm = TRUE),
    gxe = if ("GxE_effect" %in% names(pred_df)) stats::var(pred_df$GxE_effect, na.rm = TRUE) else NA_real_
  )
  data.frame(
    component = c(names(component_var), "residual_lambda", "selected_lambda"),
    estimate = c(component_var, fit$lambda * fit$y_sd^2, fit$lambda),
    scale = c(rep("prediction_component_variance", length(component_var)), "response_scale", "standardized_kernel_ratio"),
    tuning_status = c(rep(NA_character_, length(component_var)), lambda_info$status, lambda_info$status),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_gp_gaussian_model <- function(pheno_object,
                                        response,
                                        gen_name,
                                        female_parent,
                                        male_parent,
                                        heter_groups = NULL,
                                        model_type = "GP",
                                        gmatrix = NULL,
                                        female_gmatrix = NULL,
                                        male_gmatrix = NULL,
                                        geno_data = NULL,
                                        female_geno_data = NULL,
                                        male_geno_data = NULL,
                                        gmatrix_method = NULL,
                                        include_sca = TRUE,
                                        lambda = "auto",
                                        lambda_grid = NULL,
                                        component_weights = NULL,
                                        gp_backend = "auto",
                                        gp_auto_policy = TRUE,
                                        gp_dtype = "float64",
                                        tuning_max_train = NULL,
                                        env_similarity = NULL,
                                        env_ids = NULL,
                                        env_covariates = NULL,
                                        reaction_norm_feature_qc = TRUE,
                                        kenv_kernel = "matern32",
                                        kenv_bandwidth = 1.0,
                                        kenv_kernel_kwargs = NULL,
                                        random_state = 123L,
                                        python_bin = NULL,
                                        project_root = NULL) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  needed_cols <- c(gen_name, female_parent, male_parent, response)
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    needed_cols <- c(needed_cols, heter_groups)
  }
  missing_cols <- setdiff(needed_cols, names(ph))
  if (length(missing_cols)) {
    stop(paste("Hybrid GP phenotype data is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (is.null(heter_groups) || !nzchar(heter_groups)) {
    if (anyDuplicated(ph[[gen_name]])) {
      stop("Hybrid GP currently requires one phenotype row per hybrid ID unless heter_groups is supplied for multi-environment data.", call. = FALSE)
    }
  } else if (anyDuplicated(gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups))) {
    stop("Hybrid GP currently requires one phenotype row per hybrid-by-environment combination.", call. = FALSE)
  }

  ph[[female_parent]] <- as.character(ph[[female_parent]])
  ph[[male_parent]] <- as.character(ph[[male_parent]])
  ph[[gen_name]] <- as.character(ph[[gen_name]])
  ph$HybridCross <- paste(ph[[female_parent]], ph[[male_parent]], sep = "__x__")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    ph[[heter_groups]] <- as.character(ph[[heter_groups]])
  }
  env_kernel <- gp_hybrid_gp_prepare_env_similarity(
    ph = ph,
    heter_groups = heter_groups,
    env_similarity = env_similarity,
    env_ids = env_ids,
    env_covariates = env_covariates,
    reaction_norm_feature_qc = reaction_norm_feature_qc,
    kenv_kernel = kenv_kernel,
    kenv_bandwidth = kenv_bandwidth,
    kenv_kernel_kwargs = kenv_kernel_kwargs
  )

  train_idx <- which(is.finite(as.numeric(ph[[response]])))
  if (length(train_idx) < 2L) {
    stop("Hybrid GP requires at least two observed training hybrids.", call. = FALSE)
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

  kernels <- gp_hybrid_gp_component_kernels(
    ph = ph,
    female_parent = female_parent,
    male_parent = male_parent,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    include_sca = include_sca,
    component_weights = component_weights,
    heter_groups = heter_groups,
    env_similarity = env_kernel$K_env
  )
  X <- gp_hybrid_gp_fixed_design(ph, heter_groups = heter_groups, train_idx = train_idx)
  y <- as.numeric(ph[[response]])
  lambda_info <- gp_hybrid_gp_select_lambda(
    y = y,
    kernels = kernels,
    X = X,
    train_idx = train_idx,
    lambda = lambda,
    lambda_grid = lambda_grid,
    random_state = random_state,
    tuning_max_train = tuning_max_train
  )
  fit <- gp_hybrid_gp_fit_backend(
    y = y,
    kernels = kernels,
    X = X,
    train_idx = train_idx,
    lambda = lambda_info$lambda,
    ph = ph,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    env_kernel_source = env_kernel$source,
    gp_backend = gp_backend,
    gp_auto_policy = gp_auto_policy,
    gp_dtype = gp_dtype,
    random_state = random_state,
    python_bin = python_bin,
    project_root = project_root
  )
  model_label <- gp_hybrid_gp_model_label(model_type)
  pred_df <- gp_hybrid_gp_prediction_frame(
    ph = ph,
    response = response,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    fit = fit,
    model_type = model_label
  )
  eff_tables <- gp_hybrid_gp_effect_tables(
    pred_df = pred_df,
    female_parent = female_parent,
    male_parent = male_parent
  )

  list(
    gp_model = list(
      model_type = model_label,
      lambda = fit$lambda,
      lambda_status = lambda_info$status,
      beta = fit$beta,
      alpha = fit$alpha,
      train_idx = fit$train_idx,
      backend_runtime = fit$backend_runtime %||% "R",
      backend_request = fit$backend_request %||% gp_backend,
      backend_device = fit$backend_device %||% "cpu",
      backend_fallback_reason = fit$backend_fallback_reason %||% NA_character_,
      gp_auto_policy = fit$gp_auto_policy %||% NULL,
      env_kernel_source = env_kernel$source,
      env_kernel_detail = env_kernel$detail,
      device_info = fit$device_info %||% NULL
    ),
    predicted_values = pred_df,
    Predicted_value = pred_df,
    female_gca_effects = eff_tables$female_gca_effects,
    male_gca_effects = eff_tables$male_gca_effects,
    sca_effects = eff_tables$sca_effects,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    env_similarity = env_kernel$K_env,
    env_kernel = env_kernel,
    component_weights = kernels$weights,
    lambda_tuning = lambda_info$detail,
    variance_components = gp_hybrid_gp_variance_components(
      pred_df = pred_df,
      kernels = kernels,
      fit = fit,
      lambda_info = lambda_info
    ),
    diagnostic_plots = gp_hybrid_gp_diagnostic_plot(pred_df, model_label = model_label)
  )
}

gp_hybrid_gp_cv_metric_table <- function(pred_df,
                                         eval_metrics,
                                         response_family = "gaussian") {
  gp_hybrid_asreml_cv_metric_table(
    pred_df = pred_df,
    eval_metrics = eval_metrics,
    response_family = response_family
  )
}

gp_hybrid_gp_cv_process <- function(pred_df,
                                    eval_df,
                                    response,
                                    model_type) {
  out <- gp_hybrid_asreml_cv_process(
    pred_df = pred_df,
    eval_df = eval_df,
    response = response,
    model_type = gp_hybrid_gp_model_label(model_type)
  )
  out$hybrid_cv_plot <- gp_hybrid_gp_cv_plot(pred_df, model_label = gp_hybrid_gp_model_label(model_type))
  out
}

gp_hybrid_gp_gaussian_cv <- function(pheno_object,
                                     response,
                                     gen_name,
                                     female_parent,
                                     male_parent,
                                     heter_groups = NULL,
                                     cross_validation_meth,
                                     eval_metrics,
                                     model_type = "GP",
                                     gmatrix = NULL,
                                     female_gmatrix = NULL,
                                     male_gmatrix = NULL,
                                     geno_data = NULL,
                                     female_geno_data = NULL,
                                     male_geno_data = NULL,
                                     gmatrix_method = NULL,
                                     include_sca = TRUE,
                                     lambda = "auto",
                                     lambda_grid = NULL,
                                     component_weights = NULL,
                                     gp_backend = "auto",
                                     gp_auto_policy = TRUE,
                                     gp_dtype = "float64",
                                     tuning_max_train = NULL,
                                     env_similarity = NULL,
                                     env_ids = NULL,
                                     env_covariates = NULL,
                                     reaction_norm_feature_qc = TRUE,
                                     kenv_kernel = "matern32",
                                     kenv_bandwidth = 1.0,
                                     kenv_kernel_kwargs = NULL,
                                     nfolds = 5L,
                                     random_state = 123L,
                                     replication = 1L,
                                     python_bin = NULL,
                                     project_root = NULL) {
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

      fit <- gp_hybrid_gp_gaussian_model(
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
        lambda = lambda,
        lambda_grid = lambda_grid,
        component_weights = component_weights,
        gp_backend = gp_backend,
        gp_auto_policy = gp_auto_policy,
        gp_dtype = gp_dtype,
        tuning_max_train = tuning_max_train,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        random_state = random_state + rep_i,
        python_bin = python_bin,
        project_root = project_root
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

      eval_df <- gp_hybrid_gp_cv_metric_table(
        pred_df = pred_test,
        eval_metrics = eval_metrics,
        response_family = "gaussian"
      )

      all_preds[[length(all_preds) + 1L]] <- pred_test
      all_eval[[length(all_eval) + 1L]] <- eval_df
      all_raw[[length(all_raw) + 1L]] <- list(
        trait = response,
        rep = rep_i,
        model = gp_hybrid_gp_model_label(model_type),
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
  processed <- gp_hybrid_gp_cv_process(
    pred_df = pred_df,
    eval_df = eval_df,
    response = response,
    model_type = model_type
  )

  list(
    predicted_values = pred_df,
    hybrid_cv_metrics = eval_df,
    cv_results_processed = processed,
    cv_results_raw = all_raw
  )
}

gp_hybrid_gp_add_trait_column <- function(dat, trait) {
  if (is.null(dat)) {
    return(NULL)
  }
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  if (!nrow(dat) && !ncol(dat)) {
    return(dat)
  }
  dat$Trait <- trait
  first_cols <- intersect(c("Trait"), names(dat))
  dat[, c(first_cols, setdiff(names(dat), first_cols)), drop = FALSE]
}

gp_hybrid_gp_bind_trait_tables <- function(outputs, responses, field) {
  tabs <- lapply(seq_along(outputs), function(i) {
    gp_hybrid_gp_add_trait_column(outputs[[i]][[field]], responses[[i]])
  })
  tabs <- Filter(function(x) !is.null(x) && (nrow(x) > 0L || ncol(x) > 0L), tabs)
  if (!length(tabs)) {
    return(data.frame())
  }
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  out
}

gp_hybrid_gp_spd_correlation <- function(C, min_eigen = 1e-6) {
  C <- as.matrix(C)
  C[!is.finite(C)] <- 0
  C <- (C + t(C)) / 2
  diag(C) <- 1
  eig <- eigen(C, symmetric = TRUE)
  vals <- pmax(eig$values, min_eigen)
  out <- eig$vectors %*% diag(vals, length(vals)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  d <- sqrt(pmax(diag(out), min_eigen))
  out <- out / tcrossprod(d)
  diag(out) <- 1
  out
}

gp_hybrid_gp_estimate_trait_correlation <- function(Y_std, responses) {
  if (ncol(Y_std) == 1L) {
    out <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(responses, responses))
    return(out)
  }
  C <- suppressWarnings(stats::cor(Y_std, use = "pairwise.complete.obs"))
  C <- as.matrix(C)
  dimnames(C) <- list(responses, responses)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C <- gp_hybrid_gp_spd_correlation(C)
  dimnames(C) <- list(responses, responses)
  C
}

gp_hybrid_gp_estimate_gxe_trait_correlation <- function(Y_std,
                                                        responses,
                                                        env = NULL,
                                                        fallback_cor = NULL) {
  if (ncol(Y_std) == 1L) {
    out <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(responses, responses))
    return(out)
  }
  if (is.null(env)) {
    return(fallback_cor %||% gp_hybrid_gp_estimate_trait_correlation(Y_std, responses))
  }
  env <- as.character(env)
  Y_res <- Y_std
  for (j in seq_len(ncol(Y_res))) {
    yj <- Y_res[, j]
    env_mean <- stats::ave(yj, env, FUN = function(x) mean(x, na.rm = TRUE))
    env_mean[!is.finite(env_mean)] <- 0
    Y_res[, j] <- yj - env_mean
  }
  C <- suppressWarnings(stats::cor(Y_res, use = "pairwise.complete.obs"))
  C <- as.matrix(C)
  if (!all(is.finite(C))) {
    return(fallback_cor %||% gp_hybrid_gp_estimate_trait_correlation(Y_std, responses))
  }
  dimnames(C) <- list(responses, responses)
  diag(C) <- 1
  C <- gp_hybrid_gp_spd_correlation(C)
  dimnames(C) <- list(responses, responses)
  C
}

gp_hybrid_gp_joint_design <- function(X, n_traits, responses) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  out <- matrix(0, nrow = n * n_traits, ncol = p * n_traits)
  x_names <- colnames(X) %||% paste0("X", seq_len(p))
  colnames(out) <- unlist(lapply(responses, function(resp) paste(resp, x_names, sep = "::")))
  for (trait_i in seq_len(n_traits)) {
    row_idx <- ((trait_i - 1L) * n + 1L):(trait_i * n)
    col_idx <- ((trait_i - 1L) * p + 1L):(trait_i * p)
    out[row_idx, col_idx] <- X
  }
  out
}

gp_hybrid_gp_joint_kernel_cross <- function(K,
                                            trait_cor,
                                            all_row,
                                            all_trait,
                                            train_row,
                                            train_trait) {
  trait_cor[all_trait, train_trait, drop = FALSE] *
    K[all_row, train_row, drop = FALSE]
}

gp_hybrid_gp_fit_joint_core <- function(y_std_vec,
                                        kernels,
                                        trait_cor,
                                        gxe_trait_cor = trait_cor,
                                        X_joint,
                                        train_obs_idx,
                                        all_row,
                                        all_trait,
                                        lambda) {
  train_obs_idx <- as.integer(train_obs_idx)
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda <= 0) {
    lambda <- 0.1
  }
  train_row <- all_row[train_obs_idx]
  train_trait <- all_trait[train_obs_idx]
  y_train <- matrix(y_std_vec[train_obs_idx], ncol = 1L)
  K_base <- kernels$base %||% (kernels$female + kernels$male + kernels$sca)
  K_gxe <- kernels$gxe %||% matrix(0, nrow = nrow(K_base), ncol = ncol(K_base))
  K_train <- trait_cor[train_trait, train_trait, drop = FALSE] *
    K_base[train_row, train_row, drop = FALSE] +
    gxe_trait_cor[train_trait, train_trait, drop = FALSE] *
      K_gxe[train_row, train_row, drop = FALSE]
  X_train <- X_joint[train_obs_idx, , drop = FALSE]
  V <- K_train + diag(lambda, length(train_obs_idx))
  Vinv_y <- gp_hybrid_gp_solve_spd(V, y_train)
  Vinv_X <- gp_hybrid_gp_solve_spd(V, X_train)
  Xt_Vinv_X <- crossprod(X_train, Vinv_X)
  Xt_Vinv_y <- crossprod(X_train, Vinv_y)
  beta <- gp_hybrid_gp_solve_spd(
    Xt_Vinv_X + diag(1e-8, ncol(Xt_Vinv_X)),
    Xt_Vinv_y
  )
  alpha <- Vinv_y - Vinv_X %*% beta

  K_base_all_train <- gp_hybrid_gp_joint_kernel_cross(
    K = K_base,
    trait_cor = trait_cor,
    all_row = all_row,
    all_trait = all_trait,
    train_row = train_row,
    train_trait = train_trait
  )
  K_gxe_all_train <- gp_hybrid_gp_joint_kernel_cross(
    K = K_gxe,
    trait_cor = gxe_trait_cor,
    all_row = all_row,
    all_trait = all_trait,
    train_row = train_row,
    train_trait = train_trait
  )
  K_total_all_train <- K_base_all_train + K_gxe_all_train
  K_female_all_train <- gp_hybrid_gp_joint_kernel_cross(
    K = kernels$female,
    trait_cor = trait_cor,
    all_row = all_row,
    all_trait = all_trait,
    train_row = train_row,
    train_trait = train_trait
  )
  K_male_all_train <- gp_hybrid_gp_joint_kernel_cross(
    K = kernels$male,
    trait_cor = trait_cor,
    all_row = all_row,
    all_trait = all_trait,
    train_row = train_row,
    train_trait = train_trait
  )
  K_sca_all_train <- gp_hybrid_gp_joint_kernel_cross(
    K = kernels$sca,
    trait_cor = trait_cor,
    all_row = all_row,
    all_trait = all_trait,
    train_row = train_row,
    train_trait = train_trait
  )

  fixed_std <- drop(X_joint %*% beta)
  female_std <- drop(K_female_all_train %*% alpha)
  male_std <- drop(K_male_all_train %*% alpha)
  sca_std <- drop(K_sca_all_train %*% alpha)
  gxe_std <- drop(K_gxe_all_train %*% alpha)
  total_std <- fixed_std + female_std + male_std + sca_std + gxe_std

  Vinv_K <- gp_hybrid_gp_solve_spd(V, t(K_total_all_train))
  k_vinv_k <- rowSums(K_total_all_train * t(Vinv_K))
  diag_total <- diag(K_base)[all_row] * trait_cor[cbind(all_trait, all_trait)] +
    diag(K_gxe)[all_row] * gxe_trait_cor[cbind(all_trait, all_trait)]
  base_var <- pmax(diag_total - k_vinv_k, 0)
  M_inv <- gp_hybrid_gp_solve_spd(
    Xt_Vinv_X + diag(1e-8, ncol(Xt_Vinv_X)),
    diag(1, ncol(Xt_Vinv_X))
  )
  x_delta <- X_joint - K_total_all_train %*% Vinv_X
  fixed_var <- rowSums((x_delta %*% M_inv) * x_delta)
  latent_var_std <- pmax(base_var + fixed_var, 0)
  observed_var_std <- pmax(latent_var_std + lambda, 0)

  list(
    lambda = lambda,
    beta = drop(beta),
    alpha = drop(alpha),
    fixed_std = fixed_std,
    female_std = female_std,
    male_std = male_std,
    sca_std = sca_std,
    gxe_std = gxe_std,
    prediction_std = total_std,
    latent_var_std = latent_var_std,
    observed_var_std = observed_var_std,
    train_obs_idx = train_obs_idx,
    backend_runtime = "R",
    backend_device = "cpu"
  )
}

gp_hybrid_gp_prepare_multi_trait_spec <- function(y_std_vec,
                                                  kernels,
                                                  trait_cor,
                                                  gxe_trait_cor = trait_cor,
                                                  X_joint,
                                                  train_obs_idx,
                                                  all_row,
                                                  all_trait,
                                                  lambda,
                                                  route,
                                                  gp_dtype = "float64",
                                                  random_state = 123L,
                                                  project_root = NULL,
                                                  dir_path = tempfile("predictpror_hybrid_mt_gp_")) {
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  kernel_refs <- list(
    female = gp_hybrid_gp_bridge_matrix_ref(kernels$female, "K_female", input_dir),
    male = gp_hybrid_gp_bridge_matrix_ref(kernels$male, "K_male", input_dir),
    sca = gp_hybrid_gp_bridge_matrix_ref(kernels$sca, "K_sca", input_dir)
  )
  if (!is.null(kernels$gxe)) {
    kernel_refs$gxe <- gp_hybrid_gp_bridge_matrix_ref(kernels$gxe, "K_gxe", input_dir)
  }
  spec_json <- file.path(input_dir, "fit_spec.json")
  route_value <- as.character(route %||% NA_character_)[1L]
  payload <- list(
    project_root = normalizePath(project_root %||% gp_project_root(), winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    inputs = list(
      y_std = gp_hybrid_gp_bridge_vector_ref(y_std_vec, "y_std", input_dir),
      X_joint = gp_hybrid_gp_bridge_matrix_ref(X_joint, "X_joint", input_dir),
      trait_cor = gp_hybrid_gp_bridge_matrix_ref(trait_cor, "trait_cor", input_dir),
      gxe_trait_cor = gp_hybrid_gp_bridge_matrix_ref(gxe_trait_cor, "gxe_trait_cor", input_dir),
      all_row = gp_hybrid_gp_bridge_vector_ref(as.integer(all_row) - 1L, "all_row", input_dir, integer = TRUE),
      all_trait = gp_hybrid_gp_bridge_vector_ref(as.integer(all_trait) - 1L, "all_trait", input_dir, integer = TRUE),
      train_obs_idx = gp_hybrid_gp_bridge_vector_ref(as.integer(train_obs_idx) - 1L, "train_obs_idx", input_dir, integer = TRUE),
      kernels = kernel_refs
    ),
    lambda_value = as.numeric(lambda)[1L],
    dtype = as.character(gp_dtype %||% "float64")[1L],
    seed = as.integer(random_state %||% 123L)
  )
  if (!is.na(route_value) && nzchar(route_value)) {
    payload$device <- route_value
  }
  list(
    payload = payload,
    spec_json = normalizePath(spec_json, winslash = "/", mustWork = FALSE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  )
}

gp_hybrid_gp_python_fit_joint_core <- function(y_std_vec,
                                               kernels,
                                               trait_cor,
                                               gxe_trait_cor = trait_cor,
                                               X_joint,
                                               train_obs_idx,
                                               all_row,
                                               all_trait,
                                               lambda,
                                               policy,
                                               gp_backend = "auto",
                                               gp_dtype = "float64",
                                               random_state = 123L,
                                               python_bin = NULL,
                                               project_root = NULL) {
  route <- if (isTRUE(policy$enabled)) policy$device$route else NA_character_
  spec <- gp_hybrid_gp_prepare_multi_trait_spec(
    y_std_vec = y_std_vec,
    kernels = kernels,
    trait_cor = trait_cor,
    gxe_trait_cor = gxe_trait_cor,
    X_joint = X_joint,
    train_obs_idx = train_obs_idx,
    all_row = all_row,
    all_trait = all_trait,
    lambda = lambda,
    route = route,
    gp_dtype = gp_dtype,
    random_state = random_state,
    project_root = project_root %||% gp_project_root()
  )
  fit <- gp_hybrid_gp_bridge_run_spec(
    spec = spec,
    command = "fit-hybrid-multi-trait-spec",
    route = route,
    python_bin = python_bin
  )
  fit$lambda <- as.numeric(fit$lambda)[1L]
  for (nm in c(
    "beta", "alpha", "fixed_std", "female_std", "male_std", "sca_std",
    "gxe_std", "prediction_std", "latent_var_std", "observed_var_std", "train_obs_idx"
  )) {
    fit[[nm]] <- as.numeric(fit[[nm]])
  }
  fit$train_obs_idx <- as.integer(round(fit$train_obs_idx))
  fit$backend_runtime <- "PythonTorch"
  fit$backend_request <- as.character(gp_backend %||% "auto")[1L]
  fit$backend_device <- as.character(fit$device_info$device %||% route %||% "auto")[1L]
  fit$gp_auto_policy <- policy
  fit
}

gp_hybrid_gp_fit_joint_backend <- function(y_std_vec,
                                           kernels,
                                           trait_cor,
                                           gxe_trait_cor = trait_cor,
                                           X_joint,
                                           train_obs_idx,
                                           all_row,
                                           all_trait,
                                           lambda,
                                           ph = NULL,
                                           female_parent = NULL,
                                           male_parent = NULL,
                                           heter_groups = NULL,
                                           env_kernel_source = NULL,
                                           gp_backend = "auto",
                                           gp_auto_policy = TRUE,
                                           gp_dtype = "float64",
                                           random_state = 123L,
                                           python_bin = NULL,
                                           project_root = NULL) {
  policy <- gp_hybrid_gp_policy(
    ph = ph %||% data.frame(row_id = seq_along(unique(all_row))),
    train_idx = unique(all_row[train_obs_idx]),
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    env_kernel_source = env_kernel_source,
    gp_backend = gp_backend,
    enabled = isTRUE(gp_auto_policy)
  )
  if (gp_hybrid_gp_use_r_backend(gp_backend)) {
    fit <- gp_hybrid_gp_fit_joint_core(
      y_std_vec = y_std_vec,
      kernels = kernels,
      trait_cor = trait_cor,
      gxe_trait_cor = gxe_trait_cor,
      X_joint = X_joint,
      train_obs_idx = train_obs_idx,
      all_row = all_row,
      all_trait = all_trait,
      lambda = lambda
    )
    fit$backend_request <- as.character(gp_backend %||% "r")[1L]
    fit$gp_auto_policy <- policy
    return(fit)
  }
  fit_py <- tryCatch(
    gp_hybrid_gp_python_fit_joint_core(
      y_std_vec = y_std_vec,
      kernels = kernels,
      trait_cor = trait_cor,
      gxe_trait_cor = gxe_trait_cor,
      X_joint = X_joint,
      train_obs_idx = train_obs_idx,
      all_row = all_row,
      all_trait = all_trait,
      lambda = lambda,
      policy = policy,
      gp_backend = gp_backend,
      gp_dtype = gp_dtype,
      random_state = random_state,
      python_bin = python_bin,
      project_root = project_root
    ),
    error = function(e) e
  )
  if (!inherits(fit_py, "error")) {
    return(fit_py)
  }
  if (gp_hybrid_gp_require_python_backend(gp_backend)) {
    stop("Joint multi-trait hybrid GP Python/Torch backend failed: ", conditionMessage(fit_py), call. = FALSE)
  }
  fit <- gp_hybrid_gp_fit_joint_core(
    y_std_vec = y_std_vec,
    kernels = kernels,
    trait_cor = trait_cor,
    gxe_trait_cor = gxe_trait_cor,
    X_joint = X_joint,
    train_obs_idx = train_obs_idx,
    all_row = all_row,
    all_trait = all_trait,
    lambda = lambda
  )
  fit$backend_request <- as.character(gp_backend %||% "auto")[1L]
  fit$backend_fallback_reason <- conditionMessage(fit_py)
  fit$gp_auto_policy <- policy
  fit
}

gp_hybrid_gp_select_lambda_joint <- function(y_std_vec,
                                             kernels,
                                             trait_cor,
                                             gxe_trait_cor = trait_cor,
                                             X_joint,
                                             train_obs_idx,
                                             all_row,
                                             all_trait,
                                             lambda = "auto",
                                             lambda_grid = NULL,
                                             random_state = 123L,
                                             tuning_max_train = NULL) {
  if (!is.character(lambda) || !identical(tolower(lambda[1]), "auto")) {
    lam <- suppressWarnings(as.numeric(lambda)[1L])
    if (is.finite(lam) && lam > 0) {
      return(list(lambda = lam, status = "fixed", detail = data.frame()))
    }
  }
  grid <- gp_hybrid_gp_lambda_grid(lambda_grid)
  tune_train_idx <- as.integer(train_obs_idx)
  max_train <- if (exists("gp_policy_tuning_limit", mode = "function")) {
    gp_policy_tuning_limit(tuning_max_train, "PREDICTPRO_HYBRID_GP_TUNING_MAX_ROWS", 4000L)
  } else {
    suppressWarnings(as.integer(tuning_max_train %||% 4000L))
  }
  if (!is.finite(max_train) || is.na(max_train) || max_train < 1L) {
    max_train <- 4000L
  }
  bounded_tuning <- length(tune_train_idx) > max_train
  if (isTRUE(bounded_tuning)) {
  gp_set_seed(as.integer(random_state))
    tune_train_idx <- sort(sample(tune_train_idx, max_train))
  }
  folds <- gp_hybrid_gp_make_folds(tune_train_idx, nfolds = 3L, seed = random_state)
  if (length(folds) < 2L) {
    chosen <- grid[ceiling(length(grid) / 2L)]
    return(list(lambda = chosen, status = "auto_skipped_small_training_set", detail = data.frame()))
  }
  detail <- do.call(rbind, lapply(grid, function(lam) {
    fold_rmse <- vapply(folds, function(val_idx) {
      cal_idx <- setdiff(train_obs_idx, val_idx)
      if (length(cal_idx) < 2L) {
        return(NA_real_)
      }
      fit <- gp_hybrid_gp_fit_joint_core(
        y_std_vec = y_std_vec,
        kernels = kernels,
        trait_cor = trait_cor,
        gxe_trait_cor = gxe_trait_cor,
        X_joint = X_joint,
        train_obs_idx = cal_idx,
        all_row = all_row,
        all_trait = all_trait,
        lambda = lam
      )
      err <- y_std_vec[val_idx] - fit$prediction_std[val_idx]
      sqrt(mean(err^2, na.rm = TRUE))
    }, numeric(1))
    data.frame(
      lambda = lam,
      mean_rmse = mean(fold_rmse, na.rm = TRUE),
      n_folds = sum(is.finite(fold_rmse)),
      stringsAsFactors = FALSE
    )
  }))
  detail <- detail[is.finite(detail$mean_rmse), , drop = FALSE]
  chosen <- if (nrow(detail)) {
    detail$lambda[which.min(detail$mean_rmse)]
  } else {
    grid[ceiling(length(grid) / 2L)]
  }
  list(
    lambda = chosen,
    status = if (isTRUE(bounded_tuning)) "auto_cv_bounded_training_subset" else "auto_cv",
    detail = detail
  )
}

gp_hybrid_gp_joint_trait_tables <- function(trait_cov,
                                            trait_cor,
                                            responses,
                                            gxe_trait_cov = NULL,
                                            gxe_trait_cor = NULL) {
  pairs <- expand.grid(
    Trait1 = responses,
    Trait2 = responses,
    stringsAsFactors = FALSE
  )
  idx1 <- match(pairs$Trait1, responses)
  idx2 <- match(pairs$Trait2, responses)
  out <- data.frame(
    Trait1 = pairs$Trait1,
    Trait2 = pairs$Trait2,
    Genetic_covariance = trait_cov[cbind(idx1, idx2)],
    Genetic_correlation = trait_cor[cbind(idx1, idx2)],
    stringsAsFactors = FALSE
  )
  if (!is.null(gxe_trait_cov) && !is.null(gxe_trait_cor)) {
    out$GxE_covariance <- gxe_trait_cov[cbind(idx1, idx2)]
    out$GxE_correlation <- gxe_trait_cor[cbind(idx1, idx2)]
  }
  out
}

gp_hybrid_gp_joint_prediction_frame <- function(ph,
                                                responses,
                                                gen_name,
                                                female_parent,
                                                male_parent,
                                                heter_groups,
                                                fit,
                                                y_mean,
                                                y_sd,
                                                all_row,
                                                all_trait,
                                                model_type,
                                                decomposition_basis = "joint_cross_trait_hybrid_gp_additive_kernel") {
  row_order <- all_row
  trait_labels <- responses[all_trait]
  out_cols <- c(gen_name, female_parent, male_parent, "HybridCross")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    out_cols <- c(out_cols, heter_groups)
  }
  out <- ph[row_order, out_cols, drop = FALSE]
  out$Trait <- trait_labels
  out$Observed_value <- vapply(
    seq_along(row_order),
    function(i) as.numeric(ph[[trait_labels[[i]]]][row_order[[i]]]),
    numeric(1)
  )
  out$Train_Test_Label <- ifelse(is.na(out$Observed_value), "Test", "Train")
  scale <- y_sd[all_trait]
  center <- y_mean[all_trait]
  out$Prediction_intercept <- center + scale * as.numeric(fit$fixed_std)
  out$Female_GCA <- scale * as.numeric(fit$female_std)
  out$Male_GCA <- scale * as.numeric(fit$male_std)
  out$SCA_effect <- scale * as.numeric(fit$sca_std)
  out$GxE_effect <- scale * as.numeric(fit$gxe_std %||% rep(0, nrow(out)))
  out$Female_additive_contribution <- out$Female_GCA
  out$Male_additive_contribution <- out$Male_GCA
  out$Hybrid_interaction_contribution <- out$SCA_effect + out$GxE_effect
  out$Environment_interaction_contribution <- out$GxE_effect
  out$Female_GCA_like_contribution <- out$Female_GCA
  out$Male_GCA_like_contribution <- out$Male_GCA
  out$SCA_like_contribution <- out$SCA_effect
  out$GxE_like_contribution <- out$GxE_effect
  out$Predictive_decomposition_basis <- decomposition_basis
  out$Predictive_decomposition_remarks <- paste(
    "Joint hybrid GP decomposition:",
    "cross-trait covariance x additive female GCA, male GCA, SCA, and optional hybrid-by-environment kernels = predicted value"
  )
  out$Predicted_value <- center + scale * as.numeric(fit$prediction_std)
  out$Standard_error <- sqrt(pmax(as.numeric(fit$observed_var_std) * scale^2, 0))
  out$Prediction_error_variance <- pmax(as.numeric(fit$observed_var_std) * scale^2, 0)
  out$Prediction_SE_latent <- sqrt(pmax(as.numeric(fit$latent_var_std) * scale^2, 0))
  out$Prediction_Var_latent <- pmax(as.numeric(fit$latent_var_std) * scale^2, 0)
  out$hybrid_gp_model <- model_type
  out$hybrid_gp_lambda <- fit$lambda
  out$hybrid_gp_backend <- as.character(fit$backend_runtime %||% "R")[1L]
  out$hybrid_gp_device <- as.character(fit$backend_device %||% "cpu")[1L]
  first_cols <- c("Trait", gen_name, female_parent, male_parent, "HybridCross")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    first_cols <- c(first_cols, heter_groups)
  }
  out[, c(first_cols, setdiff(names(out), first_cols)), drop = FALSE]
}

gp_hybrid_gp_joint_variance_components <- function(pred_df,
                                                   fit,
                                                   lambda_info,
                                                   trait_cov,
                                                   trait_cor,
                                                   gxe_trait_cov = NULL,
                                                   gxe_trait_cor = NULL) {
  component_tab <- do.call(rbind, lapply(unique(pred_df$Trait), function(trait) {
    sub <- pred_df[pred_df$Trait == trait, , drop = FALSE]
    component_var <- c(
      female_gca = stats::var(sub$Female_GCA, na.rm = TRUE),
      male_gca = stats::var(sub$Male_GCA, na.rm = TRUE),
      sca = stats::var(sub$SCA_effect, na.rm = TRUE),
      gxe = if ("GxE_effect" %in% names(sub)) stats::var(sub$GxE_effect, na.rm = TRUE) else NA_real_
    )
    data.frame(
      Trait = trait,
      component = c(names(component_var), "residual_lambda", "selected_lambda"),
      estimate = c(component_var, fit$lambda * stats::var(sub$Observed_value, na.rm = TRUE), fit$lambda),
      scale = c(rep("prediction_component_variance", length(component_var)), "response_scale", "standardized_kernel_ratio"),
      tuning_status = c(rep(NA_character_, length(component_var)), lambda_info$status, lambda_info$status),
      stringsAsFactors = FALSE
    )
  }))
  cov_tab <- gp_hybrid_gp_joint_trait_tables(
    trait_cov = trait_cov,
    trait_cor = trait_cor,
    responses = rownames(trait_cor),
    gxe_trait_cov = gxe_trait_cov,
    gxe_trait_cor = gxe_trait_cor
  )
  cov_tab$component <- "cross_trait_genetic_covariance"
  list(
    variance_components = component_tab,
    trait_covariance = cov_tab
  )
}

gp_hybrid_gp_multi_match_observed_values <- function(pred_df,
                                                     ph,
                                                     responses,
                                                     gen_name,
                                                     heter_groups = NULL) {
  if (is.null(pred_df) || !nrow(pred_df)) {
    return(pred_df)
  }
  pred_key <- gp_hybrid_row_key(pred_df, gen_name = gen_name, heter_groups = heter_groups)
  ph_key <- gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups)
  for (resp in responses) {
    idx <- pred_df$Trait == resp
    pred_df$Observed_value[idx] <- ph[[resp]][match(pred_key[idx], ph_key)]
  }
  pred_df
}

gp_hybrid_gp_multi_response_cv_metric_table <- function(pred_df,
                                                        eval_metrics,
                                                        response_family = "gaussian") {
  if (!nrow(pred_df)) {
    return(data.frame())
  }
  split_keys <- unique(pred_df[, c("Trait", "cv_scenario", "fold", "rep"), drop = FALSE])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    idx <- pred_df$Trait == key$Trait &
      pred_df$cv_scenario == key$cv_scenario &
      pred_df$fold == key$fold &
      pred_df$rep == key$rep
    sub <- pred_df[idx, , drop = FALSE]
    row <- data.frame(
      Trait = key$Trait,
      cv_scenario = key$cv_scenario,
      fold = key$fold,
      rep = key$rep,
      stringsAsFactors = FALSE
    )
    for (m in eval_metrics) {
      row[[m]] <- safe_metric_value(
        y_true = sub$Observed_value,
        y_pred = sub$Predicted_value,
        metric = m,
        response_family = response_family
      )
    }
    row
  })
  do.call(rbind, out)
}

gp_hybrid_gp_multi_response_diagnostic_plot <- function(pred_df, model_label) {
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
    ggplot2::facet_wrap(~ Trait, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Hybrid Gaussian-process", model_label, "multi-trait observed vs predicted"),
      subtitle = "Joint cross-trait hybrid GP kernels separate female GCA, male GCA, and SCA",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_hybrid_gp_multi_response_cv_plot <- function(pred_df, model_label) {
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
    ggplot2::facet_grid(Trait ~ cv_scenario, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Hybrid GP", model_label, "multi-trait cross-validation"),
      subtitle = "Joint cross-trait hybrid GP scenarios: known parents, one new parent, both parents new",
      x = "Observed value",
      y = "Predicted value",
      color = "Scenario"
    )
}

gp_hybrid_gp_multi_response_summary <- function(predicted_object,
                                                response,
                                                female_parent,
                                                male_parent,
                                                model_type) {
  pred_df <- as.data.frame(predicted_object, stringsAsFactors = FALSE)
  responses <- unique(as.character(response))
  if (!nrow(pred_df)) {
    return(list(
      summary_statistics = data.frame(stat = character(), summary = character(), Trait = character(), stringsAsFactors = FALSE),
      hybrid_component_summary = data.frame(),
      hybrid_train_test_summary = data.frame()
    ))
  }
  per_trait <- lapply(responses, function(resp) {
    sub <- pred_df[pred_df$Trait == resp, , drop = FALSE]
    out <- summary_statistics_hybrid(
      predicted_object = sub,
      response = resp,
      female_parent = female_parent,
      male_parent = male_parent,
      mode = "hybrid_gp",
      model_type = model_type
    )
    lapply(out, function(tab) gp_hybrid_gp_add_trait_column(tab, resp))
  })
  bind_field <- function(field) {
    tabs <- Filter(
      function(x) !is.null(x) && (nrow(x) > 0L || ncol(x) > 0L),
      lapply(per_trait, `[[`, field)
    )
    if (!length(tabs)) {
      return(data.frame())
    }
    out <- do.call(rbind, tabs)
    rownames(out) <- NULL
    out
  }
  summary_tab <- bind_field("summary_statistics")
  summary_tab <- rbind(
    data.frame(
      Trait = "all",
      stat = c("Multi_Response_Mode", "Trait_Count", "Traits"),
      summary = c("joint_cross_trait_hybrid_gp", length(responses), paste(responses, collapse = ";")),
      stringsAsFactors = FALSE
    ),
    summary_tab
  )
  list(
    summary_statistics = summary_tab,
    hybrid_component_summary = bind_field("hybrid_component_summary"),
    hybrid_train_test_summary = bind_field("hybrid_train_test_summary")
  )
}

gp_hybrid_gp_multi_response_model <- function(pheno_object,
                                              response,
                                              gen_name,
                                              female_parent,
                                              male_parent,
                                              heter_groups = NULL,
                                              model_type = "GP",
                                              gmatrix = NULL,
                                              female_gmatrix = NULL,
                                              male_gmatrix = NULL,
                                              geno_data = NULL,
                                              female_geno_data = NULL,
                                              male_geno_data = NULL,
                                              gmatrix_method = NULL,
                                              include_sca = TRUE,
                                              lambda = "auto",
                                              lambda_grid = NULL,
                                              component_weights = NULL,
                                              gp_backend = "auto",
                                              gp_auto_policy = TRUE,
                                              gp_dtype = "float64",
                                              tuning_max_train = NULL,
                                              env_similarity = NULL,
                                              env_ids = NULL,
                                              env_covariates = NULL,
                                              reaction_norm_feature_qc = TRUE,
                                              kenv_kernel = "matern32",
                                              kenv_bandwidth = 1.0,
                                              kenv_kernel_kwargs = NULL,
                                              random_state = 123L,
                                              python_bin = NULL,
                                              project_root = NULL) {
  responses <- unique(as.character(response))
  responses <- responses[nzchar(responses)]
  if (!length(responses)) {
    stop("Hybrid GP multi-response prediction requires at least one response column.", call. = FALSE)
  }
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  needed_cols <- c(gen_name, female_parent, male_parent, responses)
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    needed_cols <- c(needed_cols, heter_groups)
  }
  missing_cols <- setdiff(needed_cols, names(ph))
  if (length(missing_cols)) {
    stop(paste("Hybrid GP phenotype data is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (is.null(heter_groups) || !nzchar(heter_groups)) {
    if (anyDuplicated(ph[[gen_name]])) {
      stop("Hybrid GP currently requires one phenotype row per hybrid ID unless heter_groups is supplied for multi-environment data.", call. = FALSE)
    }
  } else if (anyDuplicated(gp_hybrid_row_key(ph, gen_name = gen_name, heter_groups = heter_groups))) {
    stop("Hybrid GP currently requires one phenotype row per hybrid-by-environment combination.", call. = FALSE)
  }

  ph[[female_parent]] <- as.character(ph[[female_parent]])
  ph[[male_parent]] <- as.character(ph[[male_parent]])
  ph[[gen_name]] <- as.character(ph[[gen_name]])
  ph$HybridCross <- paste(ph[[female_parent]], ph[[male_parent]], sep = "__x__")
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    ph[[heter_groups]] <- as.character(ph[[heter_groups]])
  }
  env_kernel <- gp_hybrid_gp_prepare_env_similarity(
    ph = ph,
    heter_groups = heter_groups,
    env_similarity = env_similarity,
    env_ids = env_ids,
    env_covariates = env_covariates,
    reaction_norm_feature_qc = reaction_norm_feature_qc,
    kenv_kernel = kenv_kernel,
    kenv_bandwidth = kenv_bandwidth,
    kenv_kernel_kwargs = kenv_kernel_kwargs
  )

  Y <- as.matrix(ph[, responses, drop = FALSE])
  storage.mode(Y) <- "double"
  obs_counts <- colSums(is.finite(Y))
  if (any(obs_counts < 2L)) {
    stop("Joint multi-trait hybrid GP requires at least two observed records for every response.", call. = FALSE)
  }
  y_mean <- vapply(seq_along(responses), function(i) mean(Y[, i], na.rm = TRUE), numeric(1))
  y_sd <- vapply(seq_along(responses), function(i) stats::sd(Y[, i], na.rm = TRUE), numeric(1))
  y_sd[!is.finite(y_sd) | y_sd <= 0] <- 1
  names(y_mean) <- names(y_sd) <- responses
  Y_std <- sweep(sweep(Y, 2, y_mean, "-"), 2, y_sd, "/")
  trait_cor <- gp_hybrid_gp_estimate_trait_correlation(Y_std, responses)
  trait_cov <- diag(y_sd, length(y_sd)) %*% trait_cor %*% diag(y_sd, length(y_sd))
  dimnames(trait_cov) <- list(responses, responses)
  gxe_trait_cor <- gp_hybrid_gp_estimate_gxe_trait_correlation(
    Y_std,
    responses = responses,
    env = if (!is.null(heter_groups) && nzchar(heter_groups)) ph[[heter_groups]] else NULL,
    fallback_cor = trait_cor
  )
  gxe_trait_cov <- diag(y_sd, length(y_sd)) %*% gxe_trait_cor %*% diag(y_sd, length(y_sd))
  dimnames(gxe_trait_cov) <- list(responses, responses)

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
  kernels <- gp_hybrid_gp_component_kernels(
    ph = ph,
    female_parent = female_parent,
    male_parent = male_parent,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    include_sca = include_sca,
    component_weights = component_weights,
    heter_groups = heter_groups,
    env_similarity = env_kernel$K_env
  )
  X <- gp_hybrid_gp_fixed_design(
    ph,
    heter_groups = heter_groups,
    train_idx = which(rowSums(is.finite(Y)) > 0L)
  )
  n <- nrow(ph)
  n_traits <- length(responses)
  all_row <- rep(seq_len(n), times = n_traits)
  all_trait <- rep(seq_len(n_traits), each = n)
  y_std_vec <- as.vector(Y_std)
  train_obs_idx <- which(is.finite(y_std_vec))
  X_joint <- gp_hybrid_gp_joint_design(X, n_traits = n_traits, responses = responses)
  lambda_info <- gp_hybrid_gp_select_lambda_joint(
    y_std_vec = y_std_vec,
    kernels = kernels,
    trait_cor = trait_cor,
    gxe_trait_cor = gxe_trait_cor,
    X_joint = X_joint,
    train_obs_idx = train_obs_idx,
    all_row = all_row,
    all_trait = all_trait,
    lambda = lambda,
    lambda_grid = lambda_grid,
    random_state = random_state,
    tuning_max_train = tuning_max_train
  )
  fit <- gp_hybrid_gp_fit_joint_backend(
    y_std_vec = y_std_vec,
    kernels = kernels,
    trait_cor = trait_cor,
    gxe_trait_cor = gxe_trait_cor,
    X_joint = X_joint,
    train_obs_idx = train_obs_idx,
    all_row = all_row,
    all_trait = all_trait,
    lambda = lambda_info$lambda,
    ph = ph,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    env_kernel_source = env_kernel$source,
    gp_backend = gp_backend,
    gp_auto_policy = gp_auto_policy,
    gp_dtype = gp_dtype,
    random_state = random_state,
    python_bin = python_bin,
    project_root = project_root
  )
  model_label <- gp_hybrid_gp_model_label(model_type)
  pred_df <- gp_hybrid_gp_joint_prediction_frame(
    ph = ph,
    responses = responses,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    fit = fit,
    y_mean = y_mean,
    y_sd = y_sd,
    all_row = all_row,
    all_trait = all_trait,
    model_type = model_label
  )
  female_eff <- stats::aggregate(
    pred_df$Female_GCA,
    by = list(Trait = pred_df$Trait, parent = pred_df[[female_parent]]),
    FUN = mean,
    na.rm = TRUE
  )
  names(female_eff) <- c("Trait", female_parent, "Female_GCA")
  male_eff <- stats::aggregate(
    pred_df$Male_GCA,
    by = list(Trait = pred_df$Trait, parent = pred_df[[male_parent]]),
    FUN = mean,
    na.rm = TRUE
  )
  names(male_eff) <- c("Trait", male_parent, "Male_GCA")
  sca_eff <- stats::aggregate(
    pred_df$SCA_effect,
    by = list(Trait = pred_df$Trait, HybridCross = pred_df$HybridCross),
    FUN = mean,
    na.rm = TRUE
  )
  names(sca_eff) <- c("Trait", "HybridCross", "SCA_effect")
  vc_tables <- gp_hybrid_gp_joint_variance_components(
    pred_df = pred_df,
    fit = fit,
    lambda_info = lambda_info,
    trait_cov = trait_cov,
    trait_cor = trait_cor,
    gxe_trait_cov = gxe_trait_cov,
    gxe_trait_cor = gxe_trait_cor
  )
  lambda_tuning <- lambda_info$detail
  if (nrow(lambda_tuning)) {
    lambda_tuning$Trait <- "joint"
    lambda_tuning <- lambda_tuning[, c("Trait", setdiff(names(lambda_tuning), "Trait")), drop = FALSE]
  }
  policy <- fit$gp_auto_policy %||% gp_hybrid_gp_policy(
    ph = ph,
    train_idx = unique(all_row[train_obs_idx]),
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    env_kernel_source = env_kernel$source,
    gp_backend = gp_backend,
    enabled = isTRUE(gp_auto_policy)
  )

  list(
    gp_model = list(
      model_type = model_label,
      trait_mode = "joint_cross_trait_hybrid_gp",
      traits = responses,
      lambda = fit$lambda,
      lambda_status = lambda_info$status,
      trait_covariance = trait_cov,
      trait_correlation = trait_cor,
      gxe_trait_covariance = gxe_trait_cov,
      gxe_trait_correlation = gxe_trait_cor,
      trait_correlation_table = vc_tables$trait_covariance,
      beta = fit$beta,
      alpha = fit$alpha,
      train_obs_idx = fit$train_obs_idx,
      backend_runtime = fit$backend_runtime %||% "R",
      backend_request = fit$backend_request %||% gp_backend,
      backend_device = fit$backend_device %||% "cpu",
      backend_fallback_reason = fit$backend_fallback_reason %||% NA_character_,
      env_kernel_source = env_kernel$source,
      env_kernel_detail = env_kernel$detail,
      gp_auto_policy = policy
    ),
    predicted_values = pred_df,
    Predicted_value = pred_df,
    female_gca_effects = female_eff,
    male_gca_effects = male_eff,
    sca_effects = sca_eff,
    female_gmatrix = female_g,
    male_gmatrix = male_g,
    sca_gmatrix = sca_g,
    env_similarity = env_kernel$K_env,
    env_kernel = env_kernel,
    component_weights = kernels$weights,
    lambda_tuning = lambda_tuning,
    variance_components = vc_tables$variance_components,
    trait_covariance = vc_tables$trait_covariance,
    trait_correlation = trait_cor,
    gxe_trait_covariance = gxe_trait_cov,
    gxe_trait_correlation = gxe_trait_cor,
    diagnostic_plots = gp_hybrid_gp_multi_response_diagnostic_plot(pred_df, model_label = model_label)
  )
}

gp_hybrid_gp_multi_response_cv <- function(pheno_object,
                                           response,
                                           gen_name,
                                           female_parent,
                                           male_parent,
                                           heter_groups = NULL,
                                           cross_validation_meth,
                                           eval_metrics,
                                           model_type = "GP",
                                           gmatrix = NULL,
                                           female_gmatrix = NULL,
                                           male_gmatrix = NULL,
                                           geno_data = NULL,
                                           female_geno_data = NULL,
                                           male_geno_data = NULL,
                                           gmatrix_method = NULL,
                                           include_sca = TRUE,
                                           lambda = "auto",
                                           lambda_grid = NULL,
                                           component_weights = NULL,
                                           gp_backend = "auto",
                                           gp_auto_policy = TRUE,
                                           gp_dtype = "float64",
                                           tuning_max_train = NULL,
                                           env_similarity = NULL,
                                           env_ids = NULL,
                                           env_covariates = NULL,
                                           reaction_norm_feature_qc = TRUE,
                                           kenv_kernel = "matern32",
                                           kenv_bandwidth = 1.0,
                                           kenv_kernel_kwargs = NULL,
                                           nfolds = 5L,
                                           random_state = 123L,
                                           replication = 1L,
                                           python_bin = NULL,
                                           project_root = NULL) {
  responses <- unique(as.character(response))
  responses <- responses[nzchar(responses)]
  if (!length(responses)) {
    stop("Hybrid GP multi-response CV requires at least one response column.", call. = FALSE)
  }
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  all_preds <- list()
  all_eval <- list()
  all_raw <- list()

  for (rep_i in seq_len(replication)) {
    scenarios <- gp_hybrid_build_cv_scenarios(
      pheno_object = ph,
      response = responses[[1L]],
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
      ph_cv <- gp_hybrid_cv_mask_responses(ph, responses, scenario_rows$train)

      fit <- gp_hybrid_gp_multi_response_model(
        pheno_object = ph_cv,
        response = responses,
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
        lambda = lambda,
        lambda_grid = lambda_grid,
        component_weights = component_weights,
        gp_backend = gp_backend,
        gp_auto_policy = gp_auto_policy,
        gp_dtype = gp_dtype,
        tuning_max_train = tuning_max_train,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        random_state = random_state + rep_i,
        python_bin = python_bin,
        project_root = project_root
      )

      pred_all <- fit$predicted_values
      test_keys <- gp_hybrid_row_key(
        ph[tst, , drop = FALSE],
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      prediction_keys <- gp_hybrid_row_key(
        pred_all,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      pred_test <- pred_all[prediction_keys %in% test_keys, , drop = FALSE]
      pred_test <- gp_hybrid_gp_multi_match_observed_values(
        pred_df = pred_test,
        ph = ph,
        responses = responses,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      pred_test <- pred_test[is.finite(pred_test$Observed_value), , drop = FALSE]
      pred_test$cv_scenario <- sc$scenario
      pred_test$fold <- sc$fold
      pred_test$rep <- rep_i

      eval_df <- gp_hybrid_gp_multi_response_cv_metric_table(
        pred_df = pred_test,
        eval_metrics = eval_metrics,
        response_family = "gaussian"
      )

      all_preds[[length(all_preds) + 1L]] <- pred_test
      all_eval[[length(all_eval) + 1L]] <- eval_df
      all_raw[[length(all_raw) + 1L]] <- list(
        trait = paste(responses, collapse = ";"),
        rep = rep_i,
        model = gp_hybrid_gp_model_label(model_type),
        eval_metrics_reps = eval_df,
        ypred_cv_Reps_all = pred_test,
        trait_covariance = fit$trait_covariance,
        cv_info = list(
          method = cross_validation_meth,
          scenario = sc$scenario,
          fold = sc$fold,
          n_test = nrow(pred_test),
          n_train = sum(vapply(responses, function(resp) sum(!is.na(ph_cv[[resp]])), numeric(1))),
          n_test_rows = nrow(pred_test)
        )
      )
    }
  }

  pred_df <- if (length(all_preds)) do.call(rbind, all_preds) else data.frame()
  metric_detail <- if (length(all_eval)) do.call(rbind, all_eval) else data.frame()
  metric_summary <- if (is.null(metric_detail) || !nrow(metric_detail)) {
    data.frame()
  } else {
    numeric_cols <- names(metric_detail)[vapply(metric_detail, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, "rep")
    stats::aggregate(
      metric_detail[, numeric_cols, drop = FALSE],
      by = list(Trait = metric_detail$Trait, cv_scenario = metric_detail$cv_scenario),
      FUN = mean,
      na.rm = TRUE
    )
  }
  if (nrow(metric_summary)) {
    metric_summary$model <- gp_hybrid_gp_model_label(model_type)
    metric_summary <- metric_summary[, c(
      "Trait", "cv_scenario", "model",
      setdiff(names(metric_summary), c("Trait", "cv_scenario", "model"))
    ), drop = FALSE]
  }
  counts <- if (!nrow(pred_df)) {
    data.frame()
  } else {
    stats::aggregate(
      pred_df[[gen_name]],
      by = list(Trait = pred_df$Trait, cv_scenario = pred_df$cv_scenario),
      FUN = length
    )
  }
  if (nrow(counts)) {
    names(counts)[names(counts) == "x"] <- "n_test_hybrids"
  }
  processed <- list(
    hybrid_metric_summary = metric_summary,
    hybrid_prediction_counts = counts,
    hybrid_cv_predictions = pred_df,
    hybrid_cv_plot = gp_hybrid_gp_multi_response_cv_plot(pred_df, model_label = gp_hybrid_gp_model_label(model_type))
  )
  list(
    predicted_values = pred_df,
    hybrid_cv_metrics = metric_detail,
    cv_results_processed = processed,
    cv_results_raw = all_raw,
    diagnostic_plots = processed$hybrid_cv_plot
  )
}
