gp_lowrank_supported_models <- function() c("KRR", "GP", "GP_FA", "LowRankGP")

gp_backend_resolve_model <- function(model_name = NULL) {
  model_name <- as.character(model_name %||% "KRR")[1]
  if (identical(model_name, "LowRankGP")) {
    return("KRR")
  }
  model_name
}

gp_backend_method_map <- function(model_name = NULL) {
  switch(
    gp_backend_resolve_model(model_name),
    KRR = "krr_exact",
    GP = "gp_exact",
    GP_FA = "gp_icm_fa",
    stop("Unsupported GP backend model.", call. = FALSE)
  )
}

gp_match_gp_choice <- function(value, choices, name, default = choices[[1L]]) {
  value <- value %||% default
  if (!length(value)) {
    value <- default
  }
  value <- tolower(as.character(value)[1L])
  if (is.na(value) || !nzchar(value) || !value %in% choices) {
    stop(
      "`", name, "` must be one of: ",
      paste(sprintf('"%s"', choices), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  value
}

gp_backend_schema_cols <- function() {
  list(
    single = list(
      gid = "GID",
      env = "Env",
      y = "y",
      value = "Trait"
    ),
    long = list(
      gid = "gid",
      env = "env",
      trait = "trait",
      y = "y"
    ),
    kernels = list(
      genomic = "G",
      additive = "A"
    )
  )
}

gp_backend_schema_col <- function(section, name) {
  cols <- gp_backend_schema_cols()
  if (!section %in% names(cols) || !name %in% names(cols[[section]])) {
    stop("Unknown GP backend schema column: ", section, ".", name, call. = FALSE)
  }
  cols[[section]][[name]]
}

gp_first_existing_col <- function(data, candidates) {
  candidates <- unlist(candidates, use.names = FALSE)
  candidates <- as.character(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates <- unique(candidates)
  hits <- candidates[candidates %in% names(data)]
  if (length(hits)) hits[[1L]] else NULL
}

gp_reticulate_available <- function() {
  requireNamespace("reticulate", quietly = TRUE)
}

gp_require_reticulate <- function(context = "This GP backend path") {
  if (!gp_reticulate_available()) {
    stop(
      context,
      " requires the reticulate R package. The production default GP path ",
      "uses the direct Python bridge; keep PREDICTPRO_GP_EXECUTION=direct ",
      "and leave PREDICTPRO_GP_RETICULATE_FALLBACK disabled, or install ",
      "reticulate before enabling the reticulate fallback path.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

gp_warn_large_met_without_env_kernel <- function(env_levels,
                                                 has_env_kernel,
                                                 context = "MET GP") {
  n_env <- length(unique(as.character(env_levels)))
  if (n_env > 5L && !isTRUE(has_env_kernel)) {
    warning(
      context,
      " has ", n_env, " environments and no supplied environment kernel. ",
      "For large MET, provide `env_similarity` or `env_covariates`; a ",
      "covariate-derived environment kernel replaces environment FA or ",
      "unstructured environment covariance.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

gp_policy_is_auto <- function(x) {
  is.character(x) && length(x) == 1L && identical(tolower(x), "auto")
}

gp_policy_bool_or_auto <- function(x, name, default = FALSE) {
  if (gp_policy_is_auto(x)) {
    return("auto")
  }
  if (is.null(x) || !length(x)) {
    return(isTRUE(default))
  }
  if (is.logical(x)) {
    return(isTRUE(x[[1L]]))
  }
  if (is.character(x)) {
    value <- tolower(trimws(x[[1L]]))
    if (value %in% c("true", "t", "yes", "y", "1")) {
      return(TRUE)
    }
    if (value %in% c("false", "f", "no", "n", "0")) {
      return(FALSE)
    }
  }
  stop("`", name, "` must be TRUE, FALSE, or \"auto\".", call. = FALSE)
}

gp_policy_scalar_default <- function(value, default) {
  if (length(value) != 1L) {
    return(FALSE)
  }
  if (is.character(default)) {
    return(identical(as.character(value)[1L], default))
  }
  isTRUE(all.equal(as.numeric(value)[1L], as.numeric(default)[1L], tolerance = 1e-12))
}

gp_policy_tuning_limit <- function(value, envvar, default) {
  val <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(val) || is.na(val)) {
    val <- suppressWarnings(as.numeric(Sys.getenv(envvar, unset = as.character(default)))[1L])
  }
  if (!is.finite(val) || is.na(val)) {
    val <- as.numeric(default)[1L]
  }
  as.integer(max(1, floor(val)))
}

gp_policy_device_route <- local({
  cached <- NULL

  function() {
    requested <- Sys.getenv("PREDICTPRO_GP_DEVICE", unset = "")
    if (nzchar(requested)) {
      return(list(
        route = NA_character_,
        reason = "user_device_env_respected",
        requested = requested,
        num_gpus = NA_integer_,
        gpu_usage = NA_real_,
        gpu_busy = NA
      ))
    }

    ttl <- suppressWarnings(as.numeric(Sys.getenv("PREDICTPRO_GP_DEVICE_ROUTE_CACHE_SEC", unset = "2"))[1L])
    if (!is.finite(ttl) || is.na(ttl)) {
      ttl <- 2
    }
    now <- as.numeric(Sys.time())
    if (is.list(cached) && ttl > 0 && is.finite(cached$time) && (now - cached$time) <= ttl) {
      return(cached$value)
    }

    num_gpus <- if (exists("detect_physical_num_gpus", mode = "function")) {
      tryCatch(detect_physical_num_gpus(), error = function(e) NA_integer_)
    } else {
      NA_integer_
    }
    if (!is.finite(num_gpus) || is.na(num_gpus)) {
      num_gpus <- 0L
    }
    usage <- if (exists("get_gpu_usage", mode = "function")) {
      tryCatch(get_gpu_usage(), error = function(e) NA_real_)
    } else {
      NA_real_
    }
    busy <- if (exists("gp_gpu_is_busy", mode = "function")) {
      tryCatch(gp_gpu_is_busy(usage), error = function(e) FALSE)
    } else {
      FALSE
    }
    if (as.integer(num_gpus) <= 0L) {
      route <- "cpu"
      reason <- "gp_cpu_no_gpu"
    } else if (isTRUE(busy)) {
      route <- "cpu"
      reason <- "gp_cpu_gpu_busy"
    } else {
      route <- "auto"
      reason <- "gp_gpu_auto_available"
    }
    value <- list(
      route = route,
      reason = reason,
      requested = "auto",
      num_gpus = as.integer(num_gpus),
      gpu_usage = usage,
      gpu_busy = isTRUE(busy)
    )
    if (ttl > 0) {
      cached <<- list(time = now, value = value)
    }
    value
  }
})

gp_with_gp_device_route <- function(route, expr) {
  route <- as.character(route %||% NA_character_)[1L]
  if (is.na(route) || !nzchar(route)) {
    return(force(expr))
  }
  old <- Sys.getenv("PREDICTPRO_GP_DEVICE", unset = NA_character_)
  Sys.setenv(PREDICTPRO_GP_DEVICE = route)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("PREDICTPRO_GP_DEVICE")
    } else {
      Sys.setenv(PREDICTPRO_GP_DEVICE = old)
    }
  }, add = TRUE)
  force(expr)
}

gp_policy_prediction_profile <- function(pheno_df,
                                         split,
                                         heter_groups = NULL,
                                         gid_col = NULL,
                                         env_col = NULL,
                                         has_env_kernel = FALSE,
                                         has_env_covariates = FALSE,
                                         has_factor_cache = FALSE,
                                         large_n_threshold = 50000L) {
  schema <- gp_backend_schema_cols()
  env_name <- gp_first_existing_col(
    pheno_df,
    c(env_col, heter_groups, schema$single$env, schema$long$env)
  )
  gid_name <- gp_first_existing_col(
    pheno_df,
    c(gid_col, schema$single$gid, schema$long$gid)
  )
  env <- if (!is.null(env_name)) {
    as.character(pheno_df[[env_name]])
  } else {
    rep("ENV1", nrow(pheno_df))
  }
  env[is.na(env) | !nzchar(env)] <- "ENV1"
  gids <- if (!is.null(gid_name)) {
    as.character(pheno_df[[gid_name]])
  } else {
    as.character(seq_len(nrow(pheno_df)))
  }
  train_mask <- as.logical(split$train_mask)
  test_mask <- as.logical(split$test_mask)
  if (length(train_mask) != nrow(pheno_df) || length(test_mask) != nrow(pheno_df)) {
    stop("GP policy split masks must match phenotype row count.", call. = FALSE)
  }
  train_mask[is.na(train_mask)] <- FALSE
  test_mask[is.na(test_mask)] <- FALSE
  train_envs <- unique(env[train_mask])
  test_envs <- unique(env[test_mask])
  unseen_test_envs <- setdiff(test_envs, train_envs)
  n_train <- sum(train_mask)
  n_test <- sum(test_mask)
  n_rows <- nrow(pheno_df)
  n_env <- length(unique(env))
  large_n_threshold <- as.integer(large_n_threshold)[1L]
  if (!is.finite(large_n_threshold) || is.na(large_n_threshold)) {
    large_n_threshold <- 50000L
  }
  list(
    n_rows = as.integer(n_rows),
    n_train = as.integer(n_train),
    n_test = as.integer(n_test),
    n_genotypes = as.integer(length(unique(gids))),
    n_env = as.integer(n_env),
    n_train_env = as.integer(length(train_envs)),
    n_test_env = as.integer(length(test_envs)),
    train_envs = train_envs,
    test_envs = test_envs,
    unseen_test_envs = unseen_test_envs,
    n_unseen_test_env = as.integer(length(unseen_test_envs)),
    is_met = !is.null(env_name) && n_env > 1L,
    has_env_kernel = isTRUE(has_env_kernel),
    has_env_covariates = isTRUE(has_env_covariates),
    has_factor_cache = isTRUE(has_factor_cache),
    prospective_env_prediction = !is.null(env_name) && length(unseen_test_envs) > 0L,
    large_fit = isTRUE(has_factor_cache) || n_train >= large_n_threshold || n_rows >= large_n_threshold
  )
}

gp_single_trait_auto_policy <- function(pheno_df,
                                        split,
                                        heter_groups = NULL,
                                        gid_col = NULL,
                                        env_col = NULL,
                                        env_similarity = NULL,
                                        env_covariates = NULL,
                                        gmatrix = NULL,
                                        gp_factor_cache = NULL,
                                        geno_ids = NULL,
                                        gp_backend = "auto",
                                        gp_tiered_dispatch = FALSE,
                                        gp_dtype = "float64",
                                        large_n_threshold = 50000L,
                                        operator_tol = 1e-5,
                                        operator_max_iter = 500L,
                                        operator_dtype_compute = "float32",
                                        tune_gp = FALSE,
                                        tuning_objective = "rmse",
                                        tuning_max_calibration_rows = NULL,
                                        tuning_max_validation_rows = NULL,
                                        tuning_max_genotypes = NULL,
                                        include_components = NULL,
                                        enabled = TRUE,
                                        resolve_device = TRUE) {
  profile <- gp_policy_prediction_profile(
    pheno_df,
    split,
    heter_groups = heter_groups,
    gid_col = gid_col,
    env_col = env_col,
    has_env_kernel = !is.null(env_similarity) || !is.null(env_covariates),
    has_env_covariates = !is.null(env_covariates),
    has_factor_cache = !is.null(gp_factor_cache),
    large_n_threshold = large_n_threshold
  )
  tune_mode <- gp_policy_bool_or_auto(tune_gp, "tune_gp", default = FALSE)
  requested_objective <- as.character(tuning_objective %||% "rmse")[1L]
  device <- if (isTRUE(resolve_device)) {
    gp_policy_device_route()
  } else {
    list(
      route = NA_character_,
      reason = "python_gp_device_policy_delegated",
      requested = "auto",
      num_gpus = NA_integer_,
      gpu_usage = NA_real_,
      gpu_busy = NA
    )
  }
  decisions <- character()
  if (!isTRUE(enabled)) {
    if (gp_policy_is_auto(tune_gp)) {
      tune_mode <- FALSE
    }
    if (gp_policy_is_auto(requested_objective)) {
      requested_objective <- "rmse"
    }
    return(list(
      enabled = FALSE,
      profile = profile,
      device = device,
      decisions = "auto_policy_disabled",
      gp_backend = gp_backend,
      gp_tiered_dispatch = isTRUE(gp_tiered_dispatch),
      gp_dtype = gp_dtype,
      operator_tol = operator_tol,
      operator_max_iter = operator_max_iter,
      operator_dtype_compute = operator_dtype_compute,
      tune_gp = isTRUE(tune_mode),
      tuning_objective = requested_objective,
      include_components = include_components
    ))
  }

  gp_backend_eff <- as.character(gp_backend %||% "auto")[1L]
  if (identical(tolower(gp_backend_eff), "auto") && isTRUE(profile$large_fit)) {
    gp_backend_eff <- "operator"
    decisions <- c(decisions, "large_or_factor_cache_operator_backend")
  }

  tiered_eff <- isTRUE(gp_tiered_dispatch)
  if (!tiered_eff && isTRUE(profile$is_met) && isTRUE(profile$has_env_kernel) &&
      isTRUE(profile$prospective_env_prediction) &&
      (isTRUE(profile$has_env_covariates) || isTRUE(profile$large_fit) || profile$n_env > 5L)) {
    tiered_eff <- TRUE
    decisions <- c(decisions, "prospective_met_tiered_dispatch")
  }

  gp_dtype_eff <- as.character(gp_dtype %||% "float64")[1L]
  if (isTRUE(profile$large_fit) && gp_policy_scalar_default(gp_dtype_eff, "float64")) {
    gp_dtype_eff <- "float32"
    decisions <- c(decisions, "large_fit_float32")
  }

  operator_tol_eff <- as.numeric(operator_tol)[1L]
  if (isTRUE(profile$large_fit) && gp_policy_scalar_default(operator_tol_eff, 1e-5)) {
    operator_tol_eff <- 1e-4
    decisions <- c(decisions, "large_fit_operator_tol")
  }

  operator_max_iter_eff <- as.integer(operator_max_iter)[1L]
  if (isTRUE(profile$large_fit) && gp_policy_scalar_default(operator_max_iter_eff, 500L)) {
    operator_max_iter_eff <- 300L
    decisions <- c(decisions, "large_fit_operator_max_iter")
  }

  operator_dtype_compute_eff <- as.character(operator_dtype_compute %||% "float32")[1L]
  can_dense_tune <- !is.null(gmatrix) && !isTRUE(profile$has_factor_cache)
  can_factor_cache_tune <- !is.null(gp_factor_cache) && length(geno_ids %||% character()) > 0L
  tuning_source <- if (isTRUE(can_dense_tune)) {
    "dense_gmatrix"
  } else if (isTRUE(can_factor_cache_tune)) {
    "factor_cache_subset"
  } else {
    "unavailable"
  }
  factor_limits <- list(
    max_calibration_rows = gp_policy_tuning_limit(
      tuning_max_calibration_rows,
      "PREDICTPRO_GP_TUNING_MAX_CALIBRATION_ROWS",
      4000L
    ),
    max_validation_rows = gp_policy_tuning_limit(
      tuning_max_validation_rows,
      "PREDICTPRO_GP_TUNING_MAX_VALIDATION_ROWS",
      2000L
    ),
    max_genotypes = gp_policy_tuning_limit(
      tuning_max_genotypes,
      "PREDICTPRO_GP_TUNING_MAX_GENOTYPES",
      1500L
    )
  )
  dense_limits <- list(
    max_calibration_rows = tuning_max_calibration_rows,
    max_validation_rows = tuning_max_validation_rows,
    max_genotypes = tuning_max_genotypes
  )
  auto_tune_reason <- NA_character_
  if (identical(tune_mode, "auto")) {
    known_environment_prediction <- isTRUE(profile$is_met) &&
      !isTRUE(profile$prospective_env_prediction)
    met_tune_ok <- isTRUE(profile$is_met) &&
      (isTRUE(profile$has_env_covariates) || known_environment_prediction) &&
      is.null(env_similarity) &&
      (isTRUE(can_dense_tune) || isTRUE(can_factor_cache_tune)) &&
      profile$n_train >= 200L &&
      profile$n_train_env >= 4L
    single_env_tune_ok <- !isTRUE(profile$is_met) &&
      is.null(env_similarity) &&
      (isTRUE(can_dense_tune) || isTRUE(can_factor_cache_tune)) &&
      profile$n_train >= 50L &&
      profile$n_genotypes >= 30L
    tune_mode <- isTRUE(met_tune_ok) || isTRUE(single_env_tune_ok)
    auto_tune_reason <- if (isTRUE(tune_mode)) {
      if (isTRUE(single_env_tune_ok) && isTRUE(can_dense_tune)) {
        "auto_tune_enabled_dense_single_environment_genotype"
      } else if (isTRUE(single_env_tune_ok)) {
        "auto_tune_enabled_factor_cache_subset_single_environment_genotype"
      } else if (isTRUE(can_dense_tune) && isTRUE(profile$has_env_covariates)) {
        "auto_tune_enabled_dense_met_env_covariates"
      } else if (isTRUE(can_dense_tune)) {
        "auto_tune_enabled_dense_met_known_environment"
      } else if (isTRUE(profile$has_env_covariates)) {
        "auto_tune_enabled_factor_cache_subset_met_env_covariates"
      } else {
        "auto_tune_enabled_factor_cache_subset_met_known_environment"
      }
    } else if (!isTRUE(profile$is_met) &&
               (!isTRUE(can_dense_tune) && !isTRUE(can_factor_cache_tune))) {
      "auto_tune_skipped_single_environment_requires_dense_gmatrix_or_factor_cache"
    } else if (!isTRUE(profile$is_met) &&
               (profile$n_train < 50L || profile$n_genotypes < 30L)) {
      "auto_tune_skipped_single_environment_too_few_training_rows_or_genotypes"
    } else if (!isTRUE(profile$is_met)) {
      "auto_tune_skipped_single_environment_policy_guard"
    } else if (isTRUE(profile$prospective_env_prediction) && !isTRUE(profile$has_env_covariates)) {
      "auto_tune_skipped_no_env_covariates"
    } else if (!isTRUE(can_dense_tune) && !isTRUE(can_factor_cache_tune)) {
      "auto_tune_skipped_requires_dense_gmatrix_or_factor_cache"
    } else if (profile$n_train < 200L || profile$n_train_env < 4L) {
      "auto_tune_skipped_too_few_training_rows_or_envs"
    } else {
      "auto_tune_skipped_policy_guard"
    }
    decisions <- c(decisions, auto_tune_reason)
  }

  tuning_objective_eff <- requested_objective
  if (gp_policy_is_auto(tuning_objective_eff)) {
    tuning_objective_eff <- if (isTRUE(profile$prospective_env_prediction) &&
                                isTRUE(profile$has_env_covariates)) {
      "centered_rmse"
    } else {
      "rmse"
    }
    decisions <- c(decisions, paste0("auto_tuning_objective_", tuning_objective_eff))
  }

  include_eff <- include_components
  if (is.null(include_eff)) {
    include_eff <- if (isTRUE(profile$is_met) && isTRUE(profile$has_env_kernel)) {
      c("g", "ge", "e")
    } else {
      "g"
    }
  }

  list(
    enabled = TRUE,
    profile = profile,
    device = device,
    decisions = if (length(decisions)) unique(decisions) else "manual_or_default_settings_kept",
    gp_backend = gp_backend_eff,
    gp_tiered_dispatch = tiered_eff,
    gp_dtype = gp_dtype_eff,
    operator_tol = operator_tol_eff,
    operator_max_iter = operator_max_iter_eff,
    operator_dtype_compute = operator_dtype_compute_eff,
    tune_gp = isTRUE(tune_mode),
    tune_gp_mode = tune_mode,
    tuning_source = tuning_source,
    tuning_limits = if (identical(tuning_source, "factor_cache_subset")) factor_limits else dense_limits,
    auto_tune_reason = auto_tune_reason,
    tuning_objective = tuning_objective_eff,
    include_components = include_eff
  )
}

gp_prepare_env_covariates <- function(env_covariates,
                                      source_env_col = NULL,
                                      target_env_col = gp_backend_schema_col("single", "env")) {
  if (is.null(env_covariates) || !is.data.frame(env_covariates) ||
      is.null(source_env_col) || identical(source_env_col, target_env_col)) {
    return(env_covariates)
  }
  if (!(target_env_col %in% names(env_covariates)) &&
      source_env_col %in% names(env_covariates)) {
    env_covariates <- env_covariates
    names(env_covariates)[match(source_env_col, names(env_covariates))] <- target_env_col
  }
  env_covariates
}

gp_vendored_project_root <- function() {
  installed_root <- system.file("python/gaussian_process", package = "PredictProR")
  if (nzchar(installed_root) && dir.exists(installed_root)) {
    return(normalizePath(installed_root, winslash = "/", mustWork = TRUE))
  }

  source_candidates <- c(
    "inst/python/gaussian_process",
    "python/gaussian_process",
    "../inst/python/gaussian_process"
  )
  hits <- source_candidates[file.exists(file.path(source_candidates,
    "gp_framework.py"
  ))]
  if (!length(hits)) {
    return(NULL)
  }
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

gp_project_root <- function() {
  root <- Sys.getenv("PREDICTPRO_GP_PROJECT_ROOT", unset = "")
  vendored_root <- gp_vendored_project_root()
  candidates <- c(
    root,
    vendored_root
  )
  candidates <- unique(candidates[nzchar(candidates)])
  framework_name <- "gp_framework.py"
  hits <- candidates[file.exists(file.path(candidates, framework_name))]
  if (!length(hits)) {
    stop(
      "Gaussian_process project root was not found. Set PREDICTPRO_GP_PROJECT_ROOT ",
      "to the folder containing ", framework_name, ".",
      call. = FALSE
    )
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

gp_python_candidates_bridge <- function() {
  explicit_gp <- Sys.getenv("PREDICTPRO_GP_PYTHON", unset = "")
  if (nzchar(explicit_gp) && file.exists(explicit_gp)) {
    return(normalizePath(explicit_gp, winslash = "/", mustWork = TRUE))
  }
  root <- tryCatch(gp_project_root(), error = function(e) NULL)
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = FALSE)
  workon_home <- Sys.getenv("WORKON_HOME", unset = file.path(home, ".virtualenvs"))
  local_venv <- if (!is.null(root)) {
    gp_virtualenv_python_path(root, ".venv-gp")
  } else {
    character()
  }
  preferred_gp <- tryCatch(gp_preferred_python(purpose = "gp"), error = function(e) NULL)
  candidates <- c(
    preferred_gp,
    local_venv,
    gp_local_virtualenv_python_candidates(".venv-gp", roots = getwd()),
    gp_virtualenv_python_path(workon_home, "predictgp"),
    gp_conda_python_candidates("predictgp"),
    gp_system_python_candidates()
  )
  candidates <- candidates[nzchar(candidates)]
  unique(candidates[file.exists(candidates)])
}

gp_detect_gp_python <- function() {
  cand <- gp_python_candidates_bridge()
  if (!length(cand)) {
    return(NULL)
  }
  normalizePath(cand[[1]], winslash = "/", mustWork = TRUE)
}

gp_bridge_script_path <- function() {
  pyfile <- system.file("python/gp_bridge.py", package = "PredictProR")
  if (!nzchar(pyfile)) {
    candidates <- c("inst/python/gp_bridge.py", "python/gp_bridge.py", "../inst/python/gp_bridge.py")
    hits <- candidates[file.exists(candidates)]
    if (length(hits)) {
      pyfile <- normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
    }
  }
  if (!nzchar(pyfile)) {
    stop("gp_bridge.py not found under inst/python/.", call. = FALSE)
  }
  pyfile
}

gp_bridge_default_env <- function(pheno_data, heter_groups = NULL) {
  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    return(as.character(pheno_data[[heter_groups]]))
  }
  rep("ENV1", nrow(pheno_data))
}

gp_bridge_build_pheno <- function(pheno_data, response, gen_name, heter_groups = NULL) {
  cols <- gp_backend_schema_cols()$single
  out <- data.frame(
    .gid = as.character(pheno_data[[gen_name]]),
    .env = as.character(gp_bridge_default_env(pheno_data, heter_groups = heter_groups)),
    .value = as.numeric(pheno_data[[response]]),
    stringsAsFactors = FALSE
  )
  names(out) <- unname(c(cols$gid, cols$env, cols$value))
  out
}

gp_bridge_train_test_idx <- function(pheno_df, test_set = NULL, test_set_source = NULL) {
  cols <- gp_backend_schema_cols()$single
  gids <- as.character(pheno_df[[cols$gid]])
  y <- suppressWarnings(as.numeric(pheno_df[[cols$value]]))
  test_mask <- gp_public_test_mask(
    gids = gids,
    y = y,
    test_set = test_set,
    test_set_source = test_set_source
  )
  train_mask <- !test_mask & is.finite(y)
  list(
    train_idx0 = as.integer(which(train_mask) - 1L),
    test_idx0 = as.integer(which(test_mask) - 1L),
    test_mask = test_mask,
    train_mask = train_mask
  )
}

gp_bridge_output_level <- function(return_se = FALSE, full_vc = FALSE) {
  if (isTRUE(full_vc)) {
    return("full_vc")
  }
  if (isTRUE(return_se)) {
    return("predict_with_se")
  }
  "predict_only"
}

gp_bridge_add_gp_exact_controls <- function(args, controls = NULL) {
  controls <- controls %||% list()
  gp_iters <- suppressWarnings(as.integer(controls$gp_iters %||% NA_integer_)[1L])
  if (is.finite(gp_iters) && !is.na(gp_iters) && gp_iters > 0L) {
    args$gp_iters <- gp_iters
  }
  gp_lr <- suppressWarnings(as.numeric(controls$gp_lr %||% NA_real_)[1L])
  if (is.finite(gp_lr) && !is.na(gp_lr) && gp_lr > 0) {
    args$gp_lr <- gp_lr
  }
  args
}

gp_bridge_model_params <- function(model_name,
                                   gp_backend = "auto",
                                   gp_output_level = "predict_only",
                                   gp_varcomp_mode = "reml",
                                   gp_fa_rank = 1L,
                                   gp_prediction_output = "all",
                                   gp_factor_cache = NULL,
                                   gp_iters = NULL,
                                   gp_lr = NULL,
                                   seed = 12345L) {
  model_resolved <- gp_backend_resolve_model(model_name)
  params <- list(
    method = gp_backend_method_map(model_resolved),
    backend = gp_backend,
    output_level = gp_output_level,
    prediction_output = gp_prediction_output,
    varcomp_mode = gp_varcomp_mode %||% "reml",
    fa_rank = as.integer(gp_fa_rank %||% 1L),
    grm_factor_cache_root = if (is.list(gp_factor_cache) && identical(gp_factor_cache$type, "zarr")) gp_factor_cache$root %||% NULL else NULL,
    seed = as.integer(seed %||% 12345L)
  )
  gp_bridge_add_gp_exact_controls(params, list(gp_iters = gp_iters, gp_lr = gp_lr))
}

gp_bridge_write_matrix_bin <- function(x, bin_path, meta_path) {
  x <- unname(as.matrix(x))
  con <- file(bin_path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.double(x), con, size = 8L, endian = "little")
  writeLines(c(as.character(nrow(x)), as.character(ncol(x))), meta_path, useBytes = TRUE)
  invisible(NULL)
}

gp_bridge_read_matrix_bin <- function(bin_path, meta_path) {
  dims <- as.integer(readLines(meta_path, warn = FALSE))
  if (length(dims) < 2L || anyNA(dims[1:2])) {
    stop("Invalid matrix metadata file: ", meta_path, call. = FALSE)
  }
  con <- file(bin_path, open = "rb")
  on.exit(close(con), add = TRUE)
  values <- readBin(con, what = "double", n = prod(dims[1:2]), size = 8L, endian = "little")
  if (length(values) != prod(dims[1:2])) {
    stop("Matrix binary size mismatch: ", bin_path, call. = FALSE)
  }
  matrix(values, nrow = dims[[1L]], ncol = dims[[2L]])
}

gp_bridge_write_lines <- function(x, path) {
  writeLines(as.character(x), con = path, useBytes = TRUE)
  invisible(path)
}

gp_bridge_write_kernel_bank <- function(kernel, input_dir) {
  kernels <- kernel$Ks
  if (is.null(kernels) || !length(kernels)) {
    kernels <- list(G = kernel$K)
  }
  kernel_names <- names(kernels)
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("K", seq_along(kernels))
  }
  kernel_names <- make.unique(gp_sanitize_kernel_name(kernel_names), sep = "_")
  specs <- vector("list", length(kernels))
  names(specs) <- kernel_names
  for (i in seq_along(kernels)) {
    base <- paste0("kernel_", kernel_names[[i]])
    bin_path <- file.path(input_dir, paste0(base, ".bin"))
    meta_path <- file.path(input_dir, paste0(base, "_shape.txt"))
    gp_bridge_write_matrix_bin(kernels[[i]], bin_path, meta_path)
    specs[[i]] <- list(
      bin = normalizePath(bin_path, winslash = "/", mustWork = TRUE),
      meta = normalizePath(meta_path, winslash = "/", mustWork = TRUE)
    )
  }
  list(
    primary_name = kernel_names[[1L]],
    primary = specs[[1L]],
    kernels = specs
  )
}

gp_bridge_prepare_bundle <- function(pheno_df, gmatrix, train_idx0, test_idx0, dir_path) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  bundle <- list(
    pheno_csv = file.path(dir_path, "pheno.csv"),
    kernel_bin = file.path(dir_path, "kernel.bin"),
    kernel_meta = file.path(dir_path, "kernel_shape.txt"),
    geno_ids = file.path(dir_path, "geno_ids.txt"),
    train_idx = file.path(dir_path, "train_idx.txt"),
    test_idx = file.path(dir_path, "test_idx.txt")
  )
  utils::write.csv(pheno_df, bundle$pheno_csv, row.names = FALSE, quote = TRUE)
  gp_bridge_write_matrix_bin(gmatrix, bundle$kernel_bin, bundle$kernel_meta)
  gp_bridge_write_lines(rownames(gmatrix), bundle$geno_ids)
  gp_bridge_write_lines(train_idx0, bundle$train_idx)
  gp_bridge_write_lines(test_idx0, bundle$test_idx)
  bundle
}

gp_bridge_write_json <- function(x, path) {
  writeLines(
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = TRUE, na = "null"),
    con = path,
    useBytes = TRUE
  )
  invisible(path)
}

gp_bridge_read_json <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) list()
  )
}

gp_bridge_read_json_simplified <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = TRUE),
    error = function(e) list()
  )
}

gp_bridge_jsonable_value <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (inherits(x, c("Date", "POSIXt"))) {
    return(as.character(x))
  }
  if (is.data.frame(x)) {
    return(lapply(x, gp_bridge_jsonable_value))
  }
  if (is.matrix(x) || is.array(x)) {
    return(gp_bridge_jsonable_value(as.vector(x)))
  }
  if (is.list(x)) {
    return(lapply(x, gp_bridge_jsonable_value))
  }
  if (is.numeric(x)) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    return(x)
  }
  if (is.integer(x)) {
    return(as.integer(x))
  }
  if (is.logical(x)) {
    return(as.logical(x))
  }
  if (is.character(x)) {
    return(as.character(x))
  }
  x
}

gp_bridge_jsonable_factor_cache <- function(gp_factor_cache) {
  if (is.null(gp_factor_cache)) {
    return(NULL)
  }
  cache <- as.list(gp_factor_cache)
  type <- tolower(as.character(cache$type %||% "zarr")[1L])
  out <- list(type = type)
  if (identical(type, "zarr")) {
    root_dir <- cache$root_dir %||% cache$root
    if (!is.null(root_dir) && nzchar(as.character(root_dir)[1L])) {
      out$root_dir <- normalizePath(as.character(root_dir)[1L], winslash = "/", mustWork = FALSE)
      out$root <- out$root_dir
    }
    out$num_grms <- as.integer(cache$num_grms %||% cache$n_grms %||% 1L)
    out$paths <- gp_bridge_jsonable_value(cache$paths %||% list())
    out$shapes <- gp_bridge_jsonable_value(cache$shapes %||% list())
    out$dtypes <- gp_bridge_jsonable_value(cache$dtypes %||% list())
    return(out)
  }
  if (identical(type, "memmap")) {
    out$paths <- gp_bridge_jsonable_value(cache$paths %||% list())
    out$shapes <- gp_bridge_jsonable_value(cache$shapes %||% list())
    out$dtypes <- gp_bridge_jsonable_value(cache$dtypes %||% list())
    return(out)
  }
  gp_bridge_jsonable_value(cache)
}

gp_bridge_write_factor_cache_memmap <- function(phi, dir_path) {
  cache_dir <- file.path(dir_path, "factor_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  phi <- unname(as.matrix(phi))
  storage.mode(phi) <- "double"
  phi_path <- file.path(cache_dir, "phi0.bin")
  con <- file(phi_path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.double(as.vector(t(phi))), con, size = 8L, endian = "little")
  list(
    type = "memmap",
    paths = list(normalizePath(phi_path, winslash = "/", mustWork = TRUE)),
    shapes = list(as.integer(c(nrow(phi), ncol(phi)))),
    dtypes = list("float64")
  )
}

gp_bridge_write_factor_cache_memmap_bank <- function(kernel, dir_path) {
  kernels <- kernel$Ks
  if (is.null(kernels) || !length(kernels)) {
    kernels <- list(G = kernel$K)
  }
  cache_dir <- file.path(dir_path, "factor_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- vector("list", length(kernels))
  shapes <- vector("list", length(kernels))
  dtypes <- vector("list", length(kernels))
  for (i in seq_along(kernels)) {
    phi <- unname(as.matrix(gp_public_factor_from_kernel(kernels[[i]])))
    storage.mode(phi) <- "double"
    phi_path <- file.path(cache_dir, paste0("phi", i - 1L, ".bin"))
    con <- file(phi_path, open = "wb")
    writeBin(as.double(as.vector(t(phi))), con, size = 8L, endian = "little")
    close(con)
    paths[[i]] <- normalizePath(phi_path, winslash = "/", mustWork = TRUE)
    shapes[[i]] <- as.integer(c(nrow(phi), ncol(phi)))
    dtypes[[i]] <- "float64"
  }
  list(
    type = "memmap",
    paths = paths,
    shapes = shapes,
    dtypes = dtypes
  )
}

gp_bridge_direct_factor_cache <- function(gp_factor_cache,
                                          kernel,
                                          dir_path,
                                          required = FALSE) {
  if (!is.null(gp_factor_cache)) {
    if (!is.list(gp_factor_cache)) {
      stop(
        "Direct GP subprocess execution requires `gp_factor_cache` to be a ",
        "serializable list specification, not an in-memory Python object.",
        call. = FALSE
      )
    }
    return(gp_factor_cache)
  }
  if (!isTRUE(required)) {
    return(NULL)
  }
  gp_bridge_write_factor_cache_memmap_bank(kernel, dir_path = dir_path)
}

gp_bridge_direct_execution_enabled <- function() {
  mode <- tolower(Sys.getenv("PREDICTPRO_GP_EXECUTION", unset = "direct"))
  !mode %in% c("reticulate", "embedded", "r")
}

gp_bridge_reticulate_fallback_enabled <- function() {
  value <- tolower(Sys.getenv("PREDICTPRO_GP_RETICULATE_FALLBACK", unset = "false"))
  value %in% c("1", "true", "yes", "y")
}

gp_bridge_run_cli <- function(command, args, python_bin = NULL) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop(
      "No configured Python runtime was found for Gaussian-process models.\n",
      "Run setup_predictgp_env() or set PREDICTPRO_GP_PYTHON to the dedicated GP environment.",
      call. = FALSE
    )
  }
  cmd_args <- c(gp_bridge_script_path(), command, args)
  out <- system2(
    python_bin,
    args = gp_quote_system_args(cmd_args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) {
    stop(
      paste(
        "GP Python bridge failed.",
        paste(out, collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }
  invisible(out)
}

gp_bridge_prepare_mixed_model_spec <- function(args,
                                               kernel,
                                               gp_factor_cache = NULL,
                                               project_root = NULL,
                                               gp_theta0_warm = NULL,
                                               gp_engine = NULL,
                                               dir_path = tempfile("predictpror_gp_mixed_")) {
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pheno_csv <- file.path(input_dir, "pheno.csv")
  geno_ids_path <- file.path(input_dir, "geno_ids.txt")
  train_idx_path <- file.path(input_dir, "train_idx.txt")
  test_idx_path <- file.path(input_dir, "test_idx.txt")
  utils::write.csv(args$pheno_df, pheno_csv, row.names = FALSE, quote = TRUE)
  kernel_bank <- gp_bridge_write_kernel_bank(kernel, input_dir)
  gp_bridge_write_lines(kernel$geno_ids, geno_ids_path)
  gp_bridge_write_lines(args$train_idx %||% integer(), train_idx_path)
  gp_bridge_write_lines(args$test_idx %||% integer(), test_idx_path)

  excluded <- c(
    "pheno_df", "geno_kernels", "geno_ids", "gid_col", "env_col", "y_col",
    "train_idx", "test_idx", "env_similarity", "env_covariates",
    "grm_factor_cache"
  )
  fit_args <- args[setdiff(names(args), excluded)]
  fit_args <- lapply(fit_args, gp_bridge_jsonable_value)
  # Phase 2.1 Stage B: thread the optional warm-start theta into the Python
  # fit kwargs. gp_framework.fit_mixed_model accepts theta0_warm as a named
  # parameter (added in 0.20.19); it forwards to gp_icm_fa_with_X /
  # gp_exact_with_X / build_*_ai_varcomp where theta0_warm replaces the
  # Adam-init theta0 and tightens max_iter to min(5, max_iter).
  if (!is.null(gp_theta0_warm)) {
    fit_args[["theta0_warm"]] <- as.numeric(gp_theta0_warm)
  }
  if (!is.null(gp_engine)) {
    fit_args[["gp_engine"]] <- as.character(gp_engine)[1L]
  }

  spec <- list(
    project_root = normalizePath(project_root %||% gp_project_root(), winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    pheno_csv = normalizePath(pheno_csv, winslash = "/", mustWork = TRUE),
    kernel_bin = kernel_bank$primary$bin,
    kernel_meta = kernel_bank$primary$meta,
    geno_kernels = kernel_bank$kernels,
    geno_ids = normalizePath(geno_ids_path, winslash = "/", mustWork = TRUE),
    train_idx = normalizePath(train_idx_path, winslash = "/", mustWork = TRUE),
    test_idx = normalizePath(test_idx_path, winslash = "/", mustWork = TRUE),
    gid_col = as.character(args$gid_col %||% "GID")[1L],
    env_col = as.character(args$env_col %||% "Env")[1L],
    y_col = as.character(args$y_col %||% "y")[1L],
    kernel_name = kernel_bank$primary_name,
    fit = fit_args,
    fit_args = fit_args,
    grm_factor_cache = gp_bridge_jsonable_factor_cache(gp_factor_cache %||% args$grm_factor_cache %||% NULL)
  )

  if (!is.null(args$env_similarity)) {
    env_bin <- file.path(input_dir, "env_similarity.bin")
    env_meta <- file.path(input_dir, "env_similarity_shape.txt")
    gp_bridge_write_matrix_bin(args$env_similarity, env_bin, env_meta)
    spec$env_similarity <- list(
      bin = normalizePath(env_bin, winslash = "/", mustWork = TRUE),
      meta = normalizePath(env_meta, winslash = "/", mustWork = TRUE)
    )
  }
  if (!is.null(args$env_ids)) {
    env_ids_path <- file.path(input_dir, "env_ids.txt")
    gp_bridge_write_lines(args$env_ids, env_ids_path)
    spec$env_ids <- normalizePath(env_ids_path, winslash = "/", mustWork = TRUE)
  }
  if (!is.null(args$env_covariates)) {
    env_covariates_csv <- file.path(input_dir, "env_covariates.csv")
    utils::write.csv(args$env_covariates, env_covariates_csv, row.names = FALSE, quote = TRUE)
    spec$env_covariates_csv <- normalizePath(env_covariates_csv, winslash = "/", mustWork = TRUE)
  }

  spec_json <- file.path(input_dir, "fit_spec.json")
  gp_bridge_write_json(spec, spec_json)
  list(
    spec_json = normalizePath(spec_json, winslash = "/", mustWork = TRUE),
    out_dir = spec$out_dir,
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  )
}

gp_bridge_fit_mixed_model_direct <- function(args,
                                             kernel,
                                             gp_factor_cache = NULL,
                                             python_bin = NULL,
                                             project_root = NULL,
                                             gp_theta0_warm = NULL,
                                             gp_engine = NULL) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  spec <- gp_bridge_prepare_mixed_model_spec(
    args = args,
    kernel = kernel,
    gp_factor_cache = gp_factor_cache,
    project_root = project_root,
    gp_theta0_warm = gp_theta0_warm,
    gp_engine = gp_engine
  )
  worker_ok <- FALSE
  if (gp_bridge_worker_use_enabled()) {
    worker_ok <- tryCatch({
      gp_bridge_worker_request(
        "fit-predict-spec",
        list(
          spec = spec$spec_json,
          out_dir = NULL,
          project_root = project_root
        ),
        python_bin = python_bin,
        project_root = project_root
      )
      TRUE
    }, error = function(e) {
      warning("Warm GP worker failed; falling back to subprocess: ", conditionMessage(e), call. = FALSE)
      FALSE
    })
  }
  if (!isTRUE(worker_ok)) {
    gp_bridge_run_cli("fit-predict-spec", c("--spec", spec$spec_json), python_bin = python_bin)
  }

  pred_path <- file.path(spec$out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("GP Python bridge did not write predictions.csv.", call. = FALSE)
  }
  predictions <- utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE)
  result_extra <- gp_bridge_read_json_simplified(file.path(spec$out_dir, "result.json"))
  if (!is.list(result_extra)) {
    result_extra <- list()
  }
  meta <- gp_bridge_read_json(file.path(spec$out_dir, "meta.json"))
  info <- gp_bridge_read_json(file.path(spec$out_dir, "info.json"))
  fit <- gp_bridge_read_json(file.path(spec$out_dir, "fit.json"))
  diagnostics <- gp_bridge_read_json(file.path(spec$out_dir, "diagnostics.json"))
  out_info <- c(
    list(
      execution = "python_subprocess",
      gp_execution = "python_subprocess",
      bridge_command = "fit-predict-spec",
      bridge_dir = spec$bridge_dir,
      python = normalizePath(python_bin, winslash = "/", mustWork = FALSE)
    ),
    meta,
    info
  )
  if (length(diagnostics)) {
    out_info$diagnostics <- diagnostics
  }
  out <- list(
    result = c(list(predictions = gp_plain_data_frame(predictions)), result_extra),
    predictions = gp_plain_data_frame(predictions),
    fit = if (length(fit)) fit else list(),
    info = out_info
  )
  class(out) <- c("predictpror_gp_result", class(out))
  out
}

gp_backend_cv_batch_enabled <- function() {
  mode <- tolower(Sys.getenv("PREDICTPRO_GP_BATCH_CV", unset = "true"))
  mode %in% c("1", "true", "yes", "y", "auto")
}

gp_normalize_exact_fast_cv_mode <- function(value = NULL) {
  if (is.null(value) || length(value) < 1L) {
    return(NULL)
  }
  if (is.logical(value)) {
    return(if (isTRUE(value[[1L]])) "true" else "false")
  }
  mode <- tolower(trimws(as.character(value[[1L]])))
  if (!nzchar(mode) || mode %in% c("default", "env", "environment")) {
    return(NULL)
  }
  if (mode %in% c("0", "false", "f", "no", "n", "off", "disabled", "disable")) {
    return("false")
  }
  if (mode %in% c("1", "true", "t", "yes", "y", "on", "enabled", "enable", "auto")) {
    return("true")
  }
  if (mode %in% c("force", "forced", "always")) {
    return("force")
  }
  stop(
    "`gp_exact_fast_cv` must be NULL, TRUE/FALSE, 'auto', 'off', or 'force'.",
    call. = FALSE
  )
}

gp_bridge_prepare_mixed_model_cv_spec <- function(args,
                                                  kernel,
                                                  folds,
                                                  gp_factor_cache = NULL,
                                                  project_root = NULL,
                                                  dir_path = tempfile("predictpror_gp_cv_mixed_")) {
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  fold_dir <- file.path(input_dir, "folds")
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pheno_csv <- file.path(input_dir, "pheno.csv")
  geno_ids_path <- file.path(input_dir, "geno_ids.txt")
  utils::write.csv(args$pheno_df, pheno_csv, row.names = FALSE, quote = TRUE)
  kernel_bank <- gp_bridge_write_kernel_bank(kernel, input_dir)
  gp_bridge_write_lines(kernel$geno_ids, geno_ids_path)

  fold_specs <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    fold <- folds[[i]]
    fold_id <- as.character(fold$fold_id %||% paste0("fold", i))[1L]
    train_idx <- as.integer(fold$train_idx %||% integer())
    test_idx <- as.integer(fold$test_idx %||% integer())
    train_idx_path <- file.path(fold_dir, paste0("fold_", i, "_train_idx.txt"))
    test_idx_path <- file.path(fold_dir, paste0("fold_", i, "_test_idx.txt"))
    gp_bridge_write_lines(train_idx, train_idx_path)
    gp_bridge_write_lines(test_idx, test_idx_path)
    fold_specs[[i]] <- list(
      fold_id = fold_id,
      train_idx = normalizePath(train_idx_path, winslash = "/", mustWork = TRUE),
      test_idx = normalizePath(test_idx_path, winslash = "/", mustWork = TRUE)
    )
  }

  excluded <- c(
    "pheno_df", "geno_kernels", "geno_ids", "gid_col", "env_col", "y_col",
    "train_idx", "test_idx", "folds", "env_similarity", "env_covariates",
    "env_ids", "grm_factor_cache", "gp_exact_fast_cv"
  )
  fit_args <- args[setdiff(names(args), excluded)]
  fit_args <- lapply(fit_args, gp_bridge_jsonable_value)
  gp_exact_fast_cv <- gp_normalize_exact_fast_cv_mode(args$gp_exact_fast_cv %||% NULL)

  spec <- list(
    project_root = normalizePath(project_root %||% gp_project_root(), winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    pheno_csv = normalizePath(pheno_csv, winslash = "/", mustWork = TRUE),
    kernel_bin = kernel_bank$primary$bin,
    kernel_meta = kernel_bank$primary$meta,
    geno_kernels = kernel_bank$kernels,
    geno_ids = normalizePath(geno_ids_path, winslash = "/", mustWork = TRUE),
    folds = fold_specs,
    gid_col = as.character(args$gid_col %||% "GID")[1L],
    env_col = as.character(args$env_col %||% "Env")[1L],
    y_col = as.character(args$y_col %||% "y")[1L],
    kernel_name = kernel_bank$primary_name,
    fit = fit_args,
    fit_args = fit_args,
    grm_factor_cache = gp_bridge_jsonable_factor_cache(gp_factor_cache %||% args$grm_factor_cache %||% NULL)
  )
  if (!is.null(gp_exact_fast_cv)) {
    spec$gp_exact_fast_cv <- gp_exact_fast_cv
  }

  if (!is.null(args$env_similarity)) {
    env_bin <- file.path(input_dir, "env_similarity.bin")
    env_meta <- file.path(input_dir, "env_similarity_shape.txt")
    gp_bridge_write_matrix_bin(args$env_similarity, env_bin, env_meta)
    spec$env_similarity <- list(
      bin = normalizePath(env_bin, winslash = "/", mustWork = TRUE),
      meta = normalizePath(env_meta, winslash = "/", mustWork = TRUE)
    )
  }
  if (!is.null(args$env_ids)) {
    env_ids_path <- file.path(input_dir, "env_ids.txt")
    gp_bridge_write_lines(args$env_ids, env_ids_path)
    spec$env_ids <- normalizePath(env_ids_path, winslash = "/", mustWork = TRUE)
  }
  if (!is.null(args$env_covariates)) {
    env_covariates_csv <- file.path(input_dir, "env_covariates.csv")
    utils::write.csv(args$env_covariates, env_covariates_csv, row.names = FALSE, quote = TRUE)
    spec$env_covariates_csv <- normalizePath(env_covariates_csv, winslash = "/", mustWork = TRUE)
  }

  spec_json <- file.path(input_dir, "fit_cv_spec.json")
  gp_bridge_write_json(spec, spec_json)
  list(
    spec_json = normalizePath(spec_json, winslash = "/", mustWork = TRUE),
    out_dir = spec$out_dir,
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  )
}

gp_bridge_fit_mixed_model_cv_direct <- function(args,
                                                kernel,
                                                folds,
                                                gp_factor_cache = NULL,
                                                python_bin = NULL,
                                                project_root = NULL,
                                                read_diagnostics = TRUE) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  spec <- gp_bridge_prepare_mixed_model_cv_spec(
    args = args,
    kernel = kernel,
    folds = folds,
    gp_factor_cache = gp_factor_cache,
    project_root = project_root
  )
  worker_ok <- FALSE
  if (gp_bridge_worker_use_enabled()) {
    worker_ok <- tryCatch({
      gp_bridge_worker_request(
        "fit-predict-cv-spec",
        list(
          spec = spec$spec_json,
          out_dir = NULL,
          project_root = project_root
        ),
        python_bin = python_bin,
        project_root = project_root
      )
      TRUE
    }, error = function(e) {
      warning("Warm GP worker failed for batched CV; falling back to subprocess: ", conditionMessage(e), call. = FALSE)
      FALSE
    })
  }
  if (!isTRUE(worker_ok)) {
    gp_bridge_run_cli("fit-predict-cv-spec", c("--spec", spec$spec_json), python_bin = python_bin)
  }

  pred_path <- file.path(spec$out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("GP Python bridge did not write batched CV predictions.csv.", call. = FALSE)
  }
  predictions <- utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE)
  meta <- gp_bridge_read_json(file.path(spec$out_dir, "meta.json"))
  fold_summary <- NULL
  cv_metrics <- NULL
  if (isTRUE(read_diagnostics)) {
    folds_csv <- file.path(spec$out_dir, "folds.csv")
    metrics_csv <- file.path(spec$out_dir, "metrics.csv")
    fold_summary <- if (file.exists(folds_csv)) {
      gp_plain_data_frame(utils::read.csv(folds_csv, stringsAsFactors = FALSE, check.names = FALSE))
    } else {
      NULL
    }
    cv_metrics <- if (file.exists(metrics_csv)) {
      gp_plain_data_frame(utils::read.csv(metrics_csv, stringsAsFactors = FALSE, check.names = FALSE))
    } else {
      NULL
    }
  }
  out_info <- c(
    list(
      execution = "python_subprocess",
      gp_execution = "python_subprocess",
      package_cv_batch = TRUE,
      bridge_command = "fit-predict-cv-spec",
      bridge_dir = spec$bridge_dir,
      python = normalizePath(python_bin, winslash = "/", mustWork = FALSE)
    ),
    meta
  )
  out <- list(
    result = list(predictions = gp_plain_data_frame(predictions)),
    predictions = gp_plain_data_frame(predictions),
    fit = list(folds = fold_summary, metrics = cv_metrics),
    info = out_info
  )
  class(out) <- c("predictpror_gp_result", class(out))
  out
}

gp_bridge_prepare_public_model_spec <- function(args,
                                                kernel,
                                                gp_factor_cache = NULL,
                                                project_root = NULL,
                                                command = "fit-multi-trait-spec",
                                                dir_path = tempfile("predictpror_gp_public_")) {
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pheno_csv <- file.path(input_dir, "pheno.csv")
  geno_ids_path <- file.path(input_dir, "geno_ids.txt")
  train_idx_path <- file.path(input_dir, "train_idx.txt")
  test_idx_path <- file.path(input_dir, "test_idx.txt")
  utils::write.csv(args$pheno_df, pheno_csv, row.names = FALSE, quote = TRUE)
  kernel_bank <- gp_bridge_write_kernel_bank(kernel, input_dir)
  gp_bridge_write_lines(kernel$geno_ids, geno_ids_path)
  gp_bridge_write_lines(args$train_idx %||% integer(), train_idx_path)
  gp_bridge_write_lines(args$test_idx %||% integer(), test_idx_path)

  excluded <- c(
    "pheno_df", "geno_kernels", "geno_ids", "gid_col", "env_col",
    "trait_col", "y_col", "train_idx", "test_idx", "env_similarity",
    "env_covariates", "env_ids", "grm_factor_cache"
  )
  fit_args <- args[setdiff(names(args), excluded)]
  fit_args <- lapply(fit_args, gp_bridge_jsonable_value)
  columns <- list(
    gid_col = as.character(args$gid_col %||% "gid")[1L],
    trait_col = as.character(args$trait_col %||% "trait")[1L],
    y_col = as.character(args$y_col %||% "y")[1L]
  )
  if (!is.null(args$env_col)) {
    columns$env_col <- as.character(args$env_col)[1L]
  }

  spec <- list(
    project_root = normalizePath(project_root %||% gp_project_root(), winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    pheno_csv = normalizePath(pheno_csv, winslash = "/", mustWork = TRUE),
    kernel_bin = kernel_bank$primary$bin,
    kernel_meta = kernel_bank$primary$meta,
    geno_kernels = kernel_bank$kernels,
    kernel_name = kernel_bank$primary_name,
    geno_ids = normalizePath(geno_ids_path, winslash = "/", mustWork = TRUE),
    train_idx = normalizePath(train_idx_path, winslash = "/", mustWork = TRUE),
    test_idx = normalizePath(test_idx_path, winslash = "/", mustWork = TRUE),
    columns = columns,
    fit = fit_args,
    fit_args = fit_args,
    grm_factor_cache = gp_bridge_jsonable_factor_cache(gp_factor_cache %||% args$grm_factor_cache %||% NULL)
  )

  if (!is.null(args$env_similarity)) {
    env_bin <- file.path(input_dir, "env_similarity.bin")
    env_meta <- file.path(input_dir, "env_similarity_shape.txt")
    gp_bridge_write_matrix_bin(args$env_similarity, env_bin, env_meta)
    spec$env_similarity <- list(
      bin = normalizePath(env_bin, winslash = "/", mustWork = TRUE),
      meta = normalizePath(env_meta, winslash = "/", mustWork = TRUE)
    )
  }
  if (!is.null(args$env_ids)) {
    env_ids_path <- file.path(input_dir, "env_ids.txt")
    gp_bridge_write_lines(args$env_ids, env_ids_path)
    spec$env_ids <- normalizePath(env_ids_path, winslash = "/", mustWork = TRUE)
  }
  if (!is.null(args$env_covariates)) {
    env_covariates_csv <- file.path(input_dir, "env_covariates.csv")
    utils::write.csv(args$env_covariates, env_covariates_csv, row.names = FALSE, quote = TRUE)
    spec$env_covariates_csv <- normalizePath(env_covariates_csv, winslash = "/", mustWork = TRUE)
  }

  spec_json <- file.path(input_dir, "fit_spec.json")
  gp_bridge_write_json(spec, spec_json)
  list(
    spec_json = normalizePath(spec_json, winslash = "/", mustWork = TRUE),
    out_dir = spec$out_dir,
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE),
    command = command
  )
}

gp_bridge_public_result_from_dir <- function(out_dir,
                                             bridge_dir,
                                             command,
                                             python_bin = NULL) {
  pred_path <- file.path(out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("GP Python bridge did not write predictions.csv.", call. = FALSE)
  }
  predictions <- gp_plain_data_frame(utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE))
  result_extra <- gp_bridge_read_json_simplified(file.path(out_dir, "result.json"))
  if (!is.list(result_extra)) {
    result_extra <- list()
  }
  fit <- gp_bridge_read_json_simplified(file.path(out_dir, "fit.json"))
  info <- gp_bridge_read_json_simplified(file.path(out_dir, "info.json"))
  meta <- gp_bridge_read_json_simplified(file.path(out_dir, "meta.json"))
  out_info <- c(
    list(
      execution = "python_subprocess",
      gp_execution = "python_subprocess",
      bridge_command = command,
      bridge_dir = bridge_dir,
      python = normalizePath(python_bin %||% gp_detect_gp_python(), winslash = "/", mustWork = FALSE)
    ),
    if (is.list(meta)) meta else list(),
    if (is.list(info)) info else list()
  )
  out <- list(
    result = c(list(predictions = predictions), result_extra),
    predictions = predictions,
    fit = if (is.list(fit)) fit else list(),
    info = out_info
  )
  class(out) <- c("predictpror_gp_result", class(out))
  out
}

gp_bridge_fit_public_model_direct <- function(args,
                                              kernel,
                                              gp_factor_cache = NULL,
                                              command,
                                              python_bin = NULL,
                                              project_root = NULL,
                                              require_factor_cache = FALSE) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  bridge_dir <- tempfile("predictpror_gp_public_")
  direct_cache <- gp_bridge_direct_factor_cache(
    gp_factor_cache = gp_factor_cache,
    kernel = kernel,
    dir_path = bridge_dir,
    required = require_factor_cache
  )
  spec <- gp_bridge_prepare_public_model_spec(
    args = args,
    kernel = kernel,
    gp_factor_cache = direct_cache,
    project_root = project_root,
    command = command,
    dir_path = bridge_dir
  )
  worker_ok <- FALSE
  if (gp_bridge_worker_use_enabled()) {
    worker_ok <- tryCatch({
      gp_bridge_worker_request(
        command,
        list(
          spec = spec$spec_json,
          out_dir = NULL,
          project_root = project_root
        ),
        python_bin = python_bin,
        project_root = project_root
      )
      TRUE
    }, error = function(e) {
      warning("Warm GP worker failed; falling back to subprocess: ", conditionMessage(e), call. = FALSE)
      FALSE
    })
  }
  if (!isTRUE(worker_ok)) {
    gp_bridge_run_cli(command, c("--spec", spec$spec_json), python_bin = python_bin)
  }
  gp_bridge_public_result_from_dir(
    out_dir = spec$out_dir,
    bridge_dir = spec$bridge_dir,
    command = command,
    python_bin = python_bin
  )
}

PredictProR_gp_framework_env <- new.env(parent = emptyenv())

gp_framework_module_name <- function() {
  "gp_framework"
}

gp_framework_module_key <- function(python_bin, project_root) {
  paste(
    normalizePath(python_bin, winslash = "/", mustWork = FALSE),
    normalizePath(project_root, winslash = "/", mustWork = FALSE),
    sep = "::"
  )
}

gp_get_framework_module <- function(python_bin = NULL, project_root = NULL) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop("No configured GP Python runtime found.", call. = FALSE)
  }
  gp_require_reticulate("GP framework module import")
  gp_init_python_once(python_bin)
  key <- gp_framework_module_key(python_bin, project_root)
  if (exists(key, envir = PredictProR_gp_framework_env, inherits = FALSE)) {
    return(get(key, envir = PredictProR_gp_framework_env, inherits = FALSE))
  }
  reticulate::py_run_string(sprintf(
    "import sys\n_gp_dir = r'%s'\nif _gp_dir not in sys.path:\n    sys.path.insert(0, _gp_dir)\n",
    project_root
  ))
  mod <- reticulate::import(gp_framework_module_name(), delay_load = FALSE, convert = FALSE)
  assign(key, mod, envir = PredictProR_gp_framework_env)
  mod
}

#' Prewarm the PredictGP runtime
#'
#' @param python_bin Optional Python executable to use for the GP runtime.
#' @param project_root Optional project root containing the GP Python backend.
#'
#' @return Invisibly returns \code{TRUE}.
#' @export
gp_prewarm_runtime <- function(python_bin = NULL, project_root = NULL) {
  gp_get_framework_module(
    python_bin = python_bin %||% gp_detect_gp_python(),
    project_root = project_root %||% gp_project_root()
  )
  invisible(TRUE)
}

PredictProR_gp_worker_env <- new.env(parent = emptyenv())
PredictProR_gp_worker_cleanup_registered <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# Phase 2.1 Stage C: per-session theta warm-start cache.
#
# Stores the converged variance-component theta returned by the Python
# AI-REML driver, keyed by (model_name, response, gen_name, heter_groups,
# gp_varcomp_mode, gmatrix-digest). The next compatible fit (same key)
# reuses this theta as theta0_warm, dropping AI-Newton from 8-12 cold
# iterations to 1-2 warm iterations.
#
# The cache is per-R-session (process-local). Under mirai parallel CV
# dispatch each worker gets its own copy, which is the correct isolation:
# warm-start state is local to the worker that produced it.
# ---------------------------------------------------------------------------
PredictProR_gp_warm_start_env <- new.env(parent = emptyenv())

gp_warm_start_cache_key <- function(model_name,
                                    response,
                                    gen_name,
                                    heter_groups,
                                    gp_varcomp_mode,
                                    gmatrix) {
  # Use tools::md5sum on a serialized blob of the matrix. md5sum requires a
  # file path so we write the serialized bytes to a tempfile and md5 it.
  # The tempfile is removed immediately. Same input -> same digest within
  # one session. Negligible cost vs. a full GP fit.
  if (is.null(gmatrix)) {
    matrix_hash <- "no_gmatrix"
  } else {
    tf <- tempfile("predictpror_gp_warm_hash_", fileext = ".bin")
    on.exit(if (file.exists(tf)) unlink(tf), add = TRUE)
    suppressWarnings(saveRDS(unname(as.matrix(gmatrix)), tf, ascii = FALSE,
                             compress = FALSE))
    matrix_hash <- tryCatch(
      unname(as.character(tools::md5sum(tf))),
      error = function(e) sprintf("%s-%s", nrow(gmatrix), ncol(gmatrix))
    )
    if (is.na(matrix_hash) || !nzchar(matrix_hash)) {
      matrix_hash <- sprintf("%s-%s", nrow(gmatrix), ncol(gmatrix))
    }
  }
  paste(
    as.character(model_name %||% "NA"),
    as.character(response %||% "NA"),
    as.character(gen_name %||% "NA"),
    as.character(heter_groups %||% "NA_heter"),
    as.character(gp_varcomp_mode %||% "reml"),
    matrix_hash,
    sep = "|"
  )
}

gp_warm_start_lookup <- function(key) {
  if (!is.character(key) || !nzchar(key)) return(NULL)
  if (!exists(key, envir = PredictProR_gp_warm_start_env, inherits = FALSE)) {
    return(NULL)
  }
  obj <- get(key, envir = PredictProR_gp_warm_start_env, inherits = FALSE)
  if (!is.list(obj) || is.null(obj[["theta"]])) return(NULL)
  if (!is.numeric(obj[["theta"]]) || !length(obj[["theta"]])) return(NULL)
  obj[["theta"]]
}

gp_warm_start_save <- function(key, theta, layout = NA_character_) {
  if (!is.character(key) || !nzchar(key)) return(invisible(NULL))
  if (is.null(theta)) return(invisible(NULL))
  theta_num <- suppressWarnings(as.numeric(theta))
  if (!length(theta_num) || any(!is.finite(theta_num))) {
    return(invisible(NULL))
  }
  assign(
    key,
    list(theta = theta_num, layout = as.character(layout)[1L], saved_at = Sys.time()),
    envir = PredictProR_gp_warm_start_env
  )
  invisible(NULL)
}

gp_warm_start_clear <- function() {
  rm(
    list = ls(envir = PredictProR_gp_warm_start_env, all.names = TRUE),
    envir = PredictProR_gp_warm_start_env
  )
  invisible(NULL)
}

gp_bridge_worker_enabled <- function() {
  TRUE
}

gp_bridge_worker_key <- function(python_bin, project_root) {
  paste(normalizePath(python_bin, winslash = "/", mustWork = FALSE),
        normalizePath(project_root, winslash = "/", mustWork = FALSE),
        sep = "::")
}

gp_bridge_worker_use_enabled <- function() {
  value <- tolower(Sys.getenv("PREDICTPRO_GP_WORKER", unset = "true"))
  value %in% c("1", "true", "yes", "y", "auto")
}

gp_bridge_worker_alive <- function(worker) {
  if (!inherits(worker, "gp_bridge_worker") || is.null(worker$port)) {
    return(FALSE)
  }
  deadline <- Sys.time() + 3
  repeat {
    ok <- tryCatch({
      con <- suppressWarnings(socketConnection(
        host = "127.0.0.1",
        port = as.integer(worker$port),
        blocking = TRUE,
        open = "a+b",
        timeout = 1
      ))
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      payload <- jsonlite::toJSON(list(cmd = "ping", args = list()), auto_unbox = TRUE, null = "null")
      writeBin(charToRaw(paste0(payload, "\n")), con)
      flush(con)
      line <- readLines(con, n = 1L, warn = FALSE)
      msg <- jsonlite::fromJSON(line[[1L]])
      identical(msg$status, "ok")
    }, error = function(e) FALSE)
    if (isTRUE(ok) || Sys.time() >= deadline) {
      return(isTRUE(ok))
    }
    Sys.sleep(0.25)
  }
}

gp_bridge_find_free_port <- function() {
  if (requireNamespace("parallelly", quietly = TRUE)) {
    port <- tryCatch(parallelly::freePort(), error = function(e) NA_integer_)
    if (is.finite(port) && !is.na(port)) {
      return(as.integer(port))
    }
  }
  as.integer(sample(20000:49151, 1L))
}

gp_bridge_worker_start <- function(python_bin = NULL, project_root = NULL) {
  if (!gp_bridge_worker_use_enabled() || !gp_bridge_worker_enabled()) {
    return(NULL)
  }
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  key <- gp_bridge_worker_key(python_bin, project_root)
  if (exists(key, envir = PredictProR_gp_worker_env, inherits = FALSE)) {
    proc <- get(key, envir = PredictProR_gp_worker_env, inherits = FALSE)
    if ((inherits(proc, "process") && proc$is_alive()) ||
        gp_bridge_worker_alive(proc)) {
      return(proc)
    }
  }
  port <- gp_bridge_find_free_port()
  log_file <- tempfile("gp_worker_", fileext = ".log")
  status <- suppressWarnings(system2(
    python_bin,
    args = gp_quote_system_args(c(
      gp_bridge_script_path(), "serve-socket", "--project-root", project_root,
      "--host", "127.0.0.1", "--port", as.character(port)
    )),
    wait = FALSE,
    stdout = log_file,
    stderr = log_file
  ))
  Sys.sleep(1)
  worker <- list(port = port, log_file = log_file, status = status, python = python_bin, project_root = project_root)
  class(worker) <- "gp_bridge_worker"
  assign(key, worker, envir = PredictProR_gp_worker_env)
  if (!exists("registered", envir = PredictProR_gp_worker_cleanup_registered, inherits = FALSE)) {
    reg.finalizer(
      PredictProR_gp_worker_env,
      function(e) {
        try(gp_bridge_worker_shutdown_all(), silent = TRUE)
      },
      onexit = TRUE
    )
    assign("registered", TRUE, envir = PredictProR_gp_worker_cleanup_registered)
  }
  worker
}

gp_bridge_worker_request <- function(command, args, python_bin = NULL, project_root = NULL) {
  proc <- gp_bridge_worker_start(python_bin = python_bin, project_root = project_root)
  if (is.null(proc)) {
    return(NULL)
  }
  timeout <- as.integer(Sys.getenv("PREDICTPRO_GP_WORKER_TIMEOUT", unset = "3600"))[1L]
  if (!is.finite(timeout) || is.na(timeout) || timeout < 1L) {
    timeout <- 3600L
  }
  con <- NULL
  deadline <- Sys.time() + min(30, timeout)
  repeat {
    con <- tryCatch(
      socketConnection(
        host = "127.0.0.1",
        port = proc$port,
        blocking = TRUE,
        open = "a+b",
        timeout = timeout
      ),
      error = function(e) NULL
    )
    if (!is.null(con) || Sys.time() >= deadline) {
      break
    }
    Sys.sleep(0.25)
  }
  if (is.null(con)) {
    log_txt <- if (file.exists(proc$log_file)) paste(readLines(proc$log_file, warn = FALSE), collapse = "\n") else ""
    stop("Could not connect to warm GP bridge worker.\n", log_txt, call. = FALSE)
  }
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  payload <- jsonlite::toJSON(list(cmd = command, args = args), auto_unbox = TRUE, null = "null")
  writeBin(charToRaw(paste0(payload, "\n")), con)
  flush(con)
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) {
    log_txt <- if (file.exists(proc$log_file)) paste(readLines(proc$log_file, warn = FALSE), collapse = "\n") else ""
    stop("GP bridge worker returned no response.\n", log_txt, call. = FALSE)
  }
  msg <- tryCatch(jsonlite::fromJSON(line[[1L]]), error = function(e) NULL)
  if (is.null(msg)) {
    stop("GP bridge worker returned invalid JSON.", call. = FALSE)
  }
  if (!identical(msg$status, "ok")) {
    stop("GP bridge worker failed.\n", msg$message %||% "unknown error", call. = FALSE)
  }
  invisible(msg)
}

gp_bridge_worker_shutdown_all <- function() {
  keys <- ls(envir = PredictProR_gp_worker_env, all.names = TRUE)
  for (key in keys) {
    worker <- get(key, envir = PredictProR_gp_worker_env, inherits = FALSE)
    if (!inherits(worker, "gp_bridge_worker")) {
      next
    }
    try({
      con <- suppressWarnings(socketConnection(
        host = "127.0.0.1",
        port = as.integer(worker$port),
        blocking = TRUE,
        open = "a+b",
        timeout = 2
      ))
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      payload <- jsonlite::toJSON(list(cmd = "shutdown", args = list()), auto_unbox = TRUE, null = "null")
      writeBin(charToRaw(paste0(payload, "\n")), con)
      flush(con)
      invisible(readLines(con, n = 1L, warn = FALSE))
    }, silent = TRUE)
    rm(list = key, envir = PredictProR_gp_worker_env)
  }
  invisible(TRUE)
}

gp_bridge_collect_prediction <- function(pred_df,
                                         pheno_df,
                                         response_name,
                                         gen_name,
                                         test_mask) {
  cols <- gp_backend_schema_cols()$single
  pred_df <- as.data.frame(pred_df, stringsAsFactors = FALSE)
  names(pred_df) <- sub("^Name$", cols$gid, names(pred_df))
  names(pred_df) <- sub("^Prediction$", "Predicted_value", names(pred_df))
  merged <- merge(
    pheno_df[, unname(c(cols$gid, cols$env, cols$value)), drop = FALSE],
    pred_df,
    by = intersect(unname(c(cols$gid, cols$env)), names(pred_df)),
    all.x = TRUE,
    sort = FALSE
  )
  merged <- merged[match(
    seq_len(nrow(pheno_df)),
    match(
      paste(pheno_df[[cols$gid]], pheno_df[[cols$env]]),
      paste(merged[[cols$gid]], merged[[cols$env]])
    )
  ), , drop = FALSE]
  train_test_label <- ifelse(test_mask, "Test", "Train")
  pred_var <- if ("Prediction_variance" %in% names(merged)) {
    as.numeric(merged[["Prediction_variance"]])
  } else if ("Prediction_Var_observed" %in% names(merged)) {
    as.numeric(merged[["Prediction_Var_observed"]])
  } else if ("Prediction_Var_latent" %in% names(merged)) {
    as.numeric(merged[["Prediction_Var_latent"]])
  } else {
    rep(NA_real_, nrow(merged))
  }
  pred_se <- if ("Prediction_SE" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE"]])
  } else if ("Prediction_SE_observed" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE_observed"]])
  } else if ("Prediction_SE_latent" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE_latent"]])
  } else if ("SE" %in% names(merged)) {
    as.numeric(merged[["SE"]])
  } else if (all(is.na(pred_var))) {
    rep(NA_real_, nrow(merged))
  } else {
    sqrt(pmax(pred_var, 0))
  }
  if (all(is.na(pred_var)) && !all(is.na(pred_se))) {
    pred_var <- pred_se^2
  }

  out <- gp_ml_gaussian_prediction_output(
    ids = as.character(pheno_df[[cols$gid]]),
    predicted_value = as.numeric(merged[["Predicted_value"]]),
    train_test_label = train_test_label,
    pred_se = pred_se,
    pred_variances = pred_var,
    lower_bound = as.numeric(merged[["Predicted_value"]]) - 1.96 * pred_se,
    upper_bound = as.numeric(merged[["Predicted_value"]]) + 1.96 * pred_se,
    uncertainty = 3.92 * pred_se,
    uncertainty_remarks = ifelse(is.na(pred_se), NA_character_, "model_based"),
    reference_variance = stats::var(pheno_df[[cols$value]][!test_mask], na.rm = TRUE) %||% 1,
    observed_y = pheno_df[[cols$value]][!test_mask],
    predictor_matrix = NULL,
    gen_name = gen_name
  )
  out[["Env"]] <- as.character(pheno_df[[cols$env]])
  out
}

gp_backend_light_prediction_output <- function(pred_df,
                                               pheno_df,
                                               gen_name,
                                               test_mask) {
  cols <- gp_backend_schema_cols()$single
  pred_df <- as.data.frame(pred_df, stringsAsFactors = FALSE, check.names = FALSE)
  names(pred_df) <- sub("^Name$", cols$gid, names(pred_df))
  names(pred_df) <- sub("^Prediction$", "Predicted_value", names(pred_df))

  merged <- merge(
    pheno_df[, unname(c(cols$gid, cols$env, cols$value)), drop = FALSE],
    pred_df,
    by = intersect(unname(c(cols$gid, cols$env)), names(pred_df)),
    all.x = TRUE,
    sort = FALSE
  )
  merged <- merged[match(
    seq_len(nrow(pheno_df)),
    match(
      paste(pheno_df[[cols$gid]], pheno_df[[cols$env]]),
      paste(merged[[cols$gid]], merged[[cols$env]])
    )
  ), , drop = FALSE]

  pred_var <- if ("Prediction_Var_observed" %in% names(merged)) {
    as.numeric(merged[["Prediction_Var_observed"]])
  } else if ("Prediction_Var_latent" %in% names(merged)) {
    as.numeric(merged[["Prediction_Var_latent"]])
  } else if ("Prediction_variance" %in% names(merged)) {
    as.numeric(merged[["Prediction_variance"]])
  } else {
    rep(NA_real_, nrow(merged))
  }
  pred_se <- if ("Prediction_SE_observed" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE_observed"]])
  } else if ("Prediction_SE_latent" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE_latent"]])
  } else if ("Prediction_SE" %in% names(merged)) {
    as.numeric(merged[["Prediction_SE"]])
  } else if (!all(is.na(pred_var))) {
    sqrt(pmax(pred_var, 0))
  } else {
    rep(NA_real_, nrow(merged))
  }
  if (all(is.na(pred_var)) && !all(is.na(pred_se))) {
    pred_var <- pred_se^2
  }

  out <- data.frame(
    GID = as.character(pheno_df[[cols$gid]]),
    Env = as.character(pheno_df[[cols$env]]),
    Predicted_value = as.numeric(merged[["Predicted_value"]]),
    Train_Test_Label = ifelse(test_mask, "Test", "Train"),
    Standard_error = pred_se,
    PEV = pred_var,
    Observed_value = as.numeric(pheno_df[[cols$value]]),
    Train_Test = ifelse(test_mask, "Test", "Train"),
    stringsAsFactors = FALSE
  )
  out
}

gp_bridge_total_prediction_output <- function(pred_df, confidence_level = 0.95,
                                              env_col = NULL) {
  if (is.null(pred_df) || !is.data.frame(pred_df)) {
    return(NULL)
  }
  pred_df <- as.data.frame(pred_df, stringsAsFactors = FALSE, check.names = FALSE)
  # Phase 3.26 fix: respect user-supplied env col (heter_groups) -- the
  # public single-trait formatter renames "Env" to the user's heter_groups
  # name (e.g. "YYY"), so a hardcoded "Env" lookup returns NULL and the
  # across-environment aggregate is silently dropped.
  candidates <- unique(c(
    as.character(env_col),
    "Env", "env", "Environment", "environment", "Loc", "loc", "Location", "location"
  ))
  candidates <- candidates[nzchar(candidates) & candidates %in% names(pred_df)]
  if (!length(candidates)) {
    return(NULL)
  }
  env_resolved <- candidates[[1L]]
  if (length(unique(as.character(pred_df[[env_resolved]]))) <= 1L) {
    return(NULL)
  }
  gid_col <- if ("GID" %in% names(pred_df)) "GID" else names(pred_df)[[1L]]
  split_rows <- split(seq_len(nrow(pred_df)), as.character(pred_df[[gid_col]]), drop = TRUE)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  rows <- lapply(names(split_rows), function(gid) {
    idx <- split_rows[[gid]]
    n_env <- length(idx)
    pred <- mean(as.numeric(pred_df[["Predicted_value"]][idx]), na.rm = TRUE)
    obs <- mean(as.numeric(pred_df[["Observed_value"]][idx]), na.rm = TRUE)
    if (!is.finite(obs)) {
      obs <- NA_real_
    }
    pev_vals <- suppressWarnings(as.numeric(pred_df[["PEV"]][idx]))
    pev <- if (any(is.finite(pev_vals))) {
      sum(pev_vals[is.finite(pev_vals)], na.rm = TRUE) / (n_env^2)
    } else {
      NA_real_
    }
    se <- sqrt(pmax(pev, 0))
    rel_vals <- suppressWarnings(as.numeric(pred_df[["Reliability"]][idx]))
    rel <- if (any(is.finite(rel_vals))) mean(rel_vals, na.rm = TRUE) else NA_real_
    label <- if (any(as.character(pred_df[["Train_Test_Label"]][idx]) == "Test", na.rm = TRUE)) {
      "Test"
    } else {
      "Train"
    }
    data.frame(
      GID = gid,
      Predicted_value = pred,
      Train_Test_Label = label,
      Observed_value = obs,
      Standard_error = se,
      PEV = pev,
      lower_bound = pred - z * se,
      upper_bound = pred + z * se,
      Uncertainty = 2 * z * se,
      Uncertainty_remarks = "across_environment_aggregate",
      Reliability = rel,
      Reliability_remarks = ifelse(
        is.na(rel),
        NA_character_,
        ifelse(rel >= 0.9, "Reliable", ifelse(rel >= 0.5, "Acceptable", "Unreliable"))
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gp_bridge_fit_predict <- function(model_name,
                                  pheno_data,
                                  response,
                                  gen_name,
                                  gmatrix,
                                  heter_groups = NULL,
                                  test_set = NULL,
                                  test_set_source = NULL,
                                  gp_backend = "auto",
                                  gp_output_level = "predict_only",
                                  gp_varcomp_mode = "reml",
                                  gp_fa_rank = 1L,
                                  gp_prediction_output = "all",
                                  gp_factor_cache = NULL,
                                  gp_iters = NULL,
                                  gp_lr = NULL,
                                  seed = 12345L) {
  python_bin <- gp_detect_gp_python()
  pheno_df <- gp_bridge_build_pheno(pheno_data, response = response, gen_name = gen_name, heter_groups = heter_groups)
  split <- gp_bridge_train_test_idx(pheno_df, test_set = test_set, test_set_source = test_set_source)
  cols <- gp_backend_schema_cols()$single
  params <- gp_bridge_model_params(
    model_name = model_name,
    gp_backend = gp_backend,
    gp_output_level = gp_output_level,
    gp_varcomp_mode = gp_varcomp_mode,
    gp_fa_rank = gp_fa_rank,
    gp_prediction_output = gp_prediction_output,
    gp_factor_cache = gp_factor_cache,
    gp_iters = gp_iters,
    gp_lr = gp_lr,
    seed = seed
  )

  project_root <- gp_project_root()
  bridge_dir <- tempfile("predictpror_gp_bridge_")
  out_dir <- file.path(bridge_dir, "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bundle <- gp_bridge_prepare_bundle(
    pheno_df = pheno_df,
    gmatrix = gmatrix,
    train_idx0 = split$train_idx0,
    test_idx0 = split$test_idx0,
    dir_path = file.path(bridge_dir, "input")
  )
  cli_args <- c(
    "--project-root", project_root,
    "--pheno-csv", bundle$pheno_csv,
    "--kernel-bin", bundle$kernel_bin,
    "--kernel-meta", bundle$kernel_meta,
    "--geno-ids", bundle$geno_ids,
    "--train-idx", bundle$train_idx,
    "--test-idx", bundle$test_idx,
    "--out-dir", out_dir,
    "--gid-col", cols$gid,
    "--env-col", cols$env,
    "--y-col", cols$value,
    "--method", params$method,
    "--backend", params$backend,
    "--output-level", params$output_level,
    "--prediction-output", params$prediction_output,
    "--varcomp-mode", params$varcomp_mode,
    "--fa-rank", as.character(as.integer(params$fa_rank)),
    "--seed", as.character(as.integer(params$seed))
  )
  if (!is.null(params$gp_iters)) {
    cli_args <- c(cli_args, "--gp-iters", as.character(as.integer(params$gp_iters)))
  }
  if (!is.null(params$gp_lr)) {
    cli_args <- c(cli_args, "--gp-lr", as.character(as.numeric(params$gp_lr)))
  }
  if (!is.null(params$grm_factor_cache_root) && nzchar(params$grm_factor_cache_root)) {
    cli_args <- c(cli_args, "--grm-factor-cache-root", params$grm_factor_cache_root)
  }
  gp_bridge_run_cli("fit-predict", cli_args, python_bin = python_bin)
  pred_path <- file.path(out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("GP Python bridge did not write predictions.csv.", call. = FALSE)
  }
  pred_df <- utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (isTRUE(getOption("PredictProR.force_gc", FALSE))) {
    invisible(gc(verbose = FALSE))
  }
  out_pred <- if (identical(params$output_level, "predict_only")) {
    gp_backend_light_prediction_output(
      pred_df = pred_df,
      pheno_df = pheno_df,
      gen_name = gen_name,
      test_mask = split$test_mask
    )
  } else {
    gp_bridge_collect_prediction(
      pred_df = pred_df,
      pheno_df = pheno_df,
      response_name = response,
      gen_name = gen_name,
      test_mask = split$test_mask
    )
  }

  model_parameters <- data.frame(
    stat = c("gp_model", "gp_backend", "gp_output_level", "gp_varcomp_mode"),
    summary = c(gp_backend_resolve_model(model_name), params$backend, params$output_level, params$varcomp_mode),
    stringsAsFactors = FALSE
  )
  prediction_error_var <- if ("Prediction_error_variance" %in% names(out_pred)) {
    out_pred$Prediction_error_variance
  } else if ("PEV" %in% names(out_pred)) {
    out_pred$PEV
  } else {
    rep(NA_real_, nrow(out_pred))
  }

  raw_result <- list(
    framework = gp_framework_module_name(),
    method = params$method,
    backend = params$backend,
    execution = "python_subprocess",
    bridge_dir = normalizePath(bridge_dir, winslash = "/", mustWork = FALSE)
  )

  total_pred <- gp_bridge_total_prediction_output(out_pred, env_col = heter_groups)

  result <- list(
    model_parameters = model_parameters,
    predicted_values = out_pred,
    Total_Predicted_value = total_pred,
    raw_python_result = raw_result,
    diagnostic_plots = if (identical(params$output_level, "predict_only")) {
      NULL
    } else {
      tryCatch(
        diagnostic_plot_true_prediction(
          predictions = out_pred$Predicted_value,
          standard_errors = out_pred$Standard_error,
          prediction_error_var = prediction_error_var,
          observed_y = pheno_df[[cols$value]][!split$test_mask],
          GID_names = as.character(pheno_df[[cols$gid]]),
          train_test_label = ifelse(split$test_mask, "Test", "Train"),
          model_for_CI_cal = "ML",
          system_database = TRUE
        ),
        error = function(e) NULL
      )
    }
  )
  if (is.null(total_pred)) {
    result[["Total_Predicted_value"]] <- NULL
  }
  result
}

gp_cv_pick_numeric_column <- function(x, candidates) {
  if (is.null(x) || !length(candidates)) {
    return(NULL)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  col <- candidates[candidates %in% names(x)][1L]
  if (is.na(col) || !nzchar(col)) {
    return(NULL)
  }
  suppressWarnings(as.numeric(x[[col]]))
}

gp_cv_align_public_predictions <- function(predictions,
                                           pheno_data,
                                           tst,
                                           gen_name,
                                           heter_groups = NULL) {
  pred <- as.data.frame(predictions, stringsAsFactors = FALSE)
  if (!nrow(pred)) {
    return(pred)
  }
  if (nrow(pred) == length(tst)) {
    return(pred)
  }

  gid_col <- c("GID", "gid", gen_name, "Name", "name")
  gid_col <- gid_col[gid_col %in% names(pred)][1L]
  if (is.na(gid_col) || !nzchar(gid_col)) {
    return(pred[seq_len(min(nrow(pred), length(tst))), , drop = FALSE])
  }
  pred_gid <- as.character(pred[[gid_col]])
  test_gid <- as.character(pheno_data[[gen_name]][tst])

  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    env_col <- c("Env", "env", heter_groups)
    env_col <- env_col[env_col %in% names(pred)][1L]
    if (!is.na(env_col) && nzchar(env_col)) {
      pred_key <- paste(pred_gid, as.character(pred[[env_col]]), sep = "\r")
      test_key <- paste(test_gid, as.character(pheno_data[[heter_groups]][tst]), sep = "\r")
      idx <- match(test_key, pred_key)
      if (!anyNA(idx)) {
        return(pred[idx, , drop = FALSE])
      }
    }
  }

  idx <- match(test_gid, pred_gid)
  if (!anyNA(idx)) {
    return(pred[idx, , drop = FALSE])
  }
  pred[seq_len(min(nrow(pred), length(tst))), , drop = FALSE]
}

gp_cv_point_prediction_contract <- function(params = list()) {
  if (is.null(params)) {
    params <- list()
  }
  params$gp_output_level <- "predict_only"
  params$gp_return_se <- FALSE
  params$gp_full_vc <- FALSE
  params$gp_prediction_output <- "test_only"
  params$gp_return_trait_correlations <- FALSE
  params$return_trait_correlations <- FALSE
  params
}

gp_bridge_cv_predict <- function(model_name,
                                 y,
                                 tst,
                                 additional_params = NULL) {
  params <- gp_cv_point_prediction_contract(additional_params %||% list())
  pheno_data <- params$pheno_data %||% NULL
  response <- params$response %||% NULL
  gen_name <- params$gen_name %||% NULL
  heter_groups <- params$heter_groups %||% NULL
  gmatrix <- params$gmatrix_model_ready %||% NULL
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    omic1_kernel = params$omic1_kernel_model_ready %||% NULL,
    omic2_kernel = params$omic2_kernel_model_ready %||% NULL,
    omic3_kernel = params$omic3_kernel_model_ready %||% NULL,
    kernel_list = params$kernel_list %||% NULL
  )
  if (length(kernel_inputs) > 0L) {
    gmatrix <- kernel_inputs[[1L]]
  }
  gp_input_kernel_weights <- if (length(kernel_inputs) > 0L) {
    gp_public_kernel_weights(
      list(Ks = kernel_inputs, K = kernel_inputs[[1L]]),
      kernel_weights = params$lowrank_kernel_weights %||%
        params$kernel_weights %||%
        params$w_g %||%
        NULL
    )
  } else {
    numeric()
  }
  additional_kernel_inputs <- if (length(kernel_inputs) > 1L) {
    kernel_inputs[-1L]
  } else {
    NULL
  }
  if (is.null(pheno_data) || is.null(response) || is.null(gen_name) || is.null(gmatrix)) {
    stop("GP bridge CV requires pheno_data, response, gen_name, and at least one kernel.", call. = FALSE)
  }
  ph <- pheno_data
  ph[[response]] <- as.numeric(y)
  fit <- gp_single_trait_model(
    pheno_data = ph,
    gmatrix = gmatrix,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    test_set = NULL,
    env_similarity = params$env_similarity %||% NULL,
    env_ids = params$env_ids %||% NULL,
    env_covariates = params$env_covariates %||% NULL,
    reaction_norm_feature_qc = params$reaction_norm_feature_qc %||% TRUE,
    kenv_kernel = params$kenv_kernel %||% "matern32",
    kenv_bandwidth = params$kenv_bandwidth %||% 1.0,
    kenv_kernel_kwargs = params$kenv_kernel_kwargs %||% NULL,
    kernel_list = additional_kernel_inputs,
    w_g = unname(gp_input_kernel_weights),
    observation_weights = params$weights %||% NULL,
    model_name = model_name,
    tune_gp = FALSE,
    mean_adjustment = "none",
    prediction_blend = "none",
    gp_backend = params$gp_backend %||% "auto",
    gp_output_level = "predict_only",
    gp_return_se = FALSE,
    gp_uncertainty_calibration = "none",
    gp_full_vc = FALSE,
    gp_factor_cache = params$gp_factor_cache %||% NULL,
    gp_varcomp_mode = params$gp_varcomp_mode %||% "reml",
    gp_fa_rank = params$gp_fa_rank %||% 1L,
    gp_prediction_output = params$gp_prediction_output %||% "test_only",
    gp_iters = params$gp_iters %||% NULL,
    gp_lr = params$gp_lr %||% NULL,
    seed = params$random_state %||% 12345L
  )
  pred_df <- gp_cv_align_public_predictions(
    predictions = fit$predictions,
    pheno_data = ph,
    tst = as.integer(tst),
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  out <- gp_cv_pick_numeric_column(pred_df, c("Prediction", "Predicted_value", "yhat"))
  if (is.null(out)) {
    stop("GP bridge CV did not return a recognized prediction column.", call. = FALSE)
  }
  attr(out, "gp_info") <- fit$info
  attr(out, "gp_predictions") <- pred_df
  out
}

gp_backend_cv_predict_batch <- function(model_name,
                                        y,
                                        folds,
                                        additional_params = NULL) {
  if (!gp_backend_cv_batch_enabled()) {
    stop("Batched GP CV is disabled by PREDICTPRO_GP_BATCH_CV.", call. = FALSE)
  }
  params <- gp_cv_point_prediction_contract(additional_params %||% list())
  pheno_data <- params$pheno_data %||% NULL
  response <- params$response %||% NULL
  gen_name <- params$gen_name %||% NULL
  heter_groups <- params$heter_groups %||% NULL
  gmatrix <- params$gmatrix_model_ready %||% NULL
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    omic1_kernel = params$omic1_kernel_model_ready %||% NULL,
    omic2_kernel = params$omic2_kernel_model_ready %||% NULL,
    omic3_kernel = params$omic3_kernel_model_ready %||% NULL,
    kernel_list = params$kernel_list %||% NULL
  )
  if (is.null(gmatrix) && length(kernel_inputs) > 0L) {
    gmatrix <- kernel_inputs[[1L]]
  }
  if (is.null(pheno_data) || is.null(response) || is.null(gen_name) || is.null(gmatrix)) {
    stop("Batched GP bridge CV requires pheno_data, response, gen_name, and at least one kernel.", call. = FALSE)
  }
  if (!is.list(folds) || !length(folds)) {
    stop("Batched GP bridge CV requires a non-empty list of fold test-row indices.", call. = FALSE)
  }
  n <- nrow(pheno_data)
  y <- as.numeric(y)
  if (length(y) != n) {
    stop("Batched GP bridge CV y length does not match pheno_data rows.", call. = FALSE)
  }
  folds <- lapply(folds, function(tst) {
    tst <- as.integer(tst)
    tst <- tst[is.finite(tst) & !is.na(tst)]
    tst <- unique(tst)
    if (any(tst < 1L | tst > n)) {
      stop("Batched GP bridge CV fold contains row indices outside pheno_data.", call. = FALSE)
    }
    tst
  })
  folds <- folds[lengths(folds) > 0L]
  if (!length(folds)) {
    stop("Batched GP bridge CV has no non-empty folds.", call. = FALSE)
  }

  ph <- pheno_data
  ph[[response]] <- y
  stage2_weights <- gp_resolve_stage2_precision_weights(
    weights = params$weights %||% NULL,
    pheno_data = ph,
    response = response,
    gen_name = gen_name,
    context = "GP cross-validation Stage 2 observation weights"
  )
  pheno_ids <- as.character(ph[[gen_name]])
  kernel <- gp_public_kernel_bank_inputs(
    if (length(kernel_inputs) > 0L) kernel_inputs else gmatrix,
    pheno_ids = pheno_ids,
    geno_ids = params$geno_ids %||% NULL,
    allow_factor_cache = !is.null(params$gp_factor_cache %||% NULL)
  )
  gp_input_kernel_weights <- gp_public_kernel_weights(
    kernel,
    kernel_weights = params$lowrank_kernel_weights %||%
      params$kernel_weights %||%
      params$w_g %||%
      NULL
  )
  backend_schema <- gp_backend_schema_cols()
  single_cols <- backend_schema$single
  backend_gid_col <- single_cols$gid
  backend_env_col <- single_cols$env
  backend_y_col <- single_cols$y
  pheno_df <- data.frame(
    .gid = pheno_ids,
    .env = if (!is.null(heter_groups)) as.character(ph[[heter_groups]]) else "ENV1",
    .y = y,
    stringsAsFactors = FALSE
  )
  names(pheno_df) <- unname(c(backend_gid_col, backend_env_col, backend_y_col))
  finite_y <- is.finite(y)
  all_rows <- seq_len(n)
  fold_specs <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    tst <- folds[[i]]
    train_rows <- all_rows[finite_y]
    train_rows <- setdiff(train_rows, tst)
    if (!length(train_rows)) {
      stop("Batched GP bridge CV fold ", i, " has no finite training rows.", call. = FALSE)
    }
    fold_specs[[i]] <- list(
      fold_id = paste0("fold", i),
      train_idx = as.integer(train_rows - 1L),
      test_idx = as.integer(tst - 1L)
    )
  }

  first_test <- folds[[1L]]
  first_test_mask <- rep(FALSE, n)
  first_test_mask[first_test] <- TRUE
  first_train_mask <- finite_y & !first_test_mask
  split_policy <- list(
    train_idx = as.integer(which(first_train_mask) - 1L),
    test_idx = as.integer(which(first_test_mask) - 1L),
    train_mask = first_train_mask,
    test_mask = first_test_mask
  )

  env_similarity <- params$env_similarity %||% NULL
  env_ids <- params$env_ids %||% NULL
  env_covariates <- params$env_covariates %||% NULL
  if (!is.null(env_similarity) && !is.null(env_covariates)) {
    stop("Pass exactly one of `env_similarity` or `env_covariates` for batched GP CV.", call. = FALSE)
  }
  gp_factor_cache <- params$gp_factor_cache %||% NULL
  large_n_threshold <- as.integer(params$large_n_threshold %||% 50000L)[1L]
  operator_tol <- as.numeric(params$operator_tol %||% 1e-5)[1L]
  operator_max_iter <- as.integer(params$operator_max_iter %||% 500L)[1L]
  operator_dtype_compute <- as.character(params$operator_dtype_compute %||% "float32")[1L]
  gp_dtype <- as.character(params$gp_dtype %||% "float64")[1L]
  gp_auto_policy <- params$gp_auto_policy %||% TRUE
  include_components <- params$include_components %||% NULL
  policy <- gp_single_trait_auto_policy(
    pheno_df = pheno_df,
    split = split_policy,
    heter_groups = heter_groups,
    gid_col = backend_gid_col,
    env_col = backend_env_col,
    env_similarity = env_similarity,
    env_covariates = env_covariates,
    gmatrix = kernel$K,
    gp_factor_cache = gp_factor_cache,
    geno_ids = kernel$geno_ids,
    gp_backend = params$gp_backend %||% "auto",
    gp_tiered_dispatch = params$gp_tiered_dispatch %||% FALSE,
    gp_dtype = gp_dtype,
    large_n_threshold = large_n_threshold,
    operator_tol = operator_tol,
    operator_max_iter = operator_max_iter,
    operator_dtype_compute = operator_dtype_compute,
    tune_gp = FALSE,
    tuning_objective = params$tuning_objective %||% "rmse",
    include_components = include_components,
    enabled = isTRUE(gp_auto_policy),
    resolve_device = FALSE
  )

  fixed_effects <- params$fixed_effects %||% NULL
  if (is.null(fixed_effects) && !is.null(heter_groups)) {
    all_envs <- unique(as.character(pheno_df[[backend_env_col]]))
    fold_train_has_all_envs <- vapply(fold_specs, function(fold) {
      train_rows <- as.integer(fold$train_idx) + 1L
      all(all_envs %in% unique(as.character(pheno_df[[backend_env_col]][train_rows])))
    }, logical(1L))
    fixed_effects <- if (all(fold_train_has_all_envs)) backend_env_col else character(0)
  }

  gp_se_hutchinson_probes <- as.integer(params$gp_se_hutchinson_probes %||% 32L)[1L]
  if (!is.finite(gp_se_hutchinson_probes) ||
      is.na(gp_se_hutchinson_probes) ||
      gp_se_hutchinson_probes < 1L) {
    gp_se_hutchinson_probes <- 32L
  }

  args <- list(
    pheno_df = pheno_df,
    gid_col = backend_gid_col,
    env_col = backend_env_col,
    y_col = backend_y_col,
    geno_ids = kernel$geno_ids,
    include_components = as.list(policy$include_components %||% "g"),
    fixed_effects = as.list(fixed_effects %||% character(0)),
    method = gp_backend_method_map(model_name),
    backend = policy$gp_backend,
    prediction_output = params$gp_prediction_output %||% "test_only",
    output_level = "predict_only",
    point_predictions_only = TRUE,
    return_se = FALSE,
    compute_ai_se = FALSE,
    standardize = "global",
    varcomp_mode = params$gp_varcomp_mode %||% "reml",
    krr_lam = as.numeric(params$krr_lam %||% 1e-3)[1L],
    krr_lams = params$krr_lams %||% "none",
    lam_select = as.character(params$lam_select %||% "fixed")[1L],
    dtype = policy$gp_dtype,
    seed = as.integer(params$random_state %||% 12345L),
    large_n_threshold = as.integer(large_n_threshold),
    operator_tol = as.numeric(policy$operator_tol),
    operator_max_iter = as.integer(policy$operator_max_iter),
    operator_dtype_compute = as.character(policy$operator_dtype_compute),
    n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    operator_n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    tiered_n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    tiered_dispatch = isTRUE(policy$gp_tiered_dispatch),
    tiered_verbose = isTRUE(params$gp_tiered_verbose %||% FALSE)
  )
  if (isTRUE(stage2_weights$supplied)) {
    args$obs_weights <- as.numeric(stage2_weights$precision)
    args$obs_weight_mode <- "inverse_variance"
    args$obs_weight_global_scale <- 1
  }
  args <- gp_bridge_add_gp_exact_controls(args, params)
  if (!is.null(gp_factor_cache)) {
    args$grm_factor_cache <- gp_factor_cache
  }
  if (!is.null(env_similarity)) {
    env_similarity <- as.matrix(env_similarity)
    storage.mode(env_similarity) <- "double"
    args$env_similarity <- unname(env_similarity)
    if (is.null(env_ids) && !is.null(rownames(env_similarity))) {
      env_ids <- rownames(env_similarity)
    }
    if (!is.null(env_ids)) {
      args$env_ids <- as.character(env_ids)
    }
  }
  if (!is.null(env_covariates)) {
    args$env_covariates <- gp_prepare_env_covariates(
      env_covariates,
      source_env_col = heter_groups,
      target_env_col = backend_env_col
    )
    args$reaction_norm_feature_qc <- isTRUE(params$reaction_norm_feature_qc %||% TRUE)
    args$kenv_kernel <- as.character(params$kenv_kernel %||% "matern32")[1L]
    args$kenv_bandwidth <- as.numeric(params$kenv_bandwidth %||% 1.0)[1L]
    if (!is.null(params$kenv_kernel_kwargs)) {
      args$kenv_kernel_kwargs <- params$kenv_kernel_kwargs
    }
  }
  args$w_g <- unname(gp_input_kernel_weights)
  if (!is.null(params$w_ge)) args$w_ge <- as.numeric(params$w_ge)
  if (!is.null(params$w_e)) args$w_e <- as.numeric(params$w_e)[1L]
  if (identical(gp_backend_method_map(model_name), "gp_icm_fa")) {
    args$fa_rank <- as.integer(params$gp_fa_rank %||% 1L)
  }
  gp_exact_fast_cv <- gp_normalize_exact_fast_cv_mode(params$gp_exact_fast_cv %||% NULL)
  if (!is.null(gp_exact_fast_cv)) {
    args$gp_exact_fast_cv <- gp_exact_fast_cv
  }
  keep_cv_diagnostics <- tolower(Sys.getenv("PREDICTPRO_GP_CV_DIAGNOSTICS", unset = "false")) %in%
    c("1", "true", "yes", "y")
  keep_prediction_attr <- keep_cv_diagnostics ||
    (tolower(Sys.getenv("PREDICTPRO_GP_CV_KEEP_PREDICTIONS_ATTR", unset = "false")) %in%
       c("1", "true", "yes", "y"))

  python_bin <- params$python_bin %||% gp_detect_gp_python()
  project_root <- params$project_root %||% gp_project_root()
  fit <- gp_with_gp_device_route(
    policy$device$route,
    gp_bridge_fit_mixed_model_cv_direct(
      args = args,
      kernel = kernel,
      folds = fold_specs,
      gp_factor_cache = gp_factor_cache,
      python_bin = python_bin,
      project_root = project_root,
      read_diagnostics = keep_cv_diagnostics
    )
  )
  pred_all <- as.data.frame(fit$predictions, stringsAsFactors = FALSE)
  if (!nrow(pred_all)) {
    stop("Batched GP bridge CV returned no predictions.", call. = FALSE)
  }
  out <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    tst <- folds[[i]]
    fold_id <- paste0("fold", i)
    pred_df <- pred_all
    if ("fold_id" %in% names(pred_df)) {
      pred_df <- pred_df[as.character(pred_df$fold_id) == fold_id, , drop = FALSE]
    }
    row_col <- c("row_index", "cv_row_index", "row")
    row_col <- row_col[row_col %in% names(pred_df)][1L]
    if (!is.na(row_col) && nzchar(row_col)) {
      pred_row <- suppressWarnings(as.integer(as.numeric(pred_df[[row_col]]) + 1L))
      idx <- match(tst, pred_row)
      if (!all(is.na(idx))) {
        keep <- idx[!is.na(idx)]
        pred_df <- pred_df[keep, , drop = FALSE]
      }
    } else {
      pred_df <- gp_cv_align_public_predictions(
        predictions = pred_df,
        pheno_data = ph,
        tst = tst,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
    }
    pred <- gp_cv_pick_numeric_column(pred_df, c("Prediction", "Predicted_value", "yhat"))
    if (is.null(pred)) {
      stop("Batched GP bridge CV did not return a recognized prediction column.", call. = FALSE)
    }
    attr(pred, "gp_info") <- fit$info
    if (isTRUE(keep_prediction_attr)) {
      attr(pred, "gp_predictions") <- pred_df
    }
    out[[i]] <- pred
  }
  attr(out, "gp_info") <- fit$info
  if (isTRUE(keep_prediction_attr)) {
    attr(out, "gp_predictions") <- pred_all
  }
  out
}

gp_backend_prepare_kernel_cache <- function(ids,
                                            gmatrix = NULL,
                                            omic1_kernel = NULL,
                                            omic2_kernel = NULL,
                                            omic3_kernel = NULL,
                                            kernel_list = NULL,
                                            kernel_weights = NULL) {
  list(
    ids = as.character(ids),
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    kernel_weights = kernel_weights
  )
}

gp_backend_gaussian_model <- function(model_name,
                                      pheno_data,
                                      response,
                                      gen_name,
                                      gmatrix = NULL,
                                      omic1_kernel = NULL,
                                      omic2_kernel = NULL,
                                      omic3_kernel = NULL,
                                      kernel_list = NULL,
                                      test_set = NULL,
                                      test_set_source = NULL,
                                      lowrank_eps_trace = 1e-6,
                                      lowrank_max_rank = NULL,
                                      lowrank_jitter = NULL,
                                      lowrank_noise_grid = NULL,
                                      lowrank_kernel_weights = NULL,
                                      system_database = FALSE,
                                      gp_backend_cache = NULL,
                                      heter_groups = NULL,
                                      gp_backend = "auto",
                                      gp_output_level = "predict_only",
                                      gp_varcomp_mode = "reml",
                                      gp_fa_rank = 1L,
                                      gp_prediction_output = "all",
                                      gp_factor_cache = NULL,
                                      env_similarity = NULL,
                                      env_ids = NULL,
                                      env_covariates = NULL,
                                      reaction_norm_feature_qc = TRUE,
                                      kenv_kernel = "matern32",
                                      kenv_bandwidth = 1.0,
                                      kenv_kernel_kwargs = NULL,
                                      gp_iters = NULL,
                                      gp_lr = NULL,
                                      gp_engine = "auto",
                                      gp_learn_scales = NULL,
                                      observation_weights = NULL) {
  gp_engine_env <- tolower(Sys.getenv("PREDICTPRO_GP_ENGINE", unset = ""))
  if (nzchar(gp_engine_env)) {
    gp_engine <- gp_engine_env
  }
  gp_engine <- match.arg(
    tolower(as.character(gp_engine[1L])),
    c("auto", "dense_v", "mme", "eigen")
  )
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  if (length(kernel_inputs) > 0L) {
    gmatrix <- kernel_inputs[[1L]]
  }
  gp_input_kernel_weights <- if (length(kernel_inputs) > 0L) {
    gp_public_kernel_weights(
      list(Ks = kernel_inputs, K = kernel_inputs[[1L]]),
      kernel_weights = lowrank_kernel_weights
    )
  } else {
    numeric()
  }
  additional_kernel_inputs <- if (length(kernel_inputs) > 1L) {
    kernel_inputs[-1L]
  } else {
    NULL
  }
  if (is.null(gmatrix) && length(kernel_inputs) == 0L) {
    stop("Direct GP bridge currently requires at least one relationship matrix or kernel.", call. = FALSE)
  }
  # GP_FA is not identifiable with only one environment: lambda^2 and psi
  # cannot be separated from one scalar genetic variance. Do not relabel a
  # Kernel-GBLUP fit as FA-GBLUP. Reject the unsupported request and let the
  # user select a structurally appropriate single-environment model.
  if (identical(gp_backend_method_map(model_name), "gp_icm_fa")) {
    n_env <- 1L
    if (!is.null(heter_groups) && length(heter_groups) == 1L &&
        is.character(heter_groups) && nzchar(heter_groups) &&
        heter_groups %in% names(pheno_data)) {
      env_vec <- pheno_data[[heter_groups]]
      env_vec <- env_vec[!is.na(env_vec) & nzchar(as.character(env_vec))]
      n_env <- length(unique(as.character(env_vec)))
    }
    if (n_env <= 1L) {
      stop(
        "FA-GBLUP requires multi-environment data with at least two distinct ",
        "levels in `heter_groups`; a single-environment FA covariance is not ",
        "identifiable. Use Kernel-GBLUP or Gaussian-Process-GBLUP for ",
        "single-environment prediction.",
        call. = FALSE
      )
    }
  }
  return_se <- gp_output_level %in% c("predict_with_se", "full_vc")

  # Phase 2.1 Stage C: auto-warm-start from the session theta cache. Default
  # ON; opt-out via PREDICTPRO_GP_WARM_START=0 or by setting
  # gp_warm_start = FALSE on the caller. The cache key includes the canonical
  # model name so distinct GP parameterizations do not share warm starts.
  warm_enabled <- tolower(Sys.getenv("PREDICTPRO_GP_WARM_START", unset = "1")) %in%
    c("1", "true", "yes", "on")
  warm_key <- NULL
  warm_theta <- NULL
  if (isTRUE(warm_enabled)) {
    warm_key <- gp_warm_start_cache_key(
      model_name      = gp_backend_resolve_model(model_name),
      response        = response,
      gen_name        = gen_name,
      heter_groups    = heter_groups,
      gp_varcomp_mode = gp_varcomp_mode,
      gmatrix         = gmatrix
    )
    warm_theta <- gp_warm_start_lookup(warm_key)
  }

  public_fit <- gp_single_trait_model(
    pheno_data = pheno_data,
    gmatrix = gmatrix,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    test_set = test_set,
    test_set_source = test_set_source,
    env_similarity = env_similarity,
    env_ids = env_ids,
    env_covariates = env_covariates,
    reaction_norm_feature_qc = reaction_norm_feature_qc,
    kenv_kernel = kenv_kernel,
    kenv_bandwidth = kenv_bandwidth,
    kenv_kernel_kwargs = kenv_kernel_kwargs,
    kernel_list = additional_kernel_inputs,
    tune_gp = FALSE,
    mean_adjustment = "none",
    prediction_blend = "none",
    model_name = model_name,
    gp_backend = gp_backend,
    gp_output_level = gp_output_level,
    gp_return_se = return_se,
    gp_uncertainty_calibration = "none",
    gp_full_vc = identical(gp_output_level, "full_vc"),
    gp_varcomp_mode = gp_varcomp_mode,
    gp_fa_rank = gp_fa_rank,
    gp_prediction_output = gp_prediction_output,
    gp_factor_cache = gp_factor_cache,
    gp_iters = gp_iters,
    gp_lr = gp_lr,
    gp_theta0_warm = warm_theta,
    gp_engine = gp_engine,
    gp_learn_scales = gp_learn_scales,
    w_g = unname(gp_input_kernel_weights),
    observation_weights = observation_weights
  )

  # Save the converged theta for the next compatible fit. The Python backend
  # surfaces it as gp_result$summary$theta with the layout tag identifying
  # the spec shape (fa_icm / gp_exact). We store unconditionally on a fresh
  # fit and only use it when the next fit's key matches -- so unrelated fits
  # never see each other's theta.
  fit_layout <- tryCatch(
    public_fit[["result"]][["summary"]][["theta_layout"]],
    error = function(e) NA_character_
  )
  if (isTRUE(warm_enabled) && !is.null(warm_key)) {
    fit_theta <- tryCatch(
      public_fit[["result"]][["summary"]][["theta"]],
      error = function(e) NULL
    )
    if (!is.null(fit_theta) && is.numeric(fit_theta)) {
      gp_warm_start_save(warm_key, fit_theta, layout = fit_layout)
    }
  }
  # gp_engine_used: which REML engine actually ran. Sourced from the Python
  # summary's theta_layout tag ("gp_exact" -> dense_v, "gp_exact_mme" -> mme,
  # "fa_icm" -> dense_v FA). This lets users verify their gp_engine request
  # was honoured rather than silently downgraded by NotImplementedError
  # fallback inside gp_framework.py:gp_exact_with_X.
  gp_engine_used <- if (is.null(fit_layout) || is.na(fit_layout)) {
    NA_character_
  } else if (identical(as.character(fit_layout), "gp_exact_eigen")) {
    "eigen"
  } else if (identical(as.character(fit_layout), "gp_exact_mme")) {
    "mme"
  } else if (identical(as.character(fit_layout), "gp_exact") ||
             identical(as.character(fit_layout), "fa_icm")) {
    "dense_v"
  } else {
    as.character(fit_layout)
  }
  predicted_values <- gp_model_execute_single_trait_predicted_values(
    predictions = public_fit[["predictions"]],
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    test_set = test_set,
    test_set_source = test_set_source,
    gp_result = public_fit[["result"]]
  )
  total_predicted_values <- gp_bridge_total_prediction_output(predicted_values, env_col = heter_groups)
  # Python method that ran for predictions ("krr_exact" / "gp_exact" /
  # "gp_icm_fa"). For canonical KRR (Kernel-GBLUP, Scalable-GBLUP) the
  # prediction step is multikernel_krr_posterior_with_se. In a dense,
  # unweighted single-environment fit, GCV selects lambda = sigma2_e/sigma2_g;
  # this same ratio drives prediction uncertainty and VC reporting, with no
  # separate REML post-pass. Therefore gp_engine has no effect on KRR. For
  # canonical GP (Gaussian-Process-GBLUP), the prediction step IS the AI-REML
  # loop in gp_exact_with_X, so gp_engine governs predictions and VCs.
  gp_prediction_method <- tryCatch(
    public_fit[["info"]][["method"]],
    error = function(e) NA_character_
  )
  if (is.null(gp_prediction_method)) gp_prediction_method <- NA_character_
  stage2_weighted <- !is.null(observation_weights)
  gp_varcomp_mode_effective <- as.character(
    public_fit[["result"]][["variance_component_method"]] %||%
      public_fit[["result"]][["diagnostics"]][["variance_component_method"]] %||%
      public_fit[["info"]][["varcomp_mode"]] %||%
      gp_varcomp_mode
  )[1L]
  krr_lambda <- tryCatch(
    as.numeric(public_fit[["result"]][["diagnostics"]][["lam_eff"]])[1L],
    error = function(e) NA_real_
  )
  krr_lambda_selection <- tryCatch(
    as.character(public_fit[["result"]][["diagnostics"]][["lam_select"]])[1L],
    error = function(e) NA_character_
  )
  gp_kernel_names <- names(kernel_inputs)
  if (is.null(gp_kernel_names) || any(!nzchar(gp_kernel_names))) {
    gp_kernel_names <- paste0("kernel", seq_along(kernel_inputs))
  }
  gp_kernel_count <- length(kernel_inputs)
  gp_kernel_weights <- if (gp_kernel_count > 0L) {
    paste0(
      gp_kernel_names,
      "=",
      format(unname(gp_input_kernel_weights), digits = 15L, trim = TRUE),
      collapse = ";"
    )
  } else {
    NA_character_
  }
  gp_learns_kernel_scales <- gp_backend_method_map(model_name) %in% c("gp_exact", "gp_icm_fa") &&
    if (is.null(gp_learn_scales)) identical(gp_output_level, "full_vc") else isTRUE(gp_learn_scales)
  gp_kernel_weight_mode <- if (is.null(lowrank_kernel_weights)) {
    "equal_default"
  } else {
    "user_supplied"
  }
  if (isTRUE(gp_learns_kernel_scales)) {
    gp_kernel_weight_mode <- paste0(gp_kernel_weight_mode, "_initial_then_backend_learned_scales")
  } else {
    gp_kernel_weight_mode <- paste0(gp_kernel_weight_mode, "_fixed")
  }
  gp_is_multi_environment <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    heter_groups %in% names(pheno_data) &&
    length(unique(as.character(pheno_data[[heter_groups]][!is.na(pheno_data[[heter_groups]])]))) > 1L
  kernel_variance_components <- gp_public_single_kernel_variance_contract(
    gp_result = public_fit[["result"]],
    kernel_names = gp_kernel_names,
    model_name = model_name,
    independently_estimated = isTRUE(gp_learns_kernel_scales),
    is_met = gp_is_multi_environment
  )
  gp_fit_converged <- as.logical(
    public_fit[["result"]][["summary"]][["converged"]] %||%
      public_fit[["result"]][["diagnostics"]][["converged"]] %||%
      NA
  )[1L]
  variance_component_status <- if (isTRUE(gp_learns_kernel_scales) &&
                                   identical(gp_fit_converged, FALSE)) {
    "not_promoted_nonconverged"
  } else {
    "promoted"
  }
  fitted_kernel_scales <- if (is.data.frame(kernel_variance_components) &&
                              nrow(kernel_variance_components)) {
    report_rows <- as.logical(kernel_variance_components$Independent_kernel_estimate)
    if (any(report_rows)) {
      paste0(
        kernel_variance_components$Kernel[report_rows], "=",
        format(kernel_variance_components$Components[report_rows], digits = 15L, trim = TRUE),
        collapse = ";"
      )
    } else {
      NA_character_
    }
  } else {
    NA_character_
  }
  model_parameters <- data.frame(
    stat = c("gp_model", "gp_backend", "gp_output_level", "gp_varcomp_mode",
             "gp_varcomp_mode_requested", "gp_engine_requested", "gp_engine_used",
             "gp_prediction_method", "gp_krr_lambda", "gp_krr_lambda_selection",
             "gp_kernel_count", "gp_kernel_names", "gp_kernel_input_weights",
             "gp_kernel_weight_mode", "gp_kernel_fitted_variance_scales",
             "gp_kernel_variance_estimand",
             "stage2_observation_weighting"),
    summary = c(
      gp_backend_resolve_model(model_name),
      public_fit[["info"]][["backend"]] %||% gp_backend,
      gp_output_level,
      gp_varcomp_mode_effective,
      gp_varcomp_mode,
      as.character(gp_engine[1L]),
      gp_engine_used,
      as.character(gp_prediction_method),
      if (is.finite(krr_lambda)) format(krr_lambda, digits = 15L) else NA_character_,
      krr_lambda_selection,
      as.character(gp_kernel_count),
       paste(gp_kernel_names, collapse = ";"),
       gp_kernel_weights,
       gp_kernel_weight_mode,
       fitted_kernel_scales,
       if (isTRUE(gp_learns_kernel_scales)) {
         "independently_estimated_kernel_variance_components"
       } else if (identical(gp_backend_resolve_model(model_name), "KRR")) {
         "kernel_variance_not_identifiable_from_fixed_weights"
       } else {
         "fixed_kernel_weight_decomposition_of_combined_kernel_variance"
       },
       if (stage2_weighted) "precision_fixed_residual_variance" else "not_supplied"
    ),
    stringsAsFactors = FALSE
  )
  model_parameters <- rbind(
    model_parameters,
    gp_multi_kernel_parameter_rows(
      kernel_names = gp_kernel_names,
      strategy = paste0("GP_backend_", gp_prediction_method)
    )
  )
  list(
    model_parameters = model_parameters,
    predicted_values = predicted_values,
    Total_Predicted_value = total_predicted_values,
    raw_python_result = public_fit[["info"]],
     gp_result = public_fit[["result"]],
     gp_info = public_fit[["info"]],
     kernel_variance_components = kernel_variance_components,
     variance_component_status = variance_component_status,
     diagnostic_plots = if (identical(gp_output_level, "predict_only")) {
      NULL
    } else {
      tryCatch(
        diagnostic_plot_true_prediction(
          predictions = predicted_values$Predicted_value,
          standard_errors = predicted_values$Standard_error,
          prediction_error_var = predicted_values$PEV,
          observed_y = predicted_values$Observed_value[predicted_values$Train_Test_Label == "Train"],
          GID_names = predicted_values$GID,
          train_test_label = predicted_values$Train_Test_Label,
          model_for_CI_cal = "ML",
          system_database = TRUE
        ),
        error = function(e) NULL
      )
    }
  )
}

gp_import_gp_module <- function(module_name, python_bin = NULL, project_root = NULL) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop(
      "No configured Python runtime was found for Gaussian-process models. ",
      "Run setup_predictgp_env() or set PREDICTPRO_GP_PYTHON.",
      call. = FALSE
    )
  }
  gp_require_reticulate("GP reticulate fallback import")
  gp_init_python_once(python_bin)
  reticulate::py_run_string(sprintf(
    paste(
      "import importlib, sys",
      "_gp_dir = r'%s'",
      "if _gp_dir not in sys.path:",
      "    sys.path.insert(0, _gp_dir)",
      "try:",
      "    sys.modules['mixed_model_gpu_large_scale_backend'] = importlib.import_module('gp_large_backend')",
      "except Exception:",
      "    pass",
      sep = "\n"
    ),
    project_root
  ))
  reticulate::import(module_name, delay_load = FALSE, convert = FALSE)
}

gp_public_kernel_inputs <- function(gmatrix, pheno_ids, geno_ids = NULL,
                                    allow_factor_cache = FALSE) {
  if (is.null(gmatrix)) {
    if (!isTRUE(allow_factor_cache)) {
      stop("A genomic relationship matrix must be supplied through gmatrix.", call. = FALSE)
    }
    if (is.null(geno_ids) || !length(geno_ids)) {
      stop("geno_ids must be supplied when gmatrix is omitted for a factor-cache GP fit.", call. = FALSE)
    }
    geno_ids <- as.character(geno_ids)
    missing_pheno <- setdiff(unique(as.character(pheno_ids)), geno_ids)
    if (length(missing_pheno)) {
      stop("pheno_data contains genotypes not present in geno_ids.", call. = FALSE)
    }
    K <- matrix(1, nrow = 1L, ncol = 1L)
    rownames(K) <- colnames(K) <- geno_ids[[1L]]
    return(list(K = K, geno_ids = geno_ids))
  }
  K <- as.matrix(gmatrix)
  storage.mode(K) <- "double"
  if (nrow(K) != ncol(K)) {
    stop("gmatrix must be square.", call. = FALSE)
  }

  matrix_ids <- rownames(K)
  if (is.null(geno_ids)) {
    if (!is.null(matrix_ids)) {
      geno_ids <- matrix_ids
    } else if (!is.null(colnames(K))) {
      geno_ids <- colnames(K)
    } else {
      geno_ids <- unique(as.character(pheno_ids))
    }
  }
  geno_ids <- as.character(geno_ids)
  if (length(geno_ids) != nrow(K)) {
    stop("geno_ids must have the same length as the dimensions of gmatrix.", call. = FALSE)
  }

  if (!is.null(rownames(K)) && !is.null(colnames(K))) {
    missing_ids <- setdiff(geno_ids, intersect(rownames(K), colnames(K)))
    if (length(missing_ids)) {
      stop("geno_ids contains values missing from gmatrix row/column names.", call. = FALSE)
    }
    K <- K[geno_ids, geno_ids, drop = FALSE]
  } else {
    rownames(K) <- colnames(K) <- geno_ids
  }

  missing_pheno <- setdiff(unique(as.character(pheno_ids)), geno_ids)
  if (length(missing_pheno)) {
    stop("pheno_data contains genotypes not present in gmatrix.", call. = FALSE)
  }
  K <- 0.5 * (K + t(K))
  list(K = K, geno_ids = geno_ids)
}

gp_public_kernel_bank_inputs <- function(kernels,
                                         pheno_ids,
                                         geno_ids = NULL,
                                         allow_factor_cache = FALSE) {
  if (!gp_is_kernel_list(kernels)) {
    primary <- gp_public_kernel_inputs(
      kernels,
      pheno_ids = pheno_ids,
      geno_ids = geno_ids,
      allow_factor_cache = allow_factor_cache
    )
    return(list(
      K = primary$K,
      Ks = list(G = primary$K),
      geno_ids = primary$geno_ids
    ))
  }

  kernel_names <- names(kernels)
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("K", seq_along(kernels))
  }
  kernel_names <- make.unique(gp_sanitize_kernel_name(kernel_names), sep = "_")

  primary <- gp_public_kernel_inputs(
    kernels[[1L]],
    pheno_ids = pheno_ids,
    geno_ids = geno_ids,
    allow_factor_cache = allow_factor_cache
  )
  aligned <- vector("list", length(kernels))
  names(aligned) <- kernel_names
  aligned[[1L]] <- primary$K
  if (length(kernels) > 1L) {
    for (i in 2:length(kernels)) {
      one <- gp_public_kernel_inputs(
        kernels[[i]],
        pheno_ids = pheno_ids,
        geno_ids = primary$geno_ids,
        allow_factor_cache = FALSE
      )
      aligned[[i]] <- one$K
    }
  }

  list(
    K = aligned[[1L]],
    Ks = aligned,
    geno_ids = primary$geno_ids
  )
}

gp_public_additive_kernel_bank <- function(kernel) {
  kernels <- kernel$Ks
  if (is.null(kernels) || !length(kernels)) {
    kernels <- list(A = kernel$K)
  }
  kernel_names <- names(kernels)
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("K", seq_along(kernels))
  }
  kernel_names <- make.unique(gp_sanitize_kernel_name(kernel_names), sep = "_")
  names(kernels) <- kernel_names

  additive_idx <- match("A", names(kernels), nomatch = 0L)
  if (additive_idx > 0L) {
    if (additive_idx != 1L) {
      kernels <- c(kernels[additive_idx], kernels[-additive_idx])
    }
  } else {
    names(kernels)[1L] <- "A"
  }

  kernel$Ks <- kernels
  kernel$K <- kernels[[1L]]
  kernel
}

gp_public_kernel_weights <- function(kernel, kernel_weights = NULL) {
  kernels <- kernel$Ks
  if (is.null(kernels) || !length(kernels)) {
    kernels <- list(A = kernel$K)
  }
  kernel_names <- names(kernels)
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("K", seq_along(kernels))
  }
  if (is.null(kernel_weights)) {
    weights <- rep(1, length(kernels))
  } else {
    weights <- suppressWarnings(as.numeric(kernel_weights))
    supplied_names <- names(kernel_weights)
    if (!is.null(supplied_names) && any(nzchar(supplied_names))) {
      if (any(!nzchar(supplied_names)) || anyDuplicated(supplied_names)) {
        stop("Named `kernel_weights` must have unique, non-empty names.", call. = FALSE)
      }
      missing <- setdiff(kernel_names, supplied_names)
      extra <- setdiff(supplied_names, kernel_names)
      if (length(missing) || length(extra)) {
        stop(
          "Named `kernel_weights` must match the aligned kernel names. Expected: ",
          paste(kernel_names, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      weights <- weights[match(kernel_names, supplied_names)]
    }
  }
  if (length(weights) != length(kernels)) {
    stop(
      "`kernel_weights` length must equal the number of aligned kernels (",
      length(kernels), ").",
      call. = FALSE
    )
  }
  if (any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("`kernel_weights` must be finite, non-negative, and include a positive value.", call. = FALSE)
  }
  stats::setNames(as.numeric(weights), kernel_names)
}

gp_public_kernel_reporting_names <- function(kernel_inputs, kernel) {
  kernel_names <- names(kernel_inputs)
  if (is.null(kernel_names) || length(kernel_names) != length(kernel$Ks) ||
      any(!nzchar(kernel_names))) {
    kernel_names <- names(kernel$Ks)
  }
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("kernel", seq_along(kernel$Ks))
  }
  make.unique(gp_sanitize_kernel_name(kernel_names), sep = "_")
}

gp_public_covariance_correlation <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  d <- sqrt(pmax(diag(x), 0))
  denom <- outer(d, d)
  out <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  ok <- is.finite(denom) & denom > 0
  out[ok] <- x[ok] / denom[ok]
  diag(out)[is.finite(d) & d > 0] <- 1
  out
}

gp_public_matrix_long <- function(x, kernel_name, component,
                                  kernel_weight, kernel_diagonal_mean,
                                  estimation_method, fit_converged = NA) {
  x <- as.matrix(x)
  row_labels <- rownames(x) %||% paste0("Trait", seq_len(nrow(x)))
  col_labels <- colnames(x) %||% paste0("Trait", seq_len(ncol(x)))
  grid <- expand.grid(
    Trait_row = row_labels,
    Trait_column = col_labels,
    stringsAsFactors = FALSE
  )
  grid$Kernel <- kernel_name
  grid$Component <- component
  grid$Covariance <- as.numeric(x[cbind(
    match(grid$Trait_row, row_labels),
    match(grid$Trait_column, col_labels)
  )])
  grid$Kernel_weight <- as.numeric(kernel_weight)
  grid$Kernel_diagonal_mean <- as.numeric(kernel_diagonal_mean)
  grid$Estimation_method <- estimation_method
  grid$Independent_kernel_estimate <- FALSE
  grid$Fit_converged <- as.logical(fit_converged)
  grid[, c(
    "Kernel", "Component", "Trait_row", "Trait_column", "Covariance",
    "Kernel_weight", "Kernel_diagonal_mean", "Estimation_method",
    "Independent_kernel_estimate", "Fit_converged"
  ), drop = FALSE]
}

gp_public_joint_multikernel_contract <- function(out, kernel, kernel_weights,
                                                 kernel_names,
                                                 varcomp_mode = "reml",
                                                 is_met = FALSE,
                                                 return_trait_correlations = FALSE) {
  if (!is.list(out) || !is.list(out$fit) || !length(kernel$Ks)) {
    return(out)
  }
  kernel_names <- as.character(kernel_names)
  if (length(kernel_names) != length(kernel$Ks)) {
    kernel_names <- names(kernel$Ks)
  }
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("kernel", seq_along(kernel$Ks))
  }
  kernel_names <- make.unique(gp_sanitize_kernel_name(kernel_names), sep = "_")
  weights <- suppressWarnings(as.numeric(out$fit$kernel_weights %||% kernel_weights))
  if (length(weights) != length(kernel$Ks) || any(!is.finite(weights))) {
    weights <- suppressWarnings(as.numeric(kernel_weights))
  }
  if (length(weights) != length(kernel$Ks) || any(!is.finite(weights))) {
    return(out)
  }
  names(weights) <- kernel_names
  diagonal_means <- vapply(kernel$Ks, function(one) {
    value <- mean(diag(as.matrix(one)), na.rm = TRUE)
    if (is.finite(value)) value else NA_real_
  }, numeric(1L))
  names(diagonal_means) <- kernel_names

  trait_labels <- as.character(out$info$trait_levels %||% character())
  as_trait_matrix <- function(value) {
    value <- tryCatch(as.matrix(value), error = function(e) NULL)
    if (is.null(value) || !nrow(value) || nrow(value) != ncol(value)) {
      return(NULL)
    }
    storage.mode(value) <- "double"
    if (length(trait_labels) == nrow(value)) {
      rownames(value) <- colnames(value) <- trait_labels
    }
    value
  }
  as_trait_matrix_list <- function(value) {
    if (is.null(value) || !is.list(value) || is.data.frame(value)) return(NULL)
    out_list <- lapply(value, as_trait_matrix)
    keep <- vapply(out_list, Negate(is.null), logical(1L))
    out_list <- out_list[keep]
    if (!length(out_list)) NULL else out_list
  }
  main_coefficient <- as_trait_matrix(
    out$fit$Sigma_G_response_scale %||% out$result$genetic_covariance
  )
  if (is.null(main_coefficient)) {
    return(out)
  }
  gxe_coefficient <- if (isTRUE(is_met)) {
    as_trait_matrix(out$fit$Sigma_GE_response_scale %||% out$result$gxe_covariance)
  } else {
    NULL
  }
  residual_covariance <- as_trait_matrix(
    out$fit$Sigma_eps_response_scale %||% out$result$residual_covariance
  )
  main_by_environment_coefficient <- if (isTRUE(is_met)) {
    as_trait_matrix_list(
      out$fit$Sigma_G_response_scale_by_environment %||%
        out$result$genetic_covariance_by_environment
    )
  } else {
    NULL
  }
  gxe_by_environment_coefficient <- if (isTRUE(is_met)) {
    as_trait_matrix_list(
      out$fit$Sigma_GE_response_scale_by_environment %||%
        out$result$gxe_covariance_by_environment
    )
  } else {
    NULL
  }

  fit_converged <- if (grepl("reml", tolower(as.character(varcomp_mode)[1L]))) {
    as.logical(out$fit$reml_converged %||% out$fit$converged %||% NA)[1L]
  } else {
    as.logical(out$fit$pcg_converged %||% NA)[1L]
  }
  estimation_method <- paste0(
    toupper(as.character(varcomp_mode)[1L]),
    "_shared_covariance_fixed_kernel_weight_decomposition"
  )
  variance_promotable <- !identical(fit_converged, FALSE)

  main_by_kernel <- stats::setNames(
    lapply(seq_along(weights), function(i) main_coefficient * weights[[i]]),
    kernel_names
  )
  gxe_by_kernel <- if (!is.null(gxe_coefficient)) {
    stats::setNames(
      lapply(seq_along(weights), function(i) gxe_coefficient * weights[[i]]),
      kernel_names
    )
  } else {
    NULL
  }
  if (!isTRUE(variance_promotable)) {
    main_by_kernel <- lapply(main_by_kernel, function(one) {
      one[] <- NA_real_
      one
    })
    if (length(gxe_by_kernel)) {
      gxe_by_kernel <- lapply(gxe_by_kernel, function(one) {
        one[] <- NA_real_
        one
      })
    }
    if (!is.null(residual_covariance)) residual_covariance[] <- NA_real_
  }
  total_main <- Reduce(`+`, main_by_kernel)
  total_gxe <- if (length(gxe_by_kernel)) Reduce(`+`, gxe_by_kernel) else NULL
  total_by_kernel <- stats::setNames(lapply(seq_along(weights), function(i) {
    main_by_kernel[[i]] + if (length(gxe_by_kernel)) gxe_by_kernel[[i]] else 0
  }), kernel_names)

  variance_rows <- lapply(seq_along(weights), function(i) {
    pieces <- list(data.frame(
      Kernel = kernel_names[[i]],
      Trait = rownames(main_by_kernel[[i]]),
      Component = "kernel_genetic_variance",
      Components = diag(main_by_kernel[[i]]),
      Standard_error = NA_real_,
      Estimation_method = estimation_method,
      Kernel_weight = weights[[i]],
      Kernel_diagonal_mean = diagonal_means[[i]],
      Independent_kernel_estimate = FALSE,
      Variance_estimand = "kernel_weight_times_shared_trait_covariance_diagonal",
      Fit_converged = fit_converged,
      stringsAsFactors = FALSE
    ))
    if (length(gxe_by_kernel)) {
      pieces[[2L]] <- data.frame(
        Kernel = kernel_names[[i]],
        Trait = rownames(gxe_by_kernel[[i]]),
        Component = "kernel_gxe_variance",
        Components = diag(gxe_by_kernel[[i]]),
        Standard_error = NA_real_,
        Estimation_method = estimation_method,
        Kernel_weight = weights[[i]],
        Kernel_diagonal_mean = diagonal_means[[i]],
        Independent_kernel_estimate = FALSE,
        Variance_estimand = "kernel_weight_times_shared_gxe_covariance_diagonal",
        Fit_converged = fit_converged,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, pieces)
  })
  covariance_long <- lapply(seq_along(weights), function(i) {
    pieces <- list(gp_public_matrix_long(
      main_by_kernel[[i]], kernel_names[[i]], "genetic_covariance",
      weights[[i]], diagonal_means[[i]], estimation_method, fit_converged
    ))
    if (length(gxe_by_kernel)) {
      pieces[[2L]] <- gp_public_matrix_long(
        gxe_by_kernel[[i]], kernel_names[[i]], "gxe_covariance",
        weights[[i]], diagonal_means[[i]], estimation_method, fit_converged
      )
      pieces[[3L]] <- gp_public_matrix_long(
        total_by_kernel[[i]], kernel_names[[i]], "total_genetic_covariance",
        weights[[i]], diagonal_means[[i]], estimation_method, fit_converged
      )
    }
    do.call(rbind, pieces)
  })

  out$fit$Combined_kernel_genetic_covariance_coefficient <- main_coefficient
  out$fit$Sigma_G_response_scale <- total_main
  out$result$combined_kernel_genetic_covariance_coefficient <- main_coefficient
  out$result$genetic_covariance <- total_main
  if (isTRUE(return_trait_correlations)) {
    out$result$genetic_correlation <- gp_public_covariance_correlation(total_main)
  }
  if (length(main_by_environment_coefficient)) {
    out$fit$Combined_kernel_genetic_covariance_coefficient_by_environment <-
      main_by_environment_coefficient
    main_by_environment <- lapply(
      main_by_environment_coefficient,
      function(one) {
        one <- one * sum(weights)
        if (!isTRUE(variance_promotable)) one[] <- NA_real_
        one
      }
    )
    out$fit$Sigma_G_response_scale_by_environment <- main_by_environment
    out$result$genetic_covariance_by_environment <- main_by_environment
  }
  if (!is.null(total_gxe)) {
    out$fit$Combined_kernel_gxe_covariance_coefficient <- gxe_coefficient
    out$fit$Sigma_GE_response_scale <- total_gxe
    out$result$combined_kernel_gxe_covariance_coefficient <- gxe_coefficient
    out$result$gxe_covariance <- total_gxe
    if (isTRUE(return_trait_correlations)) {
      out$result$gxe_correlation <- gp_public_covariance_correlation(total_gxe)
    }
  }
  if (length(gxe_by_environment_coefficient)) {
    out$fit$Combined_kernel_gxe_covariance_coefficient_by_environment <-
      gxe_by_environment_coefficient
    gxe_by_environment <- lapply(
      gxe_by_environment_coefficient,
      function(one) {
        one <- one * sum(weights)
        if (!isTRUE(variance_promotable)) one[] <- NA_real_
        one
      }
    )
    out$fit$Sigma_GE_response_scale_by_environment <- gxe_by_environment
    out$result$gxe_covariance_by_environment <- gxe_by_environment
  }
  if (!is.null(residual_covariance)) {
    out$fit$Sigma_eps_response_scale <- residual_covariance
    out$result$residual_covariance <- residual_covariance
    if (isTRUE(return_trait_correlations)) {
      out$result$residual_correlation <- gp_public_covariance_correlation(residual_covariance)
    }
  }
  out$result$kernel_variance_components <- do.call(rbind, variance_rows)
  out$result$Genetic_covariance_by_kernel <- total_by_kernel
  out$result$genetic_covariance_by_kernel_long <- do.call(rbind, covariance_long)
  out$result$kernel_combination <- list(
    kernel_names = kernel_names,
    kernel_weights = stats::setNames(weights, kernel_names),
    kernel_diagonal_means = stats::setNames(diagonal_means, kernel_names),
    covariance_basis = "weighted_sum_of_shared_trait_covariance_coefficients",
    independent_kernel_covariances_estimated = FALSE
  )
  out$result$variance_component_status <- if (isTRUE(variance_promotable)) {
    "promoted"
  } else {
    "not_promoted_nonconverged"
  }
  out$info$kernel_names <- kernel_names
  out$info$kernel_weights <- stats::setNames(weights, kernel_names)
  out$info$kernel_diagonal_means <- stats::setNames(diagonal_means, kernel_names)
  out$info$kernel_covariance_basis <- "weighted_sum_of_shared_trait_covariance_coefficients"
  out$info$independent_kernel_covariances_estimated <- FALSE
  out
}

gp_inferred_missing_response_test_source <- function(test_set_source) {
  source <- as.character(test_set_source %||% character())
  source <- source[!is.na(source)]
  any(source == "inferred_missing_response") ||
    any(grepl("inferred_missing_response", source, fixed = TRUE))
}

gp_public_test_mask <- function(gids, y, test_set = NULL, test_set_source = NULL) {
  missing_mask <- !is.finite(y)
  if (isTRUE(gp_inferred_missing_response_test_source(test_set_source)) &&
      any(missing_mask, na.rm = TRUE)) {
    return(missing_mask)
  }
  if (is.null(test_set)) {
    return(missing_mask)
  }
  test_ids <- unique(as.character(test_set))
  if (length(test_ids)) gids %in% test_ids else missing_mask
}

gp_public_train_test_idx <- function(pheno_df,
                                     gid_col = gp_backend_schema_col("single", "gid"),
                                     y_col = gp_backend_schema_col("single", "y"),
                                     test_set = NULL,
                                     test_set_source = NULL) {
  gids <- as.character(pheno_df[[gid_col]])
  y <- suppressWarnings(as.numeric(pheno_df[[y_col]]))
  test_mask <- gp_public_test_mask(
    gids = gids,
    y = y,
    test_set = test_set,
    test_set_source = test_set_source
  )
  train_mask <- !test_mask & is.finite(y)
  if (!any(train_mask)) {
    stop("No finite training rows are available for the GP model.", call. = FALSE)
  }
  list(
    train_idx = as.integer(which(train_mask) - 1L),
    test_idx = as.integer(which(test_mask) - 1L),
    test_mask = test_mask,
    train_mask = train_mask
  )
}

gp_public_long_traits <- function(pheno_data, response, gen_name, heter_groups = NULL) {
  cols <- gp_backend_schema_cols()$long
  if (is.null(response) || length(response) < 2L) {
    stop("At least two response columns are required for a multi-trait GP model.", call. = FALSE)
  }
  missing_cols <- setdiff(c(gen_name, response, heter_groups), names(pheno_data))
  if (length(missing_cols)) {
    stop("pheno_data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    .gid = rep(as.character(pheno_data[[gen_name]]), times = length(response)),
    .trait = rep(as.character(response), each = nrow(pheno_data)),
    .y = as.numeric(unlist(pheno_data[response], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  names(out) <- unname(c(cols$gid, cols$trait, cols$y))
  if (!is.null(heter_groups)) {
    out[[cols$env]] <- rep(as.character(pheno_data[[heter_groups]]), times = length(response))
    out <- out[unname(c(cols$gid, cols$env, cols$trait, cols$y))]
  }
  out
}

gp_public_factor_from_kernel <- function(gmatrix, eig_tol = 1e-8) {
  K <- 0.5 * (as.matrix(gmatrix) + t(as.matrix(gmatrix)))
  storage.mode(K) <- "double"
  eig <- eigen(K, symmetric = TRUE)
  vals <- pmax(as.numeric(eig$values), 0)
  tol <- max(max(vals, na.rm = TRUE) * eig_tol, .Machine$double.eps)
  keep <- which(vals > tol)
  if (!length(keep)) {
    stop("gmatrix has no positive eigenvalues after tolerance filtering.", call. = FALSE)
  }
  phi <- sweep(eig$vectors[, keep, drop = FALSE], 2L, sqrt(vals[keep]), `*`)
  rownames(phi) <- rownames(gmatrix)
  phi
}

gp_public_inmemory_factor_cache <- function(phi, python_bin = NULL, project_root = NULL) {
  gp_require_reticulate("GP in-memory factor cache")
  gp_import_gp_module("gp_framework", python_bin = python_bin, project_root = project_root)
  phis <- if (is.list(phi) && !is.data.frame(phi)) phi else list(phi)
  phis <- lapply(phis, function(one) {
    one <- unname(as.matrix(one))
    storage.mode(one) <- "double"
    one
  })
  n_rows <- vapply(phis, nrow, integer(1L))
  if (!length(phis) || any(n_rows != n_rows[[1L]])) {
    stop("All in-memory GP factor roots must have the same number of genotype rows.", call. = FALSE)
  }
  reticulate::py_run_string(
    paste(
      "import numpy as np",
      "import torch",
      "from mixed_model_gpu_large_scale_backend import GRMFactorCacheBase",
      "class PredictProRInMemoryGRMCache(GRMFactorCacheBase):",
      "    def __init__(self, phis):",
      "        self._phis = [torch.tensor(np.asarray(phi, dtype=np.float64), dtype=torch.float64) for phi in phis]",
      "        super().__init__(num_grms=len(self._phis))",
      "        self._n_geno = int(self._phis[0].shape[0])",
      "    def ranks(self):",
      "        return [int(phi.shape[1]) for phi in self._phis]",
      "    def n_geno(self):",
      "        return self._n_geno",
      "    def get_rows(self, r, row_index, device, dtype):",
      "        ri = row_index.detach().to('cpu').long()",
      "        return self._phis[int(r)].index_select(0, ri).to(device=device, dtype=dtype)",
      sep = "\n"
    )
  )
  reticulate::py$PredictProRInMemoryGRMCache(phis)
}

gp_public_factor_cache_kernel <- function(gp_factor_cache,
                                          ids,
                                          geno_ids,
                                          python_bin = NULL,
                                          project_root = NULL,
                                          dtype = "float64") {
  if (is.null(gp_factor_cache)) {
    stop("`gp_factor_cache` is required to materialize a tuning kernel.", call. = FALSE)
  }
  ids <- unique(as.character(ids))
  geno_ids <- as.character(geno_ids)
  if (!length(ids)) {
    stop("No genotype IDs were supplied for factor-cache tuning.", call. = FALSE)
  }
  idx <- match(ids, geno_ids)
  if (anyNA(idx)) {
    missing <- ids[is.na(idx)]
    stop(
      "Factor-cache tuning IDs are missing from `geno_ids`: ",
      paste(head(missing, 5L), collapse = ", "),
      if (length(missing) > 5L) " ..." else "",
      call. = FALSE
    )
  }
  if (gp_bridge_direct_execution_enabled()) {
    if (!is.list(gp_factor_cache)) {
      stop(
        "Direct GP factor-cache materialization requires `gp_factor_cache` ",
        "to be a serializable list specification.",
        call. = FALSE
      )
    }
    python_bin <- python_bin %||% gp_detect_gp_python()
    project_root <- project_root %||% gp_project_root()
    bridge_dir <- tempfile("predictpror_gp_factor_kernel_")
    input_dir <- file.path(bridge_dir, "input")
    out_dir <- file.path(bridge_dir, "out")
    dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    cache_json <- file.path(input_dir, "factor_cache.json")
    row_idx_path <- file.path(input_dir, "row_idx.txt")
    gp_bridge_write_json(gp_bridge_jsonable_factor_cache(gp_factor_cache), cache_json)
    gp_bridge_write_lines(as.integer(idx - 1L), row_idx_path)
    gp_bridge_run_cli(
      "materialize-factor-cache",
      c(
        "--project-root", project_root,
        "--factor-cache-json", normalizePath(cache_json, winslash = "/", mustWork = TRUE),
        "--row-idx", normalizePath(row_idx_path, winslash = "/", mustWork = TRUE),
        "--out-dir", normalizePath(out_dir, winslash = "/", mustWork = FALSE),
        "--dtype", as.character(dtype)[1L]
      ),
      python_bin = python_bin
    )
    K <- gp_bridge_read_matrix_bin(
      file.path(out_dir, "kernel.bin"),
      file.path(out_dir, "kernel_shape.txt")
    )
    storage.mode(K) <- "double"
    rownames(K) <- colnames(K) <- ids
    return(K)
  }
  fw <- gp_import_gp_module("gp_framework", python_bin = python_bin, project_root = project_root)
  has_helper <- tryCatch(
    reticulate::py_has_attr(fw, "materialize_geno_kernel_from_factor_cache"),
    error = function(e) FALSE
  )
  if (!isTRUE(has_helper)) {
    stop(
      "Python GP backend does not expose `materialize_geno_kernel_from_factor_cache`.",
      call. = FALSE
    )
  }
  K <- reticulate::py_to_r(fw$materialize_geno_kernel_from_factor_cache(
    grm_factor_cache = gp_factor_cache,
    row_index = as.integer(idx - 1L),
    dtype = as.character(dtype)[1L]
  ))
  K <- as.matrix(K)
  storage.mode(K) <- "double"
  rownames(K) <- colnames(K) <- ids
  K
}

gp_public_normalize_multitrait_uncertainty <- function(out) {
  normalize_pred <- function(pred) {
    if (is.null(pred)) {
      return(pred)
    }
    pred <- gp_plain_data_frame(pred)
    if ("SE" %in% names(pred) && !"Prediction_SE_observed" %in% names(pred)) {
      pred$Prediction_SE_observed <- suppressWarnings(as.numeric(pred$SE))
    }
    if ("SE_latent" %in% names(pred) && !"Prediction_SE_latent" %in% names(pred)) {
      pred$Prediction_SE_latent <- suppressWarnings(as.numeric(pred$SE_latent))
    }
    if ("PEV" %in% names(pred) && !"Prediction_Var_latent" %in% names(pred)) {
      pred$Prediction_Var_latent <- suppressWarnings(as.numeric(pred$PEV))
    }
    if (!"Prediction_Var_observed" %in% names(pred) && "Prediction_SE_observed" %in% names(pred)) {
      pred$Prediction_Var_observed <- suppressWarnings(as.numeric(pred$Prediction_SE_observed))^2
    }
    if (!"Prediction_Var_latent" %in% names(pred) && "Prediction_SE_latent" %in% names(pred)) {
      pred$Prediction_Var_latent <- suppressWarnings(as.numeric(pred$Prediction_SE_latent))^2
    }
    if (!"SE" %in% names(pred) && "Prediction_SE_observed" %in% names(pred)) {
      pred$SE <- suppressWarnings(as.numeric(pred$Prediction_SE_observed))
    }
    if (!"SE_latent" %in% names(pred) && "Prediction_SE_latent" %in% names(pred)) {
      pred$SE_latent <- suppressWarnings(as.numeric(pred$Prediction_SE_latent))
    }
    if (!"PEV" %in% names(pred) && "Prediction_Var_latent" %in% names(pred)) {
      pred$PEV <- suppressWarnings(as.numeric(pred$Prediction_Var_latent))
    }
    pred
  }
  if (!is.null(out$predictions)) {
    out$predictions <- normalize_pred(out$predictions)
  }
  if (!is.null(out$result) && !is.null(out$result$predictions)) {
    out$result$predictions <- normalize_pred(out$result$predictions)
    if (is.null(out$predictions)) {
      out$predictions <- out$result$predictions
    }
  } else if (!is.null(out$result) && !is.null(out$predictions)) {
    out$result$predictions <- out$predictions
  }
  out
}

gp_model_execute_multitrait_predicted_values <- function(predictions,
                                                         pheno_data,
                                                         response,
                                                         gen_name,
                                                         heter_groups = NULL) {
  pred <- as.data.frame(predictions, stringsAsFactors = FALSE)
  gid_col <- c("GID", "gid", gen_name, "Name", "name")
  gid_col <- gid_col[gid_col %in% names(pred)][1L]
  trait_col <- c("Trait", "trait")
  trait_col <- trait_col[trait_col %in% names(pred)][1L]
  pred_col <- c("Predicted_value", "Prediction", "prediction", "yhat")
  pred_col <- pred_col[pred_col %in% names(pred)][1L]
  if (is.na(gid_col) || is.na(trait_col) || is.na(pred_col)) {
    stop("GP multi-trait predictions are missing genotype, trait, or prediction columns.", call. = FALSE)
  }
  out <- data.frame(GID = as.character(pred[[gid_col]]), stringsAsFactors = FALSE)
  if (!is.null(heter_groups)) {
    env_col <- c(heter_groups, "Env", "env")
    env_col <- env_col[env_col %in% names(pred)][1L]
    if (!is.na(env_col) && nzchar(env_col)) {
      # Honor user-supplied heter_groups as the output column name.
      out[[as.character(heter_groups)]] <- as.character(pred[[env_col]])
    }
  }
  out$Trait <- as.character(pred[[trait_col]])
  out$Predicted_value <- suppressWarnings(as.numeric(pred[[pred_col]]))

  se <- gp_cv_pick_numeric_column(
    pred,
    c("Prediction_SE_latent", "SE_latent", "SE", "Prediction_SE", "Standard_error", "Prediction_SE_observed")
  )
  pev <- gp_cv_pick_numeric_column(
    pred,
    c("Prediction_Var_latent", "PEV", "Prediction_error_variance", "Prediction_Var_observed")
  )
  if (is.null(pev) && !is.null(se)) {
    pev <- se^2
  }
  out$Standard_error <- se %||% rep(NA_real_, nrow(out))
  out$PEV <- pev %||% rep(NA_real_, nrow(out))

  key <- paste(out$GID, out$Trait, sep = "\r")
  ph_long <- gp_public_long_traits(
    pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  env_output_col <- if (!is.null(heter_groups) && as.character(heter_groups) %in% names(out)) {
    as.character(heter_groups)
  } else if (!is.null(heter_groups) && "Env" %in% names(out)) {
    "Env"
  } else {
    NULL
  }
  if (!is.null(env_output_col)) {
    key <- paste(out$GID, out[[env_output_col]], out$Trait, sep = "\r")
    ph_key <- paste(ph_long$gid, ph_long$env, ph_long$trait, sep = "\r")
  } else {
    ph_key <- paste(ph_long$gid, ph_long$trait, sep = "\r")
  }
  obs <- ph_long$y[match(key, ph_key)]
  out$Observed_value <- suppressWarnings(as.numeric(obs))
  out$Train_Test_Label <- ifelse(is.finite(out$Observed_value), "Train", "Test")
  out$Train_Test_Label[!is.finite(out$Predicted_value)] <- NA_character_

  z <- stats::qnorm(0.975)
  can_bound <- is.finite(out$Predicted_value) & is.finite(out$Standard_error)
  out$lower_bound <- rep(NA_real_, nrow(out))
  out$upper_bound <- rep(NA_real_, nrow(out))
  out$lower_bound[can_bound] <- out$Predicted_value[can_bound] - z * out$Standard_error[can_bound]
  out$upper_bound[can_bound] <- out$Predicted_value[can_bound] + z * out$Standard_error[can_bound]
  out$Uncertainty <- out$upper_bound - out$lower_bound
  out$Uncertainty_remarks <- ifelse(
    is.finite(out$Uncertainty),
    "model_based_prediction_interval",
    NA_character_
  )

  finite_y <- is.finite(ph_long$y)
  trait_reference <- stats::setNames(rep(NA_real_, length(response)), response)
  for (trait in response) {
    trait_y <- suppressWarnings(as.numeric(ph_long$y[ph_long$trait == trait & finite_y]))
    trait_var <- suppressWarnings(stats::var(trait_y, na.rm = TRUE))
    if (is.finite(trait_var) && trait_var > 0) {
      trait_reference[[trait]] <- trait_var
    }
  }
  global_reference <- suppressWarnings(stats::var(ph_long$y[finite_y], na.rm = TRUE))
  if (!is.finite(global_reference) || global_reference <= 0) {
    global_reference <- NA_real_
  }
  reference_variance <- as.numeric(trait_reference[out$Trait])
  missing_reference <- !is.finite(reference_variance) | reference_variance <= 0
  reference_variance[missing_reference] <- global_reference
  out$Reliability <- rep(NA_real_, nrow(out))
  rel_ok <- is.finite(out$PEV) & is.finite(reference_variance) & reference_variance > 0
  out$Reliability[rel_ok] <- pmax(0, pmin(1, 1 - (out$PEV[rel_ok] / reference_variance[rel_ok])))
  out$Reliability_reference_variance <- reference_variance
  out$Reliability_remarks <- ifelse(
    is.na(out$Reliability),
    NA_character_,
    ifelse(out$Reliability >= 0.6, "Reliable", ifelse(out$Reliability <= 0.2, "Unreliable", "Acceptable"))
  )

  cols <- c(
    "GID",
    if (!is.null(env_output_col)) env_output_col,
    "Trait",
    "Predicted_value",
    "Train_Test_Label",
    "Observed_value",
    "Standard_error",
    "PEV",
    "lower_bound",
    "upper_bound",
    "Uncertainty",
    "Uncertainty_remarks",
    "Reliability",
    "Reliability_remarks"
  )
  out[, cols, drop = FALSE]
}

gp_model_execute_multitrait_prediction_wide <- function(pred_long,
                                                        gen_name,
                                                        heter_groups = NULL) {
  gid_col <- if (gen_name %in% names(pred_long)) gen_name else "GID"
  id_cols <- gid_col
  if (!is.null(heter_groups)) {
    env_col <- if (heter_groups %in% names(pred_long)) heter_groups else "Env"
    if (env_col %in% names(pred_long)) {
      id_cols <- c(id_cols, env_col)
    }
  }
  out <- stats::reshape(
    pred_long[, c(id_cols, "Trait", "Predicted_value"), drop = FALSE],
    idvar = id_cols,
    timevar = "Trait",
    direction = "wide"
  )
  names(out) <- sub("^Predicted_value\\.", "", names(out))
  out
}

gp_model_execute_multitrait_summary <- function(pred_long,
                                                response,
                                                model_type,
                                                mode) {
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
      mode,
      "gaussian",
      model_type,
      length(response),
      paste(response, collapse = ", "),
      "joint cross-trait gaussian-process true prediction",
      counts$observed_count,
      counts$predicted_count
    ),
    stringsAsFactors = FALSE
  )
}

gp_result_as_data_frame <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  out <- tryCatch(
    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(out) || !nrow(out)) {
    return(NULL)
  }
  out
}

gp_public_single_kernel_variance_contract <- function(gp_result,
                                                       kernel_names,
                                                       model_name = NULL,
                                                       independently_estimated = FALSE,
                                                       is_met = FALSE) {
  if (!is.list(gp_result)) {
    return(NULL)
  }
  summary <- gp_result_as_data_frame(
    gp_result[["var_components_summary"]] %||% gp_result[["var_components"]]
  )
  if (!is.data.frame(summary) || !nrow(summary) ||
      !all(c("kernel", "component_type", "estimate") %in% names(summary))) {
    return(NULL)
  }
  component_type <- tolower(as.character(summary[["component_type"]]))
  keep_types <- "genetic_main"
  if (isTRUE(is_met)) keep_types <- c(keep_types, "interaction")
  keep <- component_type %in% keep_types &
    tolower(as.character(summary[["kernel"]])) != "total"
  summary <- summary[keep, , drop = FALSE]
  if (!nrow(summary)) {
    return(NULL)
  }

  kernel_names <- as.character(kernel_names)
  kernel_names <- kernel_names[!is.na(kernel_names) & nzchar(kernel_names)]
  if (!length(kernel_names)) {
    kernel_names <- paste0("kernel", seq_len(sum(component_type == "genetic_main")))
  }
  map_group_names <- function(labels, type) {
    idx <- which(tolower(as.character(summary[["component_type"]])) == type)
    if (!length(idx)) return(labels)
    supplied <- as.character(summary[["kernel"]][idx])
    sanitized <- gp_sanitize_kernel_name(supplied)
    matched <- match(sanitized, gp_sanitize_kernel_name(kernel_names))
    generic <- grepl("^(G|K)[0-9]+$", supplied, ignore.case = TRUE)
    positional <- seq_along(idx)
    use_positional <- is.na(matched) | generic
    matched[use_positional & positional <= length(kernel_names)] <-
      positional[use_positional & positional <= length(kernel_names)]
    valid <- is.finite(matched) & matched >= 1L & matched <= length(kernel_names)
    labels[idx[valid]] <- kernel_names[matched[valid]]
    labels
  }
  public_kernel <- as.character(summary[["kernel"]])
  public_kernel <- map_group_names(public_kernel, "genetic_main")
  public_kernel <- map_group_names(public_kernel, "interaction")

  estimate <- suppressWarnings(as.numeric(summary[["estimate"]]))
  se_col <- c("std.error", "Standard_error", "standard_error", "se")
  se_col <- se_col[se_col %in% names(summary)][1L]
  standard_error <- if (!is.na(se_col)) {
    suppressWarnings(as.numeric(summary[[se_col]]))
  } else {
    rep(NA_real_, nrow(summary))
  }
  boundary <- rep("", nrow(summary))
  fit_converged <- as.logical(
    gp_result[["summary"]][["converged"]] %||%
      gp_result[["diagnostics"]][["converged"]] %||%
      NA
  )[1L]

  # Exact GP returns an ASReml-shaped table with independently estimated main
  # kernel scales and their AI-REML SEs. Prefer it over the rounded reporting
  # summary, preserving a genuine zero estimate while withholding its singular
  # boundary SE. KRR has one combined-kernel row, so it deliberately remains a
  # fixed-weight decomposition from the reporting summary.
  raw <- gp_result_as_data_frame(gp_result[["varcomp"]] %||% gp_result[["var_components_ai"]])
  if (is.data.frame(raw) && nrow(raw)) {
    component_col <- c("component", "Component")
    component_col <- component_col[component_col %in% names(raw)][1L]
    estimate_col <- c("estimate", "Estimate", "Components", "value")
    estimate_col <- estimate_col[estimate_col %in% names(raw)][1L]
    raw_se_col <- c("std.error", "Standard_error", "standard_error", "se")
    raw_se_col <- raw_se_col[raw_se_col %in% names(raw)][1L]
    if (!is.na(component_col) && !is.na(estimate_col)) {
      raw_component <- as.character(raw[[component_col]])
      raw_main <- grepl("^vm\\([^,]+,[^)]+\\)!var$", raw_component)
      raw_idx <- which(raw_main)
      if (length(raw_idx) == length(kernel_names)) {
        main_idx <- which(tolower(as.character(summary[["component_type"]])) == "genetic_main")
        for (j in seq_len(min(length(main_idx), length(raw_idx)))) {
          i <- main_idx[[j]]
          r <- raw_idx[[j]]
          estimate[[i]] <- suppressWarnings(as.numeric(raw[[estimate_col]][[r]]))
          if (!is.na(raw_se_col)) {
            standard_error[[i]] <- suppressWarnings(as.numeric(raw[[raw_se_col]][[r]]))
          }
          if ("bound" %in% names(raw)) {
            status <- toupper(trimws(as.character(raw[["bound"]][[r]])))
            if (status %in% c("B", "BOUND", "BOUNDARY")) boundary[[i]] <- "bound_at_zero"
          }
          if (identical(boundary[[i]], "bound_at_zero")) standard_error[[i]] <- NA_real_
        }
      }
    }
  }

  backend_method <- gp_backend_method_map(model_name)
  row_independent <- rep(FALSE, nrow(summary))
  if (isTRUE(independently_estimated)) {
    if (identical(backend_method, "gp_icm_fa")) {
      row_independent <- tolower(as.character(summary[["component_type"]])) == "interaction"
    } else {
      row_independent[] <- TRUE
    }
  }
  is_krr <- identical(gp_backend_resolve_model(model_name), "KRR")
  method <- ifelse(
    row_independent,
    "GP_REML_independent_kernel_variance",
    if (is_krr) {
      "KRR_GCV_combined_kernel_only"
    } else {
      "GP_fixed_kernel_weight_decomposition"
    }
  )
  out <- data.frame(
    Kernel = public_kernel,
    Component = ifelse(
      tolower(as.character(summary[["component_type"]])) == "interaction",
      "kernel_gxe_variance",
      "kernel_genetic_variance"
    ),
    Components = estimate,
    Standard_error = standard_error,
    Estimation_method = method,
    Independent_kernel_estimate = row_independent,
    Variance_estimand = ifelse(
      row_independent,
      "independently_estimated_random_effect_kernel_variance",
      "fixed_kernel_weight_share_of_combined_kernel_variance"
    ),
    Fit_converged = fit_converged,
    stringsAsFactors = FALSE
  )
  if (identical(fit_converged, FALSE) && any(row_independent)) {
    out$Components <- NA_real_
    out$Standard_error <- NA_real_
  }
  # One environment with the residual estimated at zero: the kernel genetic
  # variances absorb the residual and are not separable from it.
  if (!isTRUE(is_met) && !is_krr && is.data.frame(raw) && nrow(raw) &&
      any(gp_vc_residual_at_zero(raw))) {
    out$Components <- NA_real_
    out$Standard_error <- NA_real_
    out$Variance_estimand <- "not_estimable_residual_variance_at_zero"
    boundary[] <- "not_estimable"
  }
  if (is_krr) {
    # GCV fits the combined kernel; its fixed-weight split is not an estimate
    # of a separate genetic variance for each input kernel.
    out$Fixed_weight_allocation <- out$Components
    out$Component <- ifelse(
      out$Component == "kernel_gxe_variance",
      "kernel_gxe_variance_not_identifiable",
      "kernel_genetic_variance_not_identifiable"
    )
    out$Components <- NA_real_
    out$Standard_error <- NA_real_
    out$Variance_estimand <- "not_identifiable_from_fixed_kernel_weights"
  }
  if (any(nzchar(boundary))) out$Boundary <- boundary
  rownames(out) <- NULL
  out
}

gp_first_numeric_column <- function(x, candidates) {
  if (is.null(x) || !is.data.frame(x) || !length(candidates)) {
    return(NULL)
  }
  hit <- candidates[candidates %in% names(x)][1L]
  if (is.na(hit)) {
    return(NULL)
  }
  suppressWarnings(as.numeric(x[[hit]]))
}

gp_extract_gp_env_genetic_variance <- function(gp_result, n, env_values = NULL) {
  env_tables <- list(
    gp_result[["per_env"]],
    gp_result[["env_variance_summary"]]
  )
  for (tbl in env_tables) {
    df <- gp_result_as_data_frame(tbl)
    if (is.null(df)) {
      next
    }
    gvar <- gp_first_numeric_column(
      df,
      c("total_genetic_variance", "Genetic_variance", "genetic_variance", "genetic_main_variance")
    )
    if (is.null(gvar) || !any(is.finite(gvar) & gvar > 0)) {
      next
    }
    env_col <- c("Env", "env", "Environment", "environment", "level")
    env_col <- env_col[env_col %in% names(df)][1L]
    if (!is.na(env_col) && !is.null(env_values)) {
      matched <- gvar[match(as.character(env_values), as.character(df[[env_col]]))]
      if (any(is.finite(matched) & matched > 0)) {
        return(matched)
      }
    }
    finite_gvar <- gvar[is.finite(gvar) & gvar > 0]
    if (length(finite_gvar) == 1L) {
      return(rep(finite_gvar[[1L]], n))
    }
    if (length(gvar) == n) {
      return(gvar)
    }
    return(rep(mean(finite_gvar), n))
  }
  NULL
}

gp_extract_gp_component_genetic_variance <- function(gp_result, n) {
  component_tables <- list(
    gp_result[["var_components_summary"]],
    gp_result[["var_components"]],
    gp_result[["varcomp"]]
  )
  for (tbl in component_tables) {
    df <- gp_result_as_data_frame(tbl)
    if (is.null(df)) {
      next
    }
    estimate <- gp_first_numeric_column(df, c("estimate", "Estimate", "value"))
    if (is.null(estimate)) {
      next
    }
    keep <- rep(FALSE, nrow(df))
    if ("component_type" %in% names(df)) {
      component_type <- tolower(as.character(df[["component_type"]]))
      keep <- component_type %in% c("genetic_main", "interaction")
      if ("kernel" %in% names(df)) {
        kernel <- tolower(as.character(df[["kernel"]]))
        total_keep <- keep & kernel == "total"
        if (any(total_keep)) {
          keep <- total_keep
        }
      }
    } else if ("component" %in% names(df)) {
      component <- tolower(as.character(df[["component"]]))
      keep <- !grepl("!r$|resid|residual", component) &
        grepl("vm\\(|gid|geno|genetic|kernel", component)
    }
    vals <- estimate[keep & is.finite(estimate) & estimate > 0]
    if (length(vals)) {
      return(rep(sum(vals), n))
    }
  }
  NULL
}

# Genetic variance used as the reliability denominator. Where the residual was
# estimated at zero, genetic and residual variance are not separable, so the
# genetic variance (and therefore reliability) is not estimable: NA for the
# whole trait (one environment) or for the affected environments (MET).
gp_extract_gp_genetic_variance <- function(gp_result, n, env_values = NULL) {
  out <- gp_extract_gp_genetic_variance_raw(gp_result, n, env_values)
  if (!length(out) || is.null(gp_result) || !is.list(gp_result)) {
    return(out)
  }
  varcomp <- tryCatch(gp_vc_find_named_data_frame(gp_result, c("varcomp", "var_components_ai")),
                      error = function(e) NULL)
  at_zero <- if (is.data.frame(varcomp) && nrow(varcomp)) gp_vc_residual_at_zero(varcomp) else logical(0L)
  if (!any(at_zero)) {
    return(out)
  }
  single_environment <- is.null(env_values) ||
    length(unique(as.character(env_values[!is.na(env_values)]))) <= 1L
  if (isTRUE(single_environment)) {
    out[] <- NA_real_
    return(out)
  }
  comp_col <- gp_contract_first_name(names(varcomp), c("Component", "component", "term", "source"))
  labels <- if (!is.na(comp_col)) as.character(varcomp[[comp_col]]) else rownames(varcomp)
  not_estimable_env <- sub("^Env_", "", sub("!R$", "", labels[at_zero]))
  out[as.character(env_values) %in% not_estimable_env] <- NA_real_
  out
}

gp_extract_gp_genetic_variance_raw <- function(gp_result, n, env_values = NULL) {
  n <- as.integer(n %||% 0L)[1L]
  if (!is.finite(n) || is.na(n) || n < 1L) {
    return(numeric())
  }
  if (is.null(gp_result) || !is.list(gp_result)) {
    return(rep(NA_real_, n))
  }
  # For one environment, prefer the fitted component table.  The reporting
  # bundle can also contain the unit-scale kernel diagonal, which is useful
  # structural metadata but is not the fitted genetic variance and must not be
  # used as the reliability denominator.  Multi-environment fits still prefer
  # their environment-specific fitted summaries.
  single_environment <- is.null(env_values) ||
    length(unique(as.character(env_values[!is.na(env_values)]))) <= 1L
  # Keep the reliability denominator exactly aligned with the public
  # variance-components table. That formatter deliberately prefers the REML
  # varcomp total when a nearby component-summary estimate lacks its own SE,
  # while FA delta-method summaries with a finite SE remain authoritative.
  # Reading the same formatted row here prevents tiny but user-visible
  # discrepancies between Reliability_reference_variance and the exported
  # genetic_variance component.
  if (isTRUE(single_environment) &&
      exists("gp_format_gp_single_env_variance_components", mode = "function")) {
    public_vc <- tryCatch(
      gp_format_gp_single_env_variance_components(gp_result),
      error = function(e) NULL
    )
    if (is.data.frame(public_vc) && nrow(public_vc) > 0L &&
        all(c("Component", "Components") %in% names(public_vc))) {
      public_genetic <- suppressWarnings(as.numeric(
        public_vc[["Components"]][match("genetic_variance", public_vc[["Component"]])]
      ))
      if (length(public_genetic) == 1L && is.finite(public_genetic) &&
          public_genetic > 0) {
        return(rep(public_genetic, n))
      }
    }
  }
  component_gvar <- gp_extract_gp_component_genetic_variance(gp_result, n)
  env_gvar <- gp_extract_gp_env_genetic_variance(gp_result, n, env_values)
  if (isTRUE(single_environment) && !is.null(component_gvar)) {
    return(as.numeric(component_gvar))
  }
  if (!is.null(env_gvar)) {
    return(as.numeric(env_gvar))
  }
  if (!is.null(component_gvar)) {
    return(as.numeric(component_gvar))
  }
  rep(NA_real_, n)
}

gp_model_execute_single_trait_predicted_values <- function(predictions,
                                                           pheno_data,
                                                           response,
                                                           gen_name,
                                                           heter_groups = NULL,
                                                           test_set = NULL,
                                                           test_set_source = NULL,
                                                           gp_result = NULL) {
  pred <- as.data.frame(predictions, stringsAsFactors = FALSE)
  gid_col <- c("GID", "gid", gen_name, "Name", "name")
  gid_col <- gid_col[gid_col %in% names(pred)][1L]
  pred_col <- c("Predicted_value", "Prediction", "prediction", "yhat")
  pred_col <- pred_col[pred_col %in% names(pred)][1L]
  if (is.na(gid_col) || is.na(pred_col)) {
    stop("GP predictions are missing genotype or prediction columns.", call. = FALSE)
  }
  response_y_full <- suppressWarnings(as.numeric(pheno_data[[response]]))
  pheno_ids_full <- as.character(pheno_data[[gen_name]])
  test_mask_full <- gp_public_test_mask(
    gids = pheno_ids_full,
    y = response_y_full,
    test_set = test_set,
    test_set_source = test_set_source
  )
  target_rows <- if (nrow(pred) == nrow(pheno_data)) {
    seq_len(nrow(pheno_data))
  } else {
    which(test_mask_full)
  }
  if (nrow(pred) != nrow(pheno_data)) {
    pred <- gp_cv_align_public_predictions(
      predictions = pred,
      pheno_data = pheno_data,
      tst = target_rows,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  }
  se <- gp_cv_pick_numeric_column(
    pred,
    c("Prediction_SE_latent", "SE_latent", "Standard_error", "Prediction_SE", "SE", "Prediction_SE_observed")
  )
  pev <- gp_cv_pick_numeric_column(
    pred,
    c("Prediction_Var_latent", "PEV", "Prediction_error_variance", "Prediction_Var_observed")
  )
  if (is.null(pev) && !is.null(se)) {
    pev <- se^2
  }
  ids <- as.character(pred[[gid_col]])
  response_y <- response_y_full[target_rows]
  test_mask <- test_mask_full[target_rows]
  se_out <- se %||% rep(NA_real_, length(ids))
  pev_out <- pev %||% rep(NA_real_, length(ids))
  env_values <- if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    as.character(pheno_data[[heter_groups]][target_rows])
  } else {
    pred_env_col <- c("Env", "env", "Environment", "environment")
    pred_env_col <- pred_env_col[pred_env_col %in% names(pred)][1L]
    if (!is.na(pred_env_col)) as.character(pred[[pred_env_col]]) else rep("ENV1", length(ids))
  }
  genetic_variance_out <- gp_extract_gp_genetic_variance(
    gp_result = gp_result,
    n = length(ids),
    env_values = env_values
  )

  has_genetic_variance <- any(is.finite(genetic_variance_out) & genetic_variance_out > 0)
  reference_variance_out <- if (isTRUE(has_genetic_variance)) {
    genetic_variance_out
  } else {
    rep(NA_real_, length(ids))
  }
  reliability <- rep(NA_real_, length(ids))
  rel_valid <- is.finite(pev_out) & is.finite(reference_variance_out) & reference_variance_out > 0
  reliability[rel_valid] <- pmax(0, pmin(1, 1 - (pev_out[rel_valid] / reference_variance_out[rel_valid])))
  reliability_remarks <- if (isTRUE(has_genetic_variance)) {
    ifelse(
      is.na(reliability),
      ifelse(is.finite(genetic_variance_out), NA_character_, "sigma2_g not estimable"),
      ifelse(reliability >= 0.6, "Reliable", ifelse(reliability <= 0.2, "Unreliable", "Acceptable"))
    )
  } else {
    rep("sigma2_g not estimable", length(ids))
  }
  reliability_percentage <- if (any(is.finite(reliability))) {
    mean(reliability[is.finite(reliability)] > 0.2) * 100
  } else {
    NA_real_
  }

  uncertainty_remarks_out <- ifelse(
    is.finite(se_out),
    "model_based_prediction_interval",
    NA_character_
  )

  variance_component_method <- as.character(
    gp_result[["variance_component_method"]] %||%
      gp_result[["diagnostics"]][["variance_component_method"]] %||%
      "reml"
  )[1L]
  reliability_basis <- if (identical(
    variance_component_method,
    "krr_gcv_variance_ratio"
  )) {
    "krr_gcv_variance_ratio_genetic_variance"
  } else {
    "reml_genetic_variance"
  }

  out <- data.frame(
    GID = ids,
    Predicted_value = suppressWarnings(as.numeric(pred[[pred_col]])),
    Train_Test_Label = ifelse(test_mask, "Test", "Train"),
    Standard_error = se_out,
    PEV = pev_out,
    Genetic_variance = genetic_variance_out,
    Reliability_reference_variance = reference_variance_out,
    Reliability = reliability,
    Reliability_remarks = reliability_remarks,
    Reliability_percentage = rep(reliability_percentage, length(ids)),
    Reliability_basis = rep(
      if (isTRUE(has_genetic_variance)) {
        reliability_basis
      } else {
        paste0(reliability_basis, "_not_estimable")
      },
      length(ids)
    ),
    Uncertainty_remarks = uncertainty_remarks_out,
    Observed_value = response_y,
    Train_Test = ifelse(test_mask, "Test", "Train"),
    stringsAsFactors = FALSE
  )
  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    # Honor user-supplied heter_groups as the output column name.
    env_out_name <- as.character(heter_groups)
    out[[env_out_name]] <- as.character(env_values)
    out <- out[, c("GID", env_out_name, setdiff(names(out), c("GID", env_out_name))), drop = FALSE]
  }
  out
}

gp_public_tuning_limit <- function(value, env_name, default) {
  if (exists("gp_policy_tuning_limit", mode = "function")) {
    return(gp_policy_tuning_limit(value, env_name, default))
  }
  out <- as.integer(value %||% Sys.getenv(env_name, unset = as.character(default)))[1L]
  if (!is.finite(out) || is.na(out) || out < 1L) {
    out <- as.integer(default)[1L]
  }
  out
}

gp_public_multitrait_validation_split <- function(long_df,
                                                  split,
                                                  gid_col = gp_backend_schema_col("long", "gid"),
                                                  y_col = gp_backend_schema_col("long", "y"),
                                                  validation_fraction = 0.25,
                                                  max_calibration_rows = NULL,
                                                  max_validation_rows = NULL,
                                                  seed = 12345L) {
  if (!isTRUE(length(split$train_mask) == nrow(long_df))) {
    return(list(status = "skipped", reason = "training mask length does not match phenotype rows"))
  }
  finite_train <- as.logical(split$train_mask) & is.finite(as.numeric(long_df[[y_col]]))
  if (!any(finite_train)) {
    return(list(status = "skipped", reason = "no finite training rows"))
  }
  train_gids <- unique(as.character(long_df[[gid_col]][finite_train]))
  if (length(train_gids) < 8L) {
    return(list(status = "skipped", reason = "too few training genotypes for MT validation"))
  }
  validation_fraction <- suppressWarnings(as.numeric(validation_fraction)[1L])
  if (!is.finite(validation_fraction) || validation_fraction <= 0 || validation_fraction >= 1) {
    validation_fraction <- 0.25
  }
  max_calibration_rows <- gp_public_tuning_limit(
    max_calibration_rows,
    "PREDICTPRO_MT_GP_TUNING_MAX_CALIBRATION_ROWS",
    4000L
  )
  max_validation_rows <- gp_public_tuning_limit(
    max_validation_rows,
    "PREDICTPRO_MT_GP_TUNING_MAX_VALIDATION_ROWS",
    1200L
  )
  gp_set_seed(as.integer(seed)[1L])
  n_val_gids <- max(3L, ceiling(length(train_gids) * validation_fraction))
  n_val_gids <- min(n_val_gids, length(train_gids) - 3L)
  val_gids <- sample(train_gids, n_val_gids)
  val_mask <- finite_train & as.character(long_df[[gid_col]]) %in% val_gids
  while (sum(val_mask) > max_validation_rows && length(val_gids) > 3L) {
    n_val_gids <- max(3L, floor(length(val_gids) * 0.8))
    val_gids <- sample(val_gids, n_val_gids)
    val_mask <- finite_train & as.character(long_df[[gid_col]]) %in% val_gids
  }
  cal_mask <- finite_train & !val_mask
  cal_gids <- unique(as.character(long_df[[gid_col]][cal_mask]))
  while (sum(cal_mask) > max_calibration_rows && length(cal_gids) > 3L) {
    keep_gids <- sample(cal_gids, max(3L, floor(length(cal_gids) * 0.8)))
    cal_mask <- finite_train &
      !(as.character(long_df[[gid_col]]) %in% val_gids) &
      (as.character(long_df[[gid_col]]) %in% keep_gids)
    cal_gids <- unique(as.character(long_df[[gid_col]][cal_mask]))
  }
  if (sum(cal_mask) < 10L || sum(val_mask) < 4L) {
    return(list(
      status = "skipped",
      reason = paste0("too few calibration/validation rows (cal=", sum(cal_mask), ", val=", sum(val_mask), ")")
    ))
  }
  list(
    status = "ok",
    reason = NA_character_,
    calibration_mask = cal_mask,
    validation_mask = val_mask,
    calibration_idx = as.integer(which(cal_mask) - 1L),
    validation_idx = as.integer(which(val_mask) - 1L),
    n_calibration_rows = as.integer(sum(cal_mask)),
    n_validation_rows = as.integer(sum(val_mask)),
    n_calibration_genotypes = as.integer(length(unique(as.character(long_df[[gid_col]][cal_mask])))),
    n_validation_genotypes = as.integer(length(unique(as.character(long_df[[gid_col]][val_mask]))))
  )
}

gp_public_multitrait_key_cols <- function(has_env = FALSE) {
  cols <- gp_backend_schema_cols()$long
  if (isTRUE(has_env)) {
    unname(c(cols$gid, cols$env, cols$trait))
  } else {
    unname(c(cols$gid, cols$trait))
  }
}

gp_public_multitrait_truth <- function(long_df, rows, key_cols) {
  y_col <- gp_backend_schema_col("long", "y")
  rows <- as.integer(rows)
  rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(long_df)]
  truth <- long_df[rows, c(key_cols, y_col), drop = FALSE]
  names(truth)[names(truth) == y_col] <- "Observed"
  truth$Observed <- as.numeric(truth$Observed)
  truth
}

gp_public_multitrait_baseline <- function(train_df, query_df, key_cols) {
  train_df <- train_df[is.finite(as.numeric(train_df$y)), , drop = FALSE]
  query <- query_df[, key_cols, drop = FALSE]
  query$.row_id <- seq_len(nrow(query))
  global_mean <- mean(as.numeric(train_df$y), na.rm = TRUE)
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  if ("env" %in% key_cols) {
    env_trait <- stats::aggregate(y ~ env + trait, train_df, mean, na.rm = TRUE)
    names(env_trait)[names(env_trait) == "y"] <- "Prediction_Baseline"
    out <- merge(query, env_trait, by = c("env", "trait"), all.x = TRUE, sort = FALSE)
  } else {
    out <- query
    out$Prediction_Baseline <- NA_real_
  }
  trait_mean <- stats::aggregate(y ~ trait, train_df, mean, na.rm = TRUE)
  names(trait_mean)[names(trait_mean) == "y"] <- "Prediction_TraitMean"
  out <- merge(out, trait_mean, by = "trait", all.x = TRUE, sort = FALSE)
  miss <- !is.finite(as.numeric(out$Prediction_Baseline))
  out$Prediction_Baseline[miss] <- out$Prediction_TraitMean[miss]
  miss <- !is.finite(as.numeric(out$Prediction_Baseline))
  out$Prediction_Baseline[miss] <- global_mean
  out <- out[order(out$.row_id), , drop = FALSE]
  out[, c(key_cols, "Prediction_Baseline"), drop = FALSE]
}

gp_public_multitrait_join_prediction_truth <- function(predictions, truth, key_cols) {
  pred <- gp_plain_data_frame(predictions)
  if (!all(key_cols %in% names(pred))) {
    stop("Multi-trait GP predictions are missing key columns: ", paste(setdiff(key_cols, names(pred)), collapse = ", "), call. = FALSE)
  }
  pred$.pred_row_id <- seq_len(nrow(pred))
  joined <- merge(pred, truth, by = key_cols, sort = FALSE)
  joined <- joined[order(joined$.pred_row_id), , drop = FALSE]
  joined$.pred_row_id <- NULL
  joined
}

gp_public_multitrait_score <- function(joined, baseline = NULL, key_cols = NULL) {
  observed <- as.numeric(joined$Observed)
  predicted <- as.numeric(joined$Prediction)
  keep <- is.finite(observed) & is.finite(predicted)
  if (!any(keep)) {
    return(list(rmse = NA_real_, mae = NA_real_, cor = NA_real_))
  }
  out <- list(
    rmse = sqrt(mean((predicted[keep] - observed[keep])^2)),
    mae = mean(abs(predicted[keep] - observed[keep])),
    cor = gp_tuning_safe_cor(predicted[keep], observed[keep])
  )
  if (!is.null(baseline) && "Prediction_Baseline" %in% names(baseline)) {
    base_joined <- merge(joined[, c(key_cols, "Observed"), drop = FALSE], baseline, by = key_cols, sort = FALSE)
    b <- as.numeric(base_joined$Prediction_Baseline)
    o <- as.numeric(base_joined$Observed)
    bkeep <- is.finite(b) & is.finite(o)
    out$baseline_rmse <- if (any(bkeep)) sqrt(mean((b[bkeep] - o[bkeep])^2)) else NA_real_
    out$baseline_cor <- if (any(bkeep)) gp_tuning_safe_cor(b[bkeep], o[bkeep]) else NA_real_
  }
  out
}

gp_public_multitrait_scale_uncertainty_frame <- function(x, alpha) {
  var_cols <- intersect(
    c("Prediction_Var_observed", "Prediction_Var_latent", "Prediction_variance", "PEV"),
    names(x)
  )
  se_cols <- intersect(
    c("Prediction_SE_observed", "Prediction_SE_latent", "Prediction_SE", "SE", "SE_observed", "SE_latent"),
    names(x)
  )
  for (nm in var_cols) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]])) * alpha^2
  }
  for (nm in se_cols) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]])) * alpha
  }
  x
}

gp_public_multitrait_apply_blend <- function(out,
                                             baseline,
                                             key_cols,
                                             alpha,
                                             label = "blocked_genotype_validation") {
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  if (!is.finite(alpha) || alpha < 0 || alpha > 1 ||
      is.null(out$predictions) || !nrow(gp_plain_data_frame(out$predictions))) {
    return(list(result = out, info = list(applied = FALSE, reason = "invalid blend alpha or missing predictions")))
  }
  pred <- gp_plain_data_frame(out$predictions)
  pred <- gp_public_normalize_multitrait_uncertainty(list(predictions = pred))$predictions
  pred$.pred_row_id <- seq_len(nrow(pred))
  merged <- merge(pred, baseline, by = key_cols, all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$.pred_row_id), , drop = FALSE]
  apply_rows <- is.finite(as.numeric(merged$Prediction)) & is.finite(as.numeric(merged$Prediction_Baseline))
  if (!any(apply_rows)) {
    return(list(result = out, info = list(applied = FALSE, reason = "no predictions matched MT baseline")))
  }
  gp_pred <- as.numeric(merged$Prediction)
  blended <- gp_pred
  blended[apply_rows] <- alpha * gp_pred[apply_rows] + (1 - alpha) * as.numeric(merged$Prediction_Baseline[apply_rows])
  merged$Prediction_GP <- gp_pred
  merged$Prediction <- blended
  merged$Prediction_Blend_Alpha <- alpha
  merged$Prediction_Blend_Method <- label
  merged <- gp_public_multitrait_scale_uncertainty_frame(merged, alpha)
  merged$.pred_row_id <- NULL
  out$predictions <- merged
  if (!is.null(out$result)) {
    out$result$predictions <- merged
  }
  list(
    result = out,
    info = list(
      applied = TRUE,
      alpha = alpha,
      method = label,
      se_variance_note = "Prediction SE/variance scaled by alpha and treats the validation baseline as fixed."
    )
  )
}

gp_public_multitrait_uncertainty_calibration <- function(joined,
                                                         min_n = 20L,
                                                         min_scale = NULL,
                                                         max_scale = NULL,
                                                         prior_n = NULL,
                                                         prior_env = NULL) {
  pred_var_info <- gp_tuning_pick_numeric_column_info(
    joined,
    c("Prediction_Var_observed", "Prediction_variance", "Prediction_Var_latent", "PEV")
  )
  pred_var <- pred_var_info$values
  pred_se_info <- gp_tuning_pick_numeric_column_info(
    joined,
    c("Prediction_SE_observed", "SE", "SE_observed", "Prediction_SE", "Prediction_SE_latent", "SE_latent")
  )
  pred_se <- pred_se_info$values
  variance_from_se <- FALSE
  if (all(!is.finite(pred_var)) && any(is.finite(pred_se))) {
    pred_var <- pred_se^2
    pred_var_info$column <- pred_se_info$column
    variance_from_se <- TRUE
  }
  residual <- as.numeric(joined$Prediction) - as.numeric(joined$Observed)
  keep <- is.finite(residual) & is.finite(pred_var) & pred_var > 0
  min_n <- as.integer(min_n %||% 20L)[1L]
  if (!is.finite(min_n) || is.na(min_n) || min_n < 2L) {
    min_n <- 20L
  }
  if (sum(keep) < min_n) {
    return(list(
      status = "skipped",
      reason = paste0("too few validation rows with finite MT prediction variance (n=", sum(keep), ")"),
      n_validation_se = as.integer(sum(keep))
    ))
  }
  residual <- residual[keep]
  pred_var <- pred_var[keep]
  empirical_mse <- mean(residual^2)
  model_mse <- mean(pred_var)
  if (!is.finite(empirical_mse) || !is.finite(model_mse) || model_mse <= 0) {
    return(list(status = "skipped", reason = "non-finite MT validation uncertainty moments"))
  }
  n_env <- if ("env" %in% names(joined)) length(unique(as.character(joined$env[keep]))) else 1L
  regularized <- gp_tuning_regularize_uncertainty_scale(
    sqrt(empirical_mse / model_mse),
    n_validation = length(residual),
    n_env = n_env,
    min_scale = min_scale,
    max_scale = max_scale,
    prior_n = prior_n,
    prior_env = prior_env
  )
  if (!is.finite(regularized$scale)) {
    return(list(status = "skipped", reason = "MT validation uncertainty scale was not finite"))
  }
  list(
    status = "ok",
    reason = NA_character_,
    scale = as.numeric(regularized$scale),
    variance_scale = as.numeric(regularized$scale)^2,
    raw_scale = sqrt(empirical_mse / model_mse),
    global_scale = sqrt(empirical_mse / model_mse),
    shrinkage_weight = as.numeric(regularized$shrinkage_weight),
    min_scale = as.numeric(regularized$min_scale),
    max_scale = as.numeric(regularized$max_scale),
    prior_n = as.numeric(regularized$prior_n),
    prior_env = as.numeric(regularized$prior_env),
    n_validation_se = as.integer(length(residual)),
    n_validation_env_se = as.integer(n_env),
    validation_rmse_for_se = sqrt(empirical_mse),
    validation_mean_model_se = sqrt(model_mse),
    variance_column = as.character(pred_var_info$column %||% NA_character_)[1L],
    se_column = as.character(pred_se_info$column %||% NA_character_)[1L],
    uncertainty_target = gp_tuning_uncertainty_column_target(pred_var_info$column),
    variance_from_se = isTRUE(variance_from_se),
    method = "multi_trait_blocked_genotype_residual_mse_over_mean_prediction_variance"
  )
}

gp_public_multitrait_validation_adjustment <- function(long_df,
                                                       split,
                                                       key_cols,
                                                       fit_validation,
                                                       prediction_blend = "auto",
                                                       gp_uncertainty_calibration = "auto",
                                                       return_se = TRUE,
                                                       validation_fraction = 0.25,
                                                       max_calibration_rows = NULL,
                                                       max_validation_rows = NULL,
                                                       seed = 12345L) {
  long_cols <- gp_backend_schema_cols()$long
  validation_split <- gp_public_multitrait_validation_split(
    long_df = long_df,
    split = split,
    gid_col = long_cols$gid,
    y_col = long_cols$y,
    validation_fraction = validation_fraction,
    max_calibration_rows = max_calibration_rows,
    max_validation_rows = max_validation_rows,
    seed = seed
  )
  info <- list(
    requested = TRUE,
    status = validation_split$status,
    reason = validation_split$reason %||% NA_character_,
    validation_split = validation_split
  )
  if (!identical(validation_split$status, "ok")) {
    return(list(
      info = info,
      blend = list(applied = FALSE, reason = validation_split$reason %||% "validation split unavailable"),
      calibration = list(status = "skipped", reason = validation_split$reason %||% "validation split unavailable")
    ))
  }
  validation_fit <- tryCatch(
    fit_validation(validation_split),
    error = function(e) e
  )
  if (inherits(validation_fit, "error")) {
    info$status <- "failed"
    info$reason <- conditionMessage(validation_fit)
    return(list(
      info = info,
      blend = list(applied = FALSE, reason = info$reason),
      calibration = list(status = "skipped", reason = info$reason)
    ))
  }
  validation_fit <- gp_public_normalize_multitrait_uncertainty(validation_fit)
  truth <- gp_public_multitrait_truth(
    long_df,
    rows = validation_split$validation_idx + 1L,
    key_cols = key_cols
  )
  joined <- gp_public_multitrait_join_prediction_truth(validation_fit$predictions, truth, key_cols)
  baseline <- gp_public_multitrait_baseline(
    train_df = long_df[validation_split$calibration_mask, , drop = FALSE],
    query_df = long_df[validation_split$validation_idx + 1L, , drop = FALSE],
    key_cols = key_cols
  )
  metrics <- gp_public_multitrait_score(joined, baseline = baseline, key_cols = key_cols)
  joined_for_blend <- joined
  joined_for_blend$.blend_row_id <- seq_len(nrow(joined_for_blend))
  joined_for_blend <- merge(joined_for_blend, baseline, by = key_cols, sort = FALSE)
  joined_for_blend <- joined_for_blend[order(joined_for_blend$.blend_row_id), , drop = FALSE]
  blend <- gp_tuning_select_blend(
    gp_prediction = joined_for_blend$Prediction,
    baseline_prediction = joined_for_blend$Prediction_Baseline,
    observed = joined_for_blend$Observed,
    env = if ("env" %in% key_cols) joined_for_blend$env else joined_for_blend$trait,
    objective = "rmse"
  )
  relative_improvement <- if (is.finite(blend$blend_rmse) && is.finite(blend$baseline_rmse) &&
                              is.finite(metrics$rmse) && metrics$rmse > 0) {
    (metrics$rmse - blend$blend_rmse) / metrics$rmse
  } else {
    NA_real_
  }
  blend_decision <- gp_prediction_blend_auto_decision(
    prediction_blend = prediction_blend,
    blend_alpha = blend$blend_alpha,
    relative_improvement = relative_improvement,
    allow_auto = TRUE
  )
  blend_info <- c(
    list(requested = prediction_blend),
    blend,
    list(
      applied = isTRUE(blend_decision$apply),
      reason = blend_decision$reason,
      relative_improvement = relative_improvement
    )
  )
  calibration_joined <- joined
  if (isTRUE(blend_decision$apply)) {
    baseline_for_joined <- baseline[, c(key_cols, "Prediction_Baseline"), drop = FALSE]
    calibration_joined <- merge(calibration_joined, baseline_for_joined, by = key_cols, sort = FALSE)
    alpha <- as.numeric(blend$blend_alpha)[1L]
    calibration_joined$Prediction <- alpha * as.numeric(calibration_joined$Prediction) +
      (1 - alpha) * as.numeric(calibration_joined$Prediction_Baseline)
    calibration_joined <- gp_public_multitrait_scale_uncertainty_frame(calibration_joined, alpha)
  }
  calibration <- if (isTRUE(return_se) && !identical(gp_uncertainty_calibration, "none")) {
    gp_public_multitrait_uncertainty_calibration(calibration_joined)
  } else {
    list(
      status = "skipped",
      reason = if (identical(gp_uncertainty_calibration, "none")) "uncertainty_calibration_none" else "prediction_se_not_requested"
    )
  }
  info$status <- "ok"
  info$reason <- NA_character_
  info$metrics <- metrics
  info$validation_predictions <- nrow(joined)
  list(
    info = info,
    blend = blend_info,
    calibration = calibration,
    baseline_train_mask = validation_split$calibration_mask
  )
}

gp_public_result_to_r <- function(py_out) {
  gp_require_reticulate("GP reticulate result conversion")
  converted <- tryCatch(reticulate::py_to_r(py_out), error = function(e) NULL)
  py_item <- function(name) {
    if (is.list(converted) && !is.null(names(converted)) && name %in% names(converted)) {
      return(converted[[name]])
    }
    out <- tryCatch({
      if (reticulate::py_has_attr(py_out, name)) {
        reticulate::py_to_r(reticulate::py_get_attr(py_out, name))
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
    tryCatch(reticulate::py_to_r(py_out$get(name)), error = function(e) NULL)
  }
  object_item <- function(x, name) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.list(x) && !inherits(x, "python.builtin.object") && !is.null(names(x)) && name %in% names(x)) {
      return(x[[name]])
    }
    out <- tryCatch({
      if (reticulate::py_has_attr(x, name)) {
        reticulate::py_to_r(reticulate::py_get_attr(x, name))
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
    tryCatch(reticulate::py_to_r(x$get(name)), error = function(e) NULL)
  }
  is_plain_r_list <- function(x) {
    is.list(x) && !inherits(x, "python.builtin.object") && !inherits(x, "python.builtin.module")
  }
  result <- py_item("result")
  predictions <- py_item("predictions")
  result_predictions <- object_item(result, "predictions") %||% predictions
  if (!is_plain_r_list(result)) {
    result <- list(predictions = result_predictions)
  }
  result_extra_names <- c(
    "varcomp",
    "summary",
    "ai_matrix",
    "var_components",
    "var_components_summary",
    "var_components_ai",
    "env_variance_summary",
    "residual_summary",
    "heritability_summary",
    "interaction_variance_summary",
    "interaction_correlation_summary",
    "env_covariance_fa",
    "env_correlation_fa",
    "env_covariance_fa_structure",
    "fa_kernel_weights",
    "sigma2_resid_env",
    "sigma2_resid_overall",
    "variance_component_method",
    "diagnostics"
  )
  for (nm in result_extra_names) {
    val <- py_item(nm)
    if (!is.null(val) && is.null(result[[nm]])) {
      result[[nm]] <- val
    }
  }
  if (is.null(result[["per_env"]]) && !is.null(result[["env_variance_summary"]])) {
    result[["per_env"]] <- result[["env_variance_summary"]]
  }
  fit <- py_item("fit")
  if (!is_plain_r_list(fit)) {
    fit <- list()
  }
  info <- py_item("info")
  if (is.null(info)) {
    info <- py_item("_meta")
  }
  if (!is_plain_r_list(info)) {
    info <- list()
  }
  out <- list(result = result, fit = fit, info = info)
  if (!is.null(result_predictions)) {
    out$predictions <- gp_plain_data_frame(result_predictions)
  } else if (!is.null(predictions)) {
    out$predictions <- gp_plain_data_frame(predictions)
  }
  class(out) <- c("predictpror_gp_result", class(out))
  out
}

#' Fit a Single-Trait Gaussian-Process Model
#'
#' R wrapper around the vendored GP backend for single-trait, single- or
#' multi-environment Gaussian-process prediction.
#'
#' @param pheno_data Data frame containing genotype IDs, the response, and
#' optionally an environment column.
#' @param gmatrix Optional genomic relationship matrix. Required for dense
#' REML fits. For scalable \code{varcomp_mode = "mom"} fits, it may be omitted
#' when \code{gp_factor_cache} and matching \code{geno_ids} are supplied.
#' @param response Response column name.
#' @param gen_name Genotype ID column name.
#' @param heter_groups Optional environment column name.
#' @param test_set Optional genotype IDs to predict. If omitted, rows with
#' missing response values are treated as prediction rows.
#' @param test_set_source Optional internal provenance for \code{test_set}.
#' When set to \code{"inferred_missing_response"}, the backend uses missing
#' response rows as cell-level prediction rows rather than holding out every
#' record for those genotype IDs.
#' @param env_similarity Optional environment relationship matrix for
#' multi-environment reaction-norm prediction. For large MET with more than
#' five environments, supply an environment kernel from covariates or a
#' defensible similarity matrix instead of relying on an identity environment
#' kernel.
#' @param env_covariates Optional environment covariates for backend
#' environment-kernel construction. When supplied, the backend QC-aligns them
#' to the phenotype environments and converts them to \code{env_similarity}.
#' Pass exactly one of \code{env_similarity} or \code{env_covariates}.
#' @param reaction_norm_feature_qc Logical; when \code{TRUE}, use backend
#' feature QC before converting environment covariates to an environment kernel.
#' @param kenv_kernel Environment-covariate kernel family used when
#' \code{env_covariates} is supplied.
#' @param kenv_bandwidth Numeric bandwidth used by the environment kernel.
#' @param kenv_kernel_kwargs Optional named list of backend kernel arguments.
#' @param tune_gp Logical or \code{"auto"}; when \code{TRUE}, run historical
#' tuning for MET GP fits. The tuning split uses
#' only finite training rows and masks validation responses before candidate
#' fits. With \code{"auto"}, PredictPro enables tuning only for
#' multi-environment fits with enough training rows and environments for a
#' defensible internal validation. Environment covariates are required for
#' prospective/new-environment prediction, but not for known-environment
#' genotype holdout.
#' @param tuning_strategy Tuning split strategy. \code{"auto"} uses
#' blocked-environment validation for prospective/new-environment prediction
#' and blocked-genotype validation for known-environment genotype prediction.
#' @param tuning_krr_lams Positive ridge/lambda values to evaluate during
#' blocked tuning.
#' @param tuning_kenv_kernels Environment-covariate kernel families to
#' evaluate during blocked tuning.
#' @param tuning_weight_grid Optional data frame with columns \code{w_g},
#' \code{w_ge}, and \code{w_e} for component-weight tuning.
#' @param tuning_objective Validation objective used to select the final
#' tuning candidate. \code{"auto"} chooses a conservative objective from the
#' prediction task. \code{"rmse"} preserves the historical raw-RMSE
#' behavior. \code{"centered_rmse"} removes validation environment means before
#' scoring. \code{"within_env_spearman"} optimizes weighted within-environment
#' rank correlation. \code{"centered_cor"} maximizes centered Pearson
#' correlation.
#' @param mean_adjustment One of \code{"auto"}, \code{"none"}, or
#' \code{"location_offset"}. With tuning, \code{"auto"} evaluates both no
#' adjustment and a training-location mean offset. Without tuning,
#' \code{"auto"} leaves the response unchanged.
#' @param prediction_blend One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful MET tuning, \code{"auto"} only applies
#' conservative baseline-dominant environment/location mean blends for
#' prospective/new-environment prediction when they clearly improve the
#' internal validation objective. Known-environment genotype prediction keeps
#' the GP signal by default. \code{"validation"} follows the validation-selected
#' blend directly.
#' @param env_location Optional environment-location annotation. May be a
#' phenotype column name, named vector keyed by environment, or data frame with
#' environment IDs and locations. If omitted, locations are inferred from
#' environment names by removing a year token.
#' @param env_year Optional environment-year annotation. May be a phenotype
#' column name, named vector keyed by environment, or data frame with
#' environment IDs and years. If omitted, years are parsed from environment
#' names when possible.
#' @param tuning_validation_fraction Fraction of training genotypes/environments
#' targeted for blocked validation.
#' @param tuning_max_calibration_rows Optional maximum calibration rows for
#' internal GP tuning. When \code{NULL}, dense tuning keeps the legacy full
#' blocked split and factor-cache tuning uses a bounded production-safe subset.
#' @param tuning_max_validation_rows Optional maximum validation rows for
#' internal GP tuning.
#' @param tuning_max_genotypes Optional maximum unique genotypes used by the
#' dense tuning kernel materialized from a factor cache.
#' @param geno_ids Optional genotype IDs matching the rows and columns of
#' \code{gmatrix}.
#' @param kernel_list Optional named list of additional relationship/kernel
#' matrices. When supplied, GP fits pass all matrices as backend
#' \code{geno_kernels}.
#' @param kernel_weights Optional non-negative weights for the aligned genomic
#' kernels. Named weights are reordered to the aligned kernel names. The
#' default assigns weight one to every kernel.
#' @param fixed_effects Optional fixed-effect columns passed to the backend.
#' When \code{NULL}, multi-environment fits use the environment column as a
#' fixed effect if all environments are represented in the training rows.
#' @param observation_weights Optional Stage 2 precision weights. Supply a
#' positive numeric vector, a one-column table, or a column name from
#' \code{pheno_data}. The GP residual diagonal is fixed to
#' \code{1 / observation_weights}.
#' @param gp_return_se Logical; request prediction standard errors and
#' prediction error variances.
#' @param gp_uncertainty_calibration One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful GP tuning and requested SE/PEV output,
#' \code{"auto"} applies a shrinkage-regularized validation residual scale to
#' prediction SE and variance columns. \code{"validation"} forces the same
#' validation-derived calibration when available; \code{"none"} leaves model
#' uncertainty unscaled.
#' @param gp_full_vc Logical; request full variance-component output where
#' supported.
#' @param gp_se_hutchinson_probes Integer number of Hutchinson probes used for
#' scalable operator/tiered prediction SEs. Dense paths ignore this argument.
#' @param gp_varcomp_mode Variance-component estimation mode passed to the GP
#' backend. Gaussian-process fits use the requested mode; dense
#' \code{model_name = "KRR"} fits instead report the variance ratio selected by
#' \code{lam_select}, because a separate REML post-fit would not describe the
#' KRR predictor. Use \code{"mom"} only for the scalable factor-cache path.
#' @param gp_factor_cache Optional low-rank genomic factor cache for large
#' operator-backend fits. When supplied with \code{geno_ids}, \code{gmatrix}
#' may be \code{NULL}; the backend uses the factor cache instead of a dense
#' genomic relationship matrix.
#' @param gp_tiered_dispatch Logical; route MET predictions by environment
#' connectivity. This is intended for large prospective MET exercises where
#' all-new test environments should be served by the covariate-derived
#' environment kernel without fitting an unstructured environment covariance.
#' @param gp_auto_policy Logical; when \code{TRUE}, let PredictPro choose
#' conservative backend, dtype, device route, tiered dispatch, and optional
#' \code{tune_gp = "auto"} behavior from the data shape. Set to \code{FALSE}
#' for fully manual settings.
#' @param model_name GP backend model label. One of \code{"KRR"}, \code{"GP"},
#' or \code{"GP_FA"}. \code{"LowRankGP"} is not a distinct single-trait
#' estimator and is rejected; use \code{"KRR"}. It remains available through
#' the purpose-specific joint multi-trait route.
#' @param krr_lam Numeric KRR/GP ridge parameter used when
#' \code{lam_select = "fixed"}.
#' @param krr_lams Numeric candidate residual-to-genetic variance ratios used by
#' the KRR backend when \code{lam_select} requests automatic selection. The
#' default spans 0.001 through 10.
#' @param lam_select KRR lambda selection mode. The default \code{"gcv"} uses
#' deterministic exact generalized cross-validation for dense fits up to 2,048
#' training rows; \code{"fixed"} uses \code{krr_lam}.
#' @param gp_iters Optional positive integer optimizer iterations for the
#' exact GP backend. \code{NULL} preserves the Python backend default.
#' @param gp_lr Optional positive learning rate for the exact GP backend.
#' \code{NULL} preserves the Python backend default.
#' @param gp_learn_scales Optional logical override for backend scale learning.
#' When \code{NULL}, scales are learned only for full variance-component output.
#' @param gp_theta0_warm Optional numeric starting values for the GP
#' variance-component optimizer. Intended for compatible warm-started fits.
#' @param gp_engine Character; which REML engine drives GP variance-component
#' estimation. \code{"auto"} (default) picks the fastest available
#' (currently resolves to \code{"dense_v"}); \code{"dense_v"} keeps the
#' dense V-formulation (Phase 1-2); \code{"mme"} uses the sparse Mixed
#' Model Equations engine (Phase 3.4-3.7, fastest at n_geno = 100 for MET
#' but regresses at n_geno >= 300 on GBLUP because G^-1 is dense -- see
#' Phase 3.14 NEWS); \code{"eigen"} (Phase 3.15a, experimental) uses an
#' eigen-projected REML formulation that would beat ASReml at large n on
#' the standalone validation, but the production integration is incomplete
#' (Phase 3.15b NEWS) so \code{"eigen"} currently falls back to dense_v.
#' Requires \code{PREDICTPRO_GP_DISABLE_ENV_MAIN=1} for \code{"mme"} and
#' \code{"eigen"} on MET. Honors \code{PREDICTPRO_GP_ENGINE} env var.
#' @return A list with \code{predictions}, \code{result}, \code{fit}, and
#' \code{info}.
#' @export
gp_single_trait_model <- function(pheno_data,
                                  gmatrix = NULL,
                                  response,
                                  gen_name,
                                  heter_groups = NULL,
                                  test_set = NULL,
                                  test_set_source = NULL,
                                  env_similarity = NULL,
                                  env_ids = NULL,
                                  env_covariates = NULL,
                                  reaction_norm_feature_qc = TRUE,
                                  kenv_kernel = "matern32",
                                  kenv_bandwidth = 1.0,
                                  kenv_kernel_kwargs = NULL,
                                  tune_gp = FALSE,
                                  tuning_strategy = "auto",
                                  tuning_krr_lams = c(0.03, 0.1, 0.3),
                                  tuning_kenv_kernels = c("matern32", "matern52", "rbf"),
                                  tuning_weight_grid = NULL,
                                  tuning_objective = "rmse",
                                  mean_adjustment = "auto",
                                  prediction_blend = "auto",
                                  env_location = NULL,
                                  env_year = NULL,
                                  tuning_validation_fraction = 0.25,
                                  tuning_max_calibration_rows = NULL,
                                  tuning_max_validation_rows = NULL,
                                  tuning_max_genotypes = NULL,
                                  geno_ids = NULL,
                                  kernel_list = NULL,
                                  include_components = NULL,
                                  fixed_effects = NULL,
                                  observation_weights = NULL,
                                  model_name = "KRR",
                                  gp_backend = "auto",
                                  gp_output_level = "predict_only",
                                  gp_return_se = FALSE,
                                  gp_uncertainty_calibration = "auto",
                                  gp_full_vc = FALSE,
                                  gp_factor_cache = NULL,
                                  gp_tiered_dispatch = FALSE,
                                  gp_tiered_verbose = FALSE,
                                  gp_auto_policy = TRUE,
                                  gp_dtype = "float64",
                                  large_n_threshold = 50000L,
                                  operator_tol = 1e-5,
                                  operator_max_iter = 500L,
                                  operator_dtype_compute = "float32",
                                  gp_se_hutchinson_probes = 32L,
                                  gp_varcomp_mode = "reml",
                                  gp_fa_rank = 1L,
                                  gp_prediction_output = "test_only",
                                  krr_lam = 1e-3,
                                  krr_lams = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10),
                                  lam_select = "gcv",
                                  gp_iters = NULL,
                                  gp_lr = NULL,
                                  w_g = NULL,
                                  w_ge = NULL,
                                  w_e = 0,
                                  gp_learn_scales = NULL,
                                  seed = 12345L,
                                  python_bin = NULL,
                                  project_root = NULL,
                                  gp_theta0_warm = NULL,
                                  gp_engine = "auto") {
  gp_validate_single_trait_gp_model_scope(GS_model = model_name)
  if (!is.null(env_similarity) && !is.null(env_covariates)) {
    stop("Pass exactly one of `env_similarity` or `env_covariates`.", call. = FALSE)
  }
  tuning_strategy <- gp_match_gp_choice(
    tuning_strategy,
    choices = c("auto", "blocked_environment", "blocked_genotype"),
    name = "tuning_strategy",
    default = "auto"
  )
  mean_adjustment <- gp_match_gp_choice(
    mean_adjustment,
    choices = c("auto", "none", "location_offset"),
    name = "mean_adjustment",
    default = "auto"
  )
  prediction_blend <- gp_match_gp_choice(
    prediction_blend,
    choices = c("auto", "none", "validation"),
    name = "prediction_blend",
    default = "auto"
  )
  gp_uncertainty_calibration <- gp_match_gp_choice(
    gp_uncertainty_calibration,
    choices = c("auto", "none", "validation"),
    name = "gp_uncertainty_calibration",
    default = "auto"
  )
  if (!gp_policy_is_auto(tuning_objective)) {
    tuning_objective <- gp_tuning_normalize_objective(tuning_objective)
  }
  missing_cols <- setdiff(c(gen_name, response, heter_groups), names(pheno_data))
  if (length(missing_cols)) {
    stop("pheno_data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  pheno_ids <- as.character(pheno_data[[gen_name]])
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    kernel_list = kernel_list
  )
  kernel_source <- if (length(kernel_inputs) > 0L) kernel_inputs else gmatrix
  kernel <- gp_public_kernel_bank_inputs(
    kernel_source,
    pheno_ids = pheno_ids,
    geno_ids = geno_ids,
    allow_factor_cache = !is.null(gp_factor_cache)
  )
  backend_schema <- gp_backend_schema_cols()
  single_cols <- backend_schema$single
  backend_gid_col <- single_cols$gid
  backend_env_col <- single_cols$env
  backend_y_col <- single_cols$y
  pheno_df <- data.frame(
    .gid = pheno_ids,
    .env = if (!is.null(heter_groups)) as.character(pheno_data[[heter_groups]]) else "ENV1",
    .y = as.numeric(pheno_data[[response]]),
    stringsAsFactors = FALSE
  )
  names(pheno_df) <- unname(c(backend_gid_col, backend_env_col, backend_y_col))
  env_location <- gp_tuning_resolve_pheno_annotation(
    env_location,
    pheno_data = pheno_data,
    pheno_env = pheno_df[[backend_env_col]],
    target_name = "location"
  )
  env_year <- gp_tuning_resolve_pheno_annotation(
    env_year,
    pheno_data = pheno_data,
    pheno_env = pheno_df[[backend_env_col]],
    target_name = "year"
  )
  if (!is.null(heter_groups)) {
    gp_warn_large_met_without_env_kernel(
      pheno_df[[backend_env_col]],
      has_env_kernel = !is.null(env_similarity) || !is.null(env_covariates),
      context = "Single-trait MET GP"
    )
  }
  split <- gp_public_train_test_idx(
    pheno_df,
    gid_col = backend_gid_col,
    y_col = backend_y_col,
    test_set = test_set,
    test_set_source = test_set_source
  )
  stage2_weights <- gp_resolve_stage2_precision_weights(
    weights = observation_weights,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    test_set = test_set,
    test_mask = split$test_mask,
    context = "GP Stage 2 observation weights"
  )
  if (is.null(fixed_effects) && !is.null(heter_groups)) {
    train_envs <- unique(as.character(pheno_df[[backend_env_col]][split$train_mask]))
    all_envs <- unique(as.character(pheno_df[[backend_env_col]]))
    fixed_effects <- if (all(all_envs %in% train_envs)) backend_env_col else character(0)
  }
  output_level <- if (isTRUE(gp_return_se) || isTRUE(gp_full_vc)) {
    gp_bridge_output_level(gp_return_se, gp_full_vc)
  } else {
    gp_output_level
  }
  gp_se_hutchinson_probes <- as.integer(gp_se_hutchinson_probes %||% 32L)[1L]
  if (!is.finite(gp_se_hutchinson_probes) ||
      is.na(gp_se_hutchinson_probes) ||
      gp_se_hutchinson_probes < 1L) {
    gp_se_hutchinson_probes <- 32L
  }
  policy <- gp_single_trait_auto_policy(
    pheno_df = pheno_df,
    split = split,
    heter_groups = heter_groups,
    gid_col = backend_gid_col,
    env_col = backend_env_col,
    env_similarity = env_similarity,
    env_covariates = env_covariates,
    gmatrix = kernel$K,
    gp_factor_cache = gp_factor_cache,
    geno_ids = kernel$geno_ids,
    gp_backend = gp_backend,
    gp_tiered_dispatch = gp_tiered_dispatch,
    gp_dtype = gp_dtype,
    large_n_threshold = large_n_threshold,
    operator_tol = operator_tol,
    operator_max_iter = operator_max_iter,
    operator_dtype_compute = operator_dtype_compute,
    tune_gp = tune_gp,
    tuning_objective = tuning_objective,
    tuning_max_calibration_rows = tuning_max_calibration_rows,
    tuning_max_validation_rows = tuning_max_validation_rows,
    tuning_max_genotypes = tuning_max_genotypes,
    include_components = include_components,
    enabled = isTRUE(gp_auto_policy)
  )
  gp_backend <- policy$gp_backend
  gp_tiered_dispatch <- policy$gp_tiered_dispatch
  gp_dtype <- policy$gp_dtype
  operator_tol <- policy$operator_tol
  operator_max_iter <- policy$operator_max_iter
  operator_dtype_compute <- policy$operator_dtype_compute
  tune_gp_effective <- isTRUE(policy$tune_gp)
  tuning_objective <- gp_tuning_normalize_objective(policy$tuning_objective)
  tuning_strategy_effective <- if (identical(tuning_strategy, "auto")) {
    if (isTRUE(policy$profile$prospective_env_prediction)) {
      "blocked_environment"
    } else {
      "blocked_genotype"
    }
  } else {
    tuning_strategy
  }
  include_components <- policy$include_components

  tuning_info <- list(
    requested = tune_gp_effective,
    requested_value = tune_gp,
    auto_policy_reason = policy$auto_tune_reason %||% NA_character_,
    strategy = tuning_strategy_effective,
    requested_strategy = tuning_strategy,
    tuning_objective = tuning_objective,
    status = if (tune_gp_effective) {
      "not_run"
    } else if (gp_policy_is_auto(tune_gp)) {
      "skipped"
    } else {
      "not_requested"
    },
    reason = if (!tune_gp_effective && gp_policy_is_auto(tune_gp)) {
      policy$auto_tune_reason %||% NA_character_
    } else {
      NA_character_
    }
  )
  final_mean_mode <- if (identical(mean_adjustment, "location_offset")) {
    "location_offset"
  } else {
    "none"
  }
  if (tune_gp_effective) {
    single_environment_tuning <- !isTRUE(policy$profile$is_met)
    tuning_requires_env_covariates <- isTRUE(policy$profile$prospective_env_prediction)
    if ((!isTRUE(single_environment_tuning) && is.null(heter_groups)) ||
        !is.null(env_similarity) ||
        (isTRUE(tuning_requires_env_covariates) && is.null(env_covariates))) {
      tuning_info$status <- "skipped"
      tuning_info$reason <- if (isTRUE(tuning_requires_env_covariates) && is.null(env_covariates)) {
        paste(
          "GP tuning for prospective/new-environment prediction requires",
          "env_covariates and no env_similarity"
        )
      } else if (!isTRUE(single_environment_tuning) && is.null(heter_groups)) {
        paste(
          "GP tuning requires either a single-environment genotype holdout",
          "or a multi-environment single-trait GP fit"
        )
      } else {
        paste(
          "GP tuning currently requires no env_similarity; pass",
          "env_covariates for environment-kernel tuning"
        )
      }
      warning(tuning_info$reason, call. = FALSE)
    } else {
      tuning_env_covariates <- if (is.null(env_covariates)) {
        NULL
      } else {
        gp_prepare_env_covariates(
          env_covariates,
          source_env_col = heter_groups,
          target_env_col = backend_env_col
        )
      }
      tuning <- gp_tune_single_trait_met(
        pheno_df = pheno_df,
        split = split,
        kernel = kernel,
        gp_factor_cache = gp_factor_cache,
        geno_ids = kernel$geno_ids,
        env_covariates = tuning_env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        include_components = include_components,
        gp_backend = gp_backend,
        tuning_krr_lams = tuning_krr_lams,
        tuning_kenv_kernels = tuning_kenv_kernels,
        tuning_weight_grid = tuning_weight_grid,
        tuning_objective = tuning_objective,
        mean_adjustment = mean_adjustment,
        fixed_effects = fixed_effects %||% character(0),
        observation_weights = stage2_weights$precision,
        env_location = env_location,
        env_year = env_year,
        tuning_validation_fraction = tuning_validation_fraction,
        tuning_strategy = if (isTRUE(single_environment_tuning)) "blocked_genotype" else tuning_strategy_effective,
        tuning_max_calibration_rows = policy$tuning_limits$max_calibration_rows %||% tuning_max_calibration_rows,
        tuning_max_validation_rows = policy$tuning_limits$max_validation_rows %||% tuning_max_validation_rows,
        tuning_max_genotypes = policy$tuning_limits$max_genotypes %||% tuning_max_genotypes,
        calibrate_uncertainty = output_level %in% c("predict_with_se", "full_vc") &&
          !identical(gp_uncertainty_calibration, "none"),
        gp_se_hutchinson_probes = gp_se_hutchinson_probes,
        seed = seed,
        python_bin = python_bin,
        project_root = project_root
      )
      tuning_info <- c(
        list(requested = TRUE, strategy = tuning_strategy_effective, requested_strategy = tuning_strategy),
        tuning
      )
      if (identical(tuning$status, "ok")) {
        selected <- tuning$selected[1L, , drop = FALSE]
        krr_lam <- as.numeric(selected$lambda)[1L]
        krr_lams <- "none"
        lam_select <- "fixed"
        kenv_kernel <- as.character(selected$kenv_kernel)[1L]
        w_g <- as.numeric(selected$w_g)[1L]
        w_ge <- as.numeric(selected$w_ge)[1L]
        w_e <- as.numeric(selected$w_e)[1L]
        final_mean_mode <- as.character(selected$mean_mode)[1L]
      } else {
        warning("GP tuning skipped: ", tuning$reason, call. = FALSE)
      }
    }
  }

  mean_fit <- gp_tuning_apply_mean_adjustment(
    pheno_df,
    training_rows = which(split$train_mask),
    mode = final_mean_mode,
    env_location = env_location,
    env_year = env_year
  )
  pheno_fit_df <- mean_fit$pheno_df

  args <- list(
    pheno_df = pheno_fit_df,
    gid_col = backend_gid_col,
    env_col = backend_env_col,
    y_col = backend_y_col,
    geno_ids = kernel$geno_ids,
    include_components = as.list(include_components),
    fixed_effects = as.list(fixed_effects %||% character(0)),
    method = gp_backend_method_map(model_name),
    backend = gp_backend,
    prediction_output = gp_prediction_output,
    output_level = output_level,
    point_predictions_only = identical(output_level, "predict_only"),
    return_se = output_level %in% c("predict_with_se", "full_vc"),
    compute_ai_se = identical(output_level, "full_vc"),
    standardize = "global",
    varcomp_mode = gp_varcomp_mode,
    krr_lam = as.numeric(krr_lam)[1L],
    krr_lams = if (is.character(krr_lams)) krr_lams else as.numeric(krr_lams),
    lam_select = as.character(lam_select)[1L],
    dtype = as.character(gp_dtype)[1L],
    seed = as.integer(seed),
    train_idx = split$train_idx,
    test_idx = split$test_idx,
    large_n_threshold = as.integer(large_n_threshold)[1L],
    operator_tol = as.numeric(operator_tol)[1L],
    operator_max_iter = as.integer(operator_max_iter)[1L],
    operator_dtype_compute = as.character(operator_dtype_compute)[1L],
    n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    operator_n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    tiered_n_hutchinson_probes = as.integer(gp_se_hutchinson_probes),
    tiered_dispatch = isTRUE(gp_tiered_dispatch),
    tiered_verbose = isTRUE(gp_tiered_verbose),
    learn_scales = if (is.null(gp_learn_scales)) {
      identical(output_level, "full_vc")
    } else {
      isTRUE(gp_learn_scales)
    }
  )
  if (isTRUE(stage2_weights$supplied)) {
    args$obs_weights <- as.numeric(stage2_weights$precision)
    args$obs_weight_mode <- "inverse_variance"
    args$obs_weight_global_scale <- 1
  }
  args$gp_engine <- as.character(gp_engine)[1L]
  args <- gp_bridge_add_gp_exact_controls(args, list(gp_iters = gp_iters, gp_lr = gp_lr))
  if (!is.null(gp_factor_cache)) {
    args$grm_factor_cache <- gp_factor_cache
  }
  if (!is.null(env_similarity)) {
    env_similarity <- as.matrix(env_similarity)
    storage.mode(env_similarity) <- "double"
    args$env_similarity <- unname(env_similarity)
    if (is.null(env_ids) && !is.null(rownames(env_similarity))) {
      env_ids <- rownames(env_similarity)
    }
    if (!is.null(env_ids)) {
      args$env_ids <- as.character(env_ids)
    }
  }
  if (!is.null(env_covariates)) {
    args$env_covariates <- gp_prepare_env_covariates(
      env_covariates,
      source_env_col = heter_groups,
      target_env_col = backend_env_col
    )
    args$reaction_norm_feature_qc <- isTRUE(reaction_norm_feature_qc)
    args$kenv_kernel <- as.character(kenv_kernel)[1L]
    args$kenv_bandwidth <- as.numeric(kenv_bandwidth)[1L]
    if (!is.null(kenv_kernel_kwargs)) {
      args$kenv_kernel_kwargs <- kenv_kernel_kwargs
    }
  }
  if (!is.null(w_g)) args$w_g <- as.numeric(w_g)
  if (!is.null(w_ge)) args$w_ge <- as.numeric(w_ge)
  if (!is.null(w_e)) args$w_e <- as.numeric(w_e)[1L]
  if (identical(gp_backend_method_map(model_name), "gp_icm_fa")) {
    args$fa_rank <- as.integer(gp_fa_rank)
  }

  fit_reticulate <- function() {
    fw <- gp_import_gp_module("gp_framework", python_bin = python_bin, project_root = project_root)
    reticulate_args <- args
    reticulate_args$geno_kernels <- do.call(
      reticulate::dict,
      stats::setNames(lapply(kernel$Ks, unname), names(kernel$Ks))
    )
    gp_public_result_to_r(do.call(fw$fit_mixed_model, reticulate_args))
  }
  if (gp_bridge_direct_execution_enabled()) {
    out <- tryCatch(
      gp_with_gp_device_route(
        policy$device$route,
        gp_bridge_fit_mixed_model_direct(
          args = args,
          kernel = kernel,
          gp_factor_cache = gp_factor_cache,
          python_bin = python_bin,
          project_root = project_root,
          gp_theta0_warm = gp_theta0_warm,
          gp_engine = gp_engine
        )
      ),
      error = function(e) {
        if (!gp_bridge_reticulate_fallback_enabled()) {
          stop(e)
        }
        warning(
          "Direct GP Python bridge failed; falling back to reticulate because ",
          "PREDICTPRO_GP_RETICULATE_FALLBACK is enabled. Original error: ",
          conditionMessage(e),
          call. = FALSE
        )
        gp_with_gp_device_route(policy$device$route, fit_reticulate())
      }
    )
  } else {
    out <- gp_with_gp_device_route(policy$device$route, fit_reticulate())
  }
  variance_component_method <- tryCatch(
    as.character(
      out$result$variance_component_method %||%
        out$result$diagnostics$variance_component_method %||%
        out$diagnostics$variance_component_method
    )[1L],
    error = function(e) NA_character_
  )
  if (identical(variance_component_method, "krr_gcv_variance_ratio")) {
    # KRR estimates a GCV-selected residual/genetic variance ratio; it is not
    # the REML post-fit used by Gaussian-Process-GBLUP.  Surface the effective
    # method so the public variance formatter cannot relabel it as REML.
    out$info$varcomp_mode <- variance_component_method
    if (is.list(out$fit)) {
      out$fit$varcomp_mode <- variance_component_method
    }
  }
  if (identical(final_mean_mode, "location_offset")) {
    out <- gp_apply_prediction_offset(
      out,
      pheno_df = pheno_df,
      split = split,
      offsets = mean_fit$offsets,
      label = final_mean_mode
    )
  }
  blend_info <- list(
    requested = prediction_blend,
    applied = FALSE,
    reason = if (identical(prediction_blend, "none")) "prediction_blend_none" else "tuning_blend_unavailable"
  )
  if (!identical(prediction_blend, "none") &&
      identical(tuning_info$status %||% NA_character_, "ok") &&
      !is.null(tuning_info$selected) &&
      "validation_blend_alpha" %in% names(tuning_info$selected)) {
    blend_alpha <- suppressWarnings(as.numeric(tuning_info$selected$validation_blend_alpha[[1L]]))
    blend_info$validation <- list(
      alpha = blend_alpha,
      objective = tuning_info$selected$validation_blend_objective[[1L]] %||% tuning_objective,
      objective_value = suppressWarnings(as.numeric(tuning_info$selected$validation_blend_objective_value[[1L]])),
      gp_rmse = suppressWarnings(as.numeric(tuning_info$selected$validation_rmse[[1L]])),
      envmean_rmse = suppressWarnings(as.numeric(tuning_info$selected$validation_envmean_rmse[[1L]])),
      blend_rmse = suppressWarnings(as.numeric(tuning_info$selected$validation_blend_rmse[[1L]]))
    )
    if (is.finite(blend_alpha) && blend_alpha < 1 - 1e-8) {
      gp_value <- suppressWarnings(as.numeric(tuning_info$selected$validation_objective[[1L]]))
      blend_value <- suppressWarnings(as.numeric(tuning_info$selected$validation_blend_objective_value[[1L]]))
      direction <- as.character(tuning_info$tuning_objective_direction %||% "minimize")[1L]
      improvement <- if (is.finite(gp_value) && is.finite(blend_value)) {
        denom <- max(abs(gp_value), .Machine$double.eps)
        if (identical(direction, "maximize")) {
          (blend_value - gp_value) / denom
        } else {
          (gp_value - blend_value) / denom
        }
      } else {
        NA_real_
      }
      blend_decision <- gp_prediction_blend_auto_decision(
        prediction_blend = prediction_blend,
        blend_alpha = blend_alpha,
        relative_improvement = improvement,
        allow_auto = isTRUE(policy$profile$prospective_env_prediction) ||
          !isTRUE(policy$profile$is_met)
      )
      blend_info$validation$relative_improvement <- improvement
      blend_info$validation$min_relative_improvement <- blend_decision$min_relative_improvement
      blend_info$validation$max_auto_alpha <- blend_decision$max_auto_alpha
      if (isTRUE(blend_decision$apply)) {
        blended <- gp_apply_prediction_blend(
          out,
          pheno_df = pheno_df,
          split = split,
          alpha = blend_alpha,
          env_location = env_location,
          env_year = env_year
        )
        out <- blended$result
        blend_info <- c(
          list(requested = prediction_blend),
          blended$info,
          list(validation = blend_info$validation)
        )
        blend_info$reason <- blend_decision$reason
      } else {
        blend_info$applied <- FALSE
        blend_info$reason <- blend_decision$reason
        blend_info$alpha <- blend_alpha
      }
    } else {
      blend_info$applied <- FALSE
      blend_info$reason <- "validation_selected_gp_only"
      blend_info$alpha <- blend_alpha
    }
  }
  uncertainty_calibration_info <- list(
    requested = gp_uncertainty_calibration,
    applied = FALSE,
    reason = if (identical(gp_uncertainty_calibration, "none")) {
      "uncertainty_calibration_none"
    } else if (!output_level %in% c("predict_with_se", "full_vc")) {
      "prediction_se_not_requested"
    } else if (!identical(tuning_info$status %||% NA_character_, "ok")) {
      "gp_tuning_unavailable"
    } else {
      "validation_uncertainty_calibration_unavailable"
    }
  )
  if (output_level %in% c("predict_with_se", "full_vc") &&
      !identical(gp_uncertainty_calibration, "none") &&
      identical(tuning_info$status %||% NA_character_, "ok") &&
      is.list(tuning_info$uncertainty_calibration)) {
    calibrated <- gp_apply_uncertainty_calibration(
      out,
      tuning_info$uncertainty_calibration,
      label = paste0("gp_tuning_", tuning_strategy_effective)
    )
    out <- calibrated$result
    uncertainty_calibration_info <- c(
      list(requested = gp_uncertainty_calibration),
      calibrated$info
    )
  }
  out$info$gp_tuning <- tuning_info
  out$info$gp_auto_policy <- policy
  out$info$mean_adjustment <- list(
    requested = mean_adjustment,
    applied = final_mean_mode
  )
  out$info$gp_prediction_blend <- blend_info
  out$info$gp_uncertainty_calibration <- uncertainty_calibration_info
  if (is.list(out$fit)) {
    out$fit$gp_tuning <- tuning_info
    out$fit$gp_auto_policy <- policy
    out$fit$mean_adjustment <- out$info$mean_adjustment
    out$fit$gp_prediction_blend <- blend_info
    out$fit$gp_uncertainty_calibration <- uncertainty_calibration_info
  }
  out
}

#' Fit a Multi-Trait Gaussian-Process Model
#'
#' R wrapper around the single-environment multi-trait GP backend.
#'
#' @param pheno_data Data frame with one genotype row and one column per trait.
#' @param gmatrix Optional genomic relationship matrix. For scalable large-data
#' fits, it may be omitted when \code{gp_factor_cache} and matching
#' \code{geno_ids} are supplied.
#' @param response Character vector of at least two trait columns.
#' @param gen_name Genotype ID column name.
#' @param test_set Optional genotype IDs to predict.
#' @param kernel_list Optional named list of additional relationship/kernel
#' matrices. When supplied, multi-trait GP fits pass all matrices as backend
#' \code{geno_kernels}.
#' @param kernel_weights Optional non-negative weights for the aligned named
#' kernels. Prediction uses their weighted sum. The joint backend estimates one
#' shared trait covariance for that combined kernel; exported per-kernel
#' covariances are fixed-weight decompositions, not independently estimated
#' covariance matrices.
#' @param return_se Logical; request prediction standard errors where supported.
#' @param prediction_output One of \code{"test_only"} or \code{"all"}. The
#' latter returns posterior predictions for both training and test rows while
#' the fitted model continues to condition only on the training observations.
#' @param return_trait_correlations Logical; return genetic and residual trait
#' correlation matrices.
#' @param fixed_effects Optional fixed-effect columns passed to the backend.
#' @param gp_auto_policy Logical; when \code{TRUE}, record and apply the
#' conservative GP device route selected from available GPU/CPU resources.
#' @param tune_gp Logical or \code{"auto"}; when enabled, run a bounded
#' genotype-holdout validation pass for multi-trait RMSE fallback/blending and
#' SE calibration. \code{"auto"} enables it only when enough training rows and
#' genotypes are available.
#' @param prediction_blend One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful validation, \code{"auto"} applies only
#' conservative validation-selected baseline blends.
#' @param gp_uncertainty_calibration One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful validation and requested SE output,
#' scales prediction SE/variance columns using validation residuals.
#' @param tuning_validation_fraction Fraction of training genotypes held out for
#' internal validation.
#' @param tuning_max_calibration_rows Optional row cap for internal calibration
#' fits.
#' @param tuning_max_validation_rows Optional row cap for internal validation
#' predictions.
#' @param trait_fa_rank Optional factor-analytic rank when
#' \code{trait_structure = "fa"}.
#' @param residual_fa_rank Optional factor-analytic rank when
#' \code{residual_structure = "fa"}.
#' @return A list with \code{predictions}, \code{result}, \code{fit}, and
#' \code{info}. The result includes per-kernel variance and covariance tables.
#' Variance output is withheld when the variance-component fit explicitly
#' reports non-convergence.
#' @export
gp_multi_trait_model <- function(pheno_data,
                                 gmatrix = NULL,
                                 response,
                                 gen_name,
                                 test_set = NULL,
                                 geno_ids = NULL,
                                 kernel_list = NULL,
                                 kernel_weights = NULL,
                                 estimate_kernel_weights = FALSE,
                                 force_prediction_se = FALSE,
                                 varcomp_mode = "reml",
                                 gp_factor_cache = NULL,
                                 gp_auto_policy = TRUE,
                                 tune_gp = "auto",
                                 prediction_blend = "auto",
                                 gp_uncertainty_calibration = "auto",
                                 tuning_validation_fraction = 0.25,
                                 tuning_max_calibration_rows = NULL,
                                 tuning_max_validation_rows = NULL,
                                 return_se = TRUE,
                                 prediction_output = "test_only",
                                 return_trait_correlations = FALSE,
                                 fixed_effects = NULL,
                                 trait_structure = "unstructured",
                                 trait_fa_rank = NULL,
                                 residual_structure = "unstructured",
                                 residual_fa_rank = NULL,
                                 max_iter = 100L,
                                 tol_loglik = 1e-4,
                                 seed = 12345L,
                                 python_bin = NULL,
                                 project_root = NULL) {
  long_df <- gp_public_long_traits(pheno_data, response = response, gen_name = gen_name)
  backend_schema <- gp_backend_schema_cols()
  long_cols <- backend_schema$long
  varcomp_mode <- gp_match_gp_choice(
    varcomp_mode,
    choices = c("mom", "reml"),
    name = "varcomp_mode",
    default = "reml"
  )
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    kernel_list = kernel_list
  )
  if (is.null(gmatrix) && length(kernel_inputs) > 0L) {
    gmatrix <- kernel_inputs[[1L]]
  }
  if (is.null(gmatrix) && is.null(gp_factor_cache)) {
    stop(
      "A genomic relationship matrix must be supplied through `gmatrix` unless ",
      "`gp_factor_cache` and `geno_ids` are supplied for `varcomp_mode = \"mom\"`.",
      call. = FALSE
    )
  }
  if (is.null(gmatrix) && !identical(varcomp_mode, "mom")) {
    stop(
      "`gmatrix` can be omitted only when `varcomp_mode = \"mom\"` and a ",
      "`gp_factor_cache` plus matching `geno_ids` are supplied.",
      call. = FALSE
    )
  }
  kernel <- gp_public_kernel_bank_inputs(
    if (length(kernel_inputs) > 0L) kernel_inputs else gmatrix,
    pheno_ids = long_df[[long_cols$gid]],
    geno_ids = geno_ids,
    allow_factor_cache = is.null(gmatrix) && !is.null(gp_factor_cache)
  )
  kernel_reporting_names <- gp_public_kernel_reporting_names(kernel_inputs, kernel)
  weight_kernel <- kernel
  names(weight_kernel$Ks) <- kernel_reporting_names
  kernel_weights <- gp_public_kernel_weights(weight_kernel, kernel_weights)
  kernel <- gp_public_additive_kernel_bank(kernel)
  split <- gp_public_train_test_idx(
    long_df,
    gid_col = long_cols$gid,
    y_col = long_cols$y,
    test_set = test_set
  )
  prediction_blend <- gp_match_gp_choice(
    prediction_blend,
    choices = c("auto", "none", "validation"),
    name = "prediction_blend",
    default = "auto"
  )
  gp_uncertainty_calibration <- gp_match_gp_choice(
    gp_uncertainty_calibration,
    choices = c("auto", "none", "validation"),
    name = "gp_uncertainty_calibration",
    default = "auto"
  )
  prediction_output <- gp_match_gp_choice(
    prediction_output,
    choices = c("test_only", "all"),
    name = "prediction_output",
    default = "test_only"
  )
  args <- list(
    pheno_df = long_df,
    gid_col = long_cols$gid,
    trait_col = long_cols$trait,
    y_col = long_cols$y,
    geno_ids = kernel$geno_ids,
    train_idx = split$train_idx,
    test_idx = split$test_idx,
    varcomp_mode = varcomp_mode,
    fixed_effects = as.list(fixed_effects %||% character(0)),
    trait_structure = trait_structure,
    trait_fa_rank = trait_fa_rank,
    residual_structure = residual_structure,
    residual_fa_rank = residual_fa_rank,
    kernel_weights = as.numeric(kernel_weights),
    estimate_kernel_weights = isTRUE(estimate_kernel_weights) &&
      length(kernel$Ks) > 1L && identical(tolower(varcomp_mode), "reml"),
    max_iter = as.integer(max_iter),
    tol_loglik = as.numeric(tol_loglik),
    return_se = isTRUE(return_se),
    prediction_output = prediction_output,
    return_trait_correlations = isTRUE(return_trait_correlations),
    force_prediction_se = isTRUE(force_prediction_se),
    force_cov_traits = isTRUE(force_prediction_se),
    seed = as.integer(seed)
  )
  if (isTRUE(estimate_kernel_weights) && !isTRUE(args$estimate_kernel_weights)) {
    # Say so rather than silently returning fixed weights the caller did not ask for.
    warning(
      "estimate_kernel_weights = TRUE was ignored: fitting the kernel mixture ",
      "needs varcomp_mode = 'reml' and at least two kernels (got varcomp_mode = '",
      varcomp_mode, "' with ", length(kernel$Ks), " kernel(s)).",
      call. = FALSE
    )
  }
  if (!is.null(gp_factor_cache)) {
    args$grm_factor_cache <- gp_factor_cache
  }
  policy <- list(
    enabled = isTRUE(gp_auto_policy),
    model_family = "multi_trait_single_environment",
    device = gp_policy_device_route(),
    decisions = if (isTRUE(gp_auto_policy)) "device_route_selected" else "auto_policy_disabled",
    profile = list(
      n_rows = as.integer(nrow(long_df)),
      n_train = as.integer(sum(split$train_mask)),
      n_test = as.integer(sum(split$test_mask)),
      n_traits = as.integer(length(unique(as.character(long_df[[long_cols$trait]])))),
      has_factor_cache = !is.null(gp_factor_cache) || identical(tolower(varcomp_mode), "mom")
    )
  )
  require_factor_cache <- identical(tolower(varcomp_mode), "mom")
  fit_reticulate <- function(model_args = args) {
    mt_gp <- gp_import_gp_module("mt_gp", python_bin = python_bin, project_root = project_root)
    reticulate_args <- model_args
    if (isTRUE(require_factor_cache) && is.null(reticulate_args$grm_factor_cache)) {
      reticulate_args$grm_factor_cache <- gp_public_inmemory_factor_cache(
        lapply(kernel$Ks, gp_public_factor_from_kernel),
        python_bin = python_bin,
        project_root = project_root
      )
    }
    reticulate_args$geno_kernels <- do.call(
      reticulate::dict,
      stats::setNames(lapply(kernel$Ks, unname), names(kernel$Ks))
    )
    gp_public_result_to_r(do.call(mt_gp$fit_multi_trait_gp, reticulate_args))
  }
  run_model <- function(model_args = args) {
    if (gp_bridge_direct_execution_enabled()) {
      tryCatch(
        gp_with_gp_device_route(
          if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
          gp_bridge_fit_public_model_direct(
            args = model_args,
            kernel = kernel,
            gp_factor_cache = gp_factor_cache,
            command = "fit-multi-trait-spec",
            python_bin = python_bin,
            project_root = project_root,
            require_factor_cache = require_factor_cache
          )
        ),
        error = function(e) {
          if (!gp_bridge_reticulate_fallback_enabled()) {
            stop(e)
          }
          warning(
            "Direct multi-trait GP Python bridge failed; falling back to reticulate because ",
            "PREDICTPRO_GP_RETICULATE_FALLBACK is enabled. Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
          gp_with_gp_device_route(
            if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
            fit_reticulate(model_args)
          )
        }
      )
    } else {
      gp_with_gp_device_route(
        if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
        fit_reticulate(model_args)
      )
    }
  }
  tune_mode <- gp_policy_bool_or_auto(tune_gp, "tune_gp", default = FALSE)
  n_train_gid <- length(unique(as.character(long_df[[long_cols$gid]][split$train_mask])))
  auto_tune_ok <- isTRUE(gp_auto_policy) &&
    sum(split$train_mask) >= 60L &&
    n_train_gid >= 12L &&
    length(unique(as.character(long_df[[long_cols$trait]]))) >= 2L
  tune_gp_effective <- if (identical(tune_mode, "auto")) isTRUE(auto_tune_ok) else isTRUE(tune_mode)
  validation_adjustment <- list(
    info = list(
      requested = tune_gp_effective,
      requested_value = tune_gp,
      status = if (tune_gp_effective) "not_run" else if (gp_policy_is_auto(tune_gp)) "skipped" else "not_requested",
      reason = if (!tune_gp_effective && gp_policy_is_auto(tune_gp)) {
        if (isTRUE(auto_tune_ok)) NA_character_ else "auto_tune_skipped_too_few_training_rows_or_genotypes"
      } else {
        NA_character_
      }
    ),
    blend = list(requested = prediction_blend, applied = FALSE, reason = "multi_trait_validation_unavailable"),
    calibration = list(
      requested = gp_uncertainty_calibration,
      status = if (identical(gp_uncertainty_calibration, "none")) "skipped" else "not_run",
      reason = if (identical(gp_uncertainty_calibration, "none")) "uncertainty_calibration_none" else "multi_trait_validation_unavailable"
    )
  )
  if (isTRUE(tune_gp_effective)) {
    key_cols <- gp_public_multitrait_key_cols(FALSE)
    validation_adjustment <- gp_public_multitrait_validation_adjustment(
      long_df = long_df,
      split = split,
      key_cols = key_cols,
      fit_validation = function(validation_split) {
        val_args <- args
        val_args$pheno_df <- long_df
        val_args$pheno_df[[long_cols$y]][validation_split$validation_idx + 1L] <- NA_real_
        val_args$train_idx <- validation_split$calibration_idx
        val_args$test_idx <- validation_split$validation_idx
        val_args$return_trait_correlations <- FALSE
        val_args$return_se <- isTRUE(return_se) && !identical(gp_uncertainty_calibration, "none")
        val_args$prediction_output <- "test_only"
        run_model(val_args)
      },
      prediction_blend = prediction_blend,
      gp_uncertainty_calibration = gp_uncertainty_calibration,
      return_se = return_se,
      validation_fraction = tuning_validation_fraction,
      max_calibration_rows = tuning_max_calibration_rows,
      max_validation_rows = tuning_max_validation_rows,
      seed = seed
    )
  }
  out <- gp_public_normalize_multitrait_uncertainty(run_model(args))
  key_cols <- gp_public_multitrait_key_cols(FALSE)
  blend_info <- validation_adjustment$blend
  if (!identical(prediction_blend, "none") && isTRUE(blend_info$applied)) {
    final_baseline <- gp_public_multitrait_baseline(
      train_df = long_df[split$train_mask, , drop = FALSE],
      query_df = out$predictions,
      key_cols = key_cols
    )
    blended <- gp_public_multitrait_apply_blend(
      out,
      baseline = final_baseline,
      key_cols = key_cols,
      alpha = blend_info$blend_alpha %||% blend_info$alpha %||% 1,
      label = "multi_trait_blocked_genotype_validation"
    )
    out <- gp_public_normalize_multitrait_uncertainty(blended$result)
    blend_info <- c(blend_info, blended$info)
  }
  calibration_info <- c(
    list(requested = gp_uncertainty_calibration, applied = FALSE),
    validation_adjustment$calibration
  )
  if (isTRUE(return_se) &&
      !identical(gp_uncertainty_calibration, "none") &&
      identical(validation_adjustment$calibration$status %||% NA_character_, "ok")) {
    calibrated <- gp_apply_uncertainty_calibration(
      out,
      validation_adjustment$calibration,
      label = "multi_trait_blocked_genotype_validation"
    )
    out <- gp_public_normalize_multitrait_uncertainty(calibrated$result)
    calibration_info <- c(
      list(requested = gp_uncertainty_calibration),
      calibrated$info
    )
  }
  validation_adjustment$info$requested <- tune_gp_effective
  validation_adjustment$info$requested_value <- tune_gp
  out$info$gp_auto_policy <- policy
  out$info$gp_tuning <- validation_adjustment$info
  out$info$gp_prediction_blend <- blend_info
  out$info$gp_uncertainty_calibration <- calibration_info
  if (is.list(out$fit)) {
    out$fit$gp_auto_policy <- policy
    out$fit$gp_tuning <- validation_adjustment$info
    out$fit$gp_prediction_blend <- blend_info
    out$fit$gp_uncertainty_calibration <- calibration_info
  }
  out <- gp_public_joint_multikernel_contract(
    out = out,
    kernel = kernel,
    kernel_weights = kernel_weights,
    kernel_names = kernel_reporting_names,
    varcomp_mode = varcomp_mode,
    is_met = FALSE,
    return_trait_correlations = return_trait_correlations
  )
  out
}

#' Fit a Multi-Trait Multi-Environment Gaussian-Process Model
#'
#' R wrapper around the multi-trait multi-environment GP backend.
#'
#' @param pheno_data Data frame with genotype, environment, and trait columns.
#' @param gmatrix Genomic relationship matrix.
#' @param response Character vector of at least two trait columns.
#' @param gen_name Genotype ID column name.
#' @param heter_groups Environment column name.
#' @param env_similarity Optional environment relationship matrix. A
#' covariate-derived environment kernel should be supplied here when available;
#' it replaces any need for an additional environment FA or unstructured
#' environment covariance choice.
#' @param env_covariates Optional environment covariates for backend
#' environment-kernel construction. When supplied, they define the environment
#' kernel and replace any need for an additional environment FA or unstructured
#' environment covariance choice. If both \code{env_similarity} and
#' \code{env_covariates} are omitted, an identity environment matrix is used.
#' Pass exactly one of \code{env_similarity} or \code{env_covariates}.
#' @param reaction_norm_feature_qc Logical; when \code{TRUE}, use backend
#' feature QC before converting environment covariates to an environment kernel.
#' @param kenv_kernel Environment-covariate kernel family used when
#' \code{env_covariates} is supplied.
#' @param kenv_bandwidth Numeric bandwidth used by the environment kernel.
#' @param kenv_kernel_kwargs Optional named list of backend kernel arguments.
#' @param kernel_list Optional named list of additional relationship/kernel
#' matrices. When supplied, multi-trait MET GP fits pass all matrices as
#' backend \code{geno_kernels}.
#' @param kernel_weights Optional non-negative weights for the aligned genomic
#' kernels. Named weights are reordered to the aligned kernel names. The
#' default assigns weight one to every kernel, so the genomic covariance uses
#' their sum. The backend estimates shared trait and GxE covariance
#' coefficients for the combined kernel; exported per-kernel components are
#' fixed-weight decompositions, not independent variance estimates.
#' @param varcomp_mode Variance-component estimation mode passed to the
#' MT-MET GP backend. The default \code{"reml"} reports REML-derived variance
#' components; \code{"mom"} remains available for scalable factor-cache fits.
#' @param return_se Logical; request prediction standard errors.
#' @param prediction_output One of \code{"test_only"} or \code{"all"}. The
#' latter queries every genotype-environment-trait cell from the fitted
#' posterior while the model continues to condition only on training rows.
#' @param return_trait_correlations Logical; return genetic, GxE, and residual
#' trait correlation matrices.
#' @param fixed_effects Optional fixed-effect columns passed to the backend.
#' @param gp_auto_policy Logical; when \code{TRUE}, record and apply
#' conservative GP policy metadata, including GPU/CPU routing and large-MET
#' environment-kernel guidance. Environment covariates, when supplied, remain
#' the environment kernel and do not trigger a second FA or unstructured
#' environment covariance.
#' @param tune_gp Logical or \code{"auto"}; when enabled, run a bounded
#' genotype-holdout validation pass for MT-MET RMSE fallback/blending and SE
#' calibration. \code{"auto"} enables it only when enough training rows,
#' genotypes, and environments are available.
#' @param prediction_blend One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful validation, \code{"auto"} applies only
#' conservative validation-selected baseline blends.
#' @param gp_uncertainty_calibration One of \code{"auto"}, \code{"none"}, or
#' \code{"validation"}. With successful validation and requested SE output,
#' scales prediction SE/variance columns using validation residuals.
#' @param tuning_validation_fraction Fraction of training genotypes held out for
#' internal validation.
#' @param tuning_max_calibration_rows Optional row cap for internal calibration
#' fits.
#' @param tuning_max_validation_rows Optional row cap for internal validation
#' predictions.
#' @param trait_structure Cross-trait covariance structure. Defaults to
#' \code{"unstructured"}. Use \code{"fa"} only as an explicit, benchmarked
#' opt-in. This is not an environment covariance structure.
#' @param trait_fa_rank Optional factor-analytic rank when
#' \code{trait_structure = "fa"}.
#' @param gxe_trait_structure Cross-trait covariance structure for the GxE
#' component. When \code{NULL}, it follows \code{trait_structure}. This is not
#' an environment covariance structure.
#' @param gxe_trait_fa_rank Optional factor-analytic rank when
#' \code{gxe_trait_structure = "fa"}.
#' @return A list with \code{predictions}, \code{result}, \code{fit}, and
#' \code{info}. Response-scale covariance matrices are returned even when trait
#' correlations are not requested. Per-kernel genetic/GxE variance and
#' covariance decompositions are withheld if the backend explicitly reports
#' non-convergence.
#' @export
gp_multi_trait_met_model <- function(pheno_data,
                                     gmatrix = NULL,
                                     response,
                                     gen_name,
                                     heter_groups,
                                     test_set = NULL,
                                     env_similarity = NULL,
                                     env_ids = NULL,
                                     env_covariates = NULL,
                                     reaction_norm_feature_qc = TRUE,
                                     kenv_kernel = "matern32",
                                     kenv_bandwidth = 1.0,
                                     kenv_kernel_kwargs = NULL,
                                     geno_ids = NULL,
                                     kernel_list = NULL,
                                     kernel_weights = NULL,
                                     force_prediction_se = FALSE,
                                     varcomp_mode = "reml",
                                     gp_factor_cache = NULL,
                                     gp_auto_policy = TRUE,
                                     tune_gp = "auto",
                                     prediction_blend = "auto",
                                     gp_uncertainty_calibration = "auto",
                                     tuning_validation_fraction = 0.25,
                                     tuning_max_calibration_rows = NULL,
                                     tuning_max_validation_rows = NULL,
                                     return_se = TRUE,
                                     prediction_output = "test_only",
                                     return_trait_correlations = FALSE,
                                     fixed_effects = NULL,
                                     trait_structure = "unstructured",
                                     trait_fa_rank = NULL,
                                     gxe_trait_structure = NULL,
                                     gxe_trait_fa_rank = NULL,
                                     mom_pcg_tol = 1e-4,
                                     mom_pcg_max_iter = 300L,
                                     dense_reml_max_train = 3000L,
                                     reml_max_iter = 30L,
                                     reml_tol = 1e-5,
                                     seed = 12345L,
                                     python_bin = NULL,
                                     project_root = NULL) {
  if (!is.null(env_similarity) && !is.null(env_covariates)) {
    stop("Pass exactly one of `env_similarity` or `env_covariates`.", call. = FALSE)
  }
  long_df <- gp_public_long_traits(
    pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  backend_schema <- gp_backend_schema_cols()
  single_cols <- backend_schema$single
  long_cols <- backend_schema$long
  varcomp_mode <- gp_match_gp_choice(
    varcomp_mode,
    choices = c("mom", "reml", "reml_dense"),
    name = "varcomp_mode",
    default = "reml"
  )
  kernel_inputs <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    kernel_list = kernel_list
  )
  if (is.null(gmatrix) && length(kernel_inputs) > 0L) {
    gmatrix <- kernel_inputs[[1L]]
  }
  if (is.null(gmatrix) && is.null(gp_factor_cache)) {
    stop(
      "A genomic relationship matrix must be supplied through `gmatrix` unless ",
      "`gp_factor_cache` and `geno_ids` are supplied for scalable MT-MET GP fitting.",
      call. = FALSE
    )
  }
  kernel <- gp_public_kernel_bank_inputs(
    if (length(kernel_inputs) > 0L) kernel_inputs else gmatrix,
    pheno_ids = long_df[[long_cols$gid]],
    geno_ids = geno_ids,
    allow_factor_cache = is.null(gmatrix) && !is.null(gp_factor_cache)
  )
  kernel_reporting_names <- gp_public_kernel_reporting_names(kernel_inputs, kernel)
  weight_kernel <- kernel
  names(weight_kernel$Ks) <- kernel_reporting_names
  kernel_weights <- gp_public_kernel_weights(weight_kernel, kernel_weights)
  kernel <- gp_public_additive_kernel_bank(kernel)
  split <- gp_public_train_test_idx(
    long_df,
    gid_col = long_cols$gid,
    y_col = long_cols$y,
    test_set = test_set
  )
  prediction_blend <- gp_match_gp_choice(
    prediction_blend,
    choices = c("auto", "none", "validation"),
    name = "prediction_blend",
    default = "auto"
  )
  gp_uncertainty_calibration <- gp_match_gp_choice(
    gp_uncertainty_calibration,
    choices = c("auto", "none", "validation"),
    name = "gp_uncertainty_calibration",
    default = "auto"
  )
  prediction_output <- gp_match_gp_choice(
    prediction_output,
    choices = c("test_only", "all"),
    name = "prediction_output",
    default = "test_only"
  )
  supplied_env_kernel <- !is.null(env_similarity) || !is.null(env_covariates)
  gp_warn_large_met_without_env_kernel(
    long_df[[long_cols$env]],
    has_env_kernel = supplied_env_kernel,
    context = "Multi-trait MET GP"
  )
  trait_structure <- gp_match_gp_choice(
    trait_structure,
    choices = c("unstructured", "fa"),
    name = "trait_structure",
    default = "unstructured"
  )
  gxe_trait_structure <- gp_match_gp_choice(
    gxe_trait_structure,
    choices = c("unstructured", "fa"),
    name = "gxe_trait_structure",
    default = trait_structure
  )
  if (varcomp_mode %in% c("reml", "reml_dense") &&
      (!identical(trait_structure, "unstructured") ||
       !identical(gxe_trait_structure, "unstructured"))) {
    stop(
      "Dense MT-MET REML currently supports only unstructured trait and GxE ",
      "covariance. Use `varcomp_mode = \"mom\"` for FA covariance, or keep ",
      "`trait_structure = \"unstructured\"` and ",
      "`gxe_trait_structure = \"unstructured\"` with REML.",
      call. = FALSE
    )
  }
  if (is.null(env_similarity) && is.null(env_covariates)) {
    env_levels <- unique(as.character(long_df[[long_cols$env]]))
    env_similarity <- diag(length(env_levels))
    rownames(env_similarity) <- colnames(env_similarity) <- env_levels
  }
  policy_profile_df <- data.frame(
    .gid = as.character(long_df[[long_cols$gid]]),
    .env = as.character(long_df[[long_cols$env]]),
    stringsAsFactors = FALSE
  )
  names(policy_profile_df) <- unname(c(single_cols$gid, single_cols$env))
  policy_profile <- gp_policy_prediction_profile(
    policy_profile_df,
    split = split,
    heter_groups = heter_groups,
    gid_col = single_cols$gid,
    env_col = single_cols$env,
    has_env_kernel = supplied_env_kernel,
    has_env_covariates = !is.null(env_covariates),
    has_factor_cache = TRUE,
    large_n_threshold = dense_reml_max_train
  )
  policy <- list(
    enabled = isTRUE(gp_auto_policy),
    model_family = "multi_trait_multi_environment",
    device = gp_policy_device_route(),
    decisions = c(
      if (!isTRUE(gp_auto_policy)) "auto_policy_disabled" else "device_route_selected",
      if (isTRUE(policy_profile$has_env_covariates)) "env_covariates_define_environment_kernel",
      if (!isTRUE(policy_profile$has_env_kernel) && policy_profile$n_env > 5L) {
        "large_met_without_environment_kernel_warned"
      }
    ),
    profile = policy_profile
  )
  args <- list(
    pheno_df = long_df,
    gid_col = long_cols$gid,
    env_col = long_cols$env,
    trait_col = long_cols$trait,
    y_col = long_cols$y,
    geno_ids = kernel$geno_ids,
    train_idx = split$train_idx,
    test_idx = split$test_idx,
    varcomp_mode = varcomp_mode,
    fixed_effects = as.list(fixed_effects %||% character(0)),
    trait_structure = trait_structure,
    trait_fa_rank = trait_fa_rank,
    gxe_trait_structure = gxe_trait_structure,
    gxe_trait_fa_rank = gxe_trait_fa_rank,
    kernel_weights = as.numeric(kernel_weights),
    mom_pcg_tol = as.numeric(mom_pcg_tol),
    mom_pcg_max_iter = as.integer(mom_pcg_max_iter),
    dense_reml_max_train = as.integer(dense_reml_max_train),
    reml_max_iter = as.integer(reml_max_iter),
    reml_tol = as.numeric(reml_tol),
    return_se = isTRUE(return_se),
    prediction_output = prediction_output,
    return_trait_correlations = isTRUE(return_trait_correlations),
    seed = as.integer(seed)
  )
  if (!is.null(gp_factor_cache)) {
    args$grm_factor_cache <- gp_factor_cache
  }
  if (!is.null(env_similarity)) {
    env_similarity <- as.matrix(env_similarity)
    storage.mode(env_similarity) <- "double"
    args$env_similarity <- unname(env_similarity)
    if (is.null(env_ids) && !is.null(rownames(env_similarity))) {
      env_ids <- rownames(env_similarity)
    }
    if (!is.null(env_ids)) {
      args$env_ids <- as.character(env_ids)
    }
  }
  if (!is.null(env_covariates)) {
    args$env_covariates <- gp_prepare_env_covariates(
      env_covariates,
      source_env_col = heter_groups,
      target_env_col = long_cols$env
    )
    args$reaction_norm_feature_qc <- isTRUE(reaction_norm_feature_qc)
    args$kenv_kernel <- as.character(kenv_kernel)[1L]
    args$kenv_bandwidth <- as.numeric(kenv_bandwidth)[1L]
    if (!is.null(kenv_kernel_kwargs)) {
      args$kenv_kernel_kwargs <- kenv_kernel_kwargs
    }
  }
  fit_reticulate <- function(model_args = args) {
    mt_gp_met <- gp_import_gp_module("mt_gp_met", python_bin = python_bin, project_root = project_root)
    reticulate_args <- model_args
    if (is.null(reticulate_args$grm_factor_cache)) {
      reticulate_args$grm_factor_cache <- gp_public_inmemory_factor_cache(
        lapply(kernel$Ks, gp_public_factor_from_kernel),
        python_bin = python_bin,
        project_root = project_root
      )
    }
    reticulate_args$geno_kernels <- do.call(
      reticulate::dict,
      stats::setNames(lapply(kernel$Ks, unname), names(kernel$Ks))
    )
    gp_public_result_to_r(do.call(mt_gp_met$fit_multi_trait_gp_met, reticulate_args))
  }
  run_model <- function(model_args = args) {
    if (gp_bridge_direct_execution_enabled()) {
      tryCatch(
        gp_with_gp_device_route(
          if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
          gp_bridge_fit_public_model_direct(
            args = model_args,
            kernel = kernel,
            gp_factor_cache = gp_factor_cache,
            command = "fit-multi-trait-met-spec",
            python_bin = python_bin,
            project_root = project_root,
            require_factor_cache = TRUE
          )
        ),
        error = function(e) {
          if (!gp_bridge_reticulate_fallback_enabled()) {
            stop(e)
          }
          warning(
            "Direct multi-trait MET GP Python bridge failed; falling back to reticulate because ",
            "PREDICTPRO_GP_RETICULATE_FALLBACK is enabled. Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
          gp_with_gp_device_route(
            if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
            fit_reticulate(model_args)
          )
        }
      )
    } else {
      gp_with_gp_device_route(
        if (isTRUE(gp_auto_policy)) policy$device$route else NA_character_,
        fit_reticulate(model_args)
      )
    }
  }
  tune_mode <- gp_policy_bool_or_auto(tune_gp, "tune_gp", default = FALSE)
  n_train_gid <- length(unique(as.character(long_df[[long_cols$gid]][split$train_mask])))
  n_train_env <- length(unique(as.character(long_df[[long_cols$env]][split$train_mask])))
  auto_tune_ok <- isTRUE(gp_auto_policy) &&
    sum(split$train_mask) >= 80L &&
    n_train_gid >= 12L &&
    n_train_env >= 2L &&
    length(unique(as.character(long_df[[long_cols$trait]]))) >= 2L
  tune_gp_effective <- if (identical(tune_mode, "auto")) isTRUE(auto_tune_ok) else isTRUE(tune_mode)
  validation_adjustment <- list(
    info = list(
      requested = tune_gp_effective,
      requested_value = tune_gp,
      status = if (tune_gp_effective) "not_run" else if (gp_policy_is_auto(tune_gp)) "skipped" else "not_requested",
      reason = if (!tune_gp_effective && gp_policy_is_auto(tune_gp)) {
        if (isTRUE(auto_tune_ok)) NA_character_ else "auto_tune_skipped_too_few_training_rows_genotypes_or_envs"
      } else {
        NA_character_
      }
    ),
    blend = list(requested = prediction_blend, applied = FALSE, reason = "multi_trait_met_validation_unavailable"),
    calibration = list(
      requested = gp_uncertainty_calibration,
      status = if (identical(gp_uncertainty_calibration, "none")) "skipped" else "not_run",
      reason = if (identical(gp_uncertainty_calibration, "none")) "uncertainty_calibration_none" else "multi_trait_met_validation_unavailable"
    )
  )
  if (isTRUE(tune_gp_effective)) {
    key_cols <- gp_public_multitrait_key_cols(TRUE)
    validation_adjustment <- gp_public_multitrait_validation_adjustment(
      long_df = long_df,
      split = split,
      key_cols = key_cols,
      fit_validation = function(validation_split) {
        val_args <- args
        val_args$pheno_df <- long_df
        val_args$pheno_df[[long_cols$y]][validation_split$validation_idx + 1L] <- NA_real_
        val_args$train_idx <- validation_split$calibration_idx
        val_args$test_idx <- validation_split$validation_idx
        val_args$return_trait_correlations <- FALSE
        val_args$return_se <- isTRUE(return_se) && !identical(gp_uncertainty_calibration, "none")
        val_args$prediction_output <- "test_only"
        run_model(val_args)
      },
      prediction_blend = prediction_blend,
      gp_uncertainty_calibration = gp_uncertainty_calibration,
      return_se = return_se,
      validation_fraction = tuning_validation_fraction,
      max_calibration_rows = tuning_max_calibration_rows,
      max_validation_rows = tuning_max_validation_rows,
      seed = seed
    )
  }
  out <- gp_public_normalize_multitrait_uncertainty(run_model(args))
  key_cols <- gp_public_multitrait_key_cols(TRUE)
  blend_info <- validation_adjustment$blend
  if (!identical(prediction_blend, "none") && isTRUE(blend_info$applied)) {
    final_baseline <- gp_public_multitrait_baseline(
      train_df = long_df[split$train_mask, , drop = FALSE],
      query_df = out$predictions,
      key_cols = key_cols
    )
    blended <- gp_public_multitrait_apply_blend(
      out,
      baseline = final_baseline,
      key_cols = key_cols,
      alpha = blend_info$blend_alpha %||% blend_info$alpha %||% 1,
      label = "multi_trait_met_blocked_genotype_validation"
    )
    out <- gp_public_normalize_multitrait_uncertainty(blended$result)
    blend_info <- c(blend_info, blended$info)
  }
  calibration_info <- c(
    list(requested = gp_uncertainty_calibration, applied = FALSE),
    validation_adjustment$calibration
  )
  if (isTRUE(return_se) &&
      !identical(gp_uncertainty_calibration, "none") &&
      identical(validation_adjustment$calibration$status %||% NA_character_, "ok")) {
    calibrated <- gp_apply_uncertainty_calibration(
      out,
      validation_adjustment$calibration,
      label = "multi_trait_met_blocked_genotype_validation"
    )
    out <- gp_public_normalize_multitrait_uncertainty(calibrated$result)
    calibration_info <- c(
      list(requested = gp_uncertainty_calibration),
      calibrated$info
    )
  }
  validation_adjustment$info$requested <- tune_gp_effective
  validation_adjustment$info$requested_value <- tune_gp
  out$info$gp_auto_policy <- policy
  out$info$gp_tuning <- validation_adjustment$info
  out$info$gp_prediction_blend <- blend_info
  out$info$gp_uncertainty_calibration <- calibration_info
  if (is.list(out$fit)) {
    out$fit$gp_auto_policy <- policy
    out$fit$gp_tuning <- validation_adjustment$info
    out$fit$gp_prediction_blend <- blend_info
    out$fit$gp_uncertainty_calibration <- calibration_info
  }
  out <- gp_public_joint_multikernel_contract(
    out = out,
    kernel = kernel,
    kernel_weights = kernel_weights,
    kernel_names = kernel_reporting_names,
    varcomp_mode = varcomp_mode,
    is_met = TRUE,
    return_trait_correlations = return_trait_correlations
  )
  out
}

gp_lowrank_gaussian_model <- function(pheno_data,
                                      response,
                                      gen_name,
                                      gmatrix = NULL,
                                      omic1_kernel = NULL,
                                      omic2_kernel = NULL,
                                      omic3_kernel = NULL,
                                      kernel_list = NULL,
                                      test_set = NULL,
                                      lowrank_eps_trace = 1e-6,
                                      lowrank_max_rank = NULL,
                                      lowrank_jitter = NULL,
                                      lowrank_noise_grid = NULL,
                                      lowrank_kernel_weights = NULL,
                                      system_database = FALSE,
                                      gp_engine = "auto",
                                      observation_weights = NULL) {
  gp_backend_gaussian_model(
    model_name = "LowRankGP",
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    test_set = test_set,
    lowrank_eps_trace = lowrank_eps_trace,
    lowrank_max_rank = lowrank_max_rank,
    lowrank_jitter = lowrank_jitter,
    lowrank_noise_grid = lowrank_noise_grid,
    lowrank_kernel_weights = lowrank_kernel_weights,
    system_database = system_database,
    gp_engine = gp_engine,
    observation_weights = observation_weights
  )
}

gp_lowrank_cv_predict <- function(y, tst, additional_params = NULL) {
  gp_bridge_cv_predict("LowRankGP", y = y, tst = tst, additional_params = additional_params)
}

gp_backend_cv_predict <- function(model_name, y, tst, additional_params = NULL) {
  gp_bridge_cv_predict(model_name = model_name, y = y, tst = tst, additional_params = additional_params)
}

#' Build a reusable GP factor cache
#'
#' Writes a relationship/kernel matrix to disk and asks the Python GP backend to
#' build a low-rank factor cache for repeated prediction workflows.
#'
#' @param grm Relationship or kernel matrix.
#' @param out_dir Output directory for the factor cache.
#' @param eps_trace Trace approximation tolerance.
#' @param max_rank Maximum factor rank.
#' @param store_dtype Storage dtype used by the Python backend.
#' @param overwrite Logical; overwrite an existing cache when \code{TRUE}.
#'
#' @return Invisibly returns a list describing the cache type and root path.
#' @export
build_gp_factor_cache <- function(grm,
                                  out_dir,
                                  eps_trace = 1e-4,
                                  max_rank = 4096L,
                                  store_dtype = "float32",
                                  overwrite = FALSE) {
  python_bin <- gp_detect_gp_python()
  if (is.null(python_bin)) {
    stop("No configured GP Python runtime found. Run setup_predictgp_env() first.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  input_dir <- file.path(out_dir, "_input")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  matrix_bin <- file.path(input_dir, "grm.bin")
  matrix_meta <- file.path(input_dir, "grm_shape.txt")
  gp_bridge_write_matrix_bin(grm, matrix_bin, matrix_meta)
  cli_args <- c(
    "--project-root", gp_project_root(),
    "--matrix-bin", matrix_bin,
    "--matrix-meta", matrix_meta,
    "--out-dir", normalizePath(out_dir, winslash = "/", mustWork = TRUE),
    "--eps-trace", as.character(eps_trace),
    "--max-rank", as.character(as.integer(max_rank)),
    "--store-dtype", as.character(store_dtype)
  )
  if (isTRUE(overwrite)) {
    cli_args <- c(cli_args, "--overwrite")
  }
  gp_bridge_run_cli("build-factor-cache", cli_args, python_bin = python_bin)
  invisible(list(
    type = "zarr",
    root = normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  ))
}
