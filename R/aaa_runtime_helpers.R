`%||%` <- function(a, b) if (is.null(a)) b else a

# Folder that run results are exported into: options(PredictProR.output_dir),
# else the PREDICTPRO_OUTPUT_DIR environment variable, else the working
# directory (the long-standing default).
gp_output_root <- function(create = FALSE) {
  root <- getOption("PredictProR.output_dir", NULL)
  if (is.null(root) || !nzchar(root)) root <- Sys.getenv("PREDICTPRO_OUTPUT_DIR", unset = "")
  if (!nzchar(root)) return(getwd())
  root <- path.expand(as.character(root)[1L])
  if (create && !dir.exists(root)) dir.create(root, recursive = TRUE, showWarnings = FALSE)
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

gp_parse_cli_key_values <- function(x) {
  lines <- as.character(x %||% character())
  lines <- lines[nzchar(lines)]
  kv <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(kv, function(item) item[[1L]] %||% "", character(1))
  vals <- vapply(
    kv,
    function(item) if (length(item) >= 2L) paste(item[-1L], collapse = "=") else "",
    character(1)
  )
  as.list(stats::setNames(vals[nzchar(keys)], keys[nzchar(keys)]))
}

gp_cli_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) {
    return(isTRUE(default))
  }
  tolower(trimws(as.character(x[[1L]]))) %in% c("true", "1", "yes", "y")
}

gp_cli_int <- function(x, default = 0L) {
  val <- suppressWarnings(as.integer(x[[1L]] %||% default))
  if (is.na(val)) {
    return(as.integer(default))
  }
  val
}

gp_lowrank_supported_models <- function() c("KRR", "GP", "GP_FA", "LowRankGP")

gp_single_response_single_environment_gp_supported_models <- function() {
  c("KRR", "GP")
}

gp_single_trait_multi_environment_gp_supported_models <- function() {
  c("KRR", "GP", "GP_FA")
}

gp_single_trait_gp_alias_models <- function() "LowRankGP"

gp_validate_single_trait_gp_model_scope <- function(GS_model = NULL,
                                                     GS_model_cv = NULL,
                                                     multi_trait_gp = FALSE,
                                                     hybrid_gp = FALSE) {
  if (isTRUE(multi_trait_gp) || isTRUE(hybrid_gp)) {
    return(invisible(TRUE))
  }
  requested <- unique(stats::na.omit(c(GS_model, GS_model_cv)))
  requested <- gp_canonicalize_supported_model_names(requested)
  duplicate <- intersect(requested, gp_single_trait_gp_alias_models())
  if (length(duplicate)) {
    stop(
      "Scalable-GBLUP is not supported as a separate single-trait GP model ",
      "because its single-trait route resolves to Kernel-GBLUP ",
      "(`KRR` / `krr_exact`) and produces identical predictions and ",
      "uncertainty. Use Kernel-GBLUP. Scalable-GBLUP remains available for ",
      "joint multi-trait GP analysis, where it selects the distinct ",
      "operator/MoM route.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Models that do not work on polyploid dosage (ploidy > 2). GPNet: its
# deep-kernel features collapse on dosage 0..ploidy and it predicts a constant.
gp_reject_polyploid_unsupported_models <- function(models, ploidy) {
  p <- suppressWarnings(as.integer(ploidy)[1L])
  if (!length(p) || is.na(p) || p <= 2L) {
    return(invisible(TRUE))
  }
  models <- unique(as.character(stats::na.omit(unlist(models, use.names = FALSE))))
  if (!length(models)) {
    return(invisible(TRUE))
  }
  canonical <- tryCatch(gp_canonicalize_supported_model_names(models), error = function(e) models)
  bad <- models[canonical %in% "gp_dkl" | models %in% "GPNet"]
  if (length(bad)) {
    stop(
      "GPNet is not supported for polyploid data (ploidy ", p, "): its deep-kernel features ",
      "collapse on polyploid dosage and it predicts a constant. Use another deep-learning model ",
      "(e.g. DenseNeuralNet, TabAttention, TabNet) or a GP model (Kernel-GBLUP, Gaussian-Process-GBLUP).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

gp_lowrank_friendly_models <- function() {
  c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")
}

gp_single_response_single_environment_gp_friendly_models <- function() {
  gp_display_supported_model_names(
    gp_single_response_single_environment_gp_supported_models()
  )
}

gp_deep_learning_supported_models <- function() {
  c(
    "cnn", "ft_transformer", "saint", "tabnet", "node", "deepfm", "dcnv2",
    "nam", "moe", "gp_dkl", "mlp_with_attention", "mlp", "resnet"
  )
}

gp_deep_learning_friendly_models <- function() {
  c(
    "Conv1DNet", "TabTransformer", "TabAttention", "TabNet", "LightTreeNet",
    "FactorNet", "CrossNet", "NeuralAdditive", "MixtureOfExperts", "GPNet",
    "DenseAttentionNet", "DenseNeuralNet", "ResNet"
  )
}

gp_classical_ml_supported_models <- function() {
  c(
    "Xgboost", "RandomForest", "CatBoost", "LightGBM", "PartialLeastSquare",
    "SupportVectorMachine", "K-NearestNeighbors", "Lasso", "Ridge_Regression"
  )
}

# Variance-component capability is a property of the fitted model, not of the
# workflow that happens to call it. Keep this inventory central so the same
# rule is used for single-trait, multi-trait, MET, multi-omics, feature-
# selection, and hybrid output boundaries.
gp_predictive_only_variance_models <- function() {
  unique(c(
    gp_classical_ml_supported_models(),
    gp_deep_learning_supported_models(),
    gp_deep_learning_friendly_models()
  ))
}

gp_model_based_variance_models <- function() {
  unique(c(
    "GBLUP", "GBLUP_BRR", "RKHS",
    "BRR", "BayesA", "BayesB", "BayesC", "BL",
    gp_lowrank_supported_models(),
    gp_lowrank_friendly_models()
  ))
}

gp_variance_model_key <- function(x) {
  x <- tolower(trimws(as.character(x %||% character())))
  gsub("[^a-z0-9]+", "", x)
}

gp_model_variance_role <- function(model = NULL) {
  model <- as.character(model %||% character())
  model <- model[!is.na(model) & nzchar(trimws(model))]
  if (!length(model)) {
    return("unknown")
  }
  model <- gp_canonicalize_supported_model_names(model[[1L]])
  key <- gp_variance_model_key(model[[1L]])
  predictive_keys <- gp_variance_model_key(gp_predictive_only_variance_models())
  model_based_keys <- gp_variance_model_key(gp_model_based_variance_models())
  if (key %in% predictive_keys) {
    return("predictive_only")
  }
  if (key %in% model_based_keys) {
    return("model_based")
  }
  "unknown"
}

gp_canonicalize_supported_model_names <- function(x,
                                                  include_dl = TRUE,
                                                  include_gp = TRUE) {
  if (is.null(x) || length(x) == 0L) {
    return(NULL)
  }

  out <- as.character(x)
  if (isTRUE(include_dl)) {
    dl_canonical <- gp_deep_learning_supported_models()
    dl_friendly <- gp_deep_learning_friendly_models()
    if (length(dl_canonical) > 0L && length(dl_canonical) == length(dl_friendly)) {
      out <- replace_with_canonical(
        out,
        canonical_names = dl_canonical,
        friendly_names = dl_friendly
      )
    }
  }
  if (isTRUE(include_gp)) {
    gp_canonical <- gp_lowrank_supported_models()
    gp_friendly <- gp_lowrank_friendly_models()
    if (length(gp_canonical) > 0L && length(gp_canonical) == length(gp_friendly)) {
      out <- replace_with_canonical(
        out,
        canonical_names = gp_canonical,
        friendly_names = gp_friendly
      )
    }
  }
  out
}

gp_display_supported_model_names <- function(x,
                                             include_dl = TRUE,
                                             include_gp = TRUE) {
  if (is.null(x) || length(x) == 0L) {
    return(NULL)
  }

  out <- gp_canonicalize_supported_model_names(
    x,
    include_dl = include_dl,
    include_gp = include_gp
  )
  if (isTRUE(include_dl)) {
    dl_canonical <- gp_deep_learning_supported_models()
    dl_friendly <- gp_deep_learning_friendly_models()
    if (length(dl_canonical) > 0L && length(dl_canonical) == length(dl_friendly)) {
      out <- replace_with_friendly_name(
        out,
        setNames(dl_friendly, dl_canonical)
      )
    }
  }
  if (isTRUE(include_gp)) {
    gp_canonical <- gp_lowrank_supported_models()
    gp_friendly <- gp_lowrank_friendly_models()
    if (length(gp_canonical) > 0L && length(gp_canonical) == length(gp_friendly)) {
      out <- replace_with_friendly_name(
        out,
        setNames(gp_friendly, gp_canonical)
      )
    }
  }
  out
}

gp_all_model_display_names <- function(AI_valid_models = character(),
                                       bayes_valid_models = character(),
                                       bayes_gblup_valid_models = character(),
                                       gp_valid_models = NULL,
                                       asreml_model = character()) {
  gp_valid_models <- gp_valid_models %||% gp_lowrank_supported_models()
  unique(c(
    gp_display_supported_model_names(AI_valid_models),
    bayes_valid_models,
    bayes_gblup_valid_models,
    gp_display_supported_model_names(gp_valid_models),
    asreml_model
  ))
}

gp_drop_empty_model_names <- function(x) {
  x <- as.character(x %||% character())
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

gp_multi_environment_kernel_models <- function(bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
                                               asreml_model = "GBLUP",
                                               gp_valid_models = NULL) {
  gp_valid_models <- gp_valid_models %||%
    gp_single_trait_multi_environment_gp_supported_models()
  gp_valid_models <- intersect(
    gp_canonicalize_supported_model_names(gp_valid_models),
    gp_single_trait_multi_environment_gp_supported_models()
  )
  gp_drop_empty_model_names(c(
    bayes_gblup_valid_models,
    asreml_model,
    gp_valid_models
  ))
}

gp_multi_environment_kernel_display_names <- function(bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
                                                      asreml_model = "GBLUP",
                                                      gp_valid_models = NULL) {
  gp_display_supported_model_names(
    gp_multi_environment_kernel_models(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    )
  )
}

gp_multi_environment_supported_models <- function(bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
                                                  asreml_model = "GBLUP",
                                                  gp_valid_models = NULL,
                                                  include_met_ml_dl = TRUE) {
  gp_drop_empty_model_names(c(
    gp_multi_environment_kernel_models(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    ),
    if (isTRUE(include_met_ml_dl)) gp_met_supported_models() else character()
  ))
}

gp_multi_environment_supported_display_names <- function(bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
                                                         asreml_model = "GBLUP",
                                                         gp_valid_models = NULL,
                                                         include_met_ml_dl = TRUE) {
  gp_display_supported_model_names(
    gp_multi_environment_supported_models(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models,
      include_met_ml_dl = include_met_ml_dl
    )
  )
}

gp_model_friendly_lookup <- function(include_dl = TRUE, include_gp = TRUE) {
  out <- character()
  if (isTRUE(include_dl)) {
    dl_canonical <- gp_deep_learning_supported_models()
    dl_friendly <- gp_deep_learning_friendly_models()
    if (length(dl_canonical) > 0L && length(dl_canonical) == length(dl_friendly)) {
      out <- c(out, setNames(dl_friendly, dl_canonical), setNames(dl_friendly, dl_friendly))
    }
  }
  if (isTRUE(include_gp)) {
    gp_canonical <- gp_lowrank_supported_models()
    gp_friendly <- gp_lowrank_friendly_models()
    if (length(gp_canonical) > 0L && length(gp_canonical) == length(gp_friendly)) {
      out <- c(out, setNames(gp_friendly, gp_canonical), setNames(gp_friendly, gp_friendly))
    }
  }
  out
}

gp_public_model_label <- function(model) {
  label <- gp_display_supported_model_names(model)
  if (is.null(label) || length(label) == 0L) {
    return(as.character(model))
  }
  unname(label)
}

gp_relabel_public_model_parameters <- function(model_parameters,
                                               model,
                                               stat = "gp_model",
                                               canonical_stat = "gp_model_canonical") {
  if (!is.data.frame(model_parameters) ||
      !"stat" %in% names(model_parameters) ||
      !"summary" %in% names(model_parameters)) {
    return(model_parameters)
  }
  model_label <- gp_public_model_label(model)[1]
  model_canonical <- gp_canonicalize_supported_model_names(model)[1]
  stat_idx <- which(as.character(model_parameters[["stat"]]) == stat)
  if (length(stat_idx)) {
    model_parameters[["summary"]][stat_idx] <- model_label
  }
  if (!is.na(model_canonical) &&
      !identical(model_label, model_canonical) &&
      !any(as.character(model_parameters[["stat"]]) == canonical_stat)) {
    model_parameters <- rbind(
      model_parameters,
      data.frame(stat = canonical_stat, summary = model_canonical, stringsAsFactors = FALSE)
    )
  }
  model_parameters
}

gp_annotate_requested_model_parameters <- function(model_parameters, model) {
  if (!is.data.frame(model_parameters) ||
      !all(c("stat", "summary") %in% names(model_parameters))) {
    return(model_parameters)
  }
  if (!any(as.character(model_parameters[["stat"]]) == "requested_model")) {
    model_parameters <- rbind(
      model_parameters,
      data.frame(
        stat = "requested_model",
        summary = gp_public_model_label(model)[1L],
        stringsAsFactors = FALSE
      )
    )
  }
  model_parameters
}

gp_arrange_grob_safely <- function(..., width = 10, height = 6) {
  need_device <- identical(as.integer(grDevices::dev.cur()), 1L)
  tmp <- NULL
  if (isTRUE(need_device)) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(file = tmp, width = width, height = height, onefile = TRUE)
    on.exit({
      if (!identical(as.integer(grDevices::dev.cur()), 1L)) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      if (!is.null(tmp)) {
        unlink(tmp, force = TRUE)
      }
    }, add = TRUE)
  }
  gridExtra::arrangeGrob(...)
}

gp_hybrid_asreml_supported_models <- function() c("GBLUP")

gp_hybrid_bayes_supported_models <- function() c("GBLUP_BRR", "RKHS")

gp_hybrid_gp_supported_models <- function() "GP"

gp_hybrid_ml_supported_models <- function() {
  c(
    "RandomForest",
    "Ridge_Regression",
    "PartialLeastSquare",
    "Lasso",
    "SupportVectorMachine",
    "Xgboost",
    "CatBoost"
  )
}

gp_hybrid_dl_supported_models <- function() c("mlp", "ft_transformer")

gp_multitrait_asreml_supported_models <- function() c("GBLUP")

gp_multitrait_bayes_supported_models <- function() c("GBLUP_BRR", "RKHS")

gp_multitrait_gp_supported_models <- function() c("GP", "GP_FA", "LowRankGP")

gp_kernel_models_for_input_preparation <- function(base_models = character(),
                                                   multi_trait_gp = FALSE) {
  unique(c(
    as.character(base_models %||% character()),
    if (isTRUE(multi_trait_gp)) gp_multitrait_gp_supported_models() else character()
  ))
}

gp_multitrait_gp_route_spec <- function(model_name,
                                        is_met = FALSE,
                                        fa_rank = 1L) {
  model_name <- gp_canonicalize_supported_model_names(model_name)
  supported <- gp_multitrait_gp_supported_models()
  if (length(model_name) != 1L || !model_name %in% supported) {
    stop(
      "Multi-trait GP supports only models with a distinct joint-trait ",
      "implementation: ",
      paste(gp_display_supported_model_names(supported), collapse = ", "),
      ". Kernel-GBLUP is not a multi-trait model in PredictProR.",
      call. = FALSE
    )
  }

  fa_rank <- suppressWarnings(as.integer(fa_rank)[1L])
  if (identical(model_name, "GP_FA") &&
      (!is.finite(fa_rank) || is.na(fa_rank) || fa_rank < 1L)) {
    stop("FA-GBLUP requires `gp_fa_rank` to be a positive integer.", call. = FALSE)
  }

  if (identical(model_name, "GP")) {
    return(list(
      model = model_name,
      backend_method = "joint_gp_ai_reml",
      varcomp_mode = "reml",
      trait_structure = "unstructured",
      trait_fa_rank = NULL,
      gxe_trait_structure = "unstructured",
      gxe_trait_fa_rank = NULL
    ))
  }
  if (identical(model_name, "GP_FA")) {
    return(list(
      model = model_name,
      backend_method = if (isTRUE(is_met)) "joint_gp_fa_mom" else "joint_gp_fa_ai_reml",
      varcomp_mode = if (isTRUE(is_met)) "mom" else "reml",
      trait_structure = "fa",
      trait_fa_rank = fa_rank,
      gxe_trait_structure = if (isTRUE(is_met)) "fa" else NULL,
      gxe_trait_fa_rank = if (isTRUE(is_met)) fa_rank else NULL
    ))
  }
  list(
    model = model_name,
    backend_method = "joint_gp_operator_mom",
    varcomp_mode = "mom",
    trait_structure = "unstructured",
    trait_fa_rank = NULL,
    gxe_trait_structure = "unstructured",
    gxe_trait_fa_rank = NULL
  )
}

gp_multitrait_ml_supported_models <- function() {
  c(
    "RandomForest",
    "Ridge_Regression",
    "PartialLeastSquare",
    "Lasso",
    "SupportVectorMachine",
    "Xgboost",
    "CatBoost"
  )
}

gp_multitrait_dl_supported_models <- function() c("mlp", "ft_transformer")

gp_met_supported_models <- function() {
  c(
    "CatBoost", "LightGBM", "Xgboost", "RandomForest",
    "mlp", "ft_transformer", "saint", "tabnet", "moe"
  )
}

gp_is_met_ml_dl_model <- function(model) {
  !is.null(model) && any(as.character(model) %in% gp_met_supported_models())
}

PredictProR_asreml_data_registry <- new.env(parent = emptyenv())

asreml_model_data_register <- function(pheno_data) {
  key <- paste0(
    "predictpror_asreml_data_",
    format(Sys.time(), "%Y%m%d%H%M%OS6"),
    "_",
    sprintf("%06d", sample.int(999999L, 1L))
  )
  assign(key, pheno_data, envir = PredictProR_asreml_data_registry)
  key
}

asreml_model_data_lookup <- function(key) {
  get(key, envir = PredictProR_asreml_data_registry, inherits = FALSE)
}

resolve_model_name <- function(x, name_lookup) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }

  res <- name_lookup[x]
  res <- res[!is.na(res)]
  if (length(res) == 0) {
    return(NULL)
  }

  unname(res)
}

replace_with_canonical <- function(x, canonical_names, friendly_names) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }

  k <- tolower(trimws(x))
  keys <- tolower(c(canonical_names, friendly_names))
  vals <- rep(canonical_names, 2)
  lookup <- setNames(vals, keys)
  res <- lookup[k]
  res[is.na(res)] <- x[is.na(res)]

  unname(res)
}

replace_with_friendly_name <- function(x, friendly_name_lookup) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  if (is.null(friendly_name_lookup) || length(friendly_name_lookup) == 0) {
    return(unname(x))
  }

  k <- tolower(trimws(x))
  names(friendly_name_lookup) <- tolower(names(friendly_name_lookup))
  res <- friendly_name_lookup[k]
  res[is.na(res)] <- x[is.na(res)]

  unname(res)
}

gp_dl_model_names <- function() {
  c(
    "cnn", "ft_transformer", "saint", "tabnet", "node",
    "deepfm", "dcnv2", "nam", "moe", "gp_dkl",
    "mlp_with_attention", "mlp", "resnet"
  )
}

gp_python_purpose_for_models <- function(models = NULL) {
  models <- unique(stats::na.omit(as.character(models %||% character())))
  if (!length(models)) {
    return("general")
  }
  if (any(models %in% gp_lowrank_supported_models())) {
    return("gp")
  }
  if (any(models %in% gp_dl_model_names())) {
    return("dl")
  }
  "general"
}

gp_detect_python <- function() {
  gp_env <- Sys.getenv("PREDICTPRO_GP_PYTHON", unset = "")
  if (nzchar(gp_env) && file.exists(gp_env)) {
    return(gp_env)
  }
  py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(py) && file.exists(py)) {
    return(py)
  }

  preferred <- gp_preferred_python()
  if (!is.null(preferred)) {
    return(preferred)
  }

  NULL
}

gp_virtualenv_python_path <- function(root, env_name, os_type = .Platform$OS.type) {
  if (is.null(root) || !nzchar(root) || is.null(env_name) || !nzchar(env_name)) {
    return(NULL)
  }

  file.path(root, env_name, gp_python_executable_relpath(os_type))
}

gp_virtualenv_python_candidates <- function(root, env_names = NULL, os_type = .Platform$OS.type) {
  if (is.null(root) || !nzchar(root) || !dir.exists(root)) {
    return(character())
  }

  env_dirs <- tryCatch(list.dirs(root, recursive = FALSE, full.names = TRUE), error = function(e) character())
  if (length(env_names)) {
    env_dirs <- env_dirs[basename(env_dirs) %in% env_names]
  }
  if (!length(env_dirs)) {
    return(character())
  }

  exe_rel <- gp_python_executable_relpath(os_type)
  py <- file.path(env_dirs, exe_rel)
  py <- py[file.exists(py)]
  normalizePath(py, winslash = "/", mustWork = FALSE)
}

gp_local_virtualenv_python_candidates <- function(env_dir = ".venv-gp",
                                                  roots = NULL,
                                                  os_type = .Platform$OS.type) {
  roots <- roots %||% c(getwd(), normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE))
  candidates <- character()
  for (root in roots) {
    if (is.null(root) || !nzchar(root)) {
      next
    }
    candidates <- c(
      candidates,
      file.path(root, env_dir, gp_python_executable_relpath(os_type))
    )
  }
  candidates <- candidates[file.exists(candidates)]
  normalizePath(unique(candidates), winslash = "/", mustWork = FALSE)
}

gp_conda_python_candidates <- function(env_names) {
  fallback_paths <- character()
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    exe_rel <- file.path("python.exe")
    conda_roots <- c(
      file.path(home, "AppData", "Local", "r-miniconda", "envs"),
      file.path(home, "AppData", "Local", "miniconda3", "envs"),
      file.path(home, "miniconda3", "envs")
    )
  } else {
    exe_rel <- file.path("bin", "python")
    conda_roots <- c(
      file.path(home, "miniconda3", "envs"),
      file.path(home, ".conda", "envs")
    )
  }
  for (root in conda_roots) {
    for (env_name in env_names) {
      fallback_paths <- c(fallback_paths, file.path(root, env_name, exe_rel))
    }
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    existing <- fallback_paths[file.exists(fallback_paths)]
    return(normalizePath(existing, winslash = "/", mustWork = FALSE))
  }

  envs <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.null(envs) || !nrow(envs)) {
    existing <- fallback_paths[file.exists(fallback_paths)]
    return(normalizePath(existing, winslash = "/", mustWork = FALSE))
  }

  name_col <- intersect(c("name", "envname"), names(envs))
  py_col <- intersect(c("python", "python_path"), names(envs))
  if (!length(name_col) || !length(py_col)) {
    return(character())
  }

  matched <- envs[envs[[name_col[[1L]]]] %in% env_names, , drop = FALSE]
  if (!nrow(matched)) {
    return(character())
  }

  py <- matched[[py_col[[1L]]]]
  py <- py[nzchar(py)]
  if (!length(py)) {
    py <- character()
  }

  existing <- unique(c(py[file.exists(py)], fallback_paths[file.exists(fallback_paths)]))
  normalizePath(existing, winslash = "/", mustWork = FALSE)
}

gp_python_candidates <- function() {
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = FALSE)
  workon_home <- Sys.getenv("WORKON_HOME", unset = "")
  if (!nzchar(workon_home)) {
    workon_home <- file.path(home, ".virtualenvs")
  }
  env_names <- c("predictgp", "predictdl", "r-torchgwas")
  env_paths <- vapply(
    env_names,
    function(env_name) gp_virtualenv_python_path(workon_home, env_name),
    character(1)
  )
  extra_virtualenvs <- gp_virtualenv_python_candidates(workon_home)
  conda_paths <- gp_conda_python_candidates(env_names)
  system_candidates <- gp_system_python_candidates()
  if (identical(gp_platform_family(), "windows")) {
    system_candidates <- c(
      system_candidates,
      "C:/Python/python.exe",
      "C:/Program Files/Python311/python.exe",
      "C:/Program Files/Python312/python.exe"
    )
  }

  unique(c(
    Sys.getenv("PREDICTPRO_GP_PYTHON", unset = ""),
    Sys.getenv("PREDICTPRO_PYTHON", unset = ""),
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    env_paths,
    extra_virtualenvs,
    conda_paths,
    gp_local_virtualenv_python_candidates(".venv-gp"),
    system_candidates
  ))
}

gp_python_path_bonus <- function(py_bin, purpose = c("general", "ml", "dl", "gp")) {
  purpose <- match.arg(purpose)
  if (is.null(py_bin) || !nzchar(py_bin)) {
    return(0)
  }
  py_norm <- tolower(normalizePath(py_bin, winslash = "/", mustWork = FALSE))

  if (identical(purpose, "gp")) {
    if (grepl("/\\.venv-gp/", py_norm, fixed = FALSE) || grepl("/predictgp/", py_norm, fixed = FALSE)) {
      return(50)
    }
    if (grepl("/r-torchgwas/", py_norm, fixed = FALSE) || grepl("/predictdl/", py_norm, fixed = FALSE)) {
      return(-10)
    }
  }

  if (identical(purpose, "dl")) {
    if (grepl("/r-torchgwas/", py_norm, fixed = FALSE) || grepl("/predictdl/", py_norm, fixed = FALSE)) {
      return(30)
    }
  }

  if (identical(purpose, "ml")) {
    if (basename(py_norm) %in% c("python", "python3", "python.exe") ||
        grepl("/predictml/", py_norm, fixed = FALSE)) {
      return(10)
    }
  }

  0
}

gp_python_probe <- function(py_bin,
                            modules = c("numpy", "sklearn", "xgboost", "lightgbm", "catboost", "torch")) {
  if (is.null(py_bin) || !nzchar(py_bin) || !file.exists(py_bin)) {
    return(NULL)
  }
  # Each probe starts a Python process (and importing torch takes seconds),
  # so results are cached per interpreter and module set for the session.
  cache_key <- paste0(
    "python_probe::", normalizePath(py_bin, winslash = "/", mustWork = FALSE),
    "::", paste(sort(modules), collapse = ",")
  )
  if (exists(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE)) {
    return(get(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE))
  }

  py_modules <- paste(sprintf("'%s'", modules), collapse = ", ")
  torch_lines <- if ("torch" %in% modules) c(
    "try:",
    "    import torch",
    "    out['torch_version'] = getattr(torch, '__version__', None)",
    "    out['cuda_available'] = bool(torch.cuda.is_available())",
    "    out['cuda_device_count'] = int(torch.cuda.device_count()) if out['cuda_available'] else 0",
    "    out['cuda_device_name'] = torch.cuda.get_device_name(0) if out['cuda_available'] else None",
    "except Exception as e:",
    "    out['torch_error'] = str(e)"
  )
  py_code <- paste(
    c(
      "import importlib.util as u",
      "import json",
      paste0("mods = [", py_modules, "]"),
      "out = {m: bool(u.find_spec(m)) for m in mods}",
      torch_lines,
      "print(json.dumps(out))"
    ),
    collapse = "\n"
  )
  script <- tempfile(fileext = ".py")
  writeLines(py_code, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  out <- tryCatch(
    suppressWarnings(system2(
      py_bin,
      gp_quote_system_args(script),
      stdout = TRUE,
      stderr = FALSE
    )),
    error = function(e) character()
  )
  if (!length(out)) {
    return(NULL)
  }

  probe <- tryCatch(jsonlite::fromJSON(out[[1]]), error = function(e) NULL)
  if (!is.null(probe)) {
    assign(cache_key, probe, envir = PredictProR_runtime_cache)
  }
  probe
}

gp_python_module_score <- function(py_bin,
                                   purpose = c("general", "ml", "dl", "gp"),
                                   modules = c("numpy", "sklearn", "xgboost", "lightgbm", "catboost", "torch")) {
  purpose <- match.arg(purpose)
  probe <- gp_python_probe(py_bin, modules = modules)
  if (is.null(probe)) {
    return(0)
  }

  if (identical(purpose, "gp")) {
    weights <- c(numpy = 3, torch = 5)
    score <- sum(vapply(names(weights), function(m) if (isTRUE(probe[[m]])) weights[[m]] else 0, numeric(1)))
    gpytorch_ok <- tryCatch(gp_python_probe(py_bin, modules = c("gpytorch", "linear_operator", "zarr")), error = function(e) NULL)
    if (!is.null(gpytorch_ok)) {
      score <- score +
        if (isTRUE(gpytorch_ok[["gpytorch"]])) 8 else 0 +
        if (isTRUE(gpytorch_ok[["linear_operator"]])) 4 else 0 +
        if (isTRUE(gpytorch_ok[["zarr"]])) 3 else 0
    }
    if (isTRUE(probe[["cuda_available"]])) {
      score <- score + 10 + as.numeric(probe[["cuda_device_count"]] %||% 0)
    }
    return(score)
  }

  if (identical(purpose, "dl")) {
    weights <- c(numpy = 2, sklearn = 0, xgboost = 0, lightgbm = 0, catboost = 0, torch = 8)
    score <- sum(vapply(names(weights), function(m) if (isTRUE(probe[[m]])) weights[[m]] else 0, numeric(1)))
    if (isTRUE(probe[["cuda_available"]])) {
      score <- score + 20 + as.numeric(probe[["cuda_device_count"]] %||% 0)
    }
    return(score)
  }

  weights <- c(numpy = 3, sklearn = 3, xgboost = 2, lightgbm = 2, catboost = 2, torch = 2)
  sum(vapply(names(weights), function(m) if (isTRUE(probe[[m]])) weights[[m]] else 0, numeric(1)))
}

PredictProR_runtime_cache <- new.env(parent = emptyenv())

gp_preferred_python <- function(purpose = c("general", "ml", "dl", "gp")) {
  purpose <- match.arg(purpose)
  cache_key <- paste0("preferred_python::", purpose)
  if (exists(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE)) {
    cached <- get(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE)
    if (!is.null(cached) && nzchar(cached) && file.exists(cached)) {
      return(cached)
    }
  }

  candidates <- gp_python_candidates()
  candidates <- candidates[nzchar(candidates)]
  if (!length(candidates)) {
    return(NULL)
  }

  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    return(NULL)
  }
  existing <- unique(normalizePath(existing, winslash = "/", mustWork = TRUE))
  scores <- vapply(existing, gp_python_module_score, numeric(1), purpose = purpose) +
    vapply(existing, gp_python_path_bonus, numeric(1), purpose = purpose)
  best_idx <- order(scores, decreasing = TRUE, na.last = TRUE)
  best <- existing[[best_idx[[1]]]]
  assign(cache_key, best, envir = PredictProR_runtime_cache)
  best
}

gp_configure_python_runtime <- function(force = FALSE, purpose = c("general", "ml", "dl")) {
  purpose <- match.arg(purpose)
  preferred <- gp_preferred_python(purpose = purpose)
  if (is.null(preferred)) {
    return(invisible(NULL))
  }

  current <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (!force && nzchar(current)) {
    return(invisible(normalizePath(current, winslash = "/", mustWork = FALSE)))
  }

  Sys.setenv(RETICULATE_PYTHON = preferred)
  invisible(preferred)
}

gp_configure_torch_runtime_env <- function() {
  if (!nzchar(Sys.getenv("CUBLAS_WORKSPACE_CONFIG", unset = ""))) {
    Sys.setenv(CUBLAS_WORKSPACE_CONFIG = ":4096:8")
  }
  invisible(Sys.getenv("CUBLAS_WORKSPACE_CONFIG", unset = ""))
}

gp_python_runtime_status <- function(initialize = FALSE,
                                     python_path = NULL,
                                     preferred_python = NULL,
                                     purpose = c("general", "ml", "dl", "gp"),
                                     include_accelerator = TRUE) {
  purpose <- match.arg(purpose)
  reticulate_available <- requireNamespace("reticulate", quietly = TRUE)
  preferred_python <- preferred_python %||% tryCatch(gp_preferred_python(purpose = purpose), error = function(e) NULL)
  purpose_python <- switch(
    purpose,
    dl = {
      py <- Sys.getenv("PREDICTPRO_DL_PYTHON", unset = "")
      if (!nzchar(py)) {
        py <- Sys.getenv("PREDICTPRO_PYTHON_DL", unset = "")
      }
      py
    },
    gp = Sys.getenv("PREDICTPRO_GP_PYTHON", unset = ""),
    ml = Sys.getenv("PREDICTPRO_ML_PYTHON", unset = ""),
    general = ""
  )

  out <- list(
    preferred_python = preferred_python %||% NULL,
    active_python = NULL,
    reticulate_available = reticulate_available,
    reticulate_initialized = FALSE,
    python_configured = FALSE,
    python_exists = FALSE,
    python_initializable = FALSE,
    torch_available = FALSE,
    cuda_available = FALSE,
    num_gpus = 0L,
    gpu_name = NULL
  )

  if (!reticulate_available) {
    return(out)
  }

  out$reticulate_initialized <- tryCatch(
    reticulate::py_available(initialize = FALSE),
    error = function(e) FALSE
  )

  active_python <- python_path %||% purpose_python
  if (!nzchar(active_python)) {
    active_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  }
  if (!nzchar(active_python)) {
    active_python <- tryCatch(gp_detect_python(), error = function(e) NULL) %||% ""
  }
  if (!nzchar(active_python)) {
    active_python <- preferred_python %||% ""
  }

  if (!nzchar(active_python)) {
    return(out)
  }

  out$active_python <- active_python
  out$python_configured <- TRUE
  out$python_exists <- file.exists(active_python)
  out$python_initializable <- out$reticulate_initialized || isTRUE(out$python_exists)

  if (!initialize && !out$reticulate_initialized) {
    if (isTRUE(include_accelerator)) {
      probe <- tryCatch(
        gp_python_probe(active_python, modules = "torch"),
        error = function(e) NULL
      )
      if (!is.null(probe)) {
        out$torch_available <- isTRUE(probe[["torch"]])
        out$cuda_available <- isTRUE(probe[["cuda_available"]])
        out$num_gpus <- suppressWarnings(as.integer(
          probe[["cuda_device_count"]] %||% 0L
        ))
        out$gpu_name <- probe[["cuda_device_name"]] %||% NULL
      }
    }
    return(out)
  }

  if (!out$reticulate_initialized) {
    init_ok <- tryCatch({
      gp_init_python_once(active_python)
      TRUE
    }, error = function(e) FALSE)
    out$reticulate_initialized <- tryCatch(
      reticulate::py_available(initialize = FALSE),
      error = function(e) FALSE
    )
    out$python_initializable <- isTRUE(init_ok) && out$reticulate_initialized
    if (!out$reticulate_initialized) {
      return(out)
    }
  }

  out$python_initializable <- TRUE
  out$active_python <- tryCatch({
    conf <- reticulate::py_config()
    conf$python %||% active_python
  }, error = function(e) active_python)

  if (!isTRUE(include_accelerator)) {
    return(out)
  }

  out$torch_available <- tryCatch(
    reticulate::py_module_available("torch"),
    error = function(e) FALSE
  )
  if (!out$torch_available) {
    return(out)
  }

  torch <- tryCatch(reticulate::import("torch", delay_load = FALSE), error = function(e) NULL)
  if (is.null(torch)) {
    return(out)
  }

  out$cuda_available <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  if (!out$cuda_available) {
    return(out)
  }

  out$num_gpus <- tryCatch(as.integer(torch$cuda$device_count()), error = function(e) 0L)
  out$gpu_name <- tryCatch({
    if (out$num_gpus > 0L) as.character(torch$cuda$get_device_name(as.integer(0L))) else NULL
  }, error = function(e) NULL)

  out
}

gp_runtime_metadata <- function(execution_policy = NULL,
                                python_path = NULL,
                                preferred_python = NULL,
                                purpose = c("general", "ml", "dl", "gp"),
                                include_accelerator = TRUE,
                                context = "model_execute",
                                extra_fields = NULL) {
  purpose <- match.arg(purpose)
  scalarize_metadata_value <- function(x) {
    if (length(x) == 0 || is.null(x)) {
      return(NA_character_)
    }
    if (is.function(x)) {
      nm <- tryCatch(deparse(substitute(x)), error = function(e) NULL)
      cls <- paste(class(x), collapse = ",")
      return(if (!is.null(nm) && nzchar(nm)) nm else cls)
    }
    if (is.list(x) && !is.data.frame(x)) {
      return(paste(vapply(x, scalarize_metadata_value, character(1)), collapse = ","))
    }
    val <- tryCatch(x[[1]], error = function(e) x)
    if (is.function(val)) {
      return(paste(class(val), collapse = ","))
    }
    as.character(val)
  }

  runtime <- gp_python_runtime_status(
    initialize = FALSE,
    python_path = python_path,
    preferred_python = preferred_python,
    purpose = purpose,
    include_accelerator = include_accelerator
  )

  kv <- list(
    context = context,
    preferred_python = runtime$preferred_python %||% NA_character_,
    active_python = runtime$active_python %||% NA_character_,
    reticulate_available = runtime$reticulate_available,
    reticulate_initialized = runtime$reticulate_initialized,
    python_configured = runtime$python_configured,
    python_exists = runtime$python_exists,
    python_initializable = runtime$python_initializable,
    torch_available = runtime$torch_available,
    cuda_available = runtime$cuda_available,
    runtime_num_gpus = runtime$num_gpus %||% 0L,
    runtime_gpu_name = runtime$gpu_name %||% NA_character_
  )

  if (!is.null(execution_policy)) {
    kv$policy_backend <- execution_policy$backend %||% NA_character_
    kv$policy_workers <- execution_policy$workers %||% NA_integer_
    kv$policy_plan <- execution_policy$plan %||% NA_character_
    kv$policy_chunk_size <- execution_policy$chunk_size %||% NA_integer_
    kv$policy_num_gpus <- execution_policy$num_gpus %||% NA_integer_
    kv$policy_decision_reason <- execution_policy$decision_reason %||% NA_character_
    kv$policy_backend_score <- execution_policy$backend_score %||% NA_real_
    kv$policy_hidden_threads <- tryCatch(detect_hidden_threads(), error = function(e) NA_integer_)
    kv$policy_blas_threads_env <- Sys.getenv("GP_BLAS_THREADS", unset = "")
    kv$policy_user_sequential_models <- paste(execution_policy$user_sequential_models %||% character(), collapse = ",")
    kv$policy_internal_sequential_count <- execution_policy$internal_sequential_count %||% NA_integer_
    kv$policy_parallel_count <- execution_policy$parallel_count %||% NA_integer_
    kv$policy_models <- paste(execution_policy$model_ids %||% character(), collapse = ",")
    kv$policy_worker_memory_gb <- execution_policy$worker_memory_gb %||% NA_real_
    kv$policy_worker_memory_source <- execution_policy$worker_memory_source %||% NA_character_
    kv$policy_memory_budget_gb <- execution_policy$memory_budget_gb %||% NA_real_
    kv$policy_memory_worker_cap <- execution_policy$memory_worker_cap %||% NA_integer_
  }

  if (is.list(extra_fields) && length(extra_fields) > 0L) {
    kv[names(extra_fields)] <- extra_fields
  }

  data.frame(
    key = names(kv),
    value = vapply(kv, scalarize_metadata_value, character(1)),
    stringsAsFactors = FALSE
  )
}

gp_runtime_profile_new <- function(context = "model_execute") {
  now <- Sys.time()
  list(
    context = context,
    started_at = now,
    last_at = now,
    events = list()
  )
}

gp_runtime_profile_mark <- function(profile,
                                    stage,
                                    detail = NULL) {
  if (is.null(profile)) {
    profile <- gp_runtime_profile_new()
  }
  now <- Sys.time()
  started_at <- profile$started_at %||% now
  last_at <- profile$last_at %||% started_at
  row <- data.frame(
    stage = as.character(stage[[1]]),
    elapsed_seconds = as.numeric(difftime(now, last_at, units = "secs")),
    cumulative_seconds = as.numeric(difftime(now, started_at, units = "secs")),
    detail = if (is.null(detail)) NA_character_ else as.character(detail[[1]]),
    stringsAsFactors = FALSE
  )
  profile$events[[length(profile$events) + 1L]] <- row
  profile$last_at <- now
  profile
}

gp_runtime_profile_table <- function(profile) {
  rows <- profile$events %||% list()
  if (!length(rows)) {
    return(data.frame(
      stage = character(),
      elapsed_seconds = numeric(),
      cumulative_seconds = numeric(),
      detail = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gp_runtime_profile_metadata_fields <- function(profile) {
  tbl <- gp_runtime_profile_table(profile)
  if (!nrow(tbl)) {
    return(list())
  }
  fields <- list(
    profile_context = profile$context %||% NA_character_,
    profile_stage_count = nrow(tbl),
    profile_total_seconds = tbl$cumulative_seconds[[nrow(tbl)]]
  )
  for (i in seq_len(nrow(tbl))) {
    stage_key <- gsub("[^A-Za-z0-9]+", "_", tolower(tbl$stage[[i]]))
    fields[[paste0("profile_", stage_key, "_seconds")]] <- tbl$elapsed_seconds[[i]]
  }
  fields
}

gp_init_python_once <- function(py_bin) {
  if (!requireNamespace("reticulate", quietly = TRUE) || is.null(py_bin) || !nzchar(py_bin)) {
    return(invisible(FALSE))
  }

  gp_configure_torch_runtime_env()

  if (!reticulate::py_available(initialize = FALSE)) {
    current <- Sys.getenv("RETICULATE_PYTHON", unset = "")
    target <- normalizePath(py_bin, winslash = "/", mustWork = FALSE)
    if (!nzchar(current) || !identical(normalizePath(current, winslash = "/", mustWork = FALSE), target)) {
      Sys.setenv(RETICULATE_PYTHON = target)
    }
    reticulate::use_python(py_bin, required = TRUE)
    try(reticulate::py_config(), silent = TRUE)
    return(invisible(TRUE))
  }

  conf <- reticulate::py_config()
  same <- tryCatch(
    identical(normalizePath(conf$python, winslash = "/"), normalizePath(py_bin, winslash = "/")),
    error = function(e) FALSE
  )
  if (!same) {
    message(sprintf("[reticulate] Already initialized to '%s'; requested '%s'.", conf$python, py_bin))
  }

  invisible(FALSE)
}

gp_ml_random_state <- function(default = 123L) {
  seed <- getOption("random_state", default = default)
  seed <- suppressWarnings(as.integer(seed[[1]] %||% default))
  if (!is.finite(seed)) {
    seed <- as.integer(default)
  }
  seed
}

gp_ml_numeric_predictor_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  input_rownames <- rownames(x)
  if (is.data.frame(x)) {
    x <- as.data.frame(lapply(x, function(col) {
      if (is.factor(col) || is.ordered(col)) {
        col <- as.character(col)
      }
      suppressWarnings(as.numeric(col))
    }), check.names = FALSE, stringsAsFactors = FALSE)
  }
  x_mat <- as.matrix(x)
  suppressWarnings(storage.mode(x_mat) <- "double")
  if (!is.null(input_rownames) && length(input_rownames) == nrow(x_mat)) {
    rownames(x_mat) <- input_rownames
  }
  if (is.null(colnames(x_mat))) {
    colnames(x_mat) <- paste0("V", seq_len(ncol(x_mat)))
  }
  x_mat
}

gp_ml_predictor_fill_values <- function(x_mat) {
  if (is.null(x_mat) || !ncol(x_mat)) {
    return(numeric(0))
  }
  vals <- apply(x_mat, 2L, function(col) {
    finite <- col[is.finite(col)]
    if (length(finite)) {
      mean(finite)
    } else {
      0
    }
  })
  vals[!is.finite(vals)] <- 0
  vals
}

gp_ml_fill_nonfinite_predictors <- function(x_mat, fill_values) {
  if (is.null(x_mat) || !length(x_mat) || !ncol(x_mat)) {
    return(x_mat)
  }
  fill_values <- as.numeric(fill_values)
  if (length(fill_values) != ncol(x_mat)) {
    fill_values <- rep(0, ncol(x_mat))
  }
  for (j in seq_len(ncol(x_mat))) {
    bad <- !is.finite(x_mat[, j])
    if (any(bad)) {
      x_mat[bad, j] <- fill_values[[j]]
    }
  }
  x_mat
}

gp_ml_align_predictor_matrix <- function(x_mat, target_names = NULL) {
  if (is.null(target_names) || !length(target_names)) {
    return(x_mat)
  }
  target_names <- as.character(target_names)
  if (!is.null(colnames(x_mat)) &&
      !anyDuplicated(target_names) &&
      !anyDuplicated(colnames(x_mat))) {
    out <- matrix(
      NA_real_,
      nrow = nrow(x_mat),
      ncol = length(target_names),
      dimnames = list(rownames(x_mat), target_names)
    )
    common <- intersect(target_names, colnames(x_mat))
    if (length(common)) {
      out[, common] <- x_mat[, common, drop = FALSE]
      return(out)
    }
  }
  if (ncol(x_mat) == length(target_names)) {
    colnames(x_mat) <- target_names
    return(x_mat)
  }
  stop("Predictor matrix columns cannot be aligned to the training feature layout.", call. = FALSE)
}

gp_ml_preprocess_predictors <- function(x, scaling = FALSE, centering = TRUE) {
  if (is.null(x)) {
    return(list(data = NULL, scaler = NULL, removed_cols = integer(0)))
  }

  x_mat <- gp_ml_numeric_predictor_matrix(x)
  fill_values <- gp_ml_predictor_fill_values(x_mat)
  x_filled <- gp_ml_fill_nonfinite_predictors(x_mat, fill_values)
  center <- if (isTRUE(centering) && ncol(x_filled)) colMeans(x_filled) else rep(0, ncol(x_filled))
  scale <- if (isTRUE(scaling) && ncol(x_filled)) apply(x_filled, 2L, stats::sd) else rep(1, ncol(x_filled))
  scale[!is.finite(scale) | scale == 0] <- 1
  x_proc <- sweep(sweep(x_filled, 2L, center, "-"), 2L, scale, "/")
  colnames(x_proc) <- colnames(x_mat)
  rownames(x_proc) <- rownames(x_mat)

  scaler <- list(
    center = center,
    scale = scale,
    impute_values = fill_values,
    colnames = colnames(x_mat),
    centering = isTRUE(centering),
    scaling = isTRUE(scaling)
  )

  list(
    data = x_proc,
    scaler = scaler,
    removed_cols = integer(0),
    impute_values = fill_values
  )
}

gp_ml_apply_preprocessor <- function(x, prep) {
  if (is.null(x)) {
    return(NULL)
  }

  scaler <- prep$scaler %||% list()
  x_mat <- gp_ml_numeric_predictor_matrix(x)
  x_mat <- gp_ml_align_predictor_matrix(x_mat, scaler$colnames %||% colnames(x_mat))
  fill_values <- scaler$impute_values %||% gp_ml_predictor_fill_values(x_mat)
  x_filled <- gp_ml_fill_nonfinite_predictors(x_mat, fill_values)
  center <- scaler$center %||% rep(0, ncol(x_filled))
  scale <- scaler$scale %||% rep(1, ncol(x_filled))
  scale[!is.finite(scale) | scale == 0] <- 1
  x_proc <- sweep(sweep(x_filled, 2L, center, "-"), 2L, scale, "/")
  colnames(x_proc) <- scaler$colnames %||% colnames(x_mat)
  rownames(x_proc) <- rownames(x_mat)
  x_proc
}

gp_ml_bootstrap_target_metrics <- function(boot_matrix,
                                           train_test_label = NULL,
                                           observed_y = NULL) {
  pred_variances <- apply(boot_matrix, 2, var)
  pred_se <- apply(boot_matrix, 2, sd)
  pred_mean <- apply(boot_matrix, 2, mean)

  ref_idx <- rep(TRUE, length(pred_mean))
  if (!is.null(train_test_label) && length(train_test_label) == length(pred_mean)) {
    ref_idx <- train_test_label != "Test"
    if (!any(ref_idx, na.rm = TRUE)) {
      ref_idx <- rep(TRUE, length(pred_mean))
    }
  }

  reference_variance <- NA_real_
  if (!is.null(observed_y)) {
    observed_y <- as.numeric(observed_y)
    observed_y <- observed_y[is.finite(observed_y)]
    if (length(observed_y) > 1L) {
      reference_variance <- stats::var(observed_y, na.rm = TRUE)
    }
  }
  if (!is.finite(reference_variance) || reference_variance <= 0) {
    reference_variance <- stats::var(pred_mean[ref_idx], na.rm = TRUE)
  }
  if (!is.finite(reference_variance) || reference_variance <= 0) {
    reference_variance <- NA_real_
  }

  list(
    predicted_mean = pred_mean,
    standard_error = pred_se,
    prediction_error_var = pred_variances,
    reference_variance = reference_variance
  )
}

gp_ml_gaussian_stability <- function(prediction_error_var,
                                     reference_variance,
                                     high_reliability_thres = 0.9,
                                     low_reliability_thres = 0.5) {
  if (length(reference_variance) == 1L) {
    reference_variance <- rep(reference_variance, length(prediction_error_var))
  }

  normalized_uncertainty <- rep(Inf, length(prediction_error_var))
  valid_idx <- is.finite(prediction_error_var) & prediction_error_var >= 0 &
    is.finite(reference_variance) & reference_variance > 0
  normalized_uncertainty[valid_idx] <- prediction_error_var[valid_idx] / reference_variance[valid_idx]

  stability <- rep(0, length(prediction_error_var))
  stability[valid_idx] <- 1 - normalized_uncertainty[valid_idx]
  stability <- pmax(0, pmin(1, stability))

  remarks <- ifelse(
    stability >= high_reliability_thres,
    "Stable",
    ifelse(stability <= low_reliability_thres, "Unstable", "Moderately Stable")
  )

  proportion_high <- mean(stability >= high_reliability_thres)
  proportion_medium <- mean(stability > low_reliability_thres & stability < high_reliability_thres)
  proportion_low <- mean(stability <= low_reliability_thres)

  list(
    stability = stability,
    remarks = remarks,
    normalized_uncertainty = normalized_uncertainty,
    proportion_high_reliability = proportion_high,
    proportion_medium_reliability = proportion_medium,
    proportion_low_reliability = proportion_low,
    stability_percentage = (proportion_high + proportion_medium) * 100
  )
}

gp_ml_gaussian_predictive_confidence_scale <- function(prediction_error_var,
                                                        reference_variance,
                                                        high_confidence_thres = 0.9,
                                                        low_confidence_thres = 0.5) {
  if (length(reference_variance) == 1L) {
    reference_variance <- rep(reference_variance, length(prediction_error_var))
  }
  normalized_uncertainty <- rep(Inf, length(prediction_error_var))
  valid_idx <- is.finite(prediction_error_var) & prediction_error_var >= 0 &
    is.finite(reference_variance) & reference_variance > 0
  normalized_uncertainty[valid_idx] <-
    prediction_error_var[valid_idx] / reference_variance[valid_idx]
  confidence <- rep(0, length(prediction_error_var))
  confidence[valid_idx] <- 1 / (1 + normalized_uncertainty[valid_idx])
  confidence <- pmax(0, pmin(1, confidence))
  remarks <- ifelse(
    confidence >= high_confidence_thres,
    "Stable",
    ifelse(confidence <= low_confidence_thres, "Unstable", "Moderately Stable")
  )
  proportion_high <- mean(confidence >= high_confidence_thres)
  proportion_medium <- mean(
    confidence > low_confidence_thres & confidence < high_confidence_thres
  )
  proportion_low <- mean(confidence <= low_confidence_thres)
  list(
    stability = confidence,
    remarks = remarks,
    normalized_uncertainty = normalized_uncertainty,
    proportion_high_reliability = proportion_high,
    proportion_medium_reliability = proportion_medium,
    proportion_low_reliability = proportion_low,
    stability_percentage = (proportion_high + proportion_medium) * 100
  )
}

gp_ml_predictive_reliability <- function(prediction_error_var,
                                         reference_variance,
                                         high_reliability_thres = 0.9,
                                         low_reliability_thres = 0.5) {
  prediction_error_var <- suppressWarnings(as.numeric(prediction_error_var))
  n <- length(prediction_error_var)
  reference_variance <- suppressWarnings(as.numeric(reference_variance))
  if (length(reference_variance) == 1L) {
    reference_variance <- rep(reference_variance, n)
  }
  if (length(reference_variance) != n) {
    reference_variance <- rep(NA_real_, n)
  }

  valid <- is.finite(prediction_error_var) & prediction_error_var >= 0 &
    is.finite(reference_variance) & reference_variance > 0
  reliability <- rep(NA_real_, n)
  reliability[valid] <- reference_variance[valid] /
    (reference_variance[valid] + prediction_error_var[valid])

  remarks <- rep(NA_character_, n)
  remarks[valid] <- ifelse(
    reliability[valid] >= high_reliability_thres,
    "High Predictive Precision",
    ifelse(
      reliability[valid] <= low_reliability_thres,
      "Low Predictive Precision",
      "Moderate Predictive Precision"
    )
  )

  available_basis <- paste(
    "ML/DL predictive precision index = V_y/(V_y + PEV_i);",
    "PEV_i is target-specific held-out predictive MSE and V_y is the",
    "training phenotype reference variance; not genetic reliability"
  )
  unavailable_basis <- paste(
    "Reliability unavailable: finite target-specific predictive PEV is",
    "required; Prediction_stability is a separate bootstrap descriptive",
    "index; not genetic reliability and not identifiable as such for ML/DL"
  )

  list(
    reliability = reliability,
    remarks = remarks,
    variance_input = ifelse(valid, prediction_error_var, NA_real_),
    reference_variance = ifelse(valid, reference_variance, NA_real_),
    basis = ifelse(valid, available_basis, unavailable_basis)
  )
}

gp_ml_reliability_basis <- function(variance_source,
                                    reference_source = "all_final_prediction_variance") {
  variance_source <- as.character(variance_source)
  reference_source <- rep_len(as.character(reference_source), length(variance_source))
  detail <- rep(
    paste0(
      "U_i is a legacy descriptive variance fallback used only for plotting"
    ),
    length(variance_source)
  )
  in_sample <- grepl(
    "in_sample|residual_dispersion|observed_rows",
    variance_source,
    ignore.case = TRUE
  )
  detail[in_sample] <- paste0(
    "U_i is the legacy uncalibrated in-sample residual-dispersion proxy used only for plotting"
  )
  resampling <- grepl("bootstrap|fold_prediction|repeated_prediction|resampling", variance_source, ignore.case = TRUE)
  detail[resampling] <- paste0(
    "U_i is prediction-resampling SE_i^2 used only for marker-adjustment plotting"
  )
  across_environment <- grepl("across_environment", variance_source, ignore.case = TRUE)
  detail[across_environment] <- paste0(
    "U_i is the sum of row-level stability variance inputs divided by the squared ",
    "number of environments (zero cross-environment resampling-error covariance approximation)"
  )
  ref_detail <- ifelse(
    grepl("across_environment", reference_source, ignore.case = TRUE),
    "V_yhat is the across-environment aggregate fitted-prediction variance reference",
    ifelse(
      grepl("all_final|combined_train_test|yhat", reference_source, ignore.case = TRUE),
      "V_yhat is var(yhat) across all final prediction rows (Train and Test when both are present)",
      paste0("V_yhat uses fallback source ", reference_source)
    )
  )
  paste0(
    "ML/DL marker-adjustment reliability surrogate = max(0, min(1, 1 - U_i/V_yhat)); ",
    detail, "; ", ref_detail,
    "; U_i_source=", variance_source,
    "; compatibility alias of Prediction_stability; descriptive, uncalibrated, not genetic reliability"
  )
}

gp_ml_resolve_stability_input <- function(pred_variances,
                                          predicted_value,
                                          observed_aligned,
                                          reference_variance,
                                          train_test_label = NULL) {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  n <- length(predicted_value)
  observed_aligned <- suppressWarnings(as.numeric(observed_aligned))
  if (length(observed_aligned) != n) {
    observed_aligned <- rep(NA_real_, n)
  }
  reference_variance <- suppressWarnings(as.numeric(reference_variance))
  if (length(reference_variance) == 1L) {
    reference_variance <- rep(reference_variance, n)
  }
  if (length(reference_variance) != n) {
    reference_variance <- rep(NA_real_, n)
  }

  fitted_var <- stats::var(predicted_value[is.finite(predicted_value)], na.rm = TRUE)
  reference_source <- "all_final_prediction_variance"
  if (!is.finite(fitted_var) || fitted_var <= 0) {
    observed_var <- stats::var(observed_aligned[is.finite(observed_aligned)], na.rm = TRUE)
    fitted_var <- if (is.finite(observed_var) && observed_var > 0) observed_var else NA_real_
    reference_source <- "training_phenotype_variance_fallback"
  }
  supplied_ref <- reference_variance[is.finite(reference_variance) & reference_variance > 0]
  if ((!is.finite(fitted_var) || fitted_var <= 0) && length(supplied_ref)) {
    fitted_var <- supplied_ref[[1L]]
    reference_source <- "supplied_reference_variance_fallback"
  }
  if (!is.finite(fitted_var) || fitted_var <= 0) {
    fitted_var <- 1
    reference_source <- "unit_variance_fallback"
  }
  reference_variance <- rep(fitted_var, n)

  raw_variance <- suppressWarnings(as.numeric(pred_variances))
  if (length(raw_variance) != n) {
    raw_variance <- rep(NA_real_, n)
  }
  variance_input <- raw_variance
  source <- rep("bootstrap_prediction_se_squared_marker_adjustment", n)

  residual <- observed_aligned - predicted_value
  residual_mse <- mean(residual[is.finite(residual)]^2, na.rm = TRUE)
  missing <- !is.finite(variance_input) | variance_input < 0
  if (is.finite(residual_mse) && residual_mse >= 0) {
    variance_input[missing] <- residual_mse
    source[missing] <- "in_sample_residual_mse_uncalibrated"
  }

  missing <- !is.finite(variance_input) | variance_input < 0
  variance_input[missing] <- 0.25 * reference_variance[missing]
  source[missing] <- "legacy_reference_variance_fallback_uncalibrated"

  scale_values <- c(predicted_value, observed_aligned)
  scale_values <- scale_values[is.finite(scale_values)]
  scale_ref <- if (length(scale_values)) stats::median(abs(scale_values), na.rm = TRUE) else 1
  if (!is.finite(scale_ref) || scale_ref <= 0) {
    scale_ref <- 1
  }
  min_var <- max((scale_ref * 1e-4)^2, .Machine$double.eps)
  variance_input <- pmax(variance_input, min_var)

  list(
    variance_input = variance_input,
    reference_variance = reference_variance,
    source = source,
    reference_source = rep(reference_source, n),
    basis = gp_ml_reliability_basis(source, reference_source)
  )
}

gp_ml_finite_sample_quantile <- function(x, probability) {
  x <- sort(suppressWarnings(as.numeric(x)))
  x <- x[is.finite(x)]
  probability <- suppressWarnings(as.numeric(probability))[1L]
  if (!length(x) || !is.finite(probability) || probability <= 0 || probability >= 1) {
    return(NA_real_)
  }

  # Split-conformal finite-sample correction: use the
  # ceil((n + 1) * probability)-th order statistic, capped at n.
  rank <- min(length(x), ceiling((length(x) + 1L) * probability))
  x[[max(1L, rank)]]
}

gp_ml_finite_sample_lower_quantile <- function(x, probability) {
  x <- sort(suppressWarnings(as.numeric(x)))
  x <- x[is.finite(x)]
  probability <- suppressWarnings(as.numeric(probability))[1L]
  if (!length(x) || !is.finite(probability) || probability <= 0 || probability >= 1) {
    return(NA_real_)
  }
  rank <- max(1L, floor((length(x) + 1L) * probability))
  x[[min(length(x), rank)]]
}

gp_ml_heldout_residual_summary <- function(calibration,
                                           confidence_level = 0.95) {
  if (is.null(calibration)) {
    return(NULL)
  }
  residuals <- calibration$heldout_residuals %||%
    calibration$prediction_residuals %||%
    numeric()
  residuals <- suppressWarnings(as.numeric(residuals))
  residuals <- residuals[is.finite(residuals)]
  if (length(residuals) < 2L) {
    return(NULL)
  }

  mse <- mean(residuals^2)
  interval_radius <- gp_ml_finite_sample_quantile(
    abs(residuals),
    probability = confidence_level
  )
  if (!is.finite(mse) || mse < 0 || !is.finite(interval_radius) || interval_radius < 0) {
    return(NULL)
  }

  list(
    residuals = residuals,
    mse = mse,
    rmse = sqrt(mse),
    interval_radius = interval_radius,
    calibration_n = length(residuals),
    calibration_coverage = mean(abs(residuals) <= interval_radius),
    source = calibration$calibration_source %||% "heldout_cross_fitted"
  )
}

gp_ml_local_cross_fitted_target_uncertainty <- function(
    calibration,
    final_prediction,
    confidence_level = 0.95,
    prior_neighbor_strength = 5) {
  required <- c(
    "heldout_predictors", "target_predictors", "heldout_ids", "target_ids",
    "heldout_residuals", "heldout_fold_id", "target_fold_predictions"
  )
  if (is.null(calibration) || !all(required %in% names(calibration))) {
    return(NULL)
  }

  heldout_x <- tryCatch(as.matrix(calibration$heldout_predictors), error = function(e) NULL)
  target_x <- tryCatch(as.matrix(calibration$target_predictors), error = function(e) NULL)
  target_predictions <- tryCatch(
    as.matrix(calibration$target_fold_predictions),
    error = function(e) NULL
  )
  residuals <- suppressWarnings(as.numeric(calibration$heldout_residuals))
  fold_id <- suppressWarnings(as.integer(calibration$heldout_fold_id))
  heldout_ids <- as.character(calibration$heldout_ids)
  target_ids <- as.character(calibration$target_ids)
  final_prediction <- suppressWarnings(as.numeric(final_prediction))
  if (is.null(heldout_x) || is.null(target_x) || is.null(target_predictions) ||
      nrow(heldout_x) != length(residuals) ||
      nrow(heldout_x) != length(heldout_ids) ||
      nrow(target_x) != length(target_ids) ||
      nrow(target_x) != length(final_prediction) ||
      ncol(heldout_x) != ncol(target_x) ||
      ncol(target_predictions) != nrow(target_x)) {
    return(NULL)
  }
  storage.mode(heldout_x) <- "numeric"
  storage.mode(target_x) <- "numeric"

  valid <- is.finite(residuals) & is.finite(fold_id) &
    fold_id >= 1L & fold_id <= nrow(target_predictions) &
    nzchar(heldout_ids)
  if (sum(valid) < 6L) {
    return(NULL)
  }
  heldout_x <- heldout_x[valid, , drop = FALSE]
  heldout_ids <- heldout_ids[valid]
  residuals <- residuals[valid]

  finite_column <- vapply(seq_len(ncol(heldout_x)), function(j) {
    values <- heldout_x[, j]
    sum(is.finite(values)) >= 2L && is.finite(stats::sd(values, na.rm = TRUE)) &&
      stats::sd(values, na.rm = TRUE) > sqrt(.Machine$double.eps)
  }, logical(1))
  if (!any(finite_column)) {
    return(NULL)
  }
  heldout_x <- heldout_x[, finite_column, drop = FALSE]
  target_x <- target_x[, finite_column, drop = FALSE]
  for (j in seq_len(ncol(heldout_x))) {
    center <- stats::median(heldout_x[, j], na.rm = TRUE)
    spread <- stats::sd(heldout_x[, j], na.rm = TRUE)
    heldout_x[!is.finite(heldout_x[, j]), j] <- center
    target_x[!is.finite(target_x[, j]), j] <- center
    heldout_x[, j] <- (heldout_x[, j] - center) / spread
    target_x[, j] <- (target_x[, j] - center) / spread
  }

  target_instability <- apply(target_predictions, 2L, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::var(x)
  })
  finite_instability <- target_instability[
    is.finite(target_instability) & target_instability >= 0
  ]
  instability_fallback <- if (length(finite_instability)) {
    stats::median(finite_instability)
  } else {
    0
  }
  target_instability[
    !is.finite(target_instability) | target_instability < 0
  ] <- instability_fallback

  heldout_target_match <- match(heldout_ids, target_ids)
  if (anyNA(heldout_target_match)) {
    return(NULL)
  }
  heldout_instability <- target_instability[heldout_target_match]
  noise_proxy <- pmax(residuals^2 - heldout_instability, 0)
  global_noise <- mean(noise_proxy, na.rm = TRUE)
  if (!is.finite(global_noise) || global_noise < 0) {
    return(NULL)
  }

  local_noise <- rep(NA_real_, nrow(target_x))
  neighbor_count <- max(5L, as.integer(ceiling(sqrt(nrow(heldout_x)))))
  # All target-to-held-out mean squared distances in one BLAS product,
  # (|t|^2 + |h|^2 - 2 t.h) / p. The former per-target sweep() built a full
  # held-out-sized temporary for every target (~340 s at 1.8k x 10k markers,
  # and the function runs twice per fit). Values within rounding of zero are
  # set to 0 so exact duplicates keep zero distance, as before.
  p_cols <- ncol(heldout_x)
  h2 <- rowSums(heldout_x^2)
  t2 <- rowSums(target_x^2)
  all_sq <- (outer(t2, h2, "+") - 2 * tcrossprod(target_x, heldout_x)) / p_cols
  all_sq[all_sq < 1e-10 * (outer(t2, h2, "+") / p_cols + 1)] <- 0
  for (j in seq_len(nrow(target_x))) {
    squared_distance <- all_sq[j, ]
    eligible <- is.finite(squared_distance) & is.finite(noise_proxy)
    eligible <- eligible & heldout_ids != target_ids[[j]]
    candidates <- which(eligible)
    if (!length(candidates)) {
      next
    }
    candidates <- candidates[order(squared_distance[candidates])]
    candidates <- head(candidates, min(neighbor_count, length(candidates)))
    distance <- sqrt(pmax(squared_distance[candidates], 0))
    bandwidth <- stats::median(distance[is.finite(distance) & distance > 0])
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      bandwidth <- max(distance, na.rm = TRUE)
    }
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      weights <- rep(1, length(candidates))
    } else {
      weights <- exp(-0.5 * (distance / bandwidth)^2)
    }
    local_estimate <- stats::weighted.mean(noise_proxy[candidates], weights)
    effective_n <- sum(weights)^2 / sum(weights^2)
    local_weight <- effective_n / (effective_n + prior_neighbor_strength)
    local_noise[[j]] <-
      local_weight * local_estimate + (1 - local_weight) * global_noise
  }
  local_noise[!is.finite(local_noise) | local_noise < 0] <- global_noise

  raw_pev <- local_noise + target_instability
  training_raw_pev <- raw_pev[heldout_target_match]
  heldout_mse <- mean(residuals^2)
  mean_training_raw_pev <- mean(training_raw_pev, na.rm = TRUE)
  if (!is.finite(heldout_mse) || heldout_mse < 0 ||
      !is.finite(mean_training_raw_pev) || mean_training_raw_pev <= 0) {
    return(NULL)
  }
  calibration_factor <- heldout_mse / mean_training_raw_pev
  prediction_mse <- pmax(raw_pev * calibration_factor, 0)

  interval_radius <- gp_ml_finite_sample_quantile(
    abs(residuals),
    probability = confidence_level
  )
  local_interval_scale <- sqrt(prediction_mse / heldout_mse)
  lower <- final_prediction - interval_radius * local_interval_scale
  upper <- final_prediction + interval_radius * local_interval_scale
  heldout_interval_scale <- local_interval_scale[heldout_target_match]
  empirical_coverage <- mean(
    abs(residuals) <= interval_radius * heldout_interval_scale,
    na.rm = TRUE
  )

  list(
    prediction_mse = prediction_mse,
    lower = lower,
    upper = upper,
    calibration_n = length(residuals),
    source = "local_cross_fitted_residual_risk_plus_fold_instability",
    error_estimand = paste(
      "locally weighted held-out squared residual risk after removing",
      "fold-instability, plus target-specific fold prediction variance;",
      "mean-calibrated to cross-fitted predictive MSE; not mixed-model PEV"
    ),
    interval_method = paste(
      "locally scaled cross-fitted absolute-residual interval;",
      "empirical marginal calibration"
    ),
    target_fold_instability = target_instability,
    target_local_residual_risk = local_noise,
    mean_calibration_factor = calibration_factor,
    empirical_coverage = empirical_coverage
  )
}

gp_ml_cross_fitted_target_uncertainty <- function(calibration,
                                                   final_prediction,
                                                   confidence_level = 0.95) {
  if (is.null(calibration) || is.null(calibration$target_fold_predictions) ||
      is.null(calibration$heldout_fold_id)) {
    return(NULL)
  }
  target_predictions <- tryCatch(
    as.matrix(calibration$target_fold_predictions),
    error = function(e) NULL
  )
  residuals <- suppressWarnings(as.numeric(calibration$heldout_residuals))
  fold_id <- suppressWarnings(as.integer(calibration$heldout_fold_id))
  final_prediction <- suppressWarnings(as.numeric(final_prediction))
  if (is.null(target_predictions) || ncol(target_predictions) != length(final_prediction) ||
      length(residuals) != length(fold_id)) {
    return(NULL)
  }
  local_result <- gp_ml_local_cross_fitted_target_uncertainty(
    calibration = calibration,
    final_prediction = final_prediction,
    confidence_level = confidence_level
  )
  if (!is.null(local_result)) {
    return(local_result)
  }
  valid <- is.finite(residuals) & is.finite(fold_id) &
    fold_id >= 1L & fold_id <= nrow(target_predictions)
  residuals <- residuals[valid]
  fold_id <- fold_id[valid]
  if (length(residuals) < 2L) {
    return(NULL)
  }

  alpha <- 1 - confidence_level
  n_target <- length(final_prediction)
  prediction_mse <- lower <- upper <- rep(NA_real_, n_target)
  for (j in seq_len(n_target)) {
    fold_prediction <- target_predictions[fold_id, j]
    ok <- is.finite(fold_prediction) & is.finite(residuals)
    if (sum(ok) < 2L || !is.finite(final_prediction[[j]])) next
    fold_prediction <- fold_prediction[ok]
    residual_j <- residuals[ok]

    # Cross-fitted pseudo-outcomes combine the fold prediction at the target
    # with signed held-out residuals. Their squared deviation from the final
    # prediction estimates marginal predictive MSE for this target.
    pseudo_outcome <- fold_prediction + residual_j
    prediction_mse[[j]] <- mean((pseudo_outcome - final_prediction[[j]])^2)

    absolute_residual <- abs(residual_j)
    lower[[j]] <- gp_ml_finite_sample_lower_quantile(
      fold_prediction - absolute_residual,
      probability = alpha
    )
    upper[[j]] <- gp_ml_finite_sample_quantile(
      fold_prediction + absolute_residual,
      probability = confidence_level
    )
  }

  list(
    prediction_mse = prediction_mse,
    lower = lower,
    upper = upper,
    calibration_n = length(residuals)
  )
}

gp_ml_gaussian_confidence <- function(prediction_error_var,
                                      reference_variance,
                                      rank_confidence = NULL,
                                      high_confidence_thres = 0.9,
                                      low_confidence_thres = 0.5) {
  stability_res <- gp_ml_gaussian_predictive_confidence_scale(
    prediction_error_var = prediction_error_var,
    reference_variance = reference_variance,
    high_confidence_thres = high_confidence_thres,
    low_confidence_thres = low_confidence_thres
  )

  base_confidence <- stability_res$stability
  confidence <- base_confidence
  if (!is.null(rank_confidence) && length(rank_confidence) == length(confidence)) {
    confidence <- pmin(confidence, pmax(0, pmin(1, rank_confidence)))
  }
  confidence_remarks <- ifelse(
    confidence >= high_confidence_thres,
    "High Confidence",
    ifelse(confidence <= low_confidence_thres, "Low Confidence", "Moderate Confidence")
  )

  risk_score <- 1 - confidence
  risk_remarks <- ifelse(
    risk_score >= (1 - low_confidence_thres),
    "High Risk",
    ifelse(risk_score <= (1 - high_confidence_thres), "Low Risk", "Moderate Risk")
  )

  list(
    base_confidence = base_confidence,
    rank_confidence = rank_confidence,
    confidence = confidence,
    confidence_remarks = confidence_remarks,
    risk_score = risk_score,
    risk_remarks = risk_remarks
  )
}

gp_ml_bootstrap_oob_summary <- function(boot_results,
                                        observed_y,
                                        n_train,
                                        lower_prob = 0.05,
                                        upper_prob = 0.95) {
  if (is.null(boot_results) || is.null(observed_y)) {
    return(NULL)
  }

  bt <- tryCatch(as.matrix(boot_results$t), error = function(e) NULL)
  inbag <- tryCatch(boot::boot.array(boot_results), error = function(e) NULL)
  if (is.null(bt) || is.null(inbag) || !is.matrix(bt) || !is.matrix(inbag)) {
    return(NULL)
  }

  observed_y <- suppressWarnings(as.numeric(observed_y))
  n_train <- min(as.integer(n_train %||% 0L), length(observed_y), ncol(bt), ncol(inbag))
  if (!is.finite(n_train) || n_train < 2L) {
    return(NULL)
  }

  oob_mean <- rep(NA_real_, n_train)
  oob_var <- rep(NA_real_, n_train)
  oob_lower <- rep(NA_real_, n_train)
  oob_upper <- rep(NA_real_, n_train)
  oob_count <- integer(n_train)

  for (i in seq_len(n_train)) {
    oob_idx <- which(inbag[, i] == 0)
    oob_count[[i]] <- length(oob_idx)
    if (!length(oob_idx)) next
    vals <- bt[oob_idx, i]
    vals <- vals[is.finite(vals)]
    if (!length(vals)) next
    oob_mean[[i]] <- mean(vals)
    oob_var[[i]] <- if (length(vals) >= 2L) stats::var(vals) else 0
    oob_lower[[i]] <- stats::quantile(vals, probs = lower_prob, na.rm = TRUE, names = FALSE)
    oob_upper[[i]] <- stats::quantile(vals, probs = upper_prob, na.rm = TRUE, names = FALSE)
  }

  valid_mean <- is.finite(observed_y) & is.finite(oob_mean)
  valid_var <- is.finite(oob_var)
  if (sum(valid_mean) < 2L || sum(valid_var) < 2L) {
    return(NULL)
  }

  rmse <- sqrt(mean((observed_y[valid_mean] - oob_mean[valid_mean])^2, na.rm = TRUE))
  mean_pev <- mean(oob_var[valid_var], na.rm = TRUE)
  spearman <- suppressWarnings(stats::cor(
    observed_y[valid_mean],
    oob_mean[valid_mean],
    method = "spearman",
    use = "complete.obs"
  ))

  list(
    oob_mean = oob_mean,
    oob_var = oob_var,
    oob_lower = oob_lower,
    oob_upper = oob_upper,
    oob_count = oob_count,
    calibration_rmse = rmse,
    calibration_mean_pev = mean_pev,
    rank_spearman = spearman
  )
}

gp_ml_local_rank_confidence <- function(predicted_value,
                                        pred_se,
                                        groups = NULL) {
  n <- length(predicted_value)
  rank_confidence <- rep(NA_real_, n)
  nearest_margin <- rep(NA_real_, n)

  if (is.null(groups) || length(groups) != n) {
    groups <- rep("all", n)
  }
  groups <- as.character(groups)

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(predicted_value) & is.finite(pred_se))
    if (!length(idx)) next
    if (length(idx) == 1L) {
      nearest_margin[idx] <- Inf
      rank_confidence[idx] <- 1
      next
    }

    ord <- order(predicted_value[idx], na.last = NA)
    idx_ord <- idx[ord]
    vals <- predicted_value[idx_ord]

    margins <- rep(Inf, length(idx_ord))
    if (length(idx_ord) > 1L) {
      for (j in seq_along(idx_ord)) {
        diffs <- abs(vals[j] - vals[-j])
        margins[j] <- min(diffs, na.rm = TRUE)
      }
    }
    se_vals <- pmax(pred_se[idx_ord], .Machine$double.eps)
    conf_vals <- stats::pnorm(margins / se_vals)
    conf_vals[!is.finite(conf_vals)] <- 1
    rank_confidence[idx_ord] <- pmax(0, pmin(1, conf_vals))
    nearest_margin[idx_ord] <- margins
  }

  list(
    rank_confidence = rank_confidence,
    nearest_rank_margin = nearest_margin
  )
}

gp_ml_rank_context_features <- function(predicted_value,
                                        pred_se,
                                        groups = NULL,
                                        top_fraction = 0.2) {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  pred_se <- suppressWarnings(as.numeric(pred_se))
  n <- length(predicted_value)
  if (is.null(groups) || length(groups) != n) {
    groups <- rep("all", n)
  }
  groups <- as.character(groups)

  local_density <- rep(NA_real_, n)
  neighbor_mean_se <- rep(NA_real_, n)
  pred_pct <- rep(NA_real_, n)
  top_boundary_distance <- rep(NA_real_, n)

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(predicted_value) & is.finite(pred_se))
    if (!length(idx)) next
    if (length(idx) == 1L) {
      local_density[idx] <- 0
      neighbor_mean_se[idx] <- pred_se[idx]
      pred_pct[idx] <- 1
      top_boundary_distance[idx] <- abs(1 - top_fraction)
      next
    }

    ord <- order(predicted_value[idx], decreasing = TRUE, na.last = NA)
    idx_ord <- idx[ord]
    vals <- predicted_value[idx_ord]
    se_vals <- pmax(pred_se[idx_ord], .Machine$double.eps)
    pred_pct[idx_ord] <- seq_along(idx_ord) / length(idx_ord)
    top_boundary_distance[idx_ord] <- abs(pred_pct[idx_ord] - top_fraction)

    grp_sd <- stats::sd(vals, na.rm = TRUE)
    if (!is.finite(grp_sd) || grp_sd <= 0) {
      grp_sd <- stats::median(se_vals, na.rm = TRUE)
    }
    if (!is.finite(grp_sd) || grp_sd <= 0) {
      grp_sd <- 1
    }

    for (j in seq_along(idx_ord)) {
      diffs <- abs(vals[j] - vals)
      diffs[j] <- NA_real_
      radius <- max(se_vals[j], grp_sd * 0.1, na.rm = TRUE)
      in_radius <- sum(diffs <= radius, na.rm = TRUE)
      local_density[idx_ord[j]] <- in_radius / max(1, length(idx_ord) - 1L)

      neighbor_order <- order(diffs, na.last = NA)
      neighbor_order <- neighbor_order[seq_len(min(3L, length(neighbor_order)))]
      nbr_idx <- idx_ord[neighbor_order]
      nbr_idx <- nbr_idx[is.finite(pred_se[nbr_idx])]
      if (length(nbr_idx)) {
        neighbor_mean_se[idx_ord[j]] <- mean(pred_se[nbr_idx], na.rm = TRUE)
      } else {
        neighbor_mean_se[idx_ord[j]] <- pred_se[idx_ord[j]]
      }
    }
  }

  list(
    local_density = local_density,
    neighbor_mean_se = neighbor_mean_se,
    pred_pct = pred_pct,
    top_boundary_distance = top_boundary_distance
  )
}

gp_ml_rank_risk_classes <- function(predicted_value, observed_y) {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  observed_y <- suppressWarnings(as.numeric(observed_y))
  valid <- is.finite(predicted_value) & is.finite(observed_y)
  if (sum(valid) < 6L) return(NULL)

  pred_rank <- rank(-predicted_value[valid], ties.method = "average")
  actual_rank <- rank(-observed_y[valid], ties.method = "average")
  rank_shift <- abs(pred_rank - actual_rank)
  truth_q <- stats::quantile(rank_shift, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  truth_class <- ifelse(
    rank_shift <= truth_q[1], "Low Risk",
    ifelse(rank_shift <= truth_q[2], "Moderate Risk", "High Risk")
  )
  serious_failure <- as.integer(rank_shift > truth_q[2])

  list(
    valid = valid,
    rank_shift = rank_shift,
    truth_class = truth_class,
    serious_failure = serious_failure,
    low_shift_cut = truth_q[1],
    high_shift_cut = truth_q[2]
  )
}

gp_ml_fit_rank_risk_calibrator <- function(proxy_risk,
                                           serious_failure,
                                           truth_class,
                                           source_label) {
  proxy_risk <- suppressWarnings(as.numeric(proxy_risk))
  serious_failure <- as.integer(serious_failure)
  ok <- is.finite(proxy_risk) & !is.na(serious_failure)
  if (sum(ok) < 10L || length(unique(serious_failure[ok])) < 2L) {
    return(NULL)
  }

  dat <- data.frame(
    serious_failure = serious_failure[ok],
    proxy_risk = pmax(1e-6, pmin(1 - 1e-6, proxy_risk[ok]))
  )
  fit <- tryCatch(
    stats::glm(serious_failure ~ proxy_risk, data = dat, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  prob <- as.numeric(stats::predict(fit, type = "response"))
  brier <- mean((prob - dat$serious_failure)^2, na.rm = TRUE)
  cls <- truth_class[ok]
  low_cut <- suppressWarnings(stats::quantile(prob[cls == "Low Risk"], probs = 0.75, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(prob[cls == "High Risk"], probs = 0.25, na.rm = TRUE, names = FALSE))

  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    overall_q <- stats::quantile(prob, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    low_cut <- overall_q[1]
    high_cut <- overall_q[2]
  }
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    center <- mean(prob, na.rm = TRUE)
    spread <- stats::sd(prob, na.rm = TRUE)
    if (!is.finite(center)) center <- 0.5
    if (!is.finite(spread) || spread <= 0) spread <- 0.05
    low_cut <- max(0, center - spread / 2)
    high_cut <- min(1, center + spread / 2)
  }
  if (low_cut >= high_cut) {
    low_cut <- max(0, high_cut - 0.05)
  }

  pred_class <- ifelse(
    prob <= low_cut, "Low Risk",
    ifelse(prob >= high_cut, "High Risk", "Moderate Risk")
  )

  list(
    fit = fit,
    risk_low_cut = as.numeric(low_cut),
    risk_high_cut = as.numeric(high_cut),
    calibration_score = mean(pred_class == cls, na.rm = TRUE),
    calibration_brier = brier,
    calibration_source = source_label
  )
}

gp_ml_apply_rank_risk_calibrator <- function(proxy_risk, calibrator) {
  if (is.null(calibrator) || is.null(calibrator$fit)) {
    return(pmax(0, pmin(1, proxy_risk)))
  }
  newdata <- data.frame(proxy_risk = pmax(1e-6, pmin(1 - 1e-6, as.numeric(proxy_risk))))
  pmax(0, pmin(1, as.numeric(stats::predict(calibrator$fit, newdata = newdata, type = "response"))))
}

gp_ml_pairwise_rank_risk_calibration <- function(predicted_value,
                                                 observed_y,
                                                 pred_se,
                                                 groups = NULL,
                                                 source_label = "pairwise") {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  observed_y <- suppressWarnings(as.numeric(observed_y))
  pred_se <- suppressWarnings(as.numeric(pred_se))
  n <- length(predicted_value)
  if (is.null(groups) || length(groups) != n) groups <- rep("all", n)
  groups <- as.character(groups)

  rank_res <- gp_ml_local_rank_confidence(predicted_value, pred_se, groups)
  ctx_res <- gp_ml_rank_context_features(predicted_value, pred_se, groups)
  pairwise_fail <- rep(NA_integer_, n)

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(predicted_value) & is.finite(observed_y) & is.finite(pred_se))
    if (length(idx) < 2L) next
    ord <- order(predicted_value[idx], na.last = NA)
    idx_ord <- idx[ord]
    vals <- predicted_value[idx_ord]
    for (j in seq_along(idx_ord)) {
      others <- setdiff(seq_along(idx_ord), j)
      if (!length(others)) next
      nbr_local <- others[which.min(abs(vals[j] - vals[others]))]
      nbr_idx <- idx_ord[nbr_local]
      pred_sign <- sign(predicted_value[idx_ord[j]] - predicted_value[nbr_idx])
      obs_sign <- sign(observed_y[idx_ord[j]] - observed_y[nbr_idx])
      pairwise_fail[idx_ord[j]] <- as.integer(pred_sign != 0 && obs_sign != 0 && pred_sign != obs_sign)
    }
  }

  classes <- gp_ml_rank_risk_classes(predicted_value, observed_y)
  if (is.null(classes)) return(NULL)
  valid <- classes$valid & is.finite(pred_se) & !is.na(pairwise_fail) &
    is.finite(ctx_res$pred_pct) & is.finite(ctx_res$local_density) &
    is.finite(ctx_res$neighbor_mean_se) & is.finite(ctx_res$top_boundary_distance)
  if (sum(valid) < 12L) return(NULL)

  target <- as.integer(pairwise_fail[valid] > 0 | classes$serious_failure[seq_len(sum(classes$valid))] > 0)
  if (length(unique(target)) < 2L) return(NULL)

  dat <- data.frame(
    target = target,
    proxy_risk = pmax(1e-6, pmin(1 - 1e-6, 1 - rank_res$rank_confidence[valid])),
    log_se = log1p(pmax(pred_se[valid], 0)),
    log_margin = log1p(pmax(rank_res$nearest_rank_margin[valid], 0)),
    pred_pct = ctx_res$pred_pct[valid],
    log_density = log1p(pmax(ctx_res$local_density[valid], 0)),
    log_neighbor_se = log1p(pmax(ctx_res$neighbor_mean_se[valid], 0)),
    log_top_boundary = log1p(pmax(ctx_res$top_boundary_distance[valid], 0))
  )
  fit <- tryCatch(
    stats::glm(
      target ~ proxy_risk + log_se + log_margin + pred_pct +
        log_density + log_neighbor_se + log_top_boundary,
      data = dat,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  prob <- as.numeric(stats::predict(fit, type = "response"))
  truth_class <- classes$truth_class[seq_len(sum(classes$valid))][valid[classes$valid]]
  low_cut <- suppressWarnings(stats::quantile(prob[truth_class == "Low Risk"], probs = 0.75, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(prob[truth_class == "High Risk"], probs = 0.25, na.rm = TRUE, names = FALSE))
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    overall_q <- stats::quantile(prob, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    low_cut <- overall_q[1]
    high_cut <- overall_q[2]
  }
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    low_cut <- 0.33
    high_cut <- 0.66
  }
  pred_class <- ifelse(prob <= low_cut, "Low Risk", ifelse(prob >= high_cut, "High Risk", "Moderate Risk"))

  list(
    fit = fit,
    risk_low_cut = as.numeric(low_cut),
    risk_high_cut = as.numeric(high_cut),
    calibration_score = mean(pred_class == truth_class, na.rm = TRUE),
    calibration_brier = mean((prob - target)^2, na.rm = TRUE),
    calibration_source = source_label
  )
}

gp_ml_apply_pairwise_rank_risk_calibrator <- function(predicted_value,
                                                      pred_se,
                                                      groups = NULL,
                                                      calibrator) {
  if (is.null(calibrator) || is.null(calibrator$fit)) {
    return(NULL)
  }
  n <- length(predicted_value)
  if (is.null(groups) || length(groups) != n) groups <- rep("all", n)
  groups <- as.character(groups)
  rank_res <- gp_ml_local_rank_confidence(predicted_value, pred_se, groups)
  ctx_res <- gp_ml_rank_context_features(predicted_value, pred_se, groups)
  newdata <- data.frame(
    proxy_risk = pmax(1e-6, pmin(1 - 1e-6, 1 - rank_res$rank_confidence)),
    log_se = log1p(pmax(pred_se, 0)),
    log_margin = log1p(pmax(rank_res$nearest_rank_margin, 0)),
    pred_pct = ctx_res$pred_pct,
    log_density = log1p(pmax(ctx_res$local_density, 0)),
    log_neighbor_se = log1p(pmax(ctx_res$neighbor_mean_se, 0)),
    log_top_boundary = log1p(pmax(ctx_res$top_boundary_distance, 0))
  )
  out <- suppressWarnings(as.numeric(stats::predict(calibrator$fit, newdata = newdata, type = "response")))
  pmax(0, pmin(1, out))
}

gp_ml_top_selection_risk_calibration <- function(predicted_value,
                                                 observed_y,
                                                 pred_se,
                                                 groups = NULL,
                                                 top_fraction = 0.2,
                                                 rank_shift_fraction = 0.15,
                                                 source_label = "top_selection") {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  observed_y <- suppressWarnings(as.numeric(observed_y))
  pred_se <- suppressWarnings(as.numeric(pred_se))
  n <- length(predicted_value)
  if (is.null(groups) || length(groups) != n) groups <- rep("all", n)
  groups <- as.character(groups)

  rank_res <- gp_ml_local_rank_confidence(predicted_value, pred_se, groups)
  pred_pct <- rep(NA_real_, n)
  actual_pct <- rep(NA_real_, n)
  target <- rep(NA_integer_, n)

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(predicted_value) & is.finite(observed_y) & is.finite(pred_se))
    if (length(idx) < 6L) next
    pred_rank <- rank(-predicted_value[idx], ties.method = "average")
    actual_rank <- rank(-observed_y[idx], ties.method = "average")
    pred_pct[idx] <- pred_rank / length(idx)
    actual_pct[idx] <- actual_rank / length(idx)
    top_n <- max(1L, floor(length(idx) * top_fraction))
    shift_cut <- max(2L, floor(length(idx) * rank_shift_fraction))
    pred_top <- pred_rank <= top_n
    actual_top <- actual_rank <= top_n
    target[idx] <- as.integer((pred_top != actual_top) | (abs(pred_rank - actual_rank) > shift_cut))
  }

  valid <- is.finite(predicted_value) & is.finite(observed_y) & is.finite(pred_se) &
    is.finite(rank_res$nearest_rank_margin) & is.finite(pred_pct) & !is.na(target)
  if (sum(valid) < 12L || length(unique(target[valid])) < 2L) return(NULL)

  dat <- data.frame(
    target = target[valid],
    proxy_risk = pmax(1e-6, pmin(1 - 1e-6, 1 - rank_res$rank_confidence[valid])),
    log_se = log1p(pmax(pred_se[valid], 0)),
    log_margin = log1p(pmax(rank_res$nearest_rank_margin[valid], 0)),
    pred_pct = pred_pct[valid]
  )
  fit <- tryCatch(
    stats::glm(target ~ proxy_risk + log_se + log_margin + pred_pct, data = dat, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  prob <- as.numeric(stats::predict(fit, type = "response"))
  low_cut <- suppressWarnings(stats::quantile(prob[dat$target == 0], probs = 0.75, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(prob[dat$target == 1], probs = 0.25, na.rm = TRUE, names = FALSE))
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    overall_q <- stats::quantile(prob, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE)
    low_cut <- overall_q[1]
    high_cut <- overall_q[2]
  }
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    low_cut <- 0.33
    high_cut <- 0.66
  }
  pred_class <- ifelse(prob <= low_cut, "Low Risk", ifelse(prob >= high_cut, "High Risk", "Moderate Risk"))
  truth_class <- ifelse(dat$target == 0, "Low Risk", "High Risk")

  list(
    fit = fit,
    risk_low_cut = as.numeric(low_cut),
    risk_high_cut = as.numeric(high_cut),
    calibration_score = mean((pred_class == truth_class) | (pred_class == "Moderate Risk"), na.rm = TRUE),
    calibration_brier = mean((prob - dat$target)^2, na.rm = TRUE),
    calibration_source = source_label,
    top_fraction = top_fraction,
    rank_shift_fraction = rank_shift_fraction
  )
}

gp_ml_apply_top_selection_risk_calibrator <- function(predicted_value,
                                                      pred_se,
                                                      groups = NULL,
                                                      calibrator) {
  if (is.null(calibrator) || is.null(calibrator$fit)) return(NULL)
  n <- length(predicted_value)
  if (is.null(groups) || length(groups) != n) groups <- rep("all", n)
  groups <- as.character(groups)
  rank_res <- gp_ml_local_rank_confidence(predicted_value, pred_se, groups)
  pred_pct <- rep(NA_real_, n)
  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(predicted_value) & is.finite(pred_se))
    if (!length(idx)) next
    pred_rank <- rank(-predicted_value[idx], ties.method = "average")
    pred_pct[idx] <- pred_rank / length(idx)
  }
  newdata <- data.frame(
    proxy_risk = pmax(1e-6, pmin(1 - 1e-6, 1 - rank_res$rank_confidence)),
    log_se = log1p(pmax(pred_se, 0)),
    log_margin = log1p(pmax(rank_res$nearest_rank_margin, 0)),
    pred_pct = pred_pct
  )
  out <- suppressWarnings(as.numeric(stats::predict(calibrator$fit, newdata = newdata, type = "response")))
  pmax(0, pmin(1, out))
}

gp_ml_cv_rank_risk_calibration <- function(predicted_cv,
                                           observed_y,
                                           groups = NULL,
                                           fold_id = NULL,
                                           target_fold_predictions = NULL,
                                           heldout_predictors = NULL,
                                           target_predictors = NULL,
                                           heldout_ids = NULL,
                                           target_ids = NULL) {
  predicted_cv <- suppressWarnings(as.numeric(predicted_cv))
  observed_y <- suppressWarnings(as.numeric(observed_y))
  valid <- is.finite(predicted_cv) & is.finite(observed_y)
  if (sum(valid) < 6L) {
    return(NULL)
  }

  if (is.null(groups) || length(groups) != length(predicted_cv)) {
    groups <- rep("all", length(predicted_cv))
  }
  groups <- as.character(groups)
  if (is.null(fold_id) || length(fold_id) != length(predicted_cv)) {
    fold_id <- rep(NA_integer_, length(predicted_cv))
  }

  rmse_cv <- sqrt(mean((predicted_cv[valid] - observed_y[valid])^2, na.rm = TRUE))
  if (!is.finite(rmse_cv) || rmse_cv < 0) {
    return(NULL)
  }

  heldout_residuals <- observed_y[valid] - predicted_cv[valid]
  if (is.null(heldout_ids) || length(heldout_ids) != length(predicted_cv)) {
    heldout_ids <- rownames(heldout_predictors)
  }
  if (is.null(target_ids)) {
    target_ids <- rownames(target_predictors)
  }
  predictor_payload_valid <- !is.null(heldout_predictors) &&
    !is.null(target_predictors) &&
    nrow(as.matrix(heldout_predictors)) == length(predicted_cv) &&
    length(heldout_ids) == length(predicted_cv) &&
    nrow(as.matrix(target_predictors)) == length(target_ids)

  rank_res <- gp_ml_local_rank_confidence(
    predicted_value = predicted_cv[valid],
    pred_se = rep(rmse_cv, sum(valid)),
    groups = groups[valid]
  )
  classes <- gp_ml_rank_risk_classes(predicted_cv[valid], observed_y[valid])
  if (is.null(classes)) return(NULL)
  proxy_risk <- 1 - rank_res$rank_confidence
  truth_class <- classes$truth_class
  low_cut <- suppressWarnings(stats::quantile(proxy_risk[truth_class == "Low Risk"], probs = 0.75, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(proxy_risk[truth_class == "High Risk"], probs = 0.25, na.rm = TRUE, names = FALSE))

  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    overall_q <- stats::quantile(proxy_risk, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    low_cut <- overall_q[1]
    high_cut <- overall_q[2]
  }
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    center <- mean(proxy_risk, na.rm = TRUE)
    spread <- stats::sd(proxy_risk, na.rm = TRUE)
    if (!is.finite(center)) center <- 0.5
    if (!is.finite(spread) || spread <= 0) spread <- 0.05
    low_cut <- max(0, center - spread / 2)
    high_cut <- min(1, center + spread / 2)
  }
  if (low_cut >= high_cut) {
    low_cut <- max(0, high_cut - 0.05)
  }

  pred_class <- ifelse(
    proxy_risk <= low_cut, "Low Risk",
    ifelse(proxy_risk >= high_cut, "High Risk", "Moderate Risk")
  )

  list(
    risk_low_cut = as.numeric(low_cut),
    risk_high_cut = as.numeric(high_cut),
    cv_rmse = rmse_cv,
    mean_rank_shift = mean(classes$rank_shift, na.rm = TRUE),
    median_rank_shift = stats::median(classes$rank_shift, na.rm = TRUE),
    calibration_score = mean(pred_class == truth_class, na.rm = TRUE),
    calibration_source = "heldout_internal_cv",
    heldout_predictions = predicted_cv[valid],
    heldout_observed = observed_y[valid],
    heldout_residuals = heldout_residuals,
    heldout_fold_id = suppressWarnings(as.integer(fold_id[valid])),
    target_fold_predictions = target_fold_predictions,
    heldout_predictors = if (predictor_payload_valid) {
      as.matrix(heldout_predictors)[valid, , drop = FALSE]
    } else NULL,
    target_predictors = if (predictor_payload_valid) {
      as.matrix(target_predictors)
    } else NULL,
    heldout_ids = if (predictor_payload_valid) as.character(heldout_ids[valid]) else NULL,
    target_ids = if (predictor_payload_valid) as.character(target_ids) else NULL,
    heldout_mse = mean(heldout_residuals^2),
    heldout_n = length(heldout_residuals)
  )
}

gp_ml_rank_risk_calibration <- function(observed_y,
                                        oob_mean,
                                        oob_var,
                                        reference_variance) {
  observed_y <- suppressWarnings(as.numeric(observed_y))
  oob_mean <- suppressWarnings(as.numeric(oob_mean))
  oob_var <- suppressWarnings(as.numeric(oob_var))
  valid <- is.finite(observed_y) & is.finite(oob_mean) & is.finite(oob_var)
  if (sum(valid) < 6L || !is.finite(reference_variance) || reference_variance <= 0) {
    return(NULL)
  }

  rank_res <- gp_ml_local_rank_confidence(
    predicted_value = oob_mean[valid],
    pred_se = sqrt(pmax(oob_var[valid], 0)),
    groups = rep("Train", sum(valid))
  )
  base_conf <- exp(-oob_var[valid] / reference_variance)
  proxy_conf <- pmin(base_conf, rank_res$rank_confidence)
  proxy_risk <- 1 - proxy_conf
  classes <- gp_ml_rank_risk_classes(oob_mean[valid], observed_y[valid])
  if (is.null(classes)) return(NULL)
  truth_class <- classes$truth_class
  low_cut <- suppressWarnings(stats::quantile(proxy_risk[truth_class == "Low Risk"], probs = 0.75, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(proxy_risk[truth_class == "High Risk"], probs = 0.25, na.rm = TRUE, names = FALSE))

  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    overall_q <- stats::quantile(proxy_risk, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    low_cut <- overall_q[1]
    high_cut <- overall_q[2]
  }
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    center <- mean(proxy_risk, na.rm = TRUE)
    spread <- stats::sd(proxy_risk, na.rm = TRUE)
    if (!is.finite(center)) center <- 0.5
    if (!is.finite(spread) || spread <= 0) spread <- 0.05
    low_cut <- max(0, center - spread / 2)
    high_cut <- min(1, center + spread / 2)
  }
  if (low_cut >= high_cut) {
    low_cut <- max(0, high_cut - 0.05)
  }

  pred_class <- ifelse(
    proxy_risk <= low_cut, "Low Risk",
    ifelse(proxy_risk >= high_cut, "High Risk", "Moderate Risk")
  )

  list(
    risk_low_cut = as.numeric(low_cut),
    risk_high_cut = as.numeric(high_cut),
    mean_rank_shift = mean(classes$rank_shift, na.rm = TRUE),
    median_rank_shift = stats::median(classes$rank_shift, na.rm = TRUE),
    calibration_score = mean(pred_class == truth_class, na.rm = TRUE),
    calibration_source = "bootstrap_oob"
  )
}

gp_ml_align_observed_targets <- function(observed_y,
                                         train_test_label,
                                         target_length) {
  if (is.null(observed_y)) {
    return(rep(NA_real_, target_length))
  }

  observed_y <- suppressWarnings(as.numeric(observed_y))
  if (length(observed_y) == target_length) {
    return(observed_y)
  }

  aligned <- rep(NA_real_, target_length)
  if (!is.null(train_test_label) && length(train_test_label) == target_length) {
    ref_idx <- which(train_test_label != "Test")
    if (length(ref_idx) == length(observed_y)) {
      aligned[ref_idx] <- observed_y
      return(aligned)
    }
  }

  n_copy <- min(length(observed_y), target_length)
  aligned[seq_len(n_copy)] <- observed_y[seq_len(n_copy)]
  aligned
}

gp_ml_has_substantive_target_variation <- function(
    x,
    reference_rows = NULL,
    minimum_relative_range = 1e-2) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.null(reference_rows) && length(reference_rows) == length(x)) {
    x <- x[!is.na(reference_rows) & reference_rows]
  }
  x <- x[is.finite(x) & x >= 0]
  if (length(x) < 2L) {
    return(FALSE)
  }

  scale_value <- stats::median(abs(x))
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- mean(abs(x))
  }
  if (!is.finite(scale_value) || scale_value <= 0) {
    return(FALSE)
  }

  relative_range <- diff(range(x)) / scale_value
  is.finite(relative_range) && relative_range > minimum_relative_range
}

gp_ml_unavailable_target_uncertainty_estimand <- function() {
  paste(
    "target-specific predictive SE/PEV unavailable because estimated",
    "between-target variation is at or below the 1% minimum relative",
    "range; marginal cross-fitted MSE is not repeated as individual PEV"
  )
}

gp_ml_calibrate_gaussian_uncertainty <- function(predicted_value,
                                                 pred_se,
                                                 pred_variances,
                                                 lower_bound = NULL,
                                                 upper_bound = NULL,
                                                 observed_y = NULL,
                                                 train_test_label = NULL,
                                                 boot_results = NULL,
                                                 heldout_calibration = NULL,
                                                 confidence_level = 0.95,
                                                 require_heldout_calibration = FALSE) {
  n <- length(predicted_value)
  observed_aligned <- gp_ml_align_observed_targets(
    observed_y = observed_y,
    train_test_label = train_test_label,
    target_length = n
  )
  n_train <- if (!is.null(train_test_label) && length(train_test_label) == n) {
    sum(train_test_label != "Test", na.rm = TRUE)
  } else if (!is.null(observed_y)) {
    length(observed_y)
  } else 0L
  reference_rows <- rep(TRUE, n)
  if (!is.null(train_test_label) && length(train_test_label) == n &&
      any(train_test_label != "Test", na.rm = TRUE)) {
    reference_rows <- !is.na(train_test_label) & train_test_label != "Test"
  }

  oob_res <- gp_ml_bootstrap_oob_summary(
    boot_results = boot_results,
    observed_y = observed_y,
    n_train = n_train
  )

  cal_rows <- is.finite(observed_aligned) &
    is.finite(predicted_value) &
    is.finite(pred_variances)
  if (!is.null(train_test_label) && length(train_test_label) == n && any(train_test_label != "Test")) {
    cal_rows <- cal_rows & train_test_label != "Test"
  }

  calibration_factor <- 1
  calibration_rmse <- NA_real_
  calibration_mean_pev <- NA_real_
  empirical_coverage <- NA_real_
  mean_interval_width <- NA_real_
  calibration_source <- "observed_rows"
  rank_spearman <- NA_real_

  heldout_res <- gp_ml_heldout_residual_summary(
    calibration = heldout_calibration,
    confidence_level = confidence_level
  )

  if (!is.null(heldout_res)) {
    target_uncertainty_source <- "heldout_cross_fitted_residuals"
    target_res <- gp_ml_cross_fitted_target_uncertainty(
      calibration = heldout_calibration,
      final_prediction = predicted_value,
      confidence_level = confidence_level
    )
    if (!is.null(target_res)) {
      calibrated_pev <- target_res$prediction_mse
      calibrated_se <- sqrt(pmax(calibrated_pev, 0))
      calibration_factor <- if (
        is.finite(target_res$mean_calibration_factor %||% NA_real_) &&
          target_res$mean_calibration_factor >= 0
      ) {
        sqrt(target_res$mean_calibration_factor)
      } else {
        NA_real_
      }
      calibrated_lower <- target_res$lower
      calibrated_upper <- target_res$upper
      interval_method <- target_res$interval_method %||% paste(
        "cross-fitted fold-prediction absolute-residual interval;",
        "empirical marginal calibration"
      )
      error_estimand <- target_res$error_estimand %||% paste(
        "target-specific cross-fitted predictive MSE;",
        "not mixed-model PEV"
      )
      target_uncertainty_source <- target_res$source %||%
        "heldout_cross_fitted_residuals"
    } else {
      row_instability <- suppressWarnings(as.numeric(pred_variances))
      if (length(row_instability) != n) {
        row_instability <- rep(NA_real_, n)
      }
      reference_instability <- row_instability[
        reference_rows & is.finite(row_instability) & row_instability >= 0
      ]
      instability_varies <- gp_ml_has_substantive_target_variation(
        row_instability,
        reference_rows = reference_rows
      )
      mean_instability <- if (length(reference_instability)) {
        mean(reference_instability)
      } else {
        NA_real_
      }

      if (isTRUE(instability_varies) && is.finite(mean_instability) &&
          mean_instability > 0) {
        row_instability[!is.finite(row_instability) | row_instability < 0] <-
          mean_instability
        if (heldout_res$mse >= mean_instability) {
          residual_component <- heldout_res$mse - mean_instability
          calibration_factor <- 1
        } else {
          residual_component <- 0
          calibration_factor <- sqrt(heldout_res$mse / mean_instability)
        }
        calibrated_pev <- residual_component +
          row_instability * (calibration_factor^2)
        calibrated_pev <- pmax(calibrated_pev, 0)
        calibrated_se <- sqrt(calibrated_pev)
        error_estimand <- paste(
          "row-specific bootstrap prediction instability plus a held-out",
          "cross-fitted residual component, mean-calibrated to predictive MSE;",
          "not mixed-model PEV"
        )
      } else {
        calibrated_se <- rep(NA_real_, n)
        calibrated_pev <- rep(NA_real_, n)
        calibration_factor <- NA_real_
        error_estimand <- gp_ml_unavailable_target_uncertainty_estimand()
      }
      calibrated_lower <- predicted_value - heldout_res$interval_radius
      calibrated_upper <- predicted_value + heldout_res$interval_radius
      interval_method <- paste(
        "cross-fitted absolute-residual interval with finite-sample quantile correction;",
        "empirical marginal calibration"
      )
    }
    if (!gp_ml_has_substantive_target_variation(
      calibrated_se,
      reference_rows = reference_rows
    )) {
      calibrated_se <- rep(NA_real_, n)
      calibrated_pev <- rep(NA_real_, n)
      calibration_factor <- NA_real_
      error_estimand <- gp_ml_unavailable_target_uncertainty_estimand()
    }
    return(list(
      calibration_source = target_uncertainty_source %||%
        "heldout_cross_fitted_residuals",
      calibration_factor = calibration_factor,
      calibration_rmse = heldout_res$rmse,
      calibration_mean_pev = heldout_res$mse,
      empirical_coverage = target_res$empirical_coverage %||%
        heldout_res$calibration_coverage,
      mean_interval_width = mean(calibrated_upper - calibrated_lower, na.rm = TRUE),
      rank_spearman = suppressWarnings(stats::cor(
        heldout_calibration$heldout_observed,
        heldout_calibration$heldout_predictions,
        method = "spearman",
        use = "complete.obs"
      )),
      oob_summary = oob_res,
      calibrated_standard_error = calibrated_se,
      calibrated_prediction_error_var = calibrated_pev,
      calibrated_lower_bound = calibrated_lower,
      calibrated_upper_bound = calibrated_upper,
      observed_aligned = observed_aligned,
      prediction_error_estimand = error_estimand,
      prediction_interval_method = interval_method,
      prediction_interval_nominal_coverage = confidence_level,
      prediction_interval_calibration_n = heldout_res$calibration_n
    ))
  }

  if (isTRUE(require_heldout_calibration)) {
    return(list(
      calibration_source = "unavailable_no_heldout_calibration",
      calibration_factor = NA_real_,
      calibration_rmse = NA_real_,
      calibration_mean_pev = NA_real_,
      empirical_coverage = NA_real_,
      mean_interval_width = NA_real_,
      rank_spearman = NA_real_,
      oob_summary = oob_res,
      calibrated_standard_error = rep(NA_real_, n),
      calibrated_prediction_error_var = rep(NA_real_, n),
      calibrated_lower_bound = rep(NA_real_, n),
      calibrated_upper_bound = rep(NA_real_, n),
      observed_aligned = observed_aligned,
      prediction_error_estimand =
        "unavailable without held-out predictions",
      prediction_interval_method =
        "unavailable without held-out calibration",
      prediction_interval_nominal_coverage = confidence_level,
      prediction_interval_calibration_n = 0L
    ))
  }

  if (!is.null(oob_res)) {
    calibration_rmse <- oob_res$calibration_rmse
    calibration_mean_pev <- oob_res$calibration_mean_pev
    rank_spearman <- oob_res$rank_spearman
    calibration_source <- "bootstrap_oob"
    if (is.finite(calibration_rmse) &&
        is.finite(calibration_mean_pev) &&
        calibration_mean_pev > 0) {
      calibration_factor <- calibration_rmse / sqrt(calibration_mean_pev)
    }
    if (!is.null(lower_bound) && !is.null(upper_bound)) {
      oob_valid <- is.finite(oob_res$oob_mean) &
        is.finite(oob_res$oob_lower) &
        is.finite(oob_res$oob_upper) &
        is.finite(observed_y)
      if (sum(oob_valid, na.rm = TRUE) >= 1L) {
        oob_lower_cal <- oob_res$oob_mean - calibration_factor * (oob_res$oob_mean - oob_res$oob_lower)
        oob_upper_cal <- oob_res$oob_mean + calibration_factor * (oob_res$oob_upper - oob_res$oob_mean)
        empirical_coverage <- mean(
          observed_y[oob_valid] >= oob_lower_cal[oob_valid] &
            observed_y[oob_valid] <= oob_upper_cal[oob_valid],
          na.rm = TRUE
        )
        mean_interval_width <- mean(oob_upper_cal[oob_valid] - oob_lower_cal[oob_valid], na.rm = TRUE)
      }
    }
  }

  if (!is.null(oob_res) || sum(cal_rows, na.rm = TRUE) >= 2L) {
    if (is.null(oob_res)) {
    calibration_rmse <- sqrt(mean((observed_aligned[cal_rows] - predicted_value[cal_rows])^2, na.rm = TRUE))
    calibration_mean_pev <- mean(pred_variances[cal_rows], na.rm = TRUE)
    if (is.finite(calibration_rmse) &&
        is.finite(calibration_mean_pev) &&
        calibration_mean_pev > 0) {
      calibration_factor <- calibration_rmse / sqrt(calibration_mean_pev)
    }
    }
  }

  if (!is.finite(calibration_factor) || calibration_factor <= 0) {
    calibration_factor <- 1
  }

  calibrated_se <- pred_se * calibration_factor
  calibrated_pev <- pred_variances * (calibration_factor^2)

  calibrated_lower <- lower_bound
  calibrated_upper <- upper_bound
  if (!is.null(lower_bound) && !is.null(upper_bound)) {
    calibrated_lower <- predicted_value - calibration_factor * (predicted_value - lower_bound)
    calibrated_upper <- predicted_value + calibration_factor * (upper_bound - predicted_value)
    if (sum(cal_rows, na.rm = TRUE) >= 1L) {
      empirical_coverage <- mean(
        observed_aligned[cal_rows] >= calibrated_lower[cal_rows] &
          observed_aligned[cal_rows] <= calibrated_upper[cal_rows],
        na.rm = TRUE
      )
      mean_interval_width <- mean(calibrated_upper[cal_rows] - calibrated_lower[cal_rows], na.rm = TRUE)
    }
  }

  list(
    calibration_source = calibration_source,
    calibration_factor = calibration_factor,
    calibration_rmse = calibration_rmse,
    calibration_mean_pev = calibration_mean_pev,
    empirical_coverage = empirical_coverage,
    mean_interval_width = mean_interval_width,
    rank_spearman = rank_spearman,
    oob_summary = oob_res,
    calibrated_standard_error = calibrated_se,
    calibrated_prediction_error_var = calibrated_pev,
    calibrated_lower_bound = calibrated_lower,
    calibrated_upper_bound = calibrated_upper,
    observed_aligned = observed_aligned,
    prediction_error_estimand =
      "bootstrap prediction instability scaled to observed error; not mixed-model PEV",
    prediction_interval_method =
      "scaled bootstrap interval; nominal coverage not guaranteed",
    prediction_interval_nominal_coverage = confidence_level,
    prediction_interval_calibration_n = sum(cal_rows, na.rm = TRUE)
  )
}

gp_ml_unfamiliarity_remark <- function(unfamiliarity,
                                       low_cut = 1 / 3,
                                       high_cut = 2 / 3) {
  ifelse(
    unfamiliarity <= low_cut,
    "Familiar",
    ifelse(unfamiliarity >= high_cut, "Highly Unfamiliar", "Moderately Unfamiliar")
  )
}

# Leading principal-component scores of `x` fitted on the training rows
# (uncentred: x is already standardised), for all rows. With more columns than
# training rows -- the usual marker matrix -- they come from the n x n Gram
# matrix: Xt Xt' = U D^2 U', train scores = U D, all scores = (X Xt') U / D.
# This equals prcomp()'s scores up to the sign of each component (irrelevant
# to the Euclidean distances computed from them) without a full p-wide SVD,
# which took ~48 s and several GB for 1,451 x 10,346 markers.
gp_ml_pca_scores <- function(x, train_idx, rank) {
  xt <- x[train_idx, , drop = FALSE]
  if (ncol(xt) <= nrow(xt)) {
    pca <- stats::prcomp(xt, center = FALSE, scale. = FALSE, rank. = rank)
    return(list(train = pca$x[, seq_len(rank), drop = FALSE],
                all = x %*% pca$rotation[, seq_len(rank), drop = FALSE]))
  }
  # One product serves both: the Gram matrix is its training-row block.
  cross <- tcrossprod(x, xt)
  eg <- eigen(cross[train_idx, , drop = FALSE], symmetric = TRUE)
  d <- sqrt(pmax(eg$values[seq_len(rank)], 0))
  keep <- d > max(d) * 1e-10
  u <- eg$vectors[, seq_len(rank), drop = FALSE][, keep, drop = FALSE]
  d <- d[keep]
  list(train = sweep(u, 2L, d, `*`),
       all = sweep(cross %*% u, 2L, d, `/`))
}

gp_ml_predictor_unfamiliarity <- function(predictor_matrix,
                                          train_test_label = NULL,
                                          max_components = 20L,
                                          k_neighbors = 5L) {
  if (is.null(predictor_matrix)) {
    return(NULL)
  }

  x <- tryCatch(as.matrix(predictor_matrix), error = function(e) NULL)
  if (is.null(x) || !is.matrix(x) || !nrow(x) || !ncol(x)) {
    return(NULL)
  }
  storage.mode(x) <- "double"
  n <- nrow(x)
  if (is.null(train_test_label) || length(train_test_label) != n) {
    train_test_label <- rep("Train", n)
  }
  train_idx <- which(train_test_label != "Test")
  if (length(train_idx) < 5L) {
    train_idx <- seq_len(n)
  }
  if (length(train_idx) < 5L) {
    return(NULL)
  }

  train_x <- x[train_idx, , drop = FALSE]
  keep_cols <- apply(train_x, 2, function(col) any(is.finite(col)))
  if (!any(keep_cols)) {
    return(NULL)
  }
  x <- x[, keep_cols, drop = FALSE]
  train_x <- x[train_idx, , drop = FALSE]

  train_medians <- apply(train_x, 2, function(col) {
    med <- stats::median(col[is.finite(col)], na.rm = TRUE)
    if (is.finite(med)) med else 0
  })

  for (j in seq_len(ncol(x))) {
    bad <- !is.finite(x[, j])
    if (any(bad)) {
      x[bad, j] <- train_medians[[j]]
    }
  }
  train_x <- x[train_idx, , drop = FALSE]

  train_sds <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  keep_var <- is.finite(train_sds) & train_sds > 0
  if (!any(keep_var)) {
    return(NULL)
  }
  x <- x[, keep_var, drop = FALSE]
  train_x <- x[train_idx, , drop = FALSE]

  rm(train_x)
  x_scaled <- scale(x, center = colMeans(x[train_idx, , drop = FALSE], na.rm = TRUE), scale = train_sds[keep_var])
  rm(x)
  x_scaled[!is.finite(x_scaled)] <- 0

  rank_max <- min(as.integer(max_components), ncol(x_scaled), length(train_idx) - 1L)
  if (!is.finite(rank_max) || rank_max < 1L) {
    return(NULL)
  }

  scores <- tryCatch(
    gp_ml_pca_scores(x_scaled, train_idx, rank_max),
    error = function(e) NULL
  )
  if (is.null(scores)) {
    return(NULL)
  }
  train_scores <- scores$train
  all_scores <- scores$all
  rm(x_scaled, scores)

  k_train <- max(1L, min(as.integer(k_neighbors), nrow(train_scores) - 1L))
  k_test <- max(1L, min(as.integer(k_neighbors), nrow(train_scores)))
  if (k_test < 1L) {
    return(NULL)
  }

  row_dist_mean <- function(target, reference, k, self_row = NULL) {
    d <- sqrt(rowSums((reference - matrix(target, nrow = nrow(reference), ncol = length(target), byrow = TRUE))^2))
    if (!is.null(self_row) && self_row >= 1L && self_row <= length(d)) {
      d[[self_row]] <- Inf
    }
    ord <- order(d, na.last = NA)
    ord <- ord[seq_len(min(k, length(ord)))]
    mean(d[ord], na.rm = TRUE)
  }

  train_distance <- vapply(seq_len(nrow(train_scores)), function(i) {
    row_dist_mean(train_scores[i, ], train_scores, k = k_train, self_row = i)
  }, numeric(1))

  all_distance <- vapply(seq_len(nrow(all_scores)), function(i) {
    self_row <- if (i %in% train_idx) match(i, train_idx) else NULL
    row_dist_mean(all_scores[i, ], train_scores, k = if (is.null(self_row)) k_test else k_train, self_row = self_row)
  }, numeric(1))

  dist_center <- stats::median(train_distance, na.rm = TRUE)
  dist_scale <- stats::mad(train_distance, center = dist_center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(dist_scale) || dist_scale <= 0) {
    dist_scale <- stats::sd(train_distance, na.rm = TRUE)
  }
  if (!is.finite(dist_scale) || dist_scale <= 0) {
    dist_scale <- 1
  }

  unfamiliarity <- stats::pnorm((all_distance - dist_center) / dist_scale)
  unfamiliarity <- pmax(0, pmin(1, unfamiliarity))
  low_cut <- suppressWarnings(stats::quantile(unfamiliarity[train_idx], probs = 1 / 3, na.rm = TRUE, names = FALSE))
  high_cut <- suppressWarnings(stats::quantile(unfamiliarity[train_idx], probs = 2 / 3, na.rm = TRUE, names = FALSE))
  if (!is.finite(low_cut) || !is.finite(high_cut) || low_cut >= high_cut) {
    low_cut <- 1 / 3
    high_cut <- 2 / 3
  }

  list(
    unfamiliarity = unfamiliarity,
    distance = all_distance,
    training_distance = train_distance,
    low_cut = as.numeric(low_cut),
    high_cut = as.numeric(high_cut),
    remarks = gp_ml_unfamiliarity_remark(unfamiliarity, low_cut = low_cut, high_cut = high_cut),
    calibration_source = "predictor_space_novelty"
  )
}

gp_ml_gaussian_prediction_output <- function(ids,
                                             predicted_value,
                                             train_test_label,
                                             pred_se,
                                             pred_variances,
                                             lower_bound,
                                             upper_bound,
                                             uncertainty,
                                             uncertainty_remarks,
                                             reference_variance,
                                             observed_y = NULL,
                                             boot_results = NULL,
                                             risk_calibration = NULL,
                                             predictor_matrix = NULL,
                                             gen_name,
                                             high_reliability_thres = 0.9,
                                             low_reliability_thres = 0.5,
                                             confidence_level = 0.95,
                                             require_heldout_calibration = FALSE) {
  cal_res <- gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = predicted_value,
    pred_se = pred_se,
    pred_variances = pred_variances,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    observed_y = observed_y,
    train_test_label = train_test_label,
    boot_results = boot_results,
    heldout_calibration = risk_calibration,
    confidence_level = confidence_level,
    require_heldout_calibration = require_heldout_calibration
  )

  stability_input <- gp_ml_resolve_stability_input(
    pred_variances = pred_variances,
    predicted_value = predicted_value,
    observed_aligned = cal_res$observed_aligned,
    reference_variance = reference_variance,
    train_test_label = train_test_label
  )
  stability_reference_variance <- stability_input$reference_variance
  stability_res <- gp_ml_gaussian_stability(
    prediction_error_var = stability_input$variance_input,
    reference_variance = stability_reference_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  reliability_res <- gp_ml_predictive_reliability(
    prediction_error_var = cal_res$calibrated_prediction_error_var,
    reference_variance = reference_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  rank_res <- gp_ml_local_rank_confidence(
    predicted_value = predicted_value,
    pred_se = cal_res$calibrated_standard_error,
    groups = train_test_label
  )
  oob_risk_cal <- NULL
  if (!is.null(cal_res$oob_summary)) {
    oob_risk_cal <- gp_ml_rank_risk_calibration(
      observed_y = observed_y,
      oob_mean = cal_res$oob_summary$oob_mean,
      oob_var = cal_res$oob_summary$oob_var,
      reference_variance = reference_variance
    )
  }
  risk_cal <- if (!is.null(oob_risk_cal)) oob_risk_cal else risk_calibration
  pairwise_cal <- NULL
  if (!is.null(observed_y)) {
    if (!is.null(cal_res$oob_summary)) {
      pairwise_cal <- gp_ml_pairwise_rank_risk_calibration(
        predicted_value = cal_res$oob_summary$oob_mean,
        observed_y = observed_y,
        pred_se = sqrt(pmax(cal_res$oob_summary$oob_var, 0)),
        groups = rep("Train", length(observed_y)),
        source_label = "pairwise_oob"
      )
    }
  }
  if (is.null(risk_cal) && !is.null(pairwise_cal)) {
    risk_cal <- pairwise_cal
  }
  confidence_res <- gp_ml_gaussian_confidence(
    prediction_error_var = cal_res$calibrated_prediction_error_var,
    reference_variance = reference_variance,
    rank_confidence = rank_res$rank_confidence,
    high_confidence_thres = if (!is.null(risk_cal)) 1 - risk_cal$risk_low_cut else high_reliability_thres,
    low_confidence_thres = if (!is.null(risk_cal)) 1 - risk_cal$risk_high_cut else low_reliability_thres
  )
  if (!is.null(risk_cal)) {
    if (!is.null(risk_cal$fit) && grepl("^pairwise_", risk_cal$calibration_source %||% "")) {
      calibrated_risk <- gp_ml_apply_pairwise_rank_risk_calibrator(
        predicted_value = predicted_value,
        pred_se = cal_res$calibrated_standard_error,
        groups = train_test_label,
        calibrator = risk_cal
      )
      if (!is.null(calibrated_risk) && length(calibrated_risk) == length(confidence_res$risk_score)) {
        confidence_res$risk_score <- calibrated_risk
        confidence_res$confidence <- 1 - calibrated_risk
      }
    }
    confidence_res$risk_remarks <- ifelse(
      confidence_res$risk_score >= risk_cal$risk_high_cut,
      "High Risk",
      ifelse(confidence_res$risk_score <= risk_cal$risk_low_cut, "Low Risk", "Moderate Risk")
    )
    confidence_res$confidence_remarks <- ifelse(
      confidence_res$risk_remarks == "Low Risk",
      "High Confidence",
      ifelse(confidence_res$risk_remarks == "High Risk", "Low Confidence", "Moderate Confidence")
    )
  }

  risk_basis <- ifelse(
    is.finite(confidence_res$risk_score),
    "Calibrated_Final_Prediction",
    "Unavailable_No_Target_Specific_PEV"
  )

  if (!is.null(cal_res$oob_summary) &&
      !is.null(train_test_label) &&
      length(train_test_label) == length(predicted_value)) {
    train_idx <- which(train_test_label != "Test")
    oob_mean <- cal_res$oob_summary$oob_mean
    oob_var <- cal_res$oob_summary$oob_var
    if (length(train_idx) == length(oob_mean) && length(oob_var) == length(oob_mean)) {
      oob_rank_res <- gp_ml_local_rank_confidence(
        predicted_value = oob_mean,
        pred_se = sqrt(pmax(oob_var, 0)),
        groups = rep("Train", length(oob_mean))
      )
      oob_conf <- pmin(
        exp(-oob_var / reference_variance),
        oob_rank_res$rank_confidence
      )
      oob_risk <- 1 - oob_conf
      if (!is.null(risk_cal) && !is.null(risk_cal$fit) &&
          grepl("^pairwise_", risk_cal$calibration_source %||% "")) {
        pairwise_oob_risk <- gp_ml_apply_pairwise_rank_risk_calibrator(
          predicted_value = oob_mean,
          pred_se = sqrt(pmax(oob_var, 0)),
          groups = rep("Train", length(oob_mean)),
          calibrator = risk_cal
        )
        if (!is.null(pairwise_oob_risk) && length(pairwise_oob_risk) == length(oob_risk)) {
          oob_risk <- pairwise_oob_risk
        }
      }
      oob_risk <- pmax(0, pmin(1, oob_risk))
      confidence_res$risk_score[train_idx] <- oob_risk
      confidence_res$confidence[train_idx] <- 1 - oob_risk
      risk_basis[train_idx] <- "Heldout_Training_OOB"
      if (!is.null(risk_cal)) {
        confidence_res$risk_remarks[train_idx] <- ifelse(
          oob_risk >= risk_cal$risk_high_cut,
          "High Risk",
          ifelse(oob_risk <= risk_cal$risk_low_cut, "Low Risk", "Moderate Risk")
        )
        confidence_res$confidence_remarks[train_idx] <- ifelse(
          confidence_res$risk_remarks[train_idx] == "Low Risk",
          "High Confidence",
          ifelse(confidence_res$risk_remarks[train_idx] == "High Risk", "Low Confidence", "Moderate Confidence")
        )
      }
    }
  }

  unfamiliarity_res <- gp_ml_predictor_unfamiliarity(
    predictor_matrix = predictor_matrix,
    train_test_label = train_test_label
  )
  if (!is.null(unfamiliarity_res)) {
    base_risk_score <- confidence_res$risk_score
    base_confidence <- confidence_res$confidence
    base_risk_remarks <- confidence_res$risk_remarks
    base_confidence_remarks <- confidence_res$confidence_remarks
    combined_risk <- 1 - (1 - confidence_res$risk_score) * (1 - unfamiliarity_res$unfamiliarity)
    confidence_res$risk_score <- pmax(0, pmin(1, combined_risk))
    confidence_res$confidence <- 1 - confidence_res$risk_score
    if (!is.null(risk_cal)) {
      confidence_res$risk_remarks <- ifelse(
        confidence_res$risk_score >= risk_cal$risk_high_cut,
        "High Risk",
        ifelse(confidence_res$risk_score <= risk_cal$risk_low_cut, "Low Risk", "Moderate Risk")
      )
    }
    confidence_res$confidence_remarks <- ifelse(
      confidence_res$risk_remarks == "Low Risk",
      "High Confidence",
      ifelse(confidence_res$risk_remarks == "High Risk", "Low Confidence", "Moderate Confidence")
    )
  } else {
    base_risk_score <- confidence_res$risk_score
    base_confidence <- confidence_res$confidence
    base_risk_remarks <- confidence_res$risk_remarks
    base_confidence_remarks <- confidence_res$confidence_remarks
  }

  out <- data.frame(
    name = ids,
    Predicted_value = predicted_value,
    Train_Test_Label = train_test_label,
    Observed_value = cal_res$observed_aligned,
    Standard_error = cal_res$calibrated_standard_error,
    PEV = cal_res$calibrated_prediction_error_var,
    lower_bound = cal_res$calibrated_lower_bound,
    upper_bound = cal_res$calibrated_upper_bound,
    Uncertainty = cal_res$calibrated_upper_bound - cal_res$calibrated_lower_bound,
    Uncertainty_remarks = rep(
      cal_res$calibration_source %||% uncertainty_remarks,
      length(predicted_value)
    ),
    Prediction_stability = stability_res$stability,
    Prediction_stability_remarks = stability_res$remarks,
    Prediction_stability_percentage = stability_res$stability_percentage,
    Prediction_stability_reference_variance = stability_reference_variance,
    Reliability = reliability_res$reliability,
    Reliability_remarks = reliability_res$remarks,
    Reliability_variance_input = reliability_res$variance_input,
    Reliability_reference_variance = reliability_res$reference_variance,
    Reliability_basis = reliability_res$basis,
    Prediction_uncertainty_source = rep(
      cal_res$calibration_source %||% "bootstrap_prediction_spread",
      length(predicted_value)
    ),
    PEV_basis = rep(
      cal_res$prediction_error_estimand %||% "unspecified predictive-error proxy",
      length(predicted_value)
    ),
    Prediction_interval_method = rep(
      cal_res$prediction_interval_method %||% "unspecified",
      length(predicted_value)
    ),
    Prediction_interval_nominal_coverage = rep(
      cal_res$prediction_interval_nominal_coverage %||% confidence_level,
      length(predicted_value)
    ),
    Prediction_interval_calibration_n = rep(
      cal_res$prediction_interval_calibration_n %||% NA_integer_,
      length(predicted_value)
    ),
    Prediction_confidence = confidence_res$confidence,
    Prediction_confidence_remarks = confidence_res$confidence_remarks,
    Prediction_risk_score = confidence_res$risk_score,
    Prediction_risk_remarks = confidence_res$risk_remarks,
    Prediction_risk_basis = risk_basis,
    Prediction_unfamiliarity = if (!is.null(unfamiliarity_res)) unfamiliarity_res$unfamiliarity else NA_real_,
    Prediction_unfamiliarity_remarks = if (!is.null(unfamiliarity_res)) unfamiliarity_res$remarks else NA_character_,
    Prediction_unfamiliarity_basis = if (!is.null(unfamiliarity_res)) unfamiliarity_res$calibration_source else NA_character_,
    Rank_instability_risk = base_risk_score,
    Rank_instability_risk_remarks = base_risk_remarks,
    Rank_instability_risk_basis = risk_basis,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- gen_name
  attr(out, "rank_confidence") <- rank_res
  attr(out, "rank_risk_calibration") <- risk_cal
  attr(out, "rank_risk_calibration_candidates") <- list(
    selected = risk_cal,
    pairwise_candidate = pairwise_cal
  )
  attr(out, "uncertainty_calibration") <- cal_res
  attr(out, "prediction_unfamiliarity") <- unfamiliarity_res
  out
}

gp_ml_stability_labels <- function(metric, remarks = NULL, percentage = NULL) {
  mapped_remarks <- remarks
  if (!is.null(mapped_remarks)) {
    mapped_remarks <- dplyr::case_when(
      mapped_remarks == "Reliable" ~ "Stable",
      mapped_remarks == "Acceptable" ~ "Moderately Stable",
      mapped_remarks == "Unreliable" ~ "Unstable",
      TRUE ~ as.character(mapped_remarks)
    )
  }

  list(
    Prediction_stability = metric,
    Prediction_stability_remarks = mapped_remarks,
    Prediction_stability_percentage = percentage
  )
}

gp_ml_is_classification_family <- function(response_family = "gaussian", y = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  fam %in% c("binary", "ordinal", "multiclass")
}

gp_ml_binary_levels <- function(y) {
  y_no_na <- y[!is.na(y)]
  levs <- if (is.factor(y)) {
    levels(droplevels(y))
  } else if (is.logical(y) || is.character(y)) {
    unique(as.character(y_no_na))
  } else {
    as.character(sort(unique(y_no_na)))
  }
  levs <- as.character(levs)
  if (length(levs) != 2L) {
    stop("Binary ML support requires exactly two observed classes.", call. = FALSE)
  }
  levs
}

gp_ml_binary_encode <- function(y, class_levels = NULL) {
  levs <- class_levels %||% gp_ml_binary_levels(y)
  as.integer(factor(as.character(y), levels = levs)) - 1L
}

gp_ml_binary_prob_column <- function(prob_matrix, class_levels = NULL) {
  prob_matrix <- as.data.frame(prob_matrix)
  if (ncol(prob_matrix) == 1L) {
    return(as.numeric(prob_matrix[[1]]))
  }
  levs <- class_levels %||% colnames(prob_matrix)
  pos <- levs[min(2L, length(levs))]
  pos_idx <- match(pos, colnames(prob_matrix))
  if (is.na(pos_idx)) {
    pos_idx <- ncol(prob_matrix)
  }
  as.numeric(prob_matrix[[pos_idx]])
}

gp_ml_classification_prediction_output <- function(ids,
                                                   predicted_prob,
                                                   train_test_label,
                                                   pred_se,
                                                   pred_variances,
                                                   genetic_var,
                                                   lower_bound,
                                                   upper_bound,
                                                   uncertainty,
                                                   uncertainty_remarks,
                                                   class_levels,
                                                   gen_name,
                                                   high_confidence = 0.9,
                                                   low_confidence = 0.5) {
  prob <- pmax(0, pmin(1, as.numeric(predicted_prob)))
  lower <- pmax(0, pmin(1, as.numeric(lower_bound)))
  upper <- pmax(0, pmin(1, as.numeric(upper_bound)))
  cls <- gp_classification_prediction_summary(
    prob = prob,
    class_levels = class_levels,
    high_confidence = high_confidence,
    low_confidence = low_confidence
  )
  out <- data.frame(
    name = ids,
    Predicted_value = prob,
    Train_Test_Label = train_test_label,
    Standard_error = pred_se,
    PEV = pred_variances,
    Genetic_variance = rep(genetic_var, length(prob)),
    lower_bound = lower,
    upper_bound = upper,
    Uncertainty = uncertainty,
    Uncertainty_remarks = uncertainty_remarks,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, cls)
  names(out)[1] <- gen_name
  out
}

gp_prob_column_names <- function(class_levels) {
  paste0("Prob_", make.names(as.character(class_levels), unique = TRUE))
}

gp_multiclass_bootstrap_metrics <- function(boot_matrix, n_obs, class_levels, train_test_label = NULL) {
  n_classes <- length(class_levels)
  if (ncol(boot_matrix) != n_obs * n_classes) {
    stop("Multiclass bootstrap matrix shape does not match n_obs * n_classes.", call. = FALSE)
  }
  reps <- lapply(seq_len(nrow(boot_matrix)), function(i) {
    matrix(boot_matrix[i, ], nrow = n_obs, ncol = n_classes)
  })
  prob_mean <- Reduce(`+`, reps) / length(reps)
  conf_rep <- vapply(reps, function(mat) apply(mat, 1, max), numeric(n_obs))
  conf_rep <- t(conf_rep)
  pred_conf <- apply(conf_rep, 2, mean)
  pred_var <- apply(conf_rep, 2, var)
  pred_se <- apply(conf_rep, 2, sd)

  ref_idx <- rep(TRUE, n_obs)
  if (!is.null(train_test_label) && length(train_test_label) == n_obs) {
    ref_idx <- train_test_label != "Test"
    if (!any(ref_idx, na.rm = TRUE)) ref_idx <- rep(TRUE, n_obs)
  }
  genetic_var <- stats::var(pred_conf[ref_idx], na.rm = TRUE)
  if (!is.finite(genetic_var) || genetic_var < 0) genetic_var <- NA_real_

  list(
    mean_prob = prob_mean,
    predicted_confidence = pred_conf,
    prediction_error_var = pred_var,
    standard_error = pred_se,
    genetic_var = genetic_var,
    lower_conf = apply(conf_rep, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    upper_conf = apply(conf_rep, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}

gp_multiclass_prediction_output <- function(ids,
                                            prob_mean,
                                            train_test_label,
                                            pred_se,
                                            pred_variances,
                                            genetic_var,
                                            lower_bound,
                                            upper_bound,
                                            uncertainty,
                                            uncertainty_remarks,
                                            class_levels,
                                            gen_name,
                                            high_confidence = 0.9,
                                            low_confidence = 0.5) {
  if (!is.matrix(prob_mean)) {
    prob_mean <- as.matrix(prob_mean)
    if (!is.null(class_levels) && length(class_levels) > 0L && ncol(prob_mean) != length(class_levels)) {
      if (length(prob_mean) %% length(class_levels) == 0L) {
        prob_mean <- matrix(as.numeric(prob_mean), ncol = length(class_levels), byrow = FALSE)
      }
    }
  }

  if (!is.null(class_levels) && length(class_levels) > 0L && ncol(prob_mean) != length(class_levels)) {
    stop("Multiclass probability matrix width does not match class_levels.", call. = FALSE)
  }

  n_obs <- nrow(prob_mean)
  collapse_multiclass_vector <- function(x, fun = mean, default = NA_real_) {
    if (is.null(x)) {
      return(rep(default, n_obs))
    }
    if (length(x) == n_obs) {
      return(x)
    }
    if (length(x) %% n_obs == 0L) {
      mat <- matrix(x, nrow = n_obs, byrow = FALSE)
      return(apply(mat, 1, fun))
    }
    rep(default, n_obs)
  }

  collapse_multiclass_labels <- function(x, default = "Unknown") {
    if (is.null(x)) {
      return(rep(default, n_obs))
    }
    if (length(x) == n_obs) {
      return(as.character(x))
    }
    if (length(x) %% n_obs == 0L) {
      mat <- matrix(as.character(x), nrow = n_obs, byrow = FALSE)
      return(apply(mat, 1, function(v) {
        tab <- sort(table(v), decreasing = TRUE)
        names(tab)[1]
      }))
    }
    rep(default, n_obs)
  }

  uncertainty <- collapse_multiclass_vector(uncertainty, fun = mean, default = NA_real_)
  uncertainty_remarks <- collapse_multiclass_labels(uncertainty_remarks, default = "Unknown")

  cls <- gp_classification_prediction_summary(
    prob = prob_mean,
    class_levels = class_levels,
    high_confidence = high_confidence,
    low_confidence = low_confidence
  )
  out <- data.frame(
    name = ids,
    Predicted_value = cls$Predicted_class,
    Train_Test_Label = train_test_label,
    Standard_error = pred_se,
    PEV = pred_variances,
    Genetic_variance = rep(genetic_var, nrow(prob_mean)),
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    Uncertainty = uncertainty,
    Uncertainty_remarks = uncertainty_remarks,
    cls,
    stringsAsFactors = FALSE
  )
  prob_df <- as.data.frame(prob_mean, stringsAsFactors = FALSE)
  names(prob_df) <- gp_prob_column_names(class_levels)
  out <- cbind(out, prob_df)
  names(out)[1] <- gen_name
  out
}

compose_initializers <- function(...) {
  inits <- list(...)
  function() {
    for (f in inits) {
      if (is.function(f)) {
        try(f(), silent = TRUE)
      }
    }
  }
}

normalize_cv_token <- function(x) {
  if (is.null(x)) {
    return(x)
  }

  out <- as.character(x)
  keep <- !(is.na(out) | !nzchar(out))
  if (!any(keep)) {
    return(out)
  }

  z <- tolower(out[keep])
  z <- gsub("[^a-z0-9]+", "_", z)
  z <- gsub("_+", "_", z)
  z <- sub("^_", "", z)
  z <- sub("_$", "", z)
  out[keep] <- z
  out
}

gp_resolve_cv_sampling_method <- function(cross_validation_meth, sampling_method = NULL) {
  if (!is.null(sampling_method) && length(sampling_method)) {
    method <- normalize_cv_token(as.character(sampling_method[[1]]))
    if (!method %in% c("stratified", "unstratified")) {
      stop("`sampling_method` must be either 'stratified' or 'unstratified'.", call. = FALSE)
    }
    return(method)
  }

  cv_token <- normalize_cv_token(as.character(cross_validation_meth %||% ""))
  stratified_tokens <- c(
    "stratified_hold_out",
    "repeated_stratified_hold_out",
    "stratified_k_folds",
    "repeated_stratified_k_folds"
  )
  unstratified_tokens <- c(
    "hold_out",
    "repeated_hold_out",
    "k_folds",
    "repeated_k_folds",
    "leave_one_out"
  )

  if (cv_token %in% stratified_tokens) {
    return("stratified")
  }
  if (cv_token %in% unstratified_tokens) {
    return("unstratified")
  }

  sampling_method
}

gp_normalize_cv_method_name <- function(cross_validation_meth, allowed_methods) {
  if (is.null(cross_validation_meth) || !length(cross_validation_meth)) {
    return(cross_validation_meth)
  }
  if (length(cross_validation_meth) > 1L) {
    return(cross_validation_meth)
  }

  token <- normalize_cv_token(as.character(cross_validation_meth[[1]]))
  lookup <- stats::setNames(allowed_methods, normalize_cv_token(allowed_methods))
  if (token %in% names(lookup)) {
    return(unname(lookup[[token]]))
  }
  cross_validation_meth
}

map_to_canonical <- function(x, canonical, friendly) {
  if (is.null(x) || !length(x)) {
    return(character())
  }

  key <- tolower(c(canonical, friendly))
  val <- c(canonical, canonical)
  lut <- stats::setNames(val, key)
  unname(vapply(x, function(s) {
    token <- tolower(as.character(s))
    hit <- lut[token]
    if (is.na(hit)) as.character(s) else unname(hit)
  }, character(1)))
}

gp_dispatch_legacy_crossvalidation <- function(call_env, parallel_choices, dots = list()) {
  args <- as.list(call_env, all.names = TRUE)

  pheno_data <- args$pheno_data
  args$pheno_data <- NULL
  args[["..."]] <- NULL

  if (length(dots) > 0) {
    args[names(dots)] <- dots
  }

  args$parallel_mode <- match.arg(args$parallel_mode, choices = parallel_choices)
  verbose <- isTRUE(args$verbose)

  models_execute_crossval(
    pheno_data = pheno_data,
    params = args,
    verbose = verbose
  )
}

gp_bglr_save_prefix <- function(model_name = "BGLR", response = "trait", dir = tempdir()) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%OS6")
  stamp <- gsub("[^0-9A-Za-z_]+", "_", stamp)
  token <- paste0(
    model_name, "_",
    response, "_",
    stamp, "_",
    Sys.getpid(), "_",
    sprintf("%06d", sample.int(999999L, 1L))
  )
  file.path(dir, token)
}

# TRUE for no fixed formula or an intercept-only one (~1): both mean "the
# default fixed part", so routes that build their own model accept them.
gp_formula_is_intercept_only <- function(f) {
  if (is.null(f)) return(TRUE)
  if (!inherits(f, "formula")) return(FALSE)
  tt <- stats::terms(f)
  length(attr(tt, "term.labels")) == 0L && attr(tt, "intercept") == 1L
}
