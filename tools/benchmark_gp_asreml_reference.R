#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the PredictProR package root.", call. = FALSE)
}

configured_asreml_library <- Sys.getenv("PREDICTPRO_ASREML_LIB", unset = "")
configured_asreml_library_path <- character()
if (nzchar(configured_asreml_library)) {
  if (!dir.exists(configured_asreml_library)) {
    stop("PREDICTPRO_ASREML_LIB does not name an existing directory.", call. = FALSE)
  }
  configured_asreml_library_path <- normalizePath(
    configured_asreml_library,
    winslash = "/",
    mustWork = TRUE
  )
}
project_library <- file.path(root, ".r-lib")
project_library_path <- if (dir.exists(project_library)) project_library else character()
.libPaths(unique(c(
  configured_asreml_library_path,
  project_library_path,
  .libPaths()
)))
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The pkgload package is required.", call. = FALSE)
}
if (!requireNamespace("asreml", quietly = TRUE)) {
  stop("ASReml-R is not installed for this R interpreter.", call. = FALSE)
}

license_status <- tryCatch(
  asreml::asreml.license.status(),
  error = function(e) stop("ASReml-R license checkout failed: ", conditionMessage(e), call. = FALSE)
)
if (is.null(license_status$status) || is.na(license_status$status) ||
    as.integer(license_status$status) < 0L) {
  stop("ASReml-R license checkout failed.", call. = FALSE)
}

pkgload::load_all(root, quiet = TRUE)

n_genotypes <- as.integer(parse_option("n-genotypes", "36"))
n_holdout <- as.integer(parse_option("n-holdout", "6"))
seed <- as.integer(parse_option("seed", "20260813"))
out_dir <- parse_option(
  "out-dir",
  file.path("tmp", "gp_asreml_reference_benchmark", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
scenario_filter <- trimws(strsplit(parse_option("scenarios", "all"), ",", fixed = TRUE)[[1L]])

if (!is.finite(n_genotypes) || n_genotypes < 24L) {
  stop("`--n-genotypes` must be at least 24 for the covariance models.", call. = FALSE)
}
if (!is.finite(n_holdout) || n_holdout < 3L || n_holdout >= n_genotypes / 2) {
  stop("`--n-holdout` must be at least 3 and less than half of n-genotypes.", call. = FALSE)
}

scenarios <- PredictProR:::gp_asreml_benchmark_scenarios()
if (!identical(scenario_filter, "all")) {
  unknown <- setdiff(scenario_filter, scenarios$scenario)
  if (length(unknown)) {
    stop("Unknown scenario(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  scenarios <- scenarios[scenarios$scenario %in% scenario_filter, , drop = FALSE]
}

configured_gp_python <- Sys.getenv("PREDICTPRO_GP_PYTHON", unset = "")
python_candidates <- unique(c(
  configured_gp_python,
  file.path(root, ".venv-gp", "Scripts", "python.exe"),
  file.path(root, ".venv-gp", "bin", "python"),
  unname(Sys.which("python3")),
  unname(Sys.which("python"))
))
python_candidates <- python_candidates[nzchar(python_candidates)]
python_candidates <- python_candidates[file.exists(python_candidates)]
if (!length(python_candidates)) {
  stop("A configured PredictProR GP Python interpreter was not found.", call. = FALSE)
}
gp_python <- python_candidates[[1L]]

Sys.setenv(
  PREDICTPRO_GP_PYTHON = normalizePath(gp_python, winslash = "/", mustWork = TRUE),
  PREDICTPRO_GP_EXECUTION = "direct",
  PREDICTPRO_GP_RETICULATE_FALLBACK = "false",
  PREDICTPRO_GP_WORKER = "false",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

old_asreml_options <- tryCatch(asreml::asreml.options(), error = function(e) NULL)
try(asreml::asreml.options(ai.sing = TRUE, trace = FALSE), silent = TRUE)
on.exit({
  if (is.list(old_asreml_options) && "ai.sing" %in% names(old_asreml_options)) {
    try(asreml::asreml.options(ai.sing = old_asreml_options$ai.sing), silent = TRUE)
  }
}, add = TRUE)

make_kernel <- function(gids, p, one_seed) {
  set.seed(one_seed)
  x <- matrix(stats::rnorm(length(gids) * p), nrow = length(gids), ncol = p)
  x <- scale(x, center = TRUE, scale = TRUE)
  k <- tcrossprod(x) / ncol(x)
  diag(k) <- diag(k) + 0.08
  k <- k / mean(diag(k))
  dimnames(k) <- list(gids, gids)
  k
}

left_chol <- function(x) {
  t(chol((x + t(x)) / 2 + diag(1e-10, nrow(x))))
}

matrix_normal <- function(row_cov, col_cov) {
  left_chol(row_cov) %*%
    matrix(stats::rnorm(nrow(row_cov) * nrow(col_cov)), nrow(row_cov), nrow(col_cov)) %*%
    t(left_chol(col_cov))
}

named_diagonal <- function(values, labels) {
  out <- diag(as.numeric(values), nrow = length(labels), ncol = length(labels))
  dimnames(out) <- list(labels, labels)
  out
}

truth_key <- function(gid, env, trait) {
  paste(as.character(gid), as.character(env), as.character(trait), sep = "\r")
}

attach_truth <- function(predictions, truth) {
  prediction_key <- truth_key(predictions$GID, predictions$Env, predictions$Trait)
  expected_key <- truth_key(truth$GID, truth$Env, truth$Trait)
  idx <- match(prediction_key, expected_key)
  if (anyNA(idx)) {
    stop("Predictions could not be matched to the simulated truth.", call. = FALSE)
  }
  predictions$Observed_value <- truth$Observed_value[idx]
  predictions
}

standardize_predictions <- function(data, scenario, engine, family, truth,
                                    default_env, default_trait) {
  out <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    data = data,
    scenario = scenario,
    dataset = paste0("synthetic_", scenario, "_seed_", seed),
    engine = engine,
    model_family = family,
    default_env = default_env,
    default_trait = default_trait
  )
  out <- attach_truth(out, truth)
  validation <- PredictProR:::gp_validate_asreml_benchmark_predictions(out)
  if (!isTRUE(validation$valid)) {
    stop(paste(validation$errors, collapse = " "), call. = FALSE)
  }
  out
}

elapsed_call <- function(code) {
  started <- proc.time()[["elapsed"]]
  value <- force(code)
  list(value = value, elapsed_seconds = proc.time()[["elapsed"]] - started)
}

asreml_varcomp <- function(model) {
  out <- as.data.frame(
    asreml::summary.asreml(model, param = "sigma")$varcomp,
    stringsAsFactors = FALSE
  )
  if (!"component" %in% names(out)) {
    stop("ASReml variance-component table has no `component` column.", call. = FALSE)
  }
  out$component <- suppressWarnings(as.numeric(out$component))
  out
}

asreml_scalar_component <- function(model, pattern, fixed = TRUE) {
  vc <- asreml_varcomp(model)
  label <- rownames(vc)
  hit <- if (isTRUE(fixed)) grepl(pattern, label, fixed = TRUE) else grepl(pattern, label)
  values <- vc$component[hit & is.finite(vc$component)]
  if (length(values) != 1L) {
    stop("Expected one ASReml component for pattern `", pattern, "`; found ", length(values), ".", call. = FALSE)
  }
  values[[1L]]
}

asreml_us_covariance <- function(model, traits, mask) {
  vc <- asreml_varcomp(model)
  labels <- rownames(vc)
  pairs <- PredictProR:::gp_multitrait_asreml_component_trait_pair(labels, response = traits)
  matrix_out <- matrix(NA_real_, length(traits), length(traits), dimnames = list(traits, traits))
  for (i in which(mask(labels))) {
    pair <- pairs[[i]]
    if (length(pair) != 2L || !is.finite(vc$component[[i]])) next
    matrix_out[pair[[1L]], pair[[2L]]] <- vc$component[[i]]
    matrix_out[pair[[2L]], pair[[1L]]] <- vc$component[[i]]
  }
  if (anyNA(matrix_out)) {
    stop("ASReml unstructured covariance could not be reconstructed from all trait pairs.", call. = FALSE)
  }
  matrix_out
}

asreml_corgh_covariance <- function(model, traits, mask) {
  vc <- asreml_varcomp(model)
  labels <- rownames(vc)
  pairs <- PredictProR:::gp_multitrait_asreml_component_trait_pair(labels, response = traits)
  variances <- stats::setNames(rep(NA_real_, length(traits)), traits)
  correlation <- named_diagonal(rep(1, length(traits)), traits)
  for (i in which(mask(labels))) {
    pair <- pairs[[i]]
    value <- vc$component[[i]]
    if (length(pair) != 2L || !is.finite(value)) next
    if (identical(pair[[1L]], pair[[2L]])) {
      variances[[pair[[1L]]]] <- value
    } else if (grepl("cor", labels[[i]], ignore.case = TRUE)) {
      correlation[pair[[1L]], pair[[2L]]] <- value
      correlation[pair[[2L]], pair[[1L]]] <- value
    }
  }
  if (any(!is.finite(variances)) || any(variances <= 0)) {
    stop("ASReml CORGH variances could not be reconstructed for every trait.", call. = FALSE)
  }
  scale <- diag(sqrt(variances), nrow = length(traits))
  out <- scale %*% correlation %*% scale
  dimnames(out) <- list(traits, traits)
  out
}

assert_asreml_converged <- function(model, context) {
  if (!isTRUE(model$converge)) {
    stop("ASReml did not converge for ", context, ".", call. = FALSE)
  }
  invisible(TRUE)
}

asreml_ifault <- function(model) {
  if (is.null(model$ifault) || !length(model$ifault)) return(NA_integer_)
  suppressWarnings(as.integer(model$ifault[[1L]]))
}

update_asreml <- function(model, context) {
  assert_asreml_converged(model, context)
  model
}

asreml_pvals <- function(model, classify) {
  pred <- asreml::predict.asreml(model, classify = classify, sed = FALSE)
  out <- as.data.frame(pred$pvals, stringsAsFactors = FALSE)
  if (!nrow(out)) stop("ASReml returned an empty prediction table.", call. = FALSE)
  out
}

assign_inverse <- function(kernel, stem) {
  name <- paste0(stem, "_", Sys.getpid(), "_", sample.int(1e7, 1L))
  assign(
    name,
    PredictProR::compute_inverse_and_sparse(kernel, epsilon = 1e-6, inverse = TRUE),
    envir = .GlobalEnv
  )
  name
}

covariance_rows <- function(scenario, engine, component, estimate, truth) {
  estimate <- as.matrix(estimate)
  truth <- as.matrix(truth)
  diagnostics <- PredictProR:::gp_asreml_benchmark_covariance_diagnostics(estimate, truth)
  if (!isTRUE(diagnostics$positive_semidefinite) || diagnostics$symmetry_error > 1e-7) {
    stop(engine, " returned an invalid covariance for ", scenario, "/", component, ".", call. = FALSE)
  }
  cells <- expand.grid(
    row = rownames(estimate) %||% as.character(seq_len(nrow(estimate))),
    column = colnames(estimate) %||% as.character(seq_len(ncol(estimate))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  cells$scenario <- scenario
  cells$dataset <- paste0("synthetic_", scenario, "_seed_", seed)
  cells$engine <- engine
  cells$component <- component
  cells$estimate <- as.numeric(estimate)
  cells$truth <- as.numeric(truth)
  cells <- cells[c("scenario", "dataset", "engine", "component", "row", "column", "estimate", "truth")]
  diagnostics$scenario <- scenario
  diagnostics$dataset <- paste0("synthetic_", scenario, "_seed_", seed)
  diagnostics$engine <- engine
  diagnostics$component <- component
  diagnostics <- diagnostics[c("scenario", "dataset", "engine", "component", setdiff(names(diagnostics), c("scenario", "dataset", "engine", "component")))]
  list(cells = cells, diagnostics = diagnostics)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

gp_varcomp_scalar <- function(output, genetic = TRUE) {
  vc <- as.data.frame(output$result$varcomp, stringsAsFactors = FALSE)
  component_col <- if ("component" %in% names(vc)) "component" else "Component"
  estimate_col <- if ("estimate" %in% names(vc)) "estimate" else if ("Components" %in% names(vc)) "Components" else NULL
  if (is.null(estimate_col)) stop("GP variance-component estimate column was not found.", call. = FALSE)
  labels <- tolower(as.character(vc[[component_col]]))
  mask <- if (isTRUE(genetic)) grepl("vm\\(|genetic", labels) & !grepl("!r|resid", labels) else grepl("!r|resid", labels)
  values <- suppressWarnings(as.numeric(vc[[estimate_col]][mask]))
  values <- values[is.finite(values)]
  if (length(values) != 1L) stop("Expected one GP variance component; found ", length(values), ".", call. = FALSE)
  values[[1L]]
}

make_base_data <- function() {
  gids <- paste0("g", seq_len(n_genotypes))
  envs <- paste0("E", seq_len(3L))
  traits <- c("Trait1", "Trait2")
  list(
    gids = gids,
    envs = envs,
    traits = traits,
    holdout = tail(gids, n_holdout),
    k_st = make_kernel(gids, max(10L, n_genotypes %/% 3L), seed + 11L),
    k1 = make_kernel(gids, max(48L, 2L * n_genotypes), seed + 11L),
    k2 = make_kernel(gids, max(52L, 2L * n_genotypes + 7L), seed + 23L)
  )
}

base <- make_base_data()

run_single_trait_single_environment <- function() {
  scenario <- "single_trait_single_environment"
  set.seed(seed + 101L)
  sigma_g <- 1.15
  sigma_e <- 0.45
  u <- as.numeric(matrix_normal(base$k_st, matrix(sigma_g, 1L)))
  pheno_truth <- data.frame(
    GID = base$gids,
    Yield = 4.0 + u + stats::rnorm(n_genotypes, sd = sqrt(sigma_e)),
    stringsAsFactors = FALSE
  )
  truth <- data.frame(
    GID = pheno_truth$GID,
    Env = "E1",
    Trait = "Yield",
    Observed_value = pheno_truth$Yield,
    stringsAsFactors = FALSE
  )
  masked <- pheno_truth
  masked$Yield[masked$GID %in% base$holdout] <- NA_real_

  gp_timed <- elapsed_call(PredictProR::gp_single_trait_model(
    pheno_data = masked,
    gmatrix = base$k_st,
    response = "Yield",
    gen_name = "GID",
    test_set = base$holdout,
    model_name = "GP",
    gp_backend = "torch",
    gp_output_level = "full_vc",
    gp_return_se = TRUE,
    gp_full_vc = TRUE,
    gp_varcomp_mode = "reml",
    gp_iters = 60L,
    mean_adjustment = "none",
    prediction_blend = "none",
    gp_uncertainty_calibration = "none",
    tune_gp = FALSE,
    python_bin = gp_python,
    seed = seed
  ))
  gp_out <- gp_timed$value
  gp_out$predictions$Env <- "E1"
  gp_pred <- standardize_predictions(
    gp_out$predictions, scenario, "PredictProR_GP", "GP_REML", truth, "E1", "Yield"
  )

  inv_name <- assign_inverse(base$k_st, "ppr_bench_st_inv")
  on.exit(PredictProR::remove_from_global(inv_name), add = TRUE)
  asr_timed <- elapsed_call(do.call(asreml::asreml, list(
    fixed = Yield ~ 1,
    random = stats::as.formula(paste0("~ vm(GID,", inv_name, ")")),
    residual = ~ units,
    data = transform(masked, GID = factor(GID, levels = base$gids)),
    na.action = list(x = "include", y = "include"),
    maxit = 100L,
    workspace = 1e8,
    trace = FALSE,
    ai.sing = TRUE
  )))
  asr_model <- update_asreml(asr_timed$value, scenario)
  asr_raw <- asreml_pvals(asr_model, "GID")
  asr_raw <- asr_raw[as.character(asr_raw$GID) %in% base$holdout, , drop = FALSE]
  asr_pred <- standardize_predictions(
    asr_raw, scenario, "ASReml_R", "GBLUP_REML", truth, "E1", "Yield"
  )

  gp_g <- matrix(gp_varcomp_scalar(gp_out, TRUE), 1L, dimnames = list("Yield", "Yield"))
  gp_e <- matrix(gp_varcomp_scalar(gp_out, FALSE), 1L, dimnames = list("Yield", "Yield"))
  asr_g <- matrix(asreml_scalar_component(asr_model, inv_name), 1L, dimnames = list("Yield", "Yield"))
  asr_e <- matrix(asreml_scalar_component(asr_model, "!R", fixed = FALSE), 1L, dimnames = list("Yield", "Yield"))

  list(
    predictions = rbind(gp_pred, asr_pred),
    timings = data.frame(scenario, engine = c("PredictProR_GP", "ASReml_R"), elapsed_seconds = c(gp_timed$elapsed_seconds, asr_timed$elapsed_seconds)),
    convergence = data.frame(
      scenario,
      engine = c("PredictProR_GP", "ASReml_R"),
      converged = c(TRUE, isTRUE(asr_model$converge)),
      ifault = c(NA_integer_, asreml_ifault(asr_model))
    ),
    covariance = list(
      covariance_rows(scenario, "PredictProR_GP", "genetic", gp_g, matrix(sigma_g, 1L, dimnames = dimnames(gp_g))),
      covariance_rows(scenario, "PredictProR_GP", "residual", gp_e, matrix(sigma_e, 1L, dimnames = dimnames(gp_e))),
      covariance_rows(scenario, "ASReml_R", "genetic", asr_g, matrix(sigma_g, 1L, dimnames = dimnames(asr_g))),
      covariance_rows(scenario, "ASReml_R", "residual", asr_e, matrix(sigma_e, 1L, dimnames = dimnames(asr_e)))
    )
  )
}

run_single_trait_met <- function(multi_kernel = FALSE) {
  scenario <- if (isTRUE(multi_kernel)) "single_trait_met_multi_kernel" else "single_trait_met_single_kernel"
  set.seed(seed + if (multi_kernel) 203L else 201L)
  loadings <- c(E1 = 0.95, E2 = 0.68, E3 = 0.42)
  uniqueness <- c(E1 = 0.20, E2 = 0.26, E3 = 0.31)
  sigma_env <- tcrossprod(loadings) + diag(uniqueness)
  dimnames(sigma_env) <- list(base$envs, base$envs)
  sigma_resid <- c(E1 = 0.34, E2 = 0.48, E3 = 0.39)
  genetic <- if (isTRUE(multi_kernel)) {
    matrix_normal(0.70 * base$k1, sigma_env) + matrix_normal(0.30 * base$k2, sigma_env)
  } else {
    matrix_normal(base$k1, sigma_env)
  }
  pheno_truth <- expand.grid(
    GID = base$gids,
    Env = base$envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  gi <- match(pheno_truth$GID, base$gids)
  ei <- match(pheno_truth$Env, base$envs)
  pheno_truth$Yield <- c(E1 = 4.0, E2 = 4.5, E3 = 3.7)[pheno_truth$Env] +
    genetic[cbind(gi, ei)] + stats::rnorm(nrow(pheno_truth), sd = sqrt(sigma_resid[ei]))
  truth <- transform(
    pheno_truth[c("GID", "Env")],
    Trait = "Yield",
    Observed_value = pheno_truth$Yield
  )
  masked <- pheno_truth
  masked$Yield[masked$GID %in% base$holdout] <- NA_real_

  gp_timed <- elapsed_call(PredictProR::gp_single_trait_model(
    pheno_data = masked,
    gmatrix = base$k1,
    kernel_list = if (isTRUE(multi_kernel)) list(omic = base$k2) else NULL,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    test_set = base$holdout,
    fixed_effects = "Env",
    model_name = "GP_FA",
    gp_backend = "torch",
    gp_output_level = "full_vc",
    gp_return_se = TRUE,
    gp_full_vc = TRUE,
    gp_varcomp_mode = "reml",
    gp_fa_rank = 1L,
    gp_iters = 80L,
    gp_learn_scales = TRUE,
    w_g = if (isTRUE(multi_kernel)) c(0.70, 0.30) else 1,
    w_ge = if (isTRUE(multi_kernel)) c(0.70, 0.30) else 1,
    mean_adjustment = "none",
    prediction_blend = "none",
    gp_uncertainty_calibration = "none",
    tune_gp = FALSE,
    python_bin = gp_python,
    seed = seed
  ))
  gp_out <- gp_timed$value
  gp_pred <- standardize_predictions(
    gp_out$predictions, scenario, "PredictProR_GP", "GP_FA1_REML", truth, NA_character_, "Yield"
  )

  asr_timed <- elapsed_call(PredictProR:::asreml_utilis_new(
    fixed = ~ Env,
    random = ~ GID + GID:Env,
    GS_model = "GBLUP",
    response = "Yield",
    pheno_data = transform(masked, GID = factor(GID, levels = base$gids), Env = factor(Env, levels = base$envs)),
    gmatrix = base$k1,
    kernel_list = if (isTRUE(multi_kernel)) list(omic = base$k2) else NULL,
    inverse = TRUE,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "fa1",
    engine = "asreml",
    workspace = 1e8,
    pworkspace = 1e8,
    maxit = 100L
  ))
  asr_fit <- asr_timed$value
  on.exit(PredictProR::remove_from_global(asr_fit$names_in_inv_list), add = TRUE)
  asr_model <- update_asreml(asr_fit$model, scenario)
  asr_raw <- asreml_pvals(asr_model, "Env:GID")
  asr_raw <- asr_raw[as.character(asr_raw$GID) %in% base$holdout, , drop = FALSE]
  asr_pred <- standardize_predictions(
    asr_raw, scenario, "ASReml_R", "FA1_GBLUP_REML", truth, NA_character_, "Yield"
  )

  asr_covariance <- PredictProR::asreml_herit_varCov_new(
    model = asr_model,
    heter_groups = "Env",
    var_cov_str = "fa1",
    heter_resid = TRUE,
    names_in_inv_list = asr_fit$names_in_inv_list,
    inter_gen_pos = asr_fit$inter_gen_pos,
    gen_pos = asr_fit$gen_pos
  )
  gp_g <- as.matrix(gp_out$result$env_covariance_fa)
  dimnames(gp_g) <- list(base$envs, base$envs)
  gp_e <- named_diagonal(gp_out$result$sigma2_resid_env, base$envs)
  asr_g <- Reduce(`+`, lapply(asr_covariance$Covariance, as.matrix))
  dimnames(asr_g) <- list(base$envs, base$envs)
  asr_e <- named_diagonal(asr_covariance$Residual_Var, base$envs)
  truth_e <- named_diagonal(sigma_resid, base$envs)

  list(
    predictions = rbind(gp_pred, asr_pred),
    timings = data.frame(scenario, engine = c("PredictProR_GP", "ASReml_R"), elapsed_seconds = c(gp_timed$elapsed_seconds, asr_timed$elapsed_seconds)),
    convergence = data.frame(
      scenario,
      engine = c("PredictProR_GP", "ASReml_R"),
      converged = c(TRUE, isTRUE(asr_model$converge)),
      ifault = c(NA_integer_, asreml_ifault(asr_model))
    ),
    covariance = list(
      covariance_rows(scenario, "PredictProR_GP", "genetic_environment", gp_g, sigma_env),
      covariance_rows(scenario, "PredictProR_GP", "residual_environment", gp_e, truth_e),
      covariance_rows(scenario, "ASReml_R", "genetic_environment", asr_g, sigma_env),
      covariance_rows(scenario, "ASReml_R", "residual_environment", asr_e, truth_e)
    )
  )
}

run_multi_trait_single_environment <- function() {
  scenario <- "multi_trait_single_environment"
  set.seed(seed + 301L)
  sigma_g <- matrix(c(1.00, 0.55, 0.55, 0.75), 2L, dimnames = list(base$traits, base$traits))
  sigma_e <- matrix(c(0.42, 0.12, 0.12, 0.36), 2L, dimnames = list(base$traits, base$traits))
  genetic <- matrix_normal(base$k1, sigma_g)
  residual <- matrix_normal(diag(n_genotypes), sigma_e)
  values <- sweep(genetic + residual, 2L, c(4.2, 2.7), "+")
  pheno_truth <- data.frame(GID = base$gids, Trait1 = values[, 1L], Trait2 = values[, 2L])
  truth <- rbind(
    data.frame(GID = base$gids, Env = "E1", Trait = "Trait1", Observed_value = values[, 1L]),
    data.frame(GID = base$gids, Env = "E1", Trait = "Trait2", Observed_value = values[, 2L])
  )
  masked <- pheno_truth
  masked[masked$GID %in% base$holdout, base$traits] <- NA_real_

  gp_timed <- elapsed_call(PredictProR::gp_multi_trait_model(
    pheno_data = masked,
    gmatrix = base$k1,
    response = base$traits,
    gen_name = "GID",
    test_set = base$holdout,
    varcomp_mode = "reml",
    tune_gp = FALSE,
    prediction_blend = "none",
    gp_uncertainty_calibration = "none",
    return_se = TRUE,
    return_trait_correlations = TRUE,
    max_iter = 45L,
    tol_loglik = 1e-5,
    seed = seed,
    python_bin = gp_python
  ))
  gp_out <- gp_timed$value
  gp_pred <- standardize_predictions(
    gp_out$predictions, scenario, "PredictProR_GP", "MT_GP_US_REML", truth, "E1", NA_character_
  )

  asr_timed <- elapsed_call(PredictProR:::gp_multitrait_asreml_gaussian_model(
    pheno_object = masked,
    response = base$traits,
    gen_name = "GID",
    gmatrix = base$k1,
    inverse = TRUE,
    engine = "asreml",
    workspace = 1e8,
    maxit = 100L,
    ai_sing = TRUE,
    trait_covariance = "us"
  ))
  asr_out <- asr_timed$value
  asr_model <- asr_out$Asreml_model
  assert_asreml_converged(asr_model, scenario)
  asr_raw <- asr_out$predicted_values
  asr_raw <- asr_raw[asr_raw$Train_Test_Label == "Test", , drop = FALSE]
  asr_pred <- standardize_predictions(
    asr_raw, scenario, "ASReml_R", "MT_US_GBLUP_REML", truth, "E1", NA_character_
  )

  gp_g <- as.matrix(gp_out$result$genetic_covariance)
  gp_e <- as.matrix(gp_out$result$residual_covariance)
  dimnames(gp_g) <- dimnames(gp_e) <- list(base$traits, base$traits)
  asr_g <- as.matrix(asr_out$Genetic_covariance_traits)
  asr_e <- as.matrix(asr_out$Residual_covariance_traits)

  list(
    predictions = rbind(gp_pred, asr_pred),
    timings = data.frame(scenario, engine = c("PredictProR_GP", "ASReml_R"), elapsed_seconds = c(gp_timed$elapsed_seconds, asr_timed$elapsed_seconds)),
    convergence = data.frame(
      scenario,
      engine = c("PredictProR_GP", "ASReml_R"),
      converged = c(TRUE, isTRUE(asr_model$converge)),
      ifault = c(NA_integer_, asreml_ifault(asr_model))
    ),
    covariance = list(
      covariance_rows(scenario, "PredictProR_GP", "genetic_traits", gp_g, sigma_g),
      covariance_rows(scenario, "PredictProR_GP", "residual_traits", gp_e, sigma_e),
      covariance_rows(scenario, "ASReml_R", "genetic_traits", asr_g, sigma_g),
      covariance_rows(scenario, "ASReml_R", "residual_traits", asr_e, sigma_e)
    )
  )
}

run_multi_trait_multi_environment <- function() {
  scenario <- "multi_trait_multi_environment"
  set.seed(seed + 401L)
  sigma_g <- matrix(c(0.72, 0.35, 0.35, 0.55), 2L, dimnames = list(base$traits, base$traits))
  sigma_ge <- matrix(c(0.42, 0.18, 0.18, 0.34), 2L, dimnames = list(base$traits, base$traits))
  sigma_e <- matrix(c(0.34, 0.08, 0.08, 0.29), 2L, dimnames = list(base$traits, base$traits))
  main <- matrix_normal(base$k1, sigma_g)
  gxe <- lapply(base$envs, function(x) matrix_normal(base$k1, sigma_ge))
  pheno_truth <- expand.grid(
    GID = base$gids,
    Env = base$envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  values <- matrix(NA_real_, nrow(pheno_truth), 2L, dimnames = list(NULL, base$traits))
  env_means <- rbind(E1 = c(4.0, 2.5), E2 = c(4.5, 2.2), E3 = c(3.8, 2.9))
  for (i in seq_len(nrow(pheno_truth))) {
    g <- match(pheno_truth$GID[[i]], base$gids)
    e <- match(pheno_truth$Env[[i]], base$envs)
    residual <- as.numeric(matrix_normal(matrix(1, 1L), sigma_e))
    values[i, ] <- env_means[e, ] + main[g, ] + gxe[[e]][g, ] + residual
  }
  pheno_truth[base$traits] <- values
  truth <- do.call(rbind, lapply(base$traits, function(trait) {
    data.frame(
      GID = pheno_truth$GID,
      Env = pheno_truth$Env,
      Trait = trait,
      Observed_value = pheno_truth[[trait]],
      stringsAsFactors = FALSE
    )
  }))
  masked <- pheno_truth
  masked[masked$GID %in% base$holdout, base$traits] <- NA_real_
  env_identity <- named_diagonal(rep(1, length(base$envs)), base$envs)

  gp_timed <- elapsed_call(PredictProR::gp_multi_trait_met_model(
    pheno_data = masked,
    gmatrix = base$k1,
    response = base$traits,
    gen_name = "GID",
    heter_groups = "Env",
    test_set = base$holdout,
    env_similarity = env_identity,
    varcomp_mode = "reml",
    tune_gp = FALSE,
    prediction_blend = "none",
    gp_uncertainty_calibration = "none",
    return_se = TRUE,
    return_trait_correlations = TRUE,
    trait_structure = "unstructured",
    gxe_trait_structure = "unstructured",
    dense_reml_max_train = 3000L,
    reml_max_iter = 30L,
    reml_tol = 1e-5,
    seed = seed,
    python_bin = gp_python
  ))
  gp_out <- gp_timed$value
  gp_pred <- standardize_predictions(
    gp_out$predictions, scenario, "PredictProR_GP", "MT_MET_GP_US_REML", truth, NA_character_, NA_character_
  )

  keys <- as.vector(t(outer(base$gids, base$envs, paste, sep = "__")))
  kge <- kronecker(base$k1, diag(length(base$envs)))
  dimnames(kge) <- list(keys, keys)
  inv_g <- assign_inverse(base$k1, "ppr_bench_mtmet_g_inv")
  inv_ge <- assign_inverse(kge, "ppr_bench_mtmet_ge_inv")
  on.exit(PredictProR::remove_from_global(c(inv_g, inv_ge)), add = TRUE)

  long <- do.call(rbind, lapply(base$traits, function(trait) {
    data.frame(
      GID = masked$GID,
      Env = masked$Env,
      Trait = trait,
      Observed_value = masked[[trait]],
      stringsAsFactors = FALSE
    )
  }))
  long$GID <- factor(long$GID, levels = base$gids)
  long$Env <- factor(long$Env, levels = base$envs)
  long$Trait <- factor(long$Trait, levels = base$traits)
  long$GIDEnv <- factor(paste(long$GID, long$Env, sep = "__"), levels = keys)
  long <- long[order(long$GIDEnv, long$Trait), , drop = FALSE]
  rownames(long) <- NULL

  asr_timed <- elapsed_call(do.call(asreml::asreml, list(
    fixed = Observed_value ~ Trait:Env,
    random = stats::as.formula(paste0(
      "~ corgh(Trait):vm(GID,", inv_g, ") + corgh(Trait):vm(GIDEnv,", inv_ge, ")"
    )),
    residual = ~ units:corgh(Trait),
    data = long,
    asmv = "Trait",
    na.action = list(x = "include", y = "include"),
    maxit = 150L,
    workspace = 2e8,
    trace = FALSE,
    ai.sing = TRUE
  )))
  asr_model <- update_asreml(asr_timed$value, scenario)
  asr_raw <- asreml_pvals(asr_model, "Trait:Env:GID")
  asr_raw <- asr_raw[as.character(asr_raw$GID) %in% base$holdout, , drop = FALSE]
  asr_pred <- standardize_predictions(
    asr_raw, scenario, "ASReml_R", "MT_MET_US_GBLUP_REML", truth, NA_character_, NA_character_
  )

  gp_g <- as.matrix(gp_out$result$genetic_covariance)
  gp_ge <- as.matrix(gp_out$result$gxe_covariance)
  gp_e <- as.matrix(gp_out$result$residual_covariance)
  dimnames(gp_g) <- dimnames(gp_ge) <- dimnames(gp_e) <- list(base$traits, base$traits)
  asr_g <- asreml_corgh_covariance(asr_model, base$traits, function(labels) grepl(inv_g, labels, fixed = TRUE))
  asr_ge <- asreml_corgh_covariance(asr_model, base$traits, function(labels) grepl(inv_ge, labels, fixed = TRUE))
  asr_e <- asreml_corgh_covariance(asr_model, base$traits, function(labels) grepl("units", labels, ignore.case = TRUE))

  list(
    predictions = rbind(gp_pred, asr_pred),
    timings = data.frame(scenario, engine = c("PredictProR_GP", "ASReml_R"), elapsed_seconds = c(gp_timed$elapsed_seconds, asr_timed$elapsed_seconds)),
    convergence = data.frame(
      scenario,
      engine = c("PredictProR_GP", "ASReml_R"),
      converged = c(TRUE, isTRUE(asr_model$converge)),
      ifault = c(NA_integer_, asreml_ifault(asr_model))
    ),
    covariance = list(
      covariance_rows(scenario, "PredictProR_GP", "genetic_main_traits", gp_g, sigma_g),
      covariance_rows(scenario, "PredictProR_GP", "gxe_traits", gp_ge, sigma_ge),
      covariance_rows(scenario, "PredictProR_GP", "residual_traits", gp_e, sigma_e),
      covariance_rows(scenario, "ASReml_R", "genetic_main_traits", asr_g, sigma_g),
      covariance_rows(scenario, "ASReml_R", "gxe_traits", asr_ge, sigma_ge),
      covariance_rows(scenario, "ASReml_R", "residual_traits", asr_e, sigma_e)
    )
  )
}

runner <- list(
  single_trait_single_environment = run_single_trait_single_environment,
  single_trait_met_single_kernel = function() run_single_trait_met(FALSE),
  single_trait_met_multi_kernel = function() run_single_trait_met(TRUE),
  multi_trait_single_environment = run_multi_trait_single_environment,
  multi_trait_multi_environment = run_multi_trait_multi_environment
)

results <- list()
for (scenario in scenarios$scenario) {
  cat("RUNNING ", scenario, "\n", sep = "")
  results[[scenario]] <- runner[[scenario]]()
  saveRDS(results, file.path(out_dir, "checkpoint.rds"))
}

predictions <- do.call(rbind, lapply(results, `[[`, "predictions"))
rownames(predictions) <- NULL
validation <- PredictProR:::gp_validate_asreml_benchmark_predictions(predictions)
if (!isTRUE(validation$valid)) stop(paste(validation$errors, collapse = " "), call. = FALSE)

accuracy_tables <- PredictProR:::gp_asreml_benchmark_metric_tables(predictions, scenarios)
accuracy <- accuracy_tables$overall
accuracy_by_environment <- accuracy_tables$by_environment
accuracy_by_trait <- accuracy_tables$by_trait
accuracy_by_environment_trait <- accuracy_tables$by_environment_trait
accuracy_stratified <- accuracy_tables$stratified
timings <- do.call(rbind, lapply(results, `[[`, "timings"))
convergence <- do.call(rbind, lapply(results, `[[`, "convergence"))
covariance_objects <- unlist(lapply(results, `[[`, "covariance"), recursive = FALSE)
covariance_estimates <- do.call(rbind, lapply(covariance_objects, `[[`, "cells"))
covariance_diagnostics <- do.call(rbind, lapply(covariance_objects, `[[`, "diagnostics"))

prediction_parity <- do.call(rbind, lapply(split(predictions, predictions$scenario), function(one) {
  gp <- one[one$engine == "PredictProR_GP", c("GID", "Env", "Trait", "Predicted_value")]
  asr <- one[one$engine == "ASReml_R", c("GID", "Env", "Trait", "Predicted_value")]
  joined <- merge(gp, asr, by = c("GID", "Env", "Trait"), suffixes = c("_gp", "_asreml"))
  data.frame(
    scenario = one$scenario[[1L]],
    n = nrow(joined),
    prediction_RMSE = sqrt(mean((joined$Predicted_value_gp - joined$Predicted_value_asreml)^2)),
    prediction_Pearson = if (nrow(joined) >= 3L) stats::cor(joined$Predicted_value_gp, joined$Predicted_value_asreml) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(predictions, file.path(out_dir, "predictions.csv"), row.names = FALSE, na = "")
utils::write.csv(accuracy, file.path(out_dir, "accuracy_metrics.csv"), row.names = FALSE, na = "")
utils::write.csv(accuracy_by_environment, file.path(out_dir, "accuracy_by_environment.csv"), row.names = FALSE, na = "")
utils::write.csv(accuracy_by_trait, file.path(out_dir, "accuracy_by_trait.csv"), row.names = FALSE, na = "")
utils::write.csv(accuracy_by_environment_trait, file.path(out_dir, "accuracy_by_environment_trait.csv"), row.names = FALSE, na = "")
utils::write.csv(accuracy_stratified, file.path(out_dir, "accuracy_metrics_stratified.csv"), row.names = FALSE, na = "")
utils::write.csv(prediction_parity, file.path(out_dir, "engine_prediction_parity.csv"), row.names = FALSE, na = "")
utils::write.csv(timings, file.path(out_dir, "timings.csv"), row.names = FALSE, na = "")
utils::write.csv(convergence, file.path(out_dir, "convergence.csv"), row.names = FALSE, na = "")
utils::write.csv(covariance_estimates, file.path(out_dir, "covariance_estimates.csv"), row.names = FALSE, na = "")
utils::write.csv(covariance_diagnostics, file.path(out_dir, "covariance_diagnostics.csv"), row.names = FALSE, na = "")
utils::write.csv(scenarios, file.path(out_dir, "scenario_contract.csv"), row.names = FALSE, na = "")

metadata <- list(
  timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  PredictProR_version = as.character(utils::packageVersion("PredictProR")),
  R_version = R.version.string,
  asreml_version = as.character(utils::packageVersion("asreml")),
  asreml_library = normalizePath(find.package("asreml"), winslash = "/", mustWork = TRUE),
  python = normalizePath(gp_python, winslash = "/", mustWork = TRUE),
  n_genotypes = n_genotypes,
  n_holdout = n_holdout,
  seed = seed,
  sequential = TRUE,
  gp_worker_reuse = FALSE,
  scenario_contract = scenarios,
  predictions = predictions,
  accuracy = accuracy,
  accuracy_by_environment = accuracy_by_environment,
  accuracy_by_trait = accuracy_by_trait,
  accuracy_by_environment_trait = accuracy_by_environment_trait,
  accuracy_stratified = accuracy_stratified,
  prediction_parity = prediction_parity,
  timings = timings,
  convergence = convergence,
  covariance_diagnostics = covariance_diagnostics
)
saveRDS(metadata, file.path(out_dir, "benchmark_bundle.rds"))

report <- c(
  "# PredictProR GP and ASReml-R reference benchmark",
  "",
  paste0("- Dataset: independently simulated scenario data with seed ", seed),
  paste0("- Genotypes: ", n_genotypes),
  paste0("- Held-out genotypes: ", n_holdout),
  "- Execution: sequential; GP worker reuse disabled",
  "- Interpretation: correctness/reference run, not production throughput ranking",
  "",
  "All requested fits converged, prediction contracts validated, and reported covariance matrices passed symmetry and positive-semidefinite checks.",
  "Overall prediction metrics are in accuracy_metrics.csv.",
  "Environment, trait, and environment-by-trait prediction metrics are in the corresponding accuracy_by_*.csv files and the combined accuracy_metrics_stratified.csv file.",
  "See engine_prediction_parity.csv, timings.csv, and covariance_diagnostics.csv for the remaining numeric results.",
  "The multi-kernel ASReml FA model is a nested, more flexible reference; its individual FA parameters are not one-to-one with the GP shared-FA kernel weights."
)
writeLines(report, file.path(out_dir, "benchmark_report.md"), useBytes = TRUE)

cat("BENCHMARK_OK\n")
cat("OUT_DIR ", out_dir, "\n", sep = "")
print(accuracy, row.names = FALSE)
print(accuracy_by_environment, row.names = FALSE)
print(accuracy_by_trait, row.names = FALSE)
print(accuracy_by_environment_trait, row.names = FALSE)
print(prediction_parity, row.names = FALSE)
