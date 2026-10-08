test_that("MT-MET GP defaults stay on scalable unstructured covariance", {
  defaults <- formals(PredictProR::gp_multi_trait_met_model)

  expect_identical(defaults$varcomp_mode, "reml")
  expect_identical(defaults$trait_structure, "unstructured")
  expect_null(defaults$gxe_trait_structure)
  expect_identical(defaults$reaction_norm_feature_qc, TRUE)
  expect_identical(defaults$kenv_kernel, "matern32")
  expect_identical(defaults$prediction_output, "test_only")
})

test_that("single-trait MET exposes env covariate kernel routing", {
  defaults <- formals(PredictProR::gp_single_trait_model)

  expect_true("env_covariates" %in% names(defaults))
  expect_identical(defaults$reaction_norm_feature_qc, TRUE)
  expect_identical(defaults$kenv_kernel, "matern32")
  expect_identical(defaults$tune_gp, FALSE)
  expect_identical(defaults$tuning_strategy, "auto")
  expect_identical(defaults$tuning_objective, "rmse")
  expect_identical(defaults$mean_adjustment, "auto")
  expect_identical(defaults$prediction_blend, "auto")
  expect_identical(defaults$gp_auto_policy, TRUE)
  expect_true("gp_factor_cache" %in% names(defaults))
  expect_true("gp_tiered_dispatch" %in% names(defaults))
  expect_true("gp_iters" %in% names(defaults))
  expect_true("gp_lr" %in% names(defaults))
  expect_null(defaults$gp_iters)
  expect_null(defaults$gp_lr)
})

test_that("GP backend schema constants document internal column names", {
  cols <- PredictProR:::gp_backend_schema_cols()

  if (!identical(cols$single$gid, "GID")) {
    stop("single-trait genotype column must remain GID", call. = FALSE)
  }
  if (!identical(cols$single$env, "Env")) {
    stop("single-trait environment column must remain Env", call. = FALSE)
  }
  if (!identical(cols$single$y, "y")) {
    stop("single-trait response column must remain y", call. = FALSE)
  }
  if (!identical(cols$single$value, "Trait")) {
    stop("legacy bridge response column must remain Trait", call. = FALSE)
  }
  if (!identical(cols$long, list(gid = "gid", env = "env", trait = "trait", y = "y"))) {
    stop("multi-trait long schema changed unexpectedly", call. = FALSE)
  }
  expect_true(TRUE)
})

test_that("GP auto policy profile respects caller phenotype column names", {
  pheno <- data.frame(
    sample_id = c("g1", "g2", "g3", "g4"),
    location = c("E1", "E1", "E2", "E3"),
    stringsAsFactors = FALSE
  )
  split <- list(
    train_mask = c(TRUE, TRUE, FALSE, FALSE),
    test_mask = c(FALSE, FALSE, TRUE, TRUE)
  )

  profile <- PredictProR:::gp_policy_prediction_profile(
    pheno,
    split,
    heter_groups = "location",
    gid_col = "sample_id",
    env_col = "location"
  )

  if (!identical(profile$n_env, 3L)) {
    stop("policy profile did not count noncanonical environments", call. = FALSE)
  }
  if (!identical(profile$n_train_env, 1L) || !identical(profile$n_test_env, 2L)) {
    stop("policy profile train/test environment counts are wrong", call. = FALSE)
  }
  if (!identical(profile$n_genotypes, 4L)) {
    stop("policy profile did not count noncanonical genotypes", call. = FALSE)
  }
  if (!isTRUE(profile$is_met) || !isTRUE(profile$prospective_env_prediction)) {
    stop("policy profile did not detect prospective MET prediction", call. = FALSE)
  }
  if (!identical(profile$unseen_test_envs, c("E2", "E3"))) {
    stop("policy profile unseen test environments are wrong", call. = FALSE)
  }
  expect_true(TRUE)
})

test_that("GP exact optimizer controls are preserved in bridge specs", {
  ns <- asNamespace("PredictProR")
  gids <- paste0("g", seq_len(4L))
  pheno <- data.frame(
    GID = gids,
    Env = "E1",
    y = c(1.0, 1.4, 2.1, NA_real_),
    stringsAsFactors = FALSE
  )
  K <- diag(length(gids))
  rownames(K) <- colnames(K) <- gids
  kernel <- list(K = K, geno_ids = gids)
  base_args <- list(
    pheno_df = pheno,
    gid_col = "GID",
    env_col = "Env",
    y_col = "y",
    geno_ids = gids,
    method = "gp_exact",
    backend = "auto",
    prediction_output = "test_only",
    output_level = "predict_only",
    point_predictions_only = TRUE,
    return_se = FALSE,
    compute_ai_se = FALSE,
    standardize = "global",
    varcomp_mode = "mom",
    seed = 123L,
    train_idx = c(0L, 1L, 2L),
    test_idx = 3L
  )

  invalid <- ns$gp_bridge_add_gp_exact_controls(base_args, list(gp_iters = 0L, gp_lr = -0.1))
  expect_false("gp_iters" %in% names(invalid))
  expect_false("gp_lr" %in% names(invalid))

  direct_dir <- tempfile("gp_control_direct_")
  cv_dir <- tempfile("gp_control_cv_")
  on.exit(unlink(c(direct_dir, cv_dir), recursive = TRUE, force = TRUE), add = TRUE)

  direct_args <- ns$gp_bridge_add_gp_exact_controls(base_args, list(gp_iters = 9L, gp_lr = 0.02))
  direct <- ns$gp_bridge_prepare_mixed_model_spec(
    args = direct_args,
    kernel = kernel,
    project_root = getwd(),
    dir_path = direct_dir
  )
  direct_spec <- jsonlite::fromJSON(direct$spec_json, simplifyVector = FALSE)
  expect_identical(as.integer(direct_spec$fit$gp_iters), 9L)
  expect_equal(as.numeric(direct_spec$fit$gp_lr), 0.02, tolerance = 1e-12)

  cv_args <- ns$gp_bridge_add_gp_exact_controls(base_args, list(gp_iters = 11L, gp_lr = 0.015))
  cv <- ns$gp_bridge_prepare_mixed_model_cv_spec(
    args = cv_args,
    kernel = kernel,
    folds = list(list(fold_id = "fold1", train_idx = c(0L, 1L, 2L), test_idx = 3L)),
    project_root = getwd(),
    dir_path = cv_dir
  )
  cv_spec <- jsonlite::fromJSON(cv$spec_json, simplifyVector = FALSE)
  expect_identical(as.integer(cv_spec$fit$gp_iters), 11L)
  expect_equal(as.numeric(cv_spec$fit$gp_lr), 0.015, tolerance = 1e-12)
})

test_that("single-environment GP auto tuning uses genotype holdout when data are large enough", {
  n <- 64L
  pheno <- data.frame(
    GID = paste0("G", seq_len(n)),
    Env = "ENV1",
    y = seq_len(n),
    stringsAsFactors = FALSE
  )
  split <- list(
    train_mask = seq_len(n) <= 56L,
    test_mask = seq_len(n) > 56L
  )
  split$train_idx <- which(split$train_mask) - 1L
  split$test_idx <- which(split$test_mask) - 1L
  gmatrix <- diag(n)
  rownames(gmatrix) <- colnames(gmatrix) <- pheno$GID

  policy <- PredictProR:::gp_single_trait_auto_policy(
    pheno_df = pheno,
    split = split,
    heter_groups = NULL,
    gmatrix = gmatrix,
    tune_gp = "auto",
    tuning_objective = "auto"
  )
  expect_true(policy$tune_gp)
  expect_identical(policy$auto_tune_reason, "auto_tune_enabled_dense_single_environment_genotype")
  expect_true("auto_tuning_objective_rmse" %in% policy$decisions)
  expect_identical(policy$include_components, "g")

  blocked <- PredictProR:::gp_tuning_make_genotype_split(
    pheno,
    train_mask = split$train_mask,
    allow_single_env = TRUE,
    seed = 10L
  )
  expect_identical(blocked$status, "ok")
  expect_equal(unique(pheno$Env[blocked$calibration_rows]), "ENV1")
  expect_equal(unique(pheno$Env[blocked$validation_rows]), "ENV1")
})

test_that("auto prediction blending only accepts baseline-dominant validation gains", {
  old <- Sys.getenv(
    c("PREDICTPRO_GP_BLEND_MIN_IMPROVEMENT", "PREDICTPRO_GP_BLEND_MAX_AUTO_ALPHA"),
    unset = NA_character_
  )
  on.exit({
    missing <- is.na(old)
    if (any(missing)) {
      Sys.unsetenv(names(old)[missing])
    }
    if (any(!missing)) {
      do.call(Sys.setenv, as.list(old[!missing]))
    }
  }, add = TRUE)
  Sys.setenv(
    PREDICTPRO_GP_BLEND_MIN_IMPROVEMENT = "0.05",
    PREDICTPRO_GP_BLEND_MAX_AUTO_ALPHA = "0.10"
  )

  baseline_dominant <- PredictProR:::gp_prediction_blend_auto_decision("auto", 0, 0.06)
  partial_blend <- PredictProR:::gp_prediction_blend_auto_decision("auto", 0.2, 0.06)
  weak_gain <- PredictProR:::gp_prediction_blend_auto_decision("auto", 0, 0.01)
  known_env <- PredictProR:::gp_prediction_blend_auto_decision("auto", 0, 0.20, allow_auto = FALSE)
  forced <- PredictProR:::gp_prediction_blend_auto_decision("validation", 0.8, NA_real_)

  expect_true(baseline_dominant$apply)
  expect_identical(baseline_dominant$reason, "validation_blend_auto_baseline_dominant")
  expect_false(partial_blend$apply)
  expect_identical(partial_blend$reason, "validation_blend_alpha_above_auto_limit")
  expect_false(weak_gain$apply)
  expect_identical(weak_gain$reason, "validation_blend_improvement_below_threshold")
  expect_false(known_env$apply)
  expect_identical(known_env$reason, "validation_blend_auto_disabled_for_known_environment")
  expect_true(forced$apply)
  expect_identical(forced$reason, "validation_blend_forced")
})

test_that("GP environment setup prefers GPU-capable torch by default", {
  defaults <- formals(PredictProR::setup_predictgp_env)

  expect_identical(defaults$prefer_gpu, TRUE)
  expect_true("cuda" %in% names(defaults))
  expect_true("torch_version" %in% names(defaults))
  expect_true("index_url" %in% names(defaults))
})

test_that("parallel policy routes GP to GPU unless unavailable or busy", {
  available <- PredictProR:::gp_parallel_model_capability_table(
    models = "KRR",
    model_params_list = list(list()),
    num_gpus = 1L,
    gpu_usage = 10,
    gpu_busy_threshold = 85
  )
  busy <- PredictProR:::gp_parallel_model_capability_table(
    models = "KRR",
    model_params_list = list(list()),
    num_gpus = 1L,
    gpu_usage = 95,
    gpu_busy_threshold = 85
  )
  missing <- PredictProR:::gp_parallel_model_capability_table(
    models = "KRR",
    model_params_list = list(list()),
    num_gpus = 0L,
    gpu_usage = NA_real_,
    gpu_busy_threshold = 85
  )

  expect_true(available$gpu_selected)
  expect_identical(available$device_route, "auto")
  expect_match(available$reasons, "gp_gpu_selected")
  expect_false(busy$gpu_selected)
  expect_identical(busy$device_route, "cpu")
  expect_match(busy$reasons, "gp_cpu_gpu_busy")
  expect_false(missing$gpu_selected)
  expect_identical(missing$device_route, "cpu")
  expect_match(missing$reasons, "gp_cpu_no_gpu")
})

test_that("GP tuning objective aliases normalize to supported objectives", {
  expect_identical(PredictProR:::gp_tuning_normalize_objective("raw-rmse"), "rmse")
  expect_identical(PredictProR:::gp_tuning_normalize_objective("within env rank"), "within_env_spearman")
  expect_identical(PredictProR:::gp_tuning_normalize_objective("centered"), "centered_rmse")
  expect_error(
    PredictProR:::gp_tuning_normalize_objective("not_an_objective"),
    "`tuning_objective`"
  )
})

test_that("reaction-norm tuning grid records effective nonzero env weights", {
  grid <- data.frame(
    lambda = 0.1,
    kenv_kernel = "matern32",
    mean_mode = "none",
    w_g = c(1, 0.7),
    w_ge = c(0, 0.2),
    w_e = c(0, 0.3)
  )
  effective <- PredictProR:::gp_tuning_effective_reaction_norm_grid(grid)

  expect_equal(effective$w_ge, c(0.05, 0.2), tolerance = 1e-8)
  expect_equal(effective$w_e, c(0.5, 0.3), tolerance = 1e-8)
})

test_that("single-trait auto policy enables only generalizable dense MET tuning", {
  gids <- paste0("g", seq_len(60))
  train_envs <- paste0("LOC", seq_len(4), "_2023")
  test_env <- "LOC5_2024"
  pheno <- expand.grid(
    GID = gids,
    Env = c(train_envs, test_env),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pheno$y <- seq_len(nrow(pheno))
  pheno$y[pheno$Env == test_env] <- NA_real_
  split <- PredictProR:::gp_public_train_test_idx(
    pheno,
    gid_col = "GID",
    y_col = "y"
  )

  dense_policy <- PredictProR:::gp_single_trait_auto_policy(
    pheno_df = pheno,
    split = split,
    heter_groups = "Env",
    env_covariates = data.frame(Env = c(train_envs, test_env), x = seq_len(5)),
    gmatrix = diag(length(gids)),
    tune_gp = "auto",
    tuning_objective = "auto",
    large_n_threshold = 50000L
  )
  factor_policy <- PredictProR:::gp_single_trait_auto_policy(
    pheno_df = pheno,
    split = split,
    heter_groups = "Env",
    env_covariates = data.frame(Env = c(train_envs, test_env), x = seq_len(5)),
    gmatrix = NULL,
    gp_factor_cache = list(type = "memmap"),
    geno_ids = gids,
    tune_gp = "auto",
    tuning_objective = "auto",
    large_n_threshold = 50000L
  )
  prospective_without_env_covariates <- PredictProR:::gp_single_trait_auto_policy(
    pheno_df = pheno,
    split = split,
    heter_groups = "Env",
    env_covariates = NULL,
    gmatrix = diag(length(gids)),
    tune_gp = "auto",
    tuning_objective = "auto",
    large_n_threshold = 50000L
  )

  known_pheno <- pheno
  known_pheno$y <- seq_len(nrow(known_pheno))
  known_prediction_ids <- paste0("g", 56:60)
  known_pheno$y[known_pheno$GID %in% known_prediction_ids] <- NA_real_
  known_split <- PredictProR:::gp_public_train_test_idx(
    known_pheno,
    gid_col = "GID",
    y_col = "y"
  )
  known_policy <- PredictProR:::gp_single_trait_auto_policy(
    pheno_df = known_pheno,
    split = known_split,
    heter_groups = "Env",
    env_covariates = NULL,
    gmatrix = diag(length(gids)),
    tune_gp = "auto",
    tuning_objective = "auto",
    large_n_threshold = 50000L
  )

  expect_true(dense_policy$tune_gp)
  expect_true(dense_policy$gp_tiered_dispatch)
  expect_identical(dense_policy$tuning_objective, "centered_rmse")
  expect_match(paste(dense_policy$decisions, collapse = ";"), "auto_tune_enabled")

  expect_true(factor_policy$tune_gp)
  expect_identical(factor_policy$gp_backend, "operator")
  expect_identical(factor_policy$gp_dtype, "float32")
  expect_identical(factor_policy$tuning_source, "factor_cache_subset")
  expect_match(factor_policy$auto_tune_reason, "factor_cache_subset")
  expect_equal(factor_policy$tuning_limits$max_calibration_rows, 4000L)

  expect_false(prospective_without_env_covariates$tune_gp)
  expect_match(
    paste(prospective_without_env_covariates$decisions, collapse = ";"),
    "auto_tune_skipped_no_env_covariates"
  )

  expect_true(known_policy$tune_gp)
  expect_false(known_policy$gp_tiered_dispatch)
  expect_identical(known_policy$tuning_objective, "rmse")
  expect_match(
    paste(known_policy$decisions, collapse = ";"),
    "auto_tune_enabled_dense_met_known_environment"
  )
})

test_that("GP tuning selection can optimize rank-aware objectives", {
  detail <- data.frame(
    candidate_id = 1:3,
    lambda = c(0.1, 0.1, 0.3),
    kenv_kernel = c("matern32", "matern32", "rbf"),
    mean_mode = c("none", "location_offset", "location_offset"),
    validation_rmse = c(1.0, 1.1, 1.3),
    validation_centered_rmse = c(0.9, 0.7, 0.8),
    validation_centered_rmse_se = c(NA_real_, NA_real_, NA_real_),
    validation_within_env_spearman = c(0.1, 0.2, 0.6),
    validation_within_env_spearman_se = c(NA_real_, NA_real_, NA_real_),
    complexity = c(2, 3, 5),
    failed = FALSE,
    stringsAsFactors = FALSE
  )

  by_rmse <- PredictProR:::gp_tuning_select_candidate(detail, "rmse")
  by_centered <- PredictProR:::gp_tuning_select_candidate(detail, "centered_rmse")
  by_rank <- PredictProR:::gp_tuning_select_candidate(detail, "within_env_spearman")

  expect_equal(by_rmse$selected$candidate_id, 1)
  expect_equal(by_centered$selected$candidate_id, 2)
  expect_equal(by_rank$selected$candidate_id, 3)
  expect_identical(by_rank$selected$tuning_objective, "within_env_spearman")
  expect_equal(by_rank$selected$validation_objective, 0.6)
})

test_that("GP tuning selection falls back when auto centered objective is unavailable", {
  detail <- data.frame(
    candidate_id = 1:3,
    lambda = c(0.03, 0.1, 0.3),
    kenv_kernel = c("matern32", "matern32", "rbf"),
    mean_mode = c("none", "none", "location_offset"),
    validation_rmse = c(1.2, 0.8, 1.1),
    validation_centered_rmse = c(NA_real_, NA_real_, NA_real_),
    validation_centered_rmse_se = c(NA_real_, NA_real_, NA_real_),
    validation_centered_cor = c(NA_real_, NA_real_, NA_real_),
    validation_within_env_spearman = c(NA_real_, NA_real_, NA_real_),
    complexity = c(2, 2, 5),
    failed = FALSE,
    stringsAsFactors = FALSE
  )

  selected <- PredictProR:::gp_tuning_select_candidate(detail, "centered_rmse")

  expect_equal(selected$selected$candidate_id, 2)
  expect_identical(selected$selected$requested_tuning_objective, "centered_rmse")
  expect_identical(selected$selected$tuning_objective, "rmse")
  expect_match(selected$reason, "selected by rmse")
})

test_that("GP validation uncertainty calibration estimates and applies SE scale", {
  fit <- list(predictions = data.frame(
    GID = paste0("g", 1:6),
    Env = rep(c("e1", "e2"), each = 3),
    Prediction = c(1, 2, 3, 1, 2, 3),
    Prediction_SE_observed = rep(4, 6),
    Prediction_Var_observed = rep(16, 6),
    SE_latent = rep(3, 6),
    PEV = rep(9, 6),
    stringsAsFactors = FALSE
  ))
  truth <- data.frame(
    GID = paste0("g", 1:6),
    Env = rep(c("e1", "e2"), each = 3),
    Observed = c(3, 4, 5, 3, 4, 5),
    stringsAsFactors = FALSE
  )

  calibration <- PredictProR:::gp_tuning_estimate_uncertainty_calibration(
    fit,
    truth = truth,
    offsets = rep(0, 6),
    min_n = 2,
    min_scale = 0.1,
    max_scale = 10,
    prior_n = 0,
    prior_env = 1
  )
  expect_identical(calibration$status, "ok")
  expect_equal(calibration$scale, 0.5, tolerance = 1e-8)
  expect_identical(calibration$variance_column, "Prediction_Var_observed")
  expect_identical(calibration$se_column, "Prediction_SE_observed")
  expect_identical(calibration$uncertainty_target, "observed")
  expect_false(calibration$variance_from_se)

  applied <- PredictProR:::gp_apply_uncertainty_calibration(
    list(predictions = fit$predictions, result = list(predictions = fit$predictions)),
    calibration
  )
  pred <- applied$result$predictions
  expect_true(applied$info$applied)
  expect_equal(pred$Prediction, fit$predictions$Prediction)
  expect_equal(pred$Prediction_SE_observed, rep(2, 6), tolerance = 1e-8)
  expect_equal(pred$Prediction_Var_observed, rep(4, 6), tolerance = 1e-8)
  expect_equal(pred$SE_latent, rep(1.5, 6), tolerance = 1e-8)
  expect_equal(pred$PEV, rep(2.25, 6), tolerance = 1e-8)
  expect_equal(pred$Prediction_SE_Calibration_Scale, rep(0.5, 6), tolerance = 1e-8)
})

test_that("GP uncertainty calibration falls back across available variance columns", {
  fit <- list(predictions = data.frame(
    GID = paste0("g", 1:4),
    Env = rep(c("e1", "e2"), each = 2),
    Prediction = c(1, 2, 1, 2),
    Prediction_Var_observed = NA_real_,
    Prediction_variance = rep(16, 4),
    stringsAsFactors = FALSE
  ))
  truth <- data.frame(
    GID = paste0("g", 1:4),
    Env = rep(c("e1", "e2"), each = 2),
    Observed = c(3, 4, 3, 4),
    stringsAsFactors = FALSE
  )

  calibration <- PredictProR:::gp_tuning_estimate_uncertainty_calibration(
    fit,
    truth = truth,
    offsets = rep(0, 4),
    min_n = 2,
    min_scale = 0.1,
    max_scale = 10,
    prior_n = 0,
    prior_env = 1
  )

  expect_identical(calibration$status, "ok")
  expect_equal(calibration$scale, 0.5, tolerance = 1e-8)
  expect_identical(calibration$variance_column, "Prediction_variance")
  expect_identical(calibration$uncertainty_target, "observed")
})

test_that("multi-trait GP uncertainty aliases expose observed and latent SE columns", {
  out <- list(
    predictions = data.frame(
      gid = rep(c("g1", "g2"), each = 2),
      trait = rep(c("t1", "t2"), times = 2),
      Prediction = c(1, 2, 3, 4),
      SE = c(2, 3, 4, 5),
      SE_latent = c(1, 1.5, 2, 2.5),
      PEV = c(1, 2.25, 4, 6.25),
      stringsAsFactors = FALSE
    ),
    result = list()
  )

  normalized <- PredictProR:::gp_public_normalize_multitrait_uncertainty(out)
  pred <- normalized$predictions

  expect_true(all(c(
    "Prediction_SE_observed", "Prediction_SE_latent",
    "Prediction_Var_observed", "Prediction_Var_latent"
  ) %in% names(pred)))
  expect_equal(pred$Prediction_SE_observed, pred$SE, tolerance = 1e-8)
  expect_equal(pred$Prediction_SE_latent, pred$SE_latent, tolerance = 1e-8)
  expect_equal(pred$Prediction_SE_observed^2, pred$Prediction_Var_observed, tolerance = 1e-8)
  expect_equal(pred$PEV, pred$Prediction_Var_latent, tolerance = 1e-8)
})

test_that("multi-trait GP validation blend scales SE and variance", {
  out <- list(
    predictions = data.frame(
      gid = c("g1", "g1"),
      trait = c("t1", "t2"),
      Prediction = c(10, 20),
      Prediction_SE_observed = c(4, 8),
      Prediction_Var_observed = c(16, 64),
      Prediction_SE_latent = c(2, 4),
      Prediction_Var_latent = c(4, 16),
      stringsAsFactors = FALSE
    ),
    result = list()
  )
  baseline <- data.frame(
    gid = c("g1", "g1"),
    trait = c("t1", "t2"),
    Prediction_Baseline = c(2, 4),
    stringsAsFactors = FALSE
  )

  blended <- PredictProR:::gp_public_multitrait_apply_blend(
    out,
    baseline = baseline,
    key_cols = c("gid", "trait"),
    alpha = 0.25,
    label = "test_blend"
  )
  pred <- blended$result$predictions

  expect_true(blended$info$applied)
  expect_equal(pred$Prediction, c(4, 8), tolerance = 1e-8)
  expect_equal(pred$Prediction_SE_observed, c(1, 2), tolerance = 1e-8)
  expect_equal(pred$Prediction_Var_observed, c(1, 4), tolerance = 1e-8)
  expect_equal(pred$Prediction_SE_latent, c(0.5, 1), tolerance = 1e-8)
  expect_equal(pred$Prediction_Var_latent, c(0.25, 1), tolerance = 1e-8)
})

test_that("multi-trait GP validation calibration records uncertainty target", {
  joined <- data.frame(
    gid = paste0("g", 1:6),
    trait = rep(c("t1", "t2"), 3),
    Prediction = rep(1, 6),
    Observed = rep(3, 6),
    Prediction_SE_observed = rep(4, 6),
    Prediction_Var_observed = rep(16, 6),
    stringsAsFactors = FALSE
  )

  calibration <- PredictProR:::gp_public_multitrait_uncertainty_calibration(
    joined,
    min_n = 2,
    min_scale = 0.1,
    max_scale = 10,
    prior_n = 0,
    prior_env = 1
  )

  expect_identical(calibration$status, "ok")
  expect_equal(calibration$scale, 0.5, tolerance = 1e-8)
  expect_identical(calibration$variance_column, "Prediction_Var_observed")
  expect_identical(calibration$se_column, "Prediction_SE_observed")
  expect_identical(calibration$uncertainty_target, "observed")
})

test_that("blocked GP tuning split uses only finite training rows", {
  gids <- paste0("g", seq_len(6))
  envs <- c("LOC1_2021", "LOC1_2022", "LOC2_2021", "LOC2_2022")
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pheno$y <- seq_len(nrow(pheno))
  prediction_ids <- c("g5", "g6")
  pheno$y[pheno$GID %in% prediction_ids] <- NA_real_

  split <- PredictProR:::gp_public_train_test_idx(
    pheno,
    gid_col = "GID",
    y_col = "y"
  )
  blocked <- PredictProR:::gp_tuning_make_blocked_split(
    pheno,
    train_mask = split$train_mask,
    seed = 7L
  )

  expect_identical(blocked$status, "ok")
  expect_true(all(blocked$calibration_rows %in% which(split$train_mask)))
  expect_true(all(blocked$validation_rows %in% which(split$train_mask)))
  expect_true(all(is.finite(pheno$y[blocked$validation_rows])))
  expect_false(any(blocked$validation_ids %in% prediction_ids))
  expect_true(all(pheno$Env[blocked$validation_rows] %in% blocked$validation_envs))
  expect_true(all(pheno$Env[blocked$calibration_rows] %in% blocked$calibration_envs))
})

test_that("genotype-blocked GP tuning split holds out training genotypes", {
  gids <- paste0("g", seq_len(8))
  envs <- c("LOC1_2021", "LOC1_2022", "LOC2_2021")
  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pheno$y <- seq_len(nrow(pheno))
  prediction_ids <- c("g7", "g8")
  pheno$y[pheno$GID %in% prediction_ids] <- NA_real_

  split <- PredictProR:::gp_public_train_test_idx(
    pheno,
    gid_col = "GID",
    y_col = "y"
  )
  blocked <- PredictProR:::gp_tuning_make_genotype_split(
    pheno,
    train_mask = split$train_mask,
    seed = 11L
  )

  expect_identical(blocked$status, "ok")
  expect_true(all(blocked$calibration_rows %in% which(split$train_mask)))
  expect_true(all(blocked$validation_rows %in% which(split$train_mask)))
  expect_true(all(is.finite(pheno$y[blocked$validation_rows])))
  expect_false(any(blocked$validation_ids %in% prediction_ids))
  expect_false(any(blocked$calibration_ids %in% blocked$validation_ids))
})

test_that("location mean adjustment is estimated from training rows only", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Env = c("LOC1_2021", "LOC1_2021", "LOC1_2022", "LOC1_2022"),
    y = c(1, 3, 1000, NA),
    stringsAsFactors = FALSE
  )
  adjusted <- PredictProR:::gp_tuning_apply_mean_adjustment(
    pheno,
    training_rows = c(1L, 2L),
    mode = "location_offset"
  )

  expect_equal(adjusted$offsets, rep(2, 4), tolerance = 1e-8)
  expect_equal(adjusted$pheno_df$y[1:2], c(-1, 1), tolerance = 1e-8)
  expect_equal(adjusted$pheno_df$y[[3L]], 998, tolerance = 1e-8)
  expect_true(is.na(adjusted$pheno_df$y[[4L]]))
})

test_that("MT-MET dense REML rejects FA covariance before Python startup", {
  gids <- paste0("g", seq_len(4))
  pheno <- expand.grid(
    GID = gids,
    Env = c("e1", "e2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pheno$Trait1 <- seq_len(nrow(pheno))
  pheno$Trait2 <- pheno$Trait1 + 0.5

  gmatrix <- diag(length(gids))
  rownames(gmatrix) <- colnames(gmatrix) <- gids

  expect_error(
    PredictProR::gp_multi_trait_met_model(
      pheno_data = pheno,
      gmatrix = gmatrix,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      heter_groups = "Env",
      test_set = "g4",
      varcomp_mode = "reml",
      trait_structure = "fa",
      trait_fa_rank = 1L,
      gxe_trait_structure = "fa",
      gxe_trait_fa_rank = 1L
    ),
    "Dense MT-MET REML currently supports only unstructured"
  )
})

test_that("environment covariates and environment similarities are mutually exclusive", {
  expect_error(
    PredictProR::gp_single_trait_model(
      pheno_data = data.frame(),
      gmatrix = diag(1),
      response = "y",
      gen_name = "gid",
      env_similarity = diag(1),
      env_covariates = data.frame(Env = "e1", x = 1)
    ),
    "Pass exactly one"
  )
  expect_error(
    PredictProR::gp_multi_trait_met_model(
      pheno_data = data.frame(),
      gmatrix = diag(1),
      response = c("t1", "t2"),
      gen_name = "gid",
      heter_groups = "env",
      env_similarity = diag(1),
      env_covariates = data.frame(env = "e1", x = 1)
    ),
    "Pass exactly one"
  )
})

test_that("public GP kernel banks retain gmatrix plus additional kernels", {
  ids <- paste0("g", seq_len(4L))
  K1 <- diag(4L)
  K2 <- matrix(0.2, nrow = 4L, ncol = 4L) + diag(4L)
  dimnames(K1) <- dimnames(K2) <- list(ids, ids)

  collected <- PredictProR:::gp_collect_kernel_inputs(
    gmatrix = K1,
    kernel_list = list(omic = K2)
  )
  bank <- PredictProR:::gp_public_additive_kernel_bank(
    PredictProR:::gp_public_kernel_bank_inputs(collected, pheno_ids = ids)
  )

  expect_equal(length(bank$Ks), 2L)
  expect_identical(names(bank$Ks), c("A", "omic"))
  expect_equal(bank$Ks$A, K1)
  expect_equal(bank$Ks$omic, K2)
})

test_that("MT-MET kernel weights align by name and reject invalid inputs", {
  ids <- paste0("g", seq_len(3L))
  K <- diag(3L)
  dimnames(K) <- list(ids, ids)
  bank <- list(K = K, Ks = list(A = K, omic = K), geno_ids = ids)

  expect_equal(
    PredictProR:::gp_public_kernel_weights(bank, c(omic = 0.25, A = 0.75)),
    c(A = 0.75, omic = 0.25)
  )
  expect_equal(
    PredictProR:::gp_public_kernel_weights(bank),
    c(A = 1, omic = 1)
  )
  expect_error(
    PredictProR:::gp_public_kernel_weights(bank, c(A = 1)),
    "must match"
  )
  expect_error(
    PredictProR:::gp_public_kernel_weights(bank, c(0, 0)),
    "positive value"
  )
})

test_that("environment covariate id columns are normalized for backend alignment", {
  ec <- data.frame(environ = c("e1", "e2"), x = c(0.1, 0.2))

  single <- PredictProR:::gp_prepare_env_covariates(
    ec,
    source_env_col = "environ",
    target_env_col = "Env"
  )
  met <- PredictProR:::gp_prepare_env_covariates(
    ec,
    source_env_col = "environ",
    target_env_col = "env"
  )

  expect_true("Env" %in% names(single))
  expect_true("env" %in% names(met))
})

test_that("factor-cache-only MT kernel inputs do not require a dense GRM", {
  gids <- paste0("g", seq_len(4))
  pheno_ids <- rep(gids, each = 2L)

  kernel <- PredictProR:::gp_public_kernel_inputs(
    gmatrix = NULL,
    pheno_ids = pheno_ids,
    geno_ids = gids,
    allow_factor_cache = TRUE
  )
  expect_identical(kernel$geno_ids, gids)
  expect_equal(dim(kernel$K), c(1L, 1L))

  pheno <- data.frame(
    GID = gids,
    Trait1 = seq_along(gids),
    Trait2 = seq_along(gids) + 0.5,
    stringsAsFactors = FALSE
  )
  expect_error(
    PredictProR::gp_multi_trait_model(
      pheno_data = pheno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      test_set = "g4",
      geno_ids = gids,
      varcomp_mode = "reml",
      gp_factor_cache = list(type = "memmap")
    ),
    "gmatrix.*omitted only when `varcomp_mode = \"mom\"`"
  )
  expect_error(
    PredictProR::gp_multi_trait_model(
      pheno_data = pheno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      test_set = "g4",
      geno_ids = gids,
      varcomp_mode = "mom"
    ),
    "gmatrix.*unless"
  )
})

test_that("large MET warns only when no environment kernel is supplied", {
  envs <- paste0("e", seq_len(6))

  expect_warning(
    PredictProR:::gp_warn_large_met_without_env_kernel(
      envs,
      has_env_kernel = FALSE,
      context = "Multi-trait MET GP"
    ),
    "provide `env_similarity` or `env_covariates`"
  )
  expect_silent(
    PredictProR:::gp_warn_large_met_without_env_kernel(
      envs,
      has_env_kernel = TRUE,
      context = "Multi-trait MET GP"
    )
  )
})
