if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

make_gp_bridge_smoke_data <- function(n = 12L, p = 18L, seed = 20260430L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- paste0("g", seq_len(n))
  K <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
  K <- K / mean(diag(K))
  y <- as.numeric(X[, 1] - 0.7 * X[, 2] + rnorm(n, sd = 0.2))
  pheno <- data.frame(GID = rownames(X), Trait = y, stringsAsFactors = FALSE)
  list(pheno = pheno, K = K)
}

make_gp_public_smoke_data <- function(n = 10L, p = 14L, seed = 20260501L) {
  set.seed(seed)
  gids <- paste0("g", seq_len(n))
  envs <- c("E1", "E2")
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- gids
  K <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
  K <- K / mean(diag(K))
  diag(K) <- diag(K) + 1e-6
  rownames(K) <- colnames(K) <- gids

  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gid_idx <- match(pheno$GID, gids)
  env_effect <- c(E1 = -0.2, E2 = 0.35)[pheno$Env]
  signal <- X[gid_idx, 1] - 0.4 * X[gid_idx, 2]
  pheno$Trait1 <- signal + env_effect + rnorm(nrow(pheno), sd = 0.1)
  pheno$Trait2 <- 0.7 * signal + 0.3 * X[gid_idx, 3] - 0.5 * env_effect + rnorm(nrow(pheno), sd = 0.1)

  env_similarity <- matrix(c(1, 0.25, 0.25, 1), nrow = 2L)
  rownames(env_similarity) <- colnames(env_similarity) <- envs

  list(pheno = pheno, K = K, gids = gids, envs = envs, env_similarity = env_similarity)
}

expect_gp_se_variance <- function(predictions, se_col, var_col) {
  predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
  expect_true(all(c(se_col, var_col) %in% names(predictions)))
  expect_true(all(is.finite(as.numeric(predictions[[se_col]]))))
  expect_true(all(is.finite(as.numeric(predictions[[var_col]]))))
  expect_true(all(as.numeric(predictions[[se_col]]) >= 0))
  expect_true(all(as.numeric(predictions[[var_col]]) >= 0))
  expect_equal(
    as.numeric(predictions[[se_col]])^2,
    as.numeric(predictions[[var_col]]),
    tolerance = 1e-6
  )
}

expect_gp_direct_policy <- function(out, bridge_command) {
  info <- out$info
  execution <- info$gp_execution
  if (is.null(execution)) {
    execution <- info$execution
  }

  expect_identical(execution, "python_subprocess")
  expect_identical(info$bridge_command, bridge_command)
  expect_true(is.list(info$gp_auto_policy))
  expect_true(isTRUE(info$gp_auto_policy$enabled))

  device <- info$gp_auto_policy$device
  expect_true(is.list(device))
  expect_true(as.character(device$route) %in% c("auto", "cpu"))
  expect_true(as.character(device$reason) %in% c(
    "gp_gpu_auto_available",
    "gp_cpu_no_gpu",
    "gp_cpu_gpu_busy"
  ))
  if (identical(as.character(device$reason), "gp_gpu_auto_available")) {
    expect_identical(as.character(device$route), "auto")
  } else {
    expect_identical(as.character(device$route), "cpu")
  }
}

expect_gp_multitrait_uncertainty_contract <- function(predictions) {
  predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
  expect_true(all(c("SE", "SE_latent", "PEV") %in% names(predictions)))
  expect_true(all(is.finite(as.numeric(predictions$SE))))
  expect_true(all(is.finite(as.numeric(predictions$SE_latent))))
  expect_true(all(is.finite(as.numeric(predictions$PEV))))
  expect_true(all(as.numeric(predictions$SE) >= 0))
  expect_true(all(as.numeric(predictions$SE_latent) >= 0))
  expect_true(all(as.numeric(predictions$PEV) >= 0))
  expect_true(all(as.numeric(predictions$SE) >= as.numeric(predictions$SE_latent)))
  expect_equal(
    as.numeric(predictions$SE_latent)^2,
    as.numeric(predictions$PEV),
    tolerance = 1e-6
  )
}

expect_gp_correlation_matrix <- function(x, n_traits) {
  x <- as.matrix(x)
  expect_equal(dim(x), c(n_traits, n_traits))
  expect_true(all(is.finite(x)))
  expect_equal(unname(diag(x)), rep(1, n_traits), tolerance = 1e-8)
  expect_true(max(abs(x), na.rm = TRUE) <= 1 + 1e-8)
  expect_equal(unname(x), unname(t(x)), tolerance = 1e-8)
}

test_gp_python <- function() {
  local_candidates <- if (exists("gp_local_virtualenv_python_candidates", envir = asNamespace("PredictProR"), inherits = FALSE)) {
    get("gp_local_virtualenv_python_candidates", envir = asNamespace("PredictProR"))(
      ".venv-gp",
      roots = c(
        getwd(),
        normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
      )
    )
  } else {
    c(
      file.path(getwd(), ".venv-gp", "Scripts", "python.exe"),
      file.path(getwd(), ".venv-gp", "bin", "python"),
      file.path(getwd(), "..", "..", ".venv-gp", "Scripts", "python.exe"),
      file.path(getwd(), "..", "..", ".venv-gp", "bin", "python")
    )
  }
  candidates <- unique(c(
    Sys.getenv("PREDICTPRO_GP_PYTHON", unset = ""),
    local_candidates,
    unname(Sys.which("python3")),
    unname(Sys.which("python"))
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  skip_if_not(length(candidates) > 0L, "Dedicated GP Python runtime was not found.")
  normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("GP bridge model registry is available", {
  ns <- asNamespace("PredictProR")
  expect_true(exists("gp_lowrank_supported_models", envir = ns, inherits = FALSE))
  expect_true(all(c("KRR", "GP", "GP_FA", "LowRankGP") %in% get("gp_lowrank_supported_models", envir = ns)()))
  expect_identical(get("gp_backend_resolve_model", envir = ns)("LowRankGP"), "KRR")
  expect_identical(get("gp_backend_method_map", envir = ns)("KRR"), "krr_exact")
  expect_identical(get("gp_backend_method_map", envir = ns)("GP"), "gp_exact")
  expect_identical(get("gp_backend_method_map", envir = ns)("GP_FA"), "gp_icm_fa")
})

test_that("model_execute true prediction supports KRR via Python bridge", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python
  )

  dat <- make_gp_bridge_smoke_data()
  pheno <- dat$pheno
  pheno$Trait[c(3L, 9L)] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$K,
    response = "Trait",
    gen_name = "GID",
    GS_model = "KRR",
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  )

  expect_true("Trait" %in% names(out))
  expect_true(is.data.frame(out$Trait$model_results$predicted_values))
  expect_true(any(out$Trait$model_results$predicted_values$Train_Test_Label == "Test"))
})

test_that("model_execute GP return_se emits finite SE and PEV", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python
  )

  dat <- make_gp_bridge_smoke_data(seed = 20260501L)
  pheno <- dat$pheno
  pheno$Trait[c(2L, 11L)] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$K,
    response = "Trait",
    gen_name = "GID",
    GS_model = "KRR",
    gp_return_se = TRUE,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  )

  pred <- out$Trait$model_results$predicted_values
  expect_true(all(c("Standard_error", "PEV") %in% names(pred)))
  test_rows <- pred$Train_Test_Label == "Test"
  expect_true(any(test_rows))
  expect_true(all(is.finite(pred$Standard_error[test_rows])))
  expect_true(all(is.finite(pred$PEV[test_rows])))
  expect_equal(pred$Standard_error[test_rows]^2, pred$PEV[test_rows], tolerance = 1e-6)

  params <- out$Trait$model_results$model_parameters
  expect_true(any(params$stat == "gp_output_level" & params$summary == "full_vc"))
  expect_true(PredictProR::validate_prediction_output(out)$valid)
})

test_that("single-environment GP family models share the canonical public schema", {
  gp_python <- test_gp_python()
  Sys.setenv(PREDICTPRO_GP_PYTHON = gp_python)

  dat <- make_gp_bridge_smoke_data(n = 12L, p = 16L, seed = 20260812L)
  pheno <- dat$pheno
  pheno$Trait[c(2L, 11L)] <- NA_real_
  expected_columns <- PredictProR:::gp_gaussian_prediction_columns(
    include_env = FALSE,
    include_trait = FALSE
  )

  # LowRankGP (Scalable-GBLUP) is rejected as a single-trait model: its
  # single-trait route resolves to KRR and would duplicate Kernel-GBLUP.
  results <- lapply(c("KRR", "GP"), function(model) {
    fit <- PredictProR:::gp_backend_gaussian_model(
      model_name = model,
      pheno_data = pheno,
      response = "Trait",
      gen_name = "GID",
      gmatrix = dat$K,
      gp_output_level = "predict_with_se",
      gp_prediction_output = "all",
      gp_iters = 3L,
      gp_lr = 0.02,
      gp_learn_scales = FALSE
    )
    PredictProR:::gp_standardize_public_model_result(fit, gen_name = "GID")
  })

  for (result in results) {
    expect_true(all(c(
      "model_parameters", "predicted_values", "diagnostic_plots",
      "variance_components"
    ) %in% names(result)))
    expect_identical(names(result$predicted_values), expected_columns)
    expect_true(PredictProR::validate_prediction_output(result)$valid)
    expect_true(all(is.finite(result$predicted_values$Predicted_value)))
    expect_equal(
      result$predicted_values$Standard_error^2,
      result$predicted_values$PEV,
      tolerance = 1e-6
    )
  }
})

test_that("multi-kernel FA covariance, variance, and correlation outputs agree", {
  gp_python <- test_gp_python()
  Sys.setenv(PREDICTPRO_GP_PYTHON = gp_python)

  dat <- make_gp_public_smoke_data(n = 12L, p = 16L, seed = 20260813L)
  second_kernel <- 0.55 * dat$K + 0.45 * diag(nrow(dat$K))
  dimnames(second_kernel) <- dimnames(dat$K)
  pheno <- dat$pheno[, c("GID", "Env", "Trait1")]
  pheno$Trait1[pheno$GID %in% tail(dat$gids, 3L)] <- NA_real_

  fit <- PredictProR:::gp_backend_gaussian_model(
    model_name = "GP_FA",
    pheno_data = pheno,
    response = "Trait1",
    gen_name = "GID",
    heter_groups = "Env",
    gmatrix = dat$K,
    omic1_kernel = second_kernel,
    gp_output_level = "full_vc",
    gp_varcomp_mode = "reml",
    gp_fa_rank = 1L,
    gp_prediction_output = "all",
    gp_iters = 5L,
    gp_lr = 0.02,
    gp_learn_scales = TRUE
  )
  result <- PredictProR:::gp_standardize_public_model_result(
    fit,
    gen_name = "GID",
    heter_groups = "Env"
  )
  raw <- fit$gp_result

  total_cov <- as.matrix(raw$env_covariance_fa)
  total_cor <- as.matrix(raw$env_correlation_fa)
  interaction <- raw$interaction_correlation_summary[["G:Env"]]
  by_kernel <- interaction$cov_by_kernel
  summed_cov <- Reduce(`+`, lapply(by_kernel, as.matrix))
  variance_summary <- raw$interaction_variance_summary[["G:Env"]]

  expect_equal(total_cov, t(total_cov), tolerance = 1e-8)
  expect_true(min(eigen(total_cov, symmetric = TRUE, only.values = TRUE)$values) >= -1e-8)
  expect_equal(unname(summed_cov), unname(total_cov), tolerance = 1e-7)
  expect_equal(variance_summary$interaction_variance_total, diag(total_cov), tolerance = 1e-7)
  expect_equal(total_cor, stats::cov2cor(total_cov), tolerance = 1e-7)
  expect_gp_correlation_matrix(total_cor, length(dat$envs))

  expect_equal(result$Genetic_covariance_environments, total_cov, tolerance = 1e-7)
  expect_equal(result$Genetic_correlation_environments, total_cor, tolerance = 1e-7)
  expect_identical(names(result$covariance_components), PredictProR:::gp_covariance_component_columns())
  expect_identical(names(result$correlation_components), PredictProR:::gp_covariance_component_columns())

  vc <- result$variance_components
  env_vg <- vc[match(paste0("genetic_variance_", dat$envs), vc$Component), "Components"]
  cov_diag <- diag(total_cov)[match(dat$envs, as.character(interaction$levels))]
  expect_equal(as.numeric(env_vg), as.numeric(cov_diag), tolerance = 1e-7)
  h2 <- vc[grepl("^heritability_", vc$Component), "Components"]
  expect_true(all(h2[is.finite(h2)] >= 0 & h2[is.finite(h2)] <= 1))
  expect_true(PredictProR::validate_prediction_output(result, heter_groups = "Env")$valid)
})

test_that("model_execute routes multi-trait GP through direct public wrapper", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_EXECUTION = "direct",
    PREDICTPRO_GP_RETICULATE_FALLBACK = "false",
    PREDICTPRO_GP_CV_FRONTDOOR_FAST = "true"
  )
  Sys.unsetenv("PREDICTPRO_GP_DEVICE")

  dat <- make_gp_public_smoke_data(n = 10L, p = 14L, seed = 20260504L)
  pheno <- dat$pheno[dat$pheno$Env == dat$envs[[1L]], , drop = FALSE]
  test_ids <- tail(dat$gids, 2L)
  pheno[pheno$GID %in% test_ids, c("Trait1", "Trait2")] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$K,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    # Kernel-GBLUP (KRR) is single-trait only; joint multi-trait GP uses
    # Gaussian-Process-GBLUP, FA-GBLUP or Scalable-GBLUP.
    GS_model = "Gaussian-Process-GBLUP",
    multi_trait_gp = TRUE,
    gp_return_se = TRUE,
    gp_return_trait_correlations = TRUE,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  )

  pred <- out$model_results$predicted_values
  # The multi-trait GP true-prediction route predicts every row (training and
  # held-out); the held-out rows are labelled Test.
  expect_equal(nrow(pred), nrow(pheno) * 2L)
  expect_equal(sum(pred$Train_Test_Label == "Test"), length(test_ids) * 2L)
  expect_true(all(c("Standard_error", "PEV") %in% names(pred)))
  expect_true(all(is.finite(pred$Predicted_value)))
  expect_true(all(is.finite(pred$Standard_error)))
  expect_true(all(is.finite(pred$PEV)))
  expect_equal(pred$Standard_error^2, pred$PEV, tolerance = 1e-6)
  expect_true("Genetic_correlation_traits" %in% names(out$model_results))
  expect_gp_correlation_matrix(out$model_results$Genetic_correlation_traits, 2L)
  expect_true(any(out$model_results$model_parameters$stat == "bridge_command" &
    out$model_results$model_parameters$summary == "fit-multi-trait-spec"))
  expect_true(PredictProR::validate_prediction_output(out)$valid)
})

test_that("model_execute routes multi-trait MET GP through direct public wrapper", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_EXECUTION = "direct",
    PREDICTPRO_GP_RETICULATE_FALLBACK = "false"
  )

  dat <- make_gp_public_smoke_data(n = 8L, p = 12L, seed = 20260506L)
  test_ids <- tail(dat$gids, 2L)
  pheno <- dat$pheno
  pheno[pheno$GID %in% test_ids, c("Trait1", "Trait2")] <- NA_real_

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = dat$K,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    heter_groups = "Env",
    env_similarity = dat$env_similarity,
    GS_model = "Gaussian-Process-GBLUP",   # KRR is single-trait only
    multi_trait_gp = TRUE,
    gp_return_se = TRUE,
    gp_return_trait_correlations = TRUE,
    eval_metrics = c("root_mean_squared_error"),
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  )

  pred <- out$model_results$predicted_values
  # All trait x environment rows are predicted; held-out ones are labelled Test.
  expect_equal(nrow(pred), nrow(pheno) * 2L)
  expect_equal(sum(pred$Train_Test_Label == "Test"), length(test_ids) * length(dat$envs) * 2L)
  expect_true(all(is.finite(pred$Predicted_value)))
  expect_true(all(is.finite(pred$Standard_error)))
  expect_equal(pred$Standard_error^2, pred$PEV, tolerance = 1e-6)
  expect_true(all(c(
    "Genetic_correlation_traits",
    "Genetic_correlation_trait_environment"
  ) %in% names(out$model_results)))
  expect_gp_correlation_matrix(out$model_results$Genetic_correlation_traits, 2L)
  expect_gp_correlation_matrix(out$model_results$Genetic_correlation_trait_environment, 2L)
  expect_true(any(out$model_results$model_parameters$stat == "bridge_command" &
    out$model_results$model_parameters$summary == "fit-multi-trait-met-spec"))
  expect_true(PredictProR::validate_prediction_output(out)$valid)
})

test_that("model_execute MET GP CV2 uses direct GP point predictions only", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_EXECUTION = "direct",
    PREDICTPRO_GP_RETICULATE_FALLBACK = "false"
  )

  dat <- make_gp_public_smoke_data(n = 8L, p = 12L, seed = 20260505L)
  out <- PredictProR::model_execute(
    pheno_data = dat$pheno,
    gmatrix = dat$K,
    response = "Trait1",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = "KRR",
    cross_validation = TRUE,
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    gp_return_se = TRUE,
    eval_metrics = c("root_mean_squared_error"),
    cv_evaluation_only = TRUE,
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  )

  cv_pred <- out$cv_results_raw[[1L]]$ypred_cv_Reps_all
  cv_info <- out$cv_results_raw[[1L]]$cv_info
  test_rows <- cv_pred$cv_role == "test"
  expect_true(any(test_rows))
  expect_true(all(is.finite(cv_pred$yhat[test_rows])))
  expect_false(any(c("Prediction_SE", "Prediction_error_variance") %in% names(cv_pred)))
  expect_true(isTRUE(cv_info$gp_cv_frontdoor_fast_lane))
  expect_true(isTRUE(cv_info$package_cv_batch))
  expect_identical(cv_info$bridge_command, "fit-predict-cv-spec")
  meta <- as.data.frame(out$run_metadata, stringsAsFactors = FALSE)
  expect_true(any(meta$key == "gp_package_cv_batch" & meta$value == "TRUE"))
  expect_true(any(meta$key == "gp_bridge_command" & meta$value == "fit-predict-cv-spec"))
})

test_that("public GP R wrappers cover four scenarios with SE and correlations", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_EXECUTION = "direct",
    PREDICTPRO_GP_RETICULATE_FALLBACK = "false"
  )
  Sys.unsetenv("PREDICTPRO_GP_DEVICE")

  dat <- make_gp_public_smoke_data()
  test_ids <- tail(dat$gids, 2L)
  first_env <- dat$envs[[1L]]
  pheno_single_env <- dat$pheno[dat$pheno$Env == first_env, , drop = FALSE]

  single <- PredictProR::gp_single_trait_model(
    pheno_data = pheno_single_env,
    gmatrix = dat$K,
    response = "Trait1",
    gen_name = "GID",
    test_set = test_ids,
    gp_return_se = TRUE,
    python_bin = gp_python
  )
  expect_equal(nrow(single$predictions), length(test_ids))
  expect_true(all(is.finite(as.numeric(single$predictions$Prediction))))
  expect_gp_se_variance(single$predictions, "Prediction_SE_observed", "Prediction_Var_observed")
  expect_gp_direct_policy(single, "fit-predict-spec")

  multi_env <- PredictProR::gp_single_trait_model(
    pheno_data = dat$pheno,
    gmatrix = dat$K,
    response = "Trait1",
    gen_name = "GID",
    heter_groups = "Env",
    test_set = test_ids,
    env_similarity = dat$env_similarity,
    include_components = c("g", "ge", "e"),
    gp_return_se = TRUE,
    w_g = c(0.7),
    w_ge = c(0.2),
    w_e = 0.3,
    python_bin = gp_python
  )
  expect_equal(nrow(multi_env$predictions), length(test_ids) * length(dat$envs))
  expect_true(all(is.finite(as.numeric(multi_env$predictions$Prediction))))
  expect_gp_se_variance(multi_env$predictions, "Prediction_SE_observed", "Prediction_Var_observed")
  expect_gp_direct_policy(multi_env, "fit-predict-spec")

  mt_single <- PredictProR::gp_multi_trait_model(
    pheno_data = pheno_single_env,
    gmatrix = dat$K,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    test_set = test_ids,
    return_se = TRUE,
    prediction_output = "all",
    return_trait_correlations = TRUE,
    # Genetic covariances/correlations are only published for a converged
    # REML fit (otherwise status "not_promoted_nonconverged" and NA), so allow
    # enough iterations to converge.
    max_iter = 200L,
    tol_loglik = 1e-3,
    python_bin = gp_python
  )
  expect_equal(nrow(mt_single$predictions), nrow(pheno_single_env) * 2L)
  expect_true(all(is.finite(as.numeric(mt_single$predictions$Prediction))))
  expect_gp_multitrait_uncertainty_contract(mt_single$predictions)
  expect_gp_direct_policy(mt_single, "fit-multi-trait-spec")
  expect_true("genetic_correlation" %in% names(mt_single$result))
  expect_gp_correlation_matrix(mt_single$result$genetic_correlation, 2L)

  mt_multi_env <- PredictProR::gp_multi_trait_met_model(
    pheno_data = dat$pheno,
    gmatrix = dat$K,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    heter_groups = "Env",
    test_set = test_ids,
    env_similarity = dat$env_similarity,
    return_se = TRUE,
    prediction_output = "all",
    return_trait_correlations = TRUE,
    mom_pcg_max_iter = 150L,
    python_bin = gp_python
  )
  expect_equal(nrow(mt_multi_env$predictions), nrow(dat$pheno) * 2L)
  expect_true(all(is.finite(as.numeric(mt_multi_env$predictions$Prediction))))
  expect_gp_multitrait_uncertainty_contract(mt_multi_env$predictions)
  expect_gp_direct_policy(mt_multi_env, "fit-multi-trait-met-spec")
  expect_true(all(c("genetic_correlation", "gxe_correlation") %in% names(mt_multi_env$result)))
  expect_gp_correlation_matrix(mt_multi_env$result$genetic_correlation, 2L)
  expect_gp_correlation_matrix(mt_multi_env$result$gxe_correlation, 2L)

  K2 <- 0.65 * dat$K + 0.35 * diag(nrow(dat$K))
  dimnames(K2) <- dimnames(dat$K)
  mt_multi_kernel <- PredictProR::gp_multi_trait_met_model(
    pheno_data = dat$pheno,
    gmatrix = dat$K,
    kernel_list = list(omic = K2),
    kernel_weights = c(gmatrix = 0.7, omic = 0.3),   # kernels are named gmatrix/omic
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    heter_groups = "Env",
    test_set = test_ids,
    env_similarity = dat$env_similarity,
    return_se = FALSE,
    return_trait_correlations = FALSE,
    reml_max_iter = 4L,
    python_bin = gp_python
  )
  expect_equal(nrow(mt_multi_kernel$predictions), length(test_ids) * length(dat$envs) * 2L)
  expect_equal(as.numeric(mt_multi_kernel$fit$kernel_weights), c(0.7, 0.3), tolerance = 1e-8)
})

test_that("single-trait MET GP blocked tuning runs without using prediction IDs", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python
  )

  dat <- make_gp_public_smoke_data(n = 10L, p = 12L, seed = 20260502L)
  test_ids <- tail(dat$gids, 2L)
  env_covariates <- data.frame(
    Env = dat$envs,
    ec1 = c(0, 1),
    ec2 = c(1, 0),
    stringsAsFactors = FALSE
  )
  weights <- data.frame(w_g = 0.7, w_ge = 0.2, w_e = 0.3)

  tuned <- PredictProR::gp_single_trait_model(
    pheno_data = dat$pheno,
    gmatrix = dat$K,
    response = "Trait1",
    gen_name = "GID",
    heter_groups = "Env",
    test_set = test_ids,
    env_covariates = env_covariates,
    tune_gp = TRUE,
    tuning_krr_lams = 0.1,
    tuning_kenv_kernels = "matern32",
    tuning_weight_grid = weights,
    mean_adjustment = "none",
    python_bin = gp_python
  )

  expect_equal(nrow(tuned$predictions), length(test_ids) * length(dat$envs))
  expect_identical(tuned$info$gp_tuning$status, "ok")
  expect_false(any(tuned$info$gp_tuning$validation_ids %in% test_ids))
  expect_equal(tuned$info$gp_tuning$selected$lambda, 0.1, tolerance = 1e-8)
})

test_that("factor-cache MET GP tuning materializes bounded dense subset", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python
  )

  set.seed(20260503L)
  gids <- paste0("g", seq_len(12L))
  envs <- paste0("E", seq_len(4L))
  Phi <- matrix(rnorm(length(gids) * 6L), nrow = length(gids))
  Phi <- scale(Phi, center = TRUE, scale = FALSE) / sqrt(ncol(Phi))
  phi_path <- tempfile("predictpror_phi_", fileext = ".bin")
  con <- file(phi_path, open = "wb")
  writeBin(as.numeric(t(Phi)), con, size = 4L, endian = "little")
  close(con)
  factor_cache <- list(
    type = "memmap",
    paths = list(normalizePath(phi_path, winslash = "/", mustWork = TRUE)),
    shapes = list(as.list(c(nrow(Phi), ncol(Phi)))),
    dtypes = list("float32")
  )
  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gid_idx <- match(pheno$GID, gids)
  env_effect <- setNames(c(-0.2, 0.1, 0.25, 0.4), envs)[pheno$Env]
  pheno$Trait1 <- Phi[gid_idx, 1L] - 0.5 * Phi[gid_idx, 2L] + env_effect + rnorm(nrow(pheno), sd = 0.05)
  pheno$Trait1[pheno$Env == "E4"] <- NA_real_
  env_covariates <- data.frame(
    Env = envs,
    ec1 = seq_along(envs),
    ec2 = c(0.1, 0.4, 0.7, 1.0),
    stringsAsFactors = FALSE
  )

  pheno_df <- data.frame(
    GID = pheno$GID,
    Env = pheno$Env,
    y = pheno$Trait1,
    stringsAsFactors = FALSE
  )
  split <- PredictProR:::gp_public_train_test_idx(
    pheno_df,
    gid_col = "GID",
    y_col = "y"
  )
  kernel <- PredictProR:::gp_public_kernel_inputs(
    NULL,
    pheno_ids = pheno_df$GID,
    geno_ids = gids,
    allow_factor_cache = TRUE
  )

  tuned <- PredictProR:::gp_tune_single_trait_met(
    pheno_df = pheno_df,
    split = split,
    kernel = kernel,
    env_covariates = env_covariates,
    gp_factor_cache = factor_cache,
    geno_ids = gids,
    include_components = c("g", "ge", "e"),
    gp_backend = "auto",
    tuning_krr_lams = 0.1,
    tuning_kenv_kernels = "matern32",
    tuning_weight_grid = data.frame(w_g = 0.7, w_ge = 0.2, w_e = 0.3),
    tuning_max_calibration_rows = 18L,
    tuning_max_validation_rows = 8L,
    tuning_max_genotypes = 8L,
    mean_adjustment = "none",
    python_bin = gp_python
  )

  expect_identical(tuned$status, "ok")
  expect_true(isTRUE(tuned$batch_candidate_used))
  expect_identical(tuned$batch_candidate_execution, "python_subprocess")
  expect_lte(tuned$tuning_sample$n_genotypes_after, 8L)
  expect_equal(tuned$selected$lambda, 0.1, tolerance = 1e-8)
})

test_that("package GP CV batch route uses Python shared context cache", {
  gp_python <- test_gp_python()
  old_env <- Sys.getenv(
    c("PREDICTPRO_GP_PYTHON", "PREDICTPRO_GP_BATCH_CV", "PREDICTPRO_GP_FAST_CV", "PREDICTPRO_GP_EXACT_FAST_CV"),
    unset = NA_character_
  )
  on.exit({
    missing <- names(old_env)[is.na(old_env)]
    if (length(missing)) Sys.unsetenv(missing)
    present <- old_env[!is.na(old_env)]
    if (length(present)) do.call(Sys.setenv, as.list(present))
  }, add = TRUE)
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_BATCH_CV = "true",
    PREDICTPRO_GP_FAST_CV = "auto",
    PREDICTPRO_GP_EXACT_FAST_CV = "false"
  )

  set.seed(20260506L)
  gids <- paste0("g", seq_len(6L))
  envs <- paste0("e", seq_len(5L))
  X <- matrix(rnorm(length(gids) * 12L), nrow = length(gids))
  K <- tcrossprod(scale(X, center = TRUE, scale = FALSE)) / ncol(X)
  K <- K / mean(diag(K))
  diag(K) <- diag(K) + 1e-6
  rownames(K) <- colnames(K) <- gids

  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gid_idx <- match(pheno$GID, gids)
  env_idx <- match(pheno$Env, envs)
  pheno$Trait <- X[gid_idx, 1L] - 0.25 * X[gid_idx, 2L] + seq(-0.3, 0.4, length.out = length(envs))[env_idx] +
    rnorm(nrow(pheno), sd = 0.04)
  env_covariates <- data.frame(
    Env = envs,
    temp = seq(18, 28, length.out = length(envs)),
    rain = 100 + 20 * sin(seq(0, pi, length.out = length(envs))),
    soil = cos(seq(0, pi, length.out = length(envs))),
    stringsAsFactors = FALSE
  )
  folds <- list(
    which(pheno$Env == envs[[4L]]),
    which(pheno$Env == envs[[5L]])
  )

  preds <- suppressWarnings(PredictProR:::gp_backend_cv_predict_batch(
    model_name = "GP",
    y = pheno$Trait,
    folds = folds,
    additional_params = list(
      pheno_data = pheno,
      response = "Trait",
      gen_name = "GID",
      heter_groups = "Env",
      gmatrix_model_ready = K,
      env_covariates = env_covariates,
      reaction_norm_feature_qc = FALSE,
      kenv_kernel = "auto",
      include_components = c("g", "ge", "e"),
      gp_auto_policy = FALSE,
      gp_tiered_dispatch = FALSE,
      gp_output_level = "predict_only",
      gp_prediction_output = "test_only",
      gp_iters = 2L,
      gp_lr = 0.01,
      random_state = 20260506L
    )
  ))

  expect_length(preds, length(folds))
  expect_true(all(lengths(preds) == vapply(folds, length, integer(1L))))
  expect_true(all(is.finite(unlist(preds, use.names = FALSE))))
  info <- attr(preds, "gp_info", exact = TRUE)
  expect_true(is.list(info))
  expect_true(isTRUE(info$package_cv_batch))
  expect_identical(info$bridge_command, "fit-predict-cv-spec")
  expect_identical(info$cv_fast_path, "not_used")
  expect_identical(info$cv_fast_path_skip_reason, "gp_exact_fast_cv_disabled")
  expect_identical(info$cv_preprocessed_context, "used")
  expect_true(as.integer(info$cv_preprocessed_context_tensor_cache_entries) > 0L)
})

test_that("package GP exact CV fast path is R-controlled and visibly approximate", {
  gp_python <- test_gp_python()
  old_env <- Sys.getenv(
    c("PREDICTPRO_GP_PYTHON", "PREDICTPRO_GP_BATCH_CV", "PREDICTPRO_GP_FAST_CV", "PREDICTPRO_GP_EXACT_FAST_CV"),
    unset = NA_character_
  )
  on.exit({
    missing <- names(old_env)[is.na(old_env)]
    if (length(missing)) Sys.unsetenv(missing)
    present <- old_env[!is.na(old_env)]
    if (length(present)) do.call(Sys.setenv, as.list(present))
  }, add = TRUE)
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python,
    PREDICTPRO_GP_BATCH_CV = "true",
    PREDICTPRO_GP_FAST_CV = "auto",
    PREDICTPRO_GP_EXACT_FAST_CV = "false"
  )

  set.seed(20260507L)
  gids <- paste0("g", seq_len(12L))
  X <- matrix(rnorm(length(gids) * 18L), nrow = length(gids))
  K <- tcrossprod(scale(X, center = TRUE, scale = FALSE)) / ncol(X)
  K <- K / mean(diag(K))
  diag(K) <- diag(K) + 1e-6
  rownames(K) <- colnames(K) <- gids
  pheno <- data.frame(
    GID = gids,
    Trait = X[, 1L] - 0.4 * X[, 2L] + rnorm(length(gids), sd = 0.05),
    stringsAsFactors = FALSE
  )
  folds <- split(seq_len(nrow(pheno)), rep(seq_len(3L), length.out = nrow(pheno)))

  preds <- suppressWarnings(PredictProR:::gp_backend_cv_predict_batch(
    model_name = "GP",
    y = pheno$Trait,
    folds = folds,
    additional_params = list(
      pheno_data = pheno,
      response = "Trait",
      gen_name = "GID",
      gmatrix_model_ready = K,
      include_components = "g",
      gp_auto_policy = FALSE,
      gp_tiered_dispatch = FALSE,
      gp_output_level = "predict_only",
      gp_prediction_output = "test_only",
      gp_exact_fast_cv = "force",
      gp_iters = 2L,
      gp_lr = 0.01,
      random_state = 20260507L
    )
  ))

  expect_length(preds, length(folds))
  expect_true(all(is.finite(unlist(preds, use.names = FALSE))))
  info <- attr(preds, "gp_info", exact = TRUE)
  expect_true(is.list(info))
  expect_true(isTRUE(info$package_cv_batch))
  expect_identical(info$bridge_command, "fit-predict-cv-spec")
  expect_identical(info$gp_exact_fast_cv, "force")
  expect_identical(info$cv_fast_path, "gp_exact_shared_noise_block_delete")
  expect_true(isTRUE(info$cv_fast_path_approximate))
  expect_false(isTRUE(info$cv_fast_path_exact_refit))
  expect_identical(info$cv_fast_path_solver, "krr_block_delete")
  expect_true(is.finite(as.numeric(info$cv_fast_path_likelihood_noise_variance)))
  expect_equal(
    as.numeric(info$cv_fast_path_likelihood_noise_std),
    as.numeric(info$cv_fast_path_likelihood_noise_variance)
  )
})

test_that("batch KRR candidate fast path matches the single-candidate fit", {
  gp_python <- test_gp_python()
  Sys.setenv(
    PREDICTPRO_GP_PYTHON = gp_python
  )

  set.seed(20260502L)
  gids <- paste0("g", seq_len(8L))
  envs <- paste0("e", seq_len(4L))
  X <- matrix(rnorm(length(gids) * 5L), nrow = length(gids))
  K <- tcrossprod(X) / ncol(X)
  diag(K) <- diag(K) + 1e-6
  rownames(K) <- colnames(K) <- gids
  S <- matrix(
    c(1, 0.7, 0.4, 0.2,
      0.7, 1, 0.5, 0.3,
      0.4, 0.5, 1, 0.6,
      0.2, 0.3, 0.6, 1),
    nrow = 4L,
    byrow = TRUE
  )
  rownames(S) <- colnames(S) <- envs
  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gid_idx <- match(pheno$GID, gids)
  env_eff <- setNames(c(-0.2, 0.1, 0.35, -0.1), envs)[pheno$Env]
  pheno$y <- X[gid_idx, 1L] - 0.3 * X[gid_idx, 2L] + env_eff + rnorm(nrow(pheno), sd = 0.05)
  test_rows <- pheno$Env == envs[[4L]] & pheno$GID %in% gids[1:3]
  pheno$y[test_rows] <- NA_real_
  train_idx <- as.integer(which(!test_rows) - 1L)
  test_idx <- as.integer(which(test_rows) - 1L)

  candidate <- data.frame(
    candidate_id = 1L,
    lambda = 0.1,
    kenv_kernel = "linear",
    w_g = 0.8,
    w_ge = 0.1,
    w_e = 0.2,
    stringsAsFactors = FALSE
  )

  fast <- PredictProR:::gp_tuning_batch_candidates_direct(
    prepped = list(
      fit_df = pheno,
      geno_kernel = unname(K),
      ids = gids,
      training_idx = train_idx,
      validation_idx = test_idx
    ),
    candidate_df = candidate,
    env_covariates = NULL,
    env_kernel_cache = list(kernels = list(linear = unname(S)), detail = data.frame()),
    include_components = c("g", "ge", "e"),
    gp_backend = "auto",
    reaction_norm_feature_qc = TRUE,
    kenv_bandwidth = 1,
    kenv_kernel_kwargs = NULL,
    seed = 12345L,
    python_bin = gp_python,
    project_root = PredictProR:::gp_project_root()
  )
  single <- PredictProR:::gp_bridge_fit_mixed_model_direct(
    args = list(
      pheno_df = pheno,
      gid_col = "GID",
      env_col = "Env",
      y_col = "y",
      train_idx = train_idx,
      test_idx = test_idx,
      env_similarity = unname(S),
      include_components = as.list(c("g", "ge", "e")),
      fixed_effects = list(),
      method = "krr_exact",
      backend = "auto",
      output_level = "predict_only",
      point_predictions_only = TRUE,
      return_se = FALSE,
      compute_ai_se = FALSE,
      standardize = "global",
      varcomp_mode = "mom",
      krr_lam = 0.1,
      krr_lams = "none",
      lam_select = "fixed",
      dtype = "float64",
      seed = 12345L,
      w_g = 0.8,
      w_ge = 0.1,
      w_e = 0.2,
      prediction_output = "test_only"
    ),
    kernel = list(K = unname(K), geno_ids = gids),
    python_bin = gp_python,
    project_root = PredictProR:::gp_project_root()
  )

  expect_false(any(as.logical(fast$detail$failed)))
  expect_true("fast_path" %in% names(fast$detail))
  expect_identical(as.character(fast$detail$fast_path), "dense_krr_blocks")
  expect_equal(
    as.numeric(fast$predictions$Prediction),
    as.numeric(single$result$predictions$Prediction),
    tolerance = 1e-8
  )
})

test_that("joint multi-trait Scalable-GBLUP without PEV says so (not the ML/DL wording)", {
  gp_python <- test_gp_python()
  withr::local_envvar(c(PREDICTPRO_GP_PYTHON = gp_python))
  withr::local_dir(withr::local_tempdir())
  set.seed(20261007L)
  n <- 40L; p <- 120L
  G <- matrix(stats::rbinom(n * p, 2, 0.4), n, p,
              dimnames = list(paste0("g", seq_len(n)), paste0("m", seq_len(p))))
  G <- G[, apply(G, 2, stats::var) > 0]
  s <- as.numeric(scale(G[, 1:8] %*% stats::rnorm(8)))
  ph <- data.frame(GID = rownames(G), Yield = s + stats::rnorm(n, sd = 0.6),
                   Protein = 0.6 * s + stats::rnorm(n, sd = 0.6), stringsAsFactors = FALSE)
  ph$Yield[1:6] <- NA; ph$Protein[1:6] <- NA
  out <- PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = c("Yield", "Protein"),
    multi_trait_gp = TRUE, GS_model = "Scalable-GBLUP", gmatrix_method = "VanRaden",
    het_threshold = NULL, message = FALSE, verbose = FALSE
  )
  pv <- out$model_results$predicted_values
  no_rel <- !is.finite(pv$Reliability)
  if (any(no_rel)) {
    expect_false(any(grepl("ML/DL", pv$Reliability_basis[no_rel], fixed = TRUE)))
    expect_true(all(grepl("genetic model returned no prediction error variance",
                          pv$Reliability_basis[no_rel], fixed = TRUE)))
  } else {
    succeed()
  }
})