test_that("Stage 2 precision weights resolve from supported user inputs", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    predicted.value = c(1.2, 2.3, NA_real_),
    Weight = c(1, 4, NA_real_),
    stringsAsFactors = FALSE
  )

  by_name <- PredictProR:::gp_resolve_stage2_precision_weights(
    weights = "Weight",
    pheno_data = pheno,
    response = "predicted.value",
    gen_name = "GID"
  )
  expect_true(by_name$supplied)
  expect_equal(by_name$precision, c(1, 4, 1))
  expect_identical(by_name$source, "pheno_data$Weight")

  by_table <- PredictProR:::gp_resolve_stage2_precision_weights(
    weights = data.frame(Weight = c(1, 4, 9)),
    pheno_data = transform(pheno, predicted.value = c(1.2, 2.3, 3.4)),
    response = "predicted.value"
  )
  expect_equal(by_table$precision, c(1, 4, 9))

  expect_error(
    PredictProR:::gp_resolve_stage2_precision_weights(
      weights = c(1, 0, 3),
      pheno_data = transform(pheno, predicted.value = c(1.2, 2.3, 3.4)),
      response = "predicted.value"
    ),
    "finite weight greater than zero"
  )
})

test_that("model_execute weight attachment preserves row alignment", {
  pheno <- data.frame(
    GID = c("g2", "g1", "g3"),
    predicted.value = c(2, 1, 3),
    Weight = c(4, 1, 9),
    stringsAsFactors = FALSE
  )
  attached <- PredictProR:::gp_attach_stage2_precision_weights(
    pheno_data = pheno,
    weights = "Weight",
    response = "predicted.value",
    gen_name = "GID"
  )
  col <- PredictProR:::gp_stage2_weight_column()
  expect_identical(attached$weights, col)
  expect_equal(attached$pheno_data[[col]], c(4, 1, 9))

  reordered <- attached$pheno_data[order(attached$pheno_data$GID), , drop = FALSE]
  expect_equal(reordered[[col]], c(1, 4, 9))
})

test_that("GP and ASReml use inverse-variance precision semantics", {
  clause <- PredictProR:::gp_asreml_weighted_data_clause(
    data_expr = "stage2_input",
    weight_column = "Weight",
    response_family = "gaussian"
  )
  expect_match(clause, "weights = Weight", fixed = TRUE)
  expect_match(clause, "asr_gaussian(dispersion = 1)", fixed = TRUE)

  args <- list(
    pheno_df = data.frame(GID = c("g1", "g2"), Env = "E1", y = c(1, 2)),
    gid_col = "GID",
    env_col = "Env",
    y_col = "y",
    train_idx = 0:1,
    test_idx = integer(),
    obs_weights = c(1, 4),
    obs_weight_mode = "inverse_variance",
    obs_weight_global_scale = 1
  )
  kernel <- list(
    K = diag(2),
    Ks = list(A = diag(2)),
    geno_ids = c("g1", "g2")
  )
  spec <- PredictProR:::gp_bridge_prepare_mixed_model_spec(
    args = args,
    kernel = kernel,
    dir_path = tempfile("stage2_weight_spec_")
  )
  payload <- jsonlite::fromJSON(spec$spec_json, simplifyVector = TRUE)
  expect_equal(as.numeric(payload$fit$obs_weights), c(1, 4))
  expect_identical(payload$fit$obs_weight_mode, "inverse_variance")
  expect_equal(payload$fit$obs_weight_global_scale, 1)
})

test_that("BGLR receives sqrt precision rather than squared precision", {
  seen <- NULL
  local_mocked_bindings(
    gp_bayes_bglr_fit = function(...) {
      seen <<- list(...)
      list(weights = seen$weights)
    },
    gp_bglr_save_prefix = function(...) tempfile("stage2_bglr_"),
    .package = "PredictProR"
  )

  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    predicted.value = c(1, 2, 3),
    Weight = c(1, 4, 9),
    stringsAsFactors = FALSE
  )
  out <- PredictProR::bayes_mod_execute(
    pheno_data = pheno,
    response = "predicted.value",
    weights = "Weight",
    ETA = list(),
    bayes_para = list(nIter = 20L, burnIn = 5L, thin = 1L),
    response_family = "gaussian"
  )
  expect_equal(seen$weights, c(1, 2, 3))
  expect_equal(out$model$weights, c(1, 2, 3))
})

test_that("Bayesian kernel outputs expose Stage 2 weight provenance", {
  skip_if_not_installed("BGLR")

  set.seed(902)
  pheno <- expand.grid(
    GID = paste0("g", 1:4),
    Env = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- stats::rnorm(nrow(pheno))
  pheno$Yield[pheno$GID == "g4"] <- NA_real_
  pheno$Weight <- seq(1, 4, length.out = nrow(pheno))
  K <- diag(4)
  rownames(K) <- colnames(K) <- paste0("g", 1:4)

  fit <- PredictProR::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "GBLUP_BRR",
    response = "Yield",
    weights = "Weight",
    pheno_data = pheno,
    gmatrix = K,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = FALSE,
    nIter = 80L,
    burnIn = 20L,
    thin = 2L
  )

  params <- fit$bayes_result$model_parameters
  values <- stats::setNames(as.character(params$summary), params$stat)
  expect_identical(values[["bayesian_kernel_model"]], "GBLUP_BRR")
  expect_identical(
    values[["stage2_observation_weighting"]],
    "precision_via_BGLR_inverse_squared_weight"
  )
  expect_identical(values[["stage2_weight_supplied"]], "TRUE")
  expect_identical(values[["stage2_weight_source"]], "pheno_data$Weight")
  expect_identical(values[["bglr_native_weight_transform"]], "sqrt(user_precision)")
  expect_equal(as.numeric(values[["stage2_precision_min"]]), min(pheno$Weight))
  expect_equal(as.numeric(values[["stage2_precision_max"]]), max(pheno$Weight))
})

test_that("weighted unsupported routes fail instead of ignoring weights", {
  expect_error(
    PredictProR:::gp_validate_stage2_weight_route(list(
      weights = c(1, 2),
      multi_trait_gp = TRUE,
      GS_model = "GP",
      GS_model_cv = NULL,
      gp_valid_models = "GP",
      bayes_valid_models = character(),
      bayes_gblup_valid_models = character(),
      asreml_model = "GBLUP"
    )),
    "multi_trait_gp"
  )

  expect_error(
    PredictProR:::gp_validate_stage2_weight_route(list(
      weights = c(1, 2),
      GS_model = "RandomForest",
      GS_model_cv = NULL,
      gp_valid_models = c("KRR", "GP"),
      bayes_valid_models = c("BayesA", "BRR"),
      bayes_gblup_valid_models = c("RKHS", "GBLUP_BRR"),
      asreml_model = "GBLUP"
    )),
    "RandomForest"
  )
})

test_that("weighted heterogeneous RKHS fails early without changing GBLUP_BRR validation", {
  base_ctx <- list(
    weights = c(1, 2),
    multi_trait_gp = FALSE,
    multi_trait_bayes = FALSE,
    multi_trait_asreml = FALSE,
    multi_trait_ml = FALSE,
    multi_trait_dl = FALSE,
    met_ml_dl = FALSE,
    hybrid_gp = FALSE,
    hybrid_bayes = FALSE,
    hybrid_asreml = FALSE,
    hybrid_ml = FALSE,
    hybrid_dl = FALSE,
    GS_model_cv = NULL,
    gp_valid_models = c("KRR", "GP"),
    bayes_valid_models = c("BayesA", "BRR"),
    bayes_gblup_valid_models = c("RKHS", "GBLUP_BRR"),
    asreml_model = "GBLUP",
    heter_resid = TRUE,
    bayes_kernel_heter_resid = NULL
  )

  expect_error(
    PredictProR:::gp_validate_stage2_weight_route(
      utils::modifyList(base_ctx, list(GS_model = "RKHS"))
    ),
    "RKHS cannot combine Stage 2 observation weights with heterogeneous-environment mode",
    fixed = TRUE
  )
  expect_silent(
    PredictProR:::gp_validate_stage2_weight_route(
      utils::modifyList(base_ctx, list(GS_model = "GBLUP_BRR"))
    )
  )
  expect_silent(
    PredictProR:::gp_validate_stage2_weight_route(
      utils::modifyList(
        base_ctx,
        list(GS_model = "RKHS", bayes_kernel_heter_resid = FALSE)
      )
    )
  )
})
