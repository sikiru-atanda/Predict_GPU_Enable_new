fake_asreml_model <- function(random, fixed) {
  structure(
    list(
      random = random,
      fixed = fixed,
      call = list(data = quote(pheno_data))
    ),
    class = "predictpror_fake_asreml"
  )
}

register_fake_asreml_methods <- function() {
  assign(
    "summary.predictpror_fake_asreml",
    function(object, coef = FALSE, ...) list(coef.random = object$random),
    envir = .GlobalEnv
  )
  assign(
    "coef.predictpror_fake_asreml",
    function(object, ...) list(fixed = object$fixed),
    envir = .GlobalEnv
  )
  withr::defer({
    for (method_name in c("summary.predictpror_fake_asreml", "coef.predictpror_fake_asreml")) {
      if (exists(method_name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = method_name, envir = .GlobalEnv)
      }
    }
  }, teardown_env())
}

test_that("ASReml fallback matches direct point predictions for single-trait single-environment multi-kernel output", {
  register_fake_asreml_methods()

  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(12, NA),
    stringsAsFactors = FALSE
  )
  random <- data.frame(
    solution = c(1, 0.5, 2, -0.5),
    `std.error` = c(0.10, 0.05, 0.20, 0.05),
    check.names = FALSE
  )
  rownames(random) <- c(
    "vm(GID,K1)_g1",
    "vm(GID,K2)_g1",
    "vm(GID,K1)_g2",
    "vm(GID,K2)_g2"
  )
  fixed <- data.frame(effect = 10, row.names = "(Intercept)", check.names = FALSE)
  model <- fake_asreml_model(random = random, fixed = fixed)

  old <- getOption("PredictProR.asreml_predict_impl")
  withr::defer(options(PredictProR.asreml_predict_impl = old))  # restore when this test ends

  direct_predict <- function(object, classify, sed, ...) {
    list(
      pvals = data.frame(
        GID = c("g1", "g2"),
        predicted.value = c(11.5, 11.5),
        std.error = c(0.15, 0.25),
        status = c("Estimable", "Estimable"),
        stringsAsFactors = FALSE
      )
    )
  }
  options(PredictProR.asreml_predict_impl = direct_predict)
  direct <- PredictProR:::asreml_predict_or_extract(
    mod = model,
    classify = "GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    names_in_inv_list = c("K1", "K2")
  )

  options(PredictProR.asreml_predict_impl = function(...) stop("forced ASReml predict failure"))
  fallback <- PredictProR:::asreml_predict_or_extract(
    mod = model,
    classify = "GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    names_in_inv_list = c("K1", "K2")
  )

  direct_pred <- direct$prediction[order(direct$prediction$GID), , drop = FALSE]
  fallback_pred <- fallback$prediction[order(fallback$prediction$GID), , drop = FALSE]

  expect_identical(direct$mode, "predict")
  expect_identical(fallback$mode, "extract")
  expect_equal(fallback_pred$Predicted_value, direct_pred$Predicted_value)
  expect_identical(fallback_pred$Train_Test_Label, c("Train", "Test"))
  expect_true(all(c("Standard_error", "Prediction_error_variance", "Reliability") %in% names(fallback_pred)))
  expect_true(all(is.na(fallback_pred$Standard_error)))
  expect_true(all(is.na(fallback_pred$Prediction_error_variance)))
  expect_true(all(is.na(fallback_pred$PEV)))
  expect_true(all(is.na(fallback_pred$Reliability)))
  expect_true(all(grepl("ASReml predict failed", fallback_pred$Reliability_remarks)))
})

test_that("ASReml fallback returns per-environment and across-environment point predictions for multi-kernel MET output", {
  register_fake_asreml_methods()

  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(11, NA, 9, 13),
    stringsAsFactors = FALSE
  )
  random <- data.frame(
    solution = c(
      0.1, 0.2, -0.4, 0.1,
      1.0, 0.5, 2.0, -0.5,
      -1.0, 0.25, 0.5, 0.25
    ),
    `std.error` = rep(0.05, 12),
    check.names = FALSE
  )
  rownames(random) <- c(
    "vm(GID,K1)_g1",
    "vm(GID,K2)_g1",
    "vm(GID,K1)_g2",
    "vm(GID,K2)_g2",
    "Env_E1:vm(GID,K1)_g1",
    "Env_E1:vm(GID,K2)_g1",
    "Env_E2:vm(GID,K1)_g1",
    "Env_E2:vm(GID,K2)_g1",
    "Env_E1:vm(GID,K1)_g2",
    "Env_E1:vm(GID,K2)_g2",
    "Env_E2:vm(GID,K1)_g2",
    "Env_E2:vm(GID,K2)_g2"
  )
  fixed <- data.frame(
    effect = c(10, 3),
    row.names = c("(Intercept)", "Env_E2"),
    check.names = FALSE
  )
  model <- fake_asreml_model(random = random, fixed = fixed)

  expected <- data.frame(
    Env = c("E1", "E2", "E1", "E2"),
    GID = c("g1", "g1", "g2", "g2"),
    Predicted_value = c(11.8, 14.8, 8.95, 13.45),
    stringsAsFactors = FALSE
  )

  old <- getOption("PredictProR.asreml_predict_impl")
  withr::defer(options(PredictProR.asreml_predict_impl = old))  # restore when this test ends
  options(PredictProR.asreml_predict_impl = function(object, classify, sed, ...) {
    list(
      pvals = data.frame(
        Env = expected$Env,
        GID = expected$GID,
        predicted.value = expected$Predicted_value,
        std.error = rep(0.10, nrow(expected)),
        status = rep("Estimable", nrow(expected)),
        check.names = FALSE
      )
    )
  })
  direct <- PredictProR:::asreml_predict_or_extract(
    mod = model,
    classify = "Env:GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    names_in_inv_list = c("K1", "K2")
  )

  options(PredictProR.asreml_predict_impl = function(...) stop("forced ASReml predict failure"))
  fallback <- PredictProR:::asreml_predict_or_extract(
    mod = model,
    classify = "Env:GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    names_in_inv_list = c("K1", "K2")
  )

  direct_env <- direct$across_env_prediction[order(direct$across_env_prediction$GID, direct$across_env_prediction$Env), , drop = FALSE]
  fallback_env <- fallback$across_env_prediction[order(fallback$across_env_prediction$GID, fallback$across_env_prediction$Env), , drop = FALSE]
  fallback_across <- fallback$prediction[order(fallback$prediction$GID), , drop = FALSE]

  expect_identical(fallback$mode, "extract")
  expect_equal(fallback_env$Predicted_value, direct_env$Predicted_value)
  expect_equal(fallback_across$Predicted_value, c(mean(c(11.8, 14.8)), mean(c(8.95, 13.45))))
  expect_true(all(c("Standard_error", "Prediction_error_variance", "Reliability") %in% names(fallback_env)))
  expect_true(all(is.na(fallback_env$Standard_error)))
  expect_true(all(is.na(fallback_env$Prediction_error_variance)))
  expect_true(all(is.na(fallback_env$PEV)))
  expect_true(all(is.na(fallback_env$Reliability)))
  expect_true(all(grepl("ASReml predict failed", fallback_env$Reliability_remarks)))
  expect_true(all(is.na(fallback_across$Standard_error)))
  expect_true(all(grepl("ASReml predict failed", fallback_across$Reliability_remarks)))
  expect_identical(
    fallback_env$Train_Test_Label[match(paste("g1", "E2"), paste(fallback_env$GID, fallback_env$Env))],
    "Test"
  )
})

test_that("ASReml CV fallback returns point predictions without uncertainty", {
  register_fake_asreml_methods()

  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(12, 9),
    stringsAsFactors = FALSE
  )
  random <- data.frame(
    solution = c(1, 0.5, 2, -0.5),
    `std.error` = c(0.10, 0.05, 0.20, 0.05),
    check.names = FALSE
  )
  rownames(random) <- c(
    "vm(GID,K1)_g1",
    "vm(GID,K2)_g1",
    "vm(GID,K1)_g2",
    "vm(GID,K2)_g2"
  )
  fixed <- data.frame(effect = 10, row.names = "(Intercept)", check.names = FALSE)
  model <- fake_asreml_model(random = random, fixed = fixed)

  local_mocked_bindings(
    asreml_cv_model = function(pheno_dataa, response, tst, ...) {
      pheno_dataa[tst, response] <- NA
      list(model_cv = model, pheno_dataa = pheno_dataa)
    },
    .package = "PredictProR"
  )

  old <- getOption("PredictProR.asreml_predict_impl")
  withr::defer(options(PredictProR.asreml_predict_impl = old))  # restore when this test ends
  options(PredictProR.asreml_predict_impl = function(...) stop("forced ASReml predict failure"))

  pred <- PredictProR::asreml_mod_cv(
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = NULL,
    asreml_models_prep_cv = list(
      code_asr_fit = c("", "random = ~ vm(GID,K1) + vm(GID,K2)", "", ""),
      gen_pos = 1L,
      inter_gen_pos = NULL,
      names_in_inv_list = c("K1", "K2"),
      rand_term = "GID"
    ),
    response = "Yield",
    tst = 2L
  )

  expect_equal(pred, 11.5)
})

test_that("ASReml coefficient fallback refuses non-converged fitted models", {
  register_fake_asreml_methods()

  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, NA),
    stringsAsFactors = FALSE
  )
  random <- data.frame(
    solution = c(1, 2),
    `std.error` = c(0.1, 0.2),
    check.names = FALSE
  )
  rownames(random) <- c("vm(GID,K1)_g1", "vm(GID,K1)_g2")
  fixed <- data.frame(effect = 10, row.names = "(Intercept)", check.names = FALSE)
  model <- fake_asreml_model(random = random, fixed = fixed)
  model$converge <- FALSE

  old <- getOption("PredictProR.asreml_predict_impl")
  withr::defer(options(PredictProR.asreml_predict_impl = old))  # restore when this test ends
  options(PredictProR.asreml_predict_impl = function(...) stop("forced ASReml predict failure"))

  fallback <- PredictProR:::asreml_predict_or_extract(
    mod = model,
    classify = "GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    names_in_inv_list = "K1"
  )

  expect_identical(fallback$mode, "legacy_fallback")
  expect_null(fallback$prediction)
  expect_match(fallback$error, "did not converge")
})

test_that("ASReml CV never reaches the legacy fallback after a nonconverged fit", {
  pheno <- data.frame(GID = c("g1", "g2"), Yield = c(1, NA_real_))
  local_mocked_bindings(
    asreml_cv_model = function(...) {
      list(model_cv = list(converge = FALSE), pheno_dataa = pheno)
    },
    asreml_predict_or_extract = function(...) {
      list(mode = "legacy_fallback", prediction = NULL,
           across_env_prediction = NULL, error = "Insufficient pworkspace")
    },
    asreml_predict_pvals = function(...) stop("legacy predict must not run"),
    .package = "PredictProR"
  )
  expect_error(
    PredictProR::asreml_mod_cv(
      pheno_data = pheno, gen_name = "GID", heter_groups = NULL,
      asreml_models_prep_cv = list(
        code_asr_fit = character(), gen_pos = 1L,
        inter_gen_pos = NULL, names_in_inv_list = "Kinv", rand_term = "GID"
      ),
      response = "Yield", tst = 2L
    ),
    "fit did not converge"
  )
})

test_that("multi-trait ASReml fallback matches direct point predictions for multi-kernel joint output", {
  register_fake_asreml_methods()

  random <- data.frame(
    solution = c(1, 0.5, 2, -0.5, -1, 0.25, 0.75, 0.25),
    `std.error` = rep(0.05, 8),
    check.names = FALSE
  )
  rownames(random) <- c(
    "us(Trait):vm(GID,K1)_T1_g1",
    "us(Trait):vm(GID,K2)_T1_g1",
    "us(Trait):vm(GID,K1)_T1_g2",
    "us(Trait):vm(GID,K2)_T1_g2",
    "us(Trait):vm(GID,K1)_T2_g1",
    "us(Trait):vm(GID,K2)_T2_g1",
    "us(Trait):vm(GID,K1)_T2_g2",
    "us(Trait):vm(GID,K2)_T2_g2"
  )
  fixed <- data.frame(
    effect = c(10, 3),
    row.names = c("(Intercept)", "Trait_T2"),
    check.names = FALSE
  )
  model <- fake_asreml_model(random = random, fixed = fixed)
  expected <- data.frame(
    Trait = c("T1", "T1", "T2", "T2"),
    GID = c("g1", "g2", "g1", "g2"),
    Predicted_value = c(11.5, 11.5, 12.25, 14.0),
    stringsAsFactors = FALSE
  )

  direct <- PredictProR:::gp_multitrait_asreml_predict_or_extract(
    mod = model,
    predict_fn = function(object, classify, sed, ...) {
      list(
        pvals = data.frame(
          expected,
          predicted.value = expected$Predicted_value,
          std.error = rep(0.10, nrow(expected)),
          status = rep("Estimable", nrow(expected)),
          check.names = FALSE
        )[, c("Trait", "GID", "predicted.value", "std.error", "status")]
      )
    },
    classify = "Trait:GID",
    response = c("T1", "T2"),
    gen_name = "GID",
    kernel_ids = c("K1", "K2")
  )
  fallback <- PredictProR:::gp_multitrait_asreml_predict_or_extract(
    mod = model,
    predict_fn = function(...) stop("forced ASReml predict failure"),
    classify = "Trait:GID",
    response = c("T1", "T2"),
    gen_name = "GID",
    kernel_ids = c("K1", "K2")
  )

  direct_pred <- direct$prediction[order(direct$prediction$Trait, direct$prediction$GID), , drop = FALSE]
  fallback_pred <- fallback$prediction[order(fallback$prediction$Trait, fallback$prediction$GID), , drop = FALSE]

  expect_identical(direct$mode, "predict")
  expect_identical(fallback$mode, "extract")
  expect_equal(fallback_pred$Predicted_value, direct_pred$Predicted_value)
  expect_true(all(c("Standard_error", "Prediction_error_variance") %in% names(fallback_pred)))
  expect_true(all(is.na(fallback_pred$Standard_error)))
  expect_true(all(is.na(fallback_pred$Prediction_error_variance)))
  expect_true(all(grepl("ASReml predict failed", fallback_pred$Reliability_remarks)))
})

test_that("multi-trait ASReml fallback refuses aliased coefficient tables", {
  register_fake_asreml_methods()

  random <- data.frame(
    solution = c(1, 2),
    `std.error` = c(NA_real_, 0.2),
    check.names = FALSE
  )
  rownames(random) <- c(
    "Trait_T1:vm(GID,K1)_g1",
    "Trait_T1:vm(GID,K1)_g2"
  )
  fixed <- data.frame(
    effect = 10,
    row.names = "Trait_T1",
    check.names = FALSE
  )
  model <- fake_asreml_model(random = random, fixed = fixed)

  fallback <- PredictProR:::gp_multitrait_asreml_predict_or_extract(
    mod = model,
    predict_fn = function(...) stop("forced ASReml predict failure"),
    classify = "Trait:GID",
    response = "T1",
    gen_name = "GID",
    kernel_ids = "K1"
  )

  expect_identical(fallback$mode, "predict_error")
  expect_null(fallback$prediction)
  expect_match(fallback$error, "unsafe")
})

test_that("ASReml fallback reconstructs per-environment predictions for fa()/rr() MET terms", {
  # FA/RR terms are fitted as fa(Env, k):vm(GID, K); ASReml names the
  # per-environment genetic effects `fa(Env, 1)_<env>:vm(GID, K)_<gid>` and
  # the latent scores `fa(Env, 1)_Comp1:vm(GID, K)_<gid>`. The extractor used
  # to match neither, so every FA predict failure fell through to the
  # unflagged legacy path.
  register_fake_asreml_methods()
  for (struct in c("fa", "rr")) {
    random <- data.frame(
      solution = c(0.5, -0.2, 0.8, 0.1, 9, 9),
      `std.error` = rep(0.1, 6),
      check.names = FALSE
    )
    rownames(random) <- c(
      sprintf("%s(Env, 1)_E1:vm(GID, K)_g1", struct),
      sprintf("%s(Env, 1)_E2:vm(GID, K)_g1", struct),
      sprintf("%s(Env, 1)_E1:vm(GID, K)_g2", struct),
      sprintf("%s(Env, 1)_E2:vm(GID, K)_g2", struct),
      sprintf("%s(Env, 1)_Comp1:vm(GID, K)_g1", struct),
      sprintf("%s(Env, 1)_Comp1:vm(GID, K)_g2", struct)
    )
    fixed <- data.frame(effect = c(10, 0, 2), row.names = c("(Intercept)", "Env_E1", "Env_E2"))
    out <- PredictProR:::extract_asreml_prediction(
      fake_asreml_model(random = random, fixed = fixed),
      genotype_term = "GID",
      environment_levels = c("E1", "E2"),
      heter_grp = "Env"
    )
    out <- out[order(out$Genotype, out$Environment), ]
    expect_equal(out$Predicted, c(10.5, 11.8, 10.8, 12.1), info = struct)
    expect_equal(out$Environment, c("E1", "E2", "E1", "E2"), info = struct)
  }
})
test_that("ASReml coefficient rows are labelled by parsed environment, not position", {
  rn <- c(
    "vm(GID, K)_g1",
    "Env_Zeta:vm(GID, K)_g1",
    "vm(GID, K)_g2:Env_Alpha",
    "fa(Env, 1)_Mid:vm(GID, K)_g2",
    "fa(Env, 1)_Comp1:vm(GID, K)_g2",
    "corgh(Env)_Alpha:vm(GID, K)_g1"
  )
  expect_identical(
    PredictProR:::asreml_coef_row_env(rn, "Env"),
    c(NA, "Zeta", "Alpha", "Mid", "Comp1", "Alpha")
  )
})

test_that("legacy MET fallback sums kernels and adds the genotype main effect to every environment", {
  df <- data.frame(
    GID = c("g1", "g1", "g1", "g1", "g2", "g2"),
    Env = c(NA, "Zeta", "Alpha", "Zeta", NA, "Alpha"),
    BLUP = c(1, 0.5, -0.2, 0.25, 2, 0.3),
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::asreml_legacy_met_genetic_values(df, "GID", "Env", c("Zeta", "Alpha"))
  out <- out[order(out$GID, out$Env), ]
  expect_equal(out$Env, c("Alpha", "Zeta", "Alpha", "Zeta"))
  expect_equal(out$BLUP, c(0.8, 1.75, 2.3, 2))
})