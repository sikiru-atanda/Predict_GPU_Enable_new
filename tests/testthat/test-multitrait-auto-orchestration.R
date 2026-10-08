test_that("automatic multi-trait planning preserves model-specific routes", {
  expect_true("multi_trait" %in% names(formals(PredictProR::model_execute)))

  plan <- PredictProR:::gp_multitrait_auto_plan(
    models = c(
      "GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP",
      "GBLUP_BRR", "RKHS", "RandomForest", "DenseNeuralNet"
    ),
    cross_validation = TRUE,
    is_met = FALSE
  )

  expect_identical(
    plan$route,
    c(
      "multi_trait_asreml", rep("multi_trait_gp", 3L),
      rep("multi_trait_bayes", 2L), "multi_trait_ml", "multi_trait_dl"
    )
  )
  expect_identical(
    plan$model,
    c(
      "GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP",
      "GBLUP_BRR", "RKHS", "RandomForest", "DenseNeuralNet"
    )
  )
})

test_that("multi-trait Scalable-GBLUP keeps genomic relationship preparation enabled", {
  ordinary_kernel_models <- PredictProR:::gp_multi_environment_kernel_models(
    gp_valid_models = PredictProR:::gp_lowrank_supported_models()
  )
  expect_false("LowRankGP" %in% ordinary_kernel_models)

  multi_trait_kernel_models <- PredictProR:::gp_kernel_models_for_input_preparation(
    base_models = ordinary_kernel_models,
    multi_trait_gp = TRUE
  )
  expect_true(all(PredictProR:::gp_multitrait_gp_supported_models() %in% multi_trait_kernel_models))

  ordinary_input_models <- PredictProR:::gp_kernel_models_for_input_preparation(
    base_models = ordinary_kernel_models,
    multi_trait_gp = FALSE
  )
  expect_identical(ordinary_input_models, ordinary_kernel_models)
})

multitrait_auto_test_args <- function(cv_evaluation_only = TRUE) {
  pheno <- data.frame(
    GID = paste0("G", 1:6),
    Trait1 = seq_len(6),
    Trait2 = seq_len(6) + 0.5,
    stringsAsFactors = FALSE
  )
  geno <- matrix(
    rep(c(0, 1, 2), length.out = 24),
    nrow = 6,
    dimnames = list(pheno$GID, paste0("M", 1:4))
  )
  list(
    multi_trait = TRUE,
    multi_trait_asreml = FALSE,
    multi_trait_gp = FALSE,
    multi_trait_bayes = FALSE,
    multi_trait_ml = FALSE,
    multi_trait_dl = FALSE,
    hybrid_asreml = FALSE,
    hybrid_bayes = FALSE,
    hybrid_gp = FALSE,
    hybrid_ml = FALSE,
    hybrid_dl = FALSE,
    met_ml_dl = FALSE,
    pheno_data = pheno,
    pheno_data_train = NULL,
    pheno_data_test = NULL,
    geno_data = geno,
    gmatrix = NULL,
    gkernel = NULL,
    kernel_list = NULL,
    gmatrix_method = NULL,
    gen_name = "GID",
    response = c("Trait1", "Trait2"),
    response_family = "gaussian",
    heter_groups = NULL,
    GS_model = NULL,
    GS_model_cv = c("RandomForest", "Xgboost"),
    cross_validation = TRUE,
    cross_validation_meth = "Repeated_K-Folds",
    sampling_method = "unstratified",
    cv_evaluation_only = cv_evaluation_only,
    metric_for_ranking = "accuracy",
    ranking_tie_breakers = "root_mean_squared_error",
    engine = NULL,
    plot_filename = "trait"
  )
}

test_that("automatic multi-trait CV runs one protected route per model", {
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    model <- call$GS_model_cv
    score <- if (identical(model, "Xgboost")) 0.8 else 0.6
    list(
      model_results = NULL,
      cv_results_raw = list(model = model),
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = rep(c("Trait1", "Trait2"), each = 2L),
          model = model,
          metric = rep(c("accuracy", "root_mean_squared_error"), 2L),
          value = rep(c(score, 1 - score), 2L),
          stringsAsFactors = FALSE
        )
      )
    )
  }

  out <- PredictProR:::gp_model_execute_multitrait_orchestrate(
    multitrait_auto_test_args(cv_evaluation_only = TRUE),
    execute_fn = fake_execute
  )

  expect_s3_class(out, "PredictProR_multi_trait_result")
  expect_length(calls, 2L)
  expect_true(calls[[1L]]$multi_trait_ml)
  expect_false(calls[[1L]]$multi_trait_asreml)
  expect_false(calls[[1L]]$multi_trait_dl)
  expect_identical(calls[[1L]]$GS_model_cv, "RandomForest")
  expect_identical(calls[[2L]]$GS_model_cv, "Xgboost")
  expect_true(all(vapply(calls, `[[`, logical(1L), "cv_evaluation_only")))
  expect_identical(out$selected_model, "Xgboost")
  expect_equal(nrow(out$cv_results_processed$multitrait_metric_summary), 8L)
})

test_that("automatic multi-trait CV sanitizes shared inputs for every route", {
  args <- multitrait_auto_test_args(cv_evaluation_only = TRUE)
  args$GS_model_cv <- c(
    "GBLUP", "GBLUP_BRR", "Gaussian-Process-GBLUP",
    "RandomForest", "DenseNeuralNet"
  )
  args$response <- rep(c("Trait1", "Trait2"), length(args$GS_model_cv))
  args$cross_validation_meth <- "Repeated_K-Folds"
  args$sampling_method <- "stratified"
  args$message <- FALSE
  calls <- list()

  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    list(
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = call$GS_model_cv,
          metric = "accuracy",
          value = 0.5,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_length(calls, 5L)
  expect_true(all(vapply(
    calls,
    function(call) identical(call$response, c("Trait1", "Trait2")),
    logical(1L)
  )))
  expect_true(all(vapply(
    calls,
    function(call) identical(call$cross_validation_meth, "Repeated_K-Folds"),
    logical(1L)
  )))
  expect_true(all(vapply(
    calls,
    function(call) identical(call$sampling_method, "unstratified"),
    logical(1L)
  )))
})

test_that("automatic multi-trait CV method aliases are canonicalized once", {
  expect_identical(
    PredictProR:::gp_multitrait_auto_cv_method("Repeated_K-Folds"),
    "Repeated_K-Folds"
  )
  expect_identical(
    PredictProR:::gp_multitrait_auto_cv_method("repeated k folds"),
    "Repeated_K-Folds"
  )
  expect_identical(
    PredictProR:::gp_multitrait_auto_cv_method("k_folds"),
    "K-Folds"
  )
  expect_error(
    PredictProR:::gp_multitrait_auto_cv_method("Repeated_Stratified_K-Folds"),
    "received.*normalized token"
  )
  expect_identical(
    PredictProR:::gp_multitrait_auto_cv_method("cv0", is_met = TRUE),
    "CV0"
  )
  expect_identical(
    PredictProR:::gp_multitrait_auto_cv_method("repeated cv2", is_met = TRUE),
    "Repeated_CV2"
  )
  expect_error(
    PredictProR:::gp_multitrait_auto_cv_method("K-Folds", is_met = TRUE),
    "MT-MET.*CV0"
  )
})

test_that("automatic multi-trait CV can refit the best joint model", {
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    model <- if (isTRUE(call$cross_validation)) call$GS_model_cv else call$GS_model
    if (!isTRUE(call$cross_validation)) {
      return(list(model_results = list(selected_model = model)))
    }
    score <- if (identical(model, "Xgboost")) 0.9 else 0.5
    list(
      model_results = NULL,
      cv_results_raw = list(model = model),
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = model,
          metric = "accuracy",
          value = score,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  args <- multitrait_auto_test_args(cv_evaluation_only = FALSE)
  args$ranking_tie_breakers <- NULL
  out <- PredictProR:::gp_model_execute_multitrait_orchestrate(
    args,
    execute_fn = fake_execute
  )

  expect_length(calls, 3L)
  expect_identical(out$selected_model, "Xgboost")
  expect_false(calls[[3L]]$cross_validation)
  expect_identical(calls[[3L]]$GS_model, "Xgboost")
  expect_null(calls[[3L]]$GS_model_cv)
  expect_true(calls[[3L]]$multi_trait_ml)
  expect_true(all(vapply(
    calls,
    function(call) isTRUE(call$.predictpror_shared_final_prediction),
    logical(1L)
  )))
  expect_true(out$run_metadata$preprocessing_detail$requires_final_prediction)
  expect_identical(out$final_prediction$model_results$selected_model, "Xgboost")
})

test_that("automatic MT-MET CV2 refits the selected model by direct prediction", {
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    model <- if (isTRUE(call$cross_validation)) call$GS_model_cv else call$GS_model
    if (!isTRUE(call$cross_validation)) {
      return(list(model_results = list(selected_model = model)))
    }
    score <- if (identical(model, "Gaussian-Process-GBLUP")) 0.8 else 0.4
    list(
      cv_results_raw = list(model = model),
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = model,
          metric = "accuracy",
          value = score,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  args <- multitrait_auto_test_args(cv_evaluation_only = FALSE)
  base_pheno <- args$pheno_data
  args$pheno_data <- merge(
    base_pheno,
    data.frame(Env = c("E1", "E2", "E3"), .merge_key = 1L),
    by = NULL,
    sort = FALSE
  )
  args$heter_groups <- "Env"
  args$GS_model_cv <- c("GBLUP_BRR", "Gaussian-Process-GBLUP")
  args$cross_validation_meth <- "CV2"
  args$ranking_tie_breakers <- NULL

  out <- PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_length(calls, 3L)
  expect_true(all(vapply(calls[1:2], function(x) {
    identical(x$cross_validation_meth, "CV2") &&
      isTRUE(x$cross_validation) && isTRUE(x$cv_evaluation_only)
  }, logical(1L))))
  expect_identical(out$selected_model, "Gaussian-Process-GBLUP")
  expect_false(calls[[3L]]$cross_validation)
  expect_true(calls[[3L]]$multi_trait_gp)
  expect_identical(calls[[3L]]$GS_model, "Gaussian-Process-GBLUP")
  expect_identical(
    out$final_prediction$model_results$selected_model,
    "Gaussian-Process-GBLUP"
  )
})

test_that("automatic CV and final prediction share one preprocessing cache", {
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    cache <- call$.predictpror_preprocess_cache
    if (!PredictProR:::gp_shared_preprocess_cache_has(cache, "model_input_objects")) {
      PredictProR:::gp_shared_preprocess_cache_set(
        cache,
        "model_input_objects",
        list(marker = "prepared_once")
      )
    } else {
      PredictProR:::gp_shared_preprocess_cache_get(cache, "model_input_objects")
    }
    model <- if (isTRUE(call$cross_validation)) call$GS_model_cv else call$GS_model
    if (!isTRUE(call$cross_validation)) {
      return(list(model_results = list(selected_model = model)))
    }
    list(
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = model,
          metric = "accuracy",
          value = if (identical(model, "Xgboost")) 0.9 else 0.5,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  args <- multitrait_auto_test_args(cv_evaluation_only = FALSE)
  args$ranking_tie_breakers <- NULL
  out <- PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_length(calls, 3L)
  expect_true(all(vapply(
    calls[-1L],
    function(call) identical(
      call$.predictpror_preprocess_cache,
      calls[[1L]]$.predictpror_preprocess_cache
    ),
    logical(1L)
  )))
  expect_identical(
    calls[[1L]]$.predictpror_shared_models,
    c("RandomForest", "Xgboost")
  )
  expect_true(out$run_metadata$shared_preprocessing)
  expect_identical(out$run_metadata$preprocessing_stages_built, 1L)
  expect_identical(out$run_metadata$preprocessing_stage_reuse_hits, 2L)
})

test_that("shared preparation uses the union and strictest route requirements", {
  args <- multitrait_auto_test_args(cv_evaluation_only = TRUE)
  args$GS_model_cv <- c("RandomForest", "GBLUP", "RKHS")
  args$engine <- "other"
  args$gmatrix_method <- NULL
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    list(
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = call$GS_model_cv,
          metric = "accuracy",
          value = 0.5,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_length(calls, 3L)
  expect_true(all(vapply(
    calls,
    function(call) identical(
      call$.predictpror_shared_models,
      c("RandomForest", "GBLUP", "RKHS")
    ),
    logical(1L)
  )))
  expect_true(all(vapply(
    calls,
    function(call) isTRUE(call$.predictpror_shared_requires_asreml),
    logical(1L)
  )))
  expect_true(all(vapply(
    calls,
    function(call) identical(call$gmatrix_method, "Yang"),
    logical(1L)
  )))

  requested <- PredictProR:::gp_shared_preprocess_requested_models(c(
    calls[[1L]],
    list(GS_model_cv = "RandomForest")
  ))
  expect_true(all(c("RandomForest", "GBLUP", "RKHS") %in% requested))
  expect_true(PredictProR:::gp_asreml_order_required(calls[[1L]]))
})

test_that("shared preprocessing stages are immutable", {
  cache <- PredictProR:::gp_shared_preprocess_cache_new("RandomForest")
  PredictProR:::gp_shared_preprocess_cache_set(cache, "phenotype_to_model", list(x = 1L))
  expect_identical(
    PredictProR:::gp_shared_preprocess_cache_get(cache, "phenotype_to_model"),
    list(x = 1L)
  )
  expect_error(
    PredictProR:::gp_shared_preprocess_cache_set(cache, "phenotype_to_model", list(x = 2L)),
    "replace immutable"
  )
})

test_that("automatic multi-trait CV refits a selected kernel model through its original route", {
  args <- multitrait_auto_test_args(cv_evaluation_only = FALSE)
  args$GS_model_cv <- c("FA-GBLUP", "RKHS")
  args$ranking_tie_breakers <- NULL
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    model <- if (isTRUE(call$cross_validation)) call$GS_model_cv else call$GS_model
    if (!isTRUE(call$cross_validation)) {
      return(list(model_results = list(selected_model = model)))
    }
    list(
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = model,
          metric = "accuracy",
          value = if (identical(model, "RKHS")) 0.9 else 0.4,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  out <- PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_identical(out$selected_model, "RKHS")
  expect_length(calls, 3L)
  expect_true(calls[[1L]]$multi_trait_gp)
  expect_true(calls[[2L]]$multi_trait_bayes)
  expect_false(calls[[3L]]$cross_validation)
  expect_true(calls[[3L]]$multi_trait_bayes)
  expect_identical(calls[[3L]]$GS_model, "RKHS")
})

test_that("automatic orchestration owns route-specific engine selection", {
  args <- multitrait_auto_test_args(cv_evaluation_only = TRUE)
  args$GS_model_cv <- "GBLUP"
  args$engine <- "other"
  calls <- list()
  fake_execute <- function(...) {
    call <- list(...)
    calls[[length(calls) + 1L]] <<- call
    list(
      cv_results_raw = list(),
      cv_results_processed = list(
        multitrait_metric_summary = data.frame(
          trait = c("Trait1", "Trait2"),
          model = "GBLUP",
          metric = "accuracy",
          value = c(0.4, 0.5),
          stringsAsFactors = FALSE
        )
      )
    )
  }

  PredictProR:::gp_model_execute_multitrait_orchestrate(args, fake_execute)

  expect_identical(calls[[1L]]$engine, "asreml")
  expect_true(calls[[1L]]$multi_trait_asreml)
})

test_that("automatic orchestration rejects ambiguous manual route flags", {
  args <- multitrait_auto_test_args()
  args$multi_trait_ml <- TRUE
  expect_error(
    PredictProR:::gp_model_execute_multitrait_orchestrate(
      args,
      execute_fn = function(...) NULL
    ),
    "leave multi_trait_asreml"
  )
})

test_that("automatic orchestration resolves response_family before fitting", {
  args <- multitrait_auto_test_args()
  args$response_family <- "auto"
  args$pheno_data$Trait1 <- rep(c(0L, 1L), 3L)
  args$pheno_data$Trait2 <- rep(c(1L, 0L), 3L)
  calls <- 0L

  expect_error(
    PredictProR:::gp_model_execute_multitrait_orchestrate(
      args,
      execute_fn = function(...) {
        calls <<- calls + 1L
        NULL
      }
    ),
    "gaussian traits only"
  )
  expect_identical(calls, 0L)
})

test_that("multi-trait Bayesian and ASReml standards accept raw keyed genotypes", {
  pheno <- data.frame(
    GID = paste0("G", 1:4),
    Trait1 = 1:4,
    Trait2 = 5:8,
    stringsAsFactors = FALSE
  )
  geno <- matrix(
    rep(c(0, 1, 2), length.out = 12),
    nrow = 4,
    dimnames = list(pheno$GID, paste0("M", 1:3))
  )

  expect_silent(validate_multi_trait_input_standard(
    pheno_data = pheno,
    geno_data = geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    response_family = "gaussian",
    multi_trait_bayes = TRUE,
    GS_model = "RKHS"
  ))
  expect_silent(validate_multi_trait_input_standard(
    pheno_data = pheno,
    geno_data = geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    response_family = "gaussian",
    multi_trait_asreml = TRUE,
    GS_model = "GBLUP"
  ))

  expect_silent(validate_multi_trait_input_standard(
    pheno_data = pheno,
    geno_data = geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    response_family = "gaussian",
    multi_trait_bayes = TRUE,
    GS_model_cv = "GBLUP_BRR",
    cross_validation = TRUE
  ))
  expect_silent(validate_multi_trait_input_standard(
    pheno_data = pheno,
    geno_data = geno,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    response_family = "gaussian",
    multi_trait_gp = TRUE,
    GS_model_cv = "Gaussian-Process-GBLUP",
    cross_validation = TRUE
  ))
})

test_that("joint MT-MET GP and Bayesian CV require MET scenarios", {
  pheno <- data.frame(
    GID = rep(paste0("G", 1:4), each = 2L),
    Env = rep(c("E1", "E2"), 4L),
    Trait1 = seq_len(8),
    Trait2 = seq_len(8) + 0.5,
    stringsAsFactors = FALSE
  )
  geno <- matrix(
    rep(c(0, 1, 2), length.out = 12),
    nrow = 4,
    dimnames = list(paste0("G", 1:4), paste0("M", 1:3))
  )

  expect_silent(
    validate_multi_trait_input_standard(
      pheno_data = pheno,
      geno_data = geno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      heter_groups = "Env",
      response_family = "gaussian",
      multi_trait_gp = TRUE,
      GS_model_cv = "FA-GBLUP",
      cross_validation = TRUE,
      cross_validation_meth = "CV1"
    )
  )
  expect_silent(
    validate_multi_trait_input_standard(
      pheno_data = pheno,
      geno_data = geno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      heter_groups = "Env",
      response_family = "gaussian",
      multi_trait_bayes = TRUE,
      GS_model_cv = "RKHS",
      cross_validation = TRUE,
      cross_validation_meth = "Repeated_CV2"
    )
  )
  expect_error(
    validate_multi_trait_input_standard(
      pheno_data = pheno,
      geno_data = geno,
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      heter_groups = "Env",
      response_family = "gaussian",
      multi_trait_gp = TRUE,
      GS_model_cv = "FA-GBLUP",
      cross_validation = TRUE,
      cross_validation_meth = "K-Folds"
    ),
    "MT-MET.*CV0"
  )
})
