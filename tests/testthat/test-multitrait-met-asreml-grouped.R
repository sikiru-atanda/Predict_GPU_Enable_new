test_that("grouped MT-MET ASReml structure accepts only validated FA specifications", {
  expect_equal(
    PredictProR:::gp_multitrait_asreml_mtmet_structure("fa3", 12L),
    list(token = "fa3", rank = 3L)
  )
  expect_equal(
    PredictProR:::gp_multitrait_asreml_mtmet_structure("fa(2)", 6L),
    list(token = "fa2", rank = 2L)
  )
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_structure("us", 12L),
    "validated factor-analytic"
  )
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_structure("fa4", 4L),
    "smaller than"
  )
})

test_that("grouped MT-MET long data uses an explicit collision-safe group lookup", {
  ph <- data.frame(
    GID = rep(c("g_1", "g_2"), each = 2L),
    Env = rep(c("site_A", "site_B"), 2L),
    yield_trait = 1:4,
    height_trait = 5:8,
    stringsAsFactors = FALSE
  )
  bundle <- PredictProR:::gp_multitrait_asreml_mtmet_long_data(
    ph,
    response = c("yield_trait", "height_trait"),
    gen_name = "GID",
    heter_groups = "Env"
  )

  expect_equal(nrow(bundle$data), 8L)
  expect_equal(nrow(bundle$group_lookup), 4L)
  expect_equal(bundle$group_lookup$PPGroup, sprintf("G%03d", 1:4))
  expect_setequal(bundle$group_lookup$Trait, c("yield_trait", "height_trait"))
  expect_setequal(bundle$group_lookup$Env, c("site_A", "site_B"))
  expect_false(anyNA(bundle$data$PPGroup))
})

test_that("multi-trait ASReml retries prediction workspace once on a converged fit", {
  seen <- character()
  pred <- PredictProR:::gp_multitrait_asreml_predict_pvals(
    mod = list(converge = TRUE, ifault = 0L),
    predict_fn = function(object, classify, sed, pworkspace) {
      seen <<- c(seen, as.character(pworkspace))
      if (!identical(pworkspace, "1gb")) {
        stop("Insufficient workspace available when calculating predictions")
      }
      list(pvals = data.frame(GID = "g1", predicted.value = 2))
    },
    classify = "Trait:GID", pworkspace = 1e06
  )
  expect_identical(seen, c("1e+06", "1gb"))
  expect_equal(pred$predicted.value, 2)
})

test_that("grouped MT-MET ASReml recovers a complete grid by PPGroup chunks", {
  ids <- c("g1", "g2")
  groups <- c("G001", "G002")
  bundle <- list(genotype_levels = ids,
                 group_lookup = data.frame(PPGroup = groups))
  seen <- list()
  predict_mock <- function(object, classify, sed, pworkspace,
                           levels = list()) {
    seen[[length(seen) + 1L]] <<- list(memory = pworkspace, levels = levels)
    if (!length(levels)) {
      stop("Insufficient workspace available when calculating predictions")
    }
    group <- levels$PPGroup
    list(pvals = data.frame(
      GID = ids, PPGroup = group,
      predicted.value = if (identical(group, "G001")) c(1, 2) else c(3, 4),
      std.error = rep(0.1, 2L)
    ))
  }
  out <- PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
    mod = list(converge = TRUE, ifault = 0L),
    predict_fn = predict_mock, long_bundle = bundle,
    gen_name = "GID", pworkspace = 1e06
  )
  expect_identical(out$mode, "predict_group_chunks")
  expect_equal(out$pvals$predicted.value, c(1, 2, 3, 4))
  expect_equal(length(seen), 4L)
  expect_identical(seen[[2L]]$memory, "1gb")
  expect_identical(seen[[3L]]$levels, list(PPGroup = "G001"))
})

test_that("grouped MT-MET chunk assembly preserves the public prediction table", {
  ph <- expand.grid(GID = c("g1", "g2"), Env = c("E1", "E2"),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ph$Yield <- c(1, 2, 3, NA_real_)
  ph$Height <- c(4, 5, 6, NA_real_)
  bundle <- PredictProR:::gp_multitrait_asreml_mtmet_long_data(
    ph, response = c("Yield", "Height"),
    gen_name = "GID", heter_groups = "Env"
  )
  bundle$genotype_levels <- c("g1", "g2")
  grid <- merge(data.frame(GID = bundle$genotype_levels),
                bundle$group_lookup["PPGroup"], by = NULL)
  grid$predicted.value <- match(grid$GID, bundle$genotype_levels) +
    10 * match(grid$PPGroup, bundle$group_lookup$PPGroup)
  grid$std.error <- 0.2
  direct_fn <- function(...) list(pvals = grid)
  chunk_fn <- function(object, classify, sed, pworkspace,
                       levels = list()) {
    if (!length(levels)) stop("Insufficient workspace available")
    list(pvals = grid[grid$PPGroup == levels$PPGroup, , drop = FALSE])
  }
  direct <- PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
    list(converge = TRUE), direct_fn, bundle, "GID", 1e06
  )
  chunk <- PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
    list(converge = TRUE), chunk_fn, bundle, "GID", 1e06
  )
  standardize <- function(result) {
    PredictProR:::gp_multitrait_asreml_mtmet_prediction(
      result$pvals, bundle, "GID"
    )
  }
  expect_identical(direct$mode, "predict")
  expect_identical(chunk$mode, "predict_group_chunks")
  expect_equal(standardize(chunk), standardize(direct))
})

test_that("grouped MT-MET ASReml refuses unsafe or incomplete chunk recovery", {
  bundle <- list(genotype_levels = c("g1", "g2"),
                 group_lookup = data.frame(PPGroup = c("G001", "G002")))
  fail_predict <- function(...) stop("Insufficient workspace available when calculating predictions")
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
      mod = list(converge = FALSE, ifault = 0L),
      predict_fn = fail_predict, long_bundle = bundle,
      gen_name = "GID", pworkspace = 1e06
    ),
    "did not converge"
  )
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
      mod = list(), predict_fn = fail_predict, long_bundle = bundle,
      gen_name = "GID", pworkspace = 1e06
    ),
    "did not converge"
  )
  incomplete_predict <- function(object, classify, sed, pworkspace,
                                 levels = list()) {
    if (!length(levels)) stop("Insufficient workspace available")
    ids <- if (identical(levels$PPGroup, "G001")) c("g1", "g2") else "g1"
    list(pvals = data.frame(GID = ids, PPGroup = levels$PPGroup,
                            predicted.value = rep(1, length(ids))))
  }
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_predict_or_chunk(
      mod = list(converge = TRUE, ifault = 0L),
      predict_fn = incomplete_predict, long_bundle = bundle,
      gen_name = "GID", pworkspace = 1e06
    ),
    "each genotype-by-group cell exactly once"
  )
})

test_that("automatic MT-MET inventory routes every joint covariance family for CV and prediction", {
  cv_inventory <- PredictProR:::gp_multitrait_auto_inventory(
    cross_validation = TRUE,
    is_met = TRUE
  )
  prediction_inventory <- PredictProR:::gp_multitrait_auto_inventory(
    cross_validation = FALSE,
    is_met = TRUE
  )

  expect_equal(names(cv_inventory), c("multi_trait_asreml", "multi_trait_gp", "multi_trait_bayes"))
  expect_equal(cv_inventory$multi_trait_asreml, "GBLUP")
  expect_equal(cv_inventory$multi_trait_gp, c("GP", "GP_FA", "LowRankGP"))
  expect_equal(cv_inventory$multi_trait_bayes, c("GBLUP_BRR", "RKHS"))
  expect_true("GBLUP" %in% prediction_inventory$multi_trait_asreml)
  expect_true(all(c("multi_trait_gp", "multi_trait_bayes") %in% names(prediction_inventory)))
})

test_that("specialized MT-MET indexing preserves the user environment column", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), 4L),
    loc = rep(c("site_A", "site_B"), each = 4L),
    Trait = rep(c("Yield", "Height"), each = 2L, times = 2L),
    Predicted_value = seq_len(8),
    Train_Test_Label = rep(c("Train", "Test"), 4L),
    Observed_value = c(1, NA, 3, NA, 5, NA, 7, NA),
    Standard_error = rep(0.2, 8L),
    stringsAsFactors = FALSE
  )
  vc <- unique(pred[c("Trait", "loc")])
  names(vc)[names(vc) == "loc"] <- "Env"
  vc <- vc[rep(seq_len(nrow(vc)), each = 3L), , drop = FALSE]
  vc$Component <- rep(c("genetic_variance", "residual_variance", "heritability"), nrow(vc) / 3L)
  vc$Components <- rep(c(2, 3, 0.4), nrow(vc) / 3L)
  vc$Standard_error <- 0.1
  handled <- list(model_results = list(
    predicted_values = pred,
    variance_components = vc
  ))

  indexed <- PredictProR:::gp_index_specialized_model_execute_result(
    handled_result = handled,
    model_label = "GBLUP",
    trait_name = "multi_trait",
    gen_name = "GID",
    heter_groups = "loc"
  )

  expect_true("loc" %in% names(indexed$model_results$predicted_values))
  expect_identical(indexed$model_results$predicted_values$loc, pred$loc)
  expect_setequal(unique(indexed$model_results$predicted_values$loc), c("site_A", "site_B"))
})

test_that("named FA and residual reconstruction is invariant to raw row order", {
  group_lookup <- data.frame(
    PPGroup = sprintf("G%03d", 1:4),
    Trait = rep(c("Yield", "Height"), 2L),
    Env = rep(c("E1", "E2"), each = 2L),
    stringsAsFactors = FALSE
  )
  genetic_names <- unlist(lapply(c("var", "fa1"), function(parameter) {
    paste0(
      "fa(PPGroup, 1):vm(GID, Ginv)!",
      group_lookup$PPGroup,
      "!",
      parameter
    )
  }))
  residual_names <- c(
    "PPEnv_E1!R",
    "PPEnv_E1!PPTrait!Height:!PPTrait!Yield.cor",
    "PPEnv_E1!PPTrait_Yield",
    "PPEnv_E1!PPTrait_Height",
    "PPEnv_E2!R",
    "PPEnv_E2!PPTrait!Height:!PPTrait!Yield.cor",
    "PPEnv_E2!PPTrait_Yield",
    "PPEnv_E2!PPTrait_Height"
  )
  component_names <- c(genetic_names, residual_names)
  values <- c(
    1, 2, 3, 4,
    0.5, 0.6, 0.7, 0.8,
    1, 0.2, 5, 8,
    1, -0.1, 6, 9
  )
  permutation <- c(9:16, 5:8, 1:4)
  component_names <- component_names[permutation]
  values <- values[permutation]
  raw <- data.frame(
    component = values,
    std.error = rep(0.1, length(values)),
    z.ratio = values / 0.1,
    bound = ifelse(grepl("!R$", component_names), "F", "P"),
    `%ch` = 0.1,
    row.names = component_names,
    check.names = FALSE
  )
  model <- structure(
    list(varcomp = raw, converge = TRUE, loglik = -10),
    class = "predictpror_mock_mtmet_asreml"
  )
  summary_mock <- function(object, ...) {
    list(varcomp = object$varcomp)
  }
  vpredict_mock <- function(object, xform, ...) {
    parameter_env <- list2env(
      stats::setNames(as.list(object$varcomp$component), paste0("V", seq_len(nrow(object$varcomp)))),
      parent = baseenv()
    )
    data.frame(
      Estimate = eval(xform[[3L]], envir = parameter_env),
      SE = 0.1,
      row.names = as.character(xform[[2L]])
    )
  }

  result <- PredictProR:::gp_multitrait_asreml_mtmet_variance(
    model,
    group_lookup = group_lookup,
    fa_rank = 1L,
    vpredict_fn = vpredict_mock,
    summary_fn = summary_mock
  )

  expect_equal(nrow(result$variance_components), 12L)
  expect_equal(
    result$Genetic_covariance_trait_environment["G001", "G001"],
    1 + 0.5^2
  )
  expect_equal(
    result$Genetic_covariance_trait_environment["G001", "G004"],
    0.5 * 0.8
  )
  expect_equal(unname(diag(result$Residual_covariance_by_environment$E1)), c(5, 8))
  expect_false(any(grepl("!R$", result$variance_components$Component)))
  expect_true(all(is.finite(result$variance_components$Components)))
})

test_that("grouped MT-MET ASReml sums named independent kernel covariances", {
  group_lookup <- data.frame(
    PPGroup = sprintf("G%03d", 1:4),
    Trait = rep(c("Yield", "Height"), 2L),
    Env = rep(c("E1", "E2"), each = 2L),
    stringsAsFactors = FALSE
  )
  make_genetic <- function(inverse_name, specific, loading) {
    c(
      stats::setNames(specific, paste0("fa(PPGroup, 1):vm(GID,", inverse_name, ")!", group_lookup$PPGroup, "!var")),
      stats::setNames(loading, paste0("fa(PPGroup, 1):vm(GID,", inverse_name, ")!", group_lookup$PPGroup, "!fa1"))
    )
  }
  genetic <- c(
    make_genetic("K1inv", c(1, 2, 3, 4), c(0.5, 0.6, 0.7, 0.8)),
    make_genetic("K2inv", c(0.5, 0.6, 0.7, 0.8), c(0.2, 0.3, 0.4, 0.5))
  )
  residual <- c(
    "PPEnv_E1!R" = 1,
    "PPEnv_E1!PPTrait!Height:!PPTrait!Yield.cor" = 0.2,
    "PPEnv_E1!PPTrait_Yield" = 5,
    "PPEnv_E1!PPTrait_Height" = 8,
    "PPEnv_E2!R" = 1,
    "PPEnv_E2!PPTrait!Height:!PPTrait!Yield.cor" = -0.1,
    "PPEnv_E2!PPTrait_Yield" = 6,
    "PPEnv_E2!PPTrait_Height" = 9
  )
  values <- c(genetic, residual)
  raw <- data.frame(
    component = unname(values),
    std.error = rep(0.1, length(values)),
    z.ratio = unname(values) / 0.1,
    bound = ifelse(grepl("!R$", names(values)), "F", "P"),
    `%ch` = 0.1,
    row.names = names(values),
    check.names = FALSE
  )
  model <- structure(
    list(varcomp = raw, converge = TRUE, loglik = -10),
    class = "predictpror_mock_mtmet_asreml"
  )
  summary_mock <- function(object, ...) list(varcomp = object$varcomp)
  vpredict_mock <- function(object, xform, ...) {
    parameter_env <- list2env(
      stats::setNames(as.list(object$varcomp$component), paste0("V", seq_len(nrow(object$varcomp)))),
      parent = baseenv()
    )
    data.frame(
      Estimate = eval(xform[[3L]], envir = parameter_env),
      SE = 0.1,
      row.names = as.character(xform[[2L]])
    )
  }
  specs <- data.frame(
    Kernel = c("genomic", "transcriptomic"),
    inverse_name = c("K1inv", "K2inv"),
    stringsAsFactors = FALSE
  )

  result <- PredictProR:::gp_multitrait_asreml_mtmet_variance_multikernel(
    model,
    group_lookup = group_lookup,
    fa_rank = 1L,
    kernel_specs = specs,
    vpredict_fn = vpredict_mock,
    summary_fn = summary_mock
  )

  expect_equal(
    result$Genetic_covariance_trait_environment["G001", "G001"],
    (1 + 0.5^2) + (0.5 + 0.2^2)
  )
  expect_equal(
    result$Genetic_covariance_trait_environment["G001", "G004"],
    0.5 * 0.8 + 0.2 * 0.5
  )
  expect_setequal(names(result$Genetic_covariance_by_kernel), specs$Kernel)
  expect_equal(nrow(result$kernel_variance_components), 8L)
  expect_equal(nrow(result$variance_components), 12L)
  expect_true(all(is.finite(result$variance_components$Components)))
  expect_true(all(is.finite(result$variance_components$Standard_error)))
})

test_that("public result standardization retains multi-kernel ASReml evidence", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Env = c("E1", "E1"),
    Trait = c("Yield", "Yield"),
    Predicted_value = c(1.1, 1.2),
    Train_Test_Label = c("Train", "Test"),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    Genetic_variance = c(0.7, 0.7),
    Reliability = c(0.98, 0.95),
    Reliability_remarks = c("High", "High"),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Trait = "Yield", Env = "E1", Component = "genetic_variance",
    Components = 0.7, Standard_error = 0.1,
    Estimation_method = "ASReml_REML", stringsAsFactors = FALSE
  )
  kvc <- data.frame(
    Kernel = c("genomic", "transcriptomic"), Trait = "Yield", Env = "E1",
    Component = "genetic_variance", Components = c(0.4, 0.3),
    Standard_error = c(0.08, 0.06), Estimation_method = "ASReml_REML",
    stringsAsFactors = FALSE
  )
  cov_by_kernel <- list(genomic = matrix(0.4), transcriptomic = matrix(0.3))

  public <- PredictProR:::gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      variance_components = vc,
      kernel_variance_components = kvc,
      Genetic_covariance_by_kernel = cov_by_kernel,
      genetic_covariance_by_kernel_long = kvc,
      asreml_convergence_trace = data.frame(converged = TRUE),
      model_notes = data.frame(note = "two kernels")
    ),
    gen_name = "GID",
    heter_groups = "Env",
    model = "GBLUP",
    response_family = "gaussian"
  )

  expect_identical(public$kernel_variance_components, kvc)
  expect_identical(public$Genetic_covariance_by_kernel, cov_by_kernel)
  expect_identical(public$genetic_covariance_by_kernel_long, kvc)
  expect_true(is.data.frame(public$asreml_convergence_trace))
  expect_true(is.data.frame(public$model_notes))
})

test_that("grouped MT-MET ASReml CV uses shared CV0 CV1 and CV2 assignments", {
  ph <- expand.grid(
    GID = paste0("G", seq_len(6)),
    Env = paste0("E", seq_len(3)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ph$Yield <- seq_len(nrow(ph)) / 3
  ph$Height <- seq_len(nrow(ph)) / 5 + 2
  K <- diag(6)
  rownames(K) <- colnames(K) <- paste0("G", seq_len(6))
  COP <- K + 0.1

  testthat::local_mocked_bindings(
    gp_multitrait_asreml_mtmet_gaussian_model = function(pheno_object,
                                                          response,
                                                          gen_name,
                                                          heter_groups,
                                                          ...) {
      bundle <- PredictProR:::gp_multitrait_asreml_mtmet_long_data(
        pheno_object,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      pred <- merge(
        data.frame(GID = unique(as.character(pheno_object[[gen_name]])), stringsAsFactors = FALSE),
        bundle$group_lookup,
        by = NULL
      )
      names(pred)[names(pred) == "GID"] <- gen_name
      pred$Predicted_value <- seq_len(nrow(pred)) / 10
      list(predicted_values = pred)
    },
    .package = "PredictProR"
  )

  for (method in c("CV0", "CV1", "CV2")) {
    out <- PredictProR:::gp_multitrait_asreml_mtmet_gaussian_cv(
      pheno_object = ph,
      response = c("Yield", "Height"),
      gen_name = "GID",
      heter_groups = "Env",
      gmatrix = K,
      kernel_list = list(COP = COP),
      cross_validation_meth = method,
      nfolds = if (identical(method, "CV0")) 3L else 2L,
      replication = 1L,
      eval_metrics = "mean_squared_error"
    )
    expect_identical(out$cv_results_processed$cv_scheme$method, tolower(method))
    expect_identical(out$cv_results_processed$cv_scheme$kernel_count, 2L)
    expect_setequal(out$cv_results_processed$kernel_configuration$Kernel,
                    c("gmatrix", "COP"))
    expect_equal(
      nrow(out$cv_results_processed$multitrait_oof_predictions),
      nrow(ph) * 2L
    )
  }
  expect_error(
    PredictProR:::gp_multitrait_asreml_mtmet_gaussian_cv(
      pheno_object = ph,
      response = c("Yield", "Height"),
      gen_name = "GID",
      heter_groups = "Env",
      gmatrix = K,
      cross_validation_meth = "K-Folds",
      nfolds = 2L
    ),
    "MT-MET.*CV0"
  )
})
