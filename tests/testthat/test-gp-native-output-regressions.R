test_that("direct GP routing passes the primary kernel exactly once", {
  ids <- paste0("g", 1:4)
  make_kernel <- function(scale) {
    out <- diag(scale, length(ids))
    rownames(out) <- colnames(out) <- ids
    out
  }
  gmatrix <- make_kernel(1)
  omic1 <- make_kernel(2)
  omic2 <- make_kernel(3)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    gp_single_trait_model = function(...) {
      captured$args <- list(...)
      stop("captured GP kernel arguments", call. = FALSE)
    },
    .package = "PredictProR"
  )

  expect_error(
    PredictProR:::gp_backend_gaussian_model(
      model_name = "KRR",
      pheno_data = data.frame(GID = ids, Yield = seq_along(ids)),
      response = "Yield",
      gen_name = "GID",
      gmatrix = gmatrix,
      omic1_kernel = omic1,
      omic2_kernel = omic2,
      kernel_list = list(
        gmatrix_model_ready = gmatrix,
        omic1_kernel_model_ready = omic1,
        omic2_kernel_model_ready = omic2
      ),
      lowrank_kernel_weights = c(
        gmatrix = 0.25,
        omic1_kernel = 0.5,
        omic2_kernel = 1.5
      )
    ),
    "captured GP kernel arguments"
  )

  expect_equal(captured$args$gmatrix, gmatrix)
  expect_identical(names(captured$args$kernel_list), c("omic1_kernel", "omic2_kernel"))
  expect_equal(captured$args$kernel_list[[1L]], omic1)
  expect_equal(captured$args$kernel_list[[2L]], omic2)
  expect_equal(captured$args$w_g, c(0.25, 0.5, 1.5))
})

test_that("model_execute fast GP route forwards a user kernel_list", {
  ids <- paste0("g", 1:4)
  pheno <- data.frame(GID = ids, Yield = c(1, 2, 3, NA_real_))
  make_kernel <- function(scale) {
    out <- diag(scale, length(ids))
    rownames(out) <- colnames(out) <- ids
    out
  }
  gmatrix <- make_kernel(1)
  extras <- list(ExtraA = make_kernel(2), ExtraB = make_kernel(3))
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    gp_backend_gaussian_model = function(kernel_list = NULL,
                                         lowrank_kernel_weights = NULL,
                                         ...) {
      captured$kernel_list <- kernel_list
      captured$kernel_weights <- lowrank_kernel_weights
      list(
        model_parameters = data.frame(stat = "gp_model", summary = "KRR"),
        predicted_values = data.frame(
          GID = ids,
          Predicted_value = seq_along(ids),
          Train_Test_Label = c("Train", "Train", "Train", "Test"),
          Observed_value = pheno$Yield,
          Standard_error = rep(0.2, length(ids)),
          PEV = rep(0.04, length(ids)),
          stringsAsFactors = FALSE
        ),
        variance_components = data.frame(
          Component = "predictive_variance",
          Components = 1,
          Standard_error = NA_real_
        ),
        diagnostic_plots = NULL
      )
    },
    summary_statistics_AI = function(...) {
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    results_handling = function(res_model_output = NULL, ...) res_model_output,
    .package = "PredictProR"
  )

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = gmatrix,
    kernel_list = extras,
    lowrank_kernel_weights = c(gmatrix = 0.2, ExtraA = 0.3, ExtraB = 0.5),
    response = "Yield",
    gen_name = "GID",
    GS_model = "Kernel-GBLUP",
    system_database = FALSE,
    message = FALSE
  )

  expect_true(is.list(out))
  expect_identical(names(captured$kernel_list), names(extras))
  expect_equal(captured$kernel_list, extras)
  expect_equal(
    captured$kernel_weights,
    c(gmatrix = 0.2, ExtraA = 0.3, ExtraB = 0.5)
  )
})

test_that("GP native parameters identify every fitted kernel", {
  ids <- paste0("g", 1:4)
  make_kernel <- function(scale) {
    out <- diag(scale, length(ids))
    rownames(out) <- colnames(out) <- ids
    out
  }
  gmatrix <- make_kernel(1)
  extra_a <- make_kernel(2)
  extra_b <- make_kernel(3)
  predictions <- data.frame(
    Name = ids,
    Prediction = seq_along(ids),
    Prediction_SE_latent = rep(0.2, length(ids)),
    Prediction_Var_latent = rep(0.04, length(ids)),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    gp_single_trait_model = function(...) {
      list(
        predictions = predictions,
        result = list(
          variance_component_method = "krr_gcv_variance_ratio",
          diagnostics = list(lam_eff = 1, lam_select = "gcv")
        ),
        info = list(backend = "mock", method = "krr_exact")
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_backend_gaussian_model(
    model_name = "KRR",
    pheno_data = data.frame(GID = ids, Yield = seq_along(ids)),
    response = "Yield",
    gen_name = "GID",
    gmatrix = gmatrix,
    kernel_list = list(ExtraA = extra_a, ExtraB = extra_b),
    lowrank_kernel_weights = c(gmatrix = 0.2, ExtraA = 0.3, ExtraB = 0.5)
  )
  params <- stats::setNames(out$model_parameters$summary, out$model_parameters$stat)

  expect_identical(unname(params[["gp_kernel_count"]]), "3")
  expect_identical(
    strsplit(unname(params[["gp_kernel_names"]]), ";", fixed = TRUE)[[1L]],
    c("gmatrix", "ExtraA", "ExtraB")
  )
  expect_identical(
    unname(params[["gp_kernel_input_weights"]]),
    "gmatrix=0.2;ExtraA=0.3;ExtraB=0.5"
  )
  expect_identical(unname(params[["gp_kernel_weight_mode"]]), "user_supplied_fixed")
})

test_that("SVM defaults to dimension-aware backend gamma", {
  expect_null(eval(formals(PredictProR::AI_svm)$sigma_value))
  expect_null(eval(formals(PredictProR::AI_svm)$gamma_value))
  expect_null(eval(formals(PredictProR::model_execute)$sigma_value))
})

test_that("single-environment GP reliability uses the fitted component scale", {
  gp_result <- list(
    # This is deliberately a different structural/kernel scale.  It must not
    # override the fitted variance component in a one-environment analysis.
    env_variance_summary = data.frame(
      level = "ENV1",
      total_genetic_variance = 0.88
    ),
    varcomp = data.frame(
      component = c("vm(GID,K)!var", "Residual!R"),
      estimate = c(0.38, 0.30),
      std.error = c(0.19, 0.14),
      stringsAsFactors = FALSE
    ),
    variance_component_method = "reml"
  )

  expect_equal(
    PredictProR:::gp_extract_gp_genetic_variance(
      gp_result,
      n = 3L,
      env_values = rep("ENV1", 3L)
    ),
    rep(0.38, 3L)
  )
})

test_that("single-trait KRR publishes its GCV variance-ratio provenance", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, NA_real_),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    Name = pheno$GID,
    Prediction = c(1.1, 1.8, 1.5),
    Prediction_SE_latent = c(0.2, 0.3, 0.4),
    Prediction_Var_latent = c(0.04, 0.09, 0.16),
    stringsAsFactors = FALSE
  )
  gp_result <- list(
    variance_component_method = "krr_gcv_variance_ratio",
    diagnostics = list(
      variance_component_method = "krr_gcv_variance_ratio",
      lam_eff = 1,
      lam_select = "gcv"
    ),
    varcomp = data.frame(
      component = c("vm(GID,KRR_kernel)!var", "Residual!R"),
      estimate = c(0.4, 0.4),
      stringsAsFactors = FALSE
    )
  )

  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = predictions,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    gp_result = gp_result
  )

  expect_equal(out$Observed_value, pheno$Yield)
  expect_equal(out$Reliability_reference_variance, rep(0.4, 3L))
  expect_equal(out$Reliability, pmax(0, 1 - out$PEV / 0.4))
  expect_identical(
    unique(out$Reliability_basis),
    "krr_gcv_variance_ratio_genetic_variance"
  )
  expect_identical(
    PredictProR:::gp_vc_gp_estimation_method(
      list(gp_result = gp_result)
    ),
    "KRR_GCV_variance_ratio"
  )
})

test_that("single-trait KRR defaults tune a nondegenerate lambda grid", {
  defaults <- formals(PredictProR::gp_single_trait_model)
  expect_identical(eval(defaults$lam_select), "gcv")
  expect_equal(
    eval(defaults$krr_lams),
    c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  )
})

test_that("native result export excludes GP bridge internals", {
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  case_dir <- tempfile("gp-native-output-")
  dir.create(case_dir)
  setwd(case_dir)

  pred <- data.frame(
    GID = c("g1", "g2"),
    Env = NA_character_,
    Trait = "Yield",
    Predicted_value = c(1.1, 1.2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.0, NA_real_),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    Prediction_error_variance = c(0.01, 0.04),
    Genetic_variance = c(0.4, 0.4),
    Reliability_reference_variance = c(0.4, 0.4),
    Reliability = c(0.975, 0.9),
    Reliability_remarks = c("Reliable", "Reliable"),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(0.4, 0.3, 4 / 7),
    Standard_error = c(0.1, 0.08, 0.12),
    stringsAsFactors = FALSE,
    row.names = c("genetic_variance", "residual_variance", "heritability")
  )

  out <- PredictProR:::results_handling(
    GS_model = "Gaussian-Process-GBLUP",
    res_model_output = list(
      predicted_values = pred,
      variance_components = vc,
      model_parameters = data.frame(stat = "model", summary = "Gaussian-Process-GBLUP"),
      raw_python_result = list(debug = data.frame(value = 1)),
      raw_gp_result = list(debug = data.frame(value = 2)),
      raw_gp_info = list(debug = data.frame(value = 3)),
      gp_result = list(debug = data.frame(value = 4)),
      gp_info = list(debug = data.frame(value = 5)),
      gp_fit = list(debug = data.frame(value = 6))
    ),
    system_database = FALSE,
    plot_filename = "gp_native_contract"
  )

  files <- list.files(case_dir, recursive = TRUE)
  expect_true(any(basename(files) == "Predicted_Value.csv"))
  expect_true(any(basename(files) == "variance_components.csv"))
  expect_false(any(grepl("raw_python|raw_gp|gp_result|gp_info|gp_fit", files, ignore.case = TRUE)))
  expect_false(any(c("raw_python_result", "raw_gp_result", "raw_gp_info", "gp_result", "gp_info", "gp_fit") %in%
                     names(out$model_results)))
  expect_identical(out$export_status, "Successful")
})

test_that("multi-trait GP covariance payload yields per-trait public variance components", {
  sigma_g <- matrix(c(0.8, 0.15, 0.15, 0.5), nrow = 2L)
  sigma_e <- matrix(c(0.2, 0.04, 0.04, 0.3), nrow = 2L)
  source <- list(
    gp_result = list(),
    gp_info = list(varcomp_mode = "reml", trait_levels = c("Yield", "Protein")),
    gp_fit = list(
      Sigma_G_response_scale = sigma_g,
      Sigma_eps_response_scale = sigma_e,
      # Standardized-scale aliases must not override response-scale matrices.
      Sigma_G = sigma_g * 100,
      Sigma_eps = sigma_e * 100
    )
  )

  vc <- PredictProR:::gp_public_variance_components(source)

  expect_identical(
    names(vc),
    c("Trait", "Component", "Components", "Standard_error")
  )
  expect_identical(unique(vc$Trait), c("Yield", "Protein"))
  expect_equal(vc$Components[vc$Component == "genetic_variance"], c(0.8, 0.5))
  expect_equal(vc$Components[vc$Component == "residual_variance"], c(0.2, 0.3))
  expect_equal(vc$Components[vc$Component == "heritability"], c(0.8, 0.625))
  expect_true(all(is.na(vc$Standard_error)))
  expect_false(any(vc$Component == "reml_variance_components_not_returned"))
})

test_that("export names are case-exact for every engine's list key", {
  f <- PredictProR:::format_model_table_export_name
  for (k in c("Predicted_value", "predicted_values", "Predicted_Values", "prediction_value")) {
    expect_identical(f(k), "Predicted_Value")
  }
  expect_identical(f("Variance_components"), "variance_components")
  expect_identical(f("variance_components"), "variance_components")
  expect_identical(f("Connectivity_summary"), "Connectivity_summary")
})

test_that("native writer omits empty placeholders and writes one clean summary table", {
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  case_dir <- tempfile("native-writer-contract-")
  dir.create(case_dir)
  setwd(case_dir)

  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 1.2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1, NA_real_),
    Standard_error = c(0.1, 0.2),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Trait = "Yield",
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(0.4, 0.3, 4 / 7),
    Standard_error = NA_real_,
    stringsAsFactors = FALSE
  )
  summary_table <- data.frame(
    stat = c("mode", "model"),
    summary = c("multi_trait", "mock"),
    stringsAsFactors = FALSE
  )

  PredictProR:::results_handling(
    GS_model = "DenseNeuralNet",
    res_model_output = list(
      predicted_values = pred,
      model_parameters = data.frame(stat = character(), summary = character()),
      Variance_components = vc,
      trained_model = NULL
    ),
    res_summary_stat = summary_table,
    system_database = FALSE,
    plot_filename = "native_writer_contract"
  )

  export_dir <- list.dirs(case_dir, recursive = FALSE, full.names = TRUE)[[1L]]
  files <- list.files(export_dir)
  expect_true("summary_statistics.csv" %in% files)
  expect_false(any(c("stat.csv", "summary.csv", "model_parameters.csv") %in% files))
  expect_false(any(grepl("trained_model\\.RData$", files)))

  vc_header <- readLines(file.path(export_dir, "variance_components.csv"), n = 1L)
  expect_identical(vc_header, '"Trait","Component","Components","Standard_error"')
  summary_native <- utils::read.csv(file.path(export_dir, "summary_statistics.csv"))
  expect_identical(summary_native, summary_table)
})
