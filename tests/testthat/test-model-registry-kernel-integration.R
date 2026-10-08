test_that("early runtime helpers expose model registry functions used by model_execute", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper <- file.path(root, "R", "aaa_runtime_helpers.R")
  if (!file.exists(helper)) {
    testthat::skip("raw R sources are not available in installed-package test context")
  }

  env <- new.env(parent = globalenv())
  sys.source(helper, envir = env)

  registry_functions <- c(
    "gp_lowrank_supported_models",
    "gp_single_response_single_environment_gp_supported_models",
    "gp_single_trait_multi_environment_gp_supported_models",
    "gp_validate_single_trait_gp_model_scope",
    "gp_hybrid_asreml_supported_models",
    "gp_hybrid_bayes_supported_models",
    "gp_hybrid_gp_supported_models",
    "gp_hybrid_ml_supported_models",
    "gp_hybrid_dl_supported_models",
    "gp_multitrait_asreml_supported_models",
    "gp_multitrait_ml_supported_models",
    "gp_multitrait_dl_supported_models",
    "gp_met_supported_models",
    "gp_is_met_ml_dl_model",
    "gp_multi_environment_kernel_models",
    "gp_multi_environment_kernel_display_names",
    "gp_multi_environment_supported_models",
    "gp_multi_environment_supported_display_names",
    "gp_lowrank_friendly_models",
    "gp_canonicalize_supported_model_names",
    "gp_display_supported_model_names",
    "gp_all_model_display_names",
    "gp_model_friendly_lookup",
    "gp_public_model_label",
    "gp_relabel_public_model_parameters"
  )

  expect_true(all(vapply(registry_functions, exists, logical(1), envir = env, inherits = FALSE)))
  expect_true(all(c("KRR", "GP", "GP_FA", "LowRankGP") %in% env$gp_lowrank_supported_models()))
  expect_identical(
    env$gp_canonicalize_supported_model_names(
      c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")
    ),
    c("KRR", "GP", "GP_FA", "LowRankGP")
  )
  expect_identical(
    env$gp_canonicalize_supported_model_names(c("Kernel-GBLUP", "Xgboost")),
    c("KRR", "Xgboost")
  )
  expect_identical(
    env$gp_display_supported_model_names(c("KRR", "GP", "GP_FA", "LowRankGP")),
    c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")
  )
  displayed <- env$gp_all_model_display_names(
    AI_valid_models = c("Xgboost", env$gp_deep_learning_supported_models()),
    bayes_valid_models = c("BRR", "BayesB"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    gp_valid_models = env$gp_lowrank_supported_models(),
    asreml_model = "GBLUP"
  )
  expect_true(all(c(
    "Conv1DNet", "TabNet", "Kernel-GBLUP", "FA-GBLUP",
    "Scalable-GBLUP", "Gaussian-Process-GBLUP"
  ) %in% displayed))
  expect_false(any(c("cnn", "tabnet", "KRR", "GP_FA", "LowRankGP") %in% displayed))

  friendly_lookup <- env$gp_model_friendly_lookup()
  expect_identical(unname(friendly_lookup[c("cnn", "tabnet", "KRR", "GP_FA")]),
                   c("Conv1DNet", "TabNet", "Kernel-GBLUP", "FA-GBLUP"))
  expect_identical(env$gp_public_model_label(c("KRR", "cnn")),
                   c("Kernel-GBLUP", "Conv1DNet"))
  params <- data.frame(
    stat = c("gp_model", "gp_backend"),
    summary = c("KRR", "auto"),
    stringsAsFactors = FALSE
  )
  params <- env$gp_relabel_public_model_parameters(params, "KRR")
  expect_identical(params$summary[params$stat == "gp_model"], "Kernel-GBLUP")
  expect_identical(params$summary[params$stat == "gp_model_canonical"], "KRR")
  expect_true(all(c("GBLUP_BRR", "RKHS") %in% env$gp_hybrid_bayes_supported_models()))
  expect_true(env$gp_is_met_ml_dl_model("RandomForest"))
  expect_false(env$gp_is_met_ml_dl_model("KRR"))
  expect_identical(
    env$gp_single_trait_multi_environment_gp_supported_models(),
    c("KRR", "GP", "GP_FA")
  )
  met_kernel_models <- env$gp_multi_environment_kernel_models(
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    asreml_model = "GBLUP",
    gp_valid_models = env$gp_lowrank_supported_models()
  )
  expect_true(all(c("GBLUP_BRR", "RKHS", "GBLUP", "KRR", "GP", "GP_FA") %in% met_kernel_models))
  expect_false("LowRankGP" %in% met_kernel_models)
  met_display <- env$gp_multi_environment_supported_display_names(
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    asreml_model = "GBLUP",
    gp_valid_models = env$gp_lowrank_supported_models()
  )
  expect_true(all(c(
    "GBLUP", "GBLUP_BRR", "RKHS", "Kernel-GBLUP",
    "Gaussian-Process-GBLUP", "FA-GBLUP",
    "RandomForest", "TabNet"
  ) %in% met_display))
  expect_false("Scalable-GBLUP" %in% met_display)
})

test_that("kernel source names survive phenotype matching", {
  ids <- paste0("g", 1:4)
  K1 <- diag(4)
  K2 <- diag(4) + 0.1
  dimnames(K1) <- dimnames(K2) <- list(ids, ids)
  attr(K1, "predictpror_source_kernel") <- "Gaussian"
  attr(K2, "predictpror_source_kernel") <- "Matern32"

  matched <- gp_match_kernel_bank(
    kernels = list(first_internal = K1, second_internal = K2),
    object_pheno = data.frame(GID = ids, Yield = seq_along(ids)),
    gen_name = "GID",
    message = FALSE
  )
  bank <- gp_collect_kernel_inputs(kernel_list = matched$geno_pheno_match_data)

  expect_identical(names(bank), c("Gaussian", "Matern32"))
  expect_identical(
    vapply(
      matched$geno_pheno_match_data,
      gp_kernel_source_name,
      character(1),
      fallback = "missing"
    ),
    c(first_internal = "Gaussian", second_internal = "Matern32")
  )
})

test_that("multi-kernel parameter rows have a stable public schema", {
  rows <- gp_multi_kernel_parameter_rows(
    kernel_names = c("gmatrix", "Gaussian", "Matern32"),
    strategy = "separate_terms"
  )
  values <- stats::setNames(rows$summary, rows$stat)

  expect_identical(unname(values[["multi_kernel_count"]]), "3")
  expect_identical(
    unname(values[["multi_kernel_names"]]),
    "gmatrix;Gaussian;Matern32"
  )
  expect_identical(unname(values[["multi_kernel_strategy"]]), "separate_terms")
})

test_that("single primary kernels receive the canonical public source name", {
  ids <- paste0("g", 1:3)
  K <- diag(3L)
  dimnames(K) <- list(ids, ids)

  inputs <- gp_add_kernel_inputs(list(), "gmatrix_model_ready", K)

  expect_identical(
    gp_kernel_source_name(inputs$gmatrix_model_ready, "missing"),
    "gmatrix"
  )
  expect_identical(
    names(gp_collect_kernel_inputs(
      gmatrix = inputs$gmatrix_model_ready,
      kernel_list = inputs
    )),
    "gmatrix"
  )
})

test_that("model_execute exposes kernel_list as a user-facing kernel-bank input", {
  expect_true("kernel_list" %in% names(formals(model_execute)))
})

test_that("multiple GRM methods are preserved as a GP kernel bank", {
  geno <- matrix(
    c(
      0, 1, 2, 0,
      1, 0, 1, 2,
      2, 1, 0, 1,
      0, 2, 1, 0,
      1, 1, 2, 2
    ),
    nrow = 5,
    byrow = TRUE
  )
  rownames(geno) <- paste0("g", seq_len(nrow(geno)))
  colnames(geno) <- paste0("m", seq_len(ncol(geno)))

  grms <- grm_calculation(
    geno_clean = geno,
    method = c("VanRaden", "Yang")
  )
  kernel_inputs <- gp_add_kernel_inputs(
    kernel_inputs = list(),
    kernel_key = "gmatrix_model_ready",
    kernels = grms
  )
  bank <- gp_collect_kernel_inputs(
    gmatrix = kernel_inputs[["gmatrix_model_ready"]],
    kernel_list = kernel_inputs
  )

  expect_identical(names(grms), c("VanRaden", "Yang"))
  expect_true(all(c("gmatrix_model_ready", "gmatrix_Yang_model_ready") %in% names(kernel_inputs)))
  expect_equal(length(bank), 2L)
  expect_true(all(vapply(bank, gp_is_kernel_matrix, logical(1))))
  expect_identical(rownames(bank[[1L]]), rownames(geno))
  expect_identical(rownames(bank[[2L]]), rownames(geno))
})

test_that("multiple omics kernel methods are preserved as a GP kernel bank", {
  omic <- matrix(
    c(
      0.1, 0.4, 0.7,
      0.2, 0.5, 0.8,
      1.1, 1.4, 1.7,
      1.2, 1.5, 1.8,
      2.1, 2.4, 2.7
    ),
    nrow = 5,
    byrow = TRUE
  )
  rownames(omic) <- paste0("g", seq_len(nrow(omic)))
  colnames(omic) <- paste0("o", seq_len(ncol(omic)))

  kernels <- kernel_calculation(
    M_matrix_clean = omic,
    method = c("Linear_kernel", "Gaussian_kernel"),
    scaling = TRUE,
    message = FALSE,
    backend = "r"
  )
  kernel_inputs <- gp_add_kernel_inputs(
    kernel_inputs = list(),
    kernel_key = "omic1_kernel_model_ready",
    kernels = kernels
  )
  bank <- gp_collect_kernel_inputs(
    omic1_kernel = kernel_inputs[["omic1_kernel_model_ready"]],
    kernel_list = kernel_inputs
  )

  expect_identical(names(kernels), c("Linear_kernel", "Gaussian_kernel"))
  expect_true(all(c("omic1_kernel_model_ready", "omic1_Gaussian_kernel_model_ready") %in% names(kernel_inputs)))
  expect_equal(length(bank), 2L)
  expect_true(all(vapply(bank, gp_is_kernel_matrix, logical(1))))
  expect_identical(rownames(bank[[1L]]), rownames(omic))
  expect_identical(rownames(bank[[2L]]), rownames(omic))
})

test_that("user supplied kernel_list can be the only GP relationship source", {
  ids <- paste0("g", seq_len(4))
  K1 <- matrix(
    c(
      1.0, 0.2, 0.1, 0.0,
      0.2, 1.0, 0.3, 0.1,
      0.1, 0.3, 1.0, 0.4,
      0.0, 0.1, 0.4, 1.0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )
  K2 <- matrix(
    c(
      1.0, 0.5, 0.2, 0.1,
      0.5, 1.0, 0.4, 0.2,
      0.2, 0.4, 1.0, 0.6,
      0.1, 0.2, 0.6, 1.0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  bank <- gp_collect_kernel_inputs(
    kernel_list = list(user_additive = K1, user_similarity = K2)
  )
  public_bank <- gp_public_kernel_bank_inputs(
    kernels = bank,
    pheno_ids = ids
  )

  expect_identical(names(bank), c("user_additive", "user_similarity"))
  expect_equal(length(public_bank$Ks), 2L)
  expect_identical(names(public_bank$Ks), c("user_additive", "user_similarity"))
  expect_identical(public_bank$geno_ids, ids)
  expect_equal(public_bank$K, K1)
})

test_that("general input contract accepts kernel_list for kernel models", {
  ids <- paste0("g", seq_len(4))
  pheno <- data.frame(
    GID = ids,
    Yield = c(1, 2, 3, NA),
    stringsAsFactors = FALSE
  )
  kernel <- diag(length(ids))
  rownames(kernel) <- colnames(kernel) <- ids

  expect_no_error(validate_general_input_standard(
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    GS_model = "KRR",
    AI_valid_models = c("RandomForest"),
    bayes_valid_models = c("BRR"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    gp_valid_models = gp_lowrank_supported_models(),
    asreml_model = "GBLUP",
    kernel_list = list(user_kernel = kernel)
  ))
})

test_that("kernel_calculation exposes Matern-family and related kernels", {
  mat <- matrix(
    c(
      0.1, 0.4, 0.9,
      1.0, 1.2, 1.4,
      2.1, 2.0, 1.7,
      3.0, 3.2, 3.5
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(mat) <- paste0("id", seq_len(nrow(mat)))

  methods <- c(
    "Matern_kernel",
    "Matern32_kernel",
    "Matern52_kernel",
    "Laplacian_kernel",
    "RationalQuadratic_kernel"
  )
  kernels <- kernel_calculation(
    M_matrix_clean = mat,
    method = methods,
    scaling = TRUE,
    length_scale = 1.5,
    smoothness_parameter = 1.5,
    gamma = 2,
    message = FALSE,
    backend = "r"
  )

  expect_identical(names(kernels), methods)
  for (kernel in kernels) {
    expect_true(is.matrix(kernel))
    expect_true(all(is.finite(kernel)))
    expect_equal(kernel, t(kernel), tolerance = 1e-12)
    expect_equal(unname(diag(kernel)), rep(1, nrow(kernel)), tolerance = 1e-12)
  }

  alias_kernels <- kernel_calculation(
    M_matrix_clean = mat,
    method = c("Matern", "Matern32", "Matern52", "Laplacian", "RationalQuadratic"),
    scaling = TRUE,
    length_scale = 1.5,
    smoothness_parameter = 1.5,
    gamma = 2,
    message = FALSE,
    backend = "r"
  )
  expect_identical(names(alias_kernels), methods)
})

test_that("distance-kernel default length scale adapts to predictor dimension", {
  set.seed(1204)
  mat <- matrix(stats::rnorm(48 * 120), nrow = 48, ncol = 120)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))

  kernel <- kernel_calculation(
    M_matrix_clean = mat,
    method = "Matern32_kernel",
    scaling = TRUE,
    message = FALSE,
    backend = "r"
  )

  expected_at_median_distance <- (1 + sqrt(3)) * exp(-sqrt(3))
  expect_equal(
    stats::median(kernel[lower.tri(kernel)]),
    expected_at_median_distance,
    tolerance = 1e-8
  )
  expect_gt(stats::sd(kernel[lower.tri(kernel)]), 0.01)
})

test_that("native dense kernel backend matches R for newly exposed kernels", {
  cpp_available <- getFromNamespace("gp_kernel_dense_cpp_available", "PredictProR")
  testthat::skip_if_not(cpp_available())

  mat <- matrix(
    c(
      0.2, 0.8, 1.4,
      1.1, 1.5, 1.9,
      2.0, 2.4, 2.8,
      3.1, 3.5, 3.9
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  methods <- c(
    "Matern_kernel",
    "Matern32_kernel",
    "Matern52_kernel",
    "Laplacian_kernel",
    "RationalQuadratic_kernel"
  )

  for (method in methods) {
    r_kernel <- kernel_calculation(
      mat,
      method = method,
      scaling = TRUE,
      length_scale = 1.25,
      smoothness_parameter = 1.5,
      gamma = 2,
      message = FALSE,
      backend = "r"
    )
    cpp_kernel <- kernel_calculation(
      mat,
      method = method,
      scaling = TRUE,
      length_scale = 1.25,
      smoothness_parameter = 1.5,
      gamma = 2,
      message = FALSE,
      backend = "cpp"
    )

    expect_equal(cpp_kernel, r_kernel, tolerance = 1e-10)
  }
})
