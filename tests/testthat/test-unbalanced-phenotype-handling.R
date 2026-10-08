test_that("phenotype_precheck allows unbalanced multi-environment phenotypes", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g1", "g2"),
    Env = c("E1", "E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )

  expect_warning(
    res <- PredictProR::phenotype_precheck(
      pheno_data = pheno,
      gen_name = "GID",
      response = "Yield",
      heter_groups = "Env"
    ),
    "not identical across environments"
  )

  expect_identical(attr(res, "cleared"), "pass")
  expect_true(is.list(attr(res, "connectivity_summary")))
})

test_that("phenotype_to_model keeps trait-specific inferred testing when NA patterns differ", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    Height = c(NA, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    res <- PredictProR::phenotype_to_model(
      pheno_data = pheno,
      response = c("Yield", "Height"),
      gen_name = "GID"
    ),
    "trait-specific inferred testing sets"
  )

  expect_true(is.null(res[["test_set"]]))
  expect_identical(res[["test_set_by_trait"]][["Yield"]], "g2")
  expect_identical(res[["test_set_by_trait"]][["Height"]], "g1")
})

test_that("phenotype_to_model keeps inferred missing responses trait-specific even when patterns match", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    Height = c(10, NA, 30),
    stringsAsFactors = FALSE
  )

  res <- PredictProR::phenotype_to_model(
    pheno_data = pheno,
    response = c("Yield", "Height"),
    gen_name = "GID"
  )

  expect_true(is.null(res[["test_set"]]))
  expect_identical(res[["test_set_by_trait"]][["Yield"]], "g2")
  expect_identical(res[["test_set_by_trait"]][["Height"]], "g2")
  expect_identical(res[["test_set_by_trait_source"]], "inferred_missing_response")
})

test_that("trait-specific ML data preparation subsets train and test by response", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    Height = c(NA, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    pheno_clean <- PredictProR::phenotype_to_model(
      pheno_data = pheno,
      response = c("Yield", "Height"),
      gen_name = "GID"
    ),
    "trait-specific inferred testing sets"
  )

  geno <- matrix(1:9, nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("m1", "m2", "m3")))
  ml_dat <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = c("Yield", "Height"),
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )

  expect_true(is.null(ml_dat[["test_set"]]))
  expect_identical(sort(names(ml_dat[["test_set_by_trait"]])), c("Height", "Yield"))

  yield_dat <- PredictProR:::gp_prepare_ml_data_for_trait(
    ml_dat_res = ml_dat,
    pheno_clean = pheno_clean,
    response = "Yield",
    gen_name = "GID"
  )

  height_dat <- PredictProR:::gp_prepare_ml_data_for_trait(
    ml_dat_res = ml_dat,
    pheno_clean = pheno_clean,
    response = "Height",
    gen_name = "GID"
  )

  expect_identical(yield_dat[["test_set"]], "g2")
  expect_identical(height_dat[["test_set"]], "g1")
  expect_false("g2" %in% rownames(yield_dat[["merged_data"]][["merge_data"]]))
  expect_false("g1" %in% rownames(height_dat[["merged_data"]][["merge_data"]]))
})

test_that("explicit phenotype test set overrides missing-response inference", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, NA, 3, 4),
    Height = c(1, 2, NA, 4),
    stringsAsFactors = FALSE
  )

  res <- PredictProR::phenotype_to_model(
    pheno_data = pheno,
    response = c("Yield", "Height"),
    gen_name = "GID",
    test_set = "g4"
  )

  expect_identical(as.character(res[["test_set"]]), "g4")
  expect_true(is.null(res[["test_set_by_trait"]]))
  expect_identical(res[["test_set_source"]], "explicit")
})

test_that("explicit phenotype test set excludes missing non-test rows from trait training", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, NA, 3, 4),
    Height = c(1, 2, NA, 4),
    stringsAsFactors = FALSE
  )

  pheno_clean <- PredictProR::phenotype_to_model(
    pheno_data = pheno,
    response = c("Yield", "Height"),
    gen_name = "GID",
    test_set = "g4"
  )
  geno <- matrix(1:16, nrow = 4, dimnames = list(pheno$GID, paste0("m", 1:4)))

  ml_dat <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = c("Yield", "Height"),
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )

  yield_dat <- PredictProR:::gp_prepare_ml_data_for_trait(
    ml_dat_res = ml_dat,
    pheno_clean = pheno_clean,
    response = "Yield",
    gen_name = "GID"
  )

  expect_identical(yield_dat[["test_set"]], "g4")
  expect_false("g2" %in% yield_dat[["pheno_clean_data"]][["GID"]])
  expect_false("g2" %in% rownames(yield_dat[["merged_data"]][["merge_data"]]))
  expect_true("g4" %in% rownames(yield_dat[["merged_data_test"]]))
})

test_that("explicit train set defines complementary test set and keeps trait training observed", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, NA, 3, 4),
    Height = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  pheno_clean <- PredictProR::phenotype_to_model(
    pheno_data = pheno,
    response = c("Yield", "Height"),
    gen_name = "GID",
    train_set = c("g1", "g3")
  )
  geno <- matrix(1:16, nrow = 4, dimnames = list(pheno$GID, paste0("m", 1:4)))

  ml_dat <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = c("Yield", "Height"),
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )
  yield_dat <- PredictProR:::gp_prepare_ml_data_for_trait(
    ml_dat_res = ml_dat,
    pheno_clean = pheno_clean,
    response = "Yield",
    gen_name = "GID"
  )

  expect_identical(as.character(pheno_clean[["test_set"]]), c("g2", "g4"))
  expect_identical(pheno_clean[["test_set_source"]], "train_set")
  expect_identical(yield_dat[["test_set"]], c("g2", "g4"))
  expect_identical(rownames(yield_dat[["merged_data"]][["merge_data"]]), c("g1", "g3"))
  expect_identical(rownames(yield_dat[["merged_data_test"]]), c("g2", "g4"))
})

test_that("geno-only test set is added to trait-specific missing response test sets", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    Height = c(NA, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    pheno_clean <- PredictProR::phenotype_to_model(
      pheno_data = pheno,
      response = c("Yield", "Height"),
      gen_name = "GID"
    ),
    "trait-specific inferred testing sets"
  )

  geno <- matrix(
    c(
      0, 1, 2, 0,
      1, 2, 0, 1,
      2, 0, 1, 2,
      0, 2, 1, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3", "g4"), paste0("m", 1:4))
  )

  model_inputs <- suppressMessages(PredictProR:::gp_prepare_model_input_objects(list(
    geno_data = geno,
    train_geno_data = NULL,
    test_geno_data = NULL,
    test_set = NULL,
    pheno_clean = pheno_clean,
    train_set = NULL,
    gen_name = "GID",
    kernel_method = NULL,
    gmatrix_method = NULL,
    scale = FALSE,
    map_data = NULL,
    maf_threshold = 0.01,
    het_threshold = 0.95,
    ind_call_rate_threshold = 0.95,
    snp_call_rate_threshold = 0.95,
    impute = FALSE,
    imputation_method = "knn",
    impute_knn_k = 5,
    ld_prunning_qc = FALSE,
    qc_filtering = FALSE,
    message = FALSE,
    heter_groups = NULL,
    omic1_data = NULL,
    train_omic1_data = NULL,
    test_omic1_data = NULL,
    omic2_data = NULL,
    train_omic2_data = NULL,
    test_omic2_data = NULL,
    omic3_data = NULL,
    train_omic3_data = NULL,
    test_omic3_data = NULL,
    impute_omic = FALSE,
    na_threshold = 0.9,
    gmatrix = NULL,
    gkernel = NULL,
    kernel_list = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    pedigree_matrix = NULL,
    bending = FALSE,
    bend_value = NULL,
    blending = FALSE,
    blending_value = NULL,
    high_diag_cut_off = 1.2,
    low_diag_cut_off = 0.8,
    duplicate_cut_off = 0.99,
    rcn_cutoff = 1e-12,
    optimize_diagonal = FALSE,
    optimize_duplicate = FALSE,
    kernel_check_level = "auto",
    kernel_large_n_threshold = 5000L,
    duplicate_scan = "auto",
    duplicate_sample_size = 2000L,
    duplicate_block_size = 1024L,
    duplicate_max_pairs = 10000L,
    kernel_fix_method = "auto",
    kernel_repair_priority = "speed",
    kernel_rcn_check = "auto",
    kernel_nearpd_size_limit = 2500L,
    kernel_cpp_repair_size_limit = 3000L,
    kernel_cpp_keep_diag = TRUE,
    kernel_pd_check = "auto",
    kernel_pd_sample_size = 500L,
    kernel_sanitize = "auto",
    kernel_sanitize_value = NULL,
    GS_model = NULL,
    GS_model_cv = NULL,
    AI_valid_models = character(),
    met_ml_dl = FALSE,
    gp_valid_models = character()
  )))

  yield_test <- PredictProR:::gp_trait_test_set(model_inputs[["pheno_clean"]], "Yield", "GID")
  height_test <- PredictProR:::gp_trait_test_set(model_inputs[["pheno_clean"]], "Height", "GID")

  expect_true(setequal(yield_test, c("g2", "g4")))
  expect_true(setequal(height_test, c("g1", "g4")))
})

test_that("kernel-only test set is retained in prepared phenotype metadata", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, 2),
    stringsAsFactors = FALSE
  )
  pheno_clean <- PredictProR::phenotype_to_model(
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID"
  )
  K <- diag(3)
  rownames(K) <- colnames(K) <- c("g1", "g2", "g3")

  model_inputs <- suppressMessages(PredictProR:::gp_prepare_model_input_objects(list(
    geno_data = NULL,
    train_geno_data = NULL,
    test_geno_data = NULL,
    test_set = NULL,
    pheno_clean = pheno_clean,
    train_set = NULL,
    gen_name = "GID",
    kernel_method = NULL,
    gmatrix_method = NULL,
    scale = FALSE,
    map_data = NULL,
    maf_threshold = 0.01,
    het_threshold = 0.95,
    ind_call_rate_threshold = 0.95,
    snp_call_rate_threshold = 0.95,
    impute = FALSE,
    imputation_method = "knn",
    impute_knn_k = 5,
    ld_prunning_qc = FALSE,
    qc_filtering = FALSE,
    message = FALSE,
    heter_groups = NULL,
    omic1_data = NULL,
    train_omic1_data = NULL,
    test_omic1_data = NULL,
    omic2_data = NULL,
    train_omic2_data = NULL,
    test_omic2_data = NULL,
    omic3_data = NULL,
    train_omic3_data = NULL,
    test_omic3_data = NULL,
    impute_omic = FALSE,
    na_threshold = 0.9,
    gmatrix = K,
    gkernel = NULL,
    kernel_list = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    pedigree_matrix = NULL,
    bending = FALSE,
    bend_value = NULL,
    blending = FALSE,
    blending_value = NULL,
    high_diag_cut_off = 1.2,
    low_diag_cut_off = 0.8,
    duplicate_cut_off = 0.99,
    rcn_cutoff = 1e-12,
    optimize_diagonal = FALSE,
    optimize_duplicate = FALSE,
    kernel_check_level = "auto",
    kernel_large_n_threshold = 5000L,
    duplicate_scan = "auto",
    duplicate_sample_size = 2000L,
    duplicate_block_size = 1024L,
    duplicate_max_pairs = 10000L,
    kernel_fix_method = "auto",
    kernel_repair_priority = "speed",
    kernel_rcn_check = "auto",
    kernel_nearpd_size_limit = 2500L,
    kernel_cpp_repair_size_limit = 3000L,
    kernel_cpp_keep_diag = TRUE,
    kernel_pd_check = "auto",
    kernel_pd_sample_size = 500L,
    kernel_sanitize = "auto",
    kernel_sanitize_value = NULL,
    GS_model = "KRR",
    GS_model_cv = NULL,
    AI_valid_models = character(),
    met_ml_dl = FALSE,
    gp_valid_models = PredictProR:::gp_lowrank_supported_models()
  )))

  expect_identical(model_inputs[["test_set"]], "g3")
  expect_identical(model_inputs[["pheno_clean"]][["test_set"]], "g3")
})

test_that("general input validation evaluates repeated IDs on observed rows per response", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1, NA, 3, NA),
    Height = c(10, 11, 12, 13),
    stringsAsFactors = FALSE
  )
  geno <- matrix(
    1:8,
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("m", 1:4))
  )

  expect_silent(
    PredictProR::validate_general_input_standard(
      pheno_data = pheno,
      geno_data = geno,
      response = "Yield",
      gen_name = "GID",
      GS_model = "Ridge_Regression",
      AI_valid_models = "Ridge_Regression"
    )
  )

  expect_error(
    PredictProR::validate_general_input_standard(
      pheno_data = pheno,
      geno_data = geno,
      response = "Height",
      gen_name = "GID",
      GS_model = "Ridge_Regression",
      AI_valid_models = "Ridge_Regression"
    ),
    regexp = "Repeated genotype rows"
  )
})

test_that("heter_groups are cleared when selected response is single-environment", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1, NA, 3, NA),
    Height = c(10, 11, 12, 13),
    stringsAsFactors = FALSE
  )

  single_trait <- PredictProR:::gp_normalize_single_environment_heter_controls(
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "us",
    response = "Yield"
  )
  expect_null(single_trait$heter_groups)
  expect_false(single_trait$heter_resid)
  expect_null(single_trait$var_cov_str)
  expect_true(single_trait$single_environment)

  multi_trait_pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    Height = c(10, 11, 12),
    stringsAsFactors = FALSE
  )
  multi_trait_structure <- PredictProR:::gp_normalize_single_environment_heter_controls(
    pheno_data = multi_trait_pheno,
    gen_name = "GID",
    var_cov_str = "diag",
    response = c("Yield", "Height"),
    preserve_var_cov_str = TRUE
  )
  expect_identical(multi_trait_structure$var_cov_str, "diag")
  expect_true(multi_trait_structure$single_environment)

  met_trait <- PredictProR:::gp_normalize_single_environment_heter_controls(
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "us",
    response = "Height"
  )
  expect_identical(met_trait$heter_groups, "Env")
  expect_true(met_trait$heter_resid)
  expect_identical(met_trait$var_cov_str, "us")
  expect_false(met_trait$single_environment)
})

test_that("connectivity summary tables are materialized for output export", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g1", "g2"),
    Env = c("E1", "E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )

  expect_warning(
    pheno_checked <- PredictProR::phenotype_precheck(
      pheno_data = pheno,
      gen_name = "GID",
      response = "Yield",
      heter_groups = "Env"
    ),
    "not identical across environments"
  )

  connectivity_tables <- PredictProR:::gp_connectivity_summary_tables(pheno_checked)

  expect_true(is.data.frame(connectivity_tables[["Connectivity_summary"]]))
  expect_true(is.data.frame(connectivity_tables[["Connectivity_genotype_count_by_env"]]))
  expect_true(all(c("metric", "value") %in% names(connectivity_tables[["Connectivity_summary"]])))
  expect_true(all(c("env", "genotype_count") %in% names(connectivity_tables[["Connectivity_genotype_count_by_env"]])))
})
