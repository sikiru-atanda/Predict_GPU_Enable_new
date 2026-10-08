test_that("kernel-only ML input is accepted and reaches feature fusion", {
  ids <- paste0("g", 1:5)
  pheno <- data.frame(GID = ids, y = seq_along(ids), stringsAsFactors = FALSE)
  K1 <- diag(5)
  K2 <- diag(5) + 0.1
  dimnames(K1) <- dimnames(K2) <- list(ids, ids)

  expect_silent(PredictProR::validate_general_input_standard(
    pheno_data = pheno,
    gmatrix = K1,
    kernel_list = list(omic = K2),
    response = "y",
    gen_name = "GID",
    GS_model = "RandomForest",
    AI_valid_models = "RandomForest"
  ))

  features <- PredictProR:::gp_prepare_met_kernel_features(
    kernel_list = list(genomic = K1, omic = K2),
    var_explained = 1
  )
  expect_equal(nrow(features$feature_table), length(ids))
  expect_setequal(features$feature_summary$source, c("genomic", "omic"))

  pheno_clean <- list(pheno_clean_data = pheno, test_set_by_trait = NULL)
  merged <- PredictProR:::AI_process_ml_data_if_valid(
    model_check = TRUE,
    geno_omic_model_ready_list = list(
      kernel_features_model_ready = features$feature_table
    ),
    pheno_clean = pheno_clean,
    response = "y",
    gen_name = "GID"
  )
  expect_equal(rownames(merged$merged_data$merge_data), ids)
  expect_equal(ncol(merged$merged_data$merge_data), ncol(features$feature_table))

  processed <- PredictProR:::gp_ml_preprocess_predictors(
    as.data.frame(features$feature_table),
    scaling = FALSE,
    centering = FALSE
  )
  expect_identical(rownames(processed$data), ids)
})

test_that("multi-trait ASReml kernel bank retains every named kernel", {
  ids <- paste0("g", 1:4)
  K1 <- diag(4)
  K2 <- diag(4) + 0.2
  dimnames(K1) <- dimnames(K2) <- list(ids, ids)
  attr(K1, "predictpror_source_kernel") <- "genomic"
  attr(K2, "predictpror_source_kernel") <- "metabolomic"

  bank <- PredictProR:::gp_multitrait_asreml_kernel_bank(
    gmatrix = K1,
    kernel_list = list(metabolomic = K2),
    phenotype_ids = ids
  )
  expect_equal(length(bank$kernels), 2L)
  expect_identical(bank$keep_ids, ids)
  expect_setequal(names(bank$kernels), c("genomic", "metabolomic"))

  specs <- data.frame(
    Kernel = names(bank$kernels),
    inverse_name = c("K1inv", "K2inv"),
    stringsAsFactors = FALSE
  )
  f <- PredictProR:::gp_multitrait_asreml_random_formula(
    index_factor = "Trait",
    covariance = "us",
    gen_name = "GID",
    kernel_specs = specs
  )
  formula_text <- gsub("\\s+", "", paste(deparse(f), collapse = ""))
  expect_match(formula_text, "us\\(Trait\\):vm\\(GID,K1inv\\)")
  expect_match(formula_text, "us\\(Trait\\):vm\\(GID,K2inv\\)")
})

test_that("gmatrix and gkernel are both retained when supplied together", {
  ids <- paste0("g", 1:4)
  G <- diag(4)
  K <- diag(4) + 0.15
  dimnames(G) <- dimnames(K) <- list(ids, ids)
  attr(G, "predictpror_source_kernel") <- "additive"
  attr(K, "predictpror_source_kernel") <- "gaussian"

  bank <- PredictProR:::gp_collect_kernel_inputs(gmatrix = G, gkernel = K)
  expect_length(bank, 2L)
  expect_setequal(names(bank), c("additive", "gaussian"))
})

test_that("a marker-derived GRM and gkernel remain separate model inputs", {
  ids <- paste0("g", 1:4)
  geno <- matrix(0:3, nrow = 4, dimnames = list(ids, "marker"))
  pheno <- data.frame(GID = rep(ids, 2), Env = rep(c("a", "b"), each = 4),
                      Yield = seq_len(8), Height = seq_len(8) + 1)
  G <- diag(4)
  COP <- diag(4) + 0.15
  dimnames(G) <- dimnames(COP) <- list(ids, ids)
  local_mocked_bindings(
    process_geno_data = function(...) list(gmatrix = G, geno_model_ready = geno),
    process_omic_data = function(...) NULL,
    grm_kernel_precheck = function(grm_kernel_data, ...) grm_kernel_data,
    pheno_geno_match = function(object_pheno, object_geno, ...)
      list(geno_pheno_match_data = object_geno),
    .package = "PredictProR"
  )
  inputs <- gp_prepare_model_input_objects(list(
    pheno_clean = list(pheno_clean_data = pheno),
    geno_data = geno, gen_name = "GID", response = c("Yield", "Height"),
    response_family = "gaussian", heter_groups = "Env",
    gmatrix_method = "Yang", gmatrix = NULL, gkernel = COP,
    GS_model = "GBLUP", AI_valid_models = character(),
    message = FALSE, impute = FALSE, qc_filtering = FALSE
  ))
  bank <- inputs$gmatrix_kernel_model_ready_list
  expect_setequal(names(bank), c("gmatrix_model_ready", "gkernel_model_ready"))
  expect_equal(as.numeric(bank$gmatrix_model_ready), as.numeric(G))
  expect_equal(as.numeric(bank$gkernel_model_ready), as.numeric(COP))
  collected <- gp_collect_kernel_inputs(gmatrix = bank$gmatrix_model_ready,
                                        kernel_list = bank)
  expect_length(collected, 2L)
  expect_setequal(names(collected), c("gmatrix", "gkernel"))
})
