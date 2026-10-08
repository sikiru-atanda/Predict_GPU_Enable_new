make_guardrail_ctx <- function(pheno, gmatrix = NULL, ...) {
  extra <- list(...)
  modifyList(
    list(
      pheno_clean = list(pheno_clean_data = pheno),
      gen_name = "GID",
      heter_groups = "Env",
      engine = "asreml",
      GS_model = "GBLUP",
      GS_model_cv = NULL,
      multi_trait_asreml = FALSE,
      hybrid_asreml = FALSE,
      gmatrix = gmatrix,
      gkernel = NULL,
      omic1_kernel = NULL,
      omic2_kernel = NULL,
      omic3_kernel = NULL,
      kernel_list = NULL,
      env_similarity = NULL
    ),
    extra
  )
}

test_that("model-input guardrails order ASReml phenotypes and coerce explicit kernels", {
  pheno <- data.frame(
    GID = c("g2", "g1", "g2", "g1"),
    Env = c("E2", "E2", "E1", "E1"),
    Yield = c(4, 2, 3, 1),
    stringsAsFactors = FALSE
  )
  rownames(pheno) <- paste0("r", seq_len(nrow(pheno)))

  kernel <- data.frame(
    g2 = c("1", "0.2"),
    g1 = c("0.2", "1"),
    row.names = c("g2", "g1"),
    check.names = FALSE
  )

  out <- PredictProR:::gp_apply_model_input_guardrails(
    make_guardrail_ctx(pheno, gmatrix = kernel)
  )

  pheno_out <- out$pheno_clean$pheno_clean_data
  expect_identical(
    paste(pheno_out$Env, pheno_out$GID, sep = ":"),
    c("E1:g1", "E1:g2", "E2:g1", "E2:g2")
  )
  expect_identical(rownames(pheno_out), as.character(seq_len(nrow(pheno_out))))
  expect_true(is.matrix(out$gmatrix))
  expect_true(is.numeric(out$gmatrix))
  expect_identical(rownames(out$gmatrix), colnames(out$gmatrix))
  expect_equal(out$gmatrix["g2", "g1"], 0.2)
})

test_that("model-input guardrails normalize multiple explicit omic kernels", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, NA),
    stringsAsFactors = FALSE
  )
  kernel_a <- data.frame(
    g2 = c(1, 0.3),
    g1 = c(0.3, 1),
    row.names = c("g2", "g1"),
    check.names = FALSE
  )
  kernel_b <- data.frame(
    g1 = c("1", "0.1"),
    g2 = c("0.1", "1"),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )

  out <- PredictProR:::gp_apply_model_input_guardrails(
    make_guardrail_ctx(
      pheno,
      engine = "not_asreml",
      GS_model = "RKHS",
      gmatrix = NULL,
      omic1_kernel = kernel_a,
      kernel_list = list(extra = kernel_b)
    )
  )

  expect_true(is.matrix(out$omic1_kernel))
  expect_true(is.numeric(out$omic1_kernel))
  expect_identical(rownames(out$omic1_kernel), colnames(out$omic1_kernel))
  expect_true(is.list(out$kernel_list))
  expect_true(is.matrix(out$kernel_list$extra))
  expect_true(is.numeric(out$kernel_list$extra))
})

test_that("model-input guardrails normalize environment kernels with env_ids", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Env = c("E2", "E1"),
    Yield = c(1, 2),
    stringsAsFactors = FALSE
  )
  env_kernel <- matrix(c("1", "0.4", "0.4", "1"), nrow = 2)

  out <- PredictProR:::gp_apply_model_input_guardrails(
    make_guardrail_ctx(
      pheno,
      engine = "not_asreml",
      GS_model = "GP",
      gmatrix = NULL,
      env_similarity = env_kernel,
      env_ids = c("E1", "E2")
    )
  )

  expect_true(is.matrix(out$env_similarity))
  expect_true(is.numeric(out$env_similarity))
  expect_identical(rownames(out$env_similarity), c("E1", "E2"))
  expect_identical(colnames(out$env_similarity), c("E1", "E2"))
})

test_that("model-ready feature guardrails normalize direct genotype and omic features", {
  geno_features <- data.frame(
    m1 = c("0", "2"),
    m2 = c("1", "1"),
    row.names = c("g2", "g1"),
    check.names = FALSE
  )
  omic_features <- data.frame(
    p1 = c(0.5, 0.7),
    p2 = c(1.5, 1.7),
    row.names = c("g2", "g1"),
    check.names = FALSE
  )

  out <- PredictProR:::gp_normalize_feature_bank_guardrail(
    list(geno_model_ready = geno_features, omic1_model_ready = omic_features),
    "model_ready_feature"
  )

  expect_true(is.matrix(out$geno_model_ready))
  expect_true(is.numeric(out$geno_model_ready))
  expect_identical(rownames(out$geno_model_ready), c("g2", "g1"))
  expect_true(is.matrix(out$omic1_model_ready))
  expect_true(is.numeric(out$omic1_model_ready))
})

test_that("model-ready feature guardrails use gen_name column when row names are not assigned", {
  geno_features <- data.frame(
    GID = c("g1", "g2", "g3"),
    m1 = c("0", "1", "2"),
    m2 = c("2", "1", "0"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_normalize_feature_bank_guardrail(
    list(geno_model_ready = geno_features),
    "model_ready_feature",
    id_col = "GID"
  )

  expect_true(is.matrix(out$geno_model_ready))
  expect_true(is.numeric(out$geno_model_ready))
  expect_identical(rownames(out$geno_model_ready), c("g1", "g2", "g3"))
  expect_false("GID" %in% colnames(out$geno_model_ready))
})

test_that("model-input guardrails reject malformed explicit kernels with guidance", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, 2),
    stringsAsFactors = FALSE
  )
  bad_kernel <- data.frame(
    g1 = c("bad", "0.1"),
    g2 = c("0.1", "1"),
    row.names = c("g1", "g2"),
    check.names = FALSE
  )

  expect_error(
    PredictProR:::gp_apply_model_input_guardrails(
      make_guardrail_ctx(pheno, gmatrix = bad_kernel)
    ),
    "gmatrix.*numeric"
  )
})

test_that("model_execute applies guardrails before ordinary ASReml dispatch", {
  pheno <- data.frame(
    GID = c("g2", "g1", "g2", "g1"),
    Env = c("E2", "E2", "E1", "E1"),
    Yield = c(4, 2, 3, 1),
    stringsAsFactors = FALSE
  )
  kernel <- data.frame(
    g2 = c("1", "0.2"),
    g1 = c("0.2", "1"),
    row.names = c("g2", "g1"),
    check.names = FALSE
  )

  local_mocked_bindings(
    gp_prepare_model_input_objects = function(ctx) {
      list(
        low_call_rate_inds_removed = NULL,
        pheno_clean = ctx$pheno_clean,
        geno_res = list(),
        omic1_res = list(),
        omic2_res = list(),
        omic3_res = list(),
        geno_omic_model_ready_list = list(),
        gmatrix_kernel_model_ready_list = list(gmatrix_model_ready = ctx$gmatrix),
        test_set = NULL,
        ml_dat_res = list()
      )
    },
    gp_route_true_prediction_specialized_models = function(ctx) {
      list(
        pheno_order = paste(
          ctx$pheno_clean$pheno_clean_data$Env,
          ctx$pheno_clean$pheno_clean_data$GID,
          sep = ":"
        ),
        gmatrix_is_matrix = is.matrix(ctx$gmatrix),
        gmatrix_is_numeric = is.numeric(ctx$gmatrix),
        gmatrix_dimnames_match = identical(rownames(ctx$gmatrix), colnames(ctx$gmatrix))
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR::model_execute(
    pheno_data = pheno,
    gmatrix = kernel,
    response = "Yield",
    gen_name = "GID",
    GS_model = "GBLUP",
    engine = "asreml",
    fixed = ~ Env,
    random = ~ GID + GID:Env,
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "fa1",
    eval_metrics = "accuracy",
    system_database = TRUE,
    message = FALSE
  )

  expect_identical(out$pheno_order, c("E1:g1", "E1:g2", "E2:g1", "E2:g2"))
  expect_true(out$gmatrix_is_matrix)
  expect_true(out$gmatrix_is_numeric)
  expect_true(out$gmatrix_dimnames_match)
})
