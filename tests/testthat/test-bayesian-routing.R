test_that("model_prep_bayes_cv preserves distinct keys for BRR and GBLUP_BRR", {
  local_mocked_bindings(
    bayes_finalize_A_B_C_BL_BRR = function(...) list(
      bayes_ETA = list(ETA = list(list(tag = "BRR"))),
      bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
    ),
    bayes_finalize_RKHS_GBLUPBRR = function(GS_model = NULL, ...) list(
      bayes_ETA = list(ETA = list(list(tag = GS_model))),
      bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
    ),
    .package = "PredictProR"
  )

  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  geno <- matrix(1:9, nrow = 3, dimnames = list(pheno$GID, paste0("m", 1:3)))
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- pheno$GID

  res <- PredictProR::model_prep_bayes_cv(
    random = ~ GID,
    GS_model_cv = c("BRR", "GBLUP_BRR"),
    response = "Yield",
    gen_name = "GID",
    pheno_data = pheno,
    geno_data = geno,
    gmatrix = grm
  )

  expect_identical(sort(names(res)), c("BRR", "GBLUP_BRR"))
  expect_identical(res[["BRR"]][["bayes_ETA"]][["ETA"]][[1]][["tag"]], "BRR")
  expect_identical(res[["GBLUP_BRR"]][["bayes_ETA"]][["ETA"]][[1]][["tag"]], "BRR")
})

test_that("model_prep_bayes_cv filters only available Bayesian input blocks", {
  local_mocked_bindings(
    bayes_finalize_A_B_C_BL_BRR = function(geno_data = NULL, omic1_data = NULL, ...) {
      expect_identical(rownames(geno_data), c("g1", "g3"))
      expect_null(omic1_data)
      list(
        bayes_ETA = list(ETA = list(list(tag = "BRR"))),
        bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
      )
    },
    bayes_finalize_RKHS_GBLUPBRR = function(gmatrix = NULL, omic1_kernel = NULL, ...) {
      expect_identical(rownames(gmatrix), c("g1", "g3"))
      expect_identical(colnames(gmatrix), c("g1", "g3"))
      expect_null(omic1_kernel)
      list(
        bayes_ETA = list(ETA = list(list(tag = "BRR"))),
        bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
      )
    },
    .package = "PredictProR"
  )

  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  geno <- matrix(1:9, nrow = 3, dimnames = list(pheno$GID, paste0("m", 1:3)))
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- pheno$GID

  res <- PredictProR::model_prep_bayes_cv(
    random = ~ GID,
    GS_model_cv = c("BRR", "GBLUP_BRR"),
    response = "Yield",
    gen_name = "GID",
    pheno_data = pheno,
    geno_data = geno,
    gmatrix = grm,
    test_set = "g2"
  )

  expect_identical(sort(names(res)), c("BRR", "GBLUP_BRR"))
})

test_that("bayes_finalize_RKHS_GBLUPBRR accepts GBLUP_BRR alias in prep mode", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- pheno$GID

  res <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID,
    GS_model = "GBLUP_BRR",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = grm,
    gen_name = "GID",
    nIter = 20,
    burnIn = 5,
    thin = 1,
    cross_validation = TRUE
  )

  expect_identical(res$bayes_ETA$ETA[[1]]$model, "BRR")
})

test_that("Bayesian kernel prep rejects GBLUP as a Bayesian alias", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- pheno$GID

  expect_error(
    PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
      random = ~ GID,
      GS_model = "GBLUP",
      response = "Yield",
      pheno_data = pheno,
      gmatrix = grm,
      gen_name = "GID",
      nIter = 20,
      burnIn = 5,
      thin = 1,
      cross_validation = TRUE
    ),
    "Use 'GBLUP' through the ASReml/GBLUP route",
    fixed = TRUE
  )

  prep <- PredictProR::model_prep_bayes_cv(
    random = ~ GID,
    GS_model_cv = "GBLUP",
    response = "Yield",
    gen_name = "GID",
    pheno_data = pheno,
    gmatrix = grm
  )

  expect_identical(prep, list())
})

test_that("bayes_parameter_check validates MCMC counts and silent defaults", {
  expect_silent(
    res <- PredictProR::bayes_parameter_check(message = FALSE)
  )
  expect_identical(res$nIter, 16000L)
  expect_identical(res$burnIn, 2000L)
  expect_identical(res$thin, 5L)
  expect_true(attr(res, "complete_parameters"))

  expect_error(
    PredictProR::bayes_parameter_check(nIter = 10, burnIn = 10, thin = 1, message = FALSE),
    "burnIn must be smaller than nIter"
  )
  expect_error(
    PredictProR::bayes_parameter_check(nIter = 10.5, burnIn = 1, thin = 1, message = FALSE),
    "nIter must be an integer"
  )
})

test_that("fixed ETA compilation accepts absent fixed effects", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_no_error(
    eta <- PredictProR::ETA_compiler_fixed_term(
      fixed = NULL,
      pheno_data = pheno
    )
  )
  expect_identical(eta, list())

  expect_no_error(
    eta_from_model_label <- PredictProR::ETA_compiler_fixed_term(
      fixed = "FIXED",
      fixed_term_model_bayesian = "FIXED",
      pheno_data = pheno
    )
  )
  expect_identical(eta_from_model_label, list())
})

test_that("fixed ETA compilation is treatment coded against the BGLR intercept", {
  pheno <- data.frame(
    GID = paste0("g", 1:6),
    Env = factor(rep(c("E1", "E2", "E3"), each = 2L)),
    Covariate = seq_len(6),
    Constant = 1,
    stringsAsFactors = FALSE
  )
  eta_env <- PredictProR::ETA_compiler_fixed_term(
    fixed = ~ Env,
    pheno_data = pheno
  )
  expect_length(eta_env, 1L)
  expect_identical(ncol(eta_env[[1]]$X), 2L)
  expect_false("(Intercept)" %in% colnames(eta_env[[1]]$X))

  eta_numeric <- PredictProR::ETA_compiler_fixed_term(
    fixed = ~ Covariate,
    pheno_data = pheno
  )
  expect_identical(ncol(eta_numeric[[1]]$X), 1L)

  eta_constant <- PredictProR::ETA_compiler_fixed_term(
    fixed = ~ Constant,
    pheno_data = pheno
  )
  expect_identical(eta_constant, list())
})

test_that("single-environment ASReml random setup ignores heter residual controls", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_no_error(
    out <- PredictProR::random_terms_fit_new(
      random = ~ GID,
      heter_groups = "Env",
      heter_resid = TRUE,
      var_cov_str = "fa1",
      code_asr = c("asreml::asreml(fixed=trait~1", "random=~", "residual=~"),
      names_in_inv_list = "gmatrix",
      gen_name = "GID",
      pheno_data = pheno
    )
  )
  expect_false(grepl("Env|fa\\(|idv\\(", out$code_asr[[2]]))
})

test_that("ASReml fixed formula builder drops empty trailing plus terms", {
  expect_identical(
    PredictProR:::asreml_clean_fixed_code("asreml::asreml(fixed= HT~1+"),
    "asreml::asreml(fixed= HT~1"
  )
  expect_identical(
    PredictProR:::asreml_append_rhs_terms("asreml::asreml(fixed= HT~1", character()),
    "asreml::asreml(fixed= HT~1"
  )
  expect_identical(
    PredictProR:::asreml_append_rhs_terms("asreml::asreml(fixed= HT~1", c("Env", "")),
    "asreml::asreml(fixed= HT~1+Env"
  )
  expect_null(PredictProR:::asreml_normalize_formula_arg("FIXED", labels = "FIXED"))
  expect_null(PredictProR:::asreml_normalize_formula_arg("COVA", labels = "COVA"))
})

test_that("Bayesian CV prep clears invalid heter groups for single-environment data", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- pheno$GID
  seen <- list()

  local_mocked_bindings(
    bayes_finalize_RKHS_GBLUPBRR = function(heter_groups = NULL, heter_resid = NULL, ...) {
      seen$heter_groups <<- heter_groups
      seen$heter_resid <<- heter_resid
      list(
        bayes_ETA = list(ETA = list(list(tag = "RKHS"))),
        bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
      )
    },
    .package = "PredictProR"
  )

  res <- PredictProR::model_prep_bayes_cv(
    random = ~ GID,
    GS_model_cv = "RKHS",
    response = "Yield",
    gen_name = "GID",
    pheno_data = pheno,
    gmatrix = grm,
    heter_groups = "Env",
    heter_resid = TRUE
  )

  expect_true("RKHS" %in% names(res))
  expect_null(seen$heter_groups)
  expect_false(seen$heter_resid)
})

test_that("Bayesian preprocessing defaults random models and validates marker finalizer model names", {
  expect_identical(
    PredictProR::random_term_model(rand_terms = "GID", gen_name = "GID", message = FALSE),
    "BRR"
  )

  expect_error(
    PredictProR:::bayes_finalize_A_B_C_BL_BRR(
      random = ~ GID,
      GS_model = "RKHS",
      response = "Yield",
      pheno_data = data.frame(GID = "g1", Yield = 1),
      geno_data = matrix(1, nrow = 1, dimnames = list("g1", "m1")),
      gen_name = "GID",
      cross_validation = TRUE
    ),
    "Bayesian marker finalization supports GS_model values"
  )
})

test_that("multi-environment marker-regression Bayes CV fails explicitly", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g1", "g2"),
    Env = c("E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  geno <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("m1", "m2")))

  expect_error(
    PredictProR::model_prep_bayes_cv(
      random = ~ GID + GID:Env,
      GS_model_cv = "BayesA",
      response = "Yield",
      gen_name = "GID",
      pheno_data = pheno,
      geno_data = geno,
      heter_groups = "Env"
    ),
    "single-environment only"
  )
})

test_that("kernel precheck does not fail when no duplicate pairs are present", {
  K <- matrix(
    c(1.0, 0.2, 0.1,
      0.2, 1.0, 0.3,
      0.1, 0.3, 1.0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("g1", "g2", "g3"))
  )

  expect_no_error(
    res <- PredictProR:::grm_kernel_precheck(K)
  )
  expect_true(is.matrix(res))
})

test_that("cross-validation uses the exact Bayesian prep key for GBLUP_BRR", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  seen_tag <- NULL

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) lapply(seq_along(model_ids), run_task),
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, tst = NULL, additional_params = NULL, ...) {
      seen_tag <<- additional_params$ETA[[1]]$tag
      rep(mean(y[-tst], na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Yield",
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = "mean_squared_error",
      GS_model_cv = "GBLUP_BRR",
      model_prep_all_bayes_cv = list(
        BRR = list(
          bayes_ETA = list(ETA = list(list(tag = "BRR"))),
          bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
        ),
        GBLUP_BRR = list(
          bayes_ETA = list(ETA = list(list(tag = "GBLUP_BRR"))),
          bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
        )
      )
    ),
    verbose = FALSE
  )

  expect_length(res, 1)
  expect_identical(seen_tag, "GBLUP_BRR")
})

test_that("Bayesian cross-validation never initializes an embedded Python", {
  pheno <- data.frame(GID = c("g1", "g2", "g3"), Yield = c(10, 20, 30), stringsAsFactors = FALSE)
  fake_py <- tempfile("python-")
  file.create(fake_py)

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) lapply(seq_along(model_ids), run_task),
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, tst = NULL, ...) rep(mean(y[-tst]), length(tst)),
    gp_detect_python = function(...) fake_py,
    # A python3 without a shared libpython cannot be embedded; BGLR needs none.
    gp_init_python_once = function(py_bin) stop("reticulate init attempted for an R-only model"),
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Yield", gen_name = "GID", test_size = 0.5, random_state = 1,
      replication = 1L, cross_validation_meth = "hold_out", eval_metrics = "mean_squared_error",
      GS_model_cv = "BRR",
      model_prep_all_bayes_cv = list(BRR = list(
        bayes_ETA = list(ETA = list(list(tag = "BRR"))),
        bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
      ))
    ),
    verbose = FALSE
  )
  expect_length(res, 1)
})

test_that("Bayesian true prediction never initializes an embedded Python", {
  skip_if_not_installed("BGLR")
  set.seed(3)
  G <- matrix(rbinom(20 * 15, 2, 0.4), 20, dimnames = list(sprintf("g%02d", 1:20), sprintf("m%02d", 1:15)))
  ph <- data.frame(GID = rownames(G), Yield = rnorm(20))
  ph$Yield[1:3] <- NA
  withr::local_dir(withr::local_tempdir())
  local_mocked_bindings(
    gp_init_python_once = function(py_bin) stop("reticulate init attempted for an R-only model"),
    .package = "PredictProR"
  )
  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "GBLUP_BRR", nIter = 60L, burnIn = 20L, thin = 2L, qc_filtering = FALSE, impute = FALSE,
    ld_prunning_qc = FALSE, gmatrix_method = "VanRaden", parallel_mode = "sequential",
    system_database = TRUE, message = FALSE
  ))
  expect_equal(nrow(out$model_results$predicted_values), 20L)
})

test_that("Bayesian cross-validation uses trait-specific test indices under unbalanced traits", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(10, NA, 30, 40),
    Height = c(NA, 20, 30, 40),
    stringsAsFactors = FALSE
  )

  seen_tst <- list()

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) lapply(seq_along(model_ids), run_task),
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, tst = NULL, additional_params = NULL, ...) {
      seen_tst[[additional_params$bayes_trait]] <<- tst
      rep(mean(y[-tst], na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = c("Yield", "Height"),
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = "mean_squared_error",
      GS_model_cv = "BRR",
      model_prep_all_bayes_cv = list(
        BRR = list(
          bayes_ETA = list(ETA = list(list(tag = "BRR"))),
          bayes_para = list(nIter = 10L, burnIn = 2L, thin = 1L)
        )
      )
    ),
    verbose = FALSE
  )

  expect_length(res, 2)
  expect_identical(seen_tst[["Yield"]], 1L)
  expect_identical(seen_tst[["Height"]], 2L)
})

test_that("random_term_model expands a scalar GS_model cleanly across MET terms", {
  expect_no_warning(
    res <- PredictProR::random_term_model(
      rand_terms = c("GID", "GID:Env"),
      GS_model = "RKHS",
      gen_name = "GID"
    )
  )

  expect_identical(res, c("RKHS", "RKHS"))
})

test_that("Bayesian summary statistics receive heter_groups for multi-environment RKHS", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g1", "g2"),
    Env = c("E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  seen_heter_groups <- NULL

  local_mocked_bindings(
    bayes_finalize_RKHS_GBLUPBRR = function(...) {
      list(
        bayes_result = list(pred = data.frame(x = 1)),
        bayes_model = list(model = list(y = pheno$Yield, yHat = pheno$Yield, varE = 1, ETA = list()))
      )
    },
    summary_statistics_bayes = function(heter_groups = NULL, ...) {
      seen_heter_groups <<- heter_groups
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )

  ctx <- list(
    response = "Yield",
    GS_model = "RKHS",
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    rand_term_model_bayesian = NULL,
    pheno_clean = list(pheno_clean_data = pheno),
    gen_name = "GID",
    heter_groups = "Env",
    fixed = NULL,
    random = ~ GID + GID:Env,
    weights = NULL,
    fixed_term_model_bayesian = NULL,
    gmatrix_model_ready = diag(2),
    omic1_kernel_model_ready = NULL,
    omic2_kernel_model_ready = NULL,
    omic3_kernel_model_ready = NULL,
    nIter = 10L,
    burnIn = 2L,
    thin = 1L,
    omics_kernel_label = NULL,
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.7,
    low_reliability_thres = 0.4,
    n_components = 2L,
    threshold = 10,
    interval_width_low_threshold = NULL,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    eval_metrics = "mean_squared_error",
    system_database = FALSE,
    friendly_name_lookup = c(RKHS = "RKHS"),
    ml_dat_res = list(),
    msg = "",
    cross_validation = FALSE
  )
  rownames(ctx$gmatrix_model_ready) <- colnames(ctx$gmatrix_model_ready) <- c("g1", "g2")

  res <- PredictProR:::gp_run_best_model_task(ctx)

  expect_identical(seen_heter_groups, "Env")
  expect_true("summary_statistics" %in% names(res$res_summary_stat))
})

test_that("Bayesian best-model routing preserves trait-specific true-prediction labels", {
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

  seen_test_ids <- list()

  local_mocked_bindings(
    bayes_finalize_A_B_C_BL_BRR = function(pheno_data = NULL, response = NULL, ...) {
      lbl <- ifelse(is.na(pheno_data[[response]]), "Test", "Train")
      seen_test_ids[[response]] <<- as.character(pheno_data$GID[lbl == "Test"])
      list(
        bayes_result = list(
          Predicted_value = data.frame(
            GID = pheno_data$GID,
            Predicted_value = seq_len(nrow(pheno_data)),
            Train_Test_Label = lbl,
            stringsAsFactors = FALSE
          )
        ),
        bayes_model = list(model = list(y = pheno_data[[response]], yHat = rep(0, nrow(pheno_data)), varE = 1, ETA = list()))
      )
    },
    summary_statistics_bayes = function(...) list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
    .package = "PredictProR"
  )

  base_ctx <- list(
    GS_model = "BayesA",
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    AI_valid_models = character(),
    canonical_names = character(),
    pheno_clean = pheno_clean,
    ml_dat_res = list(),
    gen_name = "GID",
    msg = "",
    fixed = NULL,
    random = ~ GID,
    weights = NULL,
    fixed_term_model_bayesian = NULL,
    rand_term_model_bayesian = NULL,
    geno_model_ready = matrix(1:9, nrow = 3, dimnames = list(c("g1", "g2", "g3"), paste0("m", 1:3))),
    omic1_model_ready = NULL,
    omic2_model_ready = NULL,
    omic3_model_ready = NULL,
    nIter = 10L,
    burnIn = 2L,
    thin = 1L,
    omics_data_label = NULL,
    scaling = FALSE,
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.7,
    low_reliability_thres = 0.4,
    n_components = 2L,
    threshold = 10,
    interval_width_high_threshold = NULL,
    interval_width_low_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    eval_metrics = "mean_squared_error",
    heter_groups = NULL,
    system_database = FALSE,
    friendly_name_lookup = c(BayesA = "BayesA"),
    cross_validation = FALSE
  )

  yield_res <- PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Yield")))
  height_res <- PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Height")))

  expect_identical(unname(seen_test_ids[["Yield"]]), "g2")
  expect_identical(unname(seen_test_ids[["Height"]]), "g1")
  expect_identical(
    as.character(yield_res$res_model_output$bayes_result$Predicted_value$Train_Test_Label),
    c("Train", "Test", "Train")
  )
  expect_identical(
    as.character(height_res$res_model_output$bayes_result$Predicted_value$Train_Test_Label),
    c("Test", "Train", "Train")
  )
})

test_that("results_handling exports non-bayes_result extras alongside bayes outputs", {
  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-bayes-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  out <- PredictProR::results_handling(
    GS_model = "BayesA",
    res_model_output = list(
      bayes_result = list(
        Predicted_value = data.frame(GID = "g1", Predicted_value = 1)
      ),
      bayes_model = list(model = "dummy"),
      Connectivity_summary = data.frame(metric = "env_count", value = 2)
    ),
    res_summary_stat = list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
    system_database = FALSE,
    plot_filename = "bayes_export_test"
  )

  created_dir <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_identical(out$export_status, "Successful")
  expect_length(created_dir, 1)
  expect_true(file.exists(file.path(created_dir, "Predicted_Value.csv")))
  expect_true(file.exists(file.path(created_dir, "Connectivity_summary.csv")))
})

test_that("RKHS EBV helper uses ETA latent effects instead of fitted-response SD", {
  mod <- list(
    model = list(
      ETA = list(
        list(
          u = c(0.2, 0.5, 0.8),
          SD.u = c(0.1, 0.2, 0.3)
        ),
        list(
          u = c(1, 1, 1),
          SD.u = c(9, 9, 9)
        )
      ),
      yHat = c(10, 10, 10),
      SD.yHat = c(5, 5, 5)
    )
  )

  res <- PredictProR:::cal_coeff_ebv_pev_rel_RHKS_glub(
    mod = mod,
    gmatrix = diag(3),
    gen_name = "GID",
    var_u = 2,
    gid_name = c("g1", "g2", "g3"),
    eta_index = 1L
  )

  expect_equal(res$Estimated_breeding_value$Estimated_breeding_value, c(0.2, 0.5, 0.8))
  expect_equal(res$Standard_error, c(0.1, 0.2, 0.3))
  expect_equal(res$PEV, c(0.01, 0.04, 0.09))
})

test_that("ASReml prediction labels always expose PEV and reliability aliases", {
  pheno <- data.frame(
    NAME = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    NAME = c("g1", "g2", "g3"),
    Predicted_value = c(1.1, 2.1, 3.1),
    Standard_error = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::asreml_add_train_test_labels(
    predicted_df = pred,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "NAME"
  )

  expect_identical(out$Train_Test_Label, c("Train", "Test", "Train"))
  expect_true(all(c(
    "Prediction_error_variance",
    "PEV",
    "Reliability_reference_variance",
    "Reliability",
    "Reliability_remarks",
    "Reliability_percentage",
    "Reliability_basis"
  ) %in% names(out)))
  expect_equal(out$PEV, pred$Standard_error^2)
  expect_true(all(is.finite(out$Reliability)))

  pheno_no_missing <- within(pheno, Yield[is.na(Yield)] <- 2)
  out_no_missing <- PredictProR:::asreml_add_train_test_labels(
    predicted_df = pred,
    pheno_data = pheno_no_missing,
    response = "Yield",
    gen_name = "NAME"
  )
  expect_identical(out_no_missing$Train_Test_Label, rep("Train", 3))
  expect_true(all(c("PEV", "Reliability") %in% names(out_no_missing)))
  expect_true(all(is.finite(out_no_missing$Reliability)))
})

test_that("bayes variance components report posterior SD scale", {
  var_u_total <- c(1, 3, 5)
  res <- PredictProR:::bayes_variance_componentsnew(
    var_u_mean_omics_list = list(gmatrix = mean(var_u_total)),
    se_var_u_omics_list = list(gmatrix = stats::sd(var_u_total)),
    var_u_total = var_u_total,
    var_residual = c(2, 4, 6),
    se_var_residual = stats::sd(c(2, 4, 6))
  )

  expect_equal(
    res["heritability", "Standard_error"],
    stats::sd(var_u_total / (var_u_total + c(2, 4, 6)))
  )
})

test_that("Bayesian saved-draw indices remove burn-in on the thinned file scale", {
  expect_identical(
    PredictProR:::gp_bayes_saved_draw_indices(nIter = 120L, burnIn = 20L, thin = 2L),
    11:60
  )

  draw_dir <- tempfile("bayes-draws-")
  dir.create(draw_dir)
  path <- file.path(draw_dir, "varE.dat")
  writeLines(as.character(seq_len(60L)), path)
  out <- PredictProR:::gp_bayes_read_varE_draws(
    output_files_names = path,
    draw_indices = 11:60
  )
  expect_identical(out, as.double(11:60))
  expect_identical(
    PredictProR::process_var_u(path, posindex = 11:60, GS_model = "RKHS"),
    as.double(11:60)
  )
})

test_that("reliability thresholds respect caller-supplied cutoffs", {
  res <- PredictProR:::reliability_thresholds(
    prediction_error_var = c(0.1, 0.4, 0.8),
    genetic_var = 1,
    high_reliability_thres = 0.85,
    low_reliability_thres = 0.3
  )

  expect_identical(res$remarks, c("Reliable", "Acceptable", "Unreliable"))
})

test_that("MPIW summary computes a non-zero medium bucket when widths fall in between", {
  res <- PredictProR:::reliability_thresholds_MPIW_from_CI(
    predictions = c(1, 1, 1),
    standard_errors = c(0.1, 0.2, 0.3),
    CI_width_thresholds = c(0.33, 0.66),
    model_for_CI_cal = "RKHS"
  )

  expect_gt(res$proportion_medium_reliability, 0)
})

test_that("MET aggregated Bayesian predictions carry target-specific uncertainty summaries", {
  set.seed(5)
  pheno <- expand.grid(
    GID = paste0("g", 1:6),
    Env = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- rnorm(nrow(pheno))
  pheno$Yield[c(2, 9)] <- NA_real_

  gmatrix <- diag(6)
  rownames(gmatrix) <- colnames(gmatrix) <- paste0("g", 1:6)

  res <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2
  )

  total_pred <- res$bayes_result$Total_Predicted_value
  expect_true(is.data.frame(total_pred))
  expect_true(all(c("Predicted_value", "Standard_error", "PEV", "Reliability", "Genetic_variance") %in% names(total_pred)))
  expect_true(any(is.finite(as.double(total_pred$Predicted_value))))
  expect_true(any(is.finite(as.double(total_pred$Standard_error))))
  expect_true(any(is.finite(as.double(total_pred$PEV))))
  expect_true(any(is.finite(as.double(total_pred$Reliability))))
  expect_true(any(is.finite(as.double(total_pred$Genetic_variance))))
})

test_that("BGLR grouped residual variances are expanded to observation order", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g1", "g2"),
    Env = c("E1", "E1", "E2", "E2"),
    Yield = c(1, 2, NA, 4),
    stringsAsFactors = FALSE
  )
  X <- diag(4)
  eta <- list(list(X = X, model = "BRR"))
  beta_draws <- matrix(
    c(
      0.2, 0.4,
      0.1, 0.3,
      0.5, 0.6,
      0.8, 1.0
    ),
    nrow = 2,
    byrow = TRUE
  )
  b_file <- file.path(tempdir(), "ETA_1_b.bin")
  con <- gzfile(b_file, open = "wb")
  writeBin(as.numeric(nrow(beta_draws)), con, size = 8)
  writeBin(as.numeric(ncol(beta_draws)), con, size = 8)
  writeBin(as.vector(beta_draws), con, size = 8)
  close(con)
  on.exit(unlink(b_file), add = TRUE)

  mod <- list(
    model = list(
      y = pheno$Yield,
      yHat = c(1.1, 2.1, 3.1, 4.1),
      mu = 0,
      ETA = list(list(b = rep(0.1, 4), varB = 0.2)),
      varE = c(E1 = 0.5, E2 = 2.0),
      weights = rep(1, 4)
    )
  )

  res <- PredictProR:::bayes_compute_prediction_targets(
    mod = mod,
    ETA = list(ETA = eta),
    GS_model = "BRR",
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = "Env",
    output_files_names = b_file
  )

  expect_true(is.list(res))
  expect_true(is.data.frame(as.data.frame(res$across)))
  expect_equal(
    PredictProR:::gp_bayes_observation_residual_variance(
      varE = mod$model$varE,
      groups = pheno$Env,
      weights = mod$model$weights
    ),
    c(0.5, 0.5, 2.0, 2.0)
  )
})

test_that("Bayesian kernel finalizer routes GBLUP_BRR MET heter_resid through the Multitrait kernel-prior fit", {
  # Background: the prior implementation let GBLUP_BRR MET take BGLR's native
  # `groups` path via bayes_mod_execute. Empirically that collapsed the
  # genetic signal (h2 ~ 0.07 vs RKHS ~ 0.51 on the same K) because BGLR's
  # default sigma^2_beta prior scale does not calibrate to the eigen-sqrt
  # BRR design the same way the kernel prior does. The kernel finalizer now
  # routes GBLUP_BRR MET through bayes_multitrait_env_heter_fit (same as
  # RKHS) so both produce the same N(0, sigma^2 * K) genetic prior and
  # matching per-env variance components.
  pheno <- data.frame(
    GID = c("g1", "g2", "g1", "g2"),
    Env = c("E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, NA),
    stringsAsFactors = FALSE
  )
  grm <- diag(2)
  rownames(grm) <- colnames(grm) <- c("g1", "g2")

  seen <- list(called = FALSE, GS_model = NULL)
  local_mocked_bindings(
    bayes_multitrait_env_heter_fit = function(GS_model = NULL, ...) {
      seen$called <<- TRUE
      seen$GS_model <<- GS_model
      list(
        bayes_result = list(Predicted_value = data.frame(ok = TRUE)),
        bayes_model = list(
          model = list(y = pheno$Yield, yHat = pheno$Yield, varE = 1, ETA = list()),
          output_files_names = character()
        )
      )
    },
    .package = "PredictProR"
  )

  PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "GBLUP_BRR",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = grm,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    nIter = 20,
    burnIn = 5,
    thin = 1
  )

  expect_true(seen$called)
  # GBLUP_BRR is normalised to "BRR" inside the finalizer; the multitrait
  # fitter accepts either label and uses it only for output naming.
  expect_identical(seen$GS_model, "BRR")
})

test_that("Bayesian RKHS finalizer treats one environment with replicates as single-environment", {
  pheno <- data.frame(
    GID = rep(c("g1", "g2", "g3"), each = 2),
    Env = "E1",
    Yield = c(1.0, NA, 2.0, 2.2, 3.1, 3.0),
    stringsAsFactors = FALSE
  )
  grm <- diag(3)
  rownames(grm) <- colnames(grm) <- c("g1", "g2", "g3")

  seen_groups <- "not called"
  local_mocked_bindings(
    bayes_mod_execute = function(groups = NULL, ...) {
      seen_groups <<- groups
      list(model = list(y = pheno$Yield, yHat = pheno$Yield, varE = 1, ETA = list()), output_files_names = character())
    },
    mod_output_bayes_gbluBRR_RKHS = function(...) list(Predicted_value = data.frame(ok = TRUE)),
    .package = "PredictProR"
  )

  res <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = grm,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    nIter = 20,
    burnIn = 5,
    thin = 1
  )

  expect_true(is.list(res))
  expect_null(seen_groups)
})

test_that("single-environment RKHS output keeps variance, SE, PEV, and reliability", {
  skip_if_not_installed("BGLR")
  set.seed(52)
  gids <- paste0("g", 1:5)
  grm <- diag(length(gids))
  rownames(grm) <- colnames(grm) <- gids
  pheno <- data.frame(
    GID = rep(gids, each = 2),
    Env = "E1",
    Yield = rnorm(length(gids) * 2),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(2, 7)] <- NA_real_

  res <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = grm,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    nIter = 80,
    burnIn = 20,
    thin = 2
  )

  required_cols <- c(
    "Predicted_value", "Standard_error", "PEV",
    "Genetic_variance", "Reliability", "Reliability_percentage"
  )
  expect_true(all(required_cols %in% names(res$bayes_result$Predicted_value)))
  expect_true(all(required_cols %in% names(res$bayes_result$Total_Predicted_value)))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$Standard_error)))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$PEV)))
  expect_true(all(res$bayes_result$Predicted_value$PEV >= 0))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$Genetic_variance)))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$Reliability)))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$Standard_error)))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$PEV)))
  expect_true(all(res$bayes_result$Total_Predicted_value$PEV >= 0))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$Genetic_variance)))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$Reliability)))
  expect_identical(dim(res$bayes_result$Residual_covariance_environments), c(1L, 1L))
  expect_identical(dim(res$bayes_result$Residual_covariance_environments_SE), c(1L, 1L))
  expect_true(all(is.finite(res$bayes_result$Residual_covariance_environments)))
  expect_true(all(is.finite(res$bayes_result$Residual_covariance_environments_SE)))
})

test_that("Bayesian CV prep carries grouped residual labels for supported GBLUP_BRR ETA", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g1", "g2"),
    Env = c("E1", "E1", "E2", "E2"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  grm <- diag(2)
  rownames(grm) <- colnames(grm) <- c("g1", "g2")

  res <- PredictProR::model_prep_bayes_cv(
    random = ~ GID + GID:Env,
    GS_model_cv = "GBLUP_BRR",
    response = "Yield",
    gen_name = "GID",
    pheno_data = pheno,
    gmatrix = grm,
    heter_groups = "Env",
    heter_resid = TRUE,
    nIter = 20,
    burnIn = 5,
    thin = 1
  )

  expect_identical(as.character(res$GBLUP_BRR$bayes_groups), pheno$Env)
})

test_that("BGLR Multitrait RKHS heter residual output carries within and across SEs", {
  skip_if_not_installed("BGLR")
  set.seed(11)
  gids <- paste0("g", 1:5)
  envs <- c("E1", "E2")
  K <- diag(length(gids))
  rownames(K) <- colnames(K) <- gids
  K2 <- diag(seq(0.7, 1.3, length.out = length(gids)))
  rownames(K2) <- colnames(K2) <- gids
  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  pheno$Yield <- rnorm(nrow(pheno))
  pheno$Yield[c(2, 8)] <- NA_real_

  res <- PredictProR:::bayes_multitrait_env_heter_fit(
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    kernels = list(gmatrix = K, transcriptome = K2),
    GS_model = "RKHS",
    bayes_para = list(nIter = 80L, burnIn = 20L, thin = 2L),
    confidence_level = 0.95
  )

  expect_true(is.data.frame(res$bayes_result$Predicted_value))
  expect_true(is.data.frame(res$bayes_result$Total_Predicted_value))
  expect_true(all(c("Predicted_value", "Standard_error", "PEV") %in% names(res$bayes_result$Predicted_value)))
  expect_true(all(c("Predicted_value", "Standard_error", "PEV") %in% names(res$bayes_result$Total_Predicted_value)))
  expect_true(any(is.finite(res$bayes_result$Predicted_value$Standard_error)))
  expect_true(any(is.finite(res$bayes_result$Total_Predicted_value$Standard_error)))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$PEV)))
  expect_true(all(res$bayes_result$Predicted_value$PEV >= 0))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$Genetic_variance)))
  expect_true(all(is.finite(res$bayes_result$Predicted_value$Reliability)))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$PEV)))
  expect_true(all(res$bayes_result$Total_Predicted_value$PEV >= 0))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$Genetic_variance)))
  expect_true(all(is.finite(res$bayes_result$Total_Predicted_value$Reliability)))
  expect_true(is.data.frame(res$bayes_result$Residual_variance_by_environment))
  expect_equal(
    sort(res$bayes_result$Residual_variance_by_environment$Environment),
    sort(envs)
  )
  expect_true(all(is.finite(res$bayes_result$Residual_variance_by_environment$Residual_variance)))
  expect_true(all(res$bayes_result$Residual_variance_by_environment$Residual_variance > 0))
  expect_identical(
    names(res$bayes_result$Variance_components),
    c("Component", "Components", "Standard_error")
  )
  expect_true(all(c(
    "genetic_variance_E1", "genetic_variance_E2",
    "residual_variance_E1", "residual_variance_E2",
    "heritability_E1", "heritability_E2",
    "across_environment_genetic_variance",
    "across_environment_residual_variance",
    "across_environment_heritability"
  ) %in% res$bayes_result$Variance_components$Component))
  params <- stats::setNames(
    as.character(res$bayes_result$model_parameters$summary),
    as.character(res$bayes_result$model_parameters$stat)
  )
  expect_identical(params[["bayesian_kernel_model"]], "RKHS")
  expect_identical(params[["rkhs_environment_mode"]], "heterogeneous_environment")
  expect_identical(params[["rkhs_heterogeneous_environment_variances"]], "TRUE")
  expect_identical(params[["bayes_kernel_heter_resid_effective"]], "TRUE")
  expect_identical(params[["stage2_weight_supplied"]], "FALSE")
  expect_length(res$bayes_result$Genetic_covariance_by_kernel, 2L)
  expect_equal(
    Reduce(`+`, res$bayes_result$Genetic_covariance_by_kernel),
    res$bayes_result$Genetic_covariance_environments,
    tolerance = 1e-8
  )
  expect_true(all(res$bayes_result$kernel_variance_components$Independent_kernel_estimate))
  expect_true(all(is.finite(res$bayes_result$kernel_variance_components$Standard_error)))
})
