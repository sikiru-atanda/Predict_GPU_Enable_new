if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

configure_test_python_runtime <- function() {
  py_ml <- PredictProR:::gp_preferred_python(purpose = "ml")
  py_dl <- PredictProR:::gp_preferred_python(purpose = "dl")
  skip_if(is.null(py_ml) || !file.exists(py_ml), "No preferred ML Python runtime found for PredictProR")
  skip_if(is.null(py_dl) || !file.exists(py_dl), "No preferred DL Python runtime found for PredictProR")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(
    PREDICTPRO_PYTHON = py_ml,
    PREDICTPRO_DL_PYTHON = py_dl,
    PREDICTPRO_ML_BRIDGE_BACKEND = "cli"
  )
  invisible(list(ml = py_ml, dl = py_dl))
}

mk_runtime_data <- function(kind = c("gaussian", "binary", "multiclass"), n = 30L, p = 8L) {
  kind <- match.arg(kind)
  set.seed(20260424 + match(kind, c("gaussian", "binary", "multiclass")))
  # Responses come from a latent continuous signal; the genotype matrix carries
  # it as homozygous 0/2 dosages: model_execute() genotype QC rejects values
  # above 2 without ploidy and, for line data, removes markers with more than
  # 10% heterozygotes.
  Z <- matrix(rnorm(n * p), nrow = n)
  X <- ifelse(Z > 0, 2, 0)   # homozygous 0/2 coding, like inbred lines
  storage.mode(X) <- "double"
  rownames(X) <- paste0(substr(kind, 1, 1), seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))

  y <- switch(
    kind,
    gaussian = Z[, 1] - 0.5 * Z[, 2] + rnorm(n, sd = 0.3),
    binary = factor(ifelse(runif(n) < plogis(Z[, 1] - Z[, 2]), "yes", "no"), levels = c("no", "yes")),
    multiclass = {
      sc <- cbind(A = Z[, 1], B = Z[, 2], C = -Z[, 1] - Z[, 2]) + matrix(rnorm(n * 3, sd = 0.2), nrow = n)
      factor(colnames(sc)[max.col(sc, ties.method = "first")], levels = c("A", "B", "C"))
    }
  )

  pheno <- data.frame(GID = rownames(X), Trait = y, stringsAsFactors = FALSE)
  tst <- sample(seq_len(n), ceiling(n * 0.25))
  list(
    train_pheno = pheno[-tst, , drop = FALSE],
    test_pheno = pheno[tst, , drop = FALSE],
    train_geno = X[-tst, , drop = FALSE],
    test_geno = X[tst, , drop = FALSE],
    all_pheno = pheno,
    all_geno = X
  )
}

expect_prediction_table_shape <- function(out, expected_rows, family) {
  expect_true(is.list(out))
  expect_equal(nrow(out$predicted_values), expected_rows)
  expect_true("tune_metric" %in% out$model_parameters$stat)

  if (family == "gaussian") {
    expect_true(all(c(
      "Predicted_value", "Prediction_stability",
      "Prediction_risk_score", "Prediction_risk_basis",
      "Rank_instability_risk", "Rank_instability_risk_basis"
    ) %in% names(out$predicted_values)))
    expect_false("Genetic_variance" %in% names(out$predicted_values))
    expect_true(all(c(
      "Reliability_variance_input", "Reliability_reference_variance", "Reliability_basis",
      "Prediction_uncertainty_source"
    ) %in% names(out$predicted_values)))
    pev_available <- is.finite(out$predicted_values$PEV)
    expect_true(all(is.na(out$predicted_values$Reliability[!pev_available])))
    expect_true(all(is.finite(out$predicted_values$Prediction_stability)))
    if (any(pev_available)) {
      expect_true(all(is.finite(out$predicted_values$Reliability[pev_available])))
      expect_equal(
        out$predicted_values$Reliability_variance_input[pev_available],
        out$predicted_values$PEV[pev_available]
      )
      expect_equal(
        out$predicted_values$Reliability[pev_available],
        out$predicted_values$Reliability_reference_variance[pev_available] /
          (out$predicted_values$Reliability_reference_variance[pev_available] +
             out$predicted_values$PEV[pev_available])
      )
    }
    expect_true(all(grepl(
      "not genetic reliability",
      out$predicted_values$Reliability_basis,
      fixed = TRUE
    )))
    expect_true(all(c(
      "genetic_variance_not_identifiable",
      "residual_variance_not_identifiable",
      "prediction_variance",
      "estimated_prediction_mse"
    ) %in% out$variance_components$Component))
  } else if (family == "binary") {
    expect_true(all(c("Predicted_class", "Prediction_confidence") %in% names(out$predicted_values)))
  } else {
    expect_true(all(c("Predicted_class", "Prediction_confidence") %in% names(out$predicted_values)))
    expect_true(all(c("Prob_A", "Prob_B", "Prob_C") %in% names(out$predicted_values)))
  }
}

expect_exported_classification_contract <- function(pred, class_levels) {
  prob_cols <- paste0("Probability_", make.names(class_levels, unique = TRUE))
  expect_identical(
    names(pred),
    c(PredictProR:::gp_public_prediction_columns("multiclass"), prob_cols)
  )
  expect_true(all(pred$Train_Test_Label %in% c("Train", "Test")))
  expect_true(all(is.finite(pred$Prediction_confidence)))
  expect_true(all(is.finite(pred$Classification_uncertainty)))
  expect_true(all(is.finite(pred$Reliability)))
  expect_true(all(!is.na(pred$Reliability_remarks) & nzchar(pred$Reliability_remarks)))
  prob <- as.matrix(pred[, prob_cols, drop = FALSE])
  expect_true(all(is.finite(prob)))
  expect_true(all(prob >= 0 & prob <= 1))
  expect_equal(rowSums(prob), rep(1, nrow(pred)), tolerance = 1e-6)
  confidence <- apply(prob, 1L, max)
  expect_equal(pred$Prediction_confidence, confidence, tolerance = 1e-6)
  expect_equal(pred$Classification_uncertainty, 1 - confidence, tolerance = 1e-6)
  expect_equal(pred$Reliability, confidence, tolerance = 1e-6)
  expect_identical(
    pred$Predicted_class,
    class_levels[max.col(prob, ties.method = "first")]
  )
}

find_runtime_predicted_file <- function(export_dir) {
  files <- list.files(export_dir, recursive = FALSE, full.names = TRUE)
  pred_files <- files[
    grepl("^predicted(_value|_values)?\\.csv$", basename(files), ignore.case = TRUE) |
      grepl("^Predicted(_Value|_Values)?\\.csv$", basename(files), ignore.case = TRUE)
  ]
  if (!length(pred_files)) {
    return(NULL)
  }
  pred_files[[1L]]
}

test_that("python-backed ML runtime smoke covers gaussian, binary, and multiclass", {
  configure_test_python_runtime()

  gdat <- mk_runtime_data("gaussian")
  out_cat <- PredictProR:::AI_CatBoost(
    pheno_object = gdat$train_pheno,
    geno_omic_object = gdat$train_geno,
    geno_omic_test_object = gdat$test_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "gaussian",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    system_database = TRUE,
    catboost_paras_tunning = list(
      catboost_iterations = c(25L),
      catboost_depth = c(4L),
      catboost_learning_rate = c(0.05),
      catboost_l2_leaf_reg = c(3L)
    )
  )
  expect_prediction_table_shape(out_cat, nrow(gdat$all_geno), "gaussian")

  bdat <- mk_runtime_data("binary")
  out_xgb <- PredictProR:::AI_Xgb(
    pheno_object = bdat$train_pheno,
    geno_omic_object = bdat$train_geno,
    geno_omic_test_object = bdat$test_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "binary",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    system_database = TRUE,
    xgb_paras_tunning = list(
      Iter_tune = c(25L),
      learning_rate_tune = c(0.05),
      max_depth = c(3L),
      xgb_gamma = c(0),
      colsample_bytree = c(0.8),
      min_child_weight = c(1),
      subsample = c(0.8),
      L1_tune = c(0),
      L2_tune = c(1)
    )
  )
  expect_prediction_table_shape(out_xgb, nrow(bdat$all_geno), "binary")

  mdat <- mk_runtime_data("multiclass")
  out_lgb <- PredictProR:::AI_LightGBM(
    pheno_object = mdat$train_pheno,
    geno_omic_object = mdat$train_geno,
    geno_omic_test_object = mdat$test_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "multi-class",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    system_database = TRUE,
    lightgbm_paras_tunning = list(
      lightgbm_nrounds = c(25L),
      lightgbm_learning_rate = c(0.05),
      lightgbm_num_leaves = c(15L),
      lightgbm_feature_fraction = c(0.8),
      lightgbm_bagging_fraction = c(0.8),
      lightgbm_min_data_in_leaf = c(5L),
      lightgbm_lambda_l1 = c(0),
      lightgbm_lambda_l2 = c(0)
    )
  )
  expect_prediction_table_shape(out_lgb, nrow(mdat$all_geno), "multiclass")

  out_rf <- PredictProR:::AI_randomForest(
    pheno_object = bdat$train_pheno,
    geno_omic_object = bdat$train_geno,
    geno_omic_test_object = bdat$test_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "binary",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    system_database = TRUE,
    rf_paras_tunning = list(
      ntree = c(100L),
      mtry = c(3L),
      nodesize = c(1L),
      maxnodes = c(20L)
    )
  )
  expect_prediction_table_shape(out_rf, nrow(bdat$all_geno), "binary")
})

test_that("public tuned classification exports keep classification columns", {
  configure_test_python_runtime()

  old_wd <- getwd()
  tmp_dir <- tempfile("predictpror-runtime-export-")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  bdat <- mk_runtime_data("binary", n = 36L, p = 8L)
  pheno_bin <- bdat$all_pheno
  pheno_bin$Trait[sample(seq_len(nrow(pheno_bin)), 8L)] <- NA

  res_bin <- PredictProR::model_execute(
    pheno_data = pheno_bin,
    geno_data = bdat$all_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "binary",
    GS_model = "CatBoost",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    system_database = FALSE,
    message = FALSE,
    catboost_paras_tunning = list(
      catboost_iterations = c(25L),
      catboost_depth = c(4L),
      catboost_learning_rate = c(0.05),
      catboost_l2_leaf_reg = c(3L)
    )
  )

  expect_true(length(res_bin) >= 1L)
  export_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  expect_true(any(grepl("^Trait_", basename(export_dirs))))
  latest_bin <- export_dirs[which.max(file.info(export_dirs)$mtime)]
  pred_bin_file <- find_runtime_predicted_file(latest_bin)
  expect_false(is.null(pred_bin_file))
  pred_bin <- utils::read.csv(pred_bin_file, stringsAsFactors = FALSE)
  expect_true(all(c("Predicted_class", "Prediction_confidence", "Classification_uncertainty") %in% names(pred_bin)))
  expect_false("Prediction_stability" %in% names(pred_bin))
  expect_exported_classification_contract(pred_bin, c("no", "yes"))

  mdat <- mk_runtime_data("multiclass", n = 36L, p = 8L)
  pheno_multi <- mdat$all_pheno
  pheno_multi$Trait[sample(seq_len(nrow(pheno_multi)), 9L)] <- NA

  res_multi <- PredictProR::model_execute(
    pheno_data = pheno_multi,
    geno_data = mdat$all_geno,
    response = "Trait",
    gen_name = "GID",
    response_family = "nominal",
    GS_model = "mlp",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE,
    system_database = FALSE,
    message = FALSE,
    dpl_paras_tunning = list(
      epochs = c(1L),
      batch_size = c(4L),
      mlp_learning_rate = c(1e-3)
    ),
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    validation_split = 0.2,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    epochs = 1L,
    batch_size = 4L,
    mlp_neurons_per_layer = as.integer(c(16L, 8L))
  )

  expect_true(length(res_multi) >= 1L)
  export_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  latest_multi <- export_dirs[which.max(file.info(export_dirs)$mtime)]
  pred_multi_file <- find_runtime_predicted_file(latest_multi)
  expect_false(is.null(pred_multi_file))
  pred_multi <- utils::read.csv(pred_multi_file, stringsAsFactors = FALSE)
  expect_true(all(c("Predicted_class", "Prediction_confidence", "Classification_uncertainty") %in% names(pred_multi)))
  expect_gte(sum(grepl("^(Prob|Probability)_", names(pred_multi))), 3L)
  expect_false("Prediction_stability" %in% names(pred_multi))
  expect_exported_classification_contract(pred_multi, c("A", "B", "C"))
})
