test_that("model_execute rejects family-incompatible eval metrics before fitting", {
  pheno <- data.frame(
    GID = paste0("g", 1:6),
    Yield = c(1.2, 2.3, 3.4, 4.5, 5.6, 6.7),
    stringsAsFactors = FALSE
  )

  geno <- matrix(
    c(
      0, 1, 2, 0,
      1, 1, 2, 0,
      2, 0, 1, 1,
      0, 2, 1, 2,
      1, 0, 0, 1,
      2, 1, 1, 0
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(pheno$GID, paste0("m", 1:4))
  )

  expect_error(
    PredictProR::model_execute(
      pheno_data = pheno,
      geno_data = geno,
      response = "Yield",
      gen_name = "GID",
      response_family = "gaussian",
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      GS_model_cv = "RandomForest",
      eval_metrics = "macro_f1",
      cross_validation_meth = "K-Folds",
      nfolds = 2L,
      replication = 1L,
      system_database = TRUE,
      parallel_mode = "sequential"
    ),
    regexp = paste(
      "Unsupported eval_metrics for response_family='gaussian': macro_f1",
      "Allowed metrics: accuracy, mean_squared_error, bias, root_mean_squared_error",
      sep = ".*"
    )
  )
})
