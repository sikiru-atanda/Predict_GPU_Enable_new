test_that("feature_score_predictors returns deterministic trait-specific metadata", {
  x <- data.frame(
    m3 = c(0, 1, 0, 1, 0, 1),
    m1 = c(1, 2, 3, 4, 5, 6),
    m2 = c(6, 5, 4, 3, 2, 1),
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    y1 = c(1, 2, 3, 4, NA, NA),
    y2 = c(NA, 6, 5, 4, 3, NA),
    stringsAsFactors = FALSE
  )

  meta1 <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = c("y1", "y2"),
    gen_name = "ID",
    scoring_model = "Ridge_Regression",
    seed = 11
  )
  meta2 <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = c("y1", "y2"),
    gen_name = "ID",
    scoring_model = "Ridge_Regression",
    seed = 11
  )

  expect_s3_class(meta1, "predictpror_feature_scores")
  expect_identical(meta1, meta2)
  expect_equal(sort(unique(meta1$trait)), c("y1", "y2"))
  expect_true(all(c(
    "trait", "predictor", "source_block", "rank", "raw_score",
    "normalized_score", "scoring_model", "seed", "n_train"
  ) %in% names(meta1)))
  trait_counts <- table(meta1$trait)
  expect_equal(as.integer(trait_counts), c(3L, 3L))
  expect_equal(names(trait_counts), c("y1", "y2"))
  expect_true(all(meta1$normalized_score >= 0 & meta1$normalized_score <= 1))
  expect_equal(unique(meta1$n_train[meta1$trait == "y1"]), 4)
  expect_equal(unique(meta1$n_train[meta1$trait == "y2"]), 4)
})

test_that("feature_select_top_k applies metadata with stable predictor order", {
  x <- data.frame(
    m1 = 1:4,
    m2 = 5:8,
    m3 = 9:12,
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  metadata <- data.frame(
    trait = "y",
    predictor = c("m3", "m1", "m2"),
    source_block = "merged",
    rank = c(1L, 2L, 3L),
    raw_score = c(3, 2, 1),
    normalized_score = c(1, 0.5, 0),
    scoring_model = "Ridge_Regression",
    seed = 1L,
    n_train = 4L,
    stringsAsFactors = FALSE
  )

  selected <- feature_select_top_k(x, metadata, trait = "y", k = 2)

  expect_equal(colnames(selected), c("m3", "m1"))
  expect_equal(nrow(selected), nrow(x))
})

test_that("gp_feature_apply_to_ml_dat_res applies the same top-k columns to train and test blocks", {
  x_train <- data.frame(
    m1 = 1:4,
    m2 = 5:8,
    m3 = 9:12,
    m4 = 13:16,
    check.names = FALSE
  )
  rownames(x_train) <- paste0("g", seq_len(nrow(x_train)))
  x_test <- data.frame(
    m1 = 21:22,
    m2 = 23:24,
    m3 = 25:26,
    m4 = 27:28,
    check.names = FALSE
  )
  rownames(x_test) <- c("g5", "g6")
  metadata <- data.frame(
    trait = "Y",
    predictor = c("m4", "m2", "m3", "m1"),
    source_block = "merged",
    rank = c(1L, 2L, 3L, 4L),
    raw_score = c(4, 3, 2, 1),
    normalized_score = c(1, 0.75, 0.5, 0.25),
    scoring_model = "Ridge_Regression",
    seed = 1L,
    n_train = 4L,
    stringsAsFactors = FALSE
  )
  ml_dat_res <- list(
    merged_data = list(merge_data = x_train),
    merged_data_test = x_test,
    omic_count = 1L
  )

  selected <- gp_feature_apply_to_ml_dat_res(
    ml_dat_res = ml_dat_res,
    feature_score_metadata = metadata,
    trait = "Y",
    k = 2L
  )

  expect_equal(colnames(selected$merged_data$merge_data), c("m4", "m2"))
  expect_equal(colnames(selected$merged_data_test), c("m4", "m2"))
  expect_equal(selected$feature_k, 2L)
})

test_that("true prediction ML route receives the selected top-k marker/omics columns", {
  x <- data.frame(
    m1 = c(1, 2, 3, 4, 5, 6),
    m2 = c(6, 5, 4, 3, 2, 1),
    m3 = c(0, 1, 0, 1, 0, 1),
    m4 = c(1, 1, 2, 2, 3, 3),
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Y = c(1, 2, 3, 4, NA, NA),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    trait = "Y",
    predictor = c("m3", "m1", "m4", "m2"),
    source_block = "merged",
    rank = c(1L, 2L, 3L, 4L),
    raw_score = c(4, 3, 2, 1),
    normalized_score = c(1, 0.75, 0.5, 0.25),
    scoring_model = "Ridge_Regression",
    seed = 7L,
    n_train = 4L,
    stringsAsFactors = FALSE
  )
  captured <- new.env(parent = emptyenv())
  captured$train_cols <- NULL
  captured$test_cols <- NULL

  local_mocked_bindings(
    AI_randomForest = function(pheno_object,
                               response,
                               geno_omic_object,
                               geno_omic_test_object = NULL,
                               gen_name,
                               ...) {
      captured$train_cols <- colnames(geno_omic_object)
      captured$test_cols <- colnames(geno_omic_test_object)
      ids <- c(rownames(geno_omic_object), rownames(geno_omic_test_object))
      list(
        model_parameters = data.frame(stat = "model", summary = "RandomForest", stringsAsFactors = FALSE),
        predicted_values = data.frame(
          ID = ids,
          Predicted_value = seq_along(ids),
          Train_Test_Label = c(rep("Train", nrow(geno_omic_object)), rep("Test", nrow(geno_omic_test_object))),
          stringsAsFactors = FALSE
        ),
        diagnostic_plots = NULL,
        variance_components = gp_empty_variance_components(),
        Variance_components = gp_empty_variance_components()
      )
    },
    summary_statistics_AI = function(...) list(metric = "mocked"),
    .package = "PredictProR"
  )

  res <- gp_run_best_model_task(list(
    response = "Y",
    GS_model = "RandomForest",
    pheno_clean = list(pheno_clean_data = pheno),
    gen_name = "ID",
    ml_dat_res = list(
      pheno_clean_data = pheno,
      merged_data = list(merge_data = x),
      omic_count = 1L
    ),
    feature_score_metadata = metadata,
    feature_k = 2L,
    response_family = "gaussian",
    cross_validation = FALSE,
    met_ml_dl = FALSE,
    heter_groups = NULL,
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    AI_valid_models = "RandomForest",
    rand_term_model_bayesian = NULL,
    friendly_name_lookup = NULL,
    msg = ""
  ))

  expect_equal(captured$train_cols, c("m3", "m1"))
  expect_equal(captured$test_cols, c("m3", "m1"))
  expect_equal(
    res$res_model_output$model_parameters$summary[
      match("n_predictors_used", res$res_model_output$model_parameters$stat)
    ],
    "2"
  )
})

test_that("true prediction ML route honors explicit feature_selected columns", {
  x <- data.frame(
    m1 = c(1, 2, 3, 4, 5, 6),
    m2 = c(6, 5, 4, 3, 2, 1),
    m3 = c(0, 1, 0, 1, 0, 1),
    m4 = c(1, 1, 2, 2, 3, 3),
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Y = c(1, 2, 3, 4, NA, NA),
    stringsAsFactors = FALSE
  )
  captured <- new.env(parent = emptyenv())
  captured$train_cols <- NULL
  captured$test_cols <- NULL

  local_mocked_bindings(
    AI_randomForest = function(pheno_object,
                               response,
                               geno_omic_object,
                               geno_omic_test_object = NULL,
                               gen_name,
                               ...) {
      captured$train_cols <- colnames(geno_omic_object)
      captured$test_cols <- colnames(geno_omic_test_object)
      ids <- c(rownames(geno_omic_object), rownames(geno_omic_test_object))
      list(
        model_parameters = data.frame(stat = "model", summary = "RandomForest", stringsAsFactors = FALSE),
        predicted_values = data.frame(
          ID = ids,
          Predicted_value = seq_along(ids),
          Train_Test_Label = c(rep("Train", nrow(geno_omic_object)), rep("Test", nrow(geno_omic_test_object))),
          stringsAsFactors = FALSE
        ),
        diagnostic_plots = NULL,
        variance_components = gp_empty_variance_components(),
        Variance_components = gp_empty_variance_components()
      )
    },
    summary_statistics_AI = function(...) list(metric = "mocked"),
    .package = "PredictProR"
  )

  res <- gp_run_best_model_task(list(
    response = "Y",
    GS_model = "RandomForest",
    pheno_clean = list(pheno_clean_data = pheno),
    gen_name = "ID",
    ml_dat_res = list(
      pheno_clean_data = pheno,
      merged_data = list(merge_data = x),
      omic_count = 1L
    ),
    feature_selected = c("m4", "m2"),
    feature_score_metadata = NULL,
    feature_k = NULL,
    response_family = "gaussian",
    cross_validation = FALSE,
    met_ml_dl = FALSE,
    heter_groups = NULL,
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    AI_valid_models = "RandomForest",
    rand_term_model_bayesian = NULL,
    friendly_name_lookup = NULL,
    msg = ""
  ))

  expect_equal(captured$train_cols, c("m4", "m2"))
  expect_equal(captured$test_cols, c("m4", "m2"))
  expect_equal(
    res$res_model_output$model_parameters$summary[
      match("n_predictors_used", res$res_model_output$model_parameters$stat)
    ],
    "2"
  )
})

test_that("feature_selected can target genotype and omics sources separately", {
  selected <- list(
    geno_data = c("m4", "m2"),
    omic1_data = c("o3", "o1")
  )

  expect_equal(
    gp_feature_selected_vector(selected, source = "geno_data"),
    c("m4", "m2")
  )
  expect_equal(
    gp_feature_selected_vector(selected, source = "omic1_data"),
    c("o3", "o1")
  )
})

test_that("feature_selected filters raw marker and omics blocks before GRM and kernel preprocessing", {
  ids <- paste0("g", seq_len(5))
  geno <- data.frame(
    m1 = c(0, 1, 2, 1, 0),
    m2 = c(2, 1, 0, 1, 2),
    m3 = c(1, 1, 1, 0, 0),
    m4 = c(0, 0, 1, 1, 2),
    check.names = FALSE
  )
  rownames(geno) <- ids
  omic1 <- data.frame(
    o1 = c(1, 2, 3, 4, 5),
    o2 = c(5, 4, 3, 2, 1),
    o3 = c(0, 1, 0, 1, 0),
    check.names = FALSE
  )
  rownames(omic1) <- ids
  pheno <- data.frame(
    ID = ids,
    Y = c(1, 2, 3, NA, NA),
    stringsAsFactors = FALSE
  )

  captured <- new.env(parent = emptyenv())
  captured$geno_cols <- NULL
  captured$omic_cols <- NULL

  local_mocked_bindings(
    process_geno_data = function(geno_data = NULL,
                                 train_geno_data = NULL,
                                 test_geno_data = NULL,
                                 ...) {
      captured$geno_cols <- colnames(geno_data)
      kernel_ids <- rownames(geno_data)
      kernel <- diag(length(kernel_ids))
      dimnames(kernel) <- list(kernel_ids, kernel_ids)
      list(
        gmatrix = kernel,
        geno_model_ready = geno_data
      )
    },
    process_omic_data = function(omic_data = NULL,
                                 train_omic_data = NULL,
                                 test_omic_data = NULL,
                                 ...) {
      if (is.null(omic_data) && is.null(train_omic_data) && is.null(test_omic_data)) {
        return(NULL)
      }
      captured$omic_cols <- colnames(omic_data)
      kernel_ids <- rownames(omic_data)
      kernel <- diag(length(kernel_ids))
      dimnames(kernel) <- list(kernel_ids, kernel_ids)
      list(
        kernel = kernel,
        omic_model_ready = omic_data
      )
    },
    grm_kernel_precheck = function(grm_kernel_data, ...) grm_kernel_data,
    pheno_geno_match = function(object_pheno, object_geno, ...) {
      list(geno_pheno_match_data = object_geno)
    },
    AI_process_ml_data_if_valid = function(...) list(),
    .package = "PredictProR"
  )

  gp_prepare_model_input_objects(list(
    pheno_clean = list(pheno_clean_data = pheno),
    geno_data = geno,
    train_geno_data = NULL,
    test_geno_data = NULL,
    omic1_data = omic1,
    train_omic1_data = NULL,
    test_omic1_data = NULL,
    omic2_data = NULL,
    train_omic2_data = NULL,
    test_omic2_data = NULL,
    omic3_data = NULL,
    train_omic3_data = NULL,
    test_omic3_data = NULL,
    feature_selected = list(
      geno_data = c("m4", "m2"),
      omic1_data = c("o3", "o1")
    ),
    response = "Y",
    response_family = "gaussian",
    gen_name = "ID",
    GS_model = "KRR",
    GS_model_cv = NULL,
    AI_valid_models = character(),
    feature_scoring = FALSE,
    feature_score_metadata = NULL,
    feature_k_grid = NULL,
    feature_k = NULL,
    met_ml_dl = FALSE,
    heter_groups = NULL,
    test_set = NULL,
    train_set = NULL,
    kernel_method = "Gaussian_kernel",
    gmatrix_method = "VanRaden",
    scale = FALSE,
    map_data = NULL,
    maf_threshold = NULL,
    het_threshold = NULL,
    ind_call_rate_threshold = NULL,
    snp_call_rate_threshold = NULL,
    impute = FALSE,
    imputation_method = "mean",
    impute_knn_k = 5L,
    ld_prunning_qc = FALSE,
    geno_data_process = NULL,
    qc_filtering = FALSE,
    message = FALSE,
    impute_omic = FALSE,
    na_threshold = 0.9,
    gmatrix = NULL,
    gkernel = NULL,
    kernel_list = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    pedigree_matrix = NULL
  ))

  expect_equal(captured$geno_cols, c("m4", "m2"))
  expect_equal(captured$omic_cols, c("o3", "o1"))
})

test_that("feature_selected filters hybrid parent genotype blocks before specialized routes", {
  ids <- paste0("p", seq_len(4))
  female <- data.frame(
    m1 = c(0, 1, 2, 1),
    m2 = c(2, 1, 0, 1),
    m3 = c(1, 0, 1, 0),
    check.names = FALSE
  )
  male <- data.frame(
    m1 = c(1, 1, 2, 2),
    m2 = c(0, 1, 0, 1),
    m3 = c(2, 2, 1, 1),
    check.names = FALSE
  )
  rownames(female) <- ids
  rownames(male) <- ids

  filtered <- gp_feature_apply_raw_selection_to_inputs(list(
    female_geno_data = female,
    male_geno_data = male,
    feature_selected = list(geno_data = c("m3", "m1")),
    gen_name = "ID"
  ))

  expect_equal(colnames(filtered$female_geno_data), c("m3", "m1"))
  expect_equal(colnames(filtered$male_geno_data), c("m3", "m1"))
})

test_that("RandomForest feature scoring uses Python feature importance", {
  captured <- NULL
  local_mocked_bindings(
    gp_py_ml_feature_importance = function(model_type,
                                           X_train,
                                           y_train,
                                           response_family = "gaussian",
                                           model_params = list(),
                                           class_levels = NULL) {
      captured <<- list(
        model_type = model_type,
        X_train = X_train,
        y_train = y_train,
        response_family = response_family,
        model_params = model_params
      )
      setNames(c(0.2, 0.8, 0.1), colnames(X_train))
    },
    .package = "PredictProR"
  )

  x <- data.frame(m1 = 1:5, m2 = 5:1, m3 = c(1, 1, 0, 0, 1), check.names = FALSE)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(ID = rownames(x), Y = c(1, NA, 3, 4, 5), stringsAsFactors = FALSE)

  meta <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = "Y",
    gen_name = "ID",
    scoring_model = "RandomForest",
    ntree = 20L,
    rf_n_jobs = 1L
  )

  expect_identical(captured$model_type, "randomforest")
  expect_false(anyNA(captured$y_train))
  expect_equal(nrow(captured$X_train), 4L)
  expect_equal(captured$model_params$n_jobs, 1L)
  expect_equal(meta$predictor[order(meta$rank)][[1L]], "m2")
})

test_that("models_execute_crossval evaluates separate feature k values", {
  skip_if_no_ml_python()
  set.seed(7)
  x <- as.data.frame(matrix(rnorm(72), nrow = 12), check.names = FALSE)
  names(x) <- paste0("m", seq_len(ncol(x)))
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Y = x$m1 * 0.8 - x$m2 * 0.4 + rnorm(nrow(x), sd = 0.05),
    stringsAsFactors = FALSE
  )
  metadata <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = "Y",
    gen_name = "ID",
    scoring_model = "Ridge_Regression",
    seed = 7
  )

  res <- models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Y",
      gen_name = "ID",
      test_size = 0.25,
      random_state = 7,
      replication = 1L,
      cross_validation_meth = "K-Folds",
      nfolds = 3L,
      eval_metrics = "mean_squared_error",
      response_family = "gaussian",
      GS_model_cv = "RandomForest",
      ml_dat_res = list(
        merged_data = list(merge_data = x),
        omic_count = 1
      ),
      ntree = 20L,
      rf_n_jobs = 1L,
      scaling = FALSE,
      centering = TRUE,
      feature_score_metadata = metadata,
      feature_k_grid = c(2L, 4L),
      feature_scoring_model = "Ridge_Regression",
      feature_scoring_cv = "fixed",
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_equal(length(res), 2)
  observed_k <- sort(vapply(res, function(x) unique(x$eval_metrics_reps$feature_k), integer(1)))
  expect_equal(observed_k, c(2L, 4L))
})

test_that("models_execute_crossval passes selected top-k columns into AI fold prediction", {
  set.seed(17)
  x <- as.data.frame(matrix(rnorm(48), nrow = 12), check.names = FALSE)
  names(x) <- paste0("m", seq_len(ncol(x)))
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Y = x$m1 - 0.25 * x$m4 + rnorm(nrow(x), sd = 0.01),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    trait = "Y",
    predictor = c("m4", "m1", "m2", "m3"),
    source_block = "merged",
    rank = c(1L, 2L, 3L, 4L),
    raw_score = c(4, 3, 2, 1),
    normalized_score = c(1, 0.75, 0.5, 0.25),
    scoring_model = "Ridge_Regression",
    seed = 17L,
    n_train = nrow(x),
    stringsAsFactors = FALSE
  )
  captured <- list()

  local_mocked_bindings(
    predict_with_model = function(model, y, omics_data = NULL, tst, additional_params = list(), ...) {
      captured[[length(captured) + 1L]] <<- colnames(omics_data)
      rep(mean(y[-tst], na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Y",
      gen_name = "ID",
      test_size = 0.25,
      random_state = 17,
      replication = 1L,
      cross_validation_meth = "K-Folds",
      nfolds = 3L,
      eval_metrics = "mean_squared_error",
      response_family = "gaussian",
      GS_model_cv = "RandomForest",
      ml_dat_res = list(
        merged_data = list(merge_data = x),
        omic_count = 1
      ),
      ntree = 20L,
      rf_n_jobs = 1L,
      scaling = FALSE,
      centering = TRUE,
      feature_score_metadata = metadata,
      feature_k_grid = 2L,
      feature_scoring_model = "Ridge_Regression",
      feature_scoring_cv = "fixed",
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_equal(length(res), 1)
  expect_true(length(captured) > 0)
  expect_true(all(vapply(captured, identical, logical(1), c("m4", "m1"))))
  expect_equal(unique(res[[1]]$eval_metrics_reps$feature_k), 2L)
})

test_that("fold_internal scoring excludes outer test responses and preserves provenance", {
  set.seed(101)
  x <- as.data.frame(matrix(rnorm(60), nrow = 12), check.names = FALSE)
  names(x) <- paste0("m", seq_len(ncol(x)))
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Y = 1.5 * x$m1 - x$m4 + rnorm(nrow(x), sd = 0.05),
    stringsAsFactors = FALSE
  )
  metadata <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = "Y",
    gen_name = "ID",
    scoring_model = "Ridge_Regression",
    seed = 101L
  )

  local_mocked_bindings(
    predict_with_model = function(model, y, omics_data = NULL, tst, additional_params = list(), ...) {
      rep(mean(as.numeric(y[-tst]), na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Y",
      gen_name = "ID",
      random_state = 101L,
      replication = 1L,
      cross_validation_meth = "K-Folds",
      nfolds = 3L,
      eval_metrics = "mean_squared_error",
      response_family = "gaussian",
      GS_model_cv = "RandomForest",
      ml_dat_res = list(merged_data = list(merge_data = x), omic_count = 1L),
      feature_score_metadata = metadata,
      feature_k_grid = 2L,
      feature_scoring_model = "Ridge_Regression",
      feature_scoring_cv = "fold_internal",
      feature_scoring_seed = 101L,
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_length(res, 1L)
  fold_meta <- res[[1]]$feature_selection_metadata
  expect_true(is.data.frame(fold_meta))
  expect_equal(sort(unique(fold_meta$fold)), 1:3)
  expect_true(all(fold_meta$selection_mode == "fold_internal"))
  expect_true(all(fold_meta$n_train == 8L))
  expect_equal(length(unique(fold_meta$training_id_hash)), 3L)
  expect_true(all(tapply(fold_meta$selected, fold_meta$fold, sum) == 2L))
  expect_equal(unique(res[[1]]$eval_metrics_reps$feature_scoring_cv), "fold_internal")
})

test_that("both feature selection modes remain separate and rank from fold_internal", {
  set.seed(102)
  x <- as.data.frame(matrix(rnorm(48), nrow = 12), check.names = FALSE)
  names(x) <- paste0("m", seq_len(ncol(x)))
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(ID = rownames(x), Y = x$m1 + rnorm(12, sd = 0.1))
  metadata <- feature_score_predictors(x, pheno, "Y", "ID", seed = 102L)

  local_mocked_bindings(
    predict_with_model = function(model, y, omics_data = NULL, tst, additional_params = list(), ...) {
      rep(mean(as.numeric(y[-tst]), na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Y", gen_name = "ID", random_state = 102L, replication = 1L,
      cross_validation_meth = "K-Folds", nfolds = 3L,
      eval_metrics = "mean_squared_error", response_family = "gaussian",
      GS_model_cv = "RandomForest",
      ml_dat_res = list(merged_data = list(merge_data = x), omic_count = 1L),
      feature_score_metadata = metadata, feature_k_grid = 2L,
      feature_scoring_model = "Ridge_Regression", feature_scoring_cv = "both",
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_length(res, 2L)
  expect_equal(
    sort(vapply(res, function(z) unique(z$eval_metrics_reps$feature_scoring_cv), character(1))),
    c("fixed", "fold_internal")
  )
})

test_that("trait-specific top-k sets do not collapse across responses", {
  x <- data.frame(
    t1_signal = seq_len(10),
    t2_signal = rep(c(-2, 2), 5),
    noise = c(2, -1, 0, 3, -2, 1, -3, 2, 0, 1),
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Trait1 = x$t1_signal,
    Trait2 = x$t2_signal,
    stringsAsFactors = FALSE
  )
  meta <- feature_score_predictors(
    x, pheno, c("Trait1", "Trait2"), "ID",
    scoring_model = "Ridge_Regression", seed = 2L
  )

  expect_equal(gp_feature_selected_metadata(meta, "Trait1", 1L)$predictor, "t1_signal")
  expect_equal(gp_feature_selected_metadata(meta, "Trait2", 1L)$predictor, "t2_signal")
})

test_that("fixed single-trait selection records the exact model columns", {
  x <- data.frame(m1 = 1:8, m2 = 8:1, m3 = rep(c(0, 1), 4))
  rownames(x) <- paste0("g", 1:8)
  pheno <- data.frame(ID = rownames(x), Y = x$m2 + 0.1 * x$m3)
  meta <- feature_score_predictors(x, pheno, "Y", "ID", seed = 8L)
  dat <- list(
    merged_data = list(merge_data = x),
    merged_data_test = x[1:2, , drop = FALSE]
  )

  selected <- gp_feature_apply_to_ml_dat_res(dat, meta, "Y", 2L)

  expect_equal(selected$feature_selected, colnames(selected$merged_data$merge_data))
  expect_equal(colnames(selected$merged_data_test), selected$feature_selected)
  expect_equal(sum(selected$feature_selection_metadata$selected), 2L)
  expect_true(all(selected$feature_selection_metadata$selection_mode == "fixed"))
})

test_that("non-gaussian coefficient scoring fails instead of depending on factor codes", {
  x <- data.frame(m1 = 1:8, m2 = 8:1)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    Class = factor(rep(c("A", "B"), 4)),
    stringsAsFactors = FALSE
  )
  expect_error(
    feature_score_predictors(
      x, pheno, "Class", "ID",
      scoring_model = "Ridge_Regression",
      response_family = "binary"
    ),
    "only defined for gaussian responses"
  )
})

test_that("kernel adapter rebuilds a relationship matrix from selected source columns", {
  ids <- paste0("g", 1:6)
  geno <- matrix(
    c(0, 1, 2, 0, 1, 2, 2, 1, 0, 2, 1, 0, 0, 0, 1, 1, 2, 2),
    nrow = 6,
    dimnames = list(ids, c("m1", "m2", "m3"))
  )
  captured <- NULL
  local_mocked_bindings(
    grm_calculation = function(geno_clean, method, ...) {
      captured <<- colnames(geno_clean)
      out <- tcrossprod(scale(geno_clean, center = TRUE, scale = FALSE))
      dimnames(out) <- list(rownames(geno_clean), rownames(geno_clean))
      out
    },
    .package = "PredictProR"
  )

  bank <- gp_feature_rebuild_true_kernel_bank(
    model = "KRR",
    trait = "Yield",
    selected_by_source = list(geno_data = c("m3", "m1")),
    source_matrices = list(geno_data = geno),
    gmatrix_method = "VanRaden"
  )

  expect_equal(captured, c("m3", "m1"))
  expect_equal(dim(bank$gmatrix_model_ready), c(6L, 6L))
})

test_that("joint multi-trait view keeps trait selections and uses their union", {
  x <- data.frame(
    a = seq_len(12),
    b = rep(c(-3, 3), 6),
    c = rep(c(0, 1, 1), 4),
    check.names = FALSE
  )
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  pheno <- data.frame(
    ID = rownames(x),
    T1 = x$a,
    T2 = x$b,
    stringsAsFactors = FALSE
  )
  meta <- feature_score_predictors(x, pheno, c("T1", "T2"), "ID", seed = 12L)
  view <- gp_feature_multitrait_view(
    predictor_data = x,
    pheno_data = pheno,
    response = c("T1", "T2"),
    gen_name = "ID",
    k = 1L,
    selection_mode = "fixed",
    feature_score_metadata = meta
  )

  expect_equal(view$selected_by_trait$T1, "a")
  expect_equal(view$selected_by_trait$T2, "b")
  expect_equal(colnames(view$predictor_data), c("a", "b"))
  expect_equal(view$joint_union_k, 2L)
  expect_true(all(c("joint_union_selected", "joint_union_k") %in% names(view$metadata)))
})

test_that("joint fold_internal rankings are fitted only on outer training rows", {
  set.seed(22)
  x <- as.data.frame(matrix(rnorm(80), nrow = 16), check.names = FALSE)
  names(x) <- paste0("m", 1:5)
  rownames(x) <- paste0("g", 1:16)
  pheno <- data.frame(
    ID = rownames(x),
    T1 = x$m1 + rnorm(16, sd = 0.01),
    T2 = x$m4 + rnorm(16, sd = 0.01),
    stringsAsFactors = FALSE
  )
  view <- gp_feature_multitrait_view(
    predictor_data = x,
    pheno_data = pheno,
    response = c("T1", "T2"),
    gen_name = "ID",
    k = 2L,
    selection_mode = "fold_internal",
    test_rows = 13:16,
    scoring_model = "Ridge_Regression",
    seed = 22L,
    replication = 1L,
    fold = 4L
  )

  expect_equal(unique(view$metadata$n_train), 12L)
  expect_equal(sort(unique(view$metadata$trait)), c("T1", "T2"))
  expect_true(all(view$metadata$selection_mode == "fold_internal"))
  expect_true(view$joint_union_k >= 2L)
  expect_true(view$joint_union_k <= 4L)
})

test_that("specialized routes require one k and one selection mode per run", {
  meta <- data.frame(trait = "Y", predictor = "m1", rank = 1L)
  expect_error(
    gp_feature_specialized_cv_config(meta, feature_k_grid = c(1L, 2L)),
    "one pre-specified"
  )
  expect_error(
    gp_feature_specialized_cv_config(meta, feature_k = 1L, feature_scoring_cv = "both"),
    "run both modes separately"
  )
  expect_equal(
    gp_feature_specialized_cv_config(meta, feature_k = 1L, feature_scoring_cv = "fold_internal")$mode,
    "fold_internal"
  )
})

test_that("multi-trait ML fold_internal route records trait-wise selections", {
  set.seed(31)
  x <- as.data.frame(matrix(rnorm(72), nrow = 12), check.names = FALSE)
  names(x) <- paste0("m", 1:6)
  rownames(x) <- paste0("g", 1:12)
  pheno <- data.frame(
    ID = rownames(x),
    T1 = x$m1 + rnorm(12, sd = 0.01),
    T2 = x$m5 + rnorm(12, sd = 0.01),
    stringsAsFactors = FALSE
  )
  meta <- feature_score_predictors(x, pheno, c("T1", "T2"), "ID", seed = 31L)

  # The CV fold fits run as one batch; predict the (scaled) training mean.
  local_mocked_bindings(
    gp_py_ml_fit_predict_batch = function(jobs) {
      lapply(jobs, function(job) {
        built <- if (is.function(job$build)) job$build() else job
        rep(mean(job$y_train), nrow(built$X_test))
      })
    },
    .package = "PredictProR"
  )

  out <- gp_multitrait_ml_gaussian_cv(
    model_type = "Ridge_Regression",
    pheno_object = pheno,
    response = c("T1", "T2"),
    geno_omic_object = x,
    gen_name = "ID",
    nfolds = 3L,
    replication = 1L,
    eval_metrics = "mean_squared_error",
    random_seed = 31L,
    feature_score_metadata = meta,
    feature_k = 2L,
    feature_scoring_cv = "fold_internal",
    feature_scoring_model = "Ridge_Regression"
  )

  feature_meta <- out$cv_results[[1]]$feature_selection_metadata
  expect_true(is.data.frame(feature_meta))
  expect_equal(sort(unique(feature_meta$trait)), c("T1", "T2"))
  expect_true(all(feature_meta$selection_mode == "fold_internal"))
  expect_true(all(out$cv_results[[1]]$eval_metrics_reps$feature_k == 2L))
  expect_true(all(out$cv_results[[1]]$eval_metrics_reps$feature_scoring_cv == "fold_internal"))
})

test_that("explicit selection data frames honor the selected flag", {
  selection <- data.frame(
    trait = c("T1", "T1", "T2"),
    source_block = c("geno_data", "geno_data", "geno_data"),
    predictor = c("m1", "m2", "m3"),
    selected = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  expect_identical(
    gp_feature_selected_vector(selection, trait = "T1", source = "geno_data"),
    "m1"
  )
})

test_that("source-keyed explicit selection excludes omitted raw sources", {
  ids <- paste0("g", 1:4)
  ctx <- list(
    feature_selected = list(geno_data = c("m2")),
    gen_name = "ID",
    geno_data = data.frame(
      ID = ids,
      m1 = 1:4,
      m2 = 4:1,
      check.names = FALSE
    ),
    omic1_data = data.frame(
      ID = ids,
      o1 = 11:14,
      o2 = 14:11,
      check.names = FALSE
    )
  )

  selected <- gp_feature_apply_raw_selection_to_inputs(ctx)

  expect_identical(colnames(selected$geno_data), c("ID", "m2"))
  expect_null(selected$omic1_data)
  expect_error(
    gp_feature_apply_raw_selection_to_inputs(utils::modifyList(
      ctx,
      list(feature_selected = list(geno_data = c("m2", "not_a_marker")))
    )),
    "not_a_marker"
  )
})

test_that("explicit source adapter preserves trait-specific sets and metadata", {
  ids <- paste0("g", 1:4)
  sources <- list(
    geno_data = matrix(1:12, nrow = 4, dimnames = list(ids, c("m1", "m2", "m3"))),
    omic1_data = matrix(21:28, nrow = 4, dimnames = list(ids, c("o1", "o2")))
  )
  selected <- list(
    T1 = list(geno_data = c("m1"), omic1_data = c("o2")),
    T2 = list(geno_data = c("m3"))
  )

  t1 <- gp_feature_explicit_selected_by_source(selected, "T1", sources)
  t2 <- gp_feature_explicit_selected_by_source(selected, "T2", sources)
  meta <- gp_feature_explicit_metadata(t1, "T1", replication = 2L, prediction_model = "BayesA")

  expect_identical(t1, list(geno_data = "m1", omic1_data = "o2"))
  expect_identical(t2, list(geno_data = "m3"))
  source_first <- list(
    geno_data = list(T1 = "m1", T2 = "m3"),
    omic1_data = list(T1 = "o2")
  )
  expect_identical(
    gp_feature_explicit_selected_by_source(source_first, "T1", sources),
    t1
  )
  expect_true(all(meta$selection_mode == "explicit"))
  expect_true(all(meta$scoring_model == "user_supplied"))
  expect_true(all(meta$selected))
  expect_true(all(meta$feature_k == 2L))
})

test_that("Bayesian true prediction receives the exact explicit set for each trait", {
  ids <- paste0("g", 1:4)
  pheno <- data.frame(
    GID = ids,
    T1 = c(1, 2, 3, NA),
    T2 = c(4, 3, NA, 1),
    stringsAsFactors = FALSE
  )
  geno <- matrix(1:12, nrow = 4, dimnames = list(ids, c("m1", "m2", "m3")))
  seen <- list()
  local_mocked_bindings(
    bayes_finalize_A_B_C_BL_BRR = function(pheno_data = NULL,
                                            response = NULL,
                                            geno_data = NULL,
                                            ...) {
      seen[[response]] <<- colnames(geno_data)
      label <- ifelse(is.na(pheno_data[[response]]), "Test", "Train")
      list(
        bayes_result = list(Predicted_value = data.frame(
          GID = pheno_data$GID,
          Predicted_value = seq_len(nrow(pheno_data)),
          Train_Test_Label = label,
          stringsAsFactors = FALSE
        )),
        bayes_model = list(model = list(
          y = pheno_data[[response]], yHat = rep(0, nrow(pheno_data)), varE = 1,
          ETA = list()
        ))
      )
    },
    summary_statistics_bayes = function(...) {
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )
  base_ctx <- list(
    GS_model = "BayesA",
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    AI_valid_models = character(),
    pheno_clean = list(pheno_clean_data = pheno),
    ml_dat_res = list(),
    gen_name = "GID",
    msg = "",
    fixed = NULL,
    random = ~ GID,
    weights = NULL,
    fixed_term_model_bayesian = NULL,
    rand_term_model_bayesian = NULL,
    geno_model_ready = geno,
    omic1_model_ready = NULL,
    omic2_model_ready = NULL,
    omic3_model_ready = NULL,
    feature_selected = list(T1 = "m1", T2 = c("m2", "m3")),
    feature_score_metadata = NULL,
    feature_k = NULL,
    nIter = 10L,
    burnIn = 2L,
    thin = 1L,
    scaling = FALSE,
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.7,
    low_reliability_thres = 0.4,
    n_components = 2L,
    threshold = 10,
    eval_metrics = "mean_squared_error",
    heter_groups = NULL,
    system_database = FALSE,
    friendly_name_lookup = c(BayesA = "BayesA"),
    cross_validation = FALSE
  )

  t1 <- gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "T1")))
  t2 <- gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "T2")))

  expect_identical(seen$T1, "m1")
  expect_identical(seen$T2, c("m2", "m3"))
  expect_identical(t1$res_model_output$feature_selection_metadata$predictor, "m1")
  expect_identical(t2$res_model_output$feature_selection_metadata$predictor, c("m2", "m3"))
})

test_that("explicit CV selection reports user-supplied provenance", {
  set.seed(77)
  x <- as.data.frame(matrix(rnorm(36), nrow = 12), check.names = FALSE)
  names(x) <- c("m1", "m2", "m3")
  rownames(x) <- paste0("g", 1:12)
  pheno <- data.frame(ID = rownames(x), Y = x$m1 + rnorm(12, sd = 0.1))
  seen <- list()
  local_mocked_bindings(
    predict_with_model = function(model, y, omics_data = NULL, tst, ...) {
      seen[[length(seen) + 1L]] <<- colnames(omics_data)
      rep(mean(y[-tst]), length(tst))
    },
    .package = "PredictProR"
  )

  out <- models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Y", gen_name = "ID", random_state = 77L, replication = 1L,
      cross_validation_meth = "K-Folds", nfolds = 3L,
      eval_metrics = "mean_squared_error", response_family = "gaussian",
      GS_model_cv = "RandomForest",
      ml_dat_res = list(merged_data = list(merge_data = x), omic_count = 1L),
      selected_raw = list(Y = c("m3", "m1")),
      feature_source_matrices = list(geno_data = as.matrix(x)),
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_true(all(vapply(seen, identical, logical(1L), c("m3", "m1"))))
  expect_true(all(out[[1]]$feature_selection_metadata$selection_mode == "explicit"))
  expect_true(all(out[[1]]$feature_selection_metadata$scoring_model == "user_supplied"))
  expect_identical(unique(out[[1]]$eval_metrics_reps$feature_scoring_cv), "explicit")
})

test_that("public model standardization preserves feature-selection provenance", {
  meta <- gp_feature_explicit_metadata(
    selected_by_source = list(geno_data = c("m1", "m3"), omic1_data = "o2"),
    trait = "Yield",
    prediction_model = "BayesA"
  )
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 1.2),
    Train_Test_Label = c("Train", "Test"),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    lower_bound = c(0.9, 0.8),
    upper_bound = c(1.3, 1.6),
    Reliability = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )

  public <- gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      variance_components = gp_empty_variance_components(),
      feature_selection_metadata = meta
    ),
    gen_name = "GID"
  )

  expect_identical(public$feature_selection_metadata, meta)
})

test_that("CV k grids always include the all-predictors candidate on standard routes", {
  expect_identical(PredictProR:::gp_feature_k_grid(c(25L, 100L), 1800L), c(25L, 100L, 1800L))
  expect_identical(PredictProR:::gp_feature_k_grid(c(25L, "all"), 1800L), c(25L, 1800L))
  expect_identical(PredictProR:::gp_feature_k_grid(c(50L, 5000L), 1800L), c(50L, 1800L))
  # specialized hybrid / joint multi-trait routes keep exactly the k supplied
  expect_identical(PredictProR:::gp_feature_k_grid(25L, 1800L, include_all = FALSE), 25L)
})

test_that("BayesB feature scoring keeps BGLR files out of the working directory", {
  testthat::skip_if_not_installed("BGLR")
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  set.seed(1)
  x <- matrix(rnorm(40 * 15), 40, 15, dimnames = list(NULL, paste0("m", 1:15)))
  y <- x[, 1] + rnorm(40, sd = 0.5)
  s1 <- PredictProR:::gp_feature_score_bayesb(x, y, seed = 7L, nIter = 200L, burnIn = 50L)
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0L)
  s2 <- PredictProR:::gp_feature_score_bayesb(x, y, seed = 7L, nIter = 200L, burnIn = 50L)
  expect_identical(s1, s2)
  expect_identical(which.max(s1), 1L)   # scores follow the column order
})