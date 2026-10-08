test_that("cross-validation uses trait-specific test indices under unbalanced traits", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(10, NA, 30),
    Height = c(NA, 20, 30),
    stringsAsFactors = FALSE
  )

  omics <- matrix(
    1:9,
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("m1", "m2", "m3"))
  )

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, omics_data = NULL, tst = NULL, additional_params = NULL) {
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
      GS_model_cv = "Xgboost",
      ml_dat_res = list(
        merged_data = list(merge_data = omics),
        omic_count = 1L
      )
    ),
    verbose = FALSE
  )

  expect_length(res, 2)
  yield_res <- res[[which(vapply(res, function(x) identical(x$trait, "Yield"), logical(1)))]]
  height_res <- res[[which(vapply(res, function(x) identical(x$trait, "Height"), logical(1)))]]

  expect_identical(which(yield_res$ypred_cv_Reps_all$cv_role == "test"), 1L)
  expect_identical(which(height_res$ypred_cv_Reps_all$cv_role == "test"), 2L)
})

test_that("deep-learning CV dispatcher forwards omics_data into deep_learning_model", {
  seen <- list()

  local_mocked_bindings(
    deep_learning_model = function(y = NULL, omics_data = NULL, tst = NULL, crossval = FALSE, ...) {
      seen <<- list(
        y = y,
        omics_data = omics_data,
        tst = tst,
        crossval = crossval
      )
      rep(0, length(tst))
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, 3)
  omics <- matrix(1:9, nrow = 3)
  tst <- 1L

  res <- PredictProR:::predict_with_model(
    model = "mlp",
    y = y,
    omics_data = omics,
    tst = tst,
    additional_params = list(
      scaling = FALSE,
      centering = FALSE,
      crossval = TRUE,
      omic_count = 1L
    )
  )

  expect_identical(res, 0)
  expect_identical(seen$omics_data, omics)
  expect_identical(seen$tst, tst)
  expect_identical(seen$crossval, TRUE)
})

test_that("Xgboost CV dispatcher forwards thread and seed controls", {
  seen <- list()

  local_mocked_bindings(
    AI_xgboost_cv = function(y = NULL, omics = NULL, tst = NULL, xgb_nthread = NULL, random_state = NULL, ...) {
      seen <<- list(
        y = y,
        omics = omics,
        tst = tst,
        xgb_nthread = xgb_nthread,
        random_state = random_state
      )
      rep(0, length(tst))
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, 3)
  omics <- matrix(1:9, nrow = 3)
  tst <- 2L

  res <- PredictProR:::predict_with_model(
    model = "Xgboost",
    y = y,
    omics_data = omics,
    tst = tst,
    additional_params = list(
      scaling = FALSE,
      centering = FALSE,
      omic_count = 1L,
      xgb_nthread = 4L,
      random_state = 99L
    )
  )

  expect_identical(res, 0)
  expect_identical(seen$omics, omics)
  expect_identical(seen$tst, tst)
  expect_identical(seen$xgb_nthread, 4L)
  expect_identical(seen$random_state, 99L)
})

test_that("CatBoost and LightGBM CV dispatcher forward model-specific controls", {
  seen_cat <- list()
  seen_lgb <- list()

  local_mocked_bindings(
    AI_catboost_cv = function(y = NULL, omics = NULL, tst = NULL, catboost_thread_count = NULL, ...) {
      seen_cat <<- list(y = y, omics = omics, tst = tst, catboost_thread_count = catboost_thread_count)
      rep(0, length(tst))
    },
    AI_lightgbm_cv = function(y = NULL, omics = NULL, tst = NULL, lightgbm_nthread = NULL, ...) {
      seen_lgb <<- list(y = y, omics = omics, tst = tst, lightgbm_nthread = lightgbm_nthread)
      rep(0, length(tst))
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, 3)
  omics <- matrix(1:9, nrow = 3)
  tst <- 2L

  expect_identical(
    PredictProR:::predict_with_model(
      model = "CatBoost",
      y = y,
      omics_data = omics,
      tst = tst,
      additional_params = list(
        scaling = FALSE,
        centering = FALSE,
        omic_count = 1L,
        catboost_thread_count = 3L
      )
    ),
    0
  )
  expect_identical(
    PredictProR:::predict_with_model(
      model = "LightGBM",
      y = y,
      omics_data = omics,
      tst = tst,
      additional_params = list(
        scaling = FALSE,
        centering = FALSE,
        omic_count = 1L,
        lightgbm_nthread = 4L
      )
    ),
    0
  )

  expect_identical(seen_cat$omics, omics)
  expect_identical(seen_cat$tst, tst)
  expect_identical(seen_cat$catboost_thread_count, 3L)
  expect_identical(seen_lgb$omics, omics)
  expect_identical(seen_lgb$tst, tst)
  expect_identical(seen_lgb$lightgbm_nthread, 4L)
})

test_that("random forest bootstrap uses resampled indices", {
  seen_nrow <- NULL

  local_mocked_bindings(
    gp_py_ml_bootstrap_predict = function(data_label_geno, indices, test_geno, model_type, ...) {
      seen_nrow <<- length(indices)
      rep(0, nrow(as.data.frame(test_geno)))
    },
    .package = "PredictProR"
  )

  data_label_geno <- data.frame(
    y = c(1, 2, 3),
    m1 = c(10, 20, 30),
    m2 = c(11, 21, 31),
    row.names = c("g1", "g2", "g3")
  )

  out <- PredictProR:::train_predict_randomForest(
    data_label_geno = data_label_geno,
    indices = c(1L, 3L),
    test_geno = as.matrix(data_label_geno[c(1L, 3L), -1, drop = FALSE]),
    ntree = 10L
  )

  expect_identical(seen_nrow, 2L)
  expect_identical(unname(out), c(0, 0))
})

test_that("xgboost bootstrap helper uses the python bridge with resampled rows", {
  seen_nrow <- NULL

  local_mocked_bindings(
    gp_py_ml_bootstrap_predict = function(data_label_geno, indices, test_geno, model_type, ...) {
      seen_nrow <<- length(indices)
      rep(0, nrow(test_geno))
    },
    .package = "PredictProR"
  )

  data_label_geno <- data.frame(
    y = c(1, 2, 3),
    m1 = c(10, 20, 30),
    m2 = c(11, 21, 31),
    row.names = c("g1", "g2", "g3")
  )

  out <- PredictProR:::train_predict_xgboost(
    data_label_geno = data_label_geno,
    indices = c(1L, 3L),
    test_geno = as.matrix(data_label_geno[c(1L, 3L), -1, drop = FALSE]),
    params = list(objective = "reg:squarederror"),
    nrounds = 5L
  )

  expect_identical(seen_nrow, 2L)
  expect_identical(unname(out), c(0, 0))
})

test_that("public single-environment ML wrappers use the Python tabular path", {
  seen <- list()

  local_mocked_bindings(
    gp_python_tabular_model = function(pheno_object = NULL,
                                       geno_omic_object = NULL,
                                       geno_omic_test_object = NULL,
                                       response = NULL,
                                       gen_name = NULL,
                                       response_family = "gaussian",
                                       model_type,
                                       model_label,
                                       model_params = list(),
                                       para_tunning = FALSE,
                                       tune_param_grid = NULL,
                                       tune_folds = 5L,
                                       ...) {
      seen[[model_label]] <<- list(
        model_type = model_type,
        model_params = model_params,
        para_tunning = para_tunning,
        tune_param_grid = tune_param_grid,
        tune_folds = tune_folds
      )
      list(model_parameters = data.frame(model = model_label))
    },
    .package = "PredictProR"
  )

  pheno <- data.frame(GID = paste0("g", 1:4), Yield = c(1, 2, 3, 4))
  geno <- matrix(seq_len(16), nrow = 4)
  rownames(geno) <- pheno$GID

  PredictProR:::AI_Xgb(
    pheno_object = pheno, geno_omic_object = geno, response = "Yield",
    gen_name = "GID", para_tunning = FALSE, n_bootstrap = 2L
  )
  PredictProR:::AI_randomForest(
    pheno_object = pheno, geno_omic_object = geno, response = "Yield",
    gen_name = "GID", para_tunning = FALSE, n_bootstrap = 2L
  )
  PredictProR:::AI_svm(
    pheno_object = pheno, geno_omic_object = geno, response = "Yield",
    gen_name = "GID", para_tunning = FALSE, n_bootstrap = 2L
  )
  PredictProR:::AI_knn(
    pheno_object = pheno, geno_omic_object = geno, response = "Yield",
    gen_name = "GID", para_tunning = FALSE, n_bootstrap = 2L
  )

  expect_setequal(names(seen), c("Xgboost", "RandomForest", "SupportVectorMachine", "K-NearestNeighbors"))
  expect_identical(seen$Xgboost$model_type, "xgboost")
  expect_identical(seen$RandomForest$model_type, "randomforest")
  expect_identical(seen$SupportVectorMachine$model_type, "svm")
  expect_identical(seen$`K-NearestNeighbors`$model_type, "knn")
  expect_false(any(vapply(seen, `[[`, logical(1), "para_tunning")))
})

test_that("multi-trait RandomForest ML uses the Python bridge", {
  seen_models <- character()

  local_mocked_bindings(
    gp_configure_python_runtime = function(...) invisible(NULL),
    gp_detect_python = function(...) NULL,
    gp_py_ml_fit_predict = function(model_type,
                                    X_train,
                                    y_train,
                                    X_test,
                                    response_family,
                                    model_params,
                                    prefer_gpu,
                                    class_levels = NULL) {
      seen_models <<- c(seen_models, model_type)
      rep(0, nrow(X_test))
    },
    .package = "PredictProR"
  )

  X_train <- matrix(seq_len(20), nrow = 5)
  X_test <- matrix(seq_len(8), nrow = 2)
  y_train <- c(1, 2, 3, 4, 5)

  out <- PredictProR:::gp_multitrait_ml_fit_trait(
    model_type = "RandomForest",
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    ntree = 10L
  )

  # one fit predicts the training and the target rows
  expect_identical(seen_models, "randomforest")
  expect_length(out$train_pred, nrow(X_train))
  expect_length(out$test_pred, nrow(X_test))
})

test_that("true-prediction ML runtime passes trait-specific test objects", {
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
    1:9,
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("m1", "m2", "m3"))
  )

  ml_dat_res <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = c("Yield", "Height"),
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )

  seen_test_ids <- list()

  local_mocked_bindings(
    AI_Xgb = function(pheno_object = NULL, geno_omic_object = NULL, geno_omic_test_object = NULL, response = NULL, ...) {
      seen_test_ids[[as.character(response)]] <<- if (is.null(geno_omic_test_object)) character() else rownames(geno_omic_test_object)
      ids <- c(rownames(geno_omic_object), if (!is.null(geno_omic_test_object)) rownames(geno_omic_test_object))
      list(
        model_parameters = data.frame(stat = "dummy", summary = 1),
        predicted_values = data.frame(
          name = ids,
          Predicted_value = seq_along(ids),
          Train_Test_Label = c(rep("Train", nrow(geno_omic_object)), rep("Test", length(ids) - nrow(geno_omic_object)))
        )
      )
    },
    summary_statistics_AI = function(...) list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
    .package = "PredictProR"
  )

  base_ctx <- list(
    GS_model = "Xgboost",
    AI_valid_models = c("Xgboost"),
    pheno_clean = pheno_clean,
    ml_dat_res = ml_dat_res,
    gen_name = "GID",
    msg = "",
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    para_tunning = FALSE,
    AI_cv_nfolds = 2L,
    xgb_paras_tunning = NULL,
    resample_method_tune = "cv",
    number_of_fold_tune = 2L,
    learning_rate = 0.1,
    xgb_gamma = 0,
    xgb_lambda = 1,
    xgb_alpha = 0,
    max_depth = 2L,
    subsample = 1,
    xgb_booster = "gbtree",
    colsample_bytree = 1,
    iteration = 1L,
    xgb_rate_drop = 0,
    xgb_skip_drop = 0,
    xgb_objective = "reg:squarederror",
    xgb_sample_type = "uniform",
    xgb_normalize_type = "tree",
    early_stop_for_iteration_xgb = FALSE,
    N_feature_impo = 0L,
    CI_width_thresholds = c(0.33, 0.66),
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5,
    n_components = 2L,
    threshold = 10,
    interval_width_low_threshold = NULL,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    n_bootstrap = 2L,
    eval_metrics = "mean_squared_error",
    friendly_name_lookup = c(Xgboost = "Xgboost"),
    cross_validation = FALSE,
    optimizer_name = "adam",
    use_amp = FALSE,
    max_grad_norm = 1,
    auto_class_weights = FALSE,
    cnn_neurons_per_layer = c(8L),
    cnn_kernel_size = 3L,
    cnn_dense_layers = c(8L),
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel = 2L,
    cnn_pool_stride = 2L,
    cnn_pool_padding = 0L,
    cnn_learning_rate = 1e-3,
    cnn_separable = FALSE,
    cnn_dilations = c(1L),
    cnn_use_se = FALSE,
    cnn_norm_type = "group",
    cnn_pool_type = "conv",
    cnn_use_global_pool = FALSE,
    resnet_neurons_per_block = c(8L),
    resnet_blocks = 1L,
    resnet_learning_rate = 1e-3,
    ft_d_model = 8L,
    ft_heads = 1L,
    ft_layers = 1L,
    ft_ff_mult = 2L,
    ft_dropout = 0.1,
    ft_token_dropout = 0,
    ft_use_cls = TRUE,
    saint_d_model = 8L,
    saint_heads = 1L,
    saint_layers = 1L,
    saint_ff_mult = 2L,
    saint_dropout = 0.1,
    saint_token_dropout = 0,
    saint_use_cls = TRUE,
    use_grouping = FALSE,
    group_trigger = 2048L,
    group_method = "auto",
    init_group_size = 64L,
    max_tokens = 1024L,
    kmeans_batch = 4096L,
    kmeans_iter = 10L,
    tabnet_steps = 2L,
    tabnet_feature_dim = 8L,
    tabnet_output_dim = 8L,
    tabnet_gamma = 1.5,
    tabnet_lambda_sparse = 1e-4,
    node_trees = 2L,
    node_depth = 2L,
    deepfm_k = 4L,
    deepfm_hidden = c(8L),
    dcn_layers = 1L,
    dcn_hidden = c(8L),
    nam_hidden = c(8L),
    nam_activation = "relu",
    nam_add_linear = TRUE,
    nam_l1 = 1e-4,
    moe_n_experts = 2L,
    moe_expert_hidden = c(8L),
    moe_gate_hidden = 8L,
    moe_temperature = 1,
    moe_sparse_topk = NA,
    moe_entropy_reg = 0,
    gp_use_variational = TRUE,
    gp_num_inducing = 8L,
    gp_feature_dim = 8L,
    gp_kernel = "rbf",
    gp_ard = TRUE,
    gp_lr_mult = 0.5,
    rff_features = 16L,
    rff_lengthscale = 1,
    rff_deep_hidden = c(8L),
    model_type = "mlp",
    epochs = 1L,
    batch_size = 4L,
    dropout = 0.1,
    l2_weight_decay = 1e-4,
    l2_regularizer_dp = 1e-4,
    dropout_rate = 0.1,
    batch_norm = FALSE,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    device = NULL,
    mlp_neurons_per_layer = c(8L),
    mlp_learning_rate = 1e-3,
    final_attention = FALSE,
    attention_across_multiple_layers = FALSE,
    heteroscedastic = FALSE,
    canonical_names = character()
  )

  PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Yield")))
  PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Height")))

  expect_identical(unname(seen_test_ids[["Yield"]]), "g2")
  expect_identical(unname(seen_test_ids[["Height"]]), "g1")
})

test_that("best-model runtime attaches connectivity outputs to returned results", {
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

  local_mocked_bindings(
    AI_Xgb = function(pheno_object = NULL, geno_omic_object = NULL, geno_omic_test_object = NULL, response = NULL, ...) {
      list(
        model_parameters = data.frame(stat = "dummy", summary = 1),
        predicted_values = data.frame(
          name = rownames(geno_omic_object),
          Predicted_value = seq_len(nrow(geno_omic_object)),
          Train_Test_Label = rep("Train", nrow(geno_omic_object))
        )
      )
    },
    summary_statistics_AI = function(...) list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
    .package = "PredictProR"
  )

  ml_dat_res <- list(
    merged_data = list(
      merge_data = matrix(
        1:15,
        nrow = 5,
        dimnames = list(c("g1", "g2", "g3", "g1_rep", "g2_rep"), c("m1", "m2", "m3"))
      )
    ),
    omic_count = 1L
  )

  ctx <- list(
    response = "Yield",
    GS_model = "Xgboost",
    AI_valid_models = "Xgboost",
    bayes_valid_models = character(),
    bayes_gblup_valid_models = character(),
    canonical_names = character(),
    pheno_clean = list(pheno_clean_data = pheno_checked),
    ml_dat_res = ml_dat_res,
    gen_name = "GID",
    msg = "",
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    para_tunning = FALSE,
    AI_cv_nfolds = 2L,
    xgb_paras_tunning = NULL,
    resample_method_tune = "cv",
    number_of_fold_tune = 2L,
    learning_rate = 0.1,
    xgb_gamma = 0,
    xgb_lambda = 1,
    xgb_alpha = 0,
    max_depth = 2L,
    subsample = 1,
    xgb_booster = "gbtree",
    colsample_bytree = 1,
    iteration = 1L,
    xgb_rate_drop = 0,
    xgb_skip_drop = 0,
    xgb_objective = "reg:squarederror",
    xgb_sample_type = "uniform",
    xgb_normalize_type = "tree",
    early_stop_for_iteration_xgb = FALSE,
    N_feature_impo = 0L,
    CI_width_thresholds = c(0.33, 0.66),
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5,
    n_components = 2L,
    threshold = 10,
    interval_width_low_threshold = NULL,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    n_bootstrap = 2L,
    eval_metrics = "mean_squared_error",
    friendly_name_lookup = c(Xgboost = "Xgboost"),
    cross_validation = FALSE,
    optimizer_name = "adam",
    use_amp = FALSE,
    max_grad_norm = 1,
    auto_class_weights = FALSE,
    cnn_neurons_per_layer = c(8L),
    cnn_kernel_size = 3L,
    cnn_dense_layers = c(8L),
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel = 2L,
    cnn_pool_stride = 2L,
    cnn_pool_padding = 0L,
    cnn_learning_rate = 1e-3,
    cnn_separable = FALSE,
    cnn_dilations = c(1L),
    cnn_use_se = FALSE,
    cnn_norm_type = "group",
    cnn_pool_type = "conv",
    cnn_use_global_pool = FALSE,
    resnet_neurons_per_block = c(8L),
    resnet_blocks = 1L,
    resnet_learning_rate = 1e-3,
    ft_d_model = 8L,
    ft_heads = 1L,
    ft_layers = 1L,
    ft_ff_mult = 2L,
    ft_dropout = 0.1,
    ft_token_dropout = 0,
    ft_use_cls = TRUE,
    saint_d_model = 8L,
    saint_heads = 1L,
    saint_layers = 1L,
    saint_ff_mult = 2L,
    saint_dropout = 0.1,
    saint_token_dropout = 0,
    saint_use_cls = TRUE,
    use_grouping = FALSE,
    group_trigger = 2048L,
    group_method = "auto",
    init_group_size = 64L,
    max_tokens = 1024L,
    kmeans_batch = 4096L,
    kmeans_iter = 10L,
    tabnet_steps = 2L,
    tabnet_feature_dim = 8L,
    tabnet_output_dim = 8L,
    tabnet_gamma = 1.5,
    tabnet_lambda_sparse = 1e-4,
    node_trees = 2L,
    node_depth = 2L,
    deepfm_k = 4L,
    deepfm_hidden = c(8L),
    dcn_layers = 1L,
    dcn_hidden = c(8L),
    nam_hidden = c(8L),
    nam_activation = "relu",
    nam_add_linear = TRUE,
    nam_l1 = 1e-4,
    moe_n_experts = 2L,
    moe_expert_hidden = c(8L),
    moe_gate_hidden = 8L,
    moe_temperature = 1,
    moe_sparse_topk = NA,
    moe_entropy_reg = 0,
    gp_use_variational = TRUE,
    gp_num_inducing = 8L,
    gp_feature_dim = 8L,
    gp_kernel = "rbf",
    gp_ard = TRUE,
    gp_lr_mult = 0.5,
    rff_features = 16L,
    rff_lengthscale = 1,
    rff_deep_hidden = c(8L),
    model_type = "mlp",
    epochs = 1L,
    batch_size = 4L,
    dropout = 0.1,
    l2_weight_decay = 1e-4,
    l2_regularizer_dp = 1e-4,
    dropout_rate = 0.1,
    batch_norm = FALSE,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    device = NULL,
    mlp_neurons_per_layer = c(8L),
    mlp_learning_rate = 1e-3,
    final_attention = FALSE,
    attention_across_multiple_layers = FALSE,
    heteroscedastic = FALSE
  )

  res <- PredictProR:::gp_run_best_model_task(ctx)

  expect_true("Connectivity_summary" %in% names(res$res_model_output))
  expect_true("Connectivity_genotype_count_by_env" %in% names(res$res_model_output))
  expect_true("Connectivity_summary" %in% names(res$res_summary_stat))
  expect_true("Connectivity_genotype_count_by_env" %in% names(res$res_summary_stat))
})

make_ml_runtime_ctx <- function(pheno_clean, ml_dat_res, gs_model, canonical_names = character()) {
  list(
    response = "Yield",
    GS_model = gs_model,
    AI_valid_models = c("Xgboost", "RandomForest", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors", "Lasso", "Ridge_Regression", canonical_names),
    bayes_valid_models = character(),
    bayes_gblup_valid_models = character(),
    canonical_names = canonical_names,
    pheno_clean = pheno_clean,
    ml_dat_res = ml_dat_res,
    gen_name = "GID",
    msg = "",
    scaling = FALSE,
    centering = FALSE,
    message = FALSE,
    para_tunning = FALSE,
    AI_cv_nfolds = 2L,
    xgb_paras_tunning = NULL,
    rr_paras_tunning = NULL,
    lasso_paras_tunning = NULL,
    rf_paras_tunning = NULL,
    pls_paras_tunning = NULL,
    svm_paras_tunning = NULL,
    knn_paras_tunning = NULL,
    dpl_paras_tunning = NULL,
    resample_method_tune = "cv",
    number_of_fold_tune = 2L,
    learning_rate = 0.1,
    xgb_gamma = 0,
    xgb_lambda = 1,
    xgb_alpha = 0,
    max_depth = 2L,
    subsample = 1,
    xgb_booster = "gbtree",
    colsample_bytree = 1,
    iteration = 1L,
    xgb_rate_drop = 0,
    xgb_skip_drop = 0,
    xgb_objective = "reg:squarederror",
    xgb_sample_type = "uniform",
    xgb_normalize_type = "tree",
    early_stop_for_iteration_xgb = FALSE,
    N_feature_impo = 0L,
    ntree = 10L,
    ncomp = 2L,
    lambda_rr = 1,
    C_value = 1,
    degree_value = 3,
    scale_value = 1,
    offset_value = 0,
    k = 3L,
    CI_width_thresholds = c(0.33, 0.66),
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5,
    n_components = 2L,
    threshold = 10,
    interval_width_low_threshold = NULL,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    n_bootstrap = 2L,
    eval_metrics = "mean_squared_error",
    friendly_name_lookup = c(
      Xgboost = "Xgboost",
      RandomForest = "RandomForest",
      PartialLeastSquare = "PartialLeastSquare",
      SupportVectorMachine = "SupportVectorMachine",
      `K-NearestNeighbors` = "K-NearestNeighbors",
      Lasso = "Lasso",
      Ridge_Regression = "Ridge_Regression",
      mlp = "mlp",
      tabnet = "tabnet"
    ),
    cross_validation = FALSE,
    optimizer_name = "adam",
    use_amp = FALSE,
    max_grad_norm = 1,
    auto_class_weights = FALSE,
    cnn_neurons_per_layer = c(8L),
    cnn_kernel_size = 3L,
    cnn_dense_layers = c(8L),
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel = 2L,
    cnn_pool_stride = 2L,
    cnn_pool_padding = 0L,
    cnn_learning_rate = 1e-3,
    cnn_separable = FALSE,
    cnn_dilations = c(1L),
    cnn_use_se = FALSE,
    cnn_norm_type = "group",
    cnn_pool_type = "conv",
    cnn_use_global_pool = FALSE,
    resnet_neurons_per_block = c(8L),
    resnet_blocks = 1L,
    resnet_learning_rate = 1e-3,
    ft_d_model = 8L,
    ft_heads = 1L,
    ft_layers = 1L,
    ft_ff_mult = 2L,
    ft_dropout = 0.1,
    ft_token_dropout = 0,
    ft_use_cls = TRUE,
    saint_d_model = 8L,
    saint_heads = 1L,
    saint_layers = 1L,
    saint_ff_mult = 2L,
    saint_dropout = 0.1,
    saint_token_dropout = 0,
    saint_use_cls = TRUE,
    use_grouping = FALSE,
    group_trigger = 2048L,
    group_method = "auto",
    init_group_size = 64L,
    max_tokens = 1024L,
    kmeans_batch = 4096L,
    kmeans_iter = 10L,
    tabnet_steps = 2L,
    tabnet_feature_dim = 8L,
    tabnet_output_dim = 8L,
    tabnet_gamma = 1.5,
    tabnet_lambda_sparse = 1e-4,
    node_trees = 2L,
    node_depth = 2L,
    deepfm_k = 4L,
    deepfm_hidden = c(8L),
    dcn_layers = 1L,
    dcn_hidden = c(8L),
    nam_hidden = c(8L),
    nam_activation = "relu",
    nam_add_linear = TRUE,
    nam_l1 = 1e-4,
    moe_n_experts = 2L,
    moe_expert_hidden = c(8L),
    moe_gate_hidden = 8L,
    moe_temperature = 1,
    moe_sparse_topk = NA,
    moe_entropy_reg = 0,
    gp_use_variational = TRUE,
    gp_num_inducing = 8L,
    gp_feature_dim = 8L,
    gp_kernel = "rbf",
    gp_ard = TRUE,
    gp_lr_mult = 0.5,
    rff_features = 16L,
    rff_lengthscale = 1,
    rff_deep_hidden = c(8L),
    model_type = "mlp",
    epochs = 1L,
    batch_size = 4L,
    dropout = 0.1,
    l2_weight_decay = 1e-4,
    l2_regularizer_dp = 1e-4,
    dropout_rate = 0.1,
    batch_norm = FALSE,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    device = NULL,
    mlp_neurons_per_layer = c(8L),
    mlp_learning_rate = 1e-3,
    final_attention = FALSE,
    attention_across_multiple_layers = FALSE,
    heteroscedastic = FALSE
  )
}

test_that("true-prediction KNN receives non-Gaussian response family", {
  pheno <- data.frame(
    GID = paste0("g", 1:4),
    BinaryTrait = factor(c("no", "yes", "no", NA), levels = c("no", "yes")),
    stringsAsFactors = FALSE
  )
  pheno_clean <- list(BinaryTrait = pheno)
  geno <- matrix(
    seq_len(16),
    nrow = 4,
    dimnames = list(pheno$GID, paste0("m", 1:4))
  )
  ml_dat_res <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = "BinaryTrait",
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )
  seen_response_family <- NULL

  local_mocked_bindings(
    AI_knn = function(pheno_object = NULL,
                      response = NULL,
                      geno_omic_object = NULL,
                      geno_omic_test_object = NULL,
                      response_family = "gaussian",
                      gen_name = NULL,
                      ...) {
      seen_response_family <<- response_family
      # Size the fake output from the matrices received (the test set may be empty).
      ids <- c(rownames(geno_omic_object), rownames(geno_omic_test_object))
      n <- length(ids)
      p_yes <- rep_len(c(0.2, 0.8), n)
      list(
        model_parameters = data.frame(stat = "model", summary = "K-NearestNeighbors"),
        predicted_values = data.frame(
          GID = ids,
          Predicted_class = ifelse(p_yes > 0.5, "yes", "no"),
          Train_Test_Label = c(rep("Train", NROW(geno_omic_object)), rep("Test", NROW(geno_omic_test_object))),
          Observed_class = rep_len(c("no", "yes"), n),
          Prediction_confidence = pmax(p_yes, 1 - p_yes),
          Classification_uncertainty = pmin(p_yes, 1 - p_yes),
          Reliability = pmax(p_yes, 1 - p_yes),
          Reliability_remarks = "Moderate Confidence",
          Probability_no = 1 - p_yes,
          Probability_yes = p_yes,
          stringsAsFactors = FALSE
        ),
        variance_components = PredictProR:::gp_empty_variance_components()
      )
    },
    summary_statistics_AI = function(...) {
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )

  ctx <- make_ml_runtime_ctx(pheno_clean, ml_dat_res, "K-NearestNeighbors")
  ctx$response <- "BinaryTrait"
  ctx$response_family <- "binary"

  log <- utils::capture.output(PredictProR:::gp_run_best_model_task(ctx))

  expect_identical(seen_response_family, "binary")
  expect_false(any(grepl("model fitting failed|Error in K-NearestNeighbors", log)))
})

test_that("predict_with_model dispatches all classical ML CV wrappers with shared tst", {
  seen <- list()

  local_mocked_bindings(
    AI_randomforest_cv = function(y, omics, tst, ...) {
      seen[["RandomForest"]] <<- tst
      rep(1, length(tst))
    },
    AI_pls_cv = function(y, omics, tst, ...) {
      seen[["PartialLeastSquare"]] <<- tst
      rep(2, length(tst))
    },
    AI_ridge_regression_cv = function(y, omics, tst, ...) {
      seen[["Ridge_Regression"]] <<- tst
      rep(3, length(tst))
    },
    AI_lasso_cv = function(y, omics, tst, ...) {
      seen[["Lasso"]] <<- tst
      rep(4, length(tst))
    },
    AI_svm_cv = function(y, omics, tst, ...) {
      seen[["SupportVectorMachine"]] <<- tst
      rep(5, length(tst))
    },
    AI_knn_cv = function(y, omics, tst, ...) {
      seen[["K-NearestNeighbors"]] <<- tst
      rep(6, length(tst))
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, 3)
  omics <- matrix(1:9, nrow = 3)
  tst <- 2L
  params <- list(
    scaling = FALSE,
    centering = FALSE,
    omic_count = 1L,
    ntree = 10L,
    ncomp = 2L,
    C_value = 1,
    degree_value = 3,
    scale_value = 1,
    offset_value = 0,
    k = 3L
  )

  models <- c("RandomForest", "PartialLeastSquare", "Ridge_Regression", "Lasso", "SupportVectorMachine", "K-NearestNeighbors")
  for (model in models) {
    res <- PredictProR:::predict_with_model(
      model = model,
      y = y,
      omics_data = omics,
      tst = tst,
      additional_params = params
    )
    expect_length(res, 1)
    expect_identical(seen[[model]], tst)
  }
})

test_that("predict_with_model dispatches deep-learning CV wrappers with model_type preserved", {
  seen <- list()

  local_mocked_bindings(
    deep_learning_model = function(y, omics, tst, model_type, crossval, ...) {
      seen[[model_type]] <<- list(tst = tst, crossval = crossval)
      rep(7, length(tst))
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, 3)
  omics <- matrix(1:9, nrow = 3)
  tst <- 3L
  params <- list(
    scaling = FALSE,
    centering = FALSE,
    crossval = TRUE,
    omic_count = 1L,
    early_stop = FALSE,
    deep_learning_model = NULL,
    optimizer_name = "adam",
    use_amp = FALSE,
    max_grad_norm = 1,
    auto_class_weights = FALSE,
    cnn_neurons_per_layer = c(8L),
    cnn_kernel_size = 3L,
    cnn_dense_layers = c(8L),
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel = 2L,
    cnn_pool_stride = 2L,
    cnn_pool_padding = 0L,
    cnn_learning_rate = 1e-3,
    cnn_separable = FALSE,
    cnn_dilations = c(1L),
    cnn_use_se = FALSE,
    cnn_norm_type = "group",
    cnn_pool_type = "conv",
    cnn_use_global_pool = FALSE,
    resnet_neurons_per_block = c(8L),
    resnet_blocks = 1L,
    resnet_learning_rate = 1e-3,
    ft_d_model = 8L,
    ft_heads = 1L,
    ft_layers = 1L,
    ft_ff_mult = 2L,
    ft_dropout = 0.1,
    ft_token_dropout = 0,
    ft_use_cls = TRUE,
    saint_d_model = 8L,
    saint_heads = 1L,
    saint_layers = 1L,
    saint_ff_mult = 2L,
    saint_dropout = 0.1,
    saint_token_dropout = 0,
    saint_use_cls = TRUE,
    use_grouping = FALSE,
    group_trigger = 2048L,
    group_method = "auto",
    init_group_size = 64L,
    max_tokens = 1024L,
    kmeans_batch = 4096L,
    kmeans_iter = 10L,
    tabnet_steps = 2L,
    tabnet_feature_dim = 8L,
    tabnet_output_dim = 8L,
    tabnet_gamma = 1.5,
    tabnet_lambda_sparse = 1e-4,
    node_trees = 2L,
    node_depth = 2L,
    deepfm_k = 4L,
    deepfm_hidden = c(8L),
    dcn_layers = 1L,
    dcn_hidden = c(8L),
    nam_hidden = c(8L),
    nam_activation = "relu",
    nam_add_linear = TRUE,
    nam_l1 = 1e-4,
    moe_n_experts = 2L,
    moe_expert_hidden = c(8L),
    moe_gate_hidden = 8L,
    moe_temperature = 1,
    moe_sparse_topk = NA,
    moe_entropy_reg = 0,
    gp_use_variational = TRUE,
    gp_num_inducing = 8L,
    gp_feature_dim = 8L,
    gp_kernel = "rbf",
    gp_ard = TRUE,
    gp_lr_mult = 0.5,
    rff_features = 16L,
    rff_lengthscale = 1,
    rff_deep_hidden = c(8L),
    epochs = 1L,
    batch_size = 4L,
    dropout = 0.1,
    l2_weight_decay = 1e-4,
    l2_regularizer_dp = 1e-4,
    dropout_rate = 0.1,
    batch_norm = FALSE,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 1L,
    device = NULL,
    mlp_neurons_per_layer = c(8L),
    mlp_learning_rate = 1e-3,
    final_attention = FALSE,
    attention_across_multiple_layers = FALSE,
    heteroscedastic = FALSE
  )

  models <- c("mlp", "tabnet", "cnn")
  for (model in models) {
    res <- PredictProR:::predict_with_model(
      model = model,
      y = y,
      omics_data = omics,
      tst = tst,
      additional_params = utils::modifyList(params, list(model_type = model))
    )
    expect_length(res, 1)
    expect_identical(seen[[model]][["tst"]], tst)
    expect_true(isTRUE(seen[[model]][["crossval"]]))
  }
})

test_that("torch_fit_filtered preserves deep-learning kwargs when torch_fit_model has dots", {
  seen <- NULL

  local_mocked_bindings(
    torch_fit_model = function(X, y, ...) {
      seen <<- list(...)
      list(ok = TRUE)
    },
    .package = "PredictProR"
  )

  PredictProR:::torch_fit_filtered(list(
    X = matrix(1, nrow = 2),
    y = c(1, 2),
    model_type = "mlp",
    compile_model = FALSE,
    epochs = 1L
  ))

  expect_identical(seen$model_type, "mlp")
  expect_false(seen$compile_model)
  expect_identical(seen$epochs, 1L)
})

test_that("true-prediction deep-learning runtime passes trait-specific test objects", {
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
    1:9,
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("m1", "m2", "m3"))
  )

  ml_dat_res <- PredictProR:::ML_data_processing(
    pheno_clean = pheno_clean,
    response = c("Yield", "Height"),
    geno_clean = geno,
    gen_name = "GID",
    omic_clean = list()
  )

  seen_test_ids <- list()

  local_mocked_bindings(
    deep_learning_model = function(pheno_object = NULL, geno_omic_object = NULL, geno_omic_test_object = NULL, response = NULL, model_type = NULL, ...) {
      seen_test_ids[[as.character(response)]] <<- if (is.null(geno_omic_test_object)) character() else rownames(geno_omic_test_object)
      ids <- c(rownames(geno_omic_object), if (!is.null(geno_omic_test_object)) rownames(geno_omic_test_object))
      list(
        model_parameters = data.frame(stat = "dummy", summary = model_type),
        predicted_values = data.frame(
          name = ids,
          Predicted_value = seq_along(ids),
          Train_Test_Label = c(rep("Train", nrow(geno_omic_object)), rep("Test", length(ids) - nrow(geno_omic_object)))
        )
      )
    },
    summary_statistics_AI = function(...) list(summary_statistics = data.frame(stat = "ok", summary = "ok")),
    .package = "PredictProR"
  )

  base_ctx <- make_ml_runtime_ctx(
    pheno_clean = pheno_clean,
    ml_dat_res = ml_dat_res,
    gs_model = "mlp",
    canonical_names = c("mlp", "tabnet", "cnn")
  )

  PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Yield", GS_model = "mlp")))
  PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Height", GS_model = "tabnet")))

  expect_identical(unname(seen_test_ids[["Yield"]]), "g2")
  expect_identical(unname(seen_test_ids[["Height"]]), "g1")
})

test_that("public model_execute keeps ML task context and returns prediction_stability output", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    stringsAsFactors = FALSE
  )

  geno <- data.frame(
    GID = c("g1", "g2", "g3"),
    m1 = c(10, 20, 30),
    m2 = c(11, 21, 31),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    process_geno_data = function(...) {
      list(geno_model_ready = geno)
    },
    AI_RidgeRegression_Lasso = function(pheno_object = NULL, response = NULL, geno_omic_object = NULL, geno_omic_test_object = NULL, ...) {
      ids <- c(rownames(geno_omic_object), rownames(geno_omic_test_object))
      pred <- data.frame(
        GID = ids,
        Predicted_value = c(1.1, 2.2, 3.3),
        Train_Test_Label = c("Train", "Train", "Test"),
        Standard_error = c(0.1, 0.1, 0.2),
        PEV = c(0.01, 0.01, 0.04),
        lower_bound = c(1.0, 2.1, 3.0),
        upper_bound = c(1.2, 2.3, 3.6),
        Uncertainty = c(0.2, 0.2, 0.4),
        Uncertainty_remarks = c("Low", "Low", "Moderate"),
        Prediction_stability = c(0.95, 0.95, 0.8),
        Prediction_stability_remarks = c("Stable", "Stable", "Moderately Stable"),
        Prediction_stability_percentage = c(100, 100, 100),
        Reliability_basis = rep(
          "not identifiable for ML/DL; not genetic reliability; Prediction_stability is a descriptive index",
          3L
        ),
        stringsAsFactors = FALSE
      )
      list(
        model_parameters = data.frame(stat = "lambda", summary = 1),
        predicted_values = pred
      )
    },
    summary_statistics_AI = function(...) {
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )

  res <- PredictProR::model_execute(
    pheno_data = pheno,
    geno_data = geno,
    response = "Yield",
    gen_name = "GID",
    GS_model = "Ridge_Regression",
    para_tunning = FALSE,
    cross_validation = FALSE,
    system_database = TRUE,
    qc_filtering = FALSE,
    ld_prunning_qc = FALSE,
    verbose = FALSE,
    n_bootstrap = 2
  )

  model_result <- res$model_results_by_model[["Ridge_Regression"]][["Yield"]]
  expect_false(is.null(model_result))
  expect_true("predicted_values" %in% names(model_result))
  expect_identical(
    colnames(model_result$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns()
  )
  expect_true("Prediction_stability" %in% colnames(model_result$predicted_values))
  expect_equal(model_result$predicted_values$Prediction_stability, c(0.95, 0.95, 0.8))
  expect_equal(
    model_result$predicted_values$Reliability,
    model_result$predicted_values$Prediction_stability
  )
  expect_match(
    unique(stats::na.omit(model_result$predicted_values$Reliability_basis)),
    "not genetic reliability",
    fixed = TRUE
  )
})
