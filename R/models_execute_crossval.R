# Legacy compatibility wrapper retained while downstream callers migrate to the
# canonical cross-validation runtime in `R/main_crossvalidation_execution_logic.R`.
#
# The historical full implementation has been preserved under:
# `inst/legacy_runtime/models_execute_crossval_legacy_full.R`

full_dp_models <- setNames(
  c("Conv1DNet", "TabTransformer", "TabAttention", "TabNet", "LightTreeNet",
    "FactorNet", "CrossNet", "NeuralAdditive", "MixtureOfExperts", "GPNet",
    "DenseAttentionNet", "DenseNeuralNet", "ResNet"),
  c("cnn", "ft_transformer", "saint", "tabnet", "node",
    "deepfm", "dcnv2", "nam", "moe", "gp_dkl", "mlp_with_attention", "mlp", "resnet")
)

maximize_metrics <- c(
  "accuracy", "kendalls_tau", "spearman", "pearson", "correlation",
  "rsq", "r2", "balanced_accuracy", "precision", "recall",
  "specificity", "f1", "mcc", "macro_precision", "macro_recall",
  "macro_f1", "within_one_class_accuracy", "quadratic_weighted_kappa"
)
minimize_metrics <- c(
  "mean_squared_error", "bias", "root_mean_squared_error",
  "relative_squared_error", "mean_absolute_error", "mean_absolute_percent_error",
  "log_loss", "brier_score", "mean_absolute_error_class", "ece"
)

metric_worst_value <- function(metric, y_true = NULL) {
  metric <- tolower(metric)
  direction <- if (exists("gp_eval_metric_direction", mode = "function")) {
    gp_eval_metric_direction(metric)
  } else if (metric %in% maximize_metrics) {
    "maximize"
  } else if (metric %in% minimize_metrics) {
    "minimize"
  } else {
    NA_character_
  }
  if (identical(direction, "maximize")) {
    return(-1e12)
  }
  if (identical(direction, "minimize")) {
    v <- tryCatch(stats::var(y_true, na.rm = TRUE), error = function(e) NA_real_)
    v <- if (is.finite(v) && v > 0) v else 1
    return(1e12 * v)
  }
  1e12
}

safe_metric_value <- function(y_true, y_pred, metric, response_family = "gaussian",
                              positive_class = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_true)
  if (fam == "gaussian") {
    idx <- is.finite(y_true) & is.finite(y_pred)
  } else if (is.matrix(y_pred) || is.data.frame(y_pred)) {
    yp_mat <- as.matrix(y_pred)
    idx <- !is.na(y_true) & rowSums(!is.na(yp_mat)) > 0
  } else {
    idx <- !is.na(y_true) & !is.na(y_pred)
  }
  if (!any(idx)) {
    return(if (identical(fam, "gaussian")) metric_worst_value(metric, y_true) else NA_real_)
  }

  yt <- y_true[idx]
  yp <- if (is.matrix(y_pred) || is.data.frame(y_pred)) y_pred[idx, , drop = FALSE] else y_pred[idx]
  needs_corr <- tolower(metric) %in% c(
    "kendalls_tau", "accuracy", "spearman", "pearson", "correlation", "rsq", "r2"
  )
  if (fam == "gaussian" && length(yt) < if (needs_corr) 2L else 1L) {
    return(metric_worst_value(metric, y_true))
  }
  if (fam == "gaussian" && needs_corr && (stats::sd(yt) == 0 || stats::sd(yp) == 0)) {
    return(metric_worst_value(metric, y_true))
  }

  out <- tryCatch(
    evaluation_metrics(
      y_observed = yt,
      y_predicted = yp,
      eval_metrics = metric,
      response_family = fam,
      positive_class = positive_class
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(out)) {
    out <- if (identical(fam, "gaussian")) metric_worst_value(metric, y_true) else NA_real_
  }
  out
}

legacy_models_execute_crossval <- function(pheno_data = NULL,
                                           test_set = NULL,
                                           response = NULL,
                                           gen_name = NULL,
                                           test_size = NULL,
                                           random_state = NULL,
                                           replication = NULL,
                                           weights = NULL,
                                           selected_raw,
                                           gam_method = NULL,
                                           max_features = 50,
                                           k_value = 5,
                                           var_explained = 0.9,
                                           engine = NULL,
                                           model_prep_all_bayes_cv = NULL,
                                           asreml_models_prep_cv = NULL,
                                           ml_dat_res = NULL,
                                           heter_groups = NULL,
                                           verbose = FALSE,
                                           num_cores = NULL,
                                           nfolds = 5,
                                           cross_validation_meth = NULL,
                                           sampling_method = NULL,
                                           eval_metrics = NULL,
                                           bayes_model = NULL,
                                           GS_model_cv = NULL,
                                           scaling = FALSE,
                                           centering = TRUE,
                                           eta = 0.1,
                                           nrounds = 100,
                                           max_depth = 6,
                                           xgb_gamma = 4,
                                           subsample = 0.5,
                                           colsample_bytree = 1,
                                           xgb_alpha = 0.001,
                                           xgb_lambda = 1,
                                           min_child_weight = 1,
                                           early_stop_for_iteration_xgb = TRUE,
                                           xgb_booster = "gbtree",
                                           xgb_rate_drop = 0.1,
                                           xgb_skip_drop = 0.5,
                                           xgb_objective = "reg:squarederror",
                                           xgb_sample_type = "uniform",
                                           xgb_normalize_type = "tree",
                                           ncomp = 3,
                                           ntree = 500,
                                           k = 5,
                                           svm_kernel = "Gaussian",
                                           sigma_value = NULL,
                                           C_value = 1,
                                           degree_value = 3,
                                           scale_value = 1,
                                           offset_value = 1,
                                           num_hidden_layers = 1,
                                           neurons_per_layer = 64,
                                           learning_rate_dp = 0.001,
                                           early_stop = TRUE,
                                           crossval = TRUE,
                                           optimizer_name = "adam",
                                           use_amp = TRUE,
                                           max_grad_norm = 1.0,
                                           auto_class_weights = FALSE,
                                           cnn_neurons_per_layer = as.integer(c(64, 64, 64)),
                                           cnn_kernel_size = 3L,
                                           cnn_dense_layers = as.integer(c(256, 128, 64)),
                                           cnn_use_max_pool = FALSE,
                                           cnn_pool_kernel = 2L,
                                           cnn_pool_stride = 2L,
                                           cnn_pool_padding = 0L,
                                           cnn_learning_rate = 1e-3,
                                           cnn_separable = TRUE,
                                           cnn_dilations = c(1, 2, 4),
                                           cnn_use_se = TRUE,
                                           cnn_norm_type = "group",
                                           cnn_pool_type = "conv",
                                           cnn_use_global_pool = FALSE,
                                           resnet_neurons_per_block = as.integer(c(256, 128, 64)),
                                           resnet_blocks = 3,
                                           resnet_learning_rate = 1e-3,
                                           ft_d_model = 192L,
                                           ft_heads = 8L,
                                           ft_layers = 3L,
                                           ft_ff_mult = 4L,
                                           ft_dropout = 0.1,
                                           ft_token_dropout = 0.0,
                                           ft_use_cls = TRUE,
                                           saint_d_model = 128L,
                                           saint_heads = 8L,
                                           saint_layers = 3L,
                                           saint_ff_mult = 4L,
                                           saint_dropout = 0.1,
                                           saint_token_dropout = 0.0,
                                           saint_use_cls = TRUE,
                                           use_grouping = FALSE,
                                           group_trigger = 2048,
                                           group_method = "auto",
                                           init_group_size = 64,
                                           max_tokens = 1024,
                                           kmeans_batch = 4096,
                                           kmeans_iter = 100,
                                           tabnet_steps = 5L,
                                           tabnet_feature_dim = 64L,
                                           tabnet_output_dim = 64L,
                                           tabnet_gamma = 1.5,
                                           tabnet_lambda_sparse = 1e-4,
                                           node_trees = 8L,
                                           node_depth = 3L,
                                           deepfm_k = 16L,
                                           deepfm_hidden = as.integer(c(128, 64)),
                                           dcn_layers = 3L,
                                           dcn_hidden = as.integer(c(256, 128, 64)),
                                           nam_hidden = as.integer(c(32, 16)),
                                           nam_activation = "relu",
                                           nam_add_linear = TRUE,
                                           nam_l1 = 1e-4,
                                           moe_n_experts = 4L,
                                           moe_expert_hidden = as.integer(c(128, 64)),
                                           moe_gate_hidden = 128L,
                                           moe_temperature = 1.0,
                                           moe_sparse_topk = NA,
                                           moe_entropy_reg = 0.0,
                                           gp_use_variational = TRUE,
                                           gp_num_inducing = 256L,
                                           gp_feature_dim = 64L,
                                           gp_kernel = "rbf",
                                           gp_ard = TRUE,
                                           gp_lr_mult = 0.5,
                                           rff_features = 1024,
                                           rff_lengthscale = 1.0,
                                           rff_deep_hidden = c(128),
                                           model_type = "resnet",
                                           epochs = 10,
                                           batch_size = 64,
                                           dropout = 0.2,
                                           l2_weight_decay = 1e-4,
                                           l2_regularizer_dp = 0.001,
                                           dropout_rate = 0.5,
                                           batch_norm = TRUE,
                                           validation_split = 0.2,
                                           compile_model = FALSE,
                                           deterministic = TRUE,
                                           random_seed = 123,
                                           device = NULL,
                                           mlp_neurons_per_layer = as.integer(c(128, 64)),
                                           mlp_learning_rate = 1e-3,
                                           final_attention = TRUE,
                                           attention_across_multiple_layers = TRUE,
                                           heteroscedastic = TRUE,
                                           docker_nd_usage = FALSE,
                                           globals_max_GB = 4,
                                           parallel_mode = c("auto", "future", "sequential", "base_parallel", "foreach"),
                                           parallel_backend_prefer_fork = TRUE,
                                           sequential_models = NULL,
                                           ...) {
  gp_dispatch_legacy_crossvalidation(
    call_env = environment(),
    parallel_choices = c("auto", "future", "sequential", "base_parallel", "foreach"),
    dots = list(...)
  )
}
