# Classical single-environment ML wrappers.
# Training, prediction, and bootstrap refits route through the Python ML bridge.

train_predict_xgboost <- function(data_label_geno,
                                  indices,
                                  test_geno,
                                  params,
                                  nrounds,
                                  response_family = "gaussian",
                                  class_levels = NULL) {
  params <- utils::modifyList(params %||% list(), list(nrounds = nrounds))
  gp_py_ml_bootstrap_predict(
    data_label_geno = data_label_geno,
    indices = indices,
    test_geno = test_geno,
    model_type = "xgboost",
    response_family = response_family,
    model_params = params,
    class_levels = class_levels
  )
}

train_predict_randomForest <- function(data_label_geno,
                                       indices,
                                       test_geno,
                                       ntree = NULL,
                                       mtry = NULL,
                                       maxnodes = NULL,
                                       nodesize = NULL,
                                       rf_n_jobs = 1L,
                                       response_family = "gaussian",
                                       class_levels = NULL) {
  gp_py_ml_bootstrap_predict(
    data_label_geno = data_label_geno,
    indices = indices,
    test_geno = test_geno,
    model_type = "randomforest",
    response_family = response_family,
    model_params = list(
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    ),
    class_levels = class_levels
  )
}

knn_predict_boost <- function(data,
                              indices,
                              best_k,
                              geno_test,
                              response_family = "gaussian",
                              class_levels = NULL) {
  gp_py_ml_bootstrap_predict(
    data_label_geno = data,
    indices = indices,
    test_geno = geno_test,
    model_type = "knn",
    response_family = response_family,
    model_params = list(k = best_k),
    class_levels = class_levels
  )
}

pls_predict_boost <- function(data, indices, best_ncomp, geno_test) {
  gp_python_tabular_train_predict(
    data_label_geno = data,
    indices = indices,
    test_geno = geno_test,
    model_type = "pls",
    model_params = list(ncomp = best_ncomp),
    response_family = "gaussian"
  )
}

ridge_lasso_predict_boost <- function(data, indices, model_type, model_params, geno_test) {
  gp_python_tabular_train_predict(
    data_label_geno = data,
    indices = indices,
    test_geno = geno_test,
    model_type = model_type,
    model_params = model_params,
    response_family = "gaussian"
  )
}

gp_rf_mtry_tune_values <- function(geno_omic_object, rf_paras_tunning, mtry) {
  n_features <- ncol(geno_omic_object)
  values <- if (isTRUE(rf_paras_tunning$mtry)) {
    c(sqrt(n_features), sqrt(n_features) / 2, n_features / 3)
  } else {
    rf_paras_tunning$mtry %||% mtry %||% floor(sqrt(n_features))
  }
  values <- unique(as.integer(values))
  values[is.finite(values) & values > 0L]
}

#' Select the Best Model Based on Preferred Metrics
#'
#' @export
select_best_model <- function(model_list,
                              preferred_metrics = c("RMSE", "Accuracy", "MAE")) {
  scores <- lapply(model_list, function(model) {
    available_metrics <- intersect(names(model$results), preferred_metrics)
    if (!length(available_metrics)) {
      return(c(BestScore = NA, BestMetric = NA))
    }
    metric_scores <- vapply(available_metrics, function(metric) {
      values <- model$results[[metric]]
      if (identical(metric, "Accuracy")) {
        max(values, na.rm = TRUE)
      } else {
        min(values, na.rm = TRUE)
      }
    }, numeric(1))
    selection_scores <- metric_scores
    selection_scores[names(selection_scores) == "Accuracy"] <- -selection_scores[names(selection_scores) == "Accuracy"]
    best_metric <- names(which.min(selection_scores))[1]
    c(
      BestScore = metric_scores[[best_metric]],
      BestMetric = best_metric,
      SelectionScore = selection_scores[[best_metric]]
    )
  })

  scores <- scores[!vapply(scores, function(x) is.na(x[["BestScore"]]), logical(1))]
  if (!length(scores)) {
    return(list(BestModel = NULL, BestModelName = NA_character_, Metrics = NULL))
  }
  best_model_index <- which.min(vapply(scores, function(x) as.numeric(x[["SelectionScore"]]), numeric(1)))
  best_model_name <- names(scores)[best_model_index]
  list(
    BestModel = model_list[[best_model_name]],
    BestModelName = best_model_name,
    Metrics = scores[[best_model_index]][c("BestScore", "BestMetric")]
  )
}

#' Extreme Gradient Boosting Machine Learning Genomic Selection Pipeline
#'
#' Trains an XGBoost (gradient-boosted trees) regressor / classifier on a
#' genomic / omics feature matrix and returns predictions, uncertainty bands,
#' reliability diagnostics and feature importances, with optional hyper-parameter
#' tuning over the boosting / regularisation grid.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param response_family Response family: `"gaussian"`, `"binary"`, or
#'   `"multiclass"`. `"binomial"` is a binary alias; `"nominal"` and
#'   `"multinomial"` are multiclass aliases. Ordinal outcomes are not supported
#'   by this ML model.
#' @param AI_cv_nfolds Number of inner cross-validation folds for the
#'   hyper-parameter search.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param xgb_paras_tunning Named list of tuning grids over XGBoost
#'   hyper-parameters (iterations, learning rate, depth, drop rates, gamma,
#'   `colsample_bytree`, `min_child_weight`, subsample, L1 / L2).
#' @param resample_method_tune Resampling method for tuning (e.g. `"cv"`).
#' @param number_of_fold_tune Number of folds used by the tuning resampler.
#' @param learning_rate Boosting learning rate (`eta`) when
#'   `para_tunning = FALSE`.
#' @param max_depth Maximum tree depth.
#' @param subsample Row subsample ratio per boosting iteration.
#' @param xgb_booster Booster type (`"gbtree"`, `"gblinear"` or `"dart"`).
#' @param xgb_alpha L1 regularisation on weights.
#' @param xgb_lambda L2 regularisation on weights.
#' @param xgb_gamma Minimum loss reduction for a further split.
#' @param min_child_weight Minimum sum of instance weight in a child.
#' @param iteration Number of boosting iterations (`nrounds`).
#' @param xgb_rate_drop,xgb_skip_drop DART dropout rate and skip-dropout rate.
#' @param xgb_objective XGBoost learning objective (e.g.
#'   `"reg:squarederror"`).
#' @param xgb_sample_type,xgb_normalize_type DART sample / normalise type.
#' @param xgb_nthread Number of XGBoost threads.
#' @param N_feature_impo Number of top features to report by importance.
#' @param colsample_bytree Column-subsample ratio per tree.
#' @param early_stop_for_iteration_xgb Number of rounds without improvement to
#'   trigger early stopping (optional).
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained XGBoost model, test-set predictions,
#'   feature importances, uncertainty / reliability diagnostics and
#'   (where applicable) tuning results.
#' @export
AI_Xgb <- function(pheno_object = NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response = NULL,
                   gen_name = NULL,
                   response_family = "gaussian",
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   xgb_paras_tunning = list(
                     Iter_tune = seq(100, 500, 100),
                     learning_rate_tune = seq(0.01, 0.1, 0.01),
                     max_depth = seq(3, 15, 2),
                     rate_drop = seq(0.05, 0.5, 0.05),
                     skip_drop = seq(0.05, 1, 0.1),
                     xgb_gamma = seq(0, 1, 0.01),
                     colsample_bytree = seq(0.1, 1, 0.1),
                     min_child_weight = seq(1, 10, 2),
                     subsample = seq(0.2, 1, 0.1),
                     L2_tune = seq(0, 1, 0.01),
                     L1_tune = seq(0, 1, 0.01)
                   ),
                   resample_method_tune = "cv",
                   number_of_fold_tune = 5,
                   learning_rate = 0.01,
                   max_depth = 6,
                   subsample = 0.7,
                   xgb_booster = "gbtree",
                   xgb_alpha = 0.001,
                   xgb_lambda = 1.0,
                   xgb_gamma = 0.01,
                   min_child_weight = 1,
                   iteration = 100,
                   xgb_rate_drop = 0.1,
                   xgb_skip_drop = 0.5,
                   xgb_objective = "reg:squarederror",
                   xgb_sample_type = "uniform",
                   xgb_normalize_type = "tree",
                   xgb_nthread = 1L,
                   N_feature_impo = 10,
                   colsample_bytree = 0.7,
                   CI_width_thresholds = c(0.33, 0.66),
                   high_reliability_thres = 0.9,
                   low_reliability_thres = 0.5,
                   n_components = 20,
                   threshold = 100,
                   target = "test_set",
                   iqr_multiplier = 1.5,
                   interval_width_high_threshold = NULL,
                   interval_width_low_threshold = NULL,
                   interval_width_moderate_threshold = NULL,
                   n_bootstrap = 30,
                   early_stop_for_iteration_xgb = FALSE,
                   system_database = FALSE,
                   ...) {
  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  tune_grid <- if (isTRUE(para_tunning)) {
    list(
      nrounds = xgb_paras_tunning$Iter_tune %||% c(iteration),
      eta = xgb_paras_tunning$learning_rate_tune %||% c(learning_rate),
      max_depth = xgb_paras_tunning$max_depth %||% c(max_depth),
      xgb_gamma = xgb_paras_tunning$xgb_gamma %||% c(xgb_gamma),
      colsample_bytree = xgb_paras_tunning$colsample_bytree %||% c(colsample_bytree),
      min_child_weight = xgb_paras_tunning$min_child_weight %||% c(min_child_weight),
      subsample = xgb_paras_tunning$subsample %||% c(subsample),
      xgb_alpha = xgb_paras_tunning$L1_tune %||% c(xgb_alpha),
      xgb_lambda = xgb_paras_tunning$L2_tune %||% c(xgb_lambda)
    )
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = fam,
    model_type = "xgboost",
    model_label = "Xgboost",
    model_params = list(
      nrounds = iteration,
      eta = learning_rate,
      max_depth = max_depth,
      subsample = subsample,
      xgb_gamma = xgb_gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      xgb_alpha = xgb_alpha,
      xgb_lambda = xgb_lambda,
      xgb_booster = xgb_booster,
      xgb_nthread = xgb_nthread,
      xgb_rate_drop = xgb_rate_drop,
      xgb_skip_drop = xgb_skip_drop,
      xgb_objective = xgb_objective,
      xgb_sample_type = xgb_sample_type,
      xgb_normalize_type = xgb_normalize_type
    ),
    para_tunning = isTRUE(para_tunning),
    tune_param_grid = tune_grid,
    tune_folds = AI_cv_nfolds,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

#' Random Forest Machine Learning Genomic Selection Pipeline
#'
#' Trains a Random Forest regressor / classifier on a genomic / omics feature
#' matrix and returns predictions, uncertainty bands and reliability
#' diagnostics, with optional hyper-parameter tuning over `mtry`, `ntree`,
#' `nodesize`, `maxnodes`.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param response_family Response family: `"gaussian"`, `"binary"`, or
#'   `"multiclass"`. `"binomial"` is a binary alias; `"nominal"` and
#'   `"multinomial"` are multiclass aliases. Ordinal outcomes are not supported
#'   by this ML model.
#' @param AI_cv_nfolds Number of inner cross-validation folds for the
#'   hyper-parameter search.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param rf_paras_tunning Named list of tuning grids (`mtry`, `ntree`,
#'   `nodesize`, `maxnodes`).
#' @param ntree Number of trees when `para_tunning = FALSE`.
#' @param mtry Number of variables randomly sampled at each split (optional;
#'   default `sqrt(p)` for classification, `p/3` for regression).
#' @param maxnodes Maximum number of terminal nodes per tree (optional).
#' @param nodesize Minimum size of terminal nodes (optional).
#' @param importance Logical; compute variable-importance scores.
#' @param rf_n_jobs Integer; number of parallel workers for Random Forest.
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained Random Forest model, test-set predictions,
#'   variable importance, uncertainty / reliability diagnostics and
#'   (where applicable) tuning results.
#' @export
AI_randomForest <- function(pheno_object = NULL,
                            geno_omic_object = NULL,
                            geno_omic_test_object = NULL,
                            response = NULL,
                            gen_name = NULL,
                            response_family = "gaussian",
                            AI_cv_nfolds = 5,
                            message = TRUE,
                            scaling = TRUE,
                            centering = FALSE,
                            omic_count = NULL,
                            para_tunning = FALSE,
                            rf_paras_tunning = list(
                              mtry = TRUE,
                              ntree = c(500, 1000, 1500),
                              nodesize = c(1, 5, 10),
                              maxnodes = c(30, 50, NULL)
                            ),
                            ntree = 500,
                            mtry = NULL,
                            maxnodes = NULL,
                            nodesize = NULL,
                            importance = TRUE,
                            rf_n_jobs = 1L,
                            CI_width_thresholds = c(0.33, 0.66),
                            high_reliability_thres = 0.9,
                            low_reliability_thres = 0.5,
                            n_components = 20,
                            threshold = 100,
                            target = "test_set",
                            iqr_multiplier = 1.5,
                            interval_width_high_threshold = NULL,
                            interval_width_low_threshold = NULL,
                            interval_width_moderate_threshold = NULL,
                            n_bootstrap = 30,
                            system_database = FALSE,
                            ...) {
  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  tune_grid <- if (isTRUE(para_tunning)) {
    list(
      ntree = rf_paras_tunning$ntree %||% c(ntree),
      mtry = gp_rf_mtry_tune_values(geno_omic_object, rf_paras_tunning, mtry),
      nodesize = rf_paras_tunning$nodesize %||% c(nodesize %||% 1L),
      maxnodes = rf_paras_tunning$maxnodes %||% c(maxnodes)
    )
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = fam,
    model_type = "randomforest",
    model_label = "RandomForest",
    model_params = list(
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      importance = importance,
      rf_n_jobs = rf_n_jobs
    ),
    para_tunning = isTRUE(para_tunning),
    tune_param_grid = tune_grid,
    tune_folds = AI_cv_nfolds,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

#' Support Vector Machine Machine Learning Genomic Selection Pipeline
#'
#' Trains a support-vector regressor / classifier on a genomic / omics feature
#' matrix and returns predictions, uncertainty bands and reliability
#' diagnostics, with optional hyper-parameter tuning over the kernel and its
#' parameters.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param response_family Response family: `"gaussian"`, `"binary"`, or
#'   `"multiclass"`. `"binomial"` is a binary alias; `"nominal"` and
#'   `"multinomial"` are multiclass aliases. Ordinal outcomes are not supported
#'   by this ML model.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param AI_cv_nfolds Number of inner cross-validation folds for the
#'   hyper-parameter search.
#' @param svm_paras_tunning Named list of tuning grids (kernel, sigma, C,
#'   gamma, degree, scale, offset).
#' @param svm_kernel Kernel family when `para_tunning = FALSE`
#'   (`"Gaussian"`, `"Linear"`, `"Polynomial"`, `"Hyperbolic_tangent"`).
#' @param svm_type SVM problem type (e.g. `"eps-regression"`).
#' @param sigma_value Optional Gaussian kernel bandwidth. When both
#'   `sigma_value` and `gamma_value` are `NULL`, the backend uses its
#'   dimension-aware `gamma = "scale"` default.
#' @param C_value Cost / regularisation parameter `C`.
#' @param degree_value Polynomial-kernel degree.
#' @param scale_value Polynomial-kernel scaling factor.
#' @param offset_value Kernel offset (polynomial / tanh).
#' @param gamma_value RBF / polynomial kernel `gamma` (optional).
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained SVM model, test-set predictions,
#'   uncertainty / reliability diagnostics and (where applicable) tuning
#'   results.
#' @export
AI_svm <- function(pheno_object = NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response = NULL,
                   gen_name = NULL,
                   response_family = "gaussian",
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   AI_cv_nfolds = 5,
                   svm_paras_tunning = list(
                     kernel = c("Gaussian", "Linear", "Polynomial", "Hyperbolic_tangent"),
                     offset_value = seq(-2, 2, length.out = 5),
                     sigma = c(0.01, 0.05, 0.1),
                     C = c(1, 10, 100),
                     gamma_value = 10^seq(-4, -1, length.out = 4),
                     degree = c(3, 4),
                     scale = c(0.1, 1)
                   ),
                   svm_kernel = "Gaussian",
                   svm_type = "eps-regression",
                   sigma_value = NULL,
                   C_value = 1,
                   degree_value = 3,
                   scale_value = 1,
                   offset_value = 0,
                   gamma_value = NULL,
                   CI_width_thresholds = c(0.33, 0.66),
                   high_reliability_thres = 0.9,
                   low_reliability_thres = 0.5,
                   n_components = 20,
                   threshold = 100,
                   target = "test_set",
                   iqr_multiplier = 1.5,
                   interval_width_high_threshold = NULL,
                   interval_width_low_threshold = NULL,
                   interval_width_moderate_threshold = NULL,
                   n_bootstrap = 30,
                   system_database = FALSE,
                   ...) {
  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  tune_grid <- if (isTRUE(para_tunning)) {
    list(
      svm_kernel = svm_paras_tunning$kernel %||% c(svm_kernel),
      C_value = svm_paras_tunning$C %||% c(C_value),
      gamma_value = svm_paras_tunning$gamma_value %||% c(gamma_value %||% sigma_value),
      degree_value = svm_paras_tunning$degree %||% c(degree_value),
      offset_value = svm_paras_tunning$offset_value %||% c(offset_value)
    )
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = fam,
    model_type = "svm",
    model_label = "SupportVectorMachine",
    model_params = list(
      svm_type = svm_type,
      svm_kernel = svm_kernel,
      sigma_value = sigma_value,
      C_value = C_value,
      gamma_value = gamma_value,
      degree_value = degree_value,
      scale_value = scale_value,
      offset_value = offset_value
    ),
    para_tunning = isTRUE(para_tunning),
    tune_param_grid = tune_grid,
    tune_folds = AI_cv_nfolds,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

#' K-Nearest Neighbors Machine Learning Genomic Selection Pipeline
#'
#' Trains a k-nearest-neighbours regressor / classifier on a genomic / omics
#' feature matrix and returns predictions, uncertainty bands and reliability
#' diagnostics, with optional hyper-parameter tuning over `k`.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param response_family Response family: `"gaussian"`, `"binary"`, or
#'   `"multiclass"`. `"binomial"` is a binary alias; `"nominal"` and
#'   `"multinomial"` are multiclass aliases. Ordinal outcomes are not supported
#'   by this ML model.
#' @param AI_cv_nfolds Number of inner cross-validation folds for the
#'   hyper-parameter search.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param knn_paras_tunning Named list of tuning grids for KNN (e.g. `k`).
#' @param k Number of neighbours when `para_tunning = FALSE`.
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained model, test-set predictions, uncertainty /
#'   reliability diagnostics and (where applicable) tuning results.
#' @export
AI_knn <- function(pheno_object = NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response = NULL,
                   gen_name = NULL,
                   response_family = "gaussian",
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   knn_paras_tunning = list(k = seq(3, 21, by = 2)),
                   k = 5,
                   CI_width_thresholds = c(0.33, 0.66),
                   high_reliability_thres = 0.9,
                   low_reliability_thres = 0.5,
                   n_components = 20,
                   threshold = 100,
                   target = "test_set",
                   iqr_multiplier = 1.5,
                   interval_width_high_threshold = NULL,
                   interval_width_low_threshold = NULL,
                   interval_width_moderate_threshold = NULL,
                   n_bootstrap = 30,
                   system_database = FALSE,
                   ...) {
  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  tune_grid <- if (isTRUE(para_tunning)) {
    list(k = knn_paras_tunning$k %||% seq(3, 21, by = 2))
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = fam,
    model_type = "knn",
    model_label = "K-NearestNeighbors",
    model_params = list(k = k),
    para_tunning = isTRUE(para_tunning),
    tune_param_grid = tune_grid,
    tune_folds = AI_cv_nfolds,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

#' Partial Least Squares Machine Learning Genomic Selection Pipeline
#'
#' Trains a partial-least-squares regression on a genomic / omics feature matrix
#' and returns predictions, uncertainty bands, reliability diagnostics and
#' feature importances, with optional hyper-parameter tuning over the number of
#' components.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param ncomp Number of PLS components when `para_tunning = FALSE`.
#' @param pls_paras_tunning Named list of tuning grids for PLS (e.g. `ncomp`).
#' @param resample_method_tune Resampling method for tuning (e.g. `"cv"`).
#' @param N_feature_impo Number of top features to report by importance score.
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained PLS model, test-set predictions,
#'   uncertainty / reliability diagnostics, feature importances and
#'   (where applicable) tuning results.
#' @export
AI_pls <- function(pheno_object = NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response = NULL,
                   gen_name = NULL,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   ncomp = 3,
                   pls_paras_tunning = list(ncomp = 10),
                   resample_method_tune = "cv",
                   N_feature_impo = 10,
                   CI_width_thresholds = c(0.33, 0.66),
                   high_reliability_thres = 0.9,
                   low_reliability_thres = 0.5,
                   n_components = 20,
                   threshold = 100,
                   target = "test_set",
                   iqr_multiplier = 1.5,
                   interval_width_high_threshold = NULL,
                   interval_width_low_threshold = NULL,
                   interval_width_moderate_threshold = NULL,
                   n_bootstrap = 30,
                   system_database = FALSE,
                   ...) {
  pls_tune_values <- as.integer(pls_paras_tunning$ncomp %||% ncomp %||% 10)
  pls_tune_values <- pls_tune_values[is.finite(pls_tune_values) & pls_tune_values > 0L]
  if (!length(pls_tune_values)) {
    pls_tune_values <- 1L
  }

  model_params <- if (isTRUE(para_tunning)) {
    list()
  } else if (is.null(ncomp) || !is.numeric(ncomp)) {
    list(
      pls_auto_components = TRUE,
      pls_max_components = max(pls_tune_values)
    )
  } else {
    list(ncomp = as.integer(ncomp[[1]]))
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = "gaussian",
    model_type = "pls",
    model_label = "PartialLeastSquare",
    model_params = model_params,
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) {
      list(ncomp = unique(pls_tune_values))
    } else NULL,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

#' Ridge Regression and Lasso Machine Learning Genomic Selection Pipeline
#'
#' Trains an L1- (Lasso) or L2- (Ridge) regularised linear model on a
#' genomic / omics feature matrix and returns predictions, uncertainty bands and
#' reliability diagnostics, with optional hyper-parameter tuning over the
#' regularisation penalty `lambda`.
#'
#' @param pheno_object Phenotype object (data frame or list) carrying cleaned
#'   phenotypes.
#' @param geno_omic_object Training feature matrix (rows = individuals,
#'   columns = genomic / omics features).
#' @param geno_omic_test_object Test-set feature matrix with the same columns
#'   as `geno_omic_object`.
#' @param response Name of the response column in `pheno_object`.
#' @param gen_name Name of the genotype-ID column in `pheno_object`.
#' @param message Logical; print progress messages.
#' @param scaling Logical; scale features to unit variance before fitting.
#' @param centering Logical; mean-centre features before fitting.
#' @param omic_count Optional integer; number of omics layers contributing to
#'   the feature matrix.
#' @param AI_cv_nfolds Number of inner cross-validation folds for the
#'   hyper-parameter search.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting.
#' @param lasso_paras_tunning Named list of tuning grids (e.g. `lambda_tune`).
#' @param lambda_rr Optional explicit `lambda` value used when
#'   `para_tunning = FALSE`.
#' @param GS_model Which model to fit: `"Lasso"` (L1) or `"Ridge_Regression"`
#'   (L2). A vector of both runs each model.
#' @param CI_width_thresholds Length-2 numeric of quantile thresholds for
#'   classifying prediction-interval width as Low / Moderate / High.
#' @param high_reliability_thres Lower bound (in `[0,1]`) for High-reliability
#'   classification.
#' @param low_reliability_thres Upper bound (in `[0,1]`) for Low-reliability
#'   classification.
#' @param n_components Number of top features / components to report.
#' @param threshold Numeric threshold passed to post-prediction risk / outlier
#'   bands.
#' @param target Which rows to score in the post-prediction summary
#'   (e.g. `"test_set"`).
#' @param iqr_multiplier IQR multiplier used in outlier detection.
#' @param interval_width_high_threshold,interval_width_low_threshold,interval_width_moderate_threshold
#'   Optional numeric overrides for the High / Low / Moderate interval-width
#'   bands; otherwise derived from quantiles of `CI_width_thresholds`.
#' @param n_bootstrap Number of bootstrap resamples for prediction uncertainty
#'   and stability.
#' @param system_database Logical; when `TRUE`, write outputs to the system
#'   database path.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the trained Lasso / Ridge model(s), test-set
#'   predictions, uncertainty / reliability diagnostics and (where applicable)
#'   tuning results.
#' @export
AI_RidgeRegression_Lasso <- function(pheno_object = NULL,
                                     geno_omic_object = NULL,
                                     geno_omic_test_object = NULL,
                                     response = NULL,
                                     gen_name = NULL,
                                     message = TRUE,
                                     scaling = TRUE,
                                     centering = FALSE,
                                     omic_count = NULL,
                                     AI_cv_nfolds = 5,
                                     para_tunning = FALSE,
                                     lasso_paras_tunning = list(lambda_tune = seq(0.000001, 0.9, length.out = 100)^4),
                                     lambda_rr = NULL,
                                     GS_model = c("Lasso", "Ridge_Regression"),
                                     CI_width_thresholds = c(0.33, 0.66),
                                     high_reliability_thres = 0.9,
                                     low_reliability_thres = 0.5,
                                     n_components = 20,
                                     threshold = 100,
                                     target = "test_set",
                                     iqr_multiplier = 1.5,
                                     interval_width_high_threshold = NULL,
                                     interval_width_low_threshold = NULL,
                                     interval_width_moderate_threshold = NULL,
                                     n_bootstrap = 30,
                                     system_database = FALSE,
                                     ...) {
  GS_model <- as.character(GS_model[[1]])
  if (!GS_model %in% c("Lasso", "Ridge_Regression")) {
    stop("GS_model must be one of 'Lasso' or 'Ridge_Regression'.", call. = FALSE)
  }

  model_type <- if (identical(GS_model, "Lasso")) "lasso" else "ridge"
  model_params <- if (identical(model_type, "lasso")) {
    list()
  } else if (!is.null(lambda_rr)) {
    list(ridge_alpha = lambda_rr)
  } else {
    list()
  }

  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = "gaussian",
    model_type = model_type,
    model_label = GS_model,
    model_params = model_params,
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) {
      if (identical(model_type, "lasso")) {
        list(lasso_alpha = lasso_paras_tunning$lambda_tune)
      } else {
        list(ridge_alpha = lambda_rr %||% 10^seq(-6, 2, length.out = 25))
      }
    } else NULL,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}
