test_that("response family inference detects binary and ordinal traits", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Bin = factor(c("no", "yes", "no", "yes")),
    stringsAsFactors = FALSE
  )
  pheno$Ord <- ordered(c("low", "medium", "high", "medium"), levels = c("low", "medium", "high"))

  expect_identical(
    PredictProR:::gp_resolve_response_family_for_responses(pheno, "Bin", "auto"),
    "binary"
  )
  expect_identical(
    PredictProR:::gp_resolve_response_family_for_responses(pheno, "Ord", "auto"),
    "ordinal"
  )
})

test_that("response family inference detects binary numeric labels without treating numeric regression traits as classification", {
  pheno <- data.frame(
    GID = paste0("g", 1:6),
    Bin01 = c(0, 1, 0, 1, NA, 0),
    Bin12 = c(1, 2, 1, 2, NA, 1),
    NumericScore = c(1L, 2L, 3L, 1L, NA, 2L),
    Continuous = c(1.1, 2.3, 1.1, 4.8, NA, 5.4),
    stringsAsFactors = FALSE
  )

  expect_identical(PredictProR:::gp_resolve_response_family_for_responses(pheno, "Bin01", "auto"), "binary")
  expect_identical(PredictProR:::gp_resolve_response_family_for_responses(pheno, "Bin12", "auto"), "binary")
  expect_identical(PredictProR:::gp_resolve_response_family_for_responses(pheno, "NumericScore", "auto"), "gaussian")
  expect_identical(PredictProR:::gp_resolve_response_family_for_responses(pheno, "Continuous", "auto"), "gaussian")
})

test_that("response family inference detects multiclass labels from categorical columns", {
  pheno <- data.frame(
    GID = paste0("g", 1:6),
    Multi = factor(c("A", "B", "C", "A", NA, "B")),
    stringsAsFactors = FALSE
  )

  expect_identical(PredictProR:::gp_resolve_response_family_for_responses(pheno, "Multi", "auto"), "multiclass")
})

test_that("response family normalization handles empty and unsupported values cleanly", {
  expect_identical(PredictProR:::gp_normalize_response_family(character()), "auto")
  expect_identical(PredictProR:::gp_normalize_response_family(NA_character_), "auto")
  expect_identical(PredictProR:::gp_normalize_response_family("continuous"), "gaussian")
  expect_identical(PredictProR:::gp_normalize_response_family("categorical"), "multiclass")
  expect_error(
    PredictProR:::gp_normalize_response_family("multitask_regression"),
    "Unsupported response_family"
  )
})

test_that("evaluation_metrics supports binary and ordinal metric validation", {
  y_bin <- factor(c("no", "yes", "no", "yes"))
  p_yes <- c(0.1, 0.9, 0.3, 0.8)

  acc <- PredictProR::evaluation_metrics(
    y_observed = y_bin,
    y_predicted = p_yes,
    eval_metrics = "accuracy",
    response_family = "binary"
  )
  ll <- PredictProR::evaluation_metrics(
    y_observed = y_bin,
    y_predicted = p_yes,
    eval_metrics = "log_loss",
    response_family = "binary"
  )

  expect_equal(acc, 1)
  expect_true(is.finite(ll))
  expect_gt(ll, 0)
  ece <- PredictProR::evaluation_metrics(
    y_observed = y_bin,
    y_predicted = p_yes,
    eval_metrics = "ece",
    response_family = "binary"
  )
  expect_true(is.finite(ece))
  expect_gte(ece, 0)

  y_ord <- ordered(c("low", "medium", "high", "medium"), levels = c("low", "medium", "high"))
  pred_ord <- c("low", "high", "high", "medium")

  mae_cls <- PredictProR::evaluation_metrics(
    y_observed = y_ord,
    y_predicted = pred_ord,
    eval_metrics = "mean_absolute_error_class",
    response_family = "ordinal"
  )

  expect_equal(mae_cls, 0.25)
  expect_error(
    PredictProR::evaluation_metrics(
      y_observed = y_ord,
      y_predicted = pred_ord,
      eval_metrics = "mean_squared_error",
      response_family = "ordinal"
    ),
    "Unsupported eval_metrics"
  )

  prob_ord <- matrix(
    c(
      0.8, 0.1, 0.1,
      0.2, 0.6, 0.2,
      0.1, 0.2, 0.7,
      0.2, 0.5, 0.3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("low", "medium", "high"))
  )
  expect_true(is.finite(PredictProR::evaluation_metrics(
    y_observed = y_ord,
    y_predicted = prob_ord,
    eval_metrics = "ece",
    response_family = "ordinal"
  )))
})

test_that("model_execute rejects unsupported model and response_family combinations early", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Class = factor(c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  geno <- data.frame(
    GID = c("g1", "g2", "g3"),
    m1 = c(0, 1, 2),
    m2 = c(2, 1, 0),
    stringsAsFactors = FALSE
  )

  expect_error(
    suppressWarnings(PredictProR::model_execute(
      pheno_data = pheno,
      geno_data = geno,
      response = "Class",
      gen_name = "GID",
      GS_model = "PartialLeastSquare",
      response_family = "binary",
      system_database = TRUE,
      message = FALSE
    )),
    "does not support response_family='binary'"
  )
})

test_that("CatBoost and LightGBM support gaussian, binary, and multiclass response families", {
  expect_silent(PredictProR:::gp_validate_model_response_family("CatBoost", "gaussian"))
  expect_silent(PredictProR:::gp_validate_model_response_family("CatBoost", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("CatBoost", "multiclass"))
  expect_silent(PredictProR:::gp_validate_model_response_family("LightGBM", "gaussian"))
  expect_silent(PredictProR:::gp_validate_model_response_family("LightGBM", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("LightGBM", "multiclass"))
  expect_error(
    PredictProR:::gp_validate_model_response_family("CatBoost", "ordinal"),
    "does not support response_family='ordinal'"
  )
})

test_that("ML and DL support classification while GP and Bayesian kernels reject unsupported multiclass", {
  expect_silent(PredictProR:::gp_validate_model_response_family("RandomForest", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("RandomForest", "multiclass"))
  expect_silent(PredictProR:::gp_validate_model_response_family("mlp", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("tabnet", "multiclass"))
  expect_error(
    PredictProR:::gp_validate_model_response_family("KRR", "binary"),
    "does not support response_family='binary'"
  )
  expect_error(
    PredictProR:::gp_validate_model_response_family("RKHS", "multiclass"),
    "does not support response_family='multiclass'"
  )
})

test_that("GP display names are validated like their internal names", {
  for (m in c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")) {
    expect_silent(PredictProR:::gp_validate_model_response_family(m, "gaussian"))
    expect_error(
      PredictProR:::gp_validate_model_response_family(m, "binary"),
      paste0("Model '", m, "' does not support response_family='binary'"),
      fixed = TRUE
    )
    expect_error(
      PredictProR:::gp_validate_model_response_family(m, "multiclass"),
      "does not support response_family='multiclass'"
    )
  }
  expect_silent(PredictProR:::gp_validate_model_response_family("DenseNeuralNet", "multiclass"))
})

test_that("BGLR-backed Bayesian models allow binary and ordinal but reject multiclass", {
  expect_silent(PredictProR:::gp_validate_model_response_family("BayesA", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("BRR", "ordinal"))
  expect_silent(PredictProR:::gp_validate_model_response_family("RKHS", "binary"))
  expect_error(
    PredictProR:::gp_validate_model_response_family("BayesA", "multiclass"),
    "does not support response_family='multiclass'"
  )
})

test_that("bayes_mod_cv returns probabilities for binary and class labels for ordinal", {
  # Mock PredictProR's BGLR wrapper, which returns full-length probabilities
  # with the original class labels (real-BGLR behaviour, including the
  # observed-rows fit for NA lines, is covered in the bayes-classification tests).
  local_mocked_bindings(
    gp_bayes_bglr_fit = function(fam = "gaussian", y, response_type = "gaussian", ...) {
      if (response_type == "ordinal" && length(unique(stats::na.omit(y))) <= 2L) {
        return(list(
          probs = matrix(c(
            0.9, 0.1,
            0.8, 0.2,
            0.2, 0.8,
            0.7, 0.3
          ), nrow = 4, byrow = TRUE,
                         dimnames = list(NULL, c("no", "yes"))),
          levels = c("no", "yes")
        ))
      }
      list(
        probs = matrix(c(
          0.7, 0.2, 0.1,
          0.2, 0.6, 0.2,
          0.6, 0.3, 0.1,
          0.1, 0.2, 0.7
        ), nrow = 4, byrow = TRUE,
                       dimnames = list(NULL, c("low", "mid", "high"))),
        levels = c("low", "mid", "high")
      )
    },
    .package = "PredictProR"
  )

  y_bin <- factor(c("no", "yes", NA, NA), levels = c("no", "yes"))
  res_bin <- PredictProR:::bayes_mod_cv(
    y = y_bin,
    ETA = list(),
    weights = NULL,
    bayes_para = list(nIter = 10, burnIn = 2, thin = 1),
    tst = 3:4,
    bayes_model = "BayesA",
    bayes_trait = "Class",
    response_family = "binary"
  )
  expect_equal(as.numeric(res_bin), c(0.8, 0.3))
  expect_true(is.matrix(attr(res_bin, "probabilities")))

  y_ord <- ordered(c("low", "mid", NA, NA), levels = c("low", "mid", "high"))
  res_ord <- PredictProR:::bayes_mod_cv(
    y = y_ord,
    ETA = list(),
    weights = NULL,
    bayes_para = list(nIter = 10, burnIn = 2, thin = 1),
    tst = 3:4,
    bayes_model = "RKHS",
    bayes_trait = "Class",
    response_family = "ordinal"
  )
  expect_identical(unname(as.vector(res_ord)), c("mid", "low"))
  expect_true(is.matrix(attr(res_ord, "probabilities")))
  expect_identical(colnames(attr(res_ord, "probabilities")), c("low", "mid", "high"))
})

test_that("Bayesian classification output uses posterior class probabilities", {
  mod <- list(
    model = list(
      y = c("no", NA, "yes"),
      probs = matrix(
        c(0.8, 0.2,
          0.4, 0.6,
          0.3, 0.7),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(NULL, c("no", "yes"))
      ),
      levels = c("no", "yes")
    ),
    output_files_names = character()
  )
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Class = c("no", NA, "yes"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_bayes_classification_output(
    mod = mod,
    pheno_data = pheno,
    gen_name = "GID",
    response_family = "binary"
  )

  expect_true(all(c(
    "Predicted_value", "Prediction_confidence", "Classification_uncertainty",
    "Prob_no", "Prob_yes", "Train_Test_Label"
  ) %in% names(out$Predicted_value)))
  expect_identical(out$Predicted_value$Predicted_value, c("no", "yes", "yes"))
  expect_identical(out$Predicted_value$Train_Test_Label, c("Train", "Test", "Train"))
})

test_that("Bayesian ordinal output preserves the ordered response level order", {
  mod <- list(
    model = list(
      y = c(1L, 2L, NA_integer_),
      probs = matrix(
        c(
          0.1, 0.8, 0.1,
          0.1, 0.2, 0.7,
          0.7, 0.2, 0.1
        ),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(NULL, c("high", "low", "medium"))
      ),
      levels = c("high", "low", "medium")
    ),
    output_files_names = character()
  )
  pheno <- data.frame(GID = c("g1", "g2", "g3"), stringsAsFactors = FALSE)
  pheno$Trait <- ordered(
    c("low", "medium", NA),
    levels = c("low", "medium", "high")
  )

  out <- PredictProR:::gp_bayes_classification_output(
    mod = mod,
    pheno_data = pheno,
    gen_name = "GID",
    response = "Trait",
    response_family = "ordinal"
  )

  expect_identical(
    grep("^Prob_", names(out$Predicted_value), value = TRUE),
    c("Prob_low", "Prob_medium", "Prob_high")
  )
  expect_identical(out$Predicted_value$Predicted_value, c("low", "medium", "high"))
})

test_that("summary_statistics_bayes documents classification outputs", {
  mod <- list(
    model = list(
      y = c("no", "yes"),
      ETA = list()
    )
  )
  model_result <- list(
    Predicted_value = data.frame(
      GID = c("g1", "g2"),
      Predicted_value = c("no", "yes"),
      Prediction_confidence = c(0.8, 0.7),
      Classification_uncertainty = c(0.2, 0.3),
      Prob_no = c(0.8, 0.3),
      Prob_yes = c(0.2, 0.7),
      stringsAsFactors = FALSE
    )
  )

  out <- PredictProR:::summary_statistics_bayes(
    mod = mod,
    model_result = model_result,
    GS_model = "BayesA",
    response_family = "binary"
  )

  expect_true("Response_Family" %in% out$summary_statistics$stat)
  expect_true("Prediction_Confidence_Definition" %in% out$summary_statistics$stat)
})

test_that("deep_learning_model binary crossval returns probabilities", {
  seen <- list()

  local_mocked_bindings(
    gp_dl_bridge_fit_predict_raw = function(model_type, X_train, y_train, X_test, response_family, dl_args,
                                            class_levels = NULL) {
      seen <<- list(
        model_type = model_type,
        y_train = y_train,
        response_family = response_family,
        dl_args = dl_args,
        class_levels = class_levels
      )
      list(
        predictions = rep(0.8, nrow(X_test)),
        probabilities = matrix(c(0.2, 0.8), nrow = nrow(X_test), ncol = 2, byrow = TRUE),
        classes = c("no", "yes")
      )
    },
    .package = "PredictProR"
  )

  y <- factor(c("no", "yes", "no"), levels = c("no", "yes"))
  omics <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("m1", "m2"))
  )

  pred <- PredictProR:::deep_learning_model(
    y = y,
    omics_data = omics,
    tst = 1L,
    crossval = TRUE,
    response_family = "binary",
    model_type = "mlp",
    device = "cpu"
  )

  expect_equal(as.numeric(pred), 0.8)
  expect_identical(as.integer(seen$y_train), c(1L, 0L))
  expect_identical(seen$model_type, "mlp")
  expect_identical(seen$response_family, "binary")
  expect_identical(seen$class_levels, c("no", "yes"))
  expect_true(is.matrix(attr(pred, "probabilities")))
})

test_that("deep_learning_model multiclass crossval returns class labels", {
  seen <- list()

  local_mocked_bindings(
    gp_dl_bridge_fit_predict_raw = function(model_type, X_train, y_train, X_test, response_family, dl_args,
                                            class_levels = NULL) {
      seen <<- list(
        model_type = model_type,
        y_train = y_train,
        response_family = response_family,
        dl_args = dl_args,
        class_levels = class_levels
      )
      list(
        predictions = c("A", "B"),
        probabilities = matrix(c(3, 1, 0, 0, 2, 1), nrow = 2, byrow = TRUE),
        classes = c("A", "B", "C")
      )
    },
    .package = "PredictProR"
  )

  y <- factor(c("A", "B", "C"), levels = c("A", "B", "C"))
  omics <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("m1", "m2"))
  )

  pred <- PredictProR:::deep_learning_model(
    y = y,
    omics_data = omics,
    tst = 1:2,
    crossval = TRUE,
    response_family = "multiclass",
    model_type = "mlp",
    device = "cpu"
  )

  expect_identical(unname(as.vector(pred)), c("A", "B"))
  expect_true(is.matrix(attr(pred, "probabilities")))
  expect_identical(colnames(attr(pred, "probabilities")), c("A", "B", "C"))
  expect_identical(as.integer(seen$y_train), 2L)
  expect_identical(seen$model_type, "mlp")
  expect_identical(seen$response_family, "multiclass")
  expect_identical(seen$class_levels, c("A", "B", "C"))
})

test_that("deep_learning_model gaussian true prediction keeps train plus test rows aligned", {
  seen <- list()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed && is.integer(old_seed) && length(old_seed) > 1L) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  local_mocked_bindings(
    gp_dl_bridge_bootstrap = function(model_type, X_train, y_train, X_pred, n_bootstrap,
                                      response_family, dl_args, seed = NULL,
                                      class_levels = NULL, training_seeds = NULL,
                                      seed_aggregation = "mean", ...) {
      seen$prediction_rows <<- nrow(X_pred)
      center <- mean(y_train)
      list(
        bootstrap = matrix(
          rep(seq(center - 0.1, center + 0.1, length.out = nrow(X_pred)), each = n_bootstrap),
          nrow = n_bootstrap,
          ncol = nrow(X_pred)
        )
      )
    },
    gp_ml_predictor_unfamiliarity = function(predictor_matrix, train_test_label = NULL, ...) {
      seen$unfamiliarity_rows <<- nrow(predictor_matrix)
      list(
        unfamiliarity = rep(0.1, nrow(predictor_matrix)),
        remarks = rep("Familiar", nrow(predictor_matrix)),
        calibration_source = rep("mock", nrow(predictor_matrix))
      )
    },
    .package = "PredictProR"
  )

  x_train <- matrix(
    seq_len(12),
    nrow = 4,
    dimnames = list(paste0("tr", 1:4), paste0("m", 1:3))
  )
  x_test <- matrix(
    seq_len(6),
    nrow = 2,
    dimnames = list(paste0("te", 1:2), paste0("m", 1:3))
  )
  pheno <- data.frame(
    GID = rownames(x_train),
    Trait = c(1.2, 1.6, 2.0, 2.4),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::deep_learning_model(
    pheno_object = pheno,
    geno_omic_object = x_train,
    geno_omic_test_object = x_test,
    response = "Trait",
    gen_name = "GID",
    response_family = "gaussian",
    model_type = "mlp",
    para_tunning = FALSE,
    n_bootstrap = 3L,
    scaling = FALSE,
    centering = FALSE,
    system_database = TRUE,
    epochs = 1L,
    batch_size = 2L,
    mlp_neurons_per_layer = as.integer(c(4L))
  )

  pred <- out$predicted_values
  expect_equal(nrow(pred), nrow(x_train) + nrow(x_test))
  expect_equal(seen$prediction_rows, nrow(x_train) + nrow(x_test))
  expect_equal(seen$unfamiliarity_rows, nrow(x_train) + nrow(x_test))
  expect_identical(pred$Train_Test_Label, c(rep("Train", nrow(x_train)), rep("Test", nrow(x_test))))
  expect_true(all(c("Standard_error", "PEV", "lower_bound", "upper_bound") %in% names(pred)))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("python-backed ML CV helpers return probabilities or labels for CatBoost and LightGBM", {
  local_mocked_bindings(
    gp_py_ml_fit_predict = function(model_type = NULL, X_train = NULL, y_train = NULL, X_test = NULL,
                                    response_family = NULL, model_params = NULL, prefer_gpu = NULL,
                                    class_levels = NULL) {
      if (identical(response_family, "binary")) {
        out <- c(0.8, 0.2)
        attr(out, "probabilities") <- matrix(c(0.2, 0.8, 0.8, 0.2), ncol = 2, byrow = TRUE)
        return(out)
      }
      out <- c("A", "B")
      attr(out, "probabilities") <- matrix(c(0.7, 0.3, 0.1, 0.9), ncol = 2, byrow = TRUE)
      out
    },
    .package = "PredictProR"
  )

  y_bin <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  y_multi <- factor(c("A", "B", "A", "B"), levels = c("A", "B"))
  omics <- matrix(
    1:16,
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("m", 1:4))
  )

  res_cat_bin <- PredictProR:::AI_catboost_cv(
    y = y_bin,
    omics = omics,
    tst = c(1L, 2L),
    response_family = "binary"
  )
  expect_equal(as.numeric(res_cat_bin), c(0.8, 0.2))
  expect_true(is.matrix(attr(res_cat_bin, "probabilities")))

  res_lgb_multi <- PredictProR:::AI_lightgbm_cv(
    y = y_multi,
    omics = omics,
    tst = c(1L, 2L),
    response_family = "multiclass"
  )
  expect_identical(unname(as.vector(res_lgb_multi)), c("A", "B"))
  expect_true(is.matrix(attr(res_lgb_multi, "probabilities")))
})

test_that("multiclass cross-validation uses stored probabilities for every class metric", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Class = factor(c("A", "B", "C", "A"), levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  omics <- matrix(
    1:8,
    nrow = 4,
    dimnames = list(pheno$GID, c("m1", "m2"))
  )

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1:2),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, omics_data = NULL, tst = NULL, additional_params = NULL) {
      # Deliberately disagree with the probability argmax. CV metrics must use
      # the same probability-derived class that feeds public reliability.
      pred <- c("C", "C")
      attr(pred, "probabilities") <- matrix(
        c(0.8, 0.1, 0.1,
          0.2, 0.7, 0.1),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "B", "C"))
      )
      pred
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Class",
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = c("accuracy", "log_loss"),
      response_family = "multiclass",
      GS_model_cv = "mlp",
      ml_dat_res = list(
        merged_data = list(merge_data = omics),
        omic_count = 1L
      )
    ),
    verbose = FALSE
  )

  expect_true(is.data.frame(res[[1]]$yprob_cv_Reps_all))
  expect_true(all(c("A", "B", "C") %in% names(res[[1]]$yprob_cv_Reps_all)))
  expect_true(is.finite(res[[1]]$results_eval_metrics_reps$log_loss))
  expect_equal(res[[1]]$results_eval_metrics_reps$accuracy, 1)
})

test_that("binary cross-validation handles factor responses for deep models", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Class = factor(c("no", "yes", "no", "yes"), levels = c("no", "yes")),
    stringsAsFactors = FALSE
  )
  omics <- matrix(
    1:8,
    nrow = 4,
    dimnames = list(pheno$GID, c("m1", "m2"))
  )

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1:2),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, omics_data = NULL, tst = NULL, additional_params = NULL) {
      rep(0.8, length(tst))
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Class",
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = c("accuracy", "log_loss"),
      response_family = "binary",
      GS_model_cv = "mlp",
      ml_dat_res = list(
        merged_data = list(merge_data = omics),
        omic_count = 1L
      )
    ),
    verbose = FALSE
  )

  expect_length(res, 1)
  expect_true(is.data.frame(res[[1]]$results_eval_metrics_reps))
  expect_true(all(c("accuracy", "log_loss") %in% names(res[[1]]$results_eval_metrics_reps)))
})

test_that("multiclass cross-validation handles class-label predictions for deep models", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Class = factor(c("A", "B", "C", "A"), levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  omics <- matrix(
    1:8,
    nrow = 4,
    dimnames = list(pheno$GID, c("m1", "m2"))
  )

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1:2),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, omics_data = NULL, tst = NULL, additional_params = NULL) {
      c("A", "B")
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Class",
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = "accuracy",
      response_family = "multiclass",
      GS_model_cv = "mlp",
      ml_dat_res = list(
        merged_data = list(merge_data = omics),
        omic_count = 1L
      )
    ),
    verbose = FALSE
  )

  expect_length(res, 1)
  expect_true(is.data.frame(res[[1]]$results_eval_metrics_reps))
  expect_true("accuracy" %in% names(res[[1]]$results_eval_metrics_reps))
})

test_that("classical ML CV helpers return probabilities for binary response_family", {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  omics <- matrix(
    1:12,
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("m", 1:3))
  )

  local_mocked_bindings(
    gp_py_ml_fit_predict = function(model_type, X_train, y_train, X_test, response_family,
                                    model_params = list(), prefer_gpu = NULL, class_levels = NULL) {
      out <- rep(0.8, nrow(X_test))
      attr(out, "probabilities") <- cbind(no = rep(0.2, nrow(X_test)), yes = rep(0.8, nrow(X_test)))
      out
    },
    .package = "PredictProR"
  )

  xgb <- PredictProR:::AI_xgboost_cv(y, omics, tst = 1:2, response_family = "binary", scaling = FALSE, centering = FALSE)
  rf <- PredictProR:::AI_randomforest_cv(y, omics, tst = 1:2, response_family = "binary", scaling = FALSE, centering = FALSE)
  svm <- PredictProR:::AI_svm_cv(y, omics, tst = 1:2, response_family = "binary", scaling = FALSE, centering = FALSE)
  knn <- PredictProR:::AI_knn_cv(y, omics, tst = 1:2, response_family = "binary", scaling = FALSE, centering = FALSE)

  for (pred in list(xgb, rf, svm, knn)) {
    expect_equal(as.numeric(pred), c(0.8, 0.8))
    expect_true(is.matrix(attr(pred, "probabilities")))
  }
})

test_that("binary ML diagnostic plot builds confidence plots", {
  boot_results <- list(
    t = matrix(c(0.7, 0.8, 0.75, 0.85), nrow = 2, byrow = TRUE)
  )

  plot_obj <- PredictProR::diagnostic_plot_true_prediction(
    boot_results = boot_results,
    GID_names = c("g1", "g2"),
    CI_width_thresholds = c(0.33, 0.66),
    predictions = c(0.75, 0.8),
    standard_errors = c(0.05, 0.04),
    prediction_error_var = c(0.0025, 0.0016),
    genetic_var = 0.02,
    model_for_CI_cal = "ML",
    response_family = "binary",
    system_database = TRUE
  )

  expect_true(is.list(plot_obj))
  expect_true(all(c("predicted_vs_reliability", "prediction_inter_vs_reliability") %in% names(plot_obj)))
})
