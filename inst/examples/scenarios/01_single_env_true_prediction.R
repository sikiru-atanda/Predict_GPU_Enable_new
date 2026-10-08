# =============================================================================
# Scenario 1: one environment, one trait - TRUE PREDICTION
# =============================================================================
# Lines with NA phenotype are predicted from the lines with phenotypes.
# One full model_execute() call, run once per model: Bayesian marker models,
# Bayesian kernel models, ASReml GBLUP, classical machine learning, deep
# learning and Gaussian process models.
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_single_env(n = 80, p = 400)
pheno <- dat$pheno            # columns GID, Yield; test lines have Yield = NA
Geno.data <- dat$geno         # lines x markers, 0/1/2, row names = GID
response <- "Yield"

single_env_single_trait_models <- c(
  "BayesA", "BayesB", "BayesC", "BL", "BRR",                  # Bayesian marker regression (BGLR)
  "GBLUP_BRR", "RKHS",                                        # Bayesian kernel models (BGLR)
  "GBLUP",                                                    # ASReml-R (engine = "asreml")
  "Lasso", "Ridge_Regression", "PartialLeastSquare", "SupportVectorMachine",
  "K-NearestNeighbors", "RandomForest", "Xgboost", "CatBoost", "LightGBM",   # machine learning
  "DenseNeuralNet", "TabAttention", "TabNet", "LightTreeNet",
  "FactorNet", "CrossNet", "MixtureOfExperts", "GPNet",
  "DenseAttentionNet", "Conv1DNet", "ResNet",                 # deep learning
  "Kernel-GBLUP", "Gaussian-Process-GBLUP"                    # Gaussian process
)
# Not in this example:
#   TabTransformer - works, but on very small panels (here ~50 lines after the
#                    validation split) this large transformer is unstable and
#                    less accurate than GBLUP; use it with larger panels. It is
#                    used in the multi-environment / multi-trait / hybrid /
#                    classification workflows (scenarios 03-06).
#   NeuralAdditive - works, but fits one small network per marker: very slow
#                    with 100 epochs on 400 markers.
#   FA-GBLUP (needs 2+ environments or traits), Scalable-GBLUP (multi-trait only).

results <- list()
for (m in single_env_single_trait_models) {
  label <- paste("Single-env TP:", m)
  if (!is.null(need <- missing_software(needs_for_model(m)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,
      gmatrix = NULL,          # constructed from geno_data
      gkernel = NULL,
      kernel_list = NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Explicit single-trait / single-environment mode
      multi_trait_gp = FALSE,
      multi_trait_bayes = FALSE,
      multi_trait_asreml = FALSE,
      multi_trait_ml = FALSE,
      multi_trait_dl = FALSE,
      met_ml_dl = FALSE,
      hybrid_asreml = FALSE,
      hybrid_bayes = FALSE,
      hybrid_gp = FALSE,
      hybrid_ml = FALSE,
      hybrid_dl = FALSE,
      heter_groups = NULL,
      heter_resid = FALSE,
      bayes_kernel_heter_resid = NULL,
      var_cov_str = NULL,

      # Mixed-model formulas (Bayesian and ASReml models)
      fixed = ~1,
      random = ~GID,
      cova = NULL,
      weights = NULL,

      # Candidate model
      GS_model = m,

      # ASReml settings used by GBLUP
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",
      kernel_method = NULL,
      inverse = TRUE,
      bending = TRUE,
      bend_value = 0.01,
      blending = FALSE,
      kernel_check_level = "auto",
      kernel_fix_method = "auto",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = 0.10,
      ind_call_rate_threshold = 0.90,
      snp_call_rate_threshold = 0.90,
      impute = TRUE,
      imputation_method = "knn",
      impute_knn_k = 5L,
      ploidy = "auto",

      # LD pruning (ld_pruning = TRUE needs PLINK and file input)
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation: fit once, predict the NA lines
      cross_validation = FALSE,

      # Evaluation and ranking (used when cross_validation = TRUE)
      metric_for_ranking = "auto",
      ranking_tie_breakers = NULL,

      # This script lists every model setting so you can see what exists.
      # Values are the package defaults, except those marked "example", which
      # are made small so the script runs in minutes.

      # General hyperparameter tuning (the *_paras_tunning grids below are
      # searched only when para_tunning = TRUE)
      para_tunning = FALSE,
      AI_cv_nfolds = 5L,
      resample_method_tune = "cv",
      number_of_fold_tune = 5L,

      # XGBoost
      iteration = 30L,             # example; default 100
      xgb_booster = "gbtree",
      xgb_paras_tunning = list(
        Iter_tune = seq(100, 500, 100),
        learning_rate_tune = c(0.01, 0.1, 0.1),
        max_depth = c(3, 6, 9),
        rate_drop = c(0.1, 0.15, 0.2),
        skip_drop = c(0.4, 0.5, 0.55),
        xgb_gamma = c(0, 0.01, 0.1),
        colsample_bytree = c(0.5, 0.75, 1),
        min_child_weight = c(1, 3, 5),
        subsample = c(0.5, 0.75, 1),
        L2_tune = c(0, 0.5, 1),
        L1_tune = c(0, 0.5, 1)
      ),
      xgb_objective = "reg:squarederror",
      xgb_nthread = 1L,
      early_stop_for_iteration_xgb = FALSE,

      # Random forest
      ntree = 50L,                 # example; default 500
      nodesize = NULL,             # NULL = automatic
      importance = TRUE,
      rf_n_jobs = 1L,
      rf_paras_tunning = list(
        mtry = TRUE,
        ntree = c(500, 1000, 1500),
        nodesize = c(1, 5, 10),
        maxnodes = c(30, 50, NULL)
      ),

      # Partial least squares, Lasso, Ridge
      ncomp = 3L,
      pls_paras_tunning = list(ncomp = 10),
      lasso_paras_tunning = list(lambda_tune = seq(1e-06, 0.9, length.out = 100)^4),
      lambda_rr = NULL,            # NULL = automatic Ridge grid

      # CatBoost
      catboost_iterations = 30L,   # example; default 500

      # LightGBM
      lightgbm_nrounds = 30L,      # example; default 100
      lightgbm_learning_rate = 0.05,
      lightgbm_num_leaves = 31L,
      lightgbm_feature_fraction = 1,
      lightgbm_bagging_fraction = 1,
      lightgbm_min_data_in_leaf = 20L,
      lightgbm_lambda_l1 = 0,
      lightgbm_lambda_l2 = 0,
      lightgbm_nthread = 1L,

      # Bayesian model controls
      nIter = 600L,                # example; use 16000 for real data
      burnIn = 200L,               # example; use 2000
      thin = 5L,

      # GP controls
      gp_backend = "auto",
      gp_output_level = "predict_with_se",   # default "predict_only"; adds standard errors
      gp_return_se = TRUE,                    # default FALSE
      gp_full_vc = FALSE,
      gp_varcomp_mode = "reml",
      gp_prediction_output = "all",
      gp_iters = NULL,
      gp_lr = NULL,
      lowrank_eps_trace = 1e-6,
      lowrank_max_rank = NULL,
      lowrank_jitter = NULL,

      # Deep-learning controls (all deep-learning models)
      epochs = 100L,                 # default 10; attention models need ~100 to learn
      batch_size = 16L,            # example; default 64
      validation_split = 0.2,
      early_stop = TRUE,
      deterministic = TRUE,
      random_seed = 20260903L,
      optimizer_name = "adam",
      use_amp = TRUE,
      max_grad_norm = 1,
      device = NULL,               # NULL = GPU when available, else CPU
      dl_internal_calibration = FALSE,   # example; default TRUE (adds calibration fits)

      # DenseNeuralNet
      mlp_neurons_per_layer = as.integer(c(128, 64)),
      mlp_learning_rate = 0.001,
      dropout = 0.2,
      batch_norm = TRUE,

      # TabTransformer
      ft_d_model = 192L,
      ft_heads = 8L,
      ft_layers = 3L,
      ft_ff_mult = 4L,
      ft_dropout = 0.1,
      ft_token_dropout = 0,
      ft_use_cls = TRUE,
      use_grouping = FALSE,
      group_trigger = 2048L,
      group_method = "auto",
      init_group_size = 64L,
      max_tokens = 1024L,

      # Prediction uncertainty
      n_bootstrap = 10,            # use 30
      confidence_level = 0.95,

      # Scaling
      scaling = TRUE,
      centering = FALSE,

      # Memory and scheduling
      num_cores = NULL,
      parallel_mode = "auto",
      parallel_backend_prefer_fork = FALSE,
      random_state = 20260903L,

      # Output (system_database = FALSE also writes a results folder)
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
  results[[m]] <- res
}

# ---- Look at one result -------------------------------------------------------------
res <- results[["GBLUP_BRR"]]
if (!inherits(res, "try-error")) {
  pred <- res$model_results$predicted_values
  print(head(pred[, c("GID", "Train_Test_Label", "Observed_value", "Predicted_value", "Standard_error", "Reliability")]))
  # accuracy against the simulated truth for the test lines
  test_rows <- pred$Train_Test_Label == "Test"
  truth <- dat$truth$Yield[match(pred$GID[test_rows], dat$truth$GID)]
  cat(sprintf("GBLUP_BRR accuracy on test lines: r = %.2f\n", stats::cor(pred$Predicted_value[test_rows], truth)))
}

scenario_report()
