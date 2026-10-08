# =============================================================================
# Scenario 2: one environment - CROSS-VALIDATION
# =============================================================================
# Compare models by cross-validation (CV); optionally refit the best model and
# predict the NA lines (cv_evaluation_only = FALSE).
# cross_validation_meth for one environment:
#   "K-Folds", "Stratified_K-Folds", "Repeated_K-Folds", "Repeated_Stratified_K-Folds",
#   "Hold_Out", "Stratified_Hold_Out", "Repeated_Hold_Out",
#   "Repeated_Stratified_Hold_Out", "Leave_one_Out"
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_single_env(n = 80, p = 400)
pheno <- dat$pheno
Geno.data <- dat$geno
response <- "Yield"

# ---- 1. The same models under each CV method --------------------------------------
cv_methods <- c("K-Folds", "Stratified_K-Folds", "Repeated_K-Folds", "Hold_Out", "Stratified_Hold_Out")
cv_models <- c("BayesB", "GBLUP_BRR", "RKHS")

for (cv_method in cv_methods) {
  label <- paste("CV:", cv_method, "-", paste(cv_models, collapse = " + "))
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,
      gmatrix = NULL,
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
      var_cov_str = NULL,

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,
      weights = NULL,

      # Candidate models compared by CV
      GS_model_cv = cv_models,

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",
      kernel_method = NULL,
      bending = TRUE,

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = 0.10,
      impute = TRUE,
      imputation_method = "knn",
      ploidy = "auto",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # Cross-validation
      cross_validation = TRUE,
      cross_validation_meth = cv_method,
      sampling_method = NULL,        # implied by the method name (stratified / unstratified)
      nfolds = 3,
      replication = 2,
      test_size = 0.30,              # hold-out methods only
      random_state = 20260903L,
      cv_evaluation_only = TRUE,     # CV metrics only, no refit
      cv_generate_plots = FALSE,

      # Evaluation and ranking
      eval_metrics = c(
        "accuracy",                  # correlation for a Gaussian response
        "root_mean_squared_error",
        "mean_absolute_error",
        "bias"
      ),
      metric_for_ranking = "auto",

      # Bayesian model controls (short chains for this example)
      nIter = 600L,                  # use 16000 for real data
      burnIn = 200L,                 # use 2000
      thin = 5L,

      # Scaling
      scaling = TRUE,
      centering = FALSE,

      # Memory and scheduling
      num_cores = NULL,
      parallel_mode = "auto",
      parallel_backend_prefer_fork = FALSE,

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 2. Several model families compared, then the best model predicts the NA lines --
# Every one-environment model can be compared in one CV call (see scenario 1).
# NeuralAdditive is left out here because it is slow (one network per marker), and
# TabTransformer because it is unstable on very small panels (see scenario 1).
all_single_env_models <- c(
  "BayesA", "BayesB", "BayesC", "BL", "BRR", "GBLUP_BRR", "RKHS",   # Bayesian (BGLR)
  "GBLUP",                                                          # ASReml-R
  "Lasso", "Ridge_Regression", "PartialLeastSquare", "SupportVectorMachine",
  "K-NearestNeighbors", "RandomForest", "Xgboost", "CatBoost", "LightGBM",   # machine learning
  setdiff(dl_models, c("NeuralAdditive", "TabTransformer")),                             # deep learning
  "Kernel-GBLUP", "Gaussian-Process-GBLUP"                          # Gaussian process
)
# keep the models whose software is installed on this machine
family_models <- Filter(function(m) is.null(missing_software(needs_for_model(m))), all_single_env_models)
label <- sprintf("CV + refit best: %d models", length(family_models))
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,
    gmatrix = NULL,
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
    var_cov_str = NULL,

    # Mixed-model formulas
    fixed = ~1,
    random = ~GID,
    weights = NULL,

    # Candidate models compared by CV
    GS_model_cv = family_models,

    # ASReml settings used by GBLUP
    engine = "asreml",
    workspace = 1e8,
    pworkspace = 1e6,
    maxit = 100L,

    # Relationship / kernel construction
    gmatrix_method = "VanRaden",
    kernel_method = NULL,
    bending = TRUE,

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    maf_threshold = 0.01,
    het_threshold = 0.10,
    impute = TRUE,
    imputation_method = "knn",
    ploidy = "auto",
    ld_prunning_qc = TRUE,
    ld_pruning = FALSE,

    # Cross-validation, then refit of the best model
    cross_validation = TRUE,
    cross_validation_meth = "Stratified_K-Folds",
    sampling_method = "stratified",
    nfolds = 3,
    replication = 1,
    random_state = 20260903L,
    cv_evaluation_only = FALSE,      # refit the best model and predict the NA lines
    cv_generate_plots = FALSE,

    # Evaluation and ranking
    eval_metrics = c("accuracy", "root_mean_squared_error", "mean_squared_error",
                     "mean_absolute_error", "bias", "kendalls_tau"),
    metric_for_ranking = "accuracy",
    ranking_tie_breakers = NULL,

    # Machine-learning settings
    para_tunning = FALSE,
    ntree = 50L,                     # use ~500-1000
    iteration = 30L,
    catboost_iterations = 30L,
    lightgbm_nrounds = 30L,
    ncomp = 3L,
    xgb_nthread = 1L,
    rf_n_jobs = 1L,
    n_bootstrap = 10,

    # Bayesian model controls
    nIter = 600L,
    burnIn = 200L,
    thin = 5L,

    # GP controls
    gp_backend = "auto",

    # Deep-learning controls
    epochs = 100L,                     # attention models (TabNet, TabAttention) need ~100
    batch_size = 16L,
    deterministic = TRUE,
    random_seed = 20260903L,
    dl_internal_calibration = FALSE,

    # Scaling
    scaling = TRUE,
    centering = FALSE,

    # Memory and scheduling
    num_cores = NULL,
    parallel_mode = "auto",
    parallel_backend_prefer_fork = FALSE,

    # Output
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  ), silent = TRUE)
})
log_result(label, res, model_runtime)
if (!inherits(res, "try-error")) {
  cat("\nModel ranking by accuracy:\n")
  print(res$cv_results_processed$best_models_list$accuracy)
  cat("\nPredictions of the refitted best model:\n")
  print(head(res$model_results$predicted_values))
}

scenario_report()
