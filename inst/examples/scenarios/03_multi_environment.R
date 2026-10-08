# =============================================================================
# Scenario 3: multi-environment trials (MET), one trait
# =============================================================================
# pheno is in LONG format: one row per line x environment, with the
# environment column given as heter_groups = "Env". Lines may be missing from
# some environments (unbalanced). Every line is predicted in every environment
# (Train_Test_Label: Train / Test / Unobserved).
#   Bayesian MET   : GBLUP_BRR, RKHS with random = ~GID + GID:Env
#   ASReml MET     : GBLUP, var_cov_str = "us","corgh","corh","corv","fa1".."fa4","rr1".."rr4"
#   GP MET         : Kernel-GBLUP, Gaussian-Process-GBLUP, FA-GBLUP
#   ML/DL MET      : met_ml_dl = TRUE with RandomForest, Xgboost, CatBoost, LightGBM,
#                    DenseNeuralNet, TabTransformer, ...
#   MET CV         : "CV0" (new environment), "CV1" (new lines), "CV2" (missing cells)
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_met(n = 60, p = 300, envs = c("E1", "E2", "E3"), unbalanced = TRUE)
pheno <- dat$pheno            # GID, Env, Yield (NA = lines to predict), Weight
Geno.data <- dat$geno
response <- "Yield"
table(pheno$Env)

# ---- 1. Bayesian, GP and ML/DL MET models -------------------------------------------
met_models <- list(
  list(name = "GBLUP_BRR", model = "GBLUP_BRR", met_ml_dl = FALSE, hetero = FALSE),
  list(name = "RKHS", model = "RKHS", met_ml_dl = FALSE, hetero = FALSE),
  list(name = "RKHS, environment-specific variances", model = "RKHS", met_ml_dl = FALSE, hetero = TRUE),
  list(name = "Kernel-GBLUP", model = "Kernel-GBLUP", met_ml_dl = FALSE, hetero = FALSE),
  list(name = "Gaussian-Process-GBLUP", model = "Gaussian-Process-GBLUP", met_ml_dl = FALSE, hetero = FALSE),
  list(name = "FA-GBLUP", model = "FA-GBLUP", met_ml_dl = FALSE, hetero = FALSE),
  list(name = "RandomForest", model = "RandomForest", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "Xgboost", model = "Xgboost", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "CatBoost", model = "CatBoost", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "LightGBM", model = "LightGBM", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "DenseNeuralNet", model = "DenseNeuralNet", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "TabTransformer", model = "TabTransformer", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "TabAttention", model = "TabAttention", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "TabNet", model = "TabNet", met_ml_dl = TRUE, hetero = FALSE),
  list(name = "MixtureOfExperts", model = "MixtureOfExperts", met_ml_dl = TRUE, hetero = FALSE)
)
# Not for MET: Bayesian marker models (BayesA/B/C, BL, BRR), Scalable-GBLUP
# (same estimator as Kernel-GBLUP) and the other ML/DL models (single-environment only).

for (cfg in met_models) {
  label <- paste("MET TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,
      gmatrix = NULL,          # constructed from geno_data
      kernel_list = NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Multi-environment mode
      multi_trait_gp = FALSE,
      multi_trait_bayes = FALSE,
      multi_trait_asreml = FALSE,
      multi_trait_ml = FALSE,
      multi_trait_dl = FALSE,
      met_ml_dl = cfg$met_ml_dl,                 # TRUE for machine / deep learning MET
      hybrid_asreml = FALSE,
      hybrid_bayes = FALSE,
      hybrid_gp = FALSE,
      hybrid_ml = FALSE,
      hybrid_dl = FALSE,
      heter_groups = "Env",
      heter_resid = FALSE,
      bayes_kernel_heter_resid = cfg$hetero,     # TRUE: environment-specific variances (RKHS)
      var_cov_str = NULL,
      met_predict_all_environments = TRUE,       # every line in every environment

      # Mixed-model formulas
      fixed = ~Env,
      random = ~GID + GID:Env,
      weights = NULL,

      # Candidate model
      GS_model = cfg$model,

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

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings
      para_tunning = FALSE,
      ntree = 50L,
      n_bootstrap = 10,
      internal_cv_nfolds = 3L,

      # Bayesian model controls (short chains)
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # GP controls
      gp_backend = "auto",
      gp_output_level = "predict_with_se",
      gp_return_se = TRUE,
      gp_fa_rank = 1L,

      # Deep-learning controls
      epochs = 100L,                 # attention models (TabNet, TabAttention) need ~100
      batch_size = 16L,
      deterministic = TRUE,
      random_seed = 20260903L,
      dl_internal_calibration = FALSE,

      # Memory and scheduling
      num_cores = NULL,
      parallel_mode = "auto",
      parallel_backend_prefer_fork = FALSE,
      random_state = 20260903L,

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 2. ASReml MET: genotype x environment variance-covariance structures --------
# Unstructured-type models ("us", "corgh") estimate one genetic correlation per
# pair of environments and need enough lines; this data set has 120.
dat120 <- sim_met(n = 120, p = 300, unbalanced = TRUE, seed = 12)
pheno120 <- dat120$pheno
Geno120 <- dat120$geno

# Structures that fit with 3 environments. fa3 / fa4 / rr3 / rr4 need more
# environments than factors (e.g. fa3 needs 4+); corv = common correlation.
for (vcs in c("us", "corgh", "corh", "corv", "fa1", "fa2", "rr1", "rr2")) {
  label <- paste("MET TP: GBLUP (ASReml)", vcs)
  if (!has_asreml()) { log_skip(label, "asreml"); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno120,
      geno_data = Geno120,
      gmatrix = NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Multi-environment mode
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
      heter_groups = "Env",
      heter_resid = TRUE,                # environment-specific residual variances
      var_cov_str = vcs,                 # genetic variance-covariance between environments
      met_predict_all_environments = TRUE,

      # Mixed-model formulas
      fixed = ~Env,
      random = ~GID + GID:Env,
      weights = NULL,

      # Candidate model
      GS_model = "GBLUP",

      # ASReml settings
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",
      inverse = TRUE,
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

      # No cross-validation
      cross_validation = FALSE,

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

# ---- 3. MET cross-validation scenarios ---------------------------------------------
for (cv_method in c("CV0", "CV1", "CV2", "Repeated_CV0", "Repeated_CV1", "Repeated_CV2")) {
  label <- paste("MET CV:", cv_method, "- GBLUP_BRR + RKHS")
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,
      gmatrix = NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Multi-environment mode
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
      heter_groups = "Env",
      heter_resid = FALSE,
      bayes_kernel_heter_resid = FALSE,
      var_cov_str = NULL,

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID + GID:Env,
      weights = "Weight",

      # Candidate models compared by CV
      GS_model_cv = c("GBLUP_BRR", "RKHS"),

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",
      bending = TRUE,

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = 0.10,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # Cross-validation: CV0 = new environment, CV1 = new lines, CV2 = missing cells
      cross_validation = TRUE,
      cross_validation_meth = cv_method,
      sampling_method = "unstratified",
      nfolds = 3,
      replication = if (grepl("^Repeated", cv_method)) 2 else 1,
      random_state = 20260903L,
      cv_evaluation_only = TRUE,
      cv_generate_plots = FALSE,

      # Evaluation and ranking
      eval_metrics = c("accuracy", "root_mean_squared_error"),
      metric_for_ranking = "accuracy",

      # Bayesian model controls
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

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

# ---- 4. ASReml MET cross-validation (CV1, fa1) --------------------------------------
label <- "MET CV: CV1 - GBLUP (ASReml fa1)"
if (!has_asreml()) log_skip(label, "asreml") else {
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno120,
      geno_data = Geno120,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Multi-environment mode
      heter_groups = "Env",
      heter_resid = TRUE,
      var_cov_str = "fa1",

      # Mixed-model formulas
      fixed = ~Env,
      random = ~GID + GID:Env,

      # Candidate model
      GS_model_cv = "GBLUP",

      # ASReml settings
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",
      bending = TRUE,

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # Cross-validation
      cross_validation = TRUE,
      cross_validation_meth = "CV1",
      sampling_method = "unstratified",
      nfolds = 3,
      replication = 1,
      random_state = 20260903L,
      cv_evaluation_only = TRUE,
      cv_generate_plots = FALSE,
      eval_metrics = c("accuracy", "root_mean_squared_error"),

      # Memory and scheduling
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

scenario_report()
