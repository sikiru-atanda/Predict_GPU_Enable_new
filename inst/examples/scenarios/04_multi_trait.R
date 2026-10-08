# =============================================================================
# Scenario 4: several traits (multi-trait), one and several environments
# =============================================================================
# response = c("Yield", "Protein").
#   No multi-trait flag      : each trait fitted on its own (one GS_model per trait)
#   multi_trait_bayes = TRUE : GBLUP_BRR, RKHS (joint model, genetic covariance)
#   multi_trait_asreml = TRUE: GBLUP (engine = "asreml"), var_cov_str "us" / "corgh" / "diag"
#   multi_trait_gp = TRUE    : Gaussian-Process-GBLUP, FA-GBLUP, Scalable-GBLUP
#   multi_trait_ml = TRUE    : RandomForest, Ridge_Regression, PartialLeastSquare, Lasso,
#                              SupportVectorMachine, Xgboost, CatBoost
#   multi_trait_dl = TRUE    : DenseNeuralNet, TabTransformer
#   multi_trait = TRUE       : compare joint models by CV and refit the best
# Traits may have different NA lines (each trait predicts its own NA lines).
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_multi_trait(n = 80, p = 300)
pheno <- dat$pheno            # GID, Yield, Protein (NA = lines to predict)
Geno.data <- dat$geno
response <- c("Yield", "Protein")

# ---- 1. Separate and joint multi-trait models ----------------------------------------
mt_models <- c(
  list(list(name = "separate fits (GBLUP_BRR per trait)", flag = "none", model = c("GBLUP_BRR", "GBLUP_BRR"))),
  lapply(c("GBLUP_BRR", "RKHS"), function(m) list(name = paste("joint Bayesian", m), flag = "bayes", model = m)),
  lapply(c("us", "corgh", "diag"), function(v) list(name = paste0("joint ASReml GBLUP (", v, ")"), flag = "asreml",
                                                   model = "GBLUP", vcs = v)),
  lapply(c("Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP"), function(m) list(name = paste("joint GP", m), flag = "gp", model = m)),
  lapply(c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
         function(m) list(name = paste("joint ML", m), flag = "ml", model = m)),
  lapply(c("DenseNeuralNet", "TabTransformer"), function(m) list(name = paste("joint DL", m), flag = "dl", model = m))
)
# Not for joint multi-trait: Bayesian marker models, Kernel-GBLUP, K-NearestNeighbors,
# LightGBM and the other deep-learning models (use them trait by trait, without a flag).

for (cfg in mt_models) {
  label <- paste("Multi-trait TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model[1])))) { log_skip(label, need); next }

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

      # Multi-trait mode (exactly one flag for a joint model)
      multi_trait_gp = cfg$flag == "gp",
      multi_trait_bayes = cfg$flag == "bayes",
      multi_trait_asreml = cfg$flag == "asreml",
      multi_trait_ml = cfg$flag == "ml",
      multi_trait_dl = cfg$flag == "dl",
      met_ml_dl = FALSE,
      hybrid_asreml = FALSE,
      hybrid_bayes = FALSE,
      hybrid_gp = FALSE,
      hybrid_ml = FALSE,
      hybrid_dl = FALSE,
      heter_groups = NULL,
      heter_resid = FALSE,
      var_cov_str = if (cfg$flag == "asreml") cfg[["vcs"]] else NULL,   # "us", "corgh" or "diag"

      # Mixed-model formulas (separate fits need random; joint models set their own;
      # fixed = ~1 is the default intercept and is accepted by every route)
      fixed = ~1,
      random = if (cfg$flag == "none") ~GID else NULL,

      # Candidate model(s)
      GS_model = cfg$model,

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
      internal_cv_nfolds = 3L,

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # GP controls
      gp_backend = "auto",
      gp_return_se = TRUE,
      gp_return_trait_correlations = TRUE,

      # Deep-learning controls
      epochs = 100L,
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

# ---- 2. Multi-trait cross-validation (joint Bayesian model) --------------------------
label <- "Multi-trait CV: joint GBLUP_BRR, K-Folds"
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    response_family = "gaussian",

    # Multi-trait mode
    multi_trait_bayes = TRUE,

    # Candidate model
    GS_model_cv = "GBLUP_BRR",

    # Relationship / kernel construction
    gmatrix_method = "VanRaden",

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    impute = TRUE,
    imputation_method = "knn",
    ld_prunning_qc = TRUE,
    ld_pruning = FALSE,

    # Cross-validation (genotype-blocked K-Folds)
    cross_validation = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3,
    replication = 1,
    random_state = 20260903L,
    cv_evaluation_only = TRUE,
    cv_generate_plots = FALSE,

    # Bayesian model controls
    nIter = 600L,
    burnIn = 200L,
    thin = 5L,

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

# ---- 3. Automatic: compare joint models by CV and refit the best ----------------------
label <- "Multi-trait auto: GBLUP_BRR vs RKHS"
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    response_family = "gaussian",

    # Automatic multi-trait orchestration (leave the five multi_trait_* flags FALSE)
    multi_trait = TRUE,

    # Candidate joint models
    GS_model_cv = c("GBLUP_BRR", "RKHS"),

    # Relationship / kernel construction
    gmatrix_method = "VanRaden",

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    impute = TRUE,
    imputation_method = "knn",
    ld_prunning_qc = TRUE,
    ld_pruning = FALSE,

    # Cross-validation, then refit the best model
    cross_validation = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 3,
    replication = 1,
    random_state = 20260903L,
    cv_evaluation_only = FALSE,
    cv_generate_plots = FALSE,

    # Bayesian model controls
    nIter = 600L,
    burnIn = 200L,
    thin = 5L,

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

# ---- 4. Multi-trait, multi-environment ---------------------------------------------
met <- sim_multi_trait(n = 60, p = 300, envs = c("E1", "E2", "E3"))
pheno_met <- met$pheno
Geno_met <- met$geno

mt_met_models <- c(
  lapply(c("GBLUP_BRR", "RKHS"), function(m) list(name = paste("joint Bayesian", m), flag = "bayes", model = m)),
  lapply(c("Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP"), function(m) list(name = paste("joint GP", m), flag = "gp", model = m)),
  list(list(name = "joint ASReml GBLUP (fa1)", flag = "asreml", model = "GBLUP", vcs = "fa1"))
)
for (cfg in mt_met_models) {
  label <- paste("Multi-trait MET TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno_met,
      geno_data = Geno_met,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Multi-trait, multi-environment mode
      multi_trait_gp = cfg$flag == "gp",
      multi_trait_bayes = cfg$flag == "bayes",
      multi_trait_asreml = cfg$flag == "asreml",
      heter_groups = "Env",
      heter_resid = cfg$flag == "asreml",   # ASReml MET: environment-specific residuals
      var_cov_str = if (cfg$flag == "asreml") cfg[["vcs"]] else NULL,   # fa1, fa2 or fa3 for MT-MET ASReml

      # ASReml settings
      engine = if (cfg$flag == "asreml") "asreml" else NULL,
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Candidate model
      GS_model = cfg$model,

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Bayesian model controls
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

      # GP controls
      gp_return_se = TRUE,
      gp_return_trait_correlations = TRUE,
      gp_fa_rank = 1L,

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



