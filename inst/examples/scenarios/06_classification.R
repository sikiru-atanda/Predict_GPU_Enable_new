# =============================================================================
# Scenario 6: classification traits (binary, ordinal, multiclass)
# =============================================================================
# Store the trait as a factor (binary / multiclass) or an ordered factor
# (ordinal) and set response_family = "binary", "ordinal" or "multiclass".
#   binary     : Bayesian models, Xgboost, RandomForest, CatBoost, LightGBM,
#                SupportVectorMachine, K-NearestNeighbors, deep learning
#   ordinal    : Bayesian models
#   multiclass : the machine-learning models above and deep learning
# Predictions give the predicted class, its confidence and class probabilities.
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_classification(n = 120, p = 300)
pheno <- dat$pheno            # GID, Rust (no/yes), Score (low<mid<high), Type (A/B/C); NA = to predict
Geno.data <- dat$geno
str(pheno)

# ---- 1. True prediction for each trait type --------------------------------------
bayes_cls_models <- c("BayesA", "BayesB", "BayesC", "BL", "BRR", "GBLUP_BRR", "RKHS")
ml_cls_models <- c("SupportVectorMachine", "K-NearestNeighbors", "RandomForest", "Xgboost", "CatBoost", "LightGBM")
# Deep-learning models (00_setup.R). NeuralAdditive also works but fits one
# small network per marker: very slow with 100 epochs, so it is left out here.
dl_cls_models <- setdiff(dl_models, "NeuralAdditive")
cls_models <- c(
  lapply(c(bayes_cls_models, ml_cls_models, dl_cls_models),
         function(m) list(name = paste("binary Rust -", m), response = "Rust", family = "binary", model = m)),
  lapply(bayes_cls_models,
         function(m) list(name = paste("ordinal Score -", m), response = "Score", family = "ordinal", model = m)),
  lapply(c(ml_cls_models, dl_cls_models),
         function(m) list(name = paste("multiclass Type -", m), response = "Type", family = "multiclass", model = m))
)
# Gaussian traits only: ASReml GBLUP, Lasso, Ridge_Regression, PartialLeastSquare
# and the Gaussian process models.

for (cfg in cls_models) {
  label <- paste("Classification TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,
      gmatrix = NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = cfg$response,
      response_family = cfg$family,          # "binary", "ordinal" or "multiclass"
      positive_class = if (cfg$family == "binary") "yes" else NULL,

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

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = cfg$model,

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

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings
      para_tunning = FALSE,
      ntree = 50L,
      iteration = 30L,
      xgb_nthread = 1L,
      n_bootstrap = 10,
      auto_class_weights = TRUE,

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # Deep-learning controls
      epochs = 100L,                 # attention models (TabNet, TabAttention) need ~100
      batch_size = 16L,
      deterministic = TRUE,
      random_seed = 20260903L,
      dl_internal_calibration = FALSE,

      # Memory and scheduling
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

# ---- 2. Classification cross-validation (stratified folds keep class proportions) --
label <- "Classification CV: binary Rust - BayesB vs GBLUP_BRR"
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = "Rust",
    response_family = "binary",
    positive_class = "yes",

    # Mixed-model formulas
    fixed = ~1,
    random = ~GID,

    # Candidate models
    GS_model_cv = c("BayesB", "GBLUP_BRR"),

    # Relationship / kernel construction
    gmatrix_method = "VanRaden",

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    impute = TRUE,
    imputation_method = "knn",
    ld_prunning_qc = TRUE,
    ld_pruning = FALSE,

    # Cross-validation
    cross_validation = TRUE,
    cross_validation_meth = "Stratified_K-Folds",
    sampling_method = "stratified",
    nfolds = 3,
    replication = 1,
    random_state = 20260903L,
    cv_evaluation_only = TRUE,
    cv_generate_plots = FALSE,

    # Evaluation (NULL = the default classification metrics)
    eval_metrics = NULL,
    metric_for_ranking = "auto",

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

# ---- 3. Multi-environment classification ---------------------------------------
#   binary     : RKHS (fixed = ~Env, random = ~GID + GID:Env) and the 9 MET ML/DL models
#   multiclass : the 9 MET ML/DL models (met_ml_dl = TRUE)
#   ordinal    : RKHS
met <- sim_met(n = 60, p = 300)
pheno_met <- met$pheno
cut3 <- stats::quantile(pheno_met$Yield, c(1 / 3, 2 / 3), na.rm = TRUE)
pheno_met$Disease <- factor(ifelse(pheno_met$Yield > stats::median(pheno_met$Yield, na.rm = TRUE), "yes", "no"),
                            levels = c("no", "yes"))
pheno_met$Grade <- ordered(cut(pheno_met$Yield, c(-Inf, cut3, Inf), labels = c("low", "mid", "high")),
                           levels = c("low", "mid", "high"))
pheno_met$Group <- factor(cut(pheno_met$Yield, c(-Inf, cut3, Inf), labels = c("A", "B", "C")))
met_ml_dl_models <- c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "DenseNeuralNet", "TabTransformer",
                      "TabAttention", "TabNet", "MixtureOfExperts")
met_cls_models <- c(
  lapply(c("RKHS", met_ml_dl_models),
         function(m) list(name = paste("MET binary Disease -", m), response = "Disease", family = "binary", model = m)),
  lapply(met_ml_dl_models,
         function(m) list(name = paste("MET multiclass Group -", m), response = "Group", family = "multiclass", model = m)),
  list(list(name = "MET ordinal Grade - RKHS", response = "Grade", family = "ordinal", model = "RKHS"))
)

for (cfg in met_cls_models) {
  label <- paste("Classification TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }
  is_ml_dl <- cfg$model %in% met_ml_dl_models

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno_met,
      geno_data = met$geno,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = cfg$response,
      response_family = cfg$family,
      positive_class = if (cfg$family == "binary") "yes" else NULL,

      # Multi-environment mode
      heter_groups = "Env",
      met_ml_dl = is_ml_dl,                  # TRUE for machine / deep learning MET

      # Mixed-model formulas (RKHS: environment fixed, genotype and G x E random)
      fixed = ~Env,
      random = ~GID + GID:Env,

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

      # Machine-learning settings
      para_tunning = FALSE,
      ntree = 50L,
      n_bootstrap = 10,
      internal_cv_nfolds = 3L,

      # Bayesian model controls
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

      # Deep-learning controls
      epochs = 100L,
      batch_size = 16L,
      deterministic = TRUE,
      random_seed = 20260903L,
      dl_internal_calibration = FALSE,

      # Memory and scheduling
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

scenario_report()
