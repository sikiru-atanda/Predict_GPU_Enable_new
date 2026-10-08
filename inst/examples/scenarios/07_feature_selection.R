# =============================================================================
# Scenario 7: feature (marker / omics) selection
# =============================================================================
#  (a) Standalone: feature_score_predictors() ranks predictors on the training
#      lines; feature_select_top_k() keeps the top k; fit any model on them.
#  (b) Inside model_execute(): feature_scoring = TRUE with
#        feature_scoring_model = "Ridge_Regression" | "BayesB" | "RandomForest"
#        feature_k      = 50          (true prediction: keep the top 50)
#        feature_k_grid = c(25, 100)  (CV: compare several k; "all" is added)
#        feature_scoring_cv = "fold_internal" (score inside each training fold,
#                                              no information leak) | "fixed" | "both"
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_single_env(n = 100, p = 400)
pheno <- dat$pheno
Geno.data <- dat$geno
response <- "Yield"
train <- pheno[!is.na(pheno$Yield), ]

# ---- (a) Standalone scoring on the training lines, then BayesB on the top 50 -------
label <- "Feature selection: standalone Ridge scores + top 50 + BayesB"
scores <- feature_score_predictors(
  predictor_data = Geno.data[train$GID, ],
  pheno_data = train,
  response = response,
  gen_name = "GID",
  scoring_model = "Ridge_Regression",
  seed = 1L
)
print(head(scores[order(scores$rank), c("predictor", "normalized_score", "rank")], 5))
Geno.selected <- feature_select_top_k(
  predictor_data = Geno.data,
  feature_score_metadata = scores,
  trait = response,
  k = 50
)

model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data (the 50 selected markers)
    pheno_data = pheno,
    geno_data = Geno.selected,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    response_family = "gaussian",

    # Mixed-model formulas
    fixed = ~1,
    random = ~GID,

    # Candidate model
    GS_model = "BayesB",

    # Genotype QC and preprocessing (already selected: keep all 50)
    qc_filtering = FALSE,
    impute = FALSE,
    ld_prunning_qc = FALSE,

    # No cross-validation
    cross_validation = FALSE,

    # Bayesian model controls
    nIter = 600L,                  # use 16000
    burnIn = 200L,                 # use 2000
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

# ---- (b1) Integrated: CV over several k (scores computed inside each fold) ---------
# Scoring models: Ridge_Regression, BayesB, RandomForest. Selection works with any model
# fitted on raw markers (Bayesian, GBLUP, machine learning, deep learning, GP).
fs_candidates <- Filter(function(m) is.null(missing_software(needs_for_model(m))),
                        c("BayesB", "GBLUP_BRR", "Ridge_Regression", "RandomForest", "Kernel-GBLUP"))
for (scoring_model in c("Ridge_Regression", "BayesB", "RandomForest")) {
  if (!is.null(need <- missing_software(needs_for_model(scoring_model)))) {
    log_skip(paste("Feature selection inside CV - scored by", scoring_model), need); next
  }
  label <- paste("Feature selection inside CV: k = 25, 100, all - scored by", scoring_model)
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.data,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate models
      GS_model_cv = fs_candidates,

      # Relationship / kernel construction (kernels are rebuilt from the selected markers)
      gmatrix_method = "VanRaden",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # Feature selection
      feature_scoring = TRUE,
      feature_scoring_model = scoring_model,
      feature_scoring_cv = "fold_internal",
      feature_k_grid = c(25L, 100L),
      feature_scoring_seed = 1L,

      # Cross-validation
      cross_validation = TRUE,
      cross_validation_meth = "K-Folds",
      nfolds = 3,
      replication = 1,
      random_state = 20260903L,
      cv_evaluation_only = TRUE,
      cv_generate_plots = FALSE,
      eval_metrics = c("accuracy", "root_mean_squared_error"),

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
  }

# ---- (b2) Integrated: true prediction with the top 50 markers ------------------------
label <- "Feature selection for true prediction: top 50 - BayesB"
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    response_family = "gaussian",

    # Mixed-model formulas
    fixed = ~1,
    random = ~GID,

    # Candidate model
    GS_model = "BayesB",

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    impute = TRUE,
    imputation_method = "knn",
    ld_prunning_qc = TRUE,
    ld_pruning = FALSE,

    # Feature selection
    feature_scoring = TRUE,
    feature_scoring_model = "Ridge_Regression",
    feature_scoring_cv = "fixed",
    feature_k = 50L,
    feature_scoring_seed = 1L,

    # No cross-validation
    cross_validation = FALSE,

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

scenario_report()
