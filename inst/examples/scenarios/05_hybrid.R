# =============================================================================
# Scenario 5: hybrid prediction (female x male crosses)
# =============================================================================
# pheno: one row per hybrid with HybridID, Female, Male and the trait; untested
# hybrids have NA. Declare them with gen_name / female_parent / male_parent.
# Genotypes:
#   hybrid_asreml / hybrid_bayes / hybrid_gp : PARENT genotypes (rows = parent IDs)
#   hybrid_ml / hybrid_dl                    : HYBRID genotypes (rows = HybridID),
#                                              e.g. mid-parent (female + male) / 2
# One hybrid_* flag and one model per call.
# Hybrid CV: "Hybrid_Known_Parents", "Hybrid_One_New_Parent", "Hybrid_Both_New_Parents".
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_hybrid(n_female = 8, n_male = 6, p = 300)
pheno <- dat$pheno                    # HybridID, Female, Male, Yield (NA = hybrids to predict)
Geno.parent <- dat$parent_geno        # parents x markers
Geno.hybrid <- dat$hybrid_geno        # hybrids x markers (mid-parent)
response <- "Yield"

# ---- 1. Hybrid true prediction with each model family ----------------------------
hybrid_models <- c(
  lapply(c("GBLUP_BRR", "RKHS"), function(m) list(name = paste("Bayes", m, "(GCA + SCA)"), flag = "bayes", model = m, geno = "parent")),
  list(list(name = "ASReml GBLUP", flag = "asreml", model = "GBLUP", geno = "parent")),
  list(list(name = "Gaussian-Process-GBLUP (R backend)", flag = "gp", model = "Gaussian-Process-GBLUP", geno = "parent")),
  lapply(c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
         function(m) list(name = paste(m, "(mid-parent genotypes)"), flag = "ml", model = m, geno = "hybrid")),
  lapply(c("DenseNeuralNet", "TabTransformer"),
         function(m) list(name = paste(m, "(mid-parent genotypes)"), flag = "dl", model = m, geno = "hybrid"))
)
# Not for hybrids: Bayesian marker models, Kernel-GBLUP / FA-GBLUP / Scalable-GBLUP,
# K-NearestNeighbors, LightGBM and the other deep-learning models.

for (cfg in hybrid_models) {
  label <- paste("Hybrid TP:", cfg$name)
  need <- if (cfg$flag == "gp") NULL else missing_software(needs_for_model(cfg$model))   # hybrid GP has an R backend
  if (!is.null(need)) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = if (cfg$geno == "parent") Geno.parent else Geno.hybrid,
      gmatrix = NULL,

      # Trait and hybrid identifiers
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      response = response,
      response_family = "gaussian",
      ploidy = 2L,

      # Hybrid mode (exactly one hybrid_* flag)
      multi_trait_gp = FALSE,
      multi_trait_bayes = FALSE,
      multi_trait_asreml = FALSE,
      multi_trait_ml = FALSE,
      multi_trait_dl = FALSE,
      met_ml_dl = FALSE,
      hybrid_asreml = cfg$flag == "asreml",
      hybrid_bayes = cfg$flag == "bayes",
      hybrid_gp = cfg$flag == "gp",
      hybrid_ml = cfg$flag == "ml",
      hybrid_dl = cfg$flag == "dl",
      hybrid_include_sca = TRUE,             # specific combining ability in addition to GCA
      heter_groups = NULL,

      # Formulas: hybrid routes build their own GCA / SCA model
      fixed = NULL,
      random = NULL,

      # Candidate model
      GS_model = cfg$model,

      # ASReml settings
      engine = if (cfg$flag == "asreml") "asreml" else NULL,
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship construction (parent genotypes)
      gmatrix_method = if (cfg$geno == "parent") "VanRaden" else NULL,
      bending = TRUE,

      # Genotype QC and preprocessing. Hybrid genotypes are heterozygous by design:
      # use het_threshold = NULL for hybrid-level genotypes (hybrid_ml / hybrid_dl),
      # otherwise the inbred-line heterozygosity filter removes most markers.
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = if (cfg$geno == "hybrid") NULL else 0.10,
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

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # GP controls (hybrid GP)
      gp_backend = if (cfg$flag == "gp") "r" else "auto",
      hybrid_gp_lambda = "auto",

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

# ---- 2. Hybrid cross-validation: how well are new crosses predicted? -------------
for (cv_method in c("Hybrid_Known_Parents", "Hybrid_One_New_Parent", "Hybrid_Both_New_Parents")) {
  label <- paste("Hybrid CV:", cv_method, "- GBLUP_BRR")
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = Geno.parent,

      # Trait and hybrid identifiers
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      response = response,
      response_family = "gaussian",
      ploidy = 2L,

      # Hybrid mode
      hybrid_bayes = TRUE,
      hybrid_include_sca = TRUE,

      # Formulas: built by the hybrid route
      fixed = NULL,
      random = NULL,

      # Candidate model
      GS_model_cv = "GBLUP_BRR",

      # Relationship construction
      gmatrix_method = "VanRaden",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # Cross-validation by parent novelty
      cross_validation = TRUE,
      cross_validation_meth = cv_method,
      nfolds = 2,
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
}

scenario_report()
