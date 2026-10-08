# =============================================================================
# Scenario 8: multi-omics data and multiple kernels
# =============================================================================
# Up to three extra data layers next to the markers:
#   raw layers         : omic1_data, omic2_data, omic3_data (lines in rows, same
#                        row names as geno_data), named with omics_data_label
#   precomputed kernels: omic1_kernel .. omic3_kernel (omics_kernel_label),
#                        gkernel, or a named kernel_list = list(Gaussian = K1, ...)
# Each layer is its own random effect in Bayesian / GP / ASReml models (variance
# per layer is reported); machine / deep learning models concatenate them.
# Kernel models built from raw layers need gmatrix_method and kernel_method.
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_single_env(n = 80, p = 300)
pheno <- dat$pheno
Geno.data <- dat$geno
response <- "Yield"
Transcripts <- sim_omics(rownames(Geno.data), n_features = 120)
Metabolites <- sim_omics(rownames(Geno.data), n_features = 60, seed = 52)
# a second kernel on the same lines (Gaussian kernel of the markers)
D2 <- as.matrix(stats::dist(Geno.data))^2
Gaussian_K <- exp(-D2 / stats::median(D2))
dimnames(Gaussian_K) <- list(rownames(Geno.data), rownames(Geno.data))
COP <- tcrossprod(scale(Transcripts)) / ncol(Transcripts) + diag(1e-4, nrow(Transcripts))   # omics kernel

kernel_models <- c("GBLUP_BRR", "RKHS", "GBLUP", "Kernel-GBLUP", "Gaussian-Process-GBLUP")
marker_models <- c("BayesA", "BayesB", "BayesC", "BL", "BRR")
ml_models <- c("Lasso", "Ridge_Regression", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors",
               "RandomForest", "Xgboost", "CatBoost", "LightGBM")
omics_setups <- c(
  # raw omics layers: one effect per layer (Bayesian marker models); ML / DL concatenate the layers
  # (NeuralAdditive also works but fits one network per feature: very slow with omics;
  #  TabTransformer is unstable on very small panels, see scenario 1)
  lapply(c(marker_models, ml_models, setdiff(dl_models, c("NeuralAdditive", "TabTransformer"))), function(m)
    list(name = paste("markers + 2 raw omics layers -", m), model = m, omics = "raw2", kernel_list = NULL)),
  # raw omics turned into a kernel (kernel_method) next to the genomic relationship
  lapply(kernel_models, function(m)
    list(name = paste("markers + raw omics -> kernel -", m), model = m, omics = "raw1", kernel_list = NULL)),
  # a precomputed omics kernel (e.g. COP) next to the genomic relationship
  lapply(kernel_models, function(m)
    list(name = paste("markers + precomputed omics kernel (COP) -", m), model = m, omics = "kernel", kernel_list = NULL)),
  # several kernels on the same lines (kernel_list). Not ASReml GBLUP here: two
  # kernels from the same markers carry nearly the same information, REML
  # cannot separate their variances and stops with a singular AI matrix.
  lapply(setdiff(kernel_models, "GBLUP"), function(m)
    list(name = paste("genomic + Gaussian kernel (kernel_list) -", m), model = m, omics = "none",
         kernel_list = list(Gaussian = Gaussian_K)))
)
for (cfg in omics_setups) {
  label <- paste("Multi-omics TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }
  uses_gmatrix <- !is.null(cfg$kernel_list)

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = pheno,
      geno_data = if (uses_gmatrix) NULL else Geno.data,
      gmatrix = if (uses_gmatrix) dat$gmatrix else NULL,
      omic1_data = if (cfg$omics %in% c("raw1", "raw2")) Transcripts else NULL,
      omic2_data = if (cfg$omics == "raw2") Metabolites else NULL,
      omics_data_label = list(omic1_data = "Transcripts", omic2_data = "Metabolites", omic3_data = NULL),
      omic1_kernel = if (cfg$omics == "kernel") COP else NULL,
      omics_kernel_label = list(omic1_kernel = if (cfg$omics == "kernel") "COP" else NULL,
                                omic2_kernel = NULL, omic3_kernel = NULL),
      gkernel = NULL,
      kernel_list = cfg$kernel_list,

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

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = cfg$model,

      # ASReml settings (used by GBLUP)
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship / kernel construction
      gmatrix_method = if (uses_gmatrix) NULL else "VanRaden",
      kernel_method = if (cfg$omics == "raw1") "Linear_kernel" else NULL,   # omics -> kernel
      bending = TRUE,

      # Genotype / omics QC and preprocessing
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = 0.10,
      impute = TRUE,
      impute_omic = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings
      para_tunning = FALSE,
      n_bootstrap = 10,

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # GP controls
      gp_backend = "auto",
      gp_return_se = TRUE,
      gp_full_vc = TRUE,           # variance per kernel

      # Machine / deep-learning controls
      ntree = 50L,
      iteration = 30L,
      catboost_iterations = 30L,
      lightgbm_nrounds = 30L,
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

scenario_report()

