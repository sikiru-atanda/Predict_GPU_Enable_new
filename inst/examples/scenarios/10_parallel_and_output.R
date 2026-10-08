# =============================================================================
# Scenario 10: parallel execution and saving results to a folder
# =============================================================================
# parallel_mode: "auto" (default; PredictProR scores the backends for your
#   machine and data), "sequential", "future", "base_parallel", "foreach",
#   "mirai". num_cores caps the workers. Every worker loads the same
#   PredictProR as your session.
# Results: system_database = FALSE writes one folder per run under
#   options(PredictProR.output_dir = "...") or PREDICTPRO_OUTPUT_DIR (default:
#   the working directory); res$export_directory gives the folder.
#   Run_metadata.csv records the parallel decision (policy_decision_reason;
#   see gp_policy_decision_reason_glossary()).
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_single_env(n = 80, p = 300)
pheno <- dat$pheno
Geno.data <- dat$geno
response <- "Yield"

# ---- 1. The same cross-validation with each parallel backend ---------------------
for (mode in c("auto", "sequential", "base_parallel", "future", "mirai")) {
  label <- paste("Parallel CV:", mode)
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
      num_cores = 2L,
      globals_max_GB = 4,
      memory_budget_gb = NULL,            # NULL: from the free memory
      parallel_mode = mode,
      parallel_backend_prefer_fork = FALSE,
      sequential_models = NULL,           # e.g. c("RKHS") to force models sequential

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
  if (!inherits(res, "try-error") && !is.null(res$run_metadata)) {
    md <- res$run_metadata
    reason <- tryCatch(md$value[md$key == "policy_decision_reason"][1], error = function(e) NA)
    if (!is.null(reason) && length(reason) && !is.na(reason)) cat("   parallel decision:", reason, "\n")
  }
}

# ---- 2. Save results to a folder of your choice -----------------------------------
out_dir <- file.path(tempdir(), "my_PredictProR_results")
old <- options(PredictProR.output_dir = out_dir)
wd_before <- getwd()
label <- "Save results to a folder (system_database = FALSE)"
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

    # No cross-validation
    cross_validation = FALSE,

    # Bayesian model controls
    nIter = 600L,
    burnIn = 200L,
    thin = 5L,

    # Memory and scheduling
    parallel_mode = "auto",
    parallel_backend_prefer_fork = FALSE,

    # Output: write CSV tables and plots
    system_database = FALSE,
    plot_extension = "pdf",
    plot_width = 17,
    plot_height = 12,
    plot_units = "in",
    plot_dpi = 300,
    message = FALSE,
    verbose = FALSE
  ), silent = TRUE)
})
options(old)
log_result(label, res, model_runtime)
files <- list.files(out_dir, recursive = TRUE)
cat("\nResults folder:", out_dir, "\n")
print(head(files, 8))
log_check("Results folder holds Predicted_Value.csv", any(grepl("Predicted_Value\\.csv$", files)),
          sprintf("%d files written", length(files)))
log_check("Working directory unchanged after export", identical(normalizePath(getwd()), normalizePath(wd_before)),
          getwd())

cat("\nTuning knobs for the parallel policy:\n")
print(utils::head(gp_parallel_policy_knobs(), 5))

scenario_report()
