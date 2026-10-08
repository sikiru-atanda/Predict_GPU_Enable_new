# =============================================================================
# Scenario 12: polyploid crops (e.g. tetraploid potato, hexaploid sweet potato)
# =============================================================================
# PredictProR is ploidy aware. Genotypes are ALT-allele dosages 0..ploidy
# (tetraploid: 0, 1, 2, 3, 4). Tell the package the ploidy:
#   * numeric matrix : ploidy = 4L (or attr(Geno.data, "ploidy") <- 4L with ploidy = "auto")
#   * VCF file       : ploidy = "auto" reads it from GT calls such as 0/0/1/1
#   * HapMap         : every allele copy written out (AAAA, AAAG, AAGG, AGGG, GGGG), ploidy = 4L
#   * CSV / TXT      : dosage 0..ploidy, ploidy = 4L
# Ploidy-aware parts: QC (MAF = mean dosage / ploidy), imputation (KNN / mean
# within 0..ploidy), VanRaden / Weighted_VanRaden / Epistasis relationship
# matrices and the generic kernels.
# Important for polyploids:
#   * het_threshold = NULL. Outbred polyploids are mostly heterozygous; the
#     default het_threshold = 0.10 is an inbred-line filter and removes them.
#   * gmatrix_method = "VanRaden". "Yang" and the dominance matrices are
#     diploid-only and stop with a message.
#   * Beagle imputation is diploid-only: use imputation_method = "knn" or "mean".
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data -----------------------------------------------------------------------
dat <- sim_polyploid(n = 100, p = 400, ploidy = 4L)
pheno <- dat$pheno            # GID, Yield (NA = lines to predict)
Geno.data <- dat$geno         # lines x markers, dosage 0..4
table(Geno.data)              # all five tetraploid dosage classes
response <- "Yield"

# ---- 1. Tetraploid, one environment: every compatible model ---------------------
tetraploid_models <- c(
  "BayesA", "BayesB", "BayesC", "BL", "BRR",                  # Bayesian marker regression
  "GBLUP_BRR", "RKHS",                                        # Bayesian kernel models
  "GBLUP",                                                    # ASReml-R
  "Lasso", "Ridge_Regression", "PartialLeastSquare", "SupportVectorMachine",
  "K-NearestNeighbors", "RandomForest", "Xgboost", "CatBoost", "LightGBM",   # machine learning
  setdiff(dl_models, c("NeuralAdditive", "TabTransformer", "GPNet")),              # deep learning
  "Kernel-GBLUP", "Gaussian-Process-GBLUP"                    # Gaussian process
)
# Not here: NeuralAdditive (works, but slow: one network per marker); TabTransformer
# (works, but unstable and less accurate than GBLUP on very small panels). GPNet is not supported for polyploids: it stops with
# a message (its deep-kernel features collapse on polyploid dosage).

for (m in tetraploid_models) {
  label <- paste("Tetraploid TP:", m)
  if (!is.null(need <- missing_software(needs_for_model(m)))) { log_skip(label, need); next }

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
      ploidy = 4L,                       # tetraploid dosage 0..4

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = m,

      # ASReml settings used by GBLUP
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship / kernel construction (ploidy-aware VanRaden)
      gmatrix_method = "VanRaden",
      bending = TRUE,

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = NULL,              # outbred polyploid: no heterozygosity filter
      ind_call_rate_threshold = 0.90,
      snp_call_rate_threshold = 0.90,
      impute = TRUE,
      imputation_method = "knn",         # Beagle is diploid-only
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings (small for this example)
      para_tunning = FALSE,
      ntree = 50L,
      iteration = 30L,
      catboost_iterations = 30L,
      lightgbm_nrounds = 30L,
      n_bootstrap = 10,

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # GP controls
      gp_backend = "auto",

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

# ---- 2. Tetraploid genotype files: VCF, HapMap, CSV ------------------------------
folder <- file.path(tempdir(), "polyploid_files")
dir.create(folder, showWarnings = FALSE)
dosage <- t(unclass(Geno.data)); attr(dosage, "ploidy") <- NULL     # markers x samples
set.seed(5)
dosage[stats::runif(length(dosage)) < 0.02] <- NA                    # 2% missing calls

# VCF: GT with one allele per chromosome copy, e.g. 0/0/1/1 = dosage 2
gt_call <- function(d) vapply(d, function(k) if (is.na(k)) "./././." else
  paste(c(rep("0", 4 - k), rep("1", k)), collapse = "/"), "")
writeLines(c("##fileformat=VCFv4.3",
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(dosage)),
                   collapse = "\t"),
             vapply(seq_len(nrow(dosage)), function(i)
               paste(c("1", i * 100, rownames(dosage)[i], "A", "G", ".", "PASS", ".", "GT", gt_call(dosage[i, ])),
                     collapse = "\t"), "")),
           file.path(folder, "tetraploid.vcf"))

# HapMap: every allele copy written out, AAGG = dosage 2 (IUPAC codes such as R
# are ambiguous above diploid and are refused)
hmp_call <- function(d) vapply(d, function(k) if (is.na(k)) "NNNN" else
  paste(c(rep("A", 4 - k), rep("G", k)), collapse = ""), "")
hapmap_table <- data.frame(`rs#` = rownames(dosage), alleles = "A/G", chrom = "1",
                           pos = seq_len(nrow(dosage)) * 100, strand = "+", `assembly#` = NA, center = NA,
                           protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
                           check.names = FALSE, stringsAsFactors = FALSE)
hapmap_table <- cbind(hapmap_table, as.data.frame(apply(dosage, 2, hmp_call), stringsAsFactors = FALSE))
utils::write.table(hapmap_table, file.path(folder, "tetraploid.hmp.txt"), sep = "\t", quote = FALSE,
                   row.names = FALSE)

# CSV: VCF-style marker columns, then one dosage column per line (0..4, "." = missing)
geno_table <- data.frame(CHROM = 1, POS = seq_len(nrow(dosage)) * 100, ID = rownames(dosage), REF = "A",
                         ALT = "G", QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT", dosage,
                         check.names = FALSE)
utils::write.csv(geno_table, file.path(folder, "tetraploid.csv"), row.names = FALSE, na = ".")

file_inputs <- list(
  list(name = "VCF, ploidy read from GT calls", vcf = "tetraploid.vcf", ploidy = "auto"),
  list(name = "HapMap file", hapmap_file = "tetraploid.hmp.txt", ploidy = 4L),
  list(name = "HapMap table in memory", hapmap = hapmap_table, ploidy = 4L),
  list(name = "CSV dosage file", csv = "tetraploid.csv", ploidy = 4L)
)
for (cfg in file_inputs) {
  label <- paste("Tetraploid file input:", cfg$name, "- GBLUP_BRR")
  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data (one genotype source)
      pheno_data = pheno,
      vcf_file_name = cfg[["vcf"]],
      vcf_file_path = if (!is.null(cfg[["vcf"]])) folder else NULL,
      hapmap = cfg[["hapmap"]],
      hapmap_file_name = cfg[["hapmap_file"]],
      hapmap_file_path = if (!is.null(cfg[["hapmap_file"]])) folder else NULL,
      csv_file_name = cfg[["csv"]],
      csv_file_path = if (!is.null(cfg[["csv"]])) folder else NULL,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",
      ploidy = cfg$ploidy,

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = "GBLUP_BRR",

      # Relationship construction
      gmatrix_method = "VanRaden",

      # Genotype QC and imputation of the missing calls
      qc_filtering = TRUE,
      maf_threshold = 0.01,
      het_threshold = NULL,
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

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 3. Tetraploid cross-validation ---------------------------------------------
label <- "Tetraploid CV: K-Folds - BayesB + GBLUP_BRR + RKHS"
model_runtime <- system.time({
  res <- try(PredictProR::model_execute(
    # Input data
    pheno_data = pheno,
    geno_data = Geno.data,

    # Trait and genotype identifiers
    gen_name = "GID",
    response = response,
    response_family = "gaussian",
    ploidy = 4L,

    # Mixed-model formulas
    fixed = ~1,
    random = ~GID,

    # Candidate models
    GS_model_cv = c("BayesB", "GBLUP_BRR", "RKHS"),

    # Relationship construction
    gmatrix_method = "VanRaden",

    # Genotype QC and preprocessing
    qc_filtering = TRUE,
    het_threshold = NULL,
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
    parallel_mode = "auto",
    parallel_backend_prefer_fork = FALSE,

    # Output
    system_database = TRUE,
    message = FALSE,
    verbose = FALSE
  ), silent = TRUE)
})
log_result(label, res, model_runtime)

# ---- 4. Tetraploid multi-environment trials --------------------------------------
met <- sim_polyploid(n = 60, p = 300, ploidy = 4L, envs = c("E1", "E2", "E3"), seed = 63)
met_ml_dl_models <- c("RandomForest", "Xgboost", "CatBoost", "LightGBM", "DenseNeuralNet", "TabTransformer",
                      "TabAttention", "TabNet", "MixtureOfExperts")
met_models <- c(
  lapply(c("GBLUP_BRR", "RKHS", "Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", met_ml_dl_models),
         function(m) list(name = m, model = m)),
  list(list(name = "GBLUP (ASReml fa1)", model = "GBLUP", vcs = "fa1"))
)
for (cfg in met_models) {
  label <- paste("Tetraploid MET TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }
  is_asreml <- !is.null(cfg[["vcs"]])

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = met$pheno,
      geno_data = met$geno,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",
      ploidy = 4L,

      # Multi-environment mode
      heter_groups = "Env",
      met_ml_dl = cfg$model %in% met_ml_dl_models,
      heter_resid = is_asreml,
      var_cov_str = cfg[["vcs"]],

      # Mixed-model formulas
      fixed = ~Env,
      random = ~GID + GID:Env,

      # Candidate model
      GS_model = cfg$model,

      # ASReml settings
      engine = if (is_asreml) "asreml" else NULL,
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship construction (also used for the ML/DL kernel features)
      gmatrix_method = "VanRaden",
      bending = TRUE,

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      het_threshold = NULL,
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
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

      # GP controls
      gp_backend = "auto",

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

# ---- 5. Tetraploid, two traits analysed jointly ---------------------------------
mt <- pheno
set.seed(9)
mt$Starch <- 0.5 * mt$Yield + stats::rnorm(nrow(mt))       # correlated second trait
mt_models <- list(
  list(name = "joint Bayesian GBLUP_BRR", flag = "bayes", model = "GBLUP_BRR"),
  list(name = "joint ASReml GBLUP (us)", flag = "asreml", model = "GBLUP"),
  list(name = "joint GP Gaussian-Process-GBLUP", flag = "gp", model = "Gaussian-Process-GBLUP"),
  list(name = "joint ML RandomForest", flag = "ml", model = "RandomForest"),
  list(name = "joint DL DenseNeuralNet", flag = "dl", model = "DenseNeuralNet")
)
for (cfg in mt_models) {
  label <- paste("Tetraploid multi-trait TP:", cfg$name)
  if (!is.null(need <- missing_software(needs_for_model(cfg$model)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = mt,
      geno_data = Geno.data,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = c("Yield", "Starch"),
      response_family = "gaussian",
      ploidy = 4L,

      # Joint multi-trait mode (exactly one multi_trait_* flag)
      multi_trait_bayes = cfg$flag == "bayes",
      multi_trait_asreml = cfg$flag == "asreml",
      multi_trait_gp = cfg$flag == "gp",
      multi_trait_ml = cfg$flag == "ml",
      multi_trait_dl = cfg$flag == "dl",
      var_cov_str = if (cfg$flag == "asreml") "us" else NULL,

      # Formulas: the joint multi-trait routes build their own model
      fixed = NULL,
      random = NULL,

      # Candidate model
      GS_model = cfg$model,

      # ASReml settings
      engine = if (cfg$flag == "asreml") "asreml" else NULL,
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship construction
      gmatrix_method = "VanRaden",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      het_threshold = NULL,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings
      ntree = 50L,
      n_bootstrap = 10,

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

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 6. Hexaploid (e.g. sweet potato): dosage 0..6 --------------------------------
hex <- sim_polyploid(n = 100, p = 400, ploidy = 6L, seed = 71)
for (m in c("BayesB", "GBLUP_BRR", "RKHS", "GBLUP", "Kernel-GBLUP", "RandomForest")) {
  label <- paste("Hexaploid TP:", m)
  if (!is.null(need <- missing_software(needs_for_model(m)))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data
      pheno_data = hex$pheno,
      geno_data = hex$geno,              # carries attr(, "ploidy") = 6

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",
      ploidy = "auto",                   # read from the matrix's ploidy attribute

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = m,

      # ASReml settings used by GBLUP
      engine = "asreml",
      workspace = 1e8,
      pworkspace = 1e6,
      maxit = 100L,

      # Relationship construction
      gmatrix_method = "VanRaden",

      # Genotype QC and preprocessing
      qc_filtering = TRUE,
      het_threshold = NULL,
      impute = TRUE,
      imputation_method = "knn",
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Machine-learning settings
      ntree = 50L,
      n_bootstrap = 10,

      # Bayesian model controls
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 7. Common polyploid mistakes and the messages you get -------------------------
expect_stop <- function(label, expr, pattern) {
  res <- try(expr, silent = TRUE)
  msg <- if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else ""
  ok <- inherits(res, "try-error") && grepl(pattern, msg)
  log_check(label, ok, if (ok) substr(trimws(msg), 1, 110) else paste("unexpected:", substr(msg, 1, 100)))
}
quick <- function(...) PredictProR::model_execute(pheno_data = pheno, gen_name = "GID", response = response,
                                                  fixed = ~1, random = ~GID, GS_model = "GBLUP_BRR",
                                                  nIter = 300L, burnIn = 100L, message = FALSE, verbose = FALSE, ...)
no_attr <- Geno.data
attr(no_attr, "ploidy") <- NULL
expect_stop("Mistake: dosage 0..4 without ploidy", quick(geno_data = no_attr, ploidy = "auto",
            gmatrix_method = "VanRaden", het_threshold = NULL), "Supply ploidy")
# het_threshold = 0.10 drops markers heterozygous (0 < dosage < 4) in more than 10%
# of lines - in an outbred polyploid that is nearly every marker:
het <- colMeans(Geno.data > 0 & Geno.data < 4)
kept <- sum(het <= 0.10)
cat(sprintf("het_threshold = 0.10 would keep %d of %d tetraploid markers\n", kept, length(het)))
log_check("Mistake: inbred heterozygosity filter", kept < 0.05 * length(het),
          sprintf("would keep only %d of %d markers: use het_threshold = NULL", kept, length(het)))
expect_stop("Mistake: diploid-only Yang GRM", quick(geno_data = Geno.data, ploidy = 4L,
            gmatrix_method = "Yang", het_threshold = NULL), "diploid-only")
expect_stop("Mistake: GPNet on tetraploid data",
            PredictProR::model_execute(pheno_data = pheno, geno_data = Geno.data, gen_name = "GID",
                                       response = response, GS_model = "GPNet", ploidy = 4L,
                                       het_threshold = NULL, message = FALSE, verbose = FALSE),
            "GPNet is not supported for polyploid")
expect_stop("Mistake: Beagle on tetraploid VCF", quick(vcf_file_name = "tetraploid.vcf", vcf_file_path = folder,
            ploidy = 4L, gmatrix_method = "VanRaden", het_threshold = NULL, impute = TRUE,
            imputation_method = "beagle"), "diploid")

scenario_report()
