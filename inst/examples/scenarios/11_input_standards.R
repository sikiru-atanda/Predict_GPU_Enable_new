# =============================================================================
# Scenario 11: the PredictProR input format - print it and check your data
# =============================================================================
#   prediction_input_standard("general" | "hybrid" | "multi_trait" | "multi_environment")
#   validate_prediction_input(...)  same arguments as model_execute(); stops with
#                                   the expected format if the data are wrong
# model_execute() runs the same checks and stops early with a clear message
# instead of failing deep inside a model. Section 2 shows common mistakes.
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- 1. Print the expected format for each workflow ---------------------------------
for (task in c("general", "hybrid", "multi_trait", "multi_environment")) {
  cat("\n=====", task, "=====\n")
  utils::str(prediction_input_standard(task), max.level = 1, give.attr = FALSE, list.len = 6)
}

dat <- sim_single_env(n = 40, p = 200)
pheno <- dat$pheno
Geno.data <- dat$geno

ok <- try(validate_prediction_input(task = "general", pheno_data = pheno, geno_data = Geno.data,
                                    response = "Yield", gen_name = "GID"), silent = TRUE)
log_check("Validate correct single-environment input", !inherits(ok, "try-error"), "input format is valid")

# ---- 2. Common mistakes and the messages they produce --------------------------------
# Each call below is wrong on purpose; PredictProR should stop with a clear message.
expect_message_check <- function(label, res, pattern) {
  msg <- if (inherits(res, "try-error")) trimws(gsub("\\s+", " ", conditionMessage(attr(res, "condition")))) else "no error"
  cat(label, "\n   PredictProR said:", substr(msg, 1, 200), "\n")
  log_check(label, grepl(pattern, msg, ignore.case = TRUE), "stopped early with a clear message")
}

# (a) a phenotyped line that is missing from the genotypes
bad_pheno <- pheno
bad_pheno$GID[1] <- "NOT_IN_GENO"
res <- try(PredictProR::model_execute(
  # Input data
  pheno_data = bad_pheno,
  geno_data = Geno.data,
  # Trait and genotype identifiers
  gen_name = "GID",
  response = "Yield",
  # Mixed-model formulas
  fixed = ~1,
  random = ~GID,
  # Candidate model
  GS_model = "BayesB",
  # Genotype QC and preprocessing
  qc_filtering = TRUE,
  impute = TRUE,
  # Bayesian model controls
  nIter = 600L, burnIn = 200L, thin = 5L,
  # Output
  system_database = TRUE, message = FALSE, verbose = FALSE
), silent = TRUE)
expect_message_check("Mistake: phenotyped line missing from the genotypes", res, "Missing IDs|genotyp")

# (b) a Bayesian model without a random term
res <- try(PredictProR::model_execute(
  # Input data
  pheno_data = pheno,
  geno_data = Geno.data,
  # Trait and genotype identifiers
  gen_name = "GID",
  response = "Yield",
  # Mixed-model formulas (random term missing)
  fixed = ~1,
  random = NULL,
  # Candidate model
  GS_model = "BayesB",
  # Bayesian model controls
  nIter = 600L, burnIn = 200L, thin = 5L,
  # Output
  system_database = TRUE, message = FALSE, verbose = FALSE
), silent = TRUE)
expect_message_check("Mistake: Bayesian model without a random term", res, "random")

# (c) a kernel model from markers without gmatrix_method
res <- try(PredictProR::model_execute(
  # Input data
  pheno_data = pheno,
  geno_data = Geno.data,
  # Trait and genotype identifiers
  gen_name = "GID",
  response = "Yield",
  # Mixed-model formulas
  fixed = ~1,
  random = ~GID,
  # Candidate model
  GS_model = "GBLUP_BRR",
  # Relationship construction (missing: gmatrix_method)
  gmatrix_method = NULL,
  # Bayesian model controls
  nIter = 600L, burnIn = 200L, thin = 5L,
  # Output
  system_database = TRUE, message = FALSE, verbose = FALSE
), silent = TRUE)
expect_message_check("Mistake: kernel model from markers without gmatrix_method", res, "gmatrix")

# (d) multi-environment data (repeated lines) without heter_groups
met <- sim_met(n = 30, p = 200)
res <- try(PredictProR::model_execute(
  # Input data
  pheno_data = met$pheno,
  geno_data = met$geno,
  # Trait and genotype identifiers
  gen_name = "GID",
  response = "Yield",
  # Multi-environment mode (missing: heter_groups = "Env")
  heter_groups = NULL,
  # Mixed-model formulas
  fixed = ~1,
  random = ~GID,
  # Candidate model
  GS_model = "RKHS",
  gmatrix_method = "VanRaden",
  # Bayesian model controls
  nIter = 600L, burnIn = 200L, thin = 5L,
  # Output
  system_database = TRUE, message = FALSE, verbose = FALSE
), silent = TRUE)
expect_message_check("Mistake: repeated lines (MET data) without heter_groups", res, "heter_groups|environment")

scenario_report()
