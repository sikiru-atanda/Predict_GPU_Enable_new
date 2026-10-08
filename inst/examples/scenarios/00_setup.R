# =============================================================================
# PredictProR scenario examples - shared setup
# =============================================================================
# Every scenario script starts with:
#
#   source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))
#
# It provides
#   1. small simulated data sets in the PredictProR input format;
#   2. checks for optional software (ASReml-R, Python backends, Java), so a
#      call that needs missing software is skipped with a message;
#   3. log_result() / log_skip(), which record each call for
#      run_all_scenarios.R (you can ignore them when you adapt a script).
#
# Optional software - set these BEFORE running, if you use the Python models:
#   Sys.setenv(PREDICTPRO_PYTHON    = "path/to/python")  # classical ML (scikit-learn, xgboost, ...)
#   Sys.setenv(PREDICTPRO_GP_PYTHON = "path/to/python")  # Gaussian process (torch, gpytorch)
#   Sys.setenv(PREDICTPRO_DL_PYTHON = "path/to/python")  # deep learning (torch)
# setup_predictgp_env() and setup_predictdl_env() can create these environments.
# ASReml models need the asreml package and a valid licence; Beagle needs Java.
# =============================================================================

suppressPackageStartupMessages(library(PredictProR))

# Exported results (system_database = FALSE) go under this folder, one
# sub-folder per run. Change it, or set PREDICTPRO_OUTPUT_DIR before starting R.
if (is.null(getOption("PredictProR.output_dir")) && !nzchar(Sys.getenv("PREDICTPRO_OUTPUT_DIR"))) {
  options(PredictProR.output_dir = file.path(tempdir(), "PredictProR_scenarios"))
}

# ---- Simulated data in PredictProR format ------------------------------------

# Genotypes: lines in rows (row names = line IDs), markers in columns, coded
# 0/1/2. Inbred lines (0/2), so the default heterozygosity QC keeps markers.
sim_geno <- function(n = 80, p = 400, ids = sprintf("L%03d", seq_len(n)), seed = 1) {
  set.seed(seed)
  freq <- stats::runif(p, 0.1, 0.9)
  G <- matrix(2L * stats::rbinom(n * p, 1, rep(freq, each = n)), n, p,
              dimnames = list(ids, sprintf("m%04d", seq_len(p))))
  G[, apply(G, 2, stats::var) > 0, drop = FALSE]
}

sim_genetic_value <- function(G, n_qtl = 30, seed = 2) {
  set.seed(seed)
  qtl <- sample(ncol(G), n_qtl)
  as.numeric(scale(G[, qtl, drop = FALSE] %*% stats::rnorm(n_qtl)))
}

# VanRaden genomic relationship matrix (square, row/col names = line IDs).
sim_gmatrix <- function(G) {
  p <- colMeans(G) / 2
  Z <- sweep(G, 2, 2 * p)
  K <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  K + diag(1e-4, nrow(K))
}

# Autopolyploid genotypes (e.g. tetraploid potato, hexaploid sweet potato):
# ALT-allele dosage 0..ploidy, outbred, so most markers are heterozygous.
# Lines are full-sib families: each offspring gets ploidy/2 of the chromosome
# copies of each of its two parents (polysomic inheritance), so lines are related.
# The ploidy is stored as attr(G, "ploidy"); model_execute() also takes ploidy =.
sim_polyploid <- function(n = 80, p = 400, ploidy = 4L, n_parents = 10, h2 = 0.7, test_fraction = 0.2,
                          envs = NULL, seed = 61) {
  set.seed(seed)
  ids <- sprintf("P%03d", seq_len(n))
  freq <- stats::runif(p, 0.1, 0.9)
  # parents: ploidy chromosome copies x markers, alleles 0/1
  parents <- lapply(seq_len(n_parents), function(i)
    matrix(stats::rbinom(ploidy * p, 1, rep(freq, each = ploidy)), ploidy, p))
  half <- ploidy %/% 2L
  gamete <- function(P) vapply(seq_len(p), function(j) sum(P[sample.int(ploidy, half), j]), numeric(1))
  G <- t(vapply(seq_len(n), function(i) {
    pp <- sample.int(n_parents, 2L)
    gamete(parents[[pp[1]]]) + gamete(parents[[pp[2]]])
  }, numeric(p)))
  dimnames(G) <- list(ids, sprintf("m%04d", seq_len(p)))
  G <- G[, apply(G, 2, stats::var) > 0, drop = FALSE]
  g <- sim_genetic_value(G, seed = seed + 1)
  set.seed(seed + 2)
  ph <- if (is.null(envs)) data.frame(GID = ids, stringsAsFactors = FALSE) else
    expand.grid(GID = ids, Env = envs, stringsAsFactors = FALSE)
  env_effect <- if (is.null(envs)) 0 else seq(0, by = 2, length.out = length(envs))[match(ph$Env, envs)]
  ph$Yield <- 10 + env_effect + g[match(ph$GID, ids)] + stats::rnorm(nrow(ph), sd = sqrt((1 - h2) / h2))
  test <- sample(ids, round(test_fraction * n))
  truth <- ph
  ph$Yield[ph$GID %in% test] <- NA
  attr(G, "ploidy") <- as.integer(ploidy)
  list(pheno = ph, geno = G, ploidy = as.integer(ploidy), test = test, truth = truth)
}

# One environment, one trait: one row per line; lines to predict have NA.
sim_single_env <- function(n = 80, p = 400, h2 = 0.7, test_fraction = 0.2, seed = 1) {
  G <- sim_geno(n, p, seed = seed)
  g <- sim_genetic_value(G, seed = seed + 1)
  set.seed(seed + 2)
  y <- 10 + g + stats::rnorm(n, sd = sqrt((1 - h2) / h2))
  pheno <- data.frame(GID = rownames(G), Yield = y, stringsAsFactors = FALSE)
  test <- sample(pheno$GID, round(test_fraction * n))
  truth <- pheno
  pheno$Yield[pheno$GID %in% test] <- NA
  list(pheno = pheno, geno = G, gmatrix = sim_gmatrix(G), test = test, truth = truth)
}

# Multi-environment trial: LONG format, one row per line x environment.
# Genetic values correlate across environments (gen_cor < 1: G x E).
# unbalanced = TRUE drops some line x environment records.
sim_met <- function(n = 60, p = 300, envs = c("E1", "E2", "E3"), unbalanced = TRUE,
                    gen_cor = 0.6, test_fraction = 0.2, seed = 11) {
  G <- sim_geno(n, p, seed = seed)
  g_common <- sim_genetic_value(G, seed = seed + 1)
  g_env <- sapply(seq_along(envs), function(k) {
    sqrt(gen_cor) * g_common + sqrt(1 - gen_cor) * sim_genetic_value(G, seed = seed + 100 + k)
  })
  colnames(g_env) <- envs
  set.seed(seed + 2)
  ph <- expand.grid(GID = rownames(G), Env = envs, stringsAsFactors = FALSE)
  env_effect <- stats::setNames(seq(0, by = 2, length.out = length(envs)), envs)
  g_rec <- g_env[cbind(match(ph$GID, rownames(G)), match(ph$Env, envs))]
  ph$Yield <- 10 + env_effect[ph$Env] + g_rec + stats::rnorm(nrow(ph), sd = 0.7)
  if (unbalanced) ph <- ph[sort(sample(nrow(ph), round(0.85 * nrow(ph)))), ]
  ph$Weight <- 1
  test <- sample(unique(ph$GID), round(test_fraction * n))
  truth <- ph
  ph$Yield[ph$GID %in% test] <- NA
  rownames(ph) <- NULL
  list(pheno = ph, geno = G, gmatrix = sim_gmatrix(G), test = test, truth = truth)
}

# Two correlated traits (Yield, Protein); single environment, or MET with envs.
sim_multi_trait <- function(n = 80, p = 300, envs = NULL, test_fraction = 0.2, seed = 21) {
  G <- sim_geno(n, p, seed = seed)
  g1 <- sim_genetic_value(G, seed = seed + 1)
  g2 <- 0.6 * g1 + 0.8 * sim_genetic_value(G, seed = seed + 2)
  set.seed(seed + 3)
  ph <- if (is.null(envs)) data.frame(GID = rownames(G), stringsAsFactors = FALSE) else
    expand.grid(GID = rownames(G), Env = envs, stringsAsFactors = FALSE)
  k <- match(ph$GID, rownames(G))
  ph$Yield <- 10 + g1[k] + stats::rnorm(nrow(ph), sd = 0.8)
  ph$Protein <- 5 + g2[k] + stats::rnorm(nrow(ph), sd = 0.8)
  test <- sample(rownames(G), round(test_fraction * n))
  ph$Yield[ph$GID %in% test] <- NA
  ph$Protein[ph$GID %in% test] <- NA
  list(pheno = ph, geno = G, gmatrix = sim_gmatrix(G), test = test)
}

# Hybrids (female x male). Parent genotypes for hybrid_asreml / hybrid_bayes /
# hybrid_gp; mid-parent hybrid genotypes for hybrid_ml / hybrid_dl.
sim_hybrid <- function(n_female = 8, n_male = 6, p = 300, test_fraction = 0.2, seed = 31) {
  parents <- c(sprintf("F%02d", seq_len(n_female)), sprintf("M%02d", seq_len(n_male)))
  P <- sim_geno(length(parents), p, ids = parents, seed = seed)
  set.seed(seed + 1)
  gca <- stats::setNames(as.numeric(scale(P %*% stats::rnorm(ncol(P), sd = 0.05))), parents)
  ph <- expand.grid(Female = parents[seq_len(n_female)], Male = parents[n_female + seq_len(n_male)],
                    stringsAsFactors = FALSE)
  ph$HybridID <- paste(ph$Female, ph$Male, sep = "x")
  ph$Yield <- 50 + 2 * gca[ph$Female] + 2 * gca[ph$Male] + stats::rnorm(nrow(ph), sd = 1)
  ph <- ph[, c("HybridID", "Female", "Male", "Yield")]
  test <- sample(ph$HybridID, round(test_fraction * nrow(ph)))
  ph$Yield[ph$HybridID %in% test] <- NA
  hyb_geno <- (P[ph$Female, , drop = FALSE] + P[ph$Male, , drop = FALSE]) / 2
  rownames(hyb_geno) <- ph$HybridID
  list(pheno = ph, parent_geno = P, hybrid_geno = hyb_geno, test = test)
}

# Classification traits: binary (factor), ordinal (ordered), multiclass (factor).
sim_classification <- function(n = 120, p = 300, seed = 41) {
  G <- sim_geno(n, p, seed = seed)
  g <- sim_genetic_value(G, seed = seed + 1)
  set.seed(seed + 2)
  liab <- g + stats::rnorm(n, sd = 0.8)
  ph <- data.frame(
    GID = rownames(G),
    Rust = factor(ifelse(liab > 0, "yes", "no"), levels = c("no", "yes")),
    Score = ordered(cut(liab, c(-Inf, -0.6, 0.6, Inf), labels = c("low", "mid", "high")),
                    levels = c("low", "mid", "high")),
    Type = factor(c("A", "B", "C")[1 + (rank(liab + stats::rnorm(n, sd = 0.5)) %% 3)]),
    stringsAsFactors = FALSE
  )
  test <- sample(ph$GID, round(0.2 * n))
  ph[ph$GID %in% test, c("Rust", "Score", "Type")] <- NA
  list(pheno = ph, geno = G, gmatrix = sim_gmatrix(G), test = test)
}

# Extra omics layer (e.g. transcripts) for the same lines.
sim_omics <- function(ids, n_features = 120, seed = 51) {
  set.seed(seed)
  matrix(stats::rnorm(length(ids) * n_features), length(ids), n_features,
         dimnames = list(ids, sprintf("t%03d", seq_len(n_features))))
}

# ---- Optional software checks -------------------------------------------------
has_asreml <- function() requireNamespace("asreml", quietly = TRUE)

python_for <- function(purpose = c("ml", "gp", "dl")) {
  purpose <- match.arg(purpose)
  env <- switch(purpose,
                ml = c("PREDICTPRO_ML_PYTHON", "PREDICTPRO_PYTHON"),
                gp = "PREDICTPRO_GP_PYTHON",
                dl = c("PREDICTPRO_DL_PYTHON", "PREDICTPRO_PYTHON_DL"))
  for (v in env) if (nzchar(Sys.getenv(v))) return(Sys.getenv(v))
  found <- tryCatch(switch(purpose,
    ml = PredictProR:::gp_preferred_python(purpose = "ml"),
    gp = PredictProR:::gp_preferred_python(purpose = "gp"),
    dl = PredictProR:::gp_detect_dl_python()), error = function(e) "")
  if (is.character(found) && length(found) == 1L && nzchar(found) && file.exists(found)) found else ""
}
has_python <- function(purpose) nzchar(python_for(purpose))
has_java <- function() nzchar(Sys.getenv("PREDICTPRO_JAVA")) || nzchar(Sys.which("java"))

# Models that need each backend (used to skip a model when it is missing).
python_ml_models <- c("Ridge_Regression", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors",
                      "RandomForest", "Xgboost", "CatBoost", "LightGBM")
dl_models <- c("DenseNeuralNet", "TabTransformer", "TabAttention", "TabNet", "LightTreeNet", "FactorNet",
               "CrossNet", "NeuralAdditive", "MixtureOfExperts", "GPNet", "DenseAttentionNet", "Conv1DNet", "ResNet")
needs_for_model <- function(model) {
  if (model %in% python_ml_models) return("python_ml")
  if (model %in% c("Kernel-GBLUP", "Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")) return("python_gp")
  if (model %in% dl_models) return("python_dl")
  if (model %in% "GBLUP") return("asreml")
  character()
}
missing_software <- function(needs) {
  for (n in needs) {
    ok <- switch(n, asreml = has_asreml(), python_ml = has_python("ml"), python_gp = has_python("gp"),
                 python_dl = has_python("dl"), java = has_java(), TRUE)
    if (!isTRUE(ok)) return(n)
  }
  NULL
}

# ---- Reading results ------------------------------------------------------------
# The first data frame holding predictions (any workflow).
find_predictions <- function(res, depth = 0) {
  if (is.data.frame(res) && any(c("Predicted_value", "Predicted_class") %in% names(res))) return(res)
  if (is.list(res) && depth < 7) for (el in res) {
    out <- find_predictions(el, depth + 1)
    if (!is.null(out)) return(out)
  }
  NULL
}
first_data_frame <- function(x, depth = 0) {
  if (is.data.frame(x) && nrow(x) > 0) return(x)
  if (is.list(x) && depth < 4) for (el in x) {
    out <- first_data_frame(el, depth + 1)
    if (!is.null(out)) return(out)
  }
  NULL
}
# Cross-validation summary (metrics per trait x model).
cv_summary <- function(res) {
  cvp <- res$cv_results_processed
  if (is.null(cvp)) return(NULL)
  first_data_frame(cvp$aggregated_data_list) %||% first_data_frame(cvp$hybrid_metric_summary) %||%
    first_data_frame(cvp$multitrait_metric_summary) %||% first_data_frame(cvp$classification_metrics) %||%
    first_data_frame(cvp$overall_model_ranking)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Recording each call for run_all_scenarios.R ---------------------------------
if (!exists("scenario_log", envir = globalenv())) {
  assign("scenario_log", data.frame(scenario = character(), status = character(), seconds = numeric(),
                                    note = character(), stringsAsFactors = FALSE), envir = globalenv())
}
log_row <- function(label, status, seconds, note) {
  tab <- get("scenario_log", envir = globalenv())
  tab[nrow(tab) + 1L, ] <- list(label, status, round(seconds, 1), note)
  assign("scenario_log", tab, envir = globalenv())
  cat(sprintf("=> %-55s %s (%.1f s) %s\n", label, status, seconds, note))
}
# log_result(label, res, model_runtime): PASS when res holds predictions or a
# CV summary, FAIL when model_execute() stopped (res is a "try-error").
log_result <- function(label, res, model_runtime = NULL) {
  secs <- if (is.null(model_runtime)) NA_real_ else unname(model_runtime[["elapsed"]])
  if (inherits(res, "try-error")) {
    msg <- trimws(gsub("\\s+", " ", gsub("=+", " ", conditionMessage(attr(res, "condition")))))
    if (grepl("licen", msg, ignore.case = TRUE)) return(invisible(log_row(label, "SKIP", secs, "ASReml licence unavailable")))
    return(invisible(log_row(label, "FAIL", secs, substr(msg, 1, 160))))
  }
  keep <- if (is.list(res)) res[setdiff(names(res), c("cv_results_processed", "cv_results_raw",
                                                      "cv_results_predicted_vs_observed"))] else res
  pred <- find_predictions(keep)
  if (!is.null(pred)) {
    # a run that returns constant or non-finite predictions did not work well
    if ("Predicted_value" %in% names(pred) && is.numeric(pred$Predicted_value)) {
      v <- pred$Predicted_value
      if (!any(is.finite(v))) return(invisible(log_row(label, "FAIL", secs, "no finite predictions")))
      if (length(unique(round(v[is.finite(v)], 10))) < 2L)
        return(invisible(log_row(label, "FAIL", secs, "all predictions are the same value")))
      # nearly constant: spread under 0.1% of the observed trait's spread
      obs <- if ("Observed_value" %in% names(pred)) suppressWarnings(as.numeric(pred$Observed_value)) else NULL
      sd_obs <- if (length(obs)) stats::sd(obs, na.rm = TRUE) else NA_real_
      if (is.finite(sd_obs) && sd_obs > 0 && stats::sd(v, na.rm = TRUE) < 1e-3 * sd_obs)
        return(invisible(log_row(label, "FAIL", secs, sprintf("predictions nearly constant (sd %.2g vs observed sd %.2g)",
                                                               stats::sd(v, na.rm = TRUE), sd_obs))))
    }
    # GP models: a genetic or residual variance may be NA only when it is
    # flagged not estimable (residual estimated at zero)
    if (grepl("Kernel-GBLUP|Gaussian-Process-GBLUP|FA-GBLUP|Scalable-GBLUP", label)) {
      bad <- unexplained_na_variance(keep)
      if (length(bad)) return(invisible(log_row(label, "FAIL", secs,
                                                 paste("NA variance without a reason:", paste(bad, collapse = ", ")))))
    }
    return(invisible(log_row(label, "PASS", secs, sprintf("%d prediction rows", nrow(pred)))))
  }
  cvs <- cv_summary(res)
  if (!is.null(cvs)) return(invisible(log_row(label, "PASS", secs, sprintf("CV summary with %d rows", nrow(cvs)))))
  invisible(log_row(label, "FAIL", secs, "no prediction table or CV summary found"))
}
# genetic / residual variance rows that are NA without Boundary = "not_estimable"
unexplained_na_variance <- function(x, depth = 0) {
  if (depth > 6 || is.null(x)) return(character())
  if (is.data.frame(x)) return(character())
  out <- character()
  if (is.list(x)) {
    vc <- x[["variance_components"]]
    if (is.data.frame(vc) && all(c("Component", "Components") %in% names(vc))) {
      rows <- grepl("^(genetic|residual)_variance", vc$Component) & is.na(vc$Components)
      flag <- if ("Boundary" %in% names(vc)) vc$Boundary else rep("", nrow(vc))
      out <- c(out, vc$Component[rows & !(flag %in% "not_estimable")])
    }
    for (n in setdiff(names(x), "variance_components")) out <- c(out, unexplained_na_variance(x[[n]], depth + 1))
  }
  unique(out)
}
log_skip <- function(label, needs) invisible(log_row(label, "SKIP", 0, paste("needs", needs)))
log_check <- function(label, ok, note) invisible(log_row(label, if (isTRUE(ok)) "PASS" else "FAIL", 0, note))

scenario_report <- function() {
  tab <- get("scenario_log", envir = globalenv())
  old <- options(width = 250)
  on.exit(options(old), add = TRUE)
  cat("\n================ Scenario summary ================\n")
  print(tab, row.names = FALSE, right = FALSE)
  cat(sprintf("\nPASS %d   SKIP %d   FAIL %d\n", sum(tab$status == "PASS"), sum(tab$status == "SKIP"),
              sum(tab$status == "FAIL")))
  invisible(tab)
}
