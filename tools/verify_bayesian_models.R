.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = TRUE), .libPaths()))

verify_args <- commandArgs(trailingOnly = TRUE)
run_full_variance_audit <- "--full-variance-audit" %in% verify_args

suppressPackageStartupMessages({
  library(pkgload)
  library(BGLR)
})

pkgload::load_all(".")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    stop(msg, call. = FALSE)
  }
}

assert_has_names <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  assert_true(length(missing) == 0, paste(label, "is missing:", paste(missing, collapse = ", ")))
}

assert_prediction_table <- function(df, label, require_env = FALSE, require_uncertainty = TRUE) {
  assert_true(is.data.frame(df), paste(label, "is not a data.frame"))
  assert_true(nrow(df) > 0, paste(label, "is empty"))
  required <- c("Predicted_value", "Standard_error", "PEV", "Reliability")
  if (isTRUE(require_env)) {
    required <- c(required, "Env")
  }
  missing <- setdiff(required, names(df))
  assert_true(length(missing) == 0, paste(label, "missing columns:", paste(missing, collapse = ", ")))
  assert_true(any(is.finite(as.double(df$Predicted_value))), paste(label, "has no finite predictions"))
  if (isTRUE(require_uncertainty)) {
    assert_true(any(is.finite(as.double(df$Standard_error))), paste(label, "has no finite standard errors"))
    assert_true(any(is.finite(as.double(df$PEV))), paste(label, "has no finite PEV values"))
    assert_true(any(is.finite(as.double(df$Reliability))), paste(label, "has no finite reliability values"))
  }
}

assert_aggregated_prediction_table <- function(df, label) {
  assert_true(is.data.frame(df), paste(label, "is not a data.frame"))
  assert_true(nrow(df) > 0, paste(label, "is empty"))
  required <- c("Predicted_value", "Standard_error", "PEV", "Reliability", "Train_Test_Label", "Genetic_variance")
  missing <- setdiff(required, names(df))
  assert_true(length(missing) == 0, paste(label, "missing columns:", paste(missing, collapse = ", ")))
  assert_true(any(is.finite(as.double(df$Predicted_value))), paste(label, "has no finite predictions"))
  assert_true(any(is.finite(as.double(df$Standard_error))), paste(label, "has no finite standard errors"))
  assert_true(any(is.finite(as.double(df$PEV))), paste(label, "has no finite PEV values"))
  assert_true(any(is.finite(as.double(df$Reliability))), paste(label, "has no finite reliability values"))
  assert_true(any(is.finite(as.double(df$Genetic_variance))), paste(label, "has no finite genetic variance values"))
}

inspect_bglr_outputs <- function() {
  set.seed(7)
  X <- matrix(rnorm(60), nrow = 10, dimnames = list(paste0("g", 1:10), paste0("m", 1:6)))
  y <- rnorm(10)
  y[1] <- NA_real_

  model_specs <- list(
    BRR = list(
      eta = list(X = X, model = "BRR", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("b", "varB", "SD.b", "SD.varB"),
      file_patterns = c("_ETA_1_b.bin", "_ETA_1_varB.dat", "_mu.dat", "_varE.dat")
    ),
    BayesA = list(
      eta = list(X = X, model = "BayesA", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("b", "varB", "S", "SD.b", "SD.varB", "SD.S"),
      file_patterns = c("_ETA_1_b.bin", "_ETA_1_ScaleBayesA.dat", "_mu.dat", "_varE.dat")
    ),
    BayesB = list(
      eta = list(X = X, model = "BayesB", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("b", "varB", "probIn", "S", "SD.b", "SD.varB", "SD.probIn", "SD.S"),
      file_patterns = c("_ETA_1_b.bin", "_ETA_1_parBayesB.dat", "_mu.dat", "_varE.dat")
    ),
    BayesC = list(
      eta = list(X = X, model = "BayesC", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("b", "varB", "probIn", "SD.b", "SD.varB", "SD.probIn"),
      file_patterns = c("_ETA_1_b.bin", "_ETA_1_parBayesC.dat", "_mu.dat", "_varE.dat")
    ),
    BL = list(
      eta = list(X = X, model = "BL", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("b", "lambda", "tau2", "SD.b"),
      file_patterns = c("_ETA_1_b.bin", "_ETA_1_lambda.dat", "_mu.dat", "_varE.dat")
    ),
    RKHS = list(
      eta = list(K = tcrossprod(scale(X)), model = "RKHS", saveEffects = TRUE),
      top = c("yHat", "SD.yHat", "mu", "SD.mu", "varE", "SD.varE", "ETA"),
      eta_names = c("u", "varU", "uStar", "SD.u", "SD.varU"),
      file_patterns = c("_ETA_1_varU.dat", "_mu.dat", "_varE.dat")
    )
  )

  for (model_name in names(model_specs)) {
    spec <- model_specs[[model_name]]
    save_at <- paste0("verify_", model_name, "_")
    fit <- BGLR::BGLR(
      y = y,
      ETA = list(spec$eta),
      nIter = 120,
      burnIn = 20,
      thin = 2,
      verbose = FALSE,
      saveAt = save_at
    )

    assert_has_names(fit, spec$top, paste("BGLR top-level output for", model_name))
    assert_has_names(fit$ETA[[1]], spec$eta_names, paste("BGLR ETA output for", model_name))

    output_files <- list.files(pattern = paste0("^", save_at))
    for (pattern in spec$file_patterns) {
      assert_true(any(grepl(pattern, output_files, fixed = TRUE)),
                  paste("BGLR files for", model_name, "missing pattern", pattern))
    }

    unlink(output_files)
    cat("verify_case:", "bglr-structure", model_name, "ok\n")
  }
}

build_single_env_inputs <- function() {
  set.seed(11)
  pheno <- data.frame(
    GID = paste0("g", 1:12),
    Yield = rnorm(12),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(2, 5, 9)] <- NA_real_
  geno <- matrix(sample(0:2, 12 * 20, replace = TRUE), nrow = 12,
                 dimnames = list(pheno$GID, paste0("m", 1:20)))
  list(pheno = pheno, geno = geno)
}

verify_single_env_models <- function() {
  inputs <- build_single_env_inputs()
  models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")

  for (model_name in models) {
    res <- PredictProR:::bayes_finalize_A_B_C_BL_BRR(
      random = ~ GID,
      GS_model = model_name,
      response = "Yield",
      pheno_data = inputs$pheno,
      geno_data = inputs$geno,
      gen_name = "GID",
      nIter = 120,
      burnIn = 20,
      thin = 2
    )

    assert_true(is.list(res), paste(model_name, "did not return a list"))
    assert_has_names(res, c("bayes_result", "bayes_model"), paste(model_name, "finalize output"))
    assert_has_names(
      res$bayes_result,
      c("Predicted_value", "Variance_components", "Estimated_breeding_value", "Total_estimated_breeding_value"),
      paste(model_name, "bayes_result")
    )
    assert_prediction_table(res$bayes_result$Predicted_value, paste(model_name, "Predicted_value"))

    var_comp <- res$bayes_result$Variance_components
    assert_true(all(c("genetic_variance", "residual_variance", "heritability") %in% rownames(var_comp)),
                paste(model_name, "variance components missing required rows"))
    assert_true(any(is.finite(as.double(var_comp[, "Components"]))), paste(model_name, "variance components are not finite"))
    assert_true(any(is.finite(as.double(var_comp[, "Standard_error"]))), paste(model_name, "variance component uncertainty is not finite"))

    ebv_total <- res$bayes_result$Total_estimated_breeding_value
    assert_true(all(c("Estimated_breeding_value", "Prediction_error_variance", "Standard_error", "Reliability") %in% names(ebv_total)),
                paste(model_name, "total EBV table missing uncertainty columns"))
    assert_true(any(is.finite(as.double(ebv_total$Standard_error))), paste(model_name, "total EBV SE is not finite"))
    assert_true(any(is.finite(as.double(ebv_total$Prediction_error_variance))), paste(model_name, "total EBV PEV is not finite"))
    assert_true(any(is.finite(as.double(ebv_total$Reliability))), paste(model_name, "total EBV reliability is not finite"))

    cat("verify_case:", "single-env", model_name, "ok\n")
  }
}

verify_single_env_cv <- function() {
  inputs <- build_single_env_inputs()
  models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  tst <- which(is.na(inputs$pheno$Yield))
  y_na <- inputs$pheno$Yield

  for (model_name in models) {
    prep <- PredictProR:::bayes_finalize_A_B_C_BL_BRR(
      random = ~ GID,
      GS_model = model_name,
      response = "Yield",
      pheno_data = inputs$pheno,
      geno_data = inputs$geno,
      gen_name = "GID",
      nIter = 120,
      burnIn = 20,
      thin = 2,
      cross_validation = TRUE
    )

    preds <- PredictProR:::bayes_mod_cv(
      y = y_na,
      ETA = prep$bayes_ETA$ETA,
      weights = NULL,
      bayes_para = prep$bayes_para,
      tst = tst,
      bayes_model = model_name,
      bayes_trait = "Yield"
    )

    assert_true(length(preds) == length(tst), paste(model_name, "CV predictions length mismatch"))
    assert_true(any(is.finite(as.double(preds))), paste(model_name, "CV predictions are not finite"))
    cat("verify_case:", "single-env-cv", model_name, "ok\n")
  }
}

build_met_inputs <- function() {
  set.seed(12)
  ids <- paste0("g", 1:8)
  pheno <- expand.grid(GID = ids, Env = c("E1", "E2"), stringsAsFactors = FALSE)
  pheno$Yield <- rnorm(nrow(pheno))
  pheno$Yield[c(2, 11)] <- NA_real_
  gmatrix <- diag(length(ids))
  rownames(gmatrix) <- colnames(gmatrix) <- ids
  list(pheno = pheno, gmatrix = gmatrix)
}

verify_met_true_prediction <- function() {
  inputs <- build_met_inputs()

  rkhs_res <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = inputs$pheno,
    gmatrix = inputs$gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2
  )

  assert_prediction_table(rkhs_res$bayes_result$Predicted_value, "RKHS MET Predicted_value", require_env = TRUE)
  assert_aggregated_prediction_table(
    rkhs_res$bayes_result$Total_Predicted_value,
    "RKHS MET Total_Predicted_value"
  )
  rkhs_total_ebv <- rkhs_res$bayes_result$Total_estimated_breeding_value
  assert_true(all(c("Estimated_breeding_value", "Prediction_error_variance", "Standard_error", "Reliability") %in% names(rkhs_total_ebv)),
              "RKHS MET total EBV table missing uncertainty columns")
  assert_true(any(is.finite(as.double(rkhs_total_ebv$Standard_error))), "RKHS MET total EBV SE is not finite")
  cat("verify_case:", "met-truepred", "RKHS", "ok\n")

  prep_gblup_brr <- PredictProR::model_prep_bayes_cv(
    random = ~ GID + GID:Env,
    GS_model_cv = "GBLUP_BRR",
    response = "Yield",
    gen_name = "GID",
    pheno_data = inputs$pheno,
    gmatrix = inputs$gmatrix,
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2
  )
  assert_true("GBLUP_BRR" %in% names(prep_gblup_brr), "GBLUP_BRR prep output missing preserved key")
  cat("verify_case:", "met-prep", "GBLUP_BRR", "ok\n")
}

verify_met_cv <- function() {
  inputs <- build_met_inputs()
  tst <- which(is.na(inputs$pheno$Yield))
  y_na <- inputs$pheno$Yield

  rkhs_prep <- PredictProR:::bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = inputs$pheno,
    gmatrix = inputs$gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  rkhs_preds <- PredictProR:::bayes_mod_cv(
    y = y_na,
    ETA = rkhs_prep$bayes_ETA$ETA,
    weights = NULL,
    bayes_para = rkhs_prep$bayes_para,
    tst = tst,
    bayes_model = "RKHS",
    bayes_trait = "Yield"
  )
  assert_true(length(rkhs_preds) == length(tst), "RKHS MET CV prediction length mismatch")
  assert_true(any(is.finite(as.double(rkhs_preds))), "RKHS MET CV predictions are not finite")
  cat("verify_case:", "met-cv", "RKHS", "ok\n")

  gblup_brr_prep <- PredictProR::model_prep_bayes_cv(
    random = ~ GID + GID:Env,
    GS_model_cv = "GBLUP_BRR",
    response = "Yield",
    gen_name = "GID",
    pheno_data = inputs$pheno,
    gmatrix = inputs$gmatrix,
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2
  )

  gblup_brr_preds <- PredictProR:::bayes_mod_cv(
    y = y_na,
    ETA = gblup_brr_prep$GBLUP_BRR$bayes_ETA$ETA,
    weights = NULL,
    bayes_para = gblup_brr_prep$GBLUP_BRR$bayes_para,
    tst = tst,
    bayes_model = "BRR",
    bayes_trait = "Yield"
  )
  assert_true(length(gblup_brr_preds) == length(tst), "GBLUP_BRR MET CV prediction length mismatch")
  assert_true(any(is.finite(as.double(gblup_brr_preds))), "GBLUP_BRR MET CV predictions are not finite")
  cat("verify_case:", "met-cv", "GBLUP_BRR", "ok\n")
}

run_case_capture <- function(case_label, expr) {
  tryCatch({
    force(expr)
    list(case = case_label, ok = TRUE, message = "ok")
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("verify_case:", case_label, "failed:", msg, "\n")
    list(case = case_label, ok = FALSE, message = msg)
  })
}

case_results <- list(
  run_case_capture("bglr-structure", inspect_bglr_outputs()),
  run_case_capture("single-env-models", verify_single_env_models()),
  run_case_capture("single-env-cv", verify_single_env_cv()),
  run_case_capture("met-true-pred", verify_met_true_prediction()),
  run_case_capture("met-cv", verify_met_cv())
)

if (isTRUE(run_full_variance_audit)) {
  case_results[[length(case_results) + 1L]] <- run_case_capture(
    "full-variance-audit",
    {
      audit_script <- normalizePath(
        file.path("tools", "audit_bglr_variance_extraction.R"),
        winslash = "/",
        mustWork = TRUE
      )
      rscript_bin <- file.path(
        R.home("bin"),
        if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript"
      )
      status <- system2(rscript_bin, args = shQuote(audit_script))
      assert_true(
        identical(status, 0L),
        paste("Full BGLR variance-extraction audit failed with status", status)
      )
      cat("verify_case:", "full-variance-audit", "ok\n")
    }
  )
}

failed_cases <- vapply(case_results, function(x) !isTRUE(x$ok), logical(1))
if (any(failed_cases)) {
  cat("verify_bayesian_models: failures\n")
  for (res in case_results[failed_cases]) {
    cat("* ", res$case, ": ", res$message, "\n", sep = "")
  }
  quit(status = 1)
}

cat("verify_bayesian_models: ok\n")
