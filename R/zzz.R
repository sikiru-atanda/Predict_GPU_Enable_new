if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "#CHROM",
    ".data",
    "ALT",
    "BLUP",
    "Class",
    "Confidence",
    "Correct",
    "Env",
    "Freq",
    "Metric",
    "Observed",
    "Observed_class",
    "Observed_value",
    "PC1",
    "PC2",
    "POS",
    "Predicted",
    "Predicted_class",
    "Predicted_value",
    "Prediction_error_variance",
    "Prediction_confidence",
    "Probability",
    "REF",
    "Remarks",
    "Reliability",
    "Reliability_remarks",
    "SelectionFrequency",
    "Set",
    "Standard_error",
    "Summed_BLUP",
    "Train_Test_Label",
    "Type",
    "bin",
    "bottom_20_obs",
    "bottom_20_pred",
    "bottom_quantile",
    "category",
    "composite_score",
    "confidence_remarks",
    "cum_composite_importance",
    "cv_role",
    "cv_scenario",
    "difference",
    "feature",
    "fit_type",
    "fold",
    "gamma_value",
    "geno_data",
    "is_classification",
    "mae",
    "mean_prob",
    "mean_rank_instability_risk",
    "metric",
    "mod",
    "mod_cv",
    "model",
    "msg",
    "obs_rank",
    "observed",
    "observed_rate",
    "pheno",
    "pred_rank",
    "predicted",
    "predicted.value",
    "prob",
    "reliability_remarks",
    "rmse",
    "row_id",
    "status",
    "std.error",
    "top_20_obs",
    "top_20_pred",
    "top_quantile",
    "total",
    "lower_bound",
    "trait",
    "upper_bound",
    "variant_id",
    "y",
    "yhat"
  ))
}

.onLoad <- function(libname, pkgname) {
  options(timeout = max(getOption("timeout", 60), 300))
  if (!nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
    try(gp_configure_torch_runtime_env(), silent = TRUE)
  }
}

.onAttach <- function(libname, pkgname) {
  backend_message <- paste(
    "Python backends use configured subprocess interpreters.",
    "`setup_predictdl_env()` and `setup_predictgp_env()` are optional one-time provisioning helpers.",
    "On servers, set `PREDICTPRO_DL_PYTHON` and/or `PREDICTPRO_GP_PYTHON` before running models."
  )

  if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
    packageStartupMessage(
      "PredictProR loaded.\n",
      backend_message
    )
    return(invisible())
  }

  # Probing every candidate interpreter starts several Python processes (and
  # imports torch); in non-interactive sessions - parallel workers, Rscript
  # batch jobs - report only what is configured and resolve lazily later.
  if (!interactive()) {
    configured <- Sys.getenv(c("PREDICTPRO_PYTHON", "PREDICTPRO_GP_PYTHON", "PREDICTPRO_DL_PYTHON"), unset = "")
    configured <- configured[nzchar(configured)]
    packageStartupMessage(
      "PredictProR loaded.\n",
      backend_message,
      if (length(configured)) paste0("\n", names(configured), ": ", configured, collapse = "") else ""
    )
    return(invisible())
  }

  normalize_python <- function(x) {
    if (inherits(x, "try-error") || is.null(x) || !nzchar(x)) {
      return(NULL)
    }
    normalizePath(x, winslash = "/", mustWork = FALSE)
  }
  configured_python <- function(env_names, fallback) {
    for (env_name in env_names) {
      val <- Sys.getenv(env_name, unset = "")
      if (nzchar(val)) {
        return(normalize_python(val))
      }
    }
    normalize_python(try(fallback(), silent = TRUE))
  }

  py_general <- configured_python("PREDICTPRO_PYTHON", function() gp_preferred_python(purpose = "ml"))
  py_gp <- configured_python("PREDICTPRO_GP_PYTHON", function() gp_preferred_python(purpose = "gp"))
  py_dl <- configured_python(c("PREDICTPRO_DL_PYTHON", "PREDICTPRO_PYTHON_DL"), function() gp_detect_dl_python())
  py_dl_probe <- if (!inherits(py_dl, "try-error") && !is.null(py_dl) && nzchar(py_dl) && file.exists(py_dl)) {
    try(gp_python_probe(py_dl), silent = TRUE)
  } else {
    NULL
  }
  gp_msg <- if (!is.null(py_gp) && nzchar(py_gp)) {
    paste0("\nConfigured GP Python: ", py_gp)
  } else ""
  general_msg <- if (!is.null(py_general) && nzchar(py_general)) {
    paste0("\nConfigured ML Python: ", py_general)
  } else ""
  dl_msg <- if (!is.null(py_dl) && nzchar(py_dl)) {
    cuda_suffix <- if (!inherits(py_dl_probe, "try-error") && is.list(py_dl_probe) && isTRUE(py_dl_probe$cuda_available)) {
      paste0(" [CUDA ", as.integer(py_dl_probe$cuda_device_count %||% 0L), " GPU]")
    } else {
      " [CPU]"
    }
    paste0("\nConfigured DL Python: ", py_dl, cuda_suffix)
  } else ""
  packageStartupMessage(
    "PredictProR loaded.\n",
    backend_message,
    general_msg,
    gp_msg,
    dl_msg
  )
}
