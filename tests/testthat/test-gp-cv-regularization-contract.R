test_that("batched CV and refitting pass the same lambda-selection rule", {
  ns <- asNamespace("PredictProR")
  seen <- new.env()
  local_mocked_bindings(
    gp_bridge_direct_execution_enabled = function() TRUE,
    gp_bridge_reticulate_fallback_enabled = function() FALSE,
    gp_backend_cv_batch_enabled = function() TRUE,
    gp_with_gp_device_route = function(route, expr) force(expr),
    gp_bridge_fit_mixed_model_cv_direct = function(args, ...) {
      seen$cv <- args; stop("captured CV", call. = FALSE)
    },
    gp_bridge_fit_mixed_model_direct = function(args, ...) {
      seen$refit <- args; stop("captured refit", call. = FALSE)
    }, .package = "PredictProR")
  ids <- letters[1:8]
  ph <- data.frame(GID = ids, Yield = seq_along(ids))
  K <- diag(8); dimnames(K) <- list(ids, ids)
  for (mode in c("default", "fixed")) {
    extra <- if (mode == "fixed") list(lam_select = "fixed", krr_lam = 2,
      krr_lams = c(.1, 2, 20)) else list()
    params <- c(list(pheno_data = ph, response = "Yield", gen_name = "GID",
      gmatrix_model_ready = K, gp_auto_policy = FALSE,
      python_bin = "unused", project_root = getwd()), extra)
    expect_error(ns$gp_backend_cv_predict_batch("GP", ph$Yield,
      folds = list(7:8, 5:6), additional_params = params), "captured CV")
    args <- c(list(pheno_data = ph, response = "Yield", gen_name = "GID", gmatrix = K,
      model_name = "GP", test_set = ids[7:8],
      gp_auto_policy = FALSE, python_bin = "unused", project_root = getwd(),
      tune_gp = FALSE, prediction_blend = "none", mean_adjustment = "none"), extra)
    expect_error(do.call(PredictProR::gp_single_trait_model, args), "captured refit")
    expect_identical(unlist(seen$refit$include_components), "g")
    for (key in c("lam_select", "krr_lam", "krr_lams")) {
      expect_equal(seen$cv[[key]], seen$refit[[key]], info = paste(mode, key))
    }
    expect_equal(seen$cv$lam_select, if (mode == "fixed") "fixed" else "gcv")
  }
})
