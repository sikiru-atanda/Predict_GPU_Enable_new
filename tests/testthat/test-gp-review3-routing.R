test_that("GP reliability uses the fitted basis and preserves missing-scale reasons", {
  ph <- data.frame(GID = c("a", "b"), Env = c("x", "y"), y = c(1, NA))
  pred <- data.frame(Genotype = ph$GID, Prediction = c(1,2), Prediction_SE_latent = c(.3,.4),
    Reliability_variance_input = c(.08,.12), Reliability_reference_variance = c(.2,.2),
    Reliability_basis = "genetic_prior_variance")
  f <- PredictProR:::gp_model_execute_single_trait_predicted_values
  out <- f(pred, ph, "y", "GID", heter_groups = "Env", test_set = "b")
  expect_identical(out$Reliability_basis, rep("genetic_prior_variance",2))
  pred$Reliability_reference_variance <- NA_real_
  pred$Reliability_missing_reason <- "genetic variance not estimated on this route"
  out <- f(pred, ph, "y", "GID", heter_groups = "Env", test_set = "b")
  expect_true(all(is.na(out$Reliability)))
  expect_identical(out$Reliability_remarks, rep(pred$Reliability_missing_reason[1],2))
  expect_equal(out$PEV, out$Standard_error^2)
})

test_that("REML GP CV requests the same fitted estimator as the delivered refit", {
  ns <- asNamespace("PredictProR"); seen <- new.env()
  local_mocked_bindings(
    gp_bridge_direct_execution_enabled = function() TRUE,
    gp_bridge_reticulate_fallback_enabled = function() FALSE,
    gp_backend_cv_batch_enabled = function() TRUE,
    gp_with_gp_device_route = function(route, expr) force(expr),
    gp_bridge_fit_mixed_model_cv_direct = function(args, ...) {seen$cv <- args; stop("captured")},
    gp_bridge_fit_mixed_model_direct = function(args, ...) {seen$refit <- args; stop("captured")},
    .package = "PredictProR")
  ids <- letters[1:8]; ph <- data.frame(GID=ids, y=1:8)
  K <- diag(8); dimnames(K) <- list(ids,ids)
  for (model in c("GP", "GP_FA")) {
    pars <- list(pheno_data=ph,response="y",gen_name="GID",gmatrix_model_ready=K,
      gp_auto_policy=FALSE,python_bin="unused",project_root=getwd(),gp_varcomp_mode="reml")
    expect_error(ns$gp_backend_cv_predict_batch(model, ph$y, list(7:8,5:6), pars),"captured")
    expect_error(PredictProR::gp_single_trait_model(pheno_data=ph,response="y",gen_name="GID",
      gmatrix=K,model_name=model,test_set=ids[7:8],gp_auto_policy=FALSE,python_bin="unused",
      project_root=getwd(),tune_gp=FALSE,prediction_blend="none",mean_adjustment="none",
      gp_full_vc=TRUE),"captured")
    expect_true(seen$cv$fit_reml)
    expect_identical(unlist(seen$refit$include_components), "g")
    expect_identical(seen$cv$learn_scales, seen$refit$learn_scales)
    expect_false(seen$cv$return_se)
    ym <- ph$y; ym[7:8] <- NA_real_
    expect_error(ns$gp_bridge_cv_predict(model,ym,7:8,pars),"captured")
    expect_true(seen$refit$fit_reml)
    expect_true(seen$refit$learn_scales)
    expect_false(seen$refit$return_se)
  }
})
