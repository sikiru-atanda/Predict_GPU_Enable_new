test_that("requested SEs learn scales unless explicitly disabled", {
  seen <- new.env()
  local_mocked_bindings(
    gp_bridge_direct_execution_enabled = function() TRUE,
    gp_bridge_reticulate_fallback_enabled = function() FALSE,
    gp_with_gp_device_route = function(route, expr) force(expr),
    gp_bridge_fit_mixed_model_direct = function(args, ...) {seen$args <- args; stop("captured")},
    .package = "PredictProR")
  ids <- letters[1:8]; K <- diag(8); dimnames(K) <- list(ids,ids)
  args <- list(pheno_data=data.frame(GID=ids,y=1:8),gmatrix=K,response="y",gen_name="GID",
    test_set=ids[7:8],gp_auto_policy=FALSE,python_bin="unused",project_root=getwd(),
    tune_gp=FALSE,prediction_blend="none",mean_adjustment="none",gp_return_se=TRUE)
  expect_error(do.call(PredictProR::gp_single_trait_model,args),"captured")
  expect_true(seen$args$learn_scales)
  expect_identical(unlist(seen$args$include_components), "g")
  args$gp_learn_scales <- FALSE
  expect_error(do.call(PredictProR::gp_single_trait_model,args),"captured")
  expect_false(seen$args$learn_scales)
})

test_that("large REML CV uses the operator MoM route and announces it", {
  ns <- asNamespace("PredictProR"); seen <- new.env()
  local_mocked_bindings(
    gp_bridge_direct_execution_enabled = function() TRUE,
    gp_bridge_reticulate_fallback_enabled = function() FALSE,
    gp_backend_cv_batch_enabled = function() TRUE,
    gp_with_gp_device_route = function(route, expr) force(expr),
    gp_bridge_fit_mixed_model_cv_direct = function(args, ...) {seen$args <- args; stop("captured")},
    gp_bridge_fit_mixed_model_direct = function(args, ...) {seen$args <- args; stop("captured")},
    .package = "PredictProR")
  ids <- letters[1:8]; ph <- data.frame(GID=ids,y=1:8); K <- diag(8); dimnames(K) <- list(ids,ids)
  pars <- list(pheno_data=ph,response="y",gen_name="GID",gmatrix_model_ready=K,
    python_bin="unused",project_root=getwd(),gp_varcomp_mode="reml",large_n_threshold=5L)
  expect_message(expect_error(ns$gp_backend_cv_predict_batch("GP",ph$y,list(7:8,5:6),pars),"captured"),
    "large.*MoM")
  expect_false(seen$args$fit_reml)
  expect_false(seen$args$compute_ai_se)
  expect_identical(seen$args$varcomp_mode,"mom")
  expect_identical(seen$args$backend,"operator")
  pars$large_n_threshold <- 50000L
  pars$gp_factor_cache <- list(type="test-cache")
  expect_message(expect_error(ns$gp_backend_cv_predict_batch("GP",ph$y,list(7:8,5:6),pars),"captured"),
    "large.*MoM")
  expect_false(seen$args$fit_reml)
  expect_identical(seen$args$varcomp_mode,"mom")
  args <- list(pheno_data=ph,gmatrix=K,response="y",gen_name="GID",test_set=ids[7:8],
    gp_auto_policy=FALSE,python_bin="unused",project_root=getwd(),tune_gp=FALSE,
    prediction_blend="none",mean_adjustment="none",gp_fit_reml=TRUE,
    large_n_threshold=5L,gp_output_level="predict_only")
  expect_message(expect_error(do.call(PredictProR::gp_single_trait_model,args),"captured"),"large.*MoM")
  expect_false(seen$args$fit_reml)
  expect_identical(seen$args$varcomp_mode,"mom")
  expect_identical(seen$args$backend,"operator")
})
