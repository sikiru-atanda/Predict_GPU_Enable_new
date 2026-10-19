test_that("direct requested GP SEs match full-VC SEs for the same REML estimator", {
  skip_if(is.null(tryCatch(PredictProR:::gp_detect_gp_python(),error=function(e) NULL)),
    "No configured GP Python runtime")
  set.seed(481)
  ids <- paste0("g",1:36); Z <- matrix(rnorm(36*10),36)
  K <- tcrossprod(Z)/10 + diag(.1,36); dimnames(K) <- list(ids,ids)
  ph <- data.frame(GID=ids,y=8 + 2*Z[,1] + Z[,2] + rnorm(36,sd=.6))
  ph$y[31:36] <- NA_real_
  args <- list(pheno_data=ph,gmatrix=K,response="y",gen_name="GID",model_name="GP",
    test_set=ids[31:36],gp_auto_policy=FALSE,gp_iters=100L,tune_gp=FALSE,gp_fit_reml=TRUE,
    mean_adjustment="none",prediction_blend="none",gp_uncertainty_calibration="none")
  se <- do.call(PredictProR::gp_single_trait_model,c(args,list(gp_return_se=TRUE)))
  vc <- do.call(PredictProR::gp_single_trait_model,c(args,list(gp_full_vc=TRUE)))
  expect_length(se$predictions$Prediction_SE_latent,6L)
  expect_true(all(is.finite(se$predictions$Prediction_SE_latent)))
  expect_equal(se$predictions$Prediction_SE_latent,vc$predictions$Prediction_SE_latent,tolerance=1e-6)
  expect_true(all(is.finite(se$predictions$Reliability_reference_variance)))
  expect_true(all(se$predictions$Reliability_reference_variance > 0))
  expect_gt(sd(se$predictions$Prediction), 1e-5)
  expect_true(any(grepl("!var", vc$result$varcomp$component)))
  args$gp_fit_reml <- FALSE
  fixed <- do.call(PredictProR::gp_single_trait_model,c(args,list(gp_return_se=TRUE,gp_learn_scales=FALSE)))
  expect_true(all(is.na(fixed$predictions$Reliability_reference_variance)))
  expect_identical(unique(fixed$predictions$Reliability_missing_reason),"genetic variance not estimated on this route")
})
