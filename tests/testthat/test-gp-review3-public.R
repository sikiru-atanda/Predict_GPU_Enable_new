if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 for GP runtime checks.")
}

test_that("direct GP SE does not label a fixed scale as estimated genetic variance", {
  set.seed(1); n <- 150L; p <- 500L
  M <- matrix(rbinom(n*p, 2,.3), n,p,
    dimnames=list(sprintf("L%03d",1:n),sprintf("m%03d",1:p)))
  Z <- scale(M); G <- tcrossprod(Z)/p; diag(G) <- diag(G)+1e-4
  g <- drop(Z[,1:30] %*% rnorm(30)); g <- g/sd(g)
  ph <- data.frame(GID=rownames(M),Yield=g+rnorm(n)); ph$Yield[121:150] <- NA
  fit <- PredictProR::gp_single_trait_model(pheno_data=ph,gmatrix=G,response="Yield",
    gen_name="GID",model_name="GP",test_set=ph$GID[121:150],gp_return_se=TRUE,gp_learn_scales=FALSE)
  pv <- fit$predictions
  expect_true(all(is.na(pv$Reliability_reference_variance)))
  expect_identical(unique(pv$Reliability_missing_reason), "genetic variance not estimated on this route")
  expect_equal(pv$PEV, pv$Prediction_SE_latent^2)
})

review3_mt_data <- function() {
  set.seed(20260501L); n <- 10L; p <- 14L
  ids <- paste0("g",1:n); X <- matrix(rnorm(n*p),n,p)
  K <- tcrossprod(scale(X,center=TRUE,scale=FALSE)); K <- K/mean(diag(K))
  diag(K) <- diag(K)+1e-6; dimnames(K) <- list(ids,ids)
  ph <- expand.grid(GID=ids,Env=c("E1","E2"),stringsAsFactors=FALSE)
  j <- match(ph$GID,ids); e <- c(E1=-.2,E2=.35)[ph$Env]
  signal <- X[j,1]-.4*X[j,2]
  ph$a <- signal+e+rnorm(nrow(ph),sd=.1)
  ph$b <- .7*signal+.3*X[j,3]-.5*e+rnorm(nrow(ph),sd=.1)
  S <- matrix(c(1,.25,.25,1),2,dimnames=list(c("E1","E2"),c("E1","E2")))
  list(pheno_data=ph,gmatrix=K,response=c("a","b"),gen_name="GID",heter_groups="Env",
    test_set=tail(ids,2),env_similarity=S,return_se=FALSE,return_trait_correlations=TRUE)
}

test_that("tiny MT-MET smoke fit must converge and publish its covariance", {
  fit <- do.call(PredictProR::gp_multi_trait_met_model, review3_mt_data())
  expect_true(isTRUE(fit$fit$reml_converged))
  expect_true(fit$result$variance_component_status %in% c("promoted","promoted_boundary"))
  expect_true(all(is.finite(fit$result$genetic_covariance)))
  expect_true(all(is.finite(fit$result$genetic_correlation)))
})

test_that("a deliberately unfinished MT-MET fit withholds inference separately", {
  args <- review3_mt_data(); args$reml_max_iter <- 1L
  fit <- do.call(PredictProR::gp_multi_trait_met_model,args)
  expect_false(isTRUE(fit$fit$reml_converged))
  expect_identical(fit$result$variance_component_status,"not_promoted_nonconverged")
  expect_true(all(is.na(fit$result$genetic_covariance)))
  expect_identical(fit$result$inference_summary$inference_status,"unavailable_nonconverged")
})
