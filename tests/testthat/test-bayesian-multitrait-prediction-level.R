# BGLR::Multitrait() returns ETAHat as the posterior mean of the full linear
# predictor (intercept included). Predictions must sit on the phenotype scale,
# not at intercept + ETAHat (twice the intercept), as they did up to 0.20.223.

bayes_level_fixture <- function(seed = 5L, n = 40L, envs = NULL) {
  set.seed(seed)
  ids <- sprintf("g%02d", seq_len(n))
  X <- matrix(rnorm(n * 30), n)
  K <- tcrossprod(scale(X)) / 30 + diag(0.05, n)
  dimnames(K) <- list(ids, ids)
  g <- drop(t(chol(K)) %*% rnorm(n)) * 0.6
  if (is.null(envs)) {
    pheno <- data.frame(GID = ids, T1 = 50 + g + rnorm(n, sd = 0.5),
                        T2 = 20 + 0.5 * g + rnorm(n, sd = 0.5), stringsAsFactors = FALSE)
    pheno$T1[1:5] <- NA_real_
    pheno$T2[3:7] <- NA_real_
  } else {
    pheno <- expand.grid(GID = ids, Env = envs, stringsAsFactors = FALSE)
    offset <- c(0, 5, 10)[match(pheno$Env, envs)]
    pheno$Yield <- 30 + offset + g[match(pheno$GID, ids)] + rnorm(nrow(pheno), sd = 0.5)
    pheno$Yield[pheno$GID %in% ids[1:5]] <- NA_real_
  }
  list(pheno = pheno, K = K)
}

train_bias <- function(pred, group) {
  tr <- pred$Train_Test_Label == "Train" & is.finite(pred$Observed_value)
  tapply(pred$Observed_value[tr] - pred$Predicted_value[tr], pred[[group]][tr], mean)
}

test_that("joint multi-trait Bayesian predictions are on the phenotype scale", {
  skip_if_not_installed("BGLR")
  d <- bayes_level_fixture()
  out <- PredictProR:::bayes_multitrait_joint_fit(
    pheno_data = d$pheno, response = c("T1", "T2"), gen_name = "GID",
    kernels = list(gmatrix = d$K), GS_model = "GBLUP_BRR",
    bayes_para = list(nIter = 600L, burnIn = 200L, thin = 5L),
    verbose = FALSE, random_state = 3L
  )
  pred <- out$bayes_result$Predicted_value
  bias <- train_bias(pred, "Trait")
  expect_true(all(abs(bias) < 0.5), info = paste(names(bias), round(bias, 3), collapse = " "))
  fit <- out$bglr_multitrait_model %||% NULL
  if (!is.null(fit)) {
    expect_equal(colMeans(as.matrix(fit$ETAHat)), as.double(fit$mu), tolerance = 0.2,
                 ignore_attr = TRUE)
  }
  # the intercept's posterior uncertainty is already inside SD.ETAHat
  expect_true(all(pred$Standard_error[pred$Train_Test_Label == "Train"] < 1))
})

test_that("environment-specific-variance Bayesian MET predictions keep environment means", {
  skip_if_not_installed("BGLR")
  d <- bayes_level_fixture(seed = 9L, envs = c("E1", "E2", "E3"))
  res <- PredictProR:::bayes_multitrait_env_heter_fit(
    pheno_data = d$pheno, response = "Yield", gen_name = "GID", heter_groups = "Env",
    kernels = list(gmatrix = d$K), GS_model = "RKHS",
    bayes_para = list(nIter = 600L, burnIn = 200L, thin = 5L)
  )
  pred <- res$bayes_result$Predicted_value
  bias <- train_bias(pred, "Env")
  expect_true(all(abs(bias) < 0.5), info = paste(names(bias), round(bias, 3), collapse = " "))
  ebv <- res$bayes_result$EBV %||% res$bayes_result$Estimated_breeding_value
  if (is.data.frame(ebv) && "Estimated_breeding_value" %in% names(ebv)) {
    # breeding values are deviations: no intercept or environment mean
    expect_lt(abs(mean(ebv$Estimated_breeding_value)), 1)
  }
})
