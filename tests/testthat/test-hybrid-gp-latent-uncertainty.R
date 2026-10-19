# Independent universal-kriging equations, conditional on the supplied covariance.
hybrid_latent_oracle <- function(K, X, train, lambda) {
  V <- K[train, train, drop = FALSE] + diag(lambda, length(train))
  C <- K[, train, drop = FALSE]
  ViX <- solve(V, X[train, , drop = FALSE])
  delta <- X - C %*% ViX
  unname(diag(K - C %*% solve(V, t(C)) +
         delta %*% solve(crossprod(X[train, , drop = FALSE], ViX), t(delta))))
}

hybrid_latent_fixture <- function() {
  ph <- data.frame(HybridID = paste0("h", 1:6), Female = rep(c("f1", "f2"), 3),
                   Male = rep(c("m1", "m2", "m3"), each = 2),
                   HybridCross = paste0("h", 1:6), Env = rep(c("e1", "e2"), each = 3),
                   a = c(1, 3, NA, 4, 2, NA), b = c(9, NA, 14, 20, 12, 17))
  K <- outer(1:6, 1:6, function(i, j) 0.6^abs(i - j))
  ge <- 0.2 * K * outer(ph$Env, ph$Env, "==")
  list(ph = ph, kernels = list(female = .5*K, male = .3*K, sca = .2*K,
                               base = K, gxe = ge, total = K + ge))
}

expect_hybrid_latent_public <- function(raw, expected) {
  # Exercise the actual public result formatter, including its uncertainty fallbacks.
  public <- PredictProR:::gp_standardize_hybrid_model_result(
    list(predicted_values = raw, gp_model = list(lambda = .4)),
    hybrid_id = "HybridID", female_parent = "Female", male_parent = "Male", model = "GP"
  )$predicted_values
  expect_equal(public$Standard_error^2, expected, tolerance = 1e-7)
  expect_equal(public$PEV, public$Standard_error^2, tolerance = 1e-12)
  expect_equal(public$lower_bound, public$Predicted_value - qnorm(.975)*sqrt(expected),
               tolerance = 1e-7)
  expect_equal(public$upper_bound, public$Predicted_value + qnorm(.975)*sqrt(expected),
               tolerance = 1e-7)
  expect_true(all(grepl("latent", public$PEV_basis)))
  expect_true(all(is.na(public$Reliability)))
  expect_true(all(is.na(public$Reliability_variance_input)))
  expect_true(all(is.na(public$Reliability_reference_variance)))
  expect_true(all(grepl("^genetic_prior_variance", public$Reliability_basis)))
  expect_true(all(public$Reliability_remarks == "genetic variance not estimated on this route"))
}

test_that("single-trait hybrid SE is latent GLS uncertainty on training and test rows", {
  d <- hybrid_latent_fixture(); ph <- d$ph
  X <- model.matrix(~ Env, ph); train <- which(is.finite(ph$a)); lambda <- .4
  fit <- PredictProR:::gp_hybrid_gp_fit_core(ph$a, d$kernels, X, train, lambda)
  expected <- hybrid_latent_oracle(d$kernels$total, X, train, lambda) * var(ph$a[train])
  raw <- PredictProR:::gp_hybrid_gp_prediction_frame(
    ph, "a", "HybridID", "Female", "Male", "Env", fit, "GP")
  expect_equal(unname(fit$latent_var), expected, tolerance = 1e-7)
  expect_equal(unname(fit$observed_var - fit$latent_var), rep(lambda * var(ph$a[train]), 6))
  expect_equal(raw$Standard_error^2, expected, tolerance = 1e-7)
  expect_equal(raw$Prediction_error_variance, expected, tolerance = 1e-7)
  expect_equal(raw$Prediction_Var_latent, expected, tolerance = 1e-7)
  expect_setequal(raw$Train_Test_Label, c("Train", "Test"))
  expect_hybrid_latent_public(raw, expected)
})

test_that("joint hybrid SE uses latent variance and each trait's response scale", {
  d <- hybrid_latent_fixture(); ph <- d$ph
  row <- rep(1:6, 2); trait <- rep(1:2, each = 6)
  means <- vapply(ph[c("a", "b")], mean, numeric(1), na.rm = TRUE)
  sds <- vapply(ph[c("a", "b")], sd, numeric(1), na.rm = TRUE)
  y <- (c(ph$a, ph$b) - means[trait]) / sds[trait]
  train <- which(is.finite(y)); lambda <- .4
  cor <- matrix(c(1, .35, .35, 1), 2)
  X <- model.matrix(~ 0 + factor(trait):factor(ph$Env[row]))
  K <- cor[trait, trait] * d$kernels$total[row, row]
  fit <- PredictProR:::gp_hybrid_gp_fit_joint_core(y, d$kernels, cor,
    X_joint = X, train_obs_idx = train, all_row = row, all_trait = trait, lambda = lambda)
  expected <- hybrid_latent_oracle(K, X, train, lambda) * unname(sds[trait])^2
  raw <- PredictProR:::gp_hybrid_gp_joint_prediction_frame(
    ph, c("a", "b"), "HybridID", "Female", "Male", "Env", fit,
    means, sds, row, trait, "GP")
  expect_equal(unname(fit$latent_var_std) * unname(sds[trait])^2, expected, tolerance = 1e-7)
  expect_equal(unname(fit$observed_var_std - fit$latent_var_std), rep(lambda, 12))
  expect_equal(raw$Standard_error^2, expected, tolerance = 1e-7)
  expect_equal(raw$Prediction_error_variance, expected, tolerance = 1e-7)
  expect_equal(raw$Prediction_Var_latent, expected, tolerance = 1e-7)
  expect_setequal(raw$Train_Test_Label, c("Train", "Test"))
  expect_hybrid_latent_public(raw, expected)
})

test_that("model_execute preserves latent hybrid uncertainty in all four layouts", {
  d <- hybrid_latent_fixture()
  female_K <- diag(2); dimnames(female_K) <- list(c("f1", "f2"), c("f1", "f2"))
  male_K <- diag(3); dimnames(male_K) <- list(paste0("m", 1:3), paste0("m", 1:3))
  single <- PredictProR:::gp_hybrid_gp_prediction_frame
  joint <- PredictProR:::gp_hybrid_gp_joint_prediction_frame
  captured <- NULL
  testthat::local_mocked_bindings(
    gp_hybrid_gp_prediction_frame = function(...) { captured <<- single(...); captured },
    gp_hybrid_gp_joint_prediction_frame = function(...) { captured <<- joint(...); captured },
    .package = "PredictProR"
  )
  for (multi in c(FALSE, TRUE)) for (met in c(FALSE, TRUE)) {
    # Matching missing cells avoid the unrelated trait-specific holdout warning.
    ph <- d$ph; ph$b[is.na(ph$a)] <- NA_real_
    ph$b[2] <- 11
    run <- function() PredictProR::model_execute(
      pheno_data = ph, female_gmatrix = female_K, male_gmatrix = male_K,
      response = if (multi) c("a", "b") else "a", gen_name = "HybridID",
      heter_groups = if (met) "Env" else NULL,
      GS_model = "GP", response_family = "gaussian", hybrid_gp = TRUE,
      female_parent = "Female", male_parent = "Male", gp_backend = "r",
      hybrid_gp_lambda = .4, system_database = TRUE, message = FALSE
    )
    if (met) {
      expect_warning(result <- run(), "Low connectivity detected across Env")
    } else {
      result <- run()
    }
    contract <- PredictProR::validate_prediction_output(result)
    expect_true(contract$valid, info = paste("multi =", multi, "met =", met))
    public <- result$model_results$predicted_values
    key <- function(x, id) paste(x[[id]], if (multi) x$Trait else "a", sep = "/")
    idx <- match(key(public, "Hybrid_ID"), key(captured, "HybridID"))
    expect_false(anyNA(idx))
    expect_equal(public$PEV, captured$Prediction_Var_latent[idx], tolerance = 1e-12)
    expect_equal(public$Standard_error^2, public$PEV, tolerance = 1e-12)
    expect_true(all(is.na(public$Reliability)))
    expect_true(all(is.na(public$Reliability_variance_input)))
    expect_true(all(is.na(public$Reliability_reference_variance)))
    expect_true(all(public$Reliability_remarks == "genetic variance not estimated on this route"))
    expect_setequal(public$Train_Test_Label, c("Train", "Test"))
  }
})
