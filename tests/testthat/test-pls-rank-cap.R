test_that("PLS caps components at the rank of each fit's training data", {
  # sklearn PLSRegression with n_components above the rank of the centred
  # training matrix fits noise and extrapolates wildly. Bootstrap resamples
  # have fewer unique rows (lower rank) than the outer training set, so an
  # ncomp that was valid for the full data produced predictions around 1e11
  # inside the bootstrap, which then failed the prediction-interval contract
  # (hybrid ML PartialLeastSquare true prediction on a small panel).
  py <- PredictProR:::gp_detect_ml_python()
  skip_if(is.null(py) || !nzchar(py) || !file.exists(py), "No ML Python runtime")
  set.seed(3)
  n <- 13L; p <- 60L
  # rank-deficient design like hybrid mid-parent genotypes: 13 rows, rank 6
  basis <- matrix(sample(0:2, 7 * p, replace = TRUE), 7, p)
  x <- basis[c(1:7, sample(7, n - 7, replace = TRUE)), ]
  storage.mode(x) <- "double"
  ids <- sprintf("H%02d", seq_len(n)); rownames(x) <- ids
  # residual noise matters: surplus components fit it along near-zero
  # directions of X, which is what makes the coefficients explode
  y <- as.vector(x %*% stats::rnorm(p, sd = 0.05)) + 8 + stats::rnorm(n, sd = 0.5)
  x_test <- matrix(sample(0:2, 3 * p, replace = TRUE), 3, p,
                   dimnames = list(sprintf("T%d", 1:3), NULL))
  storage.mode(x_test) <- "double"

  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  fit <- PredictProR:::gp_python_tabular_model(
    pheno_object = data.frame(.HybridRowID = ids, .HybridResponse = y, stringsAsFactors = FALSE),
    geno_omic_object = x, geno_omic_test_object = x_test,
    response = ".HybridResponse", gen_name = ".HybridRowID", response_family = "gaussian",
    model_type = "pls", model_label = "PartialLeastSquare",
    model_params = list(ncomp = 10L), scaling = FALSE, centering = FALSE,
    n_bootstrap = 30L, system_database = TRUE
  )
  pred <- as.data.frame(fit$predicted_values)
  expect_true(all(is.finite(pred$Predicted_value)))
  # rank-safe fits stay on the scale of the response
  expect_true(all(abs(pred$Predicted_value - mean(y)) < 20 * stats::sd(y) + 1))
})
