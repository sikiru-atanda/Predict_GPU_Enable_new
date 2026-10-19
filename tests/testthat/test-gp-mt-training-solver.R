test_that("six-trait R adapter preserves predictions and Option A with the training-grid solver", {
  skip_if(is.null(tryCatch(PredictProR:::gp_detect_gp_python(), error = function(e) NULL)),
          "No configured GP Python runtime")
  withr::local_envvar(c(PREDICTPRO_GP_DEVICE = "cpu", OMP_NUM_THREADS = "1",
                      MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"))
  set.seed(861)
  n <- 20L; ids <- sprintf("g%02d", seq_len(n)); traits <- sprintf("T%02d", 1:6)
  z <- matrix(rnorm(n * 30L), n)
  k <- tcrossprod(z) / 30 + diag(.1, n)
  dimnames(k) <- list(ids, ids)
  y <- matrix(rnorm(n * 6L), n) + z[, 1] %o% seq(.2, 1.2, length.out = 6L)
  ph <- data.frame(GID = ids, y); names(ph)[-1L] <- traits
  ph[c(2L, 7L, 18L, 20L), traits] <- NA_real_
  for (partial in c(FALSE, TRUE)) {
  if (partial) for (j in seq_along(traits)) ph[3L+j,traits[j]] <- NA_real_
  args <- list(pheno_data = ph, gmatrix = k, response = traits, gen_name = "GID",
    gp_auto_policy = FALSE, tune_gp = FALSE, prediction_blend = "none",
    gp_uncertainty_calibration = "none", max_iter = 1L, prediction_output = "all")
  # One optimizer update is an adapter check, not evidence of convergence.
  se <- suppressWarnings(do.call(PredictProR::gp_multi_trait_model, c(args, list(return_se = TRUE))))
  mean <- suppressWarnings(do.call(PredictProR::gp_multi_trait_model, c(args, list(return_se = FALSE))))
  expect_equal(se$info$training_solver, if (partial) "missing_cell_correction" else "balanced_kronecker")
  expect_equal(se$info$n_training_genotypes, 16L)
  expect_equal(nrow(se$predictions), n * 6L)
  expect_true(all(is.finite(se$predictions$Prediction_SE_latent)))
  expect_equal(se$predictions$Prediction, mean$predictions$Prediction, tolerance = 1e-10)
  expect_false("Prediction_SE_latent" %in% names(mean$predictions))
  expect_null(se$result$cov_traits)
  public <- PredictProR:::gp_model_execute_multitrait_predicted_values(
    se$predictions, ph, traits, "GID")
  expect_equal(public$PEV, public$Standard_error^2, tolerance = 1e-12)
  expect_true(all(is.finite(public$Reliability_variance_input)))
  expect_true(all(is.finite(public$Reliability_reference_variance)))
  expect_true(all(grepl("^genetic_prior_variance", public$Reliability_basis)))
  expect_true(all(c("Train", "Test") %in% public$Train_Test_Label))
  }
})
