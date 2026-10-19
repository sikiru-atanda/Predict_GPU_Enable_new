test_that("public R MET MoM uses absorbed means and preserves Option A", {
  py <- tryCatch(PredictProR:::gp_detect_gp_python(), error = function(e) NULL)
  skip_if(is.null(py), "No configured GP Python runtime")
  set.seed(2801)
  ids <- paste0("g", 1:8)
  Z <- matrix(rnorm(8*6), 8)
  K <- tcrossprod(Z)/6 + diag(.2, 8)
  dimnames(K) <- list(ids, ids)
  ph <- expand.grid(GID = ids, Env = c("e1", "e2", "e3"), stringsAsFactors = FALSE)
  gi <- match(ph$GID, ids)
  ph$a <- Z[gi, 1] + rnorm(nrow(ph))
  ph$b <- 2*Z[gi, 1] + rnorm(nrow(ph))
  ph[gi > 6, c("a", "b")] <- NA_real_
  Ke <- diag(3); dimnames(Ke) <- list(c("e1", "e2", "e3"), c("e1", "e2", "e3"))
  args <- list(pheno_data = ph, gmatrix = K, response = c("a", "b"),
    gen_name = "GID", heter_groups = "Env", env_similarity = Ke,
    test_set = ids[7:8], varcomp_mode = "mom", prediction_output = "all",
    tune_gp = FALSE, prediction_blend = "none", gp_uncertainty_calibration = "none",
    mom_pcg_tol = 1e-9, mom_pcg_max_iter = 150L)
  means <- do.call(PredictProR::gp_multi_trait_met_model, c(args, list(return_se = FALSE)))
  expect_identical(means$info$mean_design_method, "absorbed_cell_means")
  expect_equal(means$fit$pcg_rhs_count, 1)
  se <- do.call(PredictProR::gp_multi_trait_met_model, c(args, list(return_se = TRUE)))
  expect_identical(se$info$mean_design_method, "absorbed_cell_means")
  expect_true(isTRUE(se$fit$pcg_all_converged))
  expect_equal(se$predictions$Prediction, means$predictions$Prediction, tolerance = 1e-9)
  expect_true(all(is.finite(se$predictions$PEV)))
  expect_equal(se$predictions$PEV, se$predictions$SE_latent^2, tolerance = 1e-9)
})
