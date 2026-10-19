# These are ordinary integration tests: a configured GP interpreter is the only
# runtime prerequisite. They must not be hidden behind the smoke-test flag.
test_that("model_execute GP inference satisfies Option A for ST, MT and MET", {
  skip_if_no_gp_python()
  wd <- tempfile("gp-contract-"); dir.create(wd); withr::local_dir(wd)
  withr::local_envvar(c(PREDICTPRO_GP_DEVICE = "cpu", OMP_NUM_THREADS = "1",
                      MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"))
  set.seed(221)
  n <- 32L
  ids <- sprintf("g%02d", seq_len(n))
  Z <- matrix(rnorm(n * 45L), n)
  K <- tcrossprod(Z) / ncol(Z) + diag(.1, n)
  dimnames(K) <- list(ids, ids)
  g <- drop(t(chol(K)) %*% rnorm(n))
  phi <- t(chol(K))
  phi_path <- file.path(wd, "factor.bin")
  con <- file(phi_path, "wb"); writeBin(as.numeric(t(phi)), con, size = 8L); close(con)
  factor_cache <- list(type = "memmap", paths = list(phi_path),
    shapes = list(as.list(dim(phi))), dtypes = list("float64"))
  for (route in c("single", "factor_cache", "multi_trait", "multi_trait_met")) {
    met <- identical(route, "multi_trait_met")
    multi <- route %in% c("multi_trait", "multi_trait_met")
    ph <- if (met) expand.grid(GID = ids, Env = c("A", "B"), stringsAsFactors = FALSE) else data.frame(GID = ids)
    gi <- match(ph$GID, ids)
    ph$T1 <- g[gi] + rnorm(nrow(ph))
    ph$T2 <- .5 * g[gi] + rnorm(nrow(ph))
    ph[gi > n - 6L, c("T1", "T2")] <- NA_real_
    args <- list(pheno_data = ph, gmatrix = K, gen_name = "GID",
      response = if (multi) c("T1", "T2") else "T1",
      response_family = "gaussian", fixed = if (met) ~Env else ~1,
      random = ~GID, GS_model = "GP", gp_return_se = TRUE,
      cross_validation = FALSE, parallel_mode = "sequential",
      system_database = FALSE, message = FALSE, verbose = FALSE)
    if (met) args$heter_groups <- "Env"
    if (multi) args$multi_trait_gp <- TRUE
    if (route == "factor_cache") {
      args$gp_factor_cache <- factor_cache
      args$gp_varcomp_mode <- "mom"
    }
    fit <- do.call(PredictProR::model_execute, args)
    expect_true(PredictProR::validate_prediction_output(fit)$valid, info = route)
    mr <- fit$model_results_by_trait[[1L]]
    if (!"predicted_values" %in% names(mr)) mr <- mr[[1L]]
    pred <- mr$predicted_values
    expect_equal(nrow(pred), nrow(ph) * if (multi) 2L else 1L, info = route)
    expect_true(all(c("Train", "Test") %in% pred$Train_Test_Label), info = route)
    expect_equal(pred$PEV, pred$Standard_error^2, tolerance = 1e-8)
    expect_true(all(is.finite(pred$Standard_error)))
    expect_true(all(grepl("^genetic_(plus_gxe_)?prior_variance", pred$Reliability_basis)))
    expect_true(all(is.finite(pred$Reliability_variance_input)))
    expect_true(all(pred$Reliability_variance_input >= 0))
    expect_true(all(pred$Reliability_variance_input <= pred$Reliability_reference_variance + 1e-6))
    expect_equal(pred$Reliability, pmax(0, pmin(1,
      1 - pred$Reliability_variance_input / pred$Reliability_reference_variance)), tolerance = 1e-8)
    if (route == "factor_cache") {
      params <- mr$model_parameters
      selection <- params$summary[params$stat == "gp_krr_lambda_selection"]
      expect_match(selection, "gcv")
      expect_true(any(grepl("gp_krr_lambda_selection", unlist(fit$run_metadata))))
    }
    if (multi) {
      expect_true(is.data.frame(mr$covariance_parameters))
      expect_true(is.data.frame(mr$correlation_parameters))
      expect_true("std.error" %in% names(mr$covariance_parameters))
      expect_true("inference_status" %in% names(mr$inference_summary))
    }
    if (route == "single") {
      expect_equal(mr$inference_summary$inference_status, "available")
    }
  }
})
