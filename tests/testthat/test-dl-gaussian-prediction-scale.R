if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

# Regression guard for the DL gaussian double scale-reversion bug: the bootstrap
# matrix was reverted to the original y-scale, then the prediction mean and the
# lower/upper bounds were reverted a SECOND time, producing values ~y_sd x too
# large (e.g. ~1200 when the response centred on ~100). True-prediction outputs
# must sit on the observed response scale.

test_that("DL gaussian true prediction stays on the observed response scale", {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  testthat::skip_if(is.null(py) || !file.exists(py), "No DL Python runtime found")
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_PYTHON = py)

  set.seed(20260526L)
  n <- 40L; p <- 24L
  # Response from a latent signal; geno_data carries it as homozygous 0/2
  # dosages (genotype QC rejects values above 2 without ploidy and removes
  # markers with more than 10% heterozygotes).
  Z <- matrix(rnorm(n * p), n, p)
  X <- ifelse(Z > 0, 2, 0)
  rownames(X) <- paste0("g", seq_len(n)); colnames(X) <- paste0("m", seq_len(p))
  sig <- Z[, 1] - 0.6 * Z[, 2]
  # Response deliberately shifted/scaled: mean ~100, sd ~11.
  y <- 100 + 10 * scale(sig)[, 1] + rnorm(n, sd = 3)
  ph <- data.frame(GID = rownames(X), Trait = y, stringsAsFactors = FALSE)
  test_rows <- c(6L, 17L, 28L, 39L)
  ph$Trait[test_rows] <- NA_real_
  obs <- y[-test_rows]; ctr <- mean(obs); sdy <- stats::sd(obs)

  out <- suppressWarnings(suppressMessages(PredictProR::model_execute(
    pheno_data = ph, geno_data = X, response = "Trait", gen_name = "GID",
    GS_model = "mlp", response_family = "gaussian",
    deterministic = TRUE, random_seed = 1L, device = "cpu",
    epochs = 30L, batch_size = 8L, mlp_neurons_per_layer = as.integer(c(16L, 8L)),
    dropout = 0.1, n_bootstrap = 10L, compile_model = FALSE, use_amp = FALSE,
    system_database = TRUE, message = FALSE
  )))

  pred <- out[[1L]]$model_results$predicted_values
  tr <- pred$Train_Test_Label == "Test"
  pv <- pred$Predicted_value[tr]

  expect_true(all(is.finite(pv)))
  # Must be on the observed scale, NOT ~y_sd x inflated (double-revert -> ~1200).
  expect_true(
    all(abs(pv - ctr) <= 6 * sdy),
    info = sprintf("test predictions {%s} far from observed centre %.1f (sd %.1f)",
                   paste(round(pv, 1), collapse = ", "), ctr, sdy)
  )
  # Interval bounds must bracket the point prediction.
  expect_true(all(pred$lower_bound[tr] <= pv & pv <= pred$upper_bound[tr]))
})
