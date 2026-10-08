test_that("model-based reliability keeps the engine's basis text when its value is unchanged", {
  # gp_met_add_gaussian_uncertainty() recomputes model-based reliability as
  # 1 - PEV / reference variance. It kept engine-specific PEV_basis and
  # uncertainty-source text but always replaced Reliability_basis with a
  # generic sentence, so e.g. BGLR hybrid predictions lost "BGLR model-implied
  # hybrid reliability" while every number was still BGLR's own.
  pev <- c(0.2, 0.3, 0.1, 0.4)
  ref <- c(1, 1, 1, 1)
  engine_rel <- 1 - pev / ref
  engine_rel[4] <- ref[4] / (ref[4] + pev[4])   # engine used a different form here
  engine_basis <- "BGLR model-implied hybrid reliability = max(0, min(1, 1 - PEV/genetic_variance))"
  pred <- data.frame(
    GID = paste0("H", 1:4),
    Predicted_value = c(0.1, 0.2, 0.3, 0.4),
    Observed_value = c(0.2, NA, 0.1, NA),
    Train_Test_Label = c("Train", "Test", "Train", "Test"),
    Standard_error = sqrt(pev),
    PEV = pev,
    Reliability = engine_rel,
    Reliability_reference_variance = ref,
    Reliability_basis = engine_basis,
    PEV_basis = "BGLR posterior target-draw variance",
    stringsAsFactors = FALSE
  )
  out <- PredictProR:::gp_met_add_gaussian_uncertainty(pred_obs = pred, gen_name = "GID", heter_groups = NULL)

  expect_equal(out$Reliability, pmax(0, pmin(1, 1 - pev / ref)))
  # rows 1-3: value unchanged -> engine provenance kept
  expect_identical(out$Reliability_basis[1:3], rep(engine_basis, 3L))
  # row 4: value recomputed -> generic basis describes the new value
  expect_false(identical(out$Reliability_basis[4], engine_basis))
  expect_true(grepl("model-specific", out$Reliability_basis[4], fixed = TRUE))
  expect_identical(out$PEV_basis, rep("BGLR posterior target-draw variance", 4L))
})
