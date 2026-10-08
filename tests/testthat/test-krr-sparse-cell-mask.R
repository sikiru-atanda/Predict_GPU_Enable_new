test_that("KRR inferred sparse MET masks remain cell-level in GP split", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1, NA, 3, 4),
    stringsAsFactors = FALSE
  )

  inferred <- PredictProR:::gp_public_train_test_idx(
    pheno,
    gid_col = "GID",
    y_col = "Yield",
    test_set = "g1",
    test_set_source = "inferred_missing_response"
  )

  explicit <- PredictProR:::gp_public_train_test_idx(
    pheno,
    gid_col = "GID",
    y_col = "Yield",
    test_set = "g1",
    test_set_source = "explicit"
  )

  expect_identical(inferred$test_mask, c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(inferred$test_idx, as.integer(1L))
  expect_identical(inferred$train_idx, as.integer(c(0L, 2L, 3L)))
  expect_identical(explicit$test_mask, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(explicit$test_idx, as.integer(c(0L, 1L)))
})

test_that("KRR true-prediction output labels inferred sparse MET masks by cell", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1, NA, 3, 4),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = pheno$GID,
    Env = pheno$Env,
    Predicted_value = c(1.1, 2.2, 3.1, 4.1),
    Prediction_SE_latent = c(0.1, 0.2, 0.1, 0.1),
    Prediction_Var_latent = c(0.01, 0.04, 0.01, 0.01),
    stringsAsFactors = FALSE
  )

  inferred <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = predictions,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    test_set = "g1",
    test_set_source = "inferred_missing_response"
  )

  explicit <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = predictions,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    test_set = "g1",
    test_set_source = "explicit"
  )

  expect_identical(inferred$Train_Test_Label, c("Train", "Test", "Train", "Train"))
  expect_identical(inferred$Train_Test, c("Train", "Test", "Train", "Train"))
  expect_identical(explicit$Train_Test_Label, c("Test", "Test", "Train", "Train"))
  expect_identical(explicit$Train_Test, c("Test", "Test", "Train", "Train"))
})
