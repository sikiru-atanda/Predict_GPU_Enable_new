test_that("prediction input standard exposes every workflow contract", {
  standards <- prediction_input_standard("all")
  expect_identical(
    names(standards),
    c("general", "hybrid", "multi_trait", "multi_environment")
  )
  expect_identical(prediction_input_standard("hybrid"), hybrid_data_standard())
})

test_that("unified input validator resolves general and multi-trait tasks", {
  general <- data.frame(GID = c("g1", "g2"), Yield = c(1, 2))
  general_result <- validate_prediction_input(
    pheno_data = general,
    response = "Yield",
    gen_name = "GID"
  )
  expect_identical(general_result$task, "general")
  expect_true(general_result$valid)

  multi_trait <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, 2, 3),
    Protein = c(10, 11, 12)
  )
  multi_result <- validate_prediction_input(
    pheno_data = multi_trait,
    response = c("Yield", "Protein"),
    gen_name = "GID"
  )
  expect_identical(multi_result$task, "multi_trait")
  expect_true(multi_result$valid)
})

test_that("unified input validator rejects task-incompatible arguments", {
  expect_error(
    validate_prediction_input(
      task = "general",
      pheno_data = data.frame(GID = "g1", Yield = 1),
      response = "Yield",
      gen_name = "GID",
      female_parent = "Female"
    ),
    "not valid for the general input contract"
  )
})
