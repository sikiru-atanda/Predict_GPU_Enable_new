test_that("unknown processed plot prefixes fall back to the original name", {
  expect_identical(
    PredictProR:::format_processed_plot_export_name("cv_results_processed_unmapped_plot"),
    "cv_results_processed_unmapped_plot"
  )
})

test_that("output directory stems retain exact public model names", {
  expect_identical(
    PredictProR:::gp_output_directory_stem("single_trait", "Kernel-GBLUP"),
    "single_trait__model-Kernel-GBLUP"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem("multi_trait_gp", "GP"),
    "multi_trait_gp__model-Gaussian-Process-GBLUP"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem("multi_trait_gp", "GP_FA"),
    "multi_trait_gp__model-FA-GBLUP"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem("met_ml_dl", "RandomForest"),
    "met_ml_dl__model-RandomForest"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem("hybrid_dl", "DenseNeuralNet"),
    "hybrid_dl__model-DenseNeuralNet"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem(
      "CV_results",
      c("GBLUP", "GBLUP_BRR", "RKHS", "DenseNeuralNet")
    ),
    "CV_results__models-GBLUP+GBLUP_BRR+RKHS+DenseNeuralNet"
  )
  expect_identical(
    PredictProR:::gp_output_directory_stem(
      "CV_results_Gaussian-Process-GBLUP",
      "GP"
    ),
    "CV_results_Gaussian-Process-GBLUP"
  )
})

test_that("exported results report their model-identifying directory", {
  tmp_dir <- withr::local_tempdir()
  withr::local_dir(tmp_dir)

  out <- PredictProR:::results_handling(
    GS_model = "GP",
    res_model_output = list(
      predicted_values = data.frame(
        GID = c("g1", "g2"),
        Predicted_value = c(1.1, 1.3),
        stringsAsFactors = FALSE
      )
    ),
    system_database = FALSE,
    plot_filename = "multi_trait_gp"
  )

  expect_identical(out$export_status, "Successful")
  expect_true(dir.exists(out$export_directory))
  expect_match(
    basename(out$export_directory),
    "^multi_trait_gp__model-Gaussian-Process-GBLUP_[0-9]{2}-[0-9]{2}-[0-9]{4}_"
  )
  expect_true(file.exists(file.path(out$export_directory, "Predicted_Value.csv")))
})
