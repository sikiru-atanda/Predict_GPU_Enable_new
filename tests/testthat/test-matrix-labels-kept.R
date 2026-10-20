# ML/DL variance tables carry `genetic_variance_not_identifiable` once per
# environment/trait. Up to 0.20.223 those rows were read as environment
# labels and replaced the prediction covariance matrices' own dimnames, so
# every Row/Column of the published tables read "not_identifiable".

test_that("a labelled matrix keeps its dimnames", {
  m <- matrix(c(1, 0.5, 0.5, 2), 2, dimnames = list(c("E1", "E2"), c("E1", "E2")))
  tbl <- PredictProR:::gp_matrix_to_component_table(
    m, "Prediction_covariance_environments", "prediction_covariance",
    env_labels = c("not_identifiable", "not_identifiable")
  )
  expect_setequal(unique(tbl$Row), c("E1", "E2"))
  expect_false(anyDuplicated(paste(tbl$Row, tbl$Column)) > 0)
})

test_that("an unlabelled matrix still takes the environment labels", {
  m <- matrix(c(1, 0.5, 0.5, 2), 2)
  tbl <- PredictProR:::gp_matrix_to_component_table(m, "G", "environment_covariance",
                                                    env_labels = c("A", "B"))
  expect_setequal(unique(tbl$Row), c("A", "B"))
})

test_that("ML prediction covariance tables keep environment labels end to end", {
  pred <- expand.grid(GID = paste0("g", 1:8), Env = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
  set.seed(1)
  pred$Predicted_value <- rnorm(nrow(pred))
  pred$Observed_value <- pred$Predicted_value + rnorm(nrow(pred), sd = 0.3)
  bundle <- PredictProR:::gp_ml_predictive_matrix_bundle(pred, "GID", "Env", error_pred = pred,
                                                         n_bootstrap = 5L)
  vc <- data.frame(Env = c("E1", "E2", "E3"), Component = "genetic_variance_not_identifiable",
                   Components = NA_real_, Standard_error = NA_real_, stringsAsFactors = FALSE)
  public <- PredictProR:::gp_attach_correlation_matrices(
    list(variance_components = vc), c(bundle, list(variance_components = vc))
  )
  for (tbl in list(public$covariance_components, public$correlation_components)) {
    expect_true(is.data.frame(tbl) && nrow(tbl) > 0)
    expect_false(any(tbl$Row == "not_identifiable"))
    for (mtx in unique(tbl$Matrix)) {
      s <- tbl[tbl$Matrix == mtx, ]
      expect_equal(nrow(s), 9L)
      expect_setequal(unique(s$Row), c("E1", "E2", "E3"))
    }
  }
  cor_tbl <- public$correlation_components
  d <- cor_tbl[cor_tbl$Row == cor_tbl$Column, "Value"]
  expect_true(all(abs(d - 1) < 1e-12))
})
