test_that("DL seed plans are deterministic and expose fit counts", {
  plan <- dl_seed_plan(
    n_bootstrap = 4L,
    dl_n_seeds = 3L,
    random_seed = 123L,
    calibration_folds = 2L
  )

  expect_identical(
    plan$seed_manifest$training_seed,
    as.integer(c(123L, 104852L, 209581L))
  )
  expect_identical(
    plan$fit_counts$component,
    c("bootstrap_true_prediction", "heldout_risk_calibration", "total")
  )
  expect_equal(plan$fit_counts$model_fits, c(12, 6, 18))
  expect_identical(plan$aggregation, "mean")
  expect_match(plan$selection_rule, "no_test_outcome_selection", fixed = TRUE)
})

test_that("single-trait DL skips calibration fits when disabled", {
  ids <- paste0("g", 1:8)
  x <- matrix(seq_len(24), nrow = 8,
              dimnames = list(ids, paste0("m", 1:3)))
  pheno <- data.frame(GID = ids, Trait = seq_len(8))
  local_mocked_bindings(
    gp_tuning_folds = function(...) stop("calibration folds must not run"),
    gp_dl_bridge_bootstrap = function(model_type, X_train, y_train, X_pred,
                                      n_bootstrap, response_family, dl_args,
                                      seed = NULL, class_levels = NULL,
                                      training_seeds = NULL,
                                      seed_aggregation = "mean", ...) {
      list(bootstrap = matrix(rep(mean(y_train), n_bootstrap * nrow(X_pred)),
                              nrow = n_bootstrap))
    },
    .package = "PredictProR"
  )
  out <- PredictProR:::deep_learning_model(
    pheno_object = pheno, geno_omic_object = x,
    response = "Trait", gen_name = "GID", response_family = "gaussian",
    model_type = "mlp", para_tunning = FALSE, n_bootstrap = 3L,
    scaling = FALSE, centering = FALSE, system_database = TRUE,
    dl_internal_calibration = FALSE, message = FALSE
  )
  pred <- out$predicted_values
  expect_equal(nrow(pred), 8L)
  expect_true(all(is.finite(pred$Predicted_value)))
  expect_true(all(is.na(pred$Standard_error)))
  expect_true(all(is.na(pred$PEV)))
  expect_true(all(pred$Prediction_uncertainty_source ==
    "unavailable_internal_calibration_disabled_by_user"))
  expect_equal(out$dl_computation_plan$model_fits, c(3, 0, 3))
})

test_that("DL internal calibration is a validated public switch", {
  expect_identical(formals(model_execute)$dl_internal_calibration, TRUE)
  expect_true(PredictProR:::gp_dl_calibration_enabled(TRUE))
  expect_false(PredictProR:::gp_dl_calibration_enabled(FALSE))
  expect_error(PredictProR:::gp_dl_calibration_enabled(NA), "TRUE or FALSE")
  expect_error(PredictProR:::gp_dl_calibration_enabled(1), "TRUE or FALSE")
})

test_that("explicit DL seeds are validated without hidden replacement", {
  plan <- dl_seed_plan(
    n_bootstrap = 2L,
    dl_seeds = c(7L, 17L, 29L),
    random_seed = 99L,
    calibration_folds = 0L
  )
  expect_identical(
    plan$seed_manifest$training_seed,
    as.integer(c(7L, 17L, 29L))
  )
  expect_true(all(plan$seed_manifest$source == "explicit"))
  expect_error(
    dl_seed_plan(dl_n_seeds = 2L, dl_seeds = c(1L, 2L, 3L)),
    "must equal length"
  )
  expect_error(dl_seed_plan(dl_seeds = c(1L, 1L)), "unique seeds")
  expect_error(dl_seed_plan(dl_seeds = c(1L, -2L)), "non-negative integer")
  zero_plan <- dl_seed_plan(
    n_bootstrap = 1L,
    dl_seeds = c(0L, 2L),
    random_seed = 0L,
    calibration_folds = 0L
  )
  expect_identical(zero_plan$seed_manifest$training_seed, c(0L, 2L))
  expect_error(
    dl_seed_plan(dl_n_seeds = 2L, dl_seed_aggregation = "best"),
    "does not select a best seed"
  )
})

test_that("one DL training seed is serialized as a JSON array", {
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  PredictProR:::gp_dl_write_seed_json(123L, path)
  text <- paste(readLines(path, warn = FALSE), collapse = "")

  expect_match(text, "^\\[")
  expect_match(text, "\\]$")
  expect_identical(as.integer(jsonlite::fromJSON(path)), 123L)
})

test_that("seed-batch aggregation fits every requested seed and averages all", {
  observed_seeds <- integer()
  testthat::local_mocked_bindings(
    gp_dl_bridge_fit_predict_batch = function(jobs) {
      vapply(jobs, function(job) {
        seed <- as.integer(job$dl_args$random_seed)
        observed_seeds <<- c(observed_seeds, seed)
        seed + c(0, 2)
      }, numeric(2), USE.NAMES = FALSE) |>
        as.data.frame() |>
        as.list()
    },
    .package = "PredictProR"
  )
  manifest <- PredictProR:::gp_dl_seed_manifest(dl_seeds = c(11L, 31L))
  jobs <- list(list(
    id = "fold_1",
    response_family = "gaussian",
    X_train = matrix(1, 2, 1),
    y_train = c(0, 1),
    X_test = matrix(1, 2, 1),
    dl_args = list()
  ))

  out <- PredictProR:::gp_dl_bridge_fit_predict_seed_batch(jobs, manifest)
  expect_identical(observed_seeds, as.integer(c(11L, 31L)))
  expect_equal(out$predictions[[1L]], c(21, 23))
  expect_length(out$per_seed[[1L]], 2L)
})

test_that("true-prediction seed tables retain scaled per-seed results", {
  payload <- list(
    seed_manifest = data.frame(
      seed_index = 1:2,
      training_seed = c(10L, 20L),
      source = "explicit"
    ),
    seed_predictions = rbind(c(-1, 0), c(1, 2)),
    seed_variability = data.frame(
      seed_variance = c(1, 4),
      seed_range = c(2, 4)
    )
  )
  tables <- PredictProR:::gp_dl_seed_output_tables(
    boot_payload = payload,
    ids = c("g1", "g2"),
    response_family = "gaussian",
    y_scaler = list(mean = 100, std = 10),
    gen_name = "GID"
  )

  expect_equal(tables$predictions$Predicted_value, c(90, 100, 110, 120))
  expect_equal(tables$variability$DL_seed_prediction_SD, c(10, 20))
  expect_equal(tables$variability$DL_seed_prediction_range, c(20, 40))
  expect_identical(tables$predictions$training_seed, c(10L, 10L, 20L, 20L))
})

test_that("Gaussian bootstrap and retained seed predictions share one response scale", {
  scaled_predictions <- rbind(c(-1, 0), c(1, 2))
  scaler <- list(mean = 100, std = 10)
  response_bootstrap <- PredictProR:::gp_dl_bootstrap_response_scale(
    scaled_predictions,
    response_family = "gaussian",
    y_scaler = scaler
  )
  payload <- list(
    seed_manifest = data.frame(
      seed_index = 1:2,
      training_seed = c(10L, 20L),
      source = "explicit"
    ),
    seed_predictions = scaled_predictions,
    seed_variability = data.frame(
      seed_variance = c(2, 2),
      seed_range = c(2, 2)
    )
  )
  tables <- PredictProR:::gp_dl_seed_output_tables(
    boot_payload = payload,
    ids = c("g1", "g2"),
    response_family = "gaussian",
    y_scaler = scaler,
    gen_name = "GID"
  )
  seed_means <- tapply(
    tables$predictions$Predicted_value,
    tables$predictions$GID,
    mean
  )

  expect_equal(response_bootstrap, scaled_predictions * 10 + 100)
  expect_equal(as.numeric(seed_means[c("g1", "g2")]), colMeans(response_bootstrap))
})

test_that("binary DL bootstrap clipping preserves bootstrap-by-target dimensions", {
  raw <- matrix(
    c(-0.2, 0.25, 0.75, 1.2, 0.4, 0.6),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("bootstrap_1", "bootstrap_2"), c("g1", "g2", "g3"))
  )
  clipped <- PredictProR:::gp_dl_bootstrap_response_scale(
    raw,
    response_family = "binary"
  )

  expect_true(is.matrix(clipped))
  expect_identical(dim(clipped), dim(raw))
  expect_identical(dimnames(clipped), dimnames(raw))
  expect_equal(clipped, matrix(
    c(0, 0.25, 0.75, 1, 0.4, 0.6),
    nrow = 2L,
    byrow = TRUE,
    dimnames = dimnames(raw)
  ))
})

test_that("model_execute exposes the DL seed controls", {
  formals_names <- names(formals(model_execute))
  expect_true(all(c(
    "dl_n_seeds", "dl_seeds", "dl_seed_aggregation"
  ) %in% formals_names))
  expect_null(formals(model_execute)$dl_n_seeds)
  expect_null(formals(model_execute)$dl_seeds)
  expect_identical(formals(model_execute)$dl_seed_aggregation, "mean")
})

test_that("public result standardization retains seed tables separately", {
  prediction <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1, 2),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.1, NA_real_),
    Standard_error = c(0.2, 0.3),
    PEV = c(0.04, 0.09),
    lower_bound = c(0.6, 1.4),
    upper_bound = c(1.4, 2.6),
    stringsAsFactors = FALSE
  )
  raw <- list(
    model_parameters = data.frame(stat = "model_type", summary = "mlp"),
    predicted_values = prediction,
    variance_components = PredictProR:::gp_empty_variance_components(),
    dl_seed_manifest = data.frame(seed_index = 1:2, training_seed = c(11L, 31L)),
    dl_computation_plan = data.frame(component = "total", model_fits = 14),
    dl_seed_predictions = data.frame(
      seed_index = c(1L, 2L), GID = "g2", Predicted_value = c(1.8, 2.2)
    ),
    dl_seed_variability = data.frame(
      GID = c("g1", "g2"), DL_seed_prediction_SD = c(0.1, 0.2)
    )
  )
  public <- PredictProR:::gp_standardize_public_model_result(
    raw,
    gen_name = "GID"
  )

  expect_true(all(c(
    "dl_seed_manifest", "dl_computation_plan",
    "dl_seed_predictions", "dl_seed_variability"
  ) %in% names(public)))
  expect_identical(
    names(public$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns()
  )
})

test_that("seed tables are written as inspectable native CSV files", {
  output_dir <- tempfile("predictpror-dl-seed-export-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  seed_outputs <- list(
    dl_seed_manifest = data.frame(seed_index = 1:2, training_seed = c(11L, 31L)),
    dl_computation_plan = data.frame(component = "total", model_fits = 14),
    dl_seed_predictions = data.frame(
      seed_index = c(1L, 2L), GID = "g1", Predicted_value = c(1.8, 2.2)
    ),
    dl_seed_variability = data.frame(
      GID = "g1", DL_seed_prediction_SD = 0.2,
      DL_seed_prediction_range = 0.4
    )
  )
  PredictProR:::write_nested_csv_tables(
    seed_outputs,
    pathout = output_dir,
    prefix = "model_results"
  )

  expected <- file.path(output_dir, paste0(
    "model_results_",
    names(seed_outputs),
    ".csv"
  ))
  expect_true(all(file.exists(expected)))
  expect_identical(
    utils::read.csv(expected[[1L]])$training_seed,
    c(11L, 31L)
  )
})
