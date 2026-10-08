test_that("gaussian DL CV hot path predicts consistently", {
  skip_if(
    !identical(Sys.getenv("PREDICTPROR_RUN_DL_BENCHMARK_TESTS"), "1"),
    "Set PREDICTPROR_RUN_DL_BENCHMARK_TESTS=1 to run DL runtime benchmark tests."
  )
  py_bin <- gp_preferred_python(purpose = "dl")
  skip_if(is.null(py_bin) || !nzchar(py_bin) || !file.exists(py_bin), "No DL Python runtime available.")

  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_DL_PYTHON = py_bin, PREDICTPRO_PYTHON = py_bin)
  dl_prewarm_runtime(py_bin)

  set.seed(123)
  n <- 24L
  p <- 12L
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  storage.mode(x) <- "double"
  colnames(x) <- paste0("M", seq_len(p))
  y <- as.numeric(0.4 * x[, 1] - 0.25 * x[, 2] + stats::rnorm(n, sd = 0.3))
  tst <- seq(1L, n, by = 4L)

  run_model <- function(model_type) {
    if (identical(model_type, "mlp")) {
      deep_learning_model(
        y = y,
        omics_data = x,
        crossval = TRUE,
        response_family = "gaussian",
        tst = tst,
        model_type = "mlp",
        epochs = 2L,
        batch_size = 8L,
        mlp_neurons_per_layer = as.integer(c(16L, 8L)),
        validation_split = 0.1,
        deterministic = TRUE,
        random_seed = 123L,
        n_bootstrap = 2L
      )
    } else {
      deep_learning_model(
        y = y,
        omics_data = x,
        crossval = TRUE,
        response_family = "gaussian",
        tst = tst,
        model_type = "ft_transformer",
        epochs = 2L,
        batch_size = 8L,
        ft_d_model = 16L,
        ft_heads = 2L,
        ft_layers = 1L,
        ft_ff_mult = 2L,
        ft_dropout = 0.1,
        ft_token_dropout = 0.1,
        ft_use_cls = TRUE,
        use_grouping = FALSE,
        validation_split = 0.1,
        deterministic = TRUE,
        random_seed = 123L,
        n_bootstrap = 2L
      )
    }
  }

  for (model_type in c("mlp", "ft_transformer")) {
    first <- system.time(pred1 <- run_model(model_type))
    second <- system.time(pred2 <- run_model(model_type))

    expect_length(pred1, length(tst))
    expect_length(pred2, length(tst))
    expect_true(all(is.finite(pred1)))
    expect_true(all(is.finite(pred2)))
    expect_lte(unname(second[["elapsed"]]), max(5, unname(first[["elapsed"]]) * 3))
  }
})
