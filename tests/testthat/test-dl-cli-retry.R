test_that("DL bootstrap CLI retries one transient subprocess failure", {
  output_dir <- tempfile("predictpror-dl-retry-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  partial_path <- file.path(output_dir, "bootstrap.csv")
  calls <- 0L
  fake_run <- function(command, args, stdout, stderr) {
    calls <<- calls + 1L
    if (calls == 1L) {
      writeLines("partial", partial_path)
      return(structure("transient crash", status = 1033L))
    }
    expect_false(file.exists(partial_path))
    character()
  }

  expect_warning(
    PredictProR:::gp_dl_run_cli(
      "bootstrap-fit-predict",
      c("--out-dir", output_dir),
      python_bin = file.path(R.home("bin"), "Rscript.exe"),
      max_attempts = 2L,
      run_fn = fake_run
    ),
    "retrying the same deterministic command once"
  )
  expect_identical(calls, 2L)
})

test_that("single-response DL fit-predict retries one transient subprocess failure", {
  # Hybrid and single-trait DL CV call fit-predict once per fold. Without the
  # bounded retry the bootstrap and multi-trait paths already have, one
  # transient native exit (status 1033) aborted the whole CV run.
  seen <- NULL
  local_mocked_bindings(
    gp_dl_run_cli = function(command, args, python_bin = NULL,
                             max_attempts = 1L, run_fn = system2) {
      seen <<- list(command = command, max_attempts = max_attempts)
      invisible(character())
    },
    gp_dl_bridge_read_cli_result = function(out_dir) list(predictions = c(0.1, 0.2)),
    .package = "PredictProR"
  )
  set.seed(1)
  PredictProR:::gp_dl_bridge_fit_predict_raw(
    model_type = "mlp",
    X_train = matrix(rnorm(20), 10, 2),
    y_train = rnorm(10),
    X_test = matrix(rnorm(4), 2, 2)
  )
  expect_identical(seen$command, "fit-predict")
  expect_identical(as.integer(seen$max_attempts), 2L)
})

test_that("DL CLI reports the final subprocess status after bounded retry", {
  calls <- 0L
  always_fails <- function(command, args, stdout, stderr) {
    calls <<- calls + 1L
    structure("backend stopped", status = 1033L)
  }

  expect_warning(
    expect_error(
      PredictProR:::gp_dl_run_cli(
        "bootstrap-fit-predict",
        character(),
        python_bin = file.path(R.home("bin"), "Rscript.exe"),
        max_attempts = 2L,
        run_fn = always_fails
      ),
      "failed with status 1033"
    ),
    "retrying the same deterministic command once"
  )
  expect_identical(calls, 2L)
})

