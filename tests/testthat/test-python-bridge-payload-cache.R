test_that("Python bridge payload cache reuses prepared payloads", {
  prefix <- paste0("CACHE_TEST_", Sys.getpid())
  calls <- 0L
  build_payload <- function(input_dir) {
    calls <<- calls + 1L
    writeLines("cached", file.path(input_dir, "payload.txt"), useBytes = TRUE)
    invisible(input_dir)
  }

  payload <- list(
    x = matrix(seq_len(12L), nrow = 3L),
    y = c(1, 2, 3)
  )
  first <- gp_bridge_payload_cache_prepare(
    prefix = prefix,
    purpose = "unit",
    key_parts = payload,
    cells = 12L,
    build_fun = build_payload,
    min_cells = 1L,
    max_entries = 2L
  )
  second <- gp_bridge_payload_cache_prepare(
    prefix = prefix,
    purpose = "unit",
    key_parts = payload,
    cells = 12L,
    build_fun = build_payload,
    min_cells = 1L,
    max_entries = 2L
  )

  expect_false(first$cache_hit)
  expect_true(second$cache_hit)
  expect_identical(normalizePath(first$input_dir, winslash = "/", mustWork = TRUE), second$input_dir)
  expect_equal(calls, 1L)
  expect_true(file.exists(file.path(second$input_dir, "payload.txt")))

  state <- gp_bridge_payload_cache_state(prefix)
  root <- get(".root", envir = state, inherits = FALSE)
  unlink(root, recursive = TRUE, force = TRUE)
})
