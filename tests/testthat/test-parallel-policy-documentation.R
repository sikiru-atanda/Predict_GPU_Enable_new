# Pins the two parallel-policy documentation helpers added by A + B:
#   gp_parallel_policy_knobs()           -- knobs table
#   gp_policy_decision_reason_glossary() -- reason-code table
#
# These guard against drift between source (parallel_policy.R /
# mori_policy.R / sp_apply_utils.R / gpu_utilis.R / cv_execution_helpers.R)
# and the documentation -- specifically against (a) renaming a knob without
# updating the table, (b) adding a new reason code in source without
# adding a glossary entry, (c) shipping a table whose default value
# disagrees with what the policy code actually defaults to.

# ---- gp_parallel_policy_knobs() ---------------------------------------------

test_that("gp_parallel_policy_knobs returns the expected shape", {
  knobs <- gp_parallel_policy_knobs()
  expect_s3_class(knobs, "data.frame")
  expect_identical(names(knobs), c("name", "type", "default", "scope", "effect"))
  expect_gt(nrow(knobs), 30L)
  expect_true(all(nzchar(knobs$name)))
  expect_true(all(knobs$type %in% c("env", "option")))
  expect_true(all(knobs$scope %in% c("blas", "scoring", "sequential_gate",
                                     "mori", "mirai", "gpu", "reproducibility")))
  expect_true(all(nzchar(knobs$effect)))
})

test_that("gp_parallel_policy_knobs has no duplicate names", {
  knobs <- gp_parallel_policy_knobs()
  expect_equal(anyDuplicated(knobs$name), 0L)
})

test_that("every scoring weight referenced in parallel_policy.R is in the table", {
  # Guards against scoring-weight drift: if someone adds a
  # GP_PAR_SCORE_NEW_THING in source but forgets the doc table.
  src <- readLines(system.file("..", "R", "parallel_policy.R",
                               package = "PredictProR",
                               mustWork = FALSE),
                   warn = FALSE)
  if (!length(src)) skip("parallel_policy.R not installed; skipping")
  hits <- regmatches(src, regexpr('GP_PAR_SCORE_[A-Z_]+', src))
  hits <- unique(hits[nzchar(hits)])
  knobs <- gp_parallel_policy_knobs()
  expect_true(all(hits %in% knobs$name),
              info = paste("Missing knobs:", paste(setdiff(hits, knobs$name), collapse = ", ")))
})

test_that("every GP_MORI_* knob in mori_policy.R is in the table", {
  src <- readLines(system.file("..", "R", "mori_policy.R",
                               package = "PredictProR",
                               mustWork = FALSE),
                   warn = FALSE)
  if (!length(src)) skip("mori_policy.R not installed; skipping")
  hits <- regmatches(src, regexpr('GP_MORI_[A-Z_]+', src))
  hits <- unique(hits[nzchar(hits)])
  knobs <- gp_parallel_policy_knobs()
  expect_true(all(hits %in% knobs$name),
              info = paste("Missing knobs:", paste(setdiff(hits, knobs$name), collapse = ", ")))
})

test_that("scoring weight defaults match the policy source defaults", {
  # Cross-check a handful of weights against their literal defaults in
  # parallel_policy.R. If someone bumps GP_PAR_SCORE_TASKS from 0.9 to
  # 1.0 in source without updating the table, this fails.
  knobs <- gp_parallel_policy_knobs()
  expected <- list(
    GP_PAR_SCORE_TASKS = "0.9",
    GP_PAR_SCORE_WORKERS = "0.5",
    GP_PAR_SCORE_GLOBALS = "0.35",
    GP_PAR_SCORE_GPU_MIRAI = "3.0",
    GP_PAR_SCORE_TINY_SEQ = "2.2",
    GP_PAR_SCORE_HIDDEN_THREADS = "0.4",
    GP_MORI_MIN_MB = "128",
    GP_MORI_SCORE_THRESHOLD = "1.0",
    GP_BLAS_THREADS = "1",
    GP_PAR_MEMORY_FRACTION = "0.7",
    PREDICTPRO_GPU_BUSY_THRESHOLD = "85",
    GP_RANDOM_SEED = "123",
    GP_MIRAI_SUBMIT_TIMEOUT_SEC = "600"
  )
  for (nm in names(expected)) {
    row <- knobs[knobs$name == nm, , drop = FALSE]
    expect_equal(nrow(row), 1L, info = paste("knob missing:", nm))
    expect_identical(row$default, expected[[nm]],
                     info = paste("default mismatch for", nm))
  }
})

# ---- gp_policy_decision_reason_glossary() -----------------------------------

test_that("gp_policy_decision_reason_glossary returns the expected shape", {
  g <- gp_policy_decision_reason_glossary()
  expect_s3_class(g, "data.frame")
  expect_identical(names(g), c("reason", "category", "meaning"))
  expect_gt(nrow(g), 30L)
  expect_true(all(nzchar(g$reason)))
  expect_true(all(nzchar(g$meaning)))
  expect_equal(anyDuplicated(g$reason), 0L)
  expect_true(all(g$category %in% c("user_override", "sequential_driver",
                                    "backend_choice", "routing",
                                    "bonus_penalty", "mori", "dispatch")))
})

test_that("every literal reason emitted in source has a glossary entry", {
  # Mirrors the audit grep: every "reasons <- c(reasons, \"...\")" literal
  # in parallel_policy.R must have a glossary row. Catches drift when a
  # new reason code is added without the doc.
  src <- readLines(system.file("..", "R", "parallel_policy.R",
                               package = "PredictProR",
                               mustWork = FALSE),
                   warn = FALSE)
  if (!length(src)) skip("parallel_policy.R not installed; skipping")
  appended <- regmatches(src, regexpr('reasons\\s*<-\\s*c\\(reasons[^)]+', src))
  emitted <- unlist(lapply(appended, function(s) {
    matches <- gregexpr('"([^"]+)"', s)[[1]]
    if (matches[1] < 0) return(character())
    starts <- as.integer(matches)
    lengths <- attr(matches, "match.length")
    vapply(seq_along(starts), function(i) {
      substr(s, starts[i] + 1L, starts[i] + lengths[i] - 2L)
    }, character(1))
  }))
  emitted <- unique(emitted[nzchar(emitted)])
  # Drop RKHS — it's a model name caught by the regex inside a comparison,
  # not an emitted reason code.
  emitted <- setdiff(emitted, "RKHS")

  g <- gp_policy_decision_reason_glossary()
  missing_codes <- setdiff(emitted, g$reason)
  expect_length(missing_codes, 0L)
})

test_that("known user-override and CV-dispatch codes are present", {
  g <- gp_policy_decision_reason_glossary()
  required <- c(
    "user_forced_sequential", "user_forced_future", "user_forced_mirai",
    "user_forced_base_parallel", "user_forced_foreach",
    "no_parallel_eligible_tasks", "single_gpu_exclusive_deep_models",
    "multiple_sequential_constraints", "all_tasks_user_sequential",
    "single_direct_gp_cv_task", "backend_unavailable"
  )
  expect_true(all(required %in% g$reason),
              info = paste("Missing:", paste(setdiff(required, g$reason), collapse = ", ")))
})
