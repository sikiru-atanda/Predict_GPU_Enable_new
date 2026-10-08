# TDD (audit finding): the reliability column silently mixes two definitions.
# Primary is the standard r2 = 1 - PEV / sigma2_g. When that degenerates
# (PEV >= sigma2_g, so the clamped values go flat), a different ratio form
# sigma2_g / (sigma2_g + PEV) is substituted. That switch must be flagged so the
# mixed definition is transparent: gp_met_empirical_reliability now attaches a
# logical `ratio_fallback` attribute marking the rows where the fallback applied.

test_that("ratio-form reliability fallback is flagged via attribute", {
  if (!exists("gp_met_empirical_reliability", mode = "function")) {
    skip("gp_met_empirical_reliability not available")
  }
  # PEV >= reference variance everywhere -> primary form clamps flat at 0 while
  # PEV varies -> the ratio fallback kicks in for all (valid) rows.
  rel <- gp_met_empirical_reliability(c(2, 3, 4), c(1, 1, 1))
  fb <- attr(rel, "ratio_fallback")
  expect_false(is.null(fb))
  expect_length(fb, 3L)
  expect_true(all(fb))
  expect_equal(as.numeric(rel), c(1 / 3, 1 / 4, 1 / 5), tolerance = 1e-8)
})

test_that("standard reliability does not flag the fallback", {
  if (!exists("gp_met_empirical_reliability", mode = "function")) {
    skip("gp_met_empirical_reliability not available")
  }
  rel <- gp_met_empirical_reliability(c(0.1, 0.2, 0.3), 1.0)
  fb <- attr(rel, "ratio_fallback")
  expect_false(is.null(fb))
  expect_false(any(fb))
})

test_that("GP prediction contract surfaces the ratio-fallback remark", {
  if (!exists("gp_format_gaussian_prediction_table", mode = "function")) {
    skip("gp_format_gaussian_prediction_table not available")
  }
  # PEV >= reference (Genetic_variance) everywhere with varying PEV -> primary
  # reliability clamps flat -> ratio fallback -> remark must be surfaced.
  df <- data.frame(
    GID = c("a", "b", "c", "d"),
    Predicted_value = c(1.0, 2.0, 3.0, 4.0),
    PEV = c(2.0, 3.0, 4.0, 5.0),
    Genetic_variance = 1.0,
    stringsAsFactors = FALSE
  )
  out <- gp_format_gaussian_prediction_table(df)
  expect_true("Reliability_remarks" %in% names(out))
  expect_true(any(grepl("ratio-form fallback used", out[["Reliability_remarks"]], fixed = TRUE)))
})
