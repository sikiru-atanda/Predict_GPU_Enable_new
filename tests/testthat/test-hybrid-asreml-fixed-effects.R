test_that("hybrid ASReml fixed effects map ASReml level names onto model.matrix columns", {
  # ASReml names a factor level `<term>_<level>`; model.matrix uses
  # `<term><level>`. A plain name match dropped every environment effect, so
  # hybrid MET predictions all sat on the reference environment's mean.
  ph <- data.frame(Env = c("East", "North", "East", "West"), stringsAsFactors = FALSE)
  terms_fixed <- stats::delete.response(stats::terms(Yield ~ as.factor(`Env`)))
  mm <- stats::model.matrix(terms_fixed, data = ph)
  fixed_vec <- c(
    "(Intercept)" = 10.9,
    "as.factor(Env)_East" = 0,
    "as.factor(Env)_North" = -2.9,
    "as.factor(Env)_West" = 1.5
  )
  beta <- PredictProR:::gp_hybrid_asreml_match_fixed_effects(mm, terms_fixed, fixed_vec)
  expect_equal(as.numeric(mm %*% beta), c(10.9, 8.0, 10.9, 12.4))
})

test_that("hybrid ASReml fixed-effect matching keeps an intercept-only model", {
  ph <- data.frame(Yield = 1:3)
  terms_fixed <- stats::delete.response(stats::terms(Yield ~ 1))
  mm <- stats::model.matrix(terms_fixed, data = ph)
  beta <- PredictProR:::gp_hybrid_asreml_match_fixed_effects(mm, terms_fixed, c("(Intercept)" = 4.2))
  expect_equal(as.numeric(mm %*% beta), rep(4.2, 3))
})
