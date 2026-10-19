# ft_learning_rate: TabTransformer's learning rate, settable from model_execute()
# and forwarded on every route (single-environment, MET ML/DL, multi-trait,
# hybrid; true prediction and cross-validation).

test_that("ft_learning_rate has one default everywhere it is declared", {
  canon_default <- PredictProR:::gp_dl_bridge_fit_params("ft_transformer", list())$learning_rate
  expect_equal(canon_default, 3e-4)
  expect_equal(eval(formals(PredictProR::model_execute)$ft_learning_rate), canon_default)
  for (fn in c("deep_learning_model", "gp_met_deep_learning_model",
               "gp_multitrait_dl_gaussian_model", "gp_multitrait_dl_gaussian_cv",
               "gp_hybrid_dl_fit_response", "gp_hybrid_dl_gaussian_model", "gp_hybrid_dl_gaussian_cv")) {
    f <- get(fn, envir = asNamespace("PredictProR"))
    expect_true("ft_learning_rate" %in% names(formals(f)), info = fn)
    expect_equal(eval(formals(f)$ft_learning_rate), canon_default, info = fn)
  }
})

test_that("an explicit ft_learning_rate reaches the TabTransformer fit parameters", {
  expect_equal(PredictProR:::gp_dl_bridge_fit_params("ft_transformer", list(ft_learning_rate = 7e-4))$learning_rate, 7e-4)
  met <- PredictProR:::gp_met_dl_model_params("ft_transformer", list(ft_learning_rate = 7e-4))
  expect_equal(met$ft_learning_rate, 7e-4)
})

test_that("CV dispatch fits the requested deep-learning model and forwards ft_learning_rate", {
  # model_execute(model_type = "resnet") is a separate, global setting: CV used to
  # pass it on, so every single-environment DL model was fitted as a ResNet
  seen <- list()
  local_mocked_bindings(deep_learning_model = function(..., model_type, ft_learning_rate) {
    seen <<- list(model_type = model_type, ft_learning_rate = ft_learning_rate)
    numeric(0)
  })
  for (m in c("ft_transformer", "mlp", "cnn", "tabnet")) {
    PredictProR::predict_with_model(m, y = 1:4, omics_data = matrix(0, 4, 2), tst = 1L,
                                    additional_params = list(model_type = "resnet", ft_learning_rate = 7e-4))
    expect_identical(seen$model_type, m)
    expect_equal(seen$ft_learning_rate, 7e-4)
  }
})
