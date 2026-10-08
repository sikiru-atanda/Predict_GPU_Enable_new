# SVM fits use a BLAS-precomputed kernel (~1000x faster at 1.4k x 10k markers).
# It must reproduce libsvm's own kernel evaluation for every kernel and task.
test_that("precomputed-kernel SVM matches libsvm's own kernel for every kernel and task", {
  skip_if_no_ml_python()
  withr::local_envvar(c(PREDICTPRO_ML_PAYLOAD_CACHE = "false"))
  set.seed(21)
  n <- 80L; p <- 30L
  X <- matrix(rnorm(n * p), n)
  X_test <- matrix(rnorm(15L * p), 15L)
  y_num <- as.numeric(X[, 1:3] %*% c(1, -0.5, 0.3) + rnorm(n, sd = 0.3))
  y_bin <- factor(ifelse(y_num > 0, "yes", "no"), levels = c("no", "yes"))
  y_multi <- factor(cut(y_num, 3, labels = c("A", "B", "C")))
  fit <- function(y, fam, params, max_n) {
    withr::with_envvar(c(PREDICTPRO_SVM_PRECOMPUTE_MAX_N = max_n), {
      PredictProR:::gp_py_ml_fit_predict(
        model_type = "svm", X_train = X, y_train = y, X_test = X_test,
        response_family = fam, model_params = params
      )
    })
  }
  for (kernel in c("Linear", "Gaussian", "Polynomial", "Sigmoid")) {
    params <- list(svm_kernel = kernel, C_value = 1, degree_value = 2, offset_value = 0.5)
    libsvm <- fit(y_num, "gaussian", params, "0")
    blas <- fit(y_num, "gaussian", params, "15000")
    expect_equal(as.numeric(blas), as.numeric(libsvm), tolerance = 1e-8, info = paste(kernel, "eps-regression"))
    nu <- utils::modifyList(params, list(svm_type = "nu-regression"))
    expect_equal(as.numeric(fit(y_num, "gaussian", nu, "15000")), as.numeric(fit(y_num, "gaussian", nu, "0")),
                 tolerance = 1e-8, info = paste(kernel, "nu-regression"))
    for (y in list(binary = y_bin, multiclass = y_multi)) {
      fam <- if (nlevels(y) == 2L) "binary" else "multiclass"
      a <- fit(y, fam, params, "0"); b <- fit(y, fam, params, "15000")
      if (is.numeric(a)) {   # binary: the prediction is the positive-class probability
        expect_equal(as.numeric(b), as.numeric(a), tolerance = 1e-6, info = paste(kernel, fam, "prediction"))
      } else {
        expect_identical(as.character(b), as.character(a), info = paste(kernel, fam, "classes"))
      }
      expect_equal(attr(b, "probabilities"), attr(a, "probabilities"), tolerance = 1e-6, info = paste(kernel, fam, "probabilities"))
    }
  }
})
