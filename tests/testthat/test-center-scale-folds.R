# The R ML models centred and scaled with caret::preProcess() and made folds
# with caret::createFolds(): loading caret took 1.33 s (61 namespaces) on each
# worker's first ML task on the Windows test machine. The replacements must
# give caret's results exactly, so existing seeds keep their folds.

center_scale_cases <- function() {
  set.seed(4)
  y <- c(stats::rnorm(23, 10, 3), NA, NA)
  X <- matrix(stats::rnorm(40 * 6), 40, 6, dimnames = list(NULL, sprintf("m%d", 1:6)))
  X[, 3] <- 2           # zero variance
  X[5, 2] <- NA
  Xi <- matrix(sample(0:2, 60, TRUE), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
  df <- data.frame(g = letters[1:12], v = stats::rnorm(12), w = 1:12, stringsAsFactors = FALSE)
  list(
    y = as.data.frame(as.matrix(y)),
    one = as.data.frame(as.matrix(5)),
    X = X,
    Xi = Xi,
    df = df
  )
}

test_that("centring and scaling give caret's means, deviations and data", {
  skip_if_not_installed("caret")
  for (nm in names(center_scale_cases())) {
    x <- center_scale_cases()[[nm]]
    for (m in list(c("center", "scale"), "center", "scale")) {
      expected_warn <- NULL
      pp <- withCallingHandlers(
        caret::preProcess(x, method = m),
        warning = function(w) { expected_warn <<- c(expected_warn, conditionMessage(w)); invokeRestart("muffleWarning") }
      )
      got_warn <- NULL
      fit <- withCallingHandlers(
        PredictProR:::gp_center_scale_fit(x, center = "center" %in% m, scale = "scale" %in% m),
        warning = function(w) { got_warn <<- c(got_warn, conditionMessage(w)); invokeRestart("muffleWarning") }
      )
      label <- paste(nm, paste(m, collapse = "+"))
      expect_identical(fit$mean, pp$mean, label = label)
      expect_identical(fit$std, pp$std, label = label)
      expect_identical(got_warn, expected_warn, label = label)
      expect_identical(PredictProR:::gp_center_scale_apply(fit, x), stats::predict(pp, x), label = label)
    }
  }
})

test_that("a scaler from caret, as in a saved model, is applied the same way", {
  skip_if_not_installed("caret")
  x <- center_scale_cases()$X
  pp <- suppressWarnings(caret::preProcess(x, method = c("center", "scale")))
  expect_identical(PredictProR:::gp_center_scale_apply(pp, x), stats::predict(pp, x))
})

test_that("folds are caret's, draw for draw", {
  skip_if_not_installed("caret")
  set.seed(9)
  ys <- list(
    small = stats::rnorm(7),
    medium = stats::rnorm(37),
    large = stats::rnorm(200),
    ties = rep(c(1, 2, 2, 2, 3), 8),
    classes = factor(sample(c("a", "b", "c"), 50, TRUE, prob = c(0.6, 0.3, 0.1))),
    rare = factor(c(rep("x", 30), "y", "y"))
  )
  for (nm in names(ys)) {
    for (k in c(3L, 5L, 10L)) {
      for (opts in list(c(TRUE, FALSE), c(TRUE, TRUE), c(FALSE, FALSE))) {
        set.seed(123)
        expected <- caret::createFolds(ys[[nm]], k = k, list = opts[1], returnTrain = opts[2])
        expected_next <- stats::runif(1)
        set.seed(123)
        got <- PredictProR:::gp_create_folds(ys[[nm]], k = k, list = opts[1], returnTrain = opts[2])
        label <- sprintf("%s k=%d list=%s train=%s", nm, k, opts[1], opts[2])
        expect_identical(got, expected, label = label)
        # the random stream continues where caret's would
        expect_identical(stats::runif(1), expected_next, label = label)
      }
    }
  }
})

test_that("no package code calls caret", {
  files <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(!length(files), "package sources not available")
  hits <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    code <- sub("#.*$", "", lines)
    idx <- grep("caret::", code, fixed = TRUE)
    if (length(idx)) paste0(basename(f), ":", idx) else character()
  }))
  expect_identical(hits, character())
})

test_that("an R ML fold fit does not load caret", {
  skip_on_cran()
  skip_if_not_installed("callr")
  lib <- PredictProR:::gp_predictpror_loaded_lib()
  skip_if(is.null(lib), "PredictProR not loaded from an installed library")
  loaded <- callr::r(function(lib) {
    .libPaths(c(lib, .libPaths()))
    # no Python here: the fit is replaced; the scaling around it is what runs
    utils::assignInNamespace("gp_py_ml_fit_predict", function(..., X_test) rep(0, nrow(X_test)), "PredictProR")
    set.seed(1)
    X <- matrix(stats::rnorm(200), 20, dimnames = list(NULL, paste0("m", 1:10)))
    PredictProR::AI_lasso_cv(y = stats::rnorm(20), omics = X, tst = 1:5)
    isNamespaceLoaded("caret")
  }, args = list(lib = lib))
  expect_false(loaded)
})
