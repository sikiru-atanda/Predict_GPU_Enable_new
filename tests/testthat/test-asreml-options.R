test_that("ASReml option builder preserves user workspace controls", {
  opts <- PredictProR:::asreml_build_options(
    workspace = "12gb",
    pworkspace = "8gb",
    maxit = 75L
  )

  expect_identical(opts$trace, FALSE)
  expect_identical(opts$workspace, "12gb")
  expect_identical(opts$pworkspace, "8gb")
  expect_identical(opts$maxit, 75L)

  opts_null <- PredictProR:::asreml_build_options(
    workspace = NULL,
    pworkspace = NULL,
    maxit = NULL
  )
  expect_identical(opts_null, list(trace = FALSE))
})

test_that("ASReml options set for a fit are restored when the caller exits", {
  skip_if_not_installed("asreml")
  opt_fn <- getFromNamespace("asreml.options", "asreml")
  keys <- c("trace", "maxit", "workspace", "pworkspace")
  before <- tryCatch(opt_fn()[keys], error = function(e) NULL)
  skip_if(is.null(before), "asreml.options() unavailable")

  inside <- NULL
  fit_like <- function() {
    PredictProR:::asreml_apply_options(workspace = 1e8, pworkspace = 1e6, maxit = 77)
    inside <<- opt_fn()[keys]
    invisible(NULL)
  }
  fit_like()

  expect_equal(inside$maxit, 77)
  expect_equal(inside$workspace, 1e8)
  expect_identical(opt_fn()[keys], before)
})

test_that("applying ASReml options twice in one caller still restores the user's values", {
  skip_if_not_installed("asreml")
  opt_fn <- getFromNamespace("asreml.options", "asreml")
  keys <- c("trace", "maxit", "workspace", "pworkspace")
  before <- tryCatch(opt_fn()[keys], error = function(e) NULL)
  skip_if(is.null(before), "asreml.options() unavailable")
  retry_like <- function() {
    PredictProR:::asreml_apply_options(workspace = 1e8, maxit = 50)
    PredictProR:::asreml_apply_options(workspace = 4e8, maxit = 50)   # workspace retry
    invisible(NULL)
  }
  retry_like()
  expect_identical(opt_fn()[keys], before)
})

test_that("ASReml workspace parsing and retry escalation", {
  words <- PredictProR:::gp_asreml_workspace_words
  expect_equal(words(NULL), 1e8)
  expect_equal(words(2e8), 2e8)
  expect_equal(words("800mb"), 800 * 1024^2 / 8)
  expect_equal(words("4gb"), 4 * 1024^3 / 8)
  nxt <- PredictProR:::gp_asreml_next_workspace_words(1e8)
  expect_true(is.null(nxt) || (nxt > 1e8 && nxt <= 4e8))
  local_mocked_bindings(ps_system_memory = function(...) list(avail = 1024^3), .package = "ps")  # 1 GB free
  expect_null(PredictProR:::gp_asreml_next_workspace_words(1e8))   # 0.6 GB cap < 0.8 GB current
})

test_that("ASReml prediction helper forwards workspace controls to predict", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, NA),
    stringsAsFactors = FALSE
  )

  seen <- NULL
  old <- getOption("PredictProR.asreml_predict_impl")
  on.exit(options(PredictProR.asreml_predict_impl = old), add = TRUE)
  options(PredictProR.asreml_predict_impl = function(object, classify, sed, ...) {
    seen <<- list(classify = classify, sed = sed, dots = list(...))
    list(
      pvals = data.frame(
        GID = c("g1", "g2"),
        predicted.value = c(1.1, 2.2),
        std.error = c(0.2, 0.3),
        status = c("Estimable", "Estimable"),
        stringsAsFactors = FALSE
      )
    )
  })

  res <- PredictProR:::asreml_predict_or_extract(
    mod = list(call = list(data = quote(pheno_data))),
    classify = "GID",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    workspace = "12gb",
    pworkspace = "8gb"
  )

  expect_identical(res$mode, "predict")
  expect_identical(seen$classify, "GID")
  expect_identical(seen$sed, FALSE)
  expect_identical(seen$dots$workspace, "12gb")
  expect_identical(seen$dots$pworkspace, "8gb")
  expect_identical(res$prediction$Train_Test_Label, c("Train", "Test"))
})

test_that("ASReml prediction helper passes the model as a local symbol", {
  seen_object_expr <- NULL
  old <- getOption("PredictProR.asreml_predict_impl")
  on.exit(options(PredictProR.asreml_predict_impl = old), add = TRUE)
  options(PredictProR.asreml_predict_impl = function(object, classify, sed, ...) {
    seen_object_expr <<- deparse(substitute(object))
    list(
      pvals = data.frame(
        GID = "g1",
        predicted.value = 1.1,
        std.error = 0.2,
        status = "Estimable",
        stringsAsFactors = FALSE
      )
    )
  })

  pvals <- PredictProR:::asreml_predict_pvals(
    mod = list(call = list(data = quote(pheno_data))),
    classify = "GID"
  )

  expect_equal(pvals$predicted.value, 1.1)
  expect_identical(seen_object_expr, "mod")
})

test_that("ASReml CV route preserves workspace controls", {
  seen_cv <- NULL
  seen_predict <- NULL
  seen_context <- NULL
  seen_data <- NULL
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1, 2),
    stringsAsFactors = FALSE
  )
  had_pheno_dataa <- exists("pheno_dataa", envir = .GlobalEnv, inherits = FALSE)
  old_pheno_dataa <- if (had_pheno_dataa) {
    get("pheno_dataa", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  if (had_pheno_dataa) {
    rm(list = "pheno_dataa", envir = .GlobalEnv)
  }
  withr::defer({
    if (had_pheno_dataa) {
      assign("pheno_dataa", old_pheno_dataa, envir = .GlobalEnv)
    } else if (exists("pheno_dataa", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "pheno_dataa", envir = .GlobalEnv)
    }
  })

  local_mocked_bindings(
    asreml_cv_model = function(pheno_dataa = NULL, workspace = NULL, pworkspace = NULL, maxit = NULL, ...) {
      seen_cv <<- list(workspace = workspace, pworkspace = pworkspace, maxit = maxit)
      list(model_cv = list(call = list(data = quote(pheno_dataa))), pheno_dataa = pheno_dataa)
    },
    asreml_predict_pvals = function(mod, classify, workspace = NULL, pworkspace = NULL) {
      seen_predict <<- list(classify = classify, workspace = workspace, pworkspace = pworkspace)
      seen_context <<- exists("pheno_dataa", envir = .GlobalEnv, inherits = FALSE)
      seen_data <<- if (seen_context) get("pheno_dataa", envir = .GlobalEnv, inherits = FALSE) else NULL
      data.frame(
        GID = c("g1", "g2"),
        predicted.value = c(1.1, 2.2),
        std.error = c(0.2, 0.3),
        status = c("Estimable", "Estimable"),
        stringsAsFactors = FALSE
      )
    },
    .package = "PredictProR"
  )

  pred <- PredictProR::asreml_mod_cv(
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = NULL,
    asreml_models_prep_cv = list(
      code_asr_fit = c("", "", "", ""),
      gen_pos = 1L,
      inter_gen_pos = NULL,
      names_in_inv_list = "gmatrix_inv",
      rand_term = "GID"
    ),
    response = "Yield",
    tst = 2L,
    workspace = "12gb",
    pworkspace = "8gb",
    maxit = 75L
  )

  expect_equal(pred, 2.2)
  expect_identical(seen_cv, list(workspace = "12gb", pworkspace = "8gb", maxit = 75L))
  expect_identical(seen_predict, list(classify = "GID", workspace = "12gb", pworkspace = "8gb"))
  expect_true(seen_context)
  expect_equal(seen_data, pheno)
  expect_false(exists("pheno_dataa", envir = .GlobalEnv, inherits = FALSE))
})

test_that("ASReml CV fits retain fold data across bounded updates", {
  body_text <- paste(deparse(body(PredictProR::asreml_cv_model)), collapse = "\n")
  expect_match(body_text, "asreml_model_data_register\\(pheno_dataa\\)")
  expect_match(body_text, "asreml_model_data_lookup")
  expect_match(body_text, "mod_cv\\$call\\$data <- parse")

  fold_data <- data.frame(
    GID = factor(c("g1", "g2")),
    Env = factor(c("E1", "E2")),
    Yield = c(1.2, NA_real_),
    stringsAsFactors = TRUE
  )
  cache_id <- PredictProR:::asreml_model_data_register(fold_data)
  cached <- PredictProR:::asreml_model_data_lookup(cache_id)
  expect_equal(cached, fold_data)
})

test_that("predict_with_model forwards ASReml workspace controls", {
  seen <- NULL
  local_mocked_bindings(
    asreml_mod_cv = function(workspace = NULL, pworkspace = NULL, maxit = NULL, ...) {
      seen <<- list(workspace = workspace, pworkspace = pworkspace, maxit = maxit)
      1
    },
    .package = "PredictProR"
  )

  pred <- PredictProR::predict_with_model(
    model = "GBLUP",
    tst = 1L,
    additional_params = list(
      asreml_models_prep_cv = list(),
      pheno_data = data.frame(GID = "g1", Yield = 1),
      response = "Yield",
      heter_groups = NULL,
      gen_name = "GID",
      workspace = "12gb",
      pworkspace = "8gb",
      maxit = 75L
    )
  )

  expect_equal(pred, 1)
  expect_identical(seen, list(workspace = "12gb", pworkspace = "8gb", maxit = 75L))
})

test_that("ASReml exported signatures expose implemented workspace controls", {
  expect_true(all(c("workspace", "pworkspace", "maxit") %in% names(formals(PredictProR::asreml_cv_model))))
  expect_true(all(c("workspace", "pworkspace", "maxit") %in% names(formals(PredictProR::asreml_mod_cv))))
  expect_true(all(c("workspace", "pworkspace", "maxit") %in% names(formals(PredictProR::asreml_mod_output_new))))
})

test_that("ASReml coefficient recovery supports singular kernels", {
  kernel <- matrix(
    c(1, 1, 1, 1),
    nrow = 2L,
    dimnames = list(c("g1", "g2"), c("g1", "g2"))
  )
  ebv <- c(2, 2)

  out <- cal_coeff_asreml(
    gmatrix = kernel,
    ebv = ebv,
    heter_groups = NULL,
    heter_grp = NULL,
    gid_name = rownames(kernel)
  )

  expect_identical(out$x_variables, c("g1", "g2"))
  expect_equal(unname(drop(kernel %*% out$coeff)), ebv, tolerance = 1e-10)
})

test_that("ASReml output names use input roles rather than label substrings", {
  ids <- c("g1", "g2")
  make_kernel <- function(source) {
    value <- diag(2L)
    dimnames(value) <- list(ids, ids)
    attr(value, "predictpror_source_kernel") <- source
    value
  }

  expect_identical(
    asreml_kernel_output_names(
      gmatrix = make_kernel("Genomic_layer_A"),
      omic1_kernel = make_kernel("Sunflower_layer_B"),
      omic2_kernel = make_kernel("Sunflower_layer_C")
    ),
    c("gmatrix", "omic1_kernel", "omic2_kernel")
  )
})

test_that("ASReml component matching does not collide on name prefixes", {
  rows <- c(
    "vm(GID, gmatrix_inv)",
    "vm(GID, kernel_list_model_ready)",
    "vm(GID, kernel_list_model_ready_Matern32)",
    "units!R"
  )

  expect_identical(
    asreml_single_kernel_varcomp_rows(
      rows, "GID", "kernel_list_model_ready"
    ),
    2L
  )
  expect_identical(
    asreml_single_kernel_varcomp_rows(
      rows, "GID", "kernel_list_model_ready_Matern32"
    ),
    3L
  )
})
