test_that("compute_inverse_and_sparse returns inverse triplet when requested", {
  kernel <- matrix(c(2, 0.5, 0.5, 1.5), nrow = 2, byrow = TRUE)
  rownames(kernel) <- colnames(kernel) <- c("A", "B")

  inv_triplet <- PredictProR:::compute_inverse_and_sparse(kernel = kernel, inverse = TRUE)
  expected_inv <- chol2inv(chol(kernel))
  expected_triplet <- PredictProR:::sparse_matrix(expected_inv)

  expect_true(is.matrix(inv_triplet))
  expect_equal(dim(inv_triplet), c(3L, 3L))
  expect_identical(attr(inv_triplet, "rowNames"), c("A", "B"))
  expect_identical(attr(inv_triplet, "colNames"), c("A", "B"))
  expect_true(isTRUE(attr(inv_triplet, "INVERSE")))
  expect_equal(dim(inv_triplet), dim(expected_triplet))
  expect_equal(as.numeric(inv_triplet), as.numeric(expected_triplet), tolerance = 1e-5)
})

test_that("multi-kernel MET builder applies the requested structure to every kernel", {
  pheno_data <- data.frame(
    GID = factor(c("g1", "g1", "g2", "g2")),
    Env = factor(c("E1", "E2", "E1", "E2")),
    Yield = c(1, 2, 3, 4)
  )

  code_asr <- c("asreml::asreml(fixed=trait~1+Env)", "random=~", "residual=~")

  res <- PredictProR:::random_terms_fit_new(
    random = ~ GID + GID:Env,
    fixed = ~ Env,
    fixed_term = "Env",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "corgh",
    code_asr = code_asr,
    names_in_inv_list = c("gmatrix_inv", "omics_inv"),
    gen_name = "GID",
    pheno_data = pheno_data
  )

  expect_match(res$code_asr[[2]], "corgh\\(Env\\):vm\\(GID,gmatrix_inv\\)")
  expect_match(res$code_asr[[2]], "corgh\\(Env\\):vm\\(GID,omics_inv\\)")
  expect_false(grepl("idv\\(Env\\):vm\\(GID,omics_inv\\)", res$code_asr[[2]]))
})

test_that("ASReml MET compound-symmetry builder keeps a model term when Env is fixed", {
  pheno_data <- data.frame(
    GID = factor(c("g1", "g1", "g2", "g2")),
    Env = factor(c("E1", "E2", "E1", "E2")),
    Yield = c(1, 2, 3, 4)
  )

  code_asr <- c("asreml::asreml(fixed=trait~1+Env)", "random=~", "residual=~")

  expect_no_error(
    res <- PredictProR:::random_terms_fit_new(
      random = ~ GID + GID:Env,
      fixed = ~ Env,
      fixed_term = "Env",
      heter_groups = "Env",
      heter_resid = FALSE,
      var_cov_str = NULL,
      code_asr = code_asr,
      names_in_inv_list = "gmatrix_inv",
      gen_name = "GID",
      pheno_data = pheno_data
    )
  )

  expect_match(res$code_asr[[2]], "idv\\(Env\\):vm\\(GID,gmatrix_inv\\)")
  expect_false(grepl("random=~\\s*$", res$code_asr[[2]]))
})

test_that("ASReml utility coerces genotype IDs to factor before building vm terms", {
  pheno_data <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  gmatrix <- diag(2)
  rownames(gmatrix) <- colnames(gmatrix) <- c("g1", "g2")

  seen_gid_factor <- FALSE

  local_mocked_bindings(
    random_terms_fit_new = function(pheno_data = NULL, gen_name = NULL, ...) {
      seen_gid_factor <<- is.factor(pheno_data[[gen_name]])
      list(
        code_asr = c(
          "asreml::asreml(fixed=trait~1+Env)",
          "random=~idv(Env):vm(GID,gmatrix_inv)",
          "residual=~id(units)"
        ),
        gen_pos = 1L,
        inter_gen_pos = 2L,
        rand_term = c("GID", "GID:Env")
      )
    },
    .package = "PredictProR"
  )

  PredictProR:::asreml_utilis_new(
    fixed = ~ Env,
    random = ~ GID + GID:Env,
    GS_model = "GBLUP",
    response = "Yield",
    pheno_data = pheno_data,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = FALSE,
    var_cov_str = NULL,
    inverse = TRUE,
    engine = "stats",
    cross_validation = TRUE
  )

  expect_true(seen_gid_factor)
})

test_that("ASReml compound-symmetry heritability handles homogeneous residual variance", {
  model <- list(
    mf = data.frame(
      Env = factor(c("E1", "E2", "E3", "E1", "E2", "E3"))
    )
  )
  class(model) <- "mock_asreml"

  vc <- data.frame(
    component = c(0.4, 0.2),
    std.error = c(0.05, 0.03),
    bound = c("P", "P"),
    row.names = c("Env:vm(GID, gmatrix_inv)!Env", "units!R"),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    asreml_varcomp_table = function(...) vc,
    .package = "PredictProR"
  )

  expect_no_error(
    out <- PredictProR:::asreml_herit_CSM_new(
      model = model,
      heter_groups = "Env",
      heter_resid = FALSE,
      names_in_inv_list = "gmatrix_inv",
      inter_gen_pos = 2L,
      gen_pos = 1L
    )
  )

  expect_equal(names(out$Residual_Var), c("E1", "E2", "E3"))
  expect_equal(unname(out$Residual_Var), rep(0.2, 3))
  expect_equal(colnames(out$Total_genetic_var), c("E1", "E2", "E3"))
  expect_true(all(is.finite(out$Heritability)))
})

test_that("task context keeps omitted ASReml variance structure as explicit NULL", {
  parent_env <- list2env(
    list(
      var_cov_str = NULL,
      var_cov_str_available = c("us", "corgh", "corh")
    ),
    parent = emptyenv()
  )
  local_env <- list2env(
    list(
      var_cov_str = NULL,
      response = "Yield"
    ),
    parent = parent_env
  )

  ctx <- PredictProR:::gp_merge_task_context(parent_env, local_env)

  expect_true("var_cov_str" %in% names(ctx))
  expect_null(ctx[["var_cov_str"]])
  expect_null(ctx$var_cov_str)
})

test_that("obsolete ASReml corgv spelling is rejected before fitting", {
  expected <- "Use `var_cov_str = \"corgh\"`"

  expect_error(
    PredictProR:::gp_reject_obsolete_asreml_structure("corgv"),
    expected,
    fixed = TRUE
  )
  expect_error(
    random_terms_fit_new(var_cov_str = "CORGV"),
    expected,
    fixed = TRUE
  )
  expect_error(
    model_execute(var_cov_str = "corgv"),
    expected,
    fixed = TRUE
  )
  expect_invisible(
    PredictProR:::gp_reject_obsolete_asreml_structure("corgh")
  )
})

test_that("ASReml CV0 rejects fixed held-out environment means before returning empty CV output", {
  expect_error(
    PredictProR:::gp_guard_asreml_cv0_fixed_environment(
      cross_validation = TRUE,
      cross_validation_meth = "CV0",
      GS_model_cv = "GBLUP",
      engine = "asreml",
      fixed = ~ Env,
      heter_groups = "Env"
    ),
    "fixed environment means are not estimable"
  )

  expect_no_error(
    PredictProR:::gp_guard_asreml_cv0_fixed_environment(
      cross_validation = TRUE,
      cross_validation_meth = "CV1",
      GS_model_cv = "GBLUP",
      engine = "asreml",
      fixed = ~ Env,
      heter_groups = "Env"
    )
  )
})
