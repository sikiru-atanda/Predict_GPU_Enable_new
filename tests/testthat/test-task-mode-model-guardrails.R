mk_task_guardrail_pheno <- function(met = FALSE) {
  base <- data.frame(
    GID = paste0("g", seq_len(5)),
    T1 = c(1, 2, NA, 4, 5),
    T2 = c(2, NA, 3, 4, 6),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(met)) {
    return(base)
  }
  out <- rbind(
    transform(base, Env = "E1"),
    transform(base, Env = "E2")
  )
  rownames(out) <- NULL
  out
}

mk_task_guardrail_kernel <- function(ids) {
  K <- diag(length(ids))
  rownames(K) <- colnames(K) <- ids
  K
}

test_that("model_execute exposes explicit multi-trait route flags", {
  expect_true("multi_trait_gp" %in% names(formals(PredictProR::model_execute)))
  expect_true("multi_trait_bayes" %in% names(formals(PredictProR::model_execute)))

  spec <- PredictProR::multi_trait_data_standard()
  expect_true("multi_trait_gp_models" %in% names(spec$current_scope))
  expect_true("multi_trait_bayes_models" %in% names(spec$current_scope))
  expect_identical(
    spec$current_scope$multi_trait_gp_models,
    c("Gaussian-Process-GBLUP", "FA-GBLUP", "Scalable-GBLUP")
  )
  expect_false("Kernel-GBLUP" %in% spec$current_scope$multi_trait_gp_models)
  expect_true(all(c("GBLUP_BRR", "RKHS") %in% spec$current_scope$multi_trait_bayes_models))
})

test_that("multi-trait GP model labels select purpose-specific backends", {
  gp <- PredictProR:::gp_multitrait_gp_route_spec("Gaussian-Process-GBLUP")
  fa <- PredictProR:::gp_multitrait_gp_route_spec("FA-GBLUP", fa_rank = 1L)
  scalable <- PredictProR:::gp_multitrait_gp_route_spec("Scalable-GBLUP")
  fa_met <- PredictProR:::gp_multitrait_gp_route_spec(
    "FA-GBLUP", is_met = TRUE, fa_rank = 1L
  )

  expect_identical(gp$backend_method, "joint_gp_ai_reml")
  expect_identical(gp$trait_structure, "unstructured")
  expect_identical(fa$backend_method, "joint_gp_fa_ai_reml")
  expect_identical(fa$trait_structure, "fa")
  expect_identical(fa$trait_fa_rank, 1L)
  expect_identical(scalable$backend_method, "joint_gp_operator_mom")
  expect_identical(scalable$varcomp_mode, "mom")
  expect_identical(fa_met$backend_method, "joint_gp_fa_mom")
  expect_identical(fa_met$gxe_trait_structure, "fa")

  expect_error(
    PredictProR:::gp_multitrait_gp_route_spec("Kernel-GBLUP"),
    "Kernel-GBLUP is not a multi-trait model"
  )
})

test_that("single-trait GP scope removes only the exact Scalable alias", {
  expect_identical(
    PredictProR:::gp_single_response_single_environment_gp_supported_models(),
    c("KRR", "GP")
  )
  expect_silent(PredictProR:::gp_validate_single_trait_gp_model_scope(
    GS_model = c("Kernel-GBLUP", "Gaussian-Process-GBLUP")
  ))
  expect_error(
    PredictProR:::gp_validate_single_trait_gp_model_scope(
      GS_model = "Scalable-GBLUP"
    ),
    "resolves to Kernel-GBLUP.*identical predictions"
  )
  expect_silent(PredictProR:::gp_validate_single_trait_gp_model_scope(
    GS_model = "Scalable-GBLUP",
    multi_trait_gp = TRUE
  ))
  expect_identical(PredictProR:::gp_hybrid_gp_supported_models(), "GP")
  expect_identical(
    PredictProR:::gp_single_trait_multi_environment_gp_supported_models(),
    c("KRR", "GP", "GP_FA")
  )
  expect_false("LowRankGP" %in%
                 PredictProR:::gp_multi_environment_kernel_models())
})

test_that("model_execute rejects the redundant Scalable single-trait choice before fit", {
  ph <- mk_task_guardrail_pheno()
  K <- mk_task_guardrail_kernel(ph$GID)

  expect_error(
    PredictProR::model_execute(
      pheno_data = ph,
      gmatrix = K,
      response = "T1",
      gen_name = "GID",
      GS_model = "Scalable-GBLUP",
      eval_metrics = "mean_squared_error",
      system_database = TRUE,
      message = FALSE
    ),
    "not supported as a separate single-trait GP model"
  )
})

test_that("ordinary multiple response remains independent per-trait mode", {
  ph <- mk_task_guardrail_pheno()
  K <- mk_task_guardrail_kernel(ph$GID)

  out <- PredictProR:::validate_multi_trait_input_standard(
    pheno_data = ph,
    gmatrix = K,
    response = c("T1", "T2"),
    gen_name = "GID",
    GS_model = "Gaussian-Process-GBLUP",
    multi_trait_gp = FALSE
  )

  expect_null(out$mode)
})

test_that("model_execute only routes joint GP when multi_trait_gp is explicit", {
  ph <- mk_task_guardrail_pheno()
  K <- mk_task_guardrail_kernel(ph$GID)
  captured <- list()

  local_mocked_bindings(
    gp_prepare_model_input_objects = function(ctx) {
      list(
        low_call_rate_inds_removed = NULL,
        pheno_clean = list(pheno_clean_data = ph, test_set = NULL),
        geno_res = list(),
        omic1_res = list(),
        omic2_res = list(),
        omic3_res = list(),
        geno_omic_model_ready_list = list(),
        gmatrix_kernel_model_ready_list = list(gmatrix_model_ready = K),
        test_set = NULL,
        ml_dat_res = list(
          pheno_clean_data = ph,
          merged_data = list(merge_data = K),
          omic_count = 0L
        )
      )
    },
    gp_route_true_prediction_specialized_models = function(ctx) {
      captured[[length(captured) + 1L]] <<- isTRUE(ctx$gp_public_multitrait_run)
      list(gp_public_multitrait_run = isTRUE(ctx$gp_public_multitrait_run))
    },
    .package = "PredictProR"
  )

  ordinary <- PredictProR::model_execute(
    pheno_data = ph,
    gmatrix = K,
    response = c("T1", "T2"),
    gen_name = "GID",
    GS_model = "Gaussian-Process-GBLUP",
    eval_metrics = "mean_squared_error",
    system_database = TRUE,
    message = FALSE
  )
  explicit <- PredictProR::model_execute(
    pheno_data = ph,
    gmatrix = K,
    response = c("T1", "T2"),
    gen_name = "GID",
    GS_model = "Gaussian-Process-GBLUP",
    multi_trait_gp = TRUE,
    eval_metrics = "mean_squared_error",
    system_database = TRUE,
    message = FALSE
  )

  expect_false(ordinary$gp_public_multitrait_run)
  expect_true(explicit$gp_public_multitrait_run)
  expect_identical(unlist(captured), c(FALSE, TRUE))
})

test_that("multi-trait GP rejects non-GP models before dispatch", {
  ph <- mk_task_guardrail_pheno()
  K <- mk_task_guardrail_kernel(ph$GID)

  expect_error(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      gmatrix = K,
      response = c("T1", "T2"),
      gen_name = "GID",
      GS_model = "RandomForest",
      multi_trait_gp = TRUE
    ),
    "Multi-trait GP.*Gaussian-Process-GBLUP.*FA-GBLUP.*Scalable-GBLUP"
  )

  expect_error(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      gmatrix = K,
      response = c("T1", "T2"),
      gen_name = "GID",
      GS_model = "Kernel-GBLUP",
      multi_trait_gp = TRUE
    ),
    "Multi-trait GP.*supports only:.*Gaussian-Process-GBLUP.*FA-GBLUP.*Scalable-GBLUP"
  )

  for (model in c("GBLUP_BRR", "RKHS")) {
    expect_error(
      PredictProR:::validate_multi_trait_input_standard(
        pheno_data = ph,
        gmatrix = K,
        response = c("T1", "T2"),
        gen_name = "GID",
        GS_model = model,
        multi_trait_gp = TRUE
      ),
      "Multi-trait GP.*Gaussian-Process-GBLUP.*FA-GBLUP.*Scalable-GBLUP"
    )
  }
})

test_that("supported multi-trait GP accepts raw genomic data or user kernels", {
  ph <- mk_task_guardrail_pheno(met = TRUE)
  ids <- unique(ph$GID)
  geno <- matrix(
    seq_len(length(ids) * 3),
    nrow = length(ids),
    dimnames = list(ids, paste0("m", 1:3))
  )
  K <- mk_task_guardrail_kernel(ids)

  expect_silent(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      geno_data = geno,
      response = c("T1", "T2"),
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = "Gaussian-Process-GBLUP",
      multi_trait_gp = TRUE
    )
  )
  expect_silent(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      gkernel = K,
      response = c("T1", "T2"),
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = "Gaussian-Process-GBLUP",
      multi_trait_gp = TRUE
    )
  )
})

test_that("multi-trait multi-environment currently allows only GP mode", {
  ph <- mk_task_guardrail_pheno(met = TRUE)
  K <- mk_task_guardrail_kernel(unique(ph$GID))

  expect_error(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      geno_data = matrix(1, nrow = length(unique(ph$GID)), ncol = 2,
                         dimnames = list(unique(ph$GID), c("m1", "m2"))),
      response = c("T1", "T2"),
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = "RandomForest",
      multi_trait_ml = TRUE
    ),
    "multi-trait multi-environment.*multi_trait_gp"
  )

  expect_silent(
    PredictProR:::validate_multi_trait_input_standard(
      pheno_data = ph,
      gmatrix = K,
      response = c("T1", "T2"),
      gen_name = "GID",
      heter_groups = "Env",
      GS_model = "Gaussian-Process-GBLUP",
      multi_trait_gp = TRUE
    )
  )

  for (model in c("GBLUP_BRR", "RKHS")) {
    expect_error(
      PredictProR:::validate_multi_trait_input_standard(
        pheno_data = ph,
        gmatrix = K,
        response = c("T1", "T2"),
        gen_name = "GID",
        heter_groups = "Env",
        GS_model = model,
        multi_trait_gp = TRUE
      ),
      "multi-trait multi-environment prediction.*Gaussian-Process-GBLUP.*FA-GBLUP.*Scalable-GBLUP"
    )
  }
})

test_that("hybrid validators enforce task-specific model lists", {
  ph <- data.frame(
    HybridID = c("H1", "H2", "H3"),
    Female = c("F1", "F1", "F2"),
    Male = c("M1", "M2", "M1"),
    Yield = c(1, NA, 3),
    Moisture = c(4, 5, NA),
    stringsAsFactors = FALSE
  )
  parent_ids <- c("F1", "F2", "M1", "M2")
  parent_X <- matrix(1, nrow = length(parent_ids), ncol = 2,
                     dimnames = list(parent_ids, c("m1", "m2")))
  hybrid_X <- matrix(1, nrow = nrow(ph), ncol = 2,
                     dimnames = list(ph$HybridID, c("m1", "m2")))

  expect_error(
    PredictProR:::validate_hybrid_input_standard(
      pheno_data = ph,
      geno_data = parent_X,
      response = "Yield",
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      hybrid_gp = TRUE,
      GS_model = "RandomForest"
    ),
    "Hybrid GP.*Gaussian-Process-GBLUP"
  )

  expect_error(
    PredictProR:::validate_hybrid_input_standard(
      pheno_data = ph,
      geno_data = hybrid_X,
      response = c("Yield", "Moisture"),
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      hybrid_ml = TRUE,
      GS_model = "RandomForest"
    ),
    "Hybrid ML.*one response"
  )
})
