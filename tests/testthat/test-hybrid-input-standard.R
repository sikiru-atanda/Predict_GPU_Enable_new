hybrid_standard_call <- function(...) {
  if (exists("hybrid_data_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("hybrid_data_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::hybrid_data_standard(...))
  }
  stop("hybrid_data_standard is unavailable for this test.", call. = FALSE)
}

hybrid_validate_call <- function(...) {
  if (exists("validate_hybrid_input_standard", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("validate_hybrid_input_standard", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::validate_hybrid_input_standard(...))
  }
  stop("validate_hybrid_input_standard is unavailable for this test.", call. = FALSE)
}

hybrid_model_execute_call <- function(...) {
  if (exists("model_execute", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("model_execute", envir = .GlobalEnv)(...))
  }
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(PredictProR::model_execute(...))
  }
  stop("model_execute is unavailable for this test.", call. = FALSE)
}

mk_hybrid_contract_data <- function() {
  ph <- data.frame(
    HybridID = c("H1", "H2", "H3"),
    Female = c("F1", "F1", "F2"),
    Male = c("M1", "M2", "M2"),
    Env = c("E1", "E1", "E2"),
    Yield = c(1.2, 2.3, 3.4),
    stringsAsFactors = FALSE
  )

  hybrid_geno <- matrix(
    c(
      0, 1, 2,
      1, 1, 0,
      2, 0, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(hybrid_geno) <- ph$HybridID
  colnames(hybrid_geno) <- paste0("m", 1:3)

  parent_geno <- matrix(
    c(
      0, 1, 1,
      1, 0, 1,
      2, 1, 0,
      1, 2, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(parent_geno) <- c("F1", "F2", "M1", "M2")
  colnames(parent_geno) <- paste0("m", 1:3)

  list(ph = ph, hybrid_geno = hybrid_geno, parent_geno = parent_geno)
}

test_that("hybrid_data_standard returns the package contract", {
  spec <- hybrid_standard_call()
  expect_true(is.list(spec))
  expect_true(all(c("phenotype", "parent_resolved_genomics", "hybrid_level_genomics", "cv_scenarios") %in% names(spec)))
  expect_true(all(c("HybridID", "Female", "Male", "Response") %in% spec$phenotype$required))
})

test_that("hybrid output standard rejects classification families", {
  expect_error(
    prediction_output_standard(response_family = "classification", task = "hybrid"),
    "supports gaussian responses only",
    fixed = TRUE
  )
})

test_that("validate_hybrid_input_standard accepts hybrid-level ML contract", {
  dat <- mk_hybrid_contract_data()
  out <- hybrid_validate_call(
    pheno_data = dat$ph,
    geno_data = dat$hybrid_geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    heter_groups = "Env",
    hybrid_ml = TRUE
  )
  expect_true(isTRUE(out$valid))
  expect_identical(out$mode, "hybrid_ml")
})

test_that("validate_hybrid_input_standard rejects missing hybrid parent columns with standard guidance", {
  dat <- mk_hybrid_contract_data()
  bad_ph <- dat$ph[, c("HybridID", "Female", "Yield"), drop = FALSE]
  expect_error(
    hybrid_validate_call(
      pheno_data = bad_ph,
      geno_data = dat$hybrid_geno,
      response = "Yield",
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      hybrid_ml = TRUE
    ),
    regexp = "Hybrid phenotype standard"
  )
})

test_that("validate_hybrid_input_standard rejects parent-resolved modes without parent-keyed inputs", {
  dat <- mk_hybrid_contract_data()
  expect_error(
    hybrid_validate_call(
      pheno_data = dat$ph,
      geno_data = dat$hybrid_geno,
      response = "Yield",
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      hybrid_asreml = TRUE
    ),
    regexp = "Parent-resolved genomic standard"
  )
})

test_that("model_execute fails early for malformed hybrid input and shows the standard", {
  dat <- mk_hybrid_contract_data()
  bad_geno <- dat$hybrid_geno
  rownames(bad_geno) <- paste0("X", rownames(bad_geno))

  expect_error(
    hybrid_model_execute_call(
      pheno_data = dat$ph,
      geno_data = bad_geno,
      response = "Yield",
      gen_name = "HybridID",
      response_family = "gaussian",
      hybrid_ml = TRUE,
      female_parent = "Female",
      male_parent = "Male",
      GS_model = "Ridge_Regression",
      system_database = TRUE,
      message = FALSE
    ),
    regexp = "Hybrid-level genomic standard"
  )
})

test_that("hybrid MET requires CV0 CV1 or CV2 rather than parent-only CV names", {
  ph <- expand.grid(
    HybridID = paste0("H", 1:4),
    Env = paste0("E", 1:3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ph$Female <- rep(c("F1", "F1", "F2", "F2"), times = 3)
  ph$Male <- rep(c("M1", "M2", "M1", "M2"), times = 3)
  ph$Yield <- seq_len(nrow(ph)) / 10
  geno <- matrix(seq_len(16), nrow = 4)
  rownames(geno) <- paste0("H", 1:4)

  expect_silent(hybrid_validate_call(
    pheno_data = ph,
    geno_data = geno,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    heter_groups = "Env",
    hybrid_ml = TRUE,
    GS_model_cv = "Ridge_Regression",
    cross_validation = TRUE,
    cross_validation_meth = "CV1"
  ))
  expect_error(
    hybrid_validate_call(
      pheno_data = ph,
      geno_data = geno,
      response = "Yield",
      gen_name = "HybridID",
      female_parent = "Female",
      male_parent = "Male",
      heter_groups = "Env",
      hybrid_ml = TRUE,
      GS_model_cv = "Ridge_Regression",
      cross_validation = TRUE,
      cross_validation_meth = "Hybrid_Known_Parents"
    ),
    "requires CV0, CV1, CV2"
  )
})

test_that("all hybrid model families share MET CV0 CV1 and CV2 row scenarios", {
  ph <- expand.grid(
    HybridID = paste0("H", 1:6),
    Env = paste0("E", 1:3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ph$Female <- rep(c("F1", "F1", "F2", "F2", "F3", "F3"), times = 3)
  ph$Male <- rep(c("M1", "M2", "M1", "M2", "M1", "M2"), times = 3)
  ph$Yield <- seq_len(nrow(ph)) / 10
  builder <- if (exists("gp_hybrid_build_cv_scenarios", inherits = TRUE)) {
    get("gp_hybrid_build_cv_scenarios", inherits = TRUE)
  } else {
    get("gp_hybrid_build_cv_scenarios", envir = asNamespace("PredictProR"))
  }

  cv0 <- builder(ph, "Yield", "HybridID", "Female", "Male", "Env", "CV0",
                 nfolds = 3L, random_state = 9L)
  expect_true(all(vapply(cv0, function(x) {
    length(unique(ph$Env[x$test_rows])) == 1L
  }, logical(1))))

  cv1 <- builder(ph, "Yield", "HybridID", "Female", "Male", "Env", "CV1",
                 nfolds = 3L, random_state = 9L)
  expect_true(all(vapply(cv1, function(x) {
    !length(intersect(ph$HybridID[x$test_rows], ph$HybridID[x$train_rows]))
  }, logical(1))))

  cv2 <- builder(ph, "Yield", "HybridID", "Female", "Male", "Env", "CV2",
                 nfolds = 3L, random_state = 9L)
  expect_true(all(vapply(cv2, function(x) {
    length(intersect(ph$HybridID[x$test_rows], ph$HybridID[x$train_rows])) > 0L
  }, logical(1))))
})
