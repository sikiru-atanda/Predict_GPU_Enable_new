reference_scale_matrix <- function(x, scaling = TRUE, centering = FALSE) {
  x <- as.matrix(x)

  if (isTRUE(scaling)) {
    x <- scale(x, center = TRUE, scale = TRUE)
  } else if (isTRUE(centering)) {
    x <- scale(x, center = TRUE, scale = FALSE)
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop("scaled matrix contains missing or non-finite values", call. = FALSE)
  }

  x
}

reference_kernel <- function(x, method, theta = 1, alpha = 0.5) {
  x <- as.matrix(x)
  linear <- tcrossprod(x) / ncol(x)
  d2 <- as.matrix(stats::dist(x))^2
  gaussian <- exp(-theta * d2 / stats::median(d2))

  switch(
    method,
    Gaussian_kernel = gaussian,
    Linear_kernel = linear,
    Composite_kernel = alpha * linear + (1 - alpha) * gaussian,
    Poly2_kernel = (linear + 1)^2,
    Poly3_kernel = (linear + 1)^3,
    Poly4_kernel = (linear + 1)^4,
    stop("unsupported reference kernel method", call. = FALSE)
  )
}

reference_grm <- function(geno, method, weight = NULL) {
  geno <- as.matrix(geno)
  freq <- colMeans(geno) / 2
  q <- 1 - freq
  centered <- scale(geno, center = TRUE, scale = FALSE)
  denom <- sum(2 * freq * (1 - freq))

  vanraden <- tcrossprod(centered) / denom
  dominance_vitezica <- function() {
    w <- matrix(0, nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
    for (k in seq_len(ncol(geno))) {
      w[geno[, k] == 0, k] <- -2 * freq[[k]]^2
      w[geno[, k] == 1, k] <- 2 * freq[[k]] * q[[k]]
      w[geno[, k] == 2, k] <- -2 * q[[k]]^2
    }
    out <- tcrossprod(w) / sum((2 * freq * q)^2)
    dimnames(out) <- list(rownames(geno), rownames(geno))
    out
  }
  dominance_su <- function() {
    h <- matrix(0, nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
    for (k in seq_len(ncol(geno))) {
      h[geno[, k] == 0, k] <- -2 * freq[[k]] * q[[k]]
      h[geno[, k] == 1, k] <- 1 - 2 * freq[[k]] * q[[k]]
      h[geno[, k] == 2, k] <- -2 * freq[[k]] * q[[k]]
    }
    out <- tcrossprod(h) / sum(2 * freq * q * (1 - 2 * freq * q))
    dimnames(out) <- list(rownames(geno), rownames(geno))
    out
  }

  out <- switch(
    method,
    VanRaden = vanraden,
    Weighted_VanRaden = {
      if (is.null(weight)) {
        stop("weight is required", call. = FALSE)
      }
      weight <- as.matrix(weight)
      if (nrow(weight) == ncol(weight)) {
        w <- diag(weight)
      } else {
        w <- weight[, 1]
      }
      tcrossprod(sweep(centered, 2, w, `*`), centered) / denom
    },
    Epistasis = vanraden * vanraden,
    Yang = {
      n_marker <- ncol(geno)
      n_individuals <- nrow(geno)
      locus <- (1 / n_marker) *
        (centered %*% (t(centered) * (1 / (2 * freq * (1 - freq)))))
      locus[lower.tri(locus, diag = TRUE)] <- 0
      locus <- locus + t(locus)

      multiplier <- geno^2 -
        t(t(geno) * (1 + 2 * freq)) +
        matrix(rep(2 * freq^2, each = n_individuals), ncol = n_marker)
      diag(locus) <- 1 +
        (1 / n_marker) * colSums(t(multiplier) * (1 / (2 * freq * (1 - freq))))
      locus
    },
    Dominance = dominance_vitezica(),
    Dominance_Vitezica = dominance_vitezica(),
    Dominance_Su = dominance_su(),
    Dominance_Heterozygosity = dominance_su(),
    stop("unsupported reference GRM method", call. = FALSE)
  )
  attr(out, "ploidy") <- 2L
  out
}

expect_named_relationship_matrix <- function(actual, ids) {
  expect_identical(rownames(actual), ids)
  expect_identical(colnames(actual), ids)
}

test_that("native dense kernel backend matches R formulas", {
  x <- matrix(
    c(
      0.25, 1.50, -0.75,
      1.25, 0.50, 0.00,
      -0.50, 2.00, 1.25,
      0.75, -1.00, 0.50
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:4), paste0("m", 1:3))
  )

  for (method in c(
    "Gaussian_kernel",
    "Linear_kernel",
    "Composite_kernel",
    "Poly2_kernel",
    "Poly3_kernel",
    "Poly4_kernel"
  )) {
    actual <- kernel_calculation(
      M_matrix_clean = x,
      scaling = FALSE,
      theta = 0.75,
      alpha = 0.35,
      method = method,
      backend = "cpp",
      message = FALSE
    )
    expected <- reference_kernel(x, method, theta = 0.75, alpha = 0.35)

    expect_equal(as.matrix(actual), expected, tolerance = 1e-10)
    expect_named_relationship_matrix(actual, rownames(x))
  }
})

test_that("native dense GRM backend matches R formulas", {
  geno <- matrix(
    c(
      0, 1, 2, 0,
      1, 2, 0, 1,
      2, 0, 1, 2,
      0, 1, 1, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:4), paste0("snp", 1:4))
  )
  weights <- matrix(
    c(1.0, 0.7, 1.4, 0.9),
    ncol = 1,
    dimnames = list(colnames(geno), "weight")
  )
  diagonal_weights <- diag(as.vector(weights[, 1]))
  rownames(diagonal_weights) <- colnames(geno)
  colnames(diagonal_weights) <- colnames(geno)

  for (method in c("VanRaden", "Yang", "Epistasis")) {
    actual <- PredictProR:::grm_calculation(
      geno_clean = geno,
      method = method,
      backend = "cpp"
    )
    expected <- reference_grm(geno, method)

    expect_equal(as.matrix(actual), expected, tolerance = 1e-10)
    expect_named_relationship_matrix(actual, rownames(geno))
  }

  for (method in c("Dominance", "Dominance_Vitezica", "Dominance_Su", "Dominance_Heterozygosity")) {
    actual <- PredictProR:::grm_calculation(
      geno_clean = geno,
      method = method,
      backend = "cpp"
    )
    expected <- reference_grm(geno, method)

    expect_equal(as.matrix(actual), expected, tolerance = 1e-10)
    expect_named_relationship_matrix(actual, rownames(geno))
  }

  actual_weighted <- PredictProR:::grm_calculation(
    geno_clean = geno,
    weight = weights,
    method = "Weighted_VanRaden",
    backend = "cpp"
  )
  expected_weighted <- reference_grm(geno, "Weighted_VanRaden", weight = weights)
  expect_equal(as.matrix(actual_weighted), expected_weighted, tolerance = 1e-10)
  expect_named_relationship_matrix(actual_weighted, rownames(geno))

  actual_diagonal <- PredictProR:::grm_calculation(
    geno_clean = geno,
    weight = diagonal_weights,
    method = "Weighted_VanRaden",
    backend = "cpp"
  )
  expected_diagonal <- reference_grm(
    geno,
    "Weighted_VanRaden",
    weight = diagonal_weights
  )
  expect_equal(as.matrix(actual_diagonal), expected_diagonal, tolerance = 1e-10)
  expect_named_relationship_matrix(actual_diagonal, rownames(geno))
})

test_that("GRM and omics kernels can be requested as named method banks", {
  geno <- matrix(
    c(
      0, 1, 2, 0,
      1, 2, 0, 1,
      2, 0, 1, 2,
      0, 1, 1, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:4), paste0("snp", 1:4))
  )
  x <- matrix(
    c(
      0.25, 1.50, -0.75,
      1.25, 0.50, 0.00,
      -0.50, 2.00, 1.25,
      0.75, -1.00, 0.50
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:4), paste0("m", 1:3))
  )

  grms <- PredictProR:::grm_calculation(
    geno_clean = geno,
    method = c("VanRaden", "Dominance_Vitezica", "Dominance_Su"),
    backend = "cpp"
  )
  expect_identical(names(grms), c("VanRaden", "Dominance_Vitezica", "Dominance_Su"))
  expect_equal(grms$VanRaden, reference_grm(geno, "VanRaden"), tolerance = 1e-10)
  expect_equal(grms$Dominance_Vitezica, reference_grm(geno, "Dominance_Vitezica"), tolerance = 1e-10)
  expect_equal(grms$Dominance_Su, reference_grm(geno, "Dominance_Su"), tolerance = 1e-10)

  kernels <- kernel_calculation(
    M_matrix_clean = x,
    scaling = FALSE,
    method = c("Linear_kernel", "Gaussian_kernel", "Poly2_kernel"),
    backend = "r",
    message = FALSE
  )
  expect_identical(names(kernels), c("Linear_kernel", "Gaussian_kernel", "Poly2_kernel"))
  expect_equal(kernels$Linear_kernel, reference_kernel(x, "Linear_kernel"), tolerance = 1e-10)
  expect_equal(kernels$Gaussian_kernel, reference_kernel(x, "Gaussian_kernel"), tolerance = 1e-10)
  expect_equal(kernels$Poly2_kernel, reference_kernel(x, "Poly2_kernel"), tolerance = 1e-10)
})

test_that("GP bridge serializes the full kernel bank for direct execution", {
  ids <- paste0("id", 1:3)
  K1 <- diag(3)
  K2 <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), nrow = 3)
  K3 <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.5, 0.2, 0.5, 1), nrow = 3)
  dimnames(K1) <- dimnames(K2) <- dimnames(K3) <- list(ids, ids)
  kernel <- PredictProR:::gp_public_kernel_bank_inputs(
    list(VanRaden = K1, Dominance_Su = K2, omic1_Gaussian_kernel = K3),
    pheno_ids = ids
  )
  spec <- PredictProR:::gp_bridge_prepare_mixed_model_spec(
    args = list(
      pheno_df = data.frame(GID = ids, Env = "E1", y = c(1, 2, NA_real_)),
      train_idx = c(0L, 1L),
      test_idx = 2L,
      gid_col = "GID",
      env_col = "Env",
      y_col = "y"
    ),
    kernel = kernel,
    project_root = getwd(),
    dir_path = tempfile("predictpror_gp_kernel_bank_")
  )
  payload <- jsonlite::fromJSON(spec$spec_json, simplifyVector = FALSE)
  expect_identical(names(payload$geno_kernels), c("VanRaden", "Dominance_Su", "omic1_Gaussian_kernel"))
  expect_true(all(file.exists(vapply(payload$geno_kernels, `[[`, character(1), "bin"))))
  expect_true(all(file.exists(vapply(payload$geno_kernels, `[[`, character(1), "meta"))))
})

test_that("backend argument is validated", {
  x <- matrix(1:9, nrow = 3)
  geno <- matrix(c(0, 1, 2, 1, 2, 0, 2, 0, 1), nrow = 3)

  expect_error(
    kernel_calculation(
      M_matrix_clean = x,
      method = "Linear_kernel",
      backend = "not_a_backend",
      message = FALSE
    ),
    "backend"
  )
  expect_error(
    PredictProR:::grm_calculation(
      geno_clean = geno,
      method = "VanRaden",
      backend = "not_a_backend"
    ),
    "backend"
  )
})

test_that("auto backend routes BLAS-backed and dominance methods deliberately", {
  withr::local_envvar(c(PREDICTPRO_KERNEL_BACKEND = NA, PREDICTPRO_GRM_BACKEND = NA))

  for (method in c(
    "Gaussian_kernel",
    "Linear_kernel",
    "Composite_kernel",
    "Poly2_kernel",
    "Poly3_kernel",
    "Poly4_kernel"
  )) {
    expect_identical(PredictProR:::gp_kernel_backend(NULL, method), "r")
  }

  for (method in c("VanRaden", "Weighted_VanRaden", "Epistasis")) {
    expect_identical(PredictProR:::gp_grm_backend(NULL, method), "r")
  }

  for (method in c("Dominance", "Dominance_Vitezica", "Dominance_Su", "Dominance_Heterozygosity")) {
    expect_identical(PredictProR:::gp_grm_backend(NULL, method), "cpp")
  }

  expect_identical(PredictProR:::gp_grm_backend(NULL, "Yang"), "r")
})

test_that("forced native kernel backend ignores parameters unused by the selected method", {
  x <- matrix(
    c(
      0.25, 1.50, -0.75,
      1.25, 0.50, 0.00,
      -0.50, 2.00, 1.25,
      0.75, -1.00, 0.50
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:4), paste0("m", 1:3))
  )

  expect_equal(
    kernel_calculation(
      x,
      scaling = FALSE,
      method = "Linear_kernel",
      theta = NA_real_,
      alpha = NA_real_,
      length_scale = NA_real_,
      backend = "cpp",
      message = FALSE
    ),
    kernel_calculation(
      x,
      scaling = FALSE,
      method = "Linear_kernel",
      theta = NA_real_,
      alpha = NA_real_,
      length_scale = NA_real_,
      backend = "r",
      message = FALSE
    )
  )

  expect_equal(
    kernel_calculation(
      x,
      scaling = FALSE,
      method = "Poly2_kernel",
      theta = NA_real_,
      alpha = NA_real_,
      length_scale = NA_real_,
      backend = "cpp",
      message = FALSE
    ),
    kernel_calculation(
      x,
      scaling = FALSE,
      method = "Poly2_kernel",
      theta = NA_real_,
      alpha = NA_real_,
      length_scale = NA_real_,
      backend = "r",
      message = FALSE
    )
  )
})

test_that("kernel calculation honors legacy scale alias used by pipeline callers", {
  x <- matrix(
    c(
      1, 10,
      2, 20,
      3, 30
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:3), c("m1", "m2"))
  )

  actual <- kernel_calculation(
      x,
      scale = FALSE,
      method = "Linear_kernel",
      backend = "r",
      message = FALSE
  )
  comparison <- all.equal(actual, reference_kernel(x, "Linear_kernel"))
  if (!isTRUE(comparison)) {
    stop("legacy scale alias was ignored", call. = FALSE)
  }
  expect_true(TRUE)
})

test_that("square GRM weight matrices must be diagonal", {
  geno <- matrix(
    c(
      0, 1,
      1, 2,
      2, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:3), c("s1", "s2"))
  )
  bad_weight <- matrix(
    c(1, 0.2, 0.2, 1),
    nrow = 2,
    dimnames = list(colnames(geno), colnames(geno))
  )

  expect_error(
    PredictProR:::grm_calculation(
      geno,
      weight = bad_weight,
      method = "Weighted_VanRaden",
      backend = "r"
    ),
    "diagonal"
  )
})

test_that("GRM validation rejects invalid genotype coding", {
  geno <- matrix(c(0, 1, 2, 1, 3, 0, 2, 0, 1), nrow = 3)

  expect_error(
    PredictProR:::grm_calculation(
      geno_clean = geno,
      method = "VanRaden",
      backend = "r",
      ploidy = 2L
    ),
    "dosage.*\\[0, 2\\]"
  )
})

test_that("kernel scaling reports constant columns", {
  x <- matrix(
    c(
      1, 2, 5,
      1, 3, 6,
      1, 4, 7
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_error(
    kernel_calculation(
      M_matrix_clean = x,
      scaling = TRUE,
      method = "Linear_kernel",
      backend = "r",
      message = FALSE
    ),
    "constant"
  )
})

test_that("GRM cleanup recomputes frequencies and resyncs weights", {
  geno <- matrix(
    c(
      0, 1, 2,
      1, 1, 0,
      2, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("id", 1:3), c("s1", "mono", "s2"))
  )
  weights <- matrix(
    c(1.25, 7.00, 0.80),
    ncol = 1,
    dimnames = list(colnames(geno), "weight")
  )

  actual <- PredictProR:::grm_calculation(
    geno_clean = geno,
    weight = weights,
    method = "Weighted_VanRaden",
    backend = "cpp"
  )

  cleaned_geno <- geno[, c("s1", "s2"), drop = FALSE]
  cleaned_weights <- weights[c("s1", "s2"), , drop = FALSE]
  expected <- reference_grm(
    cleaned_geno,
    "Weighted_VanRaden",
    weight = cleaned_weights
  )

  expect_equal(as.matrix(actual), expected, tolerance = 1e-10)
  expect_named_relationship_matrix(actual, rownames(geno))
})
