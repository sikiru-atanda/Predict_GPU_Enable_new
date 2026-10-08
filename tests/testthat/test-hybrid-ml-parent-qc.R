mk_parent_qc_data <- function() {
  set.seed(11)
  females <- sprintf("F%d", 1:5)
  males <- sprintf("M%d", 1:4)
  markers <- sprintf("m%02d", 1:20)
  fg <- matrix(sample(c(0, 2), 5 * 20, replace = TRUE), 5, 20, dimnames = list(females, markers))
  mg <- matrix(sample(c(0, 2), 4 * 20, replace = TRUE), 4, 20, dimnames = list(males, markers))
  # m01 is heterozygous in 3 of 5 females (residual heterozygosity), m02 in
  # 2 of 4 males; every other marker is fully inbred.
  fg[1:3, "m01"] <- 1
  mg[1:2, "m02"] <- 1
  ph <- expand.grid(Female = females, Male = males, stringsAsFactors = FALSE)
  ph$HybridID <- paste(ph$Female, ph$Male, sep = "_")
  ph$Yield <- rnorm(nrow(ph))
  list(pheno = ph, fg = fg, mg = mg)
}

test_that("parent het QC drops failing markers from hybrid-level features", {
  d <- mk_parent_qc_data()
  hyb <- (d$fg[d$pheno$Female, ] + d$mg[d$pheno$Male, ]) / 2
  rownames(hyb) <- d$pheno$HybridID
  hyb <- cbind(hyb, extra = 1)

  out <- PredictProR:::gp_hybrid_ml_parent_geno_prepare(
    pheno_data = d$pheno, geno_data = hyb,
    female_geno_data = d$fg, male_geno_data = d$mg,
    gen_name = "HybridID", female_parent = "Female", male_parent = "Male",
    het_threshold = 0.1, ploidy = 2L, message = FALSE
  )
  expect_setequal(colnames(out$geno_data), c(sprintf("m%02d", 3:20), "extra"))
  expect_identical(rownames(out$geno_data), rownames(hyb))
  s <- setNames(out$summary$value, out$summary$item)
  expect_equal(s[["markers_removed_parent_heterozygosity"]], "2")
  expect_equal(s[["hybrid_markers_not_in_parent_genotypes"]], "1")
  expect_equal(s[["hybrid_feature_source"]], "hybrid_geno_data")
})

test_that("parent genotypes alone build expected F1 hybrid features", {
  d <- mk_parent_qc_data()
  out <- PredictProR:::gp_hybrid_ml_parent_geno_prepare(
    pheno_data = d$pheno, geno_data = NULL,
    female_geno_data = d$fg, male_geno_data = d$mg,
    gen_name = "HybridID", female_parent = "Female", male_parent = "Male",
    het_threshold = 0.1, ploidy = 2L, message = FALSE
  )
  x <- out$geno_data
  expect_identical(rownames(x), d$pheno$HybridID)
  expect_identical(colnames(x), sprintf("m%02d", 3:20))
  expect_equal(
    unname(x["F4_M3", ]),
    unname((d$fg["F4", colnames(x)] + d$mg["M3", colnames(x)]) / 2)
  )
  expect_true(all(x %in% c(0, 1, 2)))
})

test_that("het_threshold = NULL keeps every parent marker", {
  d <- mk_parent_qc_data()
  out <- PredictProR:::gp_hybrid_ml_parent_geno_prepare(
    pheno_data = d$pheno, female_geno_data = d$fg, male_geno_data = d$mg,
    gen_name = "HybridID", female_parent = "Female", male_parent = "Male",
    het_threshold = NULL, ploidy = 2L, message = FALSE
  )
  expect_equal(ncol(out$geno_data), 20L)
})

test_that("missing parent genotypes and conflicting parent pairs are reported", {
  d <- mk_parent_qc_data()
  expect_error(
    PredictProR:::gp_hybrid_ml_parent_geno_prepare(
      pheno_data = d$pheno, female_geno_data = d$fg[-1, ], male_geno_data = d$mg,
      gen_name = "HybridID", female_parent = "Female", male_parent = "Male",
      het_threshold = 0.1, ploidy = 2L, message = FALSE
    ),
    "female: F1"
  )
  bad <- d$pheno
  bad$HybridID[2] <- bad$HybridID[1]
  expect_error(
    PredictProR:::gp_hybrid_ml_parent_geno_prepare(
      pheno_data = bad, female_geno_data = d$fg, male_geno_data = d$mg,
      gen_name = "HybridID", female_parent = "Female", male_parent = "Male",
      het_threshold = 0.1, ploidy = 2L, message = FALSE
    ),
    "one female/male parent pair"
  )
})
