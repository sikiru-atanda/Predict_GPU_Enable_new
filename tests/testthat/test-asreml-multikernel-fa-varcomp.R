# Regression: multi-kernel factor-analytic ASReml GBLUP variance extraction.
#
# Background (user probe, plain `GBLUP` + var_cov_str = "fa1" + gmatrix +
# omics kernel): random_terms_fit_new()::build_structured_vm_term gives the
# FIRST kernel the fa() structure and every LATER kernel an idv() fallback.
# The idv kernel has a SINGLE homogeneous variance component
# (`YYY:vm(GID, omics_inv)!YYY`) and NO per-env `!var`/`!fa` rows, so the FA
# branch of asreml_herit_varCov_new() built an empty `VarG` and then errored
# on `names(VarG) <- heter_grp` with
#   'names' attribute [5] must be the same length as the vector [0]
# That error was swallowed upstream and the whole variance_components table
# collapsed to a scalar genetic_variance with NA residual_variance and NA
# heritability. This test pins the fix: the idv-fallback kernel is treated as
# a constant per-env genetic variance, extraction succeeds, and the pooled
# summary carries finite residual_variance and heritability.

# A fake `summary.asreml` object: asreml_varcomp_table() returns `$varcomp`
# verbatim for summary.asreml inputs (no ASReml needed), and the function
# reads the env factor from `$mf`.
make_fake_fa_model <- function() {
  envs <- c("B2IR", "B5I", "BLHT", "F5I", "FDRIP")  # alphabetical = levels order
  spec  <- c(0.005688, 0.047154, 0.377416, 0.283259, 0.097157)  # gmatrix !var
  load1 <- c(-0.300509, -0.074809, -0.129063, -0.039859, -0.427131)  # gmatrix !fa1
  resid <- c(0.019479, 0.007620, 0.022957, 0.045708, 0.011780)  # YYY_<env>!R
  omic_idv <- 0.007140  # single idv variance for the omics kernel

  rn <- c(
    paste0("fa(YYY, 1):vm(GID, gmatrix_inv)!", envs, "!var"),
    paste0("fa(YYY, 1):vm(GID, gmatrix_inv)!", envs, "!fa1"),
    "YYY:vm(GID, omics_inv)!YYY",
    paste0("YYY_", envs, "!R")
  )
  component <- c(spec, load1, omic_idv, resid)
  varcomp <- data.frame(
    component  = component,
    std.error  = abs(component) * 0.2 + 0.01,
    z.ratio    = component / (abs(component) * 0.2 + 0.01),
    bound      = rep("P", length(component)),
    `%ch`      = rep(0.1, length(component)),
    row.names  = rn,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  mf <- data.frame(YYY = factor(rep(envs, length.out = 50L), levels = envs))
  structure(list(varcomp = varcomp, mf = mf), class = "summary.asreml")
}

test_that("multi-kernel FA extraction does not crash on the idv-fallback kernel", {
  fake <- make_fake_fa_model()
  expect_error(
    asreml_herit_varCov_new(
      model = fake,
      heter_groups = "YYY",
      var_cov_str = "fa1",
      heter_resid = TRUE,
      names_in_inv_list = c("gmatrix_inv", "omics_inv"),
      inter_gen_pos = 2L,
      gen_pos = 1L
    ),
    NA  # expect NO error (previously: 'names' attribute [5] ... vector [0])
  )
})

test_that("multi-kernel FA GBLUP reports finite residual_variance and heritability", {
  fake <- make_fake_fa_model()
  res <- asreml_herit_varCov_new(
    model = fake,
    heter_groups = "YYY",
    var_cov_str = "fa1",
    heter_resid = TRUE,
    names_in_inv_list = c("gmatrix_inv", "omics_inv"),
    inter_gen_pos = 2L,
    gen_pos = 1L
  )

  # The idv-fallback (omics) kernel contributes a constant per-env variance.
  expect_equal(nrow(res$varG_per_omics), 2L)
  expect_true(all(is.finite(res$varG_per_omics["omics_inv", ])))
  expect_equal(unname(res$varG_per_omics["omics_inv", ]),
               rep(0.007140, 5L), tolerance = 1e-6)

  # Per-env genetic variance = gmatrix FA diagonal + omics idv constant.
  envs <- c("B2IR", "B5I", "BLHT", "F5I", "FDRIP")
  expect_true(all(is.finite(res$Total_genetic_var)))
  expect_true(all(res$Total_genetic_var > 0))
  expect_true(all(is.finite(res$Residual_Var)))
  expect_true(all(is.finite(res$Heritability)))

  # Downstream formatter: pooled summary must carry finite residual + h2
  # (the symptom was NA residual_variance and NA heritability).
  vc <- asreml_variance_components(res_var_cov_h_ve = res)
  if (is.list(vc) && !is.data.frame(vc)) vc <- vc[["Variance_components"]]
  expect_identical(rownames(vc)[1:3],
                   c("mean_genetic_variance_across_environments",
                     "mean_residual_variance_across_environments",
                     "heritability_from_mean_variances"))
  expect_true(is.finite(vc["mean_genetic_variance_across_environments", "Components"]))
  expect_true(is.finite(vc["mean_residual_variance_across_environments", "Components"]))
  expect_true(is.finite(vc["heritability_from_mean_variances", "Components"]))
  expect_gt(vc["heritability_from_mean_variances", "Components"], 0)
  expect_lt(vc["heritability_from_mean_variances", "Components"], 1)
})

test_that("malformed FA output is not silently converted to a zero idv kernel", {
  fake <- make_fake_fa_model()
  keep <- !grepl("omics_inv", rownames(fake$varcomp), fixed = TRUE)
  fake$varcomp <- fake$varcomp[keep, , drop = FALSE]

  expect_error(
    asreml_herit_varCov_new(
      model = fake,
      heter_groups = "YYY",
      var_cov_str = "fa1",
      heter_resid = TRUE,
      names_in_inv_list = c("gmatrix_inv", "omics_inv"),
      inter_gen_pos = 2L,
      gen_pos = 1L
    ),
    "Could not reconstruct fa1 covariance for kernel 'omics_inv'",
    fixed = TRUE
  )
})
