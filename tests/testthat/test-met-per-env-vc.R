# Rigor: MET variance-component output should report per-environment genetic
# variance + heritability across engines (like ASReml corgh), so users can
# compare. gp_variance_components_add_per_env() appends per-env
# genetic_variance_<env> / residual_variance_<env> / heritability_<env> rows to an
# engine's variance-components table, preserving existing rows and column layout.

resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  ns <- tryCatch(asNamespace("PredictProR"), error = function(e) NULL)
  if (!is.null(ns) && exists(name, envir = ns, mode = "function")) return(get(name, envir = ns, mode = "function"))
  NULL
}

test_that("gp_variance_components_add_per_env appends per-env Vg/Ve/h2 rows, preserving the table", {
  fn <- resolve_fn("gp_variance_components_add_per_env")
  skip_if(is.null(fn), "gp_variance_components_add_per_env not available")

  vc <- data.frame(Components = c(0.5, 0.3, 0.625), Standard_error = NA_real_,
                   row.names = c("genetic_variance", "residual_variance_mean", "heritability"),
                   stringsAsFactors = FALSE)
  eg <- c(E1 = 2, E2 = 3); er <- c(E1 = 1, E2 = 1)
  out <- fn(vc, env_genetic = eg, env_residual = er)

  # original rows preserved
  expect_true(all(rownames(vc) %in% rownames(out)))
  # per-env rows present
  expect_true(all(c("genetic_variance_E1", "genetic_variance_E2",
                    "residual_variance_E1", "residual_variance_E2",
                    "heritability_E1", "heritability_E2") %in% rownames(out)))
  # values correct
  expect_equal(out["genetic_variance_E1", "Components"], 2)
  expect_equal(out["residual_variance_E2", "Components"], 1)
  expect_equal(out["heritability_E1", "Components"], 2 / 3, tolerance = 1e-8)   # 2/(2+1)
  expect_equal(out["heritability_E2", "Components"], 3 / 4, tolerance = 1e-8)   # 3/(3+1)
  # column layout preserved
  expect_identical(names(out), names(vc))
})

test_that("gp_met_varcomp_augment_per_env reconstructs per-env Vg = sum(loading^2)+specific and h2", {
  fn <- resolve_fn("gp_met_varcomp_augment_per_env")
  skip_if(is.null(fn), "gp_met_varcomp_augment_per_env not available")
  # Mirrors a real GP_FA MET varcomp table (loadings, specific vars, residuals).
  rn <- c("fa(Env,1):vm(GID,G)!B2IR!var", "fa(Env,1):vm(GID,G)!B2IR!fa1", "Env_B2IR!R",
          "fa(Env,1):vm(GID,G)!F5I!var",  "fa(Env,1):vm(GID,G)!F5I!fa1",  "Env_F5I!R")
  vc <- data.frame(Component = rn,
                   Components = c(0.05454597, 0.06949943, 0.04212732,
                                  NA,         -0.35970896, 0.03405660),
                   Standard_error = NA_real_, row.names = rn, stringsAsFactors = FALSE)
  out <- fn(vc)
  expect_true(all(c("genetic_variance_B2IR", "heritability_B2IR",
                    "genetic_variance_F5I", "heritability_F5I") %in% rownames(out)))
  vg_b2ir <- 0.06949943^2 + 0.05454597
  expect_equal(out["genetic_variance_B2IR", "Components"], vg_b2ir, tolerance = 1e-6)
  expect_equal(out["heritability_B2IR", "Components"],
               vg_b2ir / (vg_b2ir + 0.04212732), tolerance = 1e-6)
  vg_f5i <- (-0.35970896)^2 + 0                 # specific var NA -> 0 (boundary)
  expect_equal(out["genetic_variance_F5I", "Components"], vg_f5i, tolerance = 1e-6)
  expect_equal(out["heritability_F5I", "Components"],
               vg_f5i / (vg_f5i + 0.03405660), tolerance = 1e-6)
  # original FA rows preserved
  expect_true(all(rn %in% rownames(out)))
})

test_that("gp_met_varcomp_augment_per_env returns non-FA tables unchanged", {
  fn <- resolve_fn("gp_met_varcomp_augment_per_env")
  skip_if(is.null(fn), "gp_met_varcomp_augment_per_env not available")
  vc <- data.frame(Components = c(0.5, 0.3, 0.6), Standard_error = NA_real_,
                   row.names = c("genetic_variance", "residual_variance", "heritability"),
                   stringsAsFactors = FALSE)
  expect_identical(fn(vc), vc)
})

test_that("gp_met_varcomp_augment_per_env reconstructs per-env h2 for a compound-symmetry table", {
  fn <- resolve_fn("gp_met_varcomp_augment_per_env")
  skip_if(is.null(fn), "gp_met_varcomp_augment_per_env not available")
  # GP / KRR / LowRankGP MET layout: single main genetic var + (boundary) GxE +
  # per-env residuals. Per-env Vg is homogeneous = sum of non-residual !var rows.
  rn <- c("vm(GID,gmatrix)!var", "Env!var", "Env_B2IR!R", "Env_F5I!R")
  vc <- data.frame(Components = c(0.0133, NA, 0.0955, 0.1634), Standard_error = NA_real_,
                   row.names = rn, stringsAsFactors = FALSE)
  out <- fn(vc)
  expect_true(all(c("genetic_variance_B2IR", "heritability_B2IR",
                    "residual_variance_F5I", "heritability_F5I") %in% rownames(out)))
  # homogeneous genetic variance (Env!var NA -> 0)
  expect_equal(out["genetic_variance_B2IR", "Components"], 0.0133, tolerance = 1e-8)
  expect_equal(out["genetic_variance_F5I", "Components"], 0.0133, tolerance = 1e-8)
  expect_equal(out["heritability_B2IR", "Components"], 0.0133 / (0.0133 + 0.0955), tolerance = 1e-6)
  expect_equal(out["heritability_F5I", "Components"], 0.0133 / (0.0133 + 0.1634), tolerance = 1e-6)
})

test_that("gp_variance_components_add_per_env tolerates a 'Component' label column and NA genetics", {
  fn <- resolve_fn("gp_variance_components_add_per_env")
  skip_if(is.null(fn), "gp_variance_components_add_per_env not available")
  vc <- data.frame(Component = "x", Components = 1, Standard_error = NA_real_,
                   row.names = "x", stringsAsFactors = FALSE)
  out <- fn(vc, env_genetic = c(A = NA_real_, B = 0.4), env_residual = c(A = 0.1, B = 0.1))
  expect_true("Component" %in% names(out))
  expect_identical(out["genetic_variance_A", "Component"], "genetic_variance_A")  # label filled
  expect_true(is.na(out["heritability_A", "Components"]))                          # NA genetic -> NA h2
  expect_equal(out["heritability_B", "Components"], 0.8, tolerance = 1e-8)
})
