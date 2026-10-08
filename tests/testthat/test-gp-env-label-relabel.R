# Regression: GP variance-component labels must use the user-defined env
# variable name (heter_groups), not the internal hardcoded "Env".
#
# Background (user probe): the GP backend standardises the env factor to "Env",
# so varcomp rows came back labelled `Env!var`, `Env_<level>!R` and
# `fa(Env,1):vm(GID,G)!...` even when the user supplied heter_groups = "YYY".
# gp_vc_relabel_env_to_heter_groups() rewrites those tokens to the user's name.

test_that("gp_vc_relabel_env_to_heter_groups rewrites Env tokens to heter_groups", {
  vc <- data.frame(
    Component = c("vm(GID,gmatrix)!var", "Env!var", "Env_B2IR!R", "Env_F5I!R",
                  "fa(Env,1):vm(GID,G)!B2IR!var", "fa(Env,1):vm(GID,G)!B2IR!fa1",
                  "idv(Env):vm(GID,omic1)!var", "genetic_variance_B2IR"),
    Components = c(0.5, NA, 0.3, 0.4, 0.2, 0.1, 0.05, 0.2),
    Standard_error = NA_real_,
    stringsAsFactors = FALSE
  )
  rownames(vc) <- vc$Component

  out <- gp_vc_relabel_env_to_heter_groups(vc, "YYY")

  expect_true("YYY!var" %in% out$Component)
  expect_true("YYY_B2IR!R" %in% out$Component)
  expect_true("YYY_F5I!R" %in% out$Component)
  expect_true("fa(YYY,1):vm(GID,G)!B2IR!var" %in% out$Component)
  expect_true("fa(YYY,1):vm(GID,G)!B2IR!fa1" %in% out$Component)
  expect_true("idv(YYY):vm(GID,omic1)!var" %in% out$Component)

  # No "Env" token may remain in any structural label.
  expect_false(any(grepl("(^|[(_:])Env([!_,):])", out$Component)))
  expect_identical(rownames(out), out$Component)

  # Untouched: kernel rows and per-env LEVEL summary rows.
  expect_true("vm(GID,gmatrix)!var" %in% out$Component)
  expect_true("genetic_variance_B2IR" %in% out$Component)
})

test_that("gp_vc_relabel_env_to_heter_groups is a no-op for default 'Env' / NULL", {
  vc <- data.frame(
    Component = c("Env!var", "Env_B2IR!R"),
    Components = c(NA, 0.3),
    Standard_error = NA_real_,
    stringsAsFactors = FALSE
  )
  expect_identical(gp_vc_relabel_env_to_heter_groups(vc, "Env"), vc)
  expect_identical(gp_vc_relabel_env_to_heter_groups(vc, NULL), vc)
  expect_identical(gp_vc_relabel_env_to_heter_groups(vc, ""), vc)
})
