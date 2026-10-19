# The test suite never installs Python packages: R CMD check must not reach
# the network or write into the user's Python environments. Tests that need a
# Python backend skip when its packages are missing (helper-python.R); the
# installer itself is tested with mocks and an offline wheel
# (test-python-dependencies.R). Set PREDICTPRO_AUTO_INSTALL=true before the
# run to let PredictProR install what the Python tests need.
if (!nzchar(Sys.getenv("PREDICTPRO_AUTO_INSTALL"))) {
  withr::local_envvar(c(PREDICTPRO_AUTO_INSTALL = "false"), .local_envir = testthat::teardown_env())
}
# Nor does it keep interpreter checks in the user's cache folder; the disk
# cache is tested on a private file (test-python-probe-disk-cache.R).
withr::local_envvar(c(PREDICTPRO_PROBE_CACHE = "false"), .local_envir = testthat::teardown_env())
