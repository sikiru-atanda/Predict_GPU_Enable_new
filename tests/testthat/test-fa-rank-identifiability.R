# FA-GBLUP across environments: gp_fa_rank above the identifiable maximum used
# to fit silently without converging (variance components not promoted).

test_that("largest identifiable FA rank follows the Ledermann bound", {
  k_max <- vapply(c(2, 3, 4, 5, 6, 8, 10), PredictProR:::gp_fa_max_identifiable_rank, integer(1))
  expect_identical(k_max, c(0L, 1L, 1L, 2L, 3L, 4L, 6L))
})

test_that("FA-GBLUP stops when gp_fa_rank exceeds what the environments can identify", {
  ph <- data.frame(GID = rep(paste0("g", 1:6), 3), Env = rep(c("E1", "E2", "E3"), each = 6),
                   Yield = stats::rnorm(18))
  check <- function(k, n_env = 3L, model = "FA-GBLUP") {
    p <- ph[ph$Env %in% paste0("E", seq_len(n_env)), ]
    PredictProR:::gp_validate_fa_rank_for_environments(model, k, p, "Env")
  }
  expect_true(check(1L))
  expect_error(check(1L, n_env=2L), "4 parameters.*3 covariance entries")
  expect_error(check(2L), "FA-GBLUP with 3 environments supports gp_fa_rank up to 1")
  expect_true(check(5L, model = "Kernel-GBLUP"))          # other models ignore gp_fa_rank
  ph5 <- data.frame(GID = rep(paste0("g", 1:6), 5), Env = rep(paste0("E", 1:5), each = 6), Yield = 1)
  expect_true(PredictProR:::gp_validate_fa_rank_for_environments("FA-GBLUP", 2L, ph5, "Env"))
  expect_error(PredictProR:::gp_validate_fa_rank_for_environments("FA-GBLUP", 3L, ph5, "Env"),
               "supports gp_fa_rank up to 2")
  expect_error(check(0L), "positive integer")

  G <- matrix(stats::rbinom(6 * 30, 2, 0.5), 6, 30, dimnames = list(paste0("g", 1:6), paste0("m", 1:30)))
  expect_error(
    model_execute(pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield",
                  heter_groups = "Env", GS_model = "FA-GBLUP", gp_fa_rank = 2L, message = FALSE),
    "FA-GBLUP with 3 environments supports gp_fa_rank up to 1"
  )
})

test_that("trait and GxE FA ranks stop before Python on invalid dimensions", {
  local_mocked_bindings(gp_detect_gp_python=function(...) stop("PYTHON MUST NOT START"),
                       .package="PredictProR")
  for (nk in list(c(2L,1L),c(3L,2L),c(4L,2L))) {
    nt <- nk[1]; rank <- nk[2]
    ph <- data.frame(GID=paste0("g",1:6))
    traits <- paste0("T",seq_len(nt))
    for (t in traits) ph[[t]] <- seq_len(6)
    K <- diag(6); dimnames(K) <- list(ph$GID,ph$GID)
    expect_error(gp_multi_trait_model(ph,K,traits,"GID",trait_structure="fa",
      trait_fa_rank=rank),"parameters.*covariance entries")
    expect_error(gp_multi_trait_model(ph,K,traits,"GID",residual_structure="fa",
      residual_fa_rank=rank),"parameters.*covariance entries")
    met <- rbind(transform(ph,Env="E1"),transform(ph,Env="E2"))
    expect_error(gp_multi_trait_met_model(pheno_data=met,gmatrix=K,response=traits,gen_name="GID",
      heter_groups="Env",varcomp_mode="mom",trait_structure="fa",trait_fa_rank=rank),
      "parameters.*covariance entries")
    expect_error(model_execute(pheno_data=ph,gmatrix=K,response=traits,gen_name="GID",
      multi_trait_gp=TRUE,GS_model="FA-GBLUP",gp_fa_rank=rank,message=FALSE),
      "parameters.*covariance entries")
    expect_error(gp_multi_trait_met_model(pheno_data=met,gmatrix=K,response=traits,gen_name="GID",
      heter_groups="Env",varcomp_mode="mom",gxe_trait_structure="fa",gxe_trait_fa_rank=rank),
      "parameters.*covariance entries")
  }
})

test_that("FA rank counts accept five traits at rank two and label saturation honestly", {
  expect_true(PredictProR:::gp_validate_fa_rank(2L,5L,"traits","trait_fa_rank"))
  expect_true(PredictProR:::gp_validate_fa_rank(3L,6L,"traits","trait_fa_rank"))
  expect_error(PredictProR:::gp_validate_fa_rank(1.5,5L,"traits","trait_fa_rank"),"positive integer")
})

test_that("direct single-trait FA rejects two environments before Python", {
  local_mocked_bindings(gp_detect_gp_python=function(...) stop("PYTHON MUST NOT START"),
                       .package="PredictProR")
  ph <- expand.grid(GID=paste0("g",1:6),Env=c("E1","E2")); ph$Y <- seq_len(nrow(ph))
  K <- diag(6); dimnames(K) <- list(paste0("g",1:6),paste0("g",1:6))
  expect_error(gp_single_trait_model(ph,K,"Y","GID",heter_groups="Env",model_name="GP_FA"),
               "4 parameters.*3 covariance entries")
})
