.joint_explicit_fixture <- function() {
  X <- matrix(rep(c(0,1,2,0,2,1,1,2,0,2,1,0),5),12,5,
    dimnames=list(paste0("g",1:12),paste0("m",1:5)))
  X[,2] <- rev(X[,2]); X[,4] <- c(0,2,1,2,0,1,2,1,0,1,0,2)
  list(X=X,ph=data.frame(ID=rownames(X),T1=c(1,4,2,5,3,6,2,1,5,4,7,3),
    T2=c(4,2,5,1,6,3,5,7,2,4,1,6)))
}

test_that("public joint positional selection resolves once against raw columns", {
  d <- .joint_explicit_fixture(); captured <- NULL
  local_mocked_bindings(gp_python_prepare_runtimes=function(...)NULL,
    gp_route_true_prediction_specialized_models=function(ctx) {
    captured <<- ctx; stop("captured explicit joint inputs")
  },.package="PredictProR")
  common <- list(pheno_data=d$ph,geno_data=d$X,response=c("T1","T2"),
    gen_name="ID",GS_model="RandomForest",multi_trait_ml=TRUE,
    response_family="gaussian",cross_validation=FALSE,para_tunning=FALSE,
    qc_filtering=FALSE,impute=FALSE,ld_prunning_qc=FALSE,ld_pruning=FALSE,
    message=FALSE)
  for(sel in list(c(2L,4L),c(FALSE,TRUE,FALSE,TRUE,FALSE),c("m2","m4"))) {
    expect_error(do.call(model_execute,c(common,list(feature_selected=sel))),
      "captured explicit joint inputs")
    expect_identical(colnames(captured$ml_dat_res$merged_data$merge_data),c("m2","m4"))
    out <- gp_feature_attach_joint_output(list(),c("T1","T2"),
      ml_dat_res=captured$ml_dat_res,feature_selected=captured$feature_selected,
      source_map=captured$feature_source_map)
    expect_identical(out$feature_selected,c("m2","m4"))
    expect_identical(out$feature_selection_metadata$source_block,rep("geno_data",2))
  }
})
test_that("joint source-specific positions retain actual source associations", {
  d <- .joint_explicit_fixture(); captured <- NULL
  O <- d$X[,1:3]; colnames(O) <- paste0("o",1:3)
  local_mocked_bindings(gp_python_prepare_runtimes=function(...)NULL,
    gp_route_true_prediction_specialized_models=function(ctx) {
    captured <<- ctx; stop("captured source inputs")
  },.package="PredictProR")
  expect_error(model_execute(pheno_data=d$ph,geno_data=d$X,omic1_data=O,
    response=c("T1","T2"),gen_name="ID",GS_model="RandomForest",multi_trait_ml=TRUE,
    feature_selected=list(geno_data=c(2L,4L),omic1_data=1L),
    response_family="gaussian",cross_validation=FALSE,para_tunning=FALSE,
    qc_filtering=FALSE,impute=FALSE,ld_prunning_qc=FALSE,ld_pruning=FALSE,message=FALSE),
    "captured source inputs")
  expect_identical(colnames(captured$ml_dat_res$merged_data$merge_data),c("m2","m4","o1"))
  expect_identical(unname(captured$feature_source_map[c("m2","m4","o1")]),
    c("geno_data","geno_data","omic1_data"))
})
test_that("kernel-only joint declarations remain explicitly unverified names", {
  ctx <- list(feature_selected=list(geno_data=c("m2","m4")),
    gen_name="ID",multi_trait_gp=TRUE,gmatrix=diag(3))
  out <- gp_feature_apply_raw_selection_to_inputs(ctx)
  expect_identical(out$gmatrix,ctx$gmatrix)
  expect_identical(out$feature_selected,ctx$feature_selected)
  expect_identical(unname(out$feature_selected_source_map),rep("geno_data",2))
  ctx$feature_selected <- c(2L,4L)
  expect_error(gp_feature_apply_raw_selection_to_inputs(ctx),"raw predictor columns")
  ctx$feature_selected <- list(geno_data=2L,omic1_data="o1")
  expect_error(gp_feature_apply_raw_selection_to_inputs(ctx),"raw predictor columns")
})
test_that("joint raw selection rejects ambiguous names and invalid positions", {
  d <- .joint_explicit_fixture()
  ctx <- list(geno_data=d$X,gen_name="ID",multi_trait_ml=TRUE,feature_selected=c(2.5,4))
  expect_error(gp_feature_apply_raw_selection_to_inputs(ctx),"whole")
  ctx$feature_selected <- c(2L,20L)
  expect_error(gp_feature_apply_raw_selection_to_inputs(ctx),"outside")
  ctx$feature_selected <- "m2"; ctx$omic1_data <- d$X[,2:3]
  expect_error(gp_feature_apply_raw_selection_to_inputs(ctx),"unique.*predictor names")
})


test_that("new example block preserves caller RNG state including errors and absent seed", {
  outer_kind <- RNGkind()
  outer_had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  outer_seed <- if(outer_had)get(".Random.seed",envir=.GlobalEnv) else NULL
  withr::defer({
    do.call(RNGkind,as.list(outer_kind))
    if(outer_had)assign(".Random.seed",outer_seed,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)
  })
  local_mocked_bindings(model_execute=function(...) {stats::runif(1);list()},
    .package="PredictProR")
  path <- system.file("examples","scenarios","07_feature_selection.R",package="PredictProR")
  if(!nzchar(path))path <- file.path("..","..","inst","examples","scenarios","07_feature_selection.R")
  expressions <- parse(path)
  block <- expressions[[length(expressions)-1L]]
  expect_match(paste(deparse(block),collapse=" "),"saved_kind")
  env <- new.env(parent=environment())
  env$labels <- character()
  env$sim_multi_trait <- function(...) {stats::runif(1);list(pheno=data.frame(),geno=matrix(0,1,1))}
  env$missing_software <- function(...)NULL
  env$needs_for_model <- function(...)character()
  env$log_skip <- function(...)stop("unexpected skip")
  env$log_result <- function(label,...)env$labels <- c(env$labels,label)
  set.seed(575L)
  expected <- get(".Random.seed",envir=.GlobalEnv)
  eval(block,env)
  expect_identical(get(".Random.seed",envir=.GlobalEnv),expected)
  expect_identical(RNGkind(),outer_kind)
  expect_length(env$labels,2L)
  rm(".Random.seed",envir=.GlobalEnv)
  eval(block,env)
  expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
  set.seed(575L)
  env$sim_multi_trait <- function(...) {stats::runif(1);stop("fixture error")}
  expect_error(eval(block,env),"fixture error")
  expect_identical(get(".Random.seed",envir=.GlobalEnv),expected)
})


test_that("joint raw normalization is idempotent and kernel routes retain resolved names", {
  d <- .joint_explicit_fixture()
  ctx <- list(geno_data=d$X,gen_name="ID",multi_trait_gp=TRUE,feature_selected=c(2L,4L))
  once <- gp_feature_apply_raw_selection_to_inputs(ctx)
  twice <- gp_feature_apply_raw_selection_to_inputs(once)
  expect_identical(twice$geno_data,once$geno_data)
  expect_identical(twice$feature_selected_resolved,c("m2","m4"))
  captured <- NULL
  local_mocked_bindings(gp_python_prepare_runtimes=function(...)NULL,
    gp_route_true_prediction_specialized_models=function(ctx) {
      captured <<- ctx;stop("captured joint kernel inputs")
    },.package="PredictProR")
  for(flag in c("multi_trait_bayes","multi_trait_gp","multi_trait_asreml","multi_trait_dl")) {
    model <- switch(flag,multi_trait_bayes="GBLUP_BRR",multi_trait_gp="Gaussian-Process-GBLUP",
      multi_trait_asreml="GBLUP",multi_trait_dl="mlp")
    args <- list(pheno_data=d$ph,geno_data=d$X,response=c("T1","T2"),gen_name="ID",
      GS_model=model,engine="asreml",response_family="gaussian",gmatrix_method="VanRaden",
      feature_selected=c(2L,4L),qc_filtering=FALSE,impute=FALSE,
      ld_prunning_qc=FALSE,ld_pruning=FALSE,message=FALSE)
    args[[flag]] <- TRUE
    expect_error(do.call(model_execute,args),"captured joint kernel inputs",info=flag)
    expect_identical(captured$feature_selected,c("m2","m4"),info=flag)
    expect_identical(colnames(captured$geno_omic_model_ready_list$geno_model_ready),
      c("m2","m4"),info=flag)
    expect_identical(unname(captured$feature_source_map[c("m2","m4")]),
      c("geno_data","geno_data"),info=flag)
  }
})


test_that("raw joint selection preserves ploidy when selected dosages do not reveal it", {
  d <- .joint_explicit_fixture()
  # Every observed dosage is at most two, but this is an annotated tetraploid.
  attr(d$X,"ploidy") <- 4L
  attr(d$X,"genotype_coding") <- "alternate_allele_dosage"
  attr(d$X,"genotype_coding_contract") <- "alt_dosage"
  selected <- gp_feature_select_raw_input(d$X,c("m2","m4"),"geno_data","ID")
  expect_identical(gp_resolve_matrix_ploidy(selected,"auto"),4L)
  expect_identical(attr(selected,"genotype_coding_contract"),"alt_dosage")
  expected_p <- colMeans(d$X[,c("m2","m4"),drop=FALSE])/4
  Z <- sweep(d$X[,c("m2","m4"),drop=FALSE],2,4*expected_p,"-")
  expected <- Z %*% t(Z)/sum(4*expected_p*(1-expected_p))
  bank <- gp_feature_build_kernel_bank(list(geno_data=selected),gmatrix_method="VanRaden",ploidy="auto")
  expect_lt(max(abs(bank$gmatrix_model_ready-expected)),1e-12)
})
