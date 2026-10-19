test_that("shared explicit sets preserve vectors and source lists", {
  response <- c("T1","T2")
  one <- c("m2","m1")
  expect_identical(gp_feature_joint_selected(one,response),one)
  sources <- list(geno_data=one,omic1_data="o1")
  expect_identical(gp_feature_joint_selected(sources,response),sources)
  same <- list(T1=c("m2","m1"),T2=c("m1","m2"))
  expect_identical(gp_feature_joint_selected(same,response),same$T1)
  expect_error(gp_feature_joint_selected(list(T1="m1",T2="m2"),response),"one predictor set")
  expect_error(gp_feature_joint_selected(list(T1="m1"),response),"T2")
})

test_that("per-trait source associations are compared without flattening", {
  sets <- list(T1=list(geno=c("m1","m2"),omic1="o1"),
               T2=list(geno_data=c("m2","m1"),omic1_data="o1"))
  expect_identical(gp_feature_joint_selected(sets,c("T1","T2")),sets$T1)
  sets$T2 <- list(geno_data="o1",omic1_data=c("m1","m2"))
  expect_error(gp_feature_joint_selected(sets,c("T1","T2")),"one predictor set")
})

test_that("trait-bearing selection tables respect flags and source aliases", {
  d <- data.frame(trait=c(rep("T1",3),rep("T2",3)),
                  predictor=c("m1","m2","unused","m2","unused","m1"),
                  source=c(rep("geno",3),rep("geno_data",3)),
                  selected=c(TRUE,TRUE,FALSE,TRUE,FALSE,TRUE))
  out <- gp_feature_joint_selected(d,c("T1","T2"))
  expect_false("trait"%in%names(out))
  expect_identical(gp_feature_selected_vector(out),c("m1","m2"))
  d$predictor[6L] <- "different"
  expect_error(gp_feature_joint_selected(d,c("T1","T2")),"one predictor set")
})

test_that("logical per-trait masks are compared by selected positions", {
  masks <- list(T1=c(TRUE,FALSE,FALSE),T2=c(FALSE,TRUE,FALSE))
  expect_error(gp_feature_joint_selected(masks,c("T1","T2")),"one predictor set")
  masks$T2 <- masks$T1
  expect_identical(gp_feature_joint_selected(masks,c("T1","T2")),masks$T1)
})

test_that("public joint-input refusal happens before preprocessing", {
  local_mocked_bindings(gp_shared_preprocess_cache_validate=function(...) {
    stop("preprocessing reached")
  },.package="PredictProR")
  ph <- data.frame(ID=c("a","b","c","d"),T1=c(1,3,2,4),T2=c(2,1,4,3))
  X <- matrix(c(0,0,2,2,0,2,0,2),4,2,dimnames=list(ph$ID,c("m1","m2")))
  common <- list(pheno_data=ph,geno_data=X,response=c("T1","T2"),gen_name="ID",
                 GS_model="RandomForest",response_family="gaussian",message=FALSE)
  expect_error(do.call(model_execute,c(common,list(multi_trait_ml=TRUE,
    feature_selected=list(T1="m1",T2="m2")))),"one predictor set")
  expect_error(do.call(model_execute,c(common,list(multi_trait=TRUE,
    feature_selected=list(T1="m1")))),"T2")
  expect_error(do.call(model_execute,c(common,list(hybrid_gp=TRUE,
    feature_selected=list(T1="m1",T2="m2")))),"one predictor set")
  meta <- data.frame(trait=c("T1","T2"),predictor=c("m1","m2"),
                     rank=1L,raw_score=c(2,1),scoring_scope="trait")
  expect_error(do.call(model_execute,c(common,list(multi_trait_ml=TRUE,
    feature_score_metadata=meta))),"joint = TRUE")
  expect_error(do.call(model_execute,c(common,list(
    feature_selected=list(T1="m1",T2="m2")))),"preprocessing reached")
  expect_error(do.call(model_execute,c(common,list(multi_trait_ml=TRUE,
    feature_selected=list(T1="m1",T2="m1")))),"preprocessing reached")
})

test_that("ambiguous source aliases and empty trait entries are refused", {
  expect_error(gp_feature_joint_selected(list(T1=list(geno="m1",geno_data="m2"),
    T2=list(geno_data="m1")),c("T1","T2")),"Repeated source")
  expect_error(gp_feature_joint_selected(list(T1=list(geno="m1",unknown="x"),
    T2=list(geno="m1",unknown="y")),c("T1","T2")),"Unrecognized source")
  expect_error(gp_feature_joint_selected(list(T1=NULL,T2=NULL),c("T1","T2")),"empty")
})


test_that("reserved joint label delimiter is refused for explicit selection before preprocessing", {
  local_mocked_bindings(gp_shared_preprocess_cache_validate=function(...)stop("preprocessing reached"),
    .package="PredictProR")
  ph <- data.frame(ID=letters[1:4],T1=c(1,3,2,4),T2=c(2,1,4,3))
  names(ph)[2L] <- "A+B"
  X <- matrix(c(0,0,2,2,0,2,0,2),4,2,dimnames=list(ph$ID,c("m1","m2")))
  expect_error(model_execute(pheno_data=ph,geno_data=X,response=c("A+B","T2"),
    gen_name="ID",multi_trait_ml=TRUE,GS_model="RandomForest",
    feature_selected="m1",response_family="gaussian",message=FALSE),"reserve.*[+]")
})
