.joint_tp_context <- function(flag, model, met=FALSE) {
  X <- matrix(c(0,2,0,2,0,2, 0,0,2,2,0,2, 2,0,2,0,2,0),6,3,
              dimnames=list(paste0("g",1:6),c("a","b","c")))
  ph <- data.frame(ID=rownames(X),T1=c(1,3,2,4,2,5),T2=c(4,1,5,2,3,6))
  if (met) ph <- rbind(transform(ph,Env="E1"),transform(ph,Env="E2",T1=T1+3,T2=T2-2))
  meta <- gp_feature_metadata_table("T1+T2",colnames(X),setNames(rep("geno_data",3),colnames(X)),
                                    c(1,3,2),"Ridge_Regression",1L,nrow(ph),scoring_scope="joint")
  ctx <- list(response=c("T1","T2"),response_family="gaussian",gen_name="ID",
              GS_model=model,cross_validation=FALSE,engine="asreml",fixed=~1,
              random=NULL,cova=NULL,var_cov_str="us",gp_fa_rank=1L,
              nIter=200L,burnIn=100L,thin=2L,random_state=1L,random_seed=1L,
              feature_score_metadata=meta,feature_k=2L,
              feature_source_map=setNames(rep("geno_data",3),colnames(X)),
              ml_dat_res=list(merged_data=list(merge_data=X)),
              geno_omic_model_ready_list=list(geno_model_ready=X),
              gmatrix_kernel_model_ready_list=list(gmatrix_model_ready=diag(99,6)),
              pheno_clean=list(pheno_clean_data=ph),pheno_data=ph,
              gmatrix_method="VanRaden",ploidy=2L,scaling=FALSE,centering=FALSE,
              heter_groups=if(met)"Env" else NULL,system_database=TRUE,message=FALSE)
  ctx[[flag]] <- TRUE
  if (flag=="multi_trait_gp") ctx$gp_public_multitrait_run <- TRUE
  ctx
}

test_that("all five TP routes expose the same selected columns and saved selection", {
  fit_input <- NULL; built_columns <- NULL; saved <- NULL; saved_scores <- NULL
  pred <- data.frame(ID=rep(paste0("g",1:6),2),Trait=rep(c("T1","T2"),each=6),
                     Predicted_value=1,Train_Test_Label="Train")
  output <- list(predicted_values=pred,model_parameters=data.frame(stat="fit",summary="mocked"))
  K <- diag(2,6); dimnames(K)<-list(paste0("g",1:6),paste0("g",1:6))
  local_mocked_bindings(
    gp_feature_build_kernel_bank=function(selected_sources,...) {
      built_columns <<- colnames(selected_sources$geno_data)
      list(gmatrix_model_ready=K)
    },
    bayes_multitrait_joint_fit=function(kernels,...) {
      fit_input <<- kernels[[1L]]
      list(bayes_result=output)
    },
    gp_multi_trait_model=function(gmatrix,...) {
      fit_input <<- gmatrix
      list(predictions=NULL,result=list(),info=list())
    },
    gp_multitrait_asreml_gaussian_model=function(gmatrix,...) {
      fit_input <<- gmatrix; output
    },
    gp_multitrait_asreml_mtmet_gaussian_model=function(gmatrix,...) {
      fit_input <<- gmatrix; output
    },
    gp_multitrait_ml_gaussian_model=function(geno_omic_object,...) {
      fit_input <<- geno_omic_object; output
    },
    gp_multitrait_dl_gaussian_model=function(geno_omic_object,...) {
      fit_input <<- geno_omic_object; output
    },
    gp_reconcile_multitrait_model_output=function(x,...)x,
    gp_standardize_public_model_result=function(x,...)x,
    gp_model_execute_multitrait_predicted_values=function(...)pred,
    gp_model_execute_multitrait_prediction_wide=function(...)NULL,
    gp_multitrait_ml_trait_counts=function(...)NULL,
    gp_model_execute_multitrait_summary=function(...)data.frame(),
    gp_multitrait_ml_summary_statistics=function(...)data.frame(),
    gp_multitrait_asreml_summary_statistics=function(...)data.frame(),
    gp_multitrait_dl_summary_statistics=function(...)data.frame(),
    gp_multitrait_ml_diagnostic_plot=function(...)NULL,
    gp_runtime_metadata=function(...)data.frame(),
    gp_runtime_profile_metadata_fields=function(...)list(),
    gp_detect_gp_python=function(...)NULL,
    gp_known_python=function(...)NULL,
    gp_known_preferred_python=function(...)NULL,
    results_handling=function(res_model_output,feature_selected,feature_score_metadata=NULL,...) {
      saved <<- feature_selected; saved_scores <<- feature_score_metadata; res_model_output
    },
    gp_index_specialized_model_execute_result=function(handled_result,...)handled_result,
    .package="PredictProR"
  )
  flags <- c("multi_trait_bayes","multi_trait_gp","multi_trait_asreml","multi_trait_ml","multi_trait_dl")
  models <- c("GBLUP_BRR","GP","GBLUP","RandomForest","mlp")
  for (i in seq_along(flags)) {
    ctx <- .joint_tp_context(flags[i],models[i])
    fit_input <- built_columns <- saved <- NULL
    out <- gp_route_true_prediction_specialized_models(ctx)
    expect_identical(out$feature_selected,c("b","c"),info=flags[i])
    expect_identical(out$feature_k,2L,info=flags[i])
    expect_identical(saved,c("b","c"),info=flags[i])
    expect_identical(saved_scores,ctx$feature_score_metadata,info=flags[i])
    expect_equal(sum(out$feature_selection_metadata$selected),2L,info=flags[i])
    expect_identical(unique(out$feature_selection_metadata$scoring_scope),"joint",info=flags[i])
    if (i<=3) {
      expect_identical(built_columns,c("b","c"),info=flags[i])
      expect_identical(fit_input,K,info=flags[i])
    } else expect_identical(colnames(fit_input),c("b","c"),info=flags[i])
    ctx$feature_score_metadata <- NULL; ctx$feature_k <- NULL
    control <- gp_route_true_prediction_specialized_models(ctx)
    expect_false(any(c("feature_selected","feature_k","feature_selection_metadata") %in%
                       names(control)),info=flags[i])
  }
  met <- .joint_tp_context("multi_trait_asreml","GBLUP",met=TRUE)
  met$feature_score_metadata <- feature_score_predictors(
    met$ml_dat_res$merged_data$merge_data,met$pheno_clean$pheno_clean_data,
    c("T1","T2"),"ID",joint=TRUE,heter_groups="Env",
    bayes_nIter=200L,bayes_burnIn=100L,bayes_thin=2L)
  expected <- met$feature_score_metadata$predictor[1:2]
  output_met <- gp_route_true_prediction_specialized_models(met)
  expect_identical(output_met$feature_selected,expected)
  expect_identical(built_columns,expected)
  expect_identical(unique(output_met$feature_selection_metadata$covariates_adjusted),"Env")
  expect_identical(unique(output_met$feature_selection_metadata$n_train),12L)
})

test_that("joint output fields survive public indexing without affecting single-trait fields", {
  meta <- gp_feature_metadata_table("T1+T2",c("a","b"),c(a="geno_data",b="geno_data"),
                                    c(2,1),"Ridge_Regression",1L,4L,scoring_scope="joint")
  meta$selected <- c(TRUE,FALSE)
  local_mocked_bindings(gp_gaussian_diagnostic_plots=function(...)list(),
                       .package="PredictProR")
  res <- list(predicted_values=data.frame(ID=c("g1","g1"),Trait=c("T1","T2"),
                                          Predicted_value=c(1,2)),feature_selection_metadata=meta,
              feature_selected="a",feature_k=1L)
  public <- gp_standardize_public_model_result(res,gen_name="ID",model="RandomForest",
                                  response_family="gaussian")
  expect_identical(public$feature_selected,"a")
  expect_identical(public$feature_k,1L)
  res$feature_selection_metadata$scoring_scope <- "trait"
  single <- gp_standardize_public_model_result(res,gen_name="ID",model="RandomForest",
                                  response_family="gaussian")
  expect_false(any(c("feature_selected","feature_k")%in%names(single)))
})

test_that("model_execute dispatches joint scoring for every dedicated joint flag", {
  captured <- NULL
  d <- .joint_tp_context("multi_trait_ml","RandomForest")
  X <- d$ml_dat_res$merged_data$merge_data
  local_mocked_bindings(
    gp_prepare_model_input_objects=function(ctx) {
      list(low_call_rate_inds_removed=NULL,pheno_clean=ctx$pheno_clean,
           geno_res=list(),omic1_res=list(),omic2_res=list(),omic3_res=list(),
           geno_omic_model_ready_list=list(geno_model_ready=X),
           gmatrix_kernel_model_ready_list=list(gmatrix_model_ready=diag(nrow(X))),
           test_set=NULL,ml_dat_res=list(merged_data=list(merge_data=X)))
    },
    feature_score_predictors=function(...) {
      captured <<- list(...)
      stop("scoring dispatched")
    },.package="PredictProR")
  flags <- c("multi_trait_bayes","multi_trait_gp","multi_trait_asreml","multi_trait_ml","multi_trait_dl")
  models <- c("GBLUP_BRR","GP","GBLUP","RandomForest","mlp")
  for (i in seq_along(flags)) {
    args <- list(pheno_data=d$pheno_clean$pheno_clean_data,geno_data=X,
                 response=c("T1","T2"),gen_name="ID",GS_model=models[i],
                 response_family="gaussian",feature_scoring=TRUE,feature_k=2L,
                 engine="asreml",gmatrix_method="VanRaden",message=FALSE)
    args[[flags[i]]] <- TRUE
    expect_error(do.call(model_execute,args),"scoring dispatched",info=flags[i])
    expect_true(captured$joint,info=flags[i])
    expect_null(captured$fixed,info=flags[i])
  }
  args <- list(pheno_data=d$pheno_clean$pheno_clean_data,geno_data=X,
               response=c("T1","T2"),gen_name="ID",GS_model="RandomForest",
               response_family="gaussian",feature_scoring=TRUE,feature_k=2L,message=FALSE)
  expect_error(do.call(model_execute,args),"scoring dispatched")
  expect_false(captured$joint)
})

test_that("explicit shared-set outputs carry provenance without invented scores", {
  result <- gp_feature_attach_joint_output(list(predicted_values=data.frame()),
    c("T1","T2"),ml_dat_res=list(feature_selected=c("m2","m1")),
    feature_selected=c("m1","m2"))
  expect_identical(result$feature_selected,c("m2","m1"))
  expect_identical(result$feature_k,2L)
  expect_true(all(is.na(result$feature_selection_metadata$raw_score)))
  expect_true(all(result$feature_selection_metadata$selection_mode=="explicit"))
})


test_that("public joint GP dispatch remains joint with explicit and scored selection", {
  d <- .joint_tp_context("multi_trait_gp","GP")
  X <- d$ml_dat_res$merged_data$merge_data
  captured <- NULL
  local_mocked_bindings(
    gp_prepare_model_input_objects=function(ctx)list(low_call_rate_inds_removed=NULL,
      pheno_clean=ctx$pheno_clean,geno_res=list(),omic1_res=list(),omic2_res=list(),omic3_res=list(),
      geno_omic_model_ready_list=list(geno_model_ready=X),
      gmatrix_kernel_model_ready_list=list(gmatrix_model_ready=diag(nrow(X))),
      test_set=NULL,ml_dat_res=list(merged_data=list(merge_data=X))),
    feature_score_predictors=function(...)d$feature_score_metadata,
    gp_route_true_prediction_specialized_models=function(ctx) {
      captured <<- ctx; stop("captured public joint route")
    },.package="PredictProR")
  common <- list(pheno_data=d$pheno_clean$pheno_clean_data,geno_data=X,
    response=c("T1","T2"),gen_name="ID",GS_model="Gaussian-Process-GBLUP",
    multi_trait_gp=TRUE,response_family="gaussian",gmatrix_method="VanRaden",message=FALSE)
  variants <- list(list(),list(feature_selected=c("b","c")),
    list(feature_score_metadata=d$feature_score_metadata,feature_k=2L),
    list(feature_scoring=TRUE,feature_k=2L))
  for(extra in variants) {
    expect_error(do.call(model_execute,c(common,extra)),"captured public joint route")
    expect_true(captured$gp_public_multitrait_run,info=paste(names(extra),collapse=","))
  }
})
