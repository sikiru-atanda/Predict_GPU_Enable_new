.joint_cv_data <- function() {
  withr::local_seed(6201)
  X <- matrix(rnorm(12*5),12,dimnames=list(paste0("g",1:12),paste0("m",1:5)))
  ph <- data.frame(ID=rownames(X),T1=X[,1]+rnorm(12),T2=X[,2]+rnorm(12))
  rownames(ph) <- NULL
  meta <- gp_feature_metadata_table("T1+T2",colnames(X),setNames(rep("geno_data",5),colnames(X)),5:1,"Ridge_Regression",1L,12L,scoring_scope="joint")
  K <- diag(12); dimnames(K)<-list(rownames(X),rownames(X))
  list(X=X,ph=ph,meta=meta,K=K)
}

test_that("specialized CV accepts all and rejects invalid cardinalities without truncating", {
  m <- .joint_cv_data()$meta
  expect_identical(gp_feature_specialized_cv_config(m,"all")$k,"all")
  expect_identical(gp_feature_specialized_cv_config(m,8)$k,8L)
  for (k in list(0,-1,1.5,NA_real_,Inf,c(1,2),c(2,NA)))
    expect_error(gp_feature_specialized_cv_config(m,k),"feature_k")
  expect_error(gp_feature_specialized_cv_config(m,2,feature_scoring_cv="both"),"separately")
})

test_that("joint CV callers mask the intended cells before either selection mode", {
  d <- .joint_cv_data(); captured <- NULL
  local_mocked_bindings(
    gp_tuning_folds=function(...)list(1:4,5:8,9:12),
    gp_asreml_require_engine=function(...)invisible(NULL),
    gp_feature_joint_view=function(...) {captured <<- list(...); stop("joint view captured")},
    .package="PredictProR")
  for (mode in c("fixed","fold_internal")) {
    base <- list(pheno_object=d$ph,response=c("T1","T2"),gen_name="ID",
      nfolds=3L,feature_score_metadata=d$meta,feature_k=2L,feature_scoring_cv=mode)
    for (route in c("ml","asreml","dl")) {
      captured <- NULL
      args <- switch(route,ml=c(base,list(model_type="Ridge_Regression",geno_omic_object=d$X)),
        dl=c(base,list(model_type="mlp",geno_omic_object=d$X)),
        asreml=c(base,list(gmatrix=d$K,feature_predictor_data=d$X)))
      fn <- switch(route,ml=gp_multitrait_ml_gaussian_cv,
                   dl=gp_multitrait_dl_gaussian_cv,asreml=gp_multitrait_asreml_gaussian_cv)
      expect_error(do.call(fn,args),"joint view captured",info=paste(route,mode))
      expect_identical(captured$response,c("T1","T2"))
      expect_identical(captured$selection_mode,mode)
      expect_false("test_rows"%in%names(captured))
      visible <- captured$pheno_visible
      expect_true(all(is.na(visible$T1[1:4])))
      if(route=="ml") {
        expect_identical(visible$T2,d$ph$T2)
        expect_identical(captured$held_out_trait,"T1")
      } else expect_true(all(is.na(visible$T2[1:4])))
      expect_equal(visible[5:12,c("T1","T2")],d$ph[5:12,c("T1","T2")])
      ph2 <- d$ph; ph2$T1[1:4] <- ph2$T1[1:4]+1000
      if(route!="ml") ph2$T2[1:4] <- ph2$T2[1:4]-2000
      args$pheno_object <- ph2
      expect_error(do.call(fn,args),"joint view captured")
      expect_identical(captured$pheno_visible,visible)
      if(route=="ml") {
        args$pheno_object$T2[1] <- args$pheno_object$T2[1]+50
        expect_error(do.call(fn,args),"joint view captured")
        expect_false(identical(captured$pheno_visible,visible))
      }
    }
  }
})

test_that("hidden-cell changes leave fold-internal BGLR rankings unchanged", {
  d <- .joint_cv_data(); ph <- d$ph; ph[1:4,c("T1","T2")] <- NA_real_
  view <- function(x) gp_feature_joint_view(d$X,x,c("T1","T2"),"ID",2L,
    selection_mode="fold_internal",bayes_nIter=200L,bayes_burnIn=100L,bayes_thin=2L)
  a <- view(ph)
  changed <- d$ph; changed[1:4,c("T1","T2")] <- 1e6
  changed[1:4,c("T1","T2")] <- NA_real_
  b <- view(changed)
  expect_identical(a$metadata$raw_score,b$metadata$raw_score)
  expect_identical(a$selected_predictors,b$selected_predictors)
  expect_identical(unique(a$metadata$n_train),8L)
})


test_that("ML ASReml and DL folds use the selected columns and actual k", {
  d <- .joint_cv_data(); fit_columns <- list(); kernel_columns <- list(); scored <- list()
  local_mocked_bindings(
    gp_tuning_folds=function(...)list(1:4,5:8,9:12),
    gp_asreml_require_engine=function(...)invisible(NULL),
    feature_score_predictors=function(predictor_data,pheno_data,...) {
      scored[[length(scored)+1L]] <<- pheno_data
      d$meta
    },
    gp_multitrait_ml_fit_trait=function(X_train,X_test,...) {
      expect_identical(colnames(X_train),colnames(X_test))
      fit_columns[[length(fit_columns)+1L]] <<- colnames(X_train)
      list(job=list(),scaler=NULL)
    },
    gp_py_ml_fit_predict_batch=function(jobs,...)lapply(jobs,function(job) {
      built <- job$build()
      expect_identical(colnames(built$X_train),colnames(built$X_test))
      seq_len(nrow(built$X_test))/10
    }),
    gp_multitrait_ml_job_predictions=function(result,scaler)result,
    gp_feature_genomic_kernel_from_view=function(predictor_data,...) {
      kernel_columns[[length(kernel_columns)+1L]] <<- colnames(predictor_data)
      tcrossprod(predictor_data)
    },
    gp_multitrait_asreml_gaussian_model=function(pheno_object,gmatrix,...) {
      expect_true(all(is.na(pheno_object[which(rowSums(is.na(pheno_object[,c("T1","T2")]))==2),c("T1","T2")])) )
      expect_equal(dim(gmatrix),c(12L,12L))
      list(predicted_values=data.frame(ID=rep(d$ph$ID,2),Trait=rep(c("T1","T2"),each=12),Predicted_value=seq_len(24)/10))
    },
    gp_dl_bridge_fit_predict=function(X_train,X_test,y_train,...) {
      expect_identical(colnames(X_train),colnames(X_test))
      expect_equal(nrow(X_train),8L)
      expect_equal(nrow(y_train),8L)
      fit_columns[[length(fit_columns)+1L]] <<- colnames(X_train)
      matrix(seq_len(nrow(X_test)*2)/10,nrow(X_test),2)
    },.package="PredictProR")
  for(mode in c("fixed","fold_internal")) for(k in list(2L,"all",9L)) {
    expected <- colnames(d$X)[seq_len(if(identical(k,"all"))5L else min(k,5L))]
    for(route in c("ml","asreml","dl")) {
      fit_columns <- kernel_columns <- scored <- list()
      base <- list(pheno_object=d$ph,response=c("T1","T2"),gen_name="ID",nfolds=3L,
        eval_metrics="root_mean_squared_error",feature_score_metadata=d$meta,
        feature_k=k,feature_scoring_cv=mode)
      args <- switch(route,ml=c(base,list(model_type="Ridge_Regression",geno_omic_object=d$X)),
        dl=c(base,list(model_type="mlp",geno_omic_object=d$X,scaling=FALSE,centering=FALSE)),
        asreml=c(base,list(gmatrix=d$K,feature_predictor_data=d$X)))
      fn <- switch(route,ml=gp_multitrait_ml_gaussian_cv,
                   dl=gp_multitrait_dl_gaussian_cv,asreml=gp_multitrait_asreml_gaussian_cv)
      out <- do.call(fn,args)
      columns <- if(route=="asreml")kernel_columns else fit_columns
      expect_length(columns,if(route=="ml")6L else 3L)
      for(cols in columns) expect_identical(cols[grepl("^m[1-5]$",cols)],expected,info=paste(route,mode,k))
      meta <- out$cv_results_processed$feature_selection_metadata
      expect_identical(unique(meta$trait),"T1+T2")
      group <- if(route=="ml")interaction(meta$held_out_trait,meta$fold) else meta$fold
      expect_true(all(tapply(meta$selected,group,sum)==length(expected)))
      expect_identical(unique(meta$feature_k),as.integer(length(expected)))
      expect_identical(unique(out$cv_results[[1L]]$ypred_cv_Reps_all$feature_k),as.integer(length(expected)))
      expect_length(scored,if(mode=="fixed")0L else if(route=="ml")6L else 3L)
    }
  }
})


test_that("fixed rankings learned from the full panel retain hidden-outcome dependence", {
  d <- .joint_cv_data()
  score <- function(ph) feature_score_predictors(d$X,ph,c("T1","T2"),"ID",joint=TRUE,
    bayes_nIter=200L,bayes_burnIn=100L,bayes_thin=2L)
  a <- score(d$ph)
  changed <- d$ph; changed$T1[1:4] <- 100*d$X[1:4,"m5"]
  b <- score(changed)
  masked <- d$ph; masked[1:4,c("T1","T2")] <- NA_real_
  va <- gp_feature_joint_view(d$X,masked,c("T1","T2"),"ID",2L,feature_score_metadata=a)
  vb <- gp_feature_joint_view(d$X,masked,c("T1","T2"),"ID",2L,feature_score_metadata=b)
  expect_false(identical(va$metadata$raw_score,vb$metadata$raw_score))
  expect_identical(va$metadata$raw_score,a$raw_score)
  expect_identical(vb$metadata$raw_score,b$raw_score)
})
