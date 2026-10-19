.review219_scores <- function(predictors, scores, trait = "T1+T2", extra = NULL) {
  gp_feature_metadata_table(trait, predictors, setNames(rep("geno_data", length(predictors)), predictors),
    scores, "RandomForest", 1L, 12L, scoring_scope = "joint", extra = extra)
}

test_that("supplied tied ranks survive collation changes without allowing score inversions", {
  withr::local_collate("C")
  m <- .review219_scores(c("Z", "a"), c(0, 0))
  m <- m[match(c("a", "Z"), m$predictor), ]; m$rank <- 1:2
  expect_no_error(gp_feature_joint_rows(m, c("T1", "T2")))
  m$raw_score <- c(0, 1)
  expect_error(gp_feature_joint_rows(m, c("T1", "T2")), "scores")
})

test_that("new score ties use a locale-independent predictor-name order", {
  old <- Sys.getlocale("LC_COLLATE"); on.exit(Sys.setlocale("LC_COLLATE", old))
  for (locale in c("C", "English_United States.utf8", "en_US.UTF-8")) {
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) next
    m <- .review219_scores(c("a", "Z", "a_b", "ab"), rep(0, 4))
    expect_identical(m$predictor, sort(c("a", "Z", "a_b", "ab"), method = "radix"))
  }
})

test_that("constant RF predictors follow signed active scores in both scoring scopes", {
  withr::local_seed(219)
  X <- matrix(rnorm(40*4), 40, dimnames=list(paste0("g",1:40),c("m1","m2","m3","mono")))
  X[,4] <- 2
  ph <- data.frame(ID=rownames(X),T1=rnorm(40),T2=rnorm(40))
  local_mocked_bindings(gp_py_ml_feature_importance=function(...) c(-.5,-.2,-.1,0), .package="PredictProR")
  for (joint in c(FALSE, TRUE)) {
    m <- feature_score_predictors(X,ph,if(joint)c("T1","T2") else "T1","ID",scoring_model="RandomForest",joint=joint)
    expect_identical(m$predictor,c("m3","m2","m1","mono"))
    expect_equal(m$raw_score,c(-.1,-.2,-.5,0))
    expect_true(all(is.finite(m$normalized_score)))
    expect_identical(tail(m$scoring_status,1),"no_predictor_variation")
    if(joint) expect_no_error(gp_feature_joint_rows(m,c("T1","T2")))
  }
})

test_that("active zero and inactive zero are different selection priorities", {
  m <- .review219_scores(c("constant","variable"),c(0,0),extra=list(scoring_status=c("no_predictor_variation","estimated")))
  gp_feature_reset_tie_warning(); on.exit(gp_feature_reset_tie_warning())
  expect_no_warning(top <- gp_feature_top_k_rows(m,1))
  expect_identical(top$predictor,"variable")
  expect_false(isTRUE(attr(top,"tied_at_cutoff")))
})

test_that("joint matrix mismatch reports bounded examples and both counts", {
  names <- sprintf("dropped_%03d",1:40)
  m <- .review219_scores(names,40:1)
  X <- matrix(1,3,2,dimnames=list(letters[1:3],c("retained_a","retained_b")))
  msg <- tryCatch(gp_feature_joint_view(X,NULL,c("T1","T2"),"ID",1,feature_score_metadata=m),error=conditionMessage)
  expect_match(msg,"40.*ranking.*model matrix")
  expect_match(msg,"2.*model matrix.*ranking")
  expect_match(msg,"post-QC")
  expect_match(msg,"qc_filtering = FALSE",fixed=TRUE)
  expect_lte(length(regmatches(msg,gregexpr("dropped_[0-9]+",msg))[[1]]),10)
})

test_that("joint selection tables and singleton lists round-trip shared selections", {
  m <- .review219_scores(c("a","b","c"),3:1,trait="T2+T1");m$selected<-c(TRUE,TRUE,FALSE)
  selected<-gp_feature_joint_selected(m,c("T1","T2"))
  expect_identical(gp_feature_selected_vector(selected,source="geno_data"),c("a","b"))
  expect_identical(gp_feature_joint_selected(list(c("a","b")),c("T1","T2")),c("a","b"))
  expect_error(gp_feature_joint_selected(list(T1="a",T2="b"),c("T1","T2")),"one predictor set")
  m$trait<-"T1+T3"
  expect_error(gp_feature_joint_selected(m,c("T1","T2")),"trait|Trait")
})

test_that("hybrid ML and DL CV record actual all and capped predictor counts", {
  withr::local_seed(219)
  ph<-expand.grid(Female=paste0("F",1:4),Male=paste0("M",1:3),stringsAsFactors=FALSE)
  ph$ID<-paste0("h",seq_len(nrow(ph)));ph$Y<-rnorm(nrow(ph))
  X<-matrix(rnorm(nrow(ph)*3),nrow(ph),3,dimnames=list(ph$ID,c("a","b","c")))
  m<-.review219_scores(colnames(X),3:1,trait="Y")
  fit<-function(X_train,y_train,X_test,...)list(test_pred=rep(mean(y_train),nrow(X_test)))
  local_mocked_bindings(gp_hybrid_ml_fit_response=fit,gp_hybrid_dl_fit_response=fit,.package="PredictProR")
  for(route in c("ml","dl")) for(k in list("all",99L)) {
    f<-if(route=="ml")gp_hybrid_ml_gaussian_cv else gp_hybrid_dl_gaussian_cv
    expect_no_warning(z<-f(model_type=if(route=="ml")"Ridge_Regression" else "mlp",
      pheno_object=ph,response="Y",gen_name="ID",female_parent="Female",male_parent="Male",
      geno_omic_object=X,nfolds=2L,replication=1L,feature_score_metadata=m,feature_k=k,
      eval_metrics="root_mean_squared_error"))
    for(rep in z$cv_results) {
      expect_true(nrow(rep$ypred_cv_Reps_all)>0)
      expect_identical(unique(rep$ypred_cv_Reps_all$feature_k),3L)
      expect_identical(unique(rep$eval_metrics_reps$feature_k),3L)
    }
  }
})

test_that("multi-response parent-hybrid rankings stop before preprocessing", {
  local_mocked_bindings(gp_shared_preprocess_cache_validate=function(...) stop("PREPROCESSING_REACHED"),.package="PredictProR")
  ph<-data.frame(ID=letters[1:4],T1=1:4,T2=4:1,Female=rep(c("f1","f2"),2),Male=rep(c("m1","m2"),each=2))
  X<-matrix(0:7,4,2,dimnames=list(c("f1","f2","m1","m2"),c("a","b")))
  m<-rbind(.review219_scores(c("a","b"),c(2,1),"T1"),.review219_scores(c("a","b"),c(1,2),"T2"))
  expect_error(model_execute(pheno_data=ph,geno_data=X,response=c("T1","T2"),gen_name="ID",
    female_parent="Female",male_parent="Male",hybrid_gp=TRUE,GS_model="GP",
    feature_score_metadata=m,feature_k=1,response_family="gaussian",message=FALSE),"shared.*feature_selected")
})
