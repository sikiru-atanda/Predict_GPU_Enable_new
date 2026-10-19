.joint_view_fixture <- function() {
  X <- matrix(seq_len(40),8,5,dimnames=list(paste0("L",1:8),letters[1:5]))
  meta <- gp_feature_metadata_table("T2+T1",colnames(X),
            setNames(rep("geno_data",5),colnames(X)),c(1,5,4,3,2),
            "Ridge_Regression",1L,8L,scoring_scope="joint")
  list(X=X,ph=data.frame(ID=rownames(X),T1=1:8,T2=c(1,3,2,5,4,7,6,8)),meta=meta)
}

test_that("joint fixed view selects one exact shared set with reversed labels accepted", {
  d <- .joint_view_fixture()
  v <- gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",2L,
                             feature_score_metadata=d$meta,replication=2L,fold=3L,
                             held_out_trait="T1")
  expect_identical(colnames(v$predictor_data),c("b","c"))
  expect_identical(v$selected_predictors,c("b","c"))
  expect_identical(v$feature_k,2L)
  expect_equal(sum(v$metadata$selected),2)
  expect_identical(unique(v$metadata$held_out_trait),"T1")
  expect_false(any(c("selected_by_trait","joint_union_k")%in%names(v)))
  expect_false(any(c("joint_union_selected","joint_union_k")%in%names(v$metadata)))
})

test_that("joint all/large k and cutoff ties retain exact cardinality and flags", {
  d <- .joint_view_fixture()
  gp_feature_reset_tie_warning()
  expect_no_warning(v <- gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID","all",
                                               feature_score_metadata=d$meta))
  expect_identical(v$feature_k,5L)
  expect_no_warning(gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",50L,
                                          feature_score_metadata=d$meta))
  d$meta$raw_score <- c(5,4,4,2,1)
  d$meta$rank <- seq_len(nrow(d$meta))
  gp_feature_reset_tie_warning()
  expect_warning(v <- gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",2L,
                                             feature_score_metadata=d$meta),"cut-off")
  expect_equal(sum(v$metadata$selected),2L)
  expect_true(all(v$metadata$tied_at_cutoff))
})

test_that("joint ranking validation refuses ambiguous/missing/per-trait inputs", {
  d <- .joint_view_fixture()
  bad <- d$meta; bad$scoring_scope <- "trait"
  expect_error(gp_feature_joint_rows(bad,c("T1","T2")),"joint = TRUE")
  bad <- d$meta; bad$predictor[2] <- bad$predictor[1]
  expect_error(gp_feature_joint_rows(bad,c("T1","T2")),"unique")
  bad <- d$meta; bad$raw_score[1] <- NA_real_
  expect_error(gp_feature_joint_rows(bad,c("T1","T2")),"finite")
  expect_error(gp_feature_joint_rows(d$meta,c("T1","T3")),"responses|traits")
  expect_error(gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",2L,
                                     feature_score_metadata=d$meta[-1,]),"every predictor")
  bad <- d$meta; bad$predictor[1] <- "absent"
  expect_error(gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",2L,
                                     feature_score_metadata=bad),"absent")
})

test_that("joint fold view forwards only visible cells and joint scoring arguments", {
  d <- .joint_view_fixture(); d$ph$T1[1:2] <- NA_real_
  seen <- NULL
  local_mocked_bindings(feature_score_predictors=function(...) {
    seen <<- list(...)
    d$meta
  },.package="PredictProR")
  v <- gp_feature_joint_view(d$X,d$ph,c("T1","T2"),"ID",3L,
                             selection_mode="fold_internal",seed=99L,held_out_trait="T1")
  expect_true(seen$joint)
  expect_identical(seen[[2]],d$ph)
  expect_identical(seen$seed,99L)
  expect_identical(seen[[3]],c("T1","T2"))
  expect_identical(v$feature_k,3L)
})

test_that("ML/DL apply selects the same columns for training and target matrices", {
  d <- .joint_view_fixture()
  dat <- list(merged_data=list(merge_data=d$X),merged_data_test=d$X[,5:1,drop=FALSE],
              feature_selected_by_trait=list(T1="a",T2="e"))
  out <- gp_feature_apply_multitrait_to_ml_dat_res(dat,d$ph,c("T1","T2"),"ID",d$meta,2L)
  expect_identical(colnames(out$merged_data$merge_data),c("b","c"))
  expect_identical(colnames(out$merged_data_test),c("b","c"))
  expect_identical(out$feature_selected,c("b","c"))
  expect_null(out$feature_selected_by_trait)
  expect_equal(sum(out$feature_selection_metadata$selected),2)
})

test_that("joint kernel bank rebuilds from exactly the shared selected columns", {
  d <- .joint_view_fixture()
  seen <- NULL
  local_mocked_bindings(gp_feature_build_kernel_bank=function(selected_sources,...) {
    seen <<- selected_sources
    list(gmatrix_model_ready=diag(nrow(d$X)))
  },.package="PredictProR")
  out <- gp_feature_joint_kernel_bank(
    d$X,d$ph,c("T1","T2"),"ID",d$meta,2L,
    source_map=setNames(rep("geno_data",5),colnames(d$X)),
    source_matrices=list(geno_data=d$X),gmatrix_method="VanRaden")
  expect_identical(colnames(seen$geno_data),c("b","c"))
  expect_identical(out$view$selected_predictors,c("b","c"))
})
