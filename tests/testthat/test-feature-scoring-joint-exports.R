test_that("native exports retain the full joint ranking and saved selection", {
  out_root <- tempfile("joint-export-")
  dir.create(out_root)
  withr::local_options(list(PredictProR.output_dir=out_root))
  caller_wd <- getwd()
  meta <- structure(data.frame(
    trait=rep("T1+T2",3L),predictor=c("a","b","c"),source_block="geno_data",
    rank=1:3,raw_score=c(0.7,0.4,0.2),scoring_scope="joint",
    selected=c(TRUE,TRUE,FALSE),feature_k=2L,stringsAsFactors=FALSE),
    class=c("predictpror_feature_scores","data.frame"))
  scores <- meta[,setdiff(names(meta),c("selected","feature_k")),drop=FALSE]
  keep <- c("a","b")
  res <- results_handling(
    GS_model="GBLUP_BRR",res_model_output=list(
      feature_selection_metadata=meta,feature_selected=keep,feature_k=2L),
    system_database=FALSE,plot_filename="joint-export",feature_selected=keep,
    feature_score_metadata=scores)
  expect_identical(getwd(),caller_wd)
  expect_identical(res$export_status,"Successful")
  csv <- file.path(res$export_directory,"feature_selection_metadata.csv")
  expect_true(file.exists(csv))
  if(file.exists(csv)) {
    saved <- read.csv(csv,check.names=FALSE,stringsAsFactors=FALSE)
    expect_identical(names(saved),names(meta))
    expect_identical(saved$predictor,meta$predictor)
    expect_identical(saved$selected,meta$selected)
    expect_identical(saved$scoring_scope,meta$scoring_scope)
    expect_equal(saved$raw_score,meta$raw_score,tolerance=1e-14)
  }
  e <- new.env(parent=emptyenv())
  load(file.path(res$export_directory,"feature_selected.RData"),envir=e)
  expect_identical(e$feature_selected,keep)
  load(file.path(res$export_directory,"feature_score_metadata.RData"),envir=e)
  expect_identical(e$feature_score_metadata,scores)
  expect_identical(res$model_results$feature_selection_metadata,meta)
})
