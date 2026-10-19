test_that("joint selected and all kernels equal independent raw-input formulas", {
  ids <- paste0("g",1:6)
  G <- matrix(c(0,0,1,1,2,0, 0,1,2,2,1,2,
                0,0,0,1,0,2, 2,0,1,0,1,0),6,4,
    dimnames=list(ids,paste0("m",1:4)))
  O <- matrix(c(2,4,3,8,6,1, 1,1,3,4,2,8, 5,2,8,1,7,4),6,3,
    dimnames=list(ids,paste0("o",1:3)))
  # The merged predictor view may have been scaled for ML. The GRM must still
  # use the retained raw dosages, and each omics kernel its own raw source.
  merged <- cbind(scale(G),scale(O))
  sources <- setNames(c(rep("geno_data",4),rep("omic1_data",3)),colnames(merged))
  meta <- gp_feature_metadata_table("T1+T2",colnames(merged),sources,
    c(1,7,2,6,3,4,5),"Ridge_Regression",1L,6L,scoring_scope="joint")
  ph <- data.frame(ID=ids,T1=1:6,T2=c(2,4,1,6,3,5))
  for(k in list(3L,"all")) {
    result <- gp_feature_joint_kernel_bank(merged,ph,c("T1","T2"),"ID",meta,k,
      sources,list(geno_data=G,omic1_data=O),gmatrix_method="VanRaden",
      kernel_method="Linear_kernel",ploidy=2L,scaling=FALSE,centering=FALSE)
    keep <- if(identical(k,"all"))c("m2","m4","o3","o2","o1","m3","m1") else c("m2","m4","o3")
    expect_identical(result$view$selected_predictors,keep)
    g <- G[,intersect(keep,colnames(G)),drop=FALSE]
    freq <- colMeans(g)/2
    Z <- sweep(g,2,2*freq,"-")
    expected_g <- Z %*% t(Z)/sum(2*freq*(1-freq))
    observed_g <- result$kernel_bank$gmatrix_model_ready
    expect_lt(max(abs(observed_g-expected_g)),1e-12)
    expect_identical(dimnames(observed_g),list(ids,ids))
    o <- O[,intersect(keep,colnames(O)),drop=FALSE]
    expected_o <- o %*% t(o)/ncol(o)
    observed_o <- result$kernel_bank$omic1_kernel_model_ready
    expect_lt(max(abs(observed_o-expected_o)),1e-12)
    expect_identical(dimnames(observed_o),list(ids,ids))
    expect_setequal(names(result$kernel_bank),c("gmatrix_model_ready","omic1_kernel_model_ready"))
  }
})
