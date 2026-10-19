test_that("disabled auto policy retains defaults and explicit components", {
  ph <- data.frame(GID=rep(letters[1:4],2),Env=rep(c("a","b"),each=4),y=1:8)
  for (env_kernel in c(FALSE,TRUE)) for (enabled in c(FALSE,TRUE)) {
    for (explicit in list(NULL, "ge", character())) {
      out <- PredictProR:::gp_single_trait_auto_policy(ph,
        split=list(train_mask=rep(TRUE,8),test_mask=rep(FALSE,8)),
        heter_groups="Env",gid_col="GID",env_col="Env",
        env_similarity=if(env_kernel) diag(2) else NULL,
        enabled=enabled,include_components=explicit,resolve_device=FALSE)
      expected <- if(!is.null(explicit)) explicit else if(env_kernel) c("g","ge","e") else "g"
      expect_identical(out$include_components,expected)
    }
  }
})
