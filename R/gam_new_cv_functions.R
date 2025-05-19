
ND_modes_cv <- function(y,
                      selected,
                      tst = NULL,
                      CV = TRUE,
                      k_value = 5,
                      omics,
                      var_explained = 0.9,
                      gam_method = c("ND_mod1", "ND_mod2"),
                      seed = 123,
                      select = TRUE,
                      max_features = 50){

#browser()
  test_scores <- NULL
  gam_method <- match.arg(gam_method)
  if("ND_mod1"==gam_method ||"ND_mod2"==gam_method){
    geno_omic_object_use <- omics
    yy <- y
  }


  # gam_method_use = c("ND_mod1", "ND_mod2",
  #                    "ND_mod3", "ND_mod4",
  #                    "ND_mod5", "ND_mod6", "ND_mod7", "ND_mod8")

  gam_method_use = c("ND_mod1", "ND_mod2")

  gam_method_original = c("gam_svm", "gam_svmm")

  # gam_method_original = c("pls_gam", "pls_gam_select",
  #                         "pc_gam", "pc_gam_select",
  #                         "rf_gam", "rf_gam_select",
  #                         "gam_svm", "gam_svmm")


  gam_map <- setNames(gam_method_original, gam_method_use)


  get_gam_model <- function(user_input) {
    user_input <- match.arg(user_input, choices = names(gam_map))
    return(gam_map[[user_input]])
  }

  gam_mod <- get_gam_model(gam_method)

  y <-  y[-tst]
  geno_omic_object <- omics[-tst, ]
  geno_omic_test_object <- omics[tst, ]

  if (length(selected) == 0) stop("No features were selected for modeling.")
  if (any(is.na(y))) stop("Response variable y contains missing values.")

  if(any(gam_mod%in%c("pc_gam", "pc_gam_select"))){
    pc_out <-  transform_features(X=geno_omic_object,
                                  X_test = geno_omic_test_object,
                                  method = "pca",
                                  selected = selected,
                                  n_comp = max_features)

    test_scores <- pc_out$test_scores
  }

  if(any(gam_mod%in%c("pls_gam", "pls_gam_select"))){
    pls_out <-  transform_features(X=geno_omic_object,
                                   y = y,
                                   X_test = geno_omic_test_object,
                                   method = "pls",
                                   selected = selected,
                                   n_comp = max_features)

    test_scores <- pls_out$test_scores
  }

  if(any(gam_mod%in%c("rf_gam", "rf_gam_select"))){
    scaled <- scale_omic_data(geno_omic_object[, selected[1:max_features], drop =FALSE])
    geno_scaled <- scaled$scaled
    geno_scaler <- scaled$scaler
  }

  gam_model <- switch(gam_mod,
                      "pc_gam" = tryCatch({
                        run_shrinkage_gamm(pc_out$train_scores, y, select = FALSE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "pc_gam_select" = tryCatch({
                        run_shrinkage_gamm(pc_out$train_scores, y, select = TRUE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "pls_gam" = tryCatch({
                        run_shrinkage_gamm(pls_out$train_scores, y, select = FALSE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "pls_gam_select" = tryCatch({
                        run_shrinkage_gamm(pls_out$train_scores, y, select = TRUE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "rf_gam" = tryCatch({
                        run_shrinkage_rf_gamm(X= geno_scaled, y= y, X_test = geno_omic_test_object, geno_scaler = geno_scaler, select = FALSE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "rf_gam_select" = tryCatch({
                        run_shrinkage_rf_gamm(X= geno_scaled, y= y, X_test = geno_omic_test_object, geno_scaler = geno_scaler, select = TRUE, k_value = k_value)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "gam_svm" = tryCatch({
                        AI_svm_cv(omics= geno_omic_object_use[, selected[1:max_features], drop =FALSE], y= yy, tst = tst, omic_count = NULL)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),

                      "gam_svmm" = tryCatch({
                        AI_svm_cv(omics= geno_omic_object_use[, selected, drop =FALSE], y= yy, tst = tst, omic_count = NULL)
                      }, error = function(e) {
                        message(paste("Error in", gam_method, ": "), e$message)
                        NULL
                      }),
                      stop("Unsupported GAM model method.")
  )

  if(gam_mod=="gam_svm" || gam_mod=="gam_svmm"){
    return(gam_model)
  }

  # if("test_scores"%in%names(gam_model)) {
  #   test_scores <- gam_model$test_scores
  #   gam_model <- gam_model[!names(gam_model)%in%"test_scores"]
  #
  # }
  if(!"test_scores"%in%names(gam_model)) {
  gam_model[["test_scores"]] <- test_scores
  }
  preds <- predict_gam_model(gam_model, SE_fit = FALSE)
  return(preds)

}

