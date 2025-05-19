
scale_omic_data <- function(X) {
  scaler <- caret::preProcess(X, method = c("center", "scale"))
  X_scaled <- predict(scaler, X)
  rownames(X_scaled) <- rownames(X)
  list(scaled = X_scaled, scaler = scaler)
}


rf_importance_stat <- function(data, indices, response_col = "y", mtry = NULL) {
  df <- data[indices, ]
  y <- df[[response_col]]
  X <- df[setdiff(names(df), response_col)]

  if (is.null(mtry)) {
    mtry <- floor(sqrt(ncol(X)))
  }

  rf_model <- randomForest::randomForest(x = X, y = y, mtry = mtry, importance = TRUE)
  imp <- randomForest::importance(rf_model, type = 2)[, 1]
  full_names <- colnames(X)
  imp_full <- setNames(rep(0, length(full_names)), full_names)
  imp_full[names(imp)] <- imp

  return(imp_full)
}

#####
boot_rf_feature_selection <- function(X, y, R = 100, seed = 123, mtry = 500) {
  set.seed(seed)
  df <- as.data.frame(X)
  df$y <- y

  boot_obj <- boot::boot(data = df, statistic = rf_importance_stat, R = R, response_col = "y", mtry = mtry)

  importance_matrix <- boot_obj$t
  colnames(importance_matrix) <- colnames(X)

  mean_scores <- colMeans(importance_matrix)
  sd_scores <- apply(importance_matrix, 2, sd)
  cv_scores <- sd_scores / (mean_scores + 1e-8)
  composite_scores <- mean_scores / (cv_scores + 1e-8)

  importance_summary <- data.frame(
    feature = colnames(X),
    mean = mean_scores,
    sd = sd_scores,
    cv = cv_scores,
    composite_score = composite_scores,
    stringsAsFactors = FALSE
  )

  #importance_summary <- importance_summary[order(-importance_summary$composite_score), ]

  importance_summary <- importance_summary[order(-importance_summary$mean), ]



  return(list(importance_summary = importance_summary))
}

run_shrinkage_gamm <- function(X, y, k_value = 5, select = FALSE, method = "fREML") {
  feature_names <- colnames(X)
  gam_formula <- as.formula(paste(
    "y ~",
    paste(sprintf("s(%s, k = %d)", feature_names, k_value), collapse = " + ")
  ))
  gam_data <- data.frame(y = y, X)
  gam_fit <- mgcv::bam(gam_formula, data = gam_data, select = select, method = method)
  return(list(model = gam_fit))
}

###
transform_features <- function(X, y = NULL, method, selected, X_test = NULL, n_comp = 50, var_explained = 0.9) {
  X_sel <- X[, selected, drop = FALSE]
  if (method == "pca") {
    pca <- prcomp(X_sel)
    cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    n_pc <- which(cum_var >= var_explained)[1]
    train_scores <- pca$x[, 1:n_pc, drop = FALSE]
    colnames(train_scores) <- paste0("pc", 1:ncol(train_scores))
    test_scores <- if (!is.null(X_test)) predict(pca, X_test[, selected, drop = FALSE])[, 1:n_pc, drop = FALSE] else NULL
    if(!is.null(test_scores)) colnames(test_scores) <- colnames(train_scores)
  } else if (method == "pls") {
    pls_model <- pls::plsr(y ~ ., data = data.frame(y = y, X_sel), ncomp = n_comp)
    train_scores <- pls::scores(pls_model)[, 1:n_comp, drop = FALSE]
    colnames(train_scores) <- paste0("comp", 1:ncol(train_scores))
    test_scores <- if (!is.null(X_test)) scale(X_test[, selected, drop = FALSE], center = pls_model$Xmeans, scale = FALSE) %*% pls_model$loadings[, 1:n_comp] else NULL
    if(!is.null(test_scores)) colnames(test_scores) <- colnames(train_scores)
  }
  return(list(train_scores = train_scores, test_scores = test_scores))
}


run_shrinkage_rf_gamm <- function(X, y,X_test = NULL, geno_scaler = NULL, k_value = 5, select = FALSE, method = "fREML") {
  min_k <- 3

  unique_counts <- sapply(as.data.frame(X), function(x) length(unique(x)))
  valid_cols <- names(unique_counts[unique_counts >= min_k])
  if (length(valid_cols) < ncol(X)) {
    removed <- setdiff(colnames(X), valid_cols)
    cat("Removed", length(removed), "features with <3 unique values:", paste(removed, collapse = ", "), "\n")
    X <- X[, valid_cols, drop = FALSE]
    if(!is.null(X_test) && !is.null(geno_scaler)){
      X_test <- predict(geno_scaler, X_test)
    }
    test_scores <- if (!is.null(X_test)) X_test[, valid_cols, drop = FALSE] else NULL
  }

  bs <- if (max(unique_counts) <= 3) "cr" else "tp"
  if (bs == "cr") k_value <- min_k

  feature_names <- colnames(X)
  terms <- paste0("s(", feature_names, ", bs='", bs, "', k=", k_value, ")", collapse = " + ")
  gam_formula <- as.formula(paste("y ~", terms))

  gam_data <- data.frame(y = y, X)
  gam_fit <- mgcv::bam(gam_formula, data = gam_data, select = select, method = method)
  return(list(model = gam_fit, test_scores = test_scores))
}


gam_new <- function(pheno_object,
                    selected,
                    response,
                    CV = TRUE,
                    k_value = 5,
                    geno_omic_object,
                    geno_omic_test_object = NULL,
                    var_explained = 0.9,
                    gam_method = c("ND_mod1", "ND_mod2",
                                   "ND_mod3", "ND_mod4",
                                   "ND_mod5", "ND_mod6"),
                    seed = 123,
                    select = TRUE,
                    max_features = 50){


  test_scores <- NULL
  gam_method <- match.arg(gam_method)

  gam_method_use = c("ND_mod1", "ND_mod2",
                     "ND_mod3", "ND_mod4",
                     "ND_mod5", "ND_mod6")

  gam_method_original = c("pls_gam", "pls_gam_select",
                          "pc_gam", "pc_gam_select",
                          "rf_gam", "rf_gam_select")


  gam_map <- setNames(gam_method_original, gam_method_use)


  get_gam_model <- function(user_input) {
    user_input <- match.arg(user_input, choices = names(gam_map))
    return(gam_map[[user_input]])
  }

  gam_mod <- get_gam_model(gam_method)

  if(CV){
    y <-  pheno_object
    }else{
      y <-  pheno_object[, response]
    }

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

                        stop("Unsupported GAM model method.")
    )

    if("test_scores"%in%names(gam_model)) {
      test_scores <- gam_model$test_scores
      gam_model <- gam_model[!names(gam_model)%in%"test_scores"]

    }

    return((list(model = gam_model,
                 test_scores = test_scores)))

  }

  predict_gam_model <- function(model, SE_fit = TRUE) {
    if (!is.null(model)) {
      X_test <- tryCatch({
        if (!is.null(model$test_scores)) {
          model$test_scores
        } else {
          stop("test features was not provided.")
        }
      }, error = function(e) {
        message("Prediction error: ", e$message)
        return(NULL)
      })

      predictions <- tryCatch({
        predict(model$model,
                newdata = data.frame(X_test),
                type = "response", se.fit = SE_fit)
      }, error = function(e) {
        message("Error in prediction: ", e$message)
        NULL
      })
    } else {
      return(NULL)
    }
    return(predictions)
  }







