
scale_omic_data <- function(X) {
  scaler <- caret::preProcess(X, method = c("center", "scale"))
  X_scaled <- predict(scaler, X)
  rownames(X_scaled) <- rownames(X)
  list(scaled = X_scaled, scaler = scaler)
}

get_r2 <- function(observed, predicted) {
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  1 - ss_res / ss_tot
}

get_correlation <- function(observed, predicted) {
  cor(observed, predicted, use = "complete.obs")
}

get_hybrid <- function(observed, predicted) {
  r2 <- get_r2(observed, predicted)
  corr <- get_correlation(observed, predicted)
  return((r2 + corr) / 2)
}



cv_model_score <- function(X, y, metric = c("r2", "correlation", "hybrid"), k = 5, seed = 123) {
  set.seed(seed)
  folds <- caret::createFolds(y, k = k, list = TRUE)
  scores <- numeric(k)

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)

    model <-run_shrinkage_rf_gamm(X[train_idx, , drop = FALSE], y[train_idx],
                                  X_test = X[test_idx, , drop = FALSE])
    preds <- predict(model$model, newdata = data.frame(model$test_scores))
    obs <- y[test_idx]

    scores[i] <- switch(metric[1],
                        r2 = get_r2(obs, preds),
                        correlation = get_correlation(obs, preds),
                        hybrid = get_hybrid(obs, preds),
                        stop("Invalid metric")
    )
  }
  mean(scores, na.rm = TRUE)
}

pls_gam_transform <- function(X_train, y, X_test = NULL, n_comp = NULL) {
  df_train <- data.frame(y = y, X_train)

  # Train PLS with cross-validation if n_comp is not specified
  if (is.null(n_comp)) {
    pls_model <- pls::plsr(y ~ ., data = df_train, validation = "CV")
    n_comp <- which.min(pls_model$validation$PRESS)
  } else {
    pls_model <- pls::plsr(y ~ ., data = df_train, ncomp = n_comp)
  }

  # Extract training scores
  train_scores <- pls::scores(pls_model)[, 1:n_comp, drop = FALSE]
  colnames(train_scores) <- paste0("comp", 1:ncol(train_scores))

  # Project test set (if provided)
  test_scores <- NULL
  if (!is.null(X_test)) {
    # Make sure test set has same column order and names
    X_test_scaled <- scale(X_test, center = pls_model$Xmeans, scale = FALSE)
    test_scores <- as.matrix(X_test_scaled) %*% pls_model$loadings[, 1:n_comp, drop = FALSE]
    colnames(test_scores) <- colnames(train_scores)
  }

  return(list(
    train_scores = as.data.frame(train_scores),
    test_scores = if (!is.null(test_scores)) as.data.frame(test_scores) else NULL,
    n_comp = n_comp,
    pls_model = pls_model
  ))
}




# sik = pls_gam_transform(X_train =  X[train_idx, , drop = FALSE],
#                         y= y[train_idx],
#                         X_test = X[test_idx, , drop = FALSE])
#
#
#
# y = pheno_data1$B_GLUCAN


#############################################

select_elbow_pc <- function(cum_var, sensitivity = 1) {
  x <- seq_along(cum_var)
  y <- cum_var

  # Line from (1, y1) to (n, yn)
  line_vec <- c(length(x) - 1, y[length(y)] - y[1])
  line_vec <- line_vec / sqrt(sum(line_vec^2))  # normalize

  dist_to_line <- sapply(1:length(x), function(i) {
    vec <- c(i - 1, y[i] - y[1])
    proj <- sum(vec * line_vec) * line_vec
    diff <- vec - proj
    sqrt(sum(diff^2))
  })

  elbow_idx <- which.max(dist_to_line)
  elbow_idx
}


# pca <- prcomp(X, scale. = TRUE)
# cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
# n_pc <- select_elbow_pc(cum_var)

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
boot_rf_feature_selection <- function(X, y, R = 100, seed = 123, mtry = 500, percentile = 0.90,
                                      cum_importance_threshold = 0.9,
                                      k_sd = 1.5,
                                      candidate_sizes = c(5, 10, 20, 50, 100),
                                      eval_metric = "correlation",
                                      method_percentile = TRUE,
                                      method_cum_importance_thres = TRUE,
                                      method_cross_val_tune = TRUE,
                                      method_sd = TRUE) {
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
    composite_score = composite_scores
  ) |>
    dplyr::arrange(dplyr::desc(composite_score))

  # Percentile-based threshold
  percentile_cutoff <- stats::quantile(importance_summary$composite_score, percentile)

  selected_features_percentile <- importance_summary |>
    dplyr::filter(composite_score >= percentile_cutoff) |>
    dplyr::pull(feature)

  ##### sd approach
  # Select features above mean + k * sd of importance scores
  threshold <- mean(importance_summary$mean) + k_sd * sd(importance_summary$mean)
  selected_features_sd <- importance_summary$feature[importance_summary$mean >= threshold]

  ##### using cummulative important threshold
  importance_summary_cum <- data.frame(
    feature = colnames(X),
    mean = mean_scores,
    sd = sd_scores,
    cv = cv_scores,
    composite_score = composite_scores
  ) |>
    dplyr::arrange(dplyr::desc(composite_score)) |>
    dplyr::mutate(cum_composite_importance = cumsum(composite_score) / sum(composite_score))

  # Select features based on cumulative composite importance
  selected_features_cum <- importance_summary_cum |>
    dplyr::filter(cum_composite_importance <= cum_importance_threshold) |>
    dplyr::pull(feature)

  ####
  importance_summary <- importance_summary[order(-importance_summary$composite_score), ]
  results_cor <- list()
  results_r2 <- list()
  results_hybrid <- list()
  for (k in candidate_sizes) {
    kk <- floor((nrow(importance_summary)*k)/100)
    if (kk > nrow(importance_summary)) next
    selected_features_cv <- importance_summary$feature[1:kk]
    X_subset <- X[, selected_features_cv, drop = FALSE]

    perf <- cv_model_score(X = X_subset, y = y, metric = "correlation")
    results_cor[[as.character(k)]] <- list(score = perf, features = selected_features_cv)
    ##
    perf_hybrid <- cv_model_score(X = X_subset, y = y, metric = "hybrid")
    results_hybrid[[as.character(k)]] <- list(score = perf_hybrid, features = selected_features_cv)
    ##
    perf_r2 <- cv_model_score(X = X_subset, y = y, metric = "r2")
    results_r2[[as.character(k)]] <- list(score = perf_r2, features = selected_features_cv)
  }

  best_k <- names(which.max(sapply(results_cor, function(x) x$score)))
  selected_features_cor_cv <- results_cor[[best_k]]$features
  #
  best_k_hybrid <- names(which.max(sapply(results_hybrid, function(x) x$score)))
  selected_features_hybrid_cv <- results_cor[[best_k_hybrid]]$features
  ##
  best_k_r2 <- names(which.max(sapply(results_r2, function(x) x$score)))
  selected_features_r2_cv <- results_cor[[best_k_r2]]$features

  # importance_summary <- data.frame(
  #   feature = colnames(X),
  #   mean = mean_scores,
  #   sd = sd_scores,
  #   cv = cv_scores,
  #   composite_score = composite_scores,
  #   stringsAsFactors = FALSE
  # )
  # #importance_summary <- importance_summary[order(-importance_summary$composite_score), ]
  # importance_summary <- importance_summary[order(-importance_summary$mean), ]





  return(list(selected_features_sd = selected_features_sd,
              selected_features_cum = selected_features_cum,
              selected_features_percentile = selected_features_percentile,
              selected_features_cor_cv = selected_features_cor_cv,
              selected_features_r2_cv = selected_features_r2_cv,
              selected_features_hybrid_cv = selected_features_hybrid_cv))
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
transform_features <- function(X, y = NULL, method, selected, X_test = NULL, n_comp = NULL
                               #var_explained = 0.9
                               ) {
  X_sel <- X[, selected, drop = FALSE]
  if (method == "pca") {
    pca <- prcomp(X_sel, scale. = TRUE)
    cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    #
    n_pc <- select_elbow_pc(cum_var)
    #n_pc <- which(cum_var >= var_explained)[1]
    train_scores <- pca$x[, 1:n_pc, drop = FALSE]
    colnames(train_scores) <- paste0("pc", 1:ncol(train_scores))
    test_scores <- if (!is.null(X_test)) predict(pca, X_test[, selected, drop = FALSE])[, 1:n_pc, drop = FALSE] else NULL
    if(!is.null(test_scores)) colnames(test_scores) <- colnames(train_scores)
  } else if (method == "pls") {
    df_train <- data.frame(y = y, X_sel)

    # Train PLS with cross-validation if n_comp is not specified
    if (is.null(n_comp)) {
      pls_model <- pls::plsr(y ~ ., data = df_train, validation = "CV")
      n_comp <- which.min(pls_model$validation$PRESS)
    } else {
      pls_model <- pls::plsr(y ~ ., data = df_train, ncomp = n_comp)
    }

    # Extract training scores
    train_scores <- pls::scores(pls_model)[, 1:n_comp, drop = FALSE]
    colnames(train_scores) <- paste0("comp", 1:ncol(train_scores))

    # Project test set (if provided)
    test_scores <- NULL
    if (!is.null(X_test)) {
      # Make sure test set has same column order and names
      X_test_scaled <- scale(X_test[, selected, drop = FALSE], center = pls_model$Xmeans, scale = FALSE)
      test_scores <- as.matrix(X_test_scaled) %*% pls_model$loadings[, 1:n_comp, drop = FALSE]
      colnames(test_scores) <- colnames(train_scores)
    }
  }
  return(list(train_scores = train_scores, test_scores = test_scores))
}


run_shrinkage_rf_gamm <- function(X, y,X_test = NULL, geno_scaler = NULL, k_value = 5, select = FALSE, method = "fREML") {
  min_k <- 3

  unique_counts <- sapply(as.data.frame(X), function(x) length(unique(x)))
  valid_cols <- names(unique_counts[unique_counts >= min_k])
  if (length(valid_cols) < ncol(X)) {
    removed <- setdiff(colnames(X), valid_cols)
    #cat("Removed", length(removed), "features with <3 unique values:", paste(removed, collapse = ", "), "\n")
    cat("features with <3 unique values:", paste(length(removed), collapse = ", "), "\n")
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
                                    selected = selected
                                    #n_comp = max_features
                                    )

      test_scores <- pc_out$test_scores
    }

    if(any(gam_mod%in%c("pls_gam", "pls_gam_select"))){
      pls_out <-  transform_features(X=geno_omic_object,
                                     y = y,
                                     X_test = geno_omic_test_object,
                                     method = "pls",
                                     selected = selected
                                     #n_comp = max_features
                                     )

      test_scores <- pls_out$test_scores
    }

    if(any(gam_mod%in%c("rf_gam", "rf_gam_select"))){
      #scaled <- scale_omic_data(geno_omic_object[, selected[1:max_features], drop =FALSE])
      scaled <- scale_omic_data(geno_omic_object[, selected, drop =FALSE])
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


###################################



