gp_expand_param_grid <- function(param_grid) {
  if (is.null(param_grid) || !length(param_grid)) {
    stop("param_grid must be a non-empty named list.", call. = FALSE)
  }
  as.data.frame(
    do.call(expand.grid, c(param_grid, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
  )
}

gp_default_tuning_metric <- function(response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  switch(
    fam,
    gaussian = "root_mean_squared_error",
    binary = "log_loss",
    ordinal = "log_loss",
    multiclass = "log_loss",
    "root_mean_squared_error"
  )
}

gp_metric_minimize <- function(metric) {
  tolower(metric) %in% c(
    "mean_squared_error", "root_mean_squared_error", "relative_squared_error",
    "mean_absolute_error", "mean_absolute_percent_error", "bias",
    "log_loss", "brier_score", "ece", "mean_absolute_error_class"
  )
}

gp_tuning_folds <- function(y, response_family = "gaussian", nfolds = 5L, random_state = 123L) {
  fam <- gp_resolve_response_family(response_family, y = y)
  gp_set_seed(as.integer(random_state %||% 123L))
  if (identical(fam, "gaussian")) {
    caret::createFolds(y, k = as.integer(nfolds), list = TRUE, returnTrain = FALSE)
  } else {
    caret::createFolds(as.factor(y), k = as.integer(nfolds), list = TRUE, returnTrain = FALSE)
  }
}

gp_eval_metric_from_prediction <- function(y_true, pred_obj, metric, response_family) {
  pred_val <- pred_obj$pred
  prob_val <- pred_obj$prob
  metric <- tolower(metric)
  fam <- gp_resolve_response_family(response_family, y = y_true)

  if (is.null(prob_val) && metric %in% c("log_loss", "brier_score", "ece")) {
    pred_num <- suppressWarnings(as.numeric(pred_val))
    n_obs <- length(y_true)
    n_cls <- length(unique(as.character(stats::na.omit(y_true))))
    if (length(pred_num) > n_obs && n_obs > 0L && n_cls > 1L && length(pred_num) %% n_obs == 0L) {
      inferred_cls <- length(pred_num) / n_obs
      if (inferred_cls == n_cls) {
        prob_val <- matrix(pred_num, nrow = n_obs, ncol = inferred_cls, byrow = FALSE)
        colnames(prob_val) <- sort(unique(as.character(stats::na.omit(y_true))))
      }
    }
  }

  if (metric %in% c("log_loss", "brier_score") && !is.null(prob_val)) {
    if (is.matrix(prob_val) || is.data.frame(prob_val)) {
      return(evaluation_metrics(y_true, as.matrix(prob_val), metric, response_family = fam))
    }
    return(evaluation_metrics(y_true, as.numeric(prob_val), metric, response_family = fam))
  }
  evaluation_metrics(y_true, pred_val, metric, response_family = fam)
}

gp_grid_tune_cv <- function(y,
                            response_family = "gaussian",
                            param_grid,
                            predict_fun,
                            batch_predict_fun = NULL,
                            tune_metric = NULL,
                            nfolds = 5L,
                            random_state = 123L) {
  fam <- gp_resolve_response_family(response_family, y = y)
  metric <- tune_metric %||% gp_default_tuning_metric(fam)
  combos <- gp_expand_param_grid(param_grid)
  folds <- gp_tuning_folds(y = y, response_family = fam, nfolds = nfolds, random_state = random_state)
  minimize <- gp_metric_minimize(metric)

  score_prediction <- function(pred, tst) {
    tryCatch(
      gp_eval_metric_from_prediction(
        y_true = y[tst],
        pred_obj = pred,
        metric = metric,
        response_family = fam
      ),
      error = function(e) {
        pred_len <- length(pred$pred %||% NULL)
        prob_dim <- if (!is.null(pred$prob)) paste(dim(as.matrix(pred$prob)), collapse = "x") else "NULL"
        stop(
          sprintf(
            "Tuning metric evaluation failed [family=%s metric=%s n_obs=%s pred_len=%s prob_dim=%s]: %s",
            fam,
            metric,
            length(y[tst]),
            pred_len,
            prob_dim,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  }

  batch_scores <- NULL
  if (!is.null(batch_predict_fun)) {
    batch_scores <- tryCatch({
      specs <- vector("list", nrow(combos) * length(folds))
      k <- 0L
      for (i in seq_len(nrow(combos))) {
        params <- as.list(combos[i, , drop = FALSE])
        for (j in seq_along(folds)) {
          k <- k + 1L
          specs[[k]] <- list(
            param_index = i,
            fold_index = j,
            tst = folds[[j]],
            params = params
          )
        }
      }
      preds <- batch_predict_fun(specs)
      if (length(preds) != length(specs)) {
        stop(
          "Batched tuning returned ", length(preds), " predictions for ",
          length(specs), " requested jobs.",
          call. = FALSE
        )
      }
      vals <- matrix(NA_real_, nrow = nrow(combos), ncol = length(folds))
      for (k in seq_along(specs)) {
        vals[specs[[k]]$param_index, specs[[k]]$fold_index] <- score_prediction(
          pred = preds[[k]],
          tst = specs[[k]]$tst
        )
      }
      rowMeans(vals, na.rm = TRUE)
    }, error = function(e) {
      warning(
        "Batched tuning failed; falling back to serial tuning: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    })
  }

  score_one <- function(params_row) {
    params <- as.list(params_row)
    vals <- vapply(folds, function(tst) {
      pred <- predict_fun(tst = tst, params = params)
      score_prediction(pred = pred, tst = tst)
    }, numeric(1))
    mean(vals, na.rm = TRUE)
  }

  scores <- batch_scores %||% vapply(seq_len(nrow(combos)), function(i) score_one(combos[i, , drop = FALSE]), numeric(1))
  best_idx <- if (minimize) which.min(scores) else which.max(scores)

  list(
    metric = metric,
    minimize = minimize,
    best_index = best_idx,
    best_score = scores[[best_idx]],
    best_params = as.list(combos[best_idx, , drop = FALSE]),
    all_results = cbind(combos, score = scores)
  )
}
