# Centring, scaling and folds as caret makes them (caret 7.0.1, GPL >= 2:
# preProcess(method = c("center", "scale")), its predict() method and
# createFolds()), without loading caret: 1.33 s and 61 namespaces on each
# worker's first ML task on the Windows test machine. The results are
# caret's exactly, so a seed keeps its folds and its predictions.
# gp_fast_center_scale_fit() differs in the last digits (colMeans(); the
# deviations of the raw columns), so the models that used caret use these.

# Means and standard deviations of the numeric columns of a matrix or data
# frame, NA removed. The deviations are taken after centring; an NA or zero
# deviation becomes 1, with caret's warning.
gp_center_scale_fit <- function(x, center = TRUE, scale = TRUE) {
  if (!inherits(x, "matrix") && !inherits(x, "data.frame")) {
    stop("Matrices or data frames are required for preprocessing", call. = FALSE)
  }
  if (is.null(colnames(x))) stop("`x` must have column names", call. = FALSE)
  types <- if (is.matrix(x)) {
    rep(class(x[1, 1]), ncol(x))
  } else {
    unlist(lapply(x, function(v) class(v)[1]))
  }
  cols <- colnames(x)[types %in% c("integer", "numeric", "double")]
  mean_value <- NULL
  std_value <- NULL
  if (isTRUE(center)) {
    mean_value <- apply(x[, cols, drop = FALSE], 2, mean, na.rm = TRUE)
    x[, cols] <- sweep(x[, cols, drop = FALSE], 2, mean_value, "-")
  }
  if (isTRUE(scale)) {
    std_value <- apply(x[, cols, drop = FALSE], 2, stats::sd, na.rm = TRUE)
    if (any(is.na(std_value))) {
      warning(paste("Std. deviations could not be computed for:",
                    paste(names(std_value)[which(is.na(std_value))], collapse = ", ")),
              immediate. = TRUE)
      std_value[which(is.na(std_value))] <- 1
    }
    if (any(std_value == 0)) {
      warning(paste("These variables have zero variances:",
                    paste(names(std_value)[which(std_value == 0)], collapse = ", ")),
              immediate. = TRUE)
      std_value[which(std_value == 0)] <- 1
    }
  }
  structure(list(mean = mean_value, std = std_value, columns = cols), class = "gp_center_scale")
}

# Applies a gp_center_scale_fit() scaler to `newdata`; a caret preProcess
# object (as kept by an earlier version) is applied the same way.
gp_center_scale_apply <- function(scaler, newdata) {
  caret_made <- inherits(scaler, "preProcess")
  center_cols <- if (caret_made) scaler$method$center else if (!is.null(scaler$mean)) scaler$columns
  scale_cols <- if (caret_made) scaler$method$scale else if (!is.null(scaler$std)) scaler$columns
  if (length(center_cols)) {
    newdata[, center_cols] <- sweep(newdata[, center_cols, drop = FALSE], 2, scaler$mean, "-")
  }
  if (length(scale_cols)) {
    newdata[, scale_cols] <- sweep(newdata[, scale_cols, drop = FALSE], 2, scaler$std, "/")
  }
  newdata
}

# caret::createFolds(): a numeric y is cut at up to 5 quantiles, and each of
# those groups (or each class) is spread over the k folds.
gp_create_folds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {
  if (is.numeric(y)) {
    cuts <- floor(length(y) / k)
    if (cuts < 2) cuts <- 2
    if (cuts > 5) cuts <- 5
    breaks <- unique(stats::quantile(y, probs = seq(0, 1, length.out = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    num_in_class <- table(y)
    fold_vector <- vector(mode = "integer", length(y))
    for (i in seq_along(num_in_class)) {
      min_reps <- num_in_class[i] %/% k
      if (min_reps > 0) {
        spares <- num_in_class[i] %% k
        seq_vector <- rep(1:k, min_reps)
        if (spares > 0) seq_vector <- c(seq_vector, sample(1:k, spares))
        fold_vector[which(y == names(num_in_class)[i])] <- sample(seq_vector)
      } else {
        fold_vector[which(y == names(num_in_class)[i])] <- sample(1:k, size = num_in_class[i])
      }
    }
  } else {
    fold_vector <- seq(along.with = y)
  }
  if (!list) return(fold_vector)
  out <- split(seq(along.with = y), fold_vector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along.with = out))), sep = "")
  if (returnTrain) out <- lapply(out, function(data, y) y[-data], y = seq(along.with = y))
  out
}
