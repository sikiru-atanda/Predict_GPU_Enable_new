#' Title
#'
#' @param pheno_data
#' @param response
#' @param gen_name
#' @param nfolds
#' @param random_state
#' @param sampling_method
#' @param replication
#'
#' @return
#' @export
#'
#' @examples
kfolds_stratified_un <- function(
    pheno_data,
    response,
    gen_name,
    nfolds = 5,
    random_state = NULL,
    sampling_method = c("stratified", "unstratified"),
    replication = 1,
    ...
) {

  msg <- "\n==================================================\n"

  if (!is.null(random_state) && is.numeric(random_state)) {
    set.seed(random_state)
  }

  sampling_method <- match.arg(sampling_method)
  y <- as.double(pheno_data[[response]])

  if (nfolds < 2) stop(paste(msg, "Number of nfolds should be greater than one for k-fold Cross_validation"), call. = FALSE)
  if (nfolds > length(y)) stop(paste(msg, "Y variable must be greater than number of nfolds"), call. = FALSE)

  if (length(unique(pheno_data[[gen_name]])) < length(y)) {
    warning(paste(msg, "You have more than one environment. This Cross_validation method works best with one environment."))
  }

  Rep_FoldCV <- vector("list", replication)

  for (r in 1:replication) {
    folds <- vector("integer", length(y))

    if (sampling_method == "stratified") {
      if (!is.numeric(y) || length(unique(y)) <= 2) {
        # Handle binary and categorical variables
        levels_y <- levels(factor(y))
        for (lvl in levels_y) {
          idx <- which(y == lvl)
          fold_assignments <- sample(rep(1:nfolds, length.out = length(idx)))
          folds[idx] <- fold_assignments
        }
      } else {
        # Handle continuous variables using quantile-based stratification
        num_bins <- min(max(floor(length(y) / nfolds), 2), 5)
        cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)

        # Ensure cut_points are unique
        while (any(duplicated(cut_points))) {
          num_bins <- num_bins - 1
          if (num_bins < 2) {
            warning(paste(msg, "Unable to create unique bins for stratified sampling, switching to unstratified sampling."), call. = FALSE)
            folds <- sample(rep(1:nfolds, length.out = length(y)))
            Rep_FoldCV[[r]] <- folds
            names(Rep_FoldCV)[r] <- sprintf("Rep%d_%dFold_CV", r, nfolds)
            next
          }
          cut_points <- quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
        }

        y_binned <- cut(y, breaks = cut_points, include.lowest = TRUE)
        levels_y <- levels(y_binned)
        for (lvl in levels_y) {
          idx <- which(y_binned == lvl)
          fold_assignments <- sample(rep(1:nfolds, length.out = length(idx)))
          folds[idx] <- fold_assignments
        }
      }
    } else {
      # Unstratified sampling
      folds <- sample(rep(1:nfolds, length.out = length(y)))
    }

    Rep_FoldCV[[r]] <- folds
    names(Rep_FoldCV)[r] <- sprintf("Rep%d_%dFold_CV", r, nfolds)
  }

  return(Rep_FoldCV)
}



# test_folds = nfolds_stratified_un(pheno_data = pheno,
#                                   response = "Yield",
#                                   gen_name = "GID",
#                                   sampling_method = "stratified",
#                                   replication = 1,
#                                   nFolds = 5)
#
# nfolds <- as.double(strsplit(strsplit(names(test_folds)[1], c("_"))[[1]][2],"")[[1]][1])
