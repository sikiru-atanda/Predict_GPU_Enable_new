#' Title
#'
#' @param test_size
#' @param random_state
#' @param replication
#' @param pheno_data
#' @param gen_name
#' @param response
#' @param method
#' @param ...
#' @param message
#'
#' @return
#' @export
#'
#' @examples
#'
#'
train_test_split <- function (
    pheno_data = NULL,
    gen_name = NULL,
    response = NULL,
    test_size = 0.2,
    random_state = NULL,
    replication = 1,
    method = c("stratified", "unstratified"),
    message = TRUE,
    ...
){

  msg <- "\n==================================================\n"

  test_Res <- vector(mode = "list", length = replication)

  y <- as.double(pheno_data[[response]])
  ID_GIDs = as.character(unique(pheno_data[[gen_name]]))

  if (length(y) > length(ID_GIDs) && isTRUE(message)) {
    warning(paste(msg, 'You have more than one environment. This cross-validation method works best with one environment.'), call. = FALSE)
  }

  if (replication > 1 && isTRUE(message)) {
    warning(paste(msg, 'You are using more than one replication. This might take time.'), call. = FALSE)
  }

  Grp_Sample <- if (length(unique(y)) == 2) 2 else 5

  selected_method <- match.arg(method)

  if (selected_method == "stratified") {
    if (is.numeric(y)) {
      num_bins <- Grp_Sample
      cut_points <- unique(stats::quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE))

      while (any(duplicated(cut_points))) {
        num_bins <- num_bins - 1
        if (num_bins < 2) {
          warning(paste(msg, "Unable to create unique bins for stratified sampling, switching to unstratified sampling"), call. = FALSE)
          selected_method <- "unstratified"
          break
        }
        cut_points <- unique(stats::quantile(y, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE))
      }

      if (selected_method == "stratified") {
        y <- cut(y, cut_points, include.lowest = TRUE)
      }
    }

    if (!is.null(random_state) && is.numeric(random_state)) {
      set.seed(random_state)
    }

    for (r in 1:replication) {
      strata <- split(seq_len(length(y)), y)
      TST_sample <- sort(unlist(lapply(strata, function(x) sample(x, ceiling(length(x) * test_size), replace = FALSE))))
      test_Res[[r]] <- TST_sample
    }
  }

  if (selected_method == "unstratified") {
    if (selected_method == "unstratified" && isTRUE(message)) {
      warning(paste(msg, 'Stratified is a better sampling option, consider using it.'), call. = FALSE)
    }

    if (!is.null(random_state) && is.numeric(random_state)) {
      set.seed(random_state)
    }

    for (r in 1:replication) {
      TST_sample <- sample(seq_len(length(y)), floor(length(y) * test_size), replace = FALSE)
      test_Res[[r]] <- TST_sample
    }
  }

  names(test_Res) <- paste("Rep_test_Sample", gsub(" ", "0", format(seq_along(test_Res))), sep = "")

  return(test_Res)
}

