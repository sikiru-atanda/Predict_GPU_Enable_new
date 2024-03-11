#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param response
#' @param test_size
#' @param random_state
#' @param replication
#' @param sampling_method
#' @param message
#'
#' @return
#' @export
#'
#' @examples
hold_out_stratified_and_un <- function (
                                      pheno_data,
                                      gen_name,
                                      response,
                                      test_size = 0.2,
                                      random_state = 123,
                                      replication = 1,
                                      sampling_method = c("stratified", "unstratified"),
                                      message = TRUE
                                  ) {

  #msg <- sprintf("==================================================\n")

  if (!is.null(random_state) && is.numeric(random_state)) {
    set.seed(random_state)
  }

  selected_sampling_method <- match.arg(sampling_method)

  y <- pheno_data[[response]]
  ID_GIDs <- as.character(unique(pheno_data[[gen_name]]))

  if (length(y) > length(ID_GIDs) && isTRUE(message)) {
    warning("You have more than one environment. This cross-validation method works best with one environment.")
  }

  if (replication > 1 && isTRUE(message)) {
    warning("You are using more than one replication. This might take time.")
  }

  grp_sample <- if (length(unique(y)) == 2) 2 else 5

  test_res <- vector(mode = "list", length = replication)
  names(test_res) <- sprintf("Rep_test_Sample%03d", seq_len(replication))

  for (r in seq_len(replication)) {
    if (selected_sampling_method == "stratified") {
      if (is.factor(y) || length(unique(y)) <= grp_sample) {
        strata <- split(seq_len(length(y)), y)
      } else {
        cut_points <- unique(quantile(y, probs = seq(0, 1, length.out = grp_sample + 1), na.rm = TRUE))
        strata <- split(seq_len(length(y)), cut(y, cut_points, include.lowest = TRUE))
      }
      tst_sample <- unlist(lapply(strata, function(x) sample(x, ceiling(length(x) * test_size), replace = FALSE)))
    } else {
      tst_sample <- sample(seq_len(length(y)), floor(length(y) * test_size), replace = FALSE)
    }

    test_res[[r]] <- sort(tst_sample)
  }

  if (selected_sampling_method == "unstratified" && isTRUE(message)) {
    warning("Stratified is a better sampling option, consider using it.")
  }

  return(test_res)
}


