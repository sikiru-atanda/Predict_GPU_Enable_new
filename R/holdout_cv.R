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
    method = c("stratified",
               "unstratified"),
    message = TRUE,
    ...
){

  msg <- sprintf("==================================================\n")

  test_Res <- vector(mode = "list", length = replication)

  #if(length(y) < 2) {stop("y must be greater than 1")}

  y <- pheno_data[, response]

  ID_GIDs = as.character(unique(pheno_data[, gen_name]))

  if (length(y)>length(ID_GIDs)){warning(paste(msg,'You have more than one environment. \n\tThis cross-validation method works best with one environment.'),
                                         call. = FALSE)}

  if(replication>1){warning(paste(msg,'You are using more than one replication. This might takes time.'),
                            call. = FALSE)}


  ### Initialize step to create groups to sample from using quantile for
  ## Efficient sampling

  ## If Binomial
  if(length(unique(y))==2){

    Grp_Sample = 2

  } else {

    if(length(unique(y))>2)

      Grp_Sample = 5

  }

  if (method == "stratified"){
    ##  create sample quantiles within the data  with probability[0, 1]
    if(is.numeric(y)) {
      y <- cut(y,
               unique(stats::quantile(y, probs = seq(0, 1, length = Grp_Sample))),
               include.lowest = TRUE)
      # } else {
      #   if(is.character(y)){
      #   ytab <- table(y)
      #   if(any(ytab == 0)) {
      #     warning(paste("Some classes have no records (",
      #                   paste(names(ytab)[ytab  == 0], sep = "", collapse = ", "),
      #                   ") and these will be ignored"))
      #     y <- factor(as.character(y))
      #   }
      #   if(any(ytab == 1)) {
      #     warning(paste("Some classes have a single record (",
      #                   paste(names(ytab)[ytab  == 1], sep = "", collapse = ", "),
      #                   ") and these will be selected for the sample"))
      #    }
      #
      #  }

    }

    # Rand_sample <- function(DT,  test_size) {
    #     TST_size <- ceiling(nrow(DT)*test_size)
    #     test_Res <- sample(DT$ID, size = TST_size, replace = FALSE)
    #     return(test_Res)
    #
    # }

    if(!is.null(random_state) & is.numeric(random_state)){
      set.seed(random_state)
    }

    #DT = data.frame(Response = y, ID = seq_along(y), stringsAsFactors = F)

    for (r in 1:replication) {

      strata <- split(1:length(y), y)
      TST_sample <- sort(as.vector(unlist(sapply(strata, function(x) sample(x, ceiling(length(x) * test_size), replace = FALSE)))))

      # TST_sample <- plyr::dlply(DT,"Response", Rand_sample, test_size)
      #
      # TST_sample <- sort(as.vector(unlist(TST_sample)))

      #test_Res[[r]] <- as.numeric(unlist(TST_sample))

      test_Res[[r]] <- TST_sample
    }

  } else {

    if (method == "unstratified"){
      warning(paste(msg,'stratified is a better sampling option, consider using it.'),
              call. = FALSE)

      for (r in 1:replication) {
        TST_sample <- sample(length(y), floor(length(y)*test_size),
                             replace = FALSE)

        test_Res[[r]] <- TST_sample

      }

    }
  }

  names(test_Res) <- paste("Rep_test_Sample",
                           gsub(" ", "0", format(seq_along(1:length(test_Res)))),
                           sep = "")


  return(test_Res)
}
