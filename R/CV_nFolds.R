#' Title
#'
#' @param nFolds
#' @param random_state
#' @param replication
#' @param pheno_data
#' @param response
#' @param gen_name
#' @param method
#'
#' @return
#' @export
#'
#' @examples
CV_nfolds <- function(
    pheno_data = NULL,
    response = NULL,
    gen_name =NULL,
    nFolds = 5,
    random_state = NULL,
    method = c("stratified",
               "unstratified" ),
    replication = 1) {

  msg <- sprintf("==================================================\n")

  #if(length(pheno_data[, response]) < 2) stop(print(paste(msg, "y must be greater than 1")), call. = FALSE)

  y <- pheno_data[, response]

  if (nFolds ==1) stop(print(paste(msg, 'Number of nFolds should be greater than one for k-fold Cross_validation')), call. = FALSE)

  if(nFolds> length(y)) stop(print(paste(msg, "Y variable must be greater than number of nFolds")), call. = FALSE)

  ID_GIDs = as.character(unique(pheno_data[, gen_name]))

  if (length(y)>length(ID_GIDs)){warning(paste(msg,'You have more than one environment. This Cross_validation method works best with one environment.'),
                                         call. = FALSE)}



  if(!is.null(random_state) & is.numeric(random_state)){
    set.seed(random_state)
  }


  Rep_FoldCV <- vector(mode = "list", length=replication)

  # if(!is.null(list) | list==TRUE) {
  # Rep_FoldCV <- vector(mode = "list", length=replication)
  # }
  #
  # if(is.null(list) | list==FALSE) {
  #   Rep_FoldCV <- matrix()
  # }

  msg <- sprintf("==================================================\n")

  if(replication>1){warning(paste(paste(msg,'You request for', replication), 'replications this might takes some time to run all the replications.'),
                            call. = FALSE)}

  if (method == "unstratified"){
    warning(paste(msg,'stratified is a better sampling option, consider using it.'),
            call. = FALSE)}

  for(r in 1:replication) {

    if (method == "stratified"){
      ### Initialize Stratified Sampling

      if(is.numeric(y)) {
        ## Stratified the data into groups (Grp_Sample ) using quantiles
        ## Sample within those Grp_Samples.
        ## The number of Grp_Samples will depend on the
        ## ratio of the number of nFolds to the sample size.
        ##
        ## If the size of the data is too small it might be worthwhile to do
        #  unstratified

        Grp_Sample <- floor(length(y)/nFolds)
        if(Grp_Sample < 2) Grp_Sample <- 2
        if(Grp_Sample > 5) Grp_Sample <- 5
        y <- cut(y,
                 unique(stats::quantile(y, probs = seq(0, 1, length = Grp_Sample))),
                 include.lowest = TRUE)
      }

      #if(nFolds < length(y)) {
      ## reset levels so that the possible levels and
      ## the levels in the vector are the same
      y <- factor(as.character(y))
      Grp_size <- table(y)
      All_nFolds <- vector(mode = "integer", length = length(y))

      # if(!is.null(random_state) & is.numeric(random_state)){
      #   set.seed(random_state)
      # } else{
      #   set.seed(NULL)
      # }
      ## For each group/sample, allocation of the nFolds should be balance and
      ## randomized
      for(i in 1:length(Grp_size)) {
        ##how many should be on each sample/group
        actual_sampling_size_per_grp <- Grp_size[i] %/% nFolds
        ## In a situation where there is a reminder
        if(actual_sampling_size_per_grp != 0) {
          addin_sampling_size_per_grp <- Grp_size[i] %% nFolds
          actual_seqnFolds <- rep(1:nFolds, actual_sampling_size_per_grp)
          ## In a situation where there is a reminder add both actual and
          ## reminder so that actual_seqnFolds == Grp_size[i]
          if(addin_sampling_size_per_grp != 0) {
            addin_seqnFolds <- sample(1:nFolds, addin_sampling_size_per_grp)
            actual_seqnFolds <- c(actual_seqnFolds, addin_seqnFolds)

          }
          ## Randomize the created integers for all the nFolds and
          # assign for each sample/group
          All_nFolds[which(y == names(Grp_size)[i])] <- sample(actual_seqnFolds)
        } else {
          ## Where there are few records in the sample/group than unique nFolds.
          # samples are randomly assign into class/group
          All_nFolds[which(y == names(Grp_size)[i])] <- sample(1:nFolds, size = Grp_size[i])
        }
      }

      # } else {
      #
      #   All_nFolds <- seq_along(y)
      #
      # }

      # if(list==TRUE) {
      #   Res_FoldCV <- split(seq(along = y), All_nFolds)
      #   names(Res_FoldCV) <- paste("Fold", gsub(" ", "0", format(seq(along = Res_FoldCV))),
      #                       sep = "")
      #
      # } else{
      #
      #   if(list==FALSE){
      #   Res_FoldCV <- All_nFolds
      #
      #   }

      #}

      #Rep_FoldCV[[r]] = Res_FoldCV

      Rep_FoldCV[[r]] = All_nFolds

      names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
                                          paste(nFolds, 'Fold', sep = ""), sep="_"),
                                    "CV", sep = "")

    } else {

      if (method == "unstratified"){
        # warning(paste(msg,'stratified is a better sampling option, consider using it.'),
        #         call. = FALSE)


        folds = sample(1:nFolds, size = length(y), replace = T)

        Rep_FoldCV[[r]] <-  folds

        names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
                                            paste(nFolds, 'Fold', sep = ""), sep="_"),
                                      "CV", sep = "")
      }

    }


  }

  return(Rep_FoldCV)
}
