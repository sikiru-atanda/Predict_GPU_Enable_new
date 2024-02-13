
#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param heter_groups
#' @param CV
#' @param random_state
#' @param replication
#' @param nFolds
#'
#' @return
#' @export
#'
#' @examples


CV1_CV2 <- function(
    pheno_data,
    gen_name,
    heter_groups,
    CV = NULL,
    nFolds = NULL,
    random_state = NULL,
    replication = 1){

  msg <- sprintf("==================================================\n")




  if (is.null(heter_groups)){stop(message(paste(msg,"Provide the a pointer (heter_groups) to the column contaning the environments")), call. = FALSE)}
  if(CV>2){stop(message(paste(msg,"CV must be 1 or 2")), call. = FALSE)}
  if(is.null(nFolds)){stop(message(paste(msg,"Provide value the number of desired folds")), call. = FALSE)}
  ## Order the pheno_data data by gen_name and by Environment
  pheno_data = pheno_data[order(pheno_data[, gen_name]), ]
  pheno_data = pheno_data[order(pheno_data[, heter_groups]), ]

  #nEnv <- length(unique(pheno_data[, heter_groups]))

  ID_GIDs = as.character(unique(pheno_data[, gen_name]))
  Envs_ID_GIDs = as.character(pheno_data[, gen_name])



  if(length(ID_GIDs)==length(Envs_ID_GIDs)){stop(message(paste(msg,'CV1 and CV2 works when number of environment is greater than 1')), call. = FALSE)}

  if(replication>1){warning(paste(paste(msg,'You request for', replication), 'replications this might takes some time to run all the replications.'),
                            call. = FALSE)}

  Rep_FoldCV = vector(mode = "list", length = replication)

  All_nFolds <- vector(mode = "integer",  length(pheno_data[, gen_name]))

  if(!is.null(random_state) & is.numeric(random_state)){
    set.seed(random_state)
  }

  if (CV == 1) {

    for (r in 1:replication) {

      mfold <- sample(1:nFolds, size = length(ID_GIDs), replace = TRUE)

      for (i in 1:length(pheno_data[, gen_name])) {

        All_nFolds[i] <- mfold[which(ID_GIDs == pheno_data[, gen_name][i])]

      }

      Rep_FoldCV[[r]] <- All_nFolds

      names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
                                          paste(nFolds, 'Fold', sep = ""), sep="_"),
                                    paste("CV", CV, sep=""), sep = "")

    }

  }

  if (CV == 2) {


    for (r in 1:replication) {

      for (i in ID_GIDs) {

        Env_GIDs = which(Envs_ID_GIDs == i)

        Env_GIDs_size = length(Env_GIDs)

        tmpFold <- sample(1:nFolds, size =  Env_GIDs_size, replace =  Env_GIDs_size > nFolds)

        All_nFolds[Env_GIDs] <- tmpFold

      }

      Rep_FoldCV[[r]] <- All_nFolds

      names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
                                          paste(nFolds, 'Fold', sep = ""), sep="_"),
                                    paste("CV", CV, sep=""), sep = "")

    }

  }
  return(Rep_FoldCV)
}
