

#' Title
#'
#' @param pheno_data
#' @param response
#' @param weights
#' @param ETA
#' @param bayes_para
#' @param ...
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
bayes_mod_execute <- function(pheno_data = NULL,
                              response = NULL,
                              weights = NULL,
                              ETA = NULL,
                              bayes_para = NULL,
                              verbose = FALSE,
                              core = NULL,
                              ...)
  {


  ###
  # bayes_para <- bayes_parameter_check(nIter = nIter,
  #                             burnIn = burnIn,
  #                             thin = thin)

#   # #### Initializing parallel
#   if(length(response)>1){
#     if (is.null(core)){
#       cl = parallel::detectCores()
#
#       if (cl> 4){
#         # Try in parallel
#         cl <- parallel::makeCluster(4)
#       } else{
#         cl <- parallel::makeCluster(2)
#       }
#
#     } else {
#       if(!is.null(core)){
#
#         cl <- parallel::makeCluster(core)
#       }
#     }
#
#     doParallel::registerDoParallel(cl)
#
#
# ### issue to address
#   Univariate <- foreach::foreach(trait = 1:length(response),
#                                  .errorhandling='pass') %dopar% {
#
#   current_date_time = as.character(Sys.time())
#
#   files_key = gsub(" ", "", current_date_time)
#
#   ### Create key to remove all reduant files from the wkdir
#   #files_key = strsplit(files_key, "\\.")[[1]][2]
#   files_key = strsplit(files_key, "\\.")[[1]][1]
#
#   files_key = gsub("-", "", files_key)
#
#   files_key= gsub(":", "_", files_key)
#   files_key = paste(files_key, response[trait], sep = "_")
#
#   #DateTime = paste0(gsub(":", "_", DateTime), trait)
#   ### Create key to remove all reduant files from the wkdir
#   # files_key = strsplit(files_key, "-")[[1]][1]
#   #
#   # files_key = gsub("-", "", files_key)
#
#
#
#
#   if(is.null(weights)){
#     fm <- BGLR::BGLR(
#       y=pheno_data[, response[trait]],
#       ETA = ETA,
#       nIter = bayes_para$nIter,
#       burnIn =  bayes_para$burnIn,
#       thin =  bayes_para$thin,
#       verbose = FALSE,
#       saveAt =files_key)
#
#   } else{
#
#     if(!is.null(weights)){
#       fm <- BGLR::BGLR(
#         y=pheno_data[, response[trait]],
#         ETA=ETA,
#         weights = weights,
#         nIter= bayes_para$nIter,
#         burnIn= bayes_para$burnIn,
#         thin = bayes_para$thin,
#         verbose = FALSE,
#         saveAt = files_key)
#
#     }
#
#   }
#
#
#   ## Get the name of all files stored by GBLR using the current name and time the analysis was performed
#   output_files_names = list.files(pattern=files_key)
#
#   output = list(model = fm, output_files_names = output_files_names )
#   #rm(output_files_names)
#
#   names(output) <- c("model", "output_files_names")
#
#   Univariate = output
#
#
#         }
#
#   names(Univariate) <- response
#
#   output <-  Univariate
#
#   rm(Univariate)
#
#   } else {

    #current_date_time = as.character(Sys.time())
  systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  systime <- gsub("[-: ]", "_", systime)


    if(is.null(weights)){
      fm <- BGLR::BGLR(
        y=pheno_data[, response],
        ETA = ETA,
        nIter = bayes_para[["nIter"]],
        burnIn =  bayes_para[["burnIn"]],
        thin =  bayes_para[["thin"]],
        verbose = FALSE,
        saveAt =systime)

    } else{

      if(!is.null(weights)){
        fm <- BGLR::BGLR(
          y=pheno_data[, response],
          ETA=ETA,
          weights = weights,
          nIter = bayes_para[["nIter"]],
          burnIn =  bayes_para[["burnIn"]],
          thin =  bayes_para[["thin"]],
          verbose = FALSE,
          saveAt = systime)

      }

    }


    ## Get the name of all files stored by BGLR using the current name and time the analysis was performed
    output_files_names = list.files(pattern=systime)

    output = list(model = fm, output_files_names = output_files_names )
    #rm(output_files_names)

    names(output) <- c("model", "output_files_names")


  #}
  # ### remove the generated output files from the working directory
  # unlink(output_files_names)

  return(output)

}




