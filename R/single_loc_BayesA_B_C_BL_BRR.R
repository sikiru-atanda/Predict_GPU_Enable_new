

#' Title
#'
#' @param object
#' @param response
#' @param weights
#' @param ETA
#' @param bayes_para
#' @param verbose
#' @param saveAt
#'
#' @return
#' @export
#'
#' @examples
 M_matrix_bayes_mod_single_loc <- function(object = NULL,
                                         response = NULL,
                                         weights = NULL,
                                         ETA = NULL,
                                         bayes_para = NULL,
                                         verbose = FALSE,
                                         ...
                                         ){

   ###
# bayes_para <- bayes_parameter_check(nIter = nIter,
#                             burnIn = burnIn,
#                             thin = thin)

current_date_time = as.character(Sys.time())

files_key = gsub(" ", "", current_date_time)

### Create key to remove all reduant files from the wkdir
files_key = strsplit(files_key, "\\.")[[1]][2]

### Create key to remove all reduant files from the wkdir
# files_key = strsplit(files_key, "-")[[1]][1]
#
# files_key = gsub("-", "", files_key)


if(is.null(weights)){
  fm <- BGLR::BGLR(
    y=object[, response],
    ETA = ETA$ETA,
    nIter = bayes_para$nIter,
    burnIn =  bayes_para$burnIn,
    thin =  bayes_para$thin,
    verbose = FALSE,
    saveAt =files_key)

} else{

  if(!is.null(weights)){
    fm <- BGLR::BGLR(
      y=object[, response],
      ETA=ETA$ETA,
      weights = weights,
      nIter= bayes_para$nIter,
      burnIn= bayes_para$burnIn,
      thin = bayes_para$thin,
      verbose = FALSE,
      saveAt = files_key)

  }

}


## Get the name of all files stored by GBLR using the current name and time the analysis was performed
output_files_names = list.files(pattern=files_key)

output = list(model = fm, output_files_names = output_files_names )

names(output) <- c("model", "output_files_names")

return(output)

}

