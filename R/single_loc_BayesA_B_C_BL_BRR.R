

#' Fit a Bayesian marker-effect model (BayesA / B / C / BL / BRR) at a single location
#'
#' Wraps a [BGLR::BGLR] call for the marker (single-location) Bayesian models,
#' using the pre-compiled ETA list, MCMC controls in `bayes_para` and writing
#' the run's auxiliary files under a timestamped working directory.
#'
#' @param object Cleaned phenotype data frame containing the response column.
#' @param response Name of the response (trait) column in `object`.
#' @param weights Optional positive Stage 2 observation precisions, supplied as
#'   a numeric vector, one-column table, or column name in `object`. These are
#'   converted to `sqrt(weights)` before calling BGLR because BGLR defines
#'   residual variance as inverse squared native weight.
#' @param ETA Pre-compiled ETA list (from [ETA_compiler_bayes]) describing the
#'   fixed and random terms.
#' @param bayes_para Named list of BGLR MCMC controls (`nIter`, `burnIn`,
#'   `thin`) from [bayes_parameter_check].
#' @param verbose Logical; passed through to BGLR for its progress output.
#' @param core Integer; reserved compute / parallel hint.
#' @param ... Additional arguments forwarded to BGLR.
#'
#' @return The fitted BGLR model object.
#' @export
#'
#' @examples
# @importFrom foreach %dopar%
 M_matrix_bayes_mod_single_loc <- function(object = NULL,
                                         response = NULL,
                                         weights = NULL,
                                         ETA = NULL,
                                         bayes_para = NULL,
                                         verbose = FALSE,
                                         core,
                                         ...
                                         ){

   ###
# bayes_para <- bayes_parameter_check(nIter = nIter,
#                             burnIn = burnIn,
#                             thin = thin)

current_date_time = as.character(Sys.time())

files_key = gsub(" ", "", current_date_time)

### Create key to remove all reduant files from the wkdir
#files_key = strsplit(files_key, "\\.")[[1]][2]
files_key = strsplit(files_key, "\\.")[[1]][1]

files_key = gsub("-", "", files_key)

files_key= gsub(":", "_", files_key)

#DateTime = paste0(gsub(":", "_", DateTime), trait)
### Create key to remove all reduant files from the wkdir
# files_key = strsplit(files_key, "-")[[1]][1]
#
# files_key = gsub("-", "", files_key)

# #### Initializing parallel
# if(length(response)>1){
# if (is.null(core)){
#   cl = parallel::detectCores()
#
#   if (cl> 4){
#     # Try in parallel
#     cl <- parallel::makeCluster(4)
#   } else{
#     cl <- parallel::makeCluster(2)
#   }
#
# } else {
#   if(!is.null(core)){
#
#     cl <- parallel::makeCluster(core)
#   }
# }
#
# doParallel::registerDoParallel(cl)
#
# }


  bglr_weights <- gp_bglr_weights_from_precision(
    weights = weights,
    pheno_data = object,
    response = response,
    context = "BGLR marker-model Stage 2 observation weights"
  )

if(is.null(bglr_weights)){
  fm <- BGLR::BGLR(
    y=object[, response],
    ETA = ETA,
    nIter = bayes_para$nIter,
    burnIn =  bayes_para$burnIn,
    thin =  bayes_para$thin,
    verbose = FALSE,
    saveAt =files_key)

} else{

  if(!is.null(bglr_weights)){
    fm <- BGLR::BGLR(
      y=object[, response],
      ETA=ETA,
      weights = bglr_weights,
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
#rm(output_files_names)

names(output) <- c("model", "output_files_names")

# ### remove the generated output files from the working directory
# unlink(output_files_names)

return(output)

}

