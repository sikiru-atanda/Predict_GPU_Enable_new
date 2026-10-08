#' K-fold cross-validation split assignments
#'
#' Legacy wrapper around the canonical k-fold splitter.
#'
#' @param nFolds Number of folds.
#' @param random_state Optional random seed.
#' @param replication Number of repeated split replications.
#' @param pheno_data Phenotype data frame.
#' @param response Response column.
#' @param gen_name Genotype ID column.
#' @param method Either `"stratified"` or `"unstratified"`.
#'
#' @return A list of integer fold-assignment vectors, one per replication.
#' @export
CV_nfolds <- function(
    pheno_data = NULL,
    response = NULL,
    gen_name = NULL,
    nFolds = 5,
    random_state = NULL,
    method = c("stratified", "unstratified"),
    replication = 1) {
  kfolds_stratified_un(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    nfolds = nFolds,
    random_state = random_state,
    sampling_method = match.arg(method),
    replication = replication
  )
}
