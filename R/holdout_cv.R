#' Hold-out train/test split indices
#'
#' Legacy wrapper around the canonical hold-out splitter.
#'
#' @param test_size Fraction of rows assigned to the test set.
#' @param random_state Optional random seed.
#' @param replication Number of repeated split replications.
#' @param pheno_data Phenotype data frame.
#' @param gen_name Genotype ID column.
#' @param response Response column.
#' @param method Either `"stratified"` or `"unstratified"`.
#' @param ... Additional arguments forwarded to the underlying splitter;
#'   reserved for forward compatibility.
#' @param message Logical; emit splitter messages and warnings.
#'
#' @return A list of integer test-index vectors, one per replication.
#' @export
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
) {
  hold_out_stratified_and_un(
    pheno_data = pheno_data,
    gen_name = gen_name,
    response = response,
    test_size = test_size,
    random_state = random_state,
    replication = replication,
    sampling_method = match.arg(method),
    message = message,
    ...
  )
}
