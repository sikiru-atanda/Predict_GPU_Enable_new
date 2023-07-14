

#' Title
#'
#' @param object
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
geno_precheck <- function(object = NULL,
                          message = TRUE,
                          ...) {


  msg <- sprintf("==================================================\n")
  ## Check marker matrix is not null
  if (!is.null(object)) {

    ### Check marker matrix is class matrix
    # if (inherits(object, what = 'matrix')) {
    #   object <- as.matrix(object)
    # } else {
    #   stop("The object object is not a class matrix")
    # }

    if(class(object)[[1]]!="matrix"){
      if(message){
      message(paste(msg,"The data is not class matrix: we fix it." ))
      }

      object <- as.matrix(object)
    }

    # Check row and column names in object.
    if(is.null(rownames(object))){
      stop("Individual names not assigned to rows of \'object'.")
    }

    if(is.null(colnames(object))){
      stop("Marker names not assigned to columns of \'object'.")
    }

    ## Check if allele dosage are in 0, 1, 2 format
    check_geno <- which(object== -1)

    if(length(check_geno!=0)){
      if(message){
      message(paste(msg,'The allele dosages is not in 0, 1, 2. we fix it for you'))
      }#c(-1, 0, 1) + 1
      #[1] 0 1 2
      object <- apply(object, 2, function(x) x+1)}

    # check that marker matrices did not have colmun with other variables
    # such as the genotypes name other than the alleles dosage
    check_geno <- which(object!=1 & object!=0 & object!=2)

    if(length(check_geno!=0)){stop(print(paste(msg,'Marker data contains variables other than the allele dosages:\n \t \t \t 0, 1, 2')), call. = FALSE)}

    rm(check_geno)
    ### Check and remove monomorhpic markers
    object <- Remove_NA_Mono_SNP(object)

    ## Check for Na and remove
    Na_col.omit <- which((colSums(is.na(object))==0)==FALSE)

    if (length(Na_col.omit)!=0){
      #object = object[ , colSums(is.na(object))==0]

      if (message){
        message("A total of ", length(Na_col.omit),
                " SNP (s)/ marker (s) were removed from due to missing value")

      }

      object = object[ , -Na_col.omit]
    }

    #if (any(is.na(object))) {stop(print(paste(msg,'Marker data contains Missing value')), call. = FALSE)}

    # Assign appropriate class.
    class(object) <-c("matrix", "array", "geno_data")

    attr(object, "cleared") <- "pass"

  } else {

    object = NULL
  }




  return(object)

}

