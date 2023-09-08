

#' Title
#'
#' Though it is expected the user provide clean data.
#' To keep everything tidy for model fitting the function perform the
#' following objectives:
#' 1. Check if the object_geno (geno/omic data) file contain records
#' 2. Check if the object_geno is a matrix and if not we fix it for the user
#' 3. Check the row and column names
#' 4. Check that the marker/column are all in numeric and not character
#' 5. Check if the allele dosages is not in 0, 1, 2 and it not, we fix it
#' 6. Check for monomorphic marker and remove it.
#' 7. Check for NA and remove it.
#' 8. It check for MAF and Heterozygosity '
#' 9. Here metadata is not generated because it for check not actually for QC.
#'    It is actually expected to be a clean data. If you need that please refer to marker_qc_recode function
#' 10. The output is clean snp/marker data with NA,monomorphic, maf, heterozgous
#'     removed user defined some threshould value. If not default values are used.
#'    The output is assign attribute "pass" and declared class 'geno_data'
#'
#' @param object_geno snp/marker data. NA is allowed
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
geno_precheck <- function(object_geno = NULL,
                          message = TRUE,
                          maf = 0.05,
                          freqHet = 0.2,
                          ...) {


  msg <- sprintf("==================================================\n")
  ## Check marker matrix is not null
  if (!is.null(object_geno)) {

    ### Check marker matrix is class matrix
    # if (inherits(object_geno, what = 'matrix')) {
    #   object_geno <- as.matrix(object_geno)
    # } else {
    #   stop("The object_geno object_geno is not a class matrix")
    # }

    if(class(object_geno)[[1]]!="matrix"){
      if(isTRUE(message)){
      message(paste(msg,"The data is not class matrix: we fix it." ))
      }

      object_geno <- as.matrix(object_geno)
    }

    # Check row and column names in object_geno.
    if(is.null(rownames(object_geno))){
      stop("Individual names not assigned to rows of \'object_geno'.")
    }

    if(is.null(colnames(object_geno))){
      stop("Marker names not assigned to columns of \'object_geno'.")
    }

    ## Check if allele dosage are in 0, 1, 2 format
    check_geno <- which(object_geno== -1)

    if(length(check_geno!=0)){
      if(isTRUE(message)){
      message(paste(msg,'The allele dosages is not in 0, 1, 2. we fix it for you'))
      }#c(-1, 0, 1) + 1
      #[1] 0 1 2
      object_geno <- apply(object_geno, 2, function(x) x+1)}

    # check that marker matrices did not have colmun with other variables
    # such as the genotypes name other than the alleles dosage
    check_geno <- which(object_geno!=1 & object_geno!=0 & object_geno!=2)

    if(length(check_geno!=0)){stop(print(paste(msg,'Marker data contains variables other than the allele dosages:\n \t \t \t 0, 1, 2')), call. = FALSE)}

    rm(check_geno)
    ### Check and remove monomorhpic markers
    if(message){
      message("\nMonomorphic check.  \n")
    }
    object_geno <- Remove_NA_Mono_SNP(object_geno)

    ## Check for Na and remove
    Na_col.omit <- which((colSums(is.na(object_geno))==0)==FALSE)

    if (length(Na_col.omit)!=0){
      #object_geno = object_geno[ , colSums(is.na(object_geno))==0]

      if(isTRUE(message)){
        message("A total of ", length(Na_col.omit),
                " SNP (s)/ marker (s) were removed from due to missing value")

      }

      object_geno = object_geno[ , -Na_col.omit]
    }

    ###  # Minor alele frequency
    #### MAF #############
    #compute p
    #Here the allele fre. of the second allele (denoted as phat) is calculated and if phat is less than 0.5 , phat = maf, if not maf = 1-phat
    phat=colMeans(object_geno, na.rm = T)/2
    MAF=ifelse(phat<0.5,phat,1-phat) ##

    if(isTRUE(message)){
    message("\nMinor allele frequency check (MAF).  \n")
    }

    if(length(which(MAF<maf))!=0){


    object_geno =object_geno[,-which(MAF<maf)]  ##0r (MAF ??? 0.01)

    if(maf==0.05){
      if(isTRUE(message)){
        message(paste(msg,"if higher or lower MAF threshold value is desired provide one"))
      }
    }

    if(isTRUE(message)){
      #message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',ncol(object_geno) - length(which(MAF<maf))))))
      message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',length(which(MAF<maf))))))
    }

    #dim(object_geno)
    } else {
      if(isTRUE(message)){

        message(paste(msg, paste("\tNo loci with MAF below",paste("threshold =", maf), sep = " ")))
      #cat("\tNo loci with minor allele frequency below", paste("threshold =", maf),"\n")
      }
    }


    #if (any(is.na(object_geno))) {stop(print(paste(msg,'Marker data contains Missing value')), call. = FALSE)}

    # Heterozigosity
    Het <- apply(object_geno, 2, function(x) sum( x== 1, na.rm=T)/nrow(object_geno))

    if(length(which(Het<=freqHet))!=0){

      N_int = ncol(object_geno)
    object_geno = object_geno[, which(Het<=freqHet)]

    if(isTRUE(message)){
      message(paste(msg,(paste(paste('Total loci above heterozygosity', paste("threshold", freqHet)), sep = ': \t', N_int - ncol(object_geno) ))))
    }

    }

    # Assign appropriate class.
    class(object_geno) <-c("matrix", "array", "geno_data")

    attr(object_geno, "cleared") <- "pass"

  } else {

    object_geno = NULL
  }




  return(object_geno)

}

