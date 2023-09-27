

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
#' 7. Check for NA and remove it as defined by the user in the call_rate.
#' 8. Impute missing value if user defined it
#' 9. It check for MAF and Heterozygosity '
#' 10. metadata is generated for qc and filtering.
#'    Thought it is actually expected to be a clean data.
#' 11. The output is clean snp/marker data with NA,monomorphic, maf, heterozgous
#'     removed user defined some threshould value. If not default values are used.
#'    The output is assign attribute "pass" and declared class 'geno_data'
#'
#' @param object_geno snp/marker data. NA is allowed
#' @param maf minor allele frequency user defined It must be in percentage eg 0.05
#' @param freqHet heterozygosity accepted by the user. It must be in percentage eg 0.2
#' @param call_rate missing value defined by the user. It must be in percentage eg 0.4
#' @param impute If user want to impute missing value. Logical TRUE or FALSE
#' @param qc_filtering if user want to do quality control. Logical TRUE or FALSE
#' @param ...
#' @param message If user want to print out message
#'
#' @return
#' @export
#'
#' @examples
geno_precheck <- function(object_geno = NULL,
                          maf = 0.05,
                          freqHet = 0.2,
                          call_rate = 0.5,
                          impute = TRUE,
                          qc_filtering = TRUE,
                          message = TRUE,
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

      stop(print(paste(msg, "Individual names not assigned to rows of \'object_geno'.")), call. = FALSE)

    }

    if(is.null(colnames(object_geno))){
      stop(print(paste(msg, "Marker names not assigned to columns of \'object_geno'.")), call. = FALSE)

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
    nL = ncol(object_geno)
    object_geno <- Remove_NA_Mono_SNP(object_geno)
    total_mono = nL - ncol(object_geno)

  #   ## Check for Na and remove
  #   if(isTRUE(qc_filtering)){
  #     if(!is.null(call_rate)){
  #   CR <- colMeans(is.na(object_geno))
  #   #poscr <- CR <= call_rate #selected by CR
  #
  #
  #   poscr <-  which(CR<call_rate)
  #
  #   nL = ncol(object_geno)
  #
  #   object_geno <- object_geno[, poscr]
  #   ## Gather Meta data
  #   markers_callrate_removed <- (nL - ncol(object_geno))
  #
  #
  #   if(isTRUE(message)){
  #     #message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',ncol(object_geno) - length(which(MAF<maf))))))
  #     message(paste(msg,(paste(paste('Total loci above user defined call rate ', paste("threshold", call_rate)), sep = ': \t',markers_callrate_removed))))
  #   }
  #
  #     }
  # } else {
  #   ## This is mandatory irrespective even if qc_filtering is defined FALSE by the user
  #   ###
  #   Na_col.omit <- which((colSums(is.na(object_geno))==0)==FALSE)
  #
  #   if (length(Na_col.omit)!=0){
  #     #object_geno = object_geno[ , colSums(is.na(object_geno))==0]
  #
  #     if(isTRUE(message)){
  #       message("A total of ", length(Na_col.omit),
  #               " SNP (s)/ marker (s) were removed from due to missing value")
  #
  #     }
  #
  #     object_geno = object_geno[ , -Na_col.omit]
  #   }
  #
  # }
    if(isTRUE(qc_filtering)){
    ###  # Minor alele frequency
    #### MAF #############
    #compute p
    #Here the allele fre. of the second allele (denoted as phat) is calculated and if phat is less than 0.5 , phat = maf, if not maf = 1-phat
   if(!is.null(maf)){
    phat=colMeans(object_geno, na.rm = T)/2
    MAF=ifelse(phat<0.5,phat,1-phat) ##

    if(isTRUE(message)){
    message("\nMinor allele frequency check (MAF).  \n")
    }

    if(length(which(MAF<maf))!=0){

     nl = ncol(object_geno)

    object_geno =object_geno[,-which(MAF<maf)]  ##0r (MAF ??? 0.01)

    maf_markers_removed =  nl - ncol(object_geno)

    if(maf==0.05){
      if(isTRUE(message)){
        message(paste(msg,"if higher or lower MAF threshold value is desired provide one"))
      }
    }

    if(isTRUE(message)){
      #message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',ncol(object_geno) - length(which(MAF<maf))))))
      #message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',length(which(MAF<maf))))))
      message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t', maf_markers_removed))))
    }

    #dim(object_geno)
    } else {
      if(isTRUE(message)){

        message(paste(msg, paste("\tNo loci with MAF below",paste("threshold =", maf), sep = " ")))
      #cat("\tNo loci with minor allele frequency below", paste("threshold =", maf),"\n")
      }
    }
} ## End MAF

    #if (any(is.na(object_geno))) {stop(print(paste(msg,'Marker data contains Missing value')), call. = FALSE)}

    # Heterozigosity
  if(!is.null(freqHet)){
    Het <- apply(object_geno, 2, function(x) sum( x== 1, na.rm=T)/nrow(object_geno))

    if(length(which(Het<=freqHet))!=0){

      N_int = ncol(object_geno)
    object_geno = object_geno[, which(Het<=freqHet)]

    het_markers_removed = N_int - ncol(object_geno)

    if(isTRUE(message)){
      message(paste(msg,(paste(paste('Total loci above heterozygosity', paste("threshold", freqHet)), sep = ': \t', het_markers_removed))))
    }

    }
  } ## End hetero


      ### Removing NA

      if(!is.null(call_rate)){
        CR <- colMeans(is.na(object_geno))
        #poscr <- CR <= call_rate #selected by CR


        poscr <-  which(CR<call_rate)

        nL = ncol(object_geno)

        object_geno <- object_geno[, poscr]
        ## Gather Meta data
        markers_callrate_removed <- (nL - ncol(object_geno))


        if(isTRUE(message)){
          #message(paste(msg,(paste(paste('Total loci below MAF ', paste("threshold", maf)), sep = ': \t',ncol(object_geno) - length(which(MAF<maf))))))
          message(paste(msg,(paste(paste('Total loci above user defined call rate ', paste("threshold", call_rate)), sep = ': \t',markers_callrate_removed))))
        }

      }

      ## If there is missing value and the impute is set to FALSE
      if(any(is.na(object_geno)) && isFALSE(impute)){
        stop(print(paste(msg,'Marker/snp data cannot have missing value. Set impute to TRUE to impute missing value')), call. = FALSE)
        }
      if(isTRUE(impute)){
        if(any(is.na(object_geno))==T) {
          #This function will impute the missing values
          for(j in 1:ncol(object_geno)){
            if(j%%1000==0) cat("Marker=",j,"\n")
            tmp <- object_geno[,j]
            object_geno[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
          }

        }

      }



      # End

      ## End if qc_filtering is TRUE
    } else {
      ## This is mandatory irrespective even if qc_filtering is defined FALSE by the user
      ###
      Na_col.omit <- which((colSums(is.na(object_geno))==0)==FALSE)

      if (length(Na_col.omit)!=0){
        #object_geno = object_geno[ , colSums(is.na(object_geno))==0]

        if(isTRUE(message)){
          message("A total of ", length(Na_col.omit),
                  " SNP (s)/ marker (s) were removed from due to missing value")

        }

        object_geno = object_geno[ , -Na_col.omit]
      }

    }


    ## End if qc_filtering is TRUE
    # Assign appropriate class.
    class(object_geno) <-c("matrix", "array", "geno_data")

    attr(object_geno, "cleared") <- "pass"

  } else {

    stop(print(paste(msg,'Marker/snp data cannot be empyt')), call. = FALSE)
  }

  #### Aggregate all the maker data
  total_number_genotype = nrow(object_geno)
  if(is.null(call_rate)) {
    call_rate  = 0
  } else {
    call_rate = call_rate
    }

  if(!exists("markers_callrate_removed")) {
    total_markers_removed  = 0
  } else {

    total_markers_removed  =  markers_callrate_removed
    }
  #if(!exists("ind_callrate_removed ")) {ind_callrate_removed  = 0}
  if(is.null(freqHet)) {
    heterozygosity  = 0
  } else {
    heterozygosity = freqHet
    }
  if(!exists("het_markers_removed ")) {
    total_het_markers_removed  = 0
  } else {
    total_het_markers_removed = het_markers_removed
    }
  if(is.null(maf)) {
    maf  = 0
  } else {

    maf = maf
  }
  if(!exists("maf_markers_removed")) {
    total_maf_markers_removed = 0
  } else {

    total_maf_markers_removed = maf_markers_removed
  }
  nsnp = ncol(object_geno)
  output = list(object_geno,
                total_number_genotype,
                call_rate,
                total_markers_removed,
                heterozygosity,
                total_het_markers_removed,
                maf,
                total_maf_markers_removed,
                total_mono,
                nsnp)


  names(output) <- c("marker_matrix",
                    "total_number_genotype",
                    "call_rate",
                    "total_markers_removed",
                    "heterozygosity",
                    "total_het_markers_removed",
                    "maf",
                    "total_maf_markers_removed",
                    "total_monomorphic_markers_removed",
                    "total_number_of_marker_after_qc_filtering"

  )


  return(output)

}

