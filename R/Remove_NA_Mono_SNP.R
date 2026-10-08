
#' Remove NA and monomorphic SNPs from genomic data
#'
#' Filters a genomic-marker matrix by dropping markers that are entirely missing
#' or monomorphic (a single observed allele state across individuals).
#'
#' @param geno Numeric matrix or data frame of genomic markers with rows =
#'   individuals and columns = SNPs coded as ALT dosage.
#' @param ploidy One positive integer or `"auto"`. Unannotated polyploid
#'   matrices require an explicit value.
#' @param message Logical; if `TRUE`, print a summary of how many SNPs were
#'   dropped (NA-only and monomorphic).
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return The filtered `geno` matrix with NA-only and monomorphic SNPs removed.
#' @export
#'
#' @examples
Remove_NA_Mono_SNP <- function(geno = NULL,
                               message = TRUE,
                               ploidy = "auto",
                               ...
                            ){

  na.counter <- 0
  # SNP data
  #loc.list.All.NA <- array(NA,ncol(geno))
  loc.list.Mono <- array(NA,ncol(geno))
  nL <- ncol(geno)

  N_Individuals<- nrow(geno)


  msg <- ""

  if(isTRUE(message)){

    message(insight::print_color(paste(msg,(paste('Number of Individuals', sep = ': ',N_Individuals))), "blue"))
  #message(paste(msg,(paste('Number of Individuals', sep = ': ',N_Individuals))))
    message(insight::print_color(paste(msg,(paste('Number of markers', sep = ': ', nL))), "blue"))
  #message(paste(msg,(paste('Number of markers', sep = ': ', nL))))
  }

  matrix <- as.matrix(geno)
  if (!is.numeric(matrix)) {
    stop(paste(msg, "The geno data contains non-numeric values"), call. = FALSE)
  }
  lN <- colnames(geno)
  l.names <- colnames(geno)
  resolved_ploidy <- gp_resolve_matrix_ploidy(matrix, ploidy)
  qc_metrics <- geno_qc_metrics(matrix, ploidy = resolved_ploidy)
  if (length(qc_metrics$monomorphic)) {
    loc.list.Mono[qc_metrics$monomorphic] <- lN[qc_metrics$monomorphic]
    all_missing <- which(qc_metrics$marker_missing_rate == 1)
    na.counter <- length(all_missing)
  }


    # Remove NAs from list of monomorphic loci and loci with all NAs
    # loc.list.Mono <- loc.list.Mono[!is.na(loc.list.Mono)]
    #
    # loc.list.All.NA <- loc.list.All.NA[!is.na(loc.list.All.NA)]
    #
    # if (length(loc.list.Mono)!=0){
    #
    #   message(paste(msg,(paste("Alter", sep = ":\n \t", "Markers with NA Present"))))
    # message(paste(msg,(paste('Number of markers with NA removed', sep = ': \t',(nL-length(loc.list))))))
    #
    # } else {
    #
    #   message('No markers with NA in the geno data')
    # }
    #

  # remove monomorphic loc and loci with all NAs

  loc.list.Mono <- loc.list.Mono[!is.na(loc.list.Mono)]

  if(length(loc.list.Mono)!=0){

    if(isTRUE(message)){
      message(insight::print_color(paste(msg,paste("Removing monomorphic loci\n")), "blue"))
    #message(paste(msg,("Removing monomorphic loci\n")))
    }

    # Remove loci flagged for deletion
    geno <- geno[, !colnames(geno) %in% loc.list.Mono, drop = FALSE]
    if(isTRUE(message)){
      message(insight::print_color(paste(msg,(paste('Number of monomorphic loci removed', sep = ': \t',length(loc.list.Mono)))), "blue"))
    #message(paste(msg,(paste('Number of monomorphic loci removed', sep = ': \t',length(loc.list.Mono)))))
}
  } else {

    if(isTRUE(message)){
      message(insight::print_color(paste(msg,paste("No monomorphic loci to remove\n")), "blue"))
    #message(paste(msg,("No monomorphic loci to remove\n")))
    }

  }


  attr(geno, "ploidy") <- resolved_ploidy
  return(geno)

}



