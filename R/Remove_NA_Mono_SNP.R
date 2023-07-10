
#' Title
#'
#' @param geno
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
Remove_NA_Mono_SNP <- function(geno,
                            ...
                            ){

  na.counter <- 0



  # SNP data
  #loc.list.All.NA <- array(NA,ncol(geno))
  loc.list.Mono <- array(NA,ncol(geno))
  nL <- ncol(geno)

  N_Individuals<- nrow(geno)


  msg <- sprintf("==================================================\n")
  message(paste(msg,(paste('Number of Individuals', sep = ': ',N_Individuals))))

  message(paste(msg,(paste('Number of markers', sep = ': ', nL))))

  matrix <- as.matrix(geno)
  lN <- colnames(geno)
  l.names <- colnames(geno)
  for (i in 1:nL){
    row <- matrix[,i] # Row for each locus

    # if (is.na(all(row))) {
    #
    #   loc.list.All.NA[i] <- lN[i]
    #   if (all(is.na(row))){
    #     na.counter = na.counter + 1
    #   }
    # }
    if (all(row == 0, na.rm=TRUE) | all(row == 1, na.rm=TRUE) | all(row == 2, na.rm=TRUE) | all(is.na(row))){
      loc.list.Mono[i] <- lN[i]
       if (all(is.na(row))){
         na.counter = na.counter + 1
       }
    }
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

    message(paste(msg,("Removing monomorphic loci\n")))

    # Remove loci flagged for deletion
    geno <- geno[,!colnames(geno)%in%loc.list.Mono]

    message(paste(msg,(paste('Number of monomorphic loci removed', sep = ': \t',length(loc.list.Mono)))))

  } else {

    message(paste(msg,("No monomorphic loci to remove\n")))

  }


  return(geno)

}



