
#' Title
#' This check for NA, match the grm/kernel and M_matrix to pheno_data
#' It also order the M_matrix, grm/kernel to the GID in pheno_data
#'  Ordering it not important in asreml but important in other machine.
#'  NOTE ASRgenomics only works with grm/kernel not with M_matrix
#' @param object_pheno
#' @param object_geno
#' @param gen_name
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
pheno_geno_match <- function(object_pheno = NULL,
                             object_geno = NULL,
                             gen_name = NULL,
                             message = TRUE,
                             ...){

  msg <- sprintf("==================================================\n")
  #ID_geno <- rownames(object_geno)
  ID_pheno <- as.character(unique(object_pheno[, gen_name]))

  ### It is expected that all genotype/individuals with phenotypic records
  ## have SNP/omic records or present in the grm/kernel matrix
  ## However, individuals with SNP, omic, grm/kernel matrix might not have
  ## phenotypic records such individual will be considered the testing set
  ## this is not a problem with asreml but will not work for BGLR.
  ## We can fix this problem by considering those individuals as test_data,
  ## which we already have provision for.
  #### But it good to let the user know if those individuals were actually a mistake

  if(isFALSE(all(ID_pheno%in%rownames(object_geno)))){

    stop(paste(msg, 'Not all individual with phenotypic records has genotypic/omic records.'))


  } else if(isFALSE(all(rownames(object_geno)%in%ID_pheno))){

      message(paste(msg, 'Not all individual with genotypic/omic records has phenotypic records.'))

    test_set <- setdiff(rownames(object_geno), ID_pheno)

    ## if the object_geno is a grm/kernel matrix
    if(nrow(object_geno)==ncol(object_geno)){
      object_geno <- object_geno[c(ID_pheno, test_set),  c(ID_pheno, test_set)]


    } else {

      ## if the object object_geno is a M_matrix data
      if(nrow(object_geno)!=ncol(object_geno)){
        object_geno <- object_geno[c(ID_pheno, test_set),  ]

        ### Conbined them togther. kee in mind object_geno_tst do not have
        ## phenotypic record.
        #### TO DO find way to have them as NA in BGLR or through error message if
        ## the engine if BGLR

      }
    }

  } else {

    if(isFALSE(all(ID_pheno%in%rownames(object_geno)))){

      ## if the object_geno is a grm/kernel matrix
      if(nrow(object_geno)==ncol(object_geno)){
      object_geno <- object_geno[ID_pheno,  ID_pheno]

      } else {

        ## if the object object_geno is a M_matrix data
        if(nrow(object_geno)!=ncol(object_geno)){
          object_geno <- object_geno[ID_pheno,  ]

        }
      }

    }

    }


  #### Declare it also as an object for final usage
  class(object_geno) <-c("matrix", "array", "predictor_clean")

  attr(object_geno, "cleared") <- "model_ready_use"

    if(exists("test_set")){

      output <- list(object_geno, test_set)

      names(output) <- c('object_geno', 'test_set')

    } else {

      output <- list(object_geno)

      names(output) <- c('object_geno')

    }

      return(output)

}


