
#' Title
#' This function perform the following task:
#'
#'1. This check for NA, match the grm/kernel and geno/omic (M_matrix) to pheno_data
#'2. It also order the M_matrix, grm/kernel to the GID in pheno_data
#'3. Ordering is not important in asreml but important in Baysiand and machine learning models.
#'
#'The output is a list with one or two element(s):
#'1. the clean M-matrix or grm file that perfectly match with the phenotypic records
#'2. test_set that might emmanate from the object_geno
#'
#'
#'  NOTE FOR Giovanni ASRgenomics only works with grm/kernel not with M_matrix that was why it not used
#' @param object_pheno clean phenotypic data
#' @param object_geno clean snp/omic/grm/kernel data
#' @param gen_name   column name of the genotypes in the phenotypic data
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
                             test_set = NULL,
                             train_set = NULL,
                             message = TRUE,
                             ...){

  msg <- sprintf("==================================================\n")
  #ID_geno <- rownames(object_geno)
  ID_pheno <- as.character(unique(object_pheno[, gen_name]))

  ### It is expected that all genotype/individuals with phenotypic records
  ## have SNP/omic records or present in the grm/kernel matrix
  ## However, individuals with SNP, omic, grm/kernel matrix might not have
  ## phenotypic records such individual will be considered the testing set
  ## this is not a problem with asreml but will not work for Baysian and machine learning models.
  ## We can fix this problem by considering those individuals as test_data,
  ## which we already have provision for.
  #### But it good to let the user know if those individuals were actually a mistake

  if(isFALSE(all(ID_pheno%in%rownames(object_geno)))){

    ### all Individuals in phenotypic record is expected to have genotypic record
    ## If not stop.
    ## Execption to that is if user provide coefficient of pedigree for such individuals
    ## and construct H-matrix. This is not considered here for now. This will be implemented
    ## in subsequent version
    #stop(paste(msg, 'Not all individual with phenotypic records has genotypic/omic records.'))
    stop(print(paste(msg,'Not all individual with phenotypic records has genotypic/omic records.')), call. = FALSE)


  } else {

  if(isFALSE(all(rownames(object_geno)%in%ID_pheno))){

    if(isTRUE(message)){
      message(insight::print_color(paste(msg,paste('Not all individual with genotypic/omic records has phenotypic records.')), "blue"))

    }
    ### Individuals with genotypic record but with no phenotypic record are assumed to  be
    ## testing set
    test_set <- setdiff(rownames(object_geno), ID_pheno)

    if(!is.null(test_set)){
      if(is.data.frame(test_set) | is.matrix(test_set)){

        ## To be sure the user it not providing duplicate ID
        test_set <-  unique(test_set[, 1])
        #stop(message(paste(msg,"The testing set should be a dataframe with a column named similar to the gen_name provided to represent the genotypes.")), call. = FALSE)
      }else {
        ### It can be character vector or integer vector except list
        if(!is.list(test_set)){
          test_set <-  unique(test_set)

        } else {
          if(is.list(test_set)){
            stop(message(paste(msg,'The testing set cannot be a list. Should be either dataframe, matrix or a vector.')), call. = FALSE)
          }

        }

      }
    } else {

      if(!is.null(train_set) & exists(ID_pheno)){
        if(is.data.frame(train_set) | is.matrix(train_set)){

          ## To be sure the user it not providing duplicate ID
          train_set <-  unique(train_set[, 1])
          test_set = setdiff(ID_pheno, train_set)
          ## If length of the test_set is zero remove it.
          if(length(test_set)==0) rm(test_set)
          message(insight::print_color(paste(msg,paste('No test_set available.')), "red"))

          #stop(message(paste(msg,'The training set should be a dataframe with a column named similar to the gen_name provided to represent the genotypes.')), call. = FALSE)
        }else {

          ### It can be character vector or integer vector except list
          if(!is.list(train_set)){

            ### This is to be use the user if not presenting duplicate ID
            train_set <-  unique(train_set)
            test_set = setdiff(ID_pheno, train_set)
            ## If length of the test_set is zero remove it.
            if(length(test_set)==0) {rm(test_set)
            message(insight::print_color(paste(msg,paste('No test_set available.')), "red"))
            }

          } else{
            if(is.list(train_set)){
            stop(message(paste(msg,'The training set cannot be a list. Should be either dataframe, matrix or a vector.')), call. = FALSE)
            }
          }

        }
      }
    }



    ## if the object_geno is a grm/kernel matrix
    ### First the colname must be equal rowname
    if(nrow(object_geno)==ncol(object_geno)){
      object_geno <- object_geno[c(ID_pheno, test_set),  c(ID_pheno, test_set)]


    } else {

      ## if the object_geno is a M_matrix data
      if(nrow(object_geno)!=ncol(object_geno)){
        object_geno <- object_geno[c(ID_pheno, test_set),  ]

    ### Combined them together. keep in mind test_set do not have phenotypic record.
        #### TO DO find way to have them as NA in BGLR or through error message if
        ## the engine if BGLR

      }
    }

  } else{

    ##
    if(isTRUE(all(rownames(object_geno)%in%ID_pheno))){

      if(nrow(object_geno)==ncol(object_geno)){
        object_geno <- object_geno[ID_pheno, ID_pheno]


      } else {

        ## if the object_geno is a M_matrix data
        if(nrow(object_geno)!=ncol(object_geno)){
          object_geno <- object_geno[ID_pheno,   ]

          ### Combined them together. keep in mind test_set do not have phenotypic record.
          #### TO DO find way to have them as NA in BGLR or through error message if
          ## the engine if BGLR

        }
      }

    }

  }

}

  # else {
  #
  #   if(isFALSE(all(ID_pheno%in%rownames(object_geno)))){
  #
  #     ## if the object_geno is a grm/kernel matrix
  #     if(nrow(object_geno)==ncol(object_geno)){
  #     object_geno <- object_geno[ID_pheno,  ID_pheno]
  #
  #     } else {
  #
  #       ## if the object object_geno is a M_matrix data
  #       if(nrow(object_geno)!=ncol(object_geno)){
  #         object_geno <- object_geno[ID_pheno,  ]
  #
  #       }
  #     }
  #
  #   }
  #
  #   }


  #### Declare it also as an object for final usage
  #class(object_geno) <-c("matrix", "array", "predictor_clean")

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


