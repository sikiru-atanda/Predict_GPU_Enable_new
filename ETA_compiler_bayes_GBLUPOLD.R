
#' Title
#'
#' @param fixed
#' @param random
#' @param GS_model
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param pheno_data
#' @param gen_name
#' @param ...
#' @param gmatrix
#' @param gkernel
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#'
#' @return
#' @export
#'
#' @examples

ETA_compiler_bayes_GBLUPOLD <- function(
    fixed = NULL,
    random = NULL,
    GS_model = NULL,
    fixed_term_model_bayesian = NULL,
    rand_term_model_bayesian = NULL,
    pheno_data = NULL,
    gmatrix= NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    gen_name = NULL,
    ...
){

  ### Create empty list for ETA
  ETA = list()

  msg <- "\n==================================================\n"
  ### Get the random terms. Here no interaction terms in the random effect
  rand_term_no_inter <- random_terms(random = random,
                                     object = pheno_data)

  ### Assign
  rand_model <- random_term_model(rand_terms = rand_term_no_inter,
                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                  GS_model = GS_model)

  if(!is.null(fixed)){
    ETA <-  ETA_compiler_fixed_term(fixed = fixed,
                                    pheno_data = pheno_data,
                                    fixed_term_model_bayesian = fixed_term_model_bayesian)


  }

  if(length(rand_model)==1 && length(rand_term_no_inter)==1){

    if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(!is.null(gmatrix)){
        if(rand_model== "RKHS"){
      ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                 model=rand_model,
                                 saveEffects=TRUE)
        } else {

          if(rand_model== "BRR"){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }
        }

      ETA_element_name = c("gmatrix")

      } else {

        if(!is.null(gkernel)){

          if(rand_model== "RKHS"){
            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          } else {

            if(rand_model== "BRR"){
              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }
          }

          ETA_element_name = c("gkernel")
        }
      }

    } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(rand_model== "RKHS"){
        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)
      } else {

        if(rand_model== "BRR"){
          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)
        }
      }

      ETA_element_name = c("omic1_kernel")

    } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

        if(rand_model== "RKHS"){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)
        } else {

          if(rand_model== "BRR"){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }
        }


      ETA_element_name = c("omic2_kernel")

    } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){
        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)
      } else {

        if(rand_model== "BRR"){
          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)
        }
      }


      ETA_element_name = c("omic3_kernel")
      #If user provide geno_data and omic1_data

    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
        ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA_element_name = c("gkernel", "omic1_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic1_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic1_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic1_kernel")
          }
        }

      }


      #ETA_element_name = c("gmatrix", "omic1_kernel")

    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA_element_name = c("gkernel", "omic2_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic2_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic2_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic2_kernel")
          }
        }

      }

      #ETA_element_name = c("geno_data", "omic2_data")

    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA_element_name = c("gkernel", "omic3_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic3_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic3_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic3_kernel")
          }
        }

      }

      #ETA_element_name = c("geno_data", "omic3_data")

    } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)


      } else {

        if(rand_model== "BRR"){

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

        }

      }

      ETA_element_name = c("omic1_kernel", "omic2_kernel")
      #ETA_element_name = c("omic1_data", "omic2_data")
    } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)


      } else {

        if(rand_model== "BRR"){

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

        }

      }

      ETA_element_name = c("omic1_kernel", "omic3_kernel")

      #ETA_element_name = c("omic1_data", "omic3_data")
    } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)


      } else {

        if(rand_model== "BRR"){

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

        }

      }

      ETA_element_name = c("omic2_kernel", "omic3_kernel")

      #ETA_element_name = c("omic2_data", "omic3_data")
    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
          }
        }

      }

      #ETA_element_name = c("geno_data", "omic1_data", "omic2_data")

    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
          }
        }

      }
      #ETA_element_name = c("geno_data", "omic1_data", "omic3_data")
      ###
    } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        if(!is.null(gkernel)){
          ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
        } else {

          if(!is.null(gmatrix)){

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)
          }

          ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
        }

      } else {

        if(rand_model== "BRR"){
          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
          }
        }

      }

      #ETA_element_name = c("geno_data", "omic2_data", "omic3_data")

    } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      if(rand_model== "RKHS"){

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)

        ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                   model=rand_model,
                                   saveEffects=TRUE)


      } else {

        if(rand_model== "BRR"){

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

          ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                     model=rand_model,
                                     saveEffects=TRUE)

        }

      }

      ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")



    } else {

      if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

        if(rand_model== "RKHS"){

          if(!is.null(gkernel)){
            ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                       model=rand_model,
                                       saveEffects=TRUE)

            ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
          } else {

            if(!is.null(gmatrix)){

              ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)
            }

            ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
          }

        } else {

          if(rand_model== "BRR"){
            if(!is.null(gkernel)){
              ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                         model=rand_model,
                                         saveEffects=TRUE)

              ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
            } else {

              if(!is.null(gmatrix)){

                ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
                                           model=rand_model,
                                           saveEffects=TRUE)

                ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
                                           model=rand_model,
                                           saveEffects=TRUE)

                ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_kernel),
                                           model=rand_model,
                                           saveEffects=TRUE)

                ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
                                           model=rand_model,
                                           saveEffects=TRUE)
              }

              ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
            }
          }

        }

      #ETA_element_name = c("geno_data", "omic1_data", "omic2_data", "omic3_data")
    }

    }

  } else {

    stop(paste(msg, 'This works only for one random effect'))


  }

  output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)

  names(output) <- c("ETA", "pheno_data", "ETA_element_name")

  return(output)



}


