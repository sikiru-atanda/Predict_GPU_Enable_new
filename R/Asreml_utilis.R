#' Title
#'
#' @param response
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param VarCov_str
#' h@param weights
#' @param GS_model
#' @param core
#' @param weights
#' @param pheno_data
#' @param gmatrix
#' @param gkernel
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param workspace
#' @param maxit
#' @param gen_name
#' @param ...
#'
#' @return
#' @examples
#' @importFrom foreach %dopar%

asreml_utilis <- function(
    fixed = NULL,
    random = NULL,
    cova=NULL,
    GS_model = NULL,
    response = NULL,
    pheno_data = NULL,
    gmatrix = NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    inverse = NULL,
    epsilon = TRUE,
    gen_name = NULL,
    heter_groups = NULL,
    heter_resid = FALSE,
    VarCov_str = NULL,
    weights = NULL,
    core = NULL,
    workspace=1e08,
    engine = NULL,
    #pworkspace= 1e06,
    maxit = 50,
    ...
) {

  msg <- sprintf("==================================================\n")

  if(engine %in% rownames(installed.packages())){
    do.call('library', list(engine))

    # if package is not installed locally then stop
  } else {

    stop(print(paste(msg,'You need to install asreml-R to use asreml-R')), call. = FALSE)
  }

    ######

  #   if (!is.null(gmatrix) | !is.null(gkernel)){
  #       if(!is.null(gkernel)){
  #         gmatrix <- gkernel
  #         rm(gkernel); gc()
  #       }
  #     G_inv <-  chol2inv(chol(gmatrix))
  #     G_inv <- sparse_matrix(grm_kernel_data = G_inv)
  #     attr(G_inv, "INVERSE") <- TRUE
  #
  #     # attr(G_inv, "INVERSE") <- TRUEG_inv  <- gmatrix
  #     # G_inv <-  chol2inv(chol(gmatrix))
  #     # rownames(G_inv) <- rownames(gmatrix)
  #     # colnames(G_inv) <- colnames(gmatrix)
  #     # # asreml requires the Gmatrix to be inverted in sparse format
  #     # # G_inv <- solve(G_bend)
  #     # # #check the attribute rowNames
  #     # attr(gmatrix, "rowNames") <- rownames(gmatrix)
  #     # attr(gmatrix, "colNames") <- colnames(gmatrix)
  #     # attr(G_inv, "rowNames") <- rownames(gmatrix)
  #     # attr(G_inv, "colNames") <- colnames(gmatrix)
  #     # attr(G_inv, "INVERSE") <- TRUE
  #
  #     rm(gmatrix)
  #   }
  #
  #   # if (!is.null(gkernel)){
  #   #   # GK_inv <-  chol2inv(chol(gkernel))
  #   #   # GK_inv <- sparse_matrix(grm_kernel_data = GK_inv)
  #   #   GK_inv <- gkernel
  #   #   # GK_inv <-  chol2inv(chol(gkernel))
  #   #   # rownames(GK_inv) <- rownames(gkernel)
  #   #   # colnames(GK_inv) <- colnames(gkernel)
  #   #   # attr(gkernel, "rowNames") <- rownames(gkernel)
  #   #   # attr(gkernel, "colNames") <- colnames(gkernel)
  #   #   # attr(GK_inv, "rowNames") <- rownames(gkernel)
  #   #   # attr(GK_inv, "colNames") <- colnames(gkernel)
  #   #   # attr(GK_inv, "INVERSE") <- TRUE
  #   #
  #   #   rm(gkernel)
  #   # }
  #
  #   ####
  #   if (!is.null(omic1_kernel)){
  #     omic1_inv <-  chol2inv(chol(omic1_kernel))
  #     omic1_inv <- sparse_matrix(grm_kernel_data = omic1_inv)
  #     attr(omic1_inv, "INVERSE") <- TRUE
  #     #omic1_inv <- omic1_kernel
  #     # omic1_inv <<-  chol2inv(chol(omic1_kernel))
  #     # rownames(omic1_inv) <- rownames(omic1_kernel)
  #     # colnames(omic1_inv) <- colnames(omic1_kernel)
  #     # attr(omic1_kernel, "rowNames") <- rownames(omic1_kernel)
  #     # attr(omic1_kernel, "colNames") <- colnames(omic1_kernel)
  #     # attr(omic1_inv, "rowNames") <- rownames(omic1_kernel)
  #     # attr(omic1_inv, "colNames") <- colnames(omic1_kernel)
  #     # attr(omic1_inv, "INVERSE") <- TRUE
  #
  #     rm(omic1_kernel)
  #   }
  #   #############
  # if (!is.null(omic2_kernel)){
  #   omic2_inv <-  chol2inv(chol(omic2_kernel))
  #   omic2_inv <- sparse_matrix(grm_kernel_data = omic2_inv)
  #   attr(omic2_inv, "INVERSE") <- TRUE
  #   #omic2_inv <- omic2_kernel
  #   # omic2_inv <<-  chol2inv(chol(omic2_kernel))
  #   # rownames(omic2_inv) <- rownames(omic2_kernel)
  #   # colnames(omic2_inv) <- colnames(omic2_kernel)
  #   # attr(omic2_kernel, "rowNames") <- rownames(omic2_kernel)
  #   # attr(omic2_kernel, "colNames") <- colnames(omic2_kernel)
  #   # attr(omic2_inv, "rowNames") <- rownames(omic2_kernel)
  #   # attr(omic2_inv, "colNames") <- colnames(omic2_kernel)
  #   # attr(omic2_inv, "INVERSE") <- TRUE
  #
  #   rm(omic2_kernel)
  # }
  # ####
  # if (!is.null(omic3_kernel)){
  #   omic3_inv <-  chol2inv(chol(omic3_kernel))
  #   omic3_inv <- sparse_matrix(grm_kernel_data = omic3_inv)
  #   attr(omic3_inv, "INVERSE") <- TRUE
  #   #omic3_inv <- omic3_kernel
  #   # omic3_inv <<-  chol2inv(chol(omic3_kernel))
  #   # rownames(omic3_inv) <- rownames(omic3_kernel)
  #   # colnames(omic3_inv) <- colnames(omic3_kernel)
  #   # attr(omic3_kernel, "rowNames") <- rownames(omic3_kernel)
  #   # attr(omic3_kernel, "colNames") <- colnames(omic3_kernel)
  #   # attr(omic3_inv, "rowNames") <- rownames(omic3_kernel)
  #   # attr(omic3_inv, "colNames") <- colnames(omic3_kernel)
  #   # attr(omic3_inv, "INVERSE") <- TRUE
  #
  #   rm(omic3_kernel)
  # }
#   #######
#
#     if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic3_inv")) & !exists("omic3_inv"))){
#
#       #G_copy = 1
#       if(exists("G_inv")){
#
#       G_list <<- list('G')
#
#       } else {
#
#         if(exists("GK_inv")){
#
#           #G_list <<- list('GK')
#           G_list <<- list('G')
#
#         }
#
#       }
#
#
#     }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic3_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 1
#
#     G_list <<- list('omic1')
#   }
#
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 1
#
#     G_list <<- list('omic2')
#   }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 1
#
#     G_list <<- list('omic3')
#   }
#
#   #####
#   if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 2
#     if(exists("G_inv")){
#     G_list <<- list('G', 'omic1')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic1')
#         #G_list <<- list('GK', 'omic1')
#
#       }
#
#     }
#   }
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 2
#
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic2')
#
#     } else {
#       if(exists("GK_inv")){
#         #G_list <<- list('GK', 'omic2')
#         G_list <<- list('G', 'omic2')
#       }
#
#     }
#
#   }
#
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 2
#
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic3')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic3')
#
#         #G_list <<- list('GK', 'omic3')
#
#       }
#
#     }
#   }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 2
#
#     G_list <<- list('omic1', 'omic2')
#   }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 2
#
#     G_list <<- list('omic1', 'omic3')
#   }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 2
#
#     G_list <<- list('omic2', 'omic3')
#   }
#
#   if ((!exists("G_inv") | !exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 3
#
#     G_list <<- list('omic1', 'omic2', 'omic3')
#   }
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & !exists("omic3_inv"))){
#
#     #G_copy = 3
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic1', 'omic2')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic1', 'omic2')
#         #G_list <<- list('GK', 'omic1', 'omic2')
#
#       }
#
#     }
#
#
#   }
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & !exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 3
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic1', 'omic3')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic1', 'omic3')
#         #G_list <<- list('GK', 'omic1', 'omic3')
#
#       }
#
#     }
#
#
#   }
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((!exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 3
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic2', 'omic3')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic2', 'omic3')
#         #G_list <<- list('GK', 'omic2', 'omic3')
#       }
#
#     }
#
#   }
#
#   if ((exists("G_inv") | exists("GK_inv")) & ((exists("omic1_inv") & exists("omic2_inv")) & exists("omic3_inv"))){
#
#     #G_copy = 4
#     if(exists("G_inv")){
#       G_list <<- list('G', 'omic1', 'omic2', 'omic3')
#
#     } else {
#       if(exists("GK_inv")){
#         G_list <<- list('G', 'omic1', 'omic2', 'omic3')
#         #G_list <<- list('GK', 'omic1', 'omic2', 'omic3')
#       }
#
#     }
#
#
#   }



  # Check and assign values for G_inv and GK_inv
  # if (!is.null(gmatrix)) {
  #   G_inv <- gmatrix
  #   rm(gmatrix)
  # }
  #
  # if (!is.null(gkernel)) {
  #   G_inv <- gkernel
  #   rm(gkernel)
  # }
  # if (!is.null(omic1_kernel)) {
  #   omic1_inv <- omic1_kernel
  #   rm(omic1_kernel)
  # }
  # if (!is.null(omic2_kernel)) {
  #   omic2_inv <- omic2_kernel
  #   rm(omic2_kernel)
  # }
  #
  # if (!is.null(omic3_kernel)) {
  #   omic3_inv <- omic3_kernel
  #   rm(omic3_kernel)
  # }

  # # Check and assign values for omic_inv variables
  # assign_omic_inv <- function(omic_kernel, index) {
  #   if (!is.null(omic_kernel)) {
  #     assign(paste0("omic", index, "_inv"), omic_kernel, envir = .GlobalEnv)
  #     #rm(list = paste0("omic", index, "_kernel"), envir = .GlobalEnv)
  #
  #   }
  # }
  #
  #
  #
  # assign_omic_inv(omic1_kernel, 1)
  # assign_omic_inv(omic2_kernel, 2)
  # assign_omic_inv(omic3_kernel, 3)

  # Function to compute inverse and create sparse matrix
  # Regularization parameter (adjust as needed)

  compute_inverse_and_sparse <- function(kernel,
                                         epsilon = 1e-6,
                                         inverse =NULL) {
    attr(kernel, "rowNames") <- rownames(kernel)
    attr(kernel, "colNames") <- colnames(kernel)
    kernel_extra <- kernel

    if(isTRUE(inverse)){
      # Apply regularization technique when the matrix is nearly positive definite
      # but has small negative eigenvalues due to numerical precision issues
      # to make it positive definite by adding
      ## small positive constant to the diagonal elements of the matrix
      ## Two ways of acheiving that are define here
      # 1) Apply it directly to the directly to the diagonal
          # elements of the matrix is indeed a form of regularization
      # 2) Ridge regularization
      result <- tryCatch(
        {
          diag(kernel) <- diag(kernel) + epsilon
          inverse_matrix <- chol2inv(chol(kernel))

          # Return the result
          inverse_matrix
        },
        error = function(e) {
          cat("Error occurred during computation:", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      # Check if an error occurred
      if(inherits(result, "try-error")) {
        ### # Add a small positive constant to the diagonal elements
        # Modify the kernel if needed
        kernel <- kernel + epsilon * diag(nrow(kernel))
        inverse_matrix <- chol2inv(chol(kernel))
       } else {
        inverse_matrix <- result
       }

      # Create sparse matrix
      # sparse <- sparse_matrix(grm_kernel_data = inverse_matrix)
      # attr(sparse, "rowNames") <- rownames(kernel)
      # attr(sparse, "colNames") <- colnames(kernel)
      # attr(sparse, "INVERSE") <- TRUE
      attr(inverse_matrix, "rowNames") <- rownames(kernel)
      attr(inverse_matrix, "colNames") <- colnames(kernel)
      attr(inverse_matrix, "INVERSE") <- TRUE
    } else {
      sparse <- sparse_matrix(grm_kernel_data = kernel)
      attr(sparse, "rowNames") <- rownames(kernel)
      attr(sparse, "colNames") <- colnames(kernel)
      #attr(sparse, "INVERSE") <- TRUE
    }

    return(sparse)
  }

  # Process gmatrix or gkernel
  if (!is.null(gmatrix) || !is.null(gkernel)) {
    if (!is.null(gkernel)) {
      gmatrix <- gkernel
      rm(gkernel)
    }
    G_inv <- compute_inverse_and_sparse(gmatrix, inverse)
    rm(gmatrix)
  }

  # Process omic1_kernel
  if (!is.null(omic1_kernel)) {
    omic1_inv <- compute_inverse_and_sparse(omic1_kernel, inverse)
    rm(omic1_kernel)
  }

  # Process omic2_kernel
  if (!is.null(omic2_kernel)) {
    omic2_inv <- compute_inverse_and_sparse(omic2_kernel, inverse)
    rm(omic2_kernel)
  }

  # Process omic3_kernel
  if (!is.null(omic3_kernel)) {
    omic3_inv <- compute_inverse_and_sparse(omic3_kernel, inverse)
    rm(omic3_kernel)
  }


  # Determine the value of G_list based on the existence of G_inv, GK_inv, and omic_inv variables
  G_list <- character()

  if (exists("G_inv")) {
    G_list <- c(G_list, if (exists("G_inv")) "G")
    #rm(list = paste0("G", "_inv"), envir = .GlobalEnv)
  }

  for (i in 1:3) {
    if (exists(paste0("omic", i, "_inv"))) {
      G_list <- c(G_list, paste0("omic", i))
      #rm(list = paste0("omic", i, "_inv"), envir = .GlobalEnv)
    }
  }

  # Optional: Remove G_inv, GK_inv, and omic_inv variables
  # rm(list = c("G_inv", "GK_inv", paste0("omic", 1:3, "_inv")))

  G_list = as.list(G_list)


  ## Here pworkspace is not included because predict function is not done here
    asreml::asreml.options(trace=FALSE,
                           #workspace = workspace,
                           #pworkspace = pworkspace,
                           maxit = maxit)


    #### When the gen_name are present in more than one environment/location
    if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      if (!is.null(heter_resid)) {
        if (is.null(heter_groups)) {
          stop(print(paste(msg,'No column of heterogeneous groups provided.')), call. = FALSE)
        } else {

          if(!heter_groups%in%colnames(pheno_data)) {stop(print(paste(msg,'heterogenous group provided did not match column names in pheno_data')), call. = FALSE)}

          if(!all(sapply(heter_groups, function(x, pheno_data) is.factor(pheno_data[,x]),  pheno_data))) {
            pheno_data[,heter_groups] <- as.factor(pheno_data[,heter_groups])
          }
        }

      }

    }

##### THis is important for asreml inorder to update the model if need be
    if (exists("G_inv")){.GlobalEnv$G_inv <- G_inv}
    #if (exists("GK_inv")){.GlobalEnv$GK_inv <- GK_inv}
    if (exists("omic1_inv")){.GlobalEnv$omic1_inv <- omic1_inv}
    if (exists("omic2_inv")){.GlobalEnv$omic2_inv <- omic2_inv}
    if (exists("omic3_inv")){.GlobalEnv$omic3_inv <- omic3_inv}

    #a <- a + 1
    # Code Strings for all factors y= XB + UZ + e
    code_asr <- as.character()

    # colnames(pheno_data)[colnames(pheno_data) == trait] <-
    #   deparse(substitute(trait))

    code_asr[1] <- paste0(paste('asreml::asreml(fixed=', 'trait'),  '~1')
    code_asr[2] <- 'random=~'
    code_asr[3] <- 'residual=~'

    # Adding covariates (fixed)
    if (!is.null(cova)) {

      cova_term <- strsplit(as.character(cova[2]), split = "[+]")[[1]]

      if (length(cova_term )>1) {

        for (c in 1:length(cova_term )) {
          code_asr[1] <- paste(code_asr[1], cova_term[c], sep='+')
        }

      } else {
        code_asr[1] <- paste(code_asr[1], cova_term, sep='+')
        #code.asr[1] <- paste(code.asr[1], cova)

      }

    }

    # Adding fixed factors
    if (!is.null(fixed)) {
      fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
      if (length(all.vars(fixed))>1) {
        for (v in 1:length(all.vars(fixed))) {
          code_asr[1] <- paste(code_asr[1], all.vars(fixed)[v], sep='+')
        }

      } else {
        code_asr[1] <- paste(code_asr[1], all.vars(fixed), sep='+')
      }

    } else{
      if (is.null(fixed)) {
        fixed_term = NULL
        check_heter_grp_fixed = NULL
      }
    }


    # Adding random factors

    if (is.null(random)){ stop(print('provide random term'))}

   rand_out <-  random_terms_fit(random = random,
                     fixed = fixed,
                     fixed_term = fixed_term,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     VarCov_str = VarCov_str,
                     code_asr = code_asr,
                     G_list = G_list,
                     gen_name = gen_name,
                     pheno_data = pheno_data)

   code_asr_fit <-  rand_out[[1]]
   Gen_pos <-  rand_out[[2]]
   Inter_Gen_pos <-  rand_out[[3]]
   rand_term <- rand_out[[4]]

    # if (!is.null(random)){
    #   #if (length(all.vars(random))>1) {
    #   rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    #   rand_term = gsub(" ", "", rand_term)
    #
    #
    #   #rand_InterPresent = grep(":", rand_term)
    #   ## No interaction term
    #   rand_term_No_Inter = rand_term[!grepl(":", rand_term)]
    #
    #   ## Interaction term
    #   Check_rand_Inter = rand_term[grepl(":", rand_term)]
    #
    #   ## Extract other side expect the gen_name
    #   rand_term_No_Inter_No_Gen = rand_term_No_Inter[!rand_term_No_Inter%in%gen_name]
    #   ## Get the position of other terms in the random that is not gen_name
    #   Non_Gen_pos = match(rand_term_No_Inter_No_Gen, rand_term)
    #   if(anyNA(Non_Gen_pos)){Non_Gen_pos = NULL}
    #
    #   if(length(Check_rand_Inter)!=0){
    #
    #     TestPresentofGeno = grep(gen_name, Check_rand_Inter, value = TRUE)
    #     if(anyNA(TestPresentofGeno)) {TestPresentofGeno = NULL}
    #
    #     if(length(TestPresentofGeno)!=0){
    #
    #       Inter_Gen_pos = match(TestPresentofGeno, rand_term)
    #       if(anyNA(Inter_Gen_pos)){Inter_Gen_pos = NULL}
    #
    #       Non_Gen_Inter_Test =  Check_rand_Inter[!Check_rand_Inter%in%TestPresentofGeno]
    #       if(anyNA(Non_Gen_Inter_Test)){Non_Gen_Inter_Test = NULL}
    #
    #     } else {
    #
    #       Non_Gen_Inter_Test = NULL
    #
    #       Inter_Gen_pos = NULL
    #     }
    #
    #     if(length(Non_Gen_Inter_Test)!=0){
    #
    #       Inter_Non_Gen_pos = match(Non_Gen_Inter_Test, rand_term)
    #
    #     }
    #
    #     if(is.null(Inter_Gen_pos)){warning(paste(gen_name, 'missing in the interaction term. This implies you cannot estimate GxE'))}
    #
    #   } else {
    #
    #     if(length(Check_rand_Inter)==0){
    #
    #       Inter_Gen_pos = NULL
    #     }
    #   }
    #
    #   Gen_pos=  match(gen_name, rand_term)
    #   if(anyNA(Gen_pos)){Gen_pos = NULL}
    #   if(length(Gen_pos)>1){
    #     stop(message(paste(msg, "Genotype main effect cannot be present more than one time in the model")), call. = FALSE)
    #     }
    #
    #   ## Eg when GID:Env with no variance structure is the only term in the
    #   ## random effect provided by the user
    #   if(length(rand_term)==length(Check_rand_Inter) & is.null(VarCov_str)){
    #     stop(message(paste(msg, "provide variance-covariance structure")), call. = FALSE)
    #
    #   }
    #
    #
    #   if (length(Check_rand_Inter)>=1){
    #
    #     if (is.null(heter_groups)) {stop(print(paste(msg, "hetero.groups cannot be NULL")), call. = FALSE)}
    #   }
    #
    #
    #   if(!is.null(VarCov_str)){
    #
    #     if(is.null(heter_groups)){stop(print(paste(msg, "Provide heter_groups to model specified variance-covariance structure")), call. = FALSE)}
    #
    #     ### Check if the number of hetero.Grp is greater 5 or greater than 5
    #     NN = nlevels(pheno_data[, heter_groups])
    #     if (NN >=5 & isFALSE(grepl("fa", VarCov_str))){
    #
    #       msg <- sprintf("\r==================================================\n")
    #
    #       warning(paste(msg, "The number of", heter_groups, " is", NN,  "consider using factor analytic model"), immediate. = TRUE, call. =FALSE)
    #     }
    #
    #     ### This check if heterogeneous group/environment/location is present in the fixed term
    #     if (!is.null(fixed_term)){
    #       check_heter_grp_fixed  = match(heter_groups, fixed_term)
    #       if(anyNA(check_heter_grp_fixed)) {check_heter_grp_fixed = NULL}
    #     } else {
    #       check_heter_grp_fixed = NULL
    #     }
    #
    #
    #     if (!is.null(rand_term)){
    #       check_heter_grp_rand  = match(heter_groups, rand_term)
    #       if(anyNA(check_heter_grp_rand)) {check_heter_grp_rand = NULL}
    #     } else {
    #       check_heter_grp_rand = NULL
    #     }
    #
    #     #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
    #
    #     #### When user provide only the Interaction term was provided by the user
    #     #if(!is.null(Inter_Gen_pos) & is.null(Gen_pos)){
    #     if(!is.null(Inter_Gen_pos)){
    #       if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){
    #
    #         #### Add the environment to the fixed term by default if not provided by the user
    #         ## Providing it in fixed term allow for optimal model fit compared to the random term
    #
    #         # Adding fixed factors
    #         if (!is.null(fixed_term)) {
    #
    #
    #           ### Update the fixed term if not null
    #           fixed =  stats::update(fixed,
    #                                  paste("~ . +", heter_groups))
    #
    #           fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    #
    #
    #           for (v in 1:length(all.vars(fixed))) {
    #             code.asr[1] <- paste(code.asr[1], all.vars(fixed)[v], sep='+')
    #           }
    #
    #
    #
    #         } else {
    #           #fixed = heter_groups
    #           fixed_term <- heter_groups
    #           code.asr[1] <- paste(code.asr[1], fixed_term, sep='+')
    #
    #         }
    #
    #         for (i in 1:length(G_list)){
    #
    #           if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {
    #               #if(exists("G_inv") | exists("GK_inv")){
    #               # ### This add the environment to the random term even when when Env is missing in both fixed and random terms
    #               # random= stats::update(random,
    #               #                       paste(paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
    #               #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #               #                                                   sep = ":"))), "+", heter_groups))
    #
    #
    #               ### This add the environment to the fixed term when interaction term was only provided
    #               ## by the user. That is Env is missing in both fixed and random terms
    #               random = stats::update(random,
    #                                     paste("~ . +", (paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
    #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                                                           sep = ":"))))
    #
    #
    #
    #             }else{
    #
    #               if(isTRUE(grepl("fa", VarCov_str))) {
    #                 N_fa = substr(VarCov_str, 3, 100)
    #                 #if(exists("GK_inv") | exists("G_inv")){
    #                 # random= stats::update(random,
    #                 #                       paste(paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
    #                 #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                 #                                                   sep = ":"))), "+", heter_groups))
    #
    #                 random=  stats::update(random,
    #                               paste("~ . +", (paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
    #                                                     paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                                                     sep = ":"))))
    #
    #
    #               }
    #
    #
    #               if(isTRUE(grepl("rr", VarCov_str))) {
    #                 N_rr = substr(VarCov_str, 3, 100)
    #                 #if(exists("GK_inv") | exists("G_inv")){
    #                 # random= stats::update(random,
    #                 #                       paste(paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
    #                 #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                 #                                                   sep = ":"))), "+", heter_groups))
    #
    #                 random= stats::update(random,
    #                                       paste("~ . +", (paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
    #                                                             paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                                                             sep = ":"))))
    #
    #
    #               }
    #
    #             }
    #
    #
    #         } ### End
    #
    #
    #
    #       } ### when Env is missing in both fixed and random terms
    #
    #
    #       ## If user provide ENV in the fixed term and missing in the random term
    #       if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==1){
    #
    #
    #         for (i in 1:length(G_list)){
    #
    #
    #           if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {
    #               random =stats::update(random,
    #                                     paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
    #                                                         paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #             } else{
    #
    #               if(isTRUE(grepl("fa", VarCov_str))) {
    #                 N_fa = substr(VarCov_str, 3, 100)
    #
    #                 random =stats::update(random,
    #                                       paste("~ . +",paste(paste0("fa", paste0("(",paste0(heter_groups, ",", N_fa),")")),
    #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #
    #               }
    #
    #
    #               if(isTRUE(grepl("rr", VarCov_str))) {
    #                 N_rr = substr(VarCov_str, 3, 100)
    #
    #                 random =stats::update(random,
    #                                       paste("~ . +",paste(paste0("rr", paste0("(",paste0(heter_groups, ",", N_rr),")")),
    #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #               }
    #
    #
    #             }
    #
    #
    #
    #         } ### End
    #
    #
    #       } ### End  If user provide ENV in the fixed term and missing in the random term
    #
    #       ## If user provide ENV in the random term and missing in the fixed term
    #       if (length(check_heter_grp_rand)==1 & length(check_heter_grp_fixed)==0){
    #
    #
    #         for (i in 1:length(G_list)){
    #
    #           if (VarCov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {
    #
    #               random= stats::update(random,
    #                                     paste("~ . +",paste(paste0(VarCov_str, paste0("(",heter_groups,")")),
    #                                                         paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #             } else{
    #
    #               if(isTRUE(grepl("fa", VarCov_str))) {
    #                 N_fa = substr(VarCov_str, 3, 100)
    #
    #                 random= stats::update(random,
    #                                       paste("~ . +",paste(paste0("fa",paste0("(",paste0(heter_groups, ",", N_fa),")")),
    #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #               }
    #
    #
    #               if(isTRUE(grepl("rr", VarCov_str))) {
    #                 N_rr = substr(VarCov_str, 3, 100)
    #
    #                 random= stats::update(random,
    #                                       paste("~ . +",paste(paste0("rr",paste0("(",paste0(heter_groups, ",", N_rr),")")),
    #                                                           paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #               }
    #             }
    #
    #
    #         }
    #       }
    #
    #     }
    #       ### If user provide GID:Env
    #       ### use this step to drop the orginal GID:Env
    #       if(!is.null(Inter_Gen_pos) & length(Gen_pos)==0){
    #         rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #
    #         rand_termCopy = gsub(" ", "", rand_termCopy)
    #
    #         #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
    #         if(!is.null(heter_groups)){
    #           Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
    #           if(anyNA(Inter_Gen_pos_copy)){
    #             Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
    #           }
    #
    #         }
    #
    #         random = stats::formula(stats::drop.terms(stats::terms(random),Inter_Gen_pos_copy, keep.response = F))
    #
    #       } ### End
    #
    #       ### If user provide GID, GID:Env
    #       if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
    #         ### use this step to drop the original GID and GID:Env
    #         rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #
    #         rand_termCopy = gsub(" ", "", rand_termCopy)
    #
    #         #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
    #         if(!is.null(heter_groups)){
    #           Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
    #           if(anyNA(Inter_Gen_pos_copy)){
    #             Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
    #           }
    #
    #         }
    #         Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
    #
    #         random = stats::formula(stats::drop.terms(stats::terms(random), c(Gen_pos_copy,Inter_Gen_pos_copy), keep.response = F))
    #
    #       }
    #
    #       #}
    #
    #     #}### End of when only Gen_pos and Inter_Gen_pos are provided
    #     ## That is user provide GID, GID:Env
    #
    #   }
    #
    #   ##################################################
    #   ### If user provide only the gen_name in the random term as gen_name effect
    #   ######################################################
    #
    #   if(is.null(Inter_Gen_pos) & length(Gen_pos)==1){
    #     #if (length(Gen_pos)==1 & length(Inter_Gen_pos)==0){
    #     ### stats::formula must be at least  length 1 so since it just one length. The
    #     ## adjusted random term was added and the old one was removed
    #     if (exists('G_list')){
    #       for (i in 1:length(G_list)) {
    #
    #
    #         random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
    #
    #       }
    #     }
    #     ### use this step to drop the orginal GID and GID:Env
    #     rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #
    #     rand_termCopy = gsub(" ", "", rand_termCopy)
    #
    #     Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
    #
    #     random = stats::formula(stats::drop.terms(stats::terms(random), Gen_pos_copy, keep.response = F))
    #
    #   }
    #   ######################################################################
    #   #### When Variance-Covariance Structure is missing. Compound Symmetry
    #   ####################################################################
    #   if(is.null(VarCov_str) & (length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name])))){
    #
    #     ###############################
    #     if (!is.null(fixed_term)){
    #       check_heter_grp_fixed  = match(heter_groups, fixed_term)
    #     } else {
    #       check_heter_grp_fixed = NULL
    #     }
    #
    #
    #     if (!is.null(rand_term)){
    #       check_heter_grp_rand  = match(heter_groups, rand_term)
    #     } else {
    #       check_heter_grp_rand = NULL
    #     }
    #
    #     #if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv"){
    #
    #     if(!is.null(Inter_Gen_pos)){
    #       if (length(check_heter_grp_rand)==0 & length(check_heter_grp_fixed)==0){
    #
    #         #### Add the environment to the fixed term by default if not provided by the user
    #         ## Providing it in fixed term allow for optimal model fit compared to the random term
    #
    #         # Adding fixed factors
    #         if (!is.null(fixed_term )) {
    #
    #           ### Update the fixed term if not null
    #           fixed =  stats::update(fixed,
    #                                  paste("~ . +", heter_groups))
    #
    #           fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    #
    #
    #           for (v in 1:length(all.vars(fixed))) {
    #             code.asr[1] <- paste(code.asr[1], all.vars(fixed)[v], sep='+')
    #           }
    #
    #         } else {
    #
    #           fixed_term <- heter_groups
    #
    #           code.asr[1] <- paste(code.asr[1], fixed_term, sep='+')
    #
    #         }
    #
    #         ###
    #
    #         for (i in 1:length(G_list)){
    #
    #             # random =stats::update(random, paste(paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))),
    #             #                                     "+", heter_groups))
    #
    #           # random= stats::update(random,
    #           #                       paste(paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
    #           #                                                   paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #           #                                                   sep = ":"))), "+", heter_groups))
    #
    #           random= stats::update(random,
    #                                 paste("~ . +", (paste(paste0("idv", paste0("(",heter_groups,")")),
    #                                                       paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")),
    #                                                       sep = ":"))))
    #
    #         }
    #
    #       }### when Env is missing in both fixed and random terms
    #
    #
    #       ## If user provide ENV in the fixed term and missing in the random term
    #       if (length(check_heter_grp_fixed)==1 & length(check_heter_grp_rand)==0){
    #
    #         for (i in 1:length(G_list)){
    #
    #             #random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
    #
    #           random= stats::update(random,
    #                                 paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
    #                                                     paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #         }
    #
    #       }### End  If user provide ENV in the fixed term and missing in the random term
    #
    #       ## If user provide ENV in the random term and missing in the fixed term
    #       if (length(check_heter_grp_fixed)==0 & length(check_heter_grp_rand)==1){
    #
    #         if (exists('G_list')){Inter_Gen_pos_use = length(G_list)}
    #         for (i in 1:length(G_list)){
    #
    #             #random =stats::update(random, paste("~ . +",paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = ""))))
    #
    #           random= stats::update(random,
    #                                 paste("~ . +",paste(paste0("idv", paste0("(",heter_groups,")")),
    #                                                     paste(paste0('vm(', gen_name), sep = ',', paste(G_list[[i]], "_inv)", sep = "")), sep = ":")))
    #
    #
    #         }
    #
    #       }
    #
    #
    #     }
    #
    #     ### use this step to drop the orginal GID and GID:Env
    #     # if(!is.null(Gen_pos) & is.null(Inter_Gen_pos)){
    #     # rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #     #
    #     # rand_termCopy = gsub(" ", "", rand_termCopy)
    #     #
    #     # Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
    #     #
    #     # random = stats::formula(stats::drop.terms(stats::terms(random), Gen_pos_copy, keep.response = F))
    #     # }
    #
    #     ### If user provide GID, GID:Env
    #     if(!is.null(Inter_Gen_pos) & !is.null(Gen_pos)){
    #       ### use this step to drop the original GID and GID:Env
    #       rand_termCopy <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #
    #       rand_termCopy = gsub(" ", "", rand_termCopy)
    #
    #       #Inter_Gen_pos_copy = match(rand_term[Inter_Gen_pos], rand_termCopy)
    #       if(!is.null(heter_groups)){
    #         Inter_Gen_pos_copy = match(paste(gen_name,heter_groups, sep = ":"),rand_termCopy)
    #         if(anyNA(Inter_Gen_pos_copy)){
    #           Inter_Gen_pos_copy= match(paste(heter_groups,gen_name, sep = ":"),rand_termCopy)
    #         }
    #
    #       }
    #       Gen_pos_copy = match(rand_term[Gen_pos], rand_termCopy)
    #
    #       random = stats::formula(stats::drop.terms(stats::terms(random), c(Gen_pos_copy,Inter_Gen_pos_copy), keep.response = F))
    #
    #     }
    #
    #   } ## End of when variance-covariance str is not provided by user
    #
    #   #####################################################
    #   ####
    #   # Extract each random component/terms
    #   ####################################################
    #
    #   ranTerms <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    #   ranTerms = gsub(" ", "", ranTerms)
    #
    #   for (r in 1:length(ranTerms)) {
    #
    #     if(r == 1) {
    #       code.asr[2] <- paste(code.asr[2], ranTerms[r], sep='')
    #
    #     } else {
    #
    #       code.asr[2] <- paste(code.asr[2], ranTerms[r], sep='+')
    #     }
    #   }
    #
    #
    # } ## End of fixing the random terms

    # Heterogeneous errors
    if (!is.null(heter_groups)&isTRUE(heter_resid)) {

      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        #code.asr[3] <- paste(code.asr[3], paste0('dsum(~id(units)|', paste0(heter_groups, ')')), sep='')
        code_asr_fit[3] <- paste( code_asr_fit[3], paste0('dsum(~id(units)|', paste0(heter_groups, ')')), sep='')
      } else {
        warning(paste(msg,'Heterogenous residual is not possible with one environment/location. We fix it for you.'),
                call. = FALSE)
        if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
          #code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
          code_asr_fit[3] <- paste(code_asr_fit[3], 'id(units)', sep='')
        }
      }
      #code.asr[3] <- paste(code.asr[3], paste0('dsum(~idv(units)|', paste0(heter_groups, ')')), sep='')

      #code.asr[3] <- paste(code.asr[3], 'dsum(~idv(units|', 'heter_groups)', sep='')
    } else {
      #code.asr[3] <- paste(code.asr[3], 'idv(units)', sep='')
      #code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
      code_asr_fit[3] <- paste(code_asr_fit[3], 'id(units)', sep='')
    }
    # Adding individual/gen_name (random)
    #code.asr[2] <- paste(code.asr[2], 'vm(indiv,ainv)', sep='+')

    #Univariate = c()
    #for (trait in 1:length(response)) {

      #if (trait==1){
        #code.asr[1] <-  gsub("trait", response[trait], code.asr[1])
        #code.asr[1] <-  gsub("trait", response, code.asr[1])
   code_asr_fit[1] <-  gsub("trait", response, code_asr_fit[1])
      #} else {

        #code.asr[1] <-  gsub(response[trait-1], response[trait], code.asr[1])
      #}


      if(is.null(weights)){
        #code.asr[4] <- 'na.action=list(x="include",y="include"),data=pheno_data)'
        code_asr_fit[4] <- 'na.action=list(x="include",y="include"),data=pheno_data)'
      } else {
        if(!is.null(weights)){
          #code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_gaussian(dispersion = 1), data=pheno_data)'
          code_asr_fit[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_gaussian(dispersion = 1), data=pheno_data)'
        }

        # if(length(unique(response[trait]))<=10 & is.null(weights)){
        #   code.asr[4] <- 'na.action=list(x="include",y="include"), family = asr_multinomial(),  data=pheno_data)'
        # }
        #
        # if(length(unique(response[trait]))<=10 & !is.null(weights)){
        #   code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_multinomial(dispersion = 1),  data=pheno_data)'
        # }

        if(length(unique(response[trait]))==2 & !is.null(weights)){
          #code.asr[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_binomial(dispersion = 1),  data=pheno_data)'
          code_asr_fit[4] <- 'na.action=list(x="include",y="include"), weights = weights, family = asr_binomial(dispersion = 1),  data=pheno_data)'
        }

        if(length(unique(response[trait]))==2 & is.null(weights)){
          #code.asr[4] <- 'na.action=list(x="include",y="include"), family = asr_binomial(),  data=pheno_data)'
          code_asr_fit[4] <- 'na.action=list(x="include",y="include"), family = asr_binomial(),  data=pheno_data)'
        }

      }



      ####
      #code.asr[1] <- paste('mod<-', code.asr[1], sep='')
       code_asr_fit[1] <- paste('mod<-', code_asr_fit[1], sep='')
      str_mod <- paste(code_asr_fit[1],code_asr_fit[2],code_asr_fit[3],code_asr_fit[4],sep=',')
      ## Calls the current environment for evaluation
      eval(parse(text=str_mod), envir=environment())
      if (!mod$converge) { eval(parse(text='mod<-asreml::update.asreml(mod)')) }

      ###################################################
      ##### Start the process of processing the results
      ################################################


      #Univariate[[trait]] = mod

    #} ## End loop for multiple response variables

    #names(Univariate) <- response

    #names(Univariate)[trait] <- response[trait]
    #parallel::stopCluster(cl)

    #doParallel::stopImplicitCluster()

output <- list(mod,
               str_mod,
               G_list,
               Gen_pos,
               Inter_Gen_pos,
               rand_term)

names(output) <- c("model",
                   "str_mod",
                   "G_list",
                   "Gen_pos",
                   "Inter_Gen_pos",
                   "rand_term")
  #return(c(Univariate, G_list))
return(output)

}




