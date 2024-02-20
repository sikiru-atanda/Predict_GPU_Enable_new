#' Title
#'
#' @param response
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param var_cov_str
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
    var_cov_str = NULL,
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

  #if(!matrixcalc::is.positive.definite(gmatrix)) stop("GRM issue")
    ######
  # Function to compute inverse and create sparse matrix
  # Regularization parameter (adjust as needed)

  compute_inverse_and_sparse <- function(kernel = NULL,
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
      attr(sparse, "INVERSE") <- FALSE
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
                     var_cov_str = var_cov_str,
                     code_asr = code_asr,
                     G_list = G_list,
                     gen_name = gen_name,
                     pheno_data = pheno_data)

   code_asr_fit <-  rand_out[["code_asr"]]
   gen_pos <-  rand_out[["gen_pos"]]
   inter_gen_pos <-  rand_out[["inter_gen_pos"]]
   rand_term <- rand_out[["rand_term"]]
##### Ends with random term

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
               gen_pos,
               inter_gen_pos,
               rand_term)

names(output) <- c("model",
                   "str_mod",
                   "G_list",
                   "gen_pos",
                   "inter_gen_pos",
                   "rand_term")
  #return(c(Univariate, G_list))
return(output)

}




