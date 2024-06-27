#' Compute Inverse and Sparse Matrix from a Kernel Matrix
#'
#' This function computes the inverse of a given kernel matrix using a regularization
#' technique to handle nearly positive definite matrices with small negative eigenvalues.
#' It can also return a sparse representation of the kernel matrix.
#'
#' @param kernel A numeric square matrix representing the kernel from which to compute the inverse.
#' @param epsilon A small positive value (regularization parameter) added to the diagonal elements
#'        to ensure the matrix is positive definite. Default is 1e-6. adjust as needed
#' @param inverse Logical indicating whether to compute and return the inverse of the kernel matrix.
#'        If FALSE, returns a sparse matrix representation of the original kernel matrix. Default is NULL.
#'
#' @return If `inverse = TRUE`, returns a matrix which is the inverse of the input kernel matrix,
#'         with additional attributes "rowNames", "colNames", and "INVERSE" set to TRUE.
#'         If `inverse = FALSE`, returns a sparse matrix representation of the kernel matrix,
#'         with "rowNames", "colNames", and "INVERSE" set to FALSE. If an error occurs during
#'         inversion, it returns `NULL` and prints an error message.
#'
#' @examples
#' # Create a symmetric positive definite matrix
#' kernel_matrix <- matrix(c(2, -1, -1, 2), ncol = 2)
#' rownames(kernel_matrix) <- colnames(kernel_matrix) <- c("A", "B")
#'
#' # Compute the inverse
#' inverse_matrix <- compute_inverse_and_sparse(kernel = kernel_matrix, inverse = TRUE)
#'
#' # Compute a sparse representation
#' sparse_matrix <- compute_inverse_and_sparse(kernel = kernel_matrix, inverse = FALSE)
#'
#' @export

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
    #if(inherits(result, "try-error")) {
    if(is.null(result)) {
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

# Function to check if current order of rownames or colnames is the same as unique_GIDs
## in pheno_data
# is_in_correct_order <- function(names, order) {
#   if (length(names) != length(order)) {
#     return(FALSE)  # Different lengths mean they are not in the same order
#   }
#   all(names == order)
# }

#' Adjust Random Terms for ASReml Model Fitting
#'
#' This function prepares and adjusts random terms for fitting an ASReml model,
#' facilitating the modeling of genetic data with potential heterogeneity and
#' specific variance-covariance structures. It allows for the specification of
#' fixed and random effects, interaction terms, and various models of heterogeneity
#' and variance-covariance among genetic components.
#'
#' @param random A formula specifying the random effects to be included in the model.
#' @param fixed A formula specifying the fixed effects to be included in the model.
#' @param fixed_term A character vector specifying fixed terms to be considered in the model.
#' @param heter_groups A character string identifying the grouping variable for heterogeneity.
#' @param heter_resid A logical indicating if heterogeneity in residuals should be considered.
#' @param var_cov_str A character string indicating the type of variance-covariance structure.
#' @param code_asr A character vector for storing ASReml code or model specification.
#' @param names_in_inv_list A character vector of names identifying inverse matrices in the model,
#' related to different genetic components.
#' @param gen_name A character string specifying the name of the genetic factor in the model.
#' @param pheno_data A data frame containing phenotypic data used in the model.
#' @param ... Additional arguments for future use or extensions.
#'
#' @return A list containing elements critical for ASReml model specification,
#' including adjusted code for model fitting, positions of genetic and interaction terms,
#' and the defined random terms.
#'
#' @details
#' The function is particularly useful for complex genetic models requiring precise
#' specification of random effects, handling of heterogeneity across groups, and
#' incorporation of specific variance-covariance structures. It preprocesses the input
#' model terms to ensure compatibility with ASReml software requirements and facilitates
#' the inclusion of genetic interactions and heterogeneity terms.
#'
#' @examples
#' # Assuming a hypothetical ASReml model setup:
#' random_effects <- ~ Genotype + Genotype:Environment
#' fixed_effects <- ~ Environment + Treatment
#' pheno <- data.frame(Genotype = factor(rep(1:10, each = 6)),
#'                     Environment = factor(rep(1:3, times = 20)),
#'                     Treatment = factor(rep(c("Control", "Treated"), each = 30)),
#'                     Yield = rnorm(60, mean = 100, sd = 15))
#'
#' model_setup <- random_terms_fit_new(random = random_effects,
#'                                     fixed = fixed_effects,
#'                                     heter_groups = "Environment",
#'                                     gen_name = "Genotype",
#'                                     pheno_data = pheno)
#' @export

asreml_utilis_new <- function(
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
    workspace=1e08,
    engine = NULL,
    pworkspace= 1e06,
    maxit = 50,
    cross_validation = FALSE,
    ...
) {

  msg <- sprintf("==================================================\n")

  if(engine %in% rownames(installed.packages())){
    do.call('library', list(engine))

    # if package is not installed locally then stop
  } else {

    stop(print(paste(msg,'You need to install asreml-R to use asreml-R')), call. = FALSE)
  }

  datasets <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
  dataset_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
  datasets_index <- which(!sapply(datasets, is.null))
  datasets <-  datasets[datasets_index]
  dataset_names <- dataset_names[datasets_index]

  # Extract unique GIDs from pheno_data to determine the row order
  unique_GIDs <- as.character(unique(pheno_data[[gen_name]]))

  # # Reorder the rownames and colnames of each dataset based on unique_GIDs if necessary
  # datasets <- lapply(datasets, function(mat) {
  #   correct_row_order <- is_in_correct_order(rownames(mat), unique_GIDs[unique_GIDs %in% rownames(mat)])
  #   correct_col_order <- is_in_correct_order(colnames(mat), unique_GIDs[unique_GIDs %in% colnames(mat)])
  #
  #   if (!correct_row_order || !correct_col_order) {
  #     # Reorder rows and columns if either is not in the correct order
  #     ordered_indices <- unique_GIDs[unique_GIDs %in% rownames(mat)]
  #     mat <- mat[ordered_indices, ordered_indices]
  #     message("Matrix reordered based on unique_GIDs.")
  #   } else {
  #     message("Matrix is already in the correct order; no changes made.")
  #   }
  #   return(mat)
  # })


  inv_list <- list()
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]

    if (!is.null(dataset)) {
      inv_list[[dataset_names[i]]] <- compute_inverse_and_sparse(dataset, inverse)
    }


  }

  if(!"gmatrix" %in%names(inv_list)){
    if(length(inv_list)>1){
      names(inv_list) <- paste(paste("omic", seq_along(inv_list), sep = ""), "inv", sep = "_")
    } else {
      names(inv_list) <- paste("omic", "inv", sep = "_")
    }
  } else if(length(grep("omic", names(inv_list)))==0 & "gmatrix" %in%names(inv_list)) {
    geno_index <- grep("gmatrix", names(inv_list))
    names(inv_list)[geno_index] <- paste("gmatrix", "inv", sep = "_")
  } else{
    if(length(grep("omic", names(inv_list)))>=1 & "gmatrix" %in%names(inv_list)) {
      omic_index <- grep("omic", names(inv_list))
      geno_index <- grep("gmatrix", names(inv_list))
      if(length(omic_index)>1){
        names(inv_list)[omic_index] <- paste(paste("omics", seq_along(omic_index), sep = ""), "inv", sep = "_")
      } else {
        names(inv_list)[omic_index] <- paste("omics", "inv", sep = "_")
      }
      ####
      names(inv_list)[geno_index] <- paste("gmatrix", "inv", sep = "_")
    }

  }

  names_in_inv_list <- names(inv_list)
  ## Here pworkspace is not included because predict function is not done here
  # asreml::asreml.options(trace=FALSE,
  #                        workspace = workspace,
  #                        pworkspace = pworkspace,
  #                        maxit = maxit)


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
  for (i in seq_along(inv_list)) {
    assign(names(inv_list)[i], inv_list[[i]], envir = .GlobalEnv)
  }

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

  rand_out <-  random_terms_fit_new(random = random,
                                fixed = fixed,
                                fixed_term = fixed_term,
                                heter_groups = heter_groups,
                                heter_resid = heter_resid,
                                var_cov_str = var_cov_str,
                                code_asr = code_asr,
                                names_in_inv_list = names_in_inv_list,
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
  # Adding traits
  #Univariate = c()
  #for (trait in 1:length(response)) {

  #if (trait==1){
  #code.asr[1] <-  gsub("trait", response[trait], code.asr[1])
  #code.asr[1] <-  gsub("trait", response, code.asr[1])
  code_asr_fit[1] <-  gsub("trait", response, code_asr_fit[1])
  #} else {

  #code.asr[1] <-  gsub(response[trait-1], response[trait], code.asr[1])
  #}
  if(isTRUE(cross_validation)){

    output <- list(code_asr_fit = code_asr_fit,
                   names_in_inv_list = names_in_inv_list ,
                   gen_pos = gen_pos,
                   inter_gen_pos = inter_gen_pos,
                   rand_term = rand_term,
                   inv_list = inv_list)


    return(output)
  }

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


  asreml::asreml.options(trace=FALSE)
  ####
  #code.asr[1] <- paste('mod<-', code.asr[1], sep='')
  code_asr_fit[1] <- paste('mod<-', code_asr_fit[1], sep='')
  str_mod <- paste(code_asr_fit[1],code_asr_fit[2],code_asr_fit[3],code_asr_fit[4],sep=',')
  if(isTRUE(cross_validation)){
    return(str_mod)
  }
  ## Calls the current environment for evaluation
  eval(parse(text=str_mod), envir=environment())
  #if (!mod$converge) { eval(parse(text='mod<-asreml::update.asreml(mod)')) }
  # Assuming `mod` is your initial model object
  for (i in 1:3) {
    if (!mod$converge) {

      eval(parse(text='mod<-asreml::update.asreml(mod)'))
      }else {
      # Exit the loop if model has converged
      break
    }
  }


  ###### Process if the model is not stable #######
  # specific_warning_occurred <- FALSE
  # error_occurred <- FALSE
  #
  # repeat {
  #   # Attempt to update the model and capture warnings
  #   tryCatch({
  #     # Assuming 'res' is your model object and update.asreml is the function you're using
  #     mod <- asreml::update.asreml(mod)
  #
  #     # If update.asreml runs without warnings or errors, we assume the update was successful
  #   }, warning = function(w) {
  #     # Check if the warning message matches the specific warning you're concerned with
  #     if(grepl("Some components changed by more than 1% on the last iteration", w$message)) {
  #       specific_warning_occurred <- TRUE
  #       # Log the occurrence of the specific warning for debugging
  #       cat("Specific warning occurred, attempting to update the model again...\n")
  #     }
  #   }, error = function(e) {
  #     error_occurred <- TRUE
  #     # Log the error for debugging
  #     cat("An error occurred: ", e$message, "\n")
  #   })
  #
  #   # Break the loop if an error occurred or if the specific warning did not occur in this iteration
  #   if (error_occurred || !specific_warning_occurred) {
  #     break
  #   }
  #
  #   # Reset the specific warning flag for the next iteration
  #   specific_warning_occurred <- FALSE
  # }


  ###################################################
  ##### Start the process of processing the results
  ################################################

  output <- list(mod,
                 str_mod,
                 names_in_inv_list ,
                 gen_pos,
                 inter_gen_pos,
                 rand_term)

  names(output) <- c("model",
                     "str_mod",
                     "names_in_inv_list",
                     "gen_pos",
                     "inter_gen_pos",
                     "rand_term")
  #return(c(Univariate, G_list))
  return(output)

}




