
#' Title
#' grm_kernel_data is a genomic relationship matrix or relationship matrix calculated using
#' other omics data. The relationship matrix can be calculated using any method.
#' NA is not allowed in the relationship matrix
#' ####
#' Checks
#' #######
#' 1. It check if the matrix is square matrix/symmetry, if not we fix it for the user
#' 2. It check if the matrix is positive definite, if not we fix it.
#' 3. It check for NA. If present the engine will stop further analysis.
#'
#'
#'
#' @param grm_kernel_data
#' @param bending
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param message
#' @param pedigree_matrix
#' @param rcn_cutoff  #the reciprocal conditional number of the inverse
#' @param ...
grm_kernel_precheck <- function(grm_kernel_data= NULL,
                                pedigree_matrix = NULL,
                                bending = TRUE,
                                bend_value = 0.01,
                                blending = FALSE,
                                blending_value = 0.02,
                                high_diag_cut_off = 1.2,
                                low_diag_cut_off = 0.8,
                                duplicate_cut_off = 0.95,
                                rcn_cutoff = 1e-12,
                                optimize_diagonal = FALSE,
                                optimize_duplicate = FALSE,
                                message= TRUE,
                                ...){

  msg <- sprintf("==================================================\n")

  ## Check no NA is present in the matrix
if (!is.null(grm_kernel_data)){
  if(isTRUE(anyNA(grm_kernel_data))){
    stop(message(paste(msg,'NA is not allowed in the grm or kernel matrix')), call. = FALSE)
  }
  ## Check rownames is provided
  if (is.null(rownames(grm_kernel_data))){
    stop(message(paste(msg,'Rownames containing individuals in the matrix is missing')), call. = FALSE)
  }
  ## Check colnames is provided
  if (is.null(colnames(grm_kernel_data))){
    stop(message(paste(msg,'Colnames containing individuals in the matrix is missing')), call. = FALSE)
  }

  ### Check if the matrix is in class matrix if not convert to class matrix
  if (!is.matrix(grm_kernel_data)) grm_kernel_data <- as.matrix(grm_kernel_data)
  ### Check if colname and rownames in grm/kernel matrix is the same
  if (!identical(colnames(grm_kernel_data), rownames(grm_kernel_data))) {stop(print(paste(msg,'colnames did not match rownames')), call. = FALSE)}
  ### Check if the grm/kernel matrix is symmetric


  #if(!isSymmetric.matrix(grm_kernel_data)) {stop(print(paste(msg,'grm_kernel_data is not symmetric')), call. = FALSE)}
  if(!isSymmetric.matrix(grm_kernel_data)) {
    message(insight::print_color(paste(msg,paste("Relationsip Matrix is not symmetric'. We fix it.")), "blue"))

    grm_kernel_data <- Matrix::forceSymmetric(grm_kernel_data)
    grm_kernel_data <- Matrix::as.matrix(grm_kernel_data)
  }

  if(isTRUE(bending) & !is.null(bend_value)){
  if(isFALSE(matrixcalc::is.positive.definite(grm_kernel_data))){
    message(insight::print_color(paste(msg,paste("Relationsip Matrix is not positive definite. We fix it.")), "blue"))

    grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data, posd.tol= bend_value, trace=FALSE)$mat)

  }

  } else {
#### If bending is False check for user to see the matrix is not ill-conditioned for model fit
    if(isFALSE(matrixcalc::is.positive.definite(grm_kernel_data))){

      grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data, posd.tol= bend_value, trace=FALSE)$mat)
      #message(paste(msg,"Relationsip Matrix is not positive definite. Set bending = TRUE to fix it"))
      message(insight::print_color(paste(msg,paste("Relationsip Matrix is not positive definite.\n \t We fix it by bending to make the matrix stable.")), "blue"))

    }

  }

  ### This is important to check even if the user defined blending as FALSE
  if(isFALSE(blending)){
  res = grm_kernel_diagnostic_fix(grm_kernel_data = grm_kernel_data,
                                  high_diag_cut_off = high_diag_cut_off,
                                  low_diag_cut_off = low_diag_cut_off,
                                  duplicate_cut_off = duplicate_cut_off,
                                  optimize_diagonal = optimize_diagonal,
                                  optimize_duplicate = optimize_duplicate
  )

### Also implore the reciprocal conditional number as metric to decide
  ## ill-conditioned/unstable matrix
  rcn <- rcond(res$clean_matrix)
    #if("potential_off_diag_with_duplicate"%in%names(res)){
    if("potential_off_diag_with_duplicate"%in%names(res) | rcn < rcn_cutoff){
      grm_kernel_data <- res$clean_matrix
      ncol_nrow <- ncol(grm_kernel_data)
      grm_kernel_data_ <- (1-blending_value)*grm_kernel_data + blending_value*diag(x=1, nrow=ncol_nrow , ncol=ncol_nrow )

      ## Repeat the process to ascertain the matrix is stable with no duplicate
      res <- grm_kernel_diagnostic_fix(grm_kernel_data = grm_kernel_data_,
                                       high_diag_cut_off = high_diag_cut_off,
                                       low_diag_cut_off = low_diag_cut_off,
                                       duplicate_cut_off = duplicate_cut_off,
                                       optimize_diagonal = optimize_diagonal,
                                       optimize_duplicate = optimize_duplicate)

      rcn <- rcond(grm_kernel_data_)

      ## Check if the matrix is still unstable. Call the attention of the user to provide
      ## another blending_value value.
      if("potential_off_diag_with_duplicate"%in%names(res) | rcn < rcn_cutoff){
        #if("potential_off_diag_with_duplicate"%in%names(res)){
        message(paste(insight::print_color("WARNINGS\n", "blue"),
                      insight::print_color(paste(msg,paste("Matrix contain duplicate(s) or still ill-conditioned which might be potential problem.\n \t Change the blending value eg. 0.05  etc.")), "blue")))

        ncol_nrow <-  ncol(grm_kernel_data)
        grm_kernel_data <- (1-blending_value)*grm_kernel_data_ + blending_value*diag(x=1, nrow=ncol_nrow , ncol=ncol_nrow )


      } else {
        if(isTRUE(message)){
          message(paste( insight::print_color("WARNINGS\n", "blue"),
                         insight::print_color(paste(msg,paste("Matrix contain duplicate(s) which might be potential problem.\n \t We fix it by blending using an identity matrix.")), "blue")))

        }
      }


    }

  } else {
 if(isTRUE(blending) & !is.null(blending_value)){
    grm_kernel_data <- blending_stat(grm_kernel_data = grm_kernel_data,
                                     pedigree_matrix = pedigree_matrix,
                                     blending = blending,
                                     blending_value = blending_value
    )

    ## Repeat the process to ascertain the matrix is stable with no duplicate
    res <-  grm_kernel_diagnostic_fix(grm_kernel_data = grm_kernel_data,
                                      high_diag_cut_off = high_diag_cut_off,
                                      low_diag_cut_off = low_diag_cut_off,
                                      duplicate_cut_off = duplicate_cut_off,
                                      optimize_diagonal = optimize_diagonal,
                                      optimize_duplicate = optimize_duplicate)

    ## Check if the matrix is still unstable. Call the attention of the user to provide
    ## another blending_value value.
    if("potential_off_diag_with_duplicate"%in%names(res)){
      message(paste( insight::print_color("WARNINGS\n", "blue"),
                     insight::print_color(paste(msg,paste("Matrix contain duplicate(s) which might be potential problem.\n \t Change the blending value eg. 0.05  etc.")), "blue")))


    }


  }

}

}

  if(exists("res")){
  rm(res)

  }

  #### Declare it also as an grm_kernel_data for final usage
  #class(grm_kernel_data) <-c("matrix", "array", "krm_data")

  attr(grm_kernel_data, "cleared") <- "for_model_fit"
  ##############################################
  return(grm_kernel_data)

}






