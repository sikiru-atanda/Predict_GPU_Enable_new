

#' Title
#' This function calculate the relationship for using different kernel methods
#' M_matrix_clean is not expected with NA value
#' The following kernel methods are available:
#' Gaussian, Linear, Exponential, Poly2, Poly3, Poly4. The number at the end of poly
#' means polynomial order 2, 3 and 4
#'
#' To DO: in the new version, option to impute missing data will be provided
#'
#' @param M_matrix_clean snp or omics matrix data.
#' @param center
#' @param method
#' @param message
#'
#' @return
#' @export
#'
#' @examples
#'
kernel_calculation <- function(
    M_matrix_clean = NULL,
    center=TRUE,
    theta = NULL,
    # method = c("Gaussian",
    #            "Linear",
    #            "Poly2",
    #            "Poly3",
    #            "Poly4"),
    method = NULL,
    message = TRUE,
    ...){

  msg <- sprintf("==================================================\n")

  if(is.null(theta)) theta <- 1

  if(is.null(M_matrix_clean)){

    stop(print(paste(msg,' object omics/M_matrix is missing. We fix it')), call. = FALSE)


    }

  # if(all(class(M_matrix_clean)!= c("matrix", "array", "omic_matrix"))) {
  #   stop(print('Data is not class matrix'), call. = FALSE)
  # }

  if(class(M_matrix_clean)!= c("matrix")) {
    M_matrix_clean = as.matrix(M_matrix_clean)
    stop(print(paste(msg,' M_matrix is not class matrix. We fix it')), call. = FALSE)
  }

  if (any(is.na(M_matrix_clean))){
    stop(print(paste(msg,' Missing value is not expected.')), call. = FALSE)
    }

  if(isFALSE(center)){
    if(isTRUE(message)) message(paste(msg,"if data is not previously centered.It is recommend you center the data." ))
    }

  if(isTRUE(center)){

    M_matrix_clean = scale(x = M_matrix_clean,center = T,scale = F)
  }
  Gaussian_kernel <- function(M_matrix_clean, theta){

    dist<-as.matrix(stats::dist(M_matrix_clean))^2

    #GK<-exp(-dist/stats::median(dist))

     GK<-exp(-theta*dist/stats::median(dist))



    return(GK)
  }

  Exponential_kernel <- function(M_matrix_clean, theta){


    dist <- as.matrix(stats::dist(M_matrix_clean, method = "euclidian"))/sqrt(ncol(M_matrix_clean))


    EK <- exp(-theta*dist)

    #GK<-exp(-h*dist/stats::median(dist))

    return(EK)
  }

  Polynomial2_kernel = function(M_matrix_clean){

    ### matrix package drop
    PK2 <- ((Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean))+1)^2

  }

  Polynomial3_kernel = function(M_matrix_clean){

    PK3 <- ((Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean))+1)^3

  }

  Polynomial4_kernel = function(M_matrix_clean){

    PK4 <- ((Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean))+1)^4

  }


  Liner_kernel = function(M_matrix_clean){

    LK =  Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean)
  }

  switch(method,
         "Gaussian" = {
           KRM <- Gaussian_kernel(M_matrix_clean, theta)
         },
         "Linear" = {
           KRM <- Liner_kernel(M_matrix_clean)
         },
         "Exponential" = {
           KRM <- Exponential_kernel(M_matrix_clean, theta)
         },
         "Poly2" = {
           KRM <- Polynomial2_kernel(M_matrix_clean)
         },
         "Poly3" = {
           KRM <- Polynomial3_kernel(M_matrix_clean)
         },
         "Poly4" = {
           KRM <- Polynomial4_kernel(M_matrix_clean)
         },
         {
           stop(print(paste(msg,'Select method to calculate kernel ralationship matrix')), call. = FALSE)

         })


  return(KRM)


}
