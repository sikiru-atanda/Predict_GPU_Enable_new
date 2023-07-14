

#' Title
#'
#' @param M_matrix_clean
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
    # method = c("Gaussian",
    #            "Linear",
    #            "Poly2",
    #            "Poly3",
    #            "Poly4"),
    method = NULL,
    message = TRUE){

  msg <- sprintf("==================================================\n")

  if(is.null(M_matrix_clean)){stop(print('object omics is missing'), call. = FALSE)}

  if(all(class(M_matrix_clean)!= c("matrix", "array", "omic_matrix"))) {
    stop(print('Data is not class matrix'), call. = FALSE)
  }
  if (any(is.na(M_matrix_clean))){
    stop(print(paste(msg,' Missing value is not expected.')), call. = FALSE)
    }

  if(isFALSE(center)){
    if(message) message(paste(msg,"if data is not previously centered.It is recommend you center the data." ))
    }

  if(center){

    M_matrix_clean = scale(x = M_matrix_clean,center = T,scale = F)
  }
  Gaussian_kernel <- function(M_matrix_clean){

    dist<-as.matrix(stats::dist(M_matrix_clean))^2

    GK<-exp(-dist/stats::median(dist))

    #GK<-exp(-h*dist/stats::median(dist))

    return(GK)
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
           KRM <- Gaussian_kernel(M_matrix_clean)
         },
         "Linear" = {
           KRM <- Liner_kernel(M_matrix_clean)
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
