

#' Title
#'
#' @param omics.data
#' @param center
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#'
kernelMatrix <- function(
    omics.data,
    center=TRUE,
    method = c("Gaussian",
               "Linear",
               "Poly2",
               "Poly3",
               "Poly4")){

  msg <- sprintf("==================================================\n")

  if(is.null(omics.data)){stop(print('object omics is missing'), call. = FALSE)}

  if(class(omics.data)[1]!= "matrix") {
      omics.data = as.matrix(omics.data)
    }
  if (any(is.na(omics.data))){stop(print('object omics is contain missing data'), call. = FALSE)}

  if(center==FALSE){warning('if omics is not previously centered.It is recommend you center the data')}

  if(center== TRUE){

    omics.data = scale(x = omics.data,center = T,scale = F)
  }
  Gaussian_kernel <- function(omics.data){

    dist<-as.matrix(stats::dist(omics.data))^2

    GK<-exp(-dist/stats::median(dist))

    #GK<-exp(-h*dist/stats::median(dist))

    return(GK)
  }

  Polynomial2_kernel = function(omics.data){

    ### matrix package drop
    PK2 <- ((Matrix::tcrossprod(omics.data)/ncol(omics.data))+1)^2

  }

  Polynomial3_kernel = function(omics.data){

    PK3 <- ((Matrix::tcrossprod(omics.data)/ncol(omics.data))+1)^3

  }

  Polynomial4_kernel = function(omics.data){

    PK4 <- ((Matrix::tcrossprod(omics.data)/ncol(omics.data))+1)^4

  }


  Liner_kernel = function(omics.data){

   LK =  Matrix::tcrossprod(omics.data)/ncol(omics.data)
  }

  switch(method,
         "Gaussian" = {
           KRM <- Gaussian_kernel(omics.data)
         },
         "Linear" = {
           KRM <- Liner_kernel(omics.data)
         },
         "Poly2" = {
           KRM <- Polynomial2_kernel(omics.data)
         },
         "Poly3" = {
           KRM <- Polynomial3_kernel(omics.data)
         },
         "Poly4" = {
           KRM <- Polynomial4_kernel(omics.data)
         },
         {
           stop(print("Select method to calculate kernel ralationship matrix"), call. = FALSE)
         })


  return(KRM)


}
