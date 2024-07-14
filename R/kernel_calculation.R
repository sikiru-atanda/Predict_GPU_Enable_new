#' Calculate Kernel Matrix for Omics Data
#'
#' This function computes a kernel matrix for omics data using various methods, including Gaussian, Linear,
#' Composite, and Polynomial kernels. The function allows for scaling of the input matrix and customization of
#' kernel parameters.
#'
#' @param M_matrix_clean A numeric matrix representing the cleaned omics or M matrix.
#' @param scale Logical, indicating if the input matrix should be scaled. Defaults to TRUE.
#' @param theta Numeric, the theta parameter for the Gaussian and Exponential kernels. Defaults to 1 if not provided.
#' @param alpha Numeric, the mixing parameter for the Composite kernel. Defaults to 0.5.
#' @param gamma Numeric, the gamma parameter for the Anova radial basis kernel.
#' @param smoothness_parameter Numeric, the smoothness parameter for the Matérn kernel.
#' @param length_scale Numeric, the length scale parameter for the Matérn kernel.
#' @param method Character string specifying the kernel calculation method. Supported methods include "Gaussian_kernel",
#' "Linear_kernel", "Composite_kernel", "Poly2_kernel", "Poly3_kernel", and "Poly4_kernel".
#' @param message Logical, indicating if messages should be printed. Defaults to TRUE.
#' @param ... Additional arguments passed to the kernel calculation function.
#'
#' @return Returns a numeric matrix representing the calculated kernel matrix.
#'
#' @examples
#' # Example usage with a Gaussian kernel
#' M_matrix <- matrix(rnorm(100), ncol=10)
#' kernel_matrix <- kernel_calculation(M_matrix_clean = M_matrix, method = "Gaussian_kernel")
#'
#' @export
#'
kernel_calculation <- function(
    M_matrix_clean = NULL,
    scaling = FALSE,
    centering = TRUE,
    theta = NULL,
    alpha = 0.5,
    gamma = 1,
    smoothness_parameter = 1.5,
    length_scale = 2,
    method = NULL,
    message = TRUE,
    ...){

  msg <- "\n==================================================\n"

  if(is.null(theta)) theta <- 1

  if(is.null(M_matrix_clean)){

    stop(print(paste(msg,' object omics/M_matrix is missing. We fix it')), call. = FALSE)


    }

  kernel_method_avaliable <- c("Gaussian_kernel",
                               "Linear_kernel",
                               "Composite_kernel",
                               "Poly2_kernel",
                               "Poly3_kernel",
                               "Poly4_kernel")

  if(!is.null(method)){
    if (!(method %in% kernel_method_avaliable)) {
      stop("Invalid kernel method. Choose from: ",
           paste(kernel_method_avaliable, collapse = ", "), call. = FALSE)
    }
  }


  ## checking if the class attribute "matrix" is present in the vector of class attributes returned by class.
  if(!inherits(M_matrix_clean, "matrix")) {
    M_matrix_clean <-  as.matrix(M_matrix_clean)
    message(insight::print_color(paste(msg,paste("M_matrix is not class matrix. We fix it.")), "blue"))

  }

  if (any(is.na(M_matrix_clean))){
    stop(message(paste(msg,' Missing value is not expected.')), call. = FALSE)
    }

  ##Scaling: Prior to applying the Gaussian kernel,
  ## scaling (and possibly centering) the data might still be beneficial
  ## to ensure that all features contribute equally to the distance calculations.
  ## This is particularly true if the SNP data is combined with other omic data
  ## that might have different scales or units.
  if(isFALSE(scaling)){
    if(isTRUE(message)) message(insight::print_color(paste(msg,paste("If data is not previously scaled.It is recommend you scale the data.")), "blue"))
    }

  if(isTRUE(scaling)){

    M_matrix_clean = scale(x = M_matrix_clean,center = TRUE,scale = TRUE)
  }else {
    if(isTRUE(centering)){

      M_matrix_clean = scale(x = M_matrix_clean,center = TRUE,scale = FALSE)
    }
  }


  # if(isTRUE(scale) && !method%in%(c("Normalized_laplacian_kernel",
  #                                    "spectral_kernel_matrix"
  #                                    #"Matern_kernel_matrix"
  #                                    ))){
  #
  #   M_matrix_clean = scale(x = M_matrix_clean,center = FALSE,scale = TRUE)
  # }
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

    return(PK2)

  }

  Polynomial3_kernel = function(M_matrix_clean){

    PK3 <- ((Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean))+1)^3

    return(PK3)

  }

  Polynomial4_kernel = function(M_matrix_clean){

    PK4 <- ((Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean))+1)^4

    return(PK4)

  }


  Linear_kernel = function(M_matrix_clean){

    LK =  Matrix::tcrossprod(M_matrix_clean)/ncol(M_matrix_clean)

    return(LK)
  }

  # # Function to compute Composite kernel
  # Composite_kernel <- function(kernel_list, weights) {
  #   # Check if the number of kernels matches the number of weights
  #   if (length(kernel_list) != length(weights)) {
  #     stop("Number of kernels and weights must match")
  #   }
  #
  #   # Initialize an empty kernel matrix
  #   composite_matrix <- NULL
  #
  #   # Iterate through each kernel and weight, and combine them
  #   for (i in seq_along(kernel_list)) {
  #     kernel <- kernel_list[[i]]
  #     weight <- weights[i]
  #
  #     # If it's the first kernel, initialize the composite matrix
  #     if (is.null(composite_matrix)) {
  #       composite_matrix <- kernel * weight
  #     } else {
  #       # Add the weighted kernel to the composite matrix
  #       composite_matrix <- composite_matrix + kernel * weight
  #     }
  #   }
  #
  #   return(composite_matrix)
  # }
  Composite_kernel <- function(M_matrix_clean, theta, alpha = 0.5) {
    linear_term <- Linear_kernel(M_matrix_clean)
    gaussian_term <- Gaussian_kernel(M_matrix_clean, theta)
    CK <- (alpha * linear_term + (1 - alpha) * gaussian_term)

    return(CK)
  }
  #####
  # Function to compute Anova radial basis kernel matrix
  # Anova_radial_basis_kernel <- function(M_matrix_clean, gamma) {
  #   n <- nrow(M_matrix_clean)
  #
  #   # Initialize the kernel matrix
  #   anova_rbf_kernel_matrix <- matrix(0, n, n)
  #
  #   # Compute the Anova radial basis kernel matrix
  #   for (i in 1:n) {
  #     for (j in 1:n) {
  #       diff_squared <- sum((M_matrix_clean[i, ] - M_matrix_clean[j, ])^2)
  #       anova_rbf_kernel_matrix[i, j] <- exp(-gamma * diff_squared)
  #     }
  #   }
  #
  #   return(anova_rbf_kernel_matrix)
  # }
  # ###
  # matern_kernel <- function(r, nu, rho) {
  #   term1 <- 2^(1 - nu) / gamma(nu)
  #   term2 <- (sqrt(2 * nu) * r / rho)^nu
  #   term3 <- besselK(sqrt(2 * nu) * r / rho, nu)
  #
  #   return(term1 * term2 * term3)
  # }
  #
  # # Function to compute Matérn kernel matrix
  # Matern_kernel_matrix <- function(M_matrix_clean,
  #                                  smoothness_parameter,
  #                                  length_scale
  #                                  #diag_value = 1e-6
  #                                  ) {
  #   n <- nrow(M_matrix_clean)
  #   matern_matrix <- matrix(0, n, n)
  #
  #   for (i in 1:n) {
  #     for (j in 1:n) {
  #       distance <- sqrt(sum((M_matrix_clean[i,] - M_matrix_clean[j,])^2))
  #
  #       # Handle the case when distance is zero (diagonal)
  #       if (i == j) {
  #         #matern_matrix[i, j] <- diag_value
  #         matern_matrix[i, j] <- 1.00
  #       } else {
  #         matern_matrix[i, j] <- matern_kernel(r = distance,
  #                                              nu = smoothness_parameter,
  #                                              rho = length_scale)
  #       }
  #     }
  #   }
  #
  #   return(matern_matrix)
  # }
  #
  # ####
  # Spectral_kernel_matrix <- function(M_matrix_clean) {
  #   n <- nrow(M_matrix_clean)
  #   spectral_matrix <- matrix(0, n, n)
  #
  #   for (i in 1:n) {
  #     for (j in 1:n) {
  #       # Compute spectral similarity between SNP profiles
  #       similarity <- sum(M_matrix_clean[i,] * M_matrix_clean[j,]) / (sqrt(sum(M_matrix_clean[i,]^2)) * sqrt(sum(M_matrix_clean[j,]^2)))
  #
  #       # Set the spectral matrix element
  #       spectral_matrix[i, j] <- similarity
  #     }
  #   }
  #
  #   return(spectral_matrix)
  # }
  # ###########
  # ###
  # ## Forming a graph matrix from SNP data involves defining
  # # relationships between individuals based on genetic similarity
  # # or other relevant criteria.
  # jaccard_similarity <- function(x, y) {
  #   intersection <- sum(x & y)
  #   union <- sum(x | y)
  #   return(intersection / union)
  # }
  # ##
  # hamming_distance <- function(x, y) {
  #   return(sum(x != y))
  # }
  # #
  # euclidean_distance <- function(x, y) {
  #   return(sqrt(sum((x - y)^2)))
  # }
  # #
  # cosine_similarity <- function(x, y) {
  #   return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
  # }
  # ##
  # sokal_michener_similarity <- function(x, y) {
  #   a <- sum(x & y)
  #   b <- sum(x & !y)
  #   c <- sum(!x & y)
  #   d <- sum(!x & !y)
  #   return((a + d) / (a + b + c + d))
  # }
  # ##
  # dice_similarity <- function(x, y) {
  #   intersection <- sum(x & y)
  #   union <- sum(x) + sum(y)
  #   return(2 * intersection / union)
  # }
  # ###
  # # Function to construct a graph matrix based on Jaccard similarity
  # construct_graph_matrix <- function(M_matrix_clean, similarity_method = "jaccard_similarity") {
  #   n <- nrow(M_matrix_clean)
  #   graph_matrix <- matrix(0, n, n)
  #
  #   for (i in 1:n) {
  #     for (j in 1:n) {
  #       if (i != j) {
  #         if(similarity_method =="jaccard_similarity"){
  #           similarity <- jaccard_similarity(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #
  #         } else if (similarity_method =="hamming_distance"){
  #            similarity <- hamming_distance(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #         } else if (similarity_method =="euclidean_distance"){
  #           similarity <- euclidean_distance(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #         } else if (similarity_method =="cosine_similarity"){
  #           similarity <- cosine_similarity(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #         } else if (similarity_method =="sokal_michener_similarity"){
  #           similarity <- sokal_michener_similarity(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #
  #         } else {
  #           if(similarity_method =="dice_similarity"){
  #             similarity <- dice_similarity(M_matrix_clean[i, ], M_matrix_clean[j, ])
  #           }
  #         }
  #         graph_matrix[i, j] <- similarity
  #       }
  #     }
  #   }
  #
  #   return(graph_matrix)
  # }
  # ##
  # #compute symmetric normalized Laplacian kernel matrix
  # # Function to compute symmetric normalized Laplacian kernel matrix
  # Normalized_laplacian_kernel <- function(M_matrix_clean) {
  #
  #   graph_matrix <- construct_graph_matrix(M_matrix_clean)
  #   n <- nrow(graph_matrix)
  #
  #   # Compute the degree matrix
  #   degree_matrix <- diag(rowSums(graph_matrix))
  #
  #   # Compute the Laplacian matrix
  #   laplacian_matrix <- degree_matrix - graph_matrix
  #
  #   # Compute the symmetric normalized Laplacian kernel matrix using the pseudo-inverse
  #   symmetric_normalized_laplacian_kernel <- MASS::ginv(sqrt(degree_matrix)) %*% laplacian_matrix %*% MASS::ginv(sqrt(degree_matrix))
  #
  #   return(symmetric_normalized_laplacian_kernel)
  # }

  switch(method,
         "Gaussian_kernel" = {
           KRM <- Gaussian_kernel(M_matrix_clean, theta)
         },
         "Linear_kernel" = {
           KRM <- Linear_kernel(M_matrix_clean)
         },
         "Composite_kernel" = {
           KRM <- composite_kernel(M_matrix_clean, theta, alpha)
         },
         # "Anova_radial_basis_kernel" = {
         #   KRM <- Anova_radial_basis_kernel(M_matrix_clean, gamma)
         # },
         # "Exponential_kernel" = {
         #   KRM <- Exponential_kernel(M_matrix_clean, theta)
         # },
         # "Matern_kernel" = {
         #   KRM <- Matern_kernel_matrix(M_matrix_clean,
         #                               smoothness_parameter,
         #                               length_scale
         #                               )
         #
         # },
         # "spectral_kernel" = {
         #     KRM <- Spectral_kernel_matrix(M_matrix_clean)
         #   },
         # "Normalized_laplacian_kernel" = {
         #
         #   KRM <- Normalized_laplacian_kernel(M_matrix_clean)
         # },
         "Poly2_kernel" = {
           KRM <- Polynomial2_kernel(M_matrix_clean)
         },
         "Poly3_kernel" = {
           KRM <- Polynomial3_kernel(M_matrix_clean)
         },
         "Poly4_kernel" = {
           KRM <- Polynomial4_kernel(M_matrix_clean)
         },
         {
           stop(message(paste(msg,'Select method to calculate kernel ralationship matrix')), call. = FALSE)

         })


  return(KRM)


}
