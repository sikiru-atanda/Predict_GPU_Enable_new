
# Using a kernel matrix derived from genomic data as a
# relationship matrix is a way to incorporate non-linear
# genetic relationships into the GBLUP model. The Matérn kernel,
# being a flexible and smooth kernel, can capture complex
# patterns and dependencies in the genetic data.
# This flexibility might be particularly useful if there
# are non-linear relationships or interactions among the
# genetic markers that influence the traits of interest.
##### matern_kernel

# r is the Euclidean distance between two points.
# nu is the smoothness parameter of the Matérn kernel.
# rho is the length scale parameter.

matern_kernel <- function(r, nu, rho) {
  term1 <- 2^(1 - nu) / gamma(nu)
  term2 <- (sqrt(2 * nu) * r / rho)^nu
  term3 <- besselK(sqrt(2 * nu) * r / rho, nu)

  return(term1 * term2 * term3)
}

# Function to compute Matérn kernel matrix
# Function to compute Matérn kernel matrix
compute_matern_kernel_matrix <- function(data, nu, rho, diag_value = 1e-6) {
  n <- nrow(data)
  matern_matrix <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      distance <- sqrt(sum((data[i,] - data[j,])^2))

      # Handle the case when distance is zero (diagonal)
      if (i == j) {
        matern_matrix[i, j] <- diag_value
      } else {
        matern_matrix[i, j] <- matern_kernel(distance, nu, rho)
      }
    }
  }

  return(matern_matrix)
}

######################

# Example usage
set.seed(123)
X_data <- matrix(rnorm(5000), ncol = 100, nrow = 50)  # Example data

# Parameters
smoothness_parameter <- 1.5
length_scale <- 2.0

# Compute Matérn kernel matrix
matern_matrix <- compute_matern_kernel_matrix(X_data, smoothness_parameter, length_scale)

# Display the resulting Matérn kernel matrix
print(matern_matrix)


# Example usage
set.seed(123)
X_data <- matrix(rnorm(5000), ncol = 100, nrow = 50)  # Example data

# Parameters
smoothness_parameter <- 1.5
length_scale <- 2.0

# Compute Matérn kernel matrix
matern_matrix <- compute_matern_kernel_matrix(X_data, smoothness_parameter, length_scale)

# Display the resulting Matérn kernel matrix
print(matern_matrix)

matern_matrix[1:5, 1:5]
##############################
# Creating a spectral kernel matrix involves transforming
# SNP data into a matrix of pairwise spectral similarities
# In this example:
#
# snp_data is a binary matrix representing SNP data,
# where each row corresponds to an individual,
# and each column corresponds to a SNP.
#
# The compute_spectral_kernel_matrix function computes
# the pairwise spectral similarity between SNP profiles
# and constructs the spectral matrix.
#
# The similarity between two SNP profiles is calculated as
# the dot product of the two profiles normalized by the
# product of their Euclidean norms.
#
# This simple example assumes binary SNP data, and you may
# need to adjust it based on the nature of
# your SNP data (e.g., if your SNP data is not binary).
# Additionally, you might want to explore more advanced
# techniques for spectral kernel computation,
# depending on the characteristics of your data and
# the requirements of your analysis.It's important to note that
# the spectral kernel is just one approach, and there are
# various other kernel methods that might be more suitable
# for genomic prediction tasks. The choice of the kernel
# should be guided by the properties of your data and the
# characteristics of the underlying genetic architecture.

# Function to compute Spectral kernel matrix
compute_spectral_kernel_matrix <- function(snp_data) {
  n <- nrow(snp_data)
  spectral_matrix <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      # Compute spectral similarity between SNP profiles
      similarity <- sum(snp_data[i,] * snp_data[j,]) / (sqrt(sum(snp_data[i,]^2)) * sqrt(sum(snp_data[j,]^2)))

      # Set the spectral matrix element
      spectral_matrix[i, j] <- similarity
    }
  }

  return(spectral_matrix)
}

# Example usage
set.seed(123)
snp_data <- matrix(rbinom(5000, 2, 0.5), ncol = 100, nrow = 50)  # Example SNP data (binary)

# Compute Spectral kernel matrix
spectral_matrix <- compute_spectral_kernel_matrix(snp_data)

# Display the resulting Spectral kernel matrix
print(spectral_matrix[1:5, 1:5])

##############################
# Composite Kernels:
#
# combining multiple kernels into a composite kernel
# can be effective. For example, a linear kernel combined
# with a Gaussian kernel can capture both linear and
# non-linear relationships.
#
# In the composite kernel, the weights assigned to each
# individual kernel (linear and Gaussian, in this case)
# are controlled by the parameters alpha and 1 - alpha.
# The purpose of these weights is to determine the
# contribution of each individual kernel to the overall
# composite kernel.
#
# If you have three or four kernels to combine, you can
# extend the idea by introducing additional weights
# for each kernel. For example, for three kernels
# (linear, Gaussian, and polynomial), you could define the
# composite kernel as follows:
#
# composite_kernel=α×linear_term+β×gaussian_term+γ×polynomial_termcomposite_kernel=α×linear_term+β×gaussian_term+γ×polynomial_term
#
# Here, αα, ββ, and γγ are weights for the linear,
# Gaussian, and polynomial kernels, respectively. The weights should be non-negative and sum to 1 to ensure that the composite kernel remains a valid kernel.

# Composite Kernel (Linear + Gaussian)
composite_kernel <- function(X, sigma = 1, alpha = 0.5) {
  linear_term <- linear_kernel(X)
  gaussian_term <- gaussian_kernel(X, sigma)
  return(alpha * linear_term + (1 - alpha) * gaussian_term)
}
#########################

# Install and load necessary libraries

rm(list = ls())
# Toy SNP data
set.seed(123)
snp_data <- matrix(rbinom(500, 1, 0.5), ncol = 50, nrow = 10)  # Example SNP data (binary)
######
# Forming a graph matrix from SNP data involves defining
# relationships between individuals based on genetic similarity
# or other relevant criteria.
# Here's a simple example of how you might construct a graph
# matrix based on pairwise genetic similarity using the Jaccard index
# for binary SNP data:
# Function to compute Jaccard similarity between two SNP profiles
# Pros: Simple and interpretable. It captures the proportion of
        #shared alleles between individuals.
# Cons: It treats presence or absence of a SNP as equally important,
         #regardless of allele frequency.
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
}
# Measures the number of positions at which two binary sequences (SNP profiles) differ.
# It's the complement of the Jaccard similarity.
# Pros: Straightforward and easy to interpret. It directly measures
#       the number of differing positions.
# Cons: It may not be sensitive to the biological relevance of
#        the differences.
hamming_distance <- function(x, y) {
  return(sum(x != y))
}

# Euclidean Distance:
#
#   Treats SNP profiles as points in a high-dimensional space and
# computes the Euclidean distance between them.

# Pros: Considers the magnitude of differences between SNP profiles.
#        Suitable for continuous data.
# Cons: May be sensitive to scale and may not be the
#       best choice for binary data.

euclidean_distance <- function(x, y) {
  return(sqrt(sum((x - y)^2)))
}

# Cosine Similarity:
#
#   Measures the cosine of the angle between two SNP profiles,
# treating them as vectors in a binary space.

# Pros: Accounts for the angle between SNP profiles,
#        providing a measure of direction rather than magnitude.
# Cons: Sensitive to the length of SNP profiles.

cosine_similarity <- function(x, y) {
  return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
}

# Sokal-Michener Similarity:
#
#   Takes into account the number of matching and non-matching positions,
# providing a more complex measure than Jaccard.

# Pros: Takes into account the number of matching and
#       non-matching positions.
# Cons: Similar to Jaccard, but may offer different weighting.

sokal_michener_similarity <- function(x, y) {
  a <- sum(x & y)
  b <- sum(x & !y)
  c <- sum(!x & y)
  d <- sum(!x & !y)
  return((a + d) / (a + b + c + d))
}
###
# Dice Similarity:
#
#   Similar to Jaccard, but with a slightly different weighting of
# true positives.

#
#   Pros: Similar to Jaccard but provides a different
#      weighting for true positives.
# Cons: May not be suitable for all types of genomic data.

dice_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x) + sum(y)
  return(2 * intersection / union)
}

# In graph kernels, the choice of similarity
# measure often depends on the properties of the graph you want
# to construct. For instance, cosine similarity and Euclidean distance
# are often used when treating SNP profiles as vectors in a
# high-dimensional space. The most appropriate measure may also depend
# on the nature of relationships you expect in your genomic data.
# It's recommended to experiment with multiple similarity measures
# and assess their performance using validation or downstream analysis.

# Function to construct a graph matrix based on Jaccard similarity
construct_graph_matrix <- function(snp_data) {
  n <- nrow(snp_data)
  graph_matrix <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        similarity <- jaccard_similarity(snp_data[i, ], snp_data[j, ])
        graph_matrix[i, j] <- similarity
      }
    }
  }

  return(graph_matrix)
}

# Graph Kernel (Using an example graph)
# Function to compute normalized graph Laplacian kernel matrix
# Function to compute normalized Laplacian kernel matrix
compute_normalized_laplacian_kernel <- function(graph_matrix) {
  n <- nrow(graph_matrix)

  # Compute the degree matrix
  degree_matrix <- diag(rowSums(graph_matrix))

  # Compute the Laplacian matrix
  laplacian_matrix <- degree_matrix - graph_matrix

  # Compute the normalized Laplacian kernel matrix using the pseudo-inverse
  normalized_laplacian_kernel <- MASS::ginv(sqrt(degree_matrix)) %*% laplacian_matrix %*% ginv(sqrt(degree_matrix))

  return(normalized_laplacian_kernel)
}

###
# Function to compute unnormalized Laplacian kernel matrix
compute_unnormalized_laplacian_kernel <- function(graph_matrix) {
  n <- nrow(graph_matrix)

  # Compute the degree matrix
  degree_matrix <- diag(rowSums(graph_matrix))

  # Compute the Laplacian matrix
  laplacian_matrix <- degree_matrix - graph_matrix

  # Compute the unnormalized Laplacian kernel matrix
  unnormalized_laplacian_kernel <- laplacian_matrix

  return(unnormalized_laplacian_kernel)
}
#####
##############
# Function to compute symmetric normalized Laplacian kernel matrix
# Function to compute symmetric normalized Laplacian kernel matrix
compute_symmetric_normalized_laplacian_kernel <- function(graph_matrix) {
  n <- nrow(graph_matrix)

  # Compute the degree matrix
  degree_matrix <- diag(rowSums(graph_matrix))

  # Compute the Laplacian matrix
  laplacian_matrix <- degree_matrix - graph_matrix

  # Compute the symmetric normalized Laplacian kernel matrix using the pseudo-inverse
  symmetric_normalized_laplacian_kernel <- ginv(sqrt(degree_matrix)) %*% laplacian_matrix %*% ginv(sqrt(degree_matrix))

  return(symmetric_normalized_laplacian_kernel)
}
#######

# In the context of Genomic Best Linear Unbiased Prediction (GBLUP) models,
# the choice of kernel can significantly impact the performance of the model.
# The Anova radial basis kernel is just one option among several,
# and its advantages depend on the characteristics of the data and
# the underlying genetic architecture. Here are some considerations:
#
#   Capturing Non-Linear Relationships:
#   The Anova radial basis kernel, like other non-linear kernels,
#   can capture complex and non-linear relationships in the data.
#  This is advantageous when the relationship between genetic markers and
#  the trait of interest is not well-approximated by a linear model.
#
# Feature Interaction:
#   The Anova radial basis kernel is particularly useful for capturing
#   interactions among features. In genomic data, interactions between
#   genetic markers can play a crucial role in determining phenotypic
#   outcomes. The Anova kernel can model these interactions effectively.
#
# Flexibility:
#   Non-linear kernels, including the Anova radial basis kernel,
#   provide more flexibility in modeling complex genetic architectures.
#  They can capture patterns that linear kernels might miss.
#
# Parameter Tuning:
#   The Anova radial basis kernel has a hyperparameter, gamma,
#   that controls the width of the Gaussian function. Proper tuning of
#  this parameter is crucial for achieving optimal model performance.
#  It might require careful cross-validation to find the best value for
#   your specific dataset.
#
# Interpretability:
#   While non-linear kernels can capture complex relationships,
#  they might be less interpretable than linear kernels.
#  Interpretability is an important consideration, especially in
#  genomic studies where understanding the biological relevance of
#  features is crucial.
#
# It's important to note that the choice of the kernel
# is problem-dependent, and there is no one-size-fits-all solution.
# The performance of the Anova radial basis kernel should be compared
# with other kernels (linear, polynomial, etc.) through proper validation
# procedures. It's advisable to experiment with different kernels and
# evaluate their performance in terms of predictive accuracy,
# computational efficiency, and interpretability.
#
# When applying GBLUP models, it's also common to explore genomic
# relationship matrices based on various kernels and select the one that
# provides the best predictive performance for your specific dataset
# and trait of interest.

# Function to compute Anova radial basis kernel matrix
compute_anova_rbf_kernel_custom <- function(data, gamma = 1) {
  n <- nrow(data)

  # Initialize the kernel matrix
  anova_rbf_kernel_matrix <- matrix(0, n, n)

  # Compute the Anova radial basis kernel matrix
  for (i in 1:n) {
    for (j in 1:n) {
      diff_squared <- sum((data[i, ] - data[j, ])^2)
      anova_rbf_kernel_matrix[i, j] <- exp(-gamma * diff_squared)
    }
  }

  return(anova_rbf_kernel_matrix)
}

# Example Usage
set.seed(123)
snp_data <- matrix(rbinom(500, 1, 0.5), ncol = 50, nrow = 10)  # Example SNP data (binary)

# Compute Anova radial basis kernel matrix
gamma_value <- 0.1
anova_rbf_kernel_matrix_custom <- compute_anova_rbf_kernel_custom(snp_data, gamma = gamma_value)



# Example Usage
set.seed(123)
snp_data <- matrix(rbinom(500, 1, 0.5), ncol = 50, nrow = 10)  # Example SNP data (binary)
graph_matrix <- construct_graph_matrix(snp_data)


graph_matrix <- construct_graph_matrix(snp_data)

graph_kernel_result <- compute_normalized_laplacian_kernel(graph_matrix)

graph_kernel_result2 <-compute_unnormalized_laplacian_kernel(graph_matrix)

graph_kernel_result3 <-compute_symmetric_normalized_laplacian_kernel(graph_matrix)


