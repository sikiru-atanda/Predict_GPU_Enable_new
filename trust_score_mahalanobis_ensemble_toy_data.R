set.seed(123)  # For reproducibility

###################################
# 1. Create Toy Feature-Matrix Data
###################################

# Parameters
n_train <- 30  # number of training samples
n_test  <- 10  # number of test samples
p       <- 5   # number of features

# Generate training features (X_train) from a multivariate normal
mu_vec <- rep(0, p)
Sigma  <- diag(p) * 1.0  # identity covariance
X_train <- MASS::mvrnorm(n = n_train, mu = mu_vec, Sigma = Sigma)

# Generate training labels (y_train) from a linear model + noise
beta_true <- runif(p, min = -2, max = 2)  # random true coefficients
y_train <- X_train %*% beta_true + rnorm(n_train, mean = 0, sd = 1)

# Generate test features (X_test) from a similar distribution,
# but we can shift them a bit to simulate mild distribution shift
mu_vec_test <- rep(0.5, p)
X_test <- MASS::mvrnorm(n = n_test, mu = mu_vec_test, Sigma = Sigma)

# (Optional) Test labels, if you want ground truth for evaluation
y_test <- X_test %*% beta_true + rnorm(n_test, mean = 0, sd = 1)

# Create a 2D subset for demonstration of density-based OOD
# We'll just take the first 2 columns for the 2D subset
X_train_2D <- X_train[, 1:2]
X_test_2D  <- X_test[, 1:2]

# Print summary
cat("Feature-Matrix Toy Data:\n")
cat("X_train dimensions:", dim(X_train), "\n")
cat("X_test dimensions: ", dim(X_test), "\n")
cat("y_train length:    ", length(y_train), "\n")
cat("y_test length:     ", length(y_test), "\n")
###
###################################
# 2. Create Toy Genomic (SNP) Data
###################################

n_train_geno <- 8   # number of training individuals
n_test_geno  <- 4   # number of test individuals
m_markers    <- 6   # number of SNP markers

# Generate a random SNP matrix in {0,1,2}
# Let's assume minor allele freq ~ 0.3
p_minor <- 0.3

# Training SNP matrix: rows = individuals, cols = markers
SNP_train <- matrix(
  rbinom(n_train_geno * m_markers, size = 2, prob = p_minor),
  nrow = n_train_geno,
  ncol = m_markers
)

# Test SNP matrix
SNP_test <- matrix(
  rbinom(n_test_geno * m_markers, size = 2, prob = p_minor),
  nrow = n_test_geno,
  ncol = m_markers
)

# Compute the Genomic Relationship Matrix (GRM) for training data
# A simple approach: center columns and use (W W^T) / M
compute_GRM <- function(geno_matrix) {
  geno_centered <- scale(geno_matrix, center = TRUE, scale = FALSE)
  M <- ncol(geno_matrix)  # number of markers
  G <- (geno_centered %*% t(geno_centered)) / M
  G
}

G_train_train <- compute_GRM(SNP_train)

# Compute cross-GRM for train vs test.
# We'll center each separately, or use the training centering for both:
geno_train_centered <- scale(SNP_train, center = TRUE, scale = FALSE)
geno_test_centered  <- scale(SNP_test,  center = attr(geno_train_centered, "scaled:center"), 
                             scale  = FALSE)

G_train_test <- (geno_train_centered %*% t(geno_test_centered)) / m_markers

# For demonstration: we won't compute G_test_test, but you could do similarly.

# Generate phenotypes (y_train_geno) from a linear function of SNPs
beta_snp <- runif(m_markers, min = -1, max = 1)
y_train_geno <- geno_train_centered %*% beta_snp + rnorm(n_train_geno, mean = 0, sd = 1)

cat("\nGenomic Data:\n")
cat("SNP_train dimensions:", dim(SNP_train), "\n")
cat("SNP_test dimensions: ", dim(SNP_test), "\n")
cat("G_train_train dimensions:", dim(G_train_train), "\n")
cat("G_train_test dimensions: ", dim(G_train_test), "\n")
cat("y_train_geno length:    ", length(y_train_geno), "\n")

# Fit a basic linear model with SNPs as features
lm_geno <- lm(y_train_geno ~ ., data = data.frame(SNP_train, y_train_geno))
geno_preds_train <- predict(lm_geno)  # training predictions
geno_preds_test <- predict(lm_geno, data.frame(SNP_test))  # training predictions

cat("\nExample linear model on SNP data. \n")
