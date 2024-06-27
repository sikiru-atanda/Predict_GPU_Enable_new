rm(list = ls()); ls()
gc()
graphics.off()

setwd("D:/PredictProR")
load("barley_data.Rdata")

install.packages('mt')

library('mt')

library(caret)
nzv <- nearZeroVar(geno_data, saveMetrics = TRUE)
filtered_data <- geno_data[, !nzv$nzv]

corr_matrix <- cor(filtered_data)
high_corr <- findCorrelation(corr_matrix, cutoff = 0.9)
filtered_data <- filtered_data[, -high_corr]

# Assuming pheno_data1 and filtered_data are already defined
response <- pheno_data1$B_GLUCAN

library(caret)

# Number of replications
Nreps <- 5
reps <- list()

# Run the RFE multiple times and store the results
for (rep in 1:Nreps) {
  control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
  results <- rfe(filtered_data, pheno_data1$B_GLUCAN, sizes = c(1:100), rfeControl = control)
  
  reps[[rep]] <- results$optVariables
}

# Find the common features across all replications
common_features <- Reduce(intersect, reps)

# Print the common features
print(common_features)


control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
results <- rfe(geno_data, pheno_data1$HT, 
               #sizes = c(1:50), 
               rfeControl = control)
print(results)

View(results$optVariables)

length(results$optVariables)

###
### Regularization Methods (Lasso):
library(glmnet)
x <- as.matrix(filtered_data)
y <- pheno_data1$B_GLUCAN
cv_fit <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_fit$lambda.min

res_lasso <- list()
Nreps <- 5

for(rep in  1:Nreps){
lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
selected_features <- which(coef(lasso_model) != 0)

res_lasso[[rep]] <- selected_features

}

# Find the common features across all replications
common_features_lasso <- Reduce(intersect, res_lasso)

# Print the common features
print(common_features_lasso)

library(caret)
library(GA)
# Create control function for genetic algorithm
control <- gafsControl(functions = rfGA, method = "cv", number = 5)
# Perform GAFS using random forest as the model
results <- gafs(x =geno_data, y = pheno_data1$B_GLUCAN, iters = 10, gafsControl = control)
# List the chosen features
length(results$optVariables)
