
### Toy example for calculating the standard error from the matrix

mod <- lm(mpg ~ wt + am, data = mtcars) #toy example
## X is the design matrix of the predictor including the intercept
X <- model.matrix(mod)
y <- mtcars$mpg
(betas <- solve(t(X) %*% X) %*% t(X) %*% y)

## Predicted value
yhat = X %*% betas
## Residual
modres <- fm$y[-tst] - fm$yHat[-tst]

# sigma2 is the sum of squared residuals divided by  N−k

k = ncol( Geno.data)# number of predictors including the intercept
N = length(y) ## Number of response variable
## deviance(model)/residual sum of squares of the model
(sigma2 <- sum((modres)^2)/(k-1))

gpuA = gpuMatrix(SNP_data, type="double")
#Using sigma2 , I obtained the variance/covariance matrix:
(vcm <- sigma2 * solve(t(X) %*% X))

#The square root of the diagonal is the standard error of the predictor and the intercept.

se <- sqrt(diag(vcm))


### My question is the estimation of standard error for the predicted value.
?yhat


rm(list = ls()); ls()
gc()
cat('\014')
graphics.off()

x <- matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=T)
xt <- t(x)
#w <- as.vector(c(7,8,9))
# #w needs to be 3x3 so make use diag to construct w as a matrix with those values on the

w <- diag(c(7,8,9))


H = X%*% solve(t(X) %*% X) %*% t(X)

H <- X %*% solve(t(X) %*% X) %*% t(X)

xt %*% w %*% x

A <- c(1, 2, 3)
B <- c(4, 5, 6)

A%*%B

crossprod(A, B)
tcr


g_omic_ebv = g_ebv + omic1_ebv
PEV <- apply(g_omic_ebv, 1, var)

sum_EBV <- data.frame(name = omic1_EBV[, 1],
                      EBV = GEBV$GEBV+omic1_EBV$omic1_EBV,
                      stringsAsFactors = FALSE)
colnames(sum_EBV)[1] <- gen_name


sum_EBV$Std_error <- sqrt(PEV)
sum_EBV$PEV <-   PEV

##Estimate of reliabilities
sum_EBV$Reliability <-  1 - (PEV / Var_U)

sik = list(w = c(1, 2, 4),
           ade = list(my1=A, my2=B))

names(sik) = c("my", 'mys')

names(sik$ade) = c('d', 'd1')



