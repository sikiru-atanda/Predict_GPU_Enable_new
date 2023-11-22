##Computing the point estimates and standard errors with mixed models using matrices

## https://francish.net/post/computing-point-estimates-and-standard-errors-with-mixed-models/

# The standard ordinary least squares model (OLS) can be written as:
#   y = Xβ +  ε
#   where X is a design matrix, β is a column vector of coefficents, plus some random error
#   ε. Of interest is Beta which is obtained through:
#     β=(X′X)^1(X′y)
#   For example (compare results using matrices and just using the coef function):

rm(list = ls()); ls()
gc()
graphics.off()
cat('\014')

ols1 <- lm(mpg ~ wt + am, data = mtcars) #toy example
## X is the design matrix of the predictor including the intercept
X <- model.matrix(ols1)
y <- mtcars$mpg
(betas <- solve(t(X) %*% X) %*% t(X) %*% y)

#mtcars$mpg[5:8] = NA
# The variance/covariance matrix of the betas,
# Var(^β), of which the square root of the diagonal is the standard error, is:
#   σ2(X′X)^-1.σ2
# is based on the residuals and the residuals are merely the observed value less the predicted value:

## Predicted value
yhat = X %*% betas

### Residual
modres <- residuals(ols1)

### manual computation
modres2 <- y - X %*% betas #manual computation
## compare
head(modres)

yhat = X %*% betas

# σ2 is the sum of squared residuals divided by  N−k
# where k is the number of predictors (including the intercept).

sigma(ols1)^2 #just using a function
k = 3 # number of predictors including the intercept
## deviance(model)/residual sum of squares of the model
(sigma2 <- sum((modres2)^2)/(length(y)-3))

# Also, the sum of squared residuals can be computed using
# u′u  where u contains the residuals.
# (we use as.numeric to convert it into a scalar that we can use to multiply later on).

n <- nobs(ols1)
k <- ols1$rank

u <- modres2
(sigma2 <- as.numeric(t(u) %*% u / (n - k)))

## Manual computation of the s.e using matrixes
#Using σ2 , we can obtain the variance/covariance matrix:
(vcm <- sigma2 * solve(t(X) %*% X)) #manual computation

vcov(ols1) #compare when using a function

#The square root of the diagonal is the standard error. Can compare to the results using the

se <- sqrt(diag(vcm))
tstat <- betas / se
pval <-  2 * pt(-abs(tstat), df = n - k)
data.frame(betas, se, tstat, pval)

summary(ols1)$coefficients #compare results

###############################
# Using matrices to reconstruct the mixed model coefficients --------------
##########################################


library(lme4)
data(bdf, package = 'mlmRev')
mlm1 <- lmer(aritPOST ~ aritPRET + sex +
               schoolSES + (sex | schoolNR), data = bdf)

# The data are sorted by cluster already
X <- model.matrix(mlm1) #the design matrix
B <- fixef(mlm1) #to extract the coefficients from lme4-- just as a check, don't use it here
y <- mlm1@resp$y #outcome
Z <- getME(mlm1, 'Z') #sparse Z design matrix
b <- getME(mlm1, 'b') #random effects
Gname <- names(getME(mlm1, 'l_i')) #name of clustering variable
js <- table(mlm1@frame[, Gname]) #how many observation in each cluster

qq <- getME(mlm1, 'q') #columns in Z matrix
NG <- ngrps(mlm1)
nre <- getME(mlm1, 'p_i') #qq/NG --> number of random effects
inde <- cumsum(js) #number per group summed to create an index
cols <- seq(1, (NG * nre), by = nre) #dynamic columns for Z matrix.
#if a random intercept model, cols is only a sequence from 1
#to NG
G <- bdiag(VarCorr(mlm1)) #G matrix

#############

ml <- list() #empty list to store matrices

for (i in 1:NG){

  if (i == 1) { #if first cluster
    st = 1} else { #start at row 1
      st = inde[i - 1] + 1}

  end = st + js[i] - 1

  nc <- cols[i]
  ncend <- cols[i] + (nre - 1)

  Zi <- Z[st:end, nc:ncend] #depends on how many obs in a cluster and
  # how many rand effects
  ml[[i]] <- Zi %*% G %*% t(Zi) + diag(sigma(mlm1)^2, nrow = js[i]) #ZGZ' + r
}
