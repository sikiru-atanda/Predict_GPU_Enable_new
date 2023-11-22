
ols1 <- lm(mpg ~ wt + am, data = mtcars) #toy example
## X is the design matrix of the predictor including the intercept
X <- model.matrix(ols1)
y <- mtcars$mpg

beta = solve(t(X)%*% X) %*% t(X) %*% y #for simplicity


yhat <-  X %*% beta


H = X%*% solve(t(X) %*% X) %*% t(X)
#H is the hat matrix and it depends only on Xs

yhat2 = H%*%y

#you need V(yhat), that is:

vc_mY = H*var(y, na.rm = TRUE)*t(H)

se <- sqrt(diag(vc_mY))

# but V(y) = diag(s2y)
# where s2y is just the sum of squares of y divided by n-1... or just the variance.
#
# sum((y)^2)/(length(y)-1)
#
# Also, extending to residuals you have
# but V(y) = diag(s2y)
# where s2y is just the sum of squares of y divided by n-1... or just the variance.
#
# sum((y)^2)/(length(y)-1)
#
# res = y - yhat = y - H*y = [1-H]*y
# so


res = y - (H%*%y)

vcm_res = H*var(y)*t(H)

vcm_res = 1-H*var(y)*t(1-H)

#once more that H matrix is critical
