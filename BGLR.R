library(BGLR)

rm(list = ls()); gc()
# example dataset
data(wheat)

# pedigree A
A <- wheat.A

# response
y <- wheat.Y[,1]

# genotypes
M <- wheat.X

# G matrix
vars <- apply(M, 2, var)
G <- tcrossprod(scale(M, scale = FALSE)) / sum(vars)

# add something to diagional to stabilize things
diag(G) <- diag(G) + 0.01

# using L as design matrix will lead to an equivalent animal model, which is just
# nicer to handle in BGLR
L <- t(chol(G))

# fixed effects (just intercept)
X <- array(1, dim = c(length(y), 1))

# run animal model
niter = 200
burnin = 50
mod2 <- BGLR(
  y = y,
  response_type = "gaussian",
  ETA = list(
    list(X = X, model = "FIXED"),
    #list(K = G, model = "RKHS", saveEffects = TRUE)
    list(X = L, model = "BRR", saveEffects = TRUE)
  ),
  nIter = niter,
  burnIn = burnin,
  thin = 1,
  saveAt = "myModel",
  verbose = FALSE
)

(mod$SD.yHat)^2
# get the posterior (only stored in file)
uFile <- mod$ETA[[2]]$NamefileOut
uDat <- scan(uFile, what = numeric(), sep = "\n")

eDat <- scan("myModelvarE.dat",  what = numeric(), sep = "\n")

postit <- (burnin + 1) : niter

varU <- uDat[postit]
varE <- eDat[postit]
h2 <- varU / (varU + varE)

# the heritability and its standard deviation (~ standard error)
h2Estimate = mean(h2)
h2SD = sd(h2)

# now reliabilities.
# we specified to store the posterior samples of the random effects with the option "saveEffects = TRUE".
# read those samples in now
u = readBinMat('myModelETA_2_b.bin')

# the breeding values (posterior)
g = L %*% t(u)


# PEV
PEV = apply(g, 1, var)

# reliabilities
vU <- mean(varU)
REL = 1 - (PEV / vU)

AA=  Pred |>
  dplyr::group_by(gen_name) |>
  dplyr::summarise(Pred = mean(pred))

AA= as.data.frame(Pred |>
  dplyr::group_by(across(all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
  dplyr::summarise(Pred = mean(pred))
)
