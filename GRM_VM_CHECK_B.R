#############
# CHECK GRM #

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/GRAL/Bank_LMM/Queries/Sikiru")
list.files()
load("C:/Users/sgeza/OneDrive/Desktop/GRAL/Bank_LMM/Queries/Sikiru/Wdata.Rdata")

###################
# First check quality of GRM

library(ASRgenomics)
dim(GRM)

G_summary <- kinship.diagnostics(K = GRM)

# So, there are 16 records that are potentially duplicates... 
# this is the main issue

# Matrix dimension is: 281x281
# Range diagonal values: 1.213 to 2.27264
# Mean diagonal values: 1.79401
# Range off-diagonal values: -0.4246 to 2.06884
# Mean off-diagonal values: -0.0064
# There are 281 extreme diagonal values, outside < 0.8 and > 1.2
# There are 16 records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  0.95

# Aparently diagonals are fine, but they go as low as 1.2.. is this correct?
# Not sure if that is good if you are dealing with RIL or DH.. .that should be close to 2.

G_summary$list.duplicate
G_summary$plot.diag
G_summary$plot.offdiag

# In the off-diagonal plot you can clearly see a group of individuals that are duplicates
# (anything above 1 really is potentially a problem)

# Lets see if we can get the inverse of that matrix
GinvA <- G.inverse(GRM, sparse=TRUE)
# Any value smaller than 1e-5 is usually a warning of an unstable matrix
# hence, this is possible a bad one.

# Reciprocal conditional number for original matrix is: 1.9100341420263e-06
# Reciprocal conditional number for inverted matrix is: 1.72759778226934e-06
# Inverse of matrix G appears to be ill-conditioned.

head(GinvA$Ginv.sparse) # Numbers not too large but a bit large

#########
# Lets try some blending 5%

GRM.blend <- G.tuneup(G=GRM, blend=TRUE, pblend=0.05)
# Much better conditonal number here.

GinvB <- G.inverse(GRM.blend$Gb, sparse=TRUE)
# This is overal much better
head(GinvB$Ginv.sparse) # Smaller numbers
# Reciprocal conditional number for original matrix is: 0.00028701685318875
# Reciprocal conditional number for inverted matrix is: 0.000297182930509034
# Inverse of matrix G does not appear to be ill-conditioned.

#######################################
# Lets proceed to check

# BOTH TRUE but clearly very different quality
matrixcalc::is.positive.definite(GRM)
matrixcalc::is.positive.definite(GRM.blend$Gb)

######################################
# Lets proceed to fit GBLUP models

aggregate(Yield~Env, FUN=mean, data=pheno)
aggregate(Yield~Env, FUN=var, data=pheno) # Variances are very different

# Hence, heterogeneous erros incorporated

head(pheno)
pheno$GID <- as.factor(pheno$GID)
pheno$Env <- as.factor(pheno$Env)

dim(GRM)
attr(GRM, 'INVERSE') <- FALSE

# Lets try using GRM
sik1 = asreml::asreml(fixed = Yield ~ Env,
                      random = ~  vm(GID, GRM):corgh(Env),
                      residual = ~dsum(~units|Env),
                      data = pheno)
summary(sik1)$varcomp

# This one indicates singulartiy
# but the problem is on the second term, as that is not correct
sik_1B = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, GRM) + vm(GID, GRM):corgh(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)

# Using GinvA$Ginv.sparse
sik2 = asreml::asreml(fixed = Yield ~ Env,
                      random = ~  vm(GID, GinvA$Ginv.sparse):corgh(Env),
                      residual = ~dsum(~units|Env),
                      data = pheno)
summary(sik2)$varcomp # almost identical results

# Using GinvA$Ginv.sparse
sik_2 = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, G_inv) + vm(GID, GinvA$Ginv.sparse):corgh(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)

# Lets test if we can use GinvB instead
sik_2 = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, G_inv) + vm(GID, GinvB$Ginv.sparse):corgh(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)

####################
# THESE TWO MODELS SHOULD WORK!!

# It does not work...
# Now, the above model is strange... you do the interaction as this:
sik_2 = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, G_inv) + vm(GID, GinvA$Ginv.sparse):idv(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)
# So, the above does not like GinvA
# Lets try with GinvB
sik_2 = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, G_inv) + vm(GID, GinvB$Ginv.sparse):idv(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)

######################
# Lets try now providing only the GRM
dim(GRM)
attr(GRM, 'INVERSE') <- FALSE

# THe one below runs fine...
sik_2 = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, GRM) + vm(GID, GRM):idv(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)
summary(sik_2)$varcomp

# But the next one does not. (note corgh)
asreml.options(ai.sing=TRUE)
sik_2C = asreml::asreml(fixed = Yield ~ Env,
                       random = ~  vm(GID, GRM) + vm(GID, GRM):corgh(Env),
                       residual = ~dsum(~units|Env),
                       data = pheno)
summary(sik_2C)$varcomp

# And this runs fine (note diag)
asreml.options(ai.sing=TRUE)
sik_2C = asreml::asreml(fixed = Yield ~ Env,
                        random = ~  vm(GID, GRM) + vm(GID, GRM):diag(Env),
                        residual = ~dsum(~units|Env),
                        data = pheno)
summary(sik_2C)$varcomp

# Still the best model is:
sik1 = asreml::asreml(fixed = Yield ~ Env,
                      random = ~  vm(GID, GRM):corgh(Env),
                      residual = ~dsum(~units|Env),
                      data = pheno)
summary(sik1)$varcomp
# And it does not matter if you use GRM, GinvA or GinvB