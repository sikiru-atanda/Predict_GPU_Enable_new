yh#
# GxE using marker-by-environment interactions
# The following examples illustrate how to implement marker-by-environments interaction models using BGLR, for further details about these models see Lopez-Cruz et al. (2015).
#
# (1) As a random regression on markers

rm(list = ls()); gc()

library(BGLR)
data(wheat)
Y=wheat.Y # grain yield evaluated in 4 different environments
round(cor(Y),2)

X=scale(wheat.X)/sqrt(ncol(wheat.X))

y2=Y[,2]
y3=Y[,3]
y4=Y[,4]
y=c(y2,y3, y4)

X0=matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros

X_main=rbind(X,X, X)
X_1=rbind(X,X0, X0)
X_2=rbind(X0,X, X0)
X_3=rbind(X0,X0, X)

fm=BGLR( y=y,ETA=list(
  main=list(X=X_main,model='BRR'),
  int1=list(X=X_1,model='BRR'),
  int2=list(X=X_2,model='BRR'),
  int3=list(X=X_3,model='BRR')
),
nIter=6000,burnIn=1000,saveAt='GxE_',groups=rep(1:3,each=nrow(X))
)
varU_main=scan('GxE_ETA_main_varB.dat')[-c(1:200)]
varU_int1=scan('GxE_ETA_int1_varB.dat')[-c(1:200)]
varU_int2=scan('GxE_ETA_int2_varB.dat')[-c(1:200)]
varU_int3=scan('GxE_ETA_int3_varB.dat')[-c(1:200)]
varE=read.table('GxE_varE.dat',header=FALSE)[-c(1:200),]
varU1=varU_main+varU_int1
varU2=varU_main+varU_int2
varU3=varU_main+varU_int3
h2_1=varU1/(varU1+varE[,1])
h2_2=varU2/(varU2+varE[,2])
h2_3=varU2/(varU2+varE[,3])
COR=varU_main/sqrt(varU1*varU2*varU3)
mean(h2_1)
mean(h2_2)
mean(h2_3)
mean(COR)
###################
# (2) Using genomic relationships
# A model equivalent to the one presented above can be
# implemented using G-matrices (or factorizations of it)
# with off-diagnoal blocks zeroed out for interactions,
# for further detials see Lopez-Cruz et al., 2015.

# NOTE: Heterogenous model not implemented for RKHS in BGLR
library(BGLR)
data(wheat)
Y=wheat.Y # grain yield evaluated in 4 different environments
round(cor(Y),2)
y2=Y[,2]
y3=Y[,3]
y=c(y2,y3)

library(BGData)

G=getG(wheat.X,center=TRUE,scaleG=TRUE,scale=TRUE)

G= grm_calculation(geno_clean = wheat.X, method = "Yang",
)
EVD=eigen(G)
PC=EVD$vectors[,EVD$values>1e-5]
for(i in 1:ncol(PC)){ PC[,i]=EVD$vectors[,i]*sqrt(EVD$values[i]) }

XMain=rbind(PC,PC, PC)
X0=matrix(nrow=nrow(X),ncol=ncol(PC),0) # a matrix full of zeros
X1=rbind(PC,X0, X0)
X2=rbind(X0,PC, X0)
X3=rbind(X0, X0, PC)

LP=list(main=list(X=XMain,model='BRR'),
        int1=list(X=X1,model='BRR'),
        int2=list(X=X2,model='BRR'),
        int3=list(X=X3,model='BRR'))
fmGRM=BGLR(y=y,ETA=LP,nIter=12000,burnIn=2000,saveAt='GRM_',groups=rep(1:3,each=nrow(X)))

plot(fm$yHat,fmGRM$yHat)

cor(y, fm$yHat) # 0.8329287
cor(y, fmGRM$yHat) # 0.8003565

rbind(fm$varE,fmGRM$varE)


