
# save(list = c("Pheno_data",
#               "SNP_data"),
#      file = "Pheno_SNP_data.Rdata")


## Clean the environment and free memories
rm(list = ls()); ls()
gc()
cat('\014')
### set working directories
setwd("D:/PredictProR")
## Load the R data
load('Pheno_SNP_data.Rdata')

## load required library
library(BGLR)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(randomForest)
library(xgboost)

## Detect the number of core in your machine/laptop
detectCores()

cl <- parallel::makeCluster(4)

doParallel::registerDoParallel(cl)

## Extract the name of each environments

ENVs = as.character(unique(Pheno_data$Env))

trait = "Yield"

gen_name = "GID"

pheno = droplevels(Pheno_data[Pheno_data$Env%in%c("F5I"), ])

hist(pheno$Yield)

### Create environment for the linear predictors and defined model for prediction

Model = c("BRR", "BayesB")
# ETA = list(list(X= as.matrix(SNP_data),
#                 model="BayesC",
#                 saveEffects=TRUE))
NRep <- 3

ypred_cv_Reps_sik= c()

#ypred_cv_Rep <- matrix(data=NA, nrow=1, ncol=NRep)

for (i in 1:length(Model)) {


ETA = list(list(X= as.matrix(SNP_data),
                model=Model[i],
                saveEffects=TRUE))


ypred_cv_Rep <- foreach::foreach(r = 1:NRep,
                               .errorhandling='pass') %dopar% {


nfolds <- 5
len_y <- length(pheno$Yield)
folds_generated <- rep(1:nfolds, length.out = len_y)
group <- sample(folds_generated, len_y)
table(group)


ypred_cv <- matrix(data=NA, nrow=len_y, ncol=1)



nIter = 200

burnIn = 20

thin = 1



y <-pheno$Yield   # Response variable

for (g in 1:nfolds){
  ycv <- y

  for (j in 1:len_y) {
    if(group[j] == g) { ycv[j]<-NA }
  }

  fm <- BGLR::BGLR(
    y=ycv,
    ETA = ETA,
    nIter = nIter,
    burnIn =  burnIn,
    thin =  thin,
    verbose = FALSE
    #saveAt =files_key
    )




  for (j in 1:len_y) {
    #if(group[j] == g) { ypred_cv[j] <- BLUPcv[j] }
    if(group[j] == g) { ypred_cv[j] <- fm$yHat[j] }
  }
}

(PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))

ypred_cv_Rep <- PA

}






ypred_cv_Reps_sik[[i]] <- data.frame(PredACC = do.call(rbind, ypred_cv_Rep),
                           Rep = paste("Rep", 1:NRep, sep = "_"), Model = Model[i])


}

sik = do.call(rbind, ypred_cv_Reps_sik)

ypred_cv_Reps_BayesB <- data.frame(PredACC = do.call(rbind, ypred_cv_Rep),
                                Rep = paste("Rep", 1:NRep, sep = "_"), Model = "BayesB")

ypred_cv_Reps_BayesC <- data.frame(PredACC = do.call(rbind, ypred_cv_Rep),
                                   Rep = paste("Rep", 1:NRep, sep = "_"), Model = "BayesC")

All_Model_ypredCV = rbind(ypred_cv_Reps_BRR,
                          ypred_cv_Reps_BayesB,
                          ypred_cv_Reps_BayesC)


All_Model_ypredCV$Model = as.factor(All_Model_ypredCV$Model)

means <-aggregate(PredACC ~  Model, All_Model_ypredCV, mean)
means$PredACC <- round(means$PredACC, 2)

ggplot(All_Model_ypredCV, aes(x=Model, y=PredACC)) +
  geom_boxplot(aes(fill=Model))+
  #geom_violin() +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5)
  geom_jitter(shape=16, position=position_jitter(width =0.2))+
  theme_classic()+
  stat_summary(fun=mean, colour="darkred", geom="point",
               shape=18, size=3, show.legend=FALSE) +
  #geom_text(data = means, aes(label = PredACC, y = PredACC))
stat_summary(fun.y=mean, colour="red", geom="text", show_guide = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=2)))

