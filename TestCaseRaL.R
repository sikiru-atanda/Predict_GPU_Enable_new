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


### Parallelization of the number of replication
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

## Extract the name of each environments
ENVs <- as.character(unique(Pheno_data$Env))

## Declare the response variable
trait <- "Yield"

#gen_name = "GID"
## Extract only from the B2IR from the phenotypic data
pheno <- droplevels(Pheno_data[Pheno_data$Env%in%c("FDRIP"), ])

hist(pheno$Yield)
## Models to compare
Model <- c("BayesA", "BayesB", "BayesC", "BL", "RF",
          "xgboost", "Lasso","Ridge_Reg")

## Number of replication
NRep <- 10

## Create empty list for the results output
ypred_cv_Reps_all <-  c()

### Create for loop to loop through the all the models iteratively
for (i in 1:length(Model)) {

  mod = Model[i]

  ### initialize parallization of the number of replication
ypred_cv_Rep <- foreach::foreach(r = 1:NRep
                               ) %dopar% {


## Declare number of k- folds
nfolds <- 3
## length of the response variable
len_y <- length(pheno$Yield)
## Generate the folds
folds_generated <- rep(1:nfolds, length.out = len_y)
## sample/reshuffle the  generated folds
group <- sample(folds_generated, len_y)
table(group)

## Create empyt matrix to store the predicted value for all the folds
ypred_cv <- matrix(data=NA, nrow=len_y, ncol=1)

## Extract the response variable
y <-pheno$Yield   # Response variable

## Loop through the folds
for (g in 1:nfolds){
  ycv <- y

  ## set the g fold to NA as the testing set
  for (j in 1:len_y) {
    if(group[j] == g) { ycv[j]<-NA }
  }

  ## extract the position NA which is the testing set
  tst = which(is.na(ycv))

  #### For Bayes A, B, C and BL
  if(mod=="BayesA"|mod=="BayesB"|mod=="BayesC"|mod=="BL"){

    ### Create environment for the linear predictors and defined model for prediction
    ETA = list(list(X= as.matrix(SNP_data),
                    model=mod,
                    saveEffects=TRUE))

    ## Number of iteration
    nIter = 20000
    ## Number of burning
    burnIn = 5000
    ## Number of thinning
    thin = 10

    fm <- BGLR::BGLR(
      y=ycv,
      ETA = ETA,
      nIter = nIter,
      burnIn =  burnIn,
      thin =  thin,
      verbose = FALSE
      #saveAt =files_key
    )

    ## Extract the predicted values
    preds <- fm$yHat
  }


  #### Random Forest
  if(mod =="RF"){

  fit = randomForest::randomForest(x = SNP_data[-tst, ],
                                      y = ycv[-tst],
                                      ntree = 500,
                                      importance = TRUE)

  preds <- stats::predict(fit,
                          SNP_data[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)



  preds <- preds[, 1]

  }
  ########################
  ### Extreme Gradient Boosting
  ################################
  if(mod =="xgboost"){
  xgb_params <- list(
    booster = "gblinear",
    eta = 0.001,
    max_depth = 6, #This indicates how deep the built tree can be.
    #The deeper the tree, the more splits it has and it captures more
    #information about how the data. We fit a decision tree with depths
    #ranging from 1 to 32 and plot the training and test errors
    gamma = 4,
    subsample = 0.5,
    colsample_bytree = 1,
    objective = "reg:squarederror",
    eval_metric = c("rmse", "rmsle", "mape")
  )


  SNP_Xgb <- xgboost::xgb.DMatrix(data = SNP_data[-tst, ],
                                           label = y[-tst])

  fit <- xgboost::xgb.train(
    params = xgb_params,
    data =  SNP_Xgb,
    nrounds = 5000,
    verbose = 0
  )


  preds <- stats::predict(fit,
                          SNP_data[tst, ],
                          reshape = TRUE)



}




if(mod =="Ridge_Reg"){

#####################
## Ridge Regression
######################
### Fit Ridge regression using glmnet package. Alpha = 0
fit_CV<-glmnet::cv.glmnet(x=SNP_data[-tst, ],
                             y = ycv[-tst],
                             nfolds = 5,
                             alpha = 0,
                             standardize = FALSE,
                             lambda = seq(0.000001,0.9,length.out=100)^4)

fit=  glmnet::glmnet(x=SNP_data[-tst, ],
                        y = ycv[-tst],
                        alpha = 0,
                        standardize = FALSE,
                        lambda =fit_CV$lambda.min)

preds <- stats::predict(fit,
                           SNP_data[tst, ],
                           reshape = TRUE)


preds <- preds[, 1]

}

#######

#####################
## Lasso Regression
######################
if(mod =="Lasso"){
fit_CV<-glmnet::cv.glmnet(x=SNP_data[-tst, ],
                             y = ycv[-tst],
                             nfolds = 5,
                             alpha = 1,
                             standardize = FALSE,
                             lambda = seq(0.000001,0.9,length.out=100)^4)

fit=  glmnet::glmnet(x=SNP_data[-tst, ],
                        y = ycv[-tst],
                        alpha = 1,
                        standardize = FALSE,
                        lambda =fit_CV$lambda.min)

preds <- stats::predict(fit,
                        SNP_data[tst, ],
                        reshape = TRUE)


preds <- preds[, 1]

}
  ypred_cv[tst] <- preds

}

(PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))

ypred_cv_Rep <- PA

}

ypred_cv_Reps_all[[i]] <- data.frame(PredACC = do.call(rbind, ypred_cv_Rep),
                                     Rep = paste("Rep", 1:NRep, sep = "_"), Model = mod)


}

### rbind the outputs for all the replications
All_Model_ypredCV = do.call(rbind, ypred_cv_Reps_all)

### Set model into factor
All_Model_ypredCV$Model = as.factor(All_Model_ypredCV$Model)

## Find the mean for each model
means <-aggregate(PredACC ~  Model, All_Model_ypredCV, mean)

## Round it to 2 decimal place
means$PredACC <- round(means$PredACC, 2)

#Make the box plot
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


