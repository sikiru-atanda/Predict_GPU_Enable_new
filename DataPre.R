rm(list = ls()); ls()
gc()
library(lme4GS)
data(wheat599)
rm(wheat.Pedigree)
GD = unique(wheat.Pheno$GID)

set.seed(123)
GDD = GD[GD%in%sample(GD, 200, replace = F)]

wheat.Pheno = wheat.Pheno[wheat.Pheno$Rep==1, ]

#wheat.Pheno = wheat.Pheno[wheat.Pheno$Env==1, ]

wheat.Pheno$PlantHeight = rnorm(nrow(wheat.Pheno), mean = 78, sd= 35)

wheat.Pheno$DaysToMat = rnorm(nrow(wheat.Pheno), mean = 56, sd= 10)

wheat.Pheno$Mositure_Content = rnorm(nrow(wheat.Pheno), mean = 21, sd= 6)

wheat.Pheno = wheat.Pheno[wheat.Pheno$GID%in%GDD, ]


#wheat.Pheno = wheat.Pheno[which(wheat.Pheno$GID%in%sample(GD, 200, replace = F)), ]

#wheat.Pheno = wheat.Pheno[!wheat.Pheno$Env%in%c(1, 2), ]
wheat.X = wheat.X[rownames(wheat.X)%in%unique(wheat.Pheno$GID), ]
pheno = wheat.Pheno
Geno_data = wheat.X

genotype = "GID"

pheno$Env = as.factor(pheno$Env)
pheno$Rep = as.factor(pheno$Rep)
pheno$GID = as.factor(pheno$GID)


rm(list = ls()); ls()
gc()
cat('\014')
setwd("D:/PredictProR")
load('WheatPhenoGeno.Rdata')

gen_name = "GID"

pheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])

pheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR"), ])

pheno$PH = pheno$Yield

response = "Yield"

Omic2 = COP

response = c("Yield","PH")

random = ~GID
fixed = NULL
GS_model = "BRR"


#fixed = ~Env

ZE<-model.matrix(~factor(pheno$Env)-1)
dim(ZE)
dim(Geno.data)

#Environment–marker design matrix
XEM = model.matrix(~0+Geno.data:XE)
#pheno = droplevels(pheno[pheno$Env%in%c("B2IR"), ])


pheno$GID = as.character(pheno$GID)
pheno$Env = as.character(pheno$Env)
str(pheno)

response = c("Yield")
#"B2IR","F5I", "EHT", "B5I"
# [1] B5I   BLHT  B2IR  F5I   FDRIP EHT
# Levels: B2IR B5I BLHT EHT F5I FDRIP

# library(PredictiveAnalytic)
#
# GSModels <- PredictiveAnalytic::GS_mod(
#     pheno=pheno,
#     pheno_train = NULL,
#     pheno_test = NULL,
#     Gmatrix= NULL,
#     Omic2_Relationship_Matrix = NULL,
#     Omic3_Relationship_Matrix = NULL,
#     Geno_data = Geno_data,
#     Omic2_data = NULL,
#     Omic3_data = NULL,
#     train_Geno_data = NULL,
#     test_Geno_data = NULL,
#     train_Omic2_data = NULL,
#     test_Omic2_data = NULL,
#     train_Omic3_data = NULL,
#     test_Omic3_data = NULL,
#     train_set = NULL,
#     test_set = NULL,
#     Gmatrix_method = "Yang",
#     Kernel_matrix_method = NULL,
#     engine= "BGLR",
#     response="Yield",
#     genotype="GID",
#     cova=NULL,
#     fixed= NULL,
#     random=~GID+ GID:Env,
#     heter_resid=FALSE,
#     heter_groups="Env",
#     VarCov_str = NULL,
#     weights =NULL,
#     nIter=NULL,
#     burnIn=NULL,
#     thin=NULL,
#     GS_model = "BRR",
#     Fixed_term_model_Bayesian = 'FIXED',
#     Rand_term_model_Bayesian = "BRR",
#     # Cross_validation = c("Hold_Out",
#     #                     "Stratified_Hold_Out",
#     #                     "Repeated_Hold_Out",
#     #                     "Repeated_Stratified_Hold_Out",
#     #                     "K-Folds",
#     #                     "Stratified_K-Folds",
#     #                     "Repeated_K-Folds",
#     #                     "Repeated_Stratified_K-Folds",
#     #                     "Leave_one_Out",
#     #                     "CV1",
#     #                     "Repeated_CV1",
#     #                     "CV2",
#     #                     "Repeated_CV2"),
#
#     Cross_validation = NULL,
#     test_size = 0.5,
#     random_state = 123,
#     Replication = 2,
#     nFolds = 3,
#     CV = 2,
#     core = NULL)
