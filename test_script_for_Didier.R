

# single omics, multi-omics single environment ----------------------------------------

rm(list = ls()); ls()
gc()
cat('\014')
library(PredictProR)
setwd("D:/PredictProR")
load('WheatPhenoGeno.Rdata')

COP = as.matrix(COP)
gen_name = "GID"

pheno = droplevels(pheno.data[pheno.data$Env%in%c("F5I"), ])

library(PredictProR)
bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors", "Lasso", "Ridge_Regression", "deep_learning_model")
sik <- model_execute(pheno_data = pheno2,
                     geno_data = GRM,
                     #omic1_data = COP,
                     #omic2_data = COP,
                     omic3_data = COP,
                     random = ~GID,
                     gen_name = "GID",
                     response ="Yield",
                     GS_model = "BRR",
                     #iteration = 50,
                     ntree=50,
                     system_database = FALSE
)
###
bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
sik <- model_execute(pheno_data = pheno2,
                     gmatrix = GRM,
                     #omic1_kernel = COP,
                     #omic2_kernel = COP,
                     #omic3_kernel = COP,
                     random = ~GID,
                     gen_name = "GID",
                     response ="Yield",
                     GS_model = "GBLUP_BRR",
                     heter_groups = "Env",
                     system_database = FALSE
)

###
library(PredictProR)
asreml_model <- c("GBLUP")
sik <- model_execute(pheno_data = pheno2,
                     #gmatrix = GRM,
                     omic1_kernel = GRM,
                     omic2_kernel = GRM2,
                     #omic3_kernel = COP,
                     random = ~GID,
                     gen_name = "GID",
                     response ="Yield",
                     GS_model = "GBLUP",
                     engine = "asreml",
                     system_database = FALSE
)


# multi-environment, single and multi-omics -------------------------------

sik_g = geno_precheck(object_geno = Geno.data)

GRM = grm_calculation(geno_clean = sik_g$snps_matrix, method = "Yang",
)

pheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])

library(PredictProR)
## asreml model
AA <- model_execute(random = ~ GID + GID:Env ,
                    GS_model = "GBLUP",
                    pheno_data = pheno,
                    gmatrix = GRM,
                    gen_name = "GID",
                    engine = "asreml",
                    heter_resid = TRUE,
                    heter_groups = "Env",
                    var_cov_str = "corgh",
                    response = "Yield",
                    system_database = FALSE
                    )

#### Bayesian model
bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")

AA <- model_execute(random = ~ GID + GID:Env ,
                    GS_model = "GBLUP_BRR",
                    pheno_data = pheno,
                    gmatrix = GRM,
                    gen_name = "GID",
                    heter_resid = TRUE,
                    heter_groups = "Env",
                    response = "Yield",
                    system_database = FALSE
)
