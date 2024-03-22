

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
sik <- model_execute(pheno_data = pheno,
                     #geno_data = Geno.data,
                     omic1_data = COP,
                     #omic2_data = COP,
                     omic3_data = COP,
                     random = ~GID,
                     GS_model_cv = c("BRR","RandomForest", "BayesB", "Lasso"),
                     gen_name = "GID",
                     response =c("Yield","deBLUP","BLUE"),
                     cross_validation = TRUE,
                     metric_for_ranking = "accuracy",
                     eval_metrics = c("accuracy", "mean_squared_error", "bias",
                                      "root_mean_squared_error",
                                      "relative_squared_error",
                                      "mean_absolute_error",
                                      "mean_absolute_percent_error"),
                     replication = 1,
                     cross_validation_meth = "Stratified_Hold_Out",

                     #iteration = 50,
                     ntree=50,
                     system_database = FALSE
)


#sik
###
bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
sik <- model_execute(pheno_data = pheno,
                     gmatrix = GRM,
                     omic1_kernel = COP,
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
sik <- model_execute(pheno_data = pheno,
                     #gmatrix = GRM,
                     #omic1_kernel = GRM,
                     omic2_kernel = GRM,
                     #omic3_kernel = COP,
                     random = ~GID,
                     gen_name = "GID",
                     response ="Yield",
                     GS_model = "GBLUP",
                     engine = "asreml",
                     system_database = TRUE
)

plot(sik$model_results$Predicted_value$Predicted_value,
     sik$model_results$Residual_value$Residual_value)

# multi-environment, single and multi-omics -------------------------------

sik_g = geno_precheck(object_geno = Geno.data)

GRM = grm_calculation(geno_clean = sik_g$snps_matrix, method = "Yang",
)

pheno2 = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])

library(PredictProR)
## asreml model
AA <- model_execute(random = ~ GID + GID:Env ,
                    GS_model = "GBLUP",
                    pheno_data = pheno2,
                    gmatrix = GRM,
                    omic1_kernel = GRM2,
                    gen_name = "GID",
                    engine = "asreml",
                    heter_resid = TRUE,
                    heter_groups = "Env",
                    var_cov_str = "fa2",
                    response = "Yield",
                    system_database = FALSE,
                    inverse = TRUE
                    )

res = asreml::asreml(fixed = Yield~1,
                     random = ~ Env + vm(GID, GRM):rr(Env, 2),
                     residual = ~ dsum( ~ units|Env),
                     na.action = na.method(y = "include", x = "include"),
                     data = pheno
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
