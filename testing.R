
rm(list = ls())

library(PredictProR)

#rm(list = ls())

#1, 2,
rm(sik)
sik= PredictProR::model_execute(pheno_data = pheno,
                           geno_data = Geno.data,
                           gen_name = "GID",
                           #gmatrix = COP,
                           gkernel = COP,
                           gmatrix_method = "Yang",
                           #omic1_kernel = COP,
                           #omic2_kernel = COP,
                           #omic3_kernel = COP,
                           #omic1_data = COP,
                           #omic2_data = COP,
                           #omic3_data = COP,
                           response = "Yield",
                           GS_model = "BRR",
                           #rand_term_model_bayesian = "BayesB",
                           random = ~GID,
                           message = FALSE
                           )



sik$summary_statistic$model_type

sik$model_results$genomic_variance_1



sik2 <- predict_proR(object_coeff1 = sik$model_results$coefficients_1,
                  object_coeff2 = sik$model_results$coefficients_2,
                  geno)

sik$model_results$intercept


GEBV = Geno.data%*%AA

head(GEBV)

head(GEBV2$GEBV)

GEBV2 = sik$model_results$EBV_1

AA2 = sik$model_results$coefficients_1

AA = as.matrix(AA)

class(AA)

omic1_clean <-  omic_to_model(omic_data = COP,
                              )

sik = pheno_geno_match(object_pheno = pheno_clean$pheno_data[1:100,],
                       object_geno = geno_clean,
                       gen_name = "GID")

class(sik$object_geno)

dim(sik$object_geno)
sik$test_set

sik$object_geno[1:5, 1:5]
as.character(pheno_clean$pheno_data[1:100,]$GID[1:5])

response = "Yield"

weights = NULL

random = ~GID
fixed = NULL

omic1_clean <-  omic_to_model(omic_data = COP,
)


omic1_clean= omic_to_model(omic_data = COP,
                           train_omic_data = train_omic1_data,
                           test_omic_data = test_omic1_data,
                           message =TRUE)

#sik = omic_precheck(object =  COP)

pheno_clean = phenotype_to_model(pheno = pheno,
                                 response = response,
                                 gen_name = "GID")

geno_clean = geno_to_model(geno_data = Geno.data)


ETA = ETA_compiler_bayes(random = ~GID,
                         GS_model = "BRR",
                         geno_data = geno_clean,
                         #omic1_data = omic1_clean,
                         #omic3_data = omic3_clean,
                         pheno_data = pheno_clean[[1]])

ETA = ETA_compiler_bayes_GBLUP(random = ~GID,
                               GS_model = "RKHS",
                               #geno_data = geno_clean,
                               omic1_kernel = omic1_clean,
                               #omic3_data = omic3_clean,
                               pheno_data = pheno_clean[[1]])


bayes_para  = bayes_parameter_check(nIter = NULL,
                                    burnIn = NULL,
                                    thin = NULL)

mod2 = M_matrix_bayes_mod_single_loc(object = pheno_clean[[1]],
                                    response = "Yield",
                                    ETA = ETA$ETA,
                                    bayes_para = bayes_para)

mod = bayes_mod_execute(object = pheno_clean[[1]],
                                    response = "Yield",
                                    ETA = ETA$ETA,
                                    bayes_para = bayes_para)

mod2$model$yHat


res_model_output <- mod_output_bayes(mod = mod,
                                     ETA = ETA,
                                     geno_data = geno_model_ready,
                                     gen_name = gen_name,
                                     omic1_data = omic1_clean,
                                     omic2_data = NULL,
                                     omic3_data = NULL)

res_model_output$predict_value


res_model_output$coefficients

res_summary_stat <- summary_statistics_bayes(mod = mod)


mod$model$mu
geno_data = geno_clean
omic2_data = omic2_clean
omic3_data = omic3_clean

fixed_term_model_bayesian = NULL;
rand_term_model_bayesian = NULL;
GS_model = "BRR";
nIter = NULL;
burnIn = NULL;
thin = NULL;
response = "Yield";
weights= NULL;
gen_name = "GID";
fixed = NULL;
pheno_data=pheno;
pheno_train = NULL;
pheno_test = NULL;
#geno_data = Geno.data;
omic1_data = COP;
omic2_data = COP;
omic3_data = COP;
gmatrix= NULL;
omic1_matrix = NULL;
omic2_matrix = NULL;
omic3_matrix = NULL;
train_geno_data = NULL;
train_omic1_data = NULL;
train_omic2_data = NULL;
train_omic3_data = NULL;
test_geno_data = NULL;
test_omic1_data = NULL;
test_omic2_data = NULL;
test_omic3_data = NULL;
train_coefficient = NULL;
train_set = NULL;
test_set = NULL;
