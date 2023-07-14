library(PredictProR)

#rm(list = ls())
rm(sik)
sik= PredictProR::model_execute(pheno_data = pheno,
                           geno_data = Geno.data,
                           gen_name = "GID",
                           omic1_data = COP,
                           #omic2_data = COP,
                           #omic3_data = COP,
                           response = "Yield",
                           GS_model = "BRR",
                           random = ~GID,
                           message = FALSE
                           )

(sik$model_results$coefficients_1)

omic2_clean <-  omic_to_model(omic_data = COP,
                              )

sik = pheno_geno_match(object_pheno = pheno_clean$pheno_data[1:100,],
                       object_geno = geno_clean,
                       gen_name = "GID")

class(sik$object_geno)

dim(sik$object_geno)
sik$test_set

sik$object_geno[1:5, 1:5]
as.character(pheno_clean$pheno_data[1:100,]$GID[1:5])


omic3_clean <-  omic_to_model(omic_data = COP,
)


omic1_clean= omic_to_model(omic_data = COP,
                           train_omic_data = train_omic1_data,
                           test_omic_data = test_omic1_data,
                           message =TRUE)

#sik = omic_precheck(object =  COP)

pheno_clean = phenotype_to_model(pheno = pheno,
                                 response = "Yield",
                                 gen_name = "GID")

geno_clean = geno_to_model(geno_data = Geno.data)


ETA = ETA_compiler_bayes(random = ~GID,
                         GS_model = "BRR",
                         geno_data = geno_clean,
                         omic2_data = omic2_clean,
                         omic3_data = omic3_clean,
                         pheno_data = pheno_clean[[1]])


bayes_para  = bayes_parameter_check(nIter = NULL,
                                    burnIn = NULL)

mod = M_matrix_bayes_mod_single_loc(object = pheno_clean[[1]],
                                    response = "Yield",
                                    ETA = ETA$ETA,
                                    bayes_para = bayes_para)



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
