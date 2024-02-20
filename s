[1mdiff --git a/.RData b/.RData[m
[1mnew file mode 100644[m
[1mindex 0000000..9b1bae3[m
Binary files /dev/null and b/.RData differ
[1mdiff --git a/.Rbuildignore b/.Rbuildignore[m
[1mnew file mode 100644[m
[1mindex 0000000..91114bf[m
[1m--- /dev/null[m
[1m+++ b/.Rbuildignore[m
[36m@@ -0,0 +1,2 @@[m
[32m+[m[32m^.*\.Rproj$[m
[32m+[m[32m^\.Rproj\.user$[m
[1mdiff --git a/.Rhistory b/.Rhistory[m
[1mnew file mode 100644[m
[1mindex 0000000..f304ca0[m
[1m--- /dev/null[m
[1m+++ b/.Rhistory[m
[36m@@ -0,0 +1,512 @@[m
[32m+[m[32mtrain_set = NULL,[m
[32m+[m[32mtest_set = NULL,[m
[32m+[m[32mresponse='Yield',[m
[32m+[m[32mgen_name="GID")[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR"), ])[m
[32m+[m[32msik = phenotype_to_model(pheno = DT,[m
[32m+[m[32mpheno_train = NULL,[m
[32m+[m[32mpheno_test = NULL,[m
[32m+[m[32mtrain_set = NULL,[m
[32m+[m[32mtest_set = NULL,[m
[32m+[m[32mresponse='Yield',[m
[32m+[m[32mgen_name="GID")[m
[32m+[m[32msik = PredictProR::phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno.data,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = phenotype_to_model(pheno = DT,[m
[32m+[m[32mpheno_train = NULL,[m
[32m+[m[32mpheno_test = NULL,[m
[32m+[m[32mtrain_set = NULL,[m
[32m+[m[32mtest_set = NULL,[m
[32m+[m[32mresponse='Yield',[m
[32m+[m[32mgen_name="GID")[m
[32m+[m[32msik =phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = PredictProR::phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = phenotype_to_model(pheno = DT,[m
[32m+[m[32mresponse='Yield',[m
[32m+[m[32mgen_name="GID")[m
[32m+[m[32msik = phenotype_to_model(pheno = pheno.data,[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mgen_name = "GID")[m
[32m+[m[32msik = phenotype_to_model(pheno = DT,[m
[32m+[m[32mresponse='Yield',[m
[32m+[m[32mgen_name="GID")[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = PredictProR::phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mrm(phenotype_to_model, phenotype_precheck)[m
[32m+[m[32msik = phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mgc()[m
[32m+[m[32msik = phenotype_precheck(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = phenotype_precheck(pheno = DT,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32m#' Title[m
[32m+[m[32m#'[m
[32m+[m[32m#' @param pheno[m
[32m+[m[32m#' @param gen_name[m
[32m+[m[32m#' @param response[m
[32m+[m[32m#' @param ...[m
[32m+[m[32m#'[m
[32m+[m[32m#' @return[m
[32m+[m[32m#' @export[m
[32m+[m[32m#'[m
[32m+[m[32m#' @examples[m
[32m+[m[32mphenotype_precheck <- function(pheno,[m
[32m+[m[32mgen_name,[m
[32m+[m[32mresponse,[m
[32m+[m[32m...)[m
[32m+[m[32m{[m
[32m+[m[32mmsg <- sprintf("==================================================\n")[m
[32m+[m[32mif(nrow(pheno)==0) { stop(print(paste(msg, 'No phenotypic records provided.')), call. = FALSE)[m
[32m+[m[32m}[m
[32m+[m[32mif (!inherits(pheno, what = 'data.frame')) {[m
[32m+[m[32mstop(print(paste(msg,"'phenotype' must be of class 'data.frame'")), call. = FALSE)[m
[32m+[m[32m}[m
[32m+[m[32mif (!data.table::is.data.table(pheno)){[m
[32m+[m[32mpheno <- data.table::as.data.table(pheno)[m
[32m+[m[32m}[m
[32m+[m[32m### Check the provided response name correspond to the name in the data file[m
[32m+[m[32mif(!response%in%colnames(pheno)){[m
[32m+[m[32mstop(print(paste(msg,paste(paste("The specified ",  response),[m
[32m+[m[32m" did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)[m
[32m+[m[32m}[m
[32m+[m[32mif(!gen_name%in%colnames(pheno)){[m
[32m+[m[32mstop(print(paste(msg,paste(paste("The specified column",  gen_name),[m
[32m+[m[32m"in the phenotypic data did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)[m
[32m+[m[32m}[m
[32m+[m[32m### Check to ensure no NA in the column GID/name[m
[32m+[m[32mif (anyNA(pheno[, ..gen_name]) || any(pheno[, ..gen_name]==-999)){[m
[32m+[m[32mstop(print(paste(msg,paste(paste('column',  gen_name),[m
[32m+[m[32m'should not have NA/missing'))), call. = FALSE)[m
[32m+[m[32m}[m
[32m+[m[32m# Assign appropriate class.[m
[32m+[m[32mclass(pheno) <- c("data.table", "data.frame", "phenotype")[m
[32m+[m[32mattr(pheno, "cleared") <- "pass"[m
[32m+[m[32mreturn(pheno)[m
[32m+[m[32m}[m
[32m+[m[32msik = phenotype_precheck(pheno = DT,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mrm(phenotype_to_model)[m
[32m+[m[32mrm(phenotype_precheck)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = phenotype_precheck(pheno = DT,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = phenotype_precheck(pheno = DT,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = PredictProR::phenotype_precheck(pheno = DT,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom ~GID)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32msik = PredictProR::phenotype_precheck(pheno = pheno,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = phenotype_to_model(pheno = pheno.data,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32msik = geno_to_model(geno_data = Geno.data)[m
[32m+[m[32mpheno_clean = phenotype_precheck(pheno = pheno,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mclass(pheno_clean)[m
[32m+[m[32mattr(pheno_clean)[m
[32m+[m[32mattr(pheno_clean, "pass")[m
[32m+[m[32mgeno_clean = geno_to_model(geno_data = Geno.data)[m
[32m+[m[32mclass(geno_clean)[m
[32m+[m[32mattributes(geno_clean)[m
[32m+[m[32mattributes(pheno_clean)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32mGS_model[m
[32m+[m[32mGS_model = "BRR"[m
[32m+[m[32mpheno_clean[[1]][m
[32m+[m[32mpheno_clean = phenotype_to_model(pheno = pheno,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32m### Create empty list for ETA[m
[32m+[m[32mETA = list()[m
[32m+[m[32mmsg <- sprintf("==================================================\n")[m
[32m+[m[32m### Get the random terms. Here no interaction terms in the random effect[m
+rand_term_no_inter
[32m+[m[32mobject = pheno_data)[m
[32m+[m[32mpheno_data = pheno_clean[[1]][m
[32m+[m[32m### Get the random terms. Here no interaction terms in the random effect[m
+rand_term_no_inter
[32m+[m[32mobject = pheno_data)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mrm(list = ls()); ls()[m
[32m+[m[32mgc()[m
[32m+[m[32mcat('\014')[m
[32m+[m[32msetwd("D:/PredictProR/PredictProR")[m
[32m+[m[32mload('WheatPhenoGeno.Rdata')[m
[32m+[m[32mgenotype = "GID"[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR"), ])[m
[32m+[m[32mOmic2 = as.matrix(COP)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32mgetwd()[m
[32m+[m[32mpheno_clean <- phenotype_to_model([m
[32m+[m[32mpheno = pheno_data,[m
[32m+[m[32mpheno_train = pheno_train,[m
[32m+[m[32mpheno_test = pheno_test,[m
[32m+[m[32mtrain_set = train_set,[m
[32m+[m[32mtest_set = test_set,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mpheno_clean = phenotype_to_model(pheno = pheno,[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield")[m
[32m+[m[32mgeno_clean = geno_to_model(geno_data = Geno.data)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mrandom = ~GID[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mGS_model = "BRR"[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mbayes_para <-  bayes_parameter_check(nIter = nIter,[m
[32m+[m[32mburnIn = burnIn,[m
[32m+[m[32mthin = thin)[m
[32m+[m[32mnIter = NULL[m
[32m+[m[32mburnIn = NULL[m
[32m+[m[32mthin = NULL[m
[32m+[m[32mbayes_para <-  bayes_parameter_check(nIter = nIter,[m
[32m+[m[32mburnIn = burnIn,[m
[32m+[m[32mthin = thin)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mresponse = "Yield"[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mweights = NULL[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32momic1_data = NULL,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mres_summary_stat <- summary_statistics(mod = mod)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32momic1_data = NULL,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mgc()[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mBIN = mod$output_files_names[grepl("bin", mod$output_files_names)][m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mlist.files(pattern=files_key)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mBIN = mod$output_files_names[grepl("bin", mod$output_files_names)][m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mGS_model= "BayesB"[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mgetwd()[m
[32m+[m[32mcurrent_date_time = as.character(Sys.time())[m
[32m+[m[32mfiles_key = gsub(" ", "", current_date_time)[m
[32m+[m[32m### Create key to remove all reduant files from the wkdir[m
[32m+[m[32mfiles_key = strsplit(files_key, "\\.")[[1]][2][m
[32m+[m[32mfm <- BGLR::BGLR([m
[32m+[m[32my=object[, response],[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mnIter = bayes_para$nIter,[m
[32m+[m[32mburnIn =  bayes_para$burnIn,[m
[32m+[m[32mthin =  bayes_para$thin,[m
[32m+[m[32mverbose = FALSE,[m
[32m+[m[32msaveAt =files_key)[m
[32m+[m[32mobject = pheno_clean[[1]][m
[32m+[m[32mfm <- BGLR::BGLR([m
[32m+[m[32my=object[, response],[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mnIter = bayes_para$nIter,[m
[32m+[m[32mburnIn =  bayes_para$burnIn,[m
[32m+[m[32mthin =  bayes_para$thin,[m
[32m+[m[32mverbose = FALSE,[m
[32m+[m[32msaveAt =files_key)[m
[32m+[m[32m## Get the name of all files stored by GBLR using the current name and time the analysis was performed[m
[32m+[m[32moutput_files_names = list.files(pattern=files_key)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32momic1_data = NULL,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32momic1_data = NULL,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mTA$pheno_data[m
[32m+[m[32mETA$pheno_data[m
[32m+[m[32mETA$ETA[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = NULL,[m
[32m+[m[32momic1_data = omic1_clean,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mfm <- BGLR::BGLR([m
[32m+[m[32my=object[, response],[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mnIter = bayes_para$nIter,[m
[32m+[m[32mburnIn =  bayes_para$burnIn,[m
[32m+[m[32mthin =  bayes_para$thin,[m
[32m+[m[32mverbose = FALSE,[m
[32m+[m[32msaveAt =files_key)[m
[32m+[m[32m## Get the name of all files stored by GBLR using the current name and time the analysis was performed[m
[32m+[m[32moutput_files_names = list.files(pattern=files_key)[m
[32m+[m[32moutput = list(model = fm, output_files_names = files_key)[m
[32m+[m[32mnames(output) <- c("model", "output_files_names")[m
[32m+[m[32moutput$output_files_names[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32mETA  <-  ETA_compiler_bayes([m
[32m+[m[32mrandom = random,[m
[32m+[m[32mGS_model = GS_model,[m
[32m+[m[32mpheno_data = pheno_clean[[1]],[m
[32m+[m[32mgeno_data = geno_clean,[m
[32m+[m[32mgen_name = gen_name)[m
[32m+[m[32mmod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,[m
[32m+[m[32mresponse = response,[m
[32m+[m[32mweights = weights,[m
[32m+[m[32mETA = ETA$ETA,[m
[32m+[m[32mbayes_para = bayes_para,[m
[32m+[m[32mverbose = FALSE[m
[32m+[m[32m)[m
[32m+[m[32mres_model_output <- mod_output_Bayes(mod = mod,[m
[32m+[m[32mgeno_data = NULL,[m
[32m+[m[32momic1_data = omic1_clean,[m
[32m+[m[32momic2_data = NULL,[m
[32m+[m[32momic3_data = NULL[m
[32m+[m[32m)[m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mrm(list = ls()); ls()[m
[32m+[m[32mgc()[m
[32m+[m[32mcat('\014')[m
[32m+[m[32msetwd("D:/PredictProR/PredictProR")[m
[32m+[m[32mload('WheatPhenoGeno.Rdata')[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR"), ])[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[32m+[m[32mmod$output_files_names[m
[32m+[m[32mlibrary(PredictProR)[m
[32m+[m[32msik = PredictProR::model_execute(pheno_data = pheno,[m
[32m+[m[32mgeno_data = Geno.data,[m
[32m+[m[32mGS_model = "BRR",[m
[32m+[m[32mgen_name = "GID",[m
[32m+[m[32mresponse = "Yield",[m
[32m+[m[32mrandom = ~GID)[m
[1mdiff --git a/049817mu.dat b/049817mu.dat[m
[1mnew file mode 100644[m
[1mindex 0000000..58848b5[m
[1m--- /dev/null[m
[1m+++ b/049817mu.dat[m
[36m@@ -0,0 +1,40 @@[m
[32m+[m[32m4.736122[m
[32m+[m[32m4.728996[m
[32m+[m[32m4.729454[m
[32m+[m[32m4.72706[m
[32m+[m[32m4.72531[m
[32m+[m[32m4.727752[m
[32m+[m[32m4.724302[m
[32m+[m[32m4.723962[m
[32m+[m[32m4.725295[m
[32m+[m[32m4.795576[m
[32m+[m[32m4.72814[m
[32m+[m[32m4.724376[m
[32m+[m[32m4.720094[m
[32m+[m[32m4.752382[m
[32m+[m[32m4.702455[m
[32m+[m[32m4.705013[m
[32m+[m[32m4.743193[m
[32m+[m[32m4.7051[m
[32m+[m[32m4.759175[m
[32m+[m[32m4.764916[m
[32m+[m[32m4.711221[m
[32m+[m[32m4.802366[m
[32m+[m[32m4.737439[m
[32m+[m[32m4.755238[m
[32m+[m[32m4.737584[m
[32m+[m[32m4.755642[m
[32m+[m[32m4.72004[m
[32m+[m[32m4.717995[m
[32m+[m[32m4.731267[m
[32m+[m[32m4.760735[m
[32m+[m[32m4.734419[m
[32m+[m[32m4.668416[m
[32m+[m[32m4.751231[m
[32m+[m[32m4.761734[m
[32m+[m[32m4.725341[m
[32m+[m[32m4.72162[m
[32m+[m[32m4.755956[m
[32m+[m[32m4.680149[m
[32m+[m[32m4.764139[m
[32m+[m[32m4.700551[m
[1mdiff --git a/049817varE.dat b/049817varE.dat[m
[1mnew file mode 100644[m
[1mindex 0000000..5ad4c1c[m
[1m--- /dev/null[m
[1m+++ b/049817varE.dat[m
[36m@@ -0,0 +1,40 @@[m
[32m+[m[32m0.1986176[m
[32m+[m[32m0.1689798[m
[32m+[m[32m0.2034477[m
[32m+[m[32m0.2005587[m
[32m+[m[32m0.1817149[m
[32m+[m[32m0.2113557[m
[32m+[m[32m0.181505[m
[32m+[m[32m0.190501[m
[32m+[m[32m0.1548511[m
[32m+[m[32m0.2189274[m
[32m+[m[32m0.1945712[m
[32m+[m[32m0.2101595[m
[32m+[m[32m0.1914088[m
[32m+[m[32m0.1990207[m
[32m+[m[32m0.1875685[m
[32m+[m[32m0.2150308[m
[32m+[m[32m0.1682771[m
[32m+[m[32m0.2274265[m
[32m+[m[32m0.189807[m
[32m+[m[32m0.1653081[m
[32m+[m[32m0.1967827[m
[32m+[m[32m0.1958917[m
[32m+[m[32m0.1898629[m
[32m+[m[32m0.1760104[m
[32m+[m[32m0.1895345[m
[32m+[m[32m0.163163[m
[32m+[m[32m0.161379[m
[32m+[m[32m0.1866997[m
[32m+[m[32m0.2096428[m
[32m+[m[32m0.1845779[m
[32m+[m[32m0.2030841[m
[32m+[m[32m0.2101586[m
[32m+[m[32m0.1871155[m
[32m+[m[32m0.1903905[m
[32m+[m[32m0.2184098[m
[32m+[m[32m0.2022308[m
[32m+[m[32m0.1900385[m
[32m+[m[32m0.2004087[m
[32m+[m[32m0.2035058[m
[32m+[m[32m0.1896249[m
[1mdiff --git a/192548mu.dat b/192548mu.dat[m
[1mnew file mode 100644[m
[1mindex 0000000..9602c1f[m
[1m--- /dev/null[m
[1m+++ b/192548mu.dat[m
[36m@@ -0,0 +1,40 @@[m
[32m+[m[32m4.715246[m
[32m+[m[32m4.743216[m
[32m+[m[32m4.69134[m
[32m+[m[32m4.761339[m
[32m+[m[32m4.758798[m
[32m+[m[32m4.725562[m
[32m+[m[32m4.727521[m
[32m+[m[32m4.748141[m
[32m+[m[32m4.763189[m
[32m+[m[32m4.72592[m
[32m+[m[32m4.742032[m
[32m+[m[32m4.750077[m
[32m+[m[32m4.762915[m
[32m+[m[32m4.746606[m
[32m+[m[32m4.743712[m
[32m+[m[32m4.768981[m
[32m+[m[32m4.741873[m
[32m+[m[32m4.719529[m
[32m+[m[32m4.722041[m
[32m+[m[32m4.706066[m
[32m+[m[32m4.755701[m
[32m+[m[32m4.714551[m
[32m+[m[32m4.727197[m
[32m+[m[32m4.723959[m
[32m+[m[32m4.731289[m
[32m+[m[32m4.739909[m
[32m+[m[32m4.759664[m
[32m+[m[32m4.74482[m
[32m+[m[32m4.731532[m
[32m+[m[32m4.749909[m
[32m+[m[32m4.711192[m
[32m+[m[32m4.752917[m
[32m+[m[32m4.75272[m
[32m+[m[32m4.758423[m
[32m+[m[32m4.755525[m
[32m+[m[32m4.765524[m
[32m+[m[32m4.738492[m
[32m+[m[32m4.714984[m
[32m+[m[32m4.765148[m
[32m+[m[32m4.743139[m
[1mdiff --git a/192548varE.dat b/192548varE.dat[m
[1mnew file mode 100644[m
[1mindex 0000000..b7dcf4b[m
[1m--- /dev/null[m
[1m+++ b/192548varE.dat[m
[36m@@ -0,0 +1,40 @@[m
[32m+[m[32m0.1852549[m
[32m+[m[32m0.1742704[m
[32m+[m[32m0.1728966[m
[32m+[m[32m0.1859979[m
[32m+[m[32m0.1862347[m
[32m+[m[32m0.1807402[m
[32m+[m[32m0.2275[m
[32m+[m[32m0.1810424[m
[32m+[m[32m0.1806946[m
[32m+[m[32m0.1705448[m
[32m+[m[32m0.1746602[m
[32m+[m[32m0.1732439[m
[32m+[m[32m0.1977115[m
[32m+[m[32m0.2066327[m
[32m+[m[32m0.1954405[m
[32m+[m[32m0.2068878[m
[32m+[m[32m0.1998427[m
[32m+[m[32m0.1626219[m
[32m+[m[32m0.1792248[m
[32m+[m[32m0.2002974[m
[32m+[m[32m0.2271469[m
[32m+[m[32m0.1970145[m
[32m+[m[32m0.2268265[m
[32m+[m[32m0.2103631[m
[32m+[m[32m0.1800867[m
[32m+[m[32m0.1772621[m
[32m+[m[32m0.1718[m
[32m+[m[32m0.2134375[m
[32m+[m[32m0.1810974[m
[32m+[m[32m0.194792[m
[32m+[m[32m0.1815798[m
[32m+[m[32m0.2166376[m
[32m+[m[32m0.2003887[m
[32m+[m[32m0.18614[m
[32m+[m[32m0.2003415[m
[32m+[m[32m0.2046958[m
[32m+[m[32m0.1589428[m
[32m+[m[32m0.1785318[m
[32m+[m[32m0.18742[m
[32m+[m[32m0.1948942[m
[1mdiff --git a/DESCRIPTION b/DESCRIPTION[m
[1mnew file mode 100644[m
[1mindex 0000000..4263c33[m
[1m--- /dev/null[m
[1m+++ b/DESCRIPTION[m
[36m@@ -0,0 +1,24 @@[m
[32m+[m[32mPackage: PredictProR[m
[32m+[m[32mType: Package[m
[32m+[m[32mTitle: What the Package Does (Title Case)[m
[32m+[m[32mVersion: 0.1.0[m
[32m+[m[32mAuthor: Who wrote it[m
[32m+[m[32mMaintainer: The package maintainer <yourself@somewhere.net>[m
[32m+[m[32mDescription: More about what it does (maybe more than one line)[m
[32m+[m[32m    Use four spaces when indenting paragraphs within the Description.[m
[32m+[m[32mLicense: What license is it under?[m
[32m+[m[32mEncoding: UTF-8[m
[32m+[m[32mLazyData: true[m
[32m+[m[32mEnhances:[m[41m [m
[32m+[m[32m    asreml[m
[32m+[m[32mImports:[m[41m [m
[32m+[m[32m    ASRgenomics,[m
[32m+[m[32m    BGLR,[m
[32m+[m[32m    data.table,[m
[32m+[m[32m    doBy,[m
[32m+[m[32m    Matrix,[m
[32m+[m[32m    matrixcalc,[m
[32m+[m[32m    sommer,[m
[32m+[m[32m    stats,[m
[32m+[m[32m    stringr,[m
[32m+[m[32m    utils[m
[1mdiff --git a/DataPre.R b/DataPre.R[m
[1mnew file mode 100644[m
[1mindex 0000000..d523566[m
[1m--- /dev/null[m
[1m+++ b/DataPre.R[m
[36m@@ -0,0 +1,134 @@[m
[32m+[m[32mrm(list = ls()); ls()[m
[32m+[m[32mgc()[m
[32m+[m[32mlibrary(lme4GS)[m
[32m+[m[32mdata(wheat599)[m
[32m+[m[32mrm(wheat.Pedigree)[m
[32m+[m[32mGD = unique(wheat.Pheno$GID)[m
[32m+[m
[32m+[m[32mset.seed(123)[m
[32m+[m[32mGDD = GD[GD%in%sample(GD, 200, replace = F)][m
[32m+[m
[32m+[m[32mwheat.Pheno = wheat.Pheno[wheat.Pheno$Rep==1, ][m
[32m+[m
[32m+[m[32m#wheat.Pheno = wheat.Pheno[wheat.Pheno$Env==1, ][m
[32m+[m
[32m+[m[32mwheat.Pheno$PlantHeight = rnorm(nrow(wheat.Pheno), mean = 78, sd= 35)[m
[32m+[m
[32m+[m[32mwheat.Pheno$DaysToMat = rnorm(nrow(wheat.Pheno), mean = 56, sd= 10)[m
[32m+[m
[32m+[m[32mwheat.Pheno$Mositure_Content = rnorm(nrow(wheat.Pheno), mean = 21, sd= 6)[m
[32m+[m
[32m+[m[32mwheat.Pheno = wheat.Pheno[wheat.Pheno$GID%in%GDD, ][m
[32m+[m
[32m+[m
[32m+[m[32m#wheat.Pheno = wheat.Pheno[which(wheat.Pheno$GID%in%sample(GD, 200, replace = F)), ][m
[32m+[m
[32m+[m[32m#wheat.Pheno = wheat.Pheno[!wheat.Pheno$Env%in%c(1, 2), ][m
[32m+[m[32mwheat.X = wheat.X[rownames(wheat.X)%in%unique(wheat.Pheno$GID), ][m
[32m+[m[32mpheno = wheat.Pheno[m
[32m+[m[32mGeno_data = wheat.X[m
[32m+[m
[32m+[m[32mgenotype = "GID"[m
[32m+[m
[32m+[m[32mpheno$Env = as.factor(pheno$Env)[m
[32m+[m[32mpheno$Rep = as.factor(pheno$Rep)[m
[32m+[m[32mpheno$GID = as.factor(pheno$GID)[m
[32m+[m
[32m+[m
[32m+[m[32mrm(list = ls()); ls()[m
[32m+[m[32mgc()[m
[32m+[m[32mcat('\014')[m
[32m+[m[32msetwd("D:/PredictProR/PredictProR")[m
[32m+[m[32mload('WheatPhenoGeno.Rdata')[m
[32m+[m
[32m+[m[32mgenotype = "GID"[m
[32m+[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR","F5I", "B5I"), ])[m
[32m+[m
[32m+[m[32mpheno = droplevels(pheno.data[pheno.data$Env%in%c("B2IR"), ])[m
[32m+[m
[32m+[m[32mOmic2 = as.matrix(COP)[m
[32m+[m
[32m+[m[32mresponse = "Yield"[m
[32m+[m
[32m+[m[32mrandom = ~GID[m
[32m+[m
[32m+[m[32m#fixed = ~Env[m
[32m+[m
[32m+[m[32mZE<-model.matrix(~factor(pheno$Env)-1)[m
[32m+[m[32mdim(ZE)[m
[32m+[m[32mdim(Geno.data)[m
[32m+[m
[32m+[m[32m#Environment–marker design matrix[m
[32m+[m[32mXEM = model.matrix(~0+Geno.data:XE)[m
[32m+[m[32m#pheno = droplevels(pheno[pheno$Env%in%c("B2IR"), ])[m
[32m+[m
[32m+[m
[32m+[m[32mpheno$GID = as.character(pheno$GID)[m
[32m+[m[32mpheno$Env = as.character(pheno$Env)[m
[32m+[m[32mstr(pheno)[m
[32m+[m
[32m+[m[32mresponse = c("Yield")[m
[32m+[m[32m#"B2IR","F5I", "EHT", "B5I"[m
[32m+[m[32m# [1] B5I   BLHT  B2IR  F5I   FDRIP EHT[m
[32m+[m[32m# Levels: B2IR B5I BLHT EHT F5I FDRIP[m
[32m+[m
[32m+[m[32m# library(PredictiveAnalytic)[m
[32m+[m[32m#[m
[32m+[m[32m# GSModels <- PredictiveAnalytic::GS_mod([m
[32m+[m[32m#     pheno=pheno,[m
[32m+[m[32m#     pheno_train = NULL,[m
[32m+[m[32m#     pheno_test = NULL,[m
[32m+[m[32m#     Gmatrix= NULL,[m
[32m+[m[32m#     Omic2_Relationship_Matrix = NULL,[m
[32m+[m[32m#     Omic3_Relationship_Matrix = NULL,[m
[32m+[m[32m#     Geno_data = Geno_data,[m
[32m+[m[32m#     Omic2_data = NULL,[m
[32m+[m[32m#     Omic3_data = NULL,[m
[32m+[m[32m#     train_Geno_data = NULL,[m
[32m+[m[32m#     test_Geno_data = NULL,[m
[32m+[m[32m#     train_Omic2_data = NULL,[m
[32m+[m[32m#     test_Omic2_data = NULL,[m
[32m+[m[32m#     train_Omic3_data = NULL,[m
[32m+[m[32m#     test_Omic3_data = NULL,[m
[32m+[m[32m#     train_set = NULL,[m
[32m+[m[32m#     test_set = NULL,[m
[32m+[m[32m#     Gmatrix_method = "Yang",[m
[32m+[m[32m#     Kernel_matrix_method = NULL,[m
[32m+[m[32m#     engine= "BGLR",[m
[32m+[m[32m#     response="Yield",[m
[32m+[m[32m#     genotype="GID",[m
[32m+[m[32m#     cova=NULL,[m
[32m+[m[32m#     fixed= NULL,[m
[32m+[m[32m#     random=~GID+ GID:Env,[m
[32m+[m[32m#     heter_resid=FALSE
