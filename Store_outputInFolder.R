# https://www.r-bloggers.com/2021/05/working-with-files-and-folders-in-r-ultimate-guide/
## Get the existing working directory

## Get system time

storeoutput_Folder <- function(Results,...){

  Initapth = getwd()
systime = Sys.Date()
systime = gsub("-", "_", systime)
## Create path
pathout = paste(Initapth, paste("Results", systime, sep = "_"), sep = "/")
## Create alternative path if the previous one already exist. Though not likely
pathout2 = paste(Initapth, paste("Results2", systime, sep = "_"), sep = "/")

ifelse(!dir.exists(pathout), dir.create(pathout), dir.create(pathout2))

### create output folder within the working directory
dir.create(pathout, showWarnings = FALSE)

### set the working directory to the output folder
setwd(pathout)



return(c(write.csv(Result, paste("OptimizedTRN", "NTrn.Optmize", ".csv", sep = "_")),
         write.csv(COR, paste("PredACC","NTrn.Optmize", "Gmatrix.method", ".csv", sep = "_"))))

}

########
if(system_use!="database"){

Initapth = getwd()
systime = Sys.Date()
systime = gsub("-", "_", systime)
## Create path
pathout = paste(Initapth, paste("Results", systime, sep = "_"), sep = "/")
## Create alternative path if the previous one already exist. Though not likely
pathout2 = paste(Initapth, paste("Results2", systime, sep = "_"), sep = "/")

ifelse(!dir.exists(pathout), dir.create(pathout, showWarnings = FALSE),
       dir.create(pathout2, showWarnings = FALSE))

# ### create output folder within the working directory
# dir.create(pathout, showWarnings = FALSE)

### set the working directory to the output folder
setwd(pathout)


label_index <- which(sapply(AA$model_results, function(x) class(x)=='list'))

names_res <- names(AA$model_results)

names_index = grep("estimated_breeding_value", names_res, ignore.case = TRUE)

label_index <- as.double(label_index)

#label_index <- as.double(label_index[-length(label_index)])

if(!is.null(label_index)){

  c(write.csv(AA$model_results$Predicted_value, paste("Predicted_value", "csv", sep = "."), row.names = FALSE),
    write.csv(AA$model_results$Total_estimated_breeding_value, paste("Total_estimated_breeding_value", "csv", sep = "."), row.names = FALSE),
    write.csv(AA$model_results$Variance_components, paste("Variance_components", "csv", sep = ".")),
    write.csv(AA$summary_statistic$Statics_summary, paste("summary_statistic", "csv", sep = "."), row.names = FALSE),
    for (i in 1:length(label_index)) {
      ###
      if(names(AA$model_results[label_index[i]])=="M_matrix_model_ready"){
        file_names = paste0(names(AA$model_results[[label_index[i]]]), ".txt") }
      ###

      ###
      if(names(AA$model_results[label_index[i]])=="M_matrix_model_ready"){
        for (k in 1:length(names(AA$model_results[[label_index[i]]]))) {
          write.csv(AA$model_results[[i]][k], paste(names(AA$model_results[[label_index[i]]][k]), "txt", sep = "."), row.names = FALSE)
        }
      } else {

        for (k in 1:length(AA$model_results[[i]])) {

          write.csv(AA$model_results[[i]][k], paste(names(AA$model_results[[i]][k]), "csv", sep = "."), row.names = FALSE)
        }

      }

      ##
      if(names(AA$model_results[label_index[i]])=="M_matrix_model_ready"){
        # Read the 2 CSV file names from working directory
        # Zip_Files <- list.files(path = getwd(), pattern = ".txt$")
        #
        # zip(zipfile = "TestZip", files = Zip_Files, flags = " a -tzip",
        #     zip = "C:\\Program Files\\7-Zip\\7Z")

        #zip(file.path(pathout, "zipped_files.zip"), files = file.path(pathout, file_names))

        if(length(grep("Geno", file_names, ignore.case = TRUE))!=0 & length(grep("Omic", file_names, ignore.case = TRUE))!=0){
          ## # flags="-q" silence the output message. check what appropriate in Mac
          zip(file.path(pathout, "Geno_Omics_Clean_data.zip"), files = file_names, flags="-q")

        } else if (length(grep("Geno", file_names, ignore.case = TRUE))!=0){

          zip(file.path(pathout, "Geno_Clean_data.zip"), files = file_names, flags="-q")
        } else if(length(grep("Omic", file_names, ignore.case = TRUE))!=0) {
          zip(file.path(pathout, "Omics_Clean_data.zip"), files = file_names, flags="-q")
        } else {
          if (length(grep("Omic", file_names, ignore.case = TRUE))==0){
            zip(file.path(pathout, "X_variables_Clean_data.zip"), files = file_names, flags="-q")

          }

        }

        #zip(file.path(pathout, "zipped_files.zip"), files = file_names)

        unlink(file_names)
      }
      ##

    }
  )



} else {
  c(write.csv(AA$model_results$Predicted_value, paste("Predicted_value", "csv", sep = "."), row.names = FALSE),
    write.csv(AA$model_results$Total_estimated_breeding_value, paste("Total_estimated_breeding_value", "csv", sep = "."), row.names = FALSE),
    write.csv(AA$model_results$Variance_components, paste("Variance_components", "csv", sep = ".")),
    write.csv(AA$model_results$Coefficients, paste("Coefficients", "csv", sep = "."), row.names = FALSE),
    write.csv(AA$model_results[names_index], paste(names_res[names_index], "csv", sep = "."), row.names = FALSE),
    write.csv(AA$summary_statistic$Statics_summary, paste("summary_statistic", "csv", sep = "."), row.names = FALSE)




  )



}

}else{

  output <- list(res_model_output, res_summary_stat)

  names(output) <- c('model_results', 'summary_statistic')

}

MM = AA$summary_statistic$Statics_summary
