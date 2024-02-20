#' Title
#'
#' @param GS_model
#' @param res_model_output
#' @param res_summary_stat
#' @param res_plot
#' @param system_database
#'
#' @return
#' @export
#'
#' @examples
results_handling <-  function(GS_model = NULL,
                              res_model_output = NULL,
                              res_summary_stat = NULL,
                              res_plot = NULL,
                              system_database = TRUE){
  ### This part is for GBLUP_BRR
  # if(exists("GS_modeluse")){
  #   GS_model <-  GS_modeluse
  # }

  if(!exists("res_plot")){
    res_plot = NULL
  }

  ############

  zipMMatrixModelReady <- function(pathout, geno_omic_files) {
    is_geno <- grepl("Geno", geno_omic_files, ignore.case = TRUE)
    is_omic <- grepl("Omic", geno_omic_files, ignore.case = TRUE)

    zip_name <- if (any(is_geno) && any(is_omic)) {
      "Geno_Omics_Clean_data.zip"
    } else if (any(is_geno)) {
      "Geno_Clean_data.zip"
    } else if (any(is_omic)) {
      "Omic_Clean_data.zip"
    } else {
      "X_variables_Clean_data.zip"
    }

    zip(file.path(pathout, zip_name), files = geno_omic_files, flags = "-q")
  }

  saveOutput <- function(res_model_output, res_summary_stat,  pathout, GS_model) {
    for (i in 1:length(res_model_output)) {

      for(i in 1:length(res_model_output)){
        ### check which output is a list
        if(length(class(res_model_output[[i]]))==1 && !inherits(res_model_output[[i]], "list")){
          if(names(res_model_output)[i]=="Variance_components"){
            write.csv(res_model_output[[i]],
                      file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                      row.names = TRUE)

          } else if(names(res_model_output)[i]=="Asreml_model"){
            base::saveRDS(res_model_output[[i]],
                          "asreml_model.RData")

          } else if(names(res_model_output)[i]=="trained_model"){
            base::saveRDS(res_model_output[[i]],
                          "AI_trained_model.RData")

          }else{

            write.csv(res_model_output[[i]],
                      file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                      row.names = FALSE)
          }
        } else if (length(class(res_model_output[[i]])) > 1 &&
                   "asreml.predict" %in% class(res_model_output[[i]]) &&
                   "data.frame" %in% class(res_model_output[[i]])) {

          write.csv(res_model_output[[i]],
                    file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                    row.names = TRUE)


        } else if (length(class(res_model_output[[i]])) > 1 &&
                   "matrix" %in% class(res_model_output[[i]]) &&
                   "array" %in% class(res_model_output[[i]])) {

          write.csv(res_model_output[[i]],
                    file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                    row.names = TRUE)

        } else {

          if(length(class(res_model_output[[i]]))==1 && inherits(res_model_output[[i]],"list")){
            if(names(res_model_output)[i]=="M_matrix_model_ready"){
              for(j in 1:length(res_model_output[[i]])){
                write.csv(res_model_output[[i]][[j]],
                          file.path(pathout, paste(names(res_model_output[[i]][j]), "csv", sep = ".")),
                          row.names = TRUE)

              }
            } else if(names(res_model_output)[i]=="Covariance" | names(res_model_output)[i]=="Correlation"){
              for(j in 1:length(res_model_output[[i]])){
                write.csv(res_model_output[[i]][[j]],
                          file.path(pathout, paste(names(res_model_output[[i]][j]), "csv", sep = ".")),
                          row.names = TRUE)

              }

            } else if(names(res_model_output)[i]=="Asreml_model"){
              base::saveRDS(res_model_output[[i]],
                            "asreml_model.RData")

            } else if(names(res_model_output)[i]=="trained_model"){
              base::saveRDS(res_model_output[[i]],
                            "AI_trained_model.RData")

            } else {
              for(j in 1:length(res_model_output[[i]])){
                write.csv(res_model_output[[i]][[j]],
                          file.path(pathout, paste(names(res_model_output[[i]][j]), "csv", sep = ".")),
                          row.names = FALSE)

              }

            }

          }

          # if(class(res_model_output[[i]])=="asreml"){
          #     if(names(res_model_output)[i]=="Asreml_model"){
          #         base::saveRDS(res_model_output[[i]],
          #                       "asreml_model.RData")
          #     }
          #
          # }



        }


      }

      # output <- res_model_output[[i]]
      # filename <- paste(names(output), "csv", sep = ".")
      #
      # if (class(output) == "list") {
      #     for (j in seq_along(output)) {
      #         write.csv(output[[j]], file.path(pathout, filename[j]), row.names = TRUE)
      #     }
      # } else {
      #     write.csv(output, file.path(pathout, filename), row.names = !is.null(rownames(output)))
      # }
    }
    if(!is.null(res_summary_stat)){
    for(s in 1:length(res_summary_stat)){

      if(!inherits(res_summary_stat[[s]], "list")){
        write.csv(res_summary_stat[[s]],
                  file.path(pathout,paste(names(res_summary_stat)[s], "csv", sep = ".")),
                  row.names = FALSE)

      }


    }
  }
  }

  processMMatrixModelReady <- function(pathout) {
    files_in_directory <- list.files()
    geno_omic_files <- grep("ready", files_in_directory, value = TRUE)

    if (length(geno_omic_files) > 0) {
      zipMMatrixModelReady(pathout, geno_omic_files)
      unlink(geno_omic_files)
    }
  }

  saveOutputAndZip <- function(res_model_output, res_summary_stat) {
    mainDir <- getwd()
    systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
    systime <- gsub("[-: ]", "_", systime)
    subDir <- paste("output", systime, sep = "_")
    subDir2 <- paste("outputNew", systime, sep = "_")

    if (dir.exists(file.path(mainDir, subDir))) {
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
    } else {
      dir.create(file.path(mainDir, subDir2))
      setwd(file.path(mainDir, subDir2))
    }

    pathout <- getwd()

    saveOutput(res_model_output, res_summary_stat,  pathout, GS_model)
    processMMatrixModelReady(pathout)

    setwd(mainDir)

    return(NULL)
  }

  processData <- function(GS_model,
                          res_model_output,
                          res_summary_stat,
                          res_plot,
                          system_database) {

    # if (exists("GS_modeluse")) {
    #   GS_model <- GS_modeluse
    # }

    if (GS_model %in% c("Xgboost")) {
      output <- list(model_results = res_model_output,
                     summary_statistic = res_summary_stat)
      if (exists("test_set_")) {
        output$res_plot <- res_plot
      }
      return(output)

    } else if (GS_model %in% c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS", "GBLUP_BRR", "GBLUP")) {

      if (isFALSE(system_database)) {
        output <- saveOutputAndZip(res_model_output, res_summary_stat)
      } else {
        output <- list(model_results = res_model_output,
                       summary_statistic = res_summary_stat)
      }
      return(output)

    } else {

      output <- list(model_results = res_model_output, summary_statistic = res_summary_stat)
      return(output)
    }
  }



  return(processData(GS_model = GS_model,
                     res_model_output = res_model_output,
                     res_summary_stat = res_summary_stat,
                     res_plot = res_plot,
                     system_database = system_database))

  #out <-  list(res_model_output, res_summary_stat)
  #return(out)
} ## end of function

