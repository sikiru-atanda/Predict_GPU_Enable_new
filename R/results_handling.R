#' Title
#'
#' @param res_plot
#' @param plot_filename
#' @param plot_extension
#' @param plot_width
#' @param plot_height
#' @param plot_units
#' @param plot_dpi
#'
#' @return
#' @export
#'
#' @examples
save_ggplot <- function(res_plot,
                        plot_filename = "trait",
                        plot_extension = "pdf",
                        plot_width = 17,
                        plot_height = 12,
                        plot_units = "in",
                        plot_dpi = 300) {


  # save plot plot with user defined name and parameters
  if(is.null(plot_filename)) plot_filename <- "trait"
  plot_name <- paste(plot_filename,plot_extension, sep = ".")

       ggplot2::ggsave(filename = plot_name,
                       plot = res_plot,
                       width = plot_width,
                       height = plot_height,
                       units = plot_units,
                       dpi = plot_dpi)


}


#' Title
#'
#' @param pathout
#' @param geno_omic_files
#'
#' @return
#' @export
#'
#' @examples
zipMMatrixModelReady <- function(pathout,
                                 geno_omic_files) {
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



#' Title
#'
#' @param GS_model
#' @param res_model_output
#' @param res_summary_stat
#' @param res_plot
#' @param res_plot_mean
#' @param res_plot_result_diagnostic
#' @param test_diagonistic_plots
#' @param res_mod_results_cv_per_trait_model
#' @param res_plot_result_diagnostic_cv_only
#' @param cv_results_processed
#' @param geno_qc_stat
#' @param system_database
#' @param plot_filename
#' @param plot_extension
#' @param plot_width
#' @param plot_height
#' @param plot_units
#' @param plot_dpi
#'
#' @return
#' @export
#'
#' @examples
results_handling <-  function(GS_model = NULL,
                              res_model_output = NULL,
                              res_summary_stat = NULL,
                              res_plot = NULL,
                              res_plot_mean = NULL,
                              res_plot_result_diagnostic = NULL,
                              test_diagonistic_plots = NULL,
                              res_mod_results_cv_per_trait_model = NULL,
                              res_plot_result_diagnostic_cv_only = NULL,
                              cv_results_processed = NULL,
                              geno_qc_stat = NULL,
                              system_database = TRUE,
                              plot_filename = "trait",
                              #Plot_name_result_diagnostic = "result_diagnostic",
                              plot_extension = "pdf",
                              plot_width = 17,
                              plot_height = 12,
                              plot_units = "in",
                              plot_dpi = 300
                              ){

  saveOutput <- function(res_model_output,
                         res_summary_stat,
                         pathout,
                         GS_model,
                         res_plot, res_plot_mean,
                         res_plot_result_diagnostic,
                         test_diagonistic_plots,
                         res_plot_result_diagnostic_cv_only) {
    #browser()
    if(!is.null(res_model_output)){
      if ("bayes_model" %in% names(res_model_output)) {
        res_model_output <- res_model_output[names(res_model_output) != "bayes_model"]
      }

      if ("bayes_result" %in% names(res_model_output)) {
        res_model_output <- res_model_output[["bayes_result"]]
      }

      if ("diagnostic_plots" %in% names(res_model_output)) {
        res_model_output <- res_model_output[names(res_model_output) != "diagnostic_plots"]
      }
    #for (i in 1:length(res_model_output)) {

      for(i in 1:length(res_model_output)){
        ### check which output is a list
        if(length(class(res_model_output[[i]]))==1 && !inherits(res_model_output[[i]], "list")){
          if(names(res_model_output)[i]=="Variance_components"){
            write.csv(res_model_output[[i]],
                      file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                      row.names = TRUE)

          } else if(names(res_model_output)[i]=="bayes_model"){
            base::saveRDS(res_model_output[[i]],
                          paste(GS_model, "Bayes_model.RData", sep = "_"))

          } else if(names(res_model_output)[i]=="Asreml_model"){
            base::saveRDS(res_model_output[[i]],
                          "asreml_model.RData")

          } else if(names(res_model_output)[i]=="trained_model"){
            base::saveRDS(res_model_output[[i]],
                          paste(GS_model, "trained_model.RData", sep = "_"))

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
                if(!inherits(res_model_output[[i]][[j]], "gtable")){
                write.csv(res_model_output[[i]][[j]],
                          file.path(pathout, paste(names(res_model_output[[i]][j]), "csv", sep = ".")),
                          row.names = FALSE)

                }

              }

            }

          }

        }


      }

    #}

    }

    if(!is.null(res_summary_stat)){
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
  ####
    if(!is.null(res_plot)){
      save_ggplot(res_plot = res_plot,
                  plot_filename = plot_filename,
                  plot_extension = plot_extension,
                  plot_width = plot_width,
                  plot_height = plot_height,
                  plot_units = plot_units,
                  plot_dpi = plot_dpi)
    }

    if(!is.null(res_plot_mean)){
      save_ggplot(res_plot = res_plot_mean,
                  plot_filename = paste(plot_filename, "mean", sep = "_"),
                  plot_extension = plot_extension,
                  plot_width = plot_width,
                  plot_height = plot_height,
                  plot_units = plot_units,
                  plot_dpi = plot_dpi)
    }

    if(!is.null(test_diagonistic_plots)){
      if(inherits(test_diagonistic_plots, "gtable")){
      save_ggplot(res_plot = test_diagonistic_plots,
                  plot_filename = paste("test_diagonistic_plots", "GS_model", sep = "_"),
                  plot_extension = plot_extension,
                  plot_width = plot_width,
                  plot_height = plot_height,
                  plot_units = plot_units,
                  plot_dpi = plot_dpi)
      }

    }


    if(!is.null(res_plot_result_diagnostic)){

        models <-   names(res_plot_result_diagnostic$predicted_vs_observed_plots)

        for (mod in models) {

          combined_plot <- res_plot_result_diagnostic$predicted_vs_observed_plots[[mod]]

          name_plot <- paste("Cross_validation_diagonistic_plots", mod, sep = "_")

if(inherits(combined_plot, "gtable")){
          save_ggplot(res_plot = combined_plot,
                      plot_filename = name_plot,
                      plot_extension = plot_extension,
                      plot_width = plot_width,
                      plot_height = plot_height,
                      plot_units = plot_units,
                      plot_dpi = plot_dpi)

}

        }


      }

    #######
    if(!is.null(res_plot_result_diagnostic_cv_only)){
      traits <- names(res_plot_result_diagnostic_cv_only)

      for (trait in traits) {

       res_plot_cv <-  res_plot_result_diagnostic_cv_only[[trait]]
      models <-   names(res_plot_cv$predicted_vs_observed_plots)

      for (mod in models) {

        combined_plot <- res_plot_cv$predicted_vs_observed_plots[[mod]]

        name_plot <- paste(paste("Cross_validation_diagonistic_plots",  trait, sep = "_"), mod, sep = "_")

        if(inherits(combined_plot, "gtable")){
          save_ggplot(res_plot = combined_plot,
                      plot_filename = name_plot,
                      plot_extension = plot_extension,
                      plot_width = plot_width,
                      plot_height = plot_height,
                      plot_units = plot_units,
                      plot_dpi = plot_dpi)

        }

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

  saveOutputAndZip <- function(res_model_output, res_summary_stat, output_file_name,
                               #pathout,
                               GS_model, res_plot, res_plot_mean,
                               test_diagonistic_plots, res_plot_result_diagnostic,
                               res_plot_result_diagnostic_cv_only) {
    msg <- "\n==================================================\n"
    # Check if any required object is NULL
    # if ((is.null(res_model_output) && is.null(res_summary_stat)) && (is.null(test_diagonistic_plots) && is.null(res_plot_result_diagnostic))) {
    #   cat(paste(msg,paste("The output", "from", GS_model,"is NULL.", "No output will be saved.\n")))
    #   return("Failed: Required objects are NULL")
    # }

    if (is.null(res_model_output) && is.null(res_summary_stat) &&
        is.null(res_plot) && is.null(res_plot_mean) &&
        is.null(test_diagonistic_plots) && is.null(res_plot_result_diagnostic)) {

      return("Failed: Required objects are NULL")
    }

    mainDir <- getwd()
    #systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
    #systime <- gsub("[-: ]", "_", systime)
    systime <- format(Sys.time(), "%m-%d-%Y_%I-%M%p")
    subDir <- paste(output_file_name, systime, sep = "_")
    subDir2 <- paste(paste(output_file_name, "New", sep="_"), systime, sep = "_")

    if (!dir.exists(file.path(mainDir, subDir))) {
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
    } else {
      dir.create(file.path(mainDir, subDir2))
      setwd(file.path(mainDir, subDir2))
    }

    pathout <- getwd()

    saveOutput(res_model_output = res_model_output, res_summary_stat = res_summary_stat,
               pathout = pathout,
               GS_model = GS_model, res_plot = res_plot,
               res_plot_mean = res_plot_mean, test_diagonistic_plots = test_diagonistic_plots,
               res_plot_result_diagnostic = res_plot_result_diagnostic,
               res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only)
    processMMatrixModelReady(pathout)

    setwd(mainDir)

    return("Successful")
  }

  processData <- function(GS_model,
                          res_model_output,
                          res_summary_stat,
                          res_plot,
                          #pathout,
                          res_plot_mean,
                          test_diagonistic_plots,
                          res_plot_result_diagnostic,
                          res_plot_result_diagnostic_cv_only,
                          res_mod_results_cv_per_trait_model,
                          cv_results_processed,
                          system_database,
                          plot_filename,
                          plot_extension,
                          plot_width,
                          plot_height,
                          plot_units,
                          plot_dpi) {

      if (isFALSE(system_database)) {
        output <- saveOutputAndZip(res_model_output = res_model_output,
                                   res_summary_stat = res_summary_stat,
                                   output_file_name = plot_filename,
                                   #pathout = pathout,
                                   GS_model = GS_model,
                                   res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
                                   res_plot = res_plot, res_plot_mean = res_plot_mean,
                                   test_diagonistic_plots = test_diagonistic_plots,
                                   res_plot_result_diagnostic = res_plot_result_diagnostic)
      } else {
        output <- list(model_results = res_model_output,
                       summary_statistic = res_summary_stat,
                       res_plot = res_plot,
                       cv_results_processed = cv_results_processed,
                       res_plot_result_diagnostic = res_plot_result_diagnostic,
                       test_diagonistic_plots = test_diagonistic_plots,
                       res_mod_results_cv_per_trait_model = res_mod_results_cv_per_trait_model
                       )
      }
      return(output)

  }



  return(processData(GS_model = GS_model,
                     res_model_output = res_model_output,
                     res_summary_stat = res_summary_stat,
                     res_plot = res_plot,
                     #pathout = pathout,
                     res_plot_mean = res_plot_mean,
                     test_diagonistic_plots = test_diagonistic_plots,
                     res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
                     res_plot_result_diagnostic = res_plot_result_diagnostic,
                     res_mod_results_cv_per_trait_model = res_mod_results_cv_per_trait_model,
                     cv_results_processed = cv_results_processed,
                     system_database = system_database,
                     plot_filename = plot_filename,
                     #Plot_name_result_diagnostic = Plot_name_result_diagnostic,
                     plot_extension = plot_extension,
                     plot_width = plot_width,
                     plot_height = plot_height,
                     plot_units = plot_units,
                     plot_dpi = plot_dpi))

} ## end of function

