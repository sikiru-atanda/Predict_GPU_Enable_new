
working_dir_plot <- function(output_file_name = "cross_validation"){
  mainDir <- getwd()
  systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  systime <- gsub("[-: ]", "_", systime)
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

  return(list(pathout= pathout,
              mainDir = mainDir))
  # saveOutput(res_model_output, res_summary_stat,  pathout, GS_model, res_plot)
  # processMMatrixModelReady(pathout)

  #setwd(mainDir)

}
