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
