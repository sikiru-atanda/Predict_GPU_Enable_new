

#' Title
#'
#' @param object
#' @param response
#' @param weights
#' @param ETA
#' @param bayes_para
#' @param ...
#' @param verbose
#' @param heter_groups
#' @param core
#' @param test_set_val
#' @param cross_validation
#' @param eval_metrics
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
bayes_mod_execute_crossvalGUDOLD <- function(object = NULL,
                              response = NULL,
                              weights = NULL,
                              ETA = NULL,
                              heter_groups = NULL,
                              bayes_para = NULL,
                              verbose = FALSE,
                              core = NULL,
                              test_set_val =NULL,
                              cross_validation = NULL,
                              eval_metrics = c("Accuracy",
                                               "Mean_Squared_Error",
                                               "Bias",
                                               "Root_Mean_Squared_Error",
                                               "Relative_Squared_Error",
                                               "Mean_Absolute_Error",
                                               "Mean_Absolute_Percent_Error"),
                                                    ...
                                                   )
  {




    ###############################################################
    if(!cross_validation %in%c("Hold_Out",
                               "Stratified_Hold_Out",
                               "Repeated_Hold_Out",
                               "Repeated_Stratified_Hold_Out",
                               "Leave_one_Out")){
      nfolds <- as.double(strsplit(strsplit(names(test_set_val)[1], c("_"))[[1]][2],"")[[1]][1])  # Number of folds

    }

    NRep <- length(test_set_val)

    # eval_metrics = c("Accuracy", "Mean_Squared_Error",
    #             "Bias", "Mean_Squared_Error",
    #             "Relative_Squared_Error","Mean_Absolute_Error")

    if(cross_validation %in%c("K-Folds",
                              "Stratified_K-Folds",
                              "Repeated_K-Folds",
                              "Repeated_Stratified_K-Folds")){
      yHatCV <- matrix(NA, nrow =nfolds, ncol = length(eval_metrics))

      results_eval_metrics_all <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
      rownames(results_eval_metrics_all) = paste("REP",1:NRep, sep = "_")
      colnames(results_eval_metrics_all) = eval_metrics

    }

    if(cross_validation %in%c("Hold_Out",
                              "Stratified_Hold_Out",
                              "Repeated_Hold_Out",
                              "Repeated_Stratified_Hold_Out",
                              "Leave_one_Out")){
      yHatCV <- matrix(NA, nrow =1, ncol = length(eval_metrics))

      results_eval_metrics_all <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
      rownames(results_eval_metrics_all) = paste("REP",1:NRep, sep = "_")
      colnames(results_eval_metrics_all) = eval_metrics


    }

    ########################################

    if (cross_validation %in% c("CV1",
                                "CV2",
                                "Repeated_CV1",
                                "Repeated_CV2")){



      #if (Cross_validation=="CV1" | Cross_validation=="CV2"){

      ENV = as.character(unique(object[, heter_groups]))
      # yHatCV <- matrix(NA, nrow =nfolds, ncol = length(ENV))
      # colnames(yHatCV) = ENV

      yHatCV <- matrix(NA, nrow =nfolds, ncol = length(eval_metrics))
      colnames(yHatCV) = eval_metrics

      #yHatCV_MSE <- matrix(NA, nrow =nfolds, ncol = length(ENV))
      # yHatCV_all <- vector(mode = 'list', length = length(eval_metrics))
      # names(yHatCV_all) <- eval_metrics


      yHatCV_all <- vector(mode = 'list', length = length(ENV))
      names(yHatCV_all) <- ENV

      # for (c in 1:length(eval_metrics)) {
      #
      #   yHatCV_all[[c]] <- yHatCV
      # }


      for (c in 1:length(ENV)) {

        yHatCV_all[[c]] <- yHatCV
      }
      #colnames(yHatCV_MSE) = ENV


      # results_eval_metrics_ <- matrix(NA, nrow = NRep, ncol = length(ENV))
      # rownames(results_eval_metrics_) = paste("REP",1:NRep, sep = "_")
      # colnames(results_eval_metrics_) <-  ENV

      results_eval_metrics_ <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
      rownames(results_eval_metrics_) = paste("REP",1:NRep, sep = "_")
      colnames(results_eval_metrics_) <-  eval_metrics


      # results_eval_metrics_all <- vector(mode = 'list', length = length(eval_metrics))
      # names(results_eval_metrics_all) <- eval_metrics

      results_eval_metrics_all <- vector(mode = 'list', length = length(ENV))
      names(results_eval_metrics_all) <- ENV

      for (c in 1:length(ENV)) {

        results_eval_metrics_all[[c]] <- results_eval_metrics_
      }

      # for (c in 1:length(eval_metrics)) {
      #
      #   results_eval_metrics_all[[c]] <- results_eval_metrics_
      # }

      # results_eval_metrics_MSE <- matrix(NA, nrow = NRep, ncol = length(ENV))
      # rownames(results_eval_metrics_MSE) = paste("REP",1:NRep, sep = "_")
      # colnames(results_eval_metrics_MSE) <-  ENV



      # yHatCV <- matrix(NA, nrow =nfolds, ncol = length(eval_metrics))

      # results_eval_metrics <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
      # rownames(results_eval_metrics) = paste("REP",1:NRep, sep = "_")
      # colnames(results_eval_metrics) = eval_metrics

      rm(results_eval_metrics_)

    }




    current_date_time = as.character(Sys.time())

    files_key = gsub(" ", "", current_date_time)

    ### Create key to remove all reduant files from the wkdir
    #files_key = strsplit(files_key, "\\.")[[1]][2]
    files_key = strsplit(files_key, "\\.")[[1]][1]

    files_key = gsub("-", "", files_key)

    files_key= gsub(":", "_", files_key)
    files_key = paste(files_key, response, sep = "_")

    #DateTime = paste0(gsub(":", "_", DateTime), trait)
    ### Create key to remove all reduant files from the wkdir
    # files_key = strsplit(files_key, "-")[[1]][1]
    #
    # files_key = gsub("-", "", files_key)

    y <-  object[, response]

    for (k in 1:NRep) {

      if(!cross_validation %in%c("Hold_Out",
                                 "Stratified_Hold_Out",
                                 "Repeated_Hold_Out",
                                 "Repeated_Stratified_Hold_Out",
                                 "Leave_one_Out")){

        folds = test_set_val[[k]]

        for(j in 1:max(folds)){
          tst = which(folds ==j)
          yNA = y
          yNA[tst] = NA

    current_date_time = as.character(Sys.time())

    files_key = gsub(" ", "", current_date_time)

    ### Create key to remove all reduant files from the wkdir
    #files_key = strsplit(files_key, "\\.")[[1]][2]
    files_key = strsplit(files_key, "\\.")[[1]][1]

    files_key = gsub("-", "", files_key)

    files_key= gsub(":", "_", files_key)


    if(is.null(weights)){
      fm <- BGLR::BGLR(
        y=yNA,
        ETA = ETA,
        nIter = bayes_para$nIter,
        burnIn =  bayes_para$burnIn,
        thin =  bayes_para$thin,
        verbose = FALSE,
        saveAt =files_key)

    } else{

      if(!is.null(weights)){
        fm <- BGLR::BGLR(
          y=yNA,
          ETA=ETA,
          weights = weights,
          nIter= bayes_para$nIter,
          burnIn= bayes_para$burnIn,
          thin = bayes_para$thin,
          verbose = FALSE,
          saveAt = files_key)

      }

    }

    output_files_names = list.files(pattern=files_key)

    unlink(output_files_names)

    if (cross_validation%in%c("CV1", "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){
      predictions=data.frame(Env=object[tst, heter_groups],
                             Individual=object[tst, gen_name],
                             y=object[tst, response],
                             yHat=fm$yHat[tst])


      predictions=doBy::orderBy(~Env,data=predictions)

      #for (c in 1:length(eval_metrics)) {



        for (v in 1:length(eval_metrics)) {


       sik <- unlist(doBy::lapplyBy(~Env,data=predictions,
                                                     function(x){evaluation_metrics(x$yHat,x$y,
                                                                                    eval_metrics = eval_metrics[v])}))

       for (c in 1:length(ENV)) {

         yHatCV_all[[c]][j,v] <- sik[c]


        }


      }
      # yHatCV_MSE[j, ] <- unlist(doBy::lapplyBy(~Env,data=predictions,
      #                                          function(x){Evaluation_eval_metrics(x$yHat,x$y,
      #                                                                         eval_metrics = "Mean_Squared_Error")}))
      #

    } else {

      if(cross_validation %in%c("K-Folds",
                                "Stratified_K-Folds",
                                "Repeated_K-Folds",
                                "Repeated_Stratified_K-Folds")){

        for (c in 1:ncol(yHatCV)) {

          yHatCV[j, c]  =  evaluation_metrics(y_observed = y[tst],
                                              y_predicted = fm$yHat[tst],
                                              eval_metrics = eval_metrics[c])

        }


      }

    }

        }



      } else {
        tst <-  test_set_val[[k]]
        yNA = y
        yNA[tst] = NA

        if(is.null(weights)){
          fm <- BGLR::BGLR(
            y=yNA,
            ETA = ETA,
            nIter = bayes_para$nIter,
            burnIn =  bayes_para$burnIn,
            thin =  bayes_para$thin,
            verbose = FALSE,
            saveAt =files_key)

        } else{

          if(!is.null(weights)){
            fm <- BGLR::BGLR(
              y=yNA,
              ETA=ETA,
              weights = weights,
              nIter= bayes_para$nIter,
              burnIn= bayes_para$burnIn,
              thin = bayes_para$thin,
              verbose = FALSE,
              saveAt = files_key)

          }


        }


        for (c in 1:ncol(yHatCV)) {

          yHatCV[, c]  =  evaluation_metrics(y_observed = y[tst],
                                              y_predicted = fm$yHat[tst],
                                              eval_metrics = eval_metrics[c])

        }

      }## End

      if (cross_validation%in%c("CV1",
                                "CV2",
                                "Repeated_CV1",
                                "Repeated_CV2")){
        #if (Cross_validation=="CV1" | Cross_validation=="CV2"){
        #results_eval_metrics_Acc[k, ] <- apply(yHatCV_Acc,2, mean)

        #for (c in 1:length(eval_metrics)) {
        for (c in 1:length(ENV)) {

          results_eval_metrics_all[[c]][k,] <- apply(yHatCV_all[[c]],2, mean)



        }



        #results_eval_metrics[[k]] <- list(results_eval_metrics_Acc, results_eval_metrics_MSE)

      } else {

        #ACC = apply(yHatCV,2, mean)
        results_eval_metrics_all[k, ] <- apply(yHatCV,2, mean)

      }



    } ### End Replication

    # ##### Start of mean replication across replications
    if (cross_validation%in%c("CV1",
                              "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){
      #if (Cross_validation=="CV1" | Cross_validation=="CV2"){

      sik = function(x){

        # return(data.frame(metrics= names(x),
        #            value = apply(x, 2, mean)))

        return(apply(x, 2, mean))
      }

      if(NRep>1){

        results_eval_metrics_across_reps = lapply(results_eval_metrics_all, sik)



      } else {

        results_eval_metrics_across_reps = results_eval_metrics_all


      }

    } else {


      if(NRep>1){
        results_eval_metrics_across_reps = data.frame(metrics= colnames(results_eval_metrics_all),
                                                      value = apply(results_eval_metrics_all, 2, mean))

      } else {

        results_eval_metrics_across_reps = data.frame(metrics= colnames(results_eval_metrics_all),
                                                      value = t(results_eval_metrics_all))

      }
    } ### End of mean cal across replication

    ## Get the name of all files stored by GBLR using the current name and time the analysis was performed
    # output_files_names = list.files(pattern=files_key)
    #
    # output = list(model = fm, output_files_names = output_files_names )
    # #rm(output_files_names)
    #
    # names(output) <- c("model", "output_files_names")


    output = list(replication = results_eval_metrics_all, mean_across_reps = results_eval_metrics_across_reps)
    #rm(output_files_names)

    names(output) <- c("replication", "mean")

  #}
  # ### remove the generated output files from the working directory
  # unlink(output_files_names)



  return(output)

}

#################################################################################################
##################################################################################################



#   # #### Initializing parallel
#   if(length(response)>1){
#     if (is.null(core)){
#       cl = parallel::detectCores()
#
#       if (cl> 4){
#         # Try in parallel
#         cl <- parallel::makeCluster(4)
#       } else{
#         cl <- parallel::makeCluster(2)
#       }
#
#     } else {
#       if(!is.null(core)){
#
#         cl <- parallel::makeCluster(core)
#       }
#     }
#
#     doParallel::registerDoParallel(cl)
#
#
#
#   Univariate <- foreach::foreach(trait = 1:length(response),
#                                  .errorhandling='pass') %dopar% {
#
# ###############################################################
# if(!cross_validation %in%c("Hold_Out",
#                            "Stratified_Hold_Out",
#                            "Repeated_Hold_Out",
#                            "Repeated_Stratified_Hold_Out",
#                            "Leave_one_Out")){
#   nfolds <- as.double(strsplit(strsplit(names(test_set_val)[1], c("_"))[[1]][2],"")[[1]][1])  # Number of folds
#
# }
#
# NRep <- length(test_set_val)
#
#
#
#   if(cross_validation %in%c("K-Folds",
#                             "Stratified_K-Folds",
#                             "Repeated_K-Folds",
#                             "Repeated_Stratified_K-Folds")){
#     yHatCV <- matrix(NA, nrow =nfolds, ncol = length(eval_metrics))
#
#     results_eval_metrics_all <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
#     rownames(results_eval_metrics_all) = paste("REP",1:NRep, sep = "_")
#     colnames(results_eval_metrics_all) = eval_metrics
#
#   }
#
#   if(cross_validation %in%c("Hold_Out",
#                             "Stratified_Hold_Out",
#                             "Repeated_Hold_Out",
#                             "Repeated_Stratified_Hold_Out",
#                             "Leave_one_Out")){
#     yHatCV <- matrix(NA, nrow =1, ncol = length(eval_metrics))
#
#     results_eval_metrics_all <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
#     rownames(results_eval_metrics_all) = paste("REP",1:NRep, sep = "_")
#     colnames(results_eval_metrics_all) = eval_metrics
#
#
#   }
#
# ########################################
#
# if (cross_validation %in% c("CV1",
#                             "CV2",
#                             "Repeated_CV1",
#                             "Repeated_CV2")){
#
#   ##### This is to check that condition to execute CV1/CV2 is met
#   # if(!is.null(Geno_data)){
#   #   dimG <- dim(Geno_data)[1]
#   #
#   # } else {
#   #
#   #   if(exists("Geno_data2")){
#   #
#   #     dimG <- dim(Geno_data2)[1]
#   #   }
#   #
#   #   if(!is.null(Gmatrix)){
#   #
#   #     dimG <- dim(Gmatrix)[1]
#   #   }
#   #
#   #   if(exists("GK")){
#   #
#   #     dimG <- dim(GK)[1]
#   #   }
#   #
#   #   if(exists("GK_2")){
#   #
#   #     dimG <- dim(GK_2)[1]
#   #   }
#   # }
#   #
#   # if (nrow(pheno)<=dimG){stop(print(paste(msg,'Number of enviroment/location must be greater than 1 to execute CV1/Repeated_CV1 or CV2/Repeated_CV2')), call. = FALSE)}
#   # #################################
#
#   #if (Cross_validation=="CV1" | Cross_validation=="CV2"){
#
#   ENV = as.character(unique(pheno_data[[heter_groups]]))
#   yHatCV <- matrix(NA, nrow =nfolds, ncol = length(ENV))
#   colnames(yHatCV) = ENV
#   #yHatCV_MSE <- matrix(NA, nrow =nfolds, ncol = length(ENV))
#  yHatCV_all <- vector(mode = 'list', length = length(eval_metrics))
#  names(yHatCV_all) <- eval_metrics
#
#  for (c in 1:length(eval_metrics)) {
#
#    yHatCV_all[[c]] <- yHatCV
#  }
#
#   #colnames(yHatCV_MSE) = ENV
#
#   # results_eval_metrics_Acc <- matrix(NA, nrow = NRep, ncol = length(ENV))
#   # rownames(results_eval_metrics_Acc) = paste("REP",1:NRep, sep = "_")
#   # colnames(results_eval_metrics_Acc) <-  ENV
#
#  results_eval_metrics_ <- matrix(NA, nrow = NRep, ncol = length(ENV))
#  rownames(results_eval_metrics_) = paste("REP",1:NRep, sep = "_")
#  colnames(results_eval_metrics_) <-  ENV
#
#  results_eval_metrics_all <- vector(mode = 'list', length = length(eval_metrics))
#  names(results_eval_metrics_all) <- eval_metrics
#
#  for (c in 1:length(eval_metrics)) {
#
#    results_eval_metrics_all[[c]] <- results_eval_metrics_
#  }
#
#  rm(results_eval_metrics_)
#   # results_eval_metrics_MSE <- matrix(NA, nrow = NRep, ncol = length(ENV))
#   # rownames(results_eval_metrics_MSE) = paste("REP",1:NRep, sep = "_")
#   # colnames(results_eval_metrics_MSE) <-  ENV
#
#
#
#  # yHatCV <- matrix(NA, nrow =nfolds, ncol = length(eval_metrics))
#
#   # results_eval_metrics <- matrix(NA, nrow = NRep, ncol = length(eval_metrics))
#   # rownames(results_eval_metrics) = paste("REP",1:NRep, sep = "_")
#   # colnames(results_eval_metrics) = eval_metrics
#
#
# }
#
#
#
#
#   current_date_time = as.character(Sys.time())
#
#   files_key = gsub(" ", "", current_date_time)
#
#   ### Create key to remove all reduant files from the wkdir
#   #files_key = strsplit(files_key, "\\.")[[1]][2]
#   files_key = strsplit(files_key, "\\.")[[1]][1]
#
#   files_key = gsub("-", "", files_key)
#
#   files_key= gsub(":", "_", files_key)
#   files_key = paste(files_key, response[trait], sep = "_")
#
#   #DateTime = paste0(gsub(":", "_", DateTime), trait)
#   ### Create key to remove all reduant files from the wkdir
#   # files_key = strsplit(files_key, "-")[[1]][1]
#   #
#   # files_key = gsub("-", "", files_key)
#
#   y <-  object[, response[trait]]
#
#   for (k in 1:NRep) {
#
#     if(!cross_validation %in%c("Hold_Out",
#                                "Stratified_Hold_Out",
#                                "Repeated_Hold_Out",
#                                "Repeated_Stratified_Hold_Out",
#                                "Leave_one_Out")){
#
#       folds = test_set_val[[k]]
#
#       for(j in 1:max(folds)){
#         tst = which(folds ==j)
#         yNA = y
#         yNA[tst] = NA
#
#         current_date_time = as.character(Sys.time())
#
#         files_key = gsub(" ", "", current_date_time)
#
#         ### Create key to remove all reduant files from the wkdir
#         #files_key = strsplit(files_key, "\\.")[[1]][2]
#         files_key = strsplit(files_key, "\\.")[[1]][1]
#
#         files_key = gsub("-", "", files_key)
#
#         files_key= gsub(":", "_", files_key)
#
#
#         if(is.null(weights)){
#           fm <- BGLR::BGLR(
#             y=yNA,
#             ETA = ETA,
#             nIter = bayes_para$nIter,
#             burnIn =  bayes_para$burnIn,
#             thin =  bayes_para$thin,
#             verbose = FALSE,
#             saveAt =files_key)
#
#         } else{
#
#           if(!is.null(weights)){
#             fm <- BGLR::BGLR(
#               y=yNA,
#               ETA=ETA,
#               weights = weights,
#               nIter= bayes_para$nIter,
#               burnIn= bayes_para$burnIn,
#               thin = bayes_para$thin,
#               verbose = FALSE,
#               saveAt = files_key)
#
#           }
#
#         }
#
#         output_files_names = list.files(pattern=files_key)
#
#         unlink(output_files_names)
#
#         if (cross_validation%in%c("CV1", "CV2",
#                                   "Repeated_CV1",
#                                   "Repeated_CV2")){
#           predictions=data.frame(Env=object[tst, heter_groups],
#                                  Individual=object[tst, gen_name],
#                                  y=object[tst, response],
#                                  yHat=fm$yHat[tst])
#
#
#           predictions=doBy::orderBy(~Env,data=predictions)
#
#           #for (c in 1:length(eval_metrics)) {
#
#
#
#           for (v in 1:length(eval_metrics)) {
#
#
#             sik <- unlist(doBy::lapplyBy(~Env,data=predictions,
#                                          function(x){evaluation_metrics(x$yHat,x$y,
#                                                                         eval_metrics = eval_metrics[v])}))
#
#             for (c in 1:length(ENV)) {
#
#               yHatCV_all[[c]][j,v] <- sik[c]
#
#
#             }
#
#
#           }
#           # yHatCV_MSE[j, ] <- unlist(doBy::lapplyBy(~Env,data=predictions,
#           #                                          function(x){Evaluation_eval_metrics(x$yHat,x$y,
#           #                                                                         eval_metrics = "Mean_Squared_Error")}))
#           #
#
#         } else {
#
#           if(cross_validation %in%c("K-Folds",
#                                     "Stratified_K-Folds",
#                                     "Repeated_K-Folds",
#                                     "Repeated_Stratified_K-Folds")){
#
#             for (c in 1:ncol(yHatCV)) {
#
#               yHatCV[j, c]  =  evaluation_metrics(y_observed = y[tst],
#                                                   y_predicted = fm$yHat[tst],
#                                                   eval_metrics = eval_metrics[c])
#
#             }
#
#
#           }
#
#         }
#
#       }
#
#
#
#     } else {
#       tst <-  test_set_val[[k]]
#       yNA = y
#       yNA[tst] = NA
#
#       if(is.null(weights)){
#         fm <- BGLR::BGLR(
#           y=yNA,
#           ETA = ETA,
#           nIter = bayes_para$nIter,
#           burnIn =  bayes_para$burnIn,
#           thin =  bayes_para$thin,
#           verbose = FALSE,
#           saveAt =files_key)
#
#       } else{
#
#         if(!is.null(weights)){
#           fm <- BGLR::BGLR(
#             y=yNA,
#             ETA=ETA,
#             weights = weights,
#             nIter= bayes_para$nIter,
#             burnIn= bayes_para$burnIn,
#             thin = bayes_para$thin,
#             verbose = FALSE,
#             saveAt = files_key)
#
#         }
#
#
#       }
#
#
#       for (c in 1:ncol(yHatCV)) {
#
#         yHatCV[, c]  =  evaluation_metrics(y_observed = y[tst],
#                                            y_predicted = fm$yHat[tst],
#                                            eval_metrics = eval_metrics[c])
#
#       }
#
#     }## End
#
#     if (cross_validation%in%c("CV1",
#                               "CV2",
#                               "Repeated_CV1",
#                               "Repeated_CV2")){
#       #if (Cross_validation=="CV1" | Cross_validation=="CV2"){
#       #results_eval_metrics_Acc[k, ] <- apply(yHatCV_Acc,2, mean)
#
#       #for (c in 1:length(eval_metrics)) {
#       for (c in 1:length(ENV)) {
#
#         results_eval_metrics_all[[c]][k,] <- apply(yHatCV_all[[c]],2, mean)
#
#
#
#       }
#
#
#
#       #results_eval_metrics[[k]] <- list(results_eval_metrics_Acc, results_eval_metrics_MSE)
#
#     } else {
#
#       #ACC = apply(yHatCV,2, mean)
#       results_eval_metrics_all[k, ] <- apply(yHatCV,2, mean)
#
#     }
#
#
#
#   } ### End Replication
#
#   # ##### Start of mean replication across replications
#   if (cross_validation%in%c("CV1",
#                             "CV2",
#                             "Repeated_CV1",
#                             "Repeated_CV2")){
#     #if (Cross_validation=="CV1" | Cross_validation=="CV2"){
#
#     sik = function(x){
#
#       # return(data.frame(metrics= names(x),
#       #            value = apply(x, 2, mean)))
#
#       return(apply(x, 2, mean))
#     }
#
#     if(NRep>1){
#
#       results_eval_metrics_across_reps = lapply(results_eval_metrics_all, sik)
#
#
#
#     } else {
#
#       results_eval_metrics_across_reps = results_eval_metrics_all
#
#
#     }
#
#   } else {
#
#
#     if(NRep>1){
#       results_eval_metrics_across_reps = data.frame(metrics= colnames(results_eval_metrics_all),
#                                                     value = apply(results_eval_metrics_all, 2, mean))
#
#     } else {
#
#       results_eval_metrics_across_reps = data.frame(metrics= colnames(results_eval_metrics_all),
#                                                     value = t(results_eval_metrics_all))
#
#     }
#   } ### End of mean cal across replication
#
#   ## Get the name of all files stored by GBLR using the current name and time the analysis was performed
#   # output_files_names = list.files(pattern=files_key)
#   #
#   # output = list(model = fm, output_files_names = output_files_names )
#   # #rm(output_files_names)
#   #
#   # names(output) <- c("model", "output_files_names")
#
#
#   output = list(replication = results_eval_metrics_all, mean_across_reps = results_eval_metrics_across_reps)
#   #rm(output_files_names)
#
#   names(output) <- c("replication", "mean")
#
#
#   Univariate = output
#
#
#         }
#
#   names(Univariate) <- response
#
#
#   }   else {

