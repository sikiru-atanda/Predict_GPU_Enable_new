

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
#' @importFrom (magrittr,"%>%")
#'
bayes_mod_execute_crossval <- function(pheno_object = NULL,
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



  msg <- sprintf("==================================================\n")
  ###############################################################
  if(!cross_validation %in%c("Hold_Out",
                             "Stratified_Hold_Out",
                             "Repeated_Hold_Out",
                             "Repeated_Stratified_Hold_Out",
                             "Leave_one_Out")){
    nfolds <- as.double(strsplit(strsplit(names(test_set_val)[1], c("_"))[[1]][2],"")[[1]][1])  # Number of folds

  }

  NRep <- length(test_set_val)

  len_y = nrow(pheno_object)
  ##########

  if(cross_validation %in%c("K-Folds",
                            "Stratified_K-Folds",
                            "Repeated_K-Folds",
                            "Repeated_Stratified_K-Folds",
                            "Hold_Out",
                            "Stratified_Hold_Out",
                            "Repeated_Hold_Out",
                            "Repeated_Stratified_Hold_Out",
                            "Leave_one_Out"
                            )){

    ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
    colnames(ypred_cv) <- c("y", "yhat")
    ypred_cv[, 1] = pheno_object[, response]
    ##################
    ## Add rep to the column in addition to the metrics
    results_eval_metrics_reps <- matrix(NA, nrow = NRep, ncol = length(eval_metrics)+1)
    results_eval_metrics_reps[, 1] = c(1:NRep)
    rownames(results_eval_metrics_reps) = paste("REP",1:NRep, sep = "_")
    colnames(results_eval_metrics_reps) = c("Rep", eval_metrics)

    results_eval_metrics_reps = as.data.frame(results_eval_metrics_reps)

  }

  ########################################

  if (cross_validation %in% c("CV1",
                              "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){

    if(is.null(heter_groups )){
      stop(message(paste(msg,'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)

    }
    ENV = as.character(unique(pheno_object[, heter_groups]))

    ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
    colnames(ypred_cv) <- c("y", "yhat", heter_groups)
    ypred_cv[, 1] = pheno_object[, response]

    ypred_cv[, 3] = as.character(pheno_object[, heter_groups])


    # results_eval_metrics_reps <- matrix(NA, nrow = (NRep*length(ENV)), ncol = length(eval_metrics)+2) ### Add rep and Env to the cols
    # results_eval_metrics_reps[, 1] = rep(c(1:(NRep)),each = length(ENV))
    # results_eval_metrics_reps[, 2] = rep(ENV,each = NRep)

    results_eval_metrics_reps <- matrix(NA, nrow = length(ENV), ncol = length(eval_metrics)+2) ### Add rep and Env to the cols
    results_eval_metrics_reps[, 1] = rep(1, length(ENV))
    results_eval_metrics_reps[, 2] = ENV

    #rownames(results_eval_metrics_all) = paste("REP",1:(NRep*length(ENV)), sep = "_")
    colnames(results_eval_metrics_reps) = c("Rep", "Env",
                                           eval_metrics)

    results_eval_metrics_reps_use =  results_eval_metrics_reps

    ## Convert it to empyt dataframe for final storage
    results_eval_metrics_reps = data.frame()

  }

 ################

  y <-  pheno_object[, response]

  for (k in 1:NRep) {

    group <- test_set_val[[k]]

    if(!cross_validation %in%c("Hold_Out",
                               "Stratified_Hold_Out",
                               "Repeated_Hold_Out",
                               "Repeated_Stratified_Hold_Out",
                               "Leave_one_Out")){

      #folds = test_set_val[[k]]

      for(j in 1:nfolds){

        yNA <- y
        ## set the g fold to NA as the testing set
        for (g in 1:len_y) {
          if(group[g] == j) { yNA[g]<- NA }
        }

        ## extract the position NA which is the testing set
        tst = which(is.na(yNA))

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


       ypred_cv[tst, "yhat"] <- fm$yHat[tst]

       output_files_names = list.files(pattern=files_key)

       unlink(output_files_names)



      } ### Ends folds



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


      ypred_cv[tst, "yhat"] <- fm$yHat[tst]

      ## extract created files
      output_files_names = list.files(pattern=files_key)
      ## Remove all unused files generated during the Baysian process from the dir
      unlink(output_files_names)

    }## End

    if (cross_validation %in% c("CV1",
                                "CV2",
                                "Repeated_CV1",
                                "Repeated_CV2")){
      ypred_cv = as.data.frame(ypred_cv)
      ypred_cv[, 'y'] =  as.double(ypred_cv[, 'y'])
      ypred_cv[, 'yhat'] =  as.double(ypred_cv[, 'yhat'])

      results_eval_metrics_reps_use = as.data.frame(results_eval_metrics_reps_use)
      for (eva in 1:length(eval_metrics)) {

        sik <- unlist(doBy::lapplyBy(~Env,data=ypred_cv,
                                     function(x){evaluation_metrics(x$yhat,x$y,
                                                                    eval_metrics = eval_metrics[eva])}))

        for (s in 1:length(sik)) {

          results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups]%in%names(sik)[s], c("Rep", eval_metrics[eva])] = c(k, sik[s])


        }


      }

      results_eval_metrics_reps = rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)


    } else {


    if (!cross_validation%in%c("CV1",
                              "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){
    for (eva in 1:length(eval_metrics)) {

      results_eval_metrics_reps[k, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[, "y"],
                                                                          y_predicted = ypred_cv [, "yhat"],
                                                                          eval_metrics = eval_metrics[eva])

         }

      }

    }

  } ### End Replication

 ### Check if some value are not double or numeric and convert it
  if(!all(sapply(eval_metrics, function(x, results_eval_metrics_reps) is.numeric(results_eval_metrics_reps[,x]),  results_eval_metrics_reps))) {
    results_eval_metrics_reps[, eval_metrics] <-
      lapply(results_eval_metrics_reps[, eval_metrics, drop = FALSE],
             function(x) as.double(as.character(x)))

  }

  if (cross_validation%in%c("CV1",
                             "CV2",
                             "Repeated_CV1",
                             "Repeated_CV2")){
  result_mean = results_eval_metrics_reps %>%
    dplyr::group_by(ENV) %>%
    dplyr::summarise_at(.vars = names(.)[3:ncol(results_eval_metrics_reps)],
                        .funs = c(mean="mean"))

  names(result_mean)[1] <-  heter_groups

  #result_mean2 = t(result_mean)

  } else {

    if (!cross_validation%in%c("CV1",
                              "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){
      result_mean = results_eval_metrics_reps %>%
        dplyr::summarise_at(.vars = names(.)[2:ncol(results_eval_metrics_reps)],
                            .funs = c(mean="mean"))

      result_mean = data.frame(metrics = colnames(result_mean), summary = t(result_mean))

      rownames(result_mean) = NULL

    }

  }

  output = list(results_eval_metrics_reps, result_mean)

  names(output) = c("eval_metrics_reps", "across_reps_eval_metrics")


  return(output)


}



