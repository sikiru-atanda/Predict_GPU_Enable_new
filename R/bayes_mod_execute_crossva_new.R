

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
#'
# Currently, this function is running multiple response in parallel using future.
# however, I also want to run the number of replication in parallel if greater than 1.
# To optimize resources, this is the way I think about it. Run replication in parallel when response is 1,
# run response in parallel when replication is 1 and run replication and response in parallel when they are both
# greater than 1. Modify this function to reflect these strategies. Keep in mind future has special method of handling
# output.
# with this approach only response run in parallel while replication is stil in for loop. Though this might result in
# future. Is it possible, and careful also wit the handling of result for the nested future.

bayes_mod_execute_crossvall <- function(pheno_data = NULL,
                                        response = NULL,
                                        gen_name = NULL,
                                        test_size = NULL,
                                        random_state = NULL,
                                        replication = NULL,
                                        weights = NULL,
                                        ETA = NULL,
                                        heter_groups = NULL,
                                        bayes_para = NULL,
                                        verbose = FALSE,
                                        core = NULL,
                                        nfolds = NULL,
                                        cross_validation_meth = NULL,
                                        sampling_method = NULL,
                                        eval_metrics = NULL,
                                       ...
)
{

browser()
  n_trait <- length(response)
  msg <- sprintf("==================================================\n")

  holds_out_methods_avail <- c("Hold_Out",
                               "Stratified_Hold_Out",
                               "Repeated_Hold_Out",
                               "Repeated_Stratified_Hold_Out",
                               "Leave_one_Out")

  Kfolds_methods_avail <- c("K-Folds",
                             "Stratified_K-Folds",
                             "Repeated_K-Folds",
                             "Repeated_Stratified_K-Folds")

  CVs_multi_envs_methods_avail <- c("CV1",
                                    "CV2",
                                    "Repeated_CV1",
                                    "Repeated_CV2")
  ###############################################################
  # if(!cross_validation_meth %in%holds_out_methods_avail){
  #   nfolds <- as.double(strsplit(strsplit(names(test_set_val)[1], c("_"))[[1]][2],"")[[1]][1])  # Number of folds
  #
  # }
  #
  # NRep <- length(test_set_val)
  #
  # len_y <-  nrow(pheno_data)
  ##########

  # Check if response has more than one trait and parallel processing is requested
  if (n_trait > 1 && !is.null(core)) {
    future::plan("multisession", workers = core)  # Set the plan for parallel execution
  } else {
    future::plan("sequential")  # Use sequential processing for a single trait or if core is NULL
  }

  # Initialize a list to store results for each trait
  all_traits_results <- list()

  for (trait_idx in seq_along(response)) {
    trait <- response[trait_idx]

    #######
    if(cross_validation_meth %in% c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail)) {
      if(cross_validation_meth %in% holds_out_methods_avail) {
        test_set_val <- hold_out_stratified_and_un(pheno_data = pheno_data,
                                                   gen_name = gen_name,
                                                   response = trait,
                                                   test_size = test_size,
                                                   random_state = NULL,
                                                   replication = replication,  # Use rep here if your function supports per-replication processing
                                                   sampling_method = sampling_method)
      } else if(cross_validation_meth %in% Kfolds_methods_avail) {
        test_set_val <- kfolds_stratified_un(pheno_data = pheno_data,
                                             gen_name = gen_name,
                                             response = trait,
                                             test_size = test_size,
                                             nfolds = nfolds,
                                             random_state = NULL ,
                                             replication = replication,  # Use rep here if your function supports per-replication processing
                                             sampling_method = sampling_method)
      } else if(cross_validation_meth %in% CVs_multi_envs_methods_avail) {
        CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
        test_set_val <- CV1_CV2_for_multi_environment(pheno_data = pheno_data,
                                                      gen_name = gen_name,
                                                      response = trait,
                                                      test_size = test_size,
                                                      nfolds = nfolds,
                                                      CV = CV,
                                                      heter_groups = heter_groups,
                                                      random_state = NULL,
                                                      replication = replication  # Use rep here if your function supports per-replication processing
                                                      )
      }
    } else {
      stop("Unsupported cross-validation method specified. Choose from: ",
           paste(c(holds_out_methods_avail,
                   Kfolds_methods_avail,
                   CVs_multi_envs_methods_avail), collapse = ", "), call. = FALSE)
    }

    if(!cross_validation_meth %in%holds_out_methods_avail){
      nfolds <- as.double(strsplit(strsplit(names(test_set_val)[1], c("_"))[[1]][2],"")[[1]][1])  # Number of folds

    }

    #NRep <- length(test_set_val)
    NRep <- replication

    len_y <-  nrow(pheno_data)

    options(future.rng.onMisuse = "ignore")
    # Wrap the processing logic for each trait in a future call for potential parallel execution
    all_traits_results[[trait]] <- future::future({

    y <- pheno_data[, trait]

  if(cross_validation_meth %in%c(Kfolds_methods_avail,
                                 holds_out_methods_avail)){

    ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
    colnames(ypred_cv) <- c("y", "yhat")
    ypred_cv[, "y"] <- pheno_data[, trait]
    ##################
    ## Add rep to the column in addition to the metrics
    results_eval_metrics_reps <- matrix(NA, nrow = NRep, ncol = length(eval_metrics)+1)
    results_eval_metrics_reps[, 1] <-  c(1:NRep)
    rownames(results_eval_metrics_reps) <-  paste("REP",1:NRep, sep = "_")
    colnames(results_eval_metrics_reps) <-  c("Rep", eval_metrics)

    results_eval_metrics_reps <-  as.data.frame(results_eval_metrics_reps)

  }

  ########################################

  if (cross_validation_meth %in% CVs_multi_envs_methods_avail){

    if(is.null(heter_groups )){
      stop(message(paste(msg,'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)

    }
    ENV <-  as.character(unique(pheno_data[, heter_groups]))

    ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
    colnames(ypred_cv) <- c("y", "yhat", heter_groups)
    ypred_cv[, "y"] <-  pheno_data[, trait]
    ypred_cv[, heter_groups] <-  as.character(pheno_data[, heter_groups])

    results_eval_metrics_reps <- matrix(NA, nrow = length(ENV), ncol = length(eval_metrics)+2) ### Add rep and Env to the cols
    results_eval_metrics_reps[, 1] <-  rep(1, length(ENV))
    results_eval_metrics_reps[, 2] <-  ENV

    #rownames(results_eval_metrics_all) = paste("REP",1:(NRep*length(ENV)), sep = "_")
    colnames(results_eval_metrics_reps) <-  c("Rep", "Env",
                                            eval_metrics)
    results_eval_metrics_reps_use <-   results_eval_metrics_reps

    ## Convert it to empyt dataframe for final storage
    results_eval_metrics_reps <-  data.frame()

  }

  ################

  #y <-  pheno_data[, response]


  for (k in 1:NRep) {

    group <- test_set_val[[k]]

    if(!cross_validation_meth %in%c(holds_out_methods_avail)){

      #folds = test_set_val[[k]]

      for(j in 1:nfolds){

        yNA <- y
        ## set the g fold to NA as the testing set
        for (g in 1:len_y) {
          if(group[g] == j) { yNA[g]<- NA }
        }

        ## extract the position NA which is the testing set
        tst = which(is.na(yNA))

        # systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
        # systime <- gsub("[-: ]", "_", systime)
        # systime <- paste(systime, trait, sep = "_")


            fm <- BGLR::BGLR(
              y=yNA,
              ETA=ETA,
              weights = weights,
              nIter = bayes_para[["nIter"]],
              burnIn =  bayes_para[["burnIn"]],
              thin =  bayes_para[["thin"]],
              verbose = FALSE
              #saveAt =systime
              )

        ypred_cv[tst, "yhat"] <- fm$yHat[tst]

        #output_files_names <-  list.files(pattern=systime)

        #unlink(output_files_names)

      } ### Ends folds

    } else {
      tst <-  test_set_val[[k]]
      yNA <- y
      yNA[tst] <- NA


          fm <- BGLR::BGLR(
            y=yNA,
            ETA=ETA,
            weights = weights,
            nIter = bayes_para[["nIter"]],
            burnIn =  bayes_para[["burnIn"]],
            thin =  bayes_para[["thin"]],
            verbose = FALSE
            #saveAt =systime
            )


      ypred_cv[tst, "yhat"] <- fm$yHat[tst]

      ## extract created files
      #output_files_names = list.files(pattern=systime)
      ## Remove all unused files generated during the Baysian process from the dir
      #unlink(output_files_names)

    }## End

    if (cross_validation_meth %in%CVs_multi_envs_methods_avail){
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


      if (!cross_validation_meth%in%CVs_multi_envs_methods_avail){
        for (eva in 1:length(eval_metrics)) {

          results_eval_metrics_reps[k, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
                                                                                y_predicted = ypred_cv [tst, "yhat"],
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

  if (cross_validation_meth%in%CVs_multi_envs_methods_avail){
    result_mean <- results_eval_metrics_reps |>
      dplyr::group_by(Env) |>
      dplyr::summarise(across(.cols = dplyr::all_of(names(results_eval_metrics_reps)[3:ncol(results_eval_metrics_reps)]),
                              .fns = ~ mean(.x, na.rm = TRUE),
                              .names = "{.col}_mean"))

    names(result_mean)[1] <-  heter_groups

    #result_mean2 = t(result_mean)

  } else {

    if (!cross_validation_meth%in%CVs_multi_envs_methods_avail){
      result_mean <- results_eval_metrics_reps |>
        dplyr::summarise(across(.cols = names(results_eval_metrics_reps)[2:ncol(results_eval_metrics_reps)],
                                .fns = ~ mean(.x, na.rm = TRUE),
                                .names = "{.col}_mean"))


      result_mean <-  data.frame(metrics = colnames(result_mean), summary = t(result_mean))

      rownames(result_mean) <-  NULL

    }

  }

output <- list(eval_metrics_reps= results_eval_metrics_reps, across_reps_eval_metrics = result_mean)

    })

  }

  # Use future::value to retrieve results from futures
  all_traits_results <- future::value(all_traits_results)

  # After processing all traits, you might want to combine or further process results as needed

  return(all_traits_results)


}


# TT = models_execute_crossval(pheno_data = pheno,
#                                  response = c("Yield", "BLUE","deBLUP"),
#                                  gen_name = "GID",
#                                  test_size = 0.3,
#                                  random_state = 123,
#                                  replication = 3,
#                                  heter_groups = "Env",
#                                  cross_validation_meth = "Stratified_Hold_Out",
#                                  nfolds = 5,
#                                  sampling_method = "stratified",
#                                  bayes_para = res_model_prep[["bayes_para"]],
#                                  ETA = res_model_prep[["bayes_ETA"]][["ETA"]],
#                                  ml_dat_res = ml_dat_res,
#                                  GS_model_cv = GS_model_cv,
#                                  num_cores = 5,
#                                  eval_metrics = c("accuracy", "mean_squared_error", "bias",
#                                                   "root_mean_squared_error")
#                                  )
#
#
#
#
#
# fm <- BGLR::BGLR(
#   y=pheno$Yield,
#   ETA=res_model_prep[["bayes_ETA"]][["ETA"]],
#   weights = weights,
#   nIter = res_model_prep[["bayes_para"]][["nIter"]],
#   burnIn =  res_model_prep[["bayes_para"]][["burnIn"]],
#   thin =  res_model_prep[["bayes_para"]][["thin"]],
#   verbose = FALSE,
#
# )
#
# TT = hold_out_stratified_and_un(pheno_data = pheno, gen_name = "GID",
#             response = "Yield", test_size = 0.3, random_state = 123,
#             replication = 3,
#             sampling_method = "stratified")
#
# class(TT$Rep_test_Sample001)
#
# Y = pheno$Yield
#
# Y = pheno[TT$Rep_test_Sample001, ]
