

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

 ##browser()
  n_trait <- length(response)
  msg <- "\n==================================================\n"

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

  # Determine the parallel strategy based on the number of traits and replications
  if (!is.null(core) && core > 1) {
    if (n_trait > 1 && replication <= 1) {
      future::plan("multisession", workers = core)
    } else if (replication > 1 && n_trait <= 1) {
      future::plan("multisession", workers = core)
    } else if (n_trait > 1 && replication > 1) {
      # When both traits and replications are greater than 1, set up a more complex parallel strategy
      # This example uses multisession for traits; further nested parallelism for replications could be complex and requires careful management
      future::plan("multisession", workers = core)
    } else {
      future::plan("sequential")
    }
  } else {
    future::plan("sequential")
  }

  # Create a list of all combinations of response variables and replications
  # tasks <- expand.grid(response = response, replication = seq_len(replication),
  #                      model = GS_model,
  #                      stringsAsFactors = FALSE)

  tasks <- expand.grid(response = response, replication = seq_len(replication),
                       stringsAsFactors = FALSE)

  # Execute each task in parallel
  #process_task <- function(task) {
  results <- future.apply::future_lapply(seq_len(nrow(tasks)), function(i) {
    task_row <- tasks[i, ]

    trait <- as.character(task_row$response)
    rep <- as.integer(task_row$replication)

    #print(c(trait, rep))  # For diagnostic purposes

    repp = 1
    #######
    # Check if a seed was provided and calculate a new seed based on the replication number
    if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1) {
      # Ensure random_state is an integer
      base_seed <- as.integer(random_state)
      # Generate a new, valid integer seed for each replication
      # The modulo operation ensures the seed stays within the integer range
      new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
      #new_seed <- base_seed + rep
    }

    # else {
    #   # Fallback seed if random_state is not set
    #   base_seed <- 123L
    #   new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
    # }

    # Apply the new seed for replication this is important to prevent same seed is used across replication

    #######
          if(cross_validation_meth %in% c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail)) {
            if(cross_validation_meth %in% holds_out_methods_avail) {
              test_set_val <- hold_out_stratified_and_un(pheno_data = pheno_data,
                                                         gen_name = gen_name,
                                                         response = trait,
                                                         test_size = test_size,
                                                         random_state = new_seed,
                                                         replication = repp,  # Use rep here if your function supports per-replication processing
                                                         sampling_method = sampling_method)
            } else if(cross_validation_meth %in% Kfolds_methods_avail) {
              test_set_val <- kfolds_stratified_un(pheno_data = pheno_data,
                                                   gen_name = gen_name,
                                                   response = trait,
                                                   test_size = test_size,
                                                   nfolds = nfolds,
                                                   random_state = new_seed ,
                                                   replication = repp,  # Use rep here if your function supports per-replication processing
                                                   sampling_method = sampling_method)
            } else if(cross_validation_meth %in% CVs_multi_envs_methods_avail) {
              CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
              test_set_val <- CV1_CV2_for_multi_environment(pheno_data = pheno_data,
                                                            gen_name = gen_name,
                                                            response = trait,
                                                            test_size = test_size,
                                                            CV = CV,
                                                            nfolds = nfolds,
                                                            heter_groups = heter_groups,
                                                            random_state = new_seed,
                                                            replication = repp,  # Use rep here if your function supports per-replication processing
                                                            sampling_method = sampling_method)
            }
          } else {
            stop("Unsupported cross-validation method specified. Choose from: ",
                 paste(c(holds_out_methods_avail,
                         Kfolds_methods_avail,
                         CVs_multi_envs_methods_avail), collapse = ", "), call. = FALSE)
          }

    len_y <-  nrow(pheno_data)

      y <- pheno_data[, trait]

      if(cross_validation_meth %in%c(Kfolds_methods_avail,
                                     holds_out_methods_avail)){

        ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
        colnames(ypred_cv) <- c("y", "yhat")
        ypred_cv[, "y"] <- pheno_data[, trait]
        ##################
        ## Add rep to the column in addition to the metrics
        results_eval_metrics_reps <- matrix(NA, nrow = repp, ncol = length(eval_metrics)+1)
        results_eval_metrics_reps[, 1] <-  c(1:repp)
        rownames(results_eval_metrics_reps) <-  paste("REP",1:repp, sep = "_")
        colnames(results_eval_metrics_reps) <-  c("Rep", eval_metrics)

        results_eval_metrics_reps <-  as.data.frame(results_eval_metrics_reps)

      }

      ########################################

      if (cross_validation_meth %in% CVs_multi_envs_methods_avail){

        if(is.null(heter_groups )){
          stop(message(paste(msg,'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)

        }
        ENV <-  as.character(unique(pheno_data[[heter_groups]]))

        ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
        colnames(ypred_cv) <- c("y", "yhat", heter_groups)
        ypred_cv[, "y"] <-  pheno_data[, trait]
        ypred_cv[, heter_groups] <-  as.character(pheno_data[[heter_groups]])

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
      #y <-  pheno_data[[response]]

      #for (k in 1:NRep) {

        group <- test_set_val[[repp]]

        if(!cross_validation_meth %in%holds_out_methods_avail){

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
          tst <-  test_set_val[[repp]]
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
              results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups]%in%names(sik)[s], c("Rep", eval_metrics[eva])] = c(repp, sik[s])

            }
          }
          results_eval_metrics_reps = rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)
        } else {
          if (!cross_validation_meth%in%CVs_multi_envs_methods_avail){
            for (eva in 1:length(eval_metrics)) {
              results_eval_metrics_reps[repp, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
                                                                                    y_predicted = ypred_cv [tst, "yhat"],
                                                                                    eval_metrics = eval_metrics[eva])
            }

          }

        }

      #} ### End Replication

      ### Check if some value are not double or numeric and convert it
      if(!all(sapply(eval_metrics, function(x, results_eval_metrics_reps) is.numeric(results_eval_metrics_reps[,x]),  results_eval_metrics_reps))) {
        results_eval_metrics_reps[, eval_metrics] <-
          lapply(results_eval_metrics_reps[, eval_metrics, drop = FALSE],
                 function(x) as.double(as.character(x)))

      }


      list(trait = trait, rep = rep, eval_metrics_reps = results_eval_metrics_reps)
  }, future.seed = TRUE)



  future::plan("sequential")  # Reset to default plan
  return(results) # Adjust depending on how you want to return or further process the results


}
####
# # Initialize empty lists to hold the data frames
# trait_dfs <- list()
# mean_metrics_dfs <- list()
#
# # Extract unique traits from the results
# traits <- unique(sapply(TT, function(x) x$trait))
#
# # Aggregate results by trait and calculate mean metrics
# for (trait in traits) {
#   # Filter results for the current trait
#   trait_results <- Filter(function(x) x$trait == trait, TT)
#
#   # Combine eval_metrics_reps into a single data frame for the trait
#   trait_df <- do.call(rbind, lapply(trait_results, function(x) x$eval_metrics_reps))
#   trait_dfs[[trait]] <- trait_df
#
#   # Calculate mean of evaluation metrics across all replications for the trait
#   mean_metrics <- colMeans(trait_df[, -1], na.rm = TRUE)  # Assuming first column is 'Rep'
#   mean_metrics_df <- as.data.frame(t(mean_metrics))
#   colnames(mean_metrics_df) <- names(mean_metrics)
#   mean_metrics_dfs[[trait]] <- mean_metrics_df
# }
#
# # Assign names to the lists based on trait names for easier reference
# names(trait_dfs) <- traits
# names(mean_metrics_dfs) <- traits
#
# # Now, trait_dfs contains a data frame for each trait with all replications
# # mean_metrics_dfs contains a data frame for each trait with mean metrics
#
# # To access the data frame for a specific trait, e.g., "Yield", use:
# yield_df <- trait_dfs[["Yield"]]
# yield_mean_metrics <- mean_metrics_dfs[["Yield"]]
#
# # To see the structure or print the results, you can use str() or print():
# str(yield_df)
# print(yield_mean_metrics)
#
#
