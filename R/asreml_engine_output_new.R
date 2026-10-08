
#' Extract and Process ASReml Model Output
#'
#' This function extracts and processes the output from ASReml models, including coefficients, predicted values, and variance components.
#'
#' @param mod_asreml A list containing the fitted ASReml model and related information.
#' @param pheno_data A data frame containing the phenotypic data.
#' @param response A string specifying the response variable in the phenotypic data.
#' @param gmatrix A matrix of genetic data.
#' @param omic1_kernel An optional kernel matrix for the first omic data.
#' @param omic2_kernel An optional kernel matrix for the second omic data.
#' @param omic3_kernel An optional kernel matrix for the third omic data.
#' @param kernel_list Optional named list of additional relationship or kernel matrices.
#' @param omics_kernel_label A list of labels for the omic kernel matrices.
#' @param heter_groups A factor or vector specifying the heterogeneity groups.
#' @param gen_name A string specifying the column name for the genotypic data.
#' @param var_cov_str An optional string specifying the variance-covariance structure.
#' @param heter_resid An optional string specifying the heterogeneity of residuals.
#' @param workspace A workspace setting passed to ASReml fitting and prediction.
#' @param pworkspace A numeric value specifying the workspace for ASReml. Default is 1e15.
#' @param maxit An integer specifying the maximum number of iterations for ASReml. Default is 50.
#' @param ... Additional arguments to be passed to underlying functions.
#'
#' @return A list containing the processed ASReml model output, including coefficients, predicted values, and variance components.
#' @examples
#' \dontrun{
#'   # Mock preparation of ASReml model components
#'   mod_asreml <- list(
#'     model = "asreml_model",
#'     str_mod = "model_string",
#'     gen_pos = "genetic_position",
#'     inter_gen_pos = "interaction_genetic_position",
#'     names_in_inv_list = c("GID"),
#'     rand_term = "random_term"
#'   )
#'   pheno_data <- data.frame(
#'     GID = 1:10,
#'     Trait = rnorm(10),
#'     Group = factor(rep(1:2, each = 5))
#'   )
#'   response <- "Trait"
#'   gmatrix <- matrix(rnorm(100), nrow = 10)
#'   result <- asreml_mod_output_new(
#'     mod_asreml = mod_asreml,
#'     pheno_data = pheno_data,
#'     response = response,
#'     gmatrix = gmatrix
#'   )
#'   print(result)
#' }
#' @export
#'

asreml_mod_output_new <- function(
    mod_asreml = NULL,
    pheno_data = NULL,
    response = NULL,
    gmatrix = NULL,
    omic1_kernel=NULL,
    omic2_kernel=NULL,
    omic3_kernel=NULL,
    kernel_list = NULL,
    omics_kernel_label = list(omic1_kernel = NULL,
                             omic2_kernel = NULL,
                             omic3_kernel = NULL),
    heter_groups = NULL,
    gen_name = NULL,
    var_cov_str = NULL,
    heter_resid = NULL,
    workspace = NULL,
    pworkspace= 1e15,
    maxit = 50,
    ...
)
{

  msg <- ""

  ##### Print lable
  if(inherits(omics_kernel_label,'list')){
    if(!all(sapply(omics_kernel_label, function(x){ is.null(x)}))!=FALSE){
      label <-  which(sapply(omics_kernel_label, function(x) !is.null(x)))
      print_lable <-  omics_kernel_label[label]
    } else {
      print_lable <-  NULL
    }

  }else {
    if(inherits(omics_kernel_label,"character")){
      print_lable <-  omics_kernel_label
    }
    if(is.null(omics_kernel_label)){
      print_lable <-  NULL
    }
  }
  ######

  mod <-  mod_asreml[["model"]]
  #mod = mod[[1]]
  ### if every variance component is stable the update will not run
  ## by default in asreml so it safe to keep it
  #mod <-  asreml::update.asreml(mod)
  str_mod <-  mod_asreml[["str_mod"]]
  gen_pos <-  mod_asreml[["gen_pos"]]
  inter_gen_pos <-  mod_asreml[["inter_gen_pos"]]
  names_in_inv_list <-  mod_asreml[["names_in_inv_list"]]
  rand_term <-  mod_asreml[["rand_term"]]
  Res_Va_Ve_H2_COV_COR <- NULL
  met_total_genetic_variance <- NULL
  varcov_status <- NULL   # set when variance components cannot be reported

  tst <- NULL
  tst_NA <- NULL
  yy <- as.double(pheno_data[[response]])
  heter_groups_original <- heter_groups
  heter_grp <- if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    as.character(unique(pheno_data[[heter_groups]]))
  } else {
    NULL
  }


  #sik_yna <- pheno_data

  tst_NA <- which(is.na(pheno_data[, response]))

  yNA <- pheno_data[, response]

  #############################
  ### !is.null(var_cov_str) & is.null(inter_gen_pos) incase user provide var_cov_str
  ## while the data is not MT in nature
  if(!is.null(var_cov_str) & is.null(inter_gen_pos)){
    var_cov_str <-  NULL
    heter_groups <-  NULL
    heter_resid <-  NULL
  }

  #########################
  ### Set up parameter for predict function
  asreml_apply_options(
    workspace = workspace,
    pworkspace = pworkspace,
    maxit = maxit
  )

  BLUP <- summary(mod, coef=TRUE)$coef.random

  colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

  #heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

  if (!is.null(heter_groups)){
    ### It possible the user provide the heter_groups while it actually a single environment,
    ## This will check and turn it off
    if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
      heter_groups <-  NULL
    } else{
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
        all_envs_for_met <- as.character(pheno_data[[heter_groups]])
      }
    }

  }
  #ENV_Ids = as.character(unique(pheno_data[[heter_groups]]))
  ##############################################
  #### Extract Breeding values/genetic effect estimate for all omics
  ##############################################
  estimated_breeding_value_list <- asreml_extract_ebv_tables(
    blup = BLUP,
    names_in_inv_list = names_in_inv_list,
    gen_name = gen_name,
    heter_groups = heter_groups_original,
    heter_grp = heter_grp,
    var_cov_str = var_cov_str,
    inter_gen_pos = inter_gen_pos
  )

  ### For variance structure extraction
  if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

    gp_asreml_assert_trustworthy_fit(
      mod,
      context = paste0("ASReml MET GBLUP (", as.character(var_cov_str)[1L], ")")
    )
    res_herit_varCov <- tryCatch(
      asreml_herit_varCov_new(
        model = mod,
        heter_groups = heter_groups,
        var_cov_str = var_cov_str,
        heter_resid = heter_resid,
        names_in_inv_list = names_in_inv_list,
        inter_gen_pos = inter_gen_pos,
        gen_pos = gen_pos,
        response = response,
        gen_name = gen_name
      ),
      error = function(e) {
        reason <- conditionMessage(e)
        # Unconstrained structures (us, corgh, corh, ...) can converge to
        # genetic correlations that do not form a valid (positive
        # semidefinite) matrix, typically with few lines per environment.
        # The fit and its predictions exist; only the variance-covariance
        # summary is not interpretable. Keep the predictions, say so, and
        # report the variance components as unavailable.
        if (asreml_is_invalid_correlation_error(reason)) {
          varcov_status <<- paste(
            "unavailable: the estimated genetic correlations between", heter_groups,
            "levels do not form a valid (positive semidefinite) matrix under",
            paste0("var_cov_str = '", as.character(var_cov_str)[1L], "'."),
            "Predictions are reported; variance components, correlations and reliabilities are not.",
            "Use a constrained structure (e.g. fa1 or corh) or more lines per", heter_groups, "level."
          )
          warning("ASReml MET GBLUP: ", varcov_status, call. = FALSE)
          return(NULL)
        }
        stop(
          "ASReml MET variance-covariance extraction failed: ", reason,
          " Try a constrained structure (e.g. var_cov_str = 'fa1' or 'corh') or more lines per ",
          heter_groups, " level.",
          call. = FALSE
        )
      }
    )

    #VA = Res$Genetic_Var
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in seq_along(names_in_inv_list)){

      if(!is.null(res_herit_varCov)){

        Res_Va_Ve_H2_COV_COR <-  res_herit_varCov
      if(length(names_in_inv_list)>1){
        VA <-  Res_Va_Ve_H2_COV_COR[["varG_per_omics"]][[bb]]
      } else {
        if(length(names_in_inv_list)==1){
          VA <-  unlist(Res_Va_Ve_H2_COV_COR[["Total_genetic_var"]])
        }

      }
      } else{
        VA <- NULL
      }
      met_total_genetic_variance <- suppressWarnings(as.numeric(unlist(
        Res_Va_Ve_H2_COV_COR[["Total_genetic_var"]]
      )))

      if (asreml_met_reliability_variance_valid(VA, heter_grp)) {
        estimated_breeding_value_list[[bb]] <- asreml_apply_met_reliability(
          ebv_df = estimated_breeding_value_list[[bb]],
          heter_groups = heter_groups,
          heter_grp = heter_grp,
          variance_by_env = VA
        )
      } else {
        estimated_breeding_value_list[[bb]][, "Reliability"] <- NA_real_
        message(paste(
          insight::print_color("WARNINGS\n", "blue"),
          insight::print_color(
            paste(
              msg,
              paste(
                "Reliability cannot be estimated because positive finite genetic variance",
                "was not available for every", heter_groups, "level.\n"
              )
            ),
            "blue"
          )
        ))
      }

    }


  } else {
    ## Problem
    if(is.null(var_cov_str) & is.null(inter_gen_pos) ){
      vc <- asreml_varcomp_table(mod)

      #VAR_check_Pos <- which(vc$bound=='F' | vc$bound=='U' | vc$bound=="?" | vc$bound=="S")

      res_comp_check <-  tryCatch({
        VAR_check_Pos <- which(vc$bound=="?" | vc$bound=="S")
        if(length(VAR_check_Pos)>1) {

          # stop(message(paste( insight::print_color("STOP\n", "red"),
          #                insight::print_color(paste(msg, paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
          #                                                      " is unstable, refix the model")), "red"))), call. = FALSE)
          #
          stop(print(paste(msg, paste(paste("Variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
                                      " is unstable, refix the model.\n"))), call. = FALSE)

        } else {

          if(length(VAR_check_Pos)==1) {

            stop(print(paste(msg, paste(paste("Variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" " ),
                                        " is unstable, refix the model.\n"))), call. = FALSE)

          }

        }
           vc
        },
        error = function(e) {
          # Handle the error, you can print a message or take other actions
          cat(paste("Variance component  for", conditionMessage(e), "\n"))
          return(NULL)  # Return NULL or an appropriate value to indicate the failure
        }
      )

      #### Extration of the genetic variance for all omics to calcuate heritability.
      VarG_All <- vector("list", length = length(names_in_inv_list))
      names(VarG_All) <- names_in_inv_list
      VarG <- vc[grep(paste("^vm\\(", gen_name, sep = ""), rownames(vc)), drop=FALSE, ]

      for (g in seq_along(names_in_inv_list)) {
        component_rows <- asreml_single_kernel_varcomp_rows(
          row_names = rownames(VarG),
          gen_name = gen_name,
          kernel_name = names_in_inv_list[[g]]
        )
        VarG_All[[g]] <- VarG[component_rows, "component"]

      }
#### Residual variance
      VE <- vc[grep("!R", row.names(vc)), "component"]
      varG_matrix <-  unlist(VarG_All)
      H <-  matrix(NA, nrow = 1, ncol = length(VE))

      #####
      for (bb in seq_along(names_in_inv_list)){
        estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <- estimated_breeding_value_list[[bb]][, "Standard_error"]^2
        variance_component <- suppressWarnings(as.numeric(VarG_All[[bb]]))
        if (length(variance_component) != 1L ||
            !is.finite(variance_component)) {
          stop(
            "ASReml did not return one finite variance component for ",
            names_in_inv_list[[bb]], ". Matched value: ",
            paste(capture.output(str(VarG_All[[bb]])), collapse = " "),
            "; available genetic rows: ",
            paste(rownames(VarG), collapse = " | "),
            ".",
            call. = FALSE
          )
        }
        if (variance_component <= 0) {
          # A component fixed on the zero boundary is a valid fitted outcome,
          # but component-specific reliability is undefined because its
          # denominator is zero. Keep that limitation explicit.
          estimated_breeding_value_list[[bb]][, "Reliability"] <- NA_real_
        } else {
          estimated_breeding_value_list[[bb]][, "Reliability"] <- round(
            1 - estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] /
              variance_component,
            6
          )
        }

        if(length(VE)>1){
          for (i in 1:length(VE)) {
            H[, i] <- sum(unlist(VarG_All))/(sum(unlist(VarG_All))+VE[i])
          }
        } else {
          if (length(VE)==1){
            #BV$Reliability <- round(1 - BV$PEV/VA,6)
            H <- sum(unlist(VarG_All))/(sum(unlist(VarG_All))+VE)
          }
        }
      }

      VarG_All <-  as.matrix(unlist(VarG_All))

      #colnames(VarG_All) <- "Variance"
      Total_genetic_var <- sum(VarG_All[, 1])
#### Single location
      Res_Va_Ve_H2_COV_COR <-  list(Heritability = H,
                                    varG_per_omics = VarG_All,
                                    Total_genetic_var = Total_genetic_var,
                                    Residual_Var = VE)
    } ## End


  } #### End var_Covar

  ## when variance_covariance structure is not defined by the user and CS is used
  ### For variance structure extraction
  if(is.null(var_cov_str) & !is.null(inter_gen_pos)){
    Res_Va_Ve_H2_COV_COR <-  asreml_herit_CSM_new(model= mod,
                                                 heter_groups= heter_groups,
                                                 heter_resid= heter_resid,
                                                 names_in_inv_list = names_in_inv_list,
                                                 inter_gen_pos = inter_gen_pos,
                                                 gen_pos = gen_pos)

    varG_matrix <- Res_Va_Ve_H2_COV_COR[["varG_per_omics"]]
    ##
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in seq_along(names_in_inv_list)){
      if(length(names_in_inv_list)>1){
        VA <-  Res_Va_Ve_H2_COV_COR[["varG_per_omics"]][bb, ]
      } else {
        if(length(names_in_inv_list)==1){
          VA <-  Res_Va_Ve_H2_COV_COR[["Total_genetic_var"]][1, ]
        }
      }

      if (asreml_met_reliability_variance_valid(VA, heter_grp)) {
        estimated_breeding_value_list[[bb]] <- asreml_apply_met_reliability(
          ebv_df = estimated_breeding_value_list[[bb]],
          heter_groups = heter_groups,
          heter_grp = heter_grp,
          variance_by_env = VA
        )
      } else {
        estimated_breeding_value_list[[bb]][, "Reliability"] <- NA_real_
        message(paste(
          "Reliability cannot be estimated because positive finite genetic variance",
          "was not available for every", heter_groups, "level."
        ))
      }
    }
  } ## End CS
####################
  ##### For MET analysis
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    ### Extract all environments in MET
    heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
    all_envs_for_met <-  as.character(pheno_data[,heter_groups])

    #### End MET
  } else{ ##  start single enviornment results
    ##################################

    Zg <- NULL
    heter_grp <- NULL
    all_envs_for_met <- NULL
  } ##  End single enviornment results
  #############################
  ### Predicted Values
  tst <- which(is.na(yy))
  fallback_prediction_bundle <- NULL
  result_met <- NULL

  gc() ## to free memory

  result_single_loc <- asreml_predict_or_extract(
    mod = mod,
    classify = gen_name,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    names_in_inv_list = names_in_inv_list,
    workspace = workspace,
    pworkspace = pworkspace
  )

  if (identical(result_single_loc$mode, "predict_error")) {
    cat(paste("Error in prediction for single location GBLUP model. Fallback extraction will be used if available:", result_single_loc$error, "\n"))
  }

  if(!is.null(result_single_loc$prediction)){
    fallback_prediction_bundle <- result_single_loc
    predicted_value <- result_single_loc$prediction
  } else {
    if (!is.null(result_single_loc$error) && nzchar(result_single_loc$error)) {
      cat(paste("Error in extracted ASReml fallback prediction:", result_single_loc$error, "\n"))
    }
    predicted_value <- asreml_prediction_from_ebv_tables(
      ebv_tables = estimated_breeding_value_list,
      gen_name = gen_name,
      pheno_data = pheno_data,
      response = response,
      heter_groups = NULL,
      model = mod,
      model_heter_groups = heter_groups
    )
  }
  ## sik

  if(!is.null(Zg)){
    genotype_means <- pheno_data |>
      dplyr::group_by(!!rlang::sym(gen_name)) |>
      dplyr::summarise(mean_value = mean(!!rlang::sym(response), na.rm = TRUE))
    genotype_means <-  as.data.frame(genotype_means)
    # Reorder genotype_means based on name in predicted_value
    genotype_means <- genotype_means[match(predicted_value[, gen_name], genotype_means[[gen_name]]), ]
    yy <- genotype_means[, "mean_value"]
  }

  if (!is.null(result_single_loc$prediction) || !is.null(fallback_prediction_bundle$prediction)) {
    residual_value <- asreml_build_residuals(
      predicted_df = predicted_value,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = NULL
    )
  } else {
    residual_value <- NULL
  }


  if(!is.null(heter_groups)){
    if(var(predicted_value[, "Predicted_value"])==0){
      ### Check if the the across
      message(paste( insight::print_color("WARNING\n", "blue"),
                     insight::print_color(paste(msg, paste(paste('The average prediction across', heter_groups), paste('is a constant value.\n \t Check the model to change', heter_groups), 'to fixed term.\n ')), "blue")))
    }
  }
  gc()
  ################

  if (is.null(heter_groups) & is.null(var_cov_str)) inter_gen_pos <-  NULL
  if(length(gen_pos) == length(rand_term)) inter_gen_pos <- NULL

  if(!is.null(inter_gen_pos)){

   result_met <-  asreml_predict_or_extract(
     mod = mod,
     classify = rand_term[[inter_gen_pos]],
     pheno_data = pheno_data,
     response = response,
     gen_name = gen_name,
     heter_groups = heter_groups,
     names_in_inv_list = names_in_inv_list,
     workspace = workspace,
     pworkspace = pworkspace
   )

   if (identical(result_met$mode, "predict_error")) {
     cat(paste("Error in multi-environment GS prediction. Fallback extraction will be used if available:", result_met$error, "\n"))
   }

    if(is.null(result_met$across_env_prediction)) {
      if (!is.null(fallback_prediction_bundle$across_env_prediction)) {
        across_env_predicted_value <- fallback_prediction_bundle$across_env_prediction
        residual_value_met <- asreml_build_residuals(
          predicted_df = across_env_predicted_value,
          pheno_data = pheno_data,
          response = response,
          gen_name = gen_name,
          heter_groups = heter_groups
        )
      } else {
      across_env_predicted_value <- asreml_prediction_from_ebv_tables(
        ebv_tables = estimated_breeding_value_list,
        gen_name = gen_name,
        pheno_data = pheno_data,
        response = response,
        heter_groups = heter_groups_original,
        model = mod
      )
      residual_value_met <- asreml_build_residuals(
        predicted_df = across_env_predicted_value,
        pheno_data = pheno_data,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups_original
      )
      }
    } else{
      across_env_predicted_value <- result_met[["across_env_prediction"]]
      across_env_predicted_value <- asreml_add_train_test_labels(
        predicted_df = across_env_predicted_value,
        pheno_data = pheno_data,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      residual_value_met <- asreml_build_residuals(
        predicted_df = across_env_predicted_value,
        pheno_data = pheno_data,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
    }

    # ASReml's predict table initially receives generic uncertainty fields.
    # Replace that phenotypic-scale fallback with the fitted environment-
    # specific genetic reliability used by this mixed model.  This is applied
    # to the public prediction table itself (not only the auxiliary EBV table),
    # so subsequent standardization and CSV output retain the correct
    # numerator, denominator, formula and provenance.
    across_env_predicted_value <- asreml_apply_met_reliability(
      ebv_df = across_env_predicted_value,
      heter_groups = heter_groups,
      heter_grp = heter_grp,
      variance_by_env = met_total_genetic_variance
    )

  }

  active_estimated_breeding_value_list <- estimated_breeding_value_list

  ### Initialize step to calculate coefficient for each omics
### Gather all the variables in the function some might be null
  datasets <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  dataset_names <- asreml_kernel_output_names(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  coefficients_list <-   list()
  m_matrix_model_ready_list <-   list()
  sum_estimated_breeding_value <-  0

  # gp_collect_kernel_inputs() preserves argument order in both fitting and
  # output extraction. Align by that invariant and use explicit input roles for
  # public labels. Substring matching is unsafe because names such as
  # "Genomic_layer_A" contain "omic" and were previously misclassified.
  if (length(dataset_names) != length(datasets) ||
      length(datasets) != length(names_in_inv_list)) {
    stop(paste(msg, "Fitted and output kernel components do not align."),
         call. = FALSE)
  }

  ## Intialize step to calculate EBV across environment for each omics
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    if (!is.null(dataset)) {
      ## Kernel coefficients (K b = u) are an auxiliary output: they never
      ## fail an otherwise valid fit.
      coef_name <- paste("coefficient", dataset_names[i], sep = "_")
      # `[<-` with list(NULL) keeps a skipped term's slot (`[[<-` NULL would
      # drop it): the renaming below indexes one slot per kernel, and an
      # empty list failed the whole output ("'names' attribute [1] must be
      # the same length as the vector [0]" on G2F after a 22-minute fit).
      coefficients_list[coef_name] <- list(tryCatch({
        if (!is.null(Zg)) {
          # MET: expand the kernel over the EBV table's own genotype x
          # environment rows, so EBVs and kernel rows pair by key. Expanding
          # over phenotype records paired them by position: wrong unless the
          # records were a complete grid in EBV order, and a length error for
          # unbalanced MET (e.g. G2F, where hybrids are in some environments).
          ebv_tab <- active_estimated_breeding_value_list[[i]]
          env_col <- intersect(c(heter_groups_original, heter_groups), names(ebv_tab))[1L]
          max_rows <- suppressWarnings(as.integer(Sys.getenv("PREDICTPRO_ASREML_COEF_MAX_ROWS", "6000")))
          if (is.na(env_col) || !gen_name %in% names(ebv_tab)) {
            stop("the MET EBV table has no genotype/environment columns")
          }
          if (nrow(ebv_tab) > max_rows) {
            message(sprintf(
              "Skipping MET kernel coefficients for %s: %d genotype x environment rows exceed %d (PREDICTPRO_ASREML_COEF_MAX_ROWS).",
              dataset_names[i], nrow(ebv_tab), max_rows))
            NULL
          } else {
            gid_rows <- as.character(ebv_tab[[gen_name]])
            Zg_ebv <- stats::model.matrix(~ factor(gid_rows, levels = rownames(dataset)) - 1)
            ZgZg <- Zg_ebv %*% dataset %*% t(Zg_ebv)
            suppressMessages({
              ZgZg <- grm_kernel_precheck(ZgZg)
            })
            cal_coeff_asreml(gmatrix = ZgZg, ebv = ebv_tab[["BLUP"]], heter_groups = heter_groups,
                             heter_grp = as.character(ebv_tab[[env_col]]), gid_name = gid_rows)
          }
        } else {
          cal_coeff_asreml(gmatrix = dataset,
                           ebv = active_estimated_breeding_value_list[[i]][, "BLUP"],
                           heter_groups = heter_groups,
                           heter_grp = all_envs_for_met,
                           gid_name = rownames(dataset))
        }
      }, error = function(e) {
        message("ASReml kernel coefficients for ", dataset_names[i], " were not computed: ", conditionMessage(e))
        NULL
      }))

      sum_estimated_breeding_value <- sum_estimated_breeding_value + active_estimated_breeding_value_list[[i]][["BLUP"]]

      #colnames(active_estimated_breeding_value_list[[i]])[which("BLUP"%in%colnames(active_estimated_breeding_value_list[[i]]))] <- "sik" # Estimated_breeding_value

      if(dataset_names[i]=="gmatrix"){
        m_matrix_model_ready_list[[paste(gsub("gmatrix", "geno", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
      } else{
        m_matrix_model_ready_list[[paste(gsub("_kernel", "", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
      }
    }


  }
#### process the sum_ebv
  combined_dff <- dplyr::bind_rows(active_estimated_breeding_value_list)
  if(length(datasets)>1){
    summarized_blup_use <- combined_dff |>
      dplyr::group_by(!!rlang::sym(gen_name))|>
      dplyr::summarise(
        #BLUP = sum(BLUP, na.rm = TRUE),
        BLUP = mean(BLUP, na.rm = TRUE),
        Standard_error = mean(Standard_error, na.rm = TRUE),
        #Standard_error = sqrt(sum(Standard_error^2, na.rm = TRUE)),
        #Reliability = NA,
        .groups = 'drop'  # Ensure the resulting data frame is not grouped
      )

    # Calculate the Prediction Error Variance
    sum_ebv <- summarized_blup_use |>
      dplyr::mutate(Prediction_error_variance = Standard_error^2,
                    Reliability = NA)

  } else{
    if(length(datasets)==1){
      sum_ebv <- dplyr::bind_rows(active_estimated_breeding_value_list)

    }
  }
  ### process of sum_ebv ends
  variance_components <- if (!is.null(Res_Va_Ve_H2_COV_COR)) {
    asreml_variance_components(res_var_cov_h_ve = Res_Va_Ve_H2_COV_COR)
  } else {
    data.frame(
      Components = numeric(0),
      SE = numeric(0),
      Z_ratio = numeric(0),
      Bound = character(0),
      stringsAsFactors = FALSE
    )
  }

  if(inherits(variance_components, "list")){
    covariance <-  variance_components[["Covariance"]]
    correlation <-  variance_components[["Correlation"]]
    variance_components <-  variance_components[["Variance_components"]]


} else{
  if(inherits(variance_components, "data.frame")){
    covariance <-  NULL
    correlation <-  NULL

  }
}

#######################################################
  if(!"geno_model_ready" %in%names(m_matrix_model_ready_list)){
    if(length(m_matrix_model_ready_list)>1){
      names(m_matrix_model_ready_list) <- paste(paste("Omics", seq_along(m_matrix_model_ready_list), sep = ""), "model_ready", sep = "_")
      names(coefficients_list) <- paste(paste("Omics", seq_along(coefficients_list), sep = ""), "coefficient", sep = "_")
      names(active_estimated_breeding_value_list) <- paste(paste("Omics", seq_along(active_estimated_breeding_value_list), sep = ""), "estimated_breeding_value", sep = "_")
    } else {
      names(m_matrix_model_ready_list) <- paste("Omics", "model_ready", sep = "_")
      names(coefficients_list) <- "Coefficient" # paste("Omics", "coefficient", sep = "_")
      names(active_estimated_breeding_value_list) <- "Estimated_breeding_value" # paste("Omics", "estimated_breeding_value", sep = "_")

    }
  } else if(length(grep("omic", names(m_matrix_model_ready_list)))==0 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
    geno_index <- grep("geno", names(m_matrix_model_ready_list))
    names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
    names(coefficients_list)[geno_index] <- "Coefficient"  #paste("Geno", "coefficient", sep = "_")
    names(active_estimated_breeding_value_list)[geno_index] <- "Estimated_breeding_value"  # paste("Geno", "estimated_breeding_value", sep = "_")

  } else{
    if(length(grep("omic", names(m_matrix_model_ready_list)))>=1 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
      omic_index <- grep("omic", names(m_matrix_model_ready_list))
      geno_index <- grep("geno", names(m_matrix_model_ready_list))
      if(length(omic_index)>1){
        names(m_matrix_model_ready_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "coefficient", sep = "_")
        names(active_estimated_breeding_value_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "estimated_breeding_value", sep = "_")
      } else {
        names(m_matrix_model_ready_list)[omic_index] <- paste("Omics", "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste("Omics", "coefficient", sep = "_")
        names(active_estimated_breeding_value_list)[omic_index] <- paste("Omics", "estimated_breeding_value", sep = "_")

      }
      ####
      names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
      names(coefficients_list)[geno_index] <- paste("Geno", "coefficient", sep = "_")
      names(active_estimated_breeding_value_list)[geno_index] <- paste("Geno", "estimated_breeding_value", sep = "_")

    }

  }
  ######
  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)

    # The genomic layer is the gmatrix kernel (this function has no
    # geno_data; the former check stopped every labelled multi-omics fit).
    if (!is.null(gmatrix)) {
      print_lable <- c("Genomic", print_lable)
    }
    if(length(print_lable)>length(m_matrix_model_ready_list) | length(print_lable)<length(m_matrix_model_ready_list)){
      message("The number of omics labels does not match the number of data layers; default layer names were used.")



    } else {
      if(length(print_lable)==length(m_matrix_model_ready_list)){
        rownames(variance_components)[1:length(print_lable)] <- paste(print_lable, "variance", sep="_")
        ######
        names(coefficients_list) <-  paste(print_lable, "coefficient", sep="_")
        names(active_estimated_breeding_value_list) <- paste(print_lable, "estimated_breeding_value", sep="_")
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(m_matrix_model_ready_list) <- paste(print_lable, "model_ready", sep="_")
        ####
      }
    }

  } else {
    if(is.null(gmatrix)){
    if (length(m_matrix_model_ready_list) > 1L) message("No omics labels were provided; default layer names were used.")
}

  }
#####################################################################
if(exists("sum_ebv")) sum_ebv <- as.data.frame(sum_ebv)
  has_met_prediction <- !is.null(Zg) && exists("across_env_predicted_value") && !is.null(across_env_predicted_value)

  if(is.null(Zg) || !has_met_prediction){

    res <- list(Coefficients = coefficients_list,
                Asreml_model = mod,
                Estimated_breeding_value = active_estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Residual_value = residual_value,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list
    )

  } else {
      if(is.null(correlation) & is.null(covariance)){
    res <- list(Coefficients = coefficients_list,
                Asreml_model = mod,
                Estimated_breeding_value = active_estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Total_Predicted_value =  predicted_value,
                Predicted_value = across_env_predicted_value,
                Residual_value = residual_value,
                Residual_value_MET = residual_value_met,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list)

      } else {
        if(!is.null(correlation) & !is.null(covariance)){
          res <- list(Coefficients = coefficients_list,
                      Asreml_model = mod,
                      Estimated_breeding_value = active_estimated_breeding_value_list,
                      Total_estimated_breeding_value = sum_ebv,
                      Total_Predicted_value =  predicted_value,
                      Predicted_value = across_env_predicted_value,
                      Correlation = correlation,
                      Covariance = covariance,
                      Residual_value = residual_value,
                      Residual_value_MET = residual_value_met,
                      Variance_components = variance_components,
                      M_matrix_model_ready =  m_matrix_model_ready_list)

        }
      }
  }

  if (!is.null(varcov_status)) {
    res$variance_component_status <- varcov_status
  }

  if (!is.null(Res_Va_Ve_H2_COV_COR)) {
    canonical_matrix_names <- c(
      "Genetic_covariance_environments",
      "Genetic_correlation_environments",
      "Residual_covariance_environments",
      "Residual_correlation_environments"
    )
    for (matrix_name in canonical_matrix_names) {
      matrix_value <- Res_Va_Ve_H2_COV_COR[[matrix_name]]
      if (is.matrix(matrix_value)) res[[matrix_name]] <- matrix_value
    }
    for (table_name in c(
      "Residual_variance_by_environment",
      "Residual_variance_by_observation"
    )) {
      table_value <- Res_Va_Ve_H2_COV_COR[[table_name]]
      if (is.data.frame(table_value) && nrow(table_value)) {
        res[[table_name]] <- table_value
      }
    }
    if (!is.null(var_cov_str) && !is.null(inter_gen_pos)) {
      res[["asreml_variance_covariance_structure"]] <- as.character(var_cov_str)[1L]
    }
  }

  res[["model_parameters"]] <- rbind(
    gp_multi_kernel_parameter_rows(
      kernel_names = names(datasets),
      strategy = paste0(
        "ASReml_",
        as.character(var_cov_str %||% "independent")[[1L]],
        "_per_kernel_random_effects"
      )
    ),
    data.frame(
      stat = c(
        "stage2_observation_weighting",
        "stage2_weight_source",
        "asreml_variance_covariance_structure",
        "asreml_residual_structure",
        "asreml_residual_scale_parameter",
        "asreml_residual_scale_status",
        "residual_variance_basis",
        "heritability_residual_basis"
      ),
      summary = c(
        if (!is.null(mod_asreml[["stage2_weight_column"]])) {
          "precision_fixed_residual_variance_dispersion_1"
        } else {
          "not_supplied"
        },
        as.character(mod_asreml[["stage2_weight_source"]] %||% "none"),
        as.character(var_cov_str %||% "independent")[[1L]],
        if (isTRUE(heter_resid)) "heterogeneous_by_environment" else "common",
        if (!is.null(Res_Va_Ve_H2_COV_COR[["Residual_scale_parameter"]])) {
          scale_value <- Res_Va_Ve_H2_COV_COR[["Residual_scale_parameter"]]
          paste0(names(scale_value), "=", format(as.numeric(scale_value), digits = 15), collapse = ";")
        } else {
          "not_available"
        },
        if (!is.null(mod_asreml[["stage2_weight_column"]])) {
          "fixed_by_gaussian_dispersion_1"
        } else {
          "estimated_by_ASReml"
        },
        as.character(
          Res_Va_Ve_H2_COV_COR[["Residual_variance_basis"]] %||%
            "ASReml_fitted_residual_component"
        ),
        if (!is.null(mod_asreml[["stage2_weight_column"]])) {
          "per_environment_mean_observation_residual_variance"
        } else {
          "ASReml_fitted_residual_component"
        }
      ),
      stringsAsFactors = FALSE
    )
  )

    return(res)
}

# TRUE for the variance-covariance failures caused by estimated correlations
# that do not form a valid (positive semidefinite) matrix.
asreml_is_invalid_correlation_error <- function(message) {
  grepl("Reconstructed correlation matrix is invalid|A valid correlation matrix|not positive semidefinite|Correlation parameters must be in \\[-1, 1\\]",
        message, ignore.case = TRUE)
}
