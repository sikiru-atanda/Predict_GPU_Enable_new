
# Function to extract matching variance-covariance structures in asreml model
#' Title
#'
#' @param formula
#' @param var_cov_structures
#'
#' @return
#' @export
#'
#' @examples
extract_var_cov_structures <- function(formula, var_cov_structures) {
  # Use a regular expression to find matches
  pattern <- paste0("\\b(", paste(var_cov_structures, collapse = "|"), ")\\b")
  matches <- regmatches(formula, gregexpr(pattern, formula))
  unique(unlist(matches))
}


#' Bayesian cross-validation fold predictor
#'
#' Fits the configured BGLR model (univariate `BGLR::BGLR`, or the multitrait
#' kernel-prior path when `met_kernel_cv_meta` is supplied) on the training
#' rows and returns predicted values for the held-out test indices.
#'
#' @param y Numeric response vector for the fold, with the test rows already
#'   masked to `NA` by the CV scheduler.
#' @param ETA BGLR `ETA` list prepared upstream by `model_prep_bayes_cv` /
#'   `bayes_finalize_*` (fixed and random effect specifications, kernels or
#'   marker designs, and the model labels per term).
#' @param weights Optional positive Stage 2 observation precisions, supplied as
#'   a numeric vector, one-column table, or column name in `pheno_data`. For
#'   Gaussian fits these are converted to `sqrt(weights)` before calling BGLR,
#'   because BGLR defines residual variance as inverse squared native weight.
#' @param bayes_para Named list of BGLR sampler controls: `nIter`, `burnIn`,
#'   and `thin`.
#' @param tst Integer vector of row indices to predict (the held-out rows).
#' @param bayes_model Bayesian model identifier (e.g. "RKHS", "BRR").
#' @param bayes_trait Trait name used for BGLR output isolation (saveAt prefix).
#' @param pheno_data Optional phenotype data used to resolve a named Stage 2
#'   weight column and preserve row alignment during cross-validation.
#' @param response Optional response-column name in `pheno_data`.
#' @param response_family Response family for the Bayesian fit (gaussian,
#'   ordinal, binary).
#' @param groups Optional residual-variance grouping factor for Gaussian BGLR
#'   fits (per-row environment label for heterogeneous residuals).
#' @param met_kernel_cv_meta Optional list carrying the per-fold MET kernel-prior
#'   multitrait context. When non-NULL, `bayes_mod_cv` routes the fold through
#'   `bayes_mod_cv_met_kernel_predict()` (using `BGLR::Multitrait` with kernel
#'   prior) instead of the univariate `BGLR::BGLR(groups=...)` path. This
#'   keeps CV consistent with the true-prediction fix for GBLUP_BRR/RKHS MET
#'   heter_resid. Required keys: `pheno_data`, `response`, `gen_name`,
#'   `heter_groups`, `kernels`, `GS_model`, `bayes_para`.
#'
#' @return
#' @export
#'
#' @examples
bayes_mod_cv <- function(y,
                        ETA,
                        weights,
                        bayes_para,
                        tst,
                        bayes_model,
                        bayes_trait,
                        pheno_data = NULL,
                        response = NULL,
                        response_family = "gaussian",
                        groups = NULL,
                        met_kernel_cv_meta = NULL){

  fam <- gp_resolve_response_family(response_family, y = y)
  class_levels <- if (identical(fam, "gaussian")) NULL else gp_response_class_levels(y, fam)
  precision_weights <- gp_resolve_stage2_precision_weights(
    weights = weights,
    pheno_data = pheno_data %||% data.frame(.stage2_response = y),
    response = if (is.null(pheno_data)) ".stage2_response" else response,
    test_mask = is.na(y),
    context = "BGLR cross-validation Stage 2 observation weights"
  )
  if (!identical(fam, "gaussian") && isTRUE(precision_weights$supplied) &&
      any(precision_weights$precision != 1)) {
    stop(
      "BGLR ordinal/binary cross-validation fits do not support observation weights other than 1.",
      call. = FALSE
    )
  }
  if (!identical(fam, "gaussian") && !is.null(groups)) {
    stop("BGLR residual groups are supported for gaussian Bayesian cross-validation fits only.", call. = FALSE)
  }

  # MET heter_resid kernel-Bayesian CV: route through the same multitrait
  # kernel-prior sampler the true-prediction path uses (see commit 35a3bcb).
  # The univariate BGLR::BGLR(groups=...) path below collapses the genetic
  # signal for GBLUP_BRR on MET because BGLR's default sigma^2_beta prior
  # does not calibrate to the eigen-sqrt BRR design the same way the kernel
  # prior calibrates to K; this branch makes CV match true-prediction.
  if (!is.null(met_kernel_cv_meta)) {
    if (isTRUE(precision_weights$supplied)) {
      stop(
        "Weighted MET BGLR with heterogeneous residuals is unavailable because BGLR::Multitrait has no observation-weight argument. Set `bayes_kernel_heter_resid = FALSE` to use the weighted univariate BGLR route, or use weighted GP/ASReml.",
        call. = FALSE
      )
    }
    # Phase 3.18: the prep stores met_kernel_cv_meta$response as the FULL
    # user response vector. For multi-trait CV calls this is length > 1,
    # which trips "recursive indexing failed at level 2" inside
    # bayes_mod_cv_met_kernel_predict at
    # `ph[[met_kernel_cv_meta$response]] <- as.numeric(y)`. Each fold
    # predicts one trait at a time so we override to the current
    # bayes_trait. (Phase 3.19/3.20 note: this prevents the crash but
    # bayes_multitrait_env_heter_fit may still return all-NA predictions
    # for the single-trait reduction -- tracked as a deeper architectural
    # fix.)
    if (length(met_kernel_cv_meta[["response"]]) > 1L) {
      if (!is.null(bayes_trait) && nzchar(as.character(bayes_trait)[1L])) {
        met_kernel_cv_meta[["response"]] <- as.character(bayes_trait)[1L]
      } else {
        warning(sprintf(
          "bayes_mod_cv: met_kernel_cv_meta$response has %d traits but bayes_trait is unset; using first '%s'.",
          length(met_kernel_cv_meta[["response"]]),
          as.character(met_kernel_cv_meta[["response"]])[1L]
        ), call. = FALSE)
        met_kernel_cv_meta[["response"]] <- as.character(met_kernel_cv_meta[["response"]])[1L]
      }
    }
    return(bayes_mod_cv_met_kernel_predict(
      y = y, tst = tst, met_kernel_cv_meta = met_kernel_cv_meta))
  }

  save_prefix <- gp_bglr_save_prefix(
    model_name = bayes_model %||% "BGLR",
    response = bayes_trait %||% "trait"
  )

  # Gaussian CV uses only yHat, and the saved effect draws were written to disk
  # (~26 MB per BayesB fold) and deleted unread. Classification needs the draws
  # and re-enables saveEffects in gp_bayes_bglr_fit_observed_rows().
  if (identical(fam, "gaussian")) {
    ETA <- lapply(ETA, function(term) {
      if (is.list(term) && !is.null(term[["saveEffects"]])) term[["saveEffects"]] <- FALSE
      term
    })
  }

  fit <- gp_bayes_bglr_fit(
    fam = fam,
    y = gp_bayes_prepare_response(y, fam, class_levels = class_levels),
    response_type = gp_bayes_response_type(fam),
    ETA=ETA,
    weights = if (isTRUE(precision_weights$supplied)) sqrt(precision_weights$precision) else NULL,
    groups = groups,
    nIter = bayes_para[["nIter"]],
    burnIn =  bayes_para[["burnIn"]],
    thin =  bayes_para[["thin"]],
    verbose = FALSE,
    saveAt = save_prefix
  )

  ## Get the name of all files stored by BGLR using the current name and time the analysis was performed
  output_files_names <- list.files(
    path = dirname(save_prefix),
    pattern = paste0("^", basename(save_prefix)),
    full.names = TRUE
  )

  unlink(output_files_names)

  if (identical(fam, "gaussian")) {
    return(fit$yHat[tst])
  }

  prob <- as.matrix(fit$probs[tst, , drop = FALSE])
  fit_levels <- as.character(fit$levels %||% colnames(prob) %||% character())
  target_levels <- class_levels %||% fit_levels
  if (length(fit_levels) == ncol(prob) && length(intersect(fit_levels, target_levels))) {
    colnames(prob) <- fit_levels
  } else if (length(target_levels) == ncol(prob)) {
    colnames(prob) <- target_levels
  } else {
    observed_levels <- target_levels[target_levels %in% as.character(stats::na.omit(y))]
    if (length(observed_levels) == ncol(prob)) {
      colnames(prob) <- observed_levels
    }
  }
  if (!is.null(target_levels) && length(target_levels) && length(colnames(prob))) {
    missing_levels <- setdiff(target_levels, colnames(prob))
    if (length(missing_levels)) {
      add <- matrix(0, nrow = nrow(prob), ncol = length(missing_levels))
      colnames(add) <- missing_levels
      prob <- cbind(prob, add)
    }
    prob <- prob[, target_levels, drop = FALSE]
  }

  if (identical(fam, "binary")) {
    pred <- as.numeric(prob[, ncol(prob), drop = TRUE])
    attr(pred, "probabilities") <- prob
    return(pred)
  }

  pred <- gp_pred_to_class_labels(prob, class_levels = colnames(prob))
  attr(pred, "probabilities") <- prob
  pred

}

#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param heter_groups
#' @param asreml_models_prep_cv
#' @param response
#' @param tst
#' @param weights Optional positive Stage 2 observation precisions. ASReml uses
#'   them as a data-column weight with Gaussian dispersion fixed at one.
#' @param workspace Optional ASReml workspace setting.
#' @param pworkspace Optional ASReml prediction workspace setting.
#' @param maxit Maximum number of ASReml iterations.
#'
#' @return
#' @export
#'
#' @examples
asreml_mod_cv <- function(pheno_data,
                          gen_name,
                          heter_groups,
                          #var_cov_str,
                          asreml_models_prep_cv,
                          response,
                          tst,
                          weights = NULL,
                          workspace = NULL,
                          pworkspace = NULL,
                          maxit = 50){

  #browser()
asreml_tst_model_cv <- asreml_cv_model(pheno_dataa = pheno_data,
                                       response = response,
                                       gen_name = gen_name,
                                       heter_groups = heter_groups,
                                       asreml_models_prep_cv = asreml_models_prep_cv,
                                       tst = tst,
                                       weights = weights,
                                       workspace = workspace,
                                       pworkspace = pworkspace,
                                       maxit = maxit)

modm <-  asreml_tst_model_cv[["model_cv"]]

### Asreml required the pheno data in the environment to execute the predict func
pheno_dataa <-  asreml_tst_model_cv[["pheno_dataa"]]

if(is.null(modm)){
  return(NULL)
}
GIDs <- as.character(pheno_data[[gen_name]])

GID_tst <- GIDs[tst]

var_cov_str_available <- c("us","corgh",
                           "corh","corv","fa","rr")
code_asr_fit <- asreml_models_prep_cv$code_asr_fit

# Find the index where 'random' appears in the vector
random_index <- grep("random", code_asr_fit, ignore.case = TRUE)

# Check if 'random' term exists and extract it
if (length(random_index) > 0) {
  # Extract the random part of the model
  random_formula <- asreml_models_prep_cv$code_asr_fit[random_index]
  # Apply function and extract matches for the variance-covariance structure
  var_cov_str <- extract_var_cov_structures(random_formula, var_cov_str_available)

  ## This is when var_cov_str return 0 rather than NULL when
  if(length(var_cov_str)==0) var_cov_str <- NULL

} else {
  var_cov_str <- NULL
}

gen_pos <-  asreml_models_prep_cv[["gen_pos"]]
inter_gen_pos <-  asreml_models_prep_cv[["inter_gen_pos"]]
names_in_inv_list <-  asreml_models_prep_cv[["names_in_inv_list"]]
rand_term <-  asreml_models_prep_cv[["rand_term"]]
####
# Create a reference dataframe from pheno_dataa
if(!is.null(heter_groups)){
reference_order <- pheno_dataa |>
  dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups))

# predicted_value <- reference_order |>
#   dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
#   dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)


} else{
  reference_order <- pheno_dataa |>
    dplyr::select(!!dplyr::sym(gen_name))

  # predicted_value <- reference_order |>
  #   dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
  #   dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)

}

cv_heter_groups <- if (!is.null(heter_groups) && !is.null(inter_gen_pos)) heter_groups else NULL
cv_classify <- if (!is.null(cv_heter_groups)) rand_term[[inter_gen_pos]] else gen_name
cv_prediction_bundle <- asreml_predict_or_extract(
  mod = modm,
  classify = cv_classify,
  pheno_data = pheno_dataa,
  response = response,
  gen_name = gen_name,
  heter_groups = cv_heter_groups,
  names_in_inv_list = names_in_inv_list,
  workspace = workspace,
  pworkspace = pworkspace
)
cv_prediction_table <- if (!is.null(cv_heter_groups)) {
  cv_prediction_bundle[["across_env_prediction"]]
} else {
  cv_prediction_bundle[["prediction"]]
}

if (!is.null(cv_prediction_table)) {
  if (!is.null(cv_heter_groups)) {
    predicted_value <- reference_order |>
      dplyr::left_join(cv_prediction_table, by = stats::setNames(c(gen_name, cv_heter_groups), c(gen_name, cv_heter_groups))) |>
      dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(cv_heter_groups), Predicted_value)
    return(as.double(predicted_value[tst, "Predicted_value"]))
  }

  predicted_value <- reference_order |>
    dplyr::left_join(cv_prediction_table, by = stats::setNames(c(gen_name), c(gen_name))) |>
    dplyr::select(!!dplyr::sym(gen_name), Predicted_value)
  return(as.double(predicted_value[tst, "Predicted_value"]))
}

if (!isTRUE(modm$converge)) {
  stop(
    "ASReml CV fit did not converge; prediction and coefficient fallbacks were not used. ",
    "Increase maxit or simplify the fitted covariance model before interpreting CV predictions.",
    call. = FALSE
  )
}

if (!is.null(cv_prediction_bundle$error) && nzchar(cv_prediction_bundle$error)) {
  cat(paste("Error in ASReml CV prediction fallback:", cv_prediction_bundle$error, "\n"))
}

result_met <- tryCatch(
  {
    if (is.null(heter_groups)) {
      predicted_value <- asreml_run_with_model_data_context(
        mod = modm,
        pheno_data = pheno_dataa,
        expr = asreml_predict_pvals(
          mod = modm,
          classify = gen_name,
          workspace = workspace,
          pworkspace = pworkspace
        )
      )
      predicted_value <- reference_order |>
        dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
        dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)

    } else {
      predicted_value <- asreml_run_with_model_data_context(
        mod = modm,
        pheno_data = pheno_dataa,
        expr = asreml_predict_pvals(
          mod = modm,
          classify = rand_term[[inter_gen_pos]],
          workspace = workspace,
          pworkspace = pworkspace
        )
      )
      predicted_value <- reference_order |>
        dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
        dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)

    }
    predicted_value
  },
  error = function(e) {
    # Handle the error, you can print a message or take other actions
    cat(paste("Error in prediction:", conditionMessage(e), "\n"))
    return(NULL)  # Return NULL or an appropriate value to indicate the failure
  }
)


if(is.null(result_met)) {
BLUP <- summary(asreml_tst_model_cv[["model_cv"]], coef=TRUE)$coef.random

colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

#heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

if (!is.null(heter_groups)){
  ### It possible the user provide the heter_groups while it actually a single environment,
  ## This will check and turn it off
  if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
    heter_groups <-  NULL
    inter_gen_pos <-  NULL
  } else{
    if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
      all_envs_for_met <- as.character(pheno_data[[heter_groups]])
    }
  }

}

##############################################
#### Extract Breeding values/genetic effect estimate for all omics
##############################################
estimated_breeding_value_list <- list()

##
for (b in seq_along(names_in_inv_list)) {
  estimated_breeding_value_list[[names_in_inv_list[b]]] <- BLUP[grep(paste(names_in_inv_list[b],"\\)", sep = ""),rownames(BLUP)),]

}
###
### For variance structure extraction
if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

  if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
    ## Extract the number of factors
    #N_fa = substr(var_cov_str, 3, 100)

    for (bb in seq_along(names_in_inv_list)) {
      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][!rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])%in%rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]][grep('Comp',rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])),]), ]

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

      estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,3])

    }

  } else {

    if (var_cov_str %in% c("us", "corgh", "corh", "corv")) {

      for (bb in seq_along(names_in_inv_list)) {

        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

        estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,2])

      }


    }


  } ## End

  #################################
  if(!is.null(inter_gen_pos)){
    for (bb in seq_along(names_in_inv_list)) {

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(4, 1:2)]

      # Environment parsed from the coefficient name (ASReml orders levels by
      # factor level, not by data order); genotype main-effect rows get NA.
      env_parsed <- asreml_coef_row_env(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), heter_groups)
      estimated_breeding_value_list[[names_in_inv_list[bb]]][, heter_groups] <- ifelse(env_parsed %in% heter_grp, env_parsed, NA_character_)

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(1, 4, 2:3)]

      colnames(estimated_breeding_value_list[[names_in_inv_list[bb]]])[1:3] <- c(gen_name, heter_groups, "BLUP")

      #estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Standard_error"]^2
      rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]) <-  NULL


    }

  }

  # Bind all data frames into a single data frame
  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Sum kernels and add each genotype's main effect to every environment.
  summarized_blup <- asreml_legacy_met_genetic_values(combined_df, gen_name, heter_groups, heter_grp)
  names(summarized_blup)[names(summarized_blup) == "BLUP"] <- "Summed_BLUP"

  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
    dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), Summed_BLUP)

  # Lift BLUPs back to the phenotype scale by adding intercept + per-env
  # fixed offsets. Without this, the legacy fallback returns values centered
  # around zero (deviation scale), not the original phenotype scale.
  fixed_offsets <- asreml_extract_fixed_offsets(
    model = modm, heter_groups = heter_groups,
    env_levels = unique(as.character(summarized_blup[[heter_groups]])),
    pheno_data = pheno_data, response = response
  )
  summarized_blup$Summed_BLUP <- summarized_blup$Summed_BLUP +
    as.numeric(fixed_offsets[as.character(summarized_blup[[heter_groups]])])
  return(as.double(summarized_blup[tst, c("Summed_BLUP")]))


}

### For compound symmetric
if(is.null(var_cov_str) & !is.null(inter_gen_pos)){
  for (bb in seq_along(names_in_inv_list)) {
    estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
    cs_rn <- rownames(estimated_breeding_value_list[[bb]])
    estimated_breeding_value_list[[bb]][, gen_name] <- sub("^.*\\)_", "", sub(":.*$", "", sub("^[^:]*:vm\\(", "vm(", cs_rn)))
    attr(estimated_breeding_value_list[[bb]], "coef_env") <- asreml_coef_row_env(cs_rn, heter_groups)
    rownames(estimated_breeding_value_list[[bb]]) <-  NULL
  }
  #################################
  if(!is.null(inter_gen_pos)){
    for (bb in seq_along(names_in_inv_list)) {
      env_parsed <- attr(estimated_breeding_value_list[[bb]], "coef_env")
      estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(4, 1:2)]
      estimated_breeding_value_list[[bb]][, heter_groups] <- ifelse(env_parsed %in% heter_grp, env_parsed, NA_character_)
      estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(1, 4, 2:3)]
      colnames(estimated_breeding_value_list[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
      #estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[bb]][, "Standard_error"]^2
      rownames(estimated_breeding_value_list[[bb]]) <- NULL
    }
  }

  # Bind all data frames into a single data frame
  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  if(!is.null(heter_groups)){
  #combined_df <- combined_df[order(combined_df[[heter_groups]]), ]
  #colnames(combined_df)[colnames(combined_df)%in%heter_groups] <- "Env"

  # Genotype-environment cells: kernels summed, main effect added to every env.
  summarized_blup <- asreml_legacy_met_genetic_values(combined_df, gen_name, heter_groups, heter_grp)
  names(summarized_blup)[names(summarized_blup) == "BLUP"] <- "Summed_BLUP"
  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
    dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), Summed_BLUP)
  # Lift BLUPs back to the phenotype scale (CS MET fallback).
  fixed_offsets <- asreml_extract_fixed_offsets(
    model = modm, heter_groups = heter_groups,
    env_levels = unique(as.character(summarized_blup[[heter_groups]])),
    pheno_data = pheno_data, response = response
  )
  summarized_blup$Summed_BLUP <- summarized_blup$Summed_BLUP +
    as.numeric(fixed_offsets[as.character(summarized_blup[[heter_groups]])])
  return(as.double(summarized_blup[tst, c("Summed_BLUP")]))
  }


  # Group by the user-specified environment/location and sum the BLUP values
  # summarized_blup <- combined_df |>
  #   dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups)) |>
  #   dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))
  #
  # summarized_blup <- as.data.frame(summarized_blup)
  # return(summarized_blup[tst, "Summed_BLUP"])

}

if(is.null(var_cov_str) & is.null(inter_gen_pos) ){
  for (bb in seq_along(names_in_inv_list)){
    estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
    estimated_breeding_value_list[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[bb]]), "\\)_", 3)[,2])
    rownames(estimated_breeding_value_list[[bb]]) <-  NULL

    # if(is.null(estimated_breeding_value_list[[bb]])){
    #   stop("it null")
    # }
    #colnames( estimated_breeding_value_list[[bb]])
    #message(gen_name)
    #message(colnames(estimated_breeding_value_list[[bb]]))
    estimated_breeding_value_list[[bb]] <-   estimated_breeding_value_list[[bb]][, c(4, 1:2)]
    colnames(estimated_breeding_value_list[[bb]])[1:2] <- c(gen_name, "BLUP")
  }

  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Group by GID and sum the BLUP values for each GID
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE), .groups = "drop")

  summarized_blup <- as.data.frame(summarized_blup)

  #summarized_blup <- summarized_blup[as.character(summarized_blup[[gen_name]])%in%GID_tst, ]
  #summarized_blup <- summarized_blup[order(as.character(summarized_blup[[gen_name]])%in%GID_tst), ]
  #### sik
  #gen_namess <- as.character(summarized_blup[[gen_name]])
  # Get the match positions of gen_names in GID_tst
  #order_indices <- match(gen_namess, GID_tst)
  # Order the data frame based on these match positions
  #summarized_blup <- summarized_blup[order(order_indices), ]


  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name), c(gen_name))) |>
    dplyr::select(!!dplyr::sym(gen_name), Summed_BLUP)

  # Lift BLUPs back to phenotype scale (single-env CS fallback).
  intercept <- asreml_extract_fixed_offsets(
    model = modm, heter_groups = NULL, env_levels = NULL,
    pheno_data = pheno_data, response = response
  )
  summarized_blup$Summed_BLUP <- summarized_blup$Summed_BLUP + as.numeric(intercept)
  return(as.double(summarized_blup[tst, "Summed_BLUP"]))

}

} else {

  predicted_value <- result_met

  if(!is.null(heter_groups)){
    #predicted_value <- predicted_value[order(predicted_value[[heter_groups]]), ]

    # predicted_value <- predicted_value[as.character(predicted_value[[gen_name]])%in%GID_tst, ]
    # #predicted_value <-predicted_value[order(as.character(predicted_value[[gen_name]])%in%GID_tst), ]
    # #### sik
    # gen_namess <- as.character(predicted_value[[gen_name]])
    # # Get the match positions of gen_names in GID_tst
    # order_indices <- match(gen_namess, GID_tst)
    # # Order the data frame based on these match positions
    # predicted_value <- predicted_value[order(order_indices), ]

    predicted_value <- reference_order |>
      dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
      dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)

    #colnames(predicted_value)[colnames(predicted_value)%in%heter_groups] <- "Env"
    return(as.double(predicted_value[tst, c("predicted.value")]))
  }else {
    predicted_value <- reference_order |>
      dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
      dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)


    # predicted_value <- predicted_value[as.character(predicted_value[[gen_name]])%in%GID_tst, ]
    # #predicted_value <-predicted_value[order(as.character(predicted_value[[gen_name]])%in%GID_tst), ]
    # #### sik
    # gen_namess <- as.character(predicted_value[[gen_name]])
    # # Get the match positions of gen_names in GID_tst
    # order_indices <- match(gen_namess, GID_tst)
    # # Order the data frame based on these match positions
    # predicted_value <- predicted_value[order(order_indices), ]
    return(as.double(predicted_value[tst, "predicted.value"]))
  }


}



}

## Fit the response (y) centring/scaling on TRAINING rows only. The resulting
## scaler is then applied to all rows and used to back-transform predictions.
## Fitting on the full y would let the held-out test responses leak into the
## centring/scaling used for the training fit and the reverse transform, which
## optimistically biases cross-validation RMSE-type metrics. `tst` = held-out
## (test) row indices for this fold.
gp_ml_cv_train_y_scaler <- function(y, tst) {
  y_df <- as.data.frame(as.matrix(y))
  train_df <- if (length(tst)) y_df[-tst, , drop = FALSE] else y_df
  caret::preProcess(train_df, method = c("center", "scale"))
}

#' Cross-validation fit for an XGBoost genomic prediction model
#'
#' @param y
#' @param omics
#' @param tst
#' @param eta
#' @param nrounds
#' @param max_depth
#' @param scaling
#' @param centering
#' @param omic_count
#' @param gamma
#' @param subsample
#' @param colsample_bytree
#'
#' @return
#' @export
#'
#' @examples
AI_xgboost_cv <- function(y,
                          omics,
                          tst,
                          response_family = "gaussian",
                          xgb_booster = "gbtree",
                          xgb_rate_drop = 0.1,
                          xgb_skip_drop = 0.5,
                          xgb_objective = "reg:squarederror",
                          xgb_sample_type = "uniform",
                          xgb_normalize_type = "tree",
                          eta = 0.01,
                          nrounds = 100,
                          max_depth = 6,
                          scaling = FALSE,
                          centering = TRUE,
                          omic_count,
                          xgb_gamma = 0.01, ## 4
                          min_child_weight = 1,
                          subsample = 0.7,
                          colsample_bytree = 0.7,
                          xgb_alpha = 0.001, ## gblinear
                          xgb_lambda = 1.0, # gblinear,
                          early_stop_for_iteration_xgb = FALSE, ## use when training set is large
                          xgb_nthread = 1L,
                          random_state = NULL
                          ){

  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  class_levels <- if (is_regression) NULL else gp_response_class_levels(y, fam)

  nthread <- suppressWarnings(as.integer(xgb_nthread %||% 1L))
  if (!is.finite(nthread) || nthread < 1L) {
    nthread <- 1L
  }

  prep <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)
  omics <- prep$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  xgb_params <- list(
    eta = eta,
    nrounds = nrounds,
    max_depth = max_depth,
    subsample = subsample,
    xgb_gamma = xgb_gamma,
    colsample_bytree = colsample_bytree,
    min_child_weight = min_child_weight,
    xgb_alpha = xgb_alpha,
    xgb_lambda = xgb_lambda,
    xgb_booster = xgb_booster,
    xgb_rate_drop = xgb_rate_drop,
    xgb_skip_drop = xgb_skip_drop,
    xgb_sample_type = xgb_sample_type,
    xgb_normalize_type = xgb_normalize_type,
    xgb_nthread = nthread,
    random_state = suppressWarnings(as.integer(random_state %||% gp_ml_random_state()))
  )
  preds <- gp_py_ml_fit_predict(
    model_type = "xgboost",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = xgb_params,
    class_levels = class_levels,
    prefer_gpu = TRUE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds


}
######

AI_catboost_cv <- function(y,
                           omics,
                           tst,
                           response_family = "gaussian",
                           catboost_iterations = 500,
                           catboost_depth = 6,
                           catboost_learning_rate = 0.03,
                           catboost_l2_leaf_reg = 3,
                           catboost_thread_count = 1L,
                           scaling = FALSE,
                           centering = TRUE,
                           omic_count = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  class_levels <- if (is_regression) NULL else gp_response_class_levels(y, fam)

  prep <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)
  omics <- prep$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "catboost",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = list(
      catboost_iterations = catboost_iterations,
      catboost_depth = catboost_depth,
      catboost_learning_rate = catboost_learning_rate,
      catboost_l2_leaf_reg = catboost_l2_leaf_reg,
      catboost_thread_count = as.integer(catboost_thread_count %||% 1L)
    ),
    class_levels = class_levels,
    prefer_gpu = FALSE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds
}
######

AI_lightgbm_cv <- function(y,
                           omics,
                           tst,
                           response_family = "gaussian",
                           lightgbm_nrounds = 100,
                           lightgbm_learning_rate = 0.05,
                           lightgbm_num_leaves = 31,
                           lightgbm_feature_fraction = 1.0,
                           lightgbm_bagging_fraction = 1.0,
                           lightgbm_min_data_in_leaf = 20,
                           lightgbm_lambda_l1 = 0,
                           lightgbm_lambda_l2 = 0,
                           lightgbm_nthread = 1L,
                           scaling = FALSE,
                           centering = TRUE,
                           omic_count = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  class_levels <- if (is_regression) NULL else gp_response_class_levels(y, fam)

  prep <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)
  omics <- prep$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "lightgbm",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = list(
      lightgbm_nrounds = lightgbm_nrounds,
      lightgbm_learning_rate = lightgbm_learning_rate,
      lightgbm_num_leaves = lightgbm_num_leaves,
      lightgbm_feature_fraction = lightgbm_feature_fraction,
      lightgbm_bagging_fraction = lightgbm_bagging_fraction,
      lightgbm_min_data_in_leaf = lightgbm_min_data_in_leaf,
      lightgbm_lambda_l1 = lightgbm_lambda_l1,
      lightgbm_lambda_l2 = lightgbm_lambda_l2,
      lightgbm_nthread = as.integer(lightgbm_nthread %||% 1L)
    ),
    class_levels = class_levels,
    prefer_gpu = FALSE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds
}
######

#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param ncomp
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
AI_pls_cv <- function(y,
                      omics,
                      tst,
                      ncomp = 3,
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count){
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data
  y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
  y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  use_auto <- is.null(ncomp) || !is.numeric(ncomp)
  model_params <- if (use_auto) {
    list(pls_auto_components = TRUE, pls_max_components = ncol(omics))
  } else {
    list(ncomp = as.integer(ncomp[[1]]))
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "pls",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = "gaussian",
    model_params = model_params,
    prefer_gpu = FALSE
  )
  revert_scaling(as.numeric(preds), y_scaler)
}
######
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param ntree
#' @param mtry
#' @param maxnodes
#' @param nodesize
#' @param rf_n_jobs
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
AI_randomforest_cv <- function(y,
                               omics,
                               tst,
                               response_family = "gaussian",
                               ntree = 500,
                               mtry = NULL,
                               maxnodes = NULL,
                               nodesize = NULL,
                               rf_n_jobs = 1L,
                               scaling = FALSE,
                               centering = TRUE,
                               omic_count){
  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "randomforest",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = list(
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    ),
    prefer_gpu = FALSE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds


}
########
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
AI_ridge_regression_cv <- function(y,
                                   omics,
                                   tst,
                                   scaling = FALSE,
                                   centering = TRUE,
                                   omic_count){
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data
  y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
  y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  preds <- gp_py_ml_fit_predict(
    model_type = "ridge",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = "gaussian",
    model_params = list(),
    prefer_gpu = FALSE
  )
  revert_scaling(as.numeric(preds), y_scaler)
}
######
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
AI_lasso_cv <- function(y,
                        omics,
                        tst,
                        scaling = FALSE,
                        centering = TRUE,
                        omic_count){
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data
  y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
  y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  preds <- gp_py_ml_fit_predict(
    model_type = "lasso",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = "gaussian",
    model_params = list(),
    prefer_gpu = FALSE
  )
  revert_scaling(as.numeric(preds), y_scaler)
}
#####
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#' @param k
#'
#' @return
#' @export
#'
#' @examples
AI_knn_cv <- function(y,
                      omics,
                      tst,
                      response_family = "gaussian",
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count = NULL,
                      k = 5){
  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "knn",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = list(k = k),
    prefer_gpu = FALSE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds


}

######################
#####
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#' @param c
#'
#' @return
#' @export
#'
#' @examples
AI_svm_cv <- function(y,
                      omics,
                      tst,
                      response_family = "gaussian",
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count,
                      svm_kernel = "Gaussian", # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
                      sigma_value  = NULL,      # NULL uses dimension-aware gamma = "scale"
                      C_value  = 1,             # Default cost parameter
                      degree_value = 3,        # Default degree for polynomial kernel
                      scale_value  = 1,         # Default scale for polynomial kernel
                      offset_value = 0,
                      gamma_value = NULL,
                      svm_type = "eps-regression") {       # Default offset for polynomial kernel
## Scale is not used because it inherently
  fam <- gp_resolve_response_family(response_family, y = y)
  is_regression <- identical(fam, "gaussian")
  omics <- gp_ml_preprocess_predictors(omics, scaling = scaling, centering = centering)$data

  if (is_regression) {
    y_scaler <- gp_ml_cv_train_y_scaler(y, tst)
    y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]
  } else {
    y_fit <- y
  }

  preds <- gp_py_ml_fit_predict(
    model_type = "svm",
    X_train = omics[-tst, , drop = FALSE],
    y_train = y_fit[-tst],
    X_test = omics[tst, , drop = FALSE],
    response_family = fam,
    model_params = list(
      svm_kernel = svm_kernel,
      C_value = C_value,
      degree_value = degree_value,
      scale_value = scale_value,
      offset_value = offset_value,
      gamma_value = gamma_value
    ),
    prefer_gpu = FALSE
  )
  if (is_regression) {
    preds <- revert_scaling(as.numeric(preds), y_scaler)
  }
  preds
  # if(!is.null(omics)) {
  #   if(isTRUE(scaling) || isFALSE(scaling)){
  #     omics <- scale(omics, center = TRUE, scale = TRUE)
  #   }
  # }
  # if(kernel_type!="vanilladot"){
  # fit <-  kernlab::ksvm(x = omics[-tst, ],
  #                       y = y[-tst],
  #                       kernel = kernel_type,
  #                       scaled = FALSE,
  #                       type = "nu-svr",
  #                       C = C_value,
  #                       kpar = kernel_params)
  #
  # } else {
  #   if(kernel_type=="vanilladot"){
  #     AI_fit = kernlab::ksvm(x = omics[-tst, ],
  #                            y = y[-tst],
  #                            kernel = kernel_type,
  #                            scaled = FALSE,
  #                            type = "nu-svr",
  #                            C = C_value
  #     )
  #   }
  #
  # }
  #
  # preds <- kernlab::predict(fit,
  #                           omics[tst, ])
  #
  # preds <- as.data.frame(preds)
  #
  # return(preds[, 1])


}
