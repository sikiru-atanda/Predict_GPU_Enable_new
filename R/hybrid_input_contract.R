#' Describe the hybrid prediction input standard
#'
#' @return A list describing hybrid phenotype, parent-resolved genomic,
#'   hybrid-level genomic, and CV scenario requirements.
#' @export
hybrid_data_standard <- function() {
  list(
    phenotype = list(
      required = c("HybridID", "Female", "Male", "Response"),
      optional = c("Env", "Rep", "Year", "Location", "heter_groups"),
      notes = c(
        "Use one row per hybrid observation.",
        "For multi-environment hybrid prediction, include an environment column and pass it through heter_groups.",
        "Phenotypes can be provided as one response column at a time inside model_execute(response = ...)."
      )
    ),
    parent_resolved_genomics = list(
      accepted = c(
        "female_geno_data + male_geno_data",
        "female_gmatrix + male_gmatrix",
        "shared geno_data keyed by parent IDs",
        "shared gmatrix keyed by parent IDs"
      ),
      used_by = c("hybrid_asreml", "hybrid_bayes", "hybrid_gp"),
      notes = c(
        "Row names must be parent IDs matching the Female and Male phenotype columns.",
        "Parent-resolved formats are required for explicit Female_GCA + Male_GCA + SCA kernels.",
        "hybrid_ml/hybrid_dl optionally accept female_geno_data + male_geno_data for parent heterozygosity QC, or to build expected F1 hybrid features."
      )
    ),
    hybrid_level_genomics = list(
      accepted = c("geno_data keyed by HybridID"),
      used_by = c("hybrid_ml", "hybrid_dl"),
      notes = c(
        "Row names must be hybrid IDs matching the HybridID phenotype column.",
        "Hybrid-level genotype is suitable for predictive hybrid ML/DL, not parent-resolved GCA/SCA kernels."
      )
    ),
    cv_scenarios = c(
      "Hybrid_Known_Parents",
      "Hybrid_One_New_Parent",
      "Hybrid_Both_New_Parents",
      "CV0", "CV1", "CV2",
      "Repeated_CV0", "Repeated_CV1", "Repeated_CV2"
    ),
    current_scope = list(
      hybrid_asreml_models = gp_hybrid_asreml_supported_models(),
      hybrid_bayes_models = gp_hybrid_bayes_supported_models(),
      hybrid_gp_models = gp_display_supported_model_names(gp_hybrid_gp_supported_models()),
      hybrid_ml_models = gp_hybrid_ml_supported_models(),
      hybrid_dl_models = gp_display_supported_model_names(gp_hybrid_dl_supported_models())
    )
  )
}

gp_hybrid_standard_contract_lines <- function(mode = NULL) {
  spec <- hybrid_data_standard()
  mode_label <- if (is.null(mode) || !nzchar(mode)) {
    "hybrid prediction"
  } else {
    paste0(mode, " hybrid prediction")
  }
  c(
    paste("PredictProR detected that the uploaded data is not in the required", mode_label, "format."),
    "",
    "Hybrid phenotype standard:",
    "  required columns: HybridID, Female, Male, and the selected response column",
    "  optional columns: Env, Rep, Year, Location",
    "  for multi-environment hybrid prediction: provide the environment column and pass it as heter_groups",
    "",
    "Parent-resolved genomic standard for hybrid_asreml, hybrid_bayes, and hybrid_gp:",
    "  accepted inputs: female_geno_data + male_geno_data, female_gmatrix + male_gmatrix, or a shared parent-keyed geno_data/gmatrix",
    "  row names must match parent IDs in Female and Male",
    "",
    "Hybrid-level genomic standard for hybrid_ml and hybrid_dl:",
    "  accepted input: geno_data keyed by HybridID",
    "  row names must match HybridID",
    "  optional: female_geno_data + male_geno_data keyed by parent IDs (parent heterozygosity QC; expected F1 features when geno_data is absent)",
    "",
    "Supported hybrid CV scenarios:",
    paste(" ", paste(spec$cv_scenarios, collapse = ", ")),
    "",
    "Supported hybrid models:",
    paste("  hybrid_asreml:", paste(spec$current_scope$hybrid_asreml_models, collapse = ", ")),
    paste("  hybrid_bayes:", paste(spec$current_scope$hybrid_bayes_models, collapse = ", ")),
    paste("  hybrid_gp:", paste(spec$current_scope$hybrid_gp_models, collapse = ", ")),
    paste("  hybrid_ml:", paste(spec$current_scope$hybrid_ml_models, collapse = ", ")),
    paste("  hybrid_dl:", paste(spec$current_scope$hybrid_dl_models, collapse = ", ")),
    "",
    "Please reformat the uploaded data to match one of these standards before rerunning model_execute()."
  )
}

gp_hybrid_stop_contract <- function(details = NULL, mode = NULL) {
  lines <- gp_hybrid_standard_contract_lines(mode = mode)
  if (!is.null(details) && nzchar(details)) {
    lines <- c(details, "", lines)
  }
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

gp_hybrid_assert_single_path <- function(hybrid_asreml = FALSE,
                                         hybrid_bayes = FALSE,
                                         hybrid_gp = FALSE,
                                         hybrid_ml = FALSE,
                                         hybrid_dl = FALSE) {
  active <- c(
    hybrid_asreml = isTRUE(hybrid_asreml),
    hybrid_bayes = isTRUE(hybrid_bayes),
    hybrid_gp = isTRUE(hybrid_gp),
    hybrid_ml = isTRUE(hybrid_ml),
    hybrid_dl = isTRUE(hybrid_dl)
  )
  if (sum(active) > 1L) {
    gp_hybrid_stop_contract(
      details = paste(
        "Select only one hybrid prediction path at a time.",
        "Current flags:",
        paste(names(active)[active], collapse = ", ")
      ),
      mode = "general"
    )
  }
  names(active)[active]
}

gp_hybrid_pick_pheno_input <- function(pheno_data = NULL,
                                       pheno_data_train = NULL,
                                       pheno_data_test = NULL) {
  if (!is.null(pheno_data)) {
    return(as.data.frame(pheno_data, stringsAsFactors = FALSE))
  }
  if (!is.null(pheno_data_train) || !is.null(pheno_data_test)) {
    gp_hybrid_stop_contract(
      details = paste(
        "Hybrid prediction currently expects a single phenotype table in pheno_data.",
        "Do not provide hybrid inputs only through pheno_data_train/pheno_data_test."
      ),
      mode = "general"
    )
  }
  gp_hybrid_stop_contract(
    details = "Provide pheno_data for hybrid prediction.",
    mode = "general"
  )
}

gp_hybrid_validate_phenotype_contract <- function(ph,
                                                  response,
                                                  gen_name,
                                                  female_parent,
                                                  male_parent,
                                                  heter_groups = NULL,
                                                  mode = NULL) {
  if (!is.data.frame(ph)) {
    ph <- as.data.frame(ph, stringsAsFactors = FALSE)
  }
  if (is.null(gen_name) || !nzchar(gen_name) || !gen_name %in% names(ph)) {
    gen_name_label <- if (is.null(gen_name)) "NULL" else as.character(gen_name)
    gp_hybrid_stop_contract(
      details = paste(
        "The hybrid phenotype table must include the hybrid ID column declared in gen_name.",
        "Current gen_name:", gen_name_label
      ),
      mode = mode
    )
  }
  response_cols <- unique(as.character(response %||% character()))
  invalid_response <- !length(response_cols) || any(!nzchar(response_cols)) ||
    any(is.na(response_cols))
  missing_response <- setdiff(response_cols, names(ph))
  if (identical(mode, "hybrid_gp")) {
    if (invalid_response || length(missing_response)) {
      gp_hybrid_stop_contract(
        details = paste(
          "Hybrid GP prediction expects one or more gaussian response columns, and every response must exist in pheno_data.",
          "Current response:", paste(response, collapse = ", ")
        ),
        mode = mode
      )
    }
  } else if (length(response_cols) != 1L || invalid_response || length(missing_response)) {
    gp_hybrid_stop_contract(
      details = paste(
        gp_hybrid_model_scope_label(mode),
        "currently expects one response column at a time, and that response must exist in pheno_data.",
        "Current response:", paste(response, collapse = ", ")
      ),
      mode = mode
    )
  }
  req <- c(female_parent, male_parent)
  if (any(is.null(req)) || any(!nzchar(req)) || !all(req %in% names(ph))) {
    gp_hybrid_stop_contract(
      details = paste(
        "Provide valid female_parent and male_parent column names present in pheno_data.",
        "Missing columns:",
        paste(setdiff(req, names(ph)), collapse = ", ")
      ),
      mode = mode
    )
  }
  if (!is.null(heter_groups) && nzchar(heter_groups) && !heter_groups %in% names(ph)) {
    gp_hybrid_stop_contract(
      details = paste(
        "heter_groups was provided for hybrid prediction, but that column is not present in pheno_data.",
        "Current heter_groups:", heter_groups
      ),
      mode = mode
    )
  }
  ph
}

gp_hybrid_validate_parent_keyed_matrix <- function(x, arg_name, female_ids, male_ids, mode = NULL) {
  if (is.null(x)) {
    return(invisible(FALSE))
  }
  rn <- rownames(x)
  if (is.null(rn) || !length(rn)) {
    gp_hybrid_stop_contract(
      details = paste(arg_name, "must have row names keyed by parent IDs."),
      mode = mode
    )
  }
  parent_ids <- unique(c(as.character(female_ids), as.character(male_ids)))
  overlap <- sum(parent_ids %in% rn)
  if (!overlap) {
    gp_hybrid_stop_contract(
      details = paste(
        arg_name,
        "does not overlap the Female/Male parent IDs in the phenotype table.",
        "Parent-resolved hybrid prediction requires parent-keyed row names."
      ),
      mode = mode
    )
  }
  invisible(TRUE)
}

gp_hybrid_validate_hybrid_keyed_matrix <- function(x, arg_name, hybrid_ids, mode = NULL) {
  if (is.null(x)) {
    return(invisible(FALSE))
  }
  rn <- rownames(x)
  if (is.null(rn) || !length(rn)) {
    gp_hybrid_stop_contract(
      details = paste(arg_name, "must have row names keyed by HybridID."),
      mode = mode
    )
  }
  overlap <- sum(unique(as.character(hybrid_ids)) %in% rn)
  if (!overlap) {
    gp_hybrid_stop_contract(
      details = paste(
        arg_name,
        "does not overlap HybridID values in the phenotype table.",
        "Hybrid ML/DL requires hybrid-level geno_data keyed by HybridID."
      ),
      mode = mode
    )
  }
  invisible(TRUE)
}

gp_hybrid_allowed_models_for_mode <- function(mode) {
  switch(
    mode,
    hybrid_asreml = gp_hybrid_asreml_supported_models(),
    hybrid_bayes = gp_hybrid_bayes_supported_models(),
    hybrid_gp = gp_hybrid_gp_supported_models(),
    hybrid_ml = gp_hybrid_ml_supported_models(),
    hybrid_dl = gp_hybrid_dl_supported_models(),
    character()
  )
}

gp_hybrid_model_scope_label <- function(mode) {
  switch(
    mode,
    hybrid_asreml = "Hybrid ASReml-R",
    hybrid_bayes = "Hybrid Bayesian",
    hybrid_gp = "Hybrid GP",
    hybrid_ml = "Hybrid ML",
    hybrid_dl = "Hybrid DL",
    "Hybrid prediction"
  )
}

gp_validate_hybrid_model_selection <- function(mode,
                                               GS_model = NULL,
                                               GS_model_cv = NULL,
                                               cross_validation = FALSE) {
  selected_model <- gp_selected_model_for_stage(
    GS_model = GS_model,
    GS_model_cv = GS_model_cv,
    cross_validation = cross_validation
  )
  selected_model <- gp_canonicalize_supported_model_names(selected_model)
  if (!length(selected_model)) {
    return(invisible(TRUE))
  }
  if (length(selected_model) != 1L) {
    gp_hybrid_stop_contract(
      details = paste(
        gp_hybrid_model_scope_label(mode),
        "currently expects one model at a time. Provided:",
        paste(selected_model, collapse = ", ")
      ),
      mode = mode
    )
  }
  allowed <- gp_hybrid_allowed_models_for_mode(mode)
  if (!length(allowed) || !selected_model %in% allowed) {
    allowed_display <- gp_display_supported_model_names(allowed)
    gp_hybrid_stop_contract(
      details = paste(
        gp_hybrid_model_scope_label(mode),
        "supports only:",
        paste(allowed_display, collapse = ", ")
      ),
      mode = mode
    )
  }
  invisible(TRUE)
}

#' Validate the hybrid prediction input standard
#'
#' @param pheno_data Optional phenotype data frame.
#' @param pheno_data_train Optional pre-split training phenotype data.
#' @param pheno_data_test Optional pre-split testing phenotype data.
#' @param geno_data Optional shared genomic marker matrix.
#' @param gmatrix Optional shared relationship matrix.
#' @param female_gmatrix Optional female-parent relationship matrix.
#' @param male_gmatrix Optional male-parent relationship matrix.
#' @param female_geno_data Optional female-parent genomic marker matrix.
#' @param male_geno_data Optional male-parent genomic marker matrix.
#' @param response Response column names.
#' @param gen_name Hybrid identifier column name.
#' @param female_parent Female parent column name.
#' @param male_parent Male parent column name.
#' @param heter_groups Optional environment/grouping column name.
#' @param hybrid_asreml Logical indicating whether hybrid ASReml mode is active.
#' @param hybrid_bayes Logical indicating whether hybrid Bayesian mode is active.
#' @param hybrid_gp Logical indicating whether hybrid GP mode is active.
#' @param hybrid_ml Logical indicating whether hybrid ML mode is active.
#' @param hybrid_dl Logical indicating whether hybrid DL mode is active.
#' @param GS_model Optional selected true-prediction model.
#' @param GS_model_cv Optional selected cross-validation model.
#' @param cross_validation Logical indicating whether cross-validation is active.
#' @param cross_validation_meth Optional cross-validation method. Hybrid MET
#'   workflows require CV0, CV1, CV2, or a repeated variant.
#'
#' @return Invisibly returns a validation list when inputs satisfy the contract.
#' @export
validate_hybrid_input_standard <- function(pheno_data = NULL,
                                           pheno_data_train = NULL,
                                           pheno_data_test = NULL,
                                           geno_data = NULL,
                                           gmatrix = NULL,
                                           female_gmatrix = NULL,
                                           male_gmatrix = NULL,
                                           female_geno_data = NULL,
                                           male_geno_data = NULL,
                                           response = NULL,
                                           gen_name = NULL,
                                           female_parent = NULL,
                                           male_parent = NULL,
                                           heter_groups = NULL,
                                           hybrid_asreml = FALSE,
                                           hybrid_bayes = FALSE,
                                           hybrid_gp = FALSE,
                                           hybrid_ml = FALSE,
                                           hybrid_dl = FALSE,
                                           GS_model = NULL,
                                           GS_model_cv = NULL,
                                           cross_validation = FALSE,
                                           cross_validation_meth = NULL) {
  active_mode <- gp_hybrid_assert_single_path(
    hybrid_asreml = hybrid_asreml,
    hybrid_bayes = hybrid_bayes,
    hybrid_gp = hybrid_gp,
    hybrid_ml = hybrid_ml,
    hybrid_dl = hybrid_dl
  )
  if (!length(active_mode)) {
    return(invisible(list(valid = TRUE, mode = NULL, standard = hybrid_data_standard())))
  }
  gp_validate_hybrid_model_selection(
    mode = active_mode,
    GS_model = GS_model,
    GS_model_cv = GS_model_cv,
    cross_validation = cross_validation
  )

  ph <- gp_hybrid_pick_pheno_input(
    pheno_data = pheno_data,
    pheno_data_train = pheno_data_train,
    pheno_data_test = pheno_data_test
  )
  ph <- gp_hybrid_validate_phenotype_contract(
    ph = ph,
    response = response,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    heter_groups = heter_groups,
    mode = active_mode
  )

  is_met <- !is.null(heter_groups) && length(heter_groups) == 1L &&
    !is.na(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(ph) &&
    length(unique(stats::na.omit(as.character(ph[[heter_groups]])))) > 1L
  if (isTRUE(is_met) && isTRUE(cross_validation)) {
    method <- normalize_cv_token(cross_validation_meth %||% "")
    allowed <- c(
      "cv0", "cv1", "cv2",
      "repeated_cv0", "repeated_cv1", "repeated_cv2"
    )
    if (!method %in% allowed) {
      gp_hybrid_stop_contract(
        details = paste(
          "Hybrid multi-environment cross-validation requires CV0, CV1, CV2,",
          "Repeated_CV0, Repeated_CV1, or Repeated_CV2."
        ),
        mode = active_mode
      )
    }
  }

  if (active_mode %in% c("hybrid_asreml", "hybrid_bayes", "hybrid_gp")) {
    parent_inputs_present <- sum(!vapply(
      list(female_geno_data, male_geno_data, female_gmatrix, male_gmatrix, geno_data, gmatrix),
      is.null, logical(1)
    ))
    if (!parent_inputs_present) {
      gp_hybrid_stop_contract(
        details = paste(
          active_mode,
          "requires parent-resolved genomic input.",
          "Provide female_geno_data + male_geno_data, female_gmatrix + male_gmatrix, or a shared parent-keyed geno_data/gmatrix."
        ),
        mode = active_mode
      )
    }
    gp_hybrid_validate_parent_keyed_matrix(female_geno_data, "female_geno_data", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
    gp_hybrid_validate_parent_keyed_matrix(male_geno_data, "male_geno_data", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
    gp_hybrid_validate_parent_keyed_matrix(female_gmatrix, "female_gmatrix", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
    gp_hybrid_validate_parent_keyed_matrix(male_gmatrix, "male_gmatrix", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
    gp_hybrid_validate_parent_keyed_matrix(geno_data, "geno_data", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
    gp_hybrid_validate_parent_keyed_matrix(gmatrix, "gmatrix", ph[[female_parent]], ph[[male_parent]], mode = active_mode)
  }

  if (active_mode %in% c("hybrid_ml", "hybrid_dl")) {
    if (is.null(geno_data)) {
      gp_hybrid_stop_contract(
        details = paste(
          active_mode,
          "requires hybrid-level geno_data keyed by HybridID."
        ),
        mode = active_mode
      )
    }
    gp_hybrid_validate_hybrid_keyed_matrix(geno_data, "geno_data", ph[[gen_name]], mode = active_mode)
  }

  invisible(list(valid = TRUE, mode = active_mode, standard = hybrid_data_standard()))
}
