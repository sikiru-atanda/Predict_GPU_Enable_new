# Function to add missing hybrids to phenotypic data dynamically
add_extra_gid_from_geno_omic_to_pheno <- function(geno_data, pheno_data, gen_name, response_var, heter_group = NULL) {
  # Identify genotypes in the genotypic data but not in the phenotypic data
  msg <- ""

  if (!all(response_var %in% colnames(pheno_data))) {
         stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
     }

  diff_gid <- setdiff(rownames(geno_data), unique(pheno_data[[gen_name]]))
#   # Ensure response_var exist in df
#   if (!all(response_var %in% colnames(pheno_data))) {
#     stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
#   }
#
#   # Create a logical matrix indicating NA positions for response variables
#   na_matrix <- is.na(pheno_data[response_var])
#
#   # Check if all rows have the same NA pattern
#   if (!all(rowSums(na_matrix) %in% c(0, length(response_var)))) {
#     stop(paste(msg, "Not all response variable columns have NAs in the same positions."), call. = FALSE)
#   }
#   if (all(na_matrix == FALSE)){
#   diff_gid <- setdiff(rownames(geno_data), unique(pheno_data[[gen_name]]))
#   } else{
#     diff_gid <- NULL
#     return(TRUE)
# }
  # If no missing genotypes, return the original phenotypic data
  if (length(diff_gid) == 0) {
    return(pheno_data)
  }

  # Dynamically determine column names of the phenotypic data
  col_names <- colnames(pheno_data)

  # Create a base data frame for the missing genotypes
  if (!is.null(heter_group)) {

    if(heter_group %in% col_names){
    # For multi-environment data: Expand genotypes to all environments
    environments <- as.character(unique(pheno_data[[heter_group]]))
    diff_data <- expand.grid(
      gen_name = diff_gid,
      heter_group = environments,
      stringsAsFactors = FALSE
    )

    colnames(diff_data) <- c(gen_name, heter_group)
    # Add remaining columns dynamically and fill them with NA
    for (col in setdiff(col_names, names(diff_data))) {
      diff_data[[col]] <- NA
    }

    # Ensure column order matches the original phenotypic data
    diff_data <- diff_data[, col_names, drop = FALSE]

    }
  } else {
    # For single-environment data: Create a data frame with all columns as NA
    diff_data <- data.frame(matrix(NA, nrow = length(diff_gid), ncol = length(col_names)))
    colnames(diff_data) <- col_names

    # Populate the `gen_name` column with missing genotypes
    diff_data[[gen_name]] <- diff_gid

    # Populate the `response_var` column with NA for the missing genotypes
    if (all(response_var %in% col_names)) {
      diff_data[response_var] <- lapply(diff_data[response_var], function(x) NA)
    }
  }

  # Combine the phenotypic data with the new rows
  pheno_data <- rbind(pheno_data, diff_data)

  return(pheno_data)
}

# Multi-environment true prediction: add one NA record for every line x
# environment combination that has no record, so every engine predicts the
# full grid (BGLR, GP and ML/DL predict NA records; ASReml already predicted
# the grid). Columns constant within an environment (or within a line) are
# filled from that environment (line); a weights column gets 1. If a
# fixed-effect covariate cannot be filled this way the grid is not completed
# (message) rather than guessed.
gp_met_complete_grid <- function(pheno_data, gen_name, heter_groups, response,
                                 fixed = NULL, weights = NULL) {
  if (is.null(heter_groups) || !heter_groups %in% names(pheno_data)) return(pheno_data)
  ids <- as.character(pheno_data[[gen_name]])
  envs <- as.character(pheno_data[[heter_groups]])
  present <- unique(paste(ids, envs, sep = "\r"))
  grid <- expand.grid(gid = unique(ids), env = unique(envs), stringsAsFactors = FALSE)
  missing <- grid[!paste(grid$gid, grid$env, sep = "\r") %in% present, , drop = FALSE]
  if (!nrow(missing)) return(pheno_data)

  add <- pheno_data[rep(NA_integer_, nrow(missing)), , drop = FALSE]
  rownames(add) <- NULL
  add[[gen_name]] <- missing$gid
  add[[heter_groups]] <- missing$env
  for (tr in intersect(response, names(add))) add[[tr]] <- NA
  constant_within <- function(col, key) {
    vals <- split(pheno_data[[col]], key)
    all(vapply(vals, function(v) length(unique(v[!is.na(v)])) <= 1L, logical(1)))
  }
  other_cols <- setdiff(names(pheno_data), c(gen_name, heter_groups, response))
  for (col in other_cols) {
    if (!is.null(weights) && identical(col, weights)) {
      add[[col]] <- 1
    } else if (constant_within(col, envs)) {
      add[[col]] <- pheno_data[[col]][match(missing$env, envs)]
    } else if (constant_within(col, ids)) {
      add[[col]] <- pheno_data[[col]][match(missing$gid, ids)]
    }
  }
  fixed_vars <- if (inherits(fixed, "formula")) setdiff(all.vars(fixed), c(gen_name, heter_groups)) else character()
  unfilled <- fixed_vars[fixed_vars %in% names(add) & vapply(fixed_vars, function(v) v %in% names(add) && anyNA(add[[v]]), logical(1))]
  if (length(unfilled)) {
    message("Predictions are given for observed line x ", heter_groups, " combinations only: fixed effect(s) ",
            paste(unfilled, collapse = ", "), " vary within ", heter_groups, " and cannot be filled for the other combinations.")
    return(pheno_data)
  }
  out <- rbind(pheno_data, add)
  rownames(out) <- NULL
  out
}

# Train_Test_Label for a multi-environment prediction table:
#   "Train"      the line was observed in that environment;
#   "Test"       a record the user supplied with NA (the line x environment
#                was in the input) or a line with no observation anywhere;
#   "Unobserved" a combination added to complete the grid (not in the input)
#                for a line observed in another environment.
# input_keys: "<line>\r<environment>" of the input records (before the grid
# was completed); NULL when nothing was added.
gp_met_relabel_predictions <- function(pred, pheno_data, gen_name, heter_groups, response,
                                       input_keys = NULL) {
  if (!is.data.frame(pred) || !"Train_Test_Label" %in% names(pred)) return(pred)
  env_col <- intersect(c(heter_groups, "Env"), names(pred))[1L]
  gid_col <- intersect(c(gen_name, "GID"), names(pred))[1L]
  if (is.na(env_col) || is.na(gid_col) || !heter_groups %in% names(pheno_data)) return(pred)
  trait_col <- if ("Trait" %in% names(pred)) "Trait" else NULL
  label_for <- function(tr, rows) {
    y <- pheno_data[[tr]]
    seen <- !is.na(y)
    obs_keys <- unique(paste(pheno_data[[gen_name]][seen], pheno_data[[heter_groups]][seen], sep = "\r"))
    observed_lines <- unique(as.character(pheno_data[[gen_name]][seen]))
    key <- paste(pred[[gid_col]][rows], pred[[env_col]][rows], sep = "\r")
    added <- if (is.null(input_keys)) rep(FALSE, length(key)) else !key %in% input_keys
    ifelse(key %in% obs_keys, "Train",
           ifelse(added & as.character(pred[[gid_col]][rows]) %in% observed_lines, "Unobserved", "Test"))
  }
  labels <- as.character(pred[["Train_Test_Label"]])
  if (is.null(trait_col)) {
    tr <- intersect(response, names(pheno_data))[1L]
    if (is.na(tr)) return(pred)
    labels <- label_for(tr, seq_len(nrow(pred)))
  } else {
    for (tr in intersect(unique(as.character(pred[[trait_col]])), names(pheno_data))) {
      rows <- which(as.character(pred[[trait_col]]) == tr)
      labels[rows] <- label_for(tr, rows)
    }
  }
  pred[["Train_Test_Label"]] <- labels
  if ("Train_Test" %in% names(pred)) pred[["Train_Test"]] <- labels
  pred
}
