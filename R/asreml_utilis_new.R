#' Warn when an ASReml model has not converged after update retries
#'
#' Internal helper. The REML fitting paths retry [asreml::update.asreml] a fixed
#' number of times; if the model still reports `converge == FALSE`, the variance
#' components, breeding values and prediction error variances it produces are not
#' trustworthy. This emits a single informative warning in that case. An unknown
#' convergence state (`NULL`/missing `converge`, e.g. a failed fit returning
#' `NULL`) is treated as "do not warn" to avoid spurious alarms.
#'
#' @param mod A fitted asreml model (or any list with a `converge` element).
#' @param context Optional short string naming the fitting context; included in
#'   the warning message to help the user locate the offending model.
#' @return `mod`, invisibly.
#' @keywords internal
#' @noRd
gp_asreml_warn_if_not_converged <- function(mod, context = NULL) {
  conv <- tryCatch(mod$converge, error = function(e) NULL)
  if (is.logical(conv) && length(conv) == 1L && !is.na(conv) && !conv) {
    where <- if (!is.null(context) && nzchar(context)) paste0(" (", context, ")") else ""
    warning(sprintf(
      paste0("ASReml model%s did not converge after update retries; its variance ",
             "components, breeding values and prediction error variances may be ",
             "unreliable. Consider revising the model (random/residual structure, ",
             "variance-covariance structure) or supplying better starting values."),
      where), call. = FALSE)
  }
  invisible(mod)
}

#' Decide whether an ASReml fit needs another bounded update
#'
#' A fit may report `converge = TRUE` while one or more estimable variance
#' parameters still changed materially on the final iteration. The public
#' output guard rejects that state, so fitting must refine it before applying
#' the guard rather than stopping solely on `converge`.
#'
#' @param mod A fitted ASReml model.
#' @param max_percent_change Largest accepted final percent change for an
#'   estimable, non-boundary variance parameter.
#' @return A single logical value.
#' @keywords internal
#' @noRd
gp_asreml_fit_requires_refinement <- function(mod, max_percent_change = 1) {
  if (is.null(mod)) return(FALSE)
  converge <- tryCatch(mod$converge, error = function(e) NULL)
  if (is.logical(converge) && length(converge) == 1L &&
      !is.na(converge) && !converge) {
    return(TRUE)
  }
  ifault <- suppressWarnings(as.integer(
    tryCatch(mod$ifault, error = function(e) NA_integer_)
  )[1L])
  if (is.finite(ifault) && ifault != 0L) return(TRUE)

  percent_change <- suppressWarnings(as.numeric(
    tryCatch(mod$vparameters.pc, error = function(e) numeric())
  ))
  if (!length(percent_change)) return(FALSE)
  constraint <- as.character(
    tryCatch(mod$vparameters.con, error = function(e) character())
  )
  if (length(constraint) != length(percent_change)) {
    constraint <- rep(NA_character_, length(percent_change))
  }
  estimable <- is.finite(percent_change) & !constraint %in% c("B", "F")
  any(estimable & abs(percent_change) > max_percent_change)
}

#' Refine an ASReml fit until both convergence diagnostics are stable
#'
#' @param mod A fitted ASReml model.
#' @param max_updates Maximum number of update attempts.
#' @param max_percent_change Passed to
#'   [gp_asreml_fit_requires_refinement()].
#' @param update_fn Optional update function, used by unit tests; defaults to
#'   [asreml::update.asreml()].
#' @return The last fitted model.
#' @keywords internal
#' @noRd
gp_asreml_refine_fit <- function(mod,
                                 max_updates = 8L,
                                 max_percent_change = 1,
                                 update_fn = NULL) {
  max_updates <- suppressWarnings(as.integer(max_updates)[1L])
  if (!is.finite(max_updates) || max_updates < 0L) {
    stop("max_updates must be a non-negative integer.", call. = FALSE)
  }
  if (is.null(update_fn)) {
    update_fn <- function(x) asreml::update.asreml(x)
  }
  if (!is.function(update_fn)) {
    stop("update_fn must be a function.", call. = FALSE)
  }
  if (max_updates == 0L) return(mod)
  for (attempt in seq_len(max_updates)) {
    if (!gp_asreml_fit_requires_refinement(
      mod, max_percent_change = max_percent_change
    )) break
    mod <- update_fn(mod)
  }
  mod
}

#' Reject an ASReml fit whose model-implied quantities are not trustworthy
#'
#' Prediction, PEV, variance, covariance and correlation output must never be
#' produced from a failed or materially unfinished ASReml optimization. ASReml
#' can report `converge = TRUE` while still warning that a variance component
#' changed by more than one percent on the final iteration, so both states are
#' checked. Boundary parameters are allowed; the reconstructed covariance
#' matrix is validated separately by the structure-specific extractor.
#'
#' @param mod A fitted ASReml model.
#' @param context Short label included in an error message.
#' @param max_percent_change Largest accepted final percent change for a fitted
#'   non-boundary variance parameter.
#' @return `mod`, invisibly.
#' @keywords internal
#' @noRd
gp_asreml_assert_trustworthy_fit <- function(mod,
                                             context = "ASReml model",
                                             max_percent_change = 1) {
  if (is.null(mod)) {
    stop(context, " did not return a fitted model.", call. = FALSE)
  }
  converge <- tryCatch(mod$converge, error = function(e) NULL)
  if (!isTRUE(converge)) {
    stop(
      context,
      " did not converge; prediction, PEV and variance-covariance output will not be exported.",
      call. = FALSE
    )
  }
  ifault <- suppressWarnings(as.integer(tryCatch(mod$ifault, error = function(e) NA_integer_))[1L])
  if (is.finite(ifault) && ifault != 0L) {
    stop(
      context, " has ASReml ifault = ", ifault,
      "; prediction and variance-covariance output will not be exported.",
      call. = FALSE
    )
  }
  percent_change <- suppressWarnings(as.numeric(
    tryCatch(mod$vparameters.pc, error = function(e) numeric())
  ))
  constraint <- as.character(
    tryCatch(mod$vparameters.con, error = function(e) character())
  )
  if (length(constraint) != length(percent_change)) {
    constraint <- rep(NA_character_, length(percent_change))
  }
  # B/F parameters are fixed at a boundary or by the model and ASReml reports
  # no meaningful final update for them. Test only estimable finite updates.
  estimable <- is.finite(percent_change) & !constraint %in% c("B", "F")
  excessive_change <- estimable & abs(percent_change) > max_percent_change
  if (any(excessive_change)) {
    offending <- names(mod$vparameters.pc)[excessive_change]
    stop(
      context, " has variance component(s) changing by more than ",
      format(max_percent_change, trim = TRUE), "% on the final iteration",
      if (length(offending)) paste0(": ", paste(offending, collapse = ", ")) else ".",
      " Output is withheld until the model is stable.",
      call. = FALSE
    )
  }
  invisible(mod)
}

#' Delta-method standard error of heritability via vpredict
#'
#' Internal helper. Given a fitted ASReml model and the variance-component row
#' names that form the genetic numerator and the (single) residual term of a
#' heritability ratio
#'   h2 = (Vg_1 + ... + Vg_k) / (Vg_1 + ... + Vg_k + Ve),
#' this returns the delta-method standard error of h2 from
#' [asreml::vpredict]. `vpredict` indexes components as `V<i>` by their row order
#' in `summary(model)$varcomp`, so the supplied row names are matched against
#' that table to build the formula. This is exact only when the numerator and
#' denominator are sums of raw variance components that each correspond to a
#' single varcomp row (e.g. compound-symmetry or `corgh` diagonal variances).
#'
#' It is deliberately conservative: if `asreml` is unavailable, any requested
#' component cannot be matched, or `vpredict` errors (e.g. a factor-analytic
#' structure where the per-environment genetic variance is a non-linear function
#' of loadings, which this linear-ratio formula does NOT represent), it returns
#' `NA_real_` rather than a wrong number.
#'
#' @param model A fitted asreml model.
#' @param vg_rownames Character vector of `summary(model)$varcomp` row names whose
#'   components sum to the genetic variance numerator.
#' @param ve_rowname Single `summary(model)$varcomp` row name for the residual
#'   variance of this heritability ratio.
#' @return The standard error of `h2` as a length-1 numeric, or `NA_real_`.
#' @keywords internal
#' @noRd
gp_asreml_heritability_se <- function(model, vg_rownames, ve_rowname) {
  out <- tryCatch({
    if (!requireNamespace("asreml", quietly = TRUE)) return(NA_real_)
    if (is.null(model) || length(vg_rownames) < 1L || length(ve_rowname) != 1L) return(NA_real_)
    vc <- summary(model)$varcomp
    if (is.null(vc) || is.null(rownames(vc))) return(NA_real_)
    rn <- rownames(vc)
    ig <- match(vg_rownames, rn)
    ie <- match(ve_rowname, rn)
    if (anyNA(ig) || anyNA(ie)) return(NA_real_)
    num <- paste(sprintf("V%d", ig), collapse = " + ")
    den <- paste(sprintf("V%d", c(ig, ie)), collapse = " + ")
    form <- stats::as.formula(sprintf("h2 ~ (%s) / (%s)", num, den))
    res <- asreml::vpredict(model, form)
    as.numeric(res$SE[1L])
  }, error = function(e) NA_real_)
  if (length(out) != 1L || !is.finite(out)) NA_real_ else out
}

#' Per-environment heritability SE for env-structured genetic variances
#'
#' Internal helper for MET heritability with a variance-covariance structure
#' whose diagonal entries are genuine per-environment variance-component rows
#' (e.g. `corgh`, `us`). For environment `env_label` it locates that env's
#' genetic-variance row(s) (one per omic in `names_in_inv_list`, correlation/
#' covariance rows excluded) and residual row, then returns the delta-method SE
#' of h2 from [gp_asreml_heritability_se].
#'
#' Crucially it is **self-verifying**: it only returns an SE when the matched
#' components actually reproduce the point estimate's genetic numerator
#' (`target_vg`) and residual (`target_ve`). For a factor-analytic structure the
#' per-env genetic variance is `Lambda Lambda' + psi` (a non-linear function of
#' several parameters) and matches no single row, so the check fails and the
#' helper returns `NA_real_` rather than a wrong SE.
#'
#' @param model A fitted asreml model.
#' @param vc The variance-component table (`summary(model)$varcomp`), passed in to
#'   avoid recomputing it per environment.
#' @param names_in_inv_list Inverse-relationship object name(s) that appear in the
#'   genetic variance row names.
#' @param env_label The environment level whose heritability SE is wanted.
#' @param target_vg Point-estimate genetic variance numerator for this env.
#' @param target_ve Point-estimate residual variance for this env.
#' @param tol Relative tolerance for the self-check.
#' @return The h2 standard error as a length-1 numeric, or `NA_real_`.
#' @keywords internal
#' @noRd
gp_asreml_h2_se_for_env <- function(model, vc, names_in_inv_list, env_label,
                                    target_vg, target_ve, tol = 1e-4) {
  out <- tryCatch({
    if (is.null(vc) || is.null(rownames(vc))) return(NA_real_)
    rn <- rownames(vc)
    is_cov <- grepl("\\.cor$|!cor$|\\.cov$", rn)            # correlation/covariance terms
    elab <- paste0(gsub("([\\W])", "\\\\\\1", env_label, perl = TRUE))  # escape regex metachars
    # genetic variance row(s) for this env: match the inverse name AND end in the
    # env label, excluding correlation/covariance entries.
    g_rows <- unlist(lapply(names_in_inv_list, function(nm) {
      hit <- grepl(nm, rn, fixed = TRUE) & grepl(paste0(elab, "$"), rn) & !is_cov
      rn[hit]
    }))
    if (!length(g_rows)) return(NA_real_)
    # residual row for this env, or the single shared residual if homogeneous.
    r_hit <- which(grepl(paste0(elab, "!R$"), rn))
    if (!length(r_hit)) r_hit <- which(grepl("!R$", rn))
    if (length(r_hit) != 1L) return(NA_real_)
    e_row <- rn[r_hit]
    # Self-check: matched components must reproduce the point-estimate terms.
    # Tolerance is relative to the target (plus a tiny absolute floor) so the
    # guard stays sensitive even when a variance is near zero -- a factor-analytic
    # diagonal differs from the matched rows by far more than this.
    if (abs(sum(vc[g_rows, "component"]) - target_vg) > tol * abs(target_vg) + 1e-8) return(NA_real_)
    if (abs(vc[e_row, "component"] - target_ve) > tol * abs(target_ve) + 1e-8) return(NA_real_)
    gp_asreml_heritability_se(model, vg_rownames = g_rows, ve_rowname = e_row)
  }, error = function(e) NA_real_)
  if (length(out) != 1L || !is.finite(out)) NA_real_ else out
}

#' Compute Inverse and Sparse Matrix from a Kernel Matrix
#'
#' This function computes the inverse of a given kernel matrix using a regularization
#' technique to handle nearly positive definite matrices with small negative eigenvalues.
#' It can also return a sparse representation of the kernel matrix.
#'
#' @param kernel A numeric square matrix representing the kernel from which to compute the inverse.
#' @param epsilon A small positive value (regularization parameter) added to the diagonal elements
#'        to ensure the matrix is positive definite. Default is 1e-6. adjust as needed
#' @param inverse Logical indicating whether to compute and return the inverse of the kernel matrix.
#'        If FALSE, returns a sparse matrix representation of the original kernel matrix. Default is NULL.
#'
#' @return If `inverse = TRUE`, returns a matrix which is the inverse of the input kernel matrix,
#'         with additional attributes "rowNames", "colNames", and "INVERSE" set to TRUE.
#'         If `inverse = FALSE`, returns a sparse matrix representation of the kernel matrix,
#'         with "rowNames", "colNames", and "INVERSE" set to FALSE. If an error occurs during
#'         inversion, it returns `NULL` and prints an error message.
#'
#' @examples
#' # Create a symmetric positive definite matrix
#' kernel_matrix <- matrix(c(2, -1, -1, 2), ncol = 2)
#' rownames(kernel_matrix) <- colnames(kernel_matrix) <- c("A", "B")
#'
#' # Compute the inverse
#' inverse_matrix <- compute_inverse_and_sparse(kernel = kernel_matrix, inverse = TRUE)
#'
#' # Compute a sparse representation
#' sparse_matrix <- compute_inverse_and_sparse(kernel = kernel_matrix, inverse = FALSE)
#'
#' @export

compute_inverse_and_sparse <- function(kernel = NULL,
                                       epsilon = 1e-6,
                                       inverse =NULL) {
  attr(kernel, "rowNames") <- rownames(kernel)
  attr(kernel, "colNames") <- colnames(kernel)
  kernel_extra <- kernel

  if(isTRUE(inverse)){
    # Apply regularization technique when the matrix is nearly positive definite
    # but has small negative eigenvalues due to numerical precision issues
    # to make it positive definite by adding
    ## small positive constant to the diagonal elements of the matrix
    ## Two ways of acheiving that are define here
    # 1) Apply it directly to the directly to the diagonal
    # elements of the matrix is indeed a form of regularization
    # 2) Ridge regularization
    result <- tryCatch(
      {
        diag(kernel) <- diag(kernel) + epsilon
        inverse_matrix <- chol2inv(chol(kernel))

        # Return the result
        inverse_matrix
      },
      error = function(e) {
        cat("Error occurred during computation:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    # Check if an error occurred
    #if(inherits(result, "try-error")) {
    if(is.null(result)) {
      ### # Add a small positive constant to the diagonal elements
      # Modify the kernel if needed
      kernel <- kernel + epsilon * diag(nrow(kernel))
      inverse_matrix <- chol2inv(chol(kernel))
    } else {
      inverse_matrix <- result
    }

    sparse <- sparse_matrix(grm_kernel_data = inverse_matrix)
    attr(sparse, "rowNames") <- rownames(kernel)
    attr(sparse, "colNames") <- colnames(kernel)
    attr(sparse, "INVERSE") <- TRUE
  } else {
    sparse <- sparse_matrix(grm_kernel_data = kernel)
    attr(sparse, "rowNames") <- rownames(kernel)
    attr(sparse, "colNames") <- colnames(kernel)
    attr(sparse, "INVERSE") <- FALSE
  }

  return(sparse)
}

asreml_normalize_formula_arg <- function(x, labels = character()) {
  if (is.null(x)) {
    return(NULL)
  }
  if (length(labels) && is.character(x) && length(x) == 1L) {
    x_chr <- toupper(trimws(as.character(x)))
    if (!nzchar(x_chr) || is.na(x_chr) || x_chr %in% toupper(labels)) {
      return(NULL)
    }
  }
  x
}

asreml_rhs_terms <- function(x, labels = character()) {
  x <- asreml_normalize_formula_arg(x, labels = labels)
  if (is.null(x)) {
    return(character())
  }
  if (!inherits(x, "formula")) {
    msg <- ""
    stop(print(paste(msg, "The fixed/covariate term is not a class of type 'formula'. Example: fixed = ~ X + Y")),
         call. = FALSE)
  }
  terms <- all.vars(x)
  terms <- unique(trimws(as.character(terms)))
  terms[nzchar(terms)]
}

asreml_append_rhs_terms <- function(code, terms) {
  terms <- unique(trimws(as.character(terms %||% character())))
  terms <- terms[nzchar(terms)]
  if (!length(terms)) {
    return(code)
  }
  for (term in terms) {
    code <- paste(code, term, sep = "+")
  }
  code
}

asreml_clean_fixed_code <- function(code) {
  code <- gsub("\\+\\s*$", "", code)
  code <- gsub("\\+\\s*,", ",", code)
  code
}

# Function to check if current order of rownames or colnames is the same as unique_GIDs
## in pheno_data
# is_in_correct_order <- function(names, order) {
#   if (length(names) != length(order)) {
#     return(FALSE)  # Different lengths mean they are not in the same order
#   }
#   all(names == order)
# }

#' Adjust Random Terms for ASReml Model Fitting
#'
#' This function prepares and adjusts random terms for fitting an ASReml model,
#' facilitating the modeling of genetic data with potential heterogeneity and
#' specific variance-covariance structures. It allows for the specification of
#' fixed and random effects, interaction terms, and various models of heterogeneity
#' and variance-covariance among genetic components.
#'
#' @param random A formula specifying the random effects to be included in the model.
#' @param fixed A formula specifying the fixed effects to be included in the model.
#' @param fixed_term A character vector specifying fixed terms to be considered in the model.
#' @param heter_groups A character string identifying the grouping variable for heterogeneity.
#' @param heter_resid A logical indicating if heterogeneity in residuals should be considered.
#' @param var_cov_str A character string indicating the type of variance-covariance structure.
#' @param code_asr A character vector for storing ASReml code or model specification.
#' @param names_in_inv_list A character vector of names identifying inverse matrices in the model,
#' related to different genetic components.
#' @param gen_name A character string specifying the name of the genetic factor in the model.
#' @param pheno_data A data frame containing phenotypic data used in the model.
#' @param kernel_list Optional named list of additional relationship or kernel matrices.
#' @param weights Optional positive Stage 2 observation precisions, supplied as
#'   a numeric vector, one-column table, or column name in `pheno_data`. ASReml
#'   uses an internal data column and `asr_gaussian(dispersion = 1)`, so the
#'   residual variance for row `i` is exactly `1 / weights[i]`.
#' @param ... Additional arguments for future use or extensions.
#'
#' @return A list containing elements critical for ASReml model specification,
#' including adjusted code for model fitting, positions of genetic and interaction terms,
#' and the defined random terms.
#'
#' @details
#' The function is particularly useful for complex genetic models requiring precise
#' specification of random effects, handling of heterogeneity across groups, and
#' incorporation of specific variance-covariance structures. It preprocesses the input
#' model terms to ensure compatibility with ASReml software requirements and facilitates
#' the inclusion of genetic interactions and heterogeneity terms.
#'
#' @examples
#' # Assuming a hypothetical ASReml model setup:
#' random_effects <- ~ Genotype + Genotype:Environment
#' fixed_effects <- ~ Environment + Treatment
#' pheno <- data.frame(Genotype = factor(rep(1:10, each = 6)),
#'                     Environment = factor(rep(1:3, times = 20)),
#'                     Treatment = factor(rep(c("Control", "Treated"), each = 30)),
#'                     Yield = rnorm(60, mean = 100, sd = 15))
#'
#' model_setup <- random_terms_fit_new(random = random_effects,
#'                                     fixed = fixed_effects,
#'                                     heter_groups = "Environment",
#'                                     gen_name = "Genotype",
#'                                     pheno_data = pheno)
#' @export

asreml_utilis_new <- function(
    fixed = NULL,
    random = NULL,
    cova=NULL,
    GS_model = NULL,
    response = NULL,
    pheno_data = NULL,
    gmatrix = NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    kernel_list = NULL,
    inverse = NULL,
    epsilon = TRUE,
    gen_name = NULL,
    heter_groups = NULL,
    heter_resid = FALSE,
    var_cov_str = NULL,
    weights = NULL,
    workspace=1e08,
    engine = NULL,
    pworkspace= 1e06,
    maxit = 50,
    cross_validation = FALSE,
    ...
) {

  msg <- ""

  pheno_data <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  class(pheno_data) <- "data.frame"
  if (!is.null(gen_name) && gen_name %in% names(pheno_data)) {
    pheno_data[[gen_name]] <- as.factor(pheno_data[[gen_name]])
  }
  if (!is.null(heter_groups) && heter_groups %in% names(pheno_data)) {
    pheno_data[[heter_groups]] <- as.factor(pheno_data[[heter_groups]])
  }
  stage2_weights <- gp_resolve_stage2_precision_weights(
    weights = weights,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    context = "ASReml Stage 2 observation weights"
  )
  weight_column <- NULL
  if (isTRUE(stage2_weights$supplied)) {
    weight_column <- gp_stage2_weight_column()
    pheno_data[[weight_column]] <- stage2_weights$precision
  }
  heter_control <- gp_normalize_single_environment_heter_controls(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    heter_resid = heter_resid,
    var_cov_str = var_cov_str,
    response = response
  )
  heter_groups <- heter_control$heter_groups
  heter_resid <- heter_control$heter_resid
  var_cov_str <- heter_control$var_cov_str
  fixed <- asreml_normalize_formula_arg(fixed, labels = c("FIXED", "fixed_term_model", "NULL"))
  cova <- asreml_normalize_formula_arg(cova, labels = c("COVA", "covariate", "NULL"))

  if(engine %in% rownames(installed.packages())){
    do.call('library', list(engine))

    # if package is not installed locally then stop
  } else {

    stop(print(paste(msg,'You need to install asreml-R to use asreml-R')), call. = FALSE)
  }

  datasets <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    gkernel = gkernel,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  dataset_names <- names(datasets)

  # Extract unique GIDs from pheno_data to determine the row order
  unique_GIDs <- as.character(unique(pheno_data[[gen_name]]))

  # # Reorder the rownames and colnames of each dataset based on unique_GIDs if necessary
  # datasets <- lapply(datasets, function(mat) {
  #   correct_row_order <- is_in_correct_order(rownames(mat), unique_GIDs[unique_GIDs %in% rownames(mat)])
  #   correct_col_order <- is_in_correct_order(colnames(mat), unique_GIDs[unique_GIDs %in% colnames(mat)])
  #
  #   if (!correct_row_order || !correct_col_order) {
  #     # Reorder rows and columns if either is not in the correct order
  #     ordered_indices <- unique_GIDs[unique_GIDs %in% rownames(mat)]
  #     mat <- mat[ordered_indices, ordered_indices]
  #     message("Matrix reordered based on unique_GIDs.")
  #   } else {
  #     message("Matrix is already in the correct order; no changes made.")
  #   }
  #   return(mat)
  # })


  inv_list <- list()
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]

    if (!is.null(dataset)) {
      inv_list[[dataset_names[i]]] <- compute_inverse_and_sparse(
        kernel = dataset,
        inverse = inverse
      )
    }


  }

  if(!"gmatrix" %in%names(inv_list)){
    if(length(inv_list)>1){
      names(inv_list) <- paste(paste("omic", seq_along(inv_list), sep = ""), "inv", sep = "_")
    } else {
      names(inv_list) <- paste("omic", "inv", sep = "_")
    }
  } else if(length(grep("omic", names(inv_list)))==0 & "gmatrix" %in%names(inv_list)) {
    geno_index <- grep("gmatrix", names(inv_list))
    names(inv_list)[geno_index] <- paste("gmatrix", "inv", sep = "_")
  } else{
    if(length(grep("omic", names(inv_list)))>=1 & "gmatrix" %in%names(inv_list)) {
      omic_index <- grep("omic", names(inv_list))
      geno_index <- grep("gmatrix", names(inv_list))
      if(length(omic_index)>1){
        names(inv_list)[omic_index] <- paste(paste("omics", seq_along(omic_index), sep = ""), "inv", sep = "_")
      } else {
        names(inv_list)[omic_index] <- paste("omics", "inv", sep = "_")
      }
      ####
      names(inv_list)[geno_index] <- paste("gmatrix", "inv", sep = "_")
    }

  }

  names_in_inv_list <- names(inv_list)


  #### When the gen_name are present in more than one environment/location
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    if (!is.null(heter_resid)) {
      if (is.null(heter_groups)) {
        stop(print(paste(msg,'No column of heterogeneous groups provided.')), call. = FALSE)
      } else {

        if(!heter_groups%in%colnames(pheno_data)) {stop(print(paste(msg,'heterogenous group provided did not match column names in pheno_data')), call. = FALSE)}

        if(!all(sapply(heter_groups, function(x, pheno_data) is.factor(pheno_data[,x]),  pheno_data))) {
          pheno_data[,heter_groups] <- as.factor(pheno_data[,heter_groups])
        }
      }

    }

  }

  ##### THis is important for asreml inorder to update the model if need be
  for (i in seq_along(inv_list)) {
    assign(names(inv_list)[i], inv_list[[i]], envir = .GlobalEnv)
  }

  #a <- a + 1
  # Code Strings for all factors y= XB + UZ + e
  code_asr <- as.character()

  # colnames(pheno_data)[colnames(pheno_data) == trait] <-
  #   deparse(substitute(trait))

  code_asr[1] <- paste0(paste('asreml::asreml(fixed=', 'trait'),  '~1')
  code_asr[2] <- 'random=~'
  code_asr[3] <- 'residual=~'
  fixed_term <- NULL
  check_heter_grp_fixed <- NULL

  # Adding covariates (fixed)
  if (!is.null(cova)) {
    cova_term <- asreml_rhs_terms(cova, labels = c("COVA", "covariate", "NULL"))
    code_asr[1] <- asreml_append_rhs_terms(code_asr[1], cova_term)
  }

  # Adding fixed factors
  if (!is.null(fixed)) {
    fixed_term <- asreml_rhs_terms(fixed, labels = c("FIXED", "fixed_term_model", "NULL"))
    code_asr[1] <- asreml_append_rhs_terms(code_asr[1], fixed_term)

  } else{
    if (is.null(fixed)) {
      fixed_term = NULL
      check_heter_grp_fixed = NULL
    }
  }


  # Adding random factors

  if (is.null(random)){ stop(print('provide random term'))}

  rand_out <-  random_terms_fit_new(random = random,
                                fixed = fixed,
                                fixed_term = fixed_term,
                                heter_groups = heter_groups,
                                heter_resid = heter_resid,
                                var_cov_str = var_cov_str,
                                code_asr = code_asr,
                                names_in_inv_list = names_in_inv_list,
                                gen_name = gen_name,
                                pheno_data = pheno_data)

  code_asr_fit <-  rand_out[["code_asr"]]
  gen_pos <-  rand_out[["gen_pos"]]
  inter_gen_pos <-  rand_out[["inter_gen_pos"]]
  rand_term <- rand_out[["rand_term"]]
  ##### Ends with random term

  # Heterogeneous errors
  if (!is.null(heter_groups)&isTRUE(heter_resid)) {

    if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      #code.asr[3] <- paste(code.asr[3], paste0('dsum(~id(units)|', paste0(heter_groups, ')')), sep='')
      code_asr_fit[3] <- paste( code_asr_fit[3], paste0('dsum(~id(units)|', paste0(heter_groups, ')')), sep='')
    } else {
      warning(paste(msg,'Heterogenous residual is not possible with one environment/location. We fix it for you.'),
              call. = FALSE)
      if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
        #code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
        code_asr_fit[3] <- paste(code_asr_fit[3], 'id(units)', sep='')
      }
    }
    #code.asr[3] <- paste(code.asr[3], paste0('dsum(~idv(units)|', paste0(heter_groups, ')')), sep='')

    #code.asr[3] <- paste(code.asr[3], 'dsum(~idv(units|', 'heter_groups)', sep='')
  } else {
    #code.asr[3] <- paste(code.asr[3], 'idv(units)', sep='')
    #code.asr[3] <- paste(code.asr[3], 'id(units)', sep='')
    code_asr_fit[3] <- paste(code_asr_fit[3], 'id(units)', sep='')
  }
  # Adding traits
  #Univariate = c()
  #for (trait in 1:length(response)) {

  #if (trait==1){
  #code.asr[1] <-  gsub("trait", response[trait], code.asr[1])
  #code.asr[1] <-  gsub("trait", response, code.asr[1])

  if(isTRUE(cross_validation)){

    output <- list(code_asr_fit = code_asr_fit,
                   names_in_inv_list = names_in_inv_list ,
                   gen_pos = gen_pos,
                   inter_gen_pos = inter_gen_pos,
                   rand_term = rand_term,
                   inv_list = inv_list,
                   stage2_weight_column = weight_column,
                   stage2_weight_source = stage2_weights$source)


    return(output)
  }

  code_asr_fit[1] <-  gsub("trait", response, code_asr_fit[1])
  code_asr_fit[1] <- asreml_clean_fixed_code(code_asr_fit[1])
  #} else {

  #code.asr[1] <-  gsub(response[trait-1], response[trait], code.asr[1])
  #}


  data_cache_id <- asreml_model_data_register(pheno_data)
  data_expr <- paste0("PredictProR:::asreml_model_data_lookup('", data_cache_id, "')")
  code_asr_fit[4] <- gp_asreml_weighted_data_clause(
    data_expr = data_expr,
    weight_column = weight_column,
    response_family = "gaussian"
  )


  asreml_apply_options(
    workspace = workspace,
    pworkspace = pworkspace,
    maxit = maxit,
    engine = engine %||% "asreml"
  )
  ####
  #code.asr[1] <- paste('mod<-', code.asr[1], sep='')
  code_asr_fit[1] <- paste('mod<-', code_asr_fit[1], sep='')
  str_mod <- paste(code_asr_fit[1],code_asr_fit[2],code_asr_fit[3],code_asr_fit[4],sep=',')
  if(isTRUE(cross_validation)){
    return(str_mod)
  }
  ## Calls the current environment for evaluation. ASReml's default
  ## workspace (1e8 words = 800 MB) is too small for larger MET fits (e.g.
  ## 10 G2F environments); retry with a 4x larger workspace, at most twice
  ## and within 60% of available memory, before giving up.
  fit_env <- environment()
  ws_words <- gp_asreml_workspace_words(workspace)
  for (attempt in 0:2) {
    fit_error <- tryCatch({
      eval(parse(text = str_mod), envir = fit_env)
      NULL
    }, error = function(e) e)
    if (is.null(fit_error)) break
    if (attempt == 2L || !grepl("Insufficient workspace", conditionMessage(fit_error), fixed = TRUE)) stop(fit_error)
    next_words <- gp_asreml_next_workspace_words(ws_words)
    if (is.null(next_words)) stop(fit_error)
    message(sprintf("ASReml needed more than %.1f GB of workspace; retrying with %.1f GB (set `workspace` to skip this).",
                    ws_words * 8 / 1024^3, next_words * 8 / 1024^3))
    ws_words <- next_words
    asreml_apply_options(workspace = ws_words, pworkspace = pworkspace, maxit = maxit,
                         engine = engine %||% "asreml")
  }
  mod$call$data <- parse(text = data_expr)[[1]]
  #if (!mod$converge) { eval(parse(text='mod<-asreml::update.asreml(mod)')) }
  # Assuming `mod` is your initial model object
  mod <- gp_asreml_refine_fit(mod, max_updates = 8L)
  fit_context <- if (!is.null(heter_groups) && !is.null(var_cov_str)) {
    paste0("ASReml MET GBLUP (", var_cov_str, ")")
  } else {
    "ASReml single-trait GBLUP"
  }
  gp_asreml_warn_if_not_converged(mod, context = fit_context)
  gp_asreml_assert_trustworthy_fit(mod, context = fit_context)


  ###### Process if the model is not stable #######
  # specific_warning_occurred <- FALSE
  # error_occurred <- FALSE
  #
  # repeat {
  #   # Attempt to update the model and capture warnings
  #   tryCatch({
  #     # Assuming 'res' is your model object and update.asreml is the function you're using
  #     mod <- asreml::update.asreml(mod)
  #
  #     # If update.asreml runs without warnings or errors, we assume the update was successful
  #   }, warning = function(w) {
  #     # Check if the warning message matches the specific warning you're concerned with
  #     if(grepl("Some components changed by more than 1% on the last iteration", w$message)) {
  #       specific_warning_occurred <- TRUE
  #       # Log the occurrence of the specific warning for debugging
  #       cat("Specific warning occurred, attempting to update the model again...\n")
  #     }
  #   }, error = function(e) {
  #     error_occurred <- TRUE
  #     # Log the error for debugging
  #     cat("An error occurred: ", e$message, "\n")
  #   })
  #
  #   # Break the loop if an error occurred or if the specific warning did not occur in this iteration
  #   if (error_occurred || !specific_warning_occurred) {
  #     break
  #   }
  #
  #   # Reset the specific warning flag for the next iteration
  #   specific_warning_occurred <- FALSE
  # }


  ###################################################
  ##### Start the process of processing the results
  ################################################

  output <- list(model = mod,
                 str_mod = str_mod,
                 names_in_inv_list = names_in_inv_list ,
                 gen_pos = gen_pos,
                 inter_gen_pos = inter_gen_pos,
                 rand_term = rand_term,
                 stage2_weight_column = weight_column,
                 stage2_weight_source = stage2_weights$source)

  # names(output) <- c("model",
  #                    "str_mod",
  #                    "names_in_inv_list",
  #                    "gen_pos",
  #                    "inter_gen_pos",
  #                    "rand_term")
  #return(c(Univariate, G_list))
  return(output)

}




