#' Matrix square-root design for a GBLUP-via-BRR genetic term
#'
#' Internal helper. To fit GBLUP with BGLR's `BRR` model (genetic effect
#' `u = X beta`, `beta ~ N(0, s2 I)`), the design matrix `X` must satisfy
#' `X X' = K` (the genomic relationship/kernel), so that
#' `Cov(u) = s2 * X X' = s2 * K` -- true GBLUP. Passing the relationship matrix
#' itself as `X` (i.e. `X = K`) instead implies `Cov(u) = s2 * K K' = s2 * K^2`,
#' a mis-specified covariance that under-fits the genetic signal, inflates the
#' residual variance and biases heritability downward. This returns the
#' eigen square-root `X = U D^{1/2}` (positive eigenvalues only), which is the
#' standard GBLUP-equivalent BRR parameterization and matches the `RKHS` route.
#'
#' @param K A symmetric (relationship/kernel) matrix.
#' @param tol Relative eigenvalue tolerance below which components are dropped.
#' @return A numeric matrix `X` (n x r, r = number of retained eigenvalues) with
#'   `X %*% t(X)` approximately equal to `K`. Row names are preserved.
#' @keywords internal
#' @noRd
gp_bayes_brr_design_from_kernel <- function(K, tol = 1e-8) {
  K <- as.matrix(K)
  Ksym <- (K + t(K)) / 2
  ev <- eigen(Ksym, symmetric = TRUE)
  d <- ev$values
  keep <- d > max(d) * tol
  if (!any(keep)) return(K)               # degenerate fallback (should not happen)
  X <- ev$vectors[, keep, drop = FALSE] * rep(sqrt(d[keep]), each = nrow(ev$vectors))  # column scaling, not a dense diag product
  rownames(X) <- rownames(K)
  X
}

#' Record-level kernel Z K Z' for multi-environment Bayesian terms
#'
#' Internal helper. Each phenotype record is matched to its genotype's kernel
#' row by ID, so `K[idx, idx]` equals `Z K Z'` without forming Z. The former
#' `model.matrix(~factor(id) - 1)` incidence matrix ordered its columns
#' alphabetically and dropped genotypes without records: the product failed
#' ("non-conformable arguments") when the kernel held genotypes with no
#' phenotype, and silently paired records with the wrong kernel rows when the
#' kernel was not in alphabetical order.
#'
#' @param record_ids Genotype ID of each phenotype record, in record order.
#' @param kernel Genotype-level kernel with genotype IDs as row names.
#' @return The records x records kernel (row/column names "1".."n", as the
#'   former incidence-matrix product had).
#' @keywords internal
#' @noRd
gp_bayes_record_kernel <- function(record_ids, kernel) {
  kernel <- as.matrix(kernel)
  kernel_ids <- rownames(kernel)
  if (is.null(kernel_ids)) {
    stop("The kernel has no genotype row names; records cannot be matched to it.", call. = FALSE)
  }
  idx <- match(as.character(record_ids), kernel_ids)
  if (anyNA(idx)) {
    missing <- unique(as.character(record_ids)[is.na(idx)])
    stop(sprintf("%d genotype(s) with phenotype records are not in the kernel (e.g. %s).",
                 length(missing), paste(utils::head(missing, 3), collapse = ", ")), call. = FALSE)
  }
  out <- kernel[idx, idx, drop = FALSE]
  n <- length(idx)
  dimnames(out) <- list(as.character(seq_len(n)), as.character(seq_len(n)))
  out
}

#' Eigen-decomposition of a record-level kernel from the genotype kernel
#'
#' Internal helper. The record kernel `Z K Z'` (and the genotype-by-group
#' kernel `(Z K Z') * [same group]`, which is block-diagonal by group) has at
#' most one non-zero eigenvalue per genotype (per group). With record counts
#' `C = Z'Z`, its non-zero eigenpairs follow from the genotype-level matrix
#' `C^1/2 K C^1/2 = Y L Y'`: values `L`, vectors `V = Z C^-1/2 Y` (orthonormal).
#' This replaces an eigen-decomposition of the records x records matrix
#' (O(n_records^3); about an hour at 11,270 G2F records) by one of size
#' genotypes (per group). The eigenpairs are those of the record kernel, so
#' BGLR fits are unchanged up to eigenvector sign.
#'
#' @param record_ids Genotype ID of each phenotype record, in record order.
#' @param kernel Genotype-level kernel with genotype IDs as row names.
#' @param groups Optional group (environment) of each record: the kernel is
#'   restricted to records of the same group.
#' @param rel_tol Eigenvalues at or below `rel_tol * max` are dropped.
#' @return `list(vectors, values)`, values in decreasing order.
#' @keywords internal
#' @noRd
gp_bayes_record_eigen <- function(record_ids, kernel, groups = NULL, rel_tol = 1e-12) {
  kernel <- as.matrix(kernel)
  ids <- as.character(record_ids)
  n <- length(ids)
  gp_bayes_record_kernel(ids[!duplicated(ids)], kernel)  # validates IDs against the kernel
  groups <- if (is.null(groups)) rep("all", n) else as.character(groups)
  blocks <- split(seq_len(n), factor(groups, levels = unique(groups)))
  parts <- lapply(blocks, function(rows) {
    g <- ids[rows]
    ug <- unique(g)
    s <- sqrt(tabulate(match(g, ug), nbins = length(ug)))
    gi <- match(ug, rownames(kernel))
    S <- kernel[gi, gi, drop = FALSE] * outer(s, s)
    e <- eigen((S + t(S)) / 2, symmetric = TRUE)
    list(rows = rows, g = match(g, ug), s = s, vectors = e$vectors, values = e$values)
  })
  top <- max(vapply(parts, function(p) max(p$values), numeric(1)))
  vecs <- list(); vals <- list()
  for (p in parts) {
    keep <- p$values > rel_tol * top
    if (!any(keep)) next
    Y <- p$vectors[, keep, drop = FALSE] / p$s
    V <- matrix(0, n, sum(keep))
    V[p$rows, ] <- Y[p$g, , drop = FALSE]
    vecs[[length(vecs) + 1L]] <- V
    vals[[length(vals) + 1L]] <- p$values[keep]
  }
  V <- do.call(cbind, vecs)
  d <- unlist(vals, use.names = FALSE)
  ord <- order(d, decreasing = TRUE)
  list(vectors = V[, ord, drop = FALSE], values = d[ord])
}

#' BRR design `X = V D^1/2` from an eigen-decomposition
#'
#' Internal helper; same retention rule and row names as
#' `gp_bayes_brr_design_from_kernel()` (eigenvalues above `tol * max`).
#' @keywords internal
#' @noRd
gp_bayes_brr_design_from_eigen <- function(eig, tol = 1e-8) {
  d <- eig$values
  keep <- d > max(d) * tol
  X <- eig$vectors[, keep, drop = FALSE] * rep(sqrt(d[keep]), each = nrow(eig$vectors))
  rownames(X) <- as.character(seq_len(nrow(X)))
  X
}

#' Compile ETA for Bayesian Genomic Prediction Models
#'
#' This function compiles the ETA components for Bayesian genomic prediction models,
#' including RKHS, and GBLUP_BRR, based on the provided fixed and random model terms,
#' genomic or other omics relationship matrices, and phenotypic data.
#'
#' @param fixed A formula or a list of formulas specifying the fixed effects.
#' @param random A formula or a list of formulas specifying the random effects.
#' @param GS_model A character string specifying the genomic selection model to use.
#' @param fixed_term_model_bayesian Character string specifying the model for fixed terms. Default is NULL.
#' @param rand_term_model_bayesian Character vector specifying the models for each random term. Default is NULL.
#' @param pheno_data A data.frame containing the phenotypic data.
#' @param gmatrix A numeric matrix representing the genomic relationship matrix.
#' @param gkernel A numeric matrix representing a genomic kernel for RKHS models. Default is NULL.
#' @param omic1_kernel A numeric matrix representing an omics-based kernel. Default is NULL.
#' @param omic2_kernel Same as `omic1_kernel`. Default is NULL.
#' @param omic3_kernel Same as `omic1_kernel`. Default is NULL.
#' @param kernel_list Optional named list of additional relationship or kernel matrices.
#' @param gen_name A character string specifying the column name in `pheno_data` that contains the genotype identifiers.
#' @param heter_groups A character string specifying the column name in `pheno_data` for heterogeneous groups. Default is NULL.
#' @param ... Reserved for future extensions; currently ignored.
#' @return A list containing the compiled ETA components, the modified phenotypic data, and names of ETA elements.
#' @examples
#' \dontrun{
#' # Assuming pheno_data is your phenotypic dataset, gmatrix is the genomic relationship matrix:
#' result <- ETA_compiler_bayes_GBLUP(fixed = ~ fixed_effect,
#'                                    random = ~ random_effect,
#'                                    GS_model = "RKHS",
#'                                    pheno_data = pheno_data,
#'                                    gmatrix = gmatrix,
#'                                    gen_name = "GenotypeID")
#' }
#' @export

ETA_compiler_bayes_GBLUP <- function(
    fixed = NULL,
    random = NULL,
    GS_model = NULL,
    fixed_term_model_bayesian = NULL,
    rand_term_model_bayesian = NULL,
    pheno_data = NULL,
    gmatrix= NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    kernel_list = NULL,
    gen_name = NULL,
    heter_groups = NULL,
    ...
){
#browser()
  ### Create empty list for ETA

  #rm(ZE, ZEZE, Zg, K1, K2, ETA)

  non_gen_inter_test <-  NULL

  inter_gen_pos_mod <-  NULL

  non_gen_pos_mod <-  NULL

  ETA <- list()

  msg <- ""
  heter_control <- gp_normalize_single_environment_heter_controls(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    heter_resid = FALSE,
    var_cov_str = NULL
  )
  heter_groups <- heter_control$heter_groups
  fixed <- gp_bayes_normalize_fixed_argument(
    fixed = fixed,
    fixed_term_model_bayesian = fixed_term_model_bayesian
  )
  ### Get the random terms. Both no interaction and interaction terms if present in the random terms
  rand_terms <- random_terms(random = random,
                             pheno_data = pheno_data)

  ### Assign
  rand_model <- random_term_model(rand_terms = rand_terms,
                                  rand_terms_model_bayesian = rand_term_model_bayesian,
                                  GS_model = GS_model,
                                  gen_name = gen_name)

  if(!is.null(fixed)){
    ETA <-  ETA_compiler_fixed_term(fixed = fixed,
                                    pheno_data = pheno_data,
                                    fixed_term_model_bayesian = fixed_term_model_bayesian)


  }
  ################################################
  #### Check for Interaction and and non-interaction term
  ## No interaction term
  rand_terms_no_inter <- rand_terms[!grepl(":", rand_terms)]

  ## Interaction term
  check_rand_inter <- rand_terms[grepl(":", rand_terms)]

  ## Start with the No interaction terms
  if(length(rand_terms_no_inter)!=0){
    ## Check if gen_name is present and store the position
    gen_pos_mod=  match(gen_name, rand_terms)
    if(length(gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect is missing"))), call. = FALSE)}
    if(length(gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
    ## Extract other terms from the rand_terms_no_inter  expect the gen_name
    rand_terms_no_inter_no_gen <- rand_terms_no_inter[!rand_terms_no_inter%in%gen_name]
    ## Get the position of other terms (No interaction) in the random that is not gen_name
    non_gen_pos_mod <- match(rand_terms_no_inter_no_gen, rand_terms)

  } else {
    non_gen_pos_mod <-  NULL
  }

  ######## Initialize step For model adjustment for the random term with interaction

  if(length(check_rand_inter)!=0){

    test_present_of_geno <-  grep(gen_name, check_rand_inter, value = TRUE)

    if(length(test_present_of_geno)!=0 | !is.na(test_present_of_geno)){

      inter_gen_pos_mod <-  match(test_present_of_geno, rand_terms)

      non_gen_inter_test <- check_rand_inter[!check_rand_inter%in%test_present_of_geno]

    } else {

      non_gen_inter_test <-  NULL

      inter_gen_pos_mod <-  NULL
    }

    if(!is.null(non_gen_inter_test)){

      inter_non_gen_pos_mod <- match(non_gen_inter_test, rand_terms)

    }

  }

  if (length(inter_gen_pos_mod)>=1){
    if (length(pheno_data[[gen_name]]) ==length(unique(pheno_data[[gen_name]]))){
      stop(print(paste(msg, "Phenotypic data contain single environment but you specify multi-environment analysis.")), call. = FALSE)
    }

  }

  if (identical(GS_model, "RKHS") &&
      is.null(rand_term_model_bayesian) &&
      length(inter_gen_pos_mod) >= 1 &&
      length(gen_pos_mod) == 1 &&
      identical(rand_model[gen_pos_mod], "RKHS")) {
    rand_model[gen_pos_mod] <- "BRR"
    rand_model[inter_gen_pos_mod] <- "RKHS"
  }
  #########
  #### When the genotype are present in more than one environment/location
  if(length(pheno_data[[gen_name]])>length(unique(pheno_data[[gen_name]]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[[gen_name]])-1)

    if(gp_bayes_has_multiple_residual_groups(pheno_data, heter_groups)){
      ZE <- model.matrix(~factor(pheno_data[[heter_groups]])-1)
      ZEZE <-tcrossprod(ZE)

    }

  } else {

    Zg <-  NULL
    ZE <-  NULL
    ZEZE <- NULL
  }

  rand_model_copy <- rand_model
  #####
  #####################################################
  for (ra in 1:length(rand_terms)){

    rand_model <- rand_model_copy[ra]

    ### start when no genotype
    if(is.null(non_gen_pos_mod)){
      if(length(non_gen_pos_mod)!=0){
        if((ra == non_gen_pos_mod & !is.null(ZE))){
          ETA[[length(ETA) + 1]] <- list(X= ZE,
                                         model=rand_model,
                                         saveEffects=TRUE)

        }

      }

    } ### End

    datasets <- gp_collect_kernel_inputs(
      gmatrix = gmatrix,
      gkernel = gkernel,
      omic1_kernel = omic1_kernel,
      omic2_kernel = omic2_kernel,
      omic3_kernel = omic3_kernel,
      kernel_list = kernel_list
    )
    dataset_names <- names(datasets)
    ETA_element_name <- character()

    for (i in seq_along(datasets)) {
      dataset <- datasets[[i]]
      if (!is.null(dataset)) { ## this seems redundant but useful
        if(rand_model== "RKHS"){
          if((ra == gen_pos_mod ) & is.null(Zg)){
            ETA[[length(ETA) + 1]] <- list(K= as.matrix(dataset),
                                           model = rand_model,
                                           saveEffects = TRUE)
          } else if((ra == gen_pos_mod) & !is.null(Zg)){
            K1 <- gp_bayes_record_kernel(pheno_data[[gen_name]], dataset)
            # V/d: BGLR skips its records x records eigen-decomposition; K stays
            # for the variance-component and output steps.
            eig <- gp_bayes_record_eigen(pheno_data[[gen_name]], dataset)
            ETA[[length(ETA) + 1]] <- list(K= K1, V = eig$vectors, d = eig$values,
                                           model = rand_model,
                                           saveEffects = TRUE)
          } else {

            if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
              if((ra == inter_gen_pos_mod & !is.null(Zg)) & !is.null(ZEZE)){

                K1 <- gp_bayes_record_kernel(pheno_data[[gen_name]], dataset)

                K2<-K1*ZEZE
                eig <- gp_bayes_record_eigen(pheno_data[[gen_name]], dataset,
                                             groups = pheno_data[[heter_groups]])

                ETA[[length(ETA) + 1]] <- list(K= K2, V = eig$vectors, d = eig$values,
                                               model=rand_model,
                                               saveEffects=TRUE)

              }

            }
          }
          ETA_element_name <- c(ETA_element_name, dataset_names[i])
        } else{

          if(rand_model== "BRR"){
            if((ra == gen_pos_mod ) & is.null(Zg)){
              ETA[[length(ETA) + 1]] <- list(X= gp_bayes_brr_design_from_kernel(as.matrix(dataset)),
                                             model = rand_model,
                                             saveEffects = TRUE)
            } else if((ra == gen_pos_mod) & !is.null(Zg)){
              eig <- gp_bayes_record_eigen(pheno_data[[gen_name]], dataset)
              ETA[[length(ETA) + 1]] <- list(X= gp_bayes_brr_design_from_eigen(eig),
                                             model = rand_model,
                                             saveEffects = TRUE)
            } else {

              if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
                if((ra == inter_gen_pos_mod & !is.null(Zg)) & !is.null(ZEZE)){

                  eig <- gp_bayes_record_eigen(pheno_data[[gen_name]], dataset,
                                               groups = pheno_data[[heter_groups]])

                  ETA[[length(ETA) + 1]] <- list(X= gp_bayes_brr_design_from_eigen(eig),
                                                 model=rand_model,
                                                 saveEffects=TRUE)

                }

              }
            }
            ETA_element_name <- c(ETA_element_name, dataset_names[i])
          }

        }
      }
    }
  }
  output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)

  return(output)

}


# ETA_compiler_bayes_GBLUP <- function(
#     fixed = NULL,
#     random = NULL,
#     GS_model = NULL,
#     fixed_term_model_bayesian = NULL,
#     rand_term_model_bayesian = NULL,
#     pheno_data = NULL,
#     gmatrix= NULL,
#     gkernel = NULL,
#     omic1_kernel = NULL,
#     omic2_kernel = NULL,
#     omic3_kernel = NULL,
#     gen_name = NULL,
#     heter_groups = NULL,
#     ...
# ){
#
#   ### Create empty list for ETA
#
#   #rm(ZE, ZEZE, Zg, K1, K2, ETA)
#   ETA = list()
#
#   msg <- ""
#   ### Get the random terms. Both no interaction and interaction terms if present in the random terms
#   rand_terms <- random_terms(random = random,
#                              object = pheno_data)
#
#   ### Assign
#   rand_model <- random_term_model(rand_terms = rand_terms,
#                                   rand_terms_model_bayesian = rand_term_model_bayesian,
#                                   GS_model = GS_model,
#                                   gen_name = gen_name)
#
#   if(!is.null(fixed)){
#     ETA <-  ETA_compiler_fixed_term(fixed = fixed,
#                                     pheno_data = pheno_data,
#                                     fixed_term_model_bayesian = fixed_term_model_bayesian)
#
#
#   }
#  ################################################
#   #### Check for Interaction and and non-interaction term
#   ## No interaction term
#   rand_terms_no_inter = rand_terms[!grepl(":", rand_terms)]
#
#   ## Interaction term
#   check_rand_inter = rand_terms[grepl(":", rand_terms)]
#
#   ## Start with the No interaction terms
#   if(length(rand_terms_no_inter)!=0){
#     ## Check if gen_name is present and store the position
#     gen_pos_mod=  match(gen_name, rand_terms)
#     if(length(gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect is missing"))), call. = FALSE)}
#     if(length(gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
#     ## Extract other terms from the rand_terms_no_inter  expect the gen_name
#     rand_terms_no_inter_no_gen = rand_terms_no_inter[!rand_terms_no_inter%in%gen_name]
#     ## Get the position of other terms (No interaction) in the random that is not gen_name
#     non_gen_pos_mod = match(rand_terms_no_inter_no_gen, rand_terms)
#
#   }
#
#   ######## Initialize step For model adjustment for the random term with interaction
#
#   if(length(check_rand_inter)!=0){
#
#     test_present_of_geno = grep(gen_name, check_rand_inter, value = TRUE)
#
#     if(length(test_present_of_geno)!=0 | !is.na(test_present_of_geno)){
#
#       inter_gen_pos_mod = match(test_present_of_geno, rand_terms)
#
#       non_gen_inter_test =  check_rand_inter[!check_rand_inter%in%test_present_of_geno]
#
#     } else {
#
#       non_gen_inter_test = NULL
#
#       inter_gen_pos_mod = NULL
#     }
#
#     if(!is.null(non_gen_inter_test)){
#
#       inter_non_gen_pos_mod = match(non_gen_inter_test, rand_terms)
#
#     }
#
#   }
#
#   #########
#   #### When the genotype are present in more than one environment/location
#   if(length(pheno_data[[gen_name]])>length(unique(pheno_data[[gen_name]]))){
#     ### incidence matrix for main eff. of the genotypes
#     Zg<-stats::model.matrix(~factor(pheno_data[[gen_name]])-1)
#
#     if(!is.null(heter_groups)){
#       ZE <- model.matrix(~factor(pheno_data[[heter_groups]])-1)
#       ZEZE<-tcrossprod(ZE)
#
#     }
#
#   }
#
#   rand_model_copy = rand_model
#   #####
#   #####################################################
#   for (ra in 1:length(rand_terms)){
#
#     rand_model = rand_model_copy[ra]
#
#
#     ### start when no genotype
#     if(exists("non_gen_pos_mod")){
#       if(length(non_gen_pos_mod)!=0){
#         if((ra == non_gen_pos_mod & exists("ZE"))){
#           ETA[[length(ETA) + 1]] <- list(X= ZE,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#         }
#
#       }
#
#     } ### End
#
#     if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod ) & !exists("Zg")){
#           ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("gmatrix")
#           } else if((ra == gen_pos_mod) & exists("Zg")){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix")
#             } else if(((ra == gen_pos_mod & exists("Zg")))){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix")
#                 }
#
#               }
#
#
#
#               }
#
#
#
#             }
#
#
#           }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel")
#                 }
#
#
#
#                 }
#
#               }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel")
#                 }
#
#                 if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                   if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel")
#                   }
#
#                 }
#
#               }
#             }
#
#             ETA_element_name = c("gkernel")
#           }
#
#
#         }
#
#       }
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(omic1_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                   K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic1_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#       ETA_element_name = c("omic1_kernel")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(omic2_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel")
#           } else if(((ra == gen_pos_mod & exists("Zg")))){
#
#             K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic2_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#
#       ETA_element_name = c("omic2_kernel")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(omic3_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic3_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic3_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic3_kernel")
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic3_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#
#       ETA_element_name = c("omic3_kernel")
#       #If user provide geno_data and omic1_data
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#       #ETA_element_name = c("geno_data", "omic2_data")
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic2_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic2_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic3_data")
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic3_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic3_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic3_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic3_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic3_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic3_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       } ## End
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic2_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if(((ra == gen_pos_mod & !exists("Zg")))){
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#       #ETA_element_name = c("omic1_data", "omic2_data")
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#
#       #ETA_element_name = c("omic1_data", "omic3_data")
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic2_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic2_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#
#
#
#       #ETA_element_name = c("omic2_data", "omic3_data")
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic1_data", "omic2_data")
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(rand_model== "RKHS"){
#
#         if(!is.null(gkernel)){
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#         } else {
#
#           if(!is.null(gmatrix)){
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#           }
#
#           ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#         }
#
#       } else {
#
#         if(rand_model== "BRR"){
#           if(!is.null(gkernel)){
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#           } else {
#
#             if(!is.null(gmatrix)){
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#             }
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#           }
#         }
#
#       }
#       #ETA_element_name = c("geno_data", "omic1_data", "omic3_data")
#       ###
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic2_data", "omic3_data")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_2 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           K1_3 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_2,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_3,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2 <- K1*ZEZE
#
#               K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2_1<- K1_1*ZEZE
#
#               K1_2 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_2 <- K1_2*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_2 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_2,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_3,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2 <- K1*ZEZE
#
#                 K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1<- K1_1*ZEZE
#
#                 K1_2 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_2 <- K1_2*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#         }
#
#
#       }
#
#
#     } else {
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_2omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_1omic<-K1_omic*ZEZE
#                 ###
#                 K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1omic<-K1_2omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_2omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K1_1omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_1omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_2omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_2omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K1_1omic<-K1_omic*ZEZE
#                   ###
#                   K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_1omic<-K1_2omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_2omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_2omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_2omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K1_1omic<-K1_omic*ZEZE
#                   ###
#                   K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_1omic<-K1_2omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_2omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K1_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_2omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_2omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K1_1omic<-K1_omic*ZEZE
#                     ###
#                     K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_1omic<-K1_2omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_2omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K1_1omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_1omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_2omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#     }
#
#   }
#
#   output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)
#
#   names(output) <- c("ETA", "pheno_data", "ETA_element_name")
#
#   return(output)
#
#
#
# }
#
#
