
#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param heter_groups
#' @param CV
#' @param random_state
#' @param replication
#' @param nfolds
#'
#' @return
#' @export
#'
#' @examples


# CV0_CV1_CV2_for_multi_environment <- function(
#     pheno_data,
#     gen_name,            # column with genotype IDs
#     heter_groups,        # column with environments (e.g., site-year)
#     CV = c(0, 1, 2),     # 0 = LOEO, 1 = leave genotypes out, 2 = sparse within genotype
#     nfolds = NULL,       # CV0: defaults to #envs; CV1/CV2: typical 5
#     random_state = NULL,
#     replication = 1,
#     message = TRUE,
#     ...
# ) {
#   msg_bar <- "\n==================================================\n"
#   stopifnot(is.data.frame(pheno_data))
#   if (!gen_name %in% names(pheno_data)) stop("`gen_name` column not found.")
#   if (!heter_groups %in% names(pheno_data)) stop("`heter_groups` column not found.")
#   CV <- match.arg(as.character(CV), choices = c("0","1","2"))
#   CV <- as.integer(CV)
#
#   g_full <- pheno_data[[gen_name]]
#   e_full <- pheno_data[[heter_groups]]
#   N      <- nrow(pheno_data)
#
#   keep <- !(is.na(g_full) | is.na(e_full))
#   if (any(!keep) && isTRUE(message)) {
#     warning(paste0(msg_bar, "Removed ", sum(!keep),
#                    " rows with NA in key columns *for fold assignment only*; ",
#                    "these rows will get NA folds."), call. = FALSE)
#   }
#
#   g <- g_full[keep]
#   e <- e_full[keep]
#   n <- length(g)
#
#   uniq_g <- unique(g)
#   uniq_e <- unique(e)
#   G <- length(uniq_g)
#   E <- length(uniq_e)
#
#   if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state)
#
#   # Defaults for nfolds
#   if (CV == 0) {
#     if (is.null(nfolds)) nfolds <- E
#     if (nfolds < 2) stop("For CV0, `nfolds` must be >= 2 (usually E for LOEO).")
#   } else {
#     if (is.null(nfolds)) nfolds <- 5L
#     if (nfolds < 2) stop("`nfolds` must be >= 2.")
#   }
#
#   assign_balanced <- function(ids, k) {
#     ids <- sample(ids)
#     rep(seq_len(k), length.out = length(ids))
#   }
#
#   out <- vector("list", replication)
#   names(out) <- sprintf("Rep%02d_%dFold_CV%d", seq_len(replication), nfolds, CV)
#
#   for (r in seq_len(replication)) {
#     if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state + r - 1L)
#
#     folds_kept <- integer(n)
#
#     if (CV == 0) {
#       if (nfolds == E) {
#         env_order <- sample(uniq_e)
#         env2fold  <- setNames(seq_len(E), env_order)
#       } else {
#         env_order <- sample(uniq_e)
#         env2fold  <- setNames(assign_balanced(env_order, nfolds), env_order)
#       }
#       folds_kept <- unname(env2fold[as.character(e)])
#     }
#
#     if (CV == 1) {
#       geno_order <- sample(uniq_g)
#       geno2fold  <- setNames(assign_balanced(geno_order, nfolds), geno_order)
#       folds_kept <- unname(geno2fold[as.character(g)])
#     }
#
#     if (CV == 2) {
#       folds_kept <- integer(n)
#       idx_by_g <- split(seq_len(n), g)
#       for (ix in idx_by_g) {
#         k <- length(ix)
#         if (!k) next
#         sh <- ix[sample.int(k)]
#         folds_kept[sh] <- rep(seq_len(nfolds), length.out = k)
#       }
#     }
#
#     # Expand back to full length (align to pheno_data rows)
#     folds_all <- rep(NA_integer_, N)
#     folds_all[keep] <- as.integer(folds_kept)
#     # record effective nfolds in attribute (runner will read it)
#     attr(folds_all, "nfolds_eff") <- max(folds_kept, na.rm = TRUE)
#     out[[r]] <- folds_all
#   }
#
#   if (isTRUE(message)) {
#     if (CV == 0 && nfolds == E) {
#       message(paste0(msg_bar, "CV0 (LOEO): Each environment serves once as test fold (classic leave-one-environment-out)."))
#     } else if (CV == 0) {
#       message(paste0(msg_bar, "CV0 (grouped envs): Environments grouped into ", nfolds, " folds (balanced)."))
#     } else if (CV == 1) {
#       message(paste0(msg_bar, "CV1: Genotype-level K-fold (predicting entirely new genotypes)."))
#     } else if (CV == 2) {
#       message(paste0(msg_bar, "CV2: Sparse within-genotype K-fold (predicting missing G×E cells with data in other envs)."))
#     }
#   }
#
#   out
# }



CV0_CV1_CV2_for_multi_environment <- function(
    pheno_data,
    gen_name,            # column with genotype IDs
    heter_groups,        # column with environments (e.g., site-year)
    CV = c(0, 1, 2),     # 0 = LOEO, 1 = leave genotypes out, 2 = sparse within genotype
    nfolds = NULL,       # CV0: defaults to #envs; CV1/CV2: typical 5
    random_state = NULL,
    replication = 1,
    message = TRUE,
    ...
) {
  msg_bar <- "\n==================================================\n"
  verbose <- isTRUE(message)

  stopifnot(is.data.frame(pheno_data))
  if (!gen_name %in% names(pheno_data)) stop("`gen_name` column not found.")
  if (!heter_groups %in% names(pheno_data)) stop("`heter_groups` column not found.")
  CV <- match.arg(as.character(CV), choices = c("0","1","2"))
  CV <- as.integer(CV)

  g <- pheno_data[[gen_name]]
  e <- pheno_data[[heter_groups]]

  # Drop NA rows in keys
  keep <- !(is.na(g) | is.na(e))
  if (any(!keep)) {
    if (verbose) base::warning(paste0(msg_bar, "Removed ", sum(!keep), " rows with NA in key columns."), call. = FALSE)
    g <- g[keep]; e <- e[keep]
  }

  n <- length(g)
  uniq_g <- unique(g)
  uniq_e <- unique(e)
  G <- length(uniq_g)
  E <- length(uniq_e)

  if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state)

  # Defaults for nfolds (and clamp to available groups to avoid empty folds)
  if (CV == 0) {
    if (is.null(nfolds)) nfolds <- E
    if (nfolds < 2) stop("For CV0, `nfolds` must be >= 2 (usually E for LOEO).")
    if (nfolds > E) {
      if (verbose) base::message(paste0(msg_bar, "CV0: nfolds > #environments; clamping to E = ", E, "."))
      nfolds <- E
    }
  } else {
    if (is.null(nfolds)) nfolds <- 5L
    if (nfolds < 2) stop("`nfolds` must be >= 2.")
    if (CV == 1 && nfolds > G) {
      if (verbose) base::message(paste0(msg_bar, "CV1: nfolds > #genotypes; clamping to G = ", G, "."))
      nfolds <- G
    }
  }

  # Utility: balanced assignment helper (for group labels, not rows)
  assign_balanced <- function(ids, k) {
    ids <- sample(ids)                             # shuffle groups
    setNames(rep(seq_len(k), length.out = length(ids)), ids)
  }

  out <- vector("list", replication)
  names(out) <- sprintf("Rep%02d_%dFold_CV%d", seq_len(replication), nfolds, CV)

  for (r in seq_len(replication)) {
    if (!is.null(random_state) && is.numeric(random_state)) set.seed(random_state + r - 1L)

    if (CV == 0) {
      # ---- CV0: Leave-Environment-Out (or grouped envs) ----
      if (nfolds == E) {
        env_order <- sample(uniq_e)
        env2fold  <- setNames(seq_len(E), env_order)          # one env per fold
      } else {
        env2fold  <- assign_balanced(uniq_e, nfolds)          # group envs into k folds
      }
      folds <- unname(env2fold[as.character(e)])

      # Invariant: each environment should map to exactly one fold
      chk <- tapply(folds, e, function(v) length(unique(v)))
      if (any(chk > 1, na.rm = TRUE)) stop("CV0 invariant failed: some environments map to multiple folds.")
      out[[r]] <- as.integer(folds)
    }

    if (CV == 1) {
      # ---- CV1: leave-genotypes-out (all rows of a genotype share the fold) ----
      geno2fold <- assign_balanced(uniq_g, nfolds)
      folds <- unname(geno2fold[as.character(g)])

      # Invariant: each genotype should map to exactly one fold
      chk <- tapply(folds, g, function(v) length(unique(v)))
      if (any(chk > 1, na.rm = TRUE)) stop("CV1 invariant failed: some genotypes map to multiple folds.")
      out[[r]] <- as.integer(folds)
    }

    if (CV == 2) {
      # ---- CV2: within-genotype sparse assignment across environments ----
      folds <- integer(n)
      idx_by_g <- split(seq_len(n), g)
      n_singleton <- 0L

      for (idx in idx_by_g) {
        k <- length(idx)
        if (k <= 1L) {
          # Singleton genotype: make it train-only by assigning fold = 0 (never selected as test)
          folds[idx] <- 0L
          n_singleton <- n_singleton + 1L
          next
        }
        idx_shuf <- sample(idx, k)
        kk <- min(nfolds, k)                     # can’t use more folds than rows for this genotype
        per <- sample(seq_len(kk))               # start with a permuted set of folds 1..kk
        folds[idx_shuf] <- rep(per, length.out = k)  # round-robin to spread over folds
      }

      # Optional check: warn if any multi-env genotype still ends up in a single fold (should be rare)
      split_f <- split(folds, g)
      bad <- vapply(split_f, function(fv) {
        kpos <- sum(fv > 0)
        kpos > 1L && length(unique(fv[fv > 0])) == 1L
      }, logical(1))
      if (any(bad) && verbose) {
        base::warning("CV2: some multi-environment genotypes fell into a single fold; consider larger nfolds.", call. = FALSE)
      }

      if (n_singleton > 0L && verbose) {
        base::message(sprintf("%sCV2: %d genotype(s) have only one observation; assigned fold=0 (train-only).", msg_bar, n_singleton))
      }

      out[[r]] <- as.integer(folds)
    }
  }

  # Friendly notes
  if (verbose) {
    if (CV == 0 && nfolds == E) {
      base::message(paste0(msg_bar, "CV0 (LOEO): Each environment serves once as the test fold (classic leave-one-environment-out)."))
    } else if (CV == 0) {
      base::message(paste0(msg_bar, "CV0 (grouped envs): Environments grouped into ", nfolds, " folds (balanced)."))
    } else if (CV == 1) {
      base::message(paste0(msg_bar, "CV1: Genotype-level K-fold (predicting entirely new genotypes)."))
    } else if (CV == 2) {
      base::message(paste0(msg_bar, "CV2: Sparse within-genotype K-fold (predicting missing G×E cells with data in other environments)."))
    }
  }

  out
}



# CV1_CV2_for_multi_environment <- function(
#                                           pheno_data,
#                                           gen_name,
#                                           heter_groups,
#                                           CV = NULL,
#                                           nfolds = NULL,
#                                           random_state = NULL,
#                                           replication = 1,
#                                           ...){
#
#   msg <- "\n==================================================\n"
#
#   if (is.null(heter_groups)){stop(paste(msg,"Provide the a pointer (heter_groups) to the column contaning the environments."), call. = FALSE)}
#   if(CV>2){stop(message(paste(msg,"CV must be 1 or 2")), call. = FALSE)}
#   if(is.null(nfolds)){stop(paste(msg,"Provide value the number of desired folds"), call. = FALSE)}
#   ## Order the pheno_data data by gen_name and by Environment
#   pheno_data = pheno_data[order(pheno_data[[gen_name]]), ]
#   pheno_data = pheno_data[order(pheno_data[[heter_groups]]), ]
#
#   #nEnv <- length(unique(pheno_data[[heter_groups]]))
#
#   ID_GIDs = as.character(unique(pheno_data[[gen_name]]))
#   Envs_ID_GIDs = as.character(pheno_data[[gen_name]])
#
#   if(length(ID_GIDs)==length(Envs_ID_GIDs)){stop(paste(msg,'CV1 and CV2 works when number of environment is greater than 1'), call. = FALSE)}
#
#   if(replication>1){warning(paste(msg,paste('You request for', replication, 'replications this might takes some time to run all the replications.')),
#                             call. = FALSE)}
#
#   Rep_FoldCV = vector(mode = "list", length = replication)
#
#   All_nfolds <- vector(mode = "integer",  length(pheno_data[[gen_name]]))
#
#   if(!is.null(random_state) & is.numeric(random_state)){
#     set.seed(random_state)
#   }
#
#   if (CV == 1) {
#
#     for (r in 1:replication) {
#
#       mfold <- sample(1:nfolds, size = length(ID_GIDs), replace = TRUE)
#
#       for (i in 1:length(pheno_data[[gen_name]])) {
#
#         All_nfolds[i] <- mfold[which(ID_GIDs == pheno_data[[gen_name]][i])]
#
#       }
#
#       Rep_FoldCV[[r]] <- All_nfolds
#
#       names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
#                                           paste(nfolds, 'Fold', sep = ""), sep="_"),
#                                     paste("CV", CV, sep=""), sep = "")
#
#     }
#
#   }
#
#   if (CV == 2) {
#
#     for (r in 1:replication) {
#
#       for (i in ID_GIDs) {
#
#         Env_GIDs = which(Envs_ID_GIDs == i)
#
#         Env_GIDs_size = length(Env_GIDs)
#
#         tmpFold <- sample(1:nfolds, size =  Env_GIDs_size, replace =  Env_GIDs_size > nfolds)
#
#         All_nfolds[Env_GIDs] <- tmpFold
#
#       }
#
#       Rep_FoldCV[[r]] <- All_nfolds
#
#       names(Rep_FoldCV)[r] <- paste(paste(paste("Rep", r, sep = ""),
#                                           paste(nfolds, 'Fold', sep = ""), sep="_"),
#                                     paste("CV", CV, sep=""), sep = "")
#
#     }
#
#   }
#   return(Rep_FoldCV)
# }
