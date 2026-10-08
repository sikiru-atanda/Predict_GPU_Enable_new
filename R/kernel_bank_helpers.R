gp_is_kernel_matrix <- function(x) {
  is.matrix(x) && is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x))
}

gp_is_kernel_list <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    length(x) > 0L &&
    all(vapply(x, gp_is_kernel_matrix, logical(1)))
}

gp_sanitize_kernel_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[nchar(x) == 0L] <- "kernel"
  make.names(x, unique = FALSE)
}

gp_named_method_list <- function(values, methods) {
  methods <- gp_sanitize_kernel_name(methods)
  names(values) <- make.unique(methods, sep = "_")
  values
}

gp_kernel_supported_methods <- function(include_aliases = FALSE) {
  methods <- c(
    "Gaussian_kernel",
    "Linear_kernel",
    "Composite_kernel",
    "Poly2_kernel",
    "Poly3_kernel",
    "Poly4_kernel",
    "Matern_kernel",
    "Matern12_kernel",
    "Matern32_kernel",
    "Matern52_kernel",
    "Laplacian_kernel",
    "RationalQuadratic_kernel"
  )
  if (!isTRUE(include_aliases)) {
    return(methods)
  }
  c(
    methods,
    "Gaussian",
    "Linear",
    "Composite",
    "Poly2",
    "Poly3",
    "Poly4",
    "Matern",
    "Matern12",
    "Matern32",
    "Matern52",
    "Laplacian",
    "RationalQuadratic"
  )
}

gp_normalize_kernel_methods <- function(method) {
  if (is.null(method)) {
    return(NULL)
  }
  method <- as.character(method)
  aliases <- c(
    Gaussian = "Gaussian_kernel",
    Linear = "Linear_kernel",
    Composite = "Composite_kernel",
    Poly2 = "Poly2_kernel",
    Poly3 = "Poly3_kernel",
    Poly4 = "Poly4_kernel",
    Matern = "Matern_kernel",
    Matern12 = "Matern12_kernel",
    Matern32 = "Matern32_kernel",
    Matern52 = "Matern52_kernel",
    Laplacian = "Laplacian_kernel",
    RationalQuadratic = "RationalQuadratic_kernel"
  )
  idx <- method %in% names(aliases)
  method[idx] <- unname(aliases[method[idx]])
  method
}

gp_kernel_source_name <- function(x, fallback) {
  source <- attr(x, "predictpror_source_kernel", exact = TRUE)
  if (!is.null(source) && length(source) > 0L && nzchar(as.character(source[[1L]]))) {
    return(as.character(source[[1L]]))
  }
  fallback
}

gp_kernel_key_source_name <- function(kernel_key) {
  canonical <- c(
    gmatrix_model_ready = "gmatrix",
    gkernel_model_ready = "gkernel",
    omic1_kernel_model_ready = "omic1_kernel",
    omic2_kernel_model_ready = "omic2_kernel",
    omic3_kernel_model_ready = "omic3_kernel"
  )
  mapped <- unname(canonical[as.character(kernel_key)[[1L]]])
  if (length(mapped) && !is.na(mapped) && nzchar(mapped)) {
    return(mapped)
  }
  as.character(kernel_key)[[1L]]
}

# Reports the kernels a multi-trait fit actually used. Every multi-trait route
# -- ASReml, GP or Bayesian; single-environment or MET; cross-validation or true
# prediction -- returns this table so a multi-kernel run is never reported as a
# single-kernel run.
gp_multitrait_kernel_configuration <- function(kernels,
                                               role = "independent_genetic_covariance_term") {
  kernel_names <- if (is.character(kernels)) {
    kernels
  } else {
    names(kernels) %||% character()
  }
  data.frame(
    Kernel = as.character(kernel_names),
    Kernel_role = rep(as.character(role)[[1L]], length(kernel_names)),
    stringsAsFactors = FALSE
  )
}

gp_multi_kernel_parameter_rows <- function(kernel_names = NULL,
                                           kernel_count = NULL,
                                           strategy = NA_character_) {
  kernel_names <- as.character(kernel_names %||% character())
  kernel_names <- kernel_names[!is.na(kernel_names) & nzchar(kernel_names)]
  if (is.null(kernel_count)) {
    kernel_count <- length(kernel_names)
  }
  kernel_count <- as.integer(kernel_count[[1L]])
  if (!is.finite(kernel_count) || kernel_count < 0L) {
    kernel_count <- length(kernel_names)
  }
  strategy <- as.character(strategy %||% NA_character_)[[1L]]
  data.frame(
    stat = c("multi_kernel_count", "multi_kernel_names", "multi_kernel_strategy"),
    summary = c(
      as.character(kernel_count),
      if (length(kernel_names)) paste(kernel_names, collapse = ";") else NA_character_,
      strategy
    ),
    stringsAsFactors = FALSE
  )
}

gp_kernel_extra_key <- function(base_key, method) {
  method <- gp_sanitize_kernel_name(method)
  if (identical(base_key, "gmatrix_model_ready")) {
    return(paste("gmatrix", method, "model_ready", sep = "_"))
  }
  if (grepl("_kernel_model_ready$", base_key)) {
    prefix <- sub("_kernel_model_ready$", "", base_key)
    return(paste(prefix, method, "model_ready", sep = "_"))
  }
  paste(base_key, method, sep = "_")
}

gp_add_kernel_inputs <- function(kernel_inputs, kernel_key, kernels) {
  if (is.null(kernels)) {
    return(kernel_inputs)
  }
  if (gp_is_kernel_list(kernels)) {
    kernel_names <- names(kernels)
    if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
      kernel_names <- paste0("kernel", seq_along(kernels))
    }
    kernel_inputs[[kernel_key]] <- kernels[[1]]
    attr(kernel_inputs[[kernel_key]], "predictpror_source_kernel") <- kernel_names[[1]]
    if (length(kernels) > 1L) {
      # Reserve every extra key before assignment. Distinct source labels can
      # sanitize to the same key, including labels with an existing suffix.
      extra_keys <- vapply(kernel_names[-1L], function(nm) {
        gp_kernel_extra_key(kernel_key, nm)
      }, character(1L))
      extra_keys <- make.unique(c(names(kernel_inputs), extra_keys), sep = "_")[
        length(kernel_inputs) + seq_along(extra_keys)
      ]
      for (i in 2:length(kernels)) {
        extra_key <- extra_keys[[i - 1L]]
        kernel_inputs[[extra_key]] <- kernels[[i]]
        attr(kernel_inputs[[extra_key]], "predictpror_source_kernel") <- kernel_names[[i]]
      }
    }
    return(kernel_inputs)
  }
  kernel_inputs[[kernel_key]] <- kernels
  if (is.null(attr(
    kernel_inputs[[kernel_key]],
    "predictpror_source_kernel",
    exact = TRUE
  ))) {
    attr(
      kernel_inputs[[kernel_key]],
      "predictpror_source_kernel"
    ) <- gp_kernel_key_source_name(kernel_key)
  }
  kernel_inputs
}

gp_match_kernel_bank <- function(kernels,
                                 object_pheno,
                                 gen_name,
                                 test_set = NULL,
                                 train_set = NULL,
                                 heter_groups = NULL,
                                 low_call_rate_inds_removed = NULL,
                                 test_set_overrides_inferred = FALSE,
                                 message = TRUE) {
  if (!gp_is_kernel_list(kernels)) {
    return(pheno_geno_match(
      object_pheno = object_pheno,
      object_geno = kernels,
      gen_name = gen_name,
      test_set = test_set,
      train_set = train_set,
      heter_groups = heter_groups,
      low_call_rate_inds_removed = low_call_rate_inds_removed,
      test_set_overrides_inferred = test_set_overrides_inferred,
      message = message
    ))
  }

  matched <- vector("list", length(kernels))
  kernel_names <- names(kernels)
  if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
    kernel_names <- paste0("kernel", seq_along(kernels))
  }
  names(matched) <- kernel_names
  next_test_set <- test_set
  for (i in seq_along(kernels)) {
    source_name <- gp_kernel_source_name(kernels[[i]], kernel_names[[i]])
    match_res <- pheno_geno_match(
      object_pheno = object_pheno,
      object_geno = kernels[[i]],
      gen_name = gen_name,
      test_set = next_test_set,
      train_set = train_set,
      heter_groups = heter_groups,
      low_call_rate_inds_removed = low_call_rate_inds_removed,
      test_set_overrides_inferred = test_set_overrides_inferred,
      message = message
    )
    matched[[i]] <- match_res[["geno_pheno_match_data"]]
    attr(matched[[i]], "predictpror_source_kernel") <- source_name
    if (length(match_res) > 1L) {
      next_test_set <- match_res[["test_data"]]
    }
  }
  out <- list(geno_pheno_match_data = matched)
  if (!is.null(next_test_set)) {
    out[["test_data"]] <- next_test_set
  }
  out
}

gp_collect_kernel_inputs <- function(gmatrix = NULL,
                                     gkernel = NULL,
                                     omic1_kernel = NULL,
                                     omic2_kernel = NULL,
                                     omic3_kernel = NULL,
                                     kernel_list = NULL) {
  kernels <- list()
  append_kernel <- function(value, fallback, prefix = NULL) {
    if (is.null(value)) {
      return()
    }
    if (gp_is_kernel_list(value)) {
      value_names <- names(value)
      if (is.null(value_names) || any(!nzchar(value_names))) {
        value_names <- paste0("kernel", seq_along(value))
      }
      for (i in seq_along(value)) {
        nm <- gp_kernel_source_name(value[[i]], value_names[[i]])
        if (!is.null(prefix)) {
          nm <- paste(prefix, nm, sep = "_")
        }
        kernels <<- c(kernels, stats::setNames(list(value[[i]]), nm))
      }
      return()
    }
    nm <- gp_kernel_source_name(value, fallback)
    if (!is.null(prefix)) {
      nm <- paste(prefix, nm, sep = "_")
    }
    # Append by position; source names become unique after all inputs arrive.
    # Assigning by name here would overwrite an earlier kernel with that tag.
    kernels <<- c(kernels, stats::setNames(list(value), nm))
  }
  # gmatrix and gkernel are independent user inputs.  Retain both when both
  # are supplied; choosing one by precedence silently discards a covariance
  # source and violates the multi-kernel contract.
  append_kernel(gmatrix, "gmatrix")
  append_kernel(gkernel, "gkernel")
  append_kernel(omic1_kernel, "kernel", "omic1")
  append_kernel(omic2_kernel, "kernel", "omic2")
  append_kernel(omic3_kernel, "kernel", "omic3")
  kernels <- kernels[!vapply(kernels, is.null, logical(1))]

  if (gp_is_kernel_list(kernel_list)) {
    canonical <- c(
      "gmatrix_model_ready",
      "omic1_kernel_model_ready",
      "omic2_kernel_model_ready",
      "omic3_kernel_model_ready"
    )
    input_names <- names(kernel_list)
    if (is.null(input_names)) {
      input_names <- paste0("kernel", seq_along(kernel_list))
    }
    unnamed <- is.na(input_names) | !nzchar(input_names)
    input_names[unnamed] <- paste0("kernel", which(unnamed))
    keep <- !input_names %in% canonical
    # Subset by position so repeated list names retain distinct matrices.
    extras <- stats::setNames(kernel_list[keep], input_names[keep])
    extras <- extras[vapply(extras, gp_is_kernel_matrix, logical(1))]
    if (length(extras) > 0L) {
      names(extras) <- vapply(
        seq_along(extras),
        function(i) gp_kernel_source_name(extras[[i]], names(extras)[[i]]),
        character(1)
      )
      kernels <- c(kernels, extras)
    }
  }

  if (length(kernels) > 0L) {
    names(kernels) <- make.unique(gp_sanitize_kernel_name(names(kernels)), sep = "_")
  }
  kernels
}

gp_drop_test_kernel_bank <- function(kernel_list, test_ids, label = "kernel_list") {
  if (is.null(kernel_list) || is.null(test_ids)) {
    return(kernel_list)
  }
  if (!gp_is_kernel_list(kernel_list)) {
    rn <- rownames(kernel_list)
    cn <- colnames(kernel_list)
    if (is.null(rn) || is.null(cn)) {
      stop(paste(label, "must have row and column names before Bayesian CV test-set filtering."), call. = FALSE)
    }
    return(kernel_list[!rn %in% test_ids, !cn %in% test_ids, drop = FALSE])
  }
  lapply(kernel_list, function(x) {
    rn <- rownames(x)
    cn <- colnames(x)
    if (is.null(rn) || is.null(cn)) {
      stop(paste(label, "entries must have row and column names before Bayesian CV test-set filtering."), call. = FALSE)
    }
    x[!rn %in% test_ids, !cn %in% test_ids, drop = FALSE]
  })
}
