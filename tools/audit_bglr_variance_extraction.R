.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(pkgload)
  library(BGLR)
})

pkgload::load_all(".")

n_iter <- 160L
burn_in <- 40L
thin <- 2L
tolerance <- 1e-8
audit_rows <- list()

add_audit <- function(scope, input, model, component, package, direct, source) {
  if (length(package) != 1L || length(direct) != 1L) {
    stop(
      paste0(
        "Audit value must be scalar for ", scope, " / ", input, " / ", model,
        " / ", component, "; package length=", length(package),
        ", direct length=", length(direct), "."
      ),
      call. = FALSE
    )
  }
  difference <- abs(as.double(package) - as.double(direct))
  status <- if (!is.na(difference) && is.finite(difference) && difference <= tolerance) "PASS" else "FAIL"
  audit_rows[[length(audit_rows) + 1L]] <<- data.frame(
    Scope = scope,
    Input = input,
    Model = model,
    Component = component,
    Package = as.double(package),
    Direct_BGLR = as.double(direct),
    Abs_diff = difference,
    Status = status,
    Direct_source = source,
    stringsAsFactors = FALSE
  )
  invisible(status)
}

print_section <- function(title, rows) {
  cat("\n", paste(rep("=", 96L), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 96L), collapse = ""), "\n", sep = "")
  shown <- rows[, c(
    "Scope", "Input", "Model", "Component", "Package",
    "Direct_BGLR", "Abs_diff", "Status"
  )]
  shown$Package <- formatC(shown$Package, digits = 8L, format = "fg", flag = "#")
  shown$Direct_BGLR <- formatC(shown$Direct_BGLR, digits = 8L, format = "fg", flag = "#")
  shown$Abs_diff <- formatC(shown$Abs_diff, digits = 3L, format = "e")
  utils::write.table(
    shown,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

quiet_audit <- function(expr) {
  invisible(utils::capture.output(suppressMessages(suppressWarnings(force(expr)))))
}

build_inputs <- function() {
  set.seed(20260815)
  ids <- paste0("g", seq_len(10L))
  marker <- matrix(
    stats::rnorm(length(ids) * 8L),
    nrow = length(ids),
    dimnames = list(ids, paste0("m", seq_len(8L)))
  )
  omic1 <- matrix(
    stats::rnorm(length(ids) * 5L),
    nrow = length(ids),
    dimnames = list(ids, paste0("transcript", seq_len(5L)))
  )
  omic2 <- matrix(
    stats::rnorm(length(ids) * 4L),
    nrow = length(ids),
    dimnames = list(ids, paste0("metabolite", seq_len(4L)))
  )
  signal <- as.double(scale(marker[, 1L] - 0.55 * marker[, 2L] + 0.25 * omic1[, 1L]))
  single <- data.frame(
    GID = ids,
    Yield = 3 + 0.8 * signal + stats::rnorm(length(ids), sd = 0.35),
    stringsAsFactors = FALSE
  )
  single$Yield[c(3L, 9L)] <- NA_real_

  make_kernel <- function(x, ridge = 0.05) {
    z <- scale(x, center = TRUE, scale = FALSE)
    k <- tcrossprod(z)
    k <- k / mean(diag(k)) + diag(ridge, nrow(k))
    rownames(k) <- colnames(k) <- rownames(x)
    k
  }
  k_g <- make_kernel(marker)
  k_o1 <- make_kernel(omic1)
  k_o2 <- make_kernel(omic2)

  met <- expand.grid(
    GID = ids,
    Env = c("E1", "E2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  met_signal <- signal[match(met$GID, ids)]
  met$Yield <- 2 + 0.75 * met_signal + ifelse(met$Env == "E2", 0.65, 0) +
    stats::rnorm(nrow(met), sd = 0.28)

  mt <- data.frame(
    GID = ids,
    T1 = 2 + 0.80 * signal + stats::rnorm(length(ids), sd = 0.30),
    T2 = 1 - 0.55 * signal + stats::rnorm(length(ids), sd = 0.35),
    stringsAsFactors = FALSE
  )
  mt$T1[2L] <- NA_real_
  mt$T2[8L] <- NA_real_

  mtmet <- met[, c("GID", "Env"), drop = FALSE]
  mtmet$T1 <- 2 + 0.75 * met_signal + ifelse(mtmet$Env == "E2", 0.55, 0) +
    stats::rnorm(nrow(mtmet), sd = 0.30)
  mtmet$T2 <- 1 - 0.50 * met_signal + ifelse(mtmet$Env == "E2", -0.35, 0) +
    stats::rnorm(nrow(mtmet), sd = 0.32)
  mtmet$T1[c(2L, 16L)] <- NA_real_
  mtmet$T2[c(5L, 13L)] <- NA_real_

  list(
    marker = marker,
    omic1 = omic1,
    omic2 = omic2,
    single = single,
    met = met,
    mt = mt,
    mtmet = mtmet,
    k_g = k_g,
    k_o1 = k_o1,
    k_o2 = k_o2
  )
}

saved_indices <- function(para) {
  gp_bayes_saved_draw_indices(para[["nIter"]], para[["burnIn"]], para[["thin"]])
}

read_beta <- function(mod, eta_idx, draw_indices) {
  path <- mod$output_files_names[grepl(
    paste0("ETA_", eta_idx, "_b\\.bin$"),
    basename(mod$output_files_names)
  )]
  stopifnot(length(path) == 1L)
  draws <- as.matrix(BGLR::readBinMat(path))
  if (length(draw_indices) && max(draw_indices) <= nrow(draws)) {
    draws <- draws[draw_indices, , drop = FALSE]
  }
  draws
}

read_var_u <- function(mod, eta_idx, draw_indices) {
  path <- mod$output_files_names[grepl(
    paste0("ETA_", eta_idx, "_varU\\.dat$"),
    basename(mod$output_files_names)
  )]
  stopifnot(length(path) == 1L)
  process_var_u(path, draw_indices, "RKHS")
}

direct_marker_draws <- function(mod, eta, para) {
  draw_indices <- saved_indices(para)
  random_idx <- which(vapply(eta$ETA, function(x) x$model != "FIXED", logical(1L)))
  train_idx <- which(!is.na(mod$model$y))
  component <- list()
  total_effect <- NULL
  for (idx in random_idx) {
    beta <- read_beta(mod, idx, draw_indices)
    x <- as.matrix(eta$ETA[[idx]]$X)
    effect <- x %*% t(beta)
    component[[length(component) + 1L]] <- apply(
      effect[train_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE
    )
    total_effect <- if (is.null(total_effect)) effect else total_effect + effect
  }
  list(
    component = component,
    total = apply(total_effect[train_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE),
    residual = gp_bayes_read_varE_draws(mod$output_files_names, draw_indices = draw_indices)
  )
}

audit_marker_fit <- function(dat, model, multi_omics = FALSE, seed = 1L) {
  input_label <- if (isTRUE(multi_omics)) "markers + 2 omics matrices" else "markers"
  args <- list(
    random = ~ GID,
    GS_model = model,
    response = "Yield",
    pheno_data = dat$single,
    geno_data = dat$marker,
    gen_name = "GID",
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin,
    cross_validation = TRUE
  )
  if (isTRUE(multi_omics)) {
    args$omic1_data <- dat$omic1
    args$omic2_data <- dat$omic2
  }
  set.seed(seed)
  prep <- do.call(bayes_finalize_A_B_C_BL_BRR, args)
  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = model
  )
  direct <- direct_marker_draws(mod, prep$bayes_ETA, prep$bayes_para)
  out <- mod_output_bayes(
    mod = mod,
    ETA = prep$bayes_ETA,
    pheno_data = prep$bayes_ETA$pheno_data,
    geno_data = dat$marker,
    omic1_data = if (isTRUE(multi_omics)) dat$omic1 else NULL,
    omic2_data = if (isTRUE(multi_omics)) dat$omic2 else NULL,
    gen_name = "GID",
    bayes_para = prep$bayes_para,
    GS_model = model
  )
  vc <- out$Variance_components
  genetic_rows <- grep("^genetic_variance", rownames(vc))
  stopifnot(length(genetic_rows) == length(direct$component))
  for (j in seq_along(genetic_rows)) {
    add_audit(
      "Single trait", input_label, model, rownames(vc)[genetic_rows[j]],
      vc$Components[genetic_rows[j]], mean(direct$component[[j]], na.rm = TRUE),
      "posterior mean var(X beta) on training rows"
    )
  }
  if (length(direct$component) > 1L) {
    add_audit(
      "Single trait", input_label, model, "total_genetic_variance",
      vc["total_genetic_variance", "Components"], mean(direct$total, na.rm = TRUE),
      "posterior mean var(sum X beta) on training rows"
    )
  }
  add_audit(
    "Single trait", input_label, model, "residual_variance",
    vc["residual_variance", "Components"], mean(direct$residual, na.rm = TRUE),
    "posterior mean of BGLR varE chain"
  )
}

direct_kernel_draws <- function(mod, eta, para) {
  draw_indices <- saved_indices(para)
  random_idx <- which(vapply(eta$ETA, function(x) x$model != "FIXED", logical(1L)))
  train_idx <- which(!is.na(mod$model$y))
  component <- list()
  total_effect <- NULL
  all_brr <- TRUE
  for (idx in random_idx) {
    term <- eta$ETA[[idx]]
    if (identical(term$model, "BRR")) {
      beta <- read_beta(mod, idx, draw_indices)
      effect <- as.matrix(term$X) %*% t(beta)
      component[[length(component) + 1L]] <- apply(
        effect[train_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE
      )
      total_effect <- if (is.null(total_effect)) effect else total_effect + effect
    } else if (identical(term$model, "RKHS")) {
      all_brr <- FALSE
      scale_k <- mean(diag(as.matrix(term$K))[train_idx], na.rm = TRUE)
      component[[length(component) + 1L]] <- read_var_u(mod, idx, draw_indices) * scale_k
    }
  }
  total <- if (isTRUE(all_brr)) {
    apply(total_effect[train_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
  } else {
    Reduce(`+`, component)
  }
  list(
    component = component,
    total = total,
    residual = gp_bayes_read_varE_draws(mod$output_files_names, draw_indices = draw_indices)
  )
}

audit_kernel_fit <- function(dat, model, layout, seed = 1L) {
  kernel_args <- switch(
    layout,
    single_kernel = list(gmatrix = dat$k_g),
    multiple_kernels = list(gmatrix = dat$k_g, kernel_list = list(auxiliary = dat$k_o1)),
    multi_omics = list(
      gmatrix = dat$k_g,
      omic1_kernel = dat$k_o1,
      omic2_kernel = dat$k_o2
    )
  )
  args <- c(list(
    random = ~ GID,
    GS_model = model,
    response = "Yield",
    pheno_data = dat$single,
    gen_name = "GID",
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin,
    cross_validation = TRUE
  ), kernel_args)
  set.seed(seed)
  prep <- do.call(bayes_finalize_RKHS_GBLUPBRR, args)
  internal_model <- if (identical(model, "GBLUP_BRR")) "BRR" else model
  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = internal_model
  )
  direct <- direct_kernel_draws(mod, prep$bayes_ETA, prep$bayes_para)
  out <- do.call(mod_output_bayes_gbluBRR_RKHS, c(list(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = internal_model,
    gen_name = "GID",
    pheno_data = prep$bayes_ETA$pheno_data,
    bayes_para = prep$bayes_para
  ), kernel_args))
  vc <- out$Variance_components
  genetic_rows <- grep("^genetic_variance", rownames(vc))
  stopifnot(length(genetic_rows) == length(direct$component))
  for (j in seq_along(genetic_rows)) {
    add_audit(
      "Single trait", gsub("_", " ", layout), model, rownames(vc)[genetic_rows[j]],
      vc$Components[genetic_rows[j]], mean(direct$component[[j]], na.rm = TRUE),
      if (identical(internal_model, "RKHS")) {
        "posterior mean varU x mean diag(K)"
      } else {
        "posterior mean var(X beta) on training rows"
      }
    )
  }
  if (length(direct$component) > 1L) {
    add_audit(
      "Single trait", gsub("_", " ", layout), model, "total_genetic_variance",
      vc["total_genetic_variance", "Components"], mean(direct$total, na.rm = TRUE),
      if (identical(internal_model, "RKHS")) {
        "sum of independent kernel variance draws"
      } else {
        "posterior mean var(sum X beta) on training rows"
      }
    )
  }
  add_audit(
    "Single trait", gsub("_", " ", layout), model, "residual_variance",
    vc["residual_variance", "Components"], mean(direct$residual, na.rm = TRUE),
    "posterior mean of BGLR varE chain"
  )
}

direct_multitrait_covariances <- function(fit, eta_bundle) {
  by_role <- function(role) {
    nms <- names(eta_bundle$term_roles)[eta_bundle$term_roles == role]
    pieces <- lapply(nms, function(nm) {
      bayes_multitrait_eta_scale(eta_bundle$ETA[[nm]]) *
        as.matrix(fit$ETA[[nm]]$Cov$Omega)
    })
    if (length(pieces)) Reduce(`+`, pieces) else NULL
  }
  list(
    main = by_role("genetic"),
    gxe = by_role("gxe"),
    residual = as.matrix(fit$resCov$R)
  )
}

run_direct_multitrait <- function(pheno, responses, kernels, model,
                                  heter_groups = NULL, seed = 1L) {
  kernels <- bayes_multitrait_prepare_kernels(kernels)
  kernel_ids <- Reduce(intersect, lapply(kernels, rownames))
  ph <- bayes_multitrait_expand_pheno(
    pheno_data = pheno,
    response = responses,
    gen_name = "GID",
    heter_groups = heter_groups,
    kernel_ids = kernel_ids
  )
  ph <- ph[as.character(ph$GID) %in% kernel_ids, , drop = FALSE]
  internal_model <- if (identical(model, "GBLUP_BRR")) "BRR" else model
  eta_bundle <- bayes_multitrait_build_eta(
    pheno_data = ph,
    response = responses,
    gen_name = "GID",
    heter_groups = heter_groups,
    kernels = kernels,
    GS_model = internal_model,
    trait_cov_type = "UN"
  )
  prefix <- paste0("audit_direct_mt_", seed, "_")
  set.seed(seed)
  # bayes_multitrait_joint_fit() consumes one draw while creating its unique
  # BGLR save prefix. Mirror that draw so both fits receive the identical MCMC
  # random stream and can be compared at machine precision.
  sample.int(999999L, 1L)
  fit <- BGLR::Multitrait(
    y = as.matrix(ph[, responses, drop = FALSE]),
    ETA = eta_bundle$ETA,
    resCov = bayes_multitrait_cov_prior("UN", length(responses)),
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin,
    verbose = FALSE,
    saveAt = prefix
  )
  files <- list.files(pattern = paste0("^", prefix))
  on.exit(unlink(files), add = TRUE)
  list(covariance = direct_multitrait_covariances(fit, eta_bundle), kernels = kernels)
}

audit_joint_multitrait <- function(dat, layout, is_met = FALSE, model = "RKHS", seed = 1L) {
  kernels <- switch(
    layout,
    single_kernel = list(genomic = dat$k_g),
    multiple_kernels = list(genomic = dat$k_g, auxiliary = dat$k_o1),
    multi_omics = list(genomic = dat$k_g, transcriptomic = dat$k_o1, metabolomic = dat$k_o2)
  )
  pheno <- if (isTRUE(is_met)) dat$mtmet else dat$mt
  heter_groups <- if (isTRUE(is_met)) "Env" else NULL
  responses <- c("T1", "T2")
  direct <- run_direct_multitrait(
    pheno, responses, kernels, model,
    heter_groups = heter_groups,
    seed = seed
  )
  set.seed(seed)
  out <- bayes_multitrait_joint_fit(
    pheno_data = pheno,
    response = responses,
    gen_name = "GID",
    heter_groups = heter_groups,
    kernels = kernels,
    GS_model = model,
    bayes_para = list(nIter = n_iter, burnIn = burn_in, thin = thin),
    trait_cov_type = "UN",
    residual_cov_type = "UN",
    verbose = FALSE
  )
  vc <- out$bayes_result$Variance_components
  scope <- if (isTRUE(is_met)) "Joint multi-trait MET" else "Joint multi-trait"
  for (trait in responses) {
    trait_rows <- vc$Trait == trait
    for (component in c("genetic_variance", if (isTRUE(is_met)) "gxe_variance", "residual_variance")) {
      package_value <- vc$Components[trait_rows & vc$Component == component]
      if (length(package_value) != 1L) {
        cat("\nUnexpected joint variance table for ", scope, " / ", layout, " / ", model, ":\n", sep = "")
        print(vc, row.names = FALSE)
      }
      direct_value <- switch(
        component,
        genetic_variance = diag(direct$covariance$main)[match(trait, responses)],
        gxe_variance = diag(direct$covariance$gxe)[match(trait, responses)],
        residual_variance = diag(direct$covariance$residual)[match(trait, responses)]
      )
      add_audit(
        scope, gsub("_", " ", layout), model,
        paste(trait, component, sep = ":"), package_value, direct_value,
        switch(
          component,
          genetic_variance = "sum scale(K or X) x BGLR ETA Omega diagonal",
          gxe_variance = "sum scale(GxE K or X) x BGLR ETA Omega diagonal",
          residual_variance = "diagonal of BGLR resCov R"
        )
      )
    }
  }
}

direct_heter_met <- function(out, kernels) {
  fit <- out$bayes_model$model$multitrait_fit
  genetic <- Reduce(`+`, Map(function(nm, k) {
    mean(diag(k), na.rm = TRUE) * as.matrix(fit$ETA[[nm]]$Cov$Omega)
  }, names(kernels), kernels))
  residual <- as.matrix(fit$resCov$R)
  if (is.null(rownames(residual)) || is.null(colnames(residual))) {
    dimnames(residual) <- dimnames(genetic)
  }
  list(genetic = genetic, residual = residual)
}

audit_met <- function(dat, layout, model = "RKHS", seed = 1L) {
  kernels <- switch(
    layout,
    single_kernel = list(gmatrix = dat$k_g),
    multiple_kernels = list(gmatrix = dat$k_g, kernel_list = list(auxiliary = dat$k_o1)),
    multi_omics = list(gmatrix = dat$k_g, omic1_kernel = dat$k_o1, omic2_kernel = dat$k_o2)
  )
  set.seed(seed)
  out <- do.call(bayes_finalize_RKHS_GBLUPBRR, c(list(
    random = ~ GID + GID:Env,
    GS_model = model,
    response = "Yield",
    pheno_data = dat$met,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin
  ), kernels))
  kernel_values <- gp_collect_kernel_inputs(
    gmatrix = kernels$gmatrix,
    kernel_list = kernels$kernel_list,
    omic1_kernel = kernels$omic1_kernel,
    omic2_kernel = kernels$omic2_kernel
  )
  direct <- direct_heter_met(out, kernel_values)
  vc <- out$bayes_result$Variance_components
  for (env in colnames(direct$genetic)) {
    for (component in c("genetic_variance", "residual_variance")) {
      package_value <- vc$Components[vc$Env == env & vc$Component == component]
      direct_value <- if (identical(component, "genetic_variance")) {
        diag(direct$genetic)[match(env, colnames(direct$genetic))]
      } else {
        diag(direct$residual)[match(env, colnames(direct$residual))]
      }
      add_audit(
        "MET", gsub("_", " ", layout), model,
        paste(env, component, sep = ":"), package_value, direct_value,
        if (identical(component, "genetic_variance")) {
          "sum mean diag(K) x BGLR ETA Omega diagonal"
        } else {
          "diagonal of BGLR heterogeneous resCov R"
        }
      )
    }
  }
}

dat <- build_inputs()
only_joint <- identical(Sys.getenv("PREDICTPROR_AUDIT_ONLY_JOINT"), "1")

cat("BGLR variance-extraction audit\n")
cat("MCMC:", n_iter, "iterations;", burn_in, "burn-in; thinning", thin, "\n")
cat("Tolerance for package-vs-direct extraction:", tolerance, "\n")

if (!isTRUE(only_joint)) {
  marker_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  for (j in seq_along(marker_models)) {
    quiet_audit(audit_marker_fit(dat, marker_models[j], multi_omics = FALSE, seed = 100L + j))
    quiet_audit(audit_marker_fit(dat, marker_models[j], multi_omics = TRUE, seed = 120L + j))
  }

  seed <- 200L
  for (model in c("RKHS", "GBLUP_BRR")) {
    for (layout in c("single_kernel", "multiple_kernels", "multi_omics")) {
      seed <- seed + 1L
      quiet_audit(audit_kernel_fit(dat, model, layout, seed = seed))
    }
  }

  seed <- 300L
  for (model in c("RKHS", "GBLUP_BRR")) {
    for (layout in c("single_kernel", "multiple_kernels", "multi_omics")) {
      seed <- seed + 1L
      quiet_audit(audit_met(dat, layout, model, seed = seed))
    }
  }
}

seed <- 400L
for (is_met in c(FALSE, TRUE)) {
  for (model in c("RKHS", "GBLUP_BRR")) {
    for (layout in c("single_kernel", "multiple_kernels", "multi_omics")) {
      seed <- seed + 1L
      quiet_audit(audit_joint_multitrait(dat, layout, is_met = is_met, model = model, seed = seed))
    }
  }
}

audit <- do.call(rbind, audit_rows)
rownames(audit) <- NULL
for (scope in unique(audit$Scope)) {
  print_section(scope, audit[audit$Scope == scope, , drop = FALSE])
}

cat("\nSupported-scope note:\n")
cat("- Raw marker BayesA/BayesB/BayesC/BL/BRR is single-environment only.\n")
cat("- MET and joint multi-trait BGLR use RKHS or marker-derived GBLUP_BRR kernels.\n")
cat("- Multi-omics is audited both as separate marker matrices (single trait) and separate kernels.\n")

failed <- audit$Status != "PASS"
cat("\nAUDIT SUMMARY:", sum(!failed), "passed;", sum(failed), "failed;", nrow(audit), "comparisons.\n")
if (any(failed)) {
  print(audit[failed, , drop = FALSE], row.names = FALSE)
  quit(status = 1L)
}
cat("audit_bglr_variance_extraction: ok\n")
