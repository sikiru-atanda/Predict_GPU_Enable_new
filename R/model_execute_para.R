#' Title
#'
#' @param pheno_data phenotypic data, which can be contain both training and testing set. NA is allowed. Dataframe or matrix is allowed
#' @param geno_data Genomic/SNP/Marker data NA is allowed but not expected.
#' numeric 0, 1, 2 (where 0 is minor allele, 1 is hetero and 2 is major allele)
#' and -1, 0, 1 is also allowed (where -1 is minor allele, 0 is hetero and 1 is major allele).
#'  Dataframe or matrix is allowed
#'  We allowed up to 4 different omics data for model fit
#' @param omic1_data Omic data (transcriptomic, metabolic, proteomic, environment etc) NA is allowed but not expected. Dataframe or matrix is allowed
#' @param omic2_data Similar to Omic1_data
#' @param omic3_data Similar to Omic1_data
#' @param gmatrix Genomic relationship matrix, with no missing values. When
#'   several kernels are supplied through `gmatrix`, `gkernel`, the named omics
#'   kernel arguments, or `kernel_list`, PredictProR retains every kernel.
#'   Covariance models fit separate kernel terms; ML/DL models use separate
#'   named eigenfeature blocks.
#' @param kernel_list Optional named list of additional relationship or kernel
#' matrices. Entries must be square numeric matrices keyed by the same IDs as
#' `gen_name`. ASReml, GP, and kernel-Bayesian routes treat them as covariance
#' sources according to their model-specific parameterization. ML/DL routes
#' eigen-decompose each kernel and concatenate the named feature blocks with
#' any raw genomic/omics features. ML/DL variance output remains predictive,
#' not genetic variance or heritability.
#' @param train_geno_data Genomic data for training set if geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param train_omic1_data Omic data for training set
#' @param train_omic2_data Omic data for training set
#' @param train_omic3_data Omic data for training set
#' @param test_geno_data  Genomic data for testing set if geno_data is not provided or not included in the geno_data by the user.
#' In that scenario geno_data is considered training set and training_geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param test_omic1_data Omic data for testing set if omic1_data is not provided or not included in the omic1_data by the user.
#' In that scenario omic1_data is considered training set and training_geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param test_omic2_data same as test_omic1_data
#' @param test_omic3_data same as test_omic1_data
#' @param train_set Dataframe with column name of the individual in the training set. This is useful maining
#' for purpose of cross-validation exercise.
#' @param omics_data_label  lable/name of the omics data
#' @param omics_kernel_label  lable/name for the omics_kernel if any
#' @param train_omics_label
#' @param test_omics_label
#' @param test_set Dataframe with column name of the individual in the testing set. Not required
#' if pheno_data contain individuals (testing set) with no phenotypic record as NA.
#' @param gmatrix_method Character vector of genomic relationship matrices to
#' calculate from marker data. Supported methods include additive VanRaden,
#' weighted VanRaden, Yang, epistasis, and dominance relationship matrices.
#' When `met_ml_dl = TRUE` and the user supplies `geno_data` but no
#' precomputed kernel and no explicit `gmatrix_method`, the input guardrail
#' auto-sets `gmatrix_method = "Yang"` so the MET ML/DL kernel-feature
#' pipeline has a kernel to PCA over. Set this argument explicitly to
#' override the default (e.g. `"VanRaden"`).
#' @param response trait(s) of interest to the user
#' @param gen_name Column name containing individuals/genotypes
#' @param cova covariate if any.Its epected in formula i.e cova  = ~ Rain + Temp
#' @param fixed fixed terms. Its expected in formula i.e fixed = ~ name + Env
#' @param random random terms. Its expected in formula i.r random = ~ name + Env
#' @param heter_resid True or False if user want heterogeneous residual variance or not
#' @param bayes_kernel_heter_resid Optional override for the kernel-Bayes
#'   (GBLUP_BRR / RKHS) MET routing decision, independent of
#'   `heter_resid`. `NULL` (default) inherits from `heter_resid`. `TRUE`
#'   forces the env-as-trait `BGLR::Multitrait` path (heterogeneous σ²_g
#'   per env + per-env σ²_e) and requires `weights = NULL` because
#'   `BGLR::Multitrait` has no observation-weight argument. `FALSE` forces the
#'   univariate `BGLR::BGLR`
#'   path with `~GID + GID:Loc` (single σ²_g, single σ²_gxe, single σ²_e
#'   -- the BGLR-inherent homogeneous case) and is the RKHS mode that accepts
#'   Stage 2 observation weights. ASReml and GP families keep following
#'   `heter_resid` regardless.
#' @param heter_groups  Column name for Environment or location
#' @param weights Optional Stage 2 observation precision weights for Gaussian
#'   GP, ASReml, and univariate BGLR fits. Supply a positive numeric vector, the
#'   name of a numeric column in the phenotype data (for example,
#'   `weights = "Weight"`), or a one-column data frame/matrix. PredictProR uses
#'   `Var(e_i) = 1 / weights_i` for GP and ASReml (ASReml is fitted with
#'   `asr_gaussian(dispersion = 1)`) and passes `sqrt(weights_i)` to BGLR,
#'   whose native contract is `Var(e_i)` proportional to the inverse squared
#'   BGLR weight. Unsupported weighted routes fail instead of ignoring weights.
#' @param nIter  number of iteration for Bayesian models
#' @param burnIn number of burnin  for Bayesian models
#' @param thin   number of thinning for Bayesian models
#' @param GS_model GS-model for fit. The following are available
#' BayesA, BayesB, BayesC, Baysian Ridge Regression (BRR). These models only work with
#' M-matrix(genomic data and omic data) in single location.
#' Bayesian reproducing kernel Hilbert spaces regressions (RKHS),
#' Bayesian Genomic Best linear unbias estimate (BGBLUP). Both RKHS and BGBLUP can
#' fit both single and multiple location using reaction norm.
#' Genomic Best linear unbias estimate using asreml-R package with different
#' variance structure such as (FA, RR, US, CORGH, CORH, CORV)
#' for multi-location/environment.
#' Machine learning models include:
#' Extreme Gradiant Boosting, Random Forest, KNN, Lasso, Ridge Regression,
#' Partial Least Square, Support Vector Machine. All the machine learning only work
#' in single location.
#' @param fixed_term_model_bayesian model for the fixed term which is always fixed
#' @param rand_term_model_bayesian model for the random terms which can be any of the above mentioned model
#' @param message if message/warning should be displayed
#' @param gkernel Genomic kernel. If both `gmatrix` and `gkernel` are supplied,
#'   both are retained as distinct named kernel sources rather than one taking
#'   precedence over the other.
#' @param kernel_method Character vector of kernel methods to calculate for
#' omics data. Supported methods include Gaussian, linear, composite,
#' polynomial, Matern-family, Laplacian, and rational quadratic kernels.
#' @param omic1_kernel  relationship matrix using different kernel methods
#' @param omic2_kernel  relationship matrix using different kernel methods
#' @param omic3_kernel  relationship matrix using different kernel methods
#' @param pheno_data_train phenotypic data for training set. NA not allowed. Dataframe or matrix is allowed
#' @param pheno_data_test phenotypic data for the testing set. NA is allowed. Dataframe or matrix is allowed
#' @param cross_validation
#' @param cv_generate_plots Logical or `NULL`. When `NULL`, GP-only
#' cross-validation defaults to a point-prediction metric path without building
#' diagnostic plot objects. Set `TRUE` to force legacy CV plot construction or
#' `FALSE` to skip it for faster CV output.
#' @param coefficient_1 coefficient for the training set using either genomic or any omics data.
#' We allowed up to 4 omics data for model fit
#' @param coefficient_2
#' @param coefficient_3
#' @param coefficient_4
#' @param eval_metrics
#' @param response_family Response family for the target trait. One of `auto`,
#' `gaussian`, `binary`, `ordinal`, or `multiclass`.
#' @param positive_class For binary responses, the class treated as the event
#'   for precision, recall, specificity, F1, MCC, Brier score, and one-column
#'   probability vectors. The value must match one observed class. When omitted,
#'   the second factor level is used and recorded in the CV selection output.
#' @param gp_backend Gaussian-process backend route. For `hybrid_gp`, `"auto"`
#' uses the Python/Torch additive-kernel backend when available and falls back
#' to the R solver; use `"r"` to force the R solver or `"python"`/`"torch"` to
#' require the Python/Torch backend.
#' @param gp_return_se Logical. For the non-cross-validation direct
#' Gaussian-process prediction path, request prediction standard errors and
#' prediction error variances in the public prediction output. Cross-validation
#' always uses point predictions only.
#' @param gp_iters Optional positive integer optimizer iterations for the exact
#' GP backend. `NULL` preserves the Python backend default.
#' @param gp_lr Optional positive learning rate for the exact GP backend.
#' `NULL` preserves the Python backend default.
#' @param gp_learn_scales Optional logical override for GP backend scale
#' learning. `NULL` learns scales only when full variance-component output is
#' requested.
#' @param gp_engine Character; which REML engine drives GP variance-component
#' estimation. `"auto"` (default) picks the fastest available (currently
#' resolves to `"dense_v"`); `"dense_v"` keeps the dense V-formulation
#' (Phase 1-2); `"mme"` uses the sparse Mixed Model Equations engine
#' (Phase 3.4-3.7, fastest at n_geno = 100 for MET but regresses at
#' n_geno >= 300 on GBLUP because G^-1 is dense -- see Phase 3.14 NEWS);
#' `"eigen"` (Phase 3.15a, experimental) uses an eigen-projected REML
#' formulation that beats ASReml at large n on the standalone validation,
#' but production integration is incomplete (Phase 3.15b NEWS) so
#' `"eigen"` currently falls back to `"dense_v"`. Requires
#' `PREDICTPRO_GP_DISABLE_ENV_MAIN=1` for `"mme"` and `"eigen"` on MET.
#' Honors `PREDICTPRO_GP_ENGINE` env var.
#' @param gp_exact_fast_cv Optional control for exact-GP cross-validation speed.
#' `NULL` keeps the backend default, which is exact per-fold refitting unless
#' explicitly enabled by environment policy. `TRUE` or `"auto"` enables the
#' guarded shared-hyperparameter block-delete shortcut for eligible single-
#' environment GP CV; `"force"` skips the self-check and is intended for
#' benchmarking only. The shortcut is reported as approximate in result
#' metadata.
#' @param gp_force_prediction_se Logical. Joint multi-trait Gaussian-process
#'   prediction standard errors and trait correlations cost one linear solve per
#'   genotype-by-environment test cell per trait, so they are capped: past
#'   roughly 500 test cells at two traits the fit returns predictions without
#'   them and warns, rather than spending hours. Set `TRUE` to compute them
#'   anyway. Leaving this `FALSE` never fails the run; it only omits the
#'   per-cell trait covariance, and the omission is reported rather than silent.
#' @param gp_estimate_kernel_weights Logical. Joint multi-trait Gaussian-process
#'   routes combine their kernel bank into one genetic term. By default the
#'   mixture weights are fixed (equal unless supplied), so each kernel's share
#'   is an assumption rather than an estimate. Set `TRUE` to fit the relative
#'   weights by REML profile likelihood; the first kernel is held at 1 because
#'   the overall genetic scale is already carried by the trait covariance, so
#'   the fitted values are ratios. Requires `gp_varcomp_mode = "reml"`, at least
#'   two kernels, and a single-environment panel: the multi-environment joint GP
#'   backend estimates its variance components by method of moments and has no
#'   likelihood to profile, so it keeps fixed weights and warns. The fitted
#'   weights, and the log-likelihood gained over the starting weights, are
#'   reported in `model_parameters`; a near-zero gain means the data does not
#'   identify the mixture.
#' @param gp_return_trait_correlations Logical. For non-cross-validation
#' multi-trait GP and MT-MET GP paths, request genetic and GxE/residual trait
#' correlation matrices where supported by the backend.
#' @param gp_full_vc Logical. For the non-cross-validation direct
#' Gaussian-process prediction path, request full variance-component output
#' where the backend supports it. This implies `gp_return_se`.
#' @param multi_trait Logical. High-level automatic multi-trait orchestration.
#' If `TRUE`, provide one or more models in `GS_model` (true prediction) or
#' `GS_model_cv` (cross-validation) and leave all five family-specific
#' `multi_trait_*` flags `FALSE`. PredictProR resolves each model to its existing
#' protected ASReml, GP, Bayesian, ML, or DL route and executes one model per
#' child call. With cross-validation and `cv_evaluation_only = FALSE`, models
#' are ranked by the primary metric averaged across traits and the best joint
#' model is refitted for final prediction. The family-specific flags remain
#' available for backward-compatible manual single-model routing. Current
#' single-environment joint CV candidates include `GBLUP`, `GBLUP_BRR`,
#' `RKHS`, `Gaussian-Process-GBLUP`, `FA-GBLUP`, `Scalable-GBLUP`, and the
#' supported multi-trait ML/DL models. Joint GP/Bayesian CV holds out complete
#' genotypes and masks every response trait together before each model refit.
#' Automatic joint CV accepts `K-Folds` or `Repeated_K-Folds`; its protected
#' routes use unstratified genotype-level folds. Duplicate names in `response`
#' are reduced to unique traits before child models are launched. Grouped
#' multi-trait multi-environment CV is available for ASReml GBLUP; the GP and
#' Bayesian MT-MET routes remain true-prediction only.
#' @param multi_trait_gp Logical. If `TRUE`, enable the dedicated joint
#' Gaussian-process multi-trait path for true prediction or single-environment
#' genotype-blocked CV. This is different from passing multiple response
#' columns for ordinary independent per-trait execution.
#' @param multi_trait_bayes Logical. If `TRUE`, enable the dedicated joint
#' Bayesian BGLR multi-trait path for `GBLUP_BRR` or `RKHS`, including
#' single-environment genotype-blocked CV. This is different from passing
#' multiple response columns for ordinary independent per-trait execution.
#' @param multi_trait_asreml Logical. If `TRUE`, enable the dedicated unbalanced
#' Gaussian multi-trait `ASReml-R` path (`GS_model = "GBLUP"` with
#' `engine = "asreml"`). Supports single-environment and grouped MT-MET true
#' prediction plus genotype-blocked `K-Folds` or `Repeated_K-Folds` CV. For
#' automatic `multi_trait = TRUE` comparison, `cv_evaluation_only = FALSE`
#' evaluates candidates and then refits the selected model for final
#' prediction. Each fold masks all trait-environment responses for its held-out
#' genotypes. Fold fits return predictions only; biological variance components
#' are extracted from the stable final/direct fit.
#' @param hybrid_asreml Logical. If `TRUE`, enable the dedicated Gaussian
#' hybrid `ASReml-R` true-prediction path with explicit female GCA, male GCA,
#' and SCA covariance construction.
#' @param hybrid_bayes Logical. If `TRUE`, enable the dedicated Gaussian
#' hybrid Bayesian-kernel path using explicit female GCA, male GCA, and SCA
#' covariance construction for `GBLUP_BRR` and `RKHS`.
#' @param hybrid_gp Logical. If `TRUE`, enable the dedicated Gaussian hybrid
#' GP/KRR path using explicit female GCA, male GCA, and SCA additive kernels.
#' Use `response = c(...)` for joint cross-trait multi-trait hybrid GP and
#' `heter_groups` for hybrid multi-environment data.
#' @param hybrid_ml Logical. If `TRUE`, enable the dedicated Gaussian hybrid
#' machine-learning path using hybrid-level genotype features. Single-environment
#' CV can use parent-aware scenarios; hybrid MET uses CV0, CV1, or CV2.
#' @param hybrid_dl Logical. If `TRUE`, enable the dedicated Gaussian hybrid
#' deep-learning path using hybrid-level genotype features. Single-environment
#' CV can use parent-aware scenarios; hybrid MET uses CV0, CV1, or CV2.
#' @param female_parent Column name in `pheno_data` identifying the female
#' parent of each hybrid.
#' @param male_parent Column name in `pheno_data` identifying the male parent
#' of each hybrid.
#' @param female_geno_data Optional female-parent genotype matrix (row names =
#' parent IDs) for the hybrid paths. Use this when female and male parents come
#' from separate pools. With `hybrid_ml`/`hybrid_dl`, `het_threshold` is applied
#' to these inbred parents and markers that fail in either pool are dropped from
#' the hybrid features; if no hybrid-level `geno_data` is supplied, the hybrid
#' features are the expected F1 dosages `(female + male) / 2`. The parent QC
#' summary is returned as `hybrid_parent_qc`.
#' @param male_geno_data Optional male-parent genotype matrix; see
#' `female_geno_data`. If only one of the two is supplied it is used for both
#' pools.
#' @param hybrid_gp_lambda Ridge/noise ratio for `hybrid_gp`. Use `"auto"` to
#' select from `hybrid_gp_lambda_grid` by deterministic internal CV.
#' @param hybrid_gp_lambda_grid Positive lambda values considered when
#' `hybrid_gp_lambda = "auto"`.
#' @param hybrid_gp_component_weights Optional named numeric weights for
#' `female_gca`, `male_gca`, `sca`, and `gxe` kernels in the hybrid GP
#' additive/reaction-norm kernel.
#' @param env_similarity Optional environment similarity matrix for hybrid GP
#' MET. This defines the environment side of the hybrid-by-environment
#' reaction-norm kernel.
#' @param env_ids Optional environment IDs for an unnamed `env_similarity`
#' matrix.
#' @param env_covariates Optional environment covariate data. When supplied for
#' hybrid GP, PredictPro converts the covariates to an environment similarity
#' kernel after QC. Pass exactly one of `env_similarity` or `env_covariates`.
#' @param reaction_norm_feature_qc Logical; when `TRUE`, apply environment
#' covariate QC before building a reaction-norm environment kernel.
#' @param kenv_kernel Environment-covariate kernel family used when
#' `env_covariates` is supplied. Supported hybrid GP values are `matern32`,
#' `matern52`, `rbf`, and `linear`.
#' @param kenv_bandwidth Positive bandwidth for distance-based environment
#' kernels.
#' @param kenv_kernel_kwargs Optional named list of additional environment
#' kernel arguments.
#' @param multi_trait_ml Logical. If `TRUE`, enable the dedicated unbalanced
#' Gaussian multi-trait machine-learning path. Current scope is
#' `GS_model = "RandomForest"`, `"Ridge_Regression"`, or
#' `"PartialLeastSquare"`.
#' @param para_tunning Logical; run hyper-parameter tuning before fitting the
#'   classical-ML / DL model.
#' @param var_cov_str User-defined variance-covariance structure. For
#'   single-environment multi-trait ASReml models, choose \code{"us"},
#'   \code{"corgh"}, or \code{"diag"}; \code{NULL} defaults to \code{"us"}. For
#'   grouped MT-MET ASReml GBLUP, use a validated factor-analytic structure such
#'   as \code{"fa1"}, \code{"fa2"}, or \code{"fa3"}; \code{NULL} selects rank up to
#'   three according to the number of trait-environment groups.
#' @param engine if user has asreml
#' @param workspace allocate memory for asreml model fit
#' @param pworkspace allocate memory for predict function in asreml
#' @param bending this is important when the relationship matrix is not positive definitive. It fix it for the user. it has be TRUE
#' @param maxit number of iteration for asreml
#' @param pedigree_matrix
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param rcn_cutoff
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param kernel_check_level Kernel QC level. Use `"auto"` for full checks on moderate kernels
#'        and bounded checks on large kernels.
#' @param kernel_sanitize Kernel prefit sanitizer policy. `"auto"` applies bounded ridge
#'        sanitization when exact dense checks are skipped.
#' @param kernel_repair_priority Repair priority used when `kernel_fix_method = "auto"`
#'        and ridge blending cannot make the kernel positive definite. PredictPro tries
#'        the native diagonal-preserving spectral repair first when the kernel is within
#'        `kernel_cpp_repair_size_limit`, then falls back to Matrix `nearPD`.
#' @param kernel_cpp_repair_size_limit Maximum kernel dimension for exact C++ spectral
#'        SPD repair. Larger kernels use bounded sanitizer/fallback policy.
#' @param pheno_data
#' @param pheno_data_train
#' @param pheno_data_test
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param omics_data_label
#' @param gmatrix
#' @param gkernel
#' @param kernel_list Optional named list of additional relationship or kernel
#' matrices.
#' @param pedigree_matrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param omics_kernel_label
#' @param train_geno_data
#' @param train_omic1_data
#' @param train_omic2_data
#' @param train_omic3_data
#' @param train_omics_label
#' @param test_geno_data
#' @param test_omic1_data
#' @param test_omic2_data
#' @param test_omic3_data
#' @param test_omics_label
#' @param coefficient_1
#' @param coefficient_2
#' @param coefficient_3
#' @param coefficient_4
#' @param train_set
#' @param test_set
#' @param gmatrix_method
#' @param kernel_method
#' @param response
#' @param gen_name
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param var_cov_str
#' @param nIter
#' @param burnIn
#' @param thin
#' @param GS_model
#' @param eval_metrics
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param engine
#' @param workspace
#' @param pworkspace
#' @param maxit
#' @param bending
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param rcn_cutoff
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param kernel_check_level Kernel QC level. Use `"auto"` for full checks on moderate kernels
#'        and bounded checks on large kernels.
#' @param kernel_sanitize Kernel prefit sanitizer policy. `"auto"` applies bounded ridge
#'        sanitization when exact dense checks are skipped.
#' @param kernel_repair_priority Repair priority used when `kernel_fix_method = "auto"`
#'        and ridge blending cannot make the kernel positive definite. PredictPro tries
#'        the native diagonal-preserving spectral repair first when the kernel is within
#'        `kernel_cpp_repair_size_limit`, then falls back to Matrix `nearPD`.
#' @param kernel_cpp_repair_size_limit Maximum kernel dimension for exact C++ spectral
#'        SPD repair. Larger kernels use bounded sanitizer/fallback policy.
#' @param message
#' @param system_database this dictate if the output will be created in a folder or as list
#'                         the default is FALSE. Thus output will be folder.
#' @param ... Additional arguments forwarded to the dispatched engine /
#'   workflow; reserved for forward compatibility.
#' @param scale
#' @param inverse
#' @param epsilon
#' @param vcf_file_name
#' @param vcf_file_path
#' @param vcf_file
#' @param hapmap_file_name
#' @param hapmap_file_path
#' @param hapmap
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param recode_format
#' @param out_put_map
#' @param map_data
#' @param qc_filtering
#' @param xgb_paras_tunning,rf_paras_tunning,pls_paras_tunning,svm_paras_tunning,knn_paras_tunning,lasso_paras_tunning,rr_paras_tunning,dpl_paras_tunning
#'   Per-model hyper-parameter tuning grids (XGBoost, RandomForest, PLS,
#'   SVM, KNN, Lasso, Ridge-Regression and deep-learning) used when
#'   `para_tunning = TRUE`. Each is a named list of vectors / sequences.
#' @param feature_scoring Logical; enable independent predictor scoring and
#' dynamic top-k predictor selection after QC/imputation.
#' @param feature_scoring_model Predictor scoring model. One of
#' `"Ridge_Regression"`, `"BayesB"`, or `"RandomForest"`.
#' @param feature_k_grid Integer vector of top-k predictor counts to evaluate in
#' cross-validation. `NULL` uses an automatic capped grid. On standard
#' single-trait/MET routes the grid always also contains the full predictor
#' count (k = all predictors, i.e. no selection), so CV and automatic model
#' choice can conclude that selection does not help. Metrics for several
#' candidate k values describe the candidate curve; selecting k and quoting the
#' same CV score as an unbiased final performance estimate is not nested CV.
#' Specialized hybrid and joint multi-trait routes require one pre-specified k
#' per run and fail rather than pooling several candidates.
#' @param feature_k Optional top-k predictor count for final prediction. Scores
#' and selected sets are calculated separately for every trait. A genuinely
#' joint multi-trait model uses the union of its trait-specific top-k sets and
#' records both the individual sets and the union size.
#' @param feature_scoring_cv Cross-validation scoring policy. `"fixed"` reuses
#' one full-training ranking in every fold; it is intended for final fitting or
#' externally established rankings and is not leakage-free when the same
#' responses are used to assess CV performance. `"fold_internal"` recomputes
#' the ranking using only each outer training fold. `"both"` returns fixed and
#' fold-internal tasks separately for standard model routes; their metrics are
#' never averaged, and fold-internal results drive automatic model selection.
#' Specialized hybrid/joint routes require choosing one policy per run.
#' @param feature_score_metadata Optional reusable metadata from
#' `feature_score_predictors()`. Automatic selection rebuilds marker designs and
#' kernels from named raw marker/omics columns. Precomputed kernels cannot be
#' subset by feature name and therefore fail loudly unless the user supplies an
#' externally rebuilt kernel matching `feature_selected`.
#' @param feature_scoring_seed Integer seed for deterministic feature scoring.
#' @param feature_ridge_lambda Positive ridge penalty used by Gaussian ridge
#' feature scoring.
#' @param feature_bayes_nIter,feature_bayes_burnIn,feature_bayes_thin MCMC
#' controls used by Gaussian BayesB feature scoring.
#' @param rf_n_jobs Internal Python RandomForest jobs. Keep at `1` when
#' PredictProR is already parallelizing over CV/model tasks.
#'
#' @param pedigree_matrix Optional pedigree-derived relationship matrix `A`
#'   (rows/cols = genotypes), used for pedigree-based GBLUP variants.
#' @param female_gmatrix,male_gmatrix Sex-specific genomic relationship
#'   matrices (one per parental population) used by the hybrid GBLUP path.
#' @param hybrid_include_sca Logical; include a specific-combining-ability
#'   (SCA) random term in the hybrid model in addition to female and male GCA.
#' @param train_omics_label,test_omics_label Optional labels for the train /
#'   test omics blocks (mirrors `omics_data_label`).
#' @param coefficient_2,coefficient_3,coefficient_4 Optional weighting
#'   coefficients applied to omics layers 2, 3 and 4 when combining genomic +
#'   multi-omics relationship matrices.
#' @param multi_trait_dl Logical; route to the joint deep-learning multi-trait
#'   model.
#' @param met_ml_dl Logical; enable the environment-aware MET classical-ML / DL
#'   routes for multi-environment panels. PredictProR also auto-enables this
#'   flag when an advertised MET ML/DL model is requested on a valid MET panel.
#'   Every model named in the call must be listed by `met_data_standard()`;
#'   mixed supported/unsupported requests stop before any fitting instead of
#'   dropping the unsupported model. Required inputs: `pheno_data` with repeated
#'   genotype rows across environments, a single `response` column, and an
#'   environment column named via `heter_groups`. Genomic side: any of
#'   `geno_data` (a GRM will be auto-built via `gmatrix_method = "Yang"` if
#'   no kernel is supplied), `omic1_data`/`omic2_data`/`omic3_data`, a
#'   precomputed `gmatrix`/`gkernel`/`omic*_kernel`, or a `kernel_list` of
#'   additional N×N PSD kernels keyed by `gen_name`. Each kernel is eigen-
#'   decomposed per-kernel and the top components (controlled by
#'   `met_kernel_var_explained` / `met_kernel_min_ev` / `met_kernel_max_pcs`)
#'   are concatenated into the feature matrix the ML/DL model consumes. These
#'   routes report predictive uncertainty, not genetic variance, mixed-model
#'   PEV, or heritability.
#' @param met_kernel_var_explained,met_kernel_min_ev,met_kernel_max_pcs
#'   Controls for the MET kernel reduction: target variance explained,
#'   minimum eigenvalue and maximum number of principal components.
#' @param ld_prunning_qc Logical; apply LD-based marker pruning during QC
#'   (note: the historical spelling is kept for backward compatibility).
#' @param docker_nd_usage Reserved Docker / NDSU runtime hint.
#' @param eval_metrics Character vector of evaluation metrics
#'   (e.g. `"accuracy"`, `"rmse"`).
#' @param selected,max_features Optional feature-selection controls
#'   (pre-selected feature subset and maximum number of features).
#' @param scaling,centering Logical; scale to unit variance / mean-centre the
#'   feature matrix before fitting.
#' @param inverse Logical / attribute; treat the supplied kernel as an
#'   inverse relationship matrix.
#' @param epsilon Ridge added to the kernel diagonal before inversion.
#' @param bend_value Eigenvalue floor for kernel bending toward positive
#'   definiteness.
#' @param blending,blending_value Logical / numeric; blend the kernel with the
#'   identity at the given weight to stabilise an ill-conditioned matrix.
#' @param high_diag_cut_off,low_diag_cut_off Kernel-diagnostic cut-offs for
#'   abnormally high / low diagonal entries.
#' @param duplicate_cut_off Correlation threshold above which rows are flagged
#'   as duplicates in the kernel duplicate scan.
#' @param rcn_cutoff Reciprocal-condition-number cut-off above which the
#'   kernel is considered ill-conditioned.
#' @param optimize_diagonal,optimize_duplicate Logical; let the kernel
#'   QC engine optimise the diagonal / duplicate handling automatically.
#' @param kernel_large_n_threshold Row count above which large-kernel
#'   shortcuts are used in QC.
#' @param duplicate_scan Logical; run the C++ duplicate-row scanner.
#' @param duplicate_sample_size,duplicate_block_size,duplicate_max_pairs
#'   Sampling controls for the duplicate scanner (sample size, block size and
#'   maximum pairs evaluated).
#' @param kernel_fix_method Repair method when the kernel is not positive
#'   definite (`"auto"`, `"nearpd"`, `"ridge"`, ...).
#' @param kernel_rcn_check Logical; perform the reciprocal-condition-number
#'   check.
#' @param kernel_nearpd_size_limit Maximum kernel dimension for which the
#'   `nearPD` repair path is attempted.
#' @param kernel_cpp_keep_diag Logical; have the C++ repair preserve the
#'   diagonal.
#' @param kernel_pd_check Logical; run the positive-definiteness check.
#' @param kernel_pd_sample_size Sample size for the PD-check spectral probe.
#' @param kernel_sanitize_value Numeric sanitiser threshold for the kernel
#'   pre-fit policy.
#' @param vcf_file_name,vcf_file_path VCF file name and directory used when
#'   genotypes are read from a VCF rather than supplied as a matrix.
#' @param vcf_file Optional pre-parsed VCF object.
#' @param hapmap_file_name,hapmap_file_path,hapmap HapMap counterparts of the
#'   VCF arguments.
#' @param csv_file_name,csv_file_path A CSV or TXT genotype table (comma or tab
#'   separated) with the VCF marker columns `CHROM`, `POS`, `ID`, `REF`, `ALT`,
#'   `QUAL`, `FILTER`, `INFO`, `FORMAT` followed by one column per sample
#'   (sample names = `gen_name` values). It is converted to VCF with
#'   [convert_csv_to_vcf()] and then follows the VCF route: the same QC,
#'   imputation (including `imputation_method = "beagle"`) and recoding.
#' @param met_predict_all_environments Multi-environment true prediction:
#'   `TRUE` (default) predicts every line in every environment (the full line
#'   x environment grid) with every engine; `Train_Test_Label` is `"Train"`
#'   for observed line x environment records, `"Test"` for records supplied
#'   with `NA` and for lines with no observation, and `"Unobserved"` for
#'   combinations not in `pheno_data` whose line was observed elsewhere.
#'   `FALSE` predicts only the combinations present in `pheno_data`.
#' @param csv_input_coding Coding of the numeric genotypes in `csv_file_name`:
#'   `"alt_dosage"` (0, 1, 2 copies of the ALT allele; default) or
#'   `"centered_dosage"` (-1, 0, 1). VCF GT calls such as `0/1` are also read.
#' @param maf_threshold,het_threshold,ind_call_rate_threshold,snp_call_rate_threshold
#'   QC thresholds for minor-allele frequency, per-SNP heterozygosity and
#'   per-individual / per-SNP call rate. With `hybrid_ml = TRUE` or
#'   `hybrid_dl = TRUE`, hybrid-level genotypes are never filtered on
#'   heterozygosity by default, because they are heterozygous by design. When
#'   `female_geno_data`/`male_geno_data` are supplied, `het_threshold` is
#'   applied to those inbred parents instead; otherwise it applies to the hybrid
#'   rows only when supplied explicitly.
#' @param test_train_genetic_space Logical; compute the test / train genetic
#'   space overlap diagnostic.
#' @param impute,impute_omic Logical; impute missing genotype / omics values.
#' @param imputation_method,impute_knn_k,na_threshold Imputation method
#'   (`"knn"`, `"mean"`, `"median"`, `"mode"`, or `"beagle"` for raw VCF/HapMap input),
#'   KNN neighbour count and the maximum per-column NA rate above which a
#'   column is dropped.
#' @param ploidy One positive integer or `"auto"`. Raw VCF GT calls can infer
#'   a uniform ploidy. Numeric polyploid dosage matrices must supply ploidy
#'   explicitly or carry a `ploidy` attribute. Beagle remains diploid-only.
#' @param beagle_options Named list passed to [impute_genotypes_with_beagle()]
#'   when `imputation_method = "beagle"`. Typical entries are `beagle_jar`,
#'   `output_prefix`, `ref_file`, `map_file`, `nthreads`, and `java_memory`.
#'   Without `output_prefix`, Beagle's intermediate files are written to a
#'   per-run temporary folder, not the working directory.
#' @param recode_format Genotype recoding scheme: ALT dosage, dosage centered
#'   on `ploidy / 2`, or allele frequency. Legacy diploid aliases are accepted.
#' @param ld_pruning,ld_pruning_method LD-pruning toggle and algorithm
#'   (e.g. `"indep-pairwise"`).
#' @param window_size,step_size,r2_threshold LD-pruning window size, step size
#'   and `r^2` threshold.
#' @param use_kb_window Logical; interpret `window_size` in kilobases.
#' @param phased,use_founders LD-pruning options for phased data and
#'   founder-only restriction.
#' @param out_put_map Logical; return the SNP map alongside recoded genotypes.
#' @param map_data Optional marker map data frame (chromosome, position).
#' @param qc_filtering Logical; apply MAF / het / call-rate QC filters.
#' @param optimizer_name Deep-learning optimiser (e.g. `"adam"`, `"sgd"`).
#' @param use_amp Logical; use automatic mixed-precision training.
#' @param max_grad_norm Gradient-norm clip used during DL training.
#' @param auto_class_weights Logical; auto-balance class weights for
#'   classification.
#' @param internal_cv_nfolds,internal_cv_replication Folds and replication
#'   used by multi-trait and MET DL internal CV when calibration is enabled.
#' @param dl_internal_calibration Logical. Run the extra held-out calibration
#'   fits for Gaussian DL true prediction. Defaults to `TRUE`. Set to `FALSE`
#'   to skip those fits; calibrated prediction-error variance, standard errors,
#'   intervals, and reliability then remain unavailable. This does not disable
#'   requested outer cross-validation, tuning, or bootstrap fits.
#' @param cnn_neurons_per_layer,cnn_kernel_size,cnn_dense_layers,cnn_use_max_pool,cnn_pool_kernel,cnn_pool_stride,cnn_pool_padding,cnn_learning_rate,cnn_separable,cnn_dilations,cnn_use_se,cnn_norm_type,cnn_pool_type,cnn_use_global_pool
#'   CNN architecture / training hyper-parameters (depth, kernel / pool
#'   geometry, separable convolutions, dilations, squeeze-excitation,
#'   normalisation, global pooling and learning rate).
#' @param resnet_neurons_per_block,resnet_blocks,resnet_learning_rate
#'   ResNet hyper-parameters.
#' @param ft_d_model,ft_heads,ft_layers,ft_ff_mult,ft_dropout,ft_token_dropout,ft_use_cls
#'   FT-Transformer hyper-parameters (model dimension, attention heads,
#'   transformer layers, feed-forward multiplier, dropouts, CLS token).
#' @param saint_d_model,saint_heads,saint_layers,saint_ff_mult,saint_dropout,saint_token_dropout,saint_use_cls
#'   SAINT hyper-parameters (same layout as the FT-Transformer set).
#' @param use_grouping,group_trigger,group_method,init_group_size,max_tokens,kmeans_batch,kmeans_iter
#'   Feature-tokenisation / grouping controls for the tabular-Transformer
#'   models (trigger, method, group size, max tokens, mini-batch K-means
#'   parameters).
#' @param tabnet_steps,tabnet_feature_dim,tabnet_output_dim,tabnet_gamma,tabnet_lambda_sparse
#'   TabNet hyper-parameters.
#' @param node_trees,node_depth NODE (Neural Oblivious Decision Ensembles)
#'   tree count and depth.
#' @param deepfm_k,deepfm_hidden DeepFM embedding dimension and dense hidden
#'   widths.
#' @param dcn_layers,dcn_hidden Deep & Cross Network cross / dense layers.
#' @param nam_hidden,nam_activation,nam_add_linear,nam_l1 Neural Additive
#'   Model hidden widths, activation, optional linear term and L1 penalty.
#' @param moe_n_experts,moe_expert_hidden,moe_gate_hidden,moe_temperature,moe_sparse_topk,moe_entropy_reg
#'   Mixture-of-Experts hyper-parameters.
#' @param gp_use_variational,gp_num_inducing,gp_feature_dim,gp_kernel,gp_ard,gp_lr_mult,rff_features,rff_lengthscale,rff_deep_hidden
#'   Deep-GP / random-Fourier-feature hyper-parameters for the DL-GP routes.
#' @param model_type DL model family (`"mlp"`, `"cnn"`, `"resnet"`, `"ft"`,
#'   `"saint"`, `"tabnet"`, `"node"`, `"deepfm"`, `"dcn"`, `"nam"`, `"moe"`,
#'   `"gp"`, ...).
#' @param epochs,batch_size DL training duration and mini-batch size.
#' @param dropout,l2_weight_decay,l2_regularizer_dp,dropout_rate,batch_norm,validation_split,compile_model,deterministic,random_seed,device
#'   Generic DL controls (regularisation, validation split, model compile,
#'   determinism, random seed, target device).
#' @param dl_n_seeds Optional number of DL training seeds to fit within every
#'   bootstrap sample for single-trait true prediction. The default is one,
#'   preserving the historical computational cost.
#' @param dl_seeds Optional explicit unique non-negative integer DL training seeds.
#'   When omitted, `dl_n_seeds` seeds are deterministically derived from
#'   `random_seed`; all requested seeds are retained and averaged.
#' @param dl_seed_aggregation Aggregation across DL training seeds. Currently
#'   only `"mean"` is supported; PredictProR never selects a seed using unknown
#'   test-set outcomes.
#' @param mlp_neurons_per_layer,mlp_learning_rate MLP-specific size and
#'   learning rate.
#' @param final_attention,attention_across_multiple_layers,heteroscedastic
#'   Attention pooling / heteroscedastic-output toggles for the DL wrappers.
#' @param param_grid Generic hyper-parameter grid (overrides the per-model
#'   `*_paras_tunning` lists when supplied).
#' @param early_stop Logical / control for early stopping during tuning or
#'   DL training.
#' @param k k for K-Nearest Neighbours.
#' @param learning_rate,max_depth,subsample XGBoost / boosting learning rate,
#'   tree depth and row-subsample ratio.
#' @param xgb_booster XGBoost booster (`"gbtree"`, `"gblinear"`, `"dart"`).
#' @param iteration XGBoost number of boosting iterations (`nrounds`).
#' @param N_feature_impo Top-K features to report by importance.
#' @param resample_method_tune,number_of_fold_tune Resampling method and folds
#'   for hyper-parameter tuning.
#' @param min_child_weight,colsample_bytree,xgb_alpha,xgb_gamma,lambda_rr,xgb_lambda,xgb_rate_drop,xgb_skip_drop,xgb_objective,xgb_sample_type,xgb_normalize_type,xgb_nthread
#'   XGBoost / Ridge-regression hyper-parameters (min child weight, column
#'   subsample, L1 / L2, DART drop / skip, objective, sample / normalise
#'   types, threads).
#' @param catboost_iterations,catboost_depth,catboost_learning_rate,catboost_l2_leaf_reg,catboost_thread_count
#'   CatBoost hyper-parameters.
#' @param lightgbm_nrounds,lightgbm_learning_rate,lightgbm_num_leaves,lightgbm_feature_fraction,lightgbm_bagging_fraction,lightgbm_min_data_in_leaf,lightgbm_lambda_l1,lightgbm_lambda_l2,lightgbm_nthread
#'   LightGBM hyper-parameters.
#' @param ntree,nodesize,mtry,maxnodes,importance Random Forest
#'   hyper-parameters and the variable-importance toggle. `mtry = NULL` uses
#'   R `randomForest` defaults: floor(p/3) predictors per split for
#'   regression, sqrt(p) for classification.
#' @param ncomp Number of PLS components.
#' @param svm_kernel,svm_type,sigma_value,C_value,degree_value,scale_value,gamma_value,offset_value
#'   SVM kernel family / type and kernel parameters
#'   (sigma, C, degree, scale, gamma, offset).
#' @param AI_cv_nfolds,n_bootstrap,early_stop_for_iteration_xgb Classical-ML
#'   inner-CV folds, bootstrap resamples (default 30), and XGBoost-specific
#'   early-stop rounds. In ML/DL true prediction the published prediction is
#'   the model fitted on all training lines (DL: the average of `dl_n_seeds`
#'   networks trained on all lines); the `n_bootstrap` refits only describe
#'   uncertainty and run in parallel across the available cores. Gaussian
#'   standard errors and intervals are calibrated on cross-fitted held-out
#'   errors. For Gaussian RandomForest, `PREDICTPRO_RF_UNCERTAINTY=jackknife`
#'   replaces the refits with the bias-corrected infinitesimal jackknife of the
#'   full-data forest (Wager, Hastie & Efron 2014): about 10x faster, unbiased
#'   on average with ~1000 trees but noisier than the bootstrap in simulation,
#'   so the bootstrap remains the default.
#' @param CI_width_thresholds,confidence_level Prediction-interval
#'   classification quantiles and confidence level.
#' @param high_reliability_thres,low_reliability_thres Reliability-band
#'   thresholds.
#' @param abs_very_close_threshold,abs_close_threshold Closeness thresholds
#'   used in the ranking / stability scoring.
#' @param n_components,threshold,iqr_multiplier Component / risk / outlier
#'   thresholds for the post-prediction summary.
#' @param interval_width_high_threshold,interval_width_moderate_threshold,interval_width_low_threshold
#'   Optional fixed-value overrides for the prediction-interval-width
#'   classification bands.
#' @param lowrank_eps_trace,lowrank_max_rank,lowrank_jitter,lowrank_noise_grid,lowrank_kernel_weights
#'   GP backend controls (trace tolerance, maximum rank, jitter, noise grid,
#'   and per-kernel weights). Named weights must match the aligned kernel names.
#'   They are forwarded through direct prediction, cross-validation, and
#'   automatic multi-trait final refits. They are fixed for kernel-ridge and
#'   joint multi-trait GP fits and are input/starting weights when
#'   the exact GP backend learns component scales.
#' @param gp_output_level GP backend output level
#'   (`"predict_only"`, `"predict_with_se"`, `"full_vc"`).
#' @param gp_varcomp_mode GP variance-component estimation mode
#'   (e.g. `"reml"`, `"mom"`). In joint multi-trait mode the selected model
#'   determines the statistically valid route: Gaussian-Process-GBLUP uses
#'   REML, FA-GBLUP uses factor-analytic REML (or MoM for multi-environment
#'   data), and Scalable-GBLUP uses operator/MoM.
#' @param gp_fa_rank FA rank for GP_FA / FA-GBLUP.
#' @param gp_prediction_output Which GP prediction summary to return.
#' @param gp_factor_cache Optional cached factor object reused across GP
#'   fits.
#' @param cross_validation Logical; route to the cross-validation pipeline
#'   instead of true prediction.
#' @param cv_evaluation_only Logical; only return CV evaluation outputs
#'   (skip final fit on the full data). Classification responses include
#'   family-specific metrics plus out-of-fold reliability and calibration
#'   summaries under `cv_results_processed`.
#' @param GS_model_cv Model name(s) to fit during cross-validation; defaults
#'   to `GS_model` when not supplied.
#' @param nfolds Number of folds for cross-validation.
#' @param sampling_method CV sampling scheme (e.g. `"stratified"`,
#'   `"unstratified"`).
#' @param num_cores Parallel workers for CV / bootstrap. When `NULL` (default)
#'   the policy uses `parallel::detectCores(logical = TRUE) - 1` capped at the
#'   number of parallel-eligible tasks. The full set of knobs that influence
#'   worker count and backend selection (`GP_BLAS_THREADS`, the scoring
#'   weights, mori / mirai timing, etc.) is documented in
#'   `inst/parallel-policy-knobs.md` and via the
#'   \code{\link{gp_parallel_policy_knobs}} helper.
#' @param replication Number of CV replications.
#' @param test_size Fraction reserved for the held-out test set in holdout CV.
#' @param cross_validation_meth Cross-validation method. Every multi-environment
#'   workflow uses `"CV0"`, `"CV1"`, `"CV2"`, or a repeated variant. CV0 holds
#'   out environments, CV1 holds out genotypes across environments, and CV2
#'   holds out genotype-by-environment cells. Single-environment workflows use
#'   their ordinary holdout or K-fold methods.
#' @param random_state Random seed for reproducible CV / bootstrap splits.
#' @param metric_for_ranking Metric used to rank models in the CV summary.
#'   The default, `"auto"`, uses balanced accuracy for binary responses, macro
#'   F1 for multiclass responses, quadratic-weighted kappa for ordinal
#'   responses, and correlation (`accuracy`) for Gaussian responses.
#' @param ranking_tie_breakers Optional ordered character vector of secondary
#'   metrics used when candidate models tie on `metric_for_ranking`. `NULL`
#'   uses family-aware defaults; use `"none"` to disable secondary metrics.
#' @param plot_extension,plot_width,plot_height,plot_units,plot_dpi,plot_filename,Plot_name_result_diagnostic
#'   Output-plot file extension, geometry, units, DPI, base filename and the
#'   diagnostic-plot name prefix.
#' @param feature_selected Optional exact predictor set to reuse without
#'   rescoring. Supply a character vector for one common set; a list keyed by
#'   trait; a list keyed by `geno_data`, `omic1_data`, `omic2_data`, or
#'   `omic3_data`; a nested trait/source list in either orientation; or a data
#'   frame containing `predictor` plus optional `trait`, `source_block`, and
#'   logical `selected` columns. Sets are applied separately to each trait;
#'   genuinely joint multi-trait fits use their union. In a source-keyed set,
#'   an omitted source is excluded. Named raw matrices allow PredictProR to
#'   rebuild marker designs and kernels. With a user-supplied precomputed
#'   kernel, the package cannot inspect the originating columns and assumes the
#'   kernel was externally rebuilt from this exact set.
#' @param feature_scoring_seed Random seed for the feature-scoring pipeline.
#' @param feature_ridge_lambda Ridge penalty used by the feature-scoring
#'   ridge engine.
#' @param feature_bayes_nIter,feature_bayes_burnIn,feature_bayes_thin BGLR
#'   MCMC controls used by the Bayesian feature-scoring engine.
#' @param globals_max_GB Maximum exported-globals size (GB) used by the
#'   parallel-policy scorer.
#' @param worker_memory_gb Optional memory (GiB) one parallel worker needs on
#'   top of the shared data: the R process with PredictProR loaded plus, for
#'   Python-backed ML/DL/GP models, its Python child. Together with
#'   `memory_budget_gb` it caps how many workers run at once. `NULL` (default)
#'   uses `GP_PAR_WORKER_OVERHEAD_GB` if set, otherwise measures one worker on
#'   this machine once per session (only when more than one worker could run),
#'   falling back to 0.75 GiB. The value and its source are recorded in
#'   `Run_metadata` (`policy_worker_memory_gb`, `policy_worker_memory_source`).
#' @param memory_budget_gb Optional total memory (GiB) all parallel workers may
#'   use. `NULL` (default) uses `GP_PAR_MEMORY_BUDGET_GB` if set, otherwise
#'   `GP_PAR_MEMORY_FRACTION` (default 0.7) of the currently available memory.
#' @param parallel_mode Backend selection. One of `"auto"` (default;
#'   score-based — the engine picks among mirai / future / base_parallel /
#'   foreach / sequential based on object size, fanout, OS, GPU
#'   availability, and per-model gates), or one of `"mirai"`, `"future"`,
#'   `"base_parallel"`, `"foreach"`, `"sequential"` to force a backend.
#'   The chosen backend and the contributing reason codes are recorded
#'   in `policy_decision_reason` in the exported `Run_metadata.csv`;
#'   decode them with \code{\link{gp_policy_decision_reason_glossary}}.
#'   Tuning is via the env vars / options enumerated by
#'   \code{\link{gp_parallel_policy_knobs}} and described in
#'   `inst/parallel-policy-knobs.md`.
#' @param parallel_backend_prefer_fork Logical; prefer forking on Unix
#'   workers when applicable. Set `FALSE` for PSOCK workers. Forking is
#'   disabled automatically when a worker needs embedded-runtime
#'   initialization, such as a reticulate Python session.
#' @param sequential_models Optional character vector of model names that
#'   should always run sequentially (e.g. RKHS, XGBoost, CatBoost / LightGBM
#'   which already use internal threading).
#' @param verbose Logical; print verbose progress / diagnostic messages.
#'
#' @return A nested list of per-trait model results. The prediction table
#'   (`model_results$predicted_values`) follows a shared column contract that
#'   includes `Predicted_value`, `Standard_error`, `PEV`,
#'   `lower_bound`/`upper_bound`, and `Reliability`, plus model-specific
#'   variance metadata.
#'
#'   For genetic mixed models, inspect `Reliability_basis`: the reference
#'   variance is an estimated genetic variance and PEV has breeding-value
#'   semantics. ML/DL models do not estimate that decomposition, so
#'   `Genetic_variance` remains unavailable. Their plotting-compatible
#'   `Reliability` field equals `Prediction_stability` and is calculated as
#'   `max(0, min(1, 1 - U_i / V_yhat))`, where `V_yhat` is the fitted-value
#'   variance across all rows in the final prediction table (both `Train` and
#'   `Test`, when present) and `U_i` is normally prediction-resampling `SE_i^2`.
#'   `Train_Test_Label` identifies the row set but does not subset `V_yhat`.
#'   The two inputs are returned in
#'   `Reliability_reference_variance` and `Reliability_variance_input`.
#'   Held-out predictive PEV and prediction intervals remain separate.
#'   ML/DL `Reliability` is an uncalibrated marker-adjustment plotting
#'   surrogate, not quantitative-genetic reliability.
#'
#'   Multi-seed DL true prediction additionally returns `dl_seed_manifest`,
#'   `dl_computation_plan`, `dl_seed_predictions`, and
#'   `dl_seed_variability`. These remain separate from bootstrap PEV and
#'   reliability because they describe optimizer/initialization sensitivity.
#'
#'   With `multi_trait = TRUE`, the returned
#'   `PredictProR_multi_trait_result` keeps each protected child result under
#'   `model_results_by_model`, the routing decision under
#'   `orchestration_plan`, combined CV summaries and rankings under
#'   `cv_results_processed`, and any selected full-data refit under
#'   `final_prediction`.
#' @export
#'
#' @examples
#'
model_execute <- function(
    pheno_data = NULL,
    # pheno_file_name = NULL,
    # pheno_file_path= NULL,
    pheno_data_train = NULL,
    pheno_data_test = NULL,
    geno_data = NULL,
    omic1_data = NULL,
    omic2_data = NULL,
    omic3_data = NULL,
    omics_data_label = list(omic1_data = NULL,
                            omic2_data = NULL,
                            omic3_data = NULL),
    gmatrix= NULL,
    gkernel = NULL,
    kernel_list = NULL,
    pedigree_matrix = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    omics_kernel_label = list(omic1_kernel = NULL,
                              omic2_kernel = NULL,
                              omic3_kernel = NULL),
    train_geno_data = NULL,
    train_omic1_data = NULL,
    train_omic2_data = NULL,
    train_omic3_data = NULL,
    train_omics_label = list(train_omic1_data = NULL,
                             train_omic2_data = NULL,
                             train_omic3_data = NULL),
    test_geno_data = NULL,
    test_omic1_data = NULL,
    test_omic2_data = NULL,
    test_omic3_data = NULL,
    test_omics_label = list(test_omic1_data = NULL,
                            test_omic2_data = NULL,
                            test_omic3_data = NULL),
    coefficient_1 = NULL,
    coefficient_2 = NULL,
    coefficient_3 = NULL,
    coefficient_4 = NULL,
    train_set = NULL,
    test_set = NULL,
    gmatrix_method = NULL,
    kernel_method = NULL,
    response=NULL,
    response_family = "auto",
    positive_class = NULL,
    multi_trait_gp = FALSE,
    multi_trait_bayes = FALSE,
    multi_trait_asreml = FALSE,
    hybrid_asreml = FALSE,
    hybrid_bayes = FALSE,
    hybrid_gp = FALSE,
    hybrid_ml = FALSE,
    hybrid_dl = FALSE,
    female_parent = NULL,
    male_parent = NULL,
    female_gmatrix = NULL,
    male_gmatrix = NULL,
    female_geno_data = NULL,
    male_geno_data = NULL,
    hybrid_include_sca = TRUE,
    hybrid_gp_lambda = "auto",
    hybrid_gp_lambda_grid = c(0.01, 0.03, 0.1, 0.3, 1),
    hybrid_gp_component_weights = NULL,
    env_similarity = NULL,
    env_ids = NULL,
    env_covariates = NULL,
    reaction_norm_feature_qc = TRUE,
    kenv_kernel = "matern32",
    kenv_bandwidth = 1.0,
    kenv_kernel_kwargs = NULL,
    multi_trait_ml = FALSE,
    multi_trait_dl = FALSE,
    met_ml_dl = FALSE,
    met_kernel_var_explained = 0.95,
    met_kernel_min_ev = 1e-8,
    met_kernel_max_pcs = NULL,
    gen_name=NULL,
    ld_prunning_qc = TRUE,
    docker_nd_usage = FALSE,
    cova=NULL,
    fixed=NULL,
    random=NULL,
    heter_resid=FALSE,
    bayes_kernel_heter_resid = NULL,
    heter_groups=NULL,
    var_cov_str = NULL,
    weights =NULL,
    nIter=NULL,
    burnIn=NULL,
    thin=NULL,
    GS_model = NULL,
    eval_metrics = NULL,
    fixed_term_model_bayesian = 'FIXED',
    rand_term_model_bayesian = NULL,
    selected = NULL, # nd_mods
    max_features = 50,
    #core = NULL,
    engine = NULL,
    scaling = TRUE,
    centering = FALSE,
    workspace = 1e08,
    pworkspace= 1e06,
    maxit = 50,
    inverse = TRUE,
    epsilon = 1e-6,
    bending = TRUE,
    bend_value = 0.01,
    blending = FALSE,
    blending_value = 0.02,
    high_diag_cut_off = 1.2,
    low_diag_cut_off = 0.8,
    duplicate_cut_off = 0.95,
    rcn_cutoff = 1e-12,
    optimize_diagonal = FALSE,
    optimize_duplicate = FALSE,
    kernel_check_level = "auto",
    kernel_large_n_threshold = 5000L,
    duplicate_scan = "auto",
    duplicate_sample_size = 2000L,
    duplicate_block_size = 1024L,
    duplicate_max_pairs = 10000L,
    kernel_fix_method = "auto",
    kernel_repair_priority = "speed",
    kernel_rcn_check = "auto",
    kernel_nearpd_size_limit = 2500L,
    kernel_cpp_repair_size_limit = 3000L,
    kernel_cpp_keep_diag = TRUE,
    kernel_pd_check = "auto",
    kernel_pd_sample_size = 500L,
    kernel_sanitize = "auto",
    kernel_sanitize_value = NULL,
    vcf_file_name = NULL,
    vcf_file_path = NULL,
    vcf_file = NULL,
    hapmap_file_name = NULL,
    hapmap_file_path = NULL,
    hapmap = NULL,
    csv_file_name = NULL,
    csv_file_path = NULL,
    csv_input_coding = c("alt_dosage", "centered_dosage"),
    met_predict_all_environments = TRUE,
    maf_threshold = 0.01,
    het_threshold = 0.1,
    ind_call_rate_threshold = 0.9,
    snp_call_rate_threshold = 0.9,
    test_train_genetic_space = FALSE,
    impute = TRUE,
    impute_omic = TRUE,
    imputation_method = "knn",
    beagle_options = list(),
    impute_knn_k = 5,
    ploidy = "auto",
    na_threshold = 0.9,
    recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
    ld_pruning = FALSE,         # LD pruning option
    ld_pruning_method = "indep-pairwise", # LD pruning method
    window_size = 50,           # Window size for LD pruning
    step_size = 5,              # Step size for LD pruning
    r2_threshold = 0.2,         # r^2 threshold for LD pruning
    use_kb_window = TRUE,      # Use kb for window size in LD pruning
    phased = TRUE,             # Option for phased LD pruning
    use_founders = FALSE,
    out_put_map = FALSE,
    map_data = NULL,
    qc_filtering = TRUE,
    message= TRUE,
    system_database = FALSE,
    #### optimizer / AMP / class weights / grad-clip ---
    optimizer_name = "adam",
    use_amp        = TRUE,
    max_grad_norm  = 1.0,
    auto_class_weights = FALSE,
    internal_cv_nfolds = 5L,
    internal_cv_replication = 3L,
    dl_internal_calibration = TRUE,
    ##### cnn
    cnn_neurons_per_layer = as.integer(c(64, 64, 64)),
    cnn_kernel_size = 3L,
    cnn_dense_layers = as.integer(c(256, 128, 64)),
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel = 2L,
    cnn_pool_stride = 2L,
    cnn_pool_padding = 0L,
    cnn_learning_rate = 1e-3,
    cnn_separable=TRUE,
    cnn_dilations=c(1,2,4),
    cnn_use_se=TRUE,
    cnn_norm_type="group",
    cnn_pool_type="conv",
    cnn_use_global_pool=FALSE,
    ##### resnet
    resnet_neurons_per_block = as.integer(c(256, 128, 64)),
    resnet_blocks = 3,
    resnet_learning_rate = 1e-3,
    #### ft_transformer
    ft_d_model = 192L,
    ft_heads = 8L,
    ft_layers = 3L,
    ft_ff_mult = 4L,
    ft_dropout = 0.1,
    ft_token_dropout = 0.0,
    ft_use_cls = TRUE,
    #### saint
    saint_d_model = 128L,
    saint_heads = 8L,
    saint_layers = 3L,
    saint_ff_mult = 4L,
    saint_dropout = 0.1,
    saint_token_dropout = 0.0,
    saint_use_cls       = TRUE,
    ###### Grouping controls (FT/SAINT and NAM/MoE)
    use_grouping   = FALSE,
    group_trigger  = 2048,
    group_method   = "auto",
    init_group_size = 64,
    max_tokens      = 1024,
    kmeans_batch    = 4096,
    kmeans_iter     = 100,
    #### tabnet
    tabnet_steps = 5L,
    tabnet_feature_dim = 64L,
    tabnet_output_dim = 64L,
    tabnet_gamma = 1.5,
    tabnet_lambda_sparse = 1e-4,
    #### node
    node_trees = 8L,
    node_depth = 3L,
    #### deepfm
    deepfm_k = 16L,
    deepfm_hidden = as.integer(c(128, 64)),
    #### dcnv2
    dcn_layers = 3L,
    dcn_hidden = as.integer(c(256, 128, 64)),
    #### nam
    nam_hidden = as.integer(c(32, 16)),
    nam_activation = "relu",
    nam_add_linear = TRUE,
    nam_l1 = 1e-4,
    ### moe
    moe_n_experts = 4L,
    moe_expert_hidden = as.integer(c(128, 64)),
    moe_gate_hidden = 128L,
    moe_temperature = 1.0,
    moe_sparse_topk = NA,
    moe_entropy_reg = 0.0,
    ### gp_dkl/RFF knobs
    gp_use_variational = TRUE,
    gp_num_inducing = 256L,
    gp_feature_dim = 64L,
    gp_kernel = "rbf",
    gp_ard = TRUE,
    gp_lr_mult = 0.5,
    rff_features       = 1024,
    rff_lengthscale    = 1.0,
    rff_deep_hidden    = c(128),
    ##### General dp
    model_type = "resnet",
    epochs = 10,
    batch_size = 64 ,
    dropout = 0.2,
    l2_weight_decay = 1e-4,
    l2_regularizer_dp = 0.001,
    dropout_rate = 0.5,
    batch_norm = TRUE,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 123,
    dl_n_seeds = NULL,
    dl_seeds = NULL,
    dl_seed_aggregation = "mean",
    device = NULL,
    #### mlp and attention
    mlp_neurons_per_layer = as.integer(c(128, 64)),
    mlp_learning_rate = 1e-3,
    final_attention = TRUE,
    attention_across_multiple_layers = TRUE,
    heteroscedastic = TRUE,
    ###
    #attention_on_final_layer = TRUE,
    #attention_across_multiple_layers = FALSE,
    #batch_normalization = TRUE,
    #deep_learning_model = "mlp_with_attention", #"mlp", "ResNet", "cnn"
    para_tunning = FALSE,
    param_grid = NULL,
    early_stop = TRUE,
    xgb_paras_tunning= list(Iter_tune = seq(100, 500, 100), # number of boosting iterations
                            learning_rate_tune = c(0.01, 0.1, 0.1), # learning rate, low value means model is more robust to overfitting
                            max_depth = c(3, 6, 9),
                            rate_drop = c(0.1, 0.15, 0.2),
                            skip_drop = c(0.4, 0.5, 0.55),
                            xgb_gamma = c(0, 0.01, 0.1),
                            colsample_bytree = c(0.5, 0.75, 1),
                            min_child_weight = c(1, 3, 5),
                            subsample = c(0.5, 0.75, 1),
                            L2_tune = c(0, 0.5, 1), #  for linear gbL2 Regularization (Ridge Regression)
                            L1_tune = c(0, 0.5, 1)),
    rf_paras_tunning= list(mtry = TRUE,
                           ntree = c(500, 1000, 1500),
                           nodesize = c(1, 5, 10),
                           maxnodes = c(30, 50, NULL)),  # NULL means no limit),
    pls_paras_tunning= list(ncomp = 10),
    svm_paras_tunning=list(
      kernel = c("Gaussian", "Linear",
                 "Polynomial", "Hyperbolic_tangent"),
      #cost = 10^seq(-2, 2, by = 1),
      offset_value = seq(-2, 2, length.out = 5),
      sigma = c(0.01, 0.05, 0.1),
      C = c(1, 10, 100), ## for radial kernel
      gamma_value = 10^seq(-4, -1, length.out = 4),
      degree = c(3, 4),  # Default values, used only for polynomial
      scale = c(0.1, 1) # used only for polynomial
    ),
    knn_paras_tunning= list(k = seq(3, 21, by = 2)),
    k = 5,
    lasso_paras_tunning= list(lambda_tune=seq(0.000001,0.9,length.out=100)^4),
    rr_paras_tunning = NULL,
    dpl_paras_tunning = NULL,
    learning_rate = 0.01, #xgboost
    max_depth = 6, #xgboost
    subsample = 0.7, #xgboost
    xgb_booster = "gbtree", # xgboost: "gbtree" (default), "dart" (much slower), "gblinear"
    iteration = 100, #xgboost
    N_feature_impo = 10, #xgboost
    resample_method_tune = "cv", # c("cv","boot") #xgboost
    number_of_fold_tune = 5, #xgboost
    min_child_weight = 0.8, # xgboost,
    #eta = 0.001, ## xgboost
    #nrounds = 5000, ## xgboost
    colsample_bytree = 0.7, ## xgboost
    xgb_alpha = 0.001, ## xgboost linear
    xgb_gamma = 0.01, ## xgboost it acts as a regularization parameter for controlling tree complexity
    lambda_rr = NULL,
    xgb_lambda = 1.0,  # xgboost linear
    xgb_rate_drop = 0.1,
    xgb_skip_drop = 0.5,
    xgb_objective = "reg:squarederror",
    xgb_sample_type = "uniform",
    xgb_normalize_type = "tree",
    xgb_nthread = 1L,
    catboost_iterations = 500,
    catboost_depth = 6,
    catboost_learning_rate = 0.03,
    catboost_l2_leaf_reg = 3,
    catboost_thread_count = 1L,
    lightgbm_nrounds = 100,
    lightgbm_learning_rate = 0.05,
    lightgbm_num_leaves = 31,
    lightgbm_feature_fraction = 1.0,
    lightgbm_bagging_fraction = 1.0,
    lightgbm_min_data_in_leaf = 20,
    lightgbm_lambda_l1 = 0,
    lightgbm_lambda_l2 = 0,
    lightgbm_nthread = 1L,
    ntree=500, ## RF
    nodesize =NULL,
    mtry = NULL, ## RF
    maxnodes = NULL, ## RF
    importance=TRUE, ## RF
    rf_n_jobs = 1L,
    ncomp = 3, # pls
    svm_kernel = "Gaussian", #svm "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
    svm_type = "eps-regression",
    sigma_value  = NULL,      #svm: NULL uses dimension-aware gamma = "scale"
    C_value  = 1,             #svm Default cost parameter
    degree_value = 3,        #svm Default degree for polynomial kernel
    scale_value  = 1,         #svm Default scale for polynomial kernel
    gamma_value = NULL,
    offset_value = 0,
    AI_cv_nfolds = 5,
    n_bootstrap = 30,
    early_stop_for_iteration_xgb = FALSE,
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5,
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05,
    n_components = 20,
    threshold = 100,
    iqr_multiplier = 1.5,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    interval_width_low_threshold = NULL,
    lowrank_eps_trace = 1e-6,
    lowrank_max_rank = NULL,
    lowrank_jitter = NULL,
    lowrank_noise_grid = NULL,
    lowrank_kernel_weights = NULL,
    gp_backend = "auto",
    gp_output_level = "predict_only",
    gp_return_se = FALSE,
    gp_full_vc = FALSE,
    gp_return_trait_correlations = FALSE,
    gp_estimate_kernel_weights = FALSE,
    gp_force_prediction_se = FALSE,
    gp_varcomp_mode = "reml",
    gp_fa_rank = 1L,
    gp_prediction_output = "all",
    gp_factor_cache = NULL,
    gp_iters = NULL,
    gp_lr = NULL,
    gp_engine = "auto",
    gp_learn_scales = NULL,
    gp_exact_fast_cv = NULL,
    cross_validation = FALSE,
    cv_evaluation_only = FALSE,
    cv_generate_plots = NULL,
    GS_model_cv = NULL,
    nfolds = 5,
    sampling_method = NULL,
    num_cores = NULL,
    replication = 1,
    test_size = 0.3,
    cross_validation_meth = "Stratified_Hold_Out",
    random_state = 123,
    metric_for_ranking = "auto",
    ranking_tie_breakers = NULL,
    plot_extension = "pdf",
    plot_width = 17,
    plot_height = 12,
    plot_units = "in",
    plot_dpi = 300,
    plot_filename = "trait",
    Plot_name_result_diagnostic = NULL,
    feature_selected = NULL,
    feature_scoring = FALSE,
    feature_scoring_model = "Ridge_Regression",
    feature_k_grid = NULL,
    feature_k = NULL,
    feature_scoring_cv = "fixed",
    feature_score_metadata = NULL,
    feature_scoring_seed = NULL,
    feature_ridge_lambda = 1,
    feature_bayes_nIter = 1500L,
    feature_bayes_burnIn = 500L,
    feature_bayes_thin = 5L,
    ######
    globals_max_GB = 4,
    worker_memory_gb = NULL,
    memory_budget_gb = NULL,
    parallel_mode = c("auto","future","sequential","base_parallel","foreach","mirai"),
    parallel_backend_prefer_fork = TRUE,
    sequential_models = NULL,
    verbose = TRUE,
    multi_trait = FALSE,
    ...
) {

#browser()
    # once-per-call warnings (e.g. variance components not estimable)
    assign("warned", character(), envir = PredictProR_runtime_cache)
    gp_reject_obsolete_asreml_structure(var_cov_str)
    msg <- ""
    on.exit(future::plan("sequential"), add = TRUE)
    parallel_mode <- match.arg(parallel_mode)
    dl_internal_calibration <- gp_dl_calibration_enabled(dl_internal_calibration)
    for (mem_arg in c("worker_memory_gb", "memory_budget_gb")) {
      mem_val <- get(mem_arg)
      if (!is.null(mem_val) &&
          (!is.numeric(mem_val) || length(mem_val) != 1L || !is.finite(mem_val) || mem_val <= 0)) {
        stop(mem_arg, " must be NULL or one positive number of GiB.", call. = FALSE)
      }
    }
    # Hybrid ML/DL features are hybrid-level genotypes, which are heterozygous
    # by design; the inbred-line heterozygosity filter would discard most
    # informative markers. With parent genotypes, the filter runs on the
    # inbred parents and the failing markers are dropped from the hybrid
    # features; otherwise it stays off unless the user sets it explicitly.
    hybrid_parent_qc_summary <- NULL
    if (isTRUE(hybrid_ml) || isTRUE(hybrid_dl)) {
      if (!is.null(female_geno_data) || !is.null(male_geno_data)) {
        hybrid_parent_prep <- gp_hybrid_ml_parent_geno_prepare(
          pheno_data = pheno_data,
          geno_data = geno_data,
          female_geno_data = female_geno_data,
          male_geno_data = male_geno_data,
          gen_name = gen_name,
          female_parent = female_parent,
          male_parent = male_parent,
          het_threshold = het_threshold,
          ploidy = ploidy,
          message = message
        )
        geno_data <- hybrid_parent_prep$geno_data
        hybrid_parent_qc_summary <- hybrid_parent_prep$summary
        female_geno_data <- NULL
        male_geno_data <- NULL
        het_threshold <- NULL
      } else if (missing(het_threshold)) {
        het_threshold <- NULL
      }
    }


    `%||%` <- function(a, b) if (is.null(a)) b else a
    dot_args <- list(...)
    predictpror_preprocess_cache <-
      dot_args[[".predictpror_preprocess_cache"]] %||% NULL
    predictpror_shared_models <-
      dot_args[[".predictpror_shared_models"]] %||% character()
    predictpror_shared_requires_asreml <- isTRUE(
      dot_args[[".predictpror_shared_requires_asreml"]]
    )
    predictpror_shared_final_prediction <- isTRUE(
      dot_args[[".predictpror_shared_final_prediction"]]
    )
    gp_shared_preprocess_cache_validate(predictpror_preprocess_cache)
    if (isTRUE(multi_trait)) {
      formal_names <- setdiff(names(formals(model_execute)), "...")
      orchestration_args <- mget(
        formal_names,
        envir = environment(),
        inherits = FALSE
      )
      if (length(dot_args)) {
        if (is.null(names(dot_args)) || any(!nzchar(names(dot_args)))) {
          stop(
            "All arguments passed through `...` must be named when multi_trait = TRUE.",
            call. = FALSE
          )
        }
        orchestration_args[names(dot_args)] <- dot_args
      }
      return(gp_model_execute_multitrait_orchestrate(orchestration_args))
    }
    ######
    raw_preprocess_stage <- "raw_genotype_import"
    if (gp_shared_preprocess_cache_has(
      predictpror_preprocess_cache,
      raw_preprocess_stage
    )) {
      raw_preprocess <- gp_shared_preprocess_cache_get(
        predictpror_preprocess_cache,
        raw_preprocess_stage
      )
      beagle_preprocess <- raw_preprocess$beagle_preprocess
      geno_data_process <- raw_preprocess$geno_data_process
      ploidy <- raw_preprocess$ploidy
      imputation_method <- raw_preprocess$imputation_method
    } else {
      beagle_preprocess <- NULL
      # CSV / TXT genotype tables are converted to VCF once and then take the
      # VCF route (QC, Beagle or native imputation, recoding).
      if (!is.null(csv_file_name)) {
        if (!is.null(vcf_file_name) || !is.null(hapmap_file_name) || !is.null(hapmap)) {
          stop("Give one raw genotype source: csv_file_name, vcf_file_name or a HapMap input.", call. = FALSE)
        }
        csv_input <- file.path(csv_file_path %||% getwd(), csv_file_name)
        if (!file.exists(csv_input)) {
          stop("Genotype table not found: ", csv_input, call. = FALSE)
        }
        csv_coding <- match.arg(csv_input_coding, c("alt_dosage", "centered_dosage"))
        if (!is.numeric(ploidy)) {
          if (identical(tolower(trimws(as.character(imputation_method[[1L]]))), "beagle")) {
            gp_beagle_check_csv_dosage(csv_input, csv_coding)
          } else if (tryCatch({ gp_beagle_check_csv_dosage(csv_input, csv_coding); FALSE },
                              error = function(e) TRUE)) {
            stop("The genotype table has dosages above 2 (polyploid). Give the ploidy, e.g. ploidy = 4L.",
                 call. = FALSE)
          }
        }
        csv_vcf <- file.path(tempfile("predictpror-csv-"), "genotypes_from_table.vcf")
        dir.create(dirname(csv_vcf), recursive = TRUE, showWarnings = FALSE)
        convert_csv_to_vcf(
          input_csv = csv_input,
          output_vcf = csv_vcf,
          ploidy = if (is.numeric(ploidy)) as.integer(ploidy)[1L] else 2L,
          input_coding = match.arg(csv_input_coding, c("alt_dosage", "centered_dosage"))
        )
        vcf_file_name <- basename(csv_vcf)
        vcf_file_path <- dirname(csv_vcf)
      }
      requested_imputation_method <- tolower(trimws(as.character(imputation_method[[1L]])))
      if (identical(requested_imputation_method, "beagle")) {
      if (!isTRUE(impute)) {
        stop("imputation_method = 'beagle' requires impute = TRUE.", call. = FALSE)
      }
      # Beagle is diploid-only: stop early for polyploid data and point to KNN
      requested_ploidy <- if (is.numeric(ploidy)) as.integer(ploidy)[1L] else
        suppressWarnings(as.integer(attr(geno_data, "ploidy"))[1L])
      if (isTRUE(!is.na(requested_ploidy) && requested_ploidy > 2L)) {
        gp_beagle_polyploid_stop(requested_ploidy)
      }
      if (!is.list(beagle_options) || (length(beagle_options) && is.null(names(beagle_options)))) {
        stop("beagle_options must be a named list.", call. = FALSE)
      }
      raw_sources <- c(
        vcf = !is.null(vcf_file_name),
        hapmap_file = !is.null(hapmap_file_name),
        hapmap_table = !is.null(hapmap)
      )
      if (sum(raw_sources) != 1L) {
        stop(
          "Beagle imputation requires exactly one raw genotype source: VCF, HapMap file, or in-memory HapMap table.",
          if (!is.null(geno_data)) " For a numeric dosage matrix (geno_data) use imputation_method = \"knn\" or \"mean\"." else "",
          call. = FALSE
        )
      }
      if (isTRUE(raw_sources[["vcf"]])) {
        input_path <- file.path(vcf_file_path %||% getwd(), vcf_file_name)
        input_format <- "vcf"
      } else if (isTRUE(raw_sources[["hapmap_file"]])) {
        input_path <- file.path(hapmap_file_path %||% getwd(), hapmap_file_name)
        input_format <- "hapmap"
      } else {
        input_path <- hapmap
        input_format <- "hapmap"
      }
      matrix_fallback_method <- beagle_options$matrix_fallback_method %||% "mean"
      beagle_options$matrix_fallback_method <- NULL
      beagle_call <- utils::modifyList(
        list(
          input = input_path,
          input_format = input_format,
          # Beagle's VCF, logs, converted input and map are intermediate files
          # here; keep them out of the working directory unless the user sets
          # beagle_options$output_prefix.
          output_prefix = file.path(tempfile("predictpror-beagle-"), "beagle_imputed"),
          ploidy = ploidy,
          return_genotypes = FALSE,
          # call-rate QC must see the observed calls, i.e. run before Beagle
          snp_call_rate_threshold = snp_call_rate_threshold,
          ind_call_rate_threshold = ind_call_rate_threshold
        ),
        beagle_options
      )
      beagle_preprocess <- do.call(impute_genotypes_with_beagle, beagle_call)
      vcf_file_name <- basename(beagle_preprocess$output_vcf)
      vcf_file_path <- dirname(beagle_preprocess$output_vcf)
      hapmap_file_name <- NULL
      hapmap_file_path <- NULL
      hapmap <- NULL
      imputation_method <- matrix_fallback_method
      }
      geno_data_process <- NULL
      # qc_filtering = FALSE must switch the file recoders' filters off too (NULL
      # thresholds), as it does for matrix input.
      file_qc <- !isFALSE(qc_filtering)
      if (!is.null(vcf_file_name) && isTRUE(ld_pruning)) {
        stop(
          "ld_pruning = TRUE needs the PLINK QC engine, which model_execute() does not run. ",
          "Use ld_prunning_qc = TRUE (default: built-in LD pruning of the recoded matrix), or ",
          "pre-process with vcf_qc_recode(..., qc_engine = 'plink', ld_pruning = TRUE) and pass ",
          "the resulting snps_matrix as geno_data.",
          call. = FALSE
        )
      }
      if(!is.null(vcf_file_name)){
      vcf_file_path <- vcf_file_path %||% getwd()
      geno_data_process <- vcf_qc_recode(vcf_file_name = vcf_file_name,
                                 vcf_file_path = vcf_file_path,
                                 maf_threshold = if (file_qc) maf_threshold else NULL,
                                 het_threshold = if (file_qc) het_threshold else NULL,
                                 ind_call_rate_threshold = if (file_qc) ind_call_rate_threshold else NULL,
                                 snp_call_rate_threshold = if (file_qc) snp_call_rate_threshold else NULL,
                                 impute = if (is.null(beagle_preprocess)) impute else FALSE,
                                 imputation_method = imputation_method,
                                 impute_knn_k = impute_knn_k,
                                 ploidy = ploidy,
                                 recode_format = recode_format,
                                 ld_pruning = ld_pruning,         # LD pruning option
                                 ld_pruning_method = ld_pruning_method, # LD pruning method
                                 window_size = window_size,           # Window size for LD pruning
                                 step_size = step_size,              # Step size for LD pruning
                                 r2_threshold = r2_threshold,         # r^2 threshold for LD pruning
                                 use_kb_window = use_kb_window,      # Use kb for window size in LD pruning
                                 phased = TRUE,
                                 out_put_map = out_put_map,
                                 message = message)

      } else {
        if(!is.null(hapmap_file_name) || !is.null(hapmap)){
        geno_data_process <- hmp_qc_recode(hapmap_file_name = hapmap_file_name,
                                   hapmap_file_path = hapmap_file_path,
                                   hapmap = hapmap,
                                   maf_threshold = if (file_qc) maf_threshold else NULL,
                                   het_threshold = if (file_qc) het_threshold else NULL,
                                   ind_call_rate_threshold = if (file_qc) ind_call_rate_threshold else NULL,
                                   snp_call_rate_threshold = if (file_qc) snp_call_rate_threshold else NULL,
                                   impute = if (is.null(beagle_preprocess)) impute else FALSE,
                                   imputation_method = imputation_method,
                                   impute_knn_k = impute_knn_k,
                                   ploidy = ploidy,
                                   recode_format = recode_format,
                                   out_put_map = out_put_map,
                                   message = message)
        }


      }

      if (!is.null(geno_data_process) && !is.null(beagle_preprocess)) {
        geno_data_process$beagle <- beagle_preprocess
        attr(geno_data_process$snps_matrix, "beagle_provenance") <- beagle_preprocess$beagle
      }
      if (!is.null(geno_data_process$ploidy)) {
        ploidy <- geno_data_process$ploidy
      }
      gp_shared_preprocess_cache_set(
        predictpror_preprocess_cache,
        raw_preprocess_stage,
        list(
          beagle_preprocess = beagle_preprocess,
          geno_data_process = geno_data_process,
          ploidy = ploidy,
          imputation_method = imputation_method
        )
      )
    }

    # Ploidy is now final (argument, matrix attribute, or read from the file)
    gp_reject_polyploid_unsupported_models(
      models = list(GS_model, GS_model_cv),
      ploidy = if (is.numeric(ploidy)) ploidy else
        attr(geno_data, "ploidy") %||% geno_data_process[["ploidy"]]
    )

    if(isTRUE(cross_validation) && !is.null(GS_model)) GS_model <- NULL
    if (!is.null(eval_metrics)) {
      early_eval_metric_family <- tryCatch(
        gp_resolve_response_family_for_responses(
          pheno_data = pheno_data %||% pheno_data_train %||% pheno_data_test,
          response = response,
          response_family = response_family
        ),
        error = function(e) gp_normalize_response_family(response_family)
      )
      if (!identical(early_eval_metric_family, "auto")) {
        gp_validate_eval_metrics(eval_metrics = eval_metrics, response_family = early_eval_metric_family)
      }
    }

    # Define available models and variance structures
    var_cov_str_available <- c("us","corgh",
                               "corh","corv","fa1",
                               "fa2", "fa3", "fa4",
                               "rr1","rr2", "rr3", "rr4")

    AI_valid_models <- gp_classical_ml_supported_models()

    bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
    bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
    gp_valid_models <- gp_lowrank_supported_models()
    canonical_names <- gp_deep_learning_supported_models()
    friendly_names <- gp_deep_learning_friendly_models()
    GS_model_cv_display <- gp_display_supported_model_names(GS_model_cv)
    GS_model_display <- gp_display_supported_model_names(GS_model)
    sequential_models <- gp_canonicalize_supported_model_names(sequential_models)
    GS_model <- gp_canonicalize_supported_model_names(GS_model)
    GS_model_cv <- gp_canonicalize_supported_model_names(GS_model_cv)
    gp_validate_single_trait_gp_model_scope(
      GS_model = GS_model,
      GS_model_cv = GS_model_cv,
      multi_trait_gp = multi_trait_gp,
      hybrid_gp = hybrid_gp
    )
    AI_valid_models <- unique(c(AI_valid_models, gp_deep_learning_supported_models()))
    heter_control <- gp_normalize_single_environment_heter_controls(
      pheno_data = pheno_data,
      gen_name = gen_name,
      heter_groups = heter_groups,
      heter_resid = heter_resid,
      var_cov_str = var_cov_str,
      response = response,
      response_family = response_family,
      preserve_heter_groups = isTRUE(hybrid_gp),
      preserve_var_cov_str = isTRUE(multi_trait_asreml)
    )
    heter_groups <- heter_control$heter_groups
    heter_resid <- heter_control$heter_resid
    var_cov_str <- heter_control$var_cov_str
    if (isTRUE(heter_control$changed) && isTRUE(message)) {
      base::message(
        paste(
          msg,
          "heter_groups, heter_resid, and var_cov_str were ignored because the current data are single-environment for this workflow."
        )
      )
    }
    gp_public_multitrait_run <- isTRUE(multi_trait_gp) &&
      !isTRUE(cross_validation) &&
      length(response) > 1L &&
      !is.null(GS_model) &&
      length(GS_model) == 1L &&
      GS_model %in% gp_valid_models &&
      !isTRUE(multi_trait_asreml) &&
      !isTRUE(multi_trait_bayes) &&
      !isTRUE(multi_trait_ml) &&
      !isTRUE(multi_trait_dl) &&
      !isTRUE(hybrid_asreml) &&
      !isTRUE(hybrid_bayes) &&
      !isTRUE(hybrid_gp) &&
      !isTRUE(hybrid_ml) &&
      !isTRUE(hybrid_dl) &&
      !isTRUE(met_ml_dl) &&
      !isTRUE(feature_scoring) &&
      is.null(feature_score_metadata) &&
      is.null(feature_k_grid) &&
      is.null(feature_k) &&
      is.null(feature_selected)
    if (isTRUE(gp_return_se) || isTRUE(gp_full_vc)) {
      gp_output_level <- gp_bridge_output_level(
        return_se = gp_return_se,
        full_vc = gp_full_vc
      )
    }
    if (isTRUE(cross_validation)) {
      gp_output_level <- "predict_only"
      gp_return_se <- FALSE
      gp_full_vc <- FALSE
      gp_return_trait_correlations <- FALSE
      gp_prediction_output <- "test_only"
    }
    if (!isTRUE(cross_validation)) {
      gp_requested_models <- unique(stats::na.omit(c(GS_model, GS_model_cv)))
      gp_requested_models <- gp_canonicalize_supported_model_names(gp_requested_models)
      if (length(intersect(gp_requested_models, gp_valid_models)) > 0L &&
          identical(gp_output_level, "predict_only")) {
        gp_output_level <- "full_vc"
        gp_return_se <- TRUE
        gp_full_vc <- TRUE
      }
      if (length(intersect(gp_requested_models, gp_valid_models)) > 0L &&
          identical(tolower(as.character(gp_varcomp_mode)[1L]), "reml") &&
          identical(gp_output_level, "predict_with_se")) {
        gp_output_level <- "full_vc"
        gp_return_se <- TRUE
        gp_full_vc <- TRUE
      }
    }

    asreml_model <- "GBLUP"
    kernel_relationship_models <- gp_multi_environment_kernel_models(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    )
    kernel_relationship_display <- gp_multi_environment_kernel_display_names(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    )
    multi_environment_supported_display <- gp_multi_environment_supported_display_names(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    )


    if (any(c(GS_model, GS_model_cv) %in% c(bayes_valid_models,
                                            bayes_gblup_valid_models,
                                            asreml_model)) &&
        !isTRUE(multi_trait_asreml) &&
        !isTRUE(multi_trait_gp) &&
        !isTRUE(multi_trait_bayes) &&
        !isTRUE(hybrid_asreml) &&
        !isTRUE(hybrid_bayes) &&
        !isTRUE(hybrid_gp) &&
        !isTRUE(hybrid_ml) &&
        !isTRUE(hybrid_dl)){
      if(is.null(random)){
        stop(paste(msg, "Selected model(s) required random term."), call. = FALSE)
      }
    }


    requested_dl_models <- intersect(
      c(GS_model_cv, GS_model),
      gp_deep_learning_supported_models()
    )

    if (!is.null(requested_dl_models) && length(requested_dl_models) > 0){

        if ((!is.null(dl_n_seeds) || !is.null(dl_seeds)) &&
            (isTRUE(cross_validation) || isTRUE(multi_trait_dl) ||
             isTRUE(hybrid_dl) || isTRUE(met_ml_dl))) {
          stop(
            "`dl_n_seeds`/`dl_seeds` currently apply to single-trait DL true prediction (`cross_validation = FALSE`). Use `random_seed` for other DL routes.",
            call. = FALSE
          )
        }
        gp_dl_seed_manifest(
          dl_n_seeds = dl_n_seeds,
          dl_seeds = dl_seeds,
          random_seed = random_seed
        )
        gp_dl_seed_aggregation(dl_seed_aggregation)

        runtime_status <- tryCatch(gp_python_runtime_status(initialize = FALSE, purpose = "dl"), error = function(e) NULL)
        active_python <- if (is.null(runtime_status)) NULL else (runtime_status$active_python %||% NULL)
        if (is.null(runtime_status) || !isTRUE(runtime_status$python_configured)) {
          setup_instructions <- paste(
            "No configured Python runtime was found for deep-learning models.\n",
            "Set PREDICTPRO_DL_PYTHON to a Python interpreter with the required DL packages.\n",
            "You can also run PredictProR::setup_predictdl_env(prefer_gpu = TRUE, cuda = 'auto') as a one-time provisioning helper."
          )
          stop(paste(msg, setup_instructions), call. = FALSE)
        }

    }


    all_models_avail <- c(AI_valid_models, bayes_valid_models,
                          bayes_gblup_valid_models, gp_valid_models, asreml_model)
    all_models_display <- gp_all_model_display_names(
      AI_valid_models = AI_valid_models,
      bayes_valid_models = bayes_valid_models,
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      gp_valid_models = gp_valid_models,
      asreml_model = asreml_model
    )

#####
    holds_out_methods_avail <- c("Hold_Out",
                                 "Stratified_Hold_Out",
                                 "Repeated_Hold_Out",
                                 "Repeated_Stratified_Hold_Out",
                                 "Leave_one_Out")

    Kfolds_methods_avail <- c("K-Folds",
                              "Stratified_K-Folds",
                              "Repeated_K-Folds",
                              "Repeated_Stratified_K-Folds")

    CVs_multi_envs_methods_avail <- c("CV0",
                                      "CV1",
                                      "CV2",
                                      "Repeated_CV0",
                                      "Repeated_CV1",
                                      "Repeated_CV2")

    hybrid_cv_methods_avail <- gp_hybrid_cv_supported_methods()

    all_cv_methods_avail <- c(holds_out_methods_avail,
                              Kfolds_methods_avail,
                              CVs_multi_envs_methods_avail,
                              hybrid_cv_methods_avail)

    cross_validation_meth <- gp_normalize_cv_method_name(
      cross_validation_meth,
      all_cv_methods_avail
    )

    if(is.null(heter_resid)) heter_resid <- FALSE
    ### Check for executing cross_validation
    if(isTRUE(cross_validation)) {
      if(is.null(GS_model_cv) || is.null(cross_validation_meth)) {
        stop(paste(msg, "GS_model_cv and cross_validation_meth cannot be null when cross_validation is TRUE"), call. = FALSE)
      }

      if (!is.null(GS_model) & !is.null(GS_model_cv)) {
        if (!all(GS_model %in% GS_model_cv)) {
          GS_model <- NULL
        }
      }

        if(!all(GS_model_cv%in%all_models_avail)){
          stop(paste(msg,"Invalid model. Choose from: ",
               paste(all_models_display, collapse = ", ")), call. = FALSE)
        }

      if(length(cross_validation_meth)>1){
        stop(paste(msg,'use only one cross_validation method at a time'), call. = FALSE)
      }

      hybrid_cv_allowed <- c(gp_hybrid_cv_supported_methods(), tolower(gp_hybrid_cv_supported_methods()))
      if(!is.null(cross_validation_meth) &&
         !cross_validation_meth%in%all_cv_methods_avail &&
         !(isTRUE(hybrid_asreml) || isTRUE(hybrid_bayes) || isTRUE(hybrid_gp) || isTRUE(hybrid_ml) || isTRUE(hybrid_dl)) &&
         !tolower(cross_validation_meth) %in% tolower(hybrid_cv_allowed)){
        stop(paste(msg,"Invalid cross validation method. Choose from: ",
             paste(all_cv_methods_avail, collapse = ", ")), call. = FALSE)
      }

      if(is.null(eval_metrics)){
        eval_metric_family <- tryCatch(
          gp_resolve_response_family_for_responses(
            pheno_data = pheno_data %||% pheno_data_train %||% pheno_data_test,
            response = response,
            response_family = response_family
          ),
          error = function(e) {
            gp_normalize_response_family(response_family)
          }
        )
        if (identical(eval_metric_family, "auto")) {
          eval_metric_family <- "gaussian"
        }
        early_metric_selection <- gp_prepare_cv_metric_selection(
          eval_metrics = NULL,
          metric_for_ranking = metric_for_ranking,
          ranking_tie_breakers = ranking_tie_breakers,
          response_family = eval_metric_family
        )
        eval_metrics <- early_metric_selection[["eval_metrics"]]
        metric_for_ranking <- early_metric_selection[["metric_for_ranking"]]
        ranking_tie_breakers <- early_metric_selection[["ranking_tie_breakers"]]
      }

      if (!is.null(cross_validation_meth) && any(cross_validation_meth %in% all_cv_methods_avail)) {
        if (is.null(nfolds)) {
          stop(paste(
            msg, "Provide number of nfolds for:",
            paste(cross_validation_meth, collapse = ", "),
            "required for multi-environment genomic prediction."
          ), call. = FALSE)
        }
      }
      ###
      sampling_method <- gp_resolve_cv_sampling_method(cross_validation_meth, sampling_method)

    }

#######################################

    kernel_method_avaliable <- gp_kernel_supported_methods()
    kernel_method <- gp_normalize_kernel_methods(kernel_method)

    if(!is.null(kernel_method)){
      if (!all(kernel_method %in% kernel_method_avaliable)) {
        stop(paste(msg,"Invalid kernel method. Choose from: ",
             paste(kernel_method_avaliable, collapse = ", ")), call. = FALSE)
      }
    }

    gmatrix_method_available <- c("VanRaden",
                                  "Weighted_VanRaden",
                                  "Yang",
                                  "Epistasis",
                                  "Dominance",
                                  "Dominance_Vitezica",
                                  "Dominance_Su",
                                  "Dominance_Heterozygosity")

    if(!is.null(gmatrix_method)){
      if (!all(gmatrix_method %in% gmatrix_method_available)) {
        stop(paste(msg,"Invalid genomic relationship method. Choose from: ",
             paste(gmatrix_method_available, collapse = ", ")), call. = FALSE)
      }
    }
    #############################
    if(is.null(GS_model) & isFALSE(cross_validation)){
      stop(paste(msg, "Genomic prediction model is missing."), call. = FALSE)

    }else{
      if(!isFALSE(cross_validation)){
      if (!any(GS_model_cv %in% c(bayes_valid_models,
                            bayes_gblup_valid_models,
                            gp_valid_models,
                            asreml_model,
                            AI_valid_models))) {
        stop(paste(msg,"Invalid genomic prediction model. Choose from:\n", paste(c(bayes_valid_models,
                                                                      bayes_gblup_valid_models,
                                                                      gp_display_supported_model_names(gp_valid_models),
                                                                      asreml_model,
                                                                      gp_display_supported_model_names(AI_valid_models)), collapse = ", ")),
           call. = FALSE)
      }
      } else {
      if(!is.null(GS_model) & isFALSE(cross_validation)){
        if (!any(GS_model %in% c(bayes_valid_models,
                                 bayes_gblup_valid_models,
                                 gp_valid_models,
                                 asreml_model,
                                 AI_valid_models))) {
          stop(paste(msg,"Invalid genomic prediction model. Choose from:\n", paste(c(bayes_valid_models,
                                                                           bayes_gblup_valid_models,
                                                                           gp_display_supported_model_names(gp_valid_models),
                                                                           asreml_model,
                                                                           gp_display_supported_model_names(AI_valid_models)), collapse = ", ")),
               call. = FALSE)
        }
      }
    }

    }


    # Check for mandatory phenotypic data. The chain of loop is important to ease of checking
    ## and the pheno_data_train and pheno_data_test are converted to pheno_data to make life eaier for checking
    if (is.null(pheno_data)) {
      if (is.null(pheno_data_train) && is.null(pheno_data_test)) {
        stop(paste(msg, "Phenotypic data is missing."), call. = FALSE)
      } else if (!is.null(pheno_data_train) && is.null(pheno_data_test)) {
        if ((isTRUE(cross_validation) && isFALSE(cv_evaluation_only)) || (isFALSE(cross_validation) && isTRUE(cv_evaluation_only))) {
          stop(paste(msg, "Provide a dataframe of phenotypic data for the testing set with NA in the response variable(s) if the interest is to do actual prediction after cross-validation. Otherwise set cross_validation = TRUE and cv_evaluation_only = TRUE.\n"), call. = FALSE)
        } else {
          if (isTRUE(cross_validation) && isTRUE(cv_evaluation_only)) {
            if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
              pheno_data_train <- as.data.frame(pheno_data_train)
              #message(sprintf("%s The pheno data is grouped or a tibble. We fix convert to data.frame.\n", msg))
              stop(paste(msg, "Your pheno data is grouped or a tibble. Convert data.frame.\n"), call. = FALSE)
            }
            pheno_data <- pheno_data_train
          }
        }
      } else if (is.null(pheno_data_train) && !is.null(pheno_data_test)) {
        stop(paste(msg, "Provide a dataframe of phenotypic data for the training set."), call. = FALSE)
      } else {
        if (!is.null(pheno_data_train) && !is.null(pheno_data_test)) {
          if (!identical(colnames(pheno_data_train), colnames(pheno_data_test))) {
            stop(paste(msg, "Column names do not match in the pheno_data_train and pheno_data_test.\n"), call. = FALSE)
          } else {
            if (length(as.character(pheno_data_train[[gen_name]])) > length(as.character(unique(pheno_data_train[[gen_name]]))) &&
                length(as.character(pheno_data_test[[gen_name]])) > length(as.character(unique(pheno_data_test[[gen_name]])))) {
              if (is.null(heter_groups)) {
                stop(paste(msg, paste(
                  "Your phenotypic data has a multi-environment structure,",
                  "but the column containing the environment/location is missing.",
                  "Please provide heter_groups parameter.",
                  "For example: heter_groups = 'locations'.",
                  "If you have 'location' as a column name in your phenotypic data."
                )), call. = FALSE)
              }

              # Step 1: Check if the length of unique values in heter_groups matches
              train_unique <- unique(as.character(pheno_data_train[[heter_groups]]))
              test_unique <- unique(as.character(pheno_data_test[[heter_groups]]))

              if (length(train_unique) != length(test_unique)) {
                stop(paste(msg, "Number of unique values in heter_groups is not the same."), call. = FALSE)
              }

              # Step 2: Check if the names of unique values in heter_groups are identical
              if (!identical(train_unique, test_unique)) {
                stop(paste(msg, "Names of unique values in heter_groups do not match."), call. = FALSE)
              }

              if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
                pheno_data_train <- as.data.frame(pheno_data_train)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ###
              if (inherits(pheno_data_test, "tbl_df") || inherits(pheno_data_test, "grouped_df") || inherits(pheno_data_test, "tbl")) {
                pheno_data_test <- as.data.frame(pheno_data_test)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ##

              pheno_data <- rbind(pheno_data_train, pheno_data_test)

            } else if (length(as.character(pheno_data_train[[gen_name]])) == length(as.character(unique(pheno_data_train[[gen_name]]))) &&
                       length(as.character(pheno_data_test[[gen_name]])) == length(as.character(unique(pheno_data_test[[gen_name]])))) {

              if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
                pheno_data_train <- as.data.frame(pheno_data_train)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ###
              if (inherits(pheno_data_test, "tbl_df") || inherits(pheno_data_test, "grouped_df") || inherits(pheno_data_test, "tbl")) {
                pheno_data_test <- as.data.frame(pheno_data_test)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }

              pheno_data <- rbind(pheno_data_train, pheno_data_test)

            } else {
              stop(paste(msg, "The pheno_train_data and pheno_test_data do not match.\n"), call. = FALSE)
            }
          }
        }
      }
    }

    ######
    if(is.null(gen_name)) stop(paste(msg, "Provide gen_name."), call. = FALSE)
    # Check for gen_name presence
    if(!is.null(pheno_data)){
    if (!gen_name %in% colnames(pheno_data)) {
      stop(paste(msg, sprintf("The specified column '%s' in the pheno_data did not match with your data. Please check and use appropriately.", gen_name)), call. = FALSE)
    }

    }

    if (inherits(pheno_data, "tbl_df") || inherits(pheno_data, "grouped_df") || inherits(pheno_data, "tbl")) {
      pheno_data <- as.data.frame(pheno_data)
      message(paste(msg,"The pheno data is grouped or a tibble. We fix convert to data.frame"))
      #stop("Your pheno data is grouped or a tibble. Convert data.frame")
    }

    if (!is.null(weights)) {
      stage2_weight_input <- gp_attach_stage2_precision_weights(
        pheno_data = pheno_data,
        weights = weights,
        response = response,
        gen_name = gen_name,
        test_set = test_set,
        context = "model_execute Stage 2 observation weights"
      )
      pheno_data <- stage2_weight_input$pheno_data
      # Keep weights attached to phenotype rows so later filtering/reordering
      # cannot detach them from the corresponding Stage 1 estimate.
      weights <- stage2_weight_input$weights
    }

    if(!is.null(geno_data_process)){
      geno_data <- geno_data_process[["snps_matrix"]]

    }

    raw_kernel_guardrail_ctx <- gp_apply_model_input_guardrails(as.list(environment()))
    gmatrix <- raw_kernel_guardrail_ctx[["gmatrix"]]
    gkernel <- raw_kernel_guardrail_ctx[["gkernel"]]
    female_gmatrix <- raw_kernel_guardrail_ctx[["female_gmatrix"]]
    male_gmatrix <- raw_kernel_guardrail_ctx[["male_gmatrix"]]
    pedigree_matrix <- raw_kernel_guardrail_ctx[["pedigree_matrix"]]
    omic1_kernel <- raw_kernel_guardrail_ctx[["omic1_kernel"]]
    omic2_kernel <- raw_kernel_guardrail_ctx[["omic2_kernel"]]
    omic3_kernel <- raw_kernel_guardrail_ctx[["omic3_kernel"]]
    kernel_list <- raw_kernel_guardrail_ctx[["kernel_list"]]
    env_similarity <- raw_kernel_guardrail_ctx[["env_similarity"]]
    # Phase 3.16: write back auto-promoted ctx fields so downstream gates
    # (e.g. line 1547 below: !isTRUE(met_ml_dl)) see the flipped value
    # produced by gp_auto_promote_met_ml_dl_flag().
    met_ml_dl <- isTRUE(raw_kernel_guardrail_ctx[["met_ml_dl"]])
    gmatrix_method <- raw_kernel_guardrail_ctx[["gmatrix_method"]] %||% gmatrix_method

    if (!isTRUE(hybrid_asreml) &&
        !isTRUE(hybrid_bayes) &&
        !isTRUE(hybrid_gp) &&
        !isTRUE(hybrid_ml) &&
        !isTRUE(hybrid_dl) &&
        !isTRUE(multi_trait_asreml) &&
        !isTRUE(multi_trait_gp) &&
        !isTRUE(multi_trait_bayes) &&
        !isTRUE(multi_trait_ml) &&
        !isTRUE(multi_trait_dl) &&
        !isTRUE(met_ml_dl)) {
      validate_general_input_standard(
        pheno_data = pheno_data,
        pheno_data_train = pheno_data_train,
        pheno_data_test = pheno_data_test,
        geno_data = geno_data,
        omic1_data = omic1_data,
        omic2_data = omic2_data,
        omic3_data = omic3_data,
        gmatrix = gmatrix,
        gkernel = gkernel,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups,
        GS_model = GS_model,
        GS_model_cv = GS_model_cv,
        cross_validation = cross_validation,
        cross_validation_meth = cross_validation_meth,
        AI_valid_models = AI_valid_models,
        bayes_valid_models = bayes_valid_models,
        bayes_gblup_valid_models = bayes_gblup_valid_models,
        gp_valid_models = gp_valid_models,
        asreml_model = asreml_model
      )
    }

    gp_fast_models <- gp_lowrank_supported_models()
    gp_fast_kernels_available <- !is.null(gmatrix) ||
      !is.null(omic1_kernel) ||
      !is.null(omic2_kernel) ||
      !is.null(omic3_kernel) ||
      !is.null(kernel_list)
    gp_fast_single_response <- length(response) == 1L && !is.null(response)
    gp_fast_no_duplicates <- !is.null(pheno_data) &&
      !is.null(gen_name) &&
      gen_name %in% names(pheno_data) &&
      !gp_is_multi_environment_trait_panel(
        pheno_data = pheno_data,
        gen_name = gen_name,
        heter_groups = heter_groups,
        response = response,
        response_family = response_family
      )
    gp_fast_plain_run <- !isTRUE(cross_validation) &&
      !isTRUE(multi_trait_asreml) &&
      !isTRUE(multi_trait_gp) &&
      !isTRUE(multi_trait_bayes) &&
      !isTRUE(multi_trait_ml) &&
      !isTRUE(multi_trait_dl) &&
      !isTRUE(hybrid_asreml) &&
      !isTRUE(hybrid_bayes) &&
      !isTRUE(hybrid_gp) &&
      !isTRUE(hybrid_ml) &&
      !isTRUE(hybrid_dl) &&
      !isTRUE(met_ml_dl)

    if (gp_fast_plain_run &&
        gp_fast_single_response &&
        gp_fast_no_duplicates &&
        gp_fast_kernels_available &&
        !is.null(GS_model) &&
        length(GS_model) == 1L &&
        GS_model %in% gp_fast_models &&
        identical(gp_resolve_response_family_for_responses(pheno_data, response, response_family), "gaussian")) {
      run_profile <- gp_runtime_profile_new("model_execute_gp_fast")
      run_profile <- gp_runtime_profile_mark(run_profile, "prepare_gp_fast")

      gp_cache <- gp_backend_prepare_kernel_cache(
        ids = as.character(pheno_data[[gen_name]]),
        gmatrix = gmatrix,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list,
        kernel_weights = lowrank_kernel_weights
      )

      fit <- gp_backend_gaussian_model(
        model_name = GS_model,
        pheno_data = pheno_data,
        response = response,
        gen_name = gen_name,
        gmatrix = gmatrix,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list,
        heter_groups = heter_groups,
        test_set = test_set,
        lowrank_eps_trace = lowrank_eps_trace,
        lowrank_max_rank = lowrank_max_rank,
        lowrank_jitter = lowrank_jitter,
        lowrank_noise_grid = lowrank_noise_grid,
        lowrank_kernel_weights = lowrank_kernel_weights,
        system_database = system_database,
        gp_backend_cache = gp_cache,
        gp_backend = gp_backend,
        gp_output_level = gp_output_level,
        gp_varcomp_mode = gp_varcomp_mode,
        gp_fa_rank = gp_fa_rank,
        gp_prediction_output = gp_prediction_output,
        gp_factor_cache = gp_factor_cache,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        gp_iters = gp_iters,
         gp_lr = gp_lr,
         gp_engine = gp_engine,
         gp_learn_scales = gp_learn_scales,
         observation_weights = weights
       )

      run_profile <- gp_runtime_profile_mark(run_profile, "fit_gp_fast")

      GS_model_public <- gp_public_model_label(GS_model)[1]
      fit[["model_parameters"]] <- gp_relabel_public_model_parameters(
        fit[["model_parameters"]],
        GS_model
      )
      fit <- gp_standardize_model_prediction_outputs(
        res_model_output = fit,
        gen_name = gen_name
      )
      fit_public <- gp_standardize_public_model_result(
        res_model_output = fit,
        gen_name = gen_name,
        heter_groups = heter_groups,
        model = GS_model,
        response_family = "gaussian"
      )
      summary_stat <- summary_statistics_AI(
        predicted_object = fit_public[["predicted_values"]],
        pheno_object = pheno_data[!is.na(pheno_data[[response]]), , drop = FALSE],
        response = response,
        test_set = as.character(test_set %||% pheno_data[[gen_name]][is.na(pheno_data[[response]])]),
        geno_omic_object = fit[["lowrank_feature_matrix"]],
        model_parameters = fit_public[["model_parameters"]],
        eval_metrics = eval_metrics,
        response_family = "gaussian",
        GS_model = GS_model_public
      )

      run_profile <- gp_runtime_profile_mark(run_profile, "summarize_gp_fast")

      run_metadata <- gp_runtime_metadata(
        python_path = gp_detect_gp_python(),
        preferred_python = gp_detect_gp_python(),
        purpose = "gp",
        include_accelerator = FALSE,
        context = "model_execute_gp_fast",
        extra_fields = c(
          list(
            mode = "gp_fast",
            task_unit = "genotype_level",
            task_models = GS_model_public,
            task_models_canonical = gp_canonicalize_supported_model_names(GS_model)[1],
            task_traits = paste(response, collapse = ","),
            task_replications = 1L
          ),
          gp_runtime_profile_metadata_fields(run_profile)
        )
      )

      fast_output <- results_handling(
        GS_model = GS_model_public,
        res_model_output = fit_public,
        res_summary_stat = summary_stat,
        res_plot = NULL,
        res_plot_mean = NULL,
        res_plot_result_diagnostic = NULL,
        test_diagonistic_plots = fit_public[["diagnostic_plots"]],
        res_mod_results_cv_per_trait_model = NULL,
        res_plot_result_diagnostic_cv_only = NULL,
        cv_results_processed = NULL,
        cv_results_raw = NULL,
        geno_qc_stat = geno_qc_stat,
        system_database = system_database,
        plot_filename = as.character(response[[1L]]),
        plot_extension = plot_extension,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_units = plot_units,
        plot_dpi = plot_dpi,
        feature_selected = NULL,
        feature_score_metadata = NULL,
        run_metadata = run_metadata
      )

      fast_return <- stats::setNames(list(fast_output), response)
      fast_return[["model_results_by_trait"]] <- stats::setNames(list(fit_public), response)
      fast_return[["model_results_by_model"]] <- stats::setNames(
        list(stats::setNames(list(fit_public), response)),
        make.names(GS_model_public)
      )
      fast_return[["run_metadata"]] <- run_metadata
      return(fast_return)
    }

    ## Check for scenrio where user provide pheno_train and pheno_test.
    ## Corresponding train and test geno or omic data must be provided

    # Define the error message
    error_message <- "When providing both pheno_data_train and pheno_data_test, at least one of the following pairs must be provided:
                  (train_geno_data and test_geno_data),
                  (train_omic1_data and test_omic1_data),
                  (train_omic2_data and test_omic2_data),
                  (train_omic3_data and test_omic3_data)."

    # Check if GS_model or GS_model_cv contains valid models
    valid_models_using_X_variables <- c(AI_valid_models, bayes_valid_models)
    if ((!is.null(GS_model) && any(GS_model %in% valid_models_using_X_variables)) ||
        (!is.null(GS_model_cv) && any(GS_model_cv %in% valid_models_using_X_variables))) {

      # Check if pheno_data_train and pheno_data_test are provided
      if (!is.null(pheno_data_train) && !is.null(pheno_data_test)) {

        # Check if at least one of the required pairs is provided
        conditions_met <- sum(
          !is.null(train_geno_data) && !is.null(test_geno_data),
          !is.null(train_omic1_data) && !is.null(test_omic1_data),
          !is.null(train_omic2_data) && !is.null(test_omic2_data),
          !is.null(train_omic3_data) && !is.null(test_omic3_data)
        )

        # If none of the pairs are provided, stop with an error message
        if (conditions_met < 1) {
          stop(paste(msg,error_message), call. = FALSE)
        }
      }
    } ### End

    ## Check for when both train and test set data are present in single file
    ##############
    # Define conditions
    condition1 <- is.null(geno_data) && is.null(omic1_data) && is.null(omic2_data) && is.null(omic3_data)
    condition1_1 <- is.null(gmatrix_method) && is.null(kernel_method)

    condition2 <- is.null(gmatrix) && is.null(gkernel) && is.null(omic1_kernel) && is.null(omic2_kernel) && is.null(omic3_kernel) && is.null(kernel_list)
    hybrid_relationship_inputs_present <- (isTRUE(hybrid_asreml) || isTRUE(hybrid_bayes) || isTRUE(hybrid_gp)) &&
      (
        !is.null(geno_data) ||
          !is.null(female_geno_data) || !is.null(male_geno_data) ||
          !is.null(female_gmatrix) || !is.null(male_gmatrix)
      )

    ####
    # Check if all required data sets are null
    condition11 <- is.null(geno_data) && is.null(train_geno_data) && is.null(test_geno_data)
    condition12 <- is.null(omic1_data) && is.null(train_omic1_data) && is.null(test_omic1_data)
    condition13 <- is.null(omic2_data) && is.null(train_omic2_data) && is.null(test_omic2_data)
    condition14 <- is.null(omic3_data) && is.null(train_omic3_data) && is.null(test_omic3_data)

    error_message <- paste(
      msg,
      "Machine-learning and deep-learning models require marker/omics features or one or more PSD kernels; Bayesian marker-effect models require marker/omics features. Kernel Bayesian models (GBLUP_BRR and RKHS) accept kernels."
    )
    # Kernel-aware model families consume every supplied kernel. ML/DL convert
    # each kernel to a named eigenfeature block; GBLUP_BRR/RKHS retain the
    # kernels as covariance terms. Bayesian alphabet models still require raw
    # marker/omics features because their priors are defined on marker effects.
    kernel_feature_inputs_present <- !is.null(gmatrix) || !is.null(gkernel) ||
      !is.null(omic1_kernel) || !is.null(omic2_kernel) ||
      !is.null(omic3_kernel) || !is.null(kernel_list)
    kernel_only_valid_models <- unique(c(AI_valid_models, bayes_gblup_valid_models))
    kernel_only_request_is_valid <- function(models) {
      requested <- intersect(
        as.character(models),
        unique(c(AI_valid_models, bayes_valid_models))
      )
      length(requested) > 0L &&
        all(requested %in% kernel_only_valid_models) &&
        isTRUE(kernel_feature_inputs_present)
    }
    # Check if GS_model is not null and belongs to valid models
    if (!is.null(GS_model) && any(GS_model %in% c(AI_valid_models, bayes_valid_models)) && isFALSE(cross_validation)) {
      # Check if all conditions are true
      if (condition11 && condition12 && condition13 && condition14 &&
          !kernel_only_request_is_valid(GS_model)) {
        stop(paste(msg,error_message), call. = FALSE)
      }
    }

    if (!is.null(GS_model_cv) && any(GS_model_cv %in% c(AI_valid_models, bayes_valid_models))) {
      # Check if all conditions are true
      if (condition11 && condition12 && condition13 && condition14 &&
          !kernel_only_request_is_valid(GS_model_cv)) {
        stop(error_message, call. = FALSE)
      }
    }
    ########################
    ### Bayesian GBLUP further check
    if (any(c(bayes_gblup_valid_models, gp_valid_models) %in% na.omit(c(GS_model, GS_model_cv)))) {
      omics <- list(omic1_data, omic2_data, omic3_data, geno_data,
                    train_geno_data, train_omic1_data,
                    train_omic2_data, train_omic3_data,
                    test_geno_data,
                    test_omic1_data, test_omic2_data,
                    test_omic3_data)
      omics_names <- c("omic1_data", "omic2_data", "omic3_data", "geno_data",
                       "train_geno_data", "train_omic1_data",
                       "train_omic2_data", "train_omic3_data",
                       "test_geno_data",
                       "test_omic1_data", "test_omic2_data",
                       "test_omic3_data")
      names(omics) <- omics_names

      #omics_kernel <- list(omic1_kernel, omic2_kernel, omic3_kernel, gmatrix, gkernel)

      # Remove NULL elements the list
      omics <- omics[!sapply(omics, is.null)]
      # Remove NULL elements from the list
      #omics <- Filter(Negate(is.null), omics)
      #omics_kernel <-  Filter(Negate(is.null), omics_kernel)

      ######
      # Check if kernel_method is NULL
      if (any(grepl("geno_data", names(omics)))){
        if (is.null(gmatrix_method)) {
          stop(paste(msg,"Provide gmatrix method if you use genomic data to calculate relationship matrix to fit Bayesian GBLUP model.\n"), call. = FALSE)
          #stop("Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n", call. = FALSE)
        }

      }

      if (length(grep("omic", names(omics))) > 0) {
        if (is.null(kernel_method)) {
          stop(paste(msg,"Provide Kernel method if you use omics data to calculate relationship matrix to fit Bayesian GBLUP model.\n"), call. = FALSE)
          #stop("Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n", call. = FALSE)
        }

      }
    }

    ###################

    # Check for ASReml requirement for GBLUP
    ## Define error message
    error_message <- "ASReml software is required to fit GBLUP for single or multi-environment.\n"

    if (!is.null(engine)) {
      if (!"GBLUP" %in% na.omit(c(GS_model, GS_model_cv))) {
        engine <- NULL
      }
    }

    ### Check that when asreml is used required format for omic data is provided
    if(!is.null(engine)){
      # Create a list of omics data and kernel
      omics <- list(omic1_data, omic2_data, omic3_data, geno_data,
                    train_geno_data, train_omic1_data,
                    train_omic2_data, train_omic3_data,
                    test_geno_data,
                    test_omic1_data, test_omic2_data,
                    test_omic3_data)
      omics_names <- c("omic1_data", "omic2_data", "omic3_data", "geno_data",
                       "train_geno_data", "train_omic1_data",
                       "train_omic2_data", "train_omic3_data",
                       "test_geno_data",
                       "test_omic1_data", "test_omic2_data",
                       "test_omic3_data")
      names(omics) <- omics_names

      omics_kernel <- list(omic1_kernel, omic2_kernel, omic3_kernel, gmatrix, gkernel)

      # Remove NULL elements the list
      omics <- omics[!sapply(omics, is.null)]
      # Remove NULL elements from the list
      #omics <- Filter(Negate(is.null), omics)
      omics_kernel <-  Filter(Negate(is.null), omics_kernel)
      # Define a function to check if a matrix is square
      is_square_matrix <- function(mat) {
        if (!is.matrix(mat)) {
          mat <- as.matrix(mat)
        }
        return(nrow(mat) == ncol(mat))
      }

      # Apply the function to each element of the list to check if it's a square matrix
      if (length(omics) != 0) {
        square_matrices <- sapply(omics, is_square_matrix)
        if (any(square_matrices)) {
          stop(paste(msg, "Omic data or geno data is a square matrix. If this is a relationship matrix, provide it as: omic_kernel or gmatrix.\n"), call. = FALSE)
          #stop("Omic data or geno data is a square matrix. If this is a relationship matrix, provide it as: omic_kernel or gmatrix.\n", call. = FALSE)
        }
          # Check if kernel_method is NULL
          #if("geno_data"%in%names(omics)){
        if (any(grepl("geno_data", names(omics)))){
            if (is.null(gmatrix_method)) {
              stop(paste(msg,"Provide gmatrix method if you use genomic data to calculate relationship matrix to fit GBLUP model.\n"), call. = FALSE)
              #stop("Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n", call. = FALSE)
            }

          }

          if (length(grep("omic", names(omics))) > 0) {
            if (is.null(kernel_method)) {
              stop(paste(msg,"Provide Kernel method if you use omics data to calculate relationship matrix to fit GBLUP model.\n"), call. = FALSE)
              #stop("Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n", call. = FALSE)
            }

          }

        #}

      }
      ###
      if (length(omics_kernel) != 0) {
        square_matrices <- sapply(omics_kernel, is_square_matrix)
        if (!all(square_matrices)) {
          stop(paste(msg, "Omic kernel or gmatrix or gkernel should be a relationship matrix.\n"), call. = FALSE)
          #stop("Omic kernel or gmatrix or gkernel should be a relationship matrix.\n", call. = FALSE)
        }

      }
    if(isFALSE(cross_validation) & !is.null(GS_model)){
    if (any(GS_model == "GBLUP") && engine != "asreml") {
      stop(paste(msg,error_message), call. = FALSE)
     }
    } else {
      if(!is.null(GS_model_cv)){
      if ("GBLUP" %in%GS_model_cv && engine != "asreml") {
        stop(paste(msg,error_message), call. = FALSE)
      }
      }
    }

    } else {
     if(isFALSE(cross_validation) & !is.null(GS_model)){
      if(is.null(engine) & any(GS_model == "GBLUP")){
        stop(paste(msg,error_message), call. = FALSE)
      }

     } else {
       if(!is.null(GS_model_cv)){
       if(is.null(engine) & "GBLUP" %in%GS_model_cv){
         stop(paste(msg,error_message), call. = FALSE)
       }

       }
    }

    }


    # Check for multi-environment structure and required model inputs
    if (isTRUE(hybrid_asreml) || isTRUE(hybrid_bayes) || isTRUE(hybrid_gp) || isTRUE(hybrid_ml) || isTRUE(hybrid_dl)) {
      validate_hybrid_input_standard(
        pheno_data = pheno_data,
        pheno_data_train = pheno_data_train,
        pheno_data_test = pheno_data_test,
        geno_data = geno_data,
        gmatrix = gmatrix,
        female_gmatrix = female_gmatrix,
        male_gmatrix = male_gmatrix,
        female_geno_data = female_geno_data,
        male_geno_data = male_geno_data,
        response = response,
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        hybrid_asreml = hybrid_asreml,
        hybrid_bayes = hybrid_bayes,
        hybrid_gp = hybrid_gp,
        hybrid_ml = hybrid_ml,
        hybrid_dl = hybrid_dl,
        GS_model = GS_model,
        GS_model_cv = GS_model_cv,
        cross_validation = cross_validation,
        cross_validation_meth = cross_validation_meth
      )
    }
    if (isTRUE(multi_trait_asreml) || isTRUE(multi_trait_gp) || isTRUE(multi_trait_bayes) || isTRUE(multi_trait_ml) || isTRUE(multi_trait_dl)) {
      validate_multi_trait_input_standard(
        pheno_data = pheno_data,
        pheno_data_train = pheno_data_train,
        pheno_data_test = pheno_data_test,
        geno_data = geno_data,
        omic1_data = omic1_data,
        omic2_data = omic2_data,
        omic3_data = omic3_data,
        gmatrix = gmatrix,
        gkernel = gkernel,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list,
        response = response,
        gen_name = gen_name,
        response_family = response_family,
        multi_trait_asreml = multi_trait_asreml,
        multi_trait_ml = multi_trait_ml,
        multi_trait_dl = multi_trait_dl,
        multi_trait_gp = multi_trait_gp,
        multi_trait_bayes = multi_trait_bayes,
        heter_groups = heter_groups,
        GS_model = GS_model,
        GS_model_cv = GS_model_cv,
        cross_validation = cross_validation,
        cross_validation_meth = cross_validation_meth
      )
    }
    if (isTRUE(met_ml_dl)) {
      validate_met_input_standard(
        pheno_data = pheno_data,
        pheno_data_train = pheno_data_train,
        pheno_data_test = pheno_data_test,
        geno_data = geno_data,
        omic1_data = omic1_data,
        omic2_data = omic2_data,
        omic3_data = omic3_data,
        gmatrix = gmatrix,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list,
        response = response,
        response_family = response_family,
        gen_name = gen_name,
        heter_groups = heter_groups,
        met_ml_dl = met_ml_dl,
        GS_model = GS_model,
        GS_model_cv = GS_model_cv,
        cross_validation = cross_validation,
        cross_validation_meth = cross_validation_meth
      )
    }

    observed_multi_environment <- gp_is_multi_environment_trait_panel(
      pheno_data = pheno_data,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response = response,
      response_family = response_family
    )

    if (isTRUE(observed_multi_environment)){
      if ((any(!is.null(GS_model_cv) & GS_model_cv %in% bayes_valid_models) & isTRUE(cross_validation)) ||
          (any(!is.null(GS_model) & GS_model %in% bayes_valid_models))) {
        stop(paste(
          msg,
          "Your phenotypic data has a multi-environment structure,",
          "Bayesian marker-regression models are currently single-environment only.",
          "For multi-environment genomic prediction, use one of:",
          paste(multi_environment_supported_display, collapse = ", "),
          "\n"
        ), call. = FALSE)
      }

      if (is.null(heter_groups)) {
        stop(paste(msg,paste(
          "Your phenotypic data has a multi-environment structure,",
          "but the column containing the environment/location is missing.",
          "Please provide heter_groups parameter.",
          "For example: heter_groups = 'locations'.",
          "If you have location as column name in your phenotypic data.\n"
        )), call. = FALSE)
      }

        if(isTRUE(cross_validation) &&
           !isTRUE(multi_trait_asreml) &&
           !isTRUE(hybrid_asreml) && !isTRUE(hybrid_bayes) &&
           !isTRUE(hybrid_gp) &&
           !isTRUE(hybrid_ml) && !isTRUE(hybrid_dl)){
          if (!any(cross_validation_meth %in% CVs_multi_envs_methods_avail)) {
            stop(paste(msg, paste(
              "Provide any of the following as cross-validation strategy for multi-environment GS:",
              paste(CVs_multi_envs_methods_avail, collapse = ", ")
            )), call. = FALSE)
          }
      }


      ###
      if(isFALSE(cross_validation)){
        if (any(c(GS_model, GS_model_cv) %in% AI_valid_models) &&
            !isTRUE(hybrid_gp) && !isTRUE(hybrid_ml) && !isTRUE(hybrid_dl)) {
          if (!(isTRUE(met_ml_dl) && any(c(GS_model, GS_model_cv) %in% gp_met_supported_models()))) {
            stop(paste(msg, paste(
              "Your data suggest multi-environment structure. Use a supported MET genomic prediction model:",
              paste(multi_environment_supported_display, collapse = ", ")
            )), call. = FALSE)
          }
        }
      }

      ### This is important for asreml for multi-environment analysis
      # Define the error messages
      missing_var_cov_str_msg <- paste(msg, "Your data suggest multi-environment but variance-covariance structure is missing. Choose from:", paste(var_cov_str_available, collapse = ", "))
      missing_heter_resid_msg <- paste(msg, "Your data suggest multi-environment with a variance-covariance structure (var_cov_str); set heter_resid = TRUE (environment-specific residual variances).")
      invalid_var_cov_str_msg <- paste(msg, "Invalid output variance-covariance structure. Choose from:", paste(var_cov_str_available, collapse = ", "))

      # Check GS_model
      if (!is.null(GS_model) && any(GS_model %in% c("GBLUP"))||
          !is.null(GS_model_cv) && any(GS_model_cv %in% c("GBLUP"))) {
        if (!is.null(heter_groups) && isTRUE(heter_resid) && is.null(var_cov_str)) {
          stop(missing_var_cov_str_msg, call. = FALSE)
        } else if (!is.null(heter_groups) && isFALSE(heter_resid) && !is.null(var_cov_str)) {
          stop(missing_heter_resid_msg, call. = FALSE)
        } else if (!is.null(var_cov_str) && !(var_cov_str %in% var_cov_str_available)) {
          stop(invalid_var_cov_str_msg, call. = FALSE)
        }
      }

      # Check GS_model_cv for cross-validation
      if (isTRUE(cross_validation) && !is.null(GS_model_cv) && any(GS_model_cv %in% c("GBLUP"))) {
        if (!is.null(heter_groups) && isTRUE(heter_resid) && is.null(var_cov_str)) {
          stop(missing_var_cov_str_msg, call. = FALSE)
        } else if (!is.null(heter_groups) && isFALSE(heter_resid) && !is.null(var_cov_str)) {
          stop(missing_heter_resid_msg, call. = FALSE)
        } else if (!is.null(var_cov_str) && !(var_cov_str %in% var_cov_str_available)) {
          stop(invalid_var_cov_str_msg, call. = FALSE)
        }
      }

      # Check for required inputs for multi-environment relationship/kernel models
      # Check if condition1 is true
      # Define the error message
      error_messagee <- paste(
        msg,
        paste("To fit multi-environment relationship/kernel genomic prediction models",
        paste0("(", paste(kernel_relationship_display, collapse = ", "), "),"),
        "you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel),",
        "or an omics kernel. Additionally, you can provide genomic or omics data.",
        "Ensure you provide instructions on the genomic relationship matrix method",
        "or kernel method to calculate the relationship matrix.\n"
      ))

      # Check GS_model
      if (!is.null(GS_model) && any(GS_model %in% kernel_relationship_models)) {
        if (!isTRUE(hybrid_relationship_inputs_present) &&
            condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }

      # Check GS_model_cv for cross-validation
      if (isTRUE(cross_validation) && !is.null(GS_model_cv) && any(GS_model_cv %in% kernel_relationship_models)) {
        if (!isTRUE(hybrid_relationship_inputs_present) &&
            condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }


    } else {

    # Check for single environment GBLUP and Bayesian models
    if (!isTRUE(observed_multi_environment)){
 #### In case user erroneously provide this while it is a single location
      preserve_hybrid_gp_heter_groups <- isTRUE(hybrid_gp) &&
        !is.null(heter_groups) &&
        heter_groups %in% names(pheno_data) &&
        isTRUE(gp_is_multi_environment_trait_panel(
          pheno_data = pheno_data,
          gen_name = gen_name,
          heter_groups = heter_groups,
          response = response,
          response_family = response_family
        ))
      if (!isTRUE(preserve_hybrid_gp_heter_groups)) {
        heter_groups <- NULL
        heter_resid <- FALSE
        if (!isTRUE(multi_trait_asreml)) {
          var_cov_str <- NULL
        }
      }
      if (isTRUE(cross_validation) &&
          !is.null(cross_validation_meth) &&
          !isTRUE(preserve_hybrid_gp_heter_groups) &&
          cross_validation_meth %in% CVs_multi_envs_methods_avail) {
        stop(paste(msg, paste("You selected cross validation method for multi-environment",
                   "but your phenotypic data is a single environment.\n")), call. = FALSE)
      }
      # Check for required inputs for single-environment relationship/kernel models
      # Define the error message
      error_messageee <- paste(
        msg,
        paste("To fit a single-environment relationship/kernel genomic prediction model",
        paste0("(", paste(kernel_relationship_display, collapse = ", "), "),"),
        "you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel),",
        "or an omics kernel. Additionally, you can provide genomic or omics data.",
        "Ensure you provide instructions on the genomic relationship matrix method",
        "or kernel method to calculate the relationship matrix."
      ))

      # Check GS_model
      if (!isTRUE(hybrid_asreml) &&
          !isTRUE(hybrid_bayes) &&
          !isTRUE(hybrid_gp) &&
          !is.null(GS_model) &&
          any(GS_model %in% kernel_relationship_models)) {
        if (condition1 != condition1_1 && isTRUE(condition2) & isTRUE(condition11) & isTRUE(condition12) & isTRUE(condition13) & isTRUE(condition14)) {
          stop(error_messageee, call. = FALSE)
        }
      }

      # Check GS_model_cv for cross-validation
      if (!isTRUE(hybrid_asreml) &&
          !isTRUE(hybrid_bayes) &&
          !isTRUE(hybrid_gp) &&
          isTRUE(cross_validation) &&
          !is.null(GS_model_cv) &&
          any(GS_model_cv %in% kernel_relationship_models)) {
        if (condition1 != condition1_1 && isTRUE(condition2) & isTRUE(condition11) & isTRUE(condition12) & isTRUE(condition13) & isTRUE(condition14)) {
          stop(error_messageee, call. = FALSE)
        }
      }

      # Check for required inputs for Bayesian models
      if(isTRUE(cross_validation)){
        if(any(GS_model_cv%in%c(bayes_valid_models, AI_valid_models))){
          if (isTRUE(condition1) & isTRUE(condition11) & isTRUE(condition12) & isTRUE(condition13) & isTRUE(condition14) &&
              !kernel_only_request_is_valid(GS_model_cv)) {
            #print('ok')
            stop(error_message, call. = FALSE)
          }
        }
      } else{
      if (any(GS_model %in% c(bayes_valid_models, AI_valid_models))){
        if (isTRUE(condition1) & isTRUE(condition11) & isTRUE(condition12) & isTRUE(condition13) & isTRUE(condition14) &&
            !kernel_only_request_is_valid(GS_model)) {
          #print('ok')
          stop(error_message, call. = FALSE)
        }

      }

    }

    }

    }
    ##

    # Main Script Not use again remove main scripts
    # checkForASReml(engine, GS_model, GS_model_cv, cross_validation, msg)
    # validateMultiEnvironment(pheno_data, gen_name, heter_groups, heter_resid, var_cov_str, GS_model, var_cov_str_available,cross_validation, GS_model_cv, msg)
    # validateModelRequirements(GS_model, GS_model_cv, bayes_gblup_valid_models, condition1, condition1_1, condition2, cross_validation, msg)
###############
    gp_guard_asreml_cv0_fixed_environment(
      cross_validation = cross_validation,
      cross_validation_meth = cross_validation_meth,
      GS_model_cv = GS_model_cv,
      engine = engine,
      fixed = fixed,
      heter_groups = heter_groups,
      msg = msg
    )

    # Ensure response_var exist in df
    if (!all(response %in% colnames(pheno_data))) {
      stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
    }

    # Create a logical matrix indicating NA positions for response variables
    na_matrix <- is.na(pheno_data[response])

    # Check if all rows have the same NA pattern
    # Ensure response_var exist in df
    if (!all(response %in% colnames(pheno_data))) {
      stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
    }

    # Create a logical matrix indicating NA positions for response variables
    na_matrix <- is.na(pheno_data[response])

    # Check if all rows have the same NA pattern
    # Ensure response_var exist in df
    if (!all(response %in% colnames(pheno_data))) {
      stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
    }

    # Create a logical matrix indicating NA positions for response variables
    na_matrix <- is.na(pheno_data[response])

    # Missing response values are interpreted trait by trait as prediction
    # targets. Different traits can therefore have different inferred test sets.

##########
    datasets_geno_omic <- list(geno_data, omic1_data, omic2_data, omic3_data)
    dataset_names_geno_omic <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")

    geno_omic_rownames_check <- check_names_consistency(datasets_geno_omic, dataset_names_geno_omic)

    if(isFALSE(geno_omic_rownames_check)){
      stop(paste(msg, paste("The omic/and or genotypic data do not have consistent rownames.",
                            "provide genomic or omics data with consistent rownames.")), call. = FALSE)
    }

    ###
    datasets_gmatrix_omic_kernel <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
    dataset_names_gmatrix_omic_kernel <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")

    gmatrix_omic_kernel_rownames_check <- check_names_consistency(datasets_gmatrix_omic_kernel, dataset_names_gmatrix_omic_kernel)

    if(isFALSE(gmatrix_omic_kernel_rownames_check)){
      stop(paste(msg, paste("The omic_kernel/and or gmatrix data do not have onsistent row and column names.",
                            "provide  omic_kernel/and or gmatrix data with consistent row and column names.")), call. = FALSE)
    }

    datasets_index <- which(!sapply(datasets_geno_omic, is.null))
    datasets_index_kernel <- which(!sapply(datasets_gmatrix_omic_kernel, is.null))

    names(datasets_geno_omic)[datasets_index] <- dataset_names_geno_omic[datasets_index]
    names(datasets_gmatrix_omic_kernel)[datasets_index] <- dataset_names_gmatrix_omic_kernel[datasets_index]

    if(length(datasets_index)!=0 && length(datasets_index_kernel)!=0){
      datasets_index_kernel <- NULL
    }
############
# Every line in the genotype / kernel data is predicted: lines without a
# phenotype row are added as NA rows (in every environment for MET data).
# They used to be dropped whenever pheno_data had any NA, so whether a
# genotyped-only line was predicted depended on how other lines were entered.
if (!isTRUE(hybrid_asreml) && !isTRUE(hybrid_bayes) && !isTRUE(hybrid_gp)){

if (length(datasets_index) != 0) {

  if(length(rownames(datasets_geno_omic[[1]])) > length(unique(pheno_data[[gen_name]]))){
    diff_gid <- setdiff(rownames(datasets_geno_omic[[1]]),
                        unique(pheno_data[[gen_name]]))

  pheno_data <-  add_extra_gid_from_geno_omic_to_pheno(geno_data = datasets_geno_omic[[1]],
                                                       pheno_data = pheno_data,
                                                       gen_name = gen_name,
                                                       response_var = response,
                                                       heter_group = heter_groups)

  }
}

#####

if (length(datasets_index_kernel) != 0) {

  if(length(rownames(datasets_gmatrix_omic_kernel[[1]])) > length(unique(pheno_data[[gen_name]]))){
    diff_gid <- setdiff(rownames(datasets_gmatrix_omic_kernel[[1]]),
                        unique(pheno_data[[gen_name]]))

    pheno_data <-  add_extra_gid_from_geno_omic_to_pheno(geno_data = datasets_gmatrix_omic_kernel[[1]],
                                                         pheno_data = pheno_data,
                                                         gen_name = gen_name,
                                                         response_var = response,
                                                         heter_group = heter_groups)

  }
}

} else if (!isTRUE(hybrid_asreml) && !isTRUE(hybrid_bayes) && !isTRUE(hybrid_gp)){

  if (length(datasets_index) != 0) {

    if(length(rownames(datasets_geno_omic[[1]])) > length(unique(pheno_data[[gen_name]]))){
      # Keep geno/omic-only individuals so pheno_geno_match() can infer them
      # as testing IDs.
      # prediction from seeing genomic/omic-only test candidates.
    }
  }
 ########
  if (length(datasets_index_kernel) != 0) {

    if(length(rownames(datasets_gmatrix_omic_kernel[[1]])) > length(unique(pheno_data[[gen_name]]))){
      if("gmatrix"%in%names(datasets_gmatrix_omic_kernel)){
        gmatrix <- gmatrix[rownames(gmatrix)%in%unique(pheno_data[[gen_name]]), colnames(gmatrix)%in%unique(pheno_data[[gen_name]])]
      }
      ##
      if("omic1_kernel"%in%names(datasets_gmatrix_omic_kernel)){
        omic1_kernel <- omic1_kernel[rownames(omic1_kernel)%in%unique(pheno_data[[gen_name]]), ]
      }
      ##
      if("omic2_kernel"%in%names(datasets_gmatrix_omic_kernel)){
        omic2_kernel <- omic2_kernel[rownames(omic2_kernel)%in%unique(pheno_data[[gen_name]]), ]
      }
      ###
      if("omic3_kernel"%in%names(datasets_gmatrix_omic_kernel)){
        omic3_kernel <- omic3_kernel[rownames(omic3_kernel)%in%unique(pheno_data[[gen_name]]), ]
      }
    }
  }

}

# Multi-environment true prediction: predict every line in every environment
# (the full line x environment grid) with every engine. met_input_keys keeps
# the input line x environment records, to label the added rows "Unobserved".
met_input_keys <- NULL
if (isTRUE(met_predict_all_environments) && !isTRUE(cross_validation) &&
    !is.null(heter_groups) && heter_groups %in% names(pheno_data) &&
    anyDuplicated(as.character(pheno_data[[gen_name]])) > 0L &&
    !isTRUE(hybrid_asreml) && !isTRUE(hybrid_bayes) && !isTRUE(hybrid_gp) &&
    !isTRUE(hybrid_ml) && !isTRUE(hybrid_dl)) {
  met_input_keys <- unique(paste(pheno_data[[gen_name]], pheno_data[[heter_groups]], sep = "\r"))
  pheno_data <- gp_met_complete_grid(
    pheno_data = pheno_data, gen_name = gen_name, heter_groups = heter_groups,
    response = response, fixed = fixed,
    weights = if (is.character(weights) && length(weights) == 1L) weights else NULL
  )
}

#########################
### Check phenotype_to_model for details
 #    This serve as gateway between phenotype-precheck function and readiness of
 #    the phenotypic data for model fitting.

 phenotype_preprocess_stage <- "phenotype_to_model"
 if (gp_shared_preprocess_cache_has(
   predictpror_preprocess_cache,
   phenotype_preprocess_stage
 )) {
   pheno_clean <- gp_shared_preprocess_cache_get(
     predictpror_preprocess_cache,
     phenotype_preprocess_stage
   )
 } else {
   pheno_clean <- phenotype_to_model(pheno_data = pheno_data,
                                     pheno_data_train = pheno_data_train,
                                     pheno_data_test = pheno_data_test,
                                     train_set = train_set,
                                     test_set = test_set,
                                     response = response,
                                     response_family = response_family,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     random = random,
                                     fixed = fixed,
                                     type_pheno = if(!is.null(pheno_data_test)) "test_set" else NULL)
   gp_shared_preprocess_cache_set(
     predictpror_preprocess_cache,
     phenotype_preprocess_stage,
     pheno_clean
   )
 }

 if("test_set"%in%names(pheno_clean)){
   test_set <- pheno_clean[["test_set"]]
 }
 #if(length(pheno_clean)==0) stop("pheno is null")
## pheno_clean is a list that can have one or two elements
 ## One element if only pheno_data is provided
 ## Two elements if pheno_data/pheno_training and pheno_data_testing was provided as input.
 ##

 ## Check if the pheno_data in the pheno_clean is declared model fit
 #if(attr(pheno_clean[[1]], "cleared")!="for_model_fit" && all(class(pheno_clean[[1]])!=c("data.frame"))) {
 if(attr(pheno_clean[["pheno_clean_data"]], "cleared")!="for_model_fit") {
     stop(paste(msg,'pheno_data is not phenotype data'), call. = FALSE)

 }

 response_family <- gp_resolve_response_family_for_responses(
   pheno_data = pheno_clean[["pheno_clean_data"]],
   response = response,
   response_family = response_family
 )
 for (mod in unique(stats::na.omit(c(GS_model, GS_model_cv)))) {
   gp_validate_model_response_family(
     model = mod,
     response_family = response_family
   )
 }
 positive_class_config <- gp_configure_positive_class_responses(
   pheno_data = pheno_clean[["pheno_clean_data"]],
   response = response,
   response_family = response_family,
   positive_class = positive_class
 )
 pheno_clean[["pheno_clean_data"]] <- positive_class_config[["pheno_data"]]
 positive_class <- positive_class_config[["positive_class"]]
 if (isTRUE(cross_validation)) {
   cv_metric_selection <- gp_prepare_cv_metric_selection(
     eval_metrics = eval_metrics,
     metric_for_ranking = metric_for_ranking,
     ranking_tie_breakers = ranking_tie_breakers,
     response_family = response_family
   )
   eval_metrics <- cv_metric_selection[["eval_metrics"]]
   metric_for_ranking <- cv_metric_selection[["metric_for_ranking"]]
   ranking_tie_breakers <- cv_metric_selection[["ranking_tie_breakers"]]
 } else if (!is.null(eval_metrics)) {
   gp_validate_eval_metrics(
     eval_metrics = eval_metrics,
     response_family = response_family
   )
 }
 gp_validate_stage2_weight_route(as.list(environment()))

 input_guardrail_ctx <- gp_apply_model_input_guardrails(as.list(environment()))
 pheno_clean <- input_guardrail_ctx[["pheno_clean"]]
 gmatrix <- input_guardrail_ctx[["gmatrix"]]
 gkernel <- input_guardrail_ctx[["gkernel"]]
 female_gmatrix <- input_guardrail_ctx[["female_gmatrix"]]
 male_gmatrix <- input_guardrail_ctx[["male_gmatrix"]]
 pedigree_matrix <- input_guardrail_ctx[["pedigree_matrix"]]
 omic1_kernel <- input_guardrail_ctx[["omic1_kernel"]]
 omic2_kernel <- input_guardrail_ctx[["omic2_kernel"]]
 omic3_kernel <- input_guardrail_ctx[["omic3_kernel"]]
 kernel_list <- input_guardrail_ctx[["kernel_list"]]
 env_similarity <- input_guardrail_ctx[["env_similarity"]]
 # Phase 3.16: write back auto-promoted fields from the guardrail ctx so
 # the downstream CV / true-prediction dispatch reads the flipped values
 # (met_ml_dl from gp_auto_promote_met_ml_dl_flag, gmatrix_method from
 # gp_auto_promote_met_ml_dl_kernel_build). Without these, the gates at
 # main_crossvalidation_execution_logic.R:645/785/943 and
 # model_execute_pipeline_helpers.R:3358 see the un-flipped scope-local
 # variables and silently fall through to single-env paths.
 met_ml_dl <- isTRUE(input_guardrail_ctx[["met_ml_dl"]])
 gmatrix_method <- input_guardrail_ctx[["gmatrix_method"]] %||% gmatrix_method

user_defined_model <- gp_shared_preprocess_requested_models(as.list(environment()))
run_profile <- gp_runtime_profile_new("model_execute")

pre_input_route <- gp_route_pre_input_specialized_models(as.list(environment()))
if (!is.null(pre_input_route)) {
  return(pre_input_route)
}
if (sum(isTRUE(multi_trait_asreml), isTRUE(multi_trait_gp), isTRUE(multi_trait_bayes), isTRUE(multi_trait_ml), isTRUE(multi_trait_dl), isTRUE(hybrid_asreml), isTRUE(hybrid_bayes), isTRUE(hybrid_gp), isTRUE(hybrid_ml), isTRUE(hybrid_dl)) > 1L) {
  stop("multi_trait_asreml, multi_trait_gp, multi_trait_bayes, multi_trait_ml, multi_trait_dl, hybrid_asreml, hybrid_bayes, hybrid_gp, hybrid_ml, and hybrid_dl are separate paths; enable only one at a time.", call. = FALSE)
}

valid_models <- gp_kernel_models_for_input_preparation(
  base_models = kernel_relationship_models,
  multi_trait_gp = multi_trait_gp
)

# MET ML/DL builds its features from the relationship matrix, so it keeps the
# user's gmatrix_method / kernel_method.
if (!any(user_defined_model %in% valid_models) && !isTRUE(met_ml_dl)) {
      kernel_method <- NULL
      gmatrix_method <- NULL
    }
 #### Get the clean geno_data ready for model fit
 ## The geno_to_model function depend on geno_precheck function. The expected
 ## output is clean genomic data with no missing and all QC control is checked.
 ## For details check geno_to_model and geno_precheck function description
 ##
 #if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){

gp_cv_point_only_kernel_qc <- isTRUE(cross_validation) &&
  isTRUE(cv_evaluation_only) &&
  !isTRUE(predictpror_shared_final_prediction) &&
  length(user_defined_model) > 0L &&
  all(as.character(user_defined_model) %in% gp_valid_models)
if (isTRUE(gp_cv_point_only_kernel_qc)) {
  if (identical(kernel_check_level, "auto")) {
    kernel_check_level <- "light"
  }
  if (identical(kernel_rcn_check, "auto")) {
    kernel_rcn_check <- "skip"
  }
  if (identical(kernel_pd_check, "auto")) {
    kernel_pd_check <- "sample"
  }
}

model_input_stage <- "model_input_objects"
if (gp_shared_preprocess_cache_has(
  predictpror_preprocess_cache,
  model_input_stage
)) {
  model_input_objects <- gp_shared_preprocess_cache_get(
    predictpror_preprocess_cache,
    model_input_stage
  )
} else {
  model_input_objects <- gp_prepare_model_input_objects(as.list(environment()))
  gp_shared_preprocess_cache_set(
    predictpror_preprocess_cache,
    model_input_stage,
    model_input_objects
  )
}
run_profile <- gp_runtime_profile_mark(
  run_profile,
  "prepare_model_inputs",
  detail = paste("responses", length(response %||% character()))
)
low_call_rate_inds_removed <- model_input_objects[["low_call_rate_inds_removed"]]
pheno_clean <- model_input_objects[["pheno_clean"]]
geno_res <- model_input_objects[["geno_res"]]
omic1_res <- model_input_objects[["omic1_res"]]
omic2_res <- model_input_objects[["omic2_res"]]
omic3_res <- model_input_objects[["omic3_res"]]
geno_omic_model_ready_list <- model_input_objects[["geno_omic_model_ready_list"]]
gmatrix_kernel_model_ready_list <- model_input_objects[["gmatrix_kernel_model_ready_list"]]
test_set <- gp_merge_test_sets(model_input_objects[["test_set"]], pheno_clean[["test_set"]])
if (!is.null(test_set)) {
  pheno_clean[["test_set"]] <- test_set
}
ml_dat_res <- model_input_objects[["ml_dat_res"]]
if (!is.null(feature_selected) && !is.null(ml_dat_res)) {
  ml_dat_res <- gp_feature_apply_explicit_to_ml_dat_res(
    ml_dat_res = ml_dat_res,
    feature_selected = feature_selected,
    trait = NULL
  )
}

feature_scoring_cv <- match.arg(
  gp_feature_normalize_cv_policy(feature_scoring_cv %||% "fixed"),
  choices = c("fixed", "fold_internal", "both")
)
feature_scoring_active <- isTRUE(feature_scoring) ||
  !is.null(feature_score_metadata) ||
  !is.null(feature_k_grid) ||
  !is.null(feature_k)
feature_score_metadata <- feature_score_metadata %||% NULL
feature_k_cv_summary <- NULL
feature_source_map <- NULL
if (isTRUE(feature_scoring_active) && !is.null(ml_dat_res) && !is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
  feature_scoring_model <- match.arg(
    as.character(feature_scoring_model %||% "Ridge_Regression"),
    choices = c("Ridge_Regression", "BayesB", "RandomForest")
  )
  feature_scoring_seed <- as.integer(feature_scoring_seed %||% random_state %||% random_seed %||% 123L)
  feature_matrix_for_scoring <- ml_dat_res[["merged_data"]][["merge_data"]]
  feature_source_map <- gp_feature_source_map_from_inputs(
    predictor_data = feature_matrix_for_scoring,
    source_matrices = list(
      geno_data = geno_omic_model_ready_list[["geno_model_ready"]] %||% NULL,
      omic1_data = geno_omic_model_ready_list[["omic1_model_ready"]] %||% NULL,
      omic2_data = geno_omic_model_ready_list[["omic2_model_ready"]] %||% NULL,
      omic3_data = geno_omic_model_ready_list[["omic3_model_ready"]] %||% NULL
    )
  )
  if (is.null(feature_score_metadata)) {
    feature_score_metadata <- feature_score_predictors(
      predictor_data = feature_matrix_for_scoring,
      pheno_data = pheno_clean[["pheno_clean_data"]],
      response = response,
      gen_name = gen_name,
      scoring_model = feature_scoring_model,
      seed = feature_scoring_seed,
      source_block = feature_source_map,
      ridge_lambda = feature_ridge_lambda,
      bayes_nIter = feature_bayes_nIter,
      bayes_burnIn = feature_bayes_burnIn,
      bayes_thin = feature_bayes_thin,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs,
      response_family = response_family,
      selection_mode = if (isTRUE(cross_validation) && identical(feature_scoring_cv, "fold_internal")) {
        "final_refit"
      } else {
        "fixed"
      }
    )
  }
  if (isTRUE(cross_validation)) {
    specialized_feature_route <- any(vapply(
      list(hybrid_asreml, hybrid_bayes, hybrid_gp, hybrid_ml, hybrid_dl,
           multi_trait_asreml, multi_trait_gp, multi_trait_bayes, multi_trait_ml, multi_trait_dl),
      isTRUE, logical(1)
    ))
    feature_k_grid <- gp_feature_k_grid(
      feature_k_grid,
      ncol(feature_matrix_for_scoring),
      include_all = !specialized_feature_route
    )
  } else if (is.null(feature_k)) {
    feature_k <- "all"
  }
  run_profile <- gp_runtime_profile_mark(
    run_profile,
    "feature_scoring",
    detail = paste(
      "model", feature_scoring_model,
      "traits", length(unique(feature_score_metadata$trait %||% character()))
    )
  )
}

 ### Ends

 ###################Genetic space test

 # if ("test_set" %in% names(pheno_clean) && isTRUE(test_train_genetic_space)) {
 #   test_set <- pheno_clean[["test_set"]]
 #
 #   test_set <- as.character(unique(test_set))
 #
 #   datasets_index_kernel <- which(!sapply(gmatrix_kernel_model_ready_list, is.null))
 #   datasets_index_geno_omic <- which(!sapply(geno_omic_model_ready_list, is.null))
 #
 #   if (length(datasets_index_geno_omic) != 0) {
 #     M <- geno_omic_model_ready_list[datasets_index_geno_omic]
 #     M <- do.call(rbind, M)
 #
 #     sik <- tryCatch({
 #       evaluate_genetic_space(M = M,
 #                              train_ids = setdiff(rownames(M), test_set),
 #                              test_ids = test_set)
 #     }, error = function(e) {
 #       message("Error in evaluating genetic space with geno and/or omic data: ", e$message)
 #       NULL
 #     })
 #   } else if (length(datasets_index_kernel) != 0) {
 #     M <- gmatrix_kernel_model_ready_list[[1]]
 #     sik <- tryCatch({
 #       evaluate_genetic_space(M = M,
 #                              train_ids = setdiff(rownames(M), test_set),
 #                              test_ids = test_set)
 #     }, error = function(e) {
 #       message("Error in evaluating genetic space with kernel data: ", e$message)
 #       NULL
 #     })
 #   } else {
 #     sik <- NULL
 #   }
 #
 #   if (!is.null(sik)) {
 #     tryCatch({
 #       if (!is.null(sik$plot)) {
 #         print(sik$plot)
 #         genetic_space_recommendation(recommendation = sik$recommendation)
 #       }
 #     }, error = function(e) {
 #       message("Error processing 'sik': ", e$message)
 #     })
 #   }
 # }

 #######################
#################################################################
############# Cross-Validation Start

 best_models <-  NULL
 best_models_ggplot_rep <- NULL
 best_models_ggplot_mean <- NULL
 cv_results_processed <-  NULL
 cv_results_raw <- NULL
 cv_results <- NULL
 cv_run_metadata <- NULL
 model_prep_all_bayes_cv <-  NULL
 asreml_models_prep_cv <- NULL
 res_plot_result_diagnostic <-  NULL
 res_mod_results_cv_per_trait_model <-  NULL
 cv_results_predicted_vs_observed <- NULL
 res_plot_result_diagnostic_cv_only <- NULL
 diagnostic_plots <- NULL

 ####
 test_set <- gp_merge_test_sets(test_set, pheno_clean[["test_set"]])
 if (!is.null(test_set)) {
   pheno_clean[["test_set"]] <- test_set
 }

if(docker_nd_usage) sys_name <- "Windows"

if(isTRUE(cross_validation)){
  cv_special_route <- gp_route_cross_validation_specialized_models(as.list(environment()))
  if (!is.null(cv_special_route)) {
    if (is.list(cv_special_route) && !is.null(hybrid_parent_qc_summary)) {
      cv_special_route$hybrid_parent_qc <- hybrid_parent_qc_summary
    }
    return(cv_special_route)
  }
  cv_pipeline <- tryCatch(
    gp_try_cv_frontdoor_fast_lane(as.list(environment())),
    error = function(e) {
      warning(
        "GP front-door CV fast lane failed; falling back to the standard CV path: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )
  if (!is.null(cv_pipeline)) {
    run_profile <- gp_runtime_profile_mark(
      run_profile,
      "gp_cv_frontdoor_fast_lane",
      detail = paste("results", length(cv_pipeline[["cv_results"]] %||% list()))
    )
  } else {
    cv_artifacts <- gp_prepare_cross_validation_artifacts(as.list(environment()))
    run_profile <- gp_runtime_profile_mark(
      run_profile,
      "prepare_cv_artifacts",
      detail = paste("models", length(GS_model_cv %||% character()))
    )
    pheno_data <- cv_artifacts[["pheno_data"]]
    model_prep_all_bayes_cv <- cv_artifacts[["model_prep_all_bayes_cv"]]
    asreml_models_prep_cv <- cv_artifacts[["asreml_models_prep_cv"]]

    cv_pipeline <- gp_run_cross_validation_pipeline(as.list(environment()))
    run_profile <- gp_runtime_profile_mark(
      run_profile,
      "cross_validation_pipeline",
      detail = paste("results", length(cv_pipeline[["cv_results"]] %||% list()))
    )
  }
  if (is.null(cv_pipeline)) {
     return(NULL)
   }

   cv_results <- cv_pipeline[["cv_results"]]
   cv_results_processed <- cv_pipeline[["cv_results_processed"]]
   cv_results_predicted_vs_observed <- cv_pipeline[["cv_results_predicted_vs_observed"]]
   cv_run_metadata <- cv_pipeline[["run_metadata"]]
   best_models <- cv_pipeline[["best_models"]]
   best_models_ggplot_rep <- cv_pipeline[["best_models_ggplot_rep"]]
   best_models_ggplot_mean <- cv_pipeline[["best_models_ggplot_mean"]]

 }

 if(!is.null(geno_data_process)){
   geno_qc_stat <-  geno_data_process[["qc_metrics_and_summary_stat"]]
   rm(geno_data_process); gc()
 } else{
 geno_qc_stat <- if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL

 }


 if(isTRUE(cv_evaluation_only) && isTRUE(cross_validation)){

   cv_only_result <- results_handling(GS_model = GS_model_cv,
                                      cv_results_raw = cv_results,
                                      res_model_output =  NULL,
                                      res_summary_stat =  NULL,
                                      res_plot = best_models_ggplot_rep,
                                      res_plot_mean = best_models_ggplot_mean,
                                      res_plot_result_diagnostic = NULL,
                                      test_diagonistic_plots = NULL,
                                      res_plot_result_diagnostic_cv_only = cv_results_predicted_vs_observed$predicted_vs_observed_plots,
                                      res_mod_results_cv_per_trait_model = cv_results_predicted_vs_observed$mod_res_per_trait_per_model,
                                      geno_qc_stat = NULL,
                                      cv_results_processed = cv_results_processed,
                                      run_metadata = cv_run_metadata,
                                      system_database = system_database,
                                      plot_filename = "CV_results",
                                      plot_extension = plot_extension,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      plot_units = plot_units,
                                      plot_dpi = plot_dpi,
                                      feature_selected = feature_selected,
                                      feature_score_metadata = feature_score_metadata)
   cv_only_result$cv_results_predicted_vs_observed <- cv_results_predicted_vs_observed
   return(cv_only_result)

 }

 if(is.null(best_models) & isTRUE(cross_validation)){
   cat(paste(msg, "There is a problem with the cross-vlidation process. Check the data and the model.\n"))
   return(NULL)
 }

future::plan("sequential")
###############################################################

 true_prediction_route <- gp_route_true_prediction_specialized_models(as.list(environment()))
 if (!is.null(true_prediction_route)) {
   if (is.list(true_prediction_route) && !is.null(hybrid_parent_qc_summary)) {
     true_prediction_route$hybrid_parent_qc <- hybrid_parent_qc_summary
   }
   return(true_prediction_route)
 }
 cv_best_models <- best_models
 if(is.null(best_models)){
   if (length(GS_model) > 1) {
     if (length(GS_model) != length(response)) {
       stop(paste(msg, "When the number of models is more than one, the number of models should be the same as the number of traits."), call. = FALSE)
     }
   }

   best_models <- data.frame(model = GS_model, trait = response, stringsAsFactors = FALSE)
 }

  best_models <- gp_build_true_prediction_task_table(
    best_models = best_models,
    GS_model = GS_model,
    GS_model_cv = GS_model_cv,
    response = response,
   cv_results_processed = cv_results_processed,
    metric_for_ranking = metric_for_ranking,
    cross_validation = cross_validation
  )
  cv_best_models <- best_models
  n_trait <- length(unique(best_models[["trait"]] %||% character()))
  n_model <- length(unique(best_models[["model"]] %||% character()))

 sys_name <- Sys.info()["sysname"]
 if(docker_nd_usage) sys_name <- "Windows"


 # -------- TRUE PREDICTION: smart parallel, python-safe, no DP renaming --------

 # Helpers (lightweight, keep near this block)
 `%||%` <- function(a, b) if (is.null(a)) b else a

 # Canonical names so we can safely test for DP models etc.


 # (Re)canonicalize best_models$model for policy checks & DP routing
 best_models$model <- gp_canonicalize_supported_model_names(best_models$model)

 friendly_name_lookup <- gp_model_friendly_lookup()

 # Reticulate/Python init is not used for direct CLI ML/DL/GP model paths.
 gp_py_bin <- gp_detect_python()
 direct_cli_prediction_models <- unique(c(AI_valid_models, gp_valid_models, gp_parallel_r_only_models()))
 init_py <- function(model = NULL) {
   if (!is.null(model) && as.character(model) %in% direct_cli_prediction_models) {
     return(invisible(TRUE))
   }
   gp_init_python_once(gp_py_bin)
 }

 # One runner per best-model row (single source of truth)
 outer_fixed <- fixed
 # Keep the complete prepared kernel bank in the task-local context.  The
 # individual canonical matrices below are insufficient for user-supplied
 # `kernel_list` entries, and parallel/task serialization does not guarantee
 # that an unreferenced parent-frame object remains available.
 outer_kernel_bank <- gmatrix_kernel_model_ready_list
 run_one_best <- function(i) {
   task_row <- best_models[i, ]
  response <- as.character(task_row$trait)
  GS_model <- as.character(task_row$model)   # keep canonical if DP
  task_feature_k <- if ("feature_k" %in% names(task_row) && !is.na(task_row$feature_k)) {
    as.integer(task_row$feature_k)
  } else {
    NULL
  }
  feature_k <- feature_k %||% task_feature_k
  init_py(GS_model)

   # Carry outer execution settings into the task environment so downstream
   # helpers receive the same scalar configuration used by model_execute().
   engine <- engine
   gen_name <- gen_name
   pheno_clean <- pheno_clean
   fixed <- outer_fixed
   random <- random
   cova <- cova
   heter_groups <- heter_groups
   heter_resid <- heter_resid
   bayes_kernel_heter_resid <- bayes_kernel_heter_resid
   var_cov_str <- var_cov_str
   weights <- weights
   pworkspace <- pworkspace
   workspace <- workspace
   maxit <- maxit
   inverse <- inverse
   epsilon <- epsilon
   message <- message
   eval_metrics <- eval_metrics
   system_database <- system_database
   msg <- msg
   CI_width_thresholds <- CI_width_thresholds
   confidence_level <- confidence_level
   high_reliability_thres <- high_reliability_thres
   low_reliability_thres <- low_reliability_thres
   friendly_name_lookup <- friendly_name_lookup
   gmatrix_kernel_model_ready_list <- outer_kernel_bank

   # Deterministic per-task seed derived from the public random_state.  Keeping
   # the task offset preserves stable seeds within a multi-model run, while
   # allowing independent Bayesian chains and reproducible alternative runs.
   new_seed <- gp_model_task_seed(random_state, i)
   # Pin the generator: parallel workers run L'Ecuyer-CMRG, under which the
   # same seed gives a different MCMC chain than in the main process.
   old_rng_kind_best <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
   on.exit(do.call(RNGkind, as.list(old_rng_kind_best)), add = TRUE)
   gp_set_seed(new_seed)
   options(random_state = new_seed)

   # Choose CI/calculation mode (Bayes vs ML) without mutating GS_model label
   model_for_CI_cal <- if (!GS_model %in% AI_valid_models) "Bayes" else "ML"

   # Enforce single-env constraint for ML/DP models (your original condition was inverted)
   if (GS_model %in% AI_valid_models) {
     one_per_id <- (length(pheno_clean[["pheno_clean_data"]][, gen_name]) ==
                      length(unique(pheno_clean[["pheno_clean_data"]][, gen_name])))
     if (!one_per_id && !isTRUE(met_ml_dl && gp_is_met_ml_dl_model(GS_model))) {
       stop(paste(msg, GS_model, "only works for single location/environment."), call. = FALSE)
     }
   }

   # -------------------- BAYES block (A/B/C/BL/BRR) --------------------
   # bayes_valid_models       <- c("BRR","BayesA","BayesB","BayesC","BL")
   # bayes_gblup_valid_models <- c("GBLUP_BRR","RKHS")

   res_model_output <- NULL
   res_summary_stat <- NULL
   geno_model_ready <- geno_omic_model_ready_list[["geno_model_ready"]] %||% NULL
   omic1_model_ready <- geno_omic_model_ready_list[["omic1_model_ready"]] %||% NULL
   omic2_model_ready <- geno_omic_model_ready_list[["omic2_model_ready"]] %||% NULL
   omic3_model_ready <- geno_omic_model_ready_list[["omic3_model_ready"]] %||% NULL
   gmatrix_model_ready <- gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] %||% NULL
   omic1_kernel_model_ready <- gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] %||% NULL
   omic2_kernel_model_ready <- gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] %||% NULL
   omic3_kernel_model_ready <- gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] %||% NULL
    task_ctx <- gp_merge_task_context(environment(run_one_best), environment())
      task_result <- gp_run_best_model_task(task_ctx)
   task_result$trait <- response
   task_result$response <- response
   task_result$model <- GS_model
   task_result$model_public <- gp_public_model_label(GS_model)[1]
   task_result$model_key <- make.names(task_result$model_public)
   task_result$is_cv_best <- isTRUE(task_row$is_cv_best)
   task_result$geno_qc_stat <- geno_qc_stat
   task_result
 }

results <- gp_execute_best_model_tasks(
  best_models = best_models,
  run_one_best = run_one_best,
  sequential_models = sequential_models,
  canonical_names = canonical_names,
  friendly_names = friendly_names,
  policy_param_source = as.list(environment()),
  parallel_mode = parallel_mode,
  num_cores = num_cores,
  globals_max_GB = globals_max_GB,
  verbose = verbose,
  parallel_backend_prefer_fork = parallel_backend_prefer_fork,
  pheno_clean = pheno_clean,
  ml_dat_res = ml_dat_res,
  gmatrix_kernel_model_ready_list = gmatrix_kernel_model_ready_list,
  geno_omic_model_ready_list = geno_omic_model_ready_list,
  response = response,
  init_py = init_py,
  worker_memory_gb = worker_memory_gb,
  memory_budget_gb = memory_budget_gb
)
run_profile <- gp_runtime_profile_mark(
  run_profile,
  "true_prediction_tasks",
  detail = paste("tasks", nrow(best_models %||% data.frame()))
)

metadata_python_purpose <- gp_python_purpose_for_models(
  best_models$model %||% GS_model
)
metadata_python <- if (metadata_python_purpose %in% c("dl", "gp")) {
  gp_preferred_python(purpose = metadata_python_purpose) %||% gp_py_bin
} else {
  gp_py_bin
}
run_metadata <- gp_runtime_metadata(
  execution_policy = attr(results, "gp_execution_policy"),
  python_path = metadata_python,
  preferred_python = gp_preferred_python(purpose = metadata_python_purpose),
  purpose = metadata_python_purpose,
  context = "model_execute",
  extra_fields = c(
    gp_runtime_profile_metadata_fields(run_profile),
    list(
      task_unit = if (isTRUE(cross_validation) && length(GS_model_cv %||% character()) > 0L) "cv_best_model_trait" else "model_trait",
      task_models = length(unique(best_models$model %||% character())),
      task_traits = length(unique(best_models$trait %||% character())),
      cv_best_models_retained = paste(
        paste(cv_best_models$trait %||% character(), cv_best_models$model %||% character(), sep = ":"),
        collapse = ","
      )
    )
  )
)

# keep names aligned with traits
if (!is.null(best_models) && nrow(best_models) > 0) {
   if ("task_name" %in% names(best_models)) {
     names(results) <- best_models[["task_name"]]
   } else {
     names(results) <- make.unique(make.names(paste(best_models[["trait"]], best_models[["model"]], sep = "_")))
   }
 } else {
   names(results) <- response
 }
 # ---------------------------------------------------------------------------

 final_results <- gp_finalize_model_execute_results(
   results = results,
   GS_model = GS_model,
   gen_name = gen_name,
  heter_groups = heter_groups,
  pheno_data = pheno_clean[["pheno_clean_data"]],
  response_family = response_family,
   best_models_ggplot_rep = best_models_ggplot_rep,
   best_models_ggplot_mean = best_models_ggplot_mean,
   cv_results_predicted_vs_observed = cv_results_predicted_vs_observed,
   geno_qc_stat = geno_qc_stat,
   cv_results_processed = cv_results_processed,
   cv_results = cv_results,
   feature_selected = feature_selected,
   feature_score_metadata = feature_score_metadata,
   run_metadata = run_metadata,
   system_database = system_database,
   plot_extension = plot_extension,
   plot_width = plot_width,
   plot_height = plot_height,
   plot_units = plot_units,
   plot_dpi = plot_dpi
 )

 if (is.null(final_results[["model_results"]]) &&
     is.list(final_results[["model_results_by_model"]]) &&
     length(final_results[["model_results_by_model"]]) == 1L) {
   single_model_block <- final_results[["model_results_by_model"]][[1L]]
   if (is.list(single_model_block) && length(single_model_block) == 1L) {
     final_results[["model_results"]] <- single_model_block[[1L]]
   }
 }

 # Phase 3.16: post-process to add a consistent env correlation +
 # env covariance pair for every model that produced per-env predictions.
 # GP / Bayes / ASReml keep their native env covariance; ML/DL get NA
 # covariance + empirical correlation, so user code that reads
 # sik$model_results$environment_correlation or
 # sik$cv_results_processed$environment_correlation works uniformly.
 final_results <- tryCatch(
   gp_add_consistent_met_summaries(
     final_results,
     ctx = list(heter_groups = heter_groups, gen_name = gen_name,
                cross_validation = cross_validation,
                heter_resid = heter_resid,
                bayes_kernel_heter_resid = bayes_kernel_heter_resid,
                pheno_data = pheno_data)
   ),
   error = function(e) {
     warning("gp_add_consistent_met_summaries failed (non-fatal): ",
             conditionMessage(e), call. = FALSE)
     final_results
   }
 )

 return(final_results)

} ## end of function


