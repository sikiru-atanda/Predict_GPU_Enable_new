% Generated package documentation placeholder; prefer roxygen source for future edits.
\name{mod_output_bayes_gbluBRR_RKHS}
\alias{mod_output_bayes_gbluBRR_RKHS}
\title{mod_output_bayes_gbluBRR_RKHS}
\usage{
mod_output_bayes_gbluBRR_RKHS(mod = NULL, ETA = NULL, gen_name = NULL, GS_model = NULL, 
    gkernel = NULL, gmatrix = NULL, omic1_kernel = NULL, omic2_kernel = NULL, 
    omic3_kernel = NULL, kernel_list = NULL, omics_kernel_label = list(omic1_kernel = NULL, 
        omic2_kernel = NULL, omic3_kernel = NULL), pheno_data = NULL, 
    heter_groups = NULL, bayes_para = NULL, CI_width_thresholds = c(0.33, 
        0.66), confidence_level = 0.95, high_reliability_thres = 0.9, 
    low_reliability_thres = 0.5, n_components = 20, threshold = 100, 
    target = "test_set", iqr_multiplier = 1.5, interval_width_high_threshold = NULL, 
    interval_width_low_threshold = NULL, interval_width_moderate_threshold = NULL, 
    system_database = FALSE, response_family = "gaussian", response = NULL, ...)
}
\arguments{
\item{mod}{Fitted model object.}

\item{ETA}{Parameter passed to the ETA workflow option.}

\item{gen_name}{Name of the genotype or sample identifier column.}

\item{GS_model}{Genomic selection model identifier.}

\item{gkernel}{Parameter passed to the gkernel workflow option.}

\item{gmatrix}{Parameter passed to the gmatrix workflow option.}

\item{omic1_kernel}{Parameter passed to the omic1_kernel workflow option.}

\item{omic2_kernel}{Parameter passed to the omic2_kernel workflow option.}

\item{omic3_kernel}{Parameter passed to the omic3_kernel workflow option.}

\item{kernel_list}{Optional named list of additional relationship or kernel matrices.}

\item{omics_kernel_label}{Parameter passed to the omics_kernel_label workflow option.}

\item{pheno_data}{Phenotype data frame.}

\item{response}{Response column name used to preserve original ordered class
levels for Bayesian classification output.}

\item{heter_groups}{Parameter passed to the heter_groups workflow option.}

\item{bayes_para}{Parameter passed to the bayes_para workflow option.}

\item{CI_width_thresholds}{Parameter passed to the CI_width_thresholds workflow option.}

\item{confidence_level}{Parameter passed to the confidence_level workflow option.}

\item{high_reliability_thres}{Parameter passed to the high_reliability_thres workflow option.}

\item{low_reliability_thres}{Parameter passed to the low_reliability_thres workflow option.}

\item{n_components}{Parameter passed to the n_components workflow option.}

\item{threshold}{Threshold value used by the workflow.}

\item{target}{Prediction target used for diagnostics.}

\item{iqr_multiplier}{Parameter passed to the iqr_multiplier workflow option.}

\item{interval_width_high_threshold}{Parameter passed to the interval_width_high_threshold workflow option.}

\item{interval_width_low_threshold}{Parameter passed to the interval_width_low_threshold workflow option.}

\item{interval_width_moderate_threshold}{Parameter passed to the interval_width_moderate_threshold workflow option.}

\item{system_database}{Logical indicating whether system database style outputs should be used.}

\item{response_family}{Response family used by the model workflow.}

\item{...}{Additional arguments passed to downstream helpers.}

}
\description{
Documentation for \code{mod_output_bayes_gbluBRR_RKHS}. This topic is generated to keep exported package APIs documented and synchronized with the current function signature.
}
