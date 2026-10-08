# =============================================================================
# Scenario 9: raw genotype files (VCF, HapMap, CSV/TXT) and Beagle imputation
# =============================================================================
# Instead of a 0/1/2 matrix, give a genotype file. PredictProR reads it, runs
# genotype QC (call rate, MAF, heterozygosity), imputes, and recodes to 0/1/2
# once; the matrix step then skips a second QC.
#   VCF          : vcf_file_name = "x.vcf", vcf_file_path = "<folder>"
#   HapMap       : hapmap_file_name = "x.hmp.txt", hapmap_file_path = "<folder>",
#                  or an in-memory table: hapmap = <data.frame>
#   CSV/TXT table: csv_file_name = "x.csv", csv_file_path = "<folder>"
#                  (VCF marker columns CHROM, POS, ID, REF, ALT, QUAL, FILTER,
#                  INFO, FORMAT, then one column per sample, 0/1/2 dosages)
# Sample names in the file must match gen_name in pheno; samples without a
# phenotype are predicted.
# Beagle: impute = TRUE, imputation_method = "beagle" (needs Java; the Beagle
# JAR is downloaded once). HapMap and CSV/TXT are converted to VCF for Beagle.
# Standalone: vcf_qc_recode(), hmp_qc_recode(), impute_genotypes_with_beagle(),
# convert_csv_to_vcf(), convert_hapmap_to_vcf().
# =============================================================================
source(system.file("examples", "scenarios", "00_setup.R", package = "PredictProR"))

# ---- Data: write the simulated genotypes as files ----------------------------------
dat <- sim_single_env(n = 60, p = 200)
pheno <- dat$pheno
response <- "Yield"
ids <- rownames(dat$geno)
folder <- file.path(tempdir(), "PredictProR_genotype_files")
dir.create(folder, showWarnings = FALSE)

write_vcf <- function(dosage, file) {              # dosage: markers x samples
  gt <- matrix(c("0" = "0/0", "1" = "0/1", "2" = "1/1")[as.character(dosage)], nrow = nrow(dosage))
  gt[is.na(dosage)] <- "./."
  body <- vapply(seq_len(nrow(dosage)), function(j) {
    paste(c("1", j * 1000L, rownames(dosage)[j], "A", "G", ".", "PASS", ".", "GT", gt[j, ]), collapse = "\t")
  }, character(1))
  writeLines(c("##fileformat=VCFv4.3",
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(dosage)),
                     collapse = "\t"), body), file)
}
dosage <- t(dat$geno)                                   # markers x samples
set.seed(9)
dosage_missing <- dosage
dosage_missing[stats::runif(length(dosage)) < 0.05] <- NA   # 5% missing calls

write_vcf(dosage, file.path(folder, "geno.vcf"))
write_vcf(dosage_missing, file.path(folder, "geno_missing.vcf"))

# HapMap: alleles A/G; genotype codes A = 0, R = heterozygote, G = 2, N = missing
hmp_codes <- matrix(c("0" = "A", "1" = "R", "2" = "G")[as.character(dosage_missing)], nrow = nrow(dosage))
hmp_codes[is.na(dosage_missing)] <- "N"
colnames(hmp_codes) <- ids
hapmap_table <- data.frame(`rs#` = rownames(dosage), alleles = "A/G", chrom = "1",
                           pos = seq_len(nrow(dosage)) * 1000L, strand = "+", `assembly#` = NA, center = NA,
                           protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA, hmp_codes,
                           check.names = FALSE, stringsAsFactors = FALSE)
utils::write.table(hapmap_table, file.path(folder, "geno_missing.hmp.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# CSV / TXT genotype table: VCF marker columns + one 0/1/2 column per sample
geno_table <- data.frame(CHROM = "1", POS = seq_len(nrow(dosage)) * 1000L, ID = rownames(dosage),
                         REF = "A", ALT = "G", QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
                         dosage_missing, check.names = FALSE, stringsAsFactors = FALSE)
utils::write.csv(geno_table, file.path(folder, "geno_missing.csv"), row.names = FALSE, na = ".")
utils::write.table(geno_table, file.path(folder, "geno_missing.txt"), sep = "\t", row.names = FALSE,
                   quote = FALSE, na = ".")

# ---- 1. model_execute() reading each file type -----------------------------------
file_inputs <- list(
  list(name = "VCF file", vcf = "geno.vcf", impute_method = "knn", needs = character()),
  list(name = "VCF file, stricter QC (MAF 5%)", vcf = "geno.vcf", impute_method = "knn", maf = 0.05, needs = character()),
  list(name = "HapMap table in memory", hapmap = hapmap_table, impute_method = "knn", needs = character()),
  list(name = "HapMap .txt file", hapmap_file = "geno_missing.hmp.txt", impute_method = "knn", needs = character()),
  list(name = "CSV table", csv = "geno_missing.csv", impute_method = "mean", needs = character()),
  list(name = "TXT table", csv = "geno_missing.txt", impute_method = "mean", needs = character()),
  list(name = "VCF with missing calls + Beagle", vcf = "geno_missing.vcf", impute_method = "beagle", needs = "java"),
  list(name = "HapMap .txt with missing calls + Beagle", hapmap_file = "geno_missing.hmp.txt", impute_method = "beagle", needs = "java"),
  list(name = "CSV table with missing calls + Beagle", csv = "geno_missing.csv", impute_method = "beagle", needs = "java")
)

for (cfg in file_inputs) {
  label <- paste("File input:", cfg$name, "- GBLUP_BRR")
  if (!is.null(need <- missing_software(cfg[["needs"]]))) { log_skip(label, need); next }

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data: phenotypes + ONE raw genotype source
      pheno_data = pheno,
      geno_data = NULL,
      vcf_file_name = cfg[["vcf"]],
      vcf_file_path = if (!is.null(cfg[["vcf"]])) folder else NULL,
      hapmap = cfg[["hapmap"]],
      hapmap_file_name = cfg[["hapmap_file"]],
      hapmap_file_path = if (!is.null(cfg[["hapmap_file"]])) folder else NULL,
      csv_file_name = cfg[["csv"]],
      csv_file_path = if (!is.null(cfg[["csv"]])) folder else NULL,
      csv_input_coding = "alt_dosage",

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Explicit single-trait / single-environment mode
      heter_groups = NULL,

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = "GBLUP_BRR",

      # Relationship / kernel construction
      gmatrix_method = "VanRaden",

      # Genotype QC, imputation and recoding (done once, on the file)
      qc_filtering = TRUE,
      maf_threshold = if (is.null(cfg[["maf"]])) 0.01 else cfg[["maf"]],
      het_threshold = 0.10,
      ind_call_rate_threshold = 0.90,
      snp_call_rate_threshold = 0.90,
      impute = TRUE,
      imputation_method = cfg[["impute_method"]],   # "knn", "mean", "median", "mode" or "beagle"
      impute_knn_k = 5L,
      beagle_options = list(),                 # e.g. list(beagle_jar = "path/to/beagle.jar")
      recode_format = "0,1,2",
      ploidy = 2L,
      ld_prunning_qc = TRUE,
      ld_pruning = FALSE,                      # PLINK-based; use vcf_qc_recode(qc_engine = "plink")

      # No cross-validation
      cross_validation = FALSE,

      # Bayesian model controls
      nIter = 600L,                # use 16000
      burnIn = 200L,               # use 2000
      thin = 5L,

      # Memory and scheduling
      parallel_mode = "auto",
      parallel_backend_prefer_fork = FALSE,

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

# ---- 2. Standalone recoding / imputation, then model_execute() on the matrix -------
beagle_out <- file.path(folder, "beagle")
standalone <- list(
  list(name = "vcf_qc_recode() on VCF", needs = character(),
       run = function() vcf_qc_recode(vcf_file_name = "geno.vcf", vcf_file_path = folder, maf_threshold = 0.05,
                                      impute = FALSE, recode_format = "0,1,2")$snps_matrix),
  list(name = "impute_genotypes_with_beagle() on VCF", needs = "java",
       run = function() impute_genotypes_with_beagle(file.path(folder, "geno_missing.vcf"),
                                                     output_prefix = file.path(beagle_out, "vcf"))$genotypes$snps_matrix),
  list(name = "impute_genotypes_with_beagle() on HapMap .txt", needs = "java",
       run = function() impute_genotypes_with_beagle(file.path(folder, "geno_missing.hmp.txt"), input_format = "hapmap",
                                                     output_prefix = file.path(beagle_out, "hmp"))$genotypes$snps_matrix),
  list(name = "impute_genotypes_with_beagle() on CSV table", needs = "java",
       run = function() impute_genotypes_with_beagle(file.path(folder, "geno_missing.csv"),
                                                     output_prefix = file.path(beagle_out, "csv"))$genotypes$snps_matrix),
  list(name = "convert_csv_to_vcf() then vcf_qc_recode()", needs = character(),
       run = function() {
         convert_csv_to_vcf(file.path(folder, "geno_missing.csv"), file.path(folder, "from_csv.vcf"))
         vcf_qc_recode(vcf_file_name = "from_csv.vcf", vcf_file_path = folder, impute = TRUE,
                       imputation_method = "mean", recode_format = "0,1,2")$snps_matrix
       })
)

for (cfg in standalone) {
  label <- paste("Standalone:", cfg$name, "+ BayesB")
  if (!is.null(need <- missing_software(cfg[["needs"]]))) { log_skip(label, need); next }
  Geno.data <- try(cfg$run(), silent = TRUE)
  if (inherits(Geno.data, "try-error")) { log_result(label, Geno.data); next }
  cat(sprintf("%s: %d lines x %d markers, %d missing\n", cfg$name, nrow(Geno.data), ncol(Geno.data), sum(is.na(Geno.data))))

  model_runtime <- system.time({
    res <- try(PredictProR::model_execute(
      # Input data (recoded 0/1/2 matrix)
      pheno_data = pheno,
      geno_data = Geno.data,

      # Trait and genotype identifiers
      gen_name = "GID",
      response = response,
      response_family = "gaussian",

      # Mixed-model formulas
      fixed = ~1,
      random = ~GID,

      # Candidate model
      GS_model = "BayesB",

      # Genotype QC and preprocessing (already done on the file)
      qc_filtering = FALSE,
      impute = TRUE,
      imputation_method = "mean",
      ld_prunning_qc = FALSE,

      # No cross-validation
      cross_validation = FALSE,

      # Bayesian model controls
      nIter = 600L,
      burnIn = 200L,
      thin = 5L,

      # Memory and scheduling
      parallel_mode = "auto",
      parallel_backend_prefer_fork = FALSE,

      # Output
      system_database = TRUE,
      message = FALSE,
      verbose = FALSE
    ), silent = TRUE)
  })
  log_result(label, res, model_runtime)
}

scenario_report()
