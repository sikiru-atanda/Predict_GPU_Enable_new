# run_plink_qc <- function(input_file,
#                          output_name,
#                          output_format = NULL,
#                          vcf_file_path = NULL,
#                          min_version = "1.90",
#                          #plink_path = "D:/PredictProR/plink",
#                          remove_monomorphic = TRUE,
#                          maf_threshold = NULL,
#                          heterozygosity = 0.2,
#                          snp_call_rate = NULL,
#                          allow_extra_chr = FALSE,
#                          individual_call_rate = NULL,
#                          #recode = TRUE,
#                          recode_format = "0,1,2") {
#
#   # Check if PLINK is installed and accessible
#   # if (system("plink --version", ignore.stderr = TRUE) != 0) {
#   #   message("PLINK is not installed. Attempting to install PLINK...")
#     #plink_dir <- install_plink()  # Call to install PLINK
#     plink_dir <-  check_and_install_plink(min_version, vcf_file_path)
#   #   # Re-check if PLINK is now installed. If not, exit with an error message.
#   #   if (system("plink --version", ignore.stderr = TRUE) != 0) {
#   #     stop("Failed to install PLINK. Please install PLINK manually and ensure it is in your PATH.")
#   #   }
#   # }
#
#   plink_path <- paste(plink_dir, "plink.exe", sep = "/")
#   vcf_path <- paste(getwd(), input_file, sep = "/") # Make sure to replace 'your_file.vcf' with the actual file name
#   out_put_path <- paste(getwd(), output_name, sep = "/")
#   # Determine the PLINK executable path based on the operating system
#
#   if (.Platform$OS.type == "windows") {
#     plink_path <- paste(plink_dir, "plink.exe", sep = "/")
#   } else {
#     plink_path <- paste(plink_dir, "plink", sep = "/")
#   }
#
#   # Check for file extension
#   file_ext <- tools::file_ext(input_file)
#   if (file_ext == "gz") {
#     file_ext <- tools::file_ext(sub("\\.gz$", "", input_file))
#   }
#
#   # Determine PLINK input flag based on file type
#   # input_flag <- switch(file_type,
#   #                      "vcf" = "--vcf",
#   #                      "vcf.gz" = "--gzvcf",
#   #                      "bed" = "--bfile",
#   #                      stop("Unsupported file type."))
#   # Construct the base PLINK command according to the file extension
#   plink_cmd <- ""
#   if (file_ext %in% c("bed", "bim", "fam")) {
#     # Check if .bed file was provided, and infer base name for PLINK binary format
#     bed_base <- gsub(pattern = "\\.[bim|fam|bed]+$", replacement = "", x = input_file)
#     plink_cmd <- sprintf('%s --bfile "%s"', plink_path, bed_base)
#   } else if (file_ext == "vcf") {
#     plink_cmd <- sprintf('%s --vcf "%s"', plink_path, input_file)
#     #plink_cmd <- sprintf('%s --vcf "%s" --allow-extra-chr', plink_path, input_file)
#   } else {
#     stop("Unsupported file format. Please provide a .bed (with .bim and .fam) or .vcf file.")
#   }
#
#   # Start building the PLINK command
#   #plink_cmd <- paste("plink", input_flag, shQuote(input_path), "--out", shQuote(output_prefix))
#
#
#
#   # Add the --allow-extra-chr flag if specified
#   if (isTRUE(allow_extra_chr)) {
#
#     plink_cmd <- sprintf('%s --allow-extra-chr', plink_cmd)
#   }
#
#   # Initialize a default MAF threshold for removing monomorphic SNPs if no other MAF threshold is specified
#   default_maf_threshold_for_monomorphic <- 0.0001
#
#   # Update the PLINK command with the appropriate MAF threshold
#   if (isTRUE(remove_monomorphic) && is.null(maf_threshold)) {
#     # Apply the default threshold to remove monomorphic SNPs if no specific MAF threshold is provided
#     plink_cmd <- sprintf('%s --maf %f', plink_cmd, default_maf_threshold_for_monomorphic)
#   } else if (!is.null(maf_threshold) && isFALSE(remove_monomorphic)) {
#     # Use the user-specified MAF threshold, which will also inherently remove monomorphic SNPs
#     plink_cmd <- sprintf('%s --maf %f', plink_cmd, maf_threshold)
#   } else {
#     if (!is.null(maf_threshold) && isTRUE(remove_monomorphic)) {
#       # Use the user-specified MAF threshold, which will also inherently remove monomorphic SNPs
#       plink_cmd <- sprintf('%s --maf %f', plink_cmd, maf_threshold)
#     }
#   }
#
#   # Add heterozygosity limits if specified
#   if (!is.null(heterozygosity) && length(heterozygosity) == 2) {
#     plink_cmd <- sprintf('%s --hwe %f', plink_cmd, heterozygosity)
#     #plink_cmd <- paste(plink_cmd, "--hwe", paste(heterozygosity, collapse = " "))
#   }
#
#   # Add SNP call rate threshold if specified
#   if (!is.null(snp_call_rate)) {
#     plink_cmd <- sprintf('%s --geno %f', plink_cmd, snp_call_rate)
#
#     #plink_cmd <- paste(plink_cmd, "--geno", snp_call_rate)
#   }
#
#   # Add individual call rate threshold if specified
#   if (!is.null(individual_call_rate)) {
#     plink_cmd <- sprintf('%s --mind %f', plink_cmd, snp_call_rate)
#     #plink_cmd <- paste(plink_cmd, "--mind", individual_call_rate)
#   }
#
#   # Handle recoding options
#   # if (isTRUE(recode)) {
#   #   if (recode_format == "0,1,2") {
#   #     plink_cmd <- sprintf('%s --recode AD ... %f', plink_cmd, snp_call_rate)
#   #     #plink_cmd <- paste(plink_cmd, "--recode A")
#   #   }
#   # }
#
#
#   output_cmd <- switch(output_format,
#                        "bed" = "--make-bed",
#                        "vcf" = "--recode vcf",
#                        "vcf.gz" = "--recode vcf bgz",
#                        stop("Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."))
#   plink_cmd <- sprintf('%s %s --out "%s"', plink_cmd, output_cmd, output_name)
#
#   # Execute the command
#   message("Executing command: ", plink_cmd)
#   system(plink_cmd, intern = FALSE)
#   # Execute the PLINK command
#   # message("Executing PLINK command: ", plink_cmd)
#   # system(plink_cmd)
#   #
#   # message("PLINK QC and recoding complete. Output files are prefixed with '", output_name, "'")
#   #
# }
#
#
# sik = run_plink_qc(input_file = "sik.vcf.gz",
#             output_name = "todayy",
#             output_format = "vcf",
#             recode = TRUE,
#             remove_monomorphic = FALSE,
#             maf_threshold = 0.05,
#             heterozygosity = NULL,
#             snp_call_rate = 0.95,
#             allow_extra_chr = TRUE,
#             individual_call_rate = NULL)
#
#
# run_plink_qc(input_file = "FFAR_Chinese_Pea_ref_genome_filtered.vcf.gz",
#              output_name = "pea_whole_genome",
#              output_format = "vcf",
#              vcf_file_path = "D:/RICA",
#              #recode = TRUE,
#              remove_monomorphic = FALSE,
#              maf_threshold = 0.05,
#              heterozygosity = NULL,
#              snp_call_rate = 0.95,
#              allow_extra_chr = TRUE,
#              individual_call_rate = NULL)
#
#
#
# library(vcfR)
#
# setwd("D:/PredictProR")
# vcf <- read.vcfR("2021_NDSU_AYT.vcf", verbose = FALSE)
#
# vcf@gt
