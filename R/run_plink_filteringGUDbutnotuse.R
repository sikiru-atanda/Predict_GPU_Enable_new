run_plink_filteringGUDbutnotuse <- function(input_file,
                                output_prefix,
                                snp_call_rate = NULL,
                                monomorphic = TRUE,
                                maf = NULL,
                                heterozygosity = NULL,
                                output_format = "bed",
                                plink_path = "plink",
                                allow_extra_chr = FALSE) {
  # Check for file extensions explicitly
  file_ext <- tools::file_ext(input_file)
  full_ext <- tolower(tools::file_path_sans_ext(basename(input_file)))

  # Construct the base PLINK command according to the file extension
  cmd <- ""

  if (grepl("\\.bed$", full_ext, ignore.case = TRUE)) {
    bed_base <- sub("\\.bed$", "", input_file)
    cmd <- sprintf('%s --bfile %s', plink_path, bed_base)
  } else if (file_ext == "vcf" || file_ext == "gz" && grepl("\\.vcf$", full_ext, ignore.case = TRUE)) {
    cmd <- sprintf('%s --vcf %s --allow-extra-chr', plink_path, input_file)
  } else {
    stop("Unsupported file format. Please provide a .bed (with .bim and .fam) or .vcf(.gz) file.")
  }

  # Optionally allow extra chromosomes
  if (isTRUE(allow_extra_chr)) {
    cmd <- sprintf('%s --allow-extra-chr', cmd)
  }

  # Adjust for handling monomorphic SNPs
  if(isTRUE(monomorphic)){
    cmd <- sprintf('%s --min-alleles 2 --max-alleles 2', cmd)

  }

  # Add flags for SNP call rate
  if (!is.null(snp_call_rate)) {
    cmd <- sprintf('%s --geno %f', cmd, snp_call_rate)
  }

  # Add flags for Minor Allele Frequency (MAF)
  if (!is.null(maf)) {
    cmd <- sprintf('%s --maf %f', cmd, maf)
  }

  # Add flags for heterozygosity
  if (!is.null(heterozygosity)) {
    cmd <- sprintf('%s --het --mind %f', cmd, heterozygosity)
  }

  # Specify output based on desired format
  if (output_format == "bed") {
    cmd <- sprintf('%s --make-bed --out %s', cmd, output_prefix)
  } else if (output_format == "vcf") {
    cmd <- sprintf('%s --recode vcf --out %s', cmd, output_prefix)
  } else if (output_format == "vcf.gz") {
    cmd <- sprintf('%s --recode vcf bgz --out %s', cmd, output_prefix)
  } else {
    stop("Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'.")
  }

  # Execute the command
  system(cmd, intern = FALSE)
}

# check_related_files_exist <- function(bed_file) {
#   base_name <- sub("\\.bed$", "", bed_file)
#   bim_file <- paste0(base_name, ".bim")
#   fam_file <- paste0(base_name, ".fam")
#
#   if (!file.exists(bim_file) || !file.exists(fam_file)) {
#     stop("Both .bim and .fam files must exist in the same directory as the .bed file and share the same base name.")
#   }
# }



# run_plink_filtering(
# input_file = "sik.vcf.gz",
# output_prefix = "clean",
# monomorphic = TRUE,
# snp_call_rate = 0.05,
# maf = 0.01,
# heterozygosity = 0.2,
# output_format = "vcf",
# plink_path = "plink"
# #plink_executable_name = "plink"
#
# )
