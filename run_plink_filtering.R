run_plink_filtering <- function(input_file,
                                output_prefix,
                                snp_call_rate = NULL,
                                maf = NULL,
                                heterozygosity = NULL,
                                output_format = "bed",
                                allow_extra_chr = FALSE,
                                plink_executable_name = "plink") {

  # Determine the PLINK executable path based on the operating system
  if (.Platform$OS.type == "windows") {
    plink_path <- paste0(plink_executable_name, ".exe")
  } else {
    plink_path <- plink_executable_name
  }

  # Check for file extension
  file_ext <- tools::file_ext(input_file)
  if (file_ext == "gz") {
    file_ext <- tools::file_ext(sub("\\.gz$", "", input_file))
  }

  # Construct the base PLINK command according to the file extension
  cmd <- ""
  if (file_ext %in% c("bed", "bim", "fam")) {
    # Check if .bed file was provided, and infer base name for PLINK binary format
    bed_base <- gsub(pattern = "\\.[bim|fam|bed]+$", replacement = "", x = input_file)
    cmd <- sprintf('%s --bfile "%s"', plink_path, bed_base)
  } else if (file_ext == "vcf") {
    cmd <- sprintf('%s --vcf "%s" --allow-extra-chr', plink_path, input_file)
  } else {
    stop("Unsupported file format. Please provide a .bed (with .bim and .fam) or .vcf file.")
  }

  # Optionally allow extra chromosomes
  if (isTRUE(allow_extra_chr)) {
    cmd <- sprintf('%s --allow-extra-chr', cmd)
  }

  # Adjust for handling monomorphic SNPs
  if(isTRUE(monomorphic)){
    cmd <- sprintf('%s --min-alleles 2 --max-alleles 2', cmd)

  }

  # Append additional PLINK command arguments based on function parameters
  if (!is.null(snp_call_rate)) {
    cmd <- sprintf('%s --geno %f', cmd, snp_call_rate)
  }
  if (!is.null(maf)) {
    cmd <- sprintf('%s --maf %f', cmd, maf)
  }
  if (!is.null(heterozygosity)) {
    cmd <- sprintf('%s --het --mind %f', cmd, heterozygosity)
  }

  # Specify output format
  output_cmd <- switch(output_format,
                       "bed" = "--make-bed",
                       "vcf" = "--recode vcf",
                       "vcf.gz" = "--recode vcf bgz",
                       stop("Unsupported output format. Please specify 'bed', 'vcf', or 'vcf.gz'."))
  cmd <- sprintf('%s %s --out "%s"', cmd, output_cmd, output_prefix)

  # Execute the command
  message("Executing command: ", cmd)
  system(cmd, intern = FALSE)
}


# run_plink_filtering(
#   input_file = "sik.vcf.gz",
#   output_prefix = "clean",
#   #monomorphic = TRUE,
#   snp_call_rate = 0.05,
#   maf = 0.01,
#   heterozygosity = 0.2,
#   output_format = "vcf",
#   #plink_path = "plink"
#   plink_executable_name = "plink"
#
# )
