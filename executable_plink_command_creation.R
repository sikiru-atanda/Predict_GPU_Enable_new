#' Title
#'
#' @param vcf_file
#' @param out_prefix
#' @param allow_extra_chr
#' @param min_alleles
#' @param max_alleles
#' @param geno
#' @param maf
#' @param export_format
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
create_plink_command <- function(vcf_file,
                                 out_prefix,
                                 allow_extra_chr = TRUE,
                                 min_alleles = 2, max_alleles = 2,
                                 geno = 0.15, maf = 0.05,
                                 export_format = 'vcf bgz',
                                 output_file = "plink_command.sh") {
  # Construct the command with dynamic parameters
  plink_command <- sprintf(
    "plink2 --vcf %s %s --min-alleles %d --max-alleles %d --geno %.2f --maf %.2f --export %s --out %s",
    vcf_file,
    ifelse(allow_extra_chr, "--allow-extra-chr \\", ""),
    min_alleles,
    max_alleles,
    geno,
    maf,
    export_format,
    out_prefix
  )

  # Ensure the command is formatted with line breaks for readability
  plink_command_formatted <- gsub(" \\\\", " \\\\\n       ", plink_command)

  # Write the command to the specified output file
  writeLines(plink_command_formatted, output_file)

  # Optionally, make the script executable on Unix-like systems
  system(paste("chmod +x", output_file))

  cat("PLINK command script created at:", output_file, "\n")
}

# Example usage of the function
create_plink_command(
  vcf_file = "FFAR_Chinese_Pea_ref_genome_filtered.vcf.gz",
  out_prefix = "initial_qc_data",
  geno = 0.1, # Example of changing the geno parameter
  maf = 0.05,
  output_file = "plink_command.sh"
)
