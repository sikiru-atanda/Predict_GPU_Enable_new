# Define a custom function to convert SNP data from CSV to VCF
convert_csv_to_vcf <- function(input_csv, output_vcf = NULL) {

  # Required VCF fields, ensuring all are present in the data
  required_headers <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

  # Read the CSV file into a data frame
  snp_data <- tryCatch({
     data.table::fread(input_csv, header = TRUE, stringsAsFactors = FALSE,
                      skip = "#", check.names = FALSE, na.strings = c(NA, "N", "."))
  }, error = function(e) {
    # If the skip causes an error, read without skipping
   data.table::fread(input_csv, header = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE, na.strings = c(NA, "N", "."))
  })

  snp_data <-  as.data.frame(snp_data)
  data.table::setnames(snp_data, old = names(snp_data), new = gsub("^#", "", names(snp_data)))
  # Check if all required headers are present
  if (!all(required_headers %in% colnames(snp_data))) {
    stop(paste("Missing required VCF headers:",
               paste(setdiff(required_headers, colnames(snp_data)), collapse = ", "),
               ". Please provide the necessary columns."))
  }

  # Ensure the data follows the correct VCF format with 9 required fields first, followed by sample data
  snp_data <- snp_data[, c(required_headers, setdiff(colnames(snp_data), required_headers))]

  # Prepare the VCF data content as a character vector
  vcf_content <- c()

  # Add VCF file headers
  vcf_content <- c(vcf_content, "##fileformat=VCFv4.2")
  vcf_content <- c(vcf_content, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")

  # Write column headers (with 9 required fields plus sample names)
  sample_names <- colnames(snp_data)[-(1:9)]  # Sample names start after the 9th column
  vcf_header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names)
  vcf_content <- c(vcf_content, paste(vcf_header, collapse = "\t"))

  # Helper function to check if SNP is already in VCF format (e.g., "0/0", "0/1", "1/1")
  is_vcf_format <- function(snp_value) {
    return(grepl("^([0-1]/[0-1])$", snp_value))
  }

  # Helper function to recode SNPs (handles both 0, 1, 2 and -1, 0, 1 coding)
  recode_snp <- function(snp_value) {
    if (is_vcf_format(snp_value)) {
      return(snp_value)  # Already in VCF format, no need to recode
    } else if (is.na(snp_value) || snp_value == ".") {
      return("./.")  # Missing genotype
    } else if (snp_value == 0 || snp_value == -1) {
      return("0/0")  # Homozygous reference
    } else if (snp_value == 1 || snp_value == 0) {
      return("0/1")  # Heterozygous
    } else if (snp_value == 2 || snp_value == 1) {
      return("1/1")  # Homozygous alternate
    } else {
      stop("Unexpected SNP coding value")
    }
  }

  # Process each row of SNP data and format it for VCF
  for (i in 1:nrow(snp_data)) {
    chrom <- snp_data[i, "CHROM"]
    pos <- snp_data[i, "POS"]
    id <- ifelse(snp_data[i, "ID"] == "", ".", snp_data[i, "ID"])  # Handle missing IDs
    ref <- snp_data[i, "REF"]
    alt <- snp_data[i, "ALT"]
    qual <- snp_data[i, "QUAL"]
    filter <- snp_data[i, "FILTER"]
    info <- snp_data[i, "INFO"]
    format <- snp_data[i, "FORMAT"]

    # Apply the recode_snp function to genotype data for each sample
    genotypes <- sapply(snp_data[i, -(1:9)], recode_snp)

    # Combine all columns into a single VCF row
    vcf_row <- paste(chrom, pos, id, ref, alt, qual, filter, info, format, paste(genotypes, collapse = "\t"), sep = "\t")

    # Append the row to the VCF content
    vcf_content <- c(vcf_content, vcf_row)
  }

  # If output_vcf is provided, save the VCF content to a file
  if (!is.null(output_vcf)) {
    writeLines(vcf_content, output_vcf)
    cat("VCF file has been successfully created at:", output_vcf, "\n")
  } else {
    # Otherwise, return the VCF content as a character vector
    return(vcf_content)
  }
}


# toy_snp_data <- data.frame(
#   CHROM = c(1, 1, 2),
#   POS = c(1001, 1002, 2001),
#   ID = c("rs1", "rs2", "rs3"),
#   REF = c("A", "G", "T"),
#   ALT = c("G", "A", "C"),
#   QUAL = c(100, 200, 300),
#   INFO = c(".", ".", "."),
#   FORMAT = c("GT", "GT", "GT"),
#   Sample1 = c("0/0", "0/1", "1/1"),
#   Sample2 = c("0/1", "1/1", "0/0"),
#   stringsAsFactors = FALSE
# )
# toy_snp_data_recode <- data.frame(
#   CHROM = c(1, 1, 2),
#   POS = c(1001, 1002, 2001),
#   ID = c("rs1", "rs2", "rs3"),
#   REF = c("A", "G", "T"),
#   ALT = c("G", "A", "C"),
#   QUAL = c(100, 200, 300),
#   FILTER = c("PASS", "PASS", "PASS"),
#   INFO = c(".", ".", "."),
#   FORMAT = c("GT", "GT", "GT"),
#   Sample1 = c("0", "1", "2"),
#   Sample2 = c("1", "0", "."),
#   stringsAsFactors = FALSE
# )
#
# # # Write this toy data frame to a CSV file for testing
#  write.csv(toy_snp_data, "toy_snp_data.csv", row.names = FALSE)
# #
# snp_data <-  read.csv("soy21_22_AllGenotypeData.csv", header = T, sep = ",", as.is = T,
#                       stringsAsFactors = FALSE, row.names = 1)
#
# snp_data <- t(snp_data)
# snp_data <- cbind(ID = rownames(snp_data), snp_data)
#
# write.csv(snp_data, "snp_data.csv", row.names = FALSE)

# convert_csv_to_vcf(input_csv = "toy_snp_data.csv", "toy_snp_data4.vcf")
#
#
#
#
