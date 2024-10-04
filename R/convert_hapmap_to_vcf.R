

###### Function to determine REF and ALT alleles based on allele frequencies in the genotypes
determine_ref_alt <- function(alleles_str, genotypes) {
  #### Extract allele1 and allele2 from the 'alleles' string
  allele1 <- substr(alleles_str, 1, 1)
  allele2 <- substr(alleles_str, 3, 3)

  ###### Count occurrences of each allele across all sample genotypes
  allele_counts <- c(
    allele1 = sum(sapply(genotypes, function(g) grepl(allele1, g))),
    allele2 = sum(sapply(genotypes, function(g) grepl(allele2, g)))
  )

  ##### Determine the major and minor alleles
  if (allele_counts["allele1"] >= allele_counts["allele2"]) {
    ref <- allele1  # Major allele is the reference (REF)
    alt <- allele2  # Minor allele is the alternate (ALT)
  } else {
    ref <- allele2
    alt <- allele1
  }

  return(list(ref = ref, alt = alt))
}

############## Function to convert double alleles to IUPAC codes
double_code_to_IUPAC <- function(double_code) {
  iupac_map <- list(
    "AA" = "A", "TT" = "T", "CC" = "C", "GG" = "G",
    "AT" = "W", "TA" = "W", "CG" = "S", "GC" = "S",
    "AC" = "M", "CA" = "M", "GT" = "K", "TG" = "K",
    "AG" = "R", "GA" = "R", "TC" = "Y", "CT" = "Y",
    "NN" = "N"
  )

  ##### Remove potential separators like slashes or colons and convert to uppercase
  double_code <- gsub("[:/_]", "", toupper(double_code))

  return(iupac_map[[double_code]])
}

############ Function to standardize SNP format to IUPAC
standardize_snp_format <- function(snp) {
  if (is.na(snp)) {
    return(NA)
  }

  #### Remove potential separator characters and convert to uppercase
  snp <- gsub("[:/_]", "", toupper(snp))

  #### If the SNP is a double code, convert it to IUPAC code
  if (nchar(snp) == 2) {
    snp <- double_code_to_IUPAC(snp)
  }
  return(snp)
}

############ Function to ensure all SNPs in HapMap are IUPAC-compatible
IUPAC_hapmap_compatible <- function(hapmap) {

  #### Identify the SNP columns (starting from the 12th column)
  snp_cols <- names(hapmap)[12:ncol(hapmap)]

  ##### Standardize SNP format for each SNP column
  hapmap <- as.data.frame(hapmap)
  snp_data <- apply(hapmap[, colnames(hapmap) %in% snp_cols], 2, function(col) sapply(col, standardize_snp_format))
  rownames(snp_data) <- NULL
  hapmap <- cbind(hapmap[, 1:11], snp_data)
  hapmap <- data.table::as.data.table(hapmap)
  gc()

  return(hapmap)
}


###### Main function to convert HapMap to VCF format,
convert_hapmap_to_vcf <- function(input_hapmap, output_vcf = NULL) {
  #### Define heterozygous and missing values
  heterozygous <- c('R', 'Y', 'S', 'W', 'K', 'M')
  missing_values <- c(NA, "NA", "N", "NN", "B", "V", "H", "D", ".", "-")

  ### impport HapMap data
  hapmap_data <- data.table::fread(input_hapmap, sep='\t', header=TRUE,
                                   check.names = FALSE,
                                   na.strings=c(NA, "N", "NN", "B", "V", "H", "D", ".", "-"))

  ### Prepare the VCF content as a character vector
  vcf_content <- c()
  vcf_content <- c(vcf_content, "##fileformat=VCFv4.2")
  vcf_content <- c(vcf_content, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")

  sample_names <- colnames(hapmap_data)[12:ncol(hapmap_data)]
  vcf_header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names)
  vcf_content <- c(vcf_content, paste(vcf_header, collapse = "\t"))

  hapmap_data <- IUPAC_hapmap_compatible(hapmap_data)
  hapmap_data <- as.data.frame(hapmap_data)

  ### Processing each row of HapMap data
  vcf_rows <- apply(hapmap_data, 1, function(row) {
    chrom <- row["chrom"]
    pos <- row["pos"]
    id <- row["rs#"]

    ### Extract alleles from the "alleles" column
    alleles_str <- row["alleles"]

    ### Extract genotypes across samples for this SNP
    genotypes <- row[12:length(row)]

    ### Determine REF and ALT alleles based on allele frequencies
    ref_alt <- determine_ref_alt(alleles_str, genotypes)
    ref <- ref_alt$ref
    alt <- ref_alt$alt

    ### Replace heterozygous and missing markers
    genotypes[genotypes %in% heterozygous] <- "0/1"
    genotypes[genotypes %in% missing_values] <- "./."

    #### Convert to vcf format based on allele1 (REF) and allele2 (ALT)
    genotypes <- sapply(genotypes, function(genotype) {
      if (genotype == "0/1") {
        return("0/1")  # Leave heterozygous as is
      } else if (grepl(ref, genotype)) {
        return("0/0")  # Homozygous for the major allele (REF)
      } else if (grepl(alt, genotype)) {
        return("1/1")  # Homozygous for the minor allele (ALT)
      } else {
        return("./.")  # Anything else is set to missing
      }
    })

    ### Prepare the VCF format
    qual <- "."
    filter <- "PASS"
    info <- "."
    format <- "GT"

    ### Create the VCF row
    vcf_row <- paste(chrom, pos, id, ref, alt, qual, filter, info, format, paste(genotypes, collapse = "\t"), sep = "\t")

    return(vcf_row)
  })

  ### Append all the rows to the VCF content
  vcf_content <- c(vcf_content, vcf_rows)

 ##
  if (!is.null(output_vcf)) {
    writeLines(vcf_content, output_vcf)
    cat("VCF file has been successfully created at:", output_vcf, "\n")
  } else {
    return(vcf_content)
  }
}



# # Toy HapMap data frame
# set.seed(123)
#
# # Generate hapmap data
# n_markers <- 1000
# n_samples <- 200
#
#
# rs_ids <- paste0("rs", 1:n_markers)
#
#
# alleles <- paste0(
#   sample(c("A", "C", "T", "G"), n_markers, replace = TRUE),
#   "/",
#   sample(c("A", "C", "T", "G"), n_markers, replace = TRUE)
# )
#
#
# chrom <- sample(1:5, n_markers, replace = TRUE)
# pos <- sample(1:10000, n_markers)
#
#
# strand <- rep("+", n_markers)
# assembly <- center <- protLSID <- assayLSID <- panelLSID <- QCcode <- rep(".", n_markers)
#
#
# generate_genotypes <- function() {
#   sample(c("A/A", "A/G", "G/G", "T/T", "T/C", "C/C", "C/G"), n_markers, replace = TRUE)
# }
#
#
# genotypes <- replicate(n_samples, generate_genotypes())
#
#
# sample_names <- paste0("Sample", 1:n_samples)
#
#
# toy_hapmap_data <- data.frame(
#   `rs#` = rs_ids,
#   alleles = alleles,
#   chrom = chrom,
#   pos = pos,
#   strand = strand,
#   assembly = assembly,
#   center = center,
#   protLSID = protLSID,
#   assayLSID = assayLSID,
#   panelLSID = panelLSID,
#   QCcode = QCcode,
#   stringsAsFactors = FALSE
# )
#
#
# toy_hapmap_data <- cbind(toy_hapmap_data, as.data.frame(genotypes))
# colnames(toy_hapmap_data)[12:(11 + n_samples)] <- sample_names
#
# colnames(toy_hapmap_data)[1] <- "rs#"
# # Display the first few rows of the data
# head(toy_hapmap_data)
#
# library(data.table)
# # Write this toy HapMap data frame to a file for testing
# write.table(toy_hapmap_data, "toy_hapmap_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
# # Example usage:
# convert_hapmap_to_vcf("toy_hapmap_data.txt", "output_vcf_file4.vcf")
#
#

#convert_hapmap_to_vcf("HapMap_300_Aug25_Select_Chr_No_Filtering.hmp.txt", "output_vcf_file.vcf")

