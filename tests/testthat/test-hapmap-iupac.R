test_that("IUPAC HapMap conversion preserves one-row genotype shape", {
  hapmap <- data.frame(
    rs = "m1", alleles = "A/T", chrom = 1, pos = 10, strand = "+",
    assembly = NA, center = NA, protLSID = NA, assayLSID = NA,
    panelLSID = NA, QCcode = NA, X1 = "AA", X2 = "TA",
    stringsAsFactors = FALSE
  )

  out <- IUPAC_hapmap_compatible(hapmap)

  expect_s3_class(out, "data.table")
  expect_identical(dim(out), c(1L, 13L))
  expect_identical(names(out), names(hapmap))
  expect_identical(as.character(out$X1), "A")
  expect_identical(as.character(out$X2), "W")
})

test_that("IUPAC HapMap conversion rejects missing metadata columns", {
  expect_error(
    IUPAC_hapmap_compatible(data.frame(rs = "m1", X1 = "AA")),
    "11 metadata columns"
  )
})

test_that("genomic coding detection accepts vector inputs", {
  expect_identical(
    detect_genomic_coding(c(0, 1, 2, 0, 1, 2)),
    "SNP (0, 1, 2)"
  )
  expect_identical(
    detect_genomic_coding(c("0", "1", "0", "1")),
    "Presence/Absence (0, 1)"
  )
})
