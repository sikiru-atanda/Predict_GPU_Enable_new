
load("barley_data.Rdata")

if (!requireNamespace("kernlab", quietly = TRUE)) {
  install.packages("kernlab")
  library("kernlab")
} else {
  library("kernlab")
}

toy_pred_for_souradip <- function(geno_data,
                                  pheno,
                                  response
                                  ){


  fit <- ksvm(x = as.matrix(geno_data),
              y = pheno[, response])

  return(predict(fit,
                  geno_data))


}

res <- toy_pred_for_souradip(geno_data = geno_data,
                             pheno = pheno_data1,
                             response = "B_GLUCAN")
