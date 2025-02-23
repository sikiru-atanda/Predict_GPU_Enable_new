process_var_u_new <- function(file = NULL,
                              posindex = NULL,
                              GS_model = NULL,
                              geno_data = NULL,
                              y = NULL,
                              B = NULL,
                              tst = NULL) {


  if(GS_model%in%c("RKHS")){
    #var_u <- scan(file, what = numeric(), sep = "\n", quiet = TRUE)[posindex]
    varU <- scan(file, what = numeric(), sep = "\n", quiet = TRUE)
    return(list(var_u_omics = varU))


  }else {
    if(GS_model%in%c("BayesA",  "BayesB", "BayesC", "BL", "BRR")){ # lambda
      h2 <- rep(NA,nrow(B))
      varU <- h2
      varE <- h2
      if(!is.null(tst)){
      y <- y[-tst]
      for(i in 1:length(h2)){
        u <- as.matrix(geno_data[-tst, ])%*%B[i,]
        varU[i] <- var(u)
        varE[i] <- var(y-u)
        h2[i] <- varU[i]/(varU[i]+varE[i])
      }

      } else{
        y <- y
        for(i in 1:length(h2)){
          u <- as.matrix(geno_data)%*%B[i,]
          varU[i] <- var(u)
          varE[i] <- var(y-u)
          h2[i] <- varU[i]/(varU[i]+varE[i])
        }
      }

    }

    return(list(var_u_omics = varU,
                var_residual = varE,
                h2 = h2))

  }

}
