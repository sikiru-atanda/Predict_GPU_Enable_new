mme_marker = function(y = NULL, X = NULL, M = NULL, ve = NULL, va = NULL) {
  
  if (is.null(M)) M <- diag(1, nrow = length(y))
  
  Cxx <- crossprod(X, X) / ve
  Cxm <- crossprod(X, M) / ve
  Cmm <- crossprod(M, M) / ve + solve(crossprod(M, M)) * va
  
  Xy <- crossprod(X, y) / ve
  My <- crossprod(M, y) / ve
  
  lhs <- rbind(cbind(Cxx, Cxm), cbind(t(Cxm), Cmm))
  rhs <- rbind(Xy, My)
  
  C <- solve(lhs)
  sol <- crossprod(C, rhs)
  
  aii <- diag(crossprod(M, M))
  cii <- diag(C)
  
  fixed <- rep(TRUE, ncol(Cxx))
  random <- !fixed
  
  b <- sol[fixed]
  a <- sol[random]
  
  pev <- cii[random] * ve
  sep <- sqrt(pev)
  rel <- (aii - cii[random] * ve / va) * aii
  
  return(list(b = b, a = a, pev = pev, sep = sep, rel = rel))
}

# Example usage:
# Replace 'your_data' with actual data
# result_marker <- mme_marker(y = your_data$y, X = your_data$X, M = your_data$M, ve = your_data$ve, va = your_data$va)
