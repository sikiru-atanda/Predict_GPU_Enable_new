mme_marker_with_Z = function(y = NULL, X = NULL, Z = NULL, M = NULL, ve = NULL, va = NULL) {
  
  if (is.null(M)) M <- diag(1, nrow = length(y))
  if (is.null(Z)) Z <- diag(1, nrow = length(y))
  
  Cxx <- crossprod(X, X) / ve
  Cxz <- crossprod(X, Z) / ve
  Cxm <- crossprod(X, M) / ve
  Czz <- crossprod(Z, Z) / ve + solve(crossprod(Z, Z)) * va
  Cmm <- crossprod(M, M) / ve + solve(crossprod(M, M)) * va
  
  Xy <- crossprod(X, y) / ve
  Zy <- crossprod(Z, y) / ve
  My <- crossprod(M, y) / ve
  
  lhs <- rbind(cbind(Cxx, Cxz, Cxm), cbind(t(Cxz), Czz, rep(0, nrow(Czz))), cbind(t(Cxm), rep(0, nrow(Cxm)), Cmm))
  rhs <- rbind(Xy, Zy, My)
  
  C <- solve(lhs)
  sol <- crossprod(C, rhs)
  
  aii <- diag(crossprod(Z, Z))
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
# result_marker_with_Z <- mme_marker_with_Z(y = your_data$y, X = your_data$X, Z = your_data$Z, M = your_data$M, ve = your_data$ve, va = your_data$va)
