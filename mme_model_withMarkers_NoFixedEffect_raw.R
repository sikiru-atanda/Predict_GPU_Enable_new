mme_random_effects_only = function(y = NULL, Z = NULL, M = NULL, ve = NULL, va = NULL) {
  
  if (is.null(M)) M <- diag(1, nrow = length(y))
  if (is.null(Z)) Z <- diag(1, nrow = length(y))
  
  Czz <- crossprod(Z, Z) / ve + solve(crossprod(Z, Z)) * va
  Cmm <- crossprod(M, M) / ve + solve(crossprod(M, M)) * va
  
  Zy <- crossprod(Z, y) / ve
  My <- crossprod(M, y) / ve
  
  lhs <- rbind(cbind(Czz, rep(0, nrow(Czz))), cbind(rep(0, nrow(Cmm)), Cmm))
  rhs <- rbind(Zy, My)
  
  C <- solve(lhs)
  sol <- crossprod(C, rhs)
  
  aii <- diag(crossprod(Z, Z))
  cii <- diag(C)
  
  random <- !rep(TRUE, ncol(Czz))
  
  a <- sol[random]
  
  pev <- cii[random] * ve
  sep <- sqrt(pev)
  rel <- (aii - cii[random] * ve / va) * aii
  
  return(list(a = a, pev = pev, sep = sep, rel = rel))
}

# Example usage:
# Replace 'your_data' with actual data
# result_random_effects_only <- mme_random_effects_only(y = your_data$y, Z = your_data$Z, M = your_data$M, ve = your_data$ve, va = your_data$va)
