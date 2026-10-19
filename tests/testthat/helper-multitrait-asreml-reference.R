# Dense references for the multi-trait ASReml reliability and FA tests.

# G0 / R0 from ASReml's vparameters names: "...!Trait_T2:T1".
mt_rel_matrix <- function(vp, pattern, traits) {
  out <- matrix(0, length(traits), length(traits))
  for (x in names(vp)[grepl(pattern, names(vp), fixed = TRUE)]) {
    last <- utils::tail(strsplit(x, "!", fixed = TRUE)[[1]], 1)
    ab <- match(strsplit(sub("Trait_", "", last, fixed = TRUE), ":", fixed = TRUE)[[1]], traits)
    if (length(ab) == 2L && !anyNA(ab)) out[ab[1], ab[2]] <- out[ab[2], ab[1]] <- vp[[x]]
  }
  out
}

# Var(g - ghat) = G - G Z' P Z G for records stacked by trait.
mt_rel_reference_pev <- function(ph, traits, K, G0, R0) {
  n <- nrow(K)
  yv <- unlist(ph[traits], use.names = FALSE); obs <- !is.na(yv)
  X <- kronecker(diag(length(traits)), matrix(1, n, 1))[obs, , drop = FALSE]
  Zm <- diag(length(traits) * n)[obs, , drop = FALSE]
  G <- kronecker(G0, K)
  V <- Zm %*% G %*% t(Zm) + kronecker(R0, diag(n))[obs, obs]
  Vi <- solve(V)
  P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X, t(X) %*% Vi)
  diag(G - G %*% t(Zm) %*% P %*% Zm %*% G)
}
