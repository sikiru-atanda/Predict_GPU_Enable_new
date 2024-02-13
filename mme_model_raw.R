#' Title
#'
#' @param y Vector of observed response variable (phenotype)
#' @param X Design matrix for fixed effects
#' @param Z Design matrix for random effects (genetic random effects)
#' @param GRM 
#' @param ve Residual variance
#' @param va Genetic variance
#'
#' @return
#' @export
#'
#' @examples
mme = function(y=NULL, X=NULL, Z=NULL,
               GRM=NULL, ve=NULL, va=NULL) {
  
  ## X and Z matrix should be generated as 
  # X <- model.matrix(BW ~ sex + reps, data=mouse)
  # Z <- model.matrix(BW ~ GID + Env, data=mouse)
  if(is.null(Z)) Z <- diag(1,nrow=length(y))
  Czz <- crossprod(Z,Z)/ve + solve(GRM*va) # Covariance matrix for random effects
  Zy <- crossprod(Z,y)/ve # Cross-product of Z and y
  
   if(!is.null(X)){
    Cxx <- crossprod(X,X)/ve # Covariance matrix for fixed effects
    Cxz <- crossprod(X,Z)/ve # Covariance matrix between fixed and random effects
    Xy <- crossprod(X,y)/ve # Cross-product of X and y
    # Left-hand side of the MME, formed by 
    # stacking matrices related to fixed and random effects
  lhs <- rbind(cbind(Cxx,Cxz),cbind(t(Cxz),Czz))
  # Right-hand side of the MME, formed by stacking vectors 
  # related to fixed and random effects
  rhs <- rbind(Xy,Zy)
  } else {
    
    lhs <- Czz
    rhs <- Zy
  }
  
  C <- solve(lhs)
  sol <- crossprod(C,rhs)
  aii <- diag(GRM)
  cii <- diag(C)
  if(!is.null(X)){
  fixed <- c(rep(TRUE,ncol(Cxx)),rep(FALSE,ncol(Czz)))
  random <- !fixed
  a <- sol[random] 
  pev <- cii[random]*ve
  sep <- sqrt(pev)
  rel <- (aii-cii[random]*ve/va)*aii
  b <- sol[fixed]
  seb <- cii[fixed]
  
  return(list(b=b,seb=seb,a=a,pev=pev,sep=sep,rel=rel))
  
  } else {
    
    a <- sol
    pev <- cii * ve
    sep <- sqrt(pev)
    rel <- (aii - cii * ve/va) * aii
    
    return(list(a=a, pev=pev, sep=sep, rel=rel))
  }
 
}
