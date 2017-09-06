#' Calculates the Jacobian of the score functions of the paramaters of parametric response functions  (a polynomial of any degree). Parameters are sorted per taxon, than with increasing degree.
#'
#' @param betas: a vector of length (deg+1)*(p+1) with regression parameters with deg the degree of the response function and the lagrangian multipliers
#' @param aX the nxp data matrix
#' @param reg a vector of regressors with the dimension n-by-v (\mathbf{C} \boldsymbol{\alpha})
#' @param thetaMat The n-by-p matrix with dispersion parameters
#' @param mumarg offset matrix of size nxp
#' @param psi a scalar, the importance parameter
#' @param v an integer, one plus the degree of the response function
#' @param p an integer, the number of taxa
#'
#' @return The jacobian, a square matrix of dimension (deg+1)*(p+1)

respFunJacMat = function(betas, aX, reg, thetaMat, muMarg, lambda2, psi, v, p) {
  NBparams = matrix(betas[seq_len(p*v)], ncol = p)
  mu = exp(reg %*% NBparams*psi) * muMarg
  Jac = matrix(0, (p+1)*v, (p+1)*v)
  did  = seq_len(p*v)
  didv = seq_len(v)
  indVec = matrix(rep(c(TRUE, rep(FALSE, v), TRUE),p), nrow = p*v, byrow = TRUE)
  Jac[did, didv][indVec] = NBparams
  Jac[lower.tri(Jac)] = Jac[upper.tri(Jac)]

  tmp = (1+X/thetaMat)*mu/(1+(mu/thetaMat)^2)*psi^2
  tmp2 =  vapply(didv, FUN.VALUE = tmp, function(x){reg[,x]*tmp})
  Jac[did, did] = tensor(reg, tmp2, 1, 1)

  diag(Jac)[did] =  diag(Jac)[did] + 2*betas[seq_len(v) + p*v]
  return(-Jac)
}