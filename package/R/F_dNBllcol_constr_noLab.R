#' Evaluate the score function of the response function ignoring taxon labels
#'
#' @param betas a vector of regression parameters with length v
#' @param X the nxp data matrix
#' @param reg a matrix of regressors of dimension nxv
#' @param thetasMat A matrix of dispersion parameters
#' @param muMarg offset matrix of dimension nxp
#' @param psi a scalar, the importance parameter
#' @param ... further arguments passed on to the jacobian
#'
#' @return The evaluation of the score functions (a vector length v)
dNBllcol_constr_noLab = function(betas, X, reg, thetasMat, muMarg, psi, ...) {
  mu = c(exp(reg %*% betas*psi)) * muMarg
  colSums(crossprod((X-mu)/(1+(mu/thetasMat)),reg) * psi)
}
