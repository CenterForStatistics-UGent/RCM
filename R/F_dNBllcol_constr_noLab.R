#' Evaluate the score function of the response function ignoring taxon labels
#'
#' @param beta a vector of regression parameters with length v (\boldsymbol{beta}_j)
#' @param X the nxp data matrix
#' @param reg a matrix of regressors of dimension nxv (\mathbf{C} \boldsymbol{\alpha})
#' @param thetas The dispersion parameters (a vector of length p)
#' @param muMarg offset matrix of dimension nxp
#' @param psi a scalar, the importance parameter
#'
#' @return The evaluation of the score functions (a vector length v)
dNBllcol_constr_noLab = function(betas, X, reg, thetas, muMarg, psi, ...) {
  mu = c(exp(reg %*% betas*psi)) * muMarg
  colSums(crossprod((X-mu)/(1+(mu/thetas)),reg) * psi)
}