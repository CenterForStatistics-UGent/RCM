#' A function to evaluate the score function for the parameters of the response function for 1 taxon at the time
#' @param betas a vector of v parameters of the response function of a single taxon (\boldsymbol{beta}_j)
#' @param X the count vector of length n
#' @param reg a n-by-v matrix of regressors(\mathbf{C} \boldsymbol{\alpha})
#' @param theta The dispersion parameter of this taxon
#' @param muMarg offset of length n
#' @param psi a scalar, the importance parameter
#'
#'  Even though this approach does not imply normalization over the parameters of all taxa, it is very fast and they can be normalized afterwards
#'
#' @return A vector of length v with the evaluation of the score functions
dNBllcol_constr = function(betas, X, reg, theta, muMarg, psi) {
  mu = exp(c(reg %*% betas)*psi) * muMarg
  crossprod((X-mu)/(1+mu/theta) , reg) * psi
}