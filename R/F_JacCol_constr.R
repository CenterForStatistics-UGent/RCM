#' Returns the jacobian of this function
#'
#' @param betas a vector of v parameters of the response function of a single taxon (\boldsymbol{beta}_j)
#' @param X the count vector of length n
#' @param reg a n-by-v matrix of regressors(\mathbf{C} \boldsymbol{\alpha})
#' @param theta The dispersion parameter of this taxon
#' @param muMarg offset of length n
#' @param psi a scalar, the importance parameter
#'
#' Even though this approach does not imply normalization over the parameters of all taxa, it is very fast and they can be normalized afterwards
#'
#' @return The jacobian, a square symmetric matrix of dimension v
JacCol_constr = function(betas, X, reg, theta, muMarg, psi) {
  mu = exp(c(reg %*% betas)* psi) * muMarg
  tmp = (1+X/theta)*mu/(1+mu/theta)^2
  -crossprod(tmp*reg,reg) * psi^2 #Don't forget to square psi!

}