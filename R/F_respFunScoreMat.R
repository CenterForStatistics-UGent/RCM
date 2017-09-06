#' Returns the derivative of the Lagrangian of the parameters of the parametric response function (a polynomial of any degree)
#'
#' @param betas: a (deg+1)xp matrix of regression parameters with deg the degree of the response function
#' @param aX: the nxp data matrix
#' @param reg: a matrix of regressors with the dimension (deg+1)xn (offset + \mathbf{C} \boldsymbol{\alpha})
#' @param thetaMat: The n-by-p matrix with dispersion parameters
#' @param mumarg: offset matrix of size nxp
#' @param lambda2: A vector of length (deg+1) with the Lagrange multipliers
#'
#' The parameters are restricted to be normalized, i.e. all squared intercepts, first order and second order parameters sum to 1
#'
#' @return The evaluation of the score functions, a vector of length (p+1)*(deg+1)
#'
respFunScoreMat = function(betas, aX, reg, thetaMat, muMarg, psi, lambda2) {
  mu = exp(crossprod(reg, betas)*psi) * muMarg
  score =  reg %*% ((aX-mu)/(1+mu/thetaMat)) *2*lambda2 *betas
  norm = rowSums(betas^2) - 1
  return(c(score, norm)) #Taxon per taxon
}