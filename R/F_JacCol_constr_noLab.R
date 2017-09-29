#' The jacobian of the estimating function for the response function without taxon labels
#'
#' @param beta a vector of regression parameters with length v (\boldsymbol{beta}_j)
#' @param X the nxp data matrix
#' @param reg a matrix of regressors of dimension nxv (\mathbf{C} \boldsymbol{\alpha})
#' @param thetas The dispersion parameters (a vector of length p)
#' @param muMarg offset matrix of dimension nxp
#' @param preFabMat a prefabricated matrix
#' @param psi a scalar, the importance parameter
#' @param n an integer, number of rows of X
#' @param v an integer, the number of parameters of the response function
#'
#' @return The jacobian (a v-by-v matrix)
JacCol_constr_noLab = function(betas, X, reg, thetas, muMarg, psi, n ,v, preFabMat) {
  mu = c(exp(reg %*% betas*psi)) * muMarg
  tmp = preFabMat*mu/(1+(mu/thetas))^2*psi^2 #Don't forget to square psi!
  -crossprod(reg, vapply(seq_len(v), FUN.VALUE = vector("numeric", n), function(x){rowSums(reg[,x]*tmp)}))
}