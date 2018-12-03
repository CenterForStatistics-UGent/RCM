#'Estimation of the parameters of a third degree GLM
#'
#' @param beta A vector of any length
#' @param X the data vector of length n
#' @param reg a nxlength(beta) regressor matrix
#' @param theta a scalar, the overdispersion
#' @param muMarg the offset of length n
#' @param ... further arguments passed on to the jacobian

#' @return A vector of the same length as beta with evaluations
#'  of the score function
dNBllcolNP = function(beta, X, reg, theta, muMarg, ...) {
  mu = exp(reg %*% beta) * muMarg
 crossprod(reg,((X-mu)/(1+mu/theta)))
}
