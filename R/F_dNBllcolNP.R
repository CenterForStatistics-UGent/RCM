#'A score function for the estimation of the parameters of a third degree GLM in a constrained RC(M) model, if the GAM fit fails
#'
#' @param beta A vector of any lenght
#' @param X the data vector of length n
#' @param reg a nxlength(beta) regressor matrix
#' @param theta a scalar, the overdispersion
#' @param muMarg the offset of length n
#' @param ... further arguments passed on to the jacobian

#' @return A vector of length 8 with evaluations of the score function
dNBllcolNP = function(beta, X, reg, theta, muMarg, ...) {
  mu = exp(reg %*% beta) * muMarg
 crossprod(reg,((X-mu)/(1+mu/theta)))
}
