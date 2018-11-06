#'A jacobian function for the estimation f the parameters of a third degree GLM in a constrained RC(M) model, if the GAM fit fails
#'
#' @param beta vector of any length
#' @param X the data vector of length n
#' @param reg a nxlength(beta) regressor matrix
#' @param theta a scalar, the overdispersion
#' @param muMarg the offset of length n
#'
#' @return A matrix of dimension 8-by-8
NBjacobianColNP = function(beta, X, reg, theta, muMarg){
  #Calculate the mean
  mu = exp(reg %*% beta) * muMarg
  #Return the Jacobian
  -crossprod(reg*c((1+X/theta)*mu/(1+mu/theta)^2), reg)
}
