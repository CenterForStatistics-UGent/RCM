#' A jacobian function of the NB for the row scores
#'
#' @param beta a vector of of length n + k +1 regression parameters to optimize
#' @param X the data matrix of dimensions nxp
#' @param reg a 1xp regressor matrix: outer product of column scores and psis
#' @param thetas nxp matrix with the dispersion parameters (converted to matrix for numeric reasons)
#' @param muMarg an nxp offset matrix
#' @param k a scalar, the dimension of the RC solution
#' @param p a scalar, the number of taxa
#' @param n a scalar, the number of samples
#' @param rowWeights a vector of length n, the weights used for the restrictions
#' @param nLambda an integer, the number of lagrangian multipliers
#' @param rMatK the lower dimension row scores

#' @return a symmetric jacobian matrix of size n+k + 1
NBjacobianRow = function(beta, X, reg, thetas, muMarg, k, n ,p, nlambda, rowWeights, nLambda, rMatK){
  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg)* muMarg

  Jac = matrix(0, nrow= n + nLambda, ncol= n + nLambda)
  #The symmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle

  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:n, n+1] = rowWeights
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[1:n, n+2] = 2 *rMat*rowWeights

  #dLag²/dr_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:n,(n+3):(n+nLambda)] = rMatK*rowWeights
  }
  #Symmetrize
  Jac = Jac + t(Jac)
  #dLag²/dr_{ik}²
  diag(Jac)[1:n] = -tcrossprod(reg^2 ,(1+X/thetas)*mu/(1+mu/thetas)^2) + 2*rowWeights*beta[n+2]
  Jac
}