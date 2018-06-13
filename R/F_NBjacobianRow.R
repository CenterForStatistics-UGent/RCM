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
#' @param preFabMat a prefab matrix, (1+X/thetas)
#' @param Jac an empty Jacobian matrix
#'
#' @return a symmetric jacobian matrix of size n+k + 1
NBjacobianRow = function(beta, X, reg, thetas, muMarg, k, n ,p, rowWeights, nLambda, rMatK, preFabMat, Jac){
  rMat = beta[seq_len(n)]
  mu = exp(rMat %*% reg)* muMarg

  #dLag²/dr_{ik}dlambda_{1k} already happened
  #dLag²/dr_{ik}dlambda_{2k}
  Jac[seq_len(n), n+2] = Jac[n+2, seq_len(n)] = 2 *rMat*rowWeights
  #dLag²/dr_{ik}dlambda_{3kk'} already happened

  #dLag²/dr_{ik}²
  diag(Jac)[seq_len(n)] = -tcrossprod(reg^2 ,preFabMat*mu/(1+mu/thetas)^2) + 2*rowWeights*beta[n+2]
  Jac
}
