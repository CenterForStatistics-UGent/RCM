#'A jacobian function for the estimation of the column scores in an unconstrained RC(M) model
#'
#'   @param beta vector of length p+1+1+(k-1): p row scores, 1 centering, one normalization  and (k-1) orhtogonality lagrangian multipliers
#' @param X the nxp data matrix
#' @param reg a nx1 regressor matrix: outer product of rowScores and psis
#' @param theta nxp matrix with the dispersion parameters (converted to matrix for numeric reasons)
#' @param muMarg the nxp offset
#' @param k an integer, the dimension of the RC solution
#' @param p an integer, the number of taxa
#' @param n an integer, the number of samples
#' @param nLambda an integer, the number of restrictions
#' @param colWeights the weights used for the restrictions
#' @param cMatK the lower dimensions of the colScores

#' @return A matrix of dimension p+1+1+(k-1) with evaluations of the Jacobian
NBjacobianCol = function(beta, X, reg, thetas, muMarg, k, n ,p, colWeights, nLambda, cMatK){
  cMat = matrix(beta[1:p], byrow=TRUE, nrow=1, ncol=p)

  #Calculate the mean
  mu = exp(reg %*% cMat) * muMarg

  Jac = matrix(0, nrow= p + nLambda, ncol=p + nLambda)
  #The symmetric jacobian matrix. The upper part is filled first, then mirror image is taken for lower triangle
  #
  #dLag²/dr_{ik}dlambda_{1k}
  Jac[1:p,(p+1)] = colWeights
  #Jac[1:(p*k),(p*k+1):((p+1)*k)] = sapply(1:k, function(K){c(rep(0,(K-1)*p),colWeights,rep(0,(k-K)*p))})
  Jac[1:p,p+2] = colWeights*2 *cMat

  tmp = (1+X/thetas)*mu/(1+mu/thetas)^2

  #dLag²/ds_{ik}dlambda_{3kk'}
  if(k>1){
    Jac[1:p,(p+3):(p+nLambda)] = t(cMatK)*colWeights
  }

  #Symmetrize
  Jac = Jac + t(Jac)
  #dLag²/dr_{ik}²
  diag(Jac)[1:p] = -crossprod(tmp, reg^2) + 2*beta[p+2]*colWeights
  Jac
}