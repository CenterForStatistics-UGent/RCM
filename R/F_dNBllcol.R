#'A score function for the estimation of the column scores
#' in an unconstrained RC(M) model
#'
#' @param beta vector of length p+1+1+(k-1): p row scores,
#' 1 centering, one normalization and (k-1) orhtogonality lagrangian multipliers
#' @param X the nxp data matrix
#' @param reg a nx1 regressor matrix: outer product of rowScores and psis
#' @param thetas nxp matrix with the dispersion parameters
#' (converted to matrix for numeric reasons)
#' @param muMarg the nxp offset
#' @param k an integer, the dimension of the RC solution
#' @param p an integer, the number of taxa
#' @param n an integer, the number of samples
#' @param nLambda an integer, the number of restrictions
#' @param colWeights the weights used for the restrictions
#' @param cMatK the lower dimensions of the colScores
#' @param ... further arguments passed on to the jacobian

#' @return A vector of length p+1+1+(k-1) with evaluations of the
#'  derivative of lagrangian
dNBllcol = function(beta, X, reg, thetas, muMarg, k, p, n, colWeights, nLambda,
                    cMatK, ...) {
  cMat = matrix(beta[seq_len(p)], byrow=TRUE, ncol=p, nrow=1)
  mu = exp(reg %*% cMat) * muMarg

  lambda1 = beta[p+1]
  #Lagrangian multiplier for centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[p+2]
  #Lagrangian multiplier for normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(p + 3):length(beta)]}
  #Lagrangian multiplier for orthogonalization restriction

  score = crossprod(reg,((X-mu)/(1+mu/thetas))) +
    colWeights*(lambda1 + lambda2*2*cMat + (lambda3 %*% cMatK))

  center = sum(colWeights*cMat)
  unitSum = sum(colWeights*cMat^2)-1
  if(k==1) {
    return(c(score, center, unitSum))
  }
  orthogons = tcrossprod(cMatK, cMat*colWeights)
  return(c(score, center, unitSum, orthogons))
}
