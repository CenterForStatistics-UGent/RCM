#' A score function of the NB for the row scores
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

#' @return A vector of length n + k +1 with evaluations of the derivative of the lagrangian
dNBllrow = function(beta, X, reg, thetas, muMarg, k, n , p, rowWeights, nLambda, rMatK) {

  rMat = matrix(beta[1:n], byrow=FALSE, ncol=1, nrow=n)
  mu = exp(rMat %*% reg)* muMarg

  lambda1 = beta[n+1] #Centering restrictions sum(abunds*r_{ik}) = 0
  lambda2 = beta[n+2] #normalization restrictions sum(abunds*r^2_{ik}) = 1
  lambda3 = if(k==1){0} else {beta[(n+3):length(beta)]}

  score = c(tcrossprod(reg, (X-mu)/(1+(mu/thetas)))) + rowWeights*(lambda1 + lambda2* 2*rMat + (rMatK %*% lambda3))
  center = sum(rMat*rowWeights)
  unitSum = sum(rMat^2*rowWeights)-1
  if(k==1){ return(c(score,center, unitSum))}
  orthogons = crossprod(rMatK, rMat*rowWeights)
  return(c(score,center, unitSum, orthogons))
}