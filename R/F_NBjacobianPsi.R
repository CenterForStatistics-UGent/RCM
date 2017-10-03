#' A jacobian function for the psi of a given dimension (importance parameters)
#'
#' @param beta a scalar, the current estimate
#' @param X the n-by-p count matrix
#' @param muMarg the nxp offset matrix
#' @param reg the regressor matrix, the outer product of current row and column scores
#' @param theta a n-by-p matrix with the dispersion parameters
#' @param preFabMat a prefab matrix, (1+X/thetas)
#'
#' @return The evaluation of the jacobian function at beta, a 1-by-1 matrix
NBjacobianPsi = function(beta, X, reg, muMarg, theta, preFabMat){
  # @param beta: the current guess for the psi parameter, a scalar
  # @param X: the nxp data matrix
  # @param muMarg: the nxp offset
  # @param reg: the regressor matrix, the outer product of current row and column scores
  # @param theta: nxp matrix with the dispersion parameters (converted to matrix for numeric reasons)

  # @return The evaluation of the Jacobian function at beta, a 1x1 matrix

  mu = muMarg * exp(reg* beta)
  matrix(-sum(reg^2*preFabMat*mu/(1+mu/theta)^2),1,1)

}