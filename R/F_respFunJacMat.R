#' Calculates the Jacobian of the score functions of the paramaters of parametric response functions  (a polynomial of any degree)
#'
#' @param beta: a (deg+1)*p matrix of regression parameters with deg the degree of the response function
#' @param X: the nxp data matrix
#' @param reg: a vector of regressors with the dimension (deg+1)*n (\mathbf{C} \boldsymbol{\alpha})
#' @param theta: The dispersion parameter
#' @param mumarg: offset matrix of size nxp
#' @param cReg: a vector of lentgh p psi_k * s_{jk}
#'
#' @return The jacobian, a square matrix of dimension (deg+1)*(p+1)

respFunJacMat = function(betas, aX, reg, theta, muMarg, cReg) {
  mu = exp(crossprod(reg, betas)*psi) * muMarg
  Jac = matrix(0, p*v, p*v)

  v=ncol(reg);p=ncol(aX)
  tmp = t(t((1+t(t(X)/theta))*mu/(1+t(t(mu)/theta))^2)*cReg^2)
  tmp2 =  vapply(seq_len(v), FUN.VALUE = tmp, function(x){reg[,x]*tmp})
  tmp3 =  vapply(seq_len(v), FUN.VALUE = betas, function(x){crossprod(reg,tmp2[,,x])})

  ind = as.logical(Reduce(f=adiag,lapply(seq_len(p),function(x){matrix(TRUE,v,v)})))
  Jac[ind] = apply(tmp3,2,identity)
  return(-Jac)
}