#' A function that returns the value of the partial derivative of the log-likelihood ratio to alpha, keeping the response functions fixed
#'
#' @param Alpha: a vector of length d + k*(2+(k-1)/2), the environmental gradient plus the lagrangian multipliers
#' @param X the n-by-p count matrix
#' @param CC a n-by-d covariate vector
#' @param responseFun a character string indicating the type of response function
#' @param psi a scalar, an importance parameter
#' @param NB_params Starting values for the NB_params
#' @param NB_params_noLab Starting values for the NB_params without label
#' @param d an integer, the number of covariate parameters
#' @param alphaK a matrix of environmental gradients of lower dimensions
#' @param k an integer, the current dimension
#' @param centMat a nLambda1s-by-d centering matrix
#' @param nLambdas an integer, number of lagrangian multipliers
#' @param nLambda1s an integer, number of centering restrictions
#' @param thetaMat a matrix of size n-by-p with estimated dispersion parameters
#' @param muMarg an n-by-p offset matrix
#' @param ncols a scalar, the number of columns of X
#'
#' @return: The value of the lagrangian and the constraining equations
dLR_nb <- function(Alpha, X, CC, responseFun = c("linear", "quadratic", "nonparametric","dynamic"), psi, NB_params, NB_params_noLab, d, alphaK, k, centMat, nLambda, nLambda1s, thetaMat, muMarg, ncols, ...){

  #Extract the parameters
  alpha = Alpha[seq_len(d)]
  lambda1s = Alpha[d+seq_len(nLambda1s)] #Multiple centering requirements now!
  lambda2 = Alpha[d+nLambda1s+1]
  lambda3 = if(k==1) {0} else {Alpha[(d+nLambda1s+2):(d+nLambda)]}

  sampleScore = CC %*% alpha #A linear combination of the environmental variables yields the sampleScore
  design = switch(responseFun,
                  linear = model.matrix(~ sampleScore),#With intercept
                  quadratic = model.matrix(~ sampleScore + I(sampleScore^2)),
                  stop("Non-parametric model not yet implemented")
  )

  mu = muMarg * exp(design %*% NB_params *psi)
  mu0 = muMarg * c(exp(design %*% NB_params_noLab*psi))
  tmp = (X-mu)/(1+mu/thetaMat)
  tmp0 = (X-mu0)/(1+mu0/thetaMat)

  lag = switch(responseFun, #The lagrangian depends on the shape of the response function
               "linear" = psi * (crossprod(CC, tmp) %*% (NB_params[2,]) - rowSums(crossprod(CC, tmp0*NB_params_noLab[2]))) ,
               "quadratic" = psi * (
                 c(crossprod(CC, tmp) %*% (NB_params[2,])) +
                   c(crossprod(CC * c(sampleScore), tmp) %*% (NB_params[3,]) * 2)  -
                   rowSums(crossprod(CC, tmp0) *NB_params_noLab[2]) -
                   rowSums(crossprod(CC * c(sampleScore), tmp0) * NB_params_noLab[3])*2),
               stop("Unknown response function provided! \n")) + #Restrictions do not depend on response function
    c(lambda1s %*% centMat) +
    lambda2 * 2 * alpha +
    if(k>1) alphaK %*% lambda3 else 0

  centerFactors = centMat %*% alpha #Includes overal centering
  size = sum(alpha^2)-1
  if(k==1) { return(c(lag, centerFactors, size))}
  ortho = crossprod(alphaK ,alpha)
  c(lag, centerFactors, size, ortho)
}