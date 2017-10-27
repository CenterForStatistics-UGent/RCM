#' A function that the components of the influence functions of the environmental gradient.
#'
#' @param rcm an rcm object
#' @param Dim the required dimension
#'
#' @return A list with components
#' \item{score}{a matrix with components of the score function}
#' \item{InvJac}{A square matrix of dimension n with the components of the Jacobian related to the alphas}
NBalphaInfl = function(rcm, Dim){
  if(length(Dim)>1) {stop("Influence of only one dimension at the time can be extratced! \n")}
  #Extract the parameters
  alpha = rcm$alpha[, Dim]
  lambdas = rcm$lambdasAlpha[seq_k(Dim)]
  nLambda1s = NROW(rcm$centMat)
  lambda1s = lambdas[seq_len(nLambda1s)] #Multiple centering requirements now!
  lambda2 = lambdas[nLambda1s+1]
  lambda3 = if(Dim==1) {0} else {lambdas[-seq_len(nLambda1s)]}

  sampleScore = rcm$covariates %*% alpha #A linear combination of the environmental variables yields the sampleScore
  mu = muMarg * extractE(rcm, seq_len(Dim))
  tmp = (X-mu)/(1+mu/thetaMat)
  responseFun = switch(responseFun, dynamic = "quadratic", responseFun)

  if(envGradEst == "LR"){
    mu0 = muMarg * c(exp(design %*% NB_params_noLab*psi))
    tmp0 = (X-mu0)/(1+mu0/thetaMat)
  }
  lag = switch(responseFun, #The lagrangian depends on the shape of the response function
               "linear" = if(envGradEst == "LR"){psi * (crossprod(CC, tmp) %*% (NB_params[2,])  - rowSums(crossprod(CC, tmp0*NB_params_noLab[2])))} else {psi * (crossprod(CC, tmp) %*% (NB_params[2,]))},
               "quadratic" = if(envGradEst == "LR"){psi * (
                 c(crossprod(CC, tmp) %*% (NB_params[2,])) +
                   c(crossprod(CC * c(sampleScore), tmp) %*% (NB_params[3,]) * 2)  -
                   rowSums(crossprod(CC, tmp0) *NB_params_noLab[2]) -
                   rowSums(crossprod(CC * c(sampleScore), tmp0) * NB_params_noLab[3])*2)} else {
                     psi * (c(crossprod(CC, tmp) %*% (NB_params[2,])) +
                              c(crossprod(CC * c(sampleScore), tmp) %*% (NB_params[3,]) * 2))
                   },
               stop("Unknown response function provided! \n")) + #Restrictions do not depend on response function
    c(lambda1s %*% centMat) +
    lambda2 * 2 * alpha +
    if(k>1) alphaK %*% lambda3 else 0
}