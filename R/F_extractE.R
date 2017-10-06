#' A function to extract a matrix of expected values for any dimension of the fit, along with a matrix of dispersions.
#'
#' @param rcm an object of class RCM
#' @param k the desired dimension. Defaults to the maximum of the fit. Choose 0 for the independence model
#'
#' The dispersions of each dimension are extracted if possible, fitted if necessary
#'
#' @return A list with components
#' \item{E}{The matrix of expected values}
#' \item{thetaMat}{The corresponding matrix of overdispersions}
extractE = function(rcm, k = rcm$k){
  #Expectations
  Eind = outer(rcm$libSizes, rcm$abunds) #Expected values under independence
  E = if (k==0){
    Eind
  } else {
    if(is.null(rcm$covariates)){
      Eind *exp(rcm$rMat[,1:k, drop = FALSE] %*% (rcm$cMat[1:k,, drop = FALSE]* rcm$psis[1:k]))
    } else if(rcm$responseFun == "nonparametric"){
      Eind*exp(apply(vapply(1:k,FUN.VALUE = Eind, function(j){rcm$psis[j]*rcm$nonParamRespFun[[j]]$taxonWise}),c(1,2),sum))
    } else {
      Eind*exp(apply(vapply(1:k, FUN.VALUE = Eind, function(j){rcm$psis[j]*getRowMat(sampleScore = rcm$covariates %*% rcm$alpha[,j], responseFun = rcm$responseFun, NB_params = rcm$NB_params[,,j])}),c(1,2),sum))
    }
  }

  #Overdispersions
  thetaMat = if(is.matrix(rcm$thetas)){
    matrix(rcm$thetas[,k], byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
    } else if(k==rcm$k) {
    matrix(rcm$thetas, byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
  } else {
    matrix(estDisp(X = rcm$X, muMarg = E), byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
  }

  list(E = E, thetaMat = thetaMat)
}