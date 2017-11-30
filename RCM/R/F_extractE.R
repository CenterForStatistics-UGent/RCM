#' A function to extract a matrix of expected values for any dimension of the fit
#'
#' @param rcm an object of class RCM
#' @param Dim the desired dimensions. Defaults to the maximum of the fit. Choose 0 for the independence model
#'
#' @return The matrix of expected values
extractE = function(rcm, Dim = seq_len(rcm$k)){
  #Expectations
  Eind = outer(rcm$libSizes, rcm$abunds) #Expected values under independence
  if (Dim[1] %in% c(0,NA)){
    Eind
  } else if(is.null(rcm$covariates)){
      Eind *exp(rcm$rMat[,Dim, drop = FALSE] %*% (rcm$cMat[Dim,, drop = FALSE]* rcm$psis[Dim]))
  } else if(rcm$responseFun == "nonparametric"){
      Eind*exp(apply(vapply(Dim,FUN.VALUE = Eind, function(j){rcm$psis[j]*rcm$nonParamRespFun[[j]]$taxonWise}),c(1,2),sum))
  } else {
      Eind*exp(apply(vapply(Dim, FUN.VALUE = Eind, function(j){rcm$psis[j]*getRowMat(sampleScore = rcm$covariates %*% rcm$alpha[,j], responseFun = rcm$responseFun, NB_params = rcm$NB_params[,,j])}),c(1,2),sum))
    }
  }