#' Return a matrix of row scores
#'
#' @param sampleScore a vector of length n with sample scores
#' @param responseFun a character string, the type of response function,
#'  either "linear" or "quadratic"
#' @param NB_params a v-by-p matrix of parameters of theresponse function
#' @param taxonCoef A vector of coefficients
#' @param spline The cubic smoothing spline
#'
#' Multiplying the old offset with the exponent matrix times
#' the importance parameter obtains the new one based on lower dimension
#'
#' @return a n-by-p matrix of scores
#' @importFrom VGAM predict
getRowMat = function(sampleScore, responseFun, NB_params, taxonCoef, spline){
  if(responseFun=="nonparametric"){
      cbind(1,sampleScore, predict(spline, x = sampleScore)$y)%*% c(taxonCoef,1)
  } else {
    buildDesign(sampleScore, responseFun) %*% NB_params
  }
}
