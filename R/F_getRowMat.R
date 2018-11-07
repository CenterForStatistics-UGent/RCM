#' A function that returns a matrix of scores resulting from a response function and an environmental gradient.
#'
#' @param sampleScore a vector of length n with sample scores
#' @param responseFun a character string, the type of response function, either "linear" or "quadratic"
#' @param NB_params a v-by-p matrix of parameters of theresponse function
#' @param nonParFit A list of non parametric fits as returned by vgam.edit
#'
#' Multiplying the old offset with the exponent matrix times the importance parameter obtains the new one based on lower dimension
#'
#' @return a n-by-p matrix of scores
getRowMat = function(sampleScore, responseFun, NB_params, nonParFit){
  if(responseFun=="nonparametric"){
    sapply(nonParFit, function(Fit){
      Fit = if(is.null(Fit$fit)) Fit else Fit$fit
      if(is.numeric(Fit)){
        cbind(1,sampleScore) %*% Fit
      } else {
      cbind(1,sampleScore, predict(Fit$spline, x = sampleScore)$y) %*% c(Fit$coef,1)}
      })
  } else {
    buildDesign(sampleScore, responseFun) %*% NB_params
  }

}
