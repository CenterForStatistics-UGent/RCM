#' Integrate the spline of an vgam object
#'
#' @param coef A vector of coefficients
#' @param spline The cubic smoothing spline
#' @param sampleScore the observed environmental scores
#' @param stop.on.error see ?integrate
#' @param ... additional arguments passed on to integrate()
#'
#' @return a scalar, the value of the integral
getInt = function(coef, spline, sampleScore, stop.on.error = FALSE,...){
  #Absolute values assure positive outcomes
  integrate(f = function(y, coef, spline){
    if(!is.null(spline)){
    abs(getRowMat(sampleScore = y, taxonCoef = coef, spline = spline,
                  responseFun = "nonparametric"))
    } else {#If GAM fails, GLM fit (or independence model)
abs(getModelMat(y, degree = length(coef)-1) %*% coef)
    }
    }, lower = min(sampleScore), upper = max(sampleScore), coef = coef,
spline = spline, stop.on.error = stop.on.error,...)$value
}
