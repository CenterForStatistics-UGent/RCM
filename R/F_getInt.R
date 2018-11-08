#' A function to integrate the spline of an vgam object between the boundaries of observed environmental scores.
#'
#' @param fitObj a fitted object, either gam or a vector of parameters
#' @param sampleScore the observed environmental scores
#' @param stop.on.error see ?integrate
#' @param class the class of the fit
#' @param ... additional arguments passed on to integrate()
#'
#' @return a scalar, the value of the integral
#' @importFrom VGAM predict
getInt = function(coef, spline, sampleScore, stop.on.error = FALSE,...){
  #Absolute values assure positive outcomes
  integrate(f = function(y, coef, spline){
    if(!is.null(spline)){
    abs(getRowMat(sampleScore = y, taxonCoef = coef, spline = spline, responseFun = "nonparametric")) #logMu = 0 for departure from uniformity
    } else {#If GAM fails, GLM fit (or independence model)
abs(getModelMat(y, degree = length(coef)-1) %*% coef)
    }
    }, lower = min(sampleScore), upper = max(sampleScore), taxonCoef = coef, splinesList = spline, stop.on.error = stop.on.error,...)$value
}
