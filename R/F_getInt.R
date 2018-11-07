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
getInt = function(fitObj, sampleScore, class, stop.on.error = FALSE,...){
  #Absolute values assure positive outcomes
  integrate(f = function(y, fitObj){
    if(class == "vgam"){
    abs(getRowMat(sampleScore = y, nonParFit = list(fit =fitObj), responseFun = "nonparametric")) #logMu = 0 for departure from uniformity
    } else if(class == "glm"){#If GAM fails, GLM fit
abs(getModelMat(y, degree = length(fitObj)-1) %*% fitObj)
    } else {#If GAM and GLM fail, fit indepedence model!
      rep.int(0L, length(y))
      }
    }, lower = min(sampleScore), upper = max(sampleScore), fitObj = fitObj, stop.on.error = stop.on.error,...)$value
}
