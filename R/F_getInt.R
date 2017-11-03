#' A function to integrate the spline of an vgam object between the boundaries of observed environmental scores.
#'
#' @param fitObj a fitted object, either gam or glm
#' @param sampleScore the observed environmental scores
#' @param stop.on.error see ?integrate
#'
#' @return a scalar, the value of the integral
getInt = function(fitObj, sampleScore, stop.on.error = FALSE){
  #Absolute values assure positive outcomes
  integrate(f = function(y, fitObj){
    if(class(fitObj)=="vgam"){
    abs(predict(fitObj, type = "link", newdata = data.frame(sampleScore = y, logMu = 0))) #logMu = 0 for departure from uniformity
    } else if(class(fitObj) == "list"){
      c(abs(model.matrix(~ y + I(y^2) + I(y^3)) %*% fitObj$coef))
    } else {stop("GLM fit failed! \n")}
    },lower = min(sampleScore), upper = max(sampleScore), fitObj = fitObj, stop.on.error = stop.on.error)$value
}