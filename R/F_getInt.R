#' A function to integrate the spline of an vgam object between the boundaries of observed environmental scores.
#'
#' @param fitObj a fitted object, either gam or glm
#' @param sampleScore the observed environmental scores
#' @param stop.on.error see ?integrate
#' @param ... additional arguments passed on to integrate()
#'
#' @return a scalar, the value of the integral
#' @importFrom VGAM predict
getInt = function(fitObj, sampleScore, stop.on.error = FALSE,...){
  #Absolute values assure positive outcomes
  integrate(f = function(y, fitObj){
    # requireNamespace("splines")
    if(class(fitObj) %in% c("vgam", "glm")){
    abs(predict(fitObj, type = "link", newdata = data.frame(sampleScore = y, logMu = 0))) #logMu = 0 for departure from uniformity
    } else {stop("GAM and GLM fits failed! \n")}
    },lower = min(sampleScore), upper = max(sampleScore), fitObj = fitObj, stop.on.error = stop.on.error,...)$value
}
