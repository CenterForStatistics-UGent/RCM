#' A function to make residual plots of the taxa with the strongest response to the environmental gradient, or with the highest run statistic
#'
#' @param RCM an RCM object
#' @param Dim an integer, which dimension?
#' @param whichTaxa a character string or a character vector, for which taxa to plot the diagnostic plots
#' @param resid the type of residuals to use, either "Deviance" or "Pearson"
#' @param numTaxa an integer, the number of taxa to plot
#' @param mfrow passed on to par(). If not supplied will be calculated based on numTaxa
#'
#'If whichTaxa is "run" or "response" the taxa with the highest run statistics or responses are plotted, numTax indicates the number. If whichTaxa is a character vector, these are interpreted as taxon names to plot. This function is mainly meant for linear response functions, but can be used for others too. The runs test statistic from the tseries package is used.
residualPlot = function(RCM, Dim = 1, whichTaxa = "response", resid = "Deviance", numTaxa = 9, mfrow = NULL){
  require(tseries)
  sampleScore = RCM$covariates %*% RCM$alpha[,Dim, drop = FALSE]
  if(resid == "Deviance"){
    resMat = getDevianceRes(RCM, seq_len(Dim))
    } else if(resid == "Pearson") {
      mu = extractE(RCM, seq_len(Dim)) #Residuals are also based on lower dimensions
      thetaMat = extractDisp(RCM, mu)
    resMat = (RCM$X-mu)/sqrt(mu+mu^2/thetaMat)
      } else {stop("Unknown residual type!")}
if(!whichTaxa %in% c("runs","response")){
  numTaxa = length(whichTaxa)
  idTaxa = whichTaxa
} else {
  sizes = if(whichTaxa =="runs") apply(resMat>0,2, function(x){tseries::runs.test(factor(x))$statistic}) else RCM$NB_params[2,,Dim] #Select taxa with longest runs or strongest responses
  idTaxa = which(sizes >= sort(sizes, decreasing = TRUE)[numTaxa])
}
  #Prepare the plotting facets
  mfrow = if(is.null(mfrow)) rep(ceiling(sqrt(numTaxa)),2) else mfrow
  parTmp = par(no.readonly = TRUE)
  par(mfrow = mfrow)
  resMat = resMat[,idTaxa]
  sapply(colnames(resMat), function(tax){
    plot(x = sampleScore, y = resMat[,tax], ylab = paste(resid,"residuals"), xlab = "Environmental score", main = tax)
  })
  par(parTmp)
}