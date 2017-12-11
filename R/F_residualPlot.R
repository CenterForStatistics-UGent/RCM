#' Make residual plots of the taxa with the strongest response to the environmental gradient, or with the highest run statistic
#'
#' @param RCM an RCM object
#' @param Dim an integer, which dimension?
#' @param whichTaxa a character string or a character vector, for which taxa to plot the diagnostic plots
#' @param resid the type of residuals to use, either "Deviance" or "Pearson"
#' @param numTaxa an integer, the number of taxa to plot
#' @param mfrow passed on to par(). If not supplied will be calculated based on numTaxa
#' @param samColour,samShape Vectors or character strings denoting the sample colour and shape respectively. If character string is provided, the variables with this name is extracted from the phyloseq object in RCM
#' @param legendLabSize size of the legend labels
#' @param legendTitleSize size of the legend title
#' @param axisLabSize size of the axis labels
#' @param axisTitleSize size of the axis title
#' @param taxTitle A boolean, should taxon title be printed
#'
#'@details If whichTaxa is "run" or "response" the taxa with the highest run statistics or steepest slopes of the response function are plotted, numTax indicates the number. If whichTaxa is a character vector, these are interpreted as taxon names to plot. This function is mainly meant for linear response functions, but can be used for others too. The runs test statistic from the tseries package is used.
#'@return Plots a ggplot2-object to output
#'@export
#'@import ggplot2
#'@import phyloseq
#'@importFrom tseries runs.test
residualPlot = function(RCM, Dim = 1, whichTaxa = "response", resid = "Deviance", numTaxa = 9, mfrow = NULL, samColour = NULL, samShape = NULL, legendLabSize = 15,  legendTitleSize = 16, axisLabSize = 14, axisTitleSize = 16, taxTitle = TRUE){
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
  sizes = if(whichTaxa =="runs") apply(resMat>0,2, function(x){runs.test(factor(x))$statistic}) else RCM$NB_params[2,,Dim] #Select taxa with longest runs or strongest responses
  idTaxa = which(sizes >= sort(sizes, decreasing = TRUE)[numTaxa])
}
  #Prepare the plotting facets
  mfrow = if(is.null(mfrow)) rep(ceiling(sqrt(numTaxa)),2) else mfrow
  parTmp = par(no.readonly = TRUE)
  par(mfrow = mfrow)
  resMat = resMat[,idTaxa, drop = FALSE]
  if(length(idTaxa)>1){
  sapply(colnames(resMat), function(tax){
    plot(x = sampleScore, y = resMat[,tax], ylab = paste(resid,"residuals"), xlab = paste("Environmental score in dimension",Dim), main = tax)
  })
  } else {
    Plot = ggplot(data.frame(x = c(sampleScore), y = c(resMat), samColour = if(is.null(samColour)) "black"  else get_variable(RCM$physeq, samColour), samShape = if(is.null(samShape)) "none"  else get_variable(RCM$physeq, samShape)), mapping = aes_string(x="x", y="y", colour = "samColour", shape = "samShape"))+ ylab(paste(resid,"residuals")) + xlab(paste("Environmental score in dimension",Dim)) + geom_point() + ggtitle(ifelse(taxTitle,colnames(resMat), ""))
    Plot = Plot + if(is.null(samShape)) guides(shape = FALSE) else scale_shape_discrete(name = samShape)
    Plot = Plot + if(is.null(samColour)) guides(colour = FALSE) else if(is.factor(get_variable(RCM$physeq, samColour))) scale_colour_discrete(name = samColour) else scale_colour_continuous(name = samColour)
    Plot = Plot +   # Enlarge most text
 theme(axis.title = element_text(size = axisTitleSize), axis.text = element_text(size = axisLabSize), legend.title = element_text(size = legendTitleSize), legend.text = element_text(size = legendLabSize))
    par(parTmp)
    return(Plot)}
  par(parTmp)
}
