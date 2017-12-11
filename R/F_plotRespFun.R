#' Plot the non-parametric response functions of a prespecified set of taxa
#'
#' @description Plots a number of response functions over the observed range of the environmental score. If no taxa are provided those who react most strongly to the environmental score are chosen.
#'
#' @param RCM and RCM object
#' @param taxa a character vector of taxa to be plotted
#' @param addSamples a boolean, should sample points be shown?
#' @param samShape a sample variable name or a vector of length equal to the number of samples, for the sample shape
#' @param samSize the size of the sample dots
#' @param Dim the dimension to be plotted
#' @param nPoints the number of points to be used to plot the lines
#' @param labSize the label size for the variables
#' @param yLocVar the y-location of the variables, recycled if necessary
#' @param yLocSam the y-location of the samples, recycled if necessary
#' @param Palette which color palette to use
#' @param addJitter should variable names be jittered to make them more readable
#' @param subdivisions the number of subdivisions for the integration
#' @param nTaxa an integer, number of taxa to plot
#' @param angle angle at which variable labels should be turned
#' @param legendLabSize size of the legend labels
#' @param legendTitleSize size of the legend title
#' @param axisLabSize size of the axis labels
#' @param axisTitleSize size of the axis title
#' @param ... Other argumens passed on to the ggplot() function
#'
#' @export
#' @import ggplot2
#' @import phyloseq
#' @importFrom stats runif
plotRespFun = function(RCM, taxa = NULL, addSamples = TRUE, samShape = NULL, samSize = 2, Dim = 1, nPoints = 1e3L, labSize = 2.5, yLocVar = NULL, yLocSam =  NULL, Palette = "Set3", addJitter = FALSE, subdivisions = 50L, nTaxa = 8L, angle = 90, legendLabSize = 15,  legendTitleSize = 16, axisLabSize = 14, axisTitleSize = 16,...){
  if(is.null(RCM$nonParamRespFun)){
    stop("This function can only be called on non-parametric response functions! \n")
  }
  names(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]]) = taxa_names(RCM$physeq)
  # A function to predict new values
  predictFun = function(taxon, x){
    predict(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]][[taxon]]$fit, newdata = data.frame(sampleScore=x, logMu = 0))}

  #The range of sample scores
sampleScoreRange = range(with(RCM, covariates %*% alpha)[,Dim])
sampleScoreSeq = x = seq(sampleScoreRange[1], sampleScoreRange[2], length.out = nPoints)
if(is.null(taxa)) { #If taxa not provided, pick the ones that react most strongly
  intsNonParam = sapply(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]], function(x){getInt(x$fit, sampleScore =c( RCM$covariates %*% RCM$alpha[,1]), subdivisions = subdivisions)}) #The integrals
  taxa = taxa_names(RCM$physeq)[intsNonParam > quantile(intsNonParam, (ntaxa(RCM$physeq)-nTaxa)/ntaxa(RCM$physeq))]
}

df = data.frame(sampleScore = sampleScoreSeq, lapply(taxa, predictFun, x = sampleScoreSeq))
names(df)[-1] = taxa
dfMolt = reshape2::melt(df, id.vars ="sampleScore", value.name = "responseFun", variable.name = "Taxon")

plot = ggplot(data = dfMolt, aes_string(x = "sampleScore", y = "responseFun", group = "Taxon", colour = "Taxon"),...) + geom_line()
plot = plot + xlab(paste("Environmental score of dimension", Dim)) + ylab("Response function")

#Also add the associated elements of the environmental gradient in the upper margin
textDf = data.frame(text = rownames(RCM$alpha), x = RCM$alpha[,Dim]*min(abs(sampleScoreRange))/max(abs(RCM$alpha[,Dim])))
textDf = textDf[order(textDf$x),]
textDf$y = (if(is.null(yLocVar)) (max(dfMolt$responseFun)+min(dfMolt$responseFun))/2 else yLocVar)
plot = plot + geom_text(data = textDf, mapping = aes_string(x = "x", y = "y", label = "text"), inherit.aes = FALSE, angle = angle, size = labSize) + scale_colour_brewer(palette = Palette)
#Finally add a dashed line for the independence model, and a straight line for the 0 environmental gradient
plot = plot + geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) + geom_vline(xintercept = 0, size = 0.15, linetype = "dotted")

#If samples required, add them too, as marks
if(addSamples){
dfSam = data.frame(x = RCM$covariates %*% RCM$alpha[,Dim])
dfSam$Shape = if(length(samShape)==1) get_variable(RCM$physeq, varName = samShape) else if(length(samShape)==ncol(RCM$X)) samShape else NULL
# dfSam$Fill = if(length(samColour)==1) get_variable(RCM$physeq, varName = samColour) else if(length(samColour)==ncol(RCM$X)) samColour else NULL
dfSam$y = (if(is.null(yLocSam)) min(dfMolt$responseFun)*0.8 else yLocSam) + if(addJitter) runif(min = -1,max =1, n = nrow(RCM$alpha))*diff(range(dfMolt$responseFun))/20 else 0
if(!is.null(samShape)){
plot = plot + geom_point(inherit.aes = FALSE, fill = NA, mapping = aes_string(x = "x", y = "y", shape = "Shape"), data = dfSam, size = samSize)
} else {
  plot = plot + geom_point(inherit.aes = FALSE, fill = NA, mapping = aes_string(x = "x", y = "y"),shape = 1, data = dfSam, size = samSize)
}
plot = plot + scale_shape_discrete(name = if(!is.null(samShape)) samShape else "", solid = FALSE)
}
#Adapt the text sizes
plot = plot + theme(axis.title = element_text(size = axisTitleSize), axis.text = element_text(size = axisLabSize), legend.title = element_text(size = legendTitleSize), legend.text = element_text(size = legendLabSize))

plot
}
