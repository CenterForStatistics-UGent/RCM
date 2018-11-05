#' Plot the non-parametric response functions of a prespecified set of taxa
#'
#' @description Plots a number of response functions over the observed range of the environmental score. If no taxa are provided those who react most strongly to the environmental score are chosen.
#'
#' @param RCM an RCM object
#' @param taxa a character vector of taxa to be plotted
#' @param type a character string, plot the response function on the log-scale ("link") or the abundance scale "response", similar to predict.glm().
#' @param logTransformYAxis a boolean, should y-axis be log transformed?
#' @param addSamples a boolean, should sample points be shown?
#' @param samSize a sample variable name or a vector of length equal to the number of samples, for the sample sizes
#' @param Dim An integer, the dimension to be plotted
#' @param nPoints the number of points to be used to plot the lines
#' @param labSize the label size for the variables
#' @param yLocVar the y-location of the variables, recycled if necessary
#' @param yLocSam the y-location of the samples, recycled if necessary
#' @param Palette which color palette to use
#' @param addJitter A boolean, should variable names be jittered to make them more readable
#' @param nTaxa an integer, number of taxa to plot
#' @param angle angle at which variable labels should be turned
#' @param legendLabSize size of the legend labels
#' @param legendTitleSize size of the legend title
#' @param axisLabSize size of the axis labels
#' @param axisTitleSize size of the axis title
#' @param lineSize size of the response function lines
#' @param ... Other argumens passed on to the ggplot() function
#'
#' @return Plots a ggplot2-object to output
#' @export
#' @import ggplot2
#' @import phyloseq
#' @importFrom stats runif
#' @seealso \code{\link{RCM}}
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' #Subset for a quick fit
#' zellerRCMnp = RCM(tmpPhy, k = 2, covariates = c("BMI","Age","Country","Diagnosis","Gender"),
#' round = TRUE, responseFun = "nonparametric")
#' plotRespFun(zellerRCMnp)
plotRespFun = function(RCM, taxa = NULL, type = "link", logTransformYAxis = FALSE, addSamples = TRUE, samSize = NULL, Dim = 1L, nPoints = 1e2L, labSize = 2.5, yLocVar = NULL, yLocSam =  NULL, Palette = "Set3", addJitter = FALSE, nTaxa = 8L, angle = 90, legendLabSize = 15,  legendTitleSize = 16, axisLabSize = 14, axisTitleSize = 16, lineSize = 0.75,...){
  if(is.null(RCM$nonParamRespFun)){
    stop("This function can only be called on non-parametric response functions! \n")
  }
  if(!type %in% c("link","response")){
    stop("Specify type = 'link' or 'response'!\n")
  }
  names(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]]) = taxa_names(RCM$physeq)
  # A function to predict new values
  predictFun = function(taxon, x){
    predict(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]][[taxon]]$fit, newdata = data.frame(sampleScore=x, logMu = 0))} #Fix me!

  #The range of sample scores
sampleScoreRange = range(RCM$covariates %*% RCM$alpha[,Dim])
sampleScoreSeq = seq(sampleScoreRange[1], sampleScoreRange[2], length.out = nPoints)
if(is.null(taxa)) { #If taxa not provided, pick the ones that react most strongly
  intsNonParam = sapply(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]], function(x){x$int}) #The integrals
  taxa = taxa_names(RCM$physeq)[intsNonParam > quantile(intsNonParam, (ntaxa(RCM$physeq)-nTaxa)/ntaxa(RCM$physeq))]
}

df = data.frame(sampleScore = sampleScoreSeq, lapply(taxa, predictFun, x = sampleScoreSeq))
names(df)[-1] = taxa
dfMolt = reshape2::melt(df, id.vars ="sampleScore", value.name = "responseFun", variable.name = "Taxon")
if(type=="response"){dfMolt$responseFun = exp(dfMolt$responseFun) +1e-9}

plot = ggplot(data = dfMolt, aes_string(x = "sampleScore", y = "responseFun", group = "Taxon", colour = "Taxon"),...) + geom_line(size = lineSize) + xlab(paste("Environmental score of dimension", Dim)) + ylab(ifelse(type=="link","Response function", "Response function on count scale")) + scale_y_continuous(trans = if (logTransformYAxis) "log10" else "identity", labels = if (logTransformYAxis) function(x) {sprintf("%.4f", x)} else identity)

#Also add the associated elements of the environmental gradient in the upper margin
textDf = data.frame(text = rownames(RCM$alpha), x = RCM$alpha[,Dim]*min(abs(sampleScoreRange))/max(abs(RCM$alpha[,Dim])))
textDf = textDf[order(textDf$x),]
textDf$y = (if(is.null(yLocVar)) (max(dfMolt$responseFun)+min(dfMolt$responseFun))/2 else yLocVar)
plot = plot + geom_text(data = textDf, mapping = aes_string(x = "x", y = "y", label = "text"), inherit.aes = FALSE, angle = angle, size = labSize) + scale_colour_brewer(palette = Palette)
#Finally add a dashed line for the independence model, and a straight line for the 0 environmental gradient
plot = plot + geom_hline(yintercept = switch(type, "link" = 0, "response" =1), linetype = "dashed", size = 0.3) + geom_vline(xintercept = 0, size = 0.15, linetype = "dotted")

#If samples required, add them too, as marks of different heights
if(addSamples){
dfSam = data.frame(x = RCM$covariates %*% RCM$alpha[,Dim])
dfSam$Size = if(length(samSize)==1) get_variable(RCM$physeq, varName = samSize) else if(length(samSize)==ncol(RCM$X)) samSize else NULL
# dfSam$Fill = if(length(samColour)==1) get_variable(RCM$physeq, varName = samColour) else if(length(samColour)==ncol(RCM$X)) samColour else NULL
dfSam$y = (if(is.null(yLocSam)) min(dfMolt$responseFun)*0.8 else yLocSam) + if(addJitter) runif(min = -1,max =1, n = nrow(RCM$alpha))*diff(range(dfMolt$responseFun))/20 else 0
if(!is.null(samSize)){
plot = plot + geom_point(inherit.aes = FALSE, fill = NA, mapping = aes_string(x = "x", y = "y", size = "Size"), data = dfSam, shape = 124)
} else {
  plot = plot + geom_point(inherit.aes = FALSE, fill = NA, mapping = aes_string(x = "x", y = "y"),shape = 124, data = dfSam, size = 1)
}
plot = plot + scale_size_discrete(name = if(!is.null(samSize)) samSize else "")
}
#Adapt the text sizes
plot = plot + theme_bw() + theme(axis.title = element_text(size = axisTitleSize), axis.text = element_text(size = axisLabSize), legend.title = element_text(size = legendTitleSize), legend.text = element_text(size = legendLabSize))

plot
}
