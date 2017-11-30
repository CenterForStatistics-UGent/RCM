#' Plot the non-parametric response functions of a prespecified set of taxa
#'
#' @description Plots a number of response functions over the observed range of the environmental score. If no taxa are provided those who react most strongly to the environmental score are chosen.
#'
#' @param RCM and RCM object
#' @param taxa a character vector of taxa to be plotted
#' @param Dim the dimension to be plotted
#' @param nPoints the number of points to be used to plot the lines
#' @param labSize the label size for the variables
#' @param yLoc the y-location of the variables
#' @param Palette which color palette to use
#' @param adJitter should variable names be jittered to make them more readable
#' @param subdivisions the number of subdivisions for the integration
#' @param nTaxa an integer, number of taxa to plot
#' @param angle angle at which variable labels should be turned
#' @param ... Other argumens passed on to the ggplot() function
#' @export
#' @import ggplot2
#' @import phyloseq
plotRespFun = function(RCM, taxa = NULL, Dim = 1, nPoints = 1e3L, labSize = 2.5, yLoc = NULL, Palette = "Set3", adJitter = FALSE, subdivisions = 100L, nTaxa = 8L, angle = 90,...){
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

plot = ggplot(data = dfMolt, aes(x = sampleScore, y = responseFun, group = Taxon, colour = Taxon),...) + geom_line()
plot = plot + xlab("Environmental score") + ylab("Response function")

#Also add the associated elements of the environmental gradient in the upper margin
textDf = data.frame(text = rownames(RCM$alpha), x = RCM$alpha[,Dim]*min(abs(sampleScoreRange))/max(abs(RCM$alpha[,Dim])))
textDf = textDf[order(textDf$x),]
textDf$y = (if(is.null(yLoc)) max(dfMolt$responseFun)*0.9 else yLoc) + if(adJitter) rep(c(-1,0,1)*diff(range(dfMolt$responseFun))/8, length.out = nrow(RCM$alpha)) else 0
plot = plot + geom_text(data = textDf, mapping = aes(x=x, y=y, label=text), inherit.aes = FALSE, angle = angle, size = labSize) + scale_colour_brewer(palette = Palette)
#Finally add a dashed line for the independence model
plot = plot + geom_hline(yintercept=0, linetype = "dashed", size = 0.3)
plot
}
