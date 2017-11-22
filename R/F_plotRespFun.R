#' A function that plots the non-parametric response functions of a prespecified number of taxa over the observed range of the environmental score.
#'
#' @param RCM and RCM object
#' @param taxa a character vector of taxa to be plotted
#' @param Dim the dimension to be plotted
#' @param nPoints the number of points to be used to plot the lines
#' @param labSize the label size for the variables
#' @param yLoc the y-location of the variables
#' @param ... Other argumens passed on to the ggplot() function
plotRespFun = function(RCM, taxa, Dim = 1, nPoints = 1e3L, labSize = 2.5, yLoc = NULL,...){
require(ggplot2)
require(reshape2)
  if(is.null(RCM$nonParamRespFun)){
    stop("This function can only be called on non-parametric response functions! \n")
  }

  names(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]]) = taxa_names(RCM$physeq)
  # A function to predict new values
  predictFun = function(taxon, x){
    predict(RCM$nonParamRespFun[[paste0("Dim", Dim)]][["taxonWise"]][[taxon]]$fit, newdata = data.frame(sampleScore=x, logMu = 0))
  }

  #The range of sample scores
sampleScoreRange = range(with(RCM, covariates %*% alpha)[,Dim])
sampleScoreSeq = x = seq(sampleScoreRange[1], sampleScoreRange[2], length.out = nPoints)
df = data.frame(sampleScore = sampleScoreSeq, lapply(taxa, predictFun, x = sampleScoreSeq))
names(df)[-1] = taxa
dfMolt = melt(df, id.vars ="sampleScore", value.name = "responseFun", variable.name = "Taxon")

plot = ggplot(data = dfMolt, aes(x = sampleScore, y = responseFun, group = Taxon, colour = Taxon),...) + geom_line()
plot = plot + xlab("Environmental score") + ylab("Response function")


#Also add the associated elements of the environmental gradient in the upper margin
yLoc = if(is.null(yLoc)) max(dfMolt$responseFun)*0.9 else yLoc
textDf = data.frame(text = rownames(RCM$alpha), x = RCM$alpha[,Dim]*min(abs(sampleScoreRange))/max(abs(RCM$alpha[,Dim])), y = yLoc)
plot = plot + geom_text(data = textDf, mapping = aes(x=x, y=y, label=text), inherit.aes = FALSE, angle = 90, size = labSize)
plot
}