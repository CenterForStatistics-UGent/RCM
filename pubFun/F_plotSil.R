#' A function to calculate and plot the silhouettes as boxplots
#'
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param las,... see ?boxplot
#' @param diamCol The color of the diamond for the mean
#' @param bordercol The colour of the borders of the boxplot
#'
plotSil = function(resListRows, groupFactor, method = NULL,
                   groupMeth = droplevels(factor(as.character(groupsMeth[names(groupsMeth) %in% names(resListRows[[1]])]), levels = levels(groupsMeth))),
                   las = 2, diamColor = "orange", bordercol = borderCol,...){
  parTmp = par(no.readonly = TRUE)
  par(mfrow=  c(1,1))
  meanSil1 = sapply(resListRows, function(x){
    colMeans(sapply(x, silhouette, clusters = groupFactor))
  })
  meanSil1df = orderDF(melt(meanSil1, value.name = "Silhouette", varnames = c("Method")), groupsMeth = groupMeth)

  if(!is.null(method)){
    meanSil1df = droplevels(subset(meanSil1df, meanSil1df$Method %in% method))
  }
  boxplot(Silhouette~Method, main = "Mean silhouettes per Monte Carlo simulation by method", data = meanSil1df, col = groupMeth, las = las, border = bordercol, ...)
  meansSil = tapply(meanSil1df$Silhouette, meanSil1df$Method, mean)
  points(meansSil, col=diamColor, pch=18, cex = 1.6)
  abline(h=0, lty = "dashed")
  addLegend(groupMeth = groupMeth, x = length(unique(meanSil1df$Method))+1, yloc = 0.5)
  par(parTmp)
}