#' A function to calculate and plot the silhouettes as boxplots
#'
#' @param silDF a silhouette dataframe, obtained from makeSimDf()
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param las,... see ?boxplot
#' @param diamCol The color of the diamond for the mean
#' @param bordercol The colour of the borders of the boxplot
#'
plotSil = function(silDF, groupMeth = droplevels(groupsMeth[factorMeth %in% levels(silDF$Method)]),
                   las = 2, diamColor = "orange", bordercol = borderCol,...){
  parTmp = par(no.readonly = TRUE)
  par(mfrow=  c(1,1), mar = c(10,4,4,6))
  boxplot(Silhouette~Method, main = "Mean silhouettes per Monte Carlo simulation by method", data = silDF, col = groupMeth, las = las, border = bordercol, ...)
  meansSil = tapply(silDF$Silhouette, silDF$Method, mean)
  points(meansSil, col=diamColor, pch=18, cex = 1.6)
  abline(h=0, lty = "dashed")
  addLegend(groupMeth = groupMeth, x = length(unique(silDF$Method))+1, yloc = 0.5)
  par(parTmp)
}