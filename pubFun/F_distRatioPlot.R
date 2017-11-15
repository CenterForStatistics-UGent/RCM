#' A function to plot the distance ratios as boxplots
#'
#' @param distDF the distance dataframe as resulting from the makeDistDf function
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param diamCol The color of the diamond for the mean
#' @param log,las,... see ?boxplot
distRatioPlot = function(distDF, groupMeth = droplevels(groupsMeth[factorMeth %in% levels(distDF$Method)]), log = FALSE, las = 2, diamCol = "orange",...) {
  parTmp = par(no.readonly = TRUE)
  par(mfrow=  c(1,1), mar = c(8,4,4,4))
  boxplot(DistanceRatio~Method, main = "Mean ratios of overal to within distance by method", data = distDF, col = groupMeth, ylim = if(log) NULL else c(min(c(0.9,tapply(distDF$DistanceRatio, distDF$Method, quantile, 0.05, na.rm = TRUE))),max(tapply(distDF$DistanceRatio, distDF$Method, quantile, 0.95, na.rm = TRUE))), las = las,log = ifelse(log,"y",""),...)
  meansDist = tapply(distDF$DistanceRatio, distDF$Method, mean)
  points(meansDist, col = diamCol, pch = 18, cex = 1.6)
  abline(h=1, lty = "dashed")
  addLegend(groupMeth = groupMeth, yloc = quantile(distDF$DistanceRatio, 0.95))
  par(parTmp)
}