#' A function to plot the distance ratios as boxplots
#'
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param log,las,... see ?boxplot
distRatioPlot = function(resListRows, groupFactor, method = NULL, groupMeth = droplevels(groupsMeth[names(groupsMeth) %in% names(resListRows[[1]])]),log = FALSE, las = 2,...) {
  dist = sapply(resListRows, function(x){
    sapply(x, function(y){distanceFun(y, clusters = groupFactor)$ratio})
  })
  meanDist1df = orderDF(melt(dist, value.name = "DistanceRatio", varnames = c("Method")))
  if(!is.null(method)){
    meanDist1df = subset(meanDist1df, meanDist1df$Method %in% method)
  }
  boxplot(DistanceRatio~Method, main = "Mean ratios of overal to within distance by method", data = meanDist1df, col = groupMeth, ylim = if(log) NULL else c(0, ceiling(max(meanDist1df$DistanceRatio))), las = las,log = ifelse(log,"y",""),...)
  meansDist = tapply(meanDist1df$DistanceRatio, meanDist1df$Method, mean)
  points(meansDist, col="red", pch=18, cex = 1.6)
  abline(h=1, lty = "dashed")
  addLegend(groupMeth = groupMeth, yloc = 5)
}