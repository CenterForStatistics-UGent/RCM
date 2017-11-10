#' A function to plot the result of this as boxplots. Horizontal boxplots have been tried but the overview over the best methods is a bit lost
#'
#' @param contr the contrast object, a result of the contrTaxaWrap() function
#' @param method a subset of biplot methods
#' @param expFac see addLegend
#' @param groupMeth a factor specifying the order of the methods
#' @param las,ylim,... passed on to the boxplot function
#' @param diamCol The color of the diamond for the mean
#' @param bordercol The colour of the borders of the boxplot
#' @param yloc the y-coordinate for the legend
#'
#' @return invisible(), plots the boxplot
contrPlot = function(contr, method = NULL, las = 2, expFac = 1.2, groupMeth = droplevels(factor(as.character(groupsMeth[names(groupsMeth) %in% names(contr)]))), ylim = c(0,max(sapply(contr, quantile, 0.95, na.rm = TRUE))), diamCol = "orange", bordercol = borderCol, yloc = 4,...) {
  moltenContr = orderDF(melt(contr, variable.name = "Method", value.name = "contrRatio", id.vars = NULL), groupsMeth = groupMeth)
  if(!is.null(method)){
    moltenContr = droplevels(subset(moltenContr, moltenContr$Method %in% method))
  }
  meansDist = tapply(moltenContr$contrRatio, moltenContr$Method, mean, na.rm = TRUE)
  at = seq_along(meansDist)*expFac
  boxplot(contrRatio~Method, main = "Mean contribution of DA taxa relative to \n the 50% most contributing non-DA taxa by method", data = moltenContr, ylab = "Mean contribution ratio", col = groupMeth, las = las, at = at, ylim = ylim,border = bordercol,...)
  points(y = meansDist, x = at, col = diamCol, pch=18, cex = 1.6)
  abline(h=1, lty = "dashed")
  addLegend(groupMeth = groupMeth, expFac = expFac, yloc = yloc)
}