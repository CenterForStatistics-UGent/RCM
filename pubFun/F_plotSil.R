#' A function to calculate and plot the silhouettes as boxplots
#'
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param las,... see ?boxplot
plotSil = function(resListRows, groupFactor, method = NULL, groupMeth = droplevels(groupsMeth[names(groupsMeth) %in% names(resListRows[[1]])]),las = 2,...){
  parTmp = par(no.readonly = TRUE)
  par(mfrow=  c(1,1))
  sil1  = lapply(resListRows, function(x){
    sapply(x, silhouette, clusters = groupFactor)
  })
  meanSil1 = sapply(sil1, function(x){
    colMeans(x)
  })
  meanSil1df = orderDF(melt(meanSil1, value.name = "Silhouette", varnames = c("Method")))

  if(!is.null(method)){
    meanSil1df = subset(meanSil1df, meanSil1df$Method %in% method)
  }
  boxplot(Silhouette~Method, main = "Mean silhouettes per Monte Carlo simulation by method", data = meanSil1df, col = groupMeth, las = las, ...)
  meansSil = tapply(meanSil1df$Silhouette, meanSil1df$Method, mean)
  points(meansSil, col="red", pch=18, cex = 1.6)
  abline(h=0, lty = "dashed")
  addLegend(groupMeth = groupMeth)
  par(parTmp)
}