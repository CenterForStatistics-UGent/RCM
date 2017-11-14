#' A function to calculate a dataframe of silhouettes
#'
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#'
#' @return a dataframe with Silhouette, Method and
makeSilDf = function(resListRows, groupFactor, method = NULL, groupMeth = droplevels(factor(as.character(groupsMeth[names(groupsMeth) %in% names(resListRows[[1]])]), levels = levels(groupsMeth)))){
  meanSil1 = sapply(resListRows, function(x){
    colMeans(sapply(x, silhouette, clusters = groupFactor))
  })
  meanSil1df = orderDF(melt(meanSil1, value.name = "Silhouette", varnames = c("Method")), groupsMeth = groupMeth)

  if(!is.null(method)){
    meanSil1df = droplevels(subset(meanSil1df, meanSil1df$Method %in% method))
  }
  meanSil1df
}