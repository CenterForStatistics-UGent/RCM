#' A function that prepares the dataframes for the distance ratio plots
#'
#' @param resListRows a list with row scores
#' @param groupFactor a grouping factor defining the true clusters
#' @param method a character vector to subset the methods plotted
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#'
#' @return a dataframe with "DistanceRatio" and "Method" variables
makeDistDf = function(resListRows, groupFactor, method = NULL, groupMeth = droplevels(factor(as.character(groupsMeth[names(groupsMeth) %in% names(resListRows[[1]])])))){
  dist = sapply(resListRows, function(x){
    sapply(x, function(y){distanceFun(y, clusters = groupFactor)$ratio})
  })
  meanDist1df = orderDF(melt(dist, value.name = "DistanceRatio", varnames = c("Method")), groupsMeth = groupMeth)
  if(!is.null(method)){
    meanDist1df = subset(meanDist1df, meanDist1df$Method %in% method)
  }
  meanDist1df
}