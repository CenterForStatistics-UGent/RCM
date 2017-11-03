#' A function to add the legend
addLegend = function(yloc = 0.2, groupMeth, x = length(groupMeth)+1, expFac = 1,...){
  legend(legend = levels(groupMeth), col = factor(levels(groupMeth), levels = levels(groupMeth), ordered=TRUE),pch = 16, x = x*expFac, y = yloc, xpd = TRUE,...)
}