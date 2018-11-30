#' This function adds orthogonal projections to a given plot
#'
#' @param RCMplot the RCMplot object
#' @param sample,species,variable names or approximate coordinates of sample, species or variable
#' @param Dims The dimensions of the solutions that have been plotted
#' @param addLabel a boolean, should the r-s-psi label be added?
#' @param labPos the position of the label. Will be calculated if not provided
#'
#' @return a modified ggplot object that contains the geom_segment object that draws the projection
#' @export
#' @import ggplot2
#' @import phyloseq
#' @seealso \code{\link{plot.RCM}}
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(100)],
#' prune_samples(sample_names(Zeller)[seq_len(50)], Zeller))
#' zellerRCM = RCM(tmpPhy, k = 2, round = TRUE)
#' zellerPlot = plot(zellerRCM, returnCoords = TRUE)
#' addOrthProjection(zellerPlot, species = c(-0.35,1.1), sample = c(1,1.2))
addOrthProjection = function(RCMplot, sample = NULL, species = NULL, variable = NULL, Dims = c(1,2), addLabel = FALSE, labPos = NULL){
  nulls = is.null(sample) + is.null(species) + is.null(variable)
  if(nulls != 1) stop("Provide two variables categories for a projection! \n")
  if(is.null(species)) stop("Species should be provided, cannot project sample onto variable vector! \n")

  dimNames = paste0("Dim",Dims)
  if(is.numeric(sample)){
  samp = which.min(colSums((t(RCMplot$samples[, dimNames])-sample)^2)) #Closest to approximate coordinate
  sampName = rownames(RCMplot$samples)[samp]
  } else {
    sampName = sample
  }

  if(is.numeric(species)){
    species = which.min(colSums((t(RCMplot$species[, paste0("end", Dims)])-species)^2)) #Closest to approximate coordinate
    speciesName = rownames(RCMplot$species)[species]
  } else {
    speciesName = species
  }

  if(is.numeric(variable)){
    variable = which.min(colSums((t(RCMplot$variables[, dimNames])-species)^2)) #Closest to approximate coordinate
    varName = rownames(RCMplot$variables)[variable]
  } else {
    varName = variable
  }

  mat1 = unlist(if(is.null(variable)) RCMplot$samples[sampName, dimNames] else RCMplot$variables[varName, dimNames])
  mat2 = unlist(RCMplot$species[speciesName, c(vapply(Dims, FUN.VALUE =character(2), function(x){c(paste0("end",x), paste0("origin",x))}))])

  RCMplot$plot = RCMplot$plot + geom_segment(inherit.aes = FALSE, mapping = aes_string(x=0, y= 0, xend = "Dim1", yend = "Dim2"), data = data.frame(t(mat1))) #The sample or variable vector
  IntCoordsXTip = (mat2["end1"] + mat2["end2"]*mat1[2]/mat1[1])/((mat1[2]/mat1[1])^2+1)
  IntCoordsYTip = IntCoordsXTip*mat1[2]/mat1[1]

  IntCoordsXStart = (mat2["origin1"] + mat2["origin2"]*mat1[2]/mat1[1])/((mat1[2]/mat1[1])^2+1)
  IntCoordsYStart = IntCoordsXStart*mat1[2]/mat1[1]

  dfTip = data.frame(x = mat2[grep(names(mat2), pattern = "end")][1], y = mat2[grep(names(mat2), pattern = "end")][2], xend = IntCoordsXTip, yend = IntCoordsYTip)
  dfStart = data.frame(x = mat2[grep(names(mat2), pattern = "origin")][1], y = mat2[grep(names(mat2),pattern =  "origin")][2], xend = IntCoordsXStart, yend = IntCoordsYStart)

  RCMplot$plot = RCMplot$plot + geom_segment(inherit.aes = FALSE, mapping = aes_string(x="x", y= "y", xend = "xend", yend = "yend"), data = dfTip, linetype = "dashed")
  RCMplot$plot = RCMplot$plot + geom_segment(inherit.aes = FALSE, mapping = aes_string(x="x", y= "y", xend = "xend", yend = "yend"), data = dfStart, linetype = "dashed")

  # Add a red line for the projection
  dfRed = data.frame(xend = IntCoordsXTip, yend = IntCoordsYTip, x = IntCoordsXStart, y = IntCoordsYStart)
  RCMplot$plot  =  RCMplot$plot + geom_segment(inherit.aes = FALSE, col = "orange", mapping = aes_string(x="x", y= "y", xend = "xend", yend = "yend"), data = dfRed)

  if(addLabel){
  #Add some annotation
    labPos = if(is.null(labPos)) apply(RCMplot$samples[, dimNames],2,min)*1.1 else labPos
    xLab = labPos[1]
    yLab = labPos[2]
    dfRed = within(dfRed, {
      xLab = xLab*2
      yLab = yLab*2
    })
  RCMplot$plot  =  RCMplot$plot + geom_segment(inherit.aes = FALSE, mapping = aes_string(x= "xLab", y= "yLab", xend = "xend", yend = "yend"), data = dfRed/2, arrow = arrow(length = unit(0.2, "cm")), size = 0.25) + annotate("text", col = "orange",label = "r~psi~s", x = xLab, y = yLab, parse = TRUE, size = 7)
  }

  RCMplot$plot
}
