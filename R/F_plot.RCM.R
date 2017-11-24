#' A function to plot RC(M) ordination result with the help of ggplot2
#'
#' @param RCMfit an RCM object
#' @param Dim A vector two dimensions to plot, defaults to the first two
#' @param samColour a character string, the variable to use for the colour of the sample dots
#' @param colLegend a character string, the legend text for the sample colour. Defaults to the name of the colour variable
#' @param samShape a character string, the variable to use for the shape of the sample dots
#' @param shapeLegend a character string, the text to use for the shapeLegend. Defaults to the name of the shape variable
#' @param samSize a scalar, the size of the sample dots
#' @param taxNum an integer, the number of taxa to be plotted
#' @param scalingFactor a scalar, a user supplied scaling factor for the taxon arrows If not supplied it will be calculated to make sample and taxon plots on the same scale
#' @param plotType a character string: which components to plots
#' @param nPoints an integer, number of points to evaluate for the ellipse
#' @param plotEllipse a boolean, whether to add the ellipses
#' @param taxaScale a scalar, by which to scale the rectangles of the quadratic taxon plot
#' @param Palette the colour palette
#' @param taxLabels a boolean, shoudl taxon labels be plotted
#' @param taxDots a boolean, should taxa be plotted as dots?
#' @param taxCol the taxon colour
#' @param taxColSingle the taxon colour if there is only one
#' @param nudge_y a scalar, the offet for the taxon labels
#' @param square A boolean, should the plot be square? This is highly preferred to honestly represent differences
#' @param xInd a scalar or a vector of length 2, specifying the indentation left and right of the plot to allow for the labels to be printed entirely. Defaults to 0.75 at every side
#' @param yInd a scalar or a vector of length 2, specifying the indentation top and bottom of the plot to allow for the labels to be printed entirely. Defaults to 0 at every side
#' @param labSize the size of the variable labels
#' @param taxRegExp a character vector indicating which taxa to plot
#' @param varNum an integer, number of variable arrows to draw
#' @param alpha a boolean, should small arrows be made transparent?
#' @param alphaRange The range of transparency
#' @param arrowSize a scalar, the size of the arrows
#' @param Influence a boolean, should the influence of the observation on the variable be plotted
#' @param inflDim an integer, the dimension for which the influence should be calculated
#' @param richSupported A character vector of supported richness measures
#' @param returnCoords a boolea, should final coordinates be returned?
#' @param varExpFactor a scalar, the factor by which to expand the variable coordinates
#' @param manExpFactorTaxa a manual expansion factor for the taxa Setting it to a high value allows you to plot the taxa around the samples
#' @param nPhyl an integer, number of phylogenetic levels to show
#' @param phylOther a character vector of phylogenetic levels to be included in the "other" group
#' @param legendSize a size for the coloured dots in the legend
#' @param noLegend a boolean indicating you do not want a legend
#' @param crossSize the size of the central cross
#' @param contCol  a character vector of length two, giving the low and high values of the continuous colour scale
#'
#' @return see the ggplot()-function
plot.RCM = function(RCMfit, Dim = c(1,2),
                    samColour = NULL, colLegend = if(Influence) paste("Influence on\n", samColour, "\n parameter \n in dimension",inflDim) else samColour, samShape = NULL, shapeLegend = samShape, samSize = 1.5,
                    taxNum = if(all(plotType=="species") || !is.null(taxRegExp)) {ncol(RCMfit$X)} else {10}, scalingFactor = NULL, plotType = c("samples","species","variables"), quadDrop = 0.995, nPoints = 1e3, plotEllipse = TRUE, taxaScale = 0.5,
                    Palette = if(!all(plotType=="species")) "Set1" else "Paired", taxLabels = !all(plotType=="species"), taxDots = FALSE, taxCol = "blue", taxColSingle = "blue", nudge_y = -0.08, square = TRUE, xInd = if(all(plotType=="samples")) c(0,0) else c(-0.75, 0.75), yInd = c(0,0), labSize = 2, taxRegExp = NULL, varNum = 15, alpha = TRUE, alphaRange = c(0.2,1), arrowSize = 0.25, Influence = FALSE, inflDim = 1, richSupported = c("Observed", "Chao1", "ACE", "Shannon", "Simpson","InvSimpson", "Fisher"), returnCoords = FALSE, varExpFactor = 10, manExpFactorTaxa = 0.975, nPhyl = 10, phylOther = c(""), legendSize = samSize, noLegend = is.null(samColour), crossSize = 4, contCol = c("orange","darkgreen"),...) {
  require(RColorBrewer)
  #Retrieve dots (will be passed on to aes())
  dotList = list(...)
  constrained = !is.null(RCMfit$covariates) #Constrained plot?
  #Extract the coordinates
  coords = extractCoord(RCMfit, Dim)
  if(constrained && RCMfit$responseFun  == "nonparametric") plotType = plotType[plotType!="species"] #For non-parametric response function we cannot plot the taxa

## SAMPLES
  if("samples" %in% plotType){
  # taxCol = NULL #No species colours if samples are plotted
   dataSam = coords$samples
  #Get the sample colours
  if(length(samColour)==1){
    dataSam$colourPlot = if(Influence) rowSums(NBalphaInfl(RCMfit, inflDim)[,,samColour]) else if(samColour == "Deviance") rowSums(getDevianceRes(RCMfit, Dim)^2) else if(samColour %in% richSupported) estimate_richness(RCMfit$physeq, measures = samColour) else get_variable(RCMfit$physeq, samColour)
  } else if(!is.null(samColour)){
    dataSam$colourPlot = samColour
  } else {dataSam$colourPlot = factor(rep(1, nrow(dataSam)))}
   if(is.character(dataSam$colourPlot)) dataSam$colourPlot = factor(dataSam$colourPlot)
  #Get the sample shapes
  if(length(samShape)==1){
    dataSam$shapePlot = get_variable(RCMfit$physeq, samShape)
    if(is.character(dataSam$shapePlot)) dataSam$shapePlot  = factor(dataSam$shapePlot )
  } else if(!is.null(samShape)){
    dataSam$shapePlot = samShape
  } else {dataSam$shapePlot=factor(rep(1, nrow(dataSam)))}
   if(is.character(dataSam$shapePlot)) dataSam$shapePlot = factor(dataSam$shapePlot)

   #Set colour palette
   if(is.null(Palette)){
     Palette = rainbow(length(unique(dataSam$colourPlot)))
   }

   plot = ggplot(dataSam, aes_string(x=names(dataSam)[1], y=names(dataSam)[2], dotList, col = "colourPlot", shape = "shapePlot")) +
     geom_point(size = samSize ) + #point size
     xlab(paste0(names(dataSam)[1],": ", paste0("psi",Dim[1]), " = ",round(RCMfit$psis[Dim[1]],1))) + #xlabel
     ylab(paste0(names(dataSam)[2],": ", paste0("psi",Dim[2]), " = ",round(RCMfit$psis[Dim[2]],1))) + #ylabel
     if(noLegend) {guides(colour=FALSE)}  #Legend

   #add legend names
   if(!is.null(colLegend) & is.factor(dataSam$colourPlot) ){
     plot = plot + scale_colour_manual(name = colLegend, values = colorRampPalette(brewer.pal(max(3,length(unique(dataSam$colourPlot))), Palette))(length(unique(dataSam$colourPlot))))
   }    else if(!is.null(colLegend) & !is.factor(dataSam$colourPlot) ){
     plot = plot + scale_colour_continuous(name = colLegend, low = contCol[1], high =  contCol[2])
   }
   if(!is.null(shapeLegend)){
     plot = plot + scale_shape_discrete(name = shapeLegend)
   } else {plot = plot + guides(shape=FALSE)}
} else {dataSam=NULL;plot = ggplot()}# END if samples %in% plotType

  ## VARIABLES
  if("variables" %in% plotType && constrained){
    #Add variable labels
    arrowLenghtsVar = rowSums(RCMfit$alpha[,Dim]^2) #All arrow lenghts
    attribs = attr(RCMfit$covariates, "assign")
    arrowLenghtsPerVar = tapply(arrowLenghtsVar, attribs, max) #Maximum per variable
    CumSum = cumsum(table(attribs)[unique(attribs)[order(arrowLenghtsPerVar, decreasing = TRUE)]]) <= varNum
    varID = attr(RCMfit$covariates, "dimnames")[[2]][attribs %in% as.numeric(names(CumSum)[CumSum])]
    varData = data.frame(RCMfit$alpha * if(!all(plotType=="variables")) 1 else varExpFactor)
    varData$label = rownames(RCMfit$alpha)
    # varID = attribs %in% unique(attribs)[idVar]
    #Include all levels from important factors, not just the long arrows
    varData = varData[varID,]
    # scalingFactorAlpha = min(abs(apply(dataSam[, paste0("Dim", Dim)],2, range)))/max(abs(varData[, paste0("Dim", Dim)]))*0.99
    if("samples" %in% plotType){
      scalingFactorAlphaTmp = apply(dataSam[, paste0("Dim", Dim)],2,range)/apply(varData[, paste0("Dim", Dim)],2,range)
      scalingFactorAlpha = min(scalingFactorAlphaTmp[scalingFactorAlphaTmp>0])*0.975
      varData[, paste0("Dim", Dim)] = varData[, paste0("Dim", Dim)]*scalingFactorAlpha
    }
    plot = plot + geom_text(data = varData, mapping = aes_string(x = names(varData)[1], y = names(varData)[2], label = "label"), inherit.aes = FALSE, size = labSize)
  } else {varData = NULL}

  ## TAXA
  if("species" %in% plotType){
  idTaxRegExp = if(!is.null(taxRegExp)){  #Filter out certain taxa
    apply(sapply(taxRegExp, grepl, ignore.case = TRUE, x = rownames(coords$species)),1,any) #Display only required taxa
  } else {rep(TRUE, ncol(RCMfit$X))}
  if(!any(idTaxRegExp)) {stop("Species not found! \n Check the dimnames of your RCMfit$X slot! \n")}
  taxFrac = min(taxNum/sum(idTaxRegExp),1)
  dataTax = coords$species[idTaxRegExp,] #Keep only selected taxa
  #Construct dataframe for taxa
    if(constrained) {
    if(RCMfit$responseFun=="linear"){
      dataTax$arrowLength = apply(dataTax[, c("slope1","slope2")],1,function(x){sqrt(sum(x^2))})
      id = dataTax$arrowLength >= quantile(dataTax$arrowLength,1-taxFrac)
      #Filter out small arrows
      dataTax = dataTax[id,]
      if(!all(plotType=="species")){
      scalingFactorTmp = apply(if("samples" %in% plotType) {dataSam[, paste0("Dim", Dim)]} else {varData[, paste0("Dim", Dim)]},2,range)/apply(dataTax[, c("end1","end2")]-dataTax[, c("origin1","origin2")],2,range)
      scalingFactor = min(scalingFactorTmp[scalingFactorTmp>0])*0.975
      dataTax = within(dataTax, { #Scale the arrows
        end1 = origin1 + slope1 * scalingFactor
        end2 = origin2 + slope2 * scalingFactor
      })
      }
    } else if (RCMfit$responseFun == "quadratic"){
      dataTax$colour = apply(coords$species[, paste0("a",Dim)],1, function(x){
        if(all(x>0)) {
          return("green")
        } else if (all(x<0)){
          return("red")
        } else if(x[1] > 0){
          return("brown")
        } else{
          return("purple")
        }
      })
      dataEllipseTmp = vapply(seq_along(taxa_names(RCMfit$physeq)), FUN.VALUE = matrix(0,nPoints,3),function(tax){
        x = coords$species[tax,]
        cbind(ellipseCoord(a = unlist(x[paste0("a", Dim)]) * RCMfit$psis[Dim], b = unlist(x[paste0("b", Dim)]) * RCMfit$psis[Dim], c = unlist(x[paste0("a", Dim)]) * RCMfit$psis[Dim], quadDrop = quadDrop, nPoints = nPoints), taxon = tax)
      })
      #Pick taxa with largest extrema, within observed values of the envrionmental scores (otherwise it is almost extrapolation)
      dataID = data.frame(
        meanPeakHeights = rowMeans(dataTax[,paste0("peak",Dim)]),
        id = seq_len(nrow(dataTax)),
        dataTax
      )
      envScores = RCMfit$covariates %*% RCMfit$alpha
      rownames(dataTax) = colnames(RCMfit$X)
      dataTax = dataTax[idTaxRegExp,] #Keep only selected taxa
      id = dataID[with(dataID,
                       order(end1 > max(envScores[,Dim[1]]) | end1 < min(envScores[,Dim[1]]),
                             end2 > max(envScores[,Dim[2]]) | end2 < min(envScores[,Dim[2]]),
                             -meanPeakHeights)),
                  ]$id[seq_len(ceiling(taxFrac*nrow(dataTax)))]
      dataTax = dataTax[id,]
      dataTax[, c("peak1","peak2")] = taxaScale *apply(dataTax[, c("peak1","peak2")], c(1,2), max, 0.0075) # Make sure a line always appears
      #Unfold into two dimensions
      dataEllipse = data.frame(
        apply(dataEllipseTmp[,,id], 2, c), colour = as.character(dataTax$colour)
      )
    } else {
      stop("No valid response function present in this RCM object!")
    }
  }else {
    dataTax$arrowLength = apply(dataTax[, c("end1","end2")],1,function(x){sqrt(sum(x^2))})
    id = dataTax$arrowLength >= quantile(dataTax$arrowLength, 1 - taxFrac)
    #Filter out small arrows
    dataTax = dataTax[id,]
    if("samples" %in% plotType){
    scalingFactorTmp = apply(dataSam[, paste0("Dim", Dim)],2,range)/apply(dataTax[, c("end1","end2")],2,range)
    scalingFactor = min(scalingFactorTmp[scalingFactorTmp>0])*manExpFactorTaxa
    #The scaling factor is the minimum of the ratios between the longest longest arrow and the longest species arrow in every direction of every dimension
    dataTax = within(dataTax, { #Scale the arrows
      end1 = end1 * scalingFactor
      end2 = end2 * scalingFactor
    })
  } # End scaling needed
  }
  dataTax$labels = sub(" ", "\n", rownames(dataTax))
if(!"samples" %in% plotType && length(taxCol)==1) colLegend = taxCol
    #Add arrows or labels
    if(length(taxCol)>1 && length(unique(taxCol))<10){
      dataTax$taxCol = Palette[c(taxCol[id])]
    } else if(taxCol=="Deviance"){
      dataTax$taxCol = colSums(getDevianceRes(RCMfit, Dim)^2)[id]
    } else if(taxCol %in% colnames(tax_table(RCMfit$physeq, errorIfNULL = FALSE))){
      dataTax$taxCol = tax_table(RCMfit$physeq)[, taxCol]
      mostCommon = names(sort(table(dataTax$taxCol), decreasing = TRUE)[seq_len(nPhyl)])
      dataTax$taxCol[(!dataTax$taxCol %in% mostCommon) | (dataTax$taxCol %in% phylOther)] = "Other"
      dataTax$taxCol = factor(dataTax$taxCol)
    }
    if((!constrained || RCMfit$responseFun=="linear") ){
      if(arrowSize > 0){
      plot <- plot + geom_segment(data=dataTax, aes_string(x='origin1', y = 'origin2', xend="end1", yend = "end2", alpha = if(alpha) "arrowLength" else NULL, colour = if("samples" %in% plotType) NULL else  "taxCol"), colour = taxColSingle, arrow=arrow(length=unit(0.1,"cm")), inherit.aes = FALSE, size = arrowSize) +  guides(alpha = FALSE)
      if(!"species" %in% plotType){
      plot = plot + if(is.factor(taxCol)) scale_colour_discrete(name = colLegend) else scale_colour_continuous(name = colLegend, low = contCol[1], high = contCol[2])
      }
      plot = plot +  if(alpha) scale_alpha_continuous(range = alphaRange)
      }
    } else if(RCMfit$responseFun=="quadratic"){ #quadratic response functions
      plot <- plot +
        geom_tile(data=dataTax, aes_string(x='end1', y='end2', fill="colour", width = "peak1", height = "peak2" ), pch = 21, show.legend=FALSE, inherit.aes = FALSE) + #The centers
        if(plotEllipse) {geom_path(inherit.aes = FALSE, data = dataEllipse, mapping = aes_string(x = "x", y = "y", group = "taxon"), colour = "grey50", show.legend = FALSE)}
    } else {
      plot <- plot + geom_point(data=dataTax, aes_string(x='end1', y='end2', fill="taxCol"), pch = 21, show.legend = length(taxCol)!=1, inherit.aes = FALSE)
    }
    if(!is.null(colLegend) & is.factor(dataTax$taxCol) & !taxDots ){
      plot = plot + scale_colour_brewer(palette = Palette, name = colLegend)
    } else if(!is.null(colLegend) & !is.factor(dataTax$taxCol) ){
      plot = plot + scale_fill_continuous(name = colLegend)
    }
if(taxLabels){
    plot <- plot +  if(is.null(dataTax$taxCol)){geom_text(data=dataTax, aes_string(x="end1", y = "end2", label = "labels", color = "taxCol"), color =  taxColSingle, show.legend=FALSE, nudge_y = nudge_y, size = labSize, inherit.aes = FALSE)
    } else {
      geom_text(data=dataTax, aes_string(x="end1", y = "end2", label = "labels", color = "taxCol"), show.legend=TRUE, nudge_y = nudge_y, size = labSize, inherit.aes = FALSE)
    }
} else if(taxDots){
  if(is.null(dataTax$taxCol)){plot <- plot + geom_point(data=dataTax, aes_string(x = "end1", y = "end2", color = "taxCol"), color =  taxColSingle, show.legend = FALSE, nudge_y = nudge_y, size = labSize, inherit.aes = FALSE)
  } else {
    plot <- plot + geom_point(data = dataTax, aes_string(x="end1", y = "end2", color = "taxCol"), show.legend = TRUE, size = labSize, inherit.aes = FALSE) + scale_colour_manual(values = c(brewer.pal(length(unique(dataTax$taxCol))-1, Palette), "Grey90"), name = colLegend) # "Other" is made grey
  }
}
  if(!"samples" %in% plotType){
    plot = plot +
      xlab(paste("Dim",Dim[1])) + #xlabel
      ylab(paste("Dim",Dim[2]))
  }
  } #END if "species" %in% plotType

  #Add cross in the centre, and enlarge legend sizes
  plot = plot + geom_point(data=data.frame(x=0,y=0), aes(x=x,y=y), size = crossSize, inherit.aes = FALSE, shape=3) + guides(size=legendSize)
  #Expand limits to show all text
  plot = if(square) squarePlot(plot, xInd = xInd, yInd = yInd) else plot
  if(returnCoords){
  list(plot = plot, samples = dataSam, species = dataTax, variables  = varData)
  } else {plot}
}
