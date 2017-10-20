#' A function to plot RC(M) ordination result with the help of ggplot2
#'
#' @param RCMfit an RCM object
#' @param Dim A vector two dimensions to plot, defaults to the first two
#' @param samColour a character string, the variable to use for the colour of the sample dots
#' @param colLegend a character string, the legend text for the sample colour. Defaults to the name of the colour variable
#' @param samShape a character string, the variable to use for the shape of the sample dots
#' @param shapeLegend a character string, the text to use for the shapeLegend. Defaults to the name of the shape variable
#' @param samSize a scalar, the size of the sample dots, defaults to 3
#' @param taxNum an integer, the number of taxa to be plotted
#' @param scalingFactor a scalar, a user supplied scaling factor for the taxon location. If not supplied it will be calculated to make sample and taxon plots on the same scale
#' @param plotType a character string: which components to plots
#' @param quadDrop a scalar, the value at which the response function has dropped to this fraction of its maximum on the count scale to be plotted as a halo
#' @param nPoints an integer, number of points to evaluate for the ellipse
#' @param plotEllipse a boolean, whether to add the ellipses
#' @param taxaScale a scalar, by which to scale the rectangles of the quadratic taxon plot
#' @param Palette the colour palette
#' @param taxLabels a boolean, shoudl taxon labels be plotted
#' @param taxCol the taxon colour
#' @param arrowCol a character string: the colour of the arrows of the environmental gradient.
#' @param nudge_y a scalar, the offet for the taxon labels
#' @param square A boolean, should the plot be square? This is highly preferred to honestly represent differences
#' @param xInd a scalar or a vector of length 2, specifying the indentation left and right of the plot to allow for the labels to be printed entirely. Defaults to 0.75 at every side
#' @param yInd a scalar or a vector of length 2, specifying the indentation top and bottom of the plot to allow for the labels to be printed entirely. Defaults to 0 at every side
#' @param labSize the size of the variable labels
#' @param taxRegExp a character vector indicating which taxa to plot
#' @param varNum an integer, number of variable arrows to draw
#' @param alpha a boolean, should small arrows be made transparent?
#'
#' @return see the ggplot()-function
plot.RCM = function(RCMfit, Dim = c(1,2),
                    samColour = NULL, colLegend = samColour, samShape = NULL, shapeLegend = samShape, samSize = 1.5,
                    taxNum = if(all(plotType=="species") || !is.null(taxRegExp)) {ncol(RCMfit$X)} else {10}, scalingFactor = NULL, plotType = c("samples","species","variables"), quadDrop = 0.995, nPoints = 1e3, plotEllipse = TRUE, taxaScale = 0.5,
                    Palette = NULL, taxLabels = !all(plotType=="species"), taxCol = "blue", arrowCol = "blue", nudge_y = -0.08, square = TRUE, xInd = c(-0.75,0.75), yInd = c(0,0), labSize = 2.5, taxRegExp = NULL, varNum = 5, alpha = TRUE,...) {
  #Retrieve dots (will be passed on to aes())
  dotList = list(...)
  constrained = !is.null(RCMfit$covariates)
  if(constrained){
    dataSam <- data.frame(RCMfit$covariates %*% RCMfit$alpha[,Dim] %*% diag(RCMfit$psis[Dim]))
  } else {
    dataSam <- data.frame(RCMfit$rMat[, Dim] %*% diag(RCMfit$psis[Dim]))
    plotType = plotType[plotType!="variables"]
  }
  names(dataSam)=paste0("Dim", Dim)

  #Get the colours
  if(length(samColour)==1){
    dataSam$colourPlot = get_variable(RCMfit$physeq, samColour)
    if(is.character(dataSam$colourPlot)) dataSam$colourPlot  = factor(dataSam$colourPlot )
  } else if(!is.null(samColour)){
    dataSam$colourPlot = samColour
    } else {dataSam$colourPlot=factor(rep(1, nrow(dataSam)))}

  #     #Set colour palette
  if(is.null(Palette)){
    Palette = rainbow(length(unique(dataSam$colourPlot)))
  }

  idTaxRegExp = if(!is.null(taxRegExp)){  #Filter out certain taxa
    apply(sapply(taxRegExp, grepl, ignore.case = TRUE, x = colnames(RCMfit$X)),1,any) #Display only required taxa
  } else {rep(TRUE, ncol(RCMfit$X))}
  if(!any(idTaxRegExp)) {stop("Species not found! \n Check the dimnames of your RCMfit$X slot! \n")}
  taxFrac = min(taxNum/sum(idTaxRegExp),1)

  #Construct dataframe for taxa
  if(constrained) {
    if(RCMfit$responseFun=="linear"){
      dataTax = data.frame(
        origin1 = -RCMfit$NB_params[1,,Dim[1]]/RCMfit$NB_params[2,,Dim[1]] ,
        origin2 = -RCMfit$NB_params[1,,Dim[2]]/RCMfit$NB_params[2,,Dim[2]] ,
        slope1 = RCMfit$NB_params[2,,Dim[1]] * RCMfit$psis[Dim[1]], #The gradient
        slope2 = RCMfit$NB_params[2,,Dim[2]] * RCMfit$psis[Dim[2]]
      )
      rownames(dataTax) = colnames(RCMfit$X)
      dataTax = dataTax[idTaxRegExp,] #Keep only selected taxa
            dataTax = within(dataTax, {
        end1 = origin1 + slope1
        end2 = origin2 + slope2
      })
      dataTax$arrowLength = apply(dataTax[, c("slope1","slope2")],1,function(x){sqrt(sum(x^2))})
      id = dataTax$arrowLength >= quantile(dataTax$arrowLength,1-taxFrac)
      #Filter out small arrows
      dataTax = dataTax[id,]
      if(is.null(scalingFactor)){
        scalingFactor = min(abs(apply(dataSam[, paste0("Dim", Dim)],2,range)))/max(abs(dataTax[, c("slope1","slope2")]))*0.925
      }
      dataTax = within(dataTax, { #Scale the arrows
        end1 = origin1 + slope1 * scalingFactor
        end2 = origin2 + slope2 * scalingFactor
      })
    } else if (RCMfit$responseFun == "quadratic"){
      dataTax = data.frame(
        apply(RCMfit$NB_params[c(2,3),,Dim],c(2,3), function(x){
          a = x[2]; b=x[1]
          -b/(2*a)
        })) #The location of the extrema
      names(dataTax) = c("end1","end2")
      dataTax$colour = apply(RCMfit$NB_params[3,,Dim],1, function(x){
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
        x = RCMfit$NB_params[,tax,Dim]
        cbind(ellipseCoord(a = x[3,] * RCMfit$psis[Dim], b = x[2,] * RCMfit$psis[Dim], c = x[3,] * RCMfit$psis[Dim], quadDrop = quadDrop, nPoints = nPoints), taxon = tax)
      })

      peakHeights = apply(RCMfit$NB_params[,,Dim],2, function(x){
        A = x[3,]; B = x[2,]; C = x[1,];
        sapply(exp(B^2 * RCMfit$psis[Dim]-4*A*C)/(4*A), function(y){max(y,1/y)}) #select largest relative departure
      })
      rownames(peakHeights) = c("peak1","peak2")
      meanPeakHeights = colMeans(peakHeights)
      dataTax = cbind(dataTax, t(peakHeights))
      #Pick taxa with largest extrema, within observed values of the envrionmental scores (otherwise it is almost extrapolation)
      dataID = data.frame(
        meanPeakHeights = meanPeakHeights,
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
      #Rescale peak heights for plotting
      dataTax[, c("peak1","peak2")] = taxaScale* rowMultiply(dataTax[, c("peak1","peak2")], apply(dataTax[, c("end1","end2")],2,function(x){max(abs(x))})/apply(dataTax[, c("peak1","peak2")], 2,max))
      dataTax[, c("peak1","peak2")] = apply(dataTax[, c("peak1","peak2")], c(1,2), max, 0.0075) # Make sure a line always appears
      #Unfold into two dimensions
      dataEllipse = data.frame(
        apply(dataEllipseTmp[,,id], 2, c), colour = as.character(dataTax$colour)
      )
    } else if(RCMfit$responseFun  == "nonparametric"){ #For non-parametric response function we cannot plot the taxa
      plotType = plotType[plotType!="species"]
    } else {
      stop("No valid response function present in this RCM object!")
    }
  } else { #If not constrained
    dataTax = data.frame(cbind(t(RCMfit$cMat[Dim,]),0,0))
    names(dataTax) = c("end1","end2", "origin1","origin2")
    rownames(dataTax) = colnames(RCMfit$X)
    dataTax = dataTax[idTaxRegExp,] #Keep only selected taxa
  }
  if(!constrained){ #Scaling needed
    dataTax$arrowLength = apply(dataTax[, c("end1","end2")],1,function(x){sqrt(sum(x^2))})
    id = dataTax$arrowLength >= quantile(dataTax$arrowLength,1 - taxFrac)
    #Filter out small arrows
    dataTax = dataTax[id,]
    if(is.null(scalingFactor)){
      scalingFactor = min(abs(apply(dataSam[, paste0("Dim", Dim)],2,range)))/max(abs(dataTax[, c("end1","end2")]))*0.925
    }
    dataTax = within(dataTax, { #Scale the arrows
      end1 = end1 * scalingFactor
      end2 = end2 * scalingFactor
    })
  } # End scaling needed
  if("species" %in% plotType) dataTax$labels = sub(" ", "\n", rownames(dataTax))

  if("samples" %in% plotType){
  plot = ggplot(dataSam, aes_string(x=names(dataSam)[1], y=names(dataSam)[2], col = "colourPlot", dotList)) +
    geom_point(size = samSize) + #point size
    xlab(paste0(names(dataSam)[1],": ", paste0("psi",Dim[1]), " = ",round(RCMfit$psis[Dim[1]],1))) + #xlabel
    ylab(paste0(names(dataSam)[2],": ", paste0("psi",Dim[2]), " = ",round(RCMfit$psis[Dim[2]],1))) + #ylabel
    if(is.null(samColour)) {guides(color=FALSE)} #Legend

  #add legend names
  if(!is.null(colLegend) & is.factor(dataSam$colourPlot) ){
    plot = plot + scale_colour_manual(name = colLegend, values = Palette)
  }    else if(!is.null(colLegend) & !is.factor(dataSam$colourPlot) ){
    plot = plot + scale_colour_continuous(name = colLegend)
  }
  if(!is.null(shapeLegend)){
    plot = plot + scale_shape_discrete(name = shapeLegend)
  }

  } else {plot = ggplot()}
  if("species" %in% plotType){
    #Add arrows or labels
    if(length(taxCol)>1 && length(unique(taxCol))<10){
      taxCol = Palette[c(taxCol[id])]
    }
    if((!constrained || RCMfit$responseFun=="linear")){
      plot <- plot + geom_segment(data=dataTax, aes_string(x='origin1', y='origin2', xend="end1", yend = "end2", alpha = "arrowLength"), colour = taxCol, arrow=arrow(length=unit(0.1,"cm")), show.legend=FALSE, inherit.aes = FALSE) + if(alpha) scale_alpha_continuous(range = c(0.25,1))
    } else if(RCMfit$responseFun=="quadratic"){ #quadratic response functions
      plot <- plot +
        geom_tile(data=dataTax, aes_string(x='end1', y='end2', fill="colour", width = "peak1", height = "peak2" ), pch = 21, show.legend=FALSE, inherit.aes = FALSE) + #The centers
        if(plotEllipse) {geom_path(inherit.aes = FALSE, data = dataEllipse, mapping = aes_string(x = "x", y = "y", group = "taxon"), colour = "grey50", show.legend = FALSE)}
    } else {
      plot <- plot + geom_point(data=dataTax, aes_string(x='end1', y='end2', fill="taxCol"), pch = 21, show.legend = length(taxCol)!=1, inherit.aes = FALSE)
      if(!is.null(colLegend) & is.factor(dataTax$taxCol) ){
        plot = plot + scale_fill_manual(name = colLegend, values = Palette)
      }    else if(!is.null(colLegend) & !is.factor(dataTax$taxCol) ){
        plot = plot + scale_fill_continuous(name = colLegend)
      }
    }
if(taxLabels){
    plot <- plot + geom_text(data=dataTax, aes_string(x="end1", y = "end2", label = "labels"),  alpha=0.75, color=taxCol, show.legend=FALSE, nudge_y = nudge_y, size = labSize)
}
  } else {}
  if("variables" %in% plotType){
    #Add variable labels
    arrowLenghtsVar = rowSums(RCMfit$alpha[,Dim]^2) #All arrow lenghts
    attribs = attr(RCMfit$covariates, "assign")
    arrowLenghtsPerVar = tapply(arrowLenghtsVar, attribs, max) #Maximum per variable
    #Fix me! Number of variables
    # CumSum = cumsum(table(attribs)[unique(attribs)[order(arrowLenghtsPerVar, decreasing = TRUE)]]) <= varNum
    # varID = colnames(RCMfit$covariates)[as.numeric(names(CumSum)[CumSum])]
    idVar = arrowLenghtsPerVar >= quantile(arrowLenghtsPerVar, 1 - varNum/length(unique(attribs)))
    varData = data.frame(RCMfit$alpha)
    varData$label = rownames(RCMfit$alpha)
    varID = attribs %in% unique(attribs)[idVar]
 #Include all levels from important factors, not just the long arrows
    varData = varData[varID,]
    scalingFactorAlpha = min(abs(apply(dataSam[, paste0("Dim", Dim)],2, range)))/max(abs(varData[, paste0("Dim", Dim)]))*0.99
    varData[, paste0("Dim", Dim)] = varData[, paste0("Dim", Dim)]*scalingFactorAlpha
    plot = plot + geom_text(data = varData, mapping = aes_string(x = names(dataSam)[1], y = names(dataSam)[2], label = "label"), inherit.aes = FALSE, size = labSize)
  }
  #Add cross in the centre
  plot = plot + geom_point(data=data.frame(x=0,y=0), aes(x=x,y=y), size=5, inherit.aes = FALSE, shape=3)
  #Expand limits to show all text
  if(square) squarePlot(plot, xInd = xInd, yInd = yInd) else plot
}
