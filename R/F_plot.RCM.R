#' Plot RC(M) ordination result with the help of ggplot2
#'
#' @param x an RCM object
#' @param ... further arguments, passed on to aes in the the ggplot() function
#' @param Dim An integer vector of length two, which dimensions to plot
#' @param plotType a character string: which components to plot.
#'  Can be any combination of 'samples','species' and 'variables'
#' @param samColour a character string, the variable to use for the colour
#' of the sample dots. Alternatively,
#' a vector equal to the number of samples in the RCM object can be supplied
#' @param taxNum an integer, the number of taxa to be plotted
#' @param taxRegExp a character vector indicating which taxa to plot.
#' Any taxa matcing this regular expression will be plotted
#' @param varNum an integehr, number of variable arrows to draw
#' @param varPlot the names of the variable arrows to plot.
#'  Overrides the varNum argument
#' @param arrowSize a scalar, the size of the arrows
#' @param Influence a boolean, should the influence of the observation
#' on the variable be plotted
#' @param inflDim an integer, the dimension for which the influence
#' should be calculated
#' @param returnCoords a boolean, should final coordinates be returned?
#' @param alpha a boolean, should small arrows be made transparent?
#' @param colLegend a character string, the legend text for the sample colour.
#'  Defaults to the name of the colour variable
#' @param samShape a character string, the variable to use for the shape
#'  of the sample dots
#' @param shapeLegend a character string, the text to use for the shapeLegend.
#'  Defaults to the name of the shape variable
#' @param samSize a scalar, the size of the sample dots
#' @param scalingFactor a scalar, a user supplied scaling factor
#'  for the taxon arrows. If not supplied it will be calculated to make sample
#'   and taxon plots on the same scale
#' @param quadDrop a number between 0 and 1. At this fraction of the peak height
#'  are the ellipses of the quadratic response functions drawn
#' @param plotEllipse a boolean, whether to add the ellipses
#' @param taxaScale a scalar, by which to scale the rectangles
#'  of the quadratic taxon plot
#' @param Palette the colour palette
#' @param taxLabels a boolean, should taxon labels be plotted?
#' @param taxDots a boolean, should taxa be plotted as dots?
#' @param taxCol the taxon colour
#' @param taxColSingle the taxon colour if there is only one
#' @param nudge_y a scalar, the offet for the taxon labels
#' @param axesFixed A boolean, should the aspect ratio of the plot
#'  (the scale between the x and y-axis) be fixed.
#'  It is highly recommended to keep this argument at TRUE
#'  for honest representation of the ordination. If set to FALSE,
#'  the plotting space will be optimally used but the plot
#'   may be deformed in the process.
#' @param aspRatio The aspect ratio of the plot when 'axesfixed' is TRUE
#' (otherwise this argument is ignored), passde on to ggplot2::coord_fixed().
#' It is highly recommended to keep this argument at 1 for honest
#'  representation of the ordination.
#' @param xInd a scalar or a vector of length 2, specifying the indentation
#'  left and right of the plot to allow for the labels to be printed entirely.
#'   Defaults to 0.75 at every side
#' @param yInd a scalar or a vector of length 2, specifying the indentation
#'  top and bottom of the plot to allow for the labels to be printed entirely.
#'   Defaults to 0 at every side
#' @param taxLabSize the size of taxon labels
#' @param varLabSize the size of the variable label
#' @param alphaRange The range of transparency
#' @param varExpFactor a scalar, the factor by which to expand
#' the variable coordinates
#' @param manExpFactorTaxa a manual expansion factor for the taxa.
#'  Setting it to a high value allows you to plot the taxa around the samples
#' @param nPhyl an integer, number of phylogenetic levels to show
#' @param phylOther a character vector of phylogenetic levels
#'  to be included in the 'other' group
#' @param legendSize a size for the coloured dots in the legend
#' @param noLegend a boolean indicating you do not want a legend
#' @param crossSize the size of the central cross
#' @param contCol  a character vector of length two, giving the low
#' and high values of the continuous colour scale
#' @param legendLabSize size of the legend labels
#' @param legendTitleSize size of the legend title
#' @param axisLabSize size of the axis labels
#' @param axisTitleSize size of the axis title
#' @param plotPsi a character vector, describing what to plot on the axis.
#'  Can be either 'psi', 'none' or 'loglik'.
#'  The latter plots the log-likelihood explained
#' @param breakChar a character string indicating how the taxon names
#'  should be broken
#'
#' @details
#' This function relies on the ggplot2 machinery to produce the plots,
#' and the result can be modified accordingly. Monoplots,
#' biplots and for constrained analysis even triplots can be produced,
#' depending on the 'plotType' argument.
#'
#' When one of either 'Observed', 'Chao1', 'ACE', 'Shannon', 'Simpson',
#' 'InvSimpson' or 'Fisher' are supplied to the 'samColour' argument,
#' the according richness measure (as calculated by phyloseq::estimate_richness)
#'  is mapped to the sample colour
#'
#' @return plots a ggplot2-object to output
#' @export
#' @import ggplot2
#' @import phyloseq
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom graphics par text
#' @importFrom RColorBrewer brewer.pal
#' @method plot RCM
#' @seealso \code{\link{RCM}},\code{\link{addOrthProjection}},
#' \code{\link{extractCoord}},\code{\link{plotRespFun}}
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' # Subset for a quick fit
#' zellerRCM = RCM(tmpPhy)
#' plot(zellerRCM)
plot.RCM = function(x,
                    ...,
                    Dim = c(1, 2),
                    plotType = c("samples", "species",
                    "variables"),
                    samColour = NULL,
                    taxNum = if (all(plotType == "species") ||
                    !is.null(taxRegExp)) {
                    ncol(x$X)
                    } else {
                    10
                    },
                    taxRegExp = NULL,
                    varNum = 15,
                    arrowSize = 0.25,
                    Influence = FALSE,
                    inflDim = 1,
                    returnCoords = FALSE,
                    alpha = TRUE,
                    varPlot = NULL,
                    colLegend = if (Influence)
                    paste0("Influence on\n",
                    samColour,
                    "\nparameter \nin dimension",
                    inflDim)
                    else
                    samColour,
                    samShape = NULL,
                    shapeLegend = samShape,
                    samSize = 2,
                    scalingFactor = NULL,
                    quadDrop = 0.995,
                    plotEllipse = TRUE,
                    taxaScale = 0.5,
                    Palette = if (!all(plotType ==
                    "species"))
                    "Set1"
                    else
                    "Paired",
                    taxLabels = !all(plotType == "species"),
                    taxDots = FALSE,
                    taxCol = "blue",
                    taxColSingle = "blue",
                    nudge_y = 0.08,
                    axesFixed = TRUE,
                    aspRatio = 1,
                    xInd = if (all(plotType == "samples"))
                    c(0,0) else
                    c(-0.75, 0.75),
                    yInd = c(0, 0),
                    taxLabSize = 4,
                    varLabSize = 3.5,
                    alphaRange = c(0.2, 1),
                    varExpFactor = 10,
                    manExpFactorTaxa = 0.975,
                    nPhyl = 10,
                    phylOther = c(""),
                    legendSize = samSize,
                    noLegend = is.null(samColour),
                    crossSize = 4,
                    contCol = c("orange", "darkgreen"),
                    legendLabSize = 15,
                    legendTitleSize = 16,
                    axisLabSize = 14,
                    axisTitleSize = 16,
                    plotPsi = "psi",
                    breakChar = "\n"
) {
    # Retrieve dots (will be passed on to aes())
    dotList = list(...)
    richSupported = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
        "InvSimpson", "Fisher")
    constrained = !is.null(x$covariates)  #Constrained plot?
    # Extract the coordinates
    coords = extractCoord(x, Dim)
    Dimnames = paste0("Dim", Dim)  # A text form of the dimensions
    if (constrained && x$responseFun == "nonparametric") {
        plotType = "variables"
        # For non-parametric response function we, can only plot the variables
        # meaningfully
    }

    ## SAMPLES
    if ("samples" %in% plotType) {
        dataSam = coords$samples
        # Get the sample colours
        if (length(samColour) == 1) {
            dataSam$colourPlot = if (Influence) {
                rowSums(NBalphaInfl(x, inflDim)[, , samColour])
            } else if (samColour == "Deviance") {
                rowSums(getDevianceRes(x, max(Dim))^2)
            } else if (samColour %in% richSupported)
                estimate_richness(x$physeq, measures = samColour)[[1]] else
                get_variable(x$physeq, samColour)
        } else if (!is.null(samColour)) {
            dataSam$colourPlot = samColour
        } else {
            dataSam$colourPlot = factor(rep(1, nrow(dataSam)))
        }
        if (is.character(dataSam$colourPlot))
            dataSam$colourPlot = factor(dataSam$colourPlot)
        # Get the sample shapes
        if (length(samShape) == 1) {
            dataSam$shapePlot = get_variable(x$physeq, samShape)
            if (is.character(dataSam$shapePlot))
                dataSam$shapePlot = factor(dataSam$shapePlot)
        } else if (!is.null(samShape)) {
            dataSam$shapePlot = samShape
        } else {
            dataSam$shapePlot = factor(rep(1, nrow(dataSam)))
        }
        if (is.character(dataSam$shapePlot))
            dataSam$shapePlot = factor(dataSam$shapePlot)

        # Set colour palette
        if (is.null(Palette)) {
            Palette = rainbow(length(unique(dataSam$colourPlot)))
        }

        plot = ggplot(dataSam, aes_string(x = names(dataSam)[1],
        y = names(dataSam)[2],
        dotList, col = "colourPlot", shape = "shapePlot")) +
        geom_point(size = samSize) +
            if (noLegend)
                {
        guides(colour = FALSE)
                }  #Legend

        # add legend names
        if (!is.null(colLegend) & is.factor(dataSam$colourPlot)) {
            plot = plot +
            scale_colour_manual(name = colLegend,
            values = colorRampPalette(brewer.pal(max(3,
                length(unique(dataSam$colourPlot))), Palette))(length(unique(
                dataSam$colourPlot))))
        } else if (!is.null(colLegend) & !is.factor(dataSam$colourPlot)) {
            plot = plot + scale_colour_continuous(name = colLegend,
            low = contCol[1],
            high = contCol[2])
        }
        if (!is.null(shapeLegend)) {
            plot = plot + scale_shape_discrete(name = shapeLegend)
        } else {
            plot = plot + guides(shape = FALSE)
        }
    } else {
        dataSam = NULL
        plot = ggplot()
    }  # END if samples %in% plotType

    ## TAXA
    if ("species" %in% plotType)
        {
            idTaxRegExp = if (!is.null(taxRegExp)) {
                # Filter out certain taxa
                apply(vapply(FUN.VALUE = logical(nrow(coords$species)),
                taxRegExp, grepl, ignore.case = TRUE,
                x = rownames(coords$species)),
                1, any)
                # Display only required taxa
            } else {
                rep(TRUE, ncol(x$X))
            }
            if (!any(idTaxRegExp)) {
                stop("Species not found! \n Check the dimnames
                of your x$X slot! \n")
            }
            taxFrac = min(taxNum/sum(idTaxRegExp), 1)
            dataTax = coords$species[idTaxRegExp, ]  #Keep only selected taxa
            # Construct dataframe for taxa
            if (constrained) {
                if (x$responseFun == "linear") {
                dataTax$arrowLength = apply(dataTax[, c("slope1", "slope2")],
                1, function(x) {
                sqrt(sum(x^2))
                })
                id = dataTax$arrowLength >= quantile(dataTax$arrowLength,
                1 - taxFrac)
                # Filter out small arrows
                dataTax = dataTax[id, ]
                if ("samples" %in% plotType) {
                scalingFactorTmp = apply(dataSam[, Dimnames], 2, range)/
                apply(dataTax[,
                c("end1", "end2")] - dataTax[, c("origin1", "origin2")],
                2, range)
                scalingFactor = min(scalingFactorTmp[scalingFactorTmp >
                0]) * 0.975
                # Scale the arrows
                dataTax[, c("end1", "end2")] = dataTax[, c("origin1",
                "origin2")] + dataTax[, c("slope1", "slope2")] *
                scalingFactor}
                } else if (x$responseFun == "quadratic") {
                dataTax$colour = apply(coords$species[, paste0("a", Dim)],
                1, function(x) {
                if (all(x > 0)) {
                return("green")
                } else if (all(x < 0)) {
                return("red")
                } else if (x[1] > 0) {
                return("brown")
                } else {
                return("purple")
                }
                })
                dataEllipseTmp = vapply(seq_along(taxa_names(x$physeq)),
                FUN.VALUE = matrix(0, 1000L, 3), function(tax) {
                x = coords$species[tax, ]
                cbind(ellipseCoord(a = unlist(x[paste0("a", Dim)]) *
                    x$psis[Dim], b = unlist(x[paste0("b", Dim)]) *
                        x$psis[Dim],
                        c = unlist(x[paste0("a", Dim)]) * x$psis[Dim],
                        quadDrop = quadDrop,
                        nPoints = 1000L), taxon = tax)
                    })
                # Pick taxa with largest extrema,
                # within observed values of the
                # envrionmental scores (otherwise it is almost extrapolation)
                dataID = data.frame(meanPeakHeights = rowMeans(dataTax[,
                    paste0("peak", Dim)]), id = seq_len(nrow(dataTax)),
                    dataTax)
                envScores = x$covariates %*% x$alpha
                rownames(dataTax) = colnames(x$X)
                dataTax = dataTax[idTaxRegExp, ]  #Keep only selected taxa
                id = dataID[order(dataID$end1 > max(envScores[, Dim[1]]) |
                    dataID$end1 < min(envScores[, Dim[1]]), dataID$end2 >
                    max(envScores[, Dim[2]]) | dataID$end2 < min(envScores[,
                    Dim[2]]), -dataID$meanPeakHeights), ]$id[
                    seq_len(ceiling(taxFrac *
                    nrow(dataTax)))]
                dataTax = dataTax[id, ]
                dataTax[, c("peak1", "peak2")] = taxaScale * apply(dataTax[,
                    c("peak1", "peak2")], c(1, 2), max, 0.0075)
                # Make sure a line always appears Unfold into two dimensions
                dataEllipse = data.frame(apply(dataEllipseTmp[, , id],
                    2, c), colour = as.character(dataTax$colour))
                } else {
                stop("No valid response function present in this RCM object!")
                }
            } else {
                dataTax$arrowLength = apply(dataTax[, c("end1", "end2")],
                1, function(x) {
                    sqrt(sum(x^2))
                })
                id = dataTax$arrowLength >= quantile(dataTax$arrowLength,
                1 - taxFrac)
                # Filter out small arrows
                dataTax = dataTax[id, ]
                if ("samples" %in% plotType)
                {
                    scalingFactorTmp = apply(dataSam[, Dimnames], 2, range)/
                    apply(dataTax[,
                    c("end1", "end2")], 2, range)
                    scalingFactor = min(scalingFactorTmp[scalingFactorTmp >
                    0]) * manExpFactorTaxa
                    # The scaling factor is the minimum of the ratios between
                    # the longest arrow and the longest species arrow in every
                    # direction of every dimension
                    dataTax[, c("end1", "end2")] = dataTax[, c("end1", "end2")]*
                    scalingFactor
                }  # End scaling needed
            }
            dataTax$labels = sub(" ", breakChar, rownames(dataTax))
            if (!"samples" %in% plotType && length(taxCol) == 1)
                colLegend = taxCol
            # Add arrows or labels
            if (length(taxCol) > 1 && length(unique(taxCol)) < 10) {
                dataTax$taxCol = Palette[c(taxCol[id])]
            } else if (taxCol == "Deviance") {
                dataTax$taxCol = colSums(getDevianceRes(x, max(Dim))^2)[id]
            } else if (taxCol %in% colnames(tax_table(x$physeq,
            errorIfNULL = FALSE))) {
                dataTax$taxCol = tax_table(x$physeq)[, taxCol]
                mostCommon = names(sort(table(dataTax$taxCol),
                decreasing = TRUE)[seq_len(nPhyl)])
                dataTax$taxCol[(!dataTax$taxCol %in% mostCommon) |
                (dataTax$taxCol %in%
                phylOther)] = "Other"
                dataTax$taxCol = factor(dataTax$taxCol)
            }
            if ((!constrained || x$responseFun == "linear")) {
                if (arrowSize > 0) {
                    if ("samples" %in% plotType | (length(taxCol) == 1 &&
                    taxCol != "Deviance")) {
                    plot <- plot + geom_segment(data = dataTax,
                            aes_string(x = "origin1",
                    y = "origin2", xend = "end1", yend = "end2",
                    alpha = if (alpha)
                        "arrowLength" else NULL), colour = taxColSingle,
                    arrow = arrow(length = unit(0.1,
                    "cm")), inherit.aes = FALSE, size = arrowSize) +
                    guides(alpha = FALSE)
                    } else {
                    plot <- plot + geom_segment(data = dataTax,
                    aes_string(x = "origin1",
                    y = "origin2", xend = "end1", yend = "end2",
                    alpha = if (alpha)
                    "arrowLength" else NULL, colour = "taxCol"),
                    arrow = arrow(length = unit(0.1,
                    "cm")), inherit.aes = FALSE, size = arrowSize) +
                    guides(alpha = FALSE)
                    }
                    if (!("samples" %in% plotType | (length(taxCol) == 1 &&
                    taxCol != "Deviance"))) {
                    plot = plot + if (is.factor(taxCol))
                        scale_colour_discrete(name = colLegend) else
                        scale_colour_continuous(name = colLegend,
                        low = contCol[1],
                        high = contCol[2])
                    }
                    plot = plot + if (alpha)
                    scale_alpha_continuous(range = alphaRange)
                }
            } else if (x$responseFun == "quadratic") {
                # quadratic response functions
                plot <- plot + geom_tile(data = dataTax, aes_string(x = "end1",
                    y = "end2", fill = "colour", width = "peak1",
                    height = "peak2"),
                    pch = 21, show.legend = FALSE, inherit.aes = FALSE) +
                    if (plotEllipse) {
                    geom_path(inherit.aes = FALSE, data = dataEllipse,
                            mapping = aes_string(x = "x",
                    y = "y", group = "taxon"), colour = "grey50",
                    show.legend = FALSE)
                    }
            } else {
                plot <- plot + geom_point(data = dataTax, aes_string(x = "end1",
                    y = "end2", fill = "taxCol"), pch = 21,
                    show.legend = length(taxCol) !=
                    1, inherit.aes = FALSE)
            }
            if (!is.null(colLegend) & is.factor(dataTax$taxCol) & !taxDots) {
                plot = plot + scale_colour_brewer(palette = Palette,
                name = colLegend)
            } else if (!is.null(colLegend) & !is.factor(dataTax$taxCol)) {
                plot = plot + scale_fill_continuous(name = colLegend)
            }
            if (taxLabels) {
                dataTax$end2b = dataTax$end2 + nudge_y * ifelse(dataTax$end2 >
                0, 1, -1)
                plot <- plot + if (is.null(dataTax$taxCol)) {
                geom_text(data = dataTax, aes_string(x = "end1", y = "end2b",
                    label = "labels", color = "taxCol"), color = taxColSingle,
                    show.legend = FALSE, size = taxLabSize, inherit.aes = FALSE)
                } else {
                geom_text(data = dataTax, aes_string(x = "end1", y = "end2",
                    label = "labels", color = "taxCol"), show.legend = TRUE,
                    nudge_y = nudge_y, size = taxLabSize, inherit.aes = FALSE)
                }
            } else if (taxDots) {
                if (is.null(dataTax$taxCol)) {
                plot <- plot + geom_point(data = dataTax,
                aes_string(x = "end1",
                    y = "end2", color = "taxCol"), color = taxColSingle,
                    show.legend = FALSE, nudge_y = nudge_y, size = taxLabSize,
                    inherit.aes = FALSE)
                } else {
                plot <- plot + geom_point(data = dataTax,
                aes_string(x = "end1",
                    y = "end2", color = "taxCol"), show.legend = TRUE,
                    size = taxLabSize,
                    inherit.aes = FALSE) + if (!is.numeric(dataTax$taxCol))
                    scale_colour_manual(values = c(brewer.pal(length(
                    unique(dataTax$taxCol)) -
                    1, Palette), "Grey90"), name = colLegend) else
                        scale_colour_continuous(name = colLegend)
                # 'Other' is made grey
                }
            }
            if (!"samples" %in% plotType) {
                # xlabel
                plot = plot + xlab(Dimnames[1]) + ylab(Dimnames[2])
            }
        }  #END if 'species' %in% plotType

    ## VARIABLES
    if ("variables" %in% plotType && constrained) {
        # Add variable labels
        if (is.null(varPlot)) {
            arrowLenghtsVar = rowSums(x$alpha[, Dim]^2)  #All arrow lenghts
            attribs = x$attribs
            arrowLenghtsPerVar = tapply(arrowLenghtsVar, attribs, max)
            # Maximum per variable
            CumSum = cumsum(table(attribs)[unique(attribs)[
            order(arrowLenghtsPerVar,
                decreasing = TRUE)]]) <= varNum
            varID = attr(x$covariates, "dimnames")[[2]][attribs %in%
            as.numeric(names(CumSum)[CumSum])]
        } else {
            varID = attr(x$covariates, "dimnames")[[2]] %in%
            unlist(lapply(varPlot,
                grep, value = TRUE, x = attr(x$covariates, "dimnames")[[2]]))
        }
        varData = data.frame(x$alpha * if (!all(plotType == "variables"))
            1 else varExpFactor)
        varData$label = rownames(x$alpha)
        # Include all levels from important factors, not just the long arrows
        varData = varData[varID, ]
        if (!all(plotType == "variables")) {
            if ("samples" %in% plotType) {
                scalingFactorAlphaTmp = apply(dataSam[, Dimnames], 2, range)/
                apply(varData[,
                Dimnames], 2, range)
                scalingFactorAlpha = min(scalingFactorAlphaTmp[
                scalingFactorAlphaTmp >
                0]) * 0.975
            } else if ("species" %in% plotType) {
            scalingFactorAlphaTmp = apply(dataTax[, c("end1", "end2")],
            2, range)/apply(varData[, Dimnames], 2, range)
                scalingFactorAlpha = max(scalingFactorAlphaTmp[
                scalingFactorAlphaTmp >
                0]) * 0.975
            }
            varData[, Dimnames] = varData[, Dimnames] * scalingFactorAlpha
        }
        plot = plot + geom_text(data = varData,
        mapping = aes_string(x = names(varData)[1],
            y = names(varData)[2], label = "label"), inherit.aes = FALSE,
            size = varLabSize)
    } else {
        varData = NULL
    }

    ## AXIS LABELS
    if (plotPsi == "psi") {
        plot = plot + xlab(bquote(psi[.(Dim[1])] == .(round(x$psis[Dim[1]],
            1)))) + ylab(bquote(psi[.(Dim[2])] == .(round(x$psis[Dim[2]],
            1))))
    } else if (plotPsi == "loglik") {
        liksTab = liks(x)
        if (length(x$confModelMat)) {
            # If filtered on confounders, print in title.
            plot = plot + ggtitle(paste0("Confounders' deviance explained: ",
            liksTab["logLikExplained", "filtered"] * 100, "%"))
        }
        plot = plot + xlab(paste0(Dimnames[1], ": ", liksTab["logLikExplained",
            Dimnames[1]] * 100, "%")) + ylab(paste0(Dimnames[2], ": ",
            liksTab["logLikExplained",
            Dimnames[2]] * 100, "%"))
    } else if (plotPsi == "inertia") {
        inertTab = inertia(x)
        if (length(x$confModelMat)) {
            # If filtered on confounders, print in title.
            plot = plot + ggtitle(paste0("Confounders' inertia explained: ",
                inertTab["inertiaExplained", "filtered"] * 100, "%"))
        }
        plot = plot + xlab(paste0(Dimnames[1], ": ",inertTab["inertiaExplained",
            Dimnames[1]] * 100, "%")) + ylab(paste0(Dimnames[2], ": ",
            inertTab["inertiaExplained",
            Dimnames[2]] * 100, "%"))
    } else if (plotPsi == "none") {
        plot = plot + xlab(paste0(names(dataSam)[1])) +
        ylab(paste0(names(dataSam)[2]))
    } else {
        stop("'plotPsi' argument unknown!\n")
    }

    # Add cross in the centre
    plot = plot + geom_point(data = data.frame(x = 0, y = 0),aes_string(x = "x",
        y = "y"), size = crossSize, inherit.aes = FALSE, shape = 3)
    # Enlarge most text
    plot = plot + theme_bw() +
    theme(axis.title = element_text(size = axisTitleSize),
        axis.text = element_text(size = axisLabSize),
        legend.title = element_text(size = legendTitleSize),
        legend.text = element_text(size = legendLabSize),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

    # Fix coordinates at a certain aspect ratio if required, and throw
    # warning if not
    if (axesFixed)
        plot = plot + coord_fixed(ratio = aspRatio)
    if (!(axesFixed & (aspRatio == 1)))
        warning("Axes not squared, plot may be deformed!\nConsider
setting aspRatio = 1 and axesFixed = TRUE.")
    # Expand limits to show all text
    plot = indentPlot(plot, xInd = xInd, yInd = yInd)
    if (returnCoords) {
        list(plot = plot, samples = dataSam, species = dataTax,
        variables = varData)
    } else {
        plot
    }
}
