#' Wrapper function for the RCM() function
#'
#' @param dat an nxp count matrix or a phyloseq object with an otu_table slot
#' @param k an integer, the number of dimensions of the RCM solution
#' @param round a boolean, whether to round to nearest integer. Defaults to
#'   FALSE.
#' @param prevCutOff a scalar, the prevalance cutoff for the trimming.
#'  Defaults to 2.5e-2
#' @param minFraction a scalar, each taxon's total abundance
#' should equal
#'  at least the number of samples n times minFraction,
#'   otherwise it is trimmed.
#'   Defaults to 10\%
#' @param rowWeights,colWeights character strings,
#' the weighting procedures for the normalization of row and column scores.
#'  Defaults to 'uniform' and 'marginal' respectively
#' @param covariates In case 'dat' is a phyloseq object,
#' the names of the sample
#'   variables to be used as covariates in the constrained analysis,
#'    or 'all' to
#'   indicate all variables to be used.
#'   In case 'dat' is a matrix, a nxf matrix
#'   or dataframe of covariates.
#'    Character variables will be converted to
#'   factors, with a warning. Defaults to NULL,
#'   in which case an unconstrained
#'   analysis is carried out.
#' @param confounders In case 'dat' is a phyloseq object,
#' the names of the sample variables to be used as confounders
#'  to be filtered
#'  out. In case 'dat' is a matrix, a nxf dataframe
#'  of confounders.
#'   Character variables will be converted to factors, with a warning.
#'   Defaults to NULL, in which case no filtering occurs.
#' @param confTrimMat,confModelMat,covModelMat,centMat Dedicated model matrices
#'  constructed based on phyloseq object.
#' @param ... Further arguments passed on to the RCM.NB() function
#'
#'@description This is a wrapper function,
#'which currently only fits the negative binomial distribution,
#'but which could easily be extended to other ones.
#'
#'@details This function should be called on a raw count matrix,
#' without rarefying or normalization to proportions.
#' This functions trims on prevalence and total abundance to avoid instability
#' of the algorithm. Covariate and confounder matrices are constructed,
#' so that everything is passed on
#'  to the workhorse function RCM.NB() as matrices.
#'@seealso \code{\link{RCM_NB}},\code{\link{plot.RCM}},
#'\code{\link{residualPlot}},\code{\link{plotRespFun}}
#'
#' @return see \code{\link{RCM_NB}}
#'
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, round = TRUE)
#'
#' @importFrom nleqslv nleqslv
#' @importFrom tensor tensor
#' @importFrom alabama constrOptim.nl
#' @import VGAM
#' @import phyloseq
#' @import methods
#' @importFrom stats as.formula contrasts dnbinom dpois glm integrate quantile
#'
#' @rdname RCM
#' @export
setGeneric("RCM", function(dat, ...) standardGeneric("RCM"))

#' @rdname RCM
#' @export
setMethod("RCM", "phyloseq", function(dat, covariates = NULL,
    confounders = NULL, ...) {
    # Safer than using @.Data, as this uses the
    # standard conversion function of the package as
    # opposed to the internal representation.
    X <- as(otu_table(dat), "matrix")
    if (taxa_are_rows(dat))
        X <- t(X)

    ## Build confounder matrix if applicable. ##
    if (!is.null(confounders)) {
        tmp = buildConfMat(confounders, dat)
        confTrimMat = tmp$confModelMatTrim
        confModelMat = tmp$confModelMat
        rm(tmp)
    } else {
        confModelMat = confTrimMat = NULL
    }

    ## Build covariate matrix if applicable##
    if (!is.null(covariates)) {
        tmp = buildCovMat(covariates, dat)
        covModelMat = tmp$covModelMat
        covariates = tmp$covariates
        # Already prepare the matrix that defines the
        # equations for centering the coefficients of the
        # dummy variables
        tmp = buildCentMat(tmp$datFrame)
        centMat = tmp$centMat
        datFrame = tmp$datFrame
        rm(tmp)

        # Remove rows with NA's, we might want to find
        # something better or leave it to the end user
        if (anyNA(datFrame)) {
            NArows = apply(datFrame, 1, anyNA)
            if (all(NArows)) {
                stop("All samples have missing covariates")
            }
            X = X[!NArows, ]
            datFrame = datFrame[!NArows, ]
            if (!is.null(confModelMat)) {
                confModelMat = confModelMat[!NArows,]
                confTrimMat = confTrimMat[!NArows,]
            }
            warning(paste("Some covariates contain missing values.
                    We removed samples \n",
                paste(which(NArows), collapse = ", "),
                "\n prior to analysis."), immediate. = TRUE)
        }
    } else {
        covModelMat = centMat = NULL
    }

    # Call generic again to pass on data.
    tmp = RCM(X, confModelMat = confModelMat, confTrimMat = confTrimMat,
        covModelMat = covModelMat, centMat = centMat,
        ...)
    tmp$physeq = prune_samples(rownames(tmp$X), prune_taxa(colnames(tmp$X),
        dat))
    return(tmp)
})

#' @rdname RCM
#' @export
setMethod("RCM", "matrix", function(dat, k = 2, round = FALSE,
    prevCutOff = 0.05, minFraction = 0.1, rowWeights = "uniform",
    colWeights = "marginal", confModelMat = NULL, confTrimMat = NULL,
    covModelMat = NULL, centMat = NULL, ...) {

    if (anyNA(dat)) {
    stop("NA values present in count matrix,
        please filter these out first!\n")
    }
    p = ncol(dat)
    n = nrow(dat)

    if (is.null(colnames(dat))) {
        colnames(dat) = seq_len(ncol(dat))
    }
    if (is.null(rownames(dat))) {
        rownames(dat) = seq_len(nrow(dat))
    }
    colNames = colnames(dat)
    rowNames = rownames(dat)
    if (round)
        {
            dat = round(dat, 0)
        }  #Round to integer

    # Check dat type
    if (!all(floor(dat) == dat)) {
        stop("Please provide integer count matrix
                (not a matrix of proportions!),
                or set 'round' to TRUE! \n")
    } else {
        dat = matrix(as.integer(dat), ncol = p, nrow = n)
    }

    colnames(dat) = colNames
    rownames(dat) = rowNames
    dat = dat[, (colMeans(dat == 0) < (1 - prevCutOff)) &
        (colSums(dat) > (n * minFraction))]
    # Remove taxa with low prevalence and those who do
    # not meet the minFraction requirement
    if (any(rowSums(dat) == 0)) {
        warning(immediate. = TRUE, paste0("Samples \n",
            paste(collapse = ", ", rownames(dat)[rowSums(dat) ==
                0]), "\n contained no more reads after trimming taxa
                with low prevalence, and were excluded from the fit"))
    }
    rowIDkeep = rowSums(dat) > 0
    dat = dat[rowIDkeep, ]  #Remove empty samples
    n = nrow(dat)
    confModelMat = confModelMat[rowIDkeep, ]
    confTrimMat = confTrimMat[rowIDkeep, ]
    attribs = attr(covModelMat, "assign")
    covModelMat = covModelMat[rowIDkeep, ]

    tic = proc.time()  #Time the calculation
    tmp = RCM_NB(dat, rowWeights = rowWeights, colWeights = colWeights,
        k = k, confModelMat = confModelMat, confTrimMat = confTrimMat,
        covModelMat = covModelMat, prevCutOff = prevCutOff,
        minFraction = minFraction, centMat = centMat,
        ...)
    tmp = within(tmp, {
        runtimeInMins = (proc.time() - tic)[1]/60  # The runtime
        k = k  #Store number of dimensions
    })
    tmp$call = match.call()
    tmp$attribs = attribs
    class(tmp) = "RCM"
    tmp
})
