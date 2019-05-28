#' Estimate the overdispersion
#'
#' @details Information between taxa is shared with empirical Bayes
#'  using the edgeR pacakage, where the time-limiting steps are programmed in C.
#'
#' @param X the data matrix of dimensions nxp
#' @param cMat a 1xp colum scores matrix
#' @param rMat a nx1 rowscores matrix, if unconstrained
#' @param muMarg an nxp offset matrix
#' @param psis a scalar, the current psi estimate
#' @param trended.dispersion a vector of length p with pre-calculated
#'  trended.dispersion estimates. They do not vary in function
#'  of the offset anyway
#' @param prior.df an integer, number of degrees of freedom of the prior
#'  for the Bayesian shrinkage
#' @param dispWeights Weights for estimating the dispersion
#' in a zero-inflated model
#' @param rowMat matrix of row scores in case of constrained ordination
#'
#' @return A vector of length p with dispersion estimates
estDisp = function(X, cMat = NULL, rMat = NULL,
    muMarg, psis, trended.dispersion = NULL,
    prior.df = 10, dispWeights = NULL, rowMat = NULL,
    allowMissingness = FALSE) {
    logMeansMat = if (!is.null(rMat)) {
        # Unconstrained
        t(rMat %*% (cMat * psis) + log(muMarg))
    } else if (is.null(rowMat)) {
        t(log(muMarg))  #Non-parametric
    } else {
        # Constrained
        t(log(muMarg) + psis * rowMat)
    }
    if (any(is.infinite(logMeansMat)))
        stop("Overflow! Try trimming more lowly
        abundant taxa prior to model fitting.
        \n See prevCutOff argument in ?RCM.")
    X = correctXMissingness(X, exp(logMeansMat), allowMissingness)

    trended.dispersion = if (is.null(trended.dispersion)) {
        edgeR::estimateGLMTrendedDisp(y = t(X),
            design = NULL, method = "bin.loess",
            offset = logMeansMat, weights = NULL)
    } else {
        trended.dispersion
    }
    trended.dispersion = if (is.list(trended.dispersion)) {
        trended.dispersion$dispersion
    } else trended.dispersion

    thetaEstsTmp <- edgeR::estimateGLMTagwiseDisp(y = t(X),
        design = NULL, prior.df = prior.df,
        offset = logMeansMat, dispersion = trended.dispersion,
        weights = dispWeights)

    thetaEsts = if (is.list(thetaEstsTmp)) {
        1/thetaEstsTmp$tagwise.dispersion
    } else {
        1/thetaEstsTmp
    }
    if (anyNA(thetaEsts)) {
        idNA = is.na(thetaEsts)
        thetaEsts[idNA] = mean(thetaEsts[!idNA])
        warning(paste(sum(idNA), "dispersion estimations did not converge!"))
    }
    return(thetas = thetaEsts)
}
