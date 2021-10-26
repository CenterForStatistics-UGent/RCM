#' A function to calculate the matrix of deviance residuals.
#'
#' @param RCM an RCM object
#' @param Dim The dimensions to use
#'
#' For the deviance residuals we use the overdispersions from the reduced model.
#' Standard dimensions used are only first and second,
#'  since these are also plotted
#'
#'@export
#'@return A matrix with deviance residuals of the same size
#' as the original data matrix
#'
#'@examples
#'data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:120],
#' prune_samples(sample_names(Zeller)[1:75], Zeller))
#' #Subset for a quick fit
#' zellerRCM = RCM(tmpPhy, k = 2, round = TRUE, prevCutOff = 0.03)
#' devRes = getDevianceRes(zellerRCM)
getDevianceRes = function(RCM, Dim = RCM$k) {
    mu = extractE(RCM, Dim)
    thetaMat = matrix(byrow = TRUE, nrow = nrow(RCM$X),
        ncol = ncol(RCM$X), data = RCM$thetas[,
            switch(as.character(Dim), `0` = "Independence",
                `0.5` = "Filtered", paste0("Dim",
                Dim))])
    getDevMat(X = correctXMissingness(RCM$X, mu, RCM$NApresent, is.na(RCM$X)), thetaMat = thetaMat,
        mu = mu)
}
