#' Calculate the log-likelihoods of all possible models
#'
#'@param rcm an object of the RCM class
#'
#'@return A table with inertias, proportion inertia explained
#' and cumulative proportion of inertia explained.
#'@export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, round = TRUE)
#' inertia(zellerRCM)
inertia = function(rcm) {
    vec = if (length(rcm$confModelMat))
        c(0, 0.5, seq_len(rcm$k)) else c(0:rcm$k)
    outnames = c("independence", if (length(rcm$confModelMat))
    "filtered" else NULL,
        paste0("Dim", seq_len(rcm$k)), "saturated")
    tmp = c(vapply(FUN.VALUE = numeric(1), vec, FUN = function(i) {
        eMat = extractE(rcm, i)
        round(sum(((eMat - rcm$X)^2/eMat)))
    }), 0)
    names(tmp) = outnames
    cumInertiaExplained = round((tmp - tmp[1])/(tmp[length(tmp)] - tmp[1]),
        3)
    out = rbind(inertia = tmp, inertiaExplained =
    c(0, diff(cumInertiaExplained)),
        cumInertiaExplained = cumInertiaExplained)
    return(out)
}
