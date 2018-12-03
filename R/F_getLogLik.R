#' Extract the logged likelihood of every count
#'
#' @param rcm an RCM object
#' @param Dim A vector of integers indicating which dimensions to take along,
#' or Inf for the saturated model, or 0 for the independence model
#'
#' @return A matrix with logged likelihood of the size of the data matrix
getLogLik = function(rcm, Dim) {
    if (Dim == Inf) {
        return(dpois(x = rcm$X, lambda = rcm$X, 
            log = TRUE))
    }
    E = extractE(rcm, if (Dim >= 1) 
        seq_len(Dim) else Dim)
    thetaMat = matrix(byrow = TRUE, nrow = nrow(rcm$X), 
        ncol = ncol(rcm$X), data = rcm$thetas[, 
            switch(as.character(Dim), `0` = "Independence", 
                `0.5` = "Filtered", paste0("Dim", 
                  Dim))])
    dnbinom(x = rcm$X, mu = E, size = thetaMat, 
        log = TRUE)
}
