#' An auxiliary function to fix infinite value in a count matrix
#'
#' @param data A matrix, possibly containing infinite values
#'
#'Taken from the SpiecEasi package
#'
#' @return The same matrix, but with infinite entries replaced by the column maximima plus one

.fixInf <- function(data) {
    if (any(is.infinite(data))) {
       data <-  apply(data, 2, function(x) {
              if (any(is.infinite(x))) {
                   x[ind<-which(is.infinite(x))] <- NA
                   x[ind] <- max(x, na.rm=TRUE)+1
                 }
                x
                })
    }
    data
}

