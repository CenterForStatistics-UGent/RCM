#' Perform the compostional dimension reduction
#'
#' @param dat a count matrix or phyloseq object
#' @param abCutOff a scalar, the abundance cut off
#' @param scale a boolean, should variables be scaled prior to PCA? The default is FALSE, as in Gloor and Reid (2016)
#'
#' The compositional approach is as follows: the zeroes are fixed somehow, log-ratio transform is applied to relative abundances and PCA is carried out on these transformed values. Our main concern is that this approach ignores the fact that we have count data that are heteroscedastic. This function is build after the supplementary files of Gloor and Reid, 2016. Only the abundance filtering is made less severe
#'
#' @return a list with components
#' \item{rowScores}{The row scores of the ordination}
#' \item{colScores}{The column scores of the ordination}
CoDa = function(dat, abCutOff = 1e-6, scale = FALSE){
  require(zCompositions)
  classDat = class(dat) #The class provided

  ##The count data##
  if (classDat=="matrix"){
    X = dat
  }else  if(classDat=="phyloseq"){
    X = if (taxa_are_rows(dat)) t(otu_table(dat)@.Data) else otu_table(dat)@.Data
  } else {stop("Please provide a matrix or a phyloseq object! \n")}

  Xzero <- as.matrix(cmultRepl(X = X, method="CZM", output="counts")) #Correct zeroes in a Bayesian way

  # convert to proportions
  Xnorm = Xzero/rowSums(Xzero)

  #Trim on abundance in any sample
  Xtrim = Xnorm[,apply(Xnorm,2,min)>abCutOff]

  # make our compositional dataset
  Xcrl <- t(apply(Xtrim, 1, function(x){log(x) - mean(log(x))}))

  #Column center
  Xsca = scale(Xcrl, center = TRUE, scale = FALSE)

  # Apply PCA
  Svd = svd(Xsca)
  rowScores = Svd$u %*% diag(Svd$d) #=Xsca %*% Svd$v #%*% diag(Svd$d)
  colScores = Svd$v #scale = 0 in biplot
  rownames(rowScores) = rownames(Xsca)
  rownames(colScores) = colnames(Xsca)

  list(rowScores = rowScores, colScores = colScores)
}