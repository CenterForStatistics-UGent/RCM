#'A simple wrapper function for the RC(M) object for phyloseq objects and matrices.
#'
#' @param dat an nxp count matrix or a phyloseq object with an otu_table slot
#' @param k an integer, the number of dimensions of the RCM solution
#' @param round a boolean, whether to round to nearest integer. Defaults to FALSE.
#' @param distribution a character string, the error distribution of the RC(M) model. Defaults to NB
#' @param prevCutOff a scalar, the prevalance cutoff for the trimming. Defaults to 2.5e-2
#' @param minFraction a scalar, each taxon's total abundance should equal at leat the number of samples n times minFraction, otherwise it is trimmed. Defaults to 10\%
#' @param rowWeights,colWeights character strings, the weighting procedures for the normalization of row and column scores. Defaults to "uniform" and "marginal" respectively
#' @param covariates In case "dat" is a phyloseq object, the names of the sample variables to be used as covariates in the constrained analysis, or "all" to indicate all variables to be used. In case dat is a matrix, a nxf matrix or dataframe of covariates. Character variables will be converted to factors, with a warning. Defaults to NULL, in which case an unconstrained analysis is carried out.
#' @param confounders In case "dat" is a phyloseq object, the names of the sample variables to be used as confounders to be filtered out. In case dat is a matrix, a nxf matrix or dataframe of confounders Character variables will be converted to factors, with a warning. Defaults to NULL, in which case no filtering occurs.
#' @param prevFit An object with a previous fit, normally from a lower dimension, that should be extended.
#' @param ... Further arguments passed on to the RCM.NB() function
#'
#'@description This is a wrapper function, which currently only fits the negative binomial distribution, but which could easily be extended to other ones.
#'
#'@details This functions trims on prevalence and total abundance to avoid instability of the algorithm. Covariate and confounder matrices are constructed, so that everything is passed on to the workhorse function RCM.NB() as matrices.
#'@seealso \code{\link{RCM_NB}},\code{\link{plot.RCM}},\code{\link{residualPlot}},\code{\link{plotRespFun}}
#'
#' @return see \code{\link{RCM_NB}}
#' @importFrom nleqslv nleqslv
#' @importFrom tensor tensor
#' @importFrom alabama constrOptim.nl
#' @import VGAM
#' @import phyloseq
#' @importFrom stats as.formula contrasts dnbinom dpois glm.fit integrate quantile
#' @export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, k = 2, round = TRUE)
RCM = function(dat, k, round=FALSE, distribution= "NB", prevCutOff = 0.025, minFraction = 0.1, rowWeights = "uniform", colWeights = "marginal", covariates = NULL, confounders = NULL, prevFit = NULL, ...){
  classDat = class(dat) #The class provided

  ##The count data##
  if (classDat=="matrix"){
    X = dat
  }else  if(classDat=="phyloseq"){
    X = if (taxa_are_rows(dat)) t(otu_table(dat)@.Data) else otu_table(dat)@.Data
  } else {stop("Please provide a matrix or a phyloseq object! \n")}
  p=ncol(X); n=nrow(X)

  if(is.null(colnames(X))){colnames(X)=1:ncol(X)}
  if(is.null(rownames(X))){rownames(X)=1:nrow(X)}
  colNames = colnames(X); rowNames =rownames(X)
  if(round) {X=round(X, 0) }#Round to integer

  #Check X type
  if(!all(floor(X)==X)){stop("Please provide integer count matrix, or set 'round' to TRUE! \n")
  } else{X=matrix(as.integer(X), ncol=p, nrow=n)}

  colnames(X)=colNames; rownames(X)=rowNames
  X=X[rowSums(X)>0, ] #Remove empty samples
  X=X[, (colMeans(X==0)<(1-prevCutOff)) & (colSums(X)>(n*minFraction))] #Remove txa with low prevalence and those who do not meet the minFraction requirement
  if (distribution %in% c("ZIP","ZINB")) X = X[rowSums(X==0)>0, colSums(X==0)>0] #For a zero-inflated model, make sure every row and column has zeroes

  ##Build confounder matrix if applicable. ##
  if(!is.null(confounders)){
tmp = buildConfMat(confounders, n, dat)
confModelMatTrim = tmp$confModelMatTrim
confModelMat = tmp$confModelMat
rm(tmp)
  } else{
    confModelMat = confModelMatTrim = NULL
  }

  ##Build covariate matrix if applicable##
  if(!is.null(covariates)){
    tmp = buildCovMat(covariates, n, dat)
    covModelMat = tmp$covModelMat
    covariates = tmp$covariates
    #Already prepare the matrix that defines the equations for centering the coefficients of the dummy variables
    tmp = buildCentMat(tmp$datFrame)
    centMat = tmp$centMat
    datFrame = tmp$datFrame
    rm(tmp)

  #Remove rows with NA's, we might want to find something better or leave it to the end user
  if(anyNA(datFrame)){
    NArows = apply(datFrame, 1, anyNA)
    if(all(NArows)){ stop("All samples have missing covariates")}
    X = X[!NArows,]
    datFrame = datFrame[!NArows,]
    if(!is.null(confModelMat)){
      confModelMat = confModelMat[!NArows,]
      confModelMatTrim = confModelMatTrim[!NArows,]
    }
    warning(paste("Some covariates contain missing values. We removed samples \n", paste(which(NArows), collapse = ", "), "\n prior to analysis." ),immediate. = TRUE)
    }
  } else {
    covModelMat = centMat = NULL
  }

  tic = proc.time() #Time the calculation
  tmp = RCM_NB(X, rowWeights = rowWeights, colWeights = colWeights, k = k,
                         confounders  = list(confounders = confModelMat, confoundersTrim = confModelMatTrim),
                         covariates = covModelMat, prevCutOff = prevCutOff, minFraction = minFraction, centMat = centMat, NBRCM = prevFit,...)
  if(classDat=="phyloseq"){tmp$physeq = prune_samples(rownames(X),prune_taxa(colnames(X),dat)) }
  tmp = within(tmp, {
    runtimeInMins = (proc.time()-tic)[1]/60 + if(is.null(prevFit)) {0} else {prevFit$runtimeInMins} # Sum the runtimes
    k = k #Store number of dimensions
  })
  class(tmp) = "RCM"
  tmp
}
