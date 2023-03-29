#' Perform a PERMANOVA analysis for group differences of a predefined cofactor using the pseudo F-statistic
#'
#' @param rcmObj an RCM object
#' @param groups a factor of length n with cluster memberships, or a name of a variable contained in the RCM object
#' @param nPerm Number of permutations in the PERMANOVA, defaults to 1e6
#' @param Dim Dimensions on which the test should be performed. Defaults to all dimensions of the fitted RCM object.

#' @return A list with components
#' \item{statistic}{The pseudo F-statistic}
#' \item{p.value}{The p-value of the PERMANOVA}
#'
#' @seealso \code{\link{RCM}}
#' @importFrom phyloseq get_variable
#' @importFrom stats dist
#' @export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, round = TRUE)
#' zellerPermanova = permanova(zellerRCM, "Diagnosis")
permanova = function(rcmObj, groups, nPerm = 1e6, Dim = seq_len(rcmObj$k)){
  stopifnot(is(rcmObj, "RCM"), length(nPerm)==1)
    if(nPerm <= 100){
        warning("Less than 100 permutations leads to low power of the permutation test!")
    }
  if(length(groups)==1){
      groups = get_variable(rcmObj$physeq, groups)
  }
  coord = extractCoord(rcmObj, Dim)$samples
  N = nrow(coord)
  if(N != length(groups)){
      stop("Length of grouping variable provided (", length(groups),
           "does not correspond to number of samples in RCM object (", N, ")")
  }
  a = length(unique(groups))
  if(a <= 1){
      stop("Provide more than one group in 'groups'.")
  }
  if(any(table(groups))==1){
      stop("Some groups contain only a single observation, no distances can be calculated.")
  }
  distSq = dist(coord)^2
  overalDist = sum(distSq)/N
  #Observed test statistic
  withinDistObs = sum(unlist(tapply(seq_len(nrow(coord)), groups, function(x){
        distSq[getDistCoord(x, N)]/length(x)
  })))
  FratioObs = (overalDist-withinDistObs)/withinDistObs * (a-1)/(N-a)
  #PERMANOVA
  withinDistPerm = vapply(integer(nPerm), FUN.VALUE = double(1), function(jj){
      sum(unlist(tapply(seq_len(nrow(coord)), sample(groups), function(x){
        distSq[getDistCoord(x, N)]/length(x)
        })))
  })
  FratioPerm = (overalDist-withinDistPerm)/withinDistPerm * (a-1)/(N-a)
  PvalPerm = mean(FratioObs < FratioPerm)
  return(list("statistic" = Fratio, "p.value" = PvalPerm))
}
