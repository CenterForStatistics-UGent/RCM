#' A wrapper function for phyloseq objects for Wilcoxon rank sum testing for differential abundance, yields the lfdr's of the taxa. We extend the original code to allow for Kruskal-Wallis test.
#'
#' @param physeq a phyloseq object
#' @param variable a character string, the name of the variable present in the phyloseq object for which to test for differential abundance
#' @param filterTaxa a boolean, should taxa be filtered so that there are no taxa with only zeroes in one of the groups defined by variable?
#'
#' @return a list with components
#' \item{lfdr}{The local false discovery rates of the tests}
#' \item{physeq}{The phyloseq object}
#' \item{variable}{The variable names}
#' \item{result}{A vector containing "up" or "down", indicating if the taxon had a higher or a lower mean rank in the first treatment group}
testSimSeq = function(physeq, variable, filterTaxa = FALSE){
  require(fdrtool)
  countMat = if(taxa_are_rows(physeq)){t(otu_table(physeq)@.Data)} else {otu_table(physeq)@.Data}
  treatment = get_variable(physeq, variable)
  if(filterTaxa){
    # Filter physeq to no all-zeroes in either of the groups
    physeq = prune_taxa(x = physeq, taxa = apply(countMat, 2, function(y){
      !any(tapply(y, treatment, function(z){
        all(z==0)
      }))
    }))
  }
  #Filter physeq for NA's in variables
  physeq = prune_samples(x = physeq, !is.na(treatment))
  treatment = factor(get_variable(physeq, variable))
  countMat = if(taxa_are_rows(physeq)){t(otu_table(physeq)@.Data)} else {otu_table(physeq)@.Data}
  pvals = apply(countMat, 2, function(x){kruskal.test(x, g=treatment)$p.value})
  result = apply(countMat, 2, function(x){
    avRank = length(x)/2 + 0.5
    rankX = rank(x)
    sapply(levels(treatment),avRank = avRank, function(y, avRank){
      ifelse(mean(rankX[treatment==y])>avRank, "up","down")
    })}) #Was taxon up or down regulated in group 1 (first level of treatment?
  lfdr = fdrtool(pvals, statistic = "pvalue", plot=FALSE, verbose=FALSE)$lfdr
  list(lfdr = lfdr, physeq = physeq, variable = variable, result = result)
}