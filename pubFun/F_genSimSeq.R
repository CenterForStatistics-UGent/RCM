#' A function to generate datasets non-parametrically based on lfdr's from Wilcoxon rank sum tests
#'
#' @param testRes A result from the testSimSeq() function
#' @param Ntaxa The total number of taxa in the final dataset. Defaults to the same as the original dataset
#' @param fracTaxa the fraction of taxa to be made differentially abundant
#'
#' @return a list with componenents
#' \item{data}{The newly generated dataset}
#' \item{treatment}{The grouping vector}
#' \item{DEtaxa}{A character vector for taxa made differenatially abundant}
genSimSeq = function(testres, Ntaxa = ntaxa(testres$physeq), fracTaxa = 0.2){
  physeq = testres$physeq
  weights = 1-testres$lfdr
  weights = weights/sum(weights) #Normalize sampling weights
  normFactors = sample_sums(physeq)
  treatment = get_variable(physeq, testres$variable)
  DEtaxa = sample(taxa_names(physeq), round(Ntaxa*fracTaxa), prob = weights) # The DE taxa are sampled from all taxa with the 1-lfdr as sampling weights. The NDE taxa are simply the remaining ones
  whichTrt = names(sort(table(treatment),decreasing=TRUE)[1]) # The most frequent treatment group, frow which we will sample the EE taxa
  trtLeft = unique(treatment)[unique(treatment) != whichTrt] #The remaining treatment groups
  countMat = if(taxa_are_rows(physeq)){t(otu_table(physeq)@.Data)} else {otu_table(physeq)@.Data}
  samIDee = sample(which(treatment == whichTrt),nrow(countMat), replace = TRUE) #NDA samples are taken from the first treatment group
  dataSimSeq = countMat[samIDee,] #EE taxa
  for (i in trtLeft){
    rowInd = which(treatment==i)
    dataSimSeq[rowInd,DEtaxa] = round(countMat[rowInd, DEtaxa] * normFactors[samIDee][rowInd]/normFactors[rowInd]) #DA taxa are taken from the other groups, adapting the sequencing depth to the most frequent one. samIDee refers to the sample indices of the new matrix
  }
  keepID = colSums(dataSimSeq)>0
  list(data = dataSimSeq[, keepID], treatment = treatment, DEtaxa = DEtaxa[DEtaxa %in% colnames(dataSimSeq[, keepID])])
}