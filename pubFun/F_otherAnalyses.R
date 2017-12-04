#' A function to perform all biplot analyses except for RC(M)
#'
#' @param RCMlist a list of RCM objects
#' @param unifrac a boolean, should phylogenetic methods be used?
#' @param a list of phylygenetic trees belonging to the phyloseq objects in RCMlist
#' @param cores an integer, number of cores to use
#'
#' @return a list of with results of the otherAnalysesSimple() function
otherAnalyses = function(RCMlist, gommList, unifrac = FALSE, treeList = NULL, cores = 1){
  if(unifrac) {
    RCMlist = mcmapply(mc.cores = cores, SIMPLIFY = FALSE,RCMlist, treeList, FUN = function(x,tree){
      x$physeq = phyloseq(otu_table(x$X, taxa_are_rows = FALSE), tree)
      sample_names(x$physeq) = seq_len(nsamples(x$physeq))
      x
    })
  }
  tmpList = mclapply(mc.cores = cores, RCMlist, otherAnalysesSimple, unifrac = unifrac)
  mapply(tmpList, gommList, FUN = function(a,b){a$gomm = b;a})
} #JSD works on relative abundances by default