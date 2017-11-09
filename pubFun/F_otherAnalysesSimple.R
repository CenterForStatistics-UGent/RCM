#' Perform the other analyses on a single RCM object
#'
#' @param x an RCM object
#' @param unifrac are phylogentic analyses required?
#'
#' @return a list with components
#' \item{RCM}{the supplied RCM object}
#' \item{CA}{Correspondence analysis}
#' \item{DCA}{Detrended correspondence analysis}
#' \item{CoDa}{Compositional data analysis}
#' \item{BC}{Bray-Curtis on absolute counts}
#' \item{}
#' \item{BCrel}{Bray-Curtis on relative abundances}

otherAnalysesSimple = function(x, unifrac){
  tmp = x$X; rownames(tmp) = seq_len(nrow(tmp))
  R = rowSums(tmp)
  E = outer(R, colSums(tmp)/sum(tmp))
  otuTab = otu_table(tmp, taxa_are_rows = FALSE)
  otuTabRel = otu_table(tmp/R, taxa_are_rows = FALSE)
  otuTabLog = otu_table(log(tmp+1), taxa_are_rows = FALSE)
  list(
    RCM = x, CA = caSVD(tmp), DCA = ordinate(otuTab, method = "DCA"), CoDa = CoDa(tmp),
    BC = ordinate(otuTab, method = "PCoA", distance = "bray"),
    BClog = ordinate(otuTabLog, method = "PCoA", distance = "bray"),
    JSD = ordinate(otuTab, method = "PCoA", distance = "jsd"),
    BCrel = ordinate(otuTabRel, method = "PCoA", distance = "bray"),
    BCrelNMDS = ordinate(otuTabRel, method = "NMDS", distance = "bray", k=3),
    Hellinger = svd(diag(1/sqrt(R)) %*% (sqrt(tmp)-sqrt(E))),
    UniFrac = if(unifrac) ordinate(x$physeq, method = "PCoA", distance = "unifrac") else NULL,
    wUniFrac = if(unifrac) ordinate(x$physeq, method = "PCoA", distance = "wunifrac") else NULL,
    DPCOA = if(unifrac) DPCoA(x$physeq) else NULL
  )
}