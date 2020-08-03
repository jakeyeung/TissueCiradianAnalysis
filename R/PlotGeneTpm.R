PlotGeneTpm <- function(dat, jgene, log2.pseudocount=FALSE, scale=FALSE){
  dat.sub <- subset(dat, gene == jgene)
  if (scale != FALSE){
    dat.sub$tpm <- dat.sub$tpm * scale
  }
  if (log2.pseudocount != FALSE){
    dat.sub$tpm <- log2(dat.sub$tpm + log2.pseudocount)
  }
  
  ggplot(dat.sub, aes(x = time, y = tpm)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jgene)
}