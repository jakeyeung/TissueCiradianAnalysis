LongToMat <- function(dat.fit.relamp.sub, value.var = "relamp"){
  M <- dcast(dat.fit.relamp.sub, gene ~ tissue, value.var = value.var)
  rownames(M) <- M$gene
  M$gene <- NULL
  return(M)
}