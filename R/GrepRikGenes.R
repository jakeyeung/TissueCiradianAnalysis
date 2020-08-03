# 2014-04-28
# Sometimes genes like X1700040L02Rik are tricky because they
# can be X1700040L02Rik or 1700040L02Rik
# Handle them by grepping them and allow switching between the
# two forms

GrepRikGenes <- function(gene.list){
  Xgenes <- gene.list[grep("^X\\w*Rik\\b", gene.list)]
  return(Xgenes)
}

RemoveX <- function(Xgenes){
  genes <- sapply(Xgenes, function(x){
    return(substr(x, 2, nchar(x)))
  })
}

AddX <- function(genes){
  Xgenes <- sapply(genes, function(string){
    return(paste0("X", string))
  })
}

FixRikGenes <- function(gene.list){
  gene.list <- as.character(gene.list)
  Xgenes <- GrepRikGenes(gene.list)
  genes <- RemoveX(Xgenes)
  gene.list[which(gene.list %in% Xgenes)] <- genes
  gene.list <- as.factor(gene.list)
  return(gene.list)
}