MergeToLong <- function(normalized.array, rna.seq.exprs, take.common.genes = FALSE){
  
  # Functions
  scripts.dir <- "/home/yeung/projects/tissue-specificity/scripts"
  funcs.dir <- "functions"
  source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))
  
  # Take common genes
  if (take.common.genes){
    array.genes <- rownames(normalized.array)
    rna.seq.genes <- rownames(rna.seq.exprs.filtered)
    
    common.genes <- intersect(array.genes, rna.seq.genes)
    
    print(paste(length(common.genes), "common genes between array and rnaseq"))
    
    normalized.array <- normalized.array[common.genes, ]
    rna.seq.exprs <- rna.seq.exprs[common.genes, ]    
  }

  # Get tissue and times
  tissues <- GetTissues(colnames(normalized.array))
  times.array <- GetTimes(colnames(normalized.array))
  times.rnaseq <- GetTimes(colnames(rna.seq.exprs))
  array.genes <- rownames(normalized.array)
  rnaseq.genes <- rownames(rna.seq.exprs)
  
  # Make into long
  long.array <- data.frame(gene=rep(array.genes, length(tissues) * length(times.array)),
                           tissue=as.factor(rep(tissues, each=(length(array.genes) * length(times.array)))),
                           time=as.numeric(rep(times.array, each=length(array.genes))),
                           experiment=as.factor("array"),
                           exprs=unlist(normalized.array))
  
  long.rnaseq <- data.frame(gene=rep(rnaseq.genes, length(tissues) * length(times.rnaseq)),
                            tissue=rep(tissues, each=(length(rnaseq.genes) * length(times.rnaseq))),
                            time=as.numeric(rep(times.rnaseq, each=length(rnaseq.genes))),
                            experiment=as.factor("rnaseq"),
                            exprs=unlist(rna.seq.exprs))
  
  dat <- rbind(long.array, long.rnaseq)
  return(dat)
}