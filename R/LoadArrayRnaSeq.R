# LoadArrayRnaSeq.R
# February 25 2015

LoadArrayRnaSeq <- function(normalized.array.path, rna.seq.path, fix.rik.xgene = FALSE){
  scripts.dir <- "/home/yeung/projects/tissue-specificity/scripts"
  funcs.dir <- "functions"
  source(file.path(scripts.dir, funcs.dir, "LoadAndHandleData.R"))
  source(file.path(scripts.dir, funcs.dir, "MergeToLong.R"))
  source(file.path(scripts.dir, funcs.dir, "GrepRikGenes.R"))
  # Define dirs -------------------------------------------------------------
  
  if (missing(normalized.array.path) & missing(rna.seq.path)){
    # define dirs
    data.dir <- "/home/yeung/projects/tissue-specificity/data"
    normalized.array.fname <- "array.adj.0.07.txt"
    normalized.array.path <- file.path(data.dir, normalized.array.fname)
    rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
    rna.seq.path <- file.path(data.dir, rna.seq.fname) 
  }
  
  # Load file ---------------------------------------------------------------
  
  normalized.array <- LoadNormalizedArray(normalized.array.path)
  rna.seq.exprs <- LoadRnaSeq(rna.seq.path)
  
  
  # Log2 transform of array and rnaseq --------------------------------------
  
  normalized.array <- log2(normalized.array + 1)
  rna.seq.exprs <- log2(rna.seq.exprs + 1)
  
  rna.seq.exprs.filtered <- rna.seq.exprs  # no filter
  
  # Take only common genes --------------------------------------------------
  
  array.genes <- rownames(normalized.array)
  rna.seq.genes <- rownames(rna.seq.exprs.filtered)
  
  common.genes <- intersect(array.genes, rna.seq.genes)
  
  print(paste(length(common.genes), "common genes between array and rnaseq"))
  
  normalized.array <- normalized.array[common.genes, ]
  rna.seq.exprs.filtered <- rna.seq.exprs.filtered[common.genes, ]
  
  # Merge data into long format ---------------------------------------------
  
  dat <- MergeToLong(normalized.array, rna.seq.exprs.filtered)  
  
  # Fix Rik genes ---------------------------------------------
  if (fix.rik.xgene){
    dat$gene <- FixRikGenes(dat$gene)
  }
  
  return(dat)
}
