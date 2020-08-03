LoadNormalizedArray <- function(normalized.array.path, remove.negs=TRUE, fix.rik.xgene = FALSE){
  # Normalizing has caused strange behaviours in very flat genes, causing
  # the resulting output of one or two genes to be negative.
  # We need to remove these from subsequent analysis.\
  # fix.rik.xgene optional
  scripts.dir <- "/home/yeung/projects/tissue-specificity/scripts"
  funcs.dir <- "functions"
  source(file.path(scripts.dir, funcs.dir, "GrepRikGenes.R"))
  
  normalized.array <- read.table(normalized.array.path)
  genes <- rownames(normalized.array)
  if (remove.negs){
    # How many have negative values? ------------------------------------------
    negs <- apply(normalized.array, 1, function(x){
      if (min(x) < 0){
        return(1)
      } else {
        return(0)
      }
    })
    
    problem.genes <- names(negs[which(negs == 1)])
    
    # Remove problem genes from analysis -------------------------------------
    filtered.genes <- genes[which(!genes %in% problem.genes)]
    normalized.array <- normalized.array[filtered.genes, ]
    
    print(paste("Removed problem genes:", problem.genes)) 
  }
  
  # optionally fix gene names
  if (fix.rik.xgene){
    rownames(normalized.array) <- FixRikGenes(rownames(normalized.array))
  }
  return(normalized.array)
}

LoadRnaSeq <- function(rna.seq.path, handle.duplicates=TRUE){
  if (missing(rna.seq.path)){
    rna.seq.path <- "data/rna_seq_deseq_counts_colnames_fixed.txt"
  }
  # tricky duplicate rownames. We'll fix that though.
  
  rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
  
  # Handle duplicate rownames: RNASEQ ------------------------------------
  
  rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)
  
  drop.cols <- c("gene")
  rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]
  return(rna.seq.exprs)
}