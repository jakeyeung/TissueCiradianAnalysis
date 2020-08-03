LoadKallistoGene <- function(inpath, gene_colname = "gene_name", log2.pseudocount=FALSE, form = "long"){
  source("/home/yeung/projects/tissue-specificity/scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  source("/home/yeung/projects/tissue-specificity/scripts/functions/GetTissueTimes.R")
  
  if (missing(inpath)){
    inpath <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
  }
  if (form == "wide"){
    dat <- read.table(inpath, header = TRUE, row.names = 1)
    tissues <- sapply(colnames(dat), function(s) strsplit(s, "_")[[1]][[1]])
    tissues <- ConvertRNASeqTissueNamesToArray(tissues)
    times <- GetTimes(colnames(dat), get_unique = FALSE)
    colnames(dat) <- paste(tissues, times, sep = '')
    return(dat)
  }
  
  dat <- read.table(inpath, header = TRUE)
  
  genes <- dat[, gene_colname]
  dat.mat <- dat[, colnames(dat)[which(colnames(dat) != gene_colname)]]
  
  if (log2.pseudocount != FALSE){
    dat.mat <- log2(dat.mat + log2.pseudocount)  # add a pseudocount
  }
  
  tissues <- sapply(colnames(dat.mat), function(s) strsplit(s, '_')[[1]][[1]])
  tissues <- ConvertRNASeqTissueNamesToArray(tissues)
  times <- GetTimes(colnames(dat.mat), get_unique=FALSE)
  
  tpm.long <- data.frame(gene = rep(genes, ncol(dat.mat)),
                         tissue = rep(tissues, each = nrow(dat.mat)),
                         time = as.numeric(rep(times, each = nrow(dat.mat))),
                         tpm = unlist(dat.mat))
}