LoadArray <- function(inpath, gene_colname = "gene", get.norm = FALSE, form = "long"){
  source("scripts/functions/GetTissueTimes.R")
  if (missing(inpath)){
    inpath <- "/home/yeung/projects/tissue-specificity/data/array_exprs_colnames_fixed.best.probe.selected.txt"
  }
  if (!is.na(gene_colname)){
    dat <- read.table(inpath, header=TRUE)
    dat.mat <- dat[, colnames(dat)[which(colnames(dat) != gene_colname)]]
    genes <- dat[[gene_colname]]
  } else {
    dat.mat <- read.table(inpath)
    genes <- rownames(dat.mat)
  }
  if (form == "wide"){
    rownames(dat) <- dat[[gene_colname]]
    dat[[gene_colname]] <- NULL
  } 
  if (form == "wide" & get.norm == FALSE){
    return(dat)
  } else if(form == "wide" & get.norm == TRUE){
    # genes <- as.character(dat[[gene_colname]])
    # dat.mat <- 2 ^ dat.mat
    return(2 ^ dat)
  }  
  tissues <- GetTissues(colnames(dat.mat), get_unique = FALSE)
  times <- GetTimes(colnames(dat.mat), get_unique = FALSE)
  
  array.long <- data.frame(gene = rep(genes, ncol(dat.mat)),
                           tissue = rep(tissues, each = nrow(dat.mat)),
                           time = as.numeric(rep(times, each = nrow(dat.mat))),
                           experiment = "array",
                           exprs = unlist(dat.mat))
#                          signal = unlist(dat.mat))
  if (get.norm == TRUE){
    array.long$signal.norm <- 2 ^ array.long$signal
  }
  return(array.long)
}
