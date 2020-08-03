FixColname <- function(tiss, time, exper){
  # Make MARA compatible colnames e.g., Adr22_array
  fixed_cname <- paste0(tiss, time, "_", exper)
  return(fixed_cname)
}

FixColnamesForMARA <- function(m){
  # Make colnames compatible with MARA
  # gene colname is "Gene.ID" (first column)
  # samples are TissueTime_arrayORrnaseq (e.g., Adr22_array)
  # 
  # Rearrrange colnames also to make put all arrays, THEN all rna-seq
  cnames <- colnames(m)
  genestr <- "Gene.ID"
  
  # Fix sample names (ignore first column which is gene)
  tissues <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[1]])
  times <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[2]])
  experiment <- sapply(cnames[2:length(cnames)], function(s) strsplit(s, "_")[[1]][[3]])
  
  sampnames <- mapply(FixColname, tissues, times, experiment)
  
  cnames.new <- c(genestr, sampnames)
  colnames(m) <- cnames.new
  
  # reorder colnames now
  array.samps <- cnames.new[grep("_array", cnames.new)]
  array.rnaseq <- cnames.new[grep("_rnaseq", cnames.new)]
  
  m <- m[, c(genestr, array.samps, array.rnaseq)]
  return(m)
}

WriteGeneListMat <- function(fits.rhyth.tiss, ka.long, outdir){
  # input: list of rhythmic genes by tissues (gene and tissue colnames), should be only ONE tissue (for dplyr)
  # dataframe of long matrix: (colnames: gene tissue time exprs experiment)
  # 
  # NOTE: expression should be centered across all samples! Otherwise it's a weird MARA model
  tiss <- unique(as.character(fits.rhyth.tiss$tissue))
  if (length(tiss) != 1){
    print('Warning: expect tissue length = 1')
  }
  fname.genelist <- paste(tiss, "genelist", sep = ".")
  fname.mat <- paste(tiss, "mat", sep=".")
  
  genes <- as.character(fits.rhyth.tiss$gene) 
  # write genelist
  sink(file = file.path(outdir, fname.genelist))
  for (g in genes){
    cat(g)
    cat("\n")
  }
  sink()
  
  # write matrix expression
  ka.sub <- subset(ka.long, gene %in% genes & tissue == tiss)
  m <- dcast(data = ka.sub, formula = gene ~ tissue + time + experiment, value.var = "exprs.centered")
  
  m <- FixColnamesForMARA(m)
  
  write.table(m, file = file.path(outdir, fname.mat), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  return(data.frame(NULL))
}

WriteGeneListMat2 <- function(tissues, genelist, dat, valvar, outdir){
  # More general function for WriteGeneListMat
  dat.sub <- subset(dat, gene %in% genelist)
  for (tiss in tissues){
    fname.genelist <- paste(tiss, "genelist", sep = ".")
    fname.mat <- paste(tiss, "mat", sep=".")
    
    # write genelist
    sink(file = file.path(outdir, fname.genelist))
    for (g in genelist){
      cat(g)
      cat("\n")
    }
    sink()
    
    # write matrix expression
    ka.sub <- subset(dat.sub, tissue == tiss)
    m <- dcast(data = ka.sub, formula = gene ~ tissue + time + experiment, value.var = valvar)
    
    m <- FixColnamesForMARA(m)
    
    write.table(m, file = file.path(outdir, fname.mat), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  # write for all
  fname.all <- paste("all", "mat", sep=".")
  m <- dcast(data = dat.sub, formula = gene ~ tissue + time + experiment, value.var = valvar)
  m <- FixColnamesForMARA(m)
  write.table(m, file = file.path(outdir, fname.all), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
}