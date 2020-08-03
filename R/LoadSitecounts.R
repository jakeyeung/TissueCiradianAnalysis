LoadSitecounts <- function(N.path, N.promoterpath, gene_ids=TRUE){
  if (missing(N.path)){
    if (gene_ids){
      N.path <- "data/sitecounts/motevo/sitecount_matrix_geneids"
    } else {
      N.path <- "data/sitecounts/motevo/sitecount_matrix"
    }
  } 
  if (missing(N.promoterpath)){
    N.promoterpath <- "data/sitecounts/motevo/mara_promoters_gene_name_association.bed"
  }
  N <- read.table(N.path, header = TRUE, sep = '\t')
  rnames <- make.names(N[, 1], unique = TRUE)
  rownames(N) <- rnames
  N <- N[, 2:ncol(N)]
  
  N.promoter <- read.table(N.promoterpath, header = FALSE, sep = '\t')
  N.promoter <- data.frame(Gene.ID=make.names(N.promoter$V7, unique = FALSE), maraid=N.promoter$V4, saeedid=N.promoter$V10)
  rownames(N.promoter) <- make.names(N.promoter$Gene.ID, unique = TRUE)
  return(list(N=N, N.promoter=N.promoter))
}

LoadSitecountsPromotersLong <- function(N.path){
  if (missing(N.path)){
    N.path <- "data/sitecounts/motevo/sitecount_matrix"
  } 
  N <- read.table(N.path, header = TRUE, sep = '\t')
  promoterid <- N$Promoter
  rnames <- make.names(N[, 1], unique = TRUE)
  rownames(N) <- rnames
  N <- N[, 2:ncol(N)]
  N.long <- data.frame(promoterid = promoterid, 
                       motif = rep(colnames(N), each = nrow(N)),
                       sitecount = unlist(N))
  return(N.long)
}

LoadEnsemblToPromoter <- function(annot.path){
  if (missing(annot.path)){
    annot.path <- "data/sitecounts/motevo/mara_promoters_gene_name_association.ensembl_first_exon_overlap.cut.bed"
  }
  annot <- read.table(annot.path, header = FALSE, sep = '\t')
  annot <- data.frame(Gene.ID=annot$V2, maraid=annot$V1, saeedid=annot$V3, ensemblid=annot$V4)
  transcriptid <- sapply(annot$ensemblid, function(s){
    # Convert gene_name=Lypla1;transcript_id=ENSMUST00000134384 to ENSMUST00000134384
    s <- as.character(s)
    s.transcriptstr <- strsplit(s, ';')[[1]][[2]]
    s.transcriptid <- strsplit(s.transcriptstr, '=')[[1]][[2]]
  })
  annot$ensemblid = transcriptid
  return(annot)
}

LoadSitecountsEncode <- function(inpath, tissue, with.ensemblid = TRUE){
  jdf <- read.table(inpath, header = TRUE, sep = "\t")
  rnames <- as.character(jdf$ensemblid)
  genes <- sapply(rnames, function(s) strsplit(s, ";")[[1]][[1]])
  if (with.ensemblid == TRUE){
    ensemblids <- sapply(rnames, function(s) strsplit(s, ";")[[1]][[2]])
  }
  rownames(jdf) <- rnames
  jdf$ensemblid <- NULL
  if (with.ensemblid == TRUE){
    jdf.long <- data.frame(gene = rep(genes, length(colnames(jdf))),
                           ensemblid = rep(ensemblids, length(colnames(jdf))),
                           motif = rep(colnames(jdf), each = length(rep(rownames(jdf)))),
                           motevo.value = unlist(jdf))
  } else {    
    jdf.long <- data.frame(gene = rep(genes, length(colnames(jdf))),
                           motif = rep(colnames(jdf), each = length(rep(rownames(jdf)))),
                           motevo.value = unlist(jdf))
  }
  jdf.long$tissue <- as.factor(rep(tissue, nrow(jdf.long)))
  return(jdf.long)
}

LoadSitecountsEncodeAll <- function(maindir, tissues, suffix="sitecounts.matrix", with.ensemblid = TRUE, rename.tissues = FALSE){
  # Load ENCODE sitecounts
  if (missing(maindir)){
    maindir <- "data/sitecounts/motevo/encode_sitecount_matrix"
  }
  if (missing(tissues)){
    tissues <- c("Liver", "Cere", "Mus", "Kidney", "Heart", "Lung")  # matches with Hogenesch names
    names(tissues) <- c("Liver", "Cerebellum", "SkeletalMuscle", "Kidney", "Heart", "Lung")  # matches with Hogenesch data
  } else {
    # make names(tissues) equal tissues, solves a bug inelegantly
    names(tissues) <- tissues
  }
  # suffix <- "sitecounts.matrix"
  fnames <- paste(names(tissues), suffix, sep=".")
  paths <- paste(maindir, fnames, sep="/")
  motevo.df.list <- list()
  for (p in paths){
    if (rename.tissues){
      tissue <- tissues[strsplit(basename(p), "\\.")[[1]][[1]]]
    } else {
      tissue <- strsplit(basename(p), "\\.")[[1]][[1]]
    }
    motevo.df.list[[tissue]] <- LoadSitecountsEncode(p, tissue, with.ensemblid)
  }
  motevo.long <- do.call("rbind", motevo.df.list)
}