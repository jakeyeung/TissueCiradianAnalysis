# Functions for handling Atger-Kidney from Nestle and Jerome

RemoveLowExprsPseudoShortGenes <- function(dat.long, ggbiotype="protein_coding", gglength = 250, jcutoff = 1, show.plot=FALSE){
  dat.sub <- subset(dat.long, !is.na(gene))
  
  # Remove pseudogenes nd short genes
  genes <- unique(as.character(dat.sub$gene))
  genes.biotype <- AnnotatePseudogenes(genes, return.original = FALSE)
  genes.length <- AnnotateTranscriptLength(genes, return.original = FALSE)
  biotype.hash <- hash(genes, genes.biotype)
  length.hash <- hash(genes, genes.length)
  dat.sub$gbiotype <- sapply(as.character(dat.sub$gene), function(g) biotype.hash[[g]])
  dat.sub$glength <- sapply(as.character(dat.sub$gene), function(g) length.hash[[g]])
  dat.sub <- subset(dat.sub, gbiotype == ggbiotype & glength > gglength)
  
  # remove lowly expressed genes
  dat.mean <- dat.sub %>%
    group_by(gene) %>%
    summarise(exprs.max = quantile(exprs, probs = 0.9))
  
  if (show.plot){
    plot(density(dat.mean$exprs.max))
    jcutoff <- 1
    abline(v=jcutoff)
  }
  
  genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)
  
  dat.sub <- subset(dat.sub, !gene %in% genes.cut)
  return(dat.sub)
}

CollapseTissueGeno <- function(dat.long, keep.tissue.col=FALSE){
  # Collapse tissue into tissuegeno, so you can use functions that expect "tissues" as conditions in colname (hack)
  if (keep.tissue.col){
    dat.long$tissue.old <- dat.long$tissue
  }
  dat.long$tissue <- factor(paste(dat.long$tissue, dat.long$geno, sep = "_"), levels = c("Liver_SV129", "Kidney_SV129", "Liver_BmalKO", "Kidney_BmalKO"))
  return(dat.long)
}

SameTimepointsLivKid <- function(dat.long){
  # Make Liver and Kidney same number of timepoints (tissue column). Should work for collapse or not collapse
  conds <- unique(as.character(dat.long$tissue))
  conds.lst = list()
  for (cond in conds){
    conds.lst[[cond]] <- unique(subset(dat.long, tissue == cond)$time)
  }
  common.times <- Reduce(intersect, conds.lst)
  return(subset(dat.long, time %in% common.times))
}

StaggeredTimepointsLivKid <- function(dat.long){
  # Take consecutive timepoints for liver and kidney (keep all liver, take kidney not in liver)
  # Make Liver and Kidney same number of timepoints (tissue column). Should work for collapse or not collapse
  conds <- unique(as.character(dat.long$tissue))
  conds.lst <- list()
  minsize <- Inf
  mintimes <- c()
  for (cond in conds){
    conds.lst[[cond]] <- unique(subset(dat.long, tissue == cond)$time)
    if (length(conds.lst[[cond]]) < minsize){
      minsize <- length(conds.lst[[cond]])
      mintimes <- conds.lst[[cond]]
    } 
  }
  dat.sub <- dat.long %>%
    group_by(tissue) %>%
    do(SetdiffTimepoints(., mintimes))
  return(dat.sub)
}

SetdiffTimepoints <- function(dat.tiss, mintimes){
  # if dat.tiss matches mintimes, do nothing
  # if dat.tiss larger than mintimes, take setdiff
  times <- unique(dat.tiss$time)
  if (all(times == mintimes)){
    return(dat.tiss)
  } else {
    times.new <- setdiff(times, mintimes)
    return(subset(dat.tiss, time %in% times.new))
  }
}

RemoveLowlyExpressedGenes <- function(dat.long, jquantile = 0.9, jcutoff = 1, show.plot=FALSE){
  dat.mean <- dat.long %>%
    group_by(gene) %>%
    summarise(exprs.max = quantile(exprs, probs = 0.9))
  
  if (show.plot){
    plot(density(dat.mean$exprs.max))
    abline(v=jcutoff)
  }
  genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)
  dat.long <- subset(dat.long, ! gene %in% genes.cut)
  return(dat.long)
}

# Functions for handling DAT_I_I liver kidney Cedric
TissueFromCname <- function(cname){
  # "X0_1 -> Kidney"
  # "X0_1_liver -> Liver"
  jsplit <- strsplit(cname, "_")[[1]]
  if (jsplit[[length(jsplit)]] == "liver"){
    tiss <- "Liver"
  } else {
    tiss <- "Kidney"
  }
  return(tiss)
} 

TimeFromCname <- function(cname, rm.first.char=TRUE){
  # X0_1 -> 0 + (1 - 1) * 24
  # X22_2 -> 22 + (2 - 1) * 24
  if (rm.first.char){
    cname <- substr(cname, start = 2, stop = nchar(cname))
  }
  time.base <- as.numeric(strsplit(cname, "_")[[1]][[1]])
  time.multiplier <- as.numeric(strsplit(cname, "_")[[1]][[2]]) - 1
  time <- time.base + time.multiplier * 24
  return(time)
}

LoadLivKid <- function(inf){
  if (missing(inf)){
    inf <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/DAT_l_l.txt"
  }
  dat.mat <- read.table(inf, header = TRUE)
  genes <- dat.mat$name
  exprs <- subset(dat.mat, select = c(-name))
  tissues <- rep(sapply(colnames(exprs),TissueFromCname, USE.NAMES = FALSE), each = nrow(exprs))
  tissues.uniq <- unique(tissues)
  times <- rep(sapply(colnames(exprs), TimeFromCname, USE.NAMES = FALSE), each = nrow(exprs))
  
  # make long
  dat <- data.frame(gene = dat.mat$name, 
                    tissue = tissues,
                    time = times,
                    exprs = unlist(exprs),
                    experiment = "rnaseq")
}