# PlotGeneAcrossTissues.R

PlotGeneAcrossTissues <- function(dat, jtitle, convert.linear = FALSE, make.pretty = FALSE, jxlab="CT", do.facet.wrap=TRUE, by.linetype=FALSE, sort.by.mean = FALSE){
  library(ggplot2)
  if (sort.by.mean){
    dat.mean <- dat %>%
      group_by(tissue) %>%
      summarise(exprs = mean(exprs)) %>%
      arrange(desc(exprs))
    dat$tissue <- factor(as.character(dat$tissue), levels = dat.mean$tissue)
  }
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "Log2 mRNA Abundance"
  }
  n.exper = length(as.character(unique(dat$experiment)))
  
  if (n.exper > 1){
    m <- ggplot(dat, aes(x = time, y = exprs,
                         group = experiment, 
                         colour = experiment)) 
  } else {
    m <- ggplot(dat, aes(x = time, y = exprs)) 
  }
    m <- m + 
      ggtitle(jtitle) + 
      ylab(label = jylab) + 
      xlab(label = jxlab)
    if (do.facet.wrap){
      m <- m + geom_point() + geom_line() + facet_wrap(~tissue)
    } else {
      if (!by.linetype){
        m <- m + geom_point(data = dat, aes(x = time, y = exprs, group = tissue, colour = tissue)) + geom_line(data = dat, aes(x = time, y = exprs, group = tissue, colour = tissue))
      } else {
        m <- m + geom_point(data = dat, aes(x = time, y = exprs, group = tissue, linetype = tissue)) + geom_line(data = dat, aes(x = time, y = exprs, group = tissue, linetype = tissue))
      }
    }
  if (make.pretty){
   m <- m + theme_bw() + theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               aspect.ratio = 1)
  }
  return(m)
}

ConvertToSingleDay <- function(dat, time.cname = "time", exprs.cname = "exprs", 
                               geno.cname = "geno", exper.cname = "experiment",
                               tissue.cname = "tissue", gene.cname = "gene", transcript.cname = NA, use.se=FALSE){
  library(lazyeval)
  if (any(dat$time < 0)){
    warning("Function cannot handle negative times")
  }
  mutate_call = lazyeval::interp(~jtime - 24 * floor( jtime / 24 ), jtime = as.name(time.cname))
  dat.sub <- dat %>%
    filter_(paste(time.cname, ">=", 24)) %>%
    mutate_(.dots = setNames(list(mutate_call), time.cname))
  
  if (!use.se){
    sem.cname <- "sem"
    sem_call = lazyeval::interp(~sd( exprs ), exprs = as.name(exprs.cname))
  } else {
    se.cname <- "se"
    sem.cname <- "se"
    sem_call = lazyeval::interp(~mean( se ), exprs = as.name(se.cname))
  }
  signal_call = lazyeval::interp(~mean ( exprs ), exprs = as.name(exprs.cname))
 
  if (is.na(transcript.cname)){
    dat.new <- bind_rows(dat.sub, subset(dat, time < 24)) %>%
      group_by_(gene.cname, geno.cname, exper.cname, tissue.cname, time.cname) %>%
      summarise_(.dots = setNames(list(sem_call, signal_call), c(sem.cname, exprs.cname)))
  } else {
    dat.new <- bind_rows(dat.sub, subset(dat, time < 24)) %>%
      group_by_(gene.cname, geno.cname, exper.cname, tissue.cname, time.cname, transcript.cname) %>%
      summarise_(.dots = setNames(list(sem_call, signal_call), c(sem.cname, exprs.cname)))
  }
  return(dat.new)
}

PlotGeneTissuesWTKO <- function(dat, timelabel="ZT", jtitle="", split.by="geno", ncols = 2, center = FALSE, convert.linear = FALSE, jsize = 24, single.day=FALSE, pretty.geno.names=FALSE){
  if (pretty.geno.names){
    dat$geno <- gsub("BmalKO", "Bmal1 KO", dat$geno)
    dat$geno <- gsub("SV129", "WT", dat$geno)
    dat$geno <- factor(dat$geno, levels = c("WT", "Bmal1 KO"))
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1
    jylab <- "TPM"
  } else {
    jylab <- "Log2 mRNA Abundance"
  }
  # split by geno or tissue
  if (center){
    dat <- dat %>%
      group_by(tissue, geno) %>%
      mutate(exprs = scale(exprs, center = TRUE, scale = FALSE))
  }
  if (single.day){
    dat <- ConvertToSingleDay(dat)
  }
  m <- ggplot(dat, aes(x = time, colour = tissue, linetype = geno, y = exprs)) + 
    geom_point() + geom_line() + xlab(timelabel) + ylab(jylab) + 
    theme_bw(jsize) + ggtitle(jtitle) + 
    theme(aspect.ratio = 1, legend.position = "bottom")
  if (split.by == "geno"){
    m <- m + facet_wrap(~geno, ncol = ncols)
  } else if (split.by == "tissue"){
    m <- m + facet_wrap(~tissue, ncol = ncols)
  } else {
    warning("Split by must be geno or tissue")
  }
  if (single.day){
    m <- m + geom_errorbar(data = dat, aes(ymin=exprs-sem, ymax=exprs+sem, colour = tissue), linetype = "solid", size=0.5, width=0.5)
  }
  return(m)
}


PlotGeneAcrossTissuesRnaseq <- function(dat, jtitle, convert.linear = FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "Log2 mRNA Abundance"
  }
  
  m <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab) +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
    theme(axis.text.x=element_text(angle=90,vjust = 0)) + theme_bw(24)
}

PlotTpmAcrossTissues <- function(dat, jtitle, log2.transform=FALSE, transcript_id = "transcript_id"){
  library(ggplot2)
  
  if (missing(jtitle)){
    jgene <- unique(dat$gene_name)
    jtranscript <- unique(dat$transcript_id)
    jtitle = paste(jgene, jtranscript)
  }
  if (log2.transform == FALSE){
    m <- ggplot(dat, aes_string(x = "time", y = "tpm", group = transcript_id, colour = transcript_id, shape = transcript_id)) 
    jylab <- "mRNA Abundance (normal scale)"
  } else {
    dat$log2tpm <- log2(dat$tpm + 0.01)
    m <- ggplot(dat, aes_string(x = "time", y = "log2tpm", group = transcript_id, colour = transcript_id, shape = transcript_id))
    jylab <- "Log2 mRNA Abundance"
  }
  m <- m + theme(legend.position="bottom") +
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab)
  return(m)
}

PlotTpmAcrossTissuesWTKO <- function(dat, jtitle, log2.transform=FALSE, transcript_id = "transcript_id", geno_id = "geno", tissue_id = "tissue", tstart = 0, tend = 48, jsize = 3, 
                                     single.day=FALSE, 
                                     pretty.geno.names=FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jgene <- unique(dat$gene_name)
    jtranscript <- unique(dat$transcript_id)
    jtitle = paste(jgene, jtranscript)
  }
  if (pretty.geno.names){
    dat$geno <- gsub("BmalKO", "Bmal1 KO", dat$geno)
    dat$geno <- gsub("SV129", "WT", dat$geno)
    dat$geno <- factor(dat$geno, levels = c("WT", "Bmal1 KO"))
  }
  if (log2.transform == FALSE){
    ystr <- "tpm"
    jylab <- "mRNA Abundance (normal scale)"
  } else {
    dat$log2tpm <- log2(dat$tpm + 0.01)
    ystr <- "log2tpm"
    jylab <- "Log2 mRNA Abundance"
  }
  if (single.day){
    dat <- ConvertToSingleDay(dat, transcript.cname = transcript_id, exprs.cname = ystr)
  }
  m <- ggplot(dat, aes_string(x = "time", y = ystr, group = geno_id, linetype = geno_id, colour = tissue_id, shape = transcript_id))
  m <- m + theme(legend.position="bottom") +
    geom_point(size = jsize) + 
    geom_line() + 
    facet_grid(transcript ~ tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab) + 
    theme_bw() + 
    scale_x_continuous(limits = c(tstart, tend), breaks = seq(tstart, tend, 12)) + 
    theme(legend.position = "bottom", aspect.ratio = 1) + 
    scale_shape_manual(values=c(15, 17, seq(0, 14)))
  if (single.day){
    sem.str <- "sem"
    m <- m + geom_errorbar(data = dat, aes_string(ymin=paste(ystr, sem.str, sep = "-"), 
                                                  ymax=paste(ystr, sem.str, sep = "+"), 
                                                  colour = tissue_id), linetype = "solid", size=0.5, width=0.5)
  }
  return(m)
}

PlotRnaseqAcrossTissues <- function(dat, jtitle){
  p <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = "Log2 mRNA Abundance") +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
    theme_bw(24) + 
    theme(axis.text.x=element_text(angle=90,vjust = 0))
  return(p)
}

Center <- function(x){
  # Running scale on dplyr doesnt work, try it manually
  return(x - mean(x))
}

PlotGeneNormalized <- function(dat, jtitle){
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  dat.norm <- dat %>%
    group_by(tissue) %>%
    mutate(exprs.scaled = Center(exprs))
  p <- ggplot(dat.norm, aes(x = time, y = exprs.scaled, group = tissue, colour = tissue, fill = tissue)) + 
    geom_line() +
    geom_point() + 
    ylab("Centered mRNA expression") +
    xlab("CT")
  return(p)
}

PlotEncodeRnaseq <- function(dat, jtitle, sort.by.tissue = TRUE, by.var = "tpm"){
  source("~/projects/tissue-specificity/scripts/functions/SortByTissue.R")
  dat <- SortByTissue(dat, by.var = by.var)
  
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  p <- ggplot(dat, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jtitle) + xlab("") + ylab("Gene expression (TPM)") + 
    theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
  return(p)
}

CalculatePeriodogramLong <- function(dat, jexperiment = "array", time.interval = NULL, remove.inf = TRUE){
  # interval between sampling points, in hours
  if (jexperiment == "array"){
    interval <- 2  # hrs
  } else if (jexperiment == "rnaseq"){
    interval <- 6  # hrs
  } else {
    interval <- jinterval
  }
  exprs <- dat$exprs
  p <- CalculatePeriodogram(exprs)
  periods <- signif(interval / p$freq, digits = 3)
  dat.var.s <- data.frame(periodogram = p$p.scaled, period = periods)
  dat.var.s$period <- factor(dat.var.s$period, 
                              levels = sort(unique(dat.var.s$period), decreasing = TRUE))
  if (remove.inf){
    dat.var.s <- subset(dat.var.s, period != Inf)
  }
  return(dat.var.s)
}

PlotPeriodogramLong <- function(dat, jexperiment = "array", time.interval = NULL, jtitle = "Plot title"){
  # expect dat to be by gene
  # Plot periodogram of the gene
  source("~/projects/tissue-specificity/scripts/functions/FourierFunctions.R")
  n.exper <- length(as.character(unique(dat$experiment)))
  if (n.exper > 1){
    dat <- subset(dat, experiment == jexperiment)
  }
  
  dat.periodogram <- dat %>%
    group_by(gene, tissue) %>%
    do(CalculatePeriodogramLong(., jexperiment, time.interval))
  ggplot(dat.periodogram, aes(x = period, y = periodogram)) + facet_wrap(~tissue) +  geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1) + 
    xlab("Periodogram") + ylab("Period [h]")
}
