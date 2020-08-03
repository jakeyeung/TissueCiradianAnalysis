GetSinusoid <- function(act.sub, x, w){
  phase <- unique(act.sub$phase)
  geno.act <- unique(act.sub$geno)
  tissue.act <- unique(act.sub$tissue)
  return(data.frame(time = x, exprs = cos(w * x - w * phase), tissue = tissue.act, geno = geno.act))
}

LoadProteomicsData <- function(inf = "/data/shared/OMICS/RNA-Seq/YeungGR2018/ProcessedData/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.OneGenePerLine.txt", 
                               as.long = TRUE){
  require(reshape2)
  require(dplyr)
  prot <- read.csv(inf, header = TRUE, sep = "\t")
  if (!as.long){
    return(prot)
  }
  # Make long 
  wt.sampnames <- paste("ZT", sprintf("%02d", seq(0, 45, 3)), ".WT", sep = "")
  ko.sampnames <- paste("ZT", sprintf("%02d", seq(0, 18, 6)), ".Bmal.KO", sep = "")
  bmalwt.sampnames <- paste("ZT", sprintf("%02d", seq(0, 18, 6)), ".Bmal.WT", sep = "")
  # fit.sampnames <- c("mean", "amp", "relamp", "phase", "pval", "qv", "amp.12h", "relamp.12h", "phase.12h", "pval.12h", "qv.12h")
  fit.sampnames <- c("mean", "amp", "relamp", "phase", "pval", "qv")
  
  idvars <- c("Gene.names", "Protein.IDs")
  
  prot.long.wt <- melt(prot, id.vars = idvars, measure.vars = wt.sampnames, variable.name = "samp", value.name = "rel.abund")
  prot.long.bmalko <- melt(prot, id.vars = idvars, measure.vars = ko.sampnames, variable.name = "samp", value.name = "rel.abund")
  prot.long.bmalwt <- melt(prot, id.vars = idvars, measure.vars = bmalwt.sampnames, variable.name = "samp", value.name = "rel.abund")
  fit.prot.wt <- subset(prot, select = c(idvars, fit.sampnames))
  # fit.prot.wt <- melt(prot, id.vars = "Gene.names", measure.vars = fit.sampnames)
  
  prot.long.wt$time <- GetTimeFromSamp(as.character(prot.long.wt$samp))
  prot.long.wt$geno <- GetGenoFromSamp(as.character(prot.long.wt$samp))
  
  prot.long.bmalko$time <- GetTimeFromSamp(as.character(prot.long.bmalko$samp))
  prot.long.bmalko$geno <- GetGenoFromSamp(as.character(prot.long.bmalko$samp))
  
  prot.long.bmalwt$time <- GetTimeFromSamp(as.character(prot.long.bmalwt$samp))
  prot.long.bmalwt$geno <- GetGenoFromSamp(as.character(prot.long.bmalwt$samp))
  
  # merge
  prot.long <- rbind(prot.long.wt, prot.long.bmalko, prot.long.bmalwt)
  prot.long$tissue <- "Liver"
  
  # change Gene.names to gene
  colnames(fit.prot.wt)[which(colnames(fit.prot.wt) == "Gene.names")] <- "gene"
  colnames(prot.long)[which(colnames(prot.long) == "Gene.names")] <- "gene"
  return(prot.long)
}

LoadPhosphoData <- function(inf = "/home/shared/nuclear_proteomics/Table_SXX_nuclear_phospho_all.OneGenePerLine.txt", as.long = TRUE){
  require(reshape2)
  require(dplyr)
  if (!as.long){
    prot <- read.csv(inf, header=TRUE, sep = "")
    return(prot)
  } else {
    return(LoadProteomicsData(inf, as.long = TRUE))
  }
}

PlotProteomics <- function(prot.long, jtitle = ""){
  # plot proteomics
  if (all(is.na(prot.long$rel.abund))) return(NA)
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  g <- ggplot(prot.long, aes(x = time, y = rel.abund, colour = geno, group = geno)) + geom_point() + geom_line() + theme_bw() + scale_color_manual(values = cbPalette)
  g <- g + ggtitle(jtitle)
  return(g)
}

GetTimeFromSamp <- Vectorize(function(samp){
  time <- strsplit(samp, "\\.")[[1]][[1]]
  # remove ZT
  time <- gsub("ZT", "", time)
  return(as.numeric(time))
}, vectorize.args = "samp")

GetGenoFromSamp <- Vectorize(function(samp){
  samp.split <- strsplit(samp, "\\.")[[1]]
  if (length(samp.split) == 2){
    geno <- samp.split[[2]]
  } else if (length(samp.split) == 3){
    geno <- paste(samp.split[2:3], collapse = "")
  }
  return(geno)
}, vectorize.args = "samp")

AddColname <- function(dat, cname, cvalue){
  if (nrow(dat) > 0){
    dat[[cname]] <- cvalue
    return(dat)
  } else {
    return(dat)
  }
}

SubsetGenoSignalTimeType <- function(dat){
  if (nrow(dat) > 0){
    return(as.data.frame(subset(dat, select = c(geno.std, signal, time, type, tissue))))
  } else {
    return(as.data.frame(dat))
  }
}

ScaleSignal <- function(dat, cname, jcenter = NA, jscale = NA){
  if (is.na(jcenter) & is.na(jscale)){
    return(scale(dat[[cname]], center = TRUE, scale = TRUE))
  } else {
    # make center jcenter, variance jscale
    scaled.dat <- scale(dat[[cname]], center = TRUE, scale = TRUE) * sqrt(jscale)
    scaled.dat <- scaled.dat + jcenter
    return(scaled.dat)
  }
}


PlotmRNAActivityProtein <- function(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot, jtiss = "Liver", dotsize = 3, themesize=24, 
                                    n.facetrows = 1, wt.prot = "Cry", line.for.protein=FALSE, act.in.sine.cos=FALSE, single.day=FALSE, by.color=FALSE){
  # plotting protein optional: can use this to plot just dat and act if prot.long is missing and gene.prot is ""
  if (!"geno" %in% colnames(act.long)){
    act.long$geno = tryCatch({
      sapply(as.character(act.long$tissue), function(tissgeno) strsplit(tissgeno, "_")[[1]][[2]])
    }, error = function(e) {
      "WT"
    })
    # sapply(as.character(act.long$tissue), function(tissgeno) strsplit(tissgeno, "_")[[1]][[2]])
  }
  act.long$tissue <- sapply(as.character(act.long$tissue), function(tissgeno) strsplit(tissgeno, "_")[[1]][[1]])
  if (jtiss == "Liver"){
    dat.sub <- subset(dat.long, gene == gene.dat & tissue == "Liver")
    act.sub <- subset(act.long, gene == gene.act & tissue %in% c("Liver"))
  } else if (jtiss == "Kidney"){
    dat.sub <- subset(dat.long, gene == gene.dat & tissue == "Kidney")
    act.sub <- subset(act.long, gene == gene.act & tissue %in% c("Kidney"))
  } else if (jtiss == "both"){
    dat.sub <- subset(dat.long, gene == gene.dat)
    act.sub <- subset(act.long, gene == gene.act)
  } else {
    warning("jtiss must be Liver or Kidney or both")
  }
  if (!(missing(prot.long) | all(is.na(prot.long)))){  # if not missing or not na
    if (wt.prot == "Cry"){
      prot.sub <- subset(prot.long, gene == gene.prot & geno %in% c("WT", "BmalKO"))
    } else if (wt.prot == "Bmal"){
      prot.sub <- subset(prot.long, gene == gene.prot & geno %in% c("BmalWT", "BmalKO"))
      prot.sub$geno <- as.factor(gsub("BmalWT", "WT", prot.sub$geno))
    } else {
      warning("wt.prot should be Cry or Bmal")
    }
  }
  # if (any(c(nrow(dat.sub) == 0, nrow(act.sub) == 0))){
  #   warning("Warning, empty dataframe. Returning empty")
  #   return(NULL)
  # }
  if (act.in.sine.cos){
    # Generate based on amplitude and phase
    x <- seq(0, 48)
    w <- 2 * pi / 24
    act.sub <- act.sub %>%
      group_by(geno, tissue) %>%
      do(GetSinusoid(., x, w))
    # print(data.frame(dat.sub))
    # phase <- unique(act.sub$phase)
    # geno.act <- unique(act.sub$geno)
    # tissue.act <- unique(act.sub$tissue)
    # act.sub <- data.frame(time = x, exprs = cos(w * x - w * phase), tissue = tissue.act, geno = geno.act)
  }
  

  
  # scale data, merge together, then plot in one figure
  dat.sub$signal <- ScaleSignal(dat.sub, "exprs")
  act.sub$signal <- ScaleSignal(act.sub, "exprs")
  if (!(missing(prot.long) | all(is.na(prot.long)))){
    if (jtiss != "both"){
      prot.sub$signal <- ScaleSignal(prot.sub, "rel.abund")
    } else {
      # center and scale based on mean and variance from Liver data of dat.sub
      jcenter <- mean(subset(dat.sub, tissue == "Liver")$signal)
      jscale <- var(subset(dat.sub, tissue == "Liver")$signal)
      prot.sub$signal <- ScaleSignal(prot.sub, "rel.abund", jcenter, jscale)
    }
  }
  
  # make two factors: WT and BmalKO
  dat.sub$geno.std <- gsub("SV129", "WT", as.character(dat.sub$geno))
  act.sub$geno.std <- gsub("SV129", "WT", as.character(act.sub$geno))
  
  # if geno is named Bmal rename to BmalKO. But we renamed to BmalWT and BmalKO
  if (!(missing(prot.long) | all(is.na(prot.long)))){
    prot.sub$geno.std <- prot.sub$geno
    # prot.sub$geno.std <- gsub("Bmal", "BmalKO", as.character(prot.sub$geno))
  }
  # add tissue info to geno.std
  if (jtiss == "both"){
    dat.sub$geno.std <- as.factor(paste(dat.sub$tissue, dat.sub$geno.std, sep = "_"))
    act.sub$geno.std <- as.factor(paste(act.sub$tissue, act.sub$geno.std, sep = "_"))
    if (!(missing(prot.long) | all(is.na(prot.long)))){
      prot.sub$geno.std <- as.factor(paste(prot.sub$tissue, prot.sub$geno.std, sep = "_"))
    }
  } 
  
  dat.sub <- AddColname(dat.sub, "type", "mRNA_Accum")
  act.sub <- AddColname(act.sub, "type", "TF_Activity")
  if (!(missing(prot.long) | all(is.na(prot.long)))){
    prot.sub <- AddColname(prot.sub, "type", "Nuclear_Prot_Accum")
  }

  dat.sub2 <- SubsetGenoSignalTimeType(dat.sub)
  act.sub2 <- SubsetGenoSignalTimeType(act.sub)
  # print(dat.sub)
  # print(dat.sub2)
  # print(act.sub)
  # print(act.sub2)
  if (!(missing(prot.long) | all(is.na(prot.long)))){
    prot.sub2 <- SubsetGenoSignalTimeType(prot.sub)
    factor.levels <- c("mRNA_Accum", "Nuclear_Prot_Accum", "TF_Activity")
    # ltypes <- c("solid", "twodash", "dotted")
    ltypes <- c("dotted", "twodash", "solid")
    # jshapes <- c(17, 4, 15)
    jshapes <- c(32, 17, 32)
    jsizes <- c(3, 5, 5)
    if (!by.color){
      jcols <- rep("black", 3)
    } else {
      # colorblind
      jcols <- c("#999999", "#E69F00", "#56B4E9")
      ltypes <- rep("solid", length(ltypes))
    }
  } else {
    # can rbind a null dataframe no problem
    prot.sub2 <- data.frame(NULL)
    factor.levels <- c("mRNA_Accum", "TF_Activity")
    # ltypes <- c("solid", "dotted")
    ltypes <- c("dotted", "solid")
    # jshapes <- c(17, 15)
    # jshapes <- c(1, 17)
    jshapes <- c(32, 32)
    jsizes <- c(2, 5)
    if (!by.color){
      jcols <- rep("black", 2)
    } else {
      # colorblind
      jcols <- c("#999999", "#56B4E9")
      ltypes <- rep("solid", length(ltypes))
    }
  }
  merged.dat <- rbind(dat.sub2,
                      act.sub2,
                      prot.sub2)
  
  if (jtiss != "both"){
    merged.dat$geno.std <- factor(as.character(merged.dat$geno.std), levels = c("WT", "BmalKO"))
  } else {
    merged.dat$geno.std <- factor(as.character(merged.dat$geno.std), levels = c("Liver_WT", "Liver_BmalKO", "Kidney_WT", "Kidney_BmalKO"))
  }
  merged.dat$type <- factor(merged.dat$type, levels = factor.levels)
  
  if (single.day){
    # 0 to 23, everything else gets shifted by 24 hours
    merged.dat.sub <- subset(merged.dat, time >= 24)
    merged.dat.sub$time <- merged.dat.sub$time - 24 * floor(merged.dat.sub$time / 24)
    merged.dat <- bind_rows(merged.dat.sub, subset(merged.dat, time < 24)) %>%
      group_by(geno.std, type, tissue, time) %>%
      summarise(sem = sd(signal),
                signal = mean(signal))
  }
  
  jtitle <- paste(unique(c(gene.dat, gene.prot, gene.act)), collapse = " ")
  if (jtiss != "both"){
    m <- ggplot(merged.dat, aes(x = time, y = signal, linetype = type, colour = type, shape = type, group = type, size = type))
    # m <- ggplot(merged.dat, aes(x = time, y = signal, linetype = type, shape = type, group = type, size = type))
  } else {
    m <- ggplot(merged.dat, aes(x = time, y = signal, linetype = type, shape = type, group = type, colour = tissue, size = type))
  }
  if (line.for.protein){
    m <- m + geom_line(data = merged.dat, size = 1)
  } else {
    m <- m + geom_line(data = subset(merged.dat, type != "Nuclear_Prot_Accum"), size = 1)
  }
  
  if (single.day){
    jbreaks <- seq(0, 24, 6)
  } else {
    jbreaks <- seq(0, 48, 12)
  }
  m <- m + 
    # geom_line(data = subset(merged.dat, type != "Nuclear_Prot_Accum"), size = 1) +
    # geom_point(size = dotsize) + 
    geom_point() + 
    facet_wrap(~geno.std, nrow = n.facetrows) + 
    theme_bw(themesize) + 
    xlab("ZT") + ylab("Abundance or Motif Activity\n(scaled)") +
    scale_linetype_manual(values = ltypes, drop=FALSE) +
    scale_shape_manual(values = jshapes, drop=FALSE) +
    scale_size_manual(values = jsizes, drop=FALSE) + 
    scale_colour_manual(values = jcols, drop=FALSE) + 
    theme(legend.position = "bottom", aspect.ratio = 1) +
    ggtitle(jtitle) + 
    scale_x_continuous(breaks = jbreaks)
  if (single.day){
    if (act.in.sine.cos){
      m <- m + geom_errorbar(data = subset(merged.dat, type != "TF_Activity"), aes(ymin=signal-sem, ymax=signal+sem), linetype = "solid", size=0.5, width=0.5, position=position_dodge(0.05))
    } else {
      m <- m + geom_errorbar(data = merged.dat, aes(ymin=signal-sem, ymax=signal+sem), linetype = "solid", size=0.5, width=0.5, position=position_dodge(0.05))
    }
  }
  return(m)
}

