NrhythToStr <- function(n.rhyth, max.n.rhyth = 8){
  if (n.rhyth < max.n.rhyth){
    return(as.character(n.rhyth))
  } else {
    return(paste0(max.n.rhyth, "+"))
  }
}

CountModels <- function(fits.best.sub){
  fits.sum <- fits.best.sub %>%
    group_by(model) %>%
    summarise(count = length(model))
  fits.sum <- fits.sum[order(fits.sum$count, decreasing = TRUE), ]
  fits.best.sub$model <- factor(as.character(fits.best.sub$model), levels = fits.sum$model)
  fits.sum <- OrderDecreasing(fits.sum, jfactor = "model", jval = "count")
  return(fits.sum)
}

PlotPolarHistogram <- function(fits.sub, countstring="Count", phase.cname = "phase.avg"){
  # countstring="Count" or "NormCount"
  fits.sub.sum <- CountModels(fits.sub)
  # bin manually
  # 24 and 0 are the same bin
  fits.sub$bin <- sapply(fits.sub[[phase.cname]], function(p){
    if (p > 0.5){
      phase.bin <- round(p, digits = 0)
    } else {
      phase.bin <- 24
    }
    return(phase.bin)
  })
  fits.sub.forhist <- fits.sub %>%
    group_by(model, bin) %>%
    summarise(Count = length(bin))
  fits.sub.forhist$bin <- fits.sub.forhist$bin - 0.5
  
  # add label to indicate number of genes
  n.genes.hash <- hash(as.character(fits.sub.sum$model), fits.sub.sum$count)
  # rename tissues
  fits.sub.forhist$model.n <- as.factor(sapply(as.character(fits.sub.forhist$model), function(m) paste0(m, " N=", n.genes.hash[[m]])))
  tiss.ordered <- sapply(as.character(fits.sub.sum$model), function(m) paste0(m, " N=", n.genes.hash[[m]]))
  fits.sub.forhist$model.n <- factor(as.character(fits.sub.forhist$model.n), levels = tiss.ordered)
  # keep original labels
  fits.sub.forhist$tissue <- sapply(as.character(fits.sub.forhist$model.n), function(jlab){
    strsplit(jlab, " ")[[1]][[1]]
  })
  fits.sub.forhist$tissue <- factor(as.character(fits.sub.forhist$tissue), levels = names(tiss.ordered))
  fits.sub.forhist <- fits.sub.forhist %>%
    group_by(model) %>%
    mutate(NormCount = Count / sum(Count))
  
  # normalized by sum
  m2 <- ggplot(fits.sub.forhist, aes_string(x = "bin", y = countstring)) +
    geom_bar(stat = "identity", width = 1) + 
    facet_wrap(~tissue) + 
    scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
    expand_limits(y = 0) +
    expand_limits(x = c(0, 24)) +
    xlab("Phase (h)") + 
    ylab("Count") +
    theme_bw() + 
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 6)) + 
    coord_polar(theta = "x")
  return(m2)
}

FilterKidneyLiverGenes <- function(param.vec, amp.min){
  # return TRUE if Kidney and Liver amplitudes are above amp.min
  liv.amp <- 0
  kid.amp <- 0  # inits
  livkid.amp <- 0
  
  liv.indx <- which(names(param.vec) == "Liver.amp")
  kid.indx <- which(names(param.vec) == "Kidney.amp")
  livkid.indx <- which(names(param.vec) == "Kidney,Liver.amp")
  if (length(liv.indx) == 1){
    liv.amp <- param.vec[liv.indx]
  }
  if (length(kid.indx) == 1){
    kid.amp <- param.vec[kid.indx]
  } 
  if (length(livkid.indx) == 1){
    livkid.amp <- param.vec[livkid.indx]
  }   
  if ((liv.amp > amp.min & kid.amp > amp.min) | livkid.amp > amp.min) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

FilterLiverGenes <- function(param.vec, amp.min){
  # return TRUE if Kidney and Liver amplitudes are above amp.min
  liv.amp <- 0
  
  liv.indx <- which(names(param.vec) == "Liver.amp")
  if (length(liv.indx) == 1){
    liv.amp <- param.vec[liv.indx]
  }
  if (liv.amp > amp.min) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

PlotGenesInModel <- function(fits.sub, dat.complex, filt.tiss = c("WFAT")){
  genes <- as.character(fits.sub$gene)
  
  #outobj <- PlotHeatmapNconds(fits, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
  s <- SvdOnComplex(subset(dat.complex, gene %in% genes & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  return(multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout))
}

GetTissuesFromCoefFit <- function(fit.coef.names, tissue.colname = "tissue"){
  # Get tissues involved by grepping everything after tissue.colname
  match.names.i <- grepl(tissue.colname, fit.coef.names) & !grepl(":", fit.coef.names)
  match.names <- fit.coef.names[match.names.i]
  tissues <- sapply(match.names, function(jname) strsplit(jname, split = tissue.colname)[[1]][[2]], USE.NAMES = FALSE)
  return(tissues)
}

ExtraParamsFromFit <- function(fit.coef, tissues){
  # From coefficients of fit, extract parameters
  # colnames should contain tissue, RNA-Seq intercept, phase, amplitude, and BIC
  cnames <- tissues
  # TODO
}

SubsetByMaxBicWeight <- function(dat){
  max.i <- which.max(dat$bicweight)
  return(dat[max.i, ])
}

PlotSvdFromGeneList <- function(dat.complex, gene.list, jcomp = 1, jlabel.n = 25, jvalue.var = "exprs.transformed", constant.amp = 4, show.plot=TRUE){
  # Plot SVD from gene list
  source("/home/yeung/projects/tissue-specificity/scripts/functions/SvdFunctions.R")
  jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
  s.custom <- SvdOnComplex(subset(dat.complex, gene %in% gene.list), value.var = jvalue.var)
  eigens.custom <- GetEigens(s.custom, period = 24, comp = jcomp, label.n = jlabel.n, eigenval = TRUE, adj.mag = TRUE, constant.amp = constant.amp)
  eigens.custom2 <- GetEigens(s.custom, period = 24, comp = (jcomp + 1), label.n = jlabel.n, eigenval = TRUE, adj.mag = TRUE, constant.amp = constant.amp)
  if (show.plot){
    multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)    
  }
  return(eigens.custom)
}

GetTopGenesFromSvd <- function(eigens.obj, top.n = 25){
  return(names(sort(eigens.obj$eigensamp, decreasing = TRUE)[1:top.n]))
}

PCbiplot <- function(PC, jsizes, x="PC1", y="PC2", max.size = 5) {
  if (missing(jsizes)) jsizes <- 1
  # PC being a prcomp object
  dat <- data.frame(obsnames=row.names(PC$x), size = jsizes, PC$x)
  plot <- ggplot(dat, aes_string(x=x, y=y)) + geom_text(alpha=.8, aes(label=obsnames, size = size)) + scale_size(range = c(0, max.size))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  return(plot)
}

PlotTissuePcaFromGenelist <- function(dat.means, genelist, fits.comb){
  # dat.means: gene, tissue, exprs.mean long format
  # genelist: list of genes
  dat.mat.comb <- dcast(subset(dat.means, gene %in% genelist), formula = gene ~ tissue, value.var = "exprs.mean")
  rownames(dat.mat.comb) <- dat.mat.comb$gene; dat.mat.comb$gene <- NULL
  
  # center
  dat.mat.comb <- t(scale(t(dat.mat.comb), center = TRUE, scale = FALSE))
  dat.mat.pca <- prcomp(dat.mat.comb, center = FALSE, scale. = FALSE)
  
  #   biplot(dat.mat.pca)
  jsizes <- fits.comb$bicweight[order(fits.comb$gene)]
  print(PCbiplot(dat.mat.pca, jsizes))
}