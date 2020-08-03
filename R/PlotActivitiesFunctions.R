PlotActivitiesWithSE.dhs <- function(dat, jtitle, sort.by.tissue = TRUE, by.var = "exprs"){
  # for DHS: we have no time, plot tissues on x axis
  source("~/projects/tissue-specificity/scripts/functions/SortByTissue.R")
  dat <- SortByTissue(dat, by.var = by.var)
  
  jgene <- unique(dat$gene)
  if (missing(jtitle)){
    jtitle <- jgene
  }
  ggplot(dat, aes(x = tissue, y = exprs)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    xlab("") +
    ylab("Activity") + 
    ggtitle(jtitle) + 
    theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
}

PlotActivitiesWithSE <- function(dat, jtitle, showSE = TRUE, jxlab = "ZT", jsize = 20){
  nexpers <- length(unique(as.character(dat$experiment)))
  jgene <- unique(dat$gene)
  if (missing(jtitle)){
    jtitle <- jgene
  }
  if (nexpers == 1){
    m <- ggplot(dat, aes(x = time, y = exprs))
  } else {
    m <- ggplot(dat, aes(x = time, y = exprs, group = experiment, colour = experiment))
  }
  if (showSE){
    m <- m + geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se))
  }
  m <- m + geom_line() + facet_wrap(~tissue, nrow = 1) + xlab(jxlab) + ylab("Activity") + ggtitle(jtitle) + theme_bw(jsize) + 
    # theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    theme(aspect.ratio=1, strip.text = element_blank())
  if (jxlab == "ZT"){
    m <- m + scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, 12))
  } else if (jxlab == "CT"){
    m <- m + scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))
  } else {
    warning("jxlab should be ZT or CT")
  }
  return(m)
}

PlotActivitiesWithSE.wtko <- function(dat, jtitle, showSE = TRUE, jxlab = "ZT", jsize = 20, jylab = "Activity", split.by="tissue", ncols = 2, single.day=FALSE){
  nexpers <- length(unique(as.character(dat$experiment)))
  jgene <- unique(dat$gene)
  if (missing(jtitle)){
    jtitle <- jgene
  }
  # if (nexpers == 1){
  #   m <- ggplot(dat, aes(x = time, y = exprs))
  # } else {
  #   m <- ggplot(dat, aes(x = time, y = exprs, group = experiment, colour = experiment))
  # }
  # if (showSE){
  #   m <- m + geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se))
  # }
  # m <- m + geom_line() + facet_wrap(~tissue, nrow = 1) + xlab(jxlab) + ylab("Activity") + ggtitle(jtitle) + theme_bw(jsize) + 
  #   # theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #   theme(aspect.ratio=1, strip.text = element_blank())
  
  if (single.day){
    source("/home/yeung/projects/tissue-specificity/scripts/functions/PlotGeneAcrossTissues.R")
    dat <- ConvertToSingleDay(dat, use.se = TRUE)
  }
  
  m <- ggplot(dat, aes(x = time, colour = tissue, linetype = geno, y = exprs)) + 
    geom_point() + geom_line() + xlab(jxlab) + ylab(jylab) + 
    theme_bw(jsize) + ggtitle(jtitle) + 
    theme(aspect.ratio = 1, legend.position = "bottom")
  if (showSE){
    m <- m + geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se))
  }
  if (split.by == "geno"){
    m <- m + facet_wrap(~geno, ncol = ncols)
  } else if (split.by == "tissue"){
    m <- m + facet_wrap(~tissue, ncol = ncols)
  } else {
    warning("Split by must be geno or tissue")
  }
 
  if (jxlab == "ZT"){
    jxmin <- 0
    jxmax <- ceiling(max(dat$time) / 24) * 24  # round to nearest 24
    m <- m + scale_x_continuous(limits = c(jxmin, jxmax), breaks = seq(jxmin, jxmax, 12))
  } else if (jxlab == "CT"){
    m <- m + scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))
  } else {
    warning("jxlab should be ZT or CT")
  }
  return(m)
}

PlotActivitiesWithSE.rnaseq <- function(dat, jtitle, showSE = TRUE){
  jgene <- unique(dat$gene)
  if (missing(jtitle)){
    jtitle <- jgene
  }
  if (showSE){
    ggplot(dat, 
           aes(x = time, y = exprs)) +
      geom_line() +
      geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
      facet_wrap(~tissue) + 
      xlab("CT") +
      ylab("Activity") + 
      ggtitle(jtitle)
  } else {
    ggplot(dat, 
           aes(x = time, y = exprs)) +
      geom_line() +
      facet_wrap(~tissue) + 
      xlab("CT") +
      ylab("Activity") + 
      ggtitle(jtitle)
  }
}

PlotMeanActivitiesWithSE <- function(dat){
  jgene <- unique(dat$gene)
  ggplot(dat,
         aes(x = tissue, y = exprs, group = experiment, colour = experiment)) + 
    geom_line() + 
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    ylab("Activity") +
    ggtitle(jgene)
}


PlotMeanActivitiesWithSE.singleintercept <- function(dat){
  jgene <- unique(dat$gene)
  ggplot(dat,
         aes(x = tissue, y = exprs)) + 
    geom_line() + 
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    ylab("Activity") +
    ggtitle(jgene)
}

ConvertArgToPhase <- function(phase.rads, omega){
  # convert phase in radians to phase in time, using omega.
  # expect phase to be between -pi to pi, convert that
  # to 0 to 24 hours.
  
  # convert range from -pi to pi to 0 to 2pi
  phase.rads[which(phase.rads < 0)] <- phase.rads[which(phase.rads < 0)] + 2 * pi
  phase <- phase.rads / omega
  return(phase)
}

PlotComplex2.act <- function(vec.complex, labels, omega = 2 * pi / 24, title = "My title", ylab = "Amplitude of activity", xlab = "Phase of activity (CT)"){
  # Convert complex to amplitude (2 * fourier amplitude) and phase, given omega.
  # then plot in polar coordinates
  # fourier amplitudes are half-amplitudes of the sine-wave
  # http://www.prosig.com/signal-processing/FourierAnalysis.pdf
  df <- data.frame(amp = Mod(vec.complex) * 2,
                   phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                   label = labels)
  m <- ggplot(data = df, aes(x = amp, y = phase, label = label)) + 
    geom_point(size = 0.5) +
    coord_polar(theta = "y") + 
    xlab(xlab) +
    ylab(ylab) +
    geom_text(aes(x = amp, y = phase, size = amp), vjust = 0) +
    ggtitle(title) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
}

PlotEigengene <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  eigengene <- svd.obj$v[, comp]
  if (rotate){
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
    # rotate eigengene by -phase ref
    eigengene <- eigengene * Conj(rotate.factor)
  }
  v.plot <- PlotComplex2(eigengene, labels = rownames(svd.obj$v), omega = omega, title = paste("Right singular value", comp))
  print(v.plot)
}

PlotEigensamp <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  if (rotate){
    eigengene <- svd.obj$v[, comp]
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
  }
  eigensamp <- svd.obj$u[, comp]
  eigensamp <- eigensamp * rotate.factor
  u.plot <- PlotComplex2(eigensamp, labels = rownames(svd.obj$u), omega = omega, title = paste("Left singular value", comp))   
  print(u.plot)
}