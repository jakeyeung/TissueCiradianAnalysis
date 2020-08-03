# from get_liver_peaks.R
IsTissSpec <- function(is.upper.vec, compare.vec, compare.vec2=NA, reverse=FALSE){
  # allow comparison for up to 2 vectors.
  if(is.na(compare.vec2)){
    if (all(is.upper.vec == compare.vec)){
      is.tiss <- TRUE
    } else {
      is.tiss <- FALSE
    }
    if (reverse){
      return(!is.tiss)
    } else {
      return(is.tiss)
    }
  } else {
    if (all(is.upper.vec == compare.vec | is.upper.vec == compare.vec2)){
      is.tiss <- TRUE
    } else {
      is.tiss <- FALSE
    }
    if (reverse){
      return(!is.tiss)
    } else {
      return(is.tiss)
    }
  }
  return(is.tiss)
}


IsSignalLower <- function(tiss, signal, cutoffs.tiss.lower){
  return(as.numeric(signal) < as.numeric(cutoffs.tiss.lower[[as.character(tiss)]]))
}

IsSignalUpper <- function(tiss, signal, cutoffs.tiss.upper){
  return(as.numeric(signal) > as.numeric(cutoffs.tiss.upper[[as.character(tiss)]]))
}

IsTissueSpecificLong <- function(dat, tissue.indx = 4, non.tissue.indx=c(1, 2, 3, 5, 6)){
  # check if tissue is true for upper and 
  # non.tissue is true for lower
  if (dat$is.upper[tissue.indx] == TRUE & all(dat$is.lower[non.tissue.indx]) == TRUE){
    return(data.frame(is.tiss.spec = TRUE))
  } else {
    return(data.frame(NULL))
    # return(data.frame(is.tiss.spec = FALSE))
  }
}

FilterSignal <- function(dat, cutoffs.lower, cutoffs.upper, tissue.indx, non.tissue.indx){
  if (missing(non.tissue.indx)){
    non.tissue.indx <- which(seq(nrow(dat)) != tissue.indx)
  }
  tiss.indx <- 5
  sig.indx <- 6
  is.lower <- apply(dat, 1, function(row){
    tiss <- row[[tiss.indx]]
    sig <- row[[sig.indx]]
    if (sig < cutoffs.lower[[tiss]]){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  is.upper <- apply(dat, 1, function(row){
    tiss <- row[[tiss.indx]]
    sig <- row[[sig.indx]]
    if (sig > cutoffs.upper[[tiss]]){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  # Is peak tissue-specific?
  if (is.upper[[tissue.indx]] == TRUE & all(is.lower[!tissue.indx] == FALSE)){
    is.tiss.spec <- TRUE
  } else {
    is.tiss.spec <- FALSE
  }
  
  # get difference
  jdiff <- dat$signal[tissue.indx] - mean(dat$signal[non.tissue.indx])
  return(data.frame(is.tiss.spec = is.tiss.spec, difference = jdiff))
}

# from motif_pairs_distances.R
CountDistFromMotif <- function(dat, jmotif){
  # check if jmotif is in dat, if not then return NULL
  if (length(which(dat$motif == jmotif)) == 0){
    return(NULL)
  }
  ref.row <- dat[which(dat$motif == jmotif), ]
  ref.poses <- mapply(function(start, end) mean(start, end), ref.row$start, ref.row$end)
  dat.out <- lapply(ref.poses, function(ref.pos){
    rel.pos <- apply(dat, MARGIN = 1, function(datrow){
      row.pos <- mean(as.numeric(datrow[2]), as.numeric(datrow[3]))
      return(ref.pos - row.pos)
    })
    return(data.frame(motif = dat$motif, rel.pos = rel.pos))
  })
  dat.out <- do.call(rbind, dat.out)
  # dat.out <- rbindlist(dat.out)
  return(dat.out)
  # return(do.call(rbind, dat.out))
}

# from sitecount_analysis_dhs_peak_gene_body.R
GetOrderedPeaks <- function(S, fits.best){
  # order genes by most rhythmic
  key <- as.character(fits.best$gene)
  ampval <- as.numeric(fits.best$amp.avg)
  ampdic <- hash(key, ampval)
  
  # rename Katnal1;Katnal1 to Katnal because it gives NULL otherwise
  S$amp <- sapply(S$gene, function(g){
    amp <- ampdic[[as.character(g)]]
    if (is.null(amp)) return(0)
    return(amp)
  }, USE.NAMES = FALSE, simplify = TRUE)
  S <- S[order(S$amp, decreasing = TRUE), ]
  return(S)
}


TakeIfSame <- function(x){
  if (length(unique(x)) > 1){
    return(NA)
  } else {
    return(mean(x))
  }
}

ReadDHSData <- function(path, tissues, cnames, normalize = TRUE, outlong = TRUE){
  if (missing(path)){
    path <- "/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat"
  }
  if (missing(tissues)){
    tissues <- c("Cere", "Heart", "Kidney", "Liver", "Lung", "Mus")
  }
  if (missing(cnames)){
    cnames <- c("chromo", "start", "end", tissues, "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
  }
  S <- read.table("/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat")
  
  colnames(S) <- cnames
  
  if (normalize){
    for (tiss in tissues){
      S[[tiss]] <- 10^6 * S[[tiss]] / sum(S[[tiss]])
    }
  }
  if (outlong){
    signal.vec <- unlist(S[, colnames(S) %in% tissues])
    S <- data.frame(chromo = S$chromo, start = S$start, end = S$end, 
                    peak = paste(paste(S$chromo, S$start, sep = ":"), S$end, sep = "-"), # chr1:7800234-7800734
                    tissue = rep(tissues, each = nrow(S)), 
                    signal = signal.vec, 
                    gene = S$gene, dist = S$dist)
  }
  return(S)
}

ReadSitecountsMotif <- function(path, cnames, show.time = FALSE){
  start <- Sys.time()
  N <- try(read.table(path), silent = TRUE)
  if(!is.data.frame(N)){
    return(NULL)
  }
  if (missing(cnames)){
    cnames <- c("chromo", "start", "end", "motif_peak", "sitecount", "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
  }
  colnames(N) <- cnames
  
  # split motif_peak into separate columns
  # RORA.p2;mm10_chr1:3001753-3002253 -> RORA.p2
  N$motif <- sapply(N$motif_peak, function(m) strsplit(as.character(m), ";")[[1]][[1]])
  # RORA.p2;mm10_chr1:3001753-3002253 -> chr1:3001753-3002253
  N$peak <- sapply(N$motif_peak, function(m) strsplit(strsplit(as.character(m), ";")[[1]][[2]], "_")[[1]][[2]])
  # chromo:start-end
  # N$peak <- mapply(function(chromo, start, end) )
  
  cnames.remove <- c("motif_peak", "chromo.gene", "start.gene", "end.gene", "blank", "strand")
  for (cname in cnames.remove){
    N[[cname]] <- NULL
  }
  if (show.time){
    print(Sys.time() - start)
  }
  return(N)
}

GetCoord <- function(peak, jget = "start"){
  # chr15:85950195-85950695 -> 85950195 or 85950695 depending on "start" or "end"
  # if chromo, return chr15
  if (jget == "start"){
    jgeti <- 1
  } else if(jget == "end"){
    jgeti <- 2
  } else if(jget == "chromo"){
    return(strsplit(peak, ":")[[1]][[1]])
  } else {
    print(paste("jget must be start or end", jget))
  }
  return(as.numeric(strsplit(strsplit(peak, ":")[[1]][[2]], "-")[[1]][[jgeti]]))
}

GetPeakDistance <- function(peak1, peak2){
  # Get distance between two peaks
  # check they are in same chromosomes
  if (GetCoord(peak1, jget = "chromo") != GetCoord(peak2, jget = "chromo")){
    print("Not on same chromosomes")
    return(NA)
  }
  start1 <- GetCoord(peak1, "start")
  start2 <- GetCoord(peak2, "start")
  end1 <- GetCoord(peak1, "end")
  end2 <- GetCoord(peak2, "end")
  
  dist.min <- min((end1 - start2), (start1 - end2))
  # handle negatives
  return(max(dist.min, 0))
}

AssignSitecount <- function(jkey, jhash){
  s <- jhash[[as.character(jkey)]]
  if (is.null(s)){
    return(0)
  } else {
    return(s)
  }
}

FindCutoffLong <- function(dat, signal.col = "signal", jlambdas = c(0.8, 0.2), jmus = c(-1.5, 0.8), 
                           log2.trans = TRUE, pseudo = 1e-2, take.frac = 1, jshow.fig = FALSE){
  # take.frac: sample to a fraction of the data
  if (take.frac < 1){
    dat <- dat[sample(length(dat[[signal.col]]), size = take.frac * nrow(dat), replace = F), ]
  } else if (take.frac > 1){
    print(paste("take.frac must be less than or equal to 1"))
  }
  if (log2.trans){
    cutoff <- FindCutoff(x = log2(dat[[signal.col]] + pseudo), lambdas = jlambdas, mus = jmus, k = 2, show.fig = jshow.fig)
  } else {
    cutoff <- FindCutoff(x = dat[[signal.col]], lambdas = jlambdas, mus = jmus, k = 2, show.fig = jshow.fig)
  }
  if (log2.trans){
    # return in normal scale
    return(data.frame(cutoff = 2^cutoff$maximum))
  } else {
    return(data.frame(cutoff = cutoff$maximum))
  }
}

CollapseDat <- function(dat, tissue, indx, non.tissue = "Flat", flat.style = "normal"){
  # dat: signal.cut is either 0, 1, or NA (neither 0 or 1)
  # 
  # flat.style: normal (other peaks need only one tissue to be present)
  # flat style: stringent (other peaks need all other tissues to have a peak, but liver no)
  # flat style: all (all tissues must have peak)
  # collapse dat into either "tissue" or "non-tissue"
  # indx <- which(dat$tissue == tissue)
  if (missing(indx)){
    indx <- which(dat$tissue == tissue)
  }
  tiss.sig <- dat$signal.cut[indx]
  others.sig <- dat$signal.cut[-indx]  # vec
  # check signal.cut has 1 in tissue and 0 in all others
  if (flat.style == "normal"){
    if (tiss.sig == 1 & max(others.sig) == 0){
      return(data.frame(peak.type = tissue))
    } else if (tiss.sig == 0 & max(others.sig) == 1){
      # only need one tissue to have a peak
      return(data.frame(peak.type = non.tissue))
    } else {
      return(data.frame())
    }
  } else if (flat.style == "stringent"){
    if (tiss.sig == 1 & max(others.sig) == 0){
      return(data.frame(peak.type = tissue))
    } else if (tiss.sig == 0 & min(others.sig) == 1){
      # need all non-liver tissues to have a peak
      return(data.frame(peak.type = non.tissue))
    } else {
      return(data.frame())
    }
  } else if (flat.style == "all"){
    if (tiss.sig == 1 & max(others.sig) == 0){
      return(data.frame(peak.type = tissue))
    } else if (tiss.sig == 1 & min(others.sig) == 1){
      # need all tissues to have a peak
      return(data.frame(peak.type = non.tissue))
    } else {
      return(data.frame())
    }
  } else {
    print(paste("flat.style either normal or stringent or all", flat.style))
  }

}
