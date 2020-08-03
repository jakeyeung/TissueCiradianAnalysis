GetTissSpecPeaks <- function(S.long, jgenes, distfilt, jcutoff, jcutoff.low, rhyth.tiss, flat.tiss){
  # get tiss spec peaks
  S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
  jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
  print(paste("number of peaks surrounding genes", length(jpeaks)))
  
  # take peaks with Liver signal greater than cutoff
  jtiss <- levels(S.sub$tissue)
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% flat.tiss)
  
  S.sub.tisspeaks <- S.sub %>%
    group_by(peak, gene) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
    filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)
  
  return(S.sub.tisspeaks)
}

FilterReadcounts <- function(dhs.dat, dhs.reps, good.samples, cutoff){
  dhs.clean <- data.frame(chr = dhs.dat$chr,
                          start = dhs.dat$start,
                          end = dhs.dat$end,
                          filtered_norm_counts = apply(dhs.reps[, good.samples], 1, mean))
  dhs.clean.filtered <- dhs.clean[which(dhs.clean$filtered_norm_counts > cutoff), ]
  plot(hist(log2(dhs.clean.filtered$filtered_norm_counts), 100))
  return(dhs.clean.filtered)
}

GetPeakwidth <- function(col1){
  # Given column 1 row 1 of autocorrelation.txt from homer,
  # get Peakwidth
  peakwidth <- strsplit(col1, split = '\\.\\.')[[1]][[4]] 
  peakwidth <- as.numeric(peakwidth)
  return(peakwidth)
}

GetFragLength <- function(col1){
  # given col1 row1 get fraglength from autocorrelation.txt
  fraglength <- as.numeric(strsplit(col1, split = '\\.\\.')[[1]][[2]])
  return(fraglength)
}

Cnames <- function(){
  # col2 is distance
  # col2 and col3 just + and - strand
  return(c("distance", "pos", "neg"))
}

Cnames.long <- function(){
  return(c("distance", "n.obs", "strand", "peakwidth", "fraglength", "samp"))
}

ReadAutocorrelation <- function(f.path, samp.name){
  # Read Autocorrelation.txt from homer output
  samp.ac <- read.table(f.path, sep='\t', header=TRUE)
  peakwidth <- GetPeakwidth(colnames(samp.ac)[1])
  fraglength <- GetFragLength(colnames(samp.ac)[1])
  colnames(samp.ac) <- Cnames()  # makes calling vectors easier
  # column 1 header contains peakwidth info
  samp.ac.long <- data.frame(distance = rep(samp.ac$distance, 2),
                             n.obs = c(samp.ac$pos, samp.ac$neg),
                             strand = c(rep("+", length(samp.ac$pos)), rep("-", length(samp.ac$neg))),
                             peakwidth = rep(peakwidth, nrow(samp.ac) * 2),
                             fraglength = rep(fraglength, nrow(samp.ac) * 2),
                             samp = rep(samp.name, nrow(samp.ac) * 2))
  cnames <- Cnames.long()
  colnames(samp.ac.long) <- cnames
  return(samp.ac.long)
}

GetTcdCnames <- function(){
  return(c(GetTcdCnames.partial(), "samp"))
}

GetTcdCnames.partial <- function(){
  return(c("tags.per.pos", "frac.pos"))
}

ReadTagCountsDistribution <- function(f.path, samp.name){
  samp.tcd <- read.table(f.path, sep='\t', header=TRUE)
  colnames(samp.tcd) <- GetTcdCnames.partial()
  samp.tcd.long <- data.frame(tags.per.pos = samp.tcd$tags.per.pos,
                              frac.pos = samp.tcd$frac.pos,
                              samp = rep(samp.name, nrow(samp.tcd)))
  return(samp.tcd.long)
}