Is.TissueSpecific <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  # Check that min pval and max amplitude passes cutoff
  # And max pval and min amplitude passes cutoff
  
  # first two checks: if there exists a rhythmic gene
  # last two checks: if there exists a non-rhythmic gene.
  if (min(pval) < min.pval & max(amp) > max.amp & 
        max(pval) > max.pval & min(amp) < min.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

Is.RhythmicAcrossTissues <- function(pval, amp, cutoff.pval = 1e-5, cutoff.amp = 0.5){
  if (max(pval) < cutoff.pval & min(amp) > cutoff.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

IsTissueSpecificDplyr <- function(dat, include.betweens = FALSE){
  gene <- dat$gene[1]
  n.trues <- length(dat$is.rhythmic[which(dat$is.rhythmic == TRUE)])
  if (include.betweens == FALSE){
    n.falses <- length(dat$is.rhythmic[which(dat$is.rhythmic == FALSE)])
  } else {
    n.falses <- length(dat$is.rhythmic[which(dat$is.rhythmic != TRUE)])
  }
  if (n.trues > 0 & n.falses > 0){
    dat$is.tiss.spec <- TRUE
    return(dat)
    #     return(data.frame(gene = gene, is.tiss.spec.rhyth = TRUE))
  } else {
    dat$is.tiss.spec <- FALSE
    return(dat)
    #     return(data.frame(gene = gene, is.tiss.spec.rhyth = FALSE))
  }
}


GetTissueSpecific <- function(dat, tissues){
  gene <- dat$gene[1]
  dat.sub <- subset(dat, is.rhythmic == TRUE)
  other.tissues <- as.character(dat.sub$tissue[which(!dat.sub$tissue %in% tissues)])
  if (length(other.tissues) != 0){
    return(data.frame(gene = gene, tiss.spec = TRUE))
  } else{
    return(data.frame(gene = gene, tiss.spec = FALSE))
  }
}

RhythTiss <- function(dat){
  # Count how many tissues are rhythmic and
  # record which tissues are rhythmic
  n.tiss <- length(dat$is.rhythmic[which(dat$is.rhythmic == TRUE)])
  if (n.tiss == 0){
    return(data.frame(n.tiss = n.tiss, tissues = NA))
  }
  dat.sub <- subset(dat, is.rhythmic == TRUE)
  tissues <- paste(dat.sub$tissue, collapse = ',')
  n.tiss <- length(dat.sub$tissue)
  return(data.frame(n.tiss = n.tiss, tissues = tissues))
}

GetPhaseDiff <- function(x, period = 24){
  phase.diff <- diff(x)
  if (abs(phase.diff) > (period / 2)){
    phase.diff <- (period - abs(phase.diff)) * -(phase.diff / abs(phase.diff))
  } 
  return(phase.diff)
}

GetRhythmic <- function(dat){
  subset(dat, is.rhythmic == TRUE)
}

GetTissueSpec <- function(dat, tissue.spec.rhyth){
  # this script is used for find_oscillating_genes.pairs.R
  tissue <- dat$tissue[1]
  rhyth.genes <- subset(tissue.spec.rhyth, tissues == tissue)$gene
  subset(dat, gene %in% rhyth.genes)
}