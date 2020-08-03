LoadKallisto <- function(path.kallisto){
  source("scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  source("scripts/functions/GetTissueTimes.R")
  
  if (missing(path.kallisto)){
    path.kallisto <- "data/alternative_exon_usage//abundance.merged.annotated.sorted.txt"
  }
  dat.kallisto <- read.table(path.kallisto, header = TRUE)
  
  # BEGIN: break matrix into a metadata and tpm estimates: makes converting to long easier
  dat.meta <- dat.kallisto[, c("chromo", "start", "end", "gene_name", "strand")]
  rownames(dat.meta) <- dat.kallisto$target_id
  
  colnames.meta <- c("chromo", "start", "end", "gene_name", "strand", "target_id")
  colnames.tpm <- colnames(dat.kallisto)[which(!colnames(dat.kallisto) %in% colnames.meta)]
  
  dat.tpm <- dat.kallisto[, colnames.tpm]
  rownames(dat.tpm) <- dat.kallisto$target_id
  # END: break matrix into a metadata and tpm estimates
  
  tissues <- sapply(colnames(dat.tpm), function(s) strsplit(s, '_')[[1]][[1]])
  tissues <- ConvertRNASeqTissueNamesToArray(tissues)
  times <- GetTimes(colnames(dat.tpm), get_unique=FALSE)
  
  tpm.long <- data.frame(transcript_id = rep(dat.kallisto$target_id, ncol(dat.tpm)),
                         chromo = rep(dat.kallisto$chromo, ncol(dat.tpm)),
                         start = rep(dat.kallisto$start, ncol(dat.tpm)), 
                         end = rep(dat.kallisto$end, ncol(dat.tpm)), 
                         gene_name = rep(dat.kallisto$gene_name, ncol(dat.tpm)),
                         strand = rep(dat.kallisto$strand, ncol(dat.tpm)),
                         tissue = rep(tissues, each = nrow(dat.tpm)),
                         time = as.numeric(rep(times, each = nrow(dat.tpm))),
                         tpm = unlist(dat.tpm))
  return(tpm.long)
}

LoadKallistoMerged <- function(path.kallisto){
  # colnames:
  # V1: chromo
  # V2: start
  # V3: end
  # V4: transcript_id (comma separated IDs)
  # V5: strand
  # V6: gene_name
  
  # define constants
  n.tissues <- 12
  n.times <- 8
  
  source("scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  if (missing(path.kallisto)){
    path.kallisto <- "data/alternative_exon_usage/abundance.merged.annotated.sorted.bedtoolsmerged.bed"
  }
  dat.kallisto <- read.table(path.kallisto, header = FALSE)
  
  # BEGIN: break matrix into a metadata and tpm estimates: makes converting to long easier
  dat.meta <- dat.kallisto[, c("V1", "V2", "V3", "V6")]
  rownames(dat.meta) <- dat.kallisto$V4
  
  colnames.meta <- c("V1", "V2", "V3", "V6", "V5", "V4")
  colnames.tpm <- colnames(dat.kallisto)[which(!colnames(dat.kallisto) %in% colnames.meta)]
  
  dat.tpm <- dat.kallisto[, colnames.tpm]
  rownames(dat.tpm) <- dat.kallisto$V4
  # END: break matrix into a metadata and tpm estimates
  
  tissues <- rep(TissueNames(), each = n.times)
  times <- rep(TimesNames(),  n.tissues)
  
  tpm.long <- data.frame(transcript_id = rep(dat.kallisto$V4, ncol(dat.tpm)),
                         chromo = rep(dat.kallisto$V1, ncol(dat.tpm)),
                         start = rep(dat.kallisto$V2, ncol(dat.tpm)), 
                         end = rep(dat.kallisto$V3, ncol(dat.tpm)), 
                         gene_name = rep(dat.kallisto$V6, ncol(dat.tpm)),
                         strand = rep(dat.kallisto$V5, ncol(dat.tpm)),
                         tissue = rep(tissues, each = nrow(dat.tpm)),
                         time = as.numeric(rep(times, each = nrow(dat.tpm))),
                         tpm = unlist(dat.tpm))
  return(tpm.long)
}