# Load pol2 activities from Wang
LoadPol2Activities <- function(){
  source("scripts/functions//BiomartFunctions.R")
  wang.obj.path <- "data/from_labmates//Elastic_TFs_Input_Pol2_CM_library_final.Rdata"  # mot_pol2 and res.pol2.sel
  load(wang.obj.path)
  pol2.obj <- list(sitecounts=mot_pol2, 
                   signal=res.pol2.sel)
  common.transcripts <- intersect(rownames(pol2.obj$signal), rownames(pol2.obj$sitecounts))
  pol2.obj$signal <- pol2.obj$signal[common.transcripts, ]
  pol2.obj$sitecounts <- pol2.obj$sitecounts[common.transcripts, ]
  common.genes <- Transcript2Gene(common.transcripts, return.original=TRUE)
  rownames(pol2.obj$signal) <- common.genes
  rownames(pol2.obj$sitecounts) <- common.genes
  return(pol2.obj)
}

LoadPol2Activities.all <- function(){
  source("scripts/functions//BiomartFunctions.R")
  pol2_signal.path <- "data/from_labmates/Pol2_TSS.txt"
  # make gene names as rownames then delete column of gene names.
  pol2_signal <- data.frame(read.table(pol2_signal.path, header=TRUE))
  rownames(pol2_signal) <- pol2_signal$TSname
  pol2_signal$TSname <- NULL
  
  genes <- Transcript2Gene(rownames(pol2_signal), return.original = TRUE)
  rname.genes <- make.names(genes, unique = TRUE)
  rownames(pol2_signal) <- rname.genes
  
  return(pol2_signal)
}