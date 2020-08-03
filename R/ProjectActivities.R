ProjectActivities <- function(in.fpath, remove.wfat=TRUE){
  source("scripts/functions/GetTissueTimes.R")
  source("scripts/functions/SvdFunctions.R")
  source("scripts/functions/ConvertLongToWide.R")
  source("scripts/functions/PlotActivitiesFunctions.R")
  # Create long dataframe ---------------------------------------------------
  
  act.all <- read.table(in.fpath)
  
  tissues <- GetTissues(colnames(act.all))
  times <- GetTimes(colnames(act.all))
  act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                         tissue = rep(tissues, each = nrow(act.all) * length(times)),
                         time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                         exprs = as.numeric(unlist(act.all)))
  
  
  # Convert to Fourier ------------------------------------------------------
  
  # split
  act.split <- split(act.long, act.long$tissue)
  
  # REMOVE WFAT
  if (remove.wfat){
    act.split$WFAT <- NULL
  }
  
  omega <- 2 * pi / 24
  act.split.proj <- lapply(act.split, function(df){
    ddply(df, .(gene), ProjectToFrequency2, omega = omega)
  })
  
  # add tissue information
  for (tissue in names(act.split.proj)){
    act.split.proj[[tissue]]$tissue <- tissue
  }
  
  act.proj <- do.call(rbind, act.split.proj)
  act.proj$phase <- ConvertArgToPhase(Arg(act.proj$exprs.transformed), omega)
  act.proj$amp <- 2 * Mod(act.proj$exprs.transformed)
  return(act.proj)
}