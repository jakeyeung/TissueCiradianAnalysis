GetActSvd <- function(act.long, pval.adj.cutoff = 0.05){
  # Which ones are most rhythmic? -------------------------------------------
  act.fit <- act.long %>%
    group_by(tissue, gene) %>%
    do(FitRhythmicWeighted(dat = .))
  
  # False discovery rate adj ------------------------------------------------
  
  act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")
  
  # Show top genes for each tissue ------------------------------------------
  
  head(arrange(data.frame(act.fit), pval), n = 50)
  
  suppressMessages(library(plyr))
  source("~/projects/tissue-specificity/scripts/functions/SvdFunctions.R")  # many script-specific functions here
  
  omega <- 2 * pi / 24
  
  start.time <- Sys.time()
  act.complex <- lapply(split(act.long, act.long$tissue), function(x){
    ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE)
  }) %>%
    do.call(rbind, .) %>%
    mutate(magnitude = Mod(exprs.transformed)) %>%
    arrange(desc(magnitude))
  # print(Sys.time() - start.time)
  # head(act.complex, n = 100)
  
  detach("package:plyr", unload=TRUE)
  
  library(dplyr)
  
  # SVD on complex matrix ---------------------------------------------------
  # optinally filter out motifs
  
  filter.motifs <- subset(act.fit, pval.adj <= pval.adj.cutoff)$gene
  # act.complex.mat <- dcast(data = act.complex, formula = gene ~ tissue, value.var = "exprs.transformed")
  act.complex.mat <- dcast(data = subset(act.complex, gene %in% filter.motifs), formula = gene ~ tissue, value.var = "exprs.transformed")
  rownames(act.complex.mat) <- act.complex.mat[, "gene"]
  act.complex.mat <- act.complex.mat[, 2:ncol(act.complex.mat)]
  
  act.svd <- svd(act.complex.mat) 
  
  # add row and colnames
  rownames(act.svd$u) <- rownames(act.complex.mat)
  rownames(act.svd$v) <- colnames(act.complex.mat)
  return(act.svd)
}