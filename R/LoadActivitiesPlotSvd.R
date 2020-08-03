# Jake Yeung
# LoadActivitiesPlotSvd
# 2016-04-12

LoadActivitiesPlotSvd <- function(indir, comp){
  
  source("scripts/functions/GetTissueTimes.R")
  source("scripts/functions/ActivitiesMergedFunctions.R")
  source("scripts/functions/PlotFunctions.R")
  source("scripts/functions/SvdFunctions.R")
  source("scripts/functions/PlotActivitiesFunctions.R")
  source("scripts/functions/LoadActivitiesLong.R")
  source("scripts/functions/RemoveP2Name.R")
  
  library(dplyr)  # problems with plyr old package
  library(reshape2)
  
  # Load --------------------------------------------------------------------
  act.long <- LoadActivitiesLong(indir)
  
  # To complex --------------------------------------------------------------
  
  act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")
  
  act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight
  act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
  act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)
  
  # SVD ---------------------------------------------------------------------
  
  # s.act <- SvdOnComplex(subset(act.complex, ! gene %in% c("RORA.p2", "NFIL3.p2")), value.var = "exprs.adj")
  s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")
  # s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")
  
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 14)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}
