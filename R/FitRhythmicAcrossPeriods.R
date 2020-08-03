FitRhythmicScanPeriods <- function(dat.long, periods, cores = 51){
  # likely a filtered dat.long set
  source("~/projects/tissue-specificity/scripts/functions/FitRhythmic.R")
  library(parallel)
  #   periods <- list(periods)
  dat.fitrhyth.period <- mclapply(periods, function(p){
    dat.fit <- dat.long %>%
      group_by(gene, tissue) %>%
      do(FitRhythmic(., T = p, get.residuals=TRUE))
    dat.fit$period = p
    return(dat.fit)
  }, mc.cores = 64)
  dat.fitrhyth <- do.call(rbind, dat.fitrhyth.period)
  return(dat.fitrhyth)
}

GetMinPeriod <- function(dat){
  return(subset(dat, mean.ssq.residuals == min(dat$mean.ssq.residuals)))
}

GetMinPeriodSsqResiduals <- function(dat){
  return(subset(dat, chi.sqr == min(dat$chi.sqr)))
}

PlotFitTwoPeriods <- function(dat.sub, period1, period2, tiss = "tissue", gen = "gene", exper = "array"){
  w1 <- 2 * pi / period1
  w2 <- 2 * pi / period2
  fit1 <- lm(exprs ~ 0 + experiment + sin(w1 * time) + cos(w1 * time), dat.sub)
  fit2 <- lm(exprs ~ 0 + experiment + sin(w2 * time) + cos(w2 * time), dat.sub)
  
  dat.sub.array <- subset(dat.sub, experiment == exper)
  
  plot(dat.sub.array$time, predict(fit1, dat.sub.array), "o", col = "blue", ylim = range(dat.sub$exprs), 
       main = paste0(tiss, " ", gen, " ", exper, " T=24h (blue) vs T=", period.min, "h (red)"))
  lines(dat.sub.array$time, predict(fit2, dat.sub.array), "o", col = "red")
  points(dat.sub.array$time, dat.sub.array$exprs, pch="*", col = "black")
}