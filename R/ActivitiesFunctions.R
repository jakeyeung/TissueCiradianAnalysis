# ActivitiesFunctions.R
# Functions for analyzing MARA activities for RNA-Seq only or Array only (merged is RNA-Seq + Array)

FitRhythmicWeighted.singleintercept <- function(df, T = 24){
  # Input: long df with exprs, time, se, experiment (array or rnaseq) for
  # a given tissue and a given gene. 
  # 
  # Fits a weighted regression model. With weights in SE
  w = 2 * pi / T  # omega
  tissue <- unique(df$tissue)
  
  sigma.sq <- df$se ^ 2
  jweights <- 1 / sigma.sq
  fit.rhyth <- lm(exprs ~ 1 + sin(w * time) + cos(w * time), data = df, weights = jweights)
  fit.flat <- lm(exprs ~ 1, data = df, weights = jweights)
  ftest <- anova(fit.flat, fit.rhyth)
  ftest.pval <- ftest[["Pr(>F)"]][[2]]
  amp <- unname(sqrt(coef(fit.rhyth)["sin(w * time)"] ^ 2 + coef(fit.rhyth)["cos(w * time)"] ^ 2))
  phase <- unname(atan2(coef(fit.rhyth)["sin(w * time)"], coef(fit.rhyth)["cos(w * time)"]))
  if (phase < 0){
    phase <- phase + 2 * pi
  }
  phase <- phase / w
  df.out <- data.frame(tissue = tissue, amp = amp, phase = phase, pval = ftest.pval)
}