AcosBsinToAmpPhase <- function(acos, bsin, T.period = 24){
  # convert a and b from y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  # to phase (modulo period) and amplitude. 
  amp <- sqrt(acos ^ 2 + bsin ^ 2)
  phase.rad <- atan2(acos, bsin)
  phase.time <- (phase.rad / w) %% T.period
  return(amp, phase.time)
}

FitRhythmic <- function(dat, T.period = 24, get.residuals=FALSE, get.se=FALSE){
  # Fit rhythmic model to long dat. Expect dat to be per tissue per gene.
  # use lapply and ddply to vectorize this code.
  # dat needs columns: tissue time experiment exprs
    
  # Get parameters from complex model: used later
  GetParamsRhythModel <- function(myfit, n.experiments, get.se=FALSE){
    model.params <- coef(myfit)
    if (n.experiments > 1){
      params.lst <- list(intercept.array = model.params[1],
                         intercept.rnaseq = model.params[2],
                         a = model.params[3],
                         b = model.params[4])
    } else {
      params.lst <- list(intercept = model.params[1],
                         a = model.params[2],
                         b = model.params[3])
    }
    if (get.se){
      se.lst <- summary(myfit)$coefficients[, "Std. Error"]
      # a is cos, b is sin
      cos.cname.i <- grep("^cos", x = names(se.lst))
      sin.cname.i <- grep("^sin", x = names(se.lst))
      params.lst$a.se <- se.lst[[cos.cname.i]]
      params.lst$b.se <- se.lst[[sin.cname.i]]
    }
    return(params.lst)
  }
  
  #   tissue <- unique(dat$tissue)
  tissue <- dat$tissue[1]
  w = 2 * pi / T.period
  n.experiments = length(unique(as.character(dat$experiment)))
  # Expect columns exprs and time
  if (n.experiments > 1){
    rhyth.fit <- lm(exprs ~ 0 + experiment + cos(w * time) + sin(w * time), data = dat)
    flat.fit <- lm(exprs ~ 0 + experiment, data = dat)
  } else {
    rhyth.fit <- lm(exprs ~ 1 + cos(w * time) + sin(w * time), data = dat)
    flat.fit <- lm(exprs ~ 1, data = dat)
  }
  compare.fit <- anova(flat.fit, rhyth.fit)
  pval <- compare.fit["Pr(>F)"][[1]][2]
  model.params <- GetParamsRhythModel(rhyth.fit, n.experiments, get.se)  # y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  amp <- sqrt(model.params$a ^ 2 + model.params$b ^ 2)
  phase.rad <- atan2(model.params$b, model.params$a)
  phase.time <- (phase.rad / w) %% T.period
  if (get.se){
    source("/home/yeung/projects/tissue-specificity/scripts/functions/CosSineFunctions.R")
    amp.se <- GetAmp.se(a = model.params$a, b = model.params$b, sig.a = model.params$a.se, sig.b = model.params$b.se, n = 2)  # peak to trough
    phase.se <- GetPhi.se(a = model.params$a, b = model.params$b, sig.a = model.params$a.se, sig.b = model.params$b.se, omega = 2 * pi / T.period)
  }
  if (get.residuals){
    # only care about residuals: do this if you want to see how the fit varies by changing period.
    ssq.residuals <- anova(rhyth.fit)["Residuals", "Sum Sq"]
    ssq.residuals.flat <- anova(flat.fit)["Residuals", "Sum Sq"]
    variance <- ssq.residuals / nrow(dat) # both rhyth and flat the same
    
    if (n.experiments > 1){
      dat.out <- data.frame(tissue = tissue, 
                            cos.part = model.params$a, sin.part = model.params$b, 
                            amp = amp, phase = phase.time, pval = pval,
                            int.array = model.params$intercept.array, int.rnaseq = model.params$intercept.rnaseq,
                            variance = variance,
                            ssq.residuals = ssq.residuals,
                            ssq.residuals.flat = ssq.residuals.flat,
                            variance = variance)
    } else {
        dat.out <- data.frame(tissue = tissue, 
                            cos.part = model.params$a, sin.part = model.params$b, 
                            amp = amp, phase = phase.time, pval = pval,
                            int = model.params$intercept,
                            variance = variance,
                            ssq.residuals = ssq.residuals,
                            ssq.residuals.flat = ssq.residuals.flat,
                            variance = variance)
    }
  } else {
    if (n.experiments > 1){
      dat.out <- data.frame(tissue = tissue, 
                            cos.part = model.params$a, sin.part = model.params$b, 
                            amp = amp, phase = phase.time, pval = pval,
                            int.array = model.params$intercept.array, int.rnaseq = model.params$intercept.rnaseq)
    } else {
      dat.out <- data.frame(tissue = tissue, 
                            cos.part = model.params$a, sin.part = model.params$b, 
                            amp = amp, phase = phase.time, pval = pval,
                            int = model.params$intercept)
    }
  }
  if (get.se){
    dat.out$cos.se = model.params$a.se
    dat.out$sin.se = model.params$b.se
    dat.out$amp.se = amp.se
    dat.out$phase.se = phase.se
  }
  return(dat.out)
}

FitRhythmicDatLong <- function(dat.long, jget.residuals=FALSE, get.se=FALSE, quiet=FALSE){
  library(parallel)
  dat.long.by_genetiss <- group_by(dat.long, gene, tissue)
  dat.long.by_genetiss.split <- split(dat.long.by_genetiss, dat.long.by_genetiss$tissue)
  if (!quiet){
    print("Finding rhythmic genes (~3 minutes)")
    start <- Sys.time()
  }
  dat.fitrhyth.split <- mclapply(dat.long.by_genetiss.split, function(jdf){
    rhyth <- jdf %>%
      group_by(gene) %>%
      do(FitRhythmic(dat = ., get.residuals = jget.residuals, get.se=get.se))
  }, mc.cores = 12)
  dat.fitrhyth <- do.call(rbind, dat.fitrhyth.split)
  if (!quiet){
    print(Sys.time() - start)
  }
  return(dat.fitrhyth)
}

FindMostRhythmic <- function(dat, colname="amp", decreasing = TRUE){
  # return first row after sorting by colname
  return(dat[order(dat[[colname]], decreasing = decreasing), ][1, ])
}

GetAmpRelToMax <- function(dat, fits.mostrhythmic){
  tissue <- as.character(dat$tissue[1])
  amp.max <- fits.mostrhythmic[tissue, ]$amp 
  dat$relamp <- dat$amp / amp.max
  return(dat)
}

GetRelamp <- function(fits, max.pval = 1e-3){
  # using FindMostRhythmic and GetAmpRelToMax together
  fits.mostrhythmic <- fits %>%
    group_by(tissue) %>%
    filter(pval <= max.pval) %>%
    #   filter(pval.adj <= max.pvaladj) %>%
    do(FindMostRhythmic(.)) %>%
    data.frame(.)
  rownames(fits.mostrhythmic) <- fits.mostrhythmic$tissue  # indexing
  
  # Get amplitude relative to max amp ---------------------------------------
  fits.relamp <- fits %>%
    group_by(tissue) %>%
    do(GetAmpRelToMax(., fits.mostrhythmic))
  return(fits.relamp)
}

GetAmpRelToGene <- function(dat, base.amp.dic){
  tissue <- as.character(dat$tissue[1])
  amp.base <- base.amp.dic[[tissue]]
  dat$relamp <- dat$amp / amp.base
  return(dat)
}

GetRelampByGene <- function(fits, by.gene = "Nr1d1"){
  library(hash)
  fits.base <- subset(fits, gene == by.gene)
  base.amp.dic <- hash(as.character(fits.base$tissue), fits.base$amp)
  fits.relamp <- fits %>%
    group_by(tissue) %>%
    do(GetAmpRelToGene(., base.amp.dic))
  return(fits.relamp)
}

IsRhythmicApply <- function(x, pval.min = 1e-4, pval.max = 0.05, relamp.max = 0.1, cutoff = 5.6, pvali = 7, relampi = 10, avg.exprsi = 9){
  # hard-coded shit... not the best but works fast, you should do quick sanity check 
  pval <- as.numeric(x[pvali])
  relamp <- as.numeric(x[relampi])
  avg.exprs <- as.numeric(x[avg.exprsi])
  if (is.nan(pval) | avg.exprs < cutoff){
    return(NA)
  }
  if (pval <= pval.min & relamp >= relamp.max){
    return(TRUE)
  } else if (pval >= pval.max & relamp <= relamp.max){
    return(FALSE)
  } else {
    return("Between")
  }
}
