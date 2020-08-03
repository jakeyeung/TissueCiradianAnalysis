# copied from TrackHubFunctions.R in 4c-seq project
#library(PhaseHSV)

HexColorToRgb <- function(hex){
  # HEX input: #FFC8DD
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  return(paste(R, G, B, sep = ","))
}

SaturationCurve <- function(x, Vmax, k, x0){
  y <- (Vmax * 2 / (1 + exp(-k * (x - x0)))) - 1
}

RotatePhase <- Vectorize(function(phase, rotate.hr=0){
  # rotate phase, in hours. If negative add 24. If > 24 minus 24.
  if (is.na(phase)) return(phase)
  phase <- phase + rotate.hr
  if (phase < 0) phase <- phase + 24
  if (phase > 24) phase <- phase - 24
  return(phase)
}, "phase")

StepCurve <- function(x, x.step){
  # Map x to step function from 0 to 1, 1 if x >= x.step
  return(ifelse(x >= x.step, 1, 0))
}

CutoffMap <- function(amp, pval, amp.cutoff, pval.cutoff){
  if (any(is.na(c(amp, pval)))){
    return(0)
  }
  if (amp > amp.cutoff & pval < pval.cutoff){
    return(1)
  } else {
    return(0)
  }
}

HexToHsv <- Vectorize(function(hex){
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  rgb2hsv(R, G, B)  
}, "hex")

PhaseAmpPvalToColor <- function(phase, amp, pval, rotate.hr=-8, amp.k = 2, pval.k = 0.25, method = "smooth",
                                black.to.white = TRUE){
  # method: "smooth" or "step"
  # rotate.phase: rotate phase by some hours
  phase <- RotatePhase(phase, rotate.hr)
  phase.col <- PhaseToHsv(phase, min.phase = 0, max.phase = 24)
  if (method == "smooth"){
    amp.col <- SaturationCurve(amp, Vmax = 1, k = amp.k, x0 = 0)
    pval.col <- SaturationCurve(-log10(pval), Vmax = 1, k = pval.k, x0 = 0)
  } else if (method == "step"){
    warning("Function not optimized. Some white spots may occur")
    amp.col <- StepCurve(amp, x.step = amp.k)
    pval.col <- StepCurve(-log10(pval), x.step = pval.k)
  } else if (method == "cutoff"){
    amp.col <- mapply(function(amp, pval) CutoffMap(amp, pval, amp.k, pval.k), amp, pval)
    pval.col <- amp.col
  }
  # convert bad phase amp pvals to 0,0,0
  bad.i <- which(is.na(phase.col) | is.na(amp.col) | is.na(pval.col))
  # make it 0 0 0
  if (length(bad.i) > 0){
    warning(paste0("Converting to 0,0,0 for:", paste0(bad.i, collapse = ",")))
    phase.col[bad.i] <- 0; amp.col[bad.i] <- 0; pval.col[bad.i] <- 0  
  }
  hsv.col <- hsv(phase.col, amp.col, pval.col)
  # convert black to white optionally
  if (black.to.white){
    hsv.col <- gsub(pattern = "#000000", replacement = "#FFFFFF", x = hsv.col)
    
  }
  return(hsv.col)
}

