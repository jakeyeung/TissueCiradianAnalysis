# 2017-04-18
# Jake Yeung

AdjustPhase <- function(phase.in, half.life, rads = FALSE, omega = 2 * pi / 24, fw.bw = "bw"){
  # Adjust phase given half.life
  gamma.mrna <- half.life / log(2)
  k.mrna <- 1 / gamma.mrna
  delay <- atan(omega / k.mrna)  # rads
  if (!rads){
    delay <- delay / omega  # rads to hours
  }
  if (fw.bw == "bw"){
    phase.adj <- phase.in - delay
  } else {
    phase.adj <- phase.in + delay
  }
  # set between 0 and 24
  # shift.fac <- 24 * floor(phase.adj / 24)
  # phase.adj <- phase.adj - shift.fac
  # if (phase.adj < 0){
  #   phase.adj <- phase.adj + 24
  # } else if (phase.adj >= 24){
  #   phase.adj <- phase.adj - 24
  # }
  return(phase.adj)
}