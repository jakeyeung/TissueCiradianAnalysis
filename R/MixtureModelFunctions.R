# MixtureModelFunctions.R
# 2015-04-27

# Example script
# library(mixtools)
# mixmdl <- normalmixEM(exprs.vec, lambda = c(0.25, 0.75), mu = c(2.5, 9), k = 2)
# plot(mixmdl,which=2)
# lines(density(exprs.vec), lty=2, lwd=2)
# 
# cutoff <- optimize(ShannonEntropyMixMdl, interval = c(2, 8), mixmdl = mixmdl, maximum = TRUE)
# cutoff <- cutoff$maximum  # cutoff = 4.883356

#library(mixtools)

FindCutoff <- function(x, lambdas, mus, k = 2, outdir = FALSE, show.fig = TRUE){
  if (outdir != FALSE){
    pdf(file.path(outdir, "mixturefit.pdf"))
  }
  mixmdl = normalmixEM(x, lambda=lambdas, mu=mus, k = k)
  cutoff <- optimize(f = ShannonEntropyMixMdl, interval = range(mixmdl$mu), mixmdl = mixmdl, tol = 0.0001, maximum = TRUE)
  jtitle <- paste0("Cutoff: ", 2^cutoff$maximum, "\n", "Cutoff (log2): ", cutoff$maximum)
  if (show.fig){
    plot(mixmdl, which=2)
    lines(density(x), lty=2, lwd=2)
    abline(v = cutoff$maximum)  # should intersect two gaussians
  }
  if (outdir != FALSE){
    dev.off()
    save(mixmdl, file = file.path(outdir, "mixmdl.Robj"))
    sink(file.path(outdir, "cutoff.txt"))
    cat(2^cutoff$maximum)
    cat("\n")
    cat(cutoff$maximum)
    sink()
  }
  return(cutoff)
}

GammaMixmdl <- function(x, mixmdl, i = 1){
  return(dgamma(x, shape = mixmdl$gamma.pars[1, i], scale = mixmdl$gamma.pars[2, i]))
}

FindCutoff2 <- function(xspace, mixmdl){
  # Find cutoff from mixmdl: use this to get Gamma functions
  # BEGIN: get intervals
  max1 <- optimize(f = GammaMixmdl, interval = range(xspace), mixmdl = mixmdl, i = 1, maximum = TRUE)$maximum
  max2 <- optimize(f = GammaMixmdl, interval = range(xspace), mixmdl = mixmdl, i = 2, maximum = TRUE)$maximum
  cutoff <- optimize(f = ShannonEntropyMixMdl, interval = range(max1, max2), mixmdl = mixmdl, tol = 0.0001, maximum = TRUE)
  return(cutoff$maximum)
}

ShannonEntropyMixMdl <- function(x, mixmdl){
  if (mixmdl$ft == "normalmixEM"){
    if (length(mixmdl$lambda) == 2){
      p1 <- PredictFromMM(mixmdl, x)
    } else if (length(mixmdl$lambda) == 3){
      p1 <- PredictFromMM3(mixmdl, x)
    }
  } else if (mixmdl$ft == "gammamixEM"){
    p1 <- PredictFromMMGamma(mixmdl, x)
  }
  p2 <- 1 - p1
  entropy <- 0
  for (p in c(p1, p2)){
    entropy <- entropy + (p * log(1 / p))
  }
  return(entropy)
}

GetDensity <- function(mixmdl, x, i, method = "norm"){
  # Get density from a single distribution, given by i
  if (method == "norm"){
    density.i <- mixmdl$lambda[i] * dnorm(x, mean = mixmdl$mu[i], sd = mixmdl$sigma[i])
  } else if (method == "gamma"){
    density.i <- mixmdl$lambda[i] * dgamma(x, shape = mixmdl$gamma.pars[1, i], scale = mixmdl$gamma.pars[2, i])
  } else {
    warning("Method != norm or gamma")
  }
}

MixtureDensity <- function(mixmdl, x, method = "norm"){
  # http://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
  # estimates density kernel from mixture model as a function of x
  # p(x) = lambda[1] n(x; mu[1], sigma[1]) + lambda[2] n(x; mu[2], sigma[2])
  p.of.x <- GetDensity(mixmdl, x, 1, method) + GetDensity(mixmdl, x, 2, method)
}

PredictFromMM <- function(mixmdl, x){
  # Given mixture model and a value of x,
  # what is probability that the value came from
  # the FIRST gaussian curve.
  # 
  # Assumes mixture model of two gaussians
  density.1 <- GetDensity(mixmdl, x, 1)
  p.of.x <- MixtureDensity(mixmdl, x)
  p1.given.x <- density.1 / p.of.x
  return(p1.given.x)
}

PredictFromMMGamma <- function(mixmdl, x){
  density.1 <- GetDensity(mixmdl, x, 1, method = "gamma")
  p.of.x <- MixtureDensity(mixmdl, x, method = "gamma")
  p1.given.x <- density.1 / p.of.x
  return(p1.given.x)
}

MixtureDensity3 <- function(mixmdl, x, indices){
  # Estimates density kernel from mixture model (two biggest lambdas) as
  # function of x
  p.of.x <- GetDensity(mixmdl, x, indices[1]) + GetDensity(mixmdl, x, indices[2])
}

PredictFromMM3 <- function(mixmdl, x){
  # Do PredictFromMM but from a mixmdl of 3 gaussians, pick 
  # the 2 gaussians of the largest lambdas
  # get probability the value came from the LARGEST gaussian curve
  top.2 <- order(mixmdl$lambda, decreasing = TRUE)[1:2]
  density.1 <- GetDensity(mixmdl, x, top.2[1])
  p.of.x <- MixtureDensity3(mixmdl, x, top.2)
  p1.given.x <- density.1 / p.of.x
  return(p1.given.x)
}

PlotGammaMixmdl <- function(mixmdl, counts, savedir = FALSE){
  # mixmdl from mixtools Gamma functions
  # counts is in log2
  shape1 = mixmdl$gamma.pars[1, 1]
  scale1 = mixmdl$gamma.pars[2, 1]
  shape2 = mixmdl$gamma.pars[1, 2]
  scale2 = mixmdl$gamma.pars[2, 2]
  xmin = min(counts)
  xmax = max(counts)
  xspace = seq(from = xmin, to = xmax, length.out = length(counts))
  # print(paste0("Shape: ", shape1, "Scale: ", scale1))
	# plot(xspace, density(counts))
  # print('hi')
  # plot(xspace, density(counts))
  g1 <- mixmdl$lambda[1] * dgamma(xspace, shape = shape1, scale = scale1)
  g2 <- mixmdl$lambda[2] * dgamma(xspace, shape = shape2, scale = scale2)
  g1g2 <- g1 + g2
  cutoff <- FindCutoff2(xspace, mixmdl)
  if (savedir != FALSE){
    pdf(file.path(savedir, "gamma_plot.pdf"))
    plot(density(counts))
  	lines(xspace, g1, type = 'l', col = 'red')
  	lines(xspace, g2, type = 'l', col = 'red')
    lines(xspace, g1g2, pch = "-", col = 'blue')
    abline(v = cutoff)
    plot(mixmdl)
    dev.off()
  } else {
    plot(density(counts))
    lines(xspace, g1, type = 'l', col = 'red')
    lines(xspace, g2, type = 'l', col = 'red')
    lines(xspace, g1g2, pch = "-", col = 'blue')
    abline(v = cutoff)
  }
  if (savedir != FALSE){
    sink(file.path(savedir, "cutoff_gamma.txt"))
    cat(2^cutoff)
    cat("\n")
    cat(cutoff)
    sink()
  }
  return(cutoff)
}

PlotGammaDist <- function(mixmdl, counts){
  # from http://stackoverflow.com/questions/17450591/fit-gamma-mixture-to-fertility-schedule-in-r
  print(paste0("Shape: ", mixmdl$gamma.pars[1], "Scale: ", mixmdl$gamma.pars[2]))
  d3 <- function(x) mixmdl$lambda[1]*dgamma(x, mixmdl$gamma.pars[1], 1/mixmdl$gamma.pars[2]) + mixmdl$lambda[2]*dgamma(x, mixmdl$gamma.pars[3], 1/mixmdl$gamma.pars[4])
  x <- seq(min(counts), max(counts), 0.001)
  plot(x, d3(x), col = 'red', pch = ".")
  # hist(counts, col="pink", add=T, freq=F, breaks=10000)
  lines(density(counts), col = 'blue')
}