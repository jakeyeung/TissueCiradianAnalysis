GetGeneModelKeys <- function(dat, tiss){
  rhyth.tiss <- gsub(pattern = ";", replacement = ",", x = dat$model)
  rhyth.tiss <- strsplit(rhyth.tiss, ",")[[1]]
  rhyth.tiss <- sapply(rhyth.tiss, function(tiss) strsplit(tiss, "_")[[1]][[1]], USE.NAMES = FALSE)
  
  is.rhyth <- tiss %in% rhyth.tiss
  return(data.frame(tissue = tiss, is.rhyth = is.rhyth))
}

GetMaxDist <- function(proms, tiss.vec){
  dist.tiss <- as.matrix(dist(proms))
  max.dist.tiss <- max(dist.tiss)
  rowcol.i <- which(dist.tiss == max.dist.tiss, arr.ind = TRUE)[1, ]  # take first max, there wil be two
  tiss <- paste(tiss.vec[rowcol.i[1]], tiss.vec[rowcol.i[2]], sep = ";")
  return(data.frame(max.dist = max.dist.tiss, tissue = tiss))
}

GetPromoterUsageMaxDist <- function(dat, jvar = "tpm_norm.avg"){
  proms.full <- GetPromoterUsage(dat, jvar = jvar, do.svd = TRUE, append.tiss = TRUE, get.means = FALSE, get.entropy = FALSE)  
  tiss.vec <- proms.full[[1]]$tissue
  proms <- subset(proms.full$dat.mat.trans, select = -c(amp, tissue))
  dist.tiss <- as.matrix(dist(proms))
  max.dist.tiss <- max(dist.tiss)
  rowcol.i <- which(dist.tiss == max.dist.tiss, arr.ind = TRUE)[1, ]  # take first max, there wil be two
  tiss <- paste(tiss.vec[rowcol.i[1]], tiss.vec[rowcol.i[2]], sep = ";")
  return(data.frame(max.dist = max.dist.tiss, tissue = tiss))
}
# plot top hits
PlotPromoter <- function(tpm.afe.avg.sub, jtitle="Title"){
  proms.full <- GetPromoterUsage(tpm.afe.avg.sub, jvar = "tpm_norm.avg", get.means = FALSE, get.entropy = FALSE)
  proms <- subset(proms.full$dat.mat.trans, select = -c(amp, tissue))
  plot(proms[, 1], proms[, 2], main = jtitle, xlab = "Promoter usage (1st component)", ylab = "Promoter usage (2nd component)")
  par(cex.axis=1, cex.lab=2, cex.main=2, cex.sub=1)
  text(proms[, 1], proms[, 2], proms.full$dat.mat.trans$tissue)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

PromoterSpacePlots.nostics <- function(tpm.gauss.sigs, jgene, jtitle, draw.ellipse = TRUE){
  proms <- tpm.gauss.sigs$proms
  n.proms <- ncol(proms)
  if (missing(jtitle)){
    jtitle <- paste(jgene, "n.proms", n.proms)
  }
  plot(proms[, 1], proms[, 2], main = jtitle, xlab = "Promoter usage (1st component)", ylab = "Promoter usage (2nd component)")
  par(cex.axis=1, cex.lab=2, cex.main=2, cex.sub=1)
  
  jcols <- vector(length = length(tpm.gauss.sigs$amp))
  jcols[which(tpm.gauss.sigs$amp == 1)] <- "blue"
  jcols[which(tpm.gauss.sigs$amp == 0)] <- "red"
  text(proms[, 1], proms[, 2], labels = tpm.gauss.sigs$tissue, col = jcols)
  # draw ellipse
  if (draw.ellipse){
    try(lines(ellipse(tpm.gauss.sigs$sig1$cov, level = 0.5, centre = tpm.gauss.sigs$sig1$center), type='l', col = "blue"), silent = T)
    try(lines(ellipse(tpm.gauss.sigs$sig2$cov, level = 0.5, centre = tpm.gauss.sigs$sig2$center), type='l', col = "red"), silent = T)
  }
}

PromoterSpacePlots <- function(tpm.afe.avg, jgene, jvar = "tpm_norm.avg", draw.ellipse = TRUE, use.weights = TRUE, transcript_id = "transcript_id"){
  tpm.test = tryCatch({
    tpm.test <- subset(tpm.afe.avg, gene_name == jgene & tissue != "WFAT")
  }, error = function(e) {
    tpm.test <- subset(tpm.afe.avg, gene == jgene)
  })
  test.svd <- GetPromoterUsage(tpm.test, jvar = jvar, do.svd = T, append.tiss = TRUE, get.means = TRUE, transcript_id = transcript_id)

  proms <- subset(test.svd$dat.mat.trans, select = -c(amp, tissue))
  amp <- test.svd$dat.mat.trans$amp
  weights1 <- (amp - min(amp)) / (max(amp) - min(amp))
  weights2 <- 1 - weights1
  
  # center1
  mu1 <- colSums(sweep(proms, MARGIN = 1, STATS = weights1, FUN = "*")) / sum(weights1)
  sig1 <- cov.wt(proms, wt = weights1)
  pi1 <- sum(weights1) / length(weights1)
  
  # center2
  mu2 <- colSums(sweep(proms, MARGIN = 1, STATS = weights2, FUN = "*")) / sum(weights2)
  sig2 <- cov.wt(proms, wt = weights2)
  pi2 <- sum(weights2) / length(weights2)
  
  plot(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3], main = paste(jgene, "Amp range:", diff(range(amp))))
  if (use.weights){
    myColoursAlpha <- sapply(weights1, function(a) add.alpha(1, alpha=a))
    text(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3], labels = test.svd$dat.mat.trans$tissue, col = myColoursAlpha)
  } else {
    jcols <- vector(length = length(amp))
    jcols[which(amp == 1)] <- "blue"
    jcols[which(amp == 0)] <- "red"
    text(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3], labels = test.svd$dat.mat.trans$tissue, col = jcols)
  }
  points(mu1[1], mu1[2], pch = "*", col = "blue", cex = 5)
  points(mu2[1], mu2[2], pch = "*", col = "red", cex = 5)
  # draw ellipse
  if (draw.ellipse){
    try(lines(ellipse(sig1$cov, level = 0.5, centre = sig1$center), type='l', col = "blue"), silent = T)
    try(lines(ellipse(sig2$cov, level = 0.5, centre = sig2$center), type='l', col = "red"), silent = T)
  }
  
  #   dist1 <- FuzzyDistance(proms, mu1, amp)
  #   dist2 <- FuzzyDistance(proms, mu2, amp)
  #   print(paste("Intracluster score", sum(dist1, dist2)))
  #   print(paste("Interscore", sum((mu2 - mu1) ^ 2)))
  #   
  #   # Gaussian distribution
  #   intraprob1 <- sum(weights1 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mean = mu1, sigma = sig1$cov))) / sum(weights1)
  #   intraprob2 <- sum(weights2 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu2, sigma = sig2$cov))) / sum(weights2)
  #   interprob1 <- sum(weights2 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu1, sig1$cov))) / sum(weights2)
  #   interprob2 <- sum(weights1 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu2, sig1$cov))) / sum(weights1)
}

HellingerDistance <- function(dat){
  # calculate distance betwene two rows of matrix
  if (nrow(dat) != 2){
    # print(paste("Number of columns:", ncol(dat)))
    warning("Hellinger Distance only computes distance between two rows.")
    # print("Dim:")
    # print(dim(dat))
    return(NA)
  }
  return(sqrt(sum(apply(dat, 2, function(colm) (sqrt(colm[1]) - sqrt(colm[2]))) ^ 2)))
}


CalculateDistance <- function(dat, dist.method = "hellinger", jvar = "tpm_norm.avg", transcript_id = "transcript", return.as.df=FALSE){
  # simpler than GaussianCenters
  # For Liver and Kidney, calculate distance (either euclidean or hellinger)
  if (length(unique(dat$tissue)) <= 1){
    ifelse(return.as.df, data.frame(NULL), NA)
  }
  proms <- GetPromoterUsage(dat, jvar = jvar, do.svd = FALSE, append.tiss = TRUE, get.means = FALSE, get.entropy = FALSE, transcript_id = transcript_id, get.prom.only = TRUE)
  if (nrow(proms) != 2){
    warning("Expecting only two rows in GetPromoterUsage() output")
  }
  if (dist.method == "euclidean"){
    proms.dist <- as.numeric(dist(as.matrix(proms), method = "euclidean"))
  } else if (dist.method == "hellinger"){
    proms.dist <- HellingerDistance(proms)
  }
  if (!return.as.df){
    return(proms.dist)
  } else {
    gene <- dat$gene[[1]]
    return(data.frame(gene = gene, proms.dist = proms.dist))
  }
}

CalculateGaussianCenters <- function(dat, jvar = "tpm_norm.avg", thres = 0.9, do.svd = TRUE, transcript_id = "transcript_id"){
  if (length(unique(dat$tissue)) <= 1){
    return(NA)
  }
  proms.full <- GetPromoterUsage(dat, jvar = jvar, do.svd = do.svd, append.tiss = TRUE, get.means = TRUE, get.entropy = FALSE, transcript_id = transcript_id)  
  proms <- subset(proms.full$dat.mat.trans, select = -c(amp, tissue))
  amp <- proms.full$dat.mat.trans$amp
  
  if (diff(range(amp)) == 0){
    return(NA)
  }
  weights1 <- (amp - min(amp)) / (max(amp) - min(amp))
  weights2 <- 1 - weights1
    
  # center1
  sig1 <- cov.wt(proms, wt = weights1)
  mu1 <- sig1$center
  pi1 <- sum(weights1) / length(weights1)
  
  # center2
  sig2 <- cov.wt(proms, wt = weights2)
  mu2 <- sig2$center
  pi2 <- sum(weights2) / length(weights2)
  
  return(list(mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, proms = proms, amp = amp, tissue = as.character(proms.full$dat.mat.trans$tissue)))
}

CalculateGaussianDists <- function(dat){
  # sigs is output from CalculateGaussianCenters
  # calculate inter and intra cluster scores
  
  sigs <- dat$sigs[[1]]
  proms <- sigs$proms
  amp <- sigs$amp
  
  # center1
  sig1 <- sigs$sig1
  mu1 <- sig1$center
  weights1 <- sig1$wt
  #   pi1 <- sum(weights1) / length(weights1)
  
  # center2
  sig2 <- sigs$sig2
  mu2 <- sig2$center
  weights2 <- sig2$wt
  #   pi2 <- sum(weights2) / length(weights2)
  
  center.dists <- sum((mu2 - mu1) ^ 2)
  amp.range <- max(amp) - min(amp)
  intrascore1 <- sum(weights1 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mean = mu1, sigma = sig1$cov))) / sum(weights1)
  intrascore2 <- sum(weights2 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu2, sigma = sig2$cov))) / sum(weights2)
  interscore1 <- sum(weights2 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu1, sig1$cov))) / sum(weights2)
  interscore2 <- sum(weights1 * apply(proms, 1, function(x) mvtnorm::dmvnorm(x, mu2, sig1$cov))) / sum(weights1)
  return(data.frame(center.dists=center.dists, amp.range = amp.range, intrascore1=intrascore1, intrascore2=intrascore2, interscore1=interscore1, interscore2=interscore2))
}

RunFuzzyDistance <- function(dat, jvar = "tpm_norm.avg", thres = 0.9, do.svd = TRUE){
  if (length(unique(dat$tissue)) <= 1){
    return(data.frame(amp.range = NA, intra.score = NA, inter.score = NA))
  }
  proms.full <- GetPromoterUsage(dat, jvar = jvar, do.svd = do.svd, append.tiss = TRUE, get.means = TRUE, get.entropy = FALSE)  
  proms <- subset(proms.full$dat.mat.trans, select = -c(amp, tissue))
  amp <- proms.full$dat.mat.trans$amp
  weights1 <- (amp - min(amp)) / (max(amp) - min(amp))
  weights2 <- 1 - weights1
  
  # center1
  mu1 <- colSums(sweep(proms, MARGIN = 1, STATS = weights1, FUN = "*")) / sum(weights1)
  
  # center2
  mu2 <- colSums(sweep(proms, MARGIN = 1, STATS = weights2, FUN = "*")) / sum(weights2)
  
  amp.range <- max(amp) - min(amp)
  dist1 <- FuzzyDistance(proms, mu1, amp)
  dist2 <- FuzzyDistance(proms, mu2, amp)
  intra.score <- sum(dist1, dist2)
  inter.score <- sum((mu2 - mu1) ^ 2)
  return(data.frame(amp.range = amp.range, intra.score = intra.score, inter.score = inter.score))
}

PromAmpScore <- function(prom.vec1, prom.vec2, amp1, amp2){
  euc.dist <- sqrt(sum((prom.vec1 - prom.vec2) ^ 2))
  amp.dist <- abs(amp1 - amp2)
  return((euc.dist * amp.dist))
}

# distance score
FuzzyDistance <- function(jpoints, jcluster, weights){
  # https://en.wikipedia.org/wiki/Fuzzy_clustering
  dist <- sweep(jpoints, MARGIN = 2, STATS = jcluster, FUN = "-")
  dist <- sweep(dist, MARGIN = 1, STATS = weights, FUN = "*")
  dist <- sum(apply(dist, MARGIN = 1, function(x) sum(x^2)))
  return(dist)
}


DoCanCor <- function(tpm.sub, xvar = "tpm_norm.avg", yvar = "amp"){
  # do canonical correlation of two matrices to find interesting
  # alternative promoter usage
  if (nrow(tpm.sub) == 0){
    return(data.frame(Cor = NA))
  }
  X <- dcast(data = tpm.sub, formula = transcript_id ~ tissue, value.var = xvar)
  rownames(X) <- X$transcript_id; X$transcript_id <- NULL
  X <- t(X)
  jtrans <- tpm.sub$transcript_id[[1]]
  Y <- data.frame(subset(tpm.sub, transcript_id == jtrans, select = c(yvar, "tissue")))
  rownames(Y) <- as.character(Y$tissue); Y$tissue <- NULL
  jcor <- cancor(X, Y, xcenter = T, ycenter = T)
  return(data.frame(Cor = jcor$cor))
}

KeepUpToThres <- function(vec, thres, min.dim = 2){
  # keep vector up to threshold return by index
  for (i in seq(length(vec))){
    if (vec[i] > thres){
      break
    }
  }
  if (i < min.dim){
    return(i:min.dim)
  }
  return(1:i)
}

GetPromoterUsage <- function(dat, jvar = "tpm_norm.avg", do.svd = TRUE, thres = 0.9, append.tiss = TRUE, get.means=TRUE, get.entropy=TRUE, 
                             transcript_id = "transcript_id", return.transcripts=FALSE, get.prom.only=FALSE){
  # get promoter usage
  if (!is.null(dat$mean)){
    jform <- paste0("tissue + amp + mean ~ ", transcript_id)
    # dat.mat <- dcast(dat, tissue + amp + mean~ transcript_id, value.var = jvar)
    dat.mat <- dcast(dat, as.formula(jform), value.var = jvar)
    dat.mat.prom <- subset(dat.mat, select = -c(tissue, amp, mean))
  } else {
    jform <- paste0("tissue + amp ~ ", transcript_id)
    # dat.mat <- dcast(dat, tissue + amp + mean~ transcript_id, value.var = jvar)
    dat.mat <- dcast(dat, as.formula(jform), value.var = jvar)
    dat.mat.prom <- subset(dat.mat, select = -c(tissue, amp))
  }
  if (do.svd){
    # reduce dim
    dat.mat.prom <- sweep(dat.mat.prom, MARGIN = 1, STATS = rowMeans(dat.mat.prom), FUN = "-")
    dat.mat.prom.s <- svd(dat.mat.prom)
    eigvals <- dat.mat.prom.s$d ^ 2 / sum(dat.mat.prom.s$d ^ 2)
    eigvals.cum <- cumsum(eigvals)
    # threshold for number of eigvals to keep
    keep <- KeepUpToThres(eigvals.cum, thres, min.dim = 2)
    if (length(dat.mat.prom.s$d) == 1){
      print("Unexpected here:")
      print(dat)
    }
    dat.mat.prom <- sweep(dat.mat.prom.s$u[, keep], MARGIN = 2, STATS = dat.mat.prom.s$d[keep], FUN = "*")
    dat.mat.trans <- data.frame(amp = dat.mat$amp, dat.mat.prom) 
    if (return.transcripts){
      dat.mat.v <- sweep(dat.mat.prom.s$v[, keep], MARGIN = 2, STATS = dat.mat.prom.s$d[keep], FUN = "*")
      # take most polarizing two transcripts
      transcripts <- colnames(subset(dat.mat, select = -c(tissue, amp)))
      dat.mat.v.mean <- apply(dat.mat.v, 1, function(row) sum(row ^ 2))
      
      n <- length(dat.mat.v.mean)
      tx1.i <- which.max(dat.mat.v.mean)[[1]]
      tx2.i <- which(dat.mat.v.mean == sort(dat.mat.v.mean, partial = n - 1)[n - 1])[[1]]
      
      dat.mat.v.keep <- dat.mat.v.mean[c(tx1.i, tx2.i)]
      tx <- transcripts[c(tx1.i, tx2.i)]
      
      dat.mat.trans <- cbind(dat.mat.trans, dat.mat.v.keep)
      dat.mat.trans$transcript <- tx
    }
  } else {
    dat.mat.trans <- data.frame(amp = dat.mat$amp, dat.mat.prom)
  }
  
  if (append.tiss){
    dat.mat.trans$tissue <- dat.mat$tissue
  }
  out <- list(dat.mat.trans = dat.mat.trans)
  if (get.means){
    # out$weights <- dat.mat$mean
    out$mean <- dat.mat$mean
  }
  if (get.entropy){
    dat.H <- apply(dat.mat.prom, 1, function(x) ShannonEntropy(x, normalize = FALSE))
    out$entropy <- dat.H
  }
  if (get.prom.only){
    return(subset(out$dat.mat.trans, select = c(-amp, -tissue)))
  }
  return(out)
}

CorrelateAmpPromMulti <- function(dat, jvar = "tpm_norm.avg", thres = 0.9, do.svd = TRUE, weighted = FALSE, eps = 1e-10){
  dat.mat.trans.lst <- GetPromoterUsage(dat, jvar = jvar, do.svd = do.svd, thres = thres, append.tiss = FALSE, get.means = TRUE, get.entropy = TRUE)
  dat.mat.trans <- dat.mat.trans.lst$dat.mat.trans
  weights <- dat.mat.trans.lst$weights
  entropy <- dat.mat.trans.lst$entropy
  if (weighted){
    # by shannon entropy
    weights <- 1 - entropy
    weights <- weights + eps  # if zero?
    fit.altprom <- lm(formula = amp ~ ., data = dat.mat.trans, weights = weights)
  } else {
    fit.altprom <- lm(formula = amp ~ ., data = dat.mat.trans)
  }
  # return just R-squared and best p-value
  rsqr <- summary(fit.altprom)$r.squared
  coefs <- summary(fit.altprom)$coefficients
  pval.best <- min(coefs[2:nrow(coefs), "Pr(>|t|)"])
  return(data.frame(rsqr = rsqr, pval.best = pval.best))
}

# Functions: kallisto -----------------------------------------------------

PlotDiagnostics <- function(dat.tpm, dat.arrayrnaseq, jgene, jtranscript){
  source('scripts/functions/PlotGeneAcrossTissues.R')
  # uses data objects from find_alternative_first_exons_kallisto.R OR filter_kallisto_for_start_sites.R
  # Plot diagnostic plots, given a gene and transcript
  # jgene <- "Sorbs1"; jtranscript="ENSMUST00000165469"  # Liver and BFAT looks rhythmic, but assigned as "not rhythmic"
  
  test.gene <- subset(dat.tpm, gene_name == jgene)
  test <- subset(test.gene, transcript_id == jtranscript)
  
  rhythmic.tissues <- paste(unique(test$tissue[which(test$is.rhythmic == TRUE)]), collapse = ",")
  notrhythmic.tissues <- paste(unique(test$tissue[which(test$is.rhythmic == FALSE)]), collapse = ",")
  
  print(PlotGeneAcrossTissues(subset(dat.arrayrnaseq, gene == jgene), 
                              jtitle = paste(jgene, "\nRhythmic Tissues:", rhythmic.tissues, "\nFlat Tissues:", notrhythmic.tissues)))
  
  print(PlotTpmAcrossTissues(test.gene, jtitle = paste(jgene, jtranscript), log2.transform=FALSE))
  print(PlotTpmAcrossTissues(test, jtitle = paste(jgene, jtranscript), log2.transform=TRUE))
  
  print(ggplot(test.gene, aes(y = tpm_normalized, x = tissue)) + 
          geom_boxplot() + 
          facet_wrap(~transcript_id) + 
          ggtitle(jgene) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
  print(ggplot(test, 
               aes(y = tpm_normalized, x = is.rhythmic)) + 
          geom_boxplot() + 
          ggtitle(paste(jgene, jtranscript)))
}

PlotDiagnostics2 <- function(dat.tpm, tpm.avg, dat.long, jgene, jtranscript){
  source('scripts/functions/PlotGeneAcrossTissues.R')
  # updated version plotting the linear trend
  dat.gene <- subset(dat.tpm, gene_name == jgene)
  dat.transcript <- subset(dat.gene, transcript_id == jtranscript)
  tpm.avg.sub <- subset(tpm.avg, transcript_id == jtranscript)
  
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene), jtitle = paste(jgene, jtranscript, sep = "\n")))
  print(PlotTpmAcrossTissues(dat.gene, jtitle = paste("linear-scale", jgene, jtranscript, sep = "\n"), log2.transform=FALSE))
  print(PlotTpmAcrossTissues(dat.transcript, jtitle = paste("log-scale", jgene, jtranscript, sep = "\n"), log2.transform=TRUE))
  print(ggplot(tpm.avg.sub, aes(y = tpm_norm.avg, x = relamp, label = tissue)) + geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm"))
}

IsTissueSpecific2 <- function(jdf, pval.min = 1e-5, pval.max = 0.05, cutoff = 4.89){
  # given list of pvals, check if it contains pvals less than pval.min and greater than pval.max.
  # check if pval contains values less than pval.min (rhythmic) and greater than pval.max (flat)
  # if so, return TRUE, otherwise FALSE
  avg.exprs <- jdf$int.rnaseq
  pvals <- jdf$pval
  pvals.filt <- pvals[which(!is.nan(pvals) & avg.exprs > cutoff)]
  if (min(pvals.filt) < pval.min & max(pvals.filt) > pval.max){
    jdf$is.tissue.spec.circ <- rep(TRUE, length(pvals))
  } else {
    jdf$is.tissue.spec.circ<- rep(FALSE, length(pvals))
  }
  return(jdf)
}

IsRhythmic2 <- function(pval, avg.exprs, pval.min = 1e-5, cutoff = 6){
  # ask if gene is rhythmic, not rhythmic or NA (lowly expressed or 0)
  # check if pval is less than pval.min
  if (is.nan(pval) | avg.exprs < cutoff){
    return(NA)
  }
  if (pval < pval.min){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

ModelRhythmicity <- function(dat, jformula=tpm_normalized ~ is.rhythmic){
  # check it contains both Rhythmic and NotRhythmic elements
  if (length(unique(dat$is.rhythmic)) != 2){
    return(data.frame(int = NA, coef = NA, pval = NA))
  }
  jfit <- lm(data = dat, formula = jformula)
  # int <- summary(jfit)$mat["(Intercept)", "Coef"]
  # jcoef <- summary(jfit)$mat["rhythmic.or.notRhythmic", "Coef"]
  int <- coef(jfit)[[1]]
  jcoef <- coef(jfit)[[2]]
  pval <- summary(jfit)$coefficients["is.rhythmicTRUE", "Pr(>|t|)"]
  return(data.frame(int = int, coef = jcoef, pval = pval))
}

CalculateFractionIsoformUsage <- function(tpm, pseudocount = 0){
  # given tpm of a sample across known isoforms, compute
  # fraction of isoform usage. 
  # 
  # Add pseudo count to increase robustness to lowly expressed genes
  tpm_normalized <- (tpm + pseudocount) / (sum(tpm + pseudocount))
}


# Functions: alternative first exons --------------------------------------


KsTestPosNeg <- function(jdf){
  # 
  x.neg <- jdf$sitecount[which(jdf$hit.or.not == "Neg")]
  x.pos <- jdf$sitecount[which(jdf$hit.or.not == "Pos")]
  jks <- ks.test(x.neg, x.pos)
  stat <- jks$statistic
  pval <- jks$p.value
  return(data.frame(statistic = stat, pval = pval))
}

FitPosNeg <- function(jdf){
  jfit <- lm(formula = sitecount ~ hit.or.not, data = jdf)
  summary(jfit)$coefficients["hit.or.notPos", "Pr(>|t|)"]
  pval <- summary(jfit)$coefficients["hit.or.notPos", "Pr(>|t|)"]
  coef <- summary(jfit)$coefficients["hit.or.notPos", "Estimate"]
  return(data.frame(coef = coef, pval = pval))
}

HitOrNot <- function(coef, pval, max.pval=1e-5, min.pval = 0.05){
  if (pval < max.pval){
    if (coef > 0){
      s = "Pos"
    } else if (coef < 0){
      s = "Neg"
    } else {
      warning("Coefficient is 0 but called a hit")
    }
  } else if (pval < min.pval){
    s = NA
  } else {
    s = "NotHit"
  }
  return(s)
}

FitDfToMatrix <- function(jdf, common.genes){
  dat.fitrhyth.filt <- data.frame(subset(jdf, gene %in% common.genes, select = c(gene, tissue, as.numeric(pval), amp)))  # faster
  rnames <- apply(dat.fitrhyth.filt, 1, function(x) paste0(x[2], '-', x[1]))
  rownames(dat.fitrhyth.filt) <- rnames; rm(rnames)
  dat.fitrhyth.filt <- subset(dat.fitrhyth.filt, select = c(pval, amp))
  dat.fitrhyth.filt <- data.matrix(dat.fitrhyth.filt)  # faster to work with matrices
  return(dat.fitrhyth.filt)
}

GetRhythmicOrNot <- function(x, fitdf){
  # Expect x to be a row from cov.normreads, with tissue and gene in 2nd and 3rd col
  tiss <- x[2]
  gene <- x[3]
  rname <- paste0(tiss, '-', gene)
  fitdf.sub = tryCatch({
    fitdf[rname, ]
  }, warning = function(w) {
    print("Warning")
    print(w)
  }, error = function(e) {
    # print(paste("Cannot access:", rname))
    return(NA)
  })  
  if (is.na(fitdf.sub[1])){
    return(NA)
  }
  pval <- fitdf.sub[1]
  amp <- fitdf.sub[2]
  annots <- RhythmicOrNot(pval, amp)
  return(annots)
}

RhythmicOrNot <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  if (pval < min.pval & amp > max.amp){
    return("Rhythmic")
  } else if (pval > max.pval & amp < min.amp){
    return("NotRhythmic")
  } else {
    return(NA)  # undecided
  }
}

FoldChangeRhyth <- function(jdf){
  # Calculate log2 fold change between "rhythmic" and "non rhythmic" genes 
}

cossim <- function(x, y){
  return(x %*% y / sqrt(x%*%x * y%*%y))
}

LoopCor <- function(m, show.which=FALSE, input.vec1=NA){
  imin <- NA
  jmin <- NA
  jcor.min <- 2  # init because pearson cor is between -1 and 1. Use a number outside of this range.
  for (i in 1:ncol(m)){
    if (!is.na(input.vec1)){
      i <- input.vec1
    }
    vec1 <- m[, i]
    for (j in i:ncol(m)){
      if (j == i){
        next  # dont need to compare between same vec
      }
      vec2 <- m[, j]
      jcor <- cossim(vec1, vec2)
      if (jcor < jcor.min){
        jcor.min <- jcor
        imin <- i
        jmin <- j
      }
    }
    if (!is.na(input.vec1)){
      break
    }
  }
  if (show.which){
    jtissues <- c(colnames(m)[imin], colnames(m)[jmin])
    print(jtissues)
    return(list(tissues = jtissues,
                min.cor = jcor.min))
  } 
  return(jcor.min)
}

GetMinCor <- function(df){
  m <- acast(data = df, transcript ~ tissue, value.var = "norm_reads")
  jcor.min <- LoopCor(m)
  return(data.frame(min.cor = jcor.min))
}

Normalize <- function(x, pseudocount = 1){
  if (pseudocount > 0){
    x <- x + pseudocount
  }
  x.norm <- x / sum(x)
  return(x.norm)
}

AvgAcrossTranscripts <- function(df){
  ddply(df, .(transcript, gene, tissue), summarise, mean_reads = mean(reads))
}

NormalizeReads <- function(df){
  # Normalize reads across transcripts
  ddply(df, .(gene), transform, norm_reads = Normalize(mean_reads))
}

GetGeneNameFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve RP23-27017.1
  gene.str <- strsplit(as.character(annot), ';')[[1]][[1]]
  gene.str <- strsplit(gene.str, '=')[[1]][[2]]
  return(gene.str)
}

GetTranscriptIDFromAnnot <- function(annot){
  # from gene_name=RP23-271O17.1;transcript_id=ENSMUST00000193812 retrieve ENSMUST...
  transcript.str <- strsplit(as.character(annot), ';')[[1]][[2]]
  transcript.str <- strsplit(transcript.str, '=')[[1]][[2]]
  return(transcript.str)
}

rowMax <- function(df){
  # Return vector of maximums from df
  return(apply(df, 1, max))
}

GetTissuesAFE <- function(x){
  # Get tissues from column names from Adr_CT22 format
  substr(x, "_")[[1]][[1]]
}

TissueMapping <- function(cov.to.rnaseq = TRUE){
  # Tissue names from coverage are slightly different from
  # tissue names from RNASeq data. Create the mapping between 
  # coverage to rnaseq tissue names
  list("Adr" = "Adr",
       "Aor" = "Aorta",
       "BFat" = "BFAT",
       "Bstm" = "BS",
       "Cer" = "Cere",
       "Hrt" = "Heart",
       "Hyp" = "Hypo",
       "Kid" = "Kidney",
       "Liv" = "Liver",
       "Lun" = "Lung",
       "Mus" = "Mus",
       "WFat" = "WFAT")
}

ConvertTissuesToMatchRnaSeq <- function(tissues){
  # make tissue names look like RNASeq column names using TissueMapping
  tissue.map <- TissueMapping()
  sapply(tissues, function(x){
    tissue.map[[x]]
  })
}

GetExprsAsVector <- function(dat, genes, tissuetime){
  # Given a vector of rownames and column names, extract
  # its corresponding element in the dat. Return as
  # a vector.
  if (length(genes) != length(tissuetime)) print("Genes and tissuetime is not same length")
  lookups <- mapply(function(x, y){
    return(dat[x, y])
  }, x = genes, y = tissuetime)
  return(lookups)
}

GetLocationFromAnnotation <- function(bed, gene_name, transcript_id){
  # Given bed, gene_name and transcript_id, return the chromo, start, end
  annot <- paste0("gene_name=", gene_name, ";transcript_id=", transcript_id)
  sub <- subset(bed, annotations == annot)
  # return as UCSC-style
  return(paste0(sub$chromosome, ":", sub$start, "-", sub$end))
}

SubsetBed <- function(bed, gene_name, transcript){
  # Subset bed based on grepping annotations from gene name
  if (missing(transcript)){
    return(bed[grepl(gene_name, bed$annotations), ])  
  } else if (missing(gene_name)){  
    return(bed[grepl(transcript, bed$annotations), ])  
  }
}

NormalizeBySum <- function(x){
  # Normalize a vector by its sum
  return(x / sum(x))
}

GetFirst <- function(x){
  return(x[1])
}

ShannonEntropy <- function(x.vec, normalize=FALSE){
  if (normalize){
    # should sum to 1
    x.vec <- x.vec / sum(x.vec)
  }
  entropy <- 0
  for (x in x.vec){
    entropy <- entropy + x * log2(1 / x)
  }
  entropy <- entropy / log2(length(x.vec))
  return(entropy)
}

MulitpleStarts <- function(df, min_dist = 500){
  # Check if df has multiple exons, ddply from bed
  if (nrow(df) <= 1){
    return(data.frame(MultiStart = FALSE))
  }
  dist <- diff(c(min(df$start), max(df$start)))
  if (dist >= min_dist){
    return(data.frame(MultiStart = TRUE))
  } else{
    return(data.frame(MultiStart = FALSE))
  }
}

SubsetMinPval <- function(jdf){
  # take row with minimum pval
  jdf <- jdf[order(jdf$pval), ]
  return(jdf[1, ])
}

FitRhythNonRhyth <- function(jdf, log2transf = FALSE){
  # check it contains both Rhythmic and NotRhythmic elements
  if (length(unique(jdf$rhythmic.or.not)) != 2){
    return(data.frame(int = NA, coef = NA, pval = NA))
  }
  if (log2transf){
    jform <- log2(norm_reads) ~ rhythmic.or.not
  } else {
    jform <- norm_reads ~ rhythmic.or.not
  }
  jfit <- lm(data = jdf, formula = jform)
  # int <- summary(jfit)$mat["(Intercept)", "Coef"]
  # jcoef <- summary(jfit)$mat["rhythmic.or.notRhythmic", "Coef"]
  int <- coef(jfit)[[1]]
  jcoef <- coef(jfit)[[2]]
  pval <- summary(jfit)$coefficients["rhythmic.or.notRhythmic", "Pr(>|t|)"]
  return(data.frame(int = int, coef = jcoef, pval = pval))
}

FitPromoterUsageToAmplitude <- function(dat){
  #   jweights <- 1 / sqrt(dat$tpm_norm.var)
  #   fit <- lm(formula = tpm_norm.avg ~ relamp, data = dat, weights = jweights)
  fit <- lm(formula = tpm_norm.avg ~ relamp, data = dat)
  int <- coef(fit)[["(Intercept)"]]
  relamp <- coef(fit)[["relamp"]]
  tpm_norm.range <- diff(range(dat$tpm_norm.avg))
  relamp.range <- diff(range(dat$relamp))
  # handle if relamp is "NA"
  if (is.na(fit$coefficients[["relamp"]])){
    return(data.frame(int = NA, relamp = NA, pval = NA, tpm_norm.range = NA, relamp.range = NA))
  }
  pval <- summary(fit)$coefficient["relamp", "Pr(>|t|)"]
  return(data.frame(int = int, relamp = relamp, pval = pval, tpm_norm.range = tpm_norm.range, relamp.range = relamp.range))
}
