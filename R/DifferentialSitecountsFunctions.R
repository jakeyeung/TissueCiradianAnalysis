PlotDifferentialSitecounts <- function(N, jgene, rhyth.tiss, flat.tiss, val = "motevo.value"){
  
  motevo.sub <- subset(N, gene == jgene)
  
  m <- dcast(data = motevo.sub, formula = motif ~ tissue, value.var = "motevo.value")
  rownames(m) <- m$motif
  m$motif <- NULL
  m <- as.matrix(m)
  # normalize by sum
  #   m <- sweep(m, 2, colSums(m),"/")
  #   m <- m[which(rowMeans(m) > 0), ]
  
  m.rhythmic <- m[, which(colnames(m) %in% rhyth.tiss)]
  
  # m.liver <- m[, "Liver"]
  # m.other <- m[, c("Cere", "Heart", "Kidney", "Lung", "Mus")]
  m.flat <- m[, which(colnames(m) %in% flat.tiss)]
  # m.flat <- m[, c("Liver", "Kidney")]
  if (!is.null(ncol(m.rhythmic))){
    m.rhythmic.avg <- rowMeans(m.rhythmic)
  } else {
    m.rhythmic.avg <- m.rhythmic
  }
  if (!is.null(ncol(m.flat))){
    m.flat.avg <- rowMeans(m.flat)
  } else {
    m.flat.avg <- m.flat
  }
  abs.diff <- m.rhythmic.avg - m.flat.avg
  abs.diff <- abs.diff[order(abs.diff, decreasing = TRUE)]
  # (head(sort(abs.diff, decreasing = TRUE), n = 100))
  
  par(mar=c(5, 15, 4.1, 2.1))
  barplot(abs.diff[1:50], names.arg = names(abs.diff[1:50]), las = 1, horiz = TRUE)  
}

MeanMat <- function(dat, cnames){
  # Take average of cnames and return
  if (length(cnames) > 1){
    dat <- rowMeans(dat[, cnames])
  } else {
    dat <- dat[, cnames]
  }
}

cancor.svd <- function(X, Y){
  # Do canonical correlation via SVD
  # http://www.nr.com/whp/notes/CanonCorrBySVD.pdf
  
  rank.X <- rankMatrix(X)[1]
  rank.Y <- rankMatrix(Y)[1]
  min.rank <- min(rank.X, rank.Y)
  print(paste("Rank:", min.rank))
  
  # Setp 1: SVD both X and Y
  s.X <- svd(X)
  s.Y <- svd(Y)
  rownames(s.X$v) <- colnames(X)
  rownames(s.X$u) <- rownames(X)
  rownames(s.Y$v) <- colnames(Y)
  rownames(s.Y$u) <- rownames(Y)
  
  # Step 2: Form u.xTu.y, then SVD THAT
  s.uxuy <- svd(t(s.X$u[, seq(min.rank)]) %*% s.Y$u[, seq(min.rank)])
  
  # Step 3: define a and b, matrices of linear combinations 
  d.X <- diag(s.X$d[seq(min.rank)])
  d.Y <- diag(s.Y$d[seq(min.rank)])
  
  
  a <- s.X$v[, seq(min.rank)] %*% solve(d.X) %*% s.uxuy$u[, seq(min.rank)]
  b <- s.Y$v[, seq(min.rank)] %*% solve(d.Y) %*% s.uxuy$v[, seq(min.rank)]
  
  # Do quick check that D == S and uTu == 1
  print("Checking D == S")
  D <- t(a) %*% t(X) %*% Y %*% b
  print(D)
  print(diag(s.uxuy$d))
  
  print("Checking uTu == 1")
  uTu <- t(a) %*% t(X) %*% X %*% a
  print(uTu)
  print("Checking vTv == 1")
  vTv <- t(b) %*% t(Y) %*% Y %*% b
  print(vTv)
  rownames(a) <- colnames(X)
  rownames(b) <- colnames(Y)
  return(list(s.X = s.X, s.Y = s.Y, a = a, b = b))
}



SubLongToMat <- function(N, jgene, jvar = "motevo.value"){
  # get MOTIF matrix
  N.sub <- subset(N, gene == jgene)
  if (nrow(N.sub) == 0) return(NA)
  
  X.motif <- dcast(data = N.sub, formula = motif ~ tissue, value.var = jvar)
  rownames(X.motif) <- X.motif$motif
  X.motif$motif <- NULL
  return(as.matrix(X.motif))
}

ScaleRemoveInfs <- function(X, center = jcenter, scale = jscale){
  # Scale and center, if scale and variances are zero, remove them from output
  X <- t(scale(t(as.matrix(X)), center = jcenter, scale = jscale))
  X <- X[complete.cases(X), ]
  return(X)
}

MatchColumns <- function(X.motif, X.exprs){
  # make X.expr and X.motif have same colnames
  if (ncol(X.motif) == ncol(X.exprs)){
    X.motif <- X.motif[, colnames(X.exprs)]
  } else {
    tiss.sub <- colnames(X.motif)
    if (length(tiss.sub) <= 2){
      print(paste(length(tiss.sub), "tissues is not enough to do analysis, skipping."))
      return(NA)
    }
    X.exprs <- X.exprs[, tiss.sub]
    print(paste(jgene, "does not have sitecounts in all tissues, using subset"))
    print(tiss.sub)
  }
  return(list(X.motif = X.motif, X.exprs = X.exprs))
}

PlotDiagnostics <- function(N, X.exprs, dat.rhyth.relamp, dat.long, jgene, jscale = TRUE, jvar = "motevo.value"){
  # get MOTIF matrix
  N.sub <- subset(N, gene == jgene)
  if (nrow(N.sub) == 0) return(NA)
  
  X.motif <- dcast(data = N.sub, formula = motif ~ tissue, value.var = jvar)
  rownames(X.motif) <- X.motif$motif
  X.motif$motif <- NULL
  
  # make X.expr and X.motif have same colnames
  if (ncol(X.motif) == ncol(X.exprs)){
    X.motif <- X.motif[, colnames(X.exprs)]
  } else {
    tiss.sub <- colnames(X.motif)
    if (length(tiss.sub) <= 2){
      print(paste(length(tiss.sub), "tissues is not enough to do analysis, skipping."))
      return(NA)
    }
    X.exprs <- X.exprs[, tiss.sub]
    print(paste(jgene, "does not have sitecounts in all tissues, using subset"))
    print(tiss.sub)
  }
  
  # center stuff
  X.exprs <- t(scale(t(as.matrix(X.exprs)), center = TRUE, scale = jscale))
  X.exprs <- X.exprs[complete.cases(X.exprs), ]
  X.motif <- t(scale(t(as.matrix(X.motif)), center = TRUE, scale = jscale))
  X.motif <- X.motif[complete.cases(X.motif), ]
  
  # Plot gene expression across tissues
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
  
  # Plot the biplot (learn how to interpret this) for covariance in X.motif
  p.motif <- prcomp(X.motif, center = TRUE, scale. = FALSE)
  biplot(p.motif, main = paste(jgene, "Motif PCA"), cex = c(0.5, 1.6), pch = 20) 
  
  #   p.exprs <- prcomp(X.exprs, center = TRUE, scale. = FALSE)
  #   biplot(p.exprs, main = paste(jgene, "Exprs PCA"), cex = c(0.5, 1.6), pch = 20)
  #   
  #   cancor.out <- cancor.svd(t(X.exprs), t(X.motif))
  #   rownames(cancor.out$a) <- rownames(X.exprs)
  #   rownames(cancor.out$b) <- rownames(X.motif)
  #   
  #   # biplot(cancor.out$a[, c(1, 2)], cancor.out$b[, c(1, 2)], main = paste(jgene, "Canonical Correlation"))
  #   plot(cancor.out$a[, 1], cancor.out$a[, 2], main = paste(jgene, "Canonical Correlation"), pch = ".")
  #   text(cancor.out$a[, 1], cancor.out$a[, 2], labels = rownames(cancor.out$a))
  #   abline(h = 0); abline(v = 0)
  #   
  #   plot(cancor.out$b[, 1], cancor.out$b[, 2], main = paste(jgene, "Canonical Correlation"), pch = ".")
  #   text(cancor.out$b[, 1], cancor.out$b[, 2], labels = rownames(cancor.out$b))
  #   abline(h = 0); abline(v = 0)
  return(p.motif)
}

GetGeneListAndOrder <- function(fpath, dat.rhyth.relamp){
  gene.list <- ReadListToVector(fname = fpath, HEADER = FALSE)
  # order by rhythmicity
  dat.rhyth.sub <- subset(dat.rhyth.relamp, gene %in% gene.list)
  
  gene.list.ordered <- dat.rhyth.sub %>%
    subset(., tissue %in% tiss) %>%
    arrange(desc(amp))
  return(gene.list.ordered$gene)
}

GetExprsMat <- function(dat.rhyth.relamp, tfs, dhs.tiss){
  X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
  rownames(X.exprs) <- X.exprs$gene
  X.exprs$gene <- NULL
  return(X.exprs)
}

GetMotifMat <- function(N, jgene){
  N.sub <- subset(N, gene == jgene)
  X.motif <- dcast(data = N.sub, formula = motif ~ tissue, value.var = "motevo.value")
  rownames(X.motif) <- X.motif$motif
  X.motif$motif <- NULL
  # replace NA with 0
  X.motif[is.na(X.motif)] <- 0
  return(X.motif)
}