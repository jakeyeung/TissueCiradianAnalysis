#library(penalizedLDA)
#library(wordcloud)
#library(dplyr)

SetUpMatForLda <- function(mat.fg, mat.bg, mat.bg2 = NULL, has.peaks=TRUE){
  if (is.null(mat.bg2)){
    mat.fgbg <- bind_rows(mat.fg, mat.bg)
  } else {
    mat.fgbg <- bind_rows(mat.fg, mat.bg, mat.bg2)
  }
  mat.fgbg[is.na(mat.fgbg)] <- 0
  if (has.peaks){
    rownames(mat.fgbg) <- paste(mat.fgbg$peak, mat.fgbg$gene, sep = ";"); mat.fgbg$peak <- NULL; mat.fgbg$gene <- NULL
  } else {
    rownames(mat.fgbg) <- make.unique(as.character(mat.fgbg$gene)); mat.fgbg$gene <- NULL
  }
  if (is.null(mat.bg2)){
    labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)))
  } else {
    labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)), rep(3, nrow(mat.bg2)))
  }
  return(list(mat.fgbg = mat.fgbg, labels = labels))
}

CrossProductCnames <- function(cnames, jsep = ";"){
#   cnames.new <- c()
#   for (cnamei in cnames){
#     for (cnamej in cnames){
#       if (cnamei == cnamej) next
#       cname.cross <- sort(c(cnamei, cnamej))
#       cname.cross <- paste(cname.cross, collapse = ";")
#       if (! cname.cross %in% cnames.new){
#         cnames.new <- c(cnames.new, cname.cross)
#       }
#     }
#   }
  cnames.new <- unlist(lapply(cnames, function(cname){
    paste(cname, cnames, sep = jsep)
  }))
  return(cnames.new)
}

SortCname <- function(cname, jsep = ";"){
  # make ZDF;ABC -> ABC;ZDF
  paste(sort(strsplit(cname, jsep)[[1]]), collapse = jsep)
}

IsDouble <- function(cname, jsep = ";"){
  # return TRUE if ABC;ABC
  strsplit(cname, jsep)[[1]][[1]] == strsplit(cname, jsep)[[1]][[2]]
}

DuplicateCnames <- function(cnames){
  # After running CrossProductCnames, check which ones are duplicates and remove them
  dupes.i <- c()
  sorted.cnames <- c()
  for (i in seq(length(cnames))){
    cname <- SortCname(cnames[i])
    if (IsDouble(cname)){
      dupes.i <- c(dupes.i, i)
      next
    }
    if (cname %in% sorted.cnames){
      dupes.i <- c(dupes.i, i)
      next
    }
    sorted.cnames <- c(sorted.cnames, cname)
  }
  return(dupes.i)
}

CrossProduct <- function(mat, remove.duplicates = TRUE){
  # Make cross product of each column with every other column
  if (is.data.frame(mat)){
    convert.to.mat <- TRUE
  } else {
    convert.to.mat <- FALSE
  }
  jmat <- apply(mat, MARGIN = 2, function(jcol){
    return(sweep(mat, 1, jcol, FUN = "*"))
  })
  if (convert.to.mat){
    # add cross colnames
    new.cnames <- CrossProductCnames(colnames(mat))
    mat <- bind_cols(jmat)
    colnames(mat) <- new.cnames
    if (remove.duplicates){
      mat[DuplicateCnames(colnames(mat))] <- list(NULL)  
    }
  }
  return(mat)
}

CrossProductTwoSets <- function(mat1, mat2){
  # cross product mat1 with mat2 and you're done
  if (all(is.data.frame(mat1), is.data.frame(mat2))) convert.to.mat <- TRUE
  
  jmat <- apply(mat1, MARGIN = 2, function(jcol){
    return(sweep(mat2, 1, jcol, FUN = "*"))
  })
  # add cnames
  new.cnames <- unlist(lapply(colnames(mat1), function(c1){
    paste(c1, colnames(mat2), sep = ";")
  }))
  jmat <- bind_cols(jmat)
  colnames(jmat) <- new.cnames
  return(jmat)
}

DoLdaPromoters <- function(fg.mat, bg.mat){
  # remove columns with 0 within-class SD
  fg.mat <- fg.mat[, apply(fg.mat, 2, sd) > 0]
  bg.mat <- bg.mat[, apply(bg.mat, 2, sd) > 0]
  #   flat.mat <- flat.mat[, apply(flat.mat, 2, sd) > 0]
  #   rhyth.mat <- rhyth.mat[, apply(flat.mat, 2, sd) > 0]
  
  fg.labs <- rep(1, nrow(fg.mat))
  bg.labs <- rep(2, nrow(bg.mat))
  
  common.cols <- intersect(colnames(fg.mat), colnames(bg.mat))
  
  fg.mat <- fg.mat[, common.cols]
  bg.mat <- bg.mat[, common.cols]
  
  merged.mat <- rbind(fg.mat, bg.mat)
  merged.labs <- c(fg.labs, bg.labs)
  
  out <- PenalizedLDA(x = merged.mat, y = merged.labs, lambda = 0.11, K = 1)
  return(out)
}

LongToMat <- function(N.long, jvar = "sitecount"){
  warning("Should use LongToMat.lda")
  mat <- dcast(N.long, formula = gene.uniq ~ motif, value.var = jvar)
  rownames(mat) <- mat$gene.uniq; mat$gene.uniq <- NULL
  return(mat)
}

LongToMat.lda <- function(N.long, jvar = "sitecount"){
  mat <- dcast(N.long, formula = gene.uniq ~ motif, value.var = jvar)
  rownames(mat) <- mat$gene.uniq; mat$gene.uniq <- NULL
  return(mat)
}

RunPenalizedLDA <- function(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, bg.tiss, bg.genes.type, sitecount.name, minamp, n.bg.genes = "all", jlambda=0.1){
  # from sitecounts_analysis_dhs_on_gene_body.redo.R
  if (jmodels != ""){
    gene.list <- as.character(subset(fits.best, model %in% jmodels & amp.avg > minamp)$gene)
  } else {
    # ignore amplitude requirement if model is flat
    gene.list <- as.character(subset(fits.best, model %in% jmodels)$gene)
  }
  
  # define flat genes
  if (bg.genes.type == "Flat"){
    gene.list.bg <- as.character(subset(fits.best, model == "")$gene)
    #  gene.list.bg <- unique(as.character(subset(dat.fit, tissue == "Liver" & gene %in% gene.list.bg & int.rnaseq > 9)$gene))
  } else if (bg.genes.type == "Rhythmic"){
    gene.list.bg <- as.character(subset(fits.best, n.rhyth >= 8 & amp.avg > minamp)$gene)
  } else if (bg.genes.type == "same"){
    gene.list.bg <- gene.list  # use same genes
  } else{
    warning("Unexpected bg.genes.type")
  }
  
  set.seed(0)
  if (n.bg.genes == "all"){
    # take all bg genes
    # print("No subsampling for bg genes")
    
  } else if (n.bg.genes == "same"){
    n.genes <- length(gene.list)
    gene.list.bg <- sample(gene.list.bg, size = n.genes)
  } else {
    # assume a numeric
    gene.list.bg <- sample(gene.list.bg, size = n.bg.genes)
  }
  
  N.sub <- subset(N.long.sum.bytiss, gene %in% gene.list & tissue %in% fg.tiss)
  N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% bg.tiss)
  # N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% fg.tiss)
  
  # sample to have same number as N.sub
  (n.frgrd <- nrow(N.sub))
  gene.list.N <- unique(N.sub$gene)
  
  # Run LDA -----------------------------------------------------------------
  
  # make a new label gene;tissue which will be my rowname
  N.sub$genetiss <- paste(N.sub$gene, N.sub$tissue, sep = ";")
  N.bkgrd.all$genetiss <- paste(N.bkgrd.all$gene, N.bkgrd.all$tissue, sep = ";")
  
  N.sub.mat <- dcast(data = N.sub, formula = genetiss ~ motif, fill = 0, value.var = sitecount.name)
  N.bkgrd.all.mat <- dcast(N.bkgrd.all, genetiss ~ motif, fill = 0, value.var = sitecount.name)
  
  common.cnames <- intersect(colnames(N.sub.mat), colnames(N.bkgrd.all.mat))
  # ignore non-common cnames
  N.sub.mat <- N.sub.mat[, common.cnames]
  N.bkgrd.all.mat <- N.bkgrd.all.mat[, common.cnames]
  
  labs <- c(rep(1, nrow(N.sub.mat)), rep(2, nrow(N.bkgrd.all.mat)))
  
  N.mat.merged <- rbind(N.sub.mat, N.bkgrd.all.mat)
  rownames(N.mat.merged) <- N.mat.merged$genetiss; N.mat.merged$genetiss <- NULL
  
  out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = jlambda, K = 1, standardized = FALSE)  
  
#   out.df <- data.frame(xproj = out$xproj, lab = out$y)
#   boxplot(xproj ~ lab, data = out.df)
  
#   discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
#   cnames <- colnames(N.mat.merged)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
#   plot(seq(length(discrim.filt)), discrim.filt, 
#        main=paste(paste(fg.tiss, collapse = ","), "vs", paste(bg.tiss, collapse=","), sitecount.name, "\n", 
#                   "ngenes fg:", length(gene.list), "ngenes bg:", length(gene.list.bg), "\n",
#                   "models", paste0(jmodels, collapse=","), "\n",
#                   "bg genes type:", bg.genes.type))
#   text(seq(length(discrim.filt)), discrim.filt, labels = cnames)
  out$fg.genes <- gene.list
  return(out)
}

PenalizedLdaLong <- function(fits.best, N.long, jmodel, jmodel.bg, jlambda, K, jvalue.var = "sitecount"){
  # fits.best: from Nconds
  # N.long: sitecounts from promoters
  # jmodel: get genes belonging to jmodel (can be vector) acts as foreground
  # jmodel.bg: background model for comparing against jmodel
  
  genes.fg <- subset(fits.best, model %in% jmodel)$gene
  genes.bg <- subset(fits.best, model %in% jmodel.bg)$gene
  
  N.fg <- subset(N.long, gene %in% genes.fg)
  N.bg <- subset(N.long, gene %in% genes.bg)
  
  N.mat.fg <- dcast(data = subset(N.long, gene %in% genes.fg), 
                    formula = gene.uniq ~ motif,
                    value.var = jvalue.var,
                    fill = 0)
  N.mat.bg <- dcast(data = subset(N.long, gene %in% genes.bg), 
                    formula = gene.uniq ~ motif,
                    value.var = jvalue.var,
                    fill = 0)
  
  N.mat.fgbg <- rbind(N.mat.fg, N.mat.bg)
  # add labels starting from 1...  to correspond to fg or bg
  labs <- c(rep(1, nrow(N.mat.fg)), rep(2, nrow(N.mat.bg)))
  
  # remove first column name move it to rownames
  rownames(N.mat.fgbg) <- N.mat.fgbg$gene.uniq; N.mat.fgbg$gene.uniq <- NULL
  
  # remove column sums == 0
  N.mat.fgbg[[names(which(colSums(N.mat.fgbg) == 0))]] <- NULL

  lda.out <- PenalizedLDA(N.mat.fgbg, labs, lambda = jlambda, K = K)
  return(lda.out)
}

BoxplotLdaOut <- function(out, jtitle = "Title", jdim = 1, axis.names=FALSE, jlab="Projection Space", horizontal=FALSE){
  out.df <- data.frame(xproj = out$xproj[, jdim], lab = out$y)
  if (axis.names[[1]] == FALSE){
    boxplot(xproj ~ lab, data = out.df, main = jtitle, names = axis.names, ylab=jlab, horizontal = horizontal)
  } else {
    if (horizontal == FALSE){
      boxplot(xproj ~ lab, data = out.df, main = jtitle, col = "lightgray", xaxt = "n",  xlab = "", ylab=jlab, horizontal = horizontal)
      axis(1, labels = FALSE)
      text(x =  seq_along(axis.names), y = par("usr")[3] - 1, srt = 45, adj = 1,
           labels = axis.names, xpd = TRUE)
    } else {
      boxplot(xproj ~ lab, data = out.df, main = jtitle, col = "lightgray", ylab = "", xlab=jlab, horizontal = horizontal, names = axis.names)
    }
    
  }
}

ScatterLdaOut2D <- function(out, jtitle = "Title"){
  out.df <- data.frame(xproj1 = out$xproj[, 1], xproj2 = out$xproj[, 2], lab = as.factor(out$y))
#   return(out.df)
  ggplot(out.df, aes(x = xproj1, y = xproj2, colour = lab)) + geom_point(alpha = 0.3)
}

PlotLdaOut <- function(out, jtitle = "Title", jcex = 1, take.n = NA, from.bottom = NA, jdim = NULL){
  library(wordcloud)
  try(detach("package:gplots", unload=TRUE), silent=TRUE)
  if (is.null(jdim)){
    discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
    cnames <- colnames(out$x)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
  } else {
    discrim.filt <- sort(out$discrim[, jdim][which(out$discrim[, jdim] != 0)], decreasing = FALSE)
    cnames <- colnames(out$x)[which(out$discrim[, jdim] != 0)][order(out$discrim[, jdim][which(out$discrim[, jdim] != 0)], decreasing = FALSE)]
  }
  # optionally filter if there are too many to plot
  if (is.na(take.n) & (is.na(from.bottom))){
    wordcloud::textplot(x = seq(length(discrim.filt)), y = discrim.filt, words = cnames, main = jtitle, xlab = "Index", ylab = "Loadings", cex = jcex)
  } else {
    if (from.bottom){
      # take from bottom (head)
      discrim.filt <- head(discrim.filt, take.n)
      cnames <- head(cnames, take.n)
    } else {
      discrim.filt <- tail(discrim.filt, take.n)
      cnames <- tail(cnames, take.n)
    }
    wordcloud::textplot(x = seq(length(discrim.filt)), y = discrim.filt, words = cnames, main = jtitle, xlab = "Index", ylab = "Loadings", cex = jcex)
  }
}

PlotLdaOut2D <- function(out, jdim = 1, jtitle = "Title", jcex = 1, jxlab="Discrim 1", jylab="Discrim 2"){
  discrim.vec1 <- out$discrim[, 1]
  discrim.vec2 <- out$discrim[, 2]
  discrim.vec <- cbind(out$discrim[, 1], out$discrim[, 2])
  # remove rows where BOTH are zero
  rows.keep <- which(apply(abs(discrim.vec), 1, max) > 0)
  discrim.filt <- discrim.vec[rows.keep, ]
  cnames <- colnames(out$x)[rows.keep]
  # optionally filter if there are too many to plot
  wordcloud::textplot(x = discrim.filt[, 1], y = discrim.filt[, 2], words = cnames, main = jtitle, xlab = jxlab, ylab = jylab, cex = jcex)
  #   plot(x = discrim.filt[, 1], y = discrim.filt[, 2],main = jtitle, xlab = "Discrim 1", ylab = "Discrim 2", cex = jcex)
}

SortLda <- function(out, jdim = NULL){
  if (is.null(jdim)){
    cnames <- colnames(out$x)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
  } else {
    cnames <- colnames(out$x)[which(out$discrim[, jdim] != 0)][order(out$discrim[, jdim][which(out$discrim[, jdim] != 0)], decreasing = FALSE)]
  }
}

PlotSeparation <- function(out, jtitle){
  boxplot(list("Foreground" = out$xproj[out$y == 1], "Background" = out$xproj[out$y == 2]), main = jtitle, xlab = "Group 1: foreground. Group 2: background", ylab = "Projection")
}