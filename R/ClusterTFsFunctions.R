AvgExprsAcrossTissues <- function(df){
  # take average exprs from df which contains exprs colname 
  # df.avg.tiss <- ddply(df, .(tissue), AvgExprs)
  df.avg.tiss <- ddply(df, .(tissue), summarise,
                       mean.exprs = mean(exprs))
  return(df.avg.tiss)
}

MaxExprsAcrossTissues <- function(df){
  # df.sub <- subset(df, experiment="rnaseq")
  df.max.tiss <- ddply(df, .(tissue), summarise,
                       max.exprs = max(exprs))
  return(df.max.tiss)
}

AllOnes <- function(v){
  # check if vector is all 1s
  return(all(v == 1))
}

FilterAllOnesOrZeros <- function(mat){
  # Filter out a binary matrix only for rows that do not contain
  # all 1s or 0s. They are not meaningful for differentiating TF
  # across tissues
  mat <- mat[which(apply(mat, 1, AllOnes) == FALSE), ]  # remove TFs with only 1s 
  mat <- mat[which(rowSums(mat) > 0), ] # remove TFs with only 0s
  return(mat)
}

Binarize <- function(mat, cutoff){
  mat[which(mat < cutoff)] <- 0
  mat[which(mat >= cutoff)] <- 1
  return(mat)
}

MaxExprsAcrossTissuesDplyr <- function(dat){
  df.max.tiss <- dat %>%
    group_by(tissue) %>%
    summarise(max.exprs = max(exprs))
}